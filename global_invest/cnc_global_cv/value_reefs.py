"""Value Reefs."""
import argparse
import glob
import logging
import os
import math
import sys

import ecoshard
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import pygeoprocessing
import numpy
import shutil
import taskgraph
import taskgraph_downloader_pnn
import tempfile

ECOSHARD_BUCKET_URL = (
    r'https://storage.googleapis.com/critical-natural-capital-ecoshards/')

GLOBAL_REEFS_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'ipbes-cv_reef_md5_5a90d55a505813b5aa9662faee351bf8.tif')

LS_POPULATION_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'lspop2017_md5_faaad64d15d0857894566199f62d422c.zip')
POVERTY_POPULATION_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'Poverty_Count_2017_clean_md5_0cc6e0187be07e760e66f759a0a1f7e8.tif')
GLOBAL_DEM_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'global_dem_md5_22c5c09ac4c4c722c844ab331b34996c.tif')

M_PER_DEGREE = 111300.0
REEF_PROT_DIST = 2000.0

# set a 1GB limit for the cache
gdal.SetCacheMax(2**30)

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(pathname)s.%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)

WORKSPACE_DIR = 'value_reefs_workspace'
ECOSHARD_DIR = os.path.join(WORKSPACE_DIR, 'ecoshard')
CHURN_DIR = os.path.join(WORKSPACE_DIR, 'churn')


def clean_convolve_2d(
        signal_raster_band_path, kernel_raster_band_path, target_raster_path,
        working_dir=None):
    """Do 2D convolution but mask out any close to 0 values to 0."""
    temp_workspace_dir = tempfile.mkdtemp(
        dir=os.path.dirname(signal_raster_band_path[0]),
        prefix='clean_convolve_2d')
    temp_target_raster = os.path.join(temp_workspace_dir, 'result.tif')
    pygeoprocessing.convolve_2d(
        signal_raster_band_path, kernel_raster_band_path,
        temp_target_raster, working_dir=working_dir)
    nodata = pygeoprocessing.get_raster_info(temp_target_raster)['nodata'][0]
    pygeoprocessing.raster_calculator(
        [(temp_target_raster, 1), (1e-6, 'raw')], set_almost_zero_to_zero,
        target_raster_path, gdal.GDT_Float32, nodata)
    try:
        shutil.rmtree(temp_workspace_dir)
    except OSError:
        LOGGER.exception('unable to remove %s' % temp_workspace_dir)


def calculate_utm_srs(lng, lat):
    """Calculate UTM SRS from the lng/lat point given.

    Parameters:
        lng (float): longitude point.
        lat (float): latitude point.

    Returns:
        osr.SpatialReference in the UTM zone that contains the point (lng, lat)

    """
    utm_code = (math.floor((lng+180)/6) % 60) + 1
    lat_code = 6 if lat > 0 else 7
    epsg_code = int('32%d%02d' % (lat_code, utm_code))
    utm_srs = osr.SpatialReference()
    utm_srs.ImportFromEPSG(epsg_code)
    return utm_srs


def set_almost_zero_to_zero(array, eps):
    result = numpy.empty_like(array)
    result[:] = numpy.where(numpy.abs(array) < eps, 0.0, array)
    return result


def create_averaging_kernel_raster(
        radius_in_pixels, kernel_filepath, normalize=True):
    """Create a flat raster kernel with a 2d radius given.

    Parameters:
        radius_in_pixels (tuple): the (x/y) distance of the averaging kernel.
        kernel_filepath (string): The path to the file on disk where this
            kernel should be stored.  If this file exists, it will be
            overwritten.

    Returns:
        None

    """
    driver = gdal.GetDriverByName('GTiff')
    LOGGER.debug(radius_in_pixels)
    kernel_raster = driver.Create(
        kernel_filepath, int(2*radius_in_pixels[0]),
        int(2*radius_in_pixels[1]), 1, gdal.GDT_Float32)

    # Make some kind of geotransform, it doesn't matter what but
    # will make GIS libraries behave better if it's all defined
    kernel_raster.SetGeoTransform([1, 0.1, 0, 1, 0, -0.1])
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    kernel_raster.SetProjection(srs.ExportToWkt())

    kernel_band = kernel_raster.GetRasterBand(1)
    kernel_band.SetNoDataValue(-9999)

    n_cols = kernel_raster.RasterXSize
    n_rows = kernel_raster.RasterYSize
    iv, jv = numpy.meshgrid(range(n_rows), range(n_cols), indexing='ij')

    cx = n_cols / 2.0
    cy = n_rows / 2.0

    kernel_array = numpy.where(
        ((cx-jv)**2 + (cy-iv)**2)**0.5 <= radius_in_pixels[0], 1.0, 0.0)
    kernel_array = numpy.where(
        ((cx-jv) / radius_in_pixels[0])**2 +
        ((cy-iv) / radius_in_pixels[1])**2 <= 1.0, 1.0, 0.0)
    LOGGER.debug(kernel_array)

    # normalize
    if normalize:
        kernel_array /= numpy.sum(kernel_array)
    kernel_band.WriteArray(kernel_array)
    kernel_band = None
    kernel_raster = None


def mask_by_height_op(pop_array, dem_array, mask_height, pop_nodata):
    """Set pop to 0 if > height."""
    result = numpy.zeros(shape=pop_array.shape)
    valid_mask = (
        (dem_array < mask_height) & ~numpy.isclose(pop_array, pop_nodata))
    result[valid_mask] = pop_array[valid_mask]
    return result


def align_raster_list(
        raster_path_list, target_directory, target_sr_wkt=None,
        align_index=0):
    """Aligns all the raster paths.

    Rasters are aligned using the pixel size of the first raster and use
    the intersection and near interpolation methods.

    Parameters:
        raster_path_list (list): list of str paths to rasters.
        target_directory (str): path to a directory to hold the aligned
            rasters.
        target_sr_wkt (str): if not None this is the projection.
        align_index (int): this index is used to indicate which raster should
            be used for aligned pixel size.

    Returns:
        list of raster paths that are aligned with intersection and near
            interpolation algorithm.

    """
    LOGGER.debug('aligning %s', raster_path_list)
    aligned_path_list = [
        os.path.join(target_directory, os.path.basename(path))
        for path in raster_path_list]
    target_pixel_size = pygeoprocessing.get_raster_info(
        raster_path_list[align_index])['pixel_size']
    LOGGER.debug('about to align: %s', str(raster_path_list))
    pygeoprocessing.align_and_resize_raster_stack(
        raster_path_list, aligned_path_list,
        ['near'] * len(raster_path_list), target_pixel_size,
        'intersection', target_sr_wkt=target_sr_wkt)
    return aligned_path_list


def calculate_reef_population_value(
        shore_sample_point_vector_path, dem_raster_path,
        reef_habitat_raster_path, population_raster_path_id_target_list,
        temp_workspace_dir):
    """Calculate population within protective range of reefs.

    Parameters:
        shore_sample_point_vector_path (str): path to a point shapefile that is
            used for referencing the points of interest on the coastline where.
        dem_raster_path (str): path to a dem used to mask population by height
            in wgs84 lat/lng projection.
        reef_habitat_raster_path (str): path to a mask raster where reef
            habitat exists.
        population_raster_path_id_list (list): list of
            (raster_path, field_id, target_path)
            tuples. The values in the raster paths will be masked where it
            overlaps with < 10m dem height and convolved within 2km. That
            result is in turn spread onto the habitat coverage at a distance
            of the protective distance of reefs. These rasters are in
            wgs84 lat/lng projection.

    Returns:
        None

    """
    wgs84_srs = osr.SpatialReference()
    wgs84_srs.ImportFromEPSG(4326)
    aligned_pop_raster_list = align_raster_list(
        [x[0] for x in population_raster_path_id_target_list] +
        [reef_habitat_raster_path, dem_raster_path], temp_workspace_dir,
        wgs84_srs.ExportToWkt(),
        align_index=len(population_raster_path_id_target_list))

    raster_info = pygeoprocessing.get_raster_info(
        aligned_pop_raster_list[0])
    target_pixel_size = raster_info['pixel_size']

    n_pixels_in_prot_dist = max(1, int(REEF_PROT_DIST / (
        M_PER_DEGREE * abs(target_pixel_size[0]))))
    kernel_radius = [n_pixels_in_prot_dist, n_pixels_in_prot_dist]
    kernel_filepath = os.path.join(
        temp_workspace_dir, 'reef_kernel.tif')
    create_averaging_kernel_raster(
        kernel_radius, kernel_filepath, normalize=False)

    buffered_point_raster_mask_path = os.path.join(
            temp_workspace_dir, 'reef_buffer_mask.tif')
    make_buffered_point_raster_mask(
        shore_sample_point_vector_path,
        aligned_pop_raster_list[0],
        temp_workspace_dir, 'reefs_all',
        REEF_PROT_DIST, buffered_point_raster_mask_path)

    for pop_index, (_, pop_id, target_path) in enumerate(
            population_raster_path_id_target_list):
        # mask to < 10m
        pop_height_masked_path = os.path.join(
            temp_workspace_dir, '%s_masked_by_10m.tif' % pop_id)

        pygeoprocessing.raster_calculator(
            [(aligned_pop_raster_list[pop_index], 1),
             (aligned_pop_raster_list[-1], 1),
             (10.0, 'raw'),  # mask to 10 meters
             (raster_info['nodata'][0], 'raw')],  # the -1 index is the dem
            mask_by_height_op, pop_height_masked_path, gdal.GDT_Float32,
            raster_info['nodata'][0])

        # spread the < 10m population out 2km
        n_pixels_in_2km = int(2000.0 / (
            M_PER_DEGREE * abs(target_pixel_size[0])))
        kernel_radius_2km = [n_pixels_in_2km, n_pixels_in_2km]
        kernel_2km_filepath = os.path.join(
            temp_workspace_dir, '2km_kernel_%s.tif' % pop_id)
        create_averaging_kernel_raster(
            kernel_radius_2km, kernel_2km_filepath, normalize=False)
        pop_sum_within_2km_path = os.path.join(
            temp_workspace_dir, 'pop_sum_within_2km_%s.tif' % pop_id)
        pygeoprocessing.convolve_2d(
            (pop_height_masked_path, 1), (kernel_2km_filepath, 1),
            pop_sum_within_2km_path)

        align_reef_habitat_raster_path = aligned_pop_raster_list[-2]
        population_hab_spread_raster_path = os.path.join(
            temp_workspace_dir, 'reef_%s_spread.tif' % (pop_id))
        clean_convolve_2d(
            (pop_sum_within_2km_path, 1), (kernel_filepath, 1),
            population_hab_spread_raster_path)
        hab_raster_info = pygeoprocessing.get_raster_info(
            align_reef_habitat_raster_path)

        # warp pop result to overlay
        clipped_pop_hab_spread_raster_path = os.path.join(
            temp_workspace_dir, 'reef_%s_spread_clipped.tif' % pop_id)
        pygeoprocessing.warp_raster(
            population_hab_spread_raster_path,
            hab_raster_info['pixel_size'],
            clipped_pop_hab_spread_raster_path,
            'near')

        hab_spread_nodata = pygeoprocessing.get_raster_info(
            clipped_pop_hab_spread_raster_path)['nodata'][0]
        hab_nodata = hab_raster_info['nodata'][0]

        pygeoprocessing.raster_calculator(
            [(clipped_pop_hab_spread_raster_path, 1),
             (align_reef_habitat_raster_path, 1),
             (buffered_point_raster_mask_path, 1),
             (hab_spread_nodata, 'raw'),
             (hab_nodata, 'raw')],
            intersect_and_mask_raster_op, target_path,
            gdal.GDT_Float32, hab_spread_nodata),

        ecoshard.build_overviews(target_path)


def calculate_reef_value(
        shore_sample_point_vector, template_raster_path,
        reef_habitat_raster_path, working_dir, target_reef_value_raster_path):
    """Calculate habitat value.

    Will create rasters in the `working_dir` directory named from the
    `habitat_fieldname_list` values containing relative importance of
    global habitat. The higher the value of a pixel the more important that
    pixel of habitat is for protection of the coastline.

    Parameters:
        shore_sample_point_vector (str): path to CV analysis vector containing
            at least the fields `Rt` and `Rt_nohab_[hab]` for all habitat
            types under consideration.
        template_raster_path (str): path to an existing raster whose size and
            shape will be used to be the base of the raster that's created
            for each habitat type.
        habitat_fieldname_list (list): list of habitat ids to analyise.
        habitat_vector_path_map (dict): maps fieldnames from
            `habitat_fieldname_list` to 3-tuples of
            (path to hab vector (str), risk val (float),
             protective distance (float)).
        working_dir (str): path to directory containing habitat back projection
            results
        target_reef_value_raster_path (str): path to raster value raster for
            that habitat.

    Returns:
        None.

    """
    temp_workspace_dir = os.path.join(working_dir, 'hab_value_churn')
    for dir_path in [working_dir, temp_workspace_dir]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    gpkg_driver = ogr.GetDriverByName('gpkg')
    shore_sample_point_vector = gdal.OpenEx(
        shore_sample_point_vector, gdal.OF_VECTOR)
    shore_sample_point_layer = shore_sample_point_vector.GetLayer()

    reef_service_id = 'Rt_habservice_reefs_all'

    buffer_habitat_path = os.path.join(
        temp_workspace_dir, 'reefs_all_buffer.gpkg')
    buffer_habitat_vector = gpkg_driver.CreateDataSource(
        buffer_habitat_path)
    wgs84_srs = osr.SpatialReference()
    wgs84_srs.ImportFromEPSG(4326)
    buffer_habitat_layer = (
        buffer_habitat_vector.CreateLayer(
            reef_service_id, wgs84_srs, ogr.wkbPolygon))
    buffer_habitat_layer.CreateField(ogr.FieldDefn(
        reef_service_id, ogr.OFTReal))
    buffer_habitat_layer_defn = buffer_habitat_layer.GetLayerDefn()

    shore_sample_point_layer.ResetReading()
    buffer_habitat_layer.StartTransaction()
    for point_index, point_feature in enumerate(shore_sample_point_layer):
        if point_index % 1000 == 0:
            LOGGER.debug(
                'point buffering is %.2f%% complete',
                point_index / shore_sample_point_layer.GetFeatureCount() *
                100.0)
        # for each point, convert to local UTM to buffer out a given
        # distance then back to wgs84
        point_geom = point_feature.GetGeometryRef()
        if point_geom.GetX() > 178 or point_geom.GetX() < -178:
            continue
        utm_srs = calculate_utm_srs(point_geom.GetX(), point_geom.GetY())
        wgs84_to_utm_transform = osr.CoordinateTransformation(
            wgs84_srs, utm_srs)
        utm_to_wgs84_transform = osr.CoordinateTransformation(
            utm_srs, wgs84_srs)
        point_geom.Transform(wgs84_to_utm_transform)
        buffer_poly_geom = point_geom.Buffer(REEF_PROT_DIST)
        buffer_poly_geom.Transform(utm_to_wgs84_transform)

        buffer_point_feature = ogr.Feature(buffer_habitat_layer_defn)
        buffer_point_feature.SetGeometry(buffer_poly_geom)
        buffer_point_feature.SetField(
            reef_service_id,
            point_feature.GetField('Rt_habservice_reefs_all'))
        buffer_habitat_layer.CreateFeature(buffer_point_feature)
        buffer_point_feature = None
        point_feature = None
        buffer_poly_geom = None
        point_geom = None

    # at this point every shore point has been buffered to the effective
    # habitat distance and the habitat service has been saved with it
    buffer_habitat_layer.CommitTransaction()
    buffer_habitat_layer = None
    buffer_habitat_vector = None
    value_coverage_raster_path = os.path.join(
        temp_workspace_dir, 'reefs_all_value_cover.tif')
    pygeoprocessing.new_raster_from_base(
        template_raster_path, value_coverage_raster_path,
        gdal.GDT_Float32, [0],
        raster_driver_creation_tuple=(
            'GTIFF', (
                'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=LZW',
                'BLOCKXSIZE=256', 'BLOCKYSIZE=256', 'SPARSE_OK=TRUE')))
    pygeoprocessing.rasterize(
        buffer_habitat_path, value_coverage_raster_path,
        option_list=[
            'ATTRIBUTE=%s' % reef_service_id,
            'MERGE_ALG=ADD'])

    value_coverage_nodata = pygeoprocessing.get_raster_info(
        value_coverage_raster_path)['nodata'][0]
    hab_nodata = pygeoprocessing.get_raster_info(
        reef_habitat_raster_path)['nodata'][0]

    wgs84_srs = osr.SpatialReference()
    wgs84_srs.ImportFromEPSG(4326)

    aligned_value_hab_raster_path_list = align_raster_list(
        [value_coverage_raster_path, reef_habitat_raster_path],
        temp_workspace_dir, target_sr_wkt=wgs84_srs.ExportToWkt())

    pygeoprocessing.raster_calculator(
        [(aligned_value_hab_raster_path_list[0], 1),
         (aligned_value_hab_raster_path_list[1], 1),
         (value_coverage_nodata, 'raw'), (hab_nodata, 'raw')],
        intersect_raster_op, target_reef_value_raster_path, gdal.GDT_Float32,
        value_coverage_nodata)

    ecoshard.build_overviews(target_reef_value_raster_path)


def intersect_and_mask_raster_op(
        array_a, array_b, mask_array, nodata_a, nodata_b):
    result = numpy.empty_like(array_a)
    result[:] = nodata_a
    valid_mask = (
        ~numpy.isclose(array_a, nodata_a) &
        ~numpy.isclose(array_b, nodata_b) &
        (mask_array == 1))
    result[valid_mask] = array_a[valid_mask]
    return result


def intersect_raster_op(array_a, array_b, nodata_a, nodata_b):
    """Only return values from a where a and b are defined."""
    result = numpy.empty_like(array_a)
    result[:] = nodata_a
    valid_mask = (
        ~numpy.isclose(array_a, nodata_a) &
        ~numpy.isclose(array_b, nodata_b))
    result[valid_mask] = array_a[valid_mask]
    return result


def make_buffered_point_raster_mask(
        shore_sample_point_vector_path, template_raster_path, workspace_dir,
        habitat_id,
        protective_distance, target_buffer_raster_path):
    gpkg_driver = ogr.GetDriverByName('GPKG')
    buffer_habitat_path = os.path.join(
        workspace_dir, '%s_buffer.gpkg' % habitat_id)
    buffer_habitat_vector = gpkg_driver.CreateDataSource(
        buffer_habitat_path)
    wgs84_srs = osr.SpatialReference()
    wgs84_srs.ImportFromEPSG(4326)
    buffer_habitat_layer = (
        buffer_habitat_vector.CreateLayer(
            habitat_id, wgs84_srs, ogr.wkbPolygon))
    buffer_habitat_layer_defn = buffer_habitat_layer.GetLayerDefn()

    shore_sample_point_vector = gdal.OpenEx(
        shore_sample_point_vector_path, gdal.OF_VECTOR)
    shore_sample_point_layer = shore_sample_point_vector.GetLayer()
    buffer_habitat_layer.StartTransaction()
    for point_index, point_feature in enumerate(shore_sample_point_layer):
        if point_index % 1000 == 0:
            LOGGER.debug(
                'point buffering is %.2f%% complete',
                point_index / shore_sample_point_layer.GetFeatureCount() *
                100.0)
        # for each point, convert to local UTM to buffer out a given
        # distance then back to wgs84
        point_geom = point_feature.GetGeometryRef()
        utm_srs = calculate_utm_srs(point_geom.GetX(), point_geom.GetY())
        wgs84_to_utm_transform = osr.CoordinateTransformation(
            wgs84_srs, utm_srs)
        utm_to_wgs84_transform = osr.CoordinateTransformation(
            utm_srs, wgs84_srs)
        point_geom.Transform(wgs84_to_utm_transform)
        buffer_poly_geom = point_geom.Buffer(protective_distance)
        buffer_poly_geom.Transform(utm_to_wgs84_transform)

        buffer_point_feature = ogr.Feature(buffer_habitat_layer_defn)
        buffer_point_feature.SetGeometry(buffer_poly_geom)
        buffer_habitat_layer.CreateFeature(buffer_point_feature)
        buffer_point_feature = None
        point_feature = None
        buffer_poly_geom = None
        point_geom = None

    # at this point every shore point has been buffered to the effective
    # habitat distance and the habitat service has been saved with it
    buffer_habitat_layer.CommitTransaction()
    buffer_habitat_layer = None
    buffer_habitat_vector = None
    pygeoprocessing.new_raster_from_base(
        template_raster_path, target_buffer_raster_path,
        gdal.GDT_Float32, [0],
        raster_driver_creation_tuple=(
            'GTIFF', (
                'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=LZW',
                'BLOCKXSIZE=256', 'BLOCKYSIZE=256', 'SPARSE_OK=TRUE')))
    pygeoprocessing.rasterize(
        buffer_habitat_path, target_buffer_raster_path, burn_values=[1])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calcualte risk from reefs')
    parser.add_argument(
        'cv_risk_vector_pattern', nargs='+',
        help='Can be a pattern to a file.')
    args = parser.parse_args()

    for dir_path in [WORKSPACE_DIR, ECOSHARD_DIR, CHURN_DIR]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    task_graph = taskgraph.TaskGraph(WORKSPACE_DIR, 2, 5.0)
    tdd_downloader = taskgraph_downloader_pnn.TaskGraphDownloader(
        ECOSHARD_DIR, task_graph)

    tdd_downloader.download_ecoshard(
        GLOBAL_REEFS_RASTER_URL, 'reefs')
    tdd_downloader.download_ecoshard(
        LS_POPULATION_RASTER_URL, 'total_pop', decompress='unzip',
        local_path='lspop2017')
    tdd_downloader.download_ecoshard(
        POVERTY_POPULATION_RASTER_URL, 'poor_pop')
    tdd_downloader.download_ecoshard(
        GLOBAL_DEM_RASTER_URL, 'global_dem')

    reef_degree_pixel_size = [0.004, -0.004]
    wgs84_srs = osr.SpatialReference()
    wgs84_srs.ImportFromEPSG(4326)
    projected_reef_raster_path = os.path.join(CHURN_DIR, 'wgs84_reefs.tif')
    project_reef_task = task_graph.add_task(
        func=pygeoprocessing.warp_raster,
        args=(
            tdd_downloader.get_path('reefs'), reef_degree_pixel_size,
            projected_reef_raster_path, 'near'),
        kwargs={'target_sr_wkt': wgs84_srs.ExportToWkt()},
        target_path_list=[projected_reef_raster_path],
        task_name='project reefs to wgs84')

    for cv_risk_vector_pattern in args.cv_risk_vector_pattern:
        for cv_vector_path in glob.glob(cv_risk_vector_pattern):
            basename = os.path.basename(os.path.splitext(cv_vector_path)[0])
            LOGGER.debug(basename)

            reefs_value_raster_path = os.path.join(
                WORKSPACE_DIR, "reefs_value_%s.tif" % basename)
            reefs_poor_pop_coverage_raster_path = os.path.join(
                WORKSPACE_DIR, "reefs_poor_pop_coverage_%s.tif" % basename)
            reefs_total_pop_coverage_raster_path = os.path.join(
                WORKSPACE_DIR, "reefs_total_pop_coverage_%s.tif" % basename)

            # task_graph.add_task(
            #     func=calculate_reef_value,
            #     args=(
            #         cv_vector_path, tdd_downloader.get_path('reefs'),
            #         projected_reef_raster_path, WORKSPACE_DIR,
            #         reefs_value_raster_path),
            #     target_path_list=[reefs_value_raster_path],
            #     dependent_task_list=[project_reef_task],
            #     task_name='calculate reef value for %s' % basename)

            task_graph.add_task(
                func=calculate_reef_population_value,
                args=(
                    cv_vector_path, tdd_downloader.get_path('global_dem'),
                    projected_reef_raster_path,
                    [(tdd_downloader.get_path('total_pop'), 'total_pop',
                      reefs_poor_pop_coverage_raster_path),
                     (tdd_downloader.get_path('poor_pop'), 'poor_pop',
                      reefs_total_pop_coverage_raster_path)], WORKSPACE_DIR),
                dependent_task_list=[project_reef_task],
                target_path_list=[
                    reefs_poor_pop_coverage_raster_path,
                    reefs_total_pop_coverage_raster_path],
                task_name='calculate reef population value')
            task_graph.join()

    task_graph.join()
    task_graph.close()
