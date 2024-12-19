"""Global CV Analysis for CNC.

Design doc is here:

https://docs.google.com/document/d/18AcJM-rXeIYgEsmqlaUwdtm7gdWLiaD6kkRkpARILlw/edit#heading=h.bbujb61ete53

"""
import argparse
import bisect
import collections
import datetime
import gzip
import logging
import math
import multiprocessing
import os
import shutil
import sys
import tempfile
import threading
import zipfile

import ecoshard
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import numpy
import pygeoprocessing
import retrying
import rtree
import shapely.geometry
import shapely.strtree
import shapely.wkt
import taskgraph

gdal.SetCacheMax(2**27)

WORKSPACE_DIR = 'global_cv_workspace'
CHURN_DIR = os.path.join(WORKSPACE_DIR, 'churn')
ECOSHARD_DIR = os.path.join(WORKSPACE_DIR, 'ecoshard')
GRID_WORKSPACE_DIR = os.path.join(WORKSPACE_DIR, 'grid_workspaces')
ECOSHARD_DIR = os.path.join(WORKSPACE_DIR, 'ecoshard')
HABITAT_VALUE_DIR = os.path.join(WORKSPACE_DIR, 'habitat_value')
TARGET_NODATA = -1
TARGET_CV_VECTOR_PATH = os.path.join(
    WORKSPACE_DIR, 'global_cv_analysis_result.gpkg')

# landcover inputs come from this file
LANDCOVER_RASTER_DATA_FILE = 'GLOBAL_CV_DATA_INPUTS.txt'

# [minx, miny, maxx, maxy].
GLOBAL_AOI_WGS84_BB = [-179, -65, 180, 77]
RELIEF_SAMPLE_DISTANCE = 5000.0
N_FETCH_RAYS = 16  # this is hardcoded because of WWIII fields
MAX_FETCH_DISTANCE = 60000
M_PER_DEGREE = 111300.0

ECOSHARD_BUCKET_URL = (
    r'https://storage.googleapis.com/critical-natural-capital-ecoshards/')

EMPTY_RASTER_URL = (
    ECOSHARD_BUCKET_URL + 'empty_md5_f8f71e20668060bda7567ca33149a45c.tif')

GLOBAL_POLYGON_URL = (
    ECOSHARD_BUCKET_URL +
    'ipbes-cv_global_polygon_simplified_geometries_'
    'md5_653118dde775057e24de52542b01eaee.gpkg')

BUFFER_VECTOR_URL = (
    ECOSHARD_BUCKET_URL +
    'buffered_global_shore_5km_md5_a68e1049c1c03673add014cd29b7b368.gpkg')
SHORE_GRID_URL = (
    ECOSHARD_BUCKET_URL +
    'shore_grid_md5_07aea173cf373474c096f1d5e3463c2f.gpkg')

GLOBAL_WWIII_GZ_URL = (
    ECOSHARD_BUCKET_URL +
    'wave_watch_iii_md5_c8bb1ce4739e0a27ee608303c217ab5b.gpkg.gz')
GLOBAL_DEM_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'global_dem_md5_22c5c09ac4c4c722c844ab331b34996c.tif')
LS_POPULATION_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'lspop2017_md5_faaad64d15d0857894566199f62d422c.zip')
POVERTY_POPULATION_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'Poverty_Count_nans_cleaned_md5_c3d4e9443997889f2706e9600e72c975.tif')
SLR_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'MSL_Map_MERGED_Global_AVISO_NoGIA_Adjust_'
    'md5_3072845759841d0b2523d00fe9518fee.tif')
GLOBAL_GEOMORPHOLOGY_VECTOR_URL = (
    ECOSHARD_BUCKET_URL +
    'geomorphology_md5_e65eff55840e7a80cfcb11fdad2d02d7.gpkg')
# The reefs raster is not in wgs84 but the script below projects it so
GLOBAL_REEFS_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'ipbes-cv_reef_wgs84_compressed_md5_96d95cc4f2c5348394eccff9e8b84e6b.tif')
GLOBAL_MANGROVES_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'ipbes-cv_mangrove_md5_0ec85cb51dab3c9ec3215783268111cc.tif')
GLOBAL_SEAGRASS_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'ipbes-cv_seagrass_md5_a9cc6d922d2e74a14f74b4107c94a0d6.tif')
GLOBAL_SALTMARSH_RASTER_URL = (
    ECOSHARD_BUCKET_URL +
    'ipbes-cv_saltmarsh_md5_203d8600fd4b6df91f53f66f2a011bcd.tif')

GLOBAL_MESOAMERICAN_BARRIER_REEF = (
    ECOSHARD_BUCKET_URL +
    'mesoamerican_barrier_reef_md5_4e8964bf9bf3d0f5fdddebd204eb3b32.tif')
GLOBAL_NEW_CALEDONIAN_BARRIER_REEF = (
    ECOSHARD_BUCKET_URL +
    'new_caledonian_barrier_reef_md5_ab1e8e54304a916101b9a53baa47c908.tif')
GLOBAL_GREAT_BARRIER_REEF = (
    ECOSHARD_BUCKET_URL +
    'great_barrier_reef_md5_68ec7a55480f0d8a6445ccfde21b48e0.tif')
GLOBAL_KEYS_BARRIER_REEF = (
    ECOSHARD_BUCKET_URL +
    'keys_barrier_reef_md5_21a5374ffd7237c6ccf629a8d78b3b51.tif')

GLOBAL_DATA_URL_MAP = {
    'geomorphology': GLOBAL_GEOMORPHOLOGY_VECTOR_URL,
    'mangroves_forest': GLOBAL_MANGROVES_RASTER_URL,
    'reefs': GLOBAL_REEFS_RASTER_URL,
    'seagrass': GLOBAL_SEAGRASS_RASTER_URL,
    'saltmarsh_wetland': GLOBAL_SALTMARSH_RASTER_URL,
    'dem': GLOBAL_DEM_RASTER_URL,
    'slr': SLR_RASTER_URL,
    'landmass': GLOBAL_POLYGON_URL,
    'shore_grid': SHORE_GRID_URL,
    'poverty_count': POVERTY_POPULATION_RASTER_URL,
    'mesoamerican_barrier_reef': GLOBAL_MESOAMERICAN_BARRIER_REEF,
    'new_caledonian_barrier_reef': GLOBAL_NEW_CALEDONIAN_BARRIER_REEF,
    'great_barrier_reef': GLOBAL_GREAT_BARRIER_REEF,
    'keys_barrier_reef': GLOBAL_KEYS_BARRIER_REEF,
    'global_wwiii_vector_path': GLOBAL_WWIII_GZ_URL,
    }

HAB_FIELDS = [
    '4_500',
    '2_2000',
    'mangroves_forest',
    'saltmarsh_wetland',
    'seagrass',
    'reefs_all',
]

REEF_FIELDS = [
    'reefs',
    'mesoamerican_barrier_reef',
    'new_caledonian_barrier_reef',
    'great_barrier_reef',
    'keys_barrier_reef',
]

# this is used to evaluate habitat value
FINAL_HAB_FIELDS = REEF_FIELDS + [
    'mangroves_forest',
    'saltmarsh_wetland',
    'seagrass',
    '4_500',
    '2_2000',
]


SEDTYPE_TO_RISK = {
    0: 5,  # unknown
    1: 5,  # sandy
    2: 1,  # unerodable
    3: 5,  # muddy
    4: 1,  # coral/mangrove
}


# This dictionary maps landcode id to (risk, dist) tuples
LULC_CODE_TO_HAB_MAP = {
    0: (0, None),
    10: (0, None),
    11: (0, None),
    12: (0, None),
    20: (0, None),
    30: (0, None),
    40: (4, 500),
    50: (1, 2000),
    60: (1, 2000),
    61: (1, 2000),
    62: (1, 2000),
    70: (1, 2000),
    71: (1, 2000),
    72: (1, 2000),
    80: (1, 2000),
    81: (1, 2000),
    82: (1, 2000),
    90: (1, 2000),
    100: (1, 2000),
    110: (2, 2000),
    120: (2, 2000),
    121: (2, 2000),
    122: (2, 2000),
    130: (2, 2000),
    140: (2, 2000),
    150: (4, 500),
    151: (4, 500),
    152: (4, 500),
    153: (4, 500),
    160: (2, 1000),
    170: (2, 1000),
    180: (2, 1000),
    190: (0, None),
    200: (0, None),
    201: (0, None),
    202: (0, None),
    210: (0, None),
    220: (0, None),
    }

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(pathname)s.%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)

STOP_SENTINEL = 'STOP'

HABITAT_VECTOR_PATH_MAP = {
    'reefs': (
        os.path.join(
            ECOSHARD_DIR, os.path.basename(GLOBAL_DATA_URL_MAP['reefs'])),
        1, 2000.0),
    'mangroves_forest': (
        os.path.join(
            ECOSHARD_DIR, os.path.basename(
                GLOBAL_DATA_URL_MAP['mangroves_forest'])),
        1, 2000.0),
    'saltmarsh_wetland': (
        os.path.join(
            ECOSHARD_DIR, os.path.basename(
                GLOBAL_DATA_URL_MAP['saltmarsh_wetland'])),
        2, 1000.0),
    'seagrass': (
        os.path.join(
            ECOSHARD_DIR, os.path.basename(GLOBAL_DATA_URL_MAP['seagrass'])),
        4, 500.0),
    'mesoamerican_barrier_reef': (
        os.path.join(
            ECOSHARD_DIR, os.path.basename(
                GLOBAL_DATA_URL_MAP['mesoamerican_barrier_reef'])),
        1, 35000.0),
    'new_caledonian_barrier_reef': (
        os.path.join(
            ECOSHARD_DIR, os.path.basename(
                GLOBAL_DATA_URL_MAP['new_caledonian_barrier_reef'])),
        1, 20000.0),
    'great_barrier_reef': (
        os.path.join(
            ECOSHARD_DIR, os.path.basename(
                GLOBAL_DATA_URL_MAP['great_barrier_reef'])),
        1, 35000.0),
    'keys_barrier_reef': (
        os.path.join(
            ECOSHARD_DIR, os.path.basename(
                GLOBAL_DATA_URL_MAP['keys_barrier_reef'])),
        1, 25000.0),
    }


def download_and_unzip(url, target_dir, target_token_path):
    """Download `url` to `target_dir` and touch `target_token_path`."""
    zipfile_path = os.path.join(target_dir, os.path.basename(url))
    LOGGER.debug('url %s, zipfile_path: %s', url, zipfile_path)
    ecoshard.download_url(url, zipfile_path)

    with zipfile.ZipFile(zipfile_path, 'r') as zip_ref:
        zip_ref.extractall(target_dir)

    with open(target_token_path, 'w') as touchfile:
        touchfile.write(f'unzipped {zipfile_path}')


def download_and_ungzip(url, target_path, buffer_size=2**20):
    """Download `url` to `target_dir` and touch `target_token_path`."""
    gzipfile_path = os.path.join(
        os.path.dirname(target_path), os.path.basename(url))
    ecoshard.download_url(url, gzipfile_path)

    with gzip.open(gzipfile_path, 'rb') as gzip_file:
        with open(target_path, 'wb') as target_file:
            while True:
                content = gzip_file.read(buffer_size)
                if content:
                    target_file.write(content)
                else:
                    break


def build_rtree(vector_path):
    """Build an rtree that can be queried for nearest neighbors.

    Parameters:
        vector_path (str): path to vector of geometry to build into
            r tree.

    Returns:
        rtree.Index object that will return shapely geometry objects with
            a field_val_map field that contains the 'fieldname'->value pairs
            from the original vector. The main object will also have a
            `field_name_type_list` field which contains original
            fieldname/field type pairs

    """
    geometry_prep_list = []
    vector = gdal.OpenEx(vector_path, gdal.OF_VECTOR)
    layer = vector.GetLayer()
    layer_defn = layer.GetLayerDefn()
    field_name_type_list = []
    for index in range(layer_defn.GetFieldCount()):
        field_name = layer_defn.GetFieldDefn(index).GetName()
        field_type = layer_defn.GetFieldDefn(index).GetType()
        field_name_type_list.append((field_name, field_type))

    LOGGER.debug('loop through features for rtree')
    for index, feature in enumerate(layer):
        feature_geom = feature.GetGeometryRef().Clone()
        feature_geom_shapely = shapely.wkb.loads(feature_geom.ExportToWkb())
        field_val_map = {}
        for field_name, _ in field_name_type_list:
            field_val_map[field_name] = (
                feature.GetField(field_name))
        geometry_prep_list.append(
            (index, feature_geom_shapely.bounds, field_val_map))
    LOGGER.debug('constructing the tree')
    r_tree = rtree.index.Index(geometry_prep_list)
    LOGGER.debug('all done')
    r_tree.field_name_type_list = field_name_type_list
    return r_tree


def build_strtree(vector_path):
    """Build an rtree that generates geom and preped geometry.

    Parameters:
        vector_path (str): path to vector of geometry to build into
            r tree.

    Returns:
        strtree.STRtree object that will return shapely geometry objects
            with a .prep field that is prepared geomtry for fast testing,
            a .geom field that is the base gdal geometry, and a field_val_map
            field that contains the 'fieldname'->value pairs from the original
            vector. The main object will also have a `field_name_type_list`
            field which contains original fieldname/field type pairs

    """
    geometry_prep_list = []
    vector = gdal.OpenEx(vector_path, gdal.OF_VECTOR)
    layer = vector.GetLayer()
    layer_defn = layer.GetLayerDefn()
    field_name_type_list = []
    for index in range(layer_defn.GetFieldCount()):
        field_name = layer_defn.GetFieldDefn(index).GetName()
        field_type = layer_defn.GetFieldDefn(index).GetType()
        field_name_type_list.append((field_name, field_type))

    LOGGER.debug('loop through features for rtree')
    for index, feature in enumerate(layer):
        feature_geom = feature.GetGeometryRef().Clone()
        feature_geom_shapely = shapely.wkb.loads(feature_geom.ExportToWkb())
        feature_geom_shapely.prep = shapely.prepared.prep(feature_geom_shapely)
        feature_geom_shapely.geom = feature_geom
        feature_geom_shapely.id = index
        feature_geom_shapely.field_val_map = {}
        for field_name, _ in field_name_type_list:
            feature_geom_shapely.field_val_map[field_name] = (
                feature.GetField(field_name))
        geometry_prep_list.append(feature_geom_shapely)
    LOGGER.debug('constructing the tree')
    r_tree = shapely.strtree.STRtree(geometry_prep_list)
    LOGGER.debug('all done')
    r_tree.field_name_type_list = field_name_type_list
    return r_tree


def cv_grid_worker(
        bb_work_queue,
        cv_point_complete_queue,
        global_landmass_vector_path,
        geomorphology_vector_path,
        slr_raster_path,
        global_dem_raster_path,
        wwiii_vector_path,
        habitat_raster_path_map,
        ):
    """Worker process to calculate CV for a grid.

    Parameters:
        bb_work_queue (multiprocessing.Queue): contains
            [minx, miny, maxx, maxy] bounding box values to be processed or
            `STOP_SENTINEL` values to indicate the worker should be terminated.
        cv_point_complete_queue (multiprocessing.Queue): this queue is used to
            pass completed CV point vectors for further processing. It will be
            terminated by `STOP_SENTINEL`.
        global_landmass_vector_path (str): path to global landmass vector. Used
            for intersecting rays to determine shoreline exposure.
        geomorphology_vector_path (str): path to line geometry geomorphology
            layer that has a field called 'SEDTYPE'
        global_dem_raster_path (str): path to a global dem/bathymetry raster.
        slr_raster_path (str): path to a sea level rise raster.
        wwiii_vector_path (str): path to wave watch III dataset that has
            fields REI_PCT[degree] and REI_V[degree] for degrees 0-360 in 22.5
            degree increments.
        habitat_raster_path_map (dict): mapt a habitat id to a
            (raster path, risk, dist(m)) tuple. These are the raster versions
            of habitats to use in Rhab.

    Returns:
        None.

    """
    LOGGER.info('build geomorphology rtree')
    geomorphology_strtree = build_strtree(geomorphology_vector_path)
    geomorphology_proj_wkt = pygeoprocessing.get_vector_info(
        geomorphology_vector_path)['projection']
    gegeomorphology_proj = osr.SpatialReference()
    gegeomorphology_proj.ImportFromWkt(geomorphology_proj_wkt)

    LOGGER.info('build landmass rtree')
    landmass_strtree = build_strtree(global_landmass_vector_path)

    LOGGER.info('build wwiii rtree')
    wwiii_rtree = build_rtree(wwiii_vector_path)

    target_pixel_size = [
        SHORE_POINT_SAMPLE_DISTANCE / 4,
        -SHORE_POINT_SAMPLE_DISTANCE / 4]

    while True:
        payload = bb_work_queue.get()
        if payload == STOP_SENTINEL:
            LOGGER.debug('stopping')
            # put it back so others can stop
            bb_work_queue.put(STOP_SENTINEL)
            break
        else:
            LOGGER.debug('running')
        # otherwise payload is the bounding box
        index, (lng_min, lat_min, lng_max, lat_max) = payload
        bounding_box_list = [lng_min, lat_min, lng_max, lat_max]
        buffered_bounding_box_list = [
            lng_min-0.1, lat_min-0.1, lng_max+0.1, lat_max+0.1]
        # create workspace
        workspace_dir = os.path.join(
            GRID_WORKSPACE_DIR, '%d_%s_%s_%s_%s' % (
                index, lng_min, lat_min, lng_max, lat_max))

        try:
            os.makedirs(workspace_dir)
        except OSError:
            pass

        # task_graph = taskgraph.TaskGraph(workspace_dir, -1)

        utm_srs = calculate_utm_srs((lng_min+lng_max)/2, (lat_min+lat_max)/2)
        wgs84_srs = osr.SpatialReference()
        wgs84_srs.ImportFromEPSG(4326)

        try:
            local_geomorphology_vector_path = os.path.join(
                workspace_dir, 'geomorphology.gpkg')
            clip_geometry(
                bounding_box_list, wgs84_srs, utm_srs,
                ogr.wkbMultiLineString, geomorphology_strtree,
                local_geomorphology_vector_path)

            shore_point_vector_path = os.path.join(
                workspace_dir, 'shore_points.gpkg')
            sample_line_to_points(
                local_geomorphology_vector_path, shore_point_vector_path,
                SHORE_POINT_SAMPLE_DISTANCE)

            local_landmass_vector_path = os.path.join(
                workspace_dir, 'landmass.gpkg')
            clip_geometry(
                buffered_bounding_box_list, wgs84_srs, utm_srs,
                ogr.wkbPolygon, landmass_strtree,
                local_landmass_vector_path)

            landmass_boundary_vector_path = os.path.join(
                workspace_dir, 'landmass_boundary.gpkg')
            vector_to_lines(
                local_landmass_vector_path, landmass_boundary_vector_path)

            local_dem_path = os.path.join(
                workspace_dir, 'dem.tif')
            clip_and_reproject_raster(
                global_dem_raster_path, local_dem_path, utm_srs.ExportToWkt(),
                bounding_box_list, RELIEF_SAMPLE_DISTANCE, 'bilinear', True,
                target_pixel_size)

            local_slr_path = os.path.join(
                workspace_dir, 'slr.tif')
            clip_and_reproject_raster(
                global_dem_raster_path, local_slr_path, utm_srs.ExportToWkt(),
                bounding_box_list, 0, 'bilinear', True,
                target_pixel_size)

            # Rrelief
            LOGGER.info('calculate relief on %s', workspace_dir)
            calculate_relief(
                shore_point_vector_path, local_dem_path, 'relief')
            LOGGER.info('calculate rhab on %s', workspace_dir)
            # Rhab
            calculate_rhab(
                shore_point_vector_path, habitat_raster_path_map, 'Rhab',
                target_pixel_size)

            # Rslr
            calculate_slr(shore_point_vector_path, local_slr_path, 'slr')

            # wind and wave power
            calculate_wind_and_wave(
                shore_point_vector_path, local_landmass_vector_path,
                landmass_boundary_vector_path,
                local_dem_path, wwiii_rtree, 'rei', 'ew')

            # Rsurge
            calculate_surge(shore_point_vector_path, local_dem_path, 'surge')

            # Rgeomorphology
            calculate_geomorphology(
                shore_point_vector_path, local_geomorphology_vector_path,
                'Rgeomorphology')

            LOGGER.info('completed %s', shore_point_vector_path)
            cv_point_complete_queue.put(shore_point_vector_path)

        except Exception:
            LOGGER.exception('error on %s, removing workspace', payload)
            retrying_rmtree(workspace_dir)


def make_shore_kernel(kernel_path):
    """Make a 3x3 raster with a 9 in the middle and 1s on the outside."""
    driver = gdal.GetDriverByName('GTiff')
    kernel_raster = driver.Create(
        kernel_path.encode('utf-8'), 3, 3, 1,
        gdal.GDT_Byte)

    # Make some kind of geotransform, it doesn't matter what but
    # will make GIS libraries behave better if it's all defined
    kernel_raster.SetGeoTransform([0, 1, 0, 0, 0, -1])
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS('WGS84')
    kernel_raster.SetProjection(srs.ExportToWkt())

    kernel_band = kernel_raster.GetRasterBand(1)
    kernel_band.SetNoDataValue(127)
    kernel_band.WriteArray(numpy.array([[1, 1, 1], [1, 9, 1], [1, 1, 1]]))


def calculate_geomorphology(
        shore_point_vector_path, geomorphology_vector_path,
        geomorphology_fieldname):
    """Sample the geomorphology vector path for the closest line to each point.

    Parameters:
        shore_point_vector_path (str):  path to a point shapefile to
            for relief point analysis.
        geomorphology_vector_path (str): path to a vector of lines that
            contains the integer field 'SEDTYPE'.
        geomorphology_fieldname (str): fieldname to add to
            `shore_point_vector_path`.

    Returns:
        None.

    """
    shore_point_vector = gdal.OpenEx(
        shore_point_vector_path, gdal.OF_VECTOR | gdal.GA_Update)
    shore_point_layer = shore_point_vector.GetLayer()
    shore_point_layer.CreateField(
        ogr.FieldDefn(geomorphology_fieldname, ogr.OFTReal))
    geomorphology_strtree = build_strtree(geomorphology_vector_path)

    shore_point_layer.StartTransaction()
    for shore_point_feature in shore_point_layer:
        shore_point_geom = shapely.wkb.loads(
            shore_point_feature.GetGeometryRef().ExportToWkb())
        min_dist = MAX_FETCH_DISTANCE
        geo_risk = 5
        for line in geomorphology_strtree.query(shore_point_geom.buffer(500)):
            cur_dist = line.distance(shore_point_geom)
            if cur_dist < min_dist:
                min_dist = cur_dist
                geo_risk = line.field_val_map['Rgeo']
        shore_point_feature.SetField(geomorphology_fieldname, geo_risk)
        shore_point_layer.SetFeature(shore_point_feature)
    shore_point_layer.CommitTransaction()
    shore_point_layer = None
    shore_point_vector = None


def calculate_surge(
        shore_point_vector_path, bathymetry_raster_path, surge_fieldname):
    """Calculate surge potential as distance to continental shelf (-150m).

    Parameters:
        base_shore_point_vector_path (string):  path to a point shapefile to
            for relief point analysis.
        global_dem_path (string): path to a DEM raster projected in wgs84.
        surge_fieldname (str): fieldname to add to `shore_point_vector_path`
        workspace_dir (string): path to a directory to make local calculations
            in
        target_surge_point_vector_path (string): path to output vector.
            after completion will a value for closest distance to continental
            shelf called 'surge'.

    Returns:
        None.

    """
    shore_point_vector = gdal.OpenEx(
        shore_point_vector_path, gdal.OF_VECTOR | gdal.GA_Update)
    shore_point_layer = shore_point_vector.GetLayer()
    shore_point_layer.CreateField(ogr.FieldDefn(surge_fieldname, ogr.OFTReal))

    shelf_nodata = 2

    bathymetry_nodata = pygeoprocessing.get_raster_info(
        bathymetry_raster_path)['nodata'][0]

    def mask_shelf(depth_array):
        valid_mask = ~numpy.isclose(depth_array, bathymetry_nodata)
        result_array = numpy.empty(
            depth_array.shape, dtype=numpy.int16)
        result_array[:] = shelf_nodata
        result_array[valid_mask] = 0
        result_array[depth_array < -150] = 1
        return result_array

    workspace_dir = os.path.dirname(shore_point_vector_path)

    shelf_mask_path = os.path.join(workspace_dir, 'shelf_mask.tif')
    pygeoprocessing.raster_calculator(
        [(bathymetry_raster_path, 1)], mask_shelf,
        shelf_mask_path, gdal.GDT_Byte, shelf_nodata)

    # convolve to find edges
    # grid shoreline from raster
    shelf_kernel_path = os.path.join(workspace_dir, 'shelf_kernel.tif')
    shelf_convoultion_raster_path = os.path.join(
        workspace_dir, 'shelf_convolution.tif')
    make_shore_kernel(shelf_kernel_path)
    pygeoprocessing.convolve_2d(
        (shelf_mask_path, 1), (shelf_kernel_path, 1),
        shelf_convoultion_raster_path, target_datatype=gdal.GDT_Byte,
        target_nodata=255)

    nodata = pygeoprocessing.get_raster_info(
        shelf_convoultion_raster_path)['nodata'][0]

    def _shelf_mask_op(shelf_convolution):
        """Mask values on land that border the continental shelf."""
        result = numpy.empty(shelf_convolution.shape, dtype=numpy.uint8)
        result[:] = nodata
        valid_mask = shelf_convolution != nodata
        # If a pixel is on land, it gets at least a 9, but if it's all on
        # land it gets an 17 (8 neighboring pixels), so we search between 9
        # and 17 to determine a shore pixel
        result[valid_mask] = numpy.where(
            (shelf_convolution[valid_mask] >= 9) &
            (shelf_convolution[valid_mask] < 17), 1, nodata)
        return result

    shelf_edge_raster_path = os.path.join(workspace_dir, 'shelf_edge.tif')
    pygeoprocessing.raster_calculator(
        [(shelf_convoultion_raster_path, 1)], _shelf_mask_op,
        shelf_edge_raster_path, gdal.GDT_Byte, nodata)

    shore_geotransform = pygeoprocessing.get_raster_info(
        shelf_edge_raster_path)['geotransform']

    shelf_rtree = rtree.index.Index()

    for offset_info, data_block in pygeoprocessing.iterblocks(
            (shelf_edge_raster_path, 1)):
        row_indexes, col_indexes = numpy.mgrid[
            offset_info['yoff']:offset_info['yoff']+offset_info['win_ysize'],
            offset_info['xoff']:offset_info['xoff']+offset_info['win_xsize']]
        valid_mask = data_block == 1
        x_coordinates = (
            shore_geotransform[0] +
            shore_geotransform[1] * (col_indexes[valid_mask] + 0.5) +
            shore_geotransform[2] * (row_indexes[valid_mask] + 0.5))
        y_coordinates = (
            shore_geotransform[3] +
            shore_geotransform[4] * (col_indexes[valid_mask] + 0.5) +
            shore_geotransform[5] * (row_indexes[valid_mask] + 0.5))

        for x_coord, y_coord in zip(x_coordinates, y_coordinates):
            shelf_rtree.insert(
                0, [x_coord, y_coord, x_coord, y_coord],
                obj=shapely.geometry.Point(x_coord, y_coord))

    shore_point_layer.StartTransaction()
    for point_feature in shore_point_layer:
        point_geometry = point_feature.GetGeometryRef()
        point_shapely = shapely.wkb.loads(point_geometry.ExportToWkb())
        nearest_point = list(shelf_rtree.nearest(
                (point_geometry.GetX(),
                 point_geometry.GetY(),
                 point_geometry.GetX(),
                 point_geometry.GetY()),
                objects='raw', num_results=1))
        if len(nearest_point) > 0:
            distance = nearest_point[0].distance(point_shapely)
            point_feature.SetField(surge_fieldname, float(distance))
        else:
            # so far away it's essentially not an issue
            point_feature.SetField(surge_fieldname, MAX_FETCH_DISTANCE)
        shore_point_layer.SetFeature(point_feature)

    shore_point_layer.CommitTransaction()
    shore_point_layer.SyncToDisk()
    shore_point_layer = None
    shore_point_vector = None


def calculate_wind_and_wave(
        shore_point_vector_path, landmass_vector_path,
        landmass_boundary_vector_path, bathymetry_raster_path,
        wwiii_rtree, wind_fieldname, wave_fieldname):
    """Calculate wind exposure for given points.

    Parameters:
        shore_point_vector_path (str): path to a point vector, this value will
            be modified to hold the total wind exposure at this point.
        landmass_vector_path (str): path to a vector indicating landmass that
            will block wind exposure.
        landmass_boundary_vector_path (str): path to a string vector containing
            the perimeter of `landmass_vector_path`.
        bathymetry_raster_path (str): path to a raster indicating bathymetry
            values. (negative is deeper).
        wwiii_rtree (str): path to an r_tree that can find the nearest point
            in lat/lng whose object has values 'REI_PCT', 'REI_V',
            'WavP_[DIR]', 'WavPPCT', 'V10PCT_[DIR]'.
        wind_fieldname (str): fieldname to add to `shore_point_vector_path` for
            wind power.
        wave_fieldname (str): fieldname to add to `shore_point_vector_path` for
            wave power.

    Returns:
        None

    """
    gpkg_driver = ogr.GetDriverByName('gpkg')
    temp_workspace_dir = tempfile.mkdtemp(
        dir=os.path.dirname(shore_point_vector_path),
        prefix='calculate_rwind_')
    temp_fetch_rays_vector = gpkg_driver.CreateDataSource(
        os.path.join(temp_workspace_dir, 'fetch_rays.gpkg'))
    layer_name = 'fetch_rays'
    shore_point_projection_wkt = pygeoprocessing.get_vector_info(
        shore_point_vector_path)['projection']
    shore_point_srs = osr.SpatialReference()
    shore_point_srs.ImportFromWkt(shore_point_projection_wkt)
    temp_fetch_rays_layer = (
        temp_fetch_rays_vector.CreateLayer(
            str(layer_name), shore_point_srs, ogr.wkbLineString))
    temp_fetch_rays_layer.CreateField(ogr.FieldDefn(
        'fetch_dist', ogr.OFTReal))
    temp_fetch_rays_layer.CreateField(ogr.FieldDefn(
        'direction', ogr.OFTReal))
    temp_fetch_rays_defn = temp_fetch_rays_layer.GetLayerDefn()

    # These WWIII fields are the only ones needed for wind & wave equations
    # Copy them to a new vector which also gets more fields added with
    # computed values.
    target_shore_point_vector = gdal.OpenEx(
        shore_point_vector_path, gdal.OF_VECTOR | gdal.GA_Update)
    target_shore_point_layer = target_shore_point_vector.GetLayer()
    target_shore_point_layer.CreateField(
        ogr.FieldDefn(wind_fieldname, ogr.OFTReal))
    target_shore_point_layer.CreateField(
        ogr.FieldDefn(wave_fieldname, ogr.OFTReal))
    for ray_index in range(N_FETCH_RAYS):
        compass_degree = int(ray_index * 360 / N_FETCH_RAYS)
        target_shore_point_layer.CreateField(
            ogr.FieldDefn('fdist_%d' % compass_degree, ogr.OFTReal))
        target_shore_point_layer.CreateField(
            ogr.FieldDefn('fdepth_%d' % compass_degree, ogr.OFTReal))

    # Iterate over every shore point
    LOGGER.info("Casting rays and extracting bathymetry values")
    bathy_raster = gdal.OpenEx(
        bathymetry_raster_path, gdal.OF_RASTER | gdal.GA_ReadOnly)
    bathy_band = bathy_raster.GetRasterBand(1)
    bathy_raster_info = pygeoprocessing.get_raster_info(bathymetry_raster_path)
    bathy_gt = bathy_raster_info['geotransform']
    bathy_inv_gt = gdal.InvGeoTransform(bathy_gt)

    landmass_vector = gdal.OpenEx(landmass_vector_path, gdal.OF_VECTOR)
    landmass_layer = landmass_vector.GetLayer()
    landmass_geom_list = [
        shapely.wkb.loads(f.GetGeometryRef().ExportToWkb())
        for f in landmass_layer]
    landmass_union_geom = shapely.ops.cascaded_union(landmass_geom_list)
    landmass_layer = None
    landmass_vector = None
    landmass_union_geom_prep = shapely.prepared.prep(landmass_union_geom)

    landmass_boundary_vector = gdal.OpenEx(
        landmass_boundary_vector_path, gdal.OF_VECTOR)
    landmass_boundary_layer = landmass_boundary_vector.GetLayer()
    landmass_boundary_geom_list = [
        shapely.wkb.loads(f.GetGeometryRef().ExportToWkb())
        for f in landmass_boundary_layer]
    landmass_boundary_union_geom = shapely.ops.cascaded_union(
        landmass_boundary_geom_list)
    landmass_boundary_layer = None
    landmass_boundary_vector = None
    landmass_boundary_union_geom_prep = shapely.prepared.prep(
        landmass_boundary_union_geom)
    landmass_boundary_strtree = build_strtree(landmass_boundary_vector_path)

    target_shore_point_layer.StartTransaction()
    temp_fetch_rays_layer.StartTransaction()

    # make a transfomer for local points to lat/lng for wwiii_rtree
    wgs84_srs = osr.SpatialReference()
    wgs84_srs.ImportFromEPSG(4326)
    base_to_target_transform = osr.CoordinateTransformation(
        shore_point_srs, wgs84_srs)

    for shore_point_feature in target_shore_point_layer:
        shore_point_geom = shore_point_feature.GetGeometryRef().Clone()
        _ = shore_point_geom.Transform(base_to_target_transform)
        wwiii_point = next(wwiii_rtree.nearest(
            (shore_point_geom.GetX(), shore_point_geom.GetY()), 1,
            objects='raw'))
        rei_value = 0.0
        height_list = []
        period_list = []
        e_local = 0.0
        e_ocean = 0.0

        shore_point_geometry = shore_point_feature.GetGeometryRef()
        shapely_point = shapely.wkb.loads(
            shore_point_geometry.ExportToWkb())
        if landmass_union_geom_prep.contains(shapely_point):
            new_point = shapely.ops.nearest_points(
                landmass_boundary_union_geom, shapely_point)[0]
            LOGGER.debug(
                'new point: %s %s', str(new_point), str(shore_point_geometry))
            shore_point_geometry = ogr.CreateGeometryFromWkb(new_point.wkb)

        # Iterate over every ray direction
        for sample_index in range(N_FETCH_RAYS):
            compass_degree = int(sample_index * 360 / N_FETCH_RAYS)
            compass_theta = float(sample_index) / N_FETCH_RAYS * 360

            # wwiii_point should be closest point to shore point

            rei_pct = wwiii_point[
                'REI_PCT%d' % int(compass_theta)]
            rei_v = wwiii_point[
                'REI_V%d' % int(compass_theta)]
            cartesian_theta = -(compass_theta - 90)

            # Determine the direction the ray will point
            delta_x = math.cos(cartesian_theta * math.pi / 180)
            delta_y = math.sin(cartesian_theta * math.pi / 180)

            # Start a ray offset from the shore point
            # so that rays start outside of the landmass.
            # Shore points are interpolated onto the coastline,
            # but floating point error results in points being just
            # barely inside/outside the landmass.
            offset = 10

            point_a_x = (
                shore_point_geometry.GetX() + delta_x * offset)
            point_a_y = (
                shore_point_geometry.GetY() + delta_y * offset)

            origin_point = shapely.geometry.Point(point_a_x, point_a_y)
            if landmass_union_geom_prep.intersects(origin_point):
                # the origin is inside the landmass, skip
                continue

            point_b_x = point_a_x + delta_x * (
                MAX_FETCH_DISTANCE)
            point_b_y = point_a_y + delta_y * (
                MAX_FETCH_DISTANCE)

            # build ray geometry so we can intersect it later
            ray_geometry = ogr.Geometry(ogr.wkbLineString)
            ray_geometry.AddPoint(point_a_x, point_a_y)
            ray_geometry.AddPoint(point_b_x, point_b_y)

            # keep a shapely version of the ray so we can do fast intersection
            # with it and the entire landmass
            ray_point_origin_shapely = shapely.geometry.Point(
                point_a_x, point_a_y)

            if not landmass_boundary_union_geom_prep.intersects(
                    ray_point_origin_shapely):
                # the origin is in ocean, so we'll get a ray length > 0.0

                # This algorithm searches for intersections, if one is found
                # the ray updates and a smaller intersection set is determined
                # by experimentation I've found this is significant, but not
                # an order of magnitude, faster than looping through all
                # original possible intersections.  Since this algorithm
                # will be run for a long time, it's worth the additional
                # complexity
                tested_indexes = set()
                while True:
                    intersection = False
                    ray_envelope = ray_geometry.GetEnvelope()
                    for landmass_line in landmass_boundary_strtree.query(
                             shapely.geometry.box(
                                *[ray_envelope[i] for i in [0, 2, 1, 3]])):
                        if landmass_line.id in tested_indexes:
                            continue
                        tested_indexes.add(landmass_line.id)
                        if ray_geometry.Intersects(landmass_line.geom):
                            intersection_point = ray_geometry.Intersection(
                                landmass_line.geom)
                            # offset the dist with smallest_feature_size
                            # update the endpoint of the ray
                            ray_geometry = ogr.Geometry(ogr.wkbLineString)
                            ray_geometry.AddPoint(point_a_x, point_a_y)
                            ray_geometry.AddPoint(
                                intersection_point.GetX(),
                                intersection_point.GetY())
                            intersection = True
                            break
                    if not intersection:
                        break

                ray_step_loc = 0.0
                bathy_values = []
                # walk along ray
                ray_shapely = shapely.wkb.loads(ray_geometry.ExportToWkb())
                while ray_step_loc < ray_shapely.length:
                    sample_point = ray_shapely.interpolate(ray_step_loc)
                    ray_step_loc += SHORE_POINT_SAMPLE_DISTANCE/4
                    pixel_x, pixel_y = [int(x) for x in gdal.ApplyGeoTransform(
                        bathy_inv_gt,
                        sample_point.coords[0][0], sample_point.coords[0][1])]
                    if (pixel_x < 0 or pixel_y < 0 or
                            pixel_x >= bathy_band.XSize or
                            pixel_y >= bathy_band.YSize):
                        continue
                    bathy_values.append(
                        bathy_band.ReadAsArray(
                            pixel_x, pixel_y, 1, 1)[0][0])

                if bathy_values:
                    avg_bathy_value = numpy.mean(bathy_values)
                else:
                    avg_bathy_value = 0.0
                # when we get here, we have the final ray geometry
                ray_feature = ogr.Feature(temp_fetch_rays_defn)
                ray_feature.SetField('fetch_dist', ray_shapely.length)
                ray_feature.SetField('direction', compass_degree)
                ray_feature.SetGeometry(ray_geometry)
                temp_fetch_rays_layer.CreateFeature(ray_feature)
                rei_value += ray_shapely.length * rei_pct * rei_v
                ray_length = ray_geometry.Length()
                ray_feature = None
                ray_geometry = None

                shore_point_feature.SetField(
                    'fdist_%d' % compass_degree, ray_length)
                shore_point_feature.SetField(
                    'fdepth_%d' % compass_degree, float(avg_bathy_value))

                velocity = wwiii_point['V10PCT_%d' % compass_degree]
                occurrence = wwiii_point['REI_PCT%d' % compass_degree]

                height = compute_wave_height(
                    velocity, ray_shapely.length, avg_bathy_value)
                height_list.append(height)
                period = compute_wave_period(
                    velocity, ray_shapely.length, avg_bathy_value)
                period_list.append(period)
                power = 0.5 * float(height)**2 * float(period)  # UG Eq. 8
                e_local += power * occurrence  # UG Eq. 9

                if intersection:
                    e_ocean += (
                        wwiii_point['WavP_%d' % compass_degree] *
                        wwiii_point['WavPPCT%d' % compass_degree])

                ray_feature = None
                ray_geometry = None
                rei_value += ray_length * rei_pct * rei_v
        shore_point_feature.SetField(wind_fieldname, rei_value)
        shore_point_feature.SetField(wave_fieldname, max(e_ocean, e_local))
        target_shore_point_layer.SetFeature(shore_point_feature)
        shore_point_geometry = None
    target_shore_point_layer.CommitTransaction()
    target_shore_point_layer.SyncToDisk()
    target_shore_point_layer = None
    target_shore_point_vector = None
    temp_fetch_rays_layer.CommitTransaction()
    temp_fetch_rays_layer.SyncToDisk()
    temp_fetch_rays_layer = None
    temp_fetch_rays_vector = None
    bathy_raster = None
    bathy_band = None
    try:
        shutil.rmtree(temp_workspace_dir)
    except Exception:
        LOGGER.exception('unable to remove %s', temp_workspace_dir)


def calculate_slr(shore_point_vector_path, slr_raster_path, target_fieldname):
    """Sample sea level rise raster and store values in shore points.

    Parameters:
        shore_point_vector_path (str): path to a vector of points in a local
            projected coordinate system. This vector will be modified by this
            function to include a new field called `target_fieldname`
            containing the weighted Rhab risk for the given point.
        slr_raster_path (str): path to a sea level rise raster indicating
            sea level rise amout in m.
        target_fieldname (str): fieldname to add to `shore_point_vector_path`
            that will contain the value of sea level rise for that point.

    Returns:
        None.

    """
    try:
        shore_point_vector = gdal.OpenEx(
            shore_point_vector_path, gdal.OF_VECTOR | gdal.GA_Update)
        shore_point_layer = shore_point_vector.GetLayer()
        slr_field = ogr.FieldDefn(target_fieldname, ogr.OFTReal)
        slr_field.SetPrecision(5)
        slr_field.SetWidth(24)
        shore_point_layer.CreateField(slr_field)
        slr_info = pygeoprocessing.get_raster_info(slr_raster_path)
        inv_gt = gdal.InvGeoTransform(slr_info['geotransform'])

        slr_raster = gdal.OpenEx(slr_raster_path, gdal.OF_RASTER)
        slr_band = slr_raster.GetRasterBand(1)

        shore_point_layer.ResetReading()
        for point_feature in shore_point_layer:
            point_geometry = point_feature.GetGeometryRef()
            point_x, point_y = point_geometry.GetX(), point_geometry.GetY()
            point_geometry = None

            pixel_x, pixel_y = [
                int(x) for x in
                gdal.ApplyGeoTransform(inv_gt, point_x, point_y)]
            if pixel_x < 0:
                pixel_x = 0
            if pixel_y < 0:
                pixel_y = 0
            if pixel_x >= slr_band.XSize:
                pixel_x = slr_band.XSize-1
            if pixel_y >= slr_band.YSize:
                pixel_y = slr_band.YSize-1
            try:
                pixel_value = slr_band.ReadAsArray(
                    xoff=pixel_x, yoff=pixel_y, win_xsize=1,
                    win_ysize=1)[0, 0]
            except Exception:
                LOGGER.exception(
                    'slr_band size %d %d', slr_band.XSize,
                    slr_band.YSize)
                raise
            point_feature.SetField(target_fieldname, float(pixel_value))
            shore_point_layer.SetFeature(point_feature)

        shore_point_layer.SyncToDisk()
        shore_point_layer = None
        shore_point_vector = None
        slr_raster = None
        slr_band = None

    except Exception:
        LOGGER.exception('error in slr calc')
        raise


def calculate_rhab(
        shore_point_vector_path, habitat_raster_path_map, target_fieldname,
        target_pixel_size):
    """Add Rhab risk to the shore point vector path.

    Parameters:
        shore_point_vector_path (str): path to a vector of points in a local
            projected coordinate system. This vector will be modified by this
            function to include a new field called `target_fieldname`
            containing the weighted Rhab risk for the given point.
        habitat_raster_path_map (dict): a dictionary mapping "hab id"s to
            (path to raster, risk, effective distance) tuples.
        target_fieldname (str): fieldname to add to `shore_point_vector_path`
            that will contain the value of Rhab calculated for that point.
        target_pixel_size (list): x/y size of clipped habitat in projected
            units

    Returns:
        None.

    """
    shore_point_vector = gdal.OpenEx(
        shore_point_vector_path, gdal.OF_VECTOR | gdal.GA_Update)
    shore_point_layer = shore_point_vector.GetLayer()

    for hab_id in habitat_raster_path_map:
        relief_field = ogr.FieldDefn(hab_id, ogr.OFTReal)
        relief_field.SetPrecision(5)
        relief_field.SetWidth(24)
        shore_point_layer.CreateField(relief_field)

    shore_point_info = pygeoprocessing.get_vector_info(shore_point_vector_path)

    shore_point_feature_risk_map = collections.defaultdict(list)

    tmp_working_dir = tempfile.mkdtemp(
        prefix='calculate_rhab_',
        dir=os.path.dirname(shore_point_vector_path))
    LOGGER.debug(tmp_working_dir)
    for hab_id, (hab_raster_path, risk_val, eff_dist) in (
                habitat_raster_path_map.items()):
        local_hab_raster_path = os.path.join(
            tmp_working_dir, '%s.tif' % str(hab_id))
        LOGGER.debug(
            'clip %s to %s', hab_raster_path, shore_point_info['bounding_box'])
        clip_and_reproject_raster(
            hab_raster_path, local_hab_raster_path,
            shore_point_info['projection'],
            shore_point_info['bounding_box'], eff_dist, 'near', False,
            target_pixel_size)

        # make a convolution kernel as wide as the distance but adapted to
        # the non-square size of the rasters
        kernel_filepath = '%s_kernel%s' % os.path.splitext(
            local_hab_raster_path)
        kernel_radius = [abs(eff_dist / x) for x in target_pixel_size]
        create_averaging_kernel_raster(
            kernel_radius, kernel_filepath, normalize=True)
        hab_effective_area_raster_path = (
            '%s_effective_hab%s' % os.path.splitext(local_hab_raster_path))
        pygeoprocessing.convolve_2d(
            (local_hab_raster_path, 1), (kernel_filepath, 1),
            hab_effective_area_raster_path, mask_nodata=False)
        gt = pygeoprocessing.get_raster_info(
            hab_effective_area_raster_path)['geotransform']
        inv_gt = gdal.InvGeoTransform(gt)

        hab_effective_raster = gdal.OpenEx(
            hab_effective_area_raster_path, gdal.OF_RASTER)
        hab_effective_band = hab_effective_raster.GetRasterBand(1)
        shore_point_layer.ResetReading()
        shore_point_layer.StartTransaction()
        for shore_feature in shore_point_layer:
            shore_geom = shore_feature.GetGeometryRef()
            pixel_x, pixel_y = [
                int(x) for x in
                gdal.ApplyGeoTransform(
                    inv_gt, shore_geom.GetX(), shore_geom.GetY())]
            if pixel_x < 0:
                pixel_x = 0
            if pixel_y < 0:
                pixel_y = 0
            if pixel_x >= hab_effective_band.XSize:
                pixel_x = hab_effective_band.XSize-1
            if pixel_y >= hab_effective_band.YSize:
                pixel_y = hab_effective_band.YSize-1

            try:
                pixel_val = hab_effective_band.ReadAsArray(
                    xoff=pixel_x, yoff=pixel_y, win_xsize=1,
                    win_ysize=1)[0, 0]
            except Exception:
                LOGGER.exception('error on pixel fetch for hab')
            if numpy.isclose(pixel_val, 0.0):
                pixel_val = 0
            # use max risk if no coverage
            shore_point_feature_risk_map[shore_feature.GetFID()].append(
                risk_val if pixel_val else 5)
            shore_feature.SetField(hab_id, risk_val if pixel_val else 5)
            shore_point_layer.SetFeature(shore_feature)
        shore_point_layer.CommitTransaction()

    shore_point_layer.CommitTransaction()
    shore_point_layer = None
    shore_point_vector = None
    hab_effective_raster = None
    hab_effective_band = None
    retrying_rmtree(tmp_working_dir)


@retrying.retry(
    wait_exponential_multiplier=100, wait_exponential_max=2000,
    stop_max_attempt_number=5)
def retrying_rmtree(dir_path):
    """Remove `dir_path` but try a few times."""
    try:
        shutil.rmtree(dir_path)
    except Exception:
        LOGGER.exception('unable to remove %s' % dir_path)
        raise


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


def clip_geometry(
        bounding_box_coords, base_srs, target_srs, ogr_geometry_type,
        global_geom_strtree, target_vector_path):
    """Clip geometry in `global_geom_strtree` to bounding box.

    Parameters:
        bounding_box_coords (list): a list of bounding box coordinates in
            the same coordinate system as the geometry in
            `global_geom_strtree`.
        target_srs (osr.SpatialReference): target spatial reference for
            creating the target vector.
        ogr_geometry_type (ogr.wkb[TYPE]): geometry type to create for the
            target vector.
        global_geom_strtree (shapely.strtree.STRtree): an rtree loaded with
            geometry to query via bounding box. Each geometry will contain
            parameters `field_val_map` and `prep` that have values to copy to
            `target_fieldname` and used to quickly query geometry. Main object
            will have `field_name_type_list` field used to describe the
            original field name/types.
        target_vector_path (str): path to vector to create that will contain
            locally projected geometry clipped to the given bounding box.

    Returns:
        None.

    """
    gpkg_driver = ogr.GetDriverByName("GPKG")
    vector = gpkg_driver.CreateDataSource(
        target_vector_path)
    layer = vector.CreateLayer(
        os.path.splitext(os.path.basename(target_vector_path))[0],
        target_srs, ogr_geometry_type)
    for field_name, field_type in global_geom_strtree.field_name_type_list:
        layer.CreateField(ogr.FieldDefn(field_name, field_type))
    layer_defn = layer.GetLayerDefn()
    base_to_target_transform = osr.CoordinateTransformation(
        base_srs, target_srs)

    bounding_box = shapely.geometry.box(*bounding_box_coords)

    possible_geom_list = global_geom_strtree.query(bounding_box)
    LOGGER.debug('possible intersections %d', len(possible_geom_list))
    if not possible_geom_list:
        layer = None
        vector = None
        raise ValueError('no data intersects this box')
    for geom in possible_geom_list:
        clipped_shapely_geom = bounding_box.intersection(geom)
        clipped_geom = ogr.CreateGeometryFromWkb(clipped_shapely_geom.wkb)
        error_code = clipped_geom.Transform(base_to_target_transform)
        if error_code:
            raise RuntimeError(error_code)
        feature = ogr.Feature(layer_defn)
        feature.SetGeometry(clipped_geom.Clone())
        for field_name, _ in global_geom_strtree.field_name_type_list:
            feature.SetField(
                field_name, geom.field_val_map[field_name])
        layer.CreateFeature(feature)


def sample_line_to_points(
        line_vector_path, target_point_path, point_step_size):
    """Sample lines in line vector to points along the path.

    Parameters:
        line_vector_path (str): path to line based vector.
        target_point_path (str): created by this function. A GPKG that is in
            the same projection as `line_vector` where points lie on those line
            segments spaced no less than `point_step_size` apart.
        point_step_size (float): step size in projected units of `line_vector`
            for points to drop along the line segments.

    Returns:
        None.

    """
    line_vector = gdal.OpenEx(line_vector_path)
    line_layer = line_vector.GetLayer()
    layer_name = os.path.splitext(os.path.basename(line_vector_path))[0]

    gpkg_driver = ogr.GetDriverByName('GPKG')
    if os.path.exists(target_point_path):
        os.remove(target_point_path)
    point_vector = gpkg_driver.CreateDataSource(target_point_path)
    point_layer = point_vector.CreateLayer(
        layer_name, line_layer.GetSpatialRef(), ogr.wkbPoint,
        ['OVERWRITE=YES'])
    point_defn = point_layer.GetLayerDefn()
    for feature in line_layer:
        current_distance = 0.0
        line_geom = feature.GetGeometryRef()
        line = shapely.wkb.loads(line_geom.ExportToWkb())
        if isinstance(line, shapely.geometry.collection.GeometryCollection):
            line_list = []
            for geom in list(line):
                LOGGER.debug(geom)
                if isinstance(geom, (
                        shapely.geometry.linestring.LineString,
                        shapely.geometry.multilinestring.MultiLineString)):
                    LOGGER.debug('appending')
                    line_list.append(geom)
            print('building: %s', line_list)
            line = shapely.geometry.MultiLineString(line_list)
        while current_distance < line.length:
            try:
                new_point = line.interpolate(current_distance)
                current_distance += point_step_size
                new_point_feature = ogr.Feature(point_defn)
                new_point_geom = ogr.CreateGeometryFromWkb(new_point.wkb)
                new_point_feature.SetGeometry(new_point_geom)
                point_layer.CreateFeature(new_point_feature)
            except Exception:
                LOGGER.exception('error on %s', line_geom)
                raise

    point_layer = None
    point_vector = None
    line_layer = None
    line_vector = None


def calculate_relief(
        shore_point_vector_path, dem_path, target_fieldname):
    """Calculate DEM relief as average coastal land area within 5km.

    Parameters:
        shore_point_vector_path (string):  path to a point shapefile to
            for relief point analysis.
        dem_path (string): path to a DEM raster projected in local coordinates.
        target_fieldname (string): this field name will be added to
            `shore_point_vector_path` and filled with Relief values.

    Returns:
        None.

    """
    try:
        shore_point_vector = gdal.OpenEx(
            shore_point_vector_path, gdal.OF_VECTOR | gdal.GA_Update)
        shore_point_layer = shore_point_vector.GetLayer()
        relief_field = ogr.FieldDefn(target_fieldname, ogr.OFTReal)
        relief_field.SetPrecision(5)
        relief_field.SetWidth(24)
        shore_point_layer.CreateField(relief_field)
        dem_info = pygeoprocessing.get_raster_info(dem_path)

        tmp_working_dir = tempfile.mkdtemp(
            prefix='calculate_relief_',
            dir=os.path.dirname(shore_point_vector_path))

        dem_nodata = dem_info['nodata'][0]

        def zero_negative_values(depth_array):
            valid_mask = depth_array != dem_nodata
            result_array = numpy.empty_like(depth_array)
            result_array[:] = dem_nodata
            result_array[valid_mask] = 0
            result_array[depth_array > 0] = depth_array[depth_array > 0]
            return result_array

        positive_dem_path = os.path.join(
            tmp_working_dir, 'positive_dem.tif')

        pygeoprocessing.raster_calculator(
            [(dem_path, 1)], zero_negative_values,
            positive_dem_path, gdal.GDT_Int16, dem_nodata)

        # convolve over a 5km radius
        dem_pixel_size = dem_info['pixel_size']
        kernel_radius = (
            abs(RELIEF_SAMPLE_DISTANCE // dem_pixel_size[0]),
            abs(RELIEF_SAMPLE_DISTANCE // dem_pixel_size[1]))

        kernel_filepath = os.path.join(
            tmp_working_dir, 'averaging_kernel.tif')
        create_averaging_kernel_raster(
            kernel_radius, kernel_filepath, normalize=True)

        relief_path = os.path.join(tmp_working_dir, 'relief.tif')
        pygeoprocessing.convolve_2d(
            (positive_dem_path, 1), (kernel_filepath, 1), relief_path)
        relief_raster = gdal.Open(relief_path)
        relief_band = relief_raster.GetRasterBand(1)

        inv_gt = gdal.InvGeoTransform(dem_info['geotransform'])

        shore_point_layer.ResetReading()
        for point_feature in shore_point_layer:
            point_geometry = point_feature.GetGeometryRef()
            point_x, point_y = point_geometry.GetX(), point_geometry.GetY()
            point_geometry = None

            pixel_x, pixel_y = [
                int(x) for x in
                gdal.ApplyGeoTransform(inv_gt, point_x, point_y)]

            if pixel_x < 0:
                pixel_x = 0
            if pixel_y < 0:
                pixel_y = 0
            if pixel_x >= relief_band.XSize:
                pixel_x = relief_band.XSize-1
            if pixel_y >= relief_band.YSize:
                pixel_y = relief_band.YSize-1

            try:
                pixel_value = relief_band.ReadAsArray(
                    xoff=pixel_x, yoff=pixel_y, win_xsize=1,
                    win_ysize=1)[0, 0]
            except Exception:
                LOGGER.exception(
                    'relief_band size %d %d', relief_band.XSize,
                    relief_band.YSize)
                raise
            # Make relief "negative" so when we histogram it for risk a
            # "higher" value will show a lower risk.
            point_feature.SetField(target_fieldname, -float(pixel_value))
            shore_point_layer.SetFeature(point_feature)

        shore_point_layer.SyncToDisk()
        shore_point_layer = None
        shore_point_vector = None
        relief_raster = None
        relief_band = None

        try:
            retrying_rmtree(tmp_working_dir)
        except OSError:
            LOGGER.warning('unable to rm %s' % tmp_working_dir)

    except Exception:
        LOGGER.exception('error in relief calc')
        raise


def clip_and_reproject_raster(
        base_raster_path, target_raster_path, target_srs_wkt,
        target_bounding_box, edge_buffer, resample_method,
        reproject_bounding_box, target_pixel_size):
    """Clip and reproject base to target raster.

    Parameters:
        base_raster_path (str): path to the raster to clip from.
        target_raster_path (str): path to target raster that is a clip from
            base projected in `target_srs_wkt` coordinate system.
        target_srs_wkt (str): spatial reference of target coordinate system in
            wkt.
        target_bounding_box (list): List of float describing target bounding
            box in base coordinate system as [minx, miny, maxx, maxy].
        edge_buffer (float): amount to extend sides of bounding box in target
            coordinate system units.
        resample_method (str): one of
            "near|bilinear|cubic|cubicspline|lanczos|mode".
        reproject_bounding_box (bool): If true, project `target_bounding_box`
            from base coordinate system to `target_srs_wkt`.
        target_pixel_size (float): desired target pixel size in projected
            coordinates.

    Returns:
        None.

    """
    base_raster_info = pygeoprocessing.get_raster_info(base_raster_path)
    bb_centroid = (
        (target_bounding_box[0]+target_bounding_box[2])/2,
        (target_bounding_box[1]+target_bounding_box[3])/2)

    if reproject_bounding_box:
        local_bounding_box = pygeoprocessing.transform_bounding_box(
            target_bounding_box, base_raster_info['projection'],
            target_srs_wkt, edge_samples=11)
    else:
        local_bounding_box = target_bounding_box
        base_srs = osr.SpatialReference()
        base_srs.ImportFromWkt(base_raster_info['projection'])
        target_srs = osr.SpatialReference()
        target_srs.ImportFromWkt(target_srs_wkt)
        target_to_base_transform = osr.CoordinateTransformation(
            target_srs, base_srs)
        point = ogr.CreateGeometryFromWkt("POINT (%f %f)" % bb_centroid)
        point.Transform(target_to_base_transform)
        bb_centroid = (point.GetX(), point.GetY())

    buffered_bounding_box = [
        local_bounding_box[0]-edge_buffer,
        local_bounding_box[1]-edge_buffer,
        local_bounding_box[2]+edge_buffer,
        local_bounding_box[3]+edge_buffer,
    ]

    # target_pixel_size = estimate_projected_pixel_size(
    #     base_raster_path, bb_centroid, target_srs_wkt)
    pygeoprocessing.warp_raster(
        base_raster_path, target_pixel_size, target_raster_path,
        resample_method, target_bb=buffered_bounding_box,
        target_sr_wkt=target_srs_wkt,
        working_dir=os.path.dirname(target_raster_path))


def clip_raster(
        base_raster_path, target_raster_path,
        target_bounding_box, edge_buffer):
    """Clip and reproject base to target raster.

    Parameters:
        base_raster_path (str): path to the raster to clip from.
        target_raster_path (str): path to target raster that is a clip from
            base projected in `target_srs_wkt` coordinate system.
        target_srs_wkt (str): spatial reference of target coordinate system in
            wkt.
        target_bounding_box (list): List of float describing target bounding
            box in base coordinate system as [minx, miny, maxx, maxy].
        edge_buffer (float): amount to extend sides of bounding box in target
            coordinate system units.
        resample_method (str): one of
            "near|bilinear|cubic|cubicspline|lanczos|mode".
        target_bounding_box (bool): If True, assumes bounding box is in
            base coordinate system and will transform it to target.

    Returns:
        None.

    """
    buffered_bounding_box = [
        target_bounding_box[0]-edge_buffer,
        target_bounding_box[1]-edge_buffer,
        target_bounding_box[2]+edge_buffer,
        target_bounding_box[3]+edge_buffer,
    ]

    base_raster = gdal.OpenEx(base_raster_path)
    gdal.Translate(
        target_raster_path, base_raster,
        format='GTiff',
        creationOptions=[
            'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=LZW',
            'BLOCKXSIZE=256', 'BLOCKYSIZE=256'],
        outputBounds=buffered_bounding_box,
        callback=pygeoprocessing._make_logger_callback(
            "Translate %.1f%% complete"))
    base_raster = None


def estimate_projected_pixel_size(
        base_raster_path, sample_point, target_srs_wkt):
    """Estimate the pixel size of raster if projected in `target_srs_wkt`.

    Parameters:
        base_raster_path (str): path to a raster in some coordinate system.
        sample_point (list): [x, y] coordinate in base coordinate system of
            point to estimate projected pixel size around.
        target_srs_wkt (str): desired target coordinate system in wkt for
            estimate pixel size.

    Returns:
        None.

    """
    base_raster_info = pygeoprocessing.get_raster_info(base_raster_path)
    base_pixel_size = base_raster_info['pixel_size']
    raster_center_pixel_bb = [
        sample_point[0] - abs(base_pixel_size[0]/2),
        sample_point[1] - abs(base_pixel_size[1]/2),
        sample_point[0] + abs(base_pixel_size[0]/2),
        sample_point[1] + abs(base_pixel_size[1]/2),
    ]
    pixel_bb = pygeoprocessing.transform_bounding_box(
        raster_center_pixel_bb, base_raster_info['projection'], target_srs_wkt)
    # x goes to the right, y goes down
    estimated_pixel_size = [
        pixel_bb[2]-pixel_bb[0],
        pixel_bb[1]-pixel_bb[3]]
    return estimated_pixel_size


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

    LOGGER.debug(f'kernel filepath: {kernel_filepath}')
    LOGGER.debug(f'kernel file exists? {os.path.exists(kernel_filepath)}')
    LOGGER.debug(
        f'kernel dir exists? '
        f'{os.path.exists(os.path.dirname(kernel_filepath))}')

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


def vector_to_lines(base_vector_path, target_line_vector_path):
    """Convert polygon vector to list of lines.

    Parameters:
        base_vector_path (str): path to polygon vector.
        target_line_vector_path (str): created by this file all polygons are
            converted to their line boundary equivalents.

    Returns:
        None.

    """
    # explode landmass into lines for easy intersection
    base_vector = gdal.OpenEx(base_vector_path, gdal.OF_VECTOR)
    base_layer = base_vector.GetLayer()

    gpkg_driver = ogr.GetDriverByName('GPKG')
    line_vector = gpkg_driver.CreateDataSource(
        target_line_vector_path)
    line_layer = line_vector.CreateLayer(
        target_line_vector_path, base_layer.GetSpatialRef(), ogr.wkbLineString)
    line_vector_defn = line_layer.GetLayerDefn()

    line_layer.StartTransaction()
    for base_feature in base_layer:
        base_shapely = shapely.wkb.loads(
            base_feature.GetGeometryRef().ExportToWkb())
        for line in geometry_to_lines(base_shapely):
            segment_feature = ogr.Feature(line_vector_defn)
            segement_geometry = ogr.Geometry(ogr.wkbLineString)
            segement_geometry.AddPoint(*line.coords[0])
            segement_geometry.AddPoint(*line.coords[1])
            segment_feature.SetGeometry(segement_geometry)
            line_layer.CreateFeature(segment_feature)
    line_layer.CommitTransaction()
    line_layer = None
    line_vector = None
    base_vector = None
    base_layer = None


def geometry_to_lines(geometry):
    """Convert a geometry object to a list of lines."""
    if geometry.type == 'Polygon':
        return polygon_to_lines(geometry)
    elif geometry.type == 'MultiPolygon':
        line_list = []
        for geom in geometry.geoms:
            line_list.extend(geometry_to_lines(geom))
        return line_list
    else:
        return []


def polygon_to_lines(geometry):
    """Return a list of shapely lines given higher order shapely geometry."""
    line_list = []
    last_point = geometry.exterior.coords[0]
    for point in geometry.exterior.coords[1::]:
        if point == last_point:
            continue
        line_list.append(shapely.geometry.LineString([last_point, point]))
        last_point = point
    line_list.append(shapely.geometry.LineString([
        last_point, geometry.exterior.coords[0]]))
    for interior in geometry.interiors:
        last_point = interior.coords[0]
        for point in interior.coords[1::]:
            if point == last_point:
                continue
            line_list.append(shapely.geometry.LineString([last_point, point]))
            last_point = point
        line_list.append(shapely.geometry.LineString([
            last_point, interior.coords[0]]))
    return line_list


def compute_wave_height(Un, Fn, dn):
    """Compute Wave Height by User Guide eq 10.

    This equation may not be suitable for wind speed values < 1 m/s
    The WWIII database tends to include some 0s, otherwise values > 2.

    Parameters:
        Un (float): wind velocity in meters per second.
        Fn (float): fetch ray length in meters.
        dn (float): water depth in negative meters.

    Returns:
        Float: Wave height in meters

    """
    if Un < 1.0:
        LOGGER.warning(
            'Found wind velocity of %.2f, '
            'using 1.0m/s in wave height calculation instead' % Un)
        Un = 1.0
    g = 9.81
    dn = -dn
    ds = g*dn/Un**2
    Fs = g*Fn/Un**2
    A = numpy.tanh(0.343*ds**1.14)
    B = numpy.tanh(4.41e-4*Fs**0.79/A)
    H_n = (0.24*Un**2/g)*(A*B)**0.572
    return H_n


def compute_wave_period(Un, Fn, dn):
    """Compute Wave Period by User Guide eq 10.

    This equation may not be suitable for wind speed values < 1 m/s
    The WWIII database tends to include some 0s, otherwise values > 2.

    Parameters:
        Un (float): wind velocity in meters per second.
        Fn (float): fetch ray length in meters.
        dn (float): water depth in negative meters.

    Returns:
        Float: Wave period in seconds

    """
    # This equation may not be suitable for wind speed values < 1 m/s
    # The WWIII database tends to include some 0s, otherwise values > 2
    if Un < 1.0:
        LOGGER.warning(
            'Found wind velocity of %.2f, '
            'using 1.0m/s in wave height calculation instead' % Un)
        Un = 1.0
    g = 9.81
    dn = -dn
    ds = g*dn/Un**2
    Fs = g*Fn/Un**2
    A = numpy.tanh(0.1*ds**2.01)
    B = numpy.tanh(2.77e-7*Fs**1.45/A)
    T_n = 7.69*Un/g*(A*B)**0.187
    return T_n


def merge_masks_op(mask_a, mask_b, nodata_a, nodata_b, target_nodata):
    result = numpy.empty(mask_a.shape, dtype=numpy.int16)
    valid_mask = (~numpy.isclose(mask_a, nodata_a) |
                  ~numpy.isclose(mask_b, nodata_b))
    result[:] = target_nodata
    result[valid_mask] = 1
    return result


def merge_cv_points(cv_vector_queue, target_cv_vector_path):
    """Merge vectors in `cv_vector_queue` into single vector.

    Parameters:
        cv_vector_queue (multiprocessing.Processing): a queue containing
            paths to CV workspace point vectors. Terminated with
            `STOP_SENTINEL`.
        target_cv_vector_path (str): path to a point vector created by this
            function.

    Returns:
        None.

    """
    gpkg_driver = ogr.GetDriverByName('GPKG')
    target_cv_vector = gpkg_driver.CreateDataSource(target_cv_vector_path)
    layer_name = os.path.basename(os.path.splitext(target_cv_vector_path)[0])
    wgs84_srs = osr.SpatialReference()
    wgs84_srs.ImportFromEPSG(4326)
    target_cv_layer = (
        target_cv_vector.CreateLayer(layer_name, wgs84_srs, ogr.wkbPoint))
    fields_to_copy = [
        'Rgeomorphology', 'surge', 'ew', 'rei', 'slr', 'relief',
        '4_500', '2_2000', 'reefs', 'mangroves_forest', 'saltmarsh_wetland',
        'seagrass', 'mesoamerican_barrier_reef', 'new_caledonian_barrier_reef',
        'great_barrier_reef', 'keys_barrier_reef',
        ]

    for field_id in fields_to_copy:
        target_cv_layer.CreateField(ogr.FieldDefn(field_id, ogr.OFTReal))
    target_cv_layer_defn = target_cv_layer.GetLayerDefn()

    target_cv_layer.StartTransaction()
    while True:
        cv_vector_path = cv_vector_queue.get()
        if cv_vector_path == STOP_SENTINEL:
            break
        cv_vector = gdal.OpenEx(cv_vector_path, gdal.OF_VECTOR)
        cv_layer = cv_vector.GetLayer()
        cv_projection = cv_layer.GetSpatialRef()
        base_to_target_transform = osr.CoordinateTransformation(
            cv_projection, wgs84_srs)

        for cv_feature in cv_layer:
            cv_geom = cv_feature.GetGeometryRef().Clone()
            _ = cv_geom.Transform(base_to_target_transform)
            target_feature = ogr.Feature(target_cv_layer_defn)
            target_feature.SetGeometry(cv_geom)
            for field_id in fields_to_copy:
                target_feature.SetField(
                    field_id, cv_feature.GetField(field_id))
            target_cv_layer.CreateFeature(target_feature)
        cv_feature = None
        cv_geom = None
        cv_layer = None
        cv_vector = None
    target_cv_layer.CommitTransaction()


def add_cv_vector_risk(cv_risk_vector_path):
    """Use existing biophysical fields in `cv_risk_vector_path to calc total R

    Parameters:
        cv_risk_vector_path (str): path to point vector that has at least
            the following fields in it:

            * surge
            * ew
            * rei
            * slr
            * relief

            Will add the following fields:
                * Rwave
                * Rwind
                * Rsurge
                * Rrelief
                * Rslr
    Returns:
        None

    """

    cv_risk_vector = gdal.OpenEx(
        cv_risk_vector_path, gdal.OF_VECTOR | gdal.GA_Update)
    cv_risk_layer = cv_risk_vector.GetLayer()

    for base_field, risk_field in [
            ('surge', 'Rsurge'), ('ew', 'Rwave'), ('rei', 'Rwind'),
            ('slr', 'Rslr'), ('relief', 'Rrelief')]:
        cv_risk_layer.CreateField(ogr.FieldDefn(risk_field, ogr.OFTReal))
        base_array = numpy.empty(shape=(cv_risk_layer.GetFeatureCount(),))
        for index, feature in enumerate(cv_risk_layer):
            base_array[index] = feature.GetField(base_field)
        nan_mask = numpy.isnan(base_array)
        max_val = numpy.max(base_array[~nan_mask])
        base_array[nan_mask] = max_val
        hist, bin_edges = numpy.histogram(base_array, bins=5)

        cv_risk_layer.ResetReading()
        cv_risk_layer.StartTransaction()
        for feature in cv_risk_layer:
            base_val = feature.GetField(base_field)
            risk = bisect.bisect_left(bin_edges, base_val)
            if risk < 1:
                risk = 1
            elif risk > 5:
                risk = 5
            feature.SetField(risk_field, risk)
            cv_risk_layer.SetFeature(feature)
        cv_risk_layer.CommitTransaction()
    cv_risk_layer.ResetReading()

    cv_risk_layer.CreateField(ogr.FieldDefn('Rt', ogr.OFTReal))
    for hab_field in HAB_FIELDS:
        cv_risk_layer.CreateField(
            ogr.FieldDefn('Rhab_%s' % hab_field, ogr.OFTReal))
        cv_risk_layer.CreateField(
            ogr.FieldDefn('Rnohab_%s' % hab_field, ogr.OFTReal))
        cv_risk_layer.CreateField(
            ogr.FieldDefn('Rt_nohab_%s' % hab_field, ogr.OFTReal))
        cv_risk_layer.CreateField(
            ogr.FieldDefn('Rt_habservice_%s' % hab_field, ogr.OFTReal))
    cv_risk_layer.CreateField(ogr.FieldDefn('Rhab_all', ogr.OFTReal))
    cv_risk_layer.CreateField(ogr.FieldDefn('Rt_nohab_all', ogr.OFTReal))

    # reefs are special
    cv_risk_layer.CreateField(ogr.FieldDefn('reefs_all', ogr.OFTReal))

    cv_risk_layer.ResetReading()
    cv_risk_layer.StartTransaction()
    for feature in cv_risk_layer:
        reef_risk_val = min(
            [feature.GetField(reef_field) for reef_field in REEF_FIELDS])
        feature.SetField('reefs_all', reef_risk_val)

        hab_val_map = {}
        for hab_field in HAB_FIELDS:
            hab_val = feature.GetField(hab_field)
            feature.SetField('Rhab_%s' % hab_field, hab_val)
            hab_val_map[hab_field] = hab_val

            # loop through every hab field but hab_field to calc Rhab_no
            risk_diff_list = []  # for (5-rk) vals
            for sub_hab_field in HAB_FIELDS:
                if sub_hab_field != hab_field:
                    risk_diff_list.append(5-feature.GetField(sub_hab_field))

            r_nohab = 4.8 - 0.5 * numpy.sqrt(
                (1.5 * max(risk_diff_list))**2 +
                numpy.sum([x**2 for x in risk_diff_list]) -
                max(risk_diff_list)**2)
            feature.SetField('Rnohab_%s' % hab_field, r_nohab)

        # Rhab
        # loop through every hab field but hab_field to calc Rhab_no
        risk_diff_list = []  # for (5-rk) vals
        for sub_hab_field in HAB_FIELDS:
            risk_diff_list.append(5-feature.GetField(sub_hab_field))

        r_nohab = 4.8 - 0.5 * numpy.sqrt(
            (1.5 * max(risk_diff_list))**2 +
            numpy.sum([x**2 for x in risk_diff_list]) -
            max(risk_diff_list)**2)
        feature.SetField('Rhab_all', r_nohab)

        # Rt
        exposure_index = 1.0
        for risk_field in [
                'Rgeomorphology', 'Rhab_all', 'Rsurge', 'Rwave', 'Rwind',
                'Rslr', 'Rrelief']:
            exposure_index *= feature.GetField(risk_field)
        exposure_index = (exposure_index)**(1./7.)
        feature.SetField('Rt', exposure_index)

        # Rt_nohaball
        nohab_exposure_index = 1.0
        for risk_field in [
                'Rgeomorphology', 'Rsurge', 'Rwave', 'Rwind', 'Rslr',
                'Rrelief']:
            nohab_exposure_index *= feature.GetField(risk_field)
        # the *5 is to get the "missing" habitat risk in there
        nohab_exposure_index = (nohab_exposure_index*5)**(1./7.)
        feature.SetField('Rt_nohab_all', nohab_exposure_index)

        for hab_field in HAB_FIELDS:
            nohab_exposure_index = 1.0
            for risk_field in [
                    'Rgeomorphology', 'Rsurge', 'Rwave', 'Rwind', 'Rslr',
                    'Rrelief', 'Rnohab_%s' % hab_field]:
                nohab_exposure_index *= feature.GetField(risk_field)
            nohab_exposure_index = (nohab_exposure_index)**(1./7.)
            feature.SetField('Rt_nohab_%s' % hab_field, nohab_exposure_index)
            # service is the difference between Rt without the habitat and
            # Rt with all habitats.
            hab_service = (nohab_exposure_index - feature.GetField('Rt'))
            if numpy.isclose(hab_service, 0.0):
                hab_service = 0.0
            feature.SetField('Rt_habservice_%s' % hab_field, hab_service)

        cv_risk_layer.SetFeature(feature)
    cv_risk_layer.CommitTransaction()


def calculate_habitat_population_value(
        shore_sample_point_vector_path, population_raster_path_id_list,
        dem_raster_path, habitat_fieldname_list, habitat_vector_path_map,
        results_dir, habitat_pop_value_token_path):
    """Calculate population within protective range of habitat.

    Parameters:
        shore_sample_point_vector_path (str): path to a point shapefile that is only
            used for referencing the points of interest on the coastline.
        population_raster_path_id_list (list): list of (raster_path, field_id)
            tuples. The values in the raster paths will be masked where it
            overlaps with < 10m dem height and convolved within 2km. That
            result is in turn spread onto the habitat coverage at a distance
            of the protective distance of that habitat. These rasters are in
            wgs84 lat/lng projection.
        dem_raster_path (str): path to a dem used to mask population by height
            in wgs84 lat/lng projection.
        habitat_fieldname_list (list): list of habitat ids to analyse.
        habitat_vector_path_map (dict): maps fieldnames from
            `habitat_fieldname_list` to 3-tuples of
            (path to hab raster (str), risk val (float),
             protective distance (float)).
        results_dir (str): path to directory containing habitat back projection
            results
        habitat_pop_value_token_path (str): path to a file to create if the
            run is successful.

    Returns:
        None

    """
    temp_workspace_dir = os.path.join(
        results_dir, 'calc_pop_coverage_churn')
    taskgraph_working_dir = os.path.join(temp_workspace_dir, 'taskgraph')
    for path in [results_dir, taskgraph_working_dir, temp_workspace_dir]:
        try:
            os.makedirs(results_dir)
        except OSError:
            pass
    LOGGER.info(
        f'starting taskgraph in calculate_habitat_population_value for '
        f'{taskgraph_working_dir}')
    task_graph = taskgraph.TaskGraph(taskgraph_working_dir, -1)

    aligned_pop_raster_list = align_raster_list(
        [x[0] for x in population_raster_path_id_list] + [dem_raster_path],
        temp_workspace_dir)

    for pop_index, (_, pop_id) in enumerate(population_raster_path_id_list):
        # mask to < 10m
        pop_height_masked_path = os.path.join(
            temp_workspace_dir, '%s_masked_by_10m.tif' % pop_id)
        raster_info = pygeoprocessing.get_raster_info(
            aligned_pop_raster_list[pop_index])

        LOGGER.info('pop height mask 10m %s' % pop_id)
        pop_height_mask_task = task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                [(aligned_pop_raster_list[pop_index], 1),
                 (aligned_pop_raster_list[-1], 1),
                 (10.0, 'raw'),  # mask to 10 meters
                 (raster_info['nodata'][0], 'raw')],  # the -1 index is the dem
                mask_by_height_op, pop_height_masked_path, gdal.GDT_Float32,
                raster_info['nodata'][0]),
            target_path_list=[pop_height_masked_path],
            task_name='pop height mask 10m %s' % pop_id)
        pop_height_mask_task.join()

        target_pixel_size = raster_info['pixel_size']
        # spread the < 10m population out 2km
        n_pixels_in_2km = int(2000.0 / (
            M_PER_DEGREE * abs(target_pixel_size[0])))
        kernel_radius_2km = [n_pixels_in_2km, n_pixels_in_2km]
        kernel_2km_filepath = os.path.join(
            temp_workspace_dir, '2km_kernel.tif')
        create_averaging_kernel_raster(
            kernel_radius_2km, kernel_2km_filepath, normalize=False)
        pop_sum_within_2km_path = os.path.join(
            temp_workspace_dir, '%s_pop_sum_within_2km.tif' % pop_id)
        LOGGER.info('pop sum w/in 2km %s' % pop_id)
        pop_sum_task = task_graph.add_task(
            func=pygeoprocessing.convolve_2d,
            args=(
                (pop_height_masked_path, 1), (kernel_2km_filepath, 1),
                pop_sum_within_2km_path),
            kwargs={'working_dir': temp_workspace_dir},
            target_path_list=[pop_sum_within_2km_path],
            task_name='pop sum w/in 2km %s' % pop_id)

        # spread the 2km pop out by the hab distance
        for habitat_id in habitat_fieldname_list:
            hab_raster_path, _, prot_distance = (
                habitat_vector_path_map[habitat_id])
            # make a kernel that goes out the distance of the protective
            # distance of habitat
            n_pixels_in_prot_dist = max(1, int(prot_distance / (
                M_PER_DEGREE * abs(target_pixel_size[0]))))
            kernel_radius = [n_pixels_in_prot_dist, n_pixels_in_prot_dist]
            kernel_filepath = os.path.join(
                temp_workspace_dir, '%s_kernel.tif' % habitat_id)
            create_averaging_kernel_raster(
                kernel_radius, kernel_filepath, normalize=False)
            population_hab_spread_raster_path = os.path.join(
                temp_workspace_dir, '%s_%s_spread.tif' % (habitat_id, pop_id))
            LOGGER.info('spread pop to hab %s %s ' % (pop_id, habitat_id))
            spread_to_hab_task = task_graph.add_task(
                func=clean_convolve_2d,
                args=(
                    (pop_sum_within_2km_path, 1), (kernel_filepath, 1),
                    population_hab_spread_raster_path),
                kwargs={'working_dir': temp_workspace_dir},
                target_path_list=[population_hab_spread_raster_path],
                dependent_task_list=[pop_sum_task],
                task_name='spread pop to hab %s %s ' % (pop_id, habitat_id))

            hab_raster_info = pygeoprocessing.get_raster_info(hab_raster_path)

            # warp pop result to overlay
            clipped_pop_hab_spread_raster_path = os.path.join(
                temp_workspace_dir, '%s_%s_spread_clipped.tif' % (
                    habitat_id, pop_id))
            LOGGER.info('spread to hab %s %s' % (pop_id, habitat_id))
            task_graph.add_task(
                func=pygeoprocessing.warp_raster,
                args=(
                    population_hab_spread_raster_path,
                    hab_raster_info['pixel_size'],
                    clipped_pop_hab_spread_raster_path,
                    'near'),
                kwargs={
                    'target_bb': hab_raster_info['bounding_box'],
                    'working_dir': temp_workspace_dir},
                target_path_list=[clipped_pop_hab_spread_raster_path],
                dependent_task_list=[spread_to_hab_task],
                task_name='spread to hab %s %s' % (pop_id, habitat_id))

            # mask the convolution by the habitat mask
            task_graph.join()
            hab_spread_nodata = pygeoprocessing.get_raster_info(
                clipped_pop_hab_spread_raster_path)['nodata'][0]
            hab_nodata = pygeoprocessing.get_raster_info(
                hab_raster_path)['nodata'][0]
            habitat_value_raster_path = os.path.join(
                results_dir, '%s_%s_coverage.tif' % (habitat_id, pop_id))

            buffered_point_raster_mask_path = os.path.join(
                temp_workspace_dir, '%s_buffer_mask.tif' % habitat_id)
            buffered_raster_task = task_graph.add_task(
                func=make_buffered_point_raster_mask,
                args=(
                    shore_sample_point_vector_path,
                    clipped_pop_hab_spread_raster_path,
                    temp_workspace_dir, habitat_id,
                    prot_distance, buffered_point_raster_mask_path),
                target_path_list=[buffered_point_raster_mask_path],
                task_name='buffered point for %s' % habitat_id)

            hab_value_task = task_graph.add_task(
                func=pygeoprocessing.raster_calculator,
                args=(
                    [(clipped_pop_hab_spread_raster_path, 1),
                     (hab_raster_path, 1),
                     (buffered_point_raster_mask_path, 1),
                     (hab_spread_nodata, 'raw'),
                     (hab_nodata, 'raw')],
                    intersect_and_mask_raster_op, habitat_value_raster_path,
                    gdal.GDT_Float32, hab_spread_nodata),
                dependent_task_list=[buffered_raster_task],
                target_path_list=[habitat_value_raster_path],
                task_name='mask result %s %s' % (pop_id, habitat_id))

            task_graph.add_task(
                func=ecoshard.build_overviews,
                args=(habitat_value_raster_path,),
                target_path_list=[habitat_value_raster_path],
                dependent_task_list=[hab_value_task],
                task_name='build overviews for %s' % habitat_value_raster_path)

    task_graph.join()
    task_graph.close()
    del task_graph
    with open(habitat_pop_value_token_path, 'w') as hab_pop_token_file:
        hab_pop_token_file.write(str(datetime.datetime.now()))


def mask_by_height_op(pop_array, dem_array, mask_height, pop_nodata):
    """Set pop to 0 if > height."""
    result = numpy.zeros(shape=pop_array.shape)
    valid_mask = (
        (dem_array < mask_height) & ~numpy.isclose(pop_array, pop_nodata))
    result[valid_mask] = pop_array[valid_mask]
    return result


def calculate_habitat_value(
        shore_sample_point_vector, template_raster_path,
        habitat_fieldname_list, habitat_vector_path_map, results_dir,
        habitat_value_token_path):
    """Calculate habitat value.

    Will create rasters in the `results_dir` directory named from the
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
        results_dir (str): path to directory containing habitat back projection
            results
        habitat_value_token_path (str): path to file to write when done

    Returns:
        None.

    """
    temp_workspace_dir = os.path.join(results_dir, 'hab_value_churn')
    for dir_path in [results_dir, temp_workspace_dir]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    gpkg_driver = ogr.GetDriverByName('gpkg')
    shore_sample_point_vector = gdal.OpenEx(
        shore_sample_point_vector, gdal.OF_VECTOR)
    shore_sample_point_layer = shore_sample_point_vector.GetLayer()

    for habitat_id in habitat_fieldname_list:
        habitat_service_id = 'Rt_habservice_%s' % habitat_id
        hab_raster_path, _, protective_distance = (
            habitat_vector_path_map[habitat_id])

        buffer_habitat_path = os.path.join(
            temp_workspace_dir, '%s_buffer.gpkg' % habitat_id)
        buffer_habitat_vector = gpkg_driver.CreateDataSource(
            buffer_habitat_path)
        wgs84_srs = osr.SpatialReference()
        wgs84_srs.ImportFromEPSG(4326)
        buffer_habitat_layer = (
            buffer_habitat_vector.CreateLayer(
                habitat_service_id, wgs84_srs, ogr.wkbPolygon))
        buffer_habitat_layer.CreateField(ogr.FieldDefn(
            habitat_service_id, ogr.OFTReal))
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
            x_val = point_geom.GetX()
            if (x_val < -179.8) or (x_val > 179.8):
                continue
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

            # reefs are special
            if habitat_id in REEF_FIELDS:
                buffer_point_feature.SetField(
                    habitat_service_id,
                    point_feature.GetField('Rt_habservice_reefs_all'))
            else:
                buffer_point_feature.SetField(
                    habitat_service_id,
                    point_feature.GetField(habitat_service_id))
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
            temp_workspace_dir, '%s_value_cover.tif' % habitat_id)
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
                'ATTRIBUTE=%s' % habitat_service_id,
                'MERGE_ALG=ADD'])

        habitat_value_raster_path = os.path.join(
            results_dir, '%s_value.tif' % habitat_id)

        value_coverage_nodata = pygeoprocessing.get_raster_info(
            value_coverage_raster_path)['nodata'][0]
        hab_nodata = pygeoprocessing.get_raster_info(
            hab_raster_path)['nodata'][0]

        aligned_value_hab_raster_path_list = align_raster_list(
            [value_coverage_raster_path, hab_raster_path],
            temp_workspace_dir)

        pygeoprocessing.raster_calculator(
            [(aligned_value_hab_raster_path_list[0], 1),
             (aligned_value_hab_raster_path_list[1], 1),
             (value_coverage_nodata, 'raw'), (hab_nodata, 'raw')],
            intersect_raster_op, habitat_value_raster_path, gdal.GDT_Float32,
            value_coverage_nodata)

        ecoshard.build_overviews(habitat_value_raster_path)

    with open(habitat_value_token_path, 'w') as habitat_value_file:
        habitat_value_file.write(str(datetime.datetime.now()))


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


def download_data(**kwargs):
    """Download all data necessary for analysis to run.

    Parameters:
        kwargs (dict): dictionary of key/url pairs of data to download that
            will also be added to the result.

    Returns:
        dictionary mapping GLOBAL_DATA_URL_MAP keys to local paths for those
        artifacts.

    """
    task_graph = taskgraph.TaskGraph(CHURN_DIR, -1)

    local_data_path_map = {}
    for zip_url in [LS_POPULATION_RASTER_URL]:
        target_token_path = os.path.join(
            CHURN_DIR, os.path.basename(os.path.splitext(zip_url)[0]))
        _ = task_graph.add_task(
            func=download_and_unzip,
            args=(zip_url, ECOSHARD_DIR, target_token_path),
            target_path_list=[target_token_path],
            task_name='download and unzip %s' % zip_url)
    local_data_path_map['population'] = os.path.join(
        ECOSHARD_DIR, os.path.basename(LS_POPULATION_RASTER_URL))

    for data_id, ecoshard_url in {**GLOBAL_DATA_URL_MAP, **kwargs}.items():
        local_ecoshard_path = os.path.join(
            ECOSHARD_DIR, os.path.basename(ecoshard_url))
        _ = task_graph.add_task(
            func=ecoshard.download_url,
            args=(ecoshard_url, local_ecoshard_path),
            target_path_list=[local_ecoshard_path],
            task_name='download %s' % local_ecoshard_path)
        local_data_path_map[data_id] = local_ecoshard_path

    local_data_path_map['global_wwiii_vector_path'] = os.path.join(
        ECOSHARD_DIR, os.path.basename(
            os.path.splitext(GLOBAL_WWIII_GZ_URL)[0]))

    _ = task_graph.add_task(
        func=download_and_ungzip,
        args=(
            GLOBAL_WWIII_GZ_URL,
            local_data_path_map['global_wwiii_vector_path']),
        target_path_list=[local_data_path_map['global_wwiii_vector_path']],
        task_name=(
            'download %s' % local_data_path_map['global_wwiii_vector_path']))

    local_data_path_map['shore_buffer_vector_path'] = os.path.join(
        ECOSHARD_DIR, os.path.basename(BUFFER_VECTOR_URL))
    _ = task_graph.add_task(
        func=ecoshard.download_url,
        args=(
            BUFFER_VECTOR_URL,
            local_data_path_map['shore_buffer_vector_path']),
        target_path_list=[local_data_path_map['shore_buffer_vector_path']],
        task_name='download global_vector')

    reef_degree_pixel_size = [0.004, -0.004]
    wgs84_srs = osr.SpatialReference()
    wgs84_srs.ImportFromEPSG(4326)
    projected_reef_raster_path = os.path.join(CHURN_DIR, 'wgs84_reefs.tif')
    project_reef_task = task_graph.add_task(
        func=pygeoprocessing.warp_raster,
        args=(
            local_data_path_map['reefs'], reef_degree_pixel_size,
            projected_reef_raster_path, 'near'),
        kwargs={'target_sr_wkt': wgs84_srs.ExportToWkt()},
        target_path_list=[projected_reef_raster_path],
        task_name='project reefs to wgs84')
    project_reef_task.join()
    local_data_path_map['reefs'] = projected_reef_raster_path

    task_graph.join()
    task_graph.close()
    del task_graph

    return local_data_path_map


def align_raster_list(raster_path_list, target_directory):
    """Aligns all the raster paths.

    Rasters are aligned using the pixel size of the first raster and use
    the intersection and near interpolation methods.

    Parameters:
        raster_path_list (list): list of str paths to rasters.
        target_directory (str): path to a directory to hold the aligned
            rasters.

    Returns:
        list of raster paths that are aligned with intersection and near
            interpolation algorithm.

    """
    if not hasattr(align_raster_list, 'task_graph_map'):
        align_raster_list.task_graph_map = {}
    if target_directory not in align_raster_list.task_graph_map:
        align_raster_list.task_graph_map[target_directory] = (
            taskgraph.TaskGraph(target_directory, -1))
    task_graph = align_raster_list.task_graph_map[target_directory]
    aligned_path_list = [
        os.path.join(target_directory, os.path.basename(path))
        for path in raster_path_list]
    target_pixel_size = pygeoprocessing.get_raster_info(
        raster_path_list[0])['pixel_size']
    LOGGER.debug('about to align: %s', str(raster_path_list))
    task_graph.add_task(
        func=pygeoprocessing.align_and_resize_raster_stack,
        args=(
            raster_path_list, aligned_path_list,
            ['near'] * len(raster_path_list), target_pixel_size,
            'intersection'),
        target_path_list=aligned_path_list)
    return aligned_path_list


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


def set_almost_zero_to_zero(array, eps):
    result = numpy.empty_like(array)
    result[:] = numpy.where(numpy.abs(array) < eps, 0.0, array)
    return result


def preprocess_habitat():
    """Merge and filter all the habitat layers.

    Returns:
        dictionary of habitat id to (raster path, risk, dist) tuples.

    """
    task_graph = taskgraph.TaskGraph(CHURN_DIR, -1)
    lulc_shore_mask_raster_path = os.path.join(
        CHURN_DIR, 'lulc_masked_by_shore.tif')
    mask_lulc_by_shore_task = task_graph.add_task(
        func=pygeoprocessing.mask_raster,
        args=(
            (local_data_path_map['lulc'], 1),
            local_data_path_map['shore_buffer_vector_path'],
            lulc_shore_mask_raster_path),
        target_path_list=[lulc_shore_mask_raster_path],
        # ignore filetime stamp modification by a read-only open
        ignore_path_list=[
            local_data_path_map['shore_buffer_vector_path']],
        task_name='mask shore')

    # each value in `risk_distance_to_lulc_code` can be lumped into one
    # type.
    risk_distance_to_lulc_code = collections.defaultdict(list)
    for lulc_code, risk_distance in LULC_CODE_TO_HAB_MAP.items():
        risk_distance_to_lulc_code[risk_distance].append(lulc_code)

    # this maps all the same type of codes together
    habitat_raster_risk_map = dict(HABITAT_VECTOR_PATH_MAP)
    for risk_distance_tuple, lulc_code_list in sorted(
            risk_distance_to_lulc_code.items()):
        if risk_distance_tuple[0] == 0:
            LOGGER.info('skipping hab tuple %s', str(risk_distance_tuple))
            continue
        # reclassify landcover map to be ones everywhere for `lulc_code_list`
        reclass_map = {}
        for lulc_code in LULC_CODE_TO_HAB_MAP:
            if lulc_code in lulc_code_list:
                reclass_map[lulc_code] = 1
            else:
                reclass_map[lulc_code] = 0

        risk_distance_mask_path = os.path.join(
            CHURN_DIR, '%s_%s_mask.tif' % risk_distance_tuple)

        _ = task_graph.add_task(
            func=pygeoprocessing.reclassify_raster,
            args=(
                (lulc_shore_mask_raster_path, 1), reclass_map,
                risk_distance_mask_path, gdal.GDT_Byte, 0),
            target_path_list=[risk_distance_mask_path],
            dependent_task_list=[mask_lulc_by_shore_task],
            task_name='map distance types %s' % str(risk_distance_tuple))
        habitat_raster_risk_map[risk_distance_tuple] = (
            risk_distance_mask_path, risk_distance_tuple[0],
            risk_distance_tuple[1])

    for vector_name, lucode_type_tuple in [
            ('mangroves_forest', (1, 2000)),
            ('saltmarsh_wetland', (2, 1000))]:
        aligned_raster_path_list = [
            os.path.join(CHURN_DIR, x) for x in [
                '%s_a.tif' % vector_name, '%s_b.tif' % vector_name]]
        habitat_raster_info = pygeoprocessing.get_raster_info(
            habitat_raster_risk_map[lucode_type_tuple][0])
        align_task = task_graph.add_task(
            func=pygeoprocessing.align_and_resize_raster_stack,
            args=(
                [habitat_raster_risk_map[lucode_type_tuple][0],
                 local_data_path_map[vector_name]],
                aligned_raster_path_list, ['near', 'near'],
                habitat_raster_info['pixel_size'], 'union'),
            target_path_list=aligned_raster_path_list,
            task_name='align %s' % vector_name)

        merged_hab_raster_path = os.path.join(
            CHURN_DIR, 'merged_%s.tif' % vector_name)
        nodata_0 = pygeoprocessing.get_raster_info(
            aligned_raster_path_list[0])['nodata'][0]
        nodata_1 = pygeoprocessing.get_raster_info(
            aligned_raster_path_list[1])['nodata'][0]
        _ = task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                ((aligned_raster_path_list[0], 1),
                 (aligned_raster_path_list[1], 1),
                 (nodata_0, 'raw'),
                 (nodata_1, 'raw'),
                 (0, 'raw')), merge_masks_op,
                merged_hab_raster_path, gdal.GDT_Int16, 0),
            target_path_list=[merged_hab_raster_path],
            dependent_task_list=[align_task],
            task_name='merge %s masks' % vector_name)

        del habitat_raster_risk_map[lucode_type_tuple]
        habitat_raster_risk_map[vector_name] = (
            merged_hab_raster_path, lucode_type_tuple[0],
            lucode_type_tuple[1])

    task_graph.join()
    task_graph.close()
    del task_graph
    LOGGER.debug(habitat_raster_risk_map)

    # convert tuple to strings for habitat risk so we can make fields for them
    tuple_list = [
        hab_index for hab_index in habitat_raster_risk_map
        if isinstance(hab_index, tuple)]
    for tuple_index in tuple_list:
        str_index = '%s_%s' % tuple_index
        habitat_raster_risk_map[str_index] = (
            habitat_raster_risk_map[tuple_index])
        del habitat_raster_risk_map[tuple_index]
    return habitat_raster_risk_map


def calculate_degree_cell_cv(
        local_data_path_map, habitat_raster_risk_map, target_cv_vector_path):
    """Process all global degree grids to calculate local hab risk.

    Paramters:
        local_data_path_map (dict): maps keys from GLOBAL_DATA_URL to the
            local filepaths.
        habitat_raster_risk_map (dict): maps habitat layers to raster masks,
            indexed by habitat id maps to tuple of
            (raster path, risk val, distance).

    Returns:
        None

    """
    for dir_path in [
            WORKSPACE_DIR, CHURN_DIR, ECOSHARD_DIR, GRID_WORKSPACE_DIR]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    shore_grid_vector = gdal.OpenEx(
        local_data_path_map['shore_grid'], gdal.OF_VECTOR)
    shore_grid_layer = shore_grid_vector.GetLayer()

    bb_work_queue = multiprocessing.Queue()
    cv_point_complete_queue = multiprocessing.Queue()

    cv_grid_worker_list = []
    for worker_id in range(int(args.n_workers)):
        cv_grid_worker_thread = multiprocessing.Process(
            target=cv_grid_worker,
            args=(
                bb_work_queue,
                cv_point_complete_queue,
                local_data_path_map['landmass'],
                local_data_path_map['geomorphology'],
                local_data_path_map['slr'],
                local_data_path_map['dem'],
                local_data_path_map['global_wwiii_vector_path'],
                habitat_raster_risk_map,
                ))
        cv_grid_worker_thread.start()
        cv_grid_worker_list.append(cv_grid_worker_thread)
        LOGGER.debug('starting worker %d', worker_id)

    for path in [
            local_data_path_map['population'],
            local_data_path_map['lulc'],
            local_data_path_map['global_wwiii_vector_path'],
            local_data_path_map['landmass'],
            local_data_path_map['shore_grid']]:
        LOGGER.info('%s: %s' % (os.path.exists(path), path))

    shore_grid_vector = gdal.OpenEx(
        local_data_path_map['shore_grid'], gdal.OF_VECTOR)
    shore_grid_layer = shore_grid_vector.GetLayer()

    for index, shore_grid_feature in enumerate(shore_grid_layer):
        shore_grid_geom = shore_grid_feature.GetGeometryRef()
        boundary_box = shapely.wkb.loads(shore_grid_geom.ExportToWkb())
        LOGGER.debug(boundary_box.bounds)
        bb_work_queue.put((index, boundary_box.bounds))

    bb_work_queue.put(STOP_SENTINEL)

    merge_cv_points_thread = threading.Thread(
        target=merge_cv_points,
        args=(cv_point_complete_queue, target_cv_vector_path))
    merge_cv_points_thread.start()

    for cv_grid_worker_thread in cv_grid_worker_list:
        cv_grid_worker_thread.join()

    # when workers are complete signal merger complete
    cv_point_complete_queue.put(STOP_SENTINEL)
    merge_cv_points_thread.join()

    LOGGER.debug('calculate cv vector risk')
    add_cv_vector_risk(target_cv_vector_path)


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
        x_val = point_geom.GetX()
        if (x_val < -179.8) or (x_val > 179.8):
            continue

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
    parser = argparse.ArgumentParser(description='Global CV analysis')
    parser.add_argument(
        'landcover_file',
        help='Path to file that lists landcover scenarios to run.')
    parser.add_argument(
        '--n_workers', default=multiprocessing.cpu_count(),
        help='Number of workers.')
    parser.add_argument(
        '--shore_point_sample_distance', type=float, default=2000.0,
        help='Distance between shore sample points in meters.')
    parser.add_argument(
        '--dasgupta_mode', action='store_true',
        help='Ignore offshore mangrove and saltmarsh')

    args = parser.parse_args()

    SHORE_POINT_SAMPLE_DISTANCE = args.shore_point_sample_distance

    for dir_path in [
            WORKSPACE_DIR, CHURN_DIR, ECOSHARD_DIR, GRID_WORKSPACE_DIR]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    if args.dasgupta_mode:
        GLOBAL_MANGROVES_RASTER_URL = EMPTY_RASTER_URL
        GLOBAL_SALTMARSH_RASTER_URL = EMPTY_RASTER_URL

    task_graph = taskgraph.TaskGraph(WORKSPACE_DIR, 0, 5.0)

    try:
        with open(args.landcover_file, 'r') as landcover_raster_file:
            landcover_url_list = landcover_raster_file.read().splitlines()
        for landcover_url in landcover_url_list:
            local_data_path_map = download_data(lulc=landcover_url)
            landcover_basename = os.path.splitext(
                os.path.basename(landcover_url))[0]
            local_workspace_dir = os.path.join(
                WORKSPACE_DIR, landcover_basename)
            local_habitat_value_dir = os.path.join(
                WORKSPACE_DIR, landcover_basename, 'value_rasters')
            for dir_path in [local_workspace_dir, local_habitat_value_dir]:
                try:
                    os.makedirs(dir_path)
                except OSError:
                    pass
            target_cv_vector_path = os.path.join(
                local_workspace_dir, '%s.gpkg' % landcover_basename)
            habitat_raster_risk_dist_map = preprocess_habitat()
            calculate_cv_vector_task = task_graph.add_task(
                func=calculate_degree_cell_cv,
                args=(
                    local_data_path_map, habitat_raster_risk_dist_map,
                    target_cv_vector_path),
                target_path_list=[target_cv_vector_path],
                task_name='calculate CV for %s' % landcover_basename)

            LOGGER.info('calculating population back projection')
            ls_population_raster_path = os.path.join(ECOSHARD_DIR, 'lspop2017')
            poor_population_raster_path = os.path.join(
                ECOSHARD_DIR, os.path.basename(POVERTY_POPULATION_RASTER_URL))
            global_dem_raster_path = os.path.join(
                ECOSHARD_DIR, os.path.basename(GLOBAL_DEM_RASTER_URL))
            habitat_pop_value_token_path = os.path.join(
                local_habitat_value_dir, 'hab_population_value.TOKEN')
            task_graph.add_task(
                func=calculate_habitat_population_value,
                args=(
                    target_cv_vector_path,
                    [(ls_population_raster_path, 'total_pop'),
                     (poor_population_raster_path, 'poor_pop')],
                    global_dem_raster_path, FINAL_HAB_FIELDS,
                    habitat_raster_risk_dist_map, local_habitat_value_dir,
                    habitat_pop_value_token_path),
                target_path_list=[habitat_pop_value_token_path],
                dependent_task_list=[calculate_cv_vector_task],
                task_name=(
                    'calculate habitat population for %s' %
                    landcover_basename))

            local_lulc_raster_path = os.path.join(
                ECOSHARD_DIR, os.path.basename(landcover_url))
            LOGGER.info('starting hab value calc')

            habitat_value_token_path = os.path.join(
                local_habitat_value_dir, 'hab_value.TOKEN')
            task_graph.add_task(
                func=calculate_habitat_value,
                args=(
                    target_cv_vector_path, local_lulc_raster_path,
                    FINAL_HAB_FIELDS, habitat_raster_risk_dist_map,
                    local_habitat_value_dir, habitat_value_token_path),
                dependent_task_list=[calculate_cv_vector_task],
                target_path_list=[habitat_value_token_path],
                task_name=(
                    'calculate habitat value for %s' % landcover_basename))

    except Exception:
        LOGGER.exception('error in main')

    task_graph.join()
    task_graph.close()
    LOGGER.info('completed successfully')
