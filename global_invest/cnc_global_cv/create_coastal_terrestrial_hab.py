"""Buffer CV habitat."""
import collections
import logging
import os
import sys

import ecoshard
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import pygeoprocessing
import taskgraph

LULC_RASTER_URL = 'https://storage.googleapis.com/ipbes-ndr-ecoshard-data/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7_md5_1254d25f937e6d9bdee5779d377c5aa4.tif'
BUFFER_VECTOR_URL = 'https://storage.googleapis.com/ecoshard-root/working-shards/buffered_global_shore_5km_md5_a68e1049c1c03673add014cd29b7b368.gpkg'


# I got this from https://docs.google.com/spreadsheets/d/1pYNWwPBqYYZ4tdJC3za\
# AZ8Z-CMOU2bZnB-ZvhKzDQlU/edit#gid=0

LULC_CODE_TO_HAB_MAP = {
    0: (0, None),
    10: (0, None),
    11: (0, None),
    12: (0, None),
    20: (0, None),
    30: (0, None),
    40: (0, None),
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

WORKSPACE_DIR = 'buffer_habitat_workspace'
ECOSHARD_DIR = os.path.join(WORKSPACE_DIR, 'ecoshard')
CHURN_DIR = os.path.join(WORKSPACE_DIR, 'churn')
LOGGING_LEVEL = logging.DEBUG

logging.basicConfig(
    level=LOGGING_LEVEL,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)


def main():
    """Entry point."""
    for dir_path in [ECOSHARD_DIR, CHURN_DIR]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    # make a consistent set of habitat types
    risk_distance_to_lulc_code = collections.defaultdict(list)
    for lulc_code, risk_distance in LULC_CODE_TO_HAB_MAP.items():
        risk_distance_to_lulc_code[risk_distance].append(lulc_code)

    # this maps all the same type of codes together
    current_code = -1
    lulc_code_to_reclass_value = {}
    hab_code_to_dist_risk = {}
    for risk_distance_tuple, lulc_code_list in sorted(
            risk_distance_to_lulc_code.items()):
        hab_code_to_dist_risk[current_code] = risk_distance_tuple
        for lulc_code in lulc_code_list:
            lulc_code_to_reclass_value[lulc_code] = current_code
        current_code += 1
    LOGGER.debug(lulc_code_to_reclass_value)

    task_graph = taskgraph.TaskGraph(WORKSPACE_DIR, -1, 5)

    buffer_vector_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(BUFFER_VECTOR_URL))
    download_buffer_task = task_graph.add_task(
        func=ecoshard.download_url,
        args=(BUFFER_VECTOR_URL, buffer_vector_path),
        target_path_list=[buffer_vector_path],
        task_name='download global_vector')

    lulc_raster_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(LULC_RASTER_URL))
    download_lulc_task = task_graph.add_task(
        func=ecoshard.download_url,
        args=(LULC_RASTER_URL, lulc_raster_path),
        target_path_list=[lulc_raster_path],
        task_name='download lulc raster')

    shore_lulc_raster_path = os.path.join(CHURN_DIR, 'shore_lulc.tif')
    mask_shore_task = task_graph.add_task(
        func=pygeoprocessing.mask_raster,
        args=((lulc_raster_path, 1), buffer_vector_path,
              shore_lulc_raster_path),
        target_path_list=[shore_lulc_raster_path],
        dependent_task_list=[download_buffer_task, download_lulc_task],
        task_name='mask shore')

    shore_hab_raster_path = os.path.join(CHURN_DIR, 'shore_hab.tif')
    reclassify_shore_task = task_graph.add_task(
        func=pygeoprocessing.reclassify_raster,
        args=(
            (shore_lulc_raster_path, 1), lulc_code_to_reclass_value,
            shore_hab_raster_path, gdal.GDT_Int16, -1),
        target_path_list=[shore_hab_raster_path],
        dependent_task_list=[mask_shore_task],
        task_name='shore habitat reclassify')

    shore_hab_vector_path = os.path.join(WORKSPACE_DIR, 'shore_hab.gpkg')
    polygonalize_task = task_graph.add_task(
        func=polygonalize_raster,
        args=(shore_hab_raster_path, shore_hab_vector_path),
        target_path_list=[shore_hab_vector_path],
        dependent_task_list=[mask_shore_task],
        task_name='polygonalize shore hab')

    shore_hab_risk_vector_path = os.path.join(
        WORKSPACE_DIR, 'shore_hab_risk_dist.gpkg')

    add_risk_and_prot_dist_task = task_graph.add_task(
        func=add_risk_dist_to_vector,
        args=(shore_hab_vector_path, hab_code_to_dist_risk,
              shore_hab_risk_vector_path),
        target_path_list=[shore_hab_risk_vector_path])

    task_graph.join()
    task_graph.close()


def add_risk_dist_to_vector(
        base_vector_path, base_id_to_dist_risk_map, target_vector_path):
    """Make a new vector with `risk` and `dist` fields.

    Parameters:
        base_vector_path (str): path to a vector that contains a field
            'pix_val'. This value is indexable into `base_id_to_disk_risk_map`.
        base_id_to_dist_risk_map (dict): maps values in `base_vector`:'pix_val'
            field to a (risk, distance) tuple to add to the target vector.
        target_vector_path (str): path to target vector that is a copy of
            the base vector with `risk` and `dist` fields added to it created
            by mapping `pix_val` through `base_id_to_disk_risk_map`.

    Returns:
        None.

    """
    gpkg_driver = ogr.GetDriverByName('GPKG')
    base_vector = ogr.Open(base_vector_path)
    target_vector = gpkg_driver.CopyDataSource(base_vector, target_vector_path)
    target_layer = target_vector.GetLayer()
    target_layer.CreateField(ogr.FieldDefn('risk', ogr.OFTReal))
    target_layer.CreateField(ogr.FieldDefn('dist', ogr.OFTReal))
    target_layer.StartTransaction()
    for feature in target_layer:
        pix_val = feature.GetField('pix_val')
        risk, dist = base_id_to_dist_risk_map[pix_val]
        feature.SetField('risk', float(risk))
        feature.SetField('dist', float(dist))
        target_layer.SetFeature(feature)
    target_layer.CommitTransaction()
    target_layer.SyncToDisk()
    target_layer = None
    target_vector = None


def polygonalize_raster(raster_path, target_vector_path):
    """Polygonalize `raster_path` to `target_vector_path` gpkg."""
    raster = gdal.OpenEx(raster_path, gdal.OF_RASTER)
    band = raster.GetRasterBand(1)
    gpkg_driver = ogr.GetDriverByName('gpkg')
    vector = gpkg_driver.CreateDataSource(target_vector_path)
    layer_name = os.path.basename(os.path.splitext(target_vector_path)[0])
    raster_srs = osr.SpatialReference(raster.GetProjection())
    target_layer = vector.CreateLayer(
        layer_name, raster_srs, ogr.wkbPolygon)
    target_layer.CreateField(ogr.FieldDefn('pix_val', ogr.OFTInteger))
    rasterize_callback = pygeoprocessing._make_logger_callback(
        "polygonalizing %.1f%% complete")
    gdal.Polygonize(band, None, target_layer, 0, callback=rasterize_callback)


def mask_raster(mask_array, value_array, value_nodata):
    """Mask the value array where mask array is 1."""
    result = value_array.copy()
    result[mask_array != 1] = value_nodata
    return result


def make_mask(base_raster_path, mask_vector_path, target_raster_path):
    """Mask vector onto target using base as the reference."""
    pygeoprocessing.new_raster_from_base(
        base_raster_path, target_raster_path, gdal.GDT_Byte, [0],
        fill_value_list=[0])
    pygeoprocessing.rasterize(
        mask_vector_path, target_raster_path, burn_values=[1])


if __name__ == '__main__':
    main()
