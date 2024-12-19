"""Create a global shore grid layer."""
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import shapely.geometry
import shapely.ops
import shapely.prepared
import shapely.wkb

SHORELINE_VECTOR_PATH = r"C:\Users\richp\Documents\code_repos\cnc\cnc_global_cv\global_cv_workspace\ecoshard\Sediments.shp"
TARGET_GRID_VECTOR_PATH = 'shore_grid.gpkg'


if __name__ == '__main__':
    shoreline_vector = gdal.OpenEx(SHORELINE_VECTOR_PATH, gdal.OF_VECTOR)
    shoreline_layer = shoreline_vector.GetLayer()
    shoreline_geom_list = []
    for shoreline_feature in shoreline_layer:
        shoreline_geom = shapely.wkb.loads(
            shoreline_feature.GetGeometryRef().ExportToWkb())
        shoreline_geom_list.append(shoreline_geom)
    shoreline_union = shapely.ops.cascaded_union(shoreline_geom_list)
    shoreline_prep = shapely.prepared.prep(shoreline_union)

    gpkg_driver = ogr.GetDriverByName('GPKG')
    wgs84_srs = osr.SpatialReference()
    wgs84_srs.ImportFromEPSG(4326)
    shore_grid_vector = gpkg_driver.CreateDataSource(TARGET_GRID_VECTOR_PATH)
    shore_grid_layer = shore_grid_vector.CreateLayer(
        'shore_grid', wgs84_srs, ogr.wkbPolygon)
    shore_grid_layer_defn = shore_grid_layer.GetLayerDefn()
    for lat in range(-80, 81):
        print(lat)
        for lng in range(-180, 180):
            degree_box = shapely.geometry.box(lng, lat, lng+1, lat+1)
            if shoreline_prep.intersects(degree_box):
                degree_geom = ogr.CreateGeometryFromWkb(degree_box.wkb)
                degree_feature = ogr.Feature(shore_grid_layer_defn)
                degree_feature.SetGeometry(degree_geom)
                shore_grid_layer.CreateFeature(degree_feature)
