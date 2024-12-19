"""Scrub geomorphology to be just the risk types."""
from osgeo import ogr
from osgeo import gdal

SEDTYPE_TO_RISK = {
    0: 5,  # unknown
    1: 5,  # sandy
    2: 1,  # unerodable
    3: 5,  # muddy
    4: 5,  # coral/mangrove
}

if __name__ == '__main__':
    base_geomorphology_path = r"C:\Users\richp\Documents\code_repos\cnc\cnc_global_cv\global_cv_workspace\ecoshard\Sediments.shp"
    target_geomorphology_path = 'geomorphology.gpkg'

    gpkg_driver = ogr.GetDriverByName('GPKG')

    base_geomorphology_vector = gdal.OpenEx(
        base_geomorphology_path, gdal.OF_VECTOR)
    base_geomorphology_layer = base_geomorphology_vector.GetLayer()
    target_geomorphology_vector = gpkg_driver.CreateDataSource(
        target_geomorphology_path)
    target_geomorphology_layer = (
        target_geomorphology_vector.CreateLayer(
            'geomorphology', base_geomorphology_layer.GetSpatialRef(),
            ogr.wkbLineString))
    target_geomorphology_layer.CreateField(ogr.FieldDefn(
        'Rgeo', ogr.OFTReal))
    target_geomorphology_defn = target_geomorphology_layer.GetLayerDefn()
    target_geomorphology_layer.StartTransaction()
    for geomorphology_feature in base_geomorphology_layer:
        target_feature = ogr.Feature(target_geomorphology_defn)
        feature_risk = SEDTYPE_TO_RISK[
            int(geomorphology_feature.GetField('SEDTYPE'))]
        if feature_risk == 0:
            continue
        geometry = geomorphology_feature.GetGeometryRef()
        target_feature.SetGeometry(geometry.Clone())
        target_feature.SetField('Rgeo', feature_risk)
        target_geomorphology_layer.CreateFeature(target_feature)
    target_geomorphology_layer.CommitTransaction()
