import os
import pygeoprocessing
from osgeo import ogr, gdal
from pathlib import Path
import pandas as pd




# Function to calculate pollination-dependent production
def calculate_pollination_dependent_production(poll_suff, crop_production, dependence_ratio):
    return poll_suff * crop_production * dependence_ratio

# Function to perform zonal statistics and save results
def calculate_zonal_statistics(raster_input_path, zone_input_path, output_table):
    with open(output_table, 'w', newline='') as results_file:
        results_file.write('Zone_ID,Count,Sum\n')

        # Run zonal statistics to calculate the sum within each zone
        zonal_stats = pygeoprocessing.zonal_statistics((raster_input_path, 1), zone_input_path, polygons_might_overlap=False)

        # Open the shapefile and access the layer
        # Get Gpkg driver
        driver = ogr.GetDriverByName('GPKG')
        # driver = ogr.GetDriverByName('ESRI Shapefile')
        target_vector = driver.Open(zone_input_path, 0)
        target_layer = target_vector.GetLayer()
        ID_field = 'ee_r264_name'
        # Iterate through each feature in the shapefile
        for feature in target_layer:
            feature_FID = feature.GetFID()
            feature_ID = feature.GetField(ID_field)

            # Retrieve zonal statistics for the current feature
            feature_stats = zonal_stats.get(feature_FID, {})
            feature_count = float(feature_stats.get('count', 0))
            feature_sum = float(feature_stats.get('sum', 0))

            # Write the results to the output file
            results_file.write(f'"{feature_ID}",{feature_count},{feature_sum}\n')

        # Close the shapefile to release resources
        target_vector = None
        target_layer = None

    print(f"Finished writing zonal stats for raster input: {raster_input_path}")
