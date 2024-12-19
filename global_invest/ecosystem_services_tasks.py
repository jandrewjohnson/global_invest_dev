import os, sys
import hazelbean as hb
import seals.seals_utils as seals_utils
import pandas as pd
import global_invest
import numpy as np
from global_invest import ecosystem_services_functions
import dask
import logging
import sys

import natcap.invest.carbon
import natcap.invest.utils
# from seals import seals_utils

def project_aoi(p):
    
    p.fine_resolution_degrees = hb.get_cell_size_from_path(p.base_year_lulc_path)
    p.coarse_resolution_degrees = hb.get_cell_size_from_path(p.region_ids_coarse_path)
    
    p.fine_resolution_arcseconds = hb.pyramid_compatible_resolution_to_arcseconds[p.fine_resolution_degrees]
    p.coarse_resolution_arcseconds = hb.pyramid_compatible_resolution_to_arcseconds[p.coarse_resolution_degrees]
    
    p.processing_resolution_arcseconds = p.coarse_resolution_arcseconds

    p.ha_per_cell_coarse_path = p.get_path(hb.ha_per_cell_ref_paths[p.coarse_resolution_arcseconds])
    p.ha_per_cell_fine_path = p.get_path(hb.ha_per_cell_ref_paths[p.fine_resolution_arcseconds])
    
    # Process p.aoi to set the regional_vector, bb, bb_exact, and aoi_ha_per_cell_paths
    if isinstance(p.aoi, str):
        if p.aoi == 'global':
            p.aoi_path = p.global_regions_vector_path
            p.aoi_label = 'global'
            p.bb_exact = hb.global_bounding_box
            p.bb = p.bb_exact

            p.aoi_ha_per_cell_coarse_path = p.ha_per_cell_coarse_path
            p.aoi_ha_per_cell_fine_path = p.ha_per_cell_fine_path
        
        elif isinstance(p.aoi, str):
            if len(p.aoi) == 3: # Then it might be an ISO3 code. For now, assume so.
                p.aoi_path = os.path.join(p.cur_dir, 'aoi_' + str(p.aoi) + '.gpkg')
                p.aoi_label = p.aoi
            else: # Then it's a path to a shapefile.
                p.aoi_path = p.aoi
                p.aoi_label = os.path.splitext(os.path.basename(p.aoi))[0]

            for current_aoi_path in hb.list_filtered_paths_nonrecursively(p.cur_dir, include_strings='aoi'):
                if current_aoi_path != p.aoi_path:
                    raise NameError('There is more than one AOI in the current directory. This means you are trying to run a project in a new area of interst in a project that was already run in a different area of interest. This is not allowed! You probably want to create a new project directory and set the p = hb.ProjectFlow(...) line to point to the new directory.')

            if not hb.path_exists(p.aoi_path):
                hb.extract_features_in_shapefile_by_attribute(p.global_regions_vector_path, p.aoi_path, 'eemarine_r566_label', p.aoi.upper())
            
            p.bb_exact = hb.spatial_projection.get_bounding_box(p.aoi_path)
            p.bb = hb.pyramids.get_pyramid_compatible_bb_from_vector_and_resolution(p.aoi_path, p.processing_resolution_arcseconds)
            
            
            # Create a PROJECT-SPECIFIC version of these clipped ones.
            p.aoi_ha_per_cell_fine_path = os.path.join(p.cur_dir, 'pyramids', 'aoi_ha_per_cell_fine.tif')
            if not hb.path_exists(p.aoi_ha_per_cell_fine_path):
                hb.create_directories(p.aoi_ha_per_cell_fine_path)
                
                #  make ha_per_cell_paths not be a dict but a project level ha_per_cell_fine_path etc
                hb.clip_raster_by_bb(p.ha_per_cell_fine_path, p.bb, p.aoi_ha_per_cell_fine_path)
            
            p.aoi_ha_per_cell_coarse_path = os.path.join(p.cur_dir, 'pyramids', 'aoi_ha_per_cell_coarse.tif')
            if not hb.path_exists(p.aoi_ha_per_cell_coarse_path):
                hb.create_directories(p.aoi_ha_per_cell_coarse_path)
                hb.clip_raster_by_bb(p.ha_per_cell_coarse_path, p.bb, p.aoi_ha_per_cell_coarse_path)
        
            
            
        else:
            p.bb_exact = hb.spatial_projection.get_bounding_box(p.aoi_path)
            p.bb = hb.pyramids.get_pyramid_compatible_bb_from_vector_and_resolution(p.aoi_path, p.processing_resolution_arcseconds)

            # Create a PROJECT-SPECIFIC version of these clipped ones.
            p.aoi_ha_per_cell_fine_path = os.path.join(p.cur_dir, 'pyramids', 'aoi_ha_per_cell_fine.tif')
            if not hb.path_exists(p.aoi_ha_per_cell_fine_path):
                hb.create_directories(p.aoi_ha_per_cell_fine_path)
                hb.clip_raster_by_bb(p.ha_per_cell_paths[p.fine_resolution_arcseconds], p.bb, p.aoi_ha_per_cell_fine_path)
            
            p.aoi_ha_per_cell_coarse_path = os.path.join(p.cur_dir, 'pyramids', 'aoi_ha_per_cell_coarse.tif')
            if not hb.path_exists(p.aoi_ha_per_cell_coarse_path):
                hb.create_directories(p.aoi_ha_per_cell_coarse_path)
                hb.clip_raster_by_bb(p.ha_per_cell_paths[p.coarse_resolution_arcseconds], p.bb, p.aoi_ha_per_cell_coarse_path)
                    
    else:
        raise NameError('Unable to interpret p.aoi.')

    
def ecosystem_services(p):
    pass # Just to generate a folder


def carbon_storage_biophysical_invest(p):
    
    # Based on the RAW InVEST launch script generated by InVEST user interface
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter(fmt=natcap.invest.utils.LOG_FMT, datefmt='%m/%d/%Y %H:%M:%S ')
    handler.setFormatter(formatter)
    logging.basicConfig(level=logging.INFO, handlers=[handler])

    # Identify customized input paths based on ProjectFlow attributes
    p.example_lulc_input_path = p.get_path("invest_sample_data", "Carbon", "lulc_current_willamette.tif")
    p.example_biophysical_table_input_path = p.get_path("invest_sample_data", "Carbon", "carbon_pools_willamette.csv")

    # Plug attributes into the InVEST args dict (model spec)
    args = {
        'calc_sequestration': False,
        'carbon_pools_path': p.example_biophysical_table_input_path,
        'discount_rate': '',
        'do_redd': False,
        'do_valuation': False,
        'lulc_cur_path': p.example_lulc_input_path,
        'lulc_cur_year': '',
        'lulc_fut_path': '',
        'lulc_fut_year': '',
        'lulc_redd_path': '',
        'n_workers': '-1',
        'price_per_metric_ton_of_c': '',
        'rate_change': '',
        'results_suffix': '',
        'workspace_dir': p.cur_dir,
    }

    # Call the invest carbon model's execute function with the args dict. 
    # This function could be replaced by whatever function you create for other global ES or GEP models.
    natcap.invest.carbon.execute(args)
    


def aoi_inputs(p):
    p.carbon_zones_path = p.get_path("global_invest", "carbon", 'carbon_zones_rasterized.tif')
    if p.aoi != 'global':
        p.aoi_carbon_zones_path = os.path.join(p.cur_dir, 'carbon_zones.tif')                    
        if not hb.path_exists(p.aoi_carbon_zones_path):
            hb.clip_raster_by_bb(p.carbon_zones_path, p.bb, p.aoi_carbon_zones_path)
        p.aoi_base_year_lulc_path = os.path.join(p.cur_dir, 'lulc.tif')                    
        if not hb.path_exists(p.aoi_base_year_lulc_path):
            hb.clip_raster_by_bb(p.base_year_lulc_path, p.bb, p.aoi_base_year_lulc_path)
    else:
        # Then it's global so create a shortcut to the base_data_dir
        hb.create_shortcut('base_data_shortcut', p.base_data_dir)
        p.aoi_carbon_zones_path = p.carbon_zones_path    
        p.aoi_base_year_lulc_path = p.base_year_lulc_path    
    
def carbon_storage_simple(p):
    """Run a single version of calculating carbon storage present from the project's base_year LULC maps."""

    # Input Paths
    p.exhaustive_carbon_table_path = p.get_path("global_invest", "carbon", "exhaustive_carbon_table.csv")


    # Calculate zonal statistics
    p.carbon_storage_Mg_per_ha_path = os.path.join(p.cur_dir, 'carbon_storage_Mg_per_ha.csv')


    carbon_Mg_per_ha_output_path = os.path.join(p.cur_dir, 'carbon_Mg_per_ha.tif')
    if not hb.path_exists(carbon_Mg_per_ha_output_path):
        ecosystem_services_functions.carbon_storage_ipcc_tier_1(p.aoi_base_year_lulc_path, p.aoi_carbon_zones_path, p.exhaustive_carbon_table_path, carbon_Mg_per_ha_output_path)

    # Calculate zonal statistics
    p.carbon_storage_csv_path = os.path.join(p.cur_dir, 'carbon_storage.csv')
    if not hb.path_exists(p.carbon_storage_csv_path):
        5

def carbon_storage_biophysical_for_seals_scenarios(p):
    """Iterate over a scenarios file to calculate carbon storage presentfrom LULC maps."""

    p.exhaustive_carbon_table_path = os.path.join(p.base_data_dir, "global_invest", "carbon", "exhaustive_carbon_table.csv")

    p.carbon_zones_path = os.path.join(p.base_data_dir, "global_invest", "carbon", 'carbon_zones_rasterized.tif')

    p.joined_carbon_table_stacked_with_seals_simplified_classes_path = os.path.join(p.cur_dir, 'joined_carbon_table_stacked_with_seals_simplified_classes.csv')

    p.carbon_storage_csv_path = os.path.join(p.cur_dir, 'carbon_storage.csv')

    # TODOO I may want to redo the run_subset_tiles functionality to instead just require a small aoi 
    # polygon (which would in principle only cover 1 tile. This is because the definition of BB gets confused
    # but also i would need to specify the difference between bb_pyramid vs bb_vector. Think about this.
    # Also note that the bb_pyramid would be determined by the processing_resolution

    # if not hasattr(p, 'scenarios_df'):
    #     p.scenarios_df = hb.get_dummy_scenarios_df()
        
    if p.run_this:      
        csvs_to_merge = []
        gpkgs_to_merge = []

        for index, row in p.scenarios_df.iterrows():

            # NOTE THAT WE CALL THE SEALS_UTILS version of this function and not the hazelbean function
            # because we are assuming the structure of a seals scenarios.csv
            seals_utils.assign_df_row_to_object_attributes(p, row)
            hb.log('Calculating carbon storage on ' + str(index) + ' of ' + str(len(p.scenarios_df)))

            if p.scenario_type == 'baseline':
                for year in p.base_years:
                    current_lulc_path = os.path.join(p.stitched_lulc_simplified_scenarios_dir  , 'lulc_' + p.lulc_src_label + '_' + p.lulc_simplification_label + '_' + p.model_label + '_' + str(year) + '.tif')
                    current_lulc_bb = hb.get_bounding_box(current_lulc_path)
                    carbon_Mg_per_ha_output_path = os.path.join(p.cur_dir, 'carbon_Mg_per_ha_' + p.model_label + '_' + str(year) + '.tif')
                    if not hb.path_exists(carbon_Mg_per_ha_output_path):
                        
                        current_carbon_zones_path = os.path.join(p.cur_dir, 'carbon_zones.tif')
                        if not hb.path_exists(current_carbon_zones_path):
                            hb.clip_raster_by_bb(p.carbon_zones_path, current_lulc_bb, current_carbon_zones_path)

                        current_ha_per_cell_path = os.path.join(p.cur_dir, 'ha_per_cell_fine.tif')
                        if not hb.path_exists(current_ha_per_cell_path):
                            hb.clip_raster_by_bb(p.aoi_ha_per_cell_fine_path, current_lulc_bb, current_ha_per_cell_path)

                        ecosystem_services_functions.carbon_storage_ipcc_tier_1(current_lulc_path, current_carbon_zones_path, p.exhaustive_carbon_table_path, carbon_Mg_per_ha_output_path)

                    vector_output_path = os.path.join(p.cur_dir, 'carbon_Mg_per_ha_' + p.model_label + '_' + str(year) + '.gpkg')
                    csv_output_path = os.path.join(p.cur_dir, 'carbon_Mg_per_ha_' + p.model_label + '_' + str(year) + '.csv')

                    csvs_to_merge.append(csv_output_path)
                    gpkgs_to_merge.append(vector_output_path)                                                  
                    
                    skip = False
                    if not skip and not hb.path_exists(vector_output_path):

                        gdf = hb.zonal_statistics(carbon_Mg_per_ha_output_path,
                                p.aoi,
                                id_column_label=None,
                                zone_ids_raster_path=None,
                                stats_to_retrieve='sums',
                                vector_columns_to_keep = 'just_id',
                                csv_output_path=csv_output_path,
                                vector_output_path=vector_output_path)


            elif p.scenario_type != 'baseline':
                for year in p.years:
                    current_lulc_path = os.path.join(p.stitched_lulc_simplified_scenarios_dir  , 'lulc_' + p.lulc_src_label + '_' + p.lulc_simplification_label + '_' + p.exogenous_label + '_' + p.climate_label + '_' + p.model_label + '_' + p.counterfactual_label + '_' + str(year) + '.tif')
                    current_lulc_bb = hb.get_bounding_box(current_lulc_path)
                    carbon_Mg_per_ha_output_path = os.path.join(p.cur_dir, 'carbon_Mg_per_ha_' + p.exogenous_label + '_' + p.climate_label + '_' + p.model_label + '_' + p.counterfactual_label + '_' + str(year) + '.tif')
                    if not hb.path_exists(carbon_Mg_per_ha_output_path):

                        current_carbon_zones_path = os.path.join(p.cur_dir, 'carbon_zones.tif')
                        if not hb.path_exists(current_carbon_zones_path):
                            hb.clip_raster_by_bb(p.carbon_zones_path, current_lulc_bb, current_carbon_zones_path)

                        current_ha_per_cell_path = os.path.join(p.cur_dir, 'ha_per_cell_fine.tif')
                        if not hb.path_exists(current_ha_per_cell_path):
                            hb.clip_raster_by_bb(p.aoi_ha_per_cell_fine_path, current_lulc_bb, current_ha_per_cell_path)


                        ecosystem_services_functions.carbon_storage_ipcc_tier_1(current_lulc_path, current_carbon_zones_path, p.exhaustive_carbon_table_path, carbon_Mg_per_ha_output_path)
                    
                    vector_output_path = os.path.join(p.cur_dir, 'carbon_Mg_per_ha_' + p.exogenous_label + '_' + p.climate_label + '_' + p.model_label + '_' + p.counterfactual_label + '_' + str(year) + '.gpkg')
                    csv_output_path = os.path.join(p.cur_dir, 'carbon_Mg_per_ha_' + p.exogenous_label + '_' + p.climate_label + '_' + p.model_label + '_' + p.counterfactual_label + '_' + str(year) + '.csv')                    
                    
                    csvs_to_merge.append(csv_output_path)
                    gpkgs_to_merge.append(vector_output_path)                                                  
                    
                    skip = False
                    if not skip and not hb.path_exists(vector_output_path):

                        gdf = hb.zonal_statistics(carbon_Mg_per_ha_output_path,
                                p.aoi,
                                id_column_label=None,
                                zone_ids_raster_path=None,
                                stats_to_retrieve='sums',
                                vector_columns_to_keep = 'just_id',
                                csv_output_path=csv_output_path,
                                vector_output_path=vector_output_path)
        skip = False
                           
        if not skip and not hb.path_exists(p.carbon_storage_csv_path):
            hb.df_merge_list_of_csv_paths(csvs_to_merge, p.carbon_storage_csv_path, on='generated_ids', column_suffix='ignore', verbose=False)


def carbon_storage_economic(p):
    """Converts carbon storage biophysical outputs into per-econ region shockfiles."""

    
    p.carbon_storage_shockfile_csv_path = os.path.join(p.cur_dir, 'carbon_storage_shockfile.csv')

    if p.run_this:              
        
        if not hb.path_exists(p.carbon_storage_shockfile_csv_path):      
            carbon_df = pd.read_csv(p.carbon_storage_csv_path)
        
            # Iterate over the non-baseline scenarios and compare them with the appropriate baseline 
            # to generate percent change.
            shock_df = pd.DataFrame(carbon_df['generated_ids'])
            for index, row in p.scenarios_df.iterrows():
                seals_utils.assign_df_row_to_object_attributes(p, row)
                hb.log('Calculating carbon storage shockfile on ' + str(index) + ' of ' + str(len(p.scenarios_df)))
        
                if p.scenario_type != 'baseline':
                    print('Analyzing', p.scenario_type, 'which has a baseline scenario of', p.baseline_reference_label)
                            
                    for year_c, year in enumerate(p.years):

                        current_scenario_label = p.exogenous_label + '_' + p.climate_label + '_' + p.model_label + '_' + p.counterfactual_label + '_' + str(year)   
                        
                        # Get the label for the previous scenario (tricky if it's the first scenario year, then we need to use the baseline)
                        if year_c == 0:
                            previous_year = p.key_base_year
                            row = p.scenarios_df.loc[p.scenarios_df['scenario_label'] == p.baseline_reference_label]
                            current_model_label = row['model_label'].values[0]
                            previous_scenario_label = current_model_label + '_' + str(previous_year)                        
                        else:
                            previous_year = p.years[year_c - 1]
                            previous_scenario_label = p.exogenous_label + '_' + p.climate_label + '_' + p.model_label + '_' + p.counterfactual_label + '_' + str(previous_year)   

                        previous_carbon_label = 'carbon_Mg_per_ha_' + previous_scenario_label + '_sums'
                        current_carbon_label = 'carbon_Mg_per_ha_' +current_scenario_label + '_sums'
                        new_label = 'carbon_shock_' + current_scenario_label


                        a = (carbon_df[current_carbon_label] - carbon_df[previous_carbon_label]) / carbon_df[previous_carbon_label]
                        carbon_df[new_label] = a
                        subset_df = carbon_df[['generated_ids', new_label]]
                        shock_df = hb.df_merge(shock_df, subset_df, on='generated_ids', verbose=False)

            shock_df.to_csv(p.carbon_storage_shockfile_csv_path, index=False)

            
       
def pollination_biophysical(p):
    """Calculate carbon storage presentfrom LULC maps."""

    p.exhaustive_carbon_table_path = os.path.join(p.base_data_dir, "global_invest", "carbon", "exhaustive_carbon_table.csv")

    p.carbon_zones_path = os.path.join(p.base_data_dir, "global_invest", "carbon", 'carbon_zones_rasterized.tif')

    p.joined_carbon_table_stacked_with_seals_simplified_classes_path = os.path.join(p.cur_dir, 'joined_carbon_table_stacked_with_seals_simplified_classes.csv')

    # TODOO I may want to redo the run_subset_tiles functionality to instead just require a small aoi 
    # polygon (which would in principle only cover 1 tile. This is because the definition of BB gets confused
    # but also i would need to specify the difference between bb_pyramid vs bb_vector. Think about this.
    # Also note that the bb_pyramid would be determined by the processing_resolution


    if p.run_this:     

        # make a blank DF to store results
        df = None

        csvs_to_merge = []
        gpkgs_to_merge = []
            
        for index, row in p.scenarios_df.iterrows():
            seals_utils.assign_df_row_to_object_attributes(p, row)
            hb.log('Calculating carbon storage on ' + str(index) + ' of ' + str(len(p.scenarios_df)))

            if p.scenario_type == 'baseline':
                for year in p.base_years:

                    current_lulc_path = os.path.join(p.stitched_lulc_simplified_scenarios_dir  , 'lulc_' + p.lulc_src_label + '_' + p.lulc_simplification_label + '_' + p.model_label + '_' + str(year) + '.tif')
                    current_lulc_bb = hb.get_bounding_box(current_lulc_path)
                    pollination_sufficiency_output_path = os.path.join(p.cur_dir, 'pollination_sufficiency_' + p.model_label + '_' + str(year) + '.tif')
                    if not hb.path_exists(pollination_sufficiency_output_path):
                        hb.log('Running global_invest_main.make_poll_suff on LULC: ' + str(current_lulc_path) + ' and saving results to ' + str(pollination_sufficiency_output_path))
                        ecosystem_services_functions.pollination_sufficiency(current_lulc_path, pollination_sufficiency_output_path)

                    # Calculate zonal statistics for pollination sufficiency
                    # df_sum_col = 'pollination_sufficiency_' + p.model_label + '_' + str(year) + '_sum'
                    # df_count_col = 'pollination_sufficiency_' + p.model_label + '_' + str(year) + '_count'

                    vector_output_path = os.path.join(p.cur_dir, hb.path_replace_extension(pollination_sufficiency_output_path, '.gpkg'))
                    csv_output_path = os.path.join(p.cur_dir, hb.path_replace_extension(pollination_sufficiency_output_path, '.csv'))

                    csvs_to_merge.append(csv_output_path)
                    gpkgs_to_merge.append(vector_output_path)                                                  
                    
                    if not hb.path_exists(vector_output_path):

                        gdf = hb.zonal_statistics(pollination_sufficiency_output_path,
                                p.aoi,
                                id_column_label=None,
                                zone_ids_raster_path=None,
                                stats_to_retrieve='sums_counts',
                                csv_output_path=csv_output_path,
                                vector_output_path=vector_output_path)

            
            elif p.scenario_type != 'baseline':
                for year in p.years:
                    current_lulc_path = os.path.join(p.stitched_lulc_simplified_scenarios_dir  , 'lulc_' + p.lulc_src_label + '_' + p.lulc_simplification_label + '_' + p.exogenous_label + '_' + p.climate_label + '_' + p.model_label + '_' + p.counterfactual_label + '_' + str(year) + '.tif')
                    current_lulc_bb = hb.get_bounding_box(current_lulc_path)
                    pollination_sufficiency_output_path = os.path.join(p.cur_dir, 'pollination_sufficiency_' + p.exogenous_label + '_' + p.climate_label + '_' + p.model_label + '_' + p.counterfactual_label + '_' + str(year) + '.tif')
                    if not hb.path_exists(pollination_sufficiency_output_path):
                        hb.log('Running global_invest_main.make_poll_suff on LULC: ' + str(current_lulc_path) + ' and saving results to ' + str(pollination_sufficiency_output_path))
                        ecosystem_services_functions.pollination_sufficiency(current_lulc_path, pollination_sufficiency_output_path)

                    vector_output_path = os.path.join(p.cur_dir, hb.path_replace_extension(pollination_sufficiency_output_path, '.gpkg'))
                    csv_output_path = os.path.join(p.cur_dir, hb.path_replace_extension(pollination_sufficiency_output_path, '.csv'))

                    csvs_to_merge.append(csv_output_path)
                    gpkgs_to_merge.append(vector_output_path)                                                  
                    
                    if not hb.path_exists(vector_output_path):

                        gdf = hb.zonal_statistics(pollination_sufficiency_output_path,
                                p.aoi,
                                id_column_label=None,
                                zone_ids_raster_path=None,
                                stats_to_retrieve='sums_counts',
                                csv_output_path=csv_output_path,
                                vector_output_path=vector_output_path)
        
        # Merge all the csvs and gpkg files#

        # NOTE This currently doesn't make sense to sum up until after i run the pollination_economic function
        skip = True
        if not skip:
            merged_csv_output_path = os.path.join(p.cur_dir, 'pollination_sufficiency.csv')
            if not hb.path_exists(merged_csv_output_path):
                hb.df_merge_list_of_csv_paths(csvs_to_merge, merged_csv_output_path, on='generated_ids', column_suffix='ignore', verbose=False)

def pollination_economic(p):
    """Convert pollination into shockfile."""

    ### Thought: think again about adding a penalty to only crediting pollination losses in bau
    p.pollination_shock_change_per_region_path = os.path.join(p.cur_dir, 'pollination_shock_change_per_region.gpkg')
    p.pollination_shock_csv_path = os.path.join(p.cur_dir, 'pollination_shock.csv')
    p.crop_data_dir = os.path.join(p.base_data_dir, "crops\earthstat\crop_production")
    
    # DATA NOT HERE YET TODO
    p.crop_prices_dir = os.path.join(p.base_data_dir, "crops\earthstat\crop_prices")
    # p.crop_prices_dir = r"C:\Users\jajohns\Files\Research\base_data\pyramids\crops\price"
    # p.pollination_biophysical_dir = r"C:\Users\jajohns\Files\Research\cge\gtap_invest\projects\feedback_with_policies\intermediate\pollination_biophysical"
    p.pollination_dependence_spreadsheet_input_path = os.path.join(p.base_data_dir, "pollination/rspb20141799supp3.xls") # Note had to fix pol.dep for cofee and greenbroadbean as it was 25 not .25


    if p.run_this:
        df = None

        for index, row in p.scenarios_df.iterrows():
            seals_utils.assign_df_row_to_object_attributes(p, row)
            hb.log('Calculating pollination_economic on ' + str(index) + ' of ' + str(len(p.scenarios_df)))

            if p.scenario_type == 'baseline':
                for year in p.base_years:

                    current_pollination_sufficiency_path = os.path.join(p.pollination_dir, 'pollination_sufficiency_' + p.model_label + '_' + str(year) + '.tif')
                    if not hb.path_exists(current_pollination_sufficiency_path):
                        raise NameError('Pollination sufficiency raster not found at ' + str(current_pollination_sufficiency_path))

                    crop_value_baseline_output_path = os.path.join(p.cur_dir, 'crop_value_' + p.model_label + '_' + str(year) + '.tif')

                    # TODOOO: Figure out how to properly load baseline and scenario to calculate shock from previous period.

                    
                    # crop_value_baseline_output_path = os.path.join(p.cur_dir, 'crop_value_' + p.model_label + '_' + str(year) + '.tif')
                    # crop_value_baseline_output_path = os.path.join(p.cur_dir, 'crop_value_' + p.model_label + '_' + str(year) + '.tif')
                    # crop_value_baseline_output_path = os.path.join(p.cur_dir, 'crop_value_' + p.model_label + '_' + str(year) + '.tif')
                    # crop_value_baseline_output_path = os.path.join(p.cur_dir, 'crop_value_' + p.model_label + '_' + str(year) + '.tif')
                    
                    p.crop_value_no_pollination_path = os.path.join(p.cur_dir, 'crop_value_no_pollination.tif')
                    p.crop_value_max_lost_path = os.path.join(p.cur_dir, 'crop_value_max_lost.tif')
                    p.crop_value_max_lost_10s_path = os.path.join(p.cur_dir, 'crop_value_max_lost_10s.tif')
                    p.crop_value_baseline_10s_path = os.path.join(p.cur_dir, 'crop_value_baseline_10s.tif')


            elif p.scenario_type != 'baseline':
                for year in p.years:
                    current_lulc_path = os.path.join(p.stitched_lulc_simplified_scenarios_dir  , 'lulc_' + p.lulc_src_label + '_' + p.lulc_simplification_label + '_' + p.exogenous_label + '_' + p.climate_label + '_' + p.model_label + '_' + p.counterfactual_label + '_' + str(year) + '.tif')

        ###########################################
        ###### Calculate base-data necessary to do conversion of biophysical to shockfile
        ###########################################

        if not all([hb.path_exists(i) for i in [p.crop_value_baseline_path,
                                                p.crop_value_no_pollination_path,
                                                p.crop_value_max_lost_path,]]):
            df = pd.read_excel(p.pollination_dependence_spreadsheet_input_path, sheet_name='Crop nutrient content')

            crop_names = list(df['Crop map file name'])[:-3] # Drop last three which were custom addons in manuscript and don't seem to have earthstat data for.
            pollination_dependence = list(df['poll.dep'])
            crop_value_baseline = np.zeros(hb.get_shape_from_dataset_path(p.ha_per_cell_300sec_path))
            crop_value_no_pollination = np.zeros(hb.get_shape_from_dataset_path(p.ha_per_cell_300sec_path))
            for c, crop_name in enumerate(crop_names):
                L.info('Calculating value yield effect from pollination for ' + str(crop_name) + ' with pollination dependence ' + str(pollination_dependence[c]))
                crop_price_path = os.path.join(p.crop_prices_dir, crop_name + '_prices_per_ton.tif')
                crop_price = hb.as_array(crop_price_path)
                crop_price = np.where(crop_price > 0, crop_price, 0.0)
                crop_yield = hb.as_array(os.path.join(p.crop_data_dir, crop_name + '_HarvAreaYield_Geotiff', crop_name + '_Production.tif'))
                crop_yield = np.where(crop_yield > 0, crop_yield, 0.0)

                crop_value_baseline += (crop_yield * crop_price)
                crop_value_no_pollination += (crop_yield * crop_price) * (1 - float(pollination_dependence[c]))

            crop_value_max_lost = crop_value_baseline - crop_value_no_pollination
            #
            # crop_value_baseline_path = os.path.join(p.cur_dir, 'crop_value_baseline.tif')
            # crop_value_no_pollination_path = os.path.join(p.cur_dir, 'crop_value_no_pollination.tif')
            # crop_value_max_lost_path = os.path.join(p.cur_dir, 'crop_value_max_lost.tif')

            hb.save_array_as_geotiff(crop_value_baseline, p.crop_value_baseline_path, p.match_300sec_path, ndv=-9999, data_type=6)
            hb.save_array_as_geotiff(crop_value_no_pollination, p.crop_value_no_pollination_path, p.match_300sec_path, ndv=-9999, data_type=6)
            hb.save_array_as_geotiff(crop_value_max_lost, p.crop_value_max_lost_path, p.match_300sec_path, ndv=-9999, data_type=6)


        ### Resample the base data to match LULC
        global_bb = hb.get_bounding_box(p.base_year_lulc_path)
        stitched_bb = hb.get_bounding_box(p.baseline_clipped_lulc_path)
        if stitched_bb != global_bb:
            current_path = os.path.join(p.cur_dir, 'crop_value_max_lost_clipped.tif')
            hb.clip_raster_by_bb(p.crop_value_max_lost_path, stitched_bb, current_path)
            p.crop_value_max_lost_path = current_path

        if not hb.path_exists(p.crop_value_baseline_10s_path):
            hb.resample_to_match(p.crop_value_baseline_path, p.baseline_clipped_lulc_path, p.crop_value_baseline_10s_path, ndv=-9999., output_data_type=6)


        if not hb.path_exists(p.crop_value_max_lost_10s_path):
            hb.resample_to_match(p.crop_value_max_lost_path, p.baseline_clipped_lulc_path, p.crop_value_max_lost_10s_path, ndv=-9999., output_data_type=6)


        ###########################################
        ###### Calculate crop_value_pollinator_adjusted.
        ###########################################

        # Incorporate the "sufficient pollination threshold" of 30%
        # TODOO Go through and systematically pull into config files to initialize model and write output summary of what were used.
        sufficient_pollination_threshold = 0.3

        ### BASELINE crop_value_pollinator_adjusted:
        policy_scenario_label = 'gtap1_baseline_' + str(p.base_year)
        current_output_excel_path = os.path.join(p.cur_dir, 'crop_value_pollinator_adjusted_' + policy_scenario_label + '_zonal_stats.xlsx')
        suff_path = os.path.join(p.pollination_biophysical_dir, 'poll_suff_ag_coverage_prop_lulc_esa_' + policy_scenario_label + '.tif')
        crop_value_pollinator_adjusted_path = os.path.join(p.cur_dir, 'crop_value_pollinator_adjusted_' + policy_scenario_label + '.tif')

        if not hb.path_exists(crop_value_pollinator_adjusted_path):
            hb.raster_calculator_af_flex([p.crop_value_baseline_10s_path, p.crop_value_max_lost_10s_path, suff_path, p.base_year_simplified_lulc_path], lambda baseline_value, max_loss, suff, lulc:
                    np.where((max_loss > 0) & (suff < sufficient_pollination_threshold) & (lulc == 2), baseline_value - max_loss * (1 - (1/sufficient_pollination_threshold) * suff),
                        np.where((max_loss > 0) & (suff >= sufficient_pollination_threshold) & (lulc == 2), baseline_value, -9999.)), output_path=crop_value_pollinator_adjusted_path)

        # Do zonal statistics on outputed raster by AEZ-REG. Note that we need sum and count for when/if we calculate mean ON GRIDCELLS WITH AG.
        if not hb.path_exists(current_output_excel_path):
            df = hb.zonal_statistics_flex(crop_value_pollinator_adjusted_path,
                                          p.gtap37_aez18_path,
                                          zone_ids_raster_path=p.zone_ids_raster_path,
                                          id_column_label='pyramid_id',
                                          zones_raster_data_type=5,
                                          values_raster_data_type=6,
                                          zones_ndv=-9999,
                                          values_ndv=-9999,
                                          all_touched=None,
                                          stats_to_retrieve='sums_counts',
                                          assert_projections_same=False, )
            generated_scenario_label = 'gtap1_baseline_2014'
            df.rename(columns={'sums': generated_scenario_label + '_sum', 'counts': generated_scenario_label + '_count'}, inplace=True)
            df.to_excel(current_output_excel_path)
        else:
            generated_scenario_label = 'gtap1_baseline_2014'
            df = pd.read_excel(current_output_excel_path, index_col=0)
            df.rename(columns={'sums': generated_scenario_label + '_sum', 'counts': generated_scenario_label + '_count'}, inplace=True)
        merged_df = df

        ### SCENARIO crop_value_pollinator_adjusted
        for luh_scenario_label in p.luh_scenario_labels:
            for scenario_year in p.scenario_years:
                for policy_scenario_label in p.policy_scenario_labels:
                    current_output_excel_path = os.path.join(p.cur_dir, 'crop_value_pollinator_adjusted_' + policy_scenario_label + '_zonal_stats.xlsx')
                    suff_path = os.path.join(p.pollination_biophysical_dir, 'poll_suff_ag_coverage_prop_gtap1_' + luh_scenario_label
                                                + '_' + str(scenario_year) + '_' + policy_scenario_label + '.tif')
                    lulc_path = os.path.join(p.stitched_lulc_esa_scenarios_dir, 'lulc_seals7_gtap1_' + luh_scenario_label
                                                + '_' + str(scenario_year) + '_' + policy_scenario_label + '.tif')

                    crop_value_pollinator_adjusted_path = os.path.join(p.cur_dir, 'crop_value_pollinator_adjusted_' + policy_scenario_label + '.tif')

                    if not hb.path_exists(crop_value_pollinator_adjusted_path):
                        hb.raster_calculator_af_flex([p.crop_value_baseline_10s_path, p.crop_value_max_lost_10s_path, suff_path, lulc_path], lambda baseline_value, max_loss, suff, lulc:
                                np.where((max_loss > 0) & (suff < sufficient_pollination_threshold) & (lulc == 2), baseline_value - max_loss * (1 - (1/sufficient_pollination_threshold) * suff),
                                        np.where((max_loss > 0) & (suff >= sufficient_pollination_threshold) & (lulc == 2), baseline_value, -9999.)), output_path=crop_value_pollinator_adjusted_path)


                        # TODOOO: Continue thinking about what the right shock is overall. Is it the average on NEW land? Or the aggregate value
                        # To isolate the effect, maybe calculate the average value of crop loss on cells that are cultivated in both scenarios? Start on a dask function that does that?

                    if not hb.path_exists(current_output_excel_path):
                        df = hb.zonal_statistics_flex(crop_value_pollinator_adjusted_path,
                                                      p.gtap37_aez18_path,
                                                      zone_ids_raster_path=p.zone_ids_raster_path,
                                                      id_column_label='pyramid_id',
                                                      zones_raster_data_type=5,
                                                      values_raster_data_type=6,
                                                      zones_ndv=-9999,
                                                      values_ndv=-9999,
                                                      all_touched=None,
                                                      stats_to_retrieve='sums_counts',
                                                      assert_projections_same=False, )

                        generated_scenario_label = 'gtap2_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label
                        df.rename(columns={'sums': generated_scenario_label + '_sum', 'counts': generated_scenario_label + '_count'}, inplace=True)
                        df.to_excel(current_output_excel_path)
                    else:
                        generated_scenario_label = 'gtap2_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label
                        df = pd.read_excel(current_output_excel_path, index_col=0)
                        df.rename(columns={'sums': generated_scenario_label + '_sum', 'counts': generated_scenario_label + '_count', generated_scenario_label + '_total': generated_scenario_label + '_sum',}, inplace=True)
                    merged_df = pd.merge(merged_df, df, how='outer', left_index=True, right_index=True)

        ###########################################
        ###### Calculate change from scenario to baseline, on and not on existing ag
        ###########################################

        baseline_policy_scenario_label = 'gtap1_baseline_' + str(p.base_year)
        baseline_crop_value_pollinator_adjusted_path = os.path.join(p.cur_dir, 'crop_value_pollinator_adjusted_' + baseline_policy_scenario_label + '.tif')
        for luh_scenario_label in p.luh_scenario_labels:
            for scenario_year in p.scenario_years:
                for policy_scenario_label in p.policy_scenario_labels:

                    # Calculate difference between scenario and BASELINE for crop value adjusted
                    bau_crop_value_pollinator_adjusted_path = os.path.join(p.cur_dir, 'crop_value_pollinator_adjusted_BAU.tif')
                    current_crop_value_pollinator_adjusted_path = os.path.join(p.cur_dir, 'crop_value_pollinator_adjusted_' + policy_scenario_label + '.tif')

                    crop_value_difference_from_baseline_path = os.path.join(p.cur_dir, 'crop_value_difference_from_baseline_' + policy_scenario_label + '.tif')
                    if not hb.path_exists(crop_value_difference_from_baseline_path):
                        hb.dask_compute([baseline_crop_value_pollinator_adjusted_path, current_crop_value_pollinator_adjusted_path], lambda x, y: y - x, crop_value_difference_from_baseline_path)

                    # Zonal stats for difference from Baseline
                    current_output_excel_path = os.path.join(p.cur_dir, 'crop_value_difference_from_baseline_' + policy_scenario_label + '_zonal_stats.xlsx')
                    if not hb.path_exists(current_output_excel_path):
                        df = hb.zonal_statistics_flex(crop_value_difference_from_baseline_path,
                                                      p.gtap37_aez18_path,
                                                      zone_ids_raster_path=p.zone_ids_raster_path,
                                                      id_column_label='pyramid_id',
                                                      zones_raster_data_type=5,
                                                      values_raster_data_type=6,
                                                      zones_ndv=-9999,
                                                      values_ndv=-9999,
                                                      all_touched=None,
                                                      stats_to_retrieve='sums_counts',
                                                      assert_projections_same=False, )
                        generated_scenario_label = 'gtap2_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label
                        df.rename(columns={'sums': generated_scenario_label + '_sum', 'counts': generated_scenario_label + '_count'}, inplace=True)
                        df.to_excel(current_output_excel_path)

                    # Calc difference between scenario and BASELINE for crop_value on grid-cells that were agri in both lulc maps.
                    lulc_path = os.path.join(p.stitched_lulc_esa_scenarios_dir, 'lulc_seals7_gtap1_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label + '.tif')
                    crop_value_difference_from_baseline_existing_ag_path = os.path.join(p.cur_dir, 'crop_value_difference_from_baseline_existing_ag_' + policy_scenario_label + '.tif')

                    def op(x, y, w, z):
                        r = dask.array.where((w == 2) & (z == 2), y - x, 0.)
                        rr = (z * 0.0) + r # HACK. Dask.array.where was returning a standard xarray rather than a rioxarray. This dumb hack makes it inherit the rioxarray parameters of z
                        return rr

                    if not hb.path_exists(crop_value_difference_from_baseline_existing_ag_path):
                        op_paths = [
                            baseline_crop_value_pollinator_adjusted_path,
                            crop_value_pollinator_adjusted_path,
                            lulc_path,
                            p.base_year_simplified_lulc_path,
                        ]

                        hb.dask_compute(op_paths, op, crop_value_difference_from_baseline_existing_ag_path)

                    # Zonal stats for difference from Baseline
                    current_output_excel_path = os.path.join(p.cur_dir, 'crop_value_difference_from_baseline_existing_ag_' + policy_scenario_label + '_zonal_stats.xlsx')
                    if not hb.path_exists(current_output_excel_path):
                        df = hb.zonal_statistics_flex(crop_value_difference_from_baseline_existing_ag_path,
                                                      p.gtap37_aez18_path,
                                                      zone_ids_raster_path=p.zone_ids_raster_path,
                                                      id_column_label='pyramid_id',
                                                      zones_raster_data_type=5,
                                                      values_raster_data_type=6,
                                                      zones_ndv=-9999,
                                                      values_ndv=-9999,
                                                      all_touched=None,
                                                      stats_to_retrieve='sums_counts',
                                                      assert_projections_same=False, )
                        generated_scenario_label = 'gtap2_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label + '_existing_ag'
                        df.rename(columns={'sums': generated_scenario_label + '_sum', 'counts': generated_scenario_label + '_count'}, inplace=True)
                        df.to_excel(current_output_excel_path)
                    else:
                        df = pd.read_excel(current_output_excel_path, index_col=0)
                    merged_df = pd.merge(merged_df, df, how='outer', left_index=True, right_index=True)

                    # Also need to compute the value on that cropland that was cropland in both
                    # crop_value_baseline_existing_ag_path = os.path.join(p.cur_dir, 'crop_value_baseline_existing_ag_' + policy_scenario_label + '.tif')

                    def op(y, w, z):
                        r = dask.array.where((w == 2) & (z == 2), y, 0.)
                        rr = (z * 0.0) + r  # HACK. Dask.array.where was returning a standard xarray rather than a rioxarray. This dumb hack makes it inherit the rioxarray parameters of z
                        return rr

                    if not hb.path_exists(crop_value_baseline_existing_ag_path):
                        op_paths = [
                            baseline_crop_value_pollinator_adjusted_path,
                            lulc_path,
                            p.base_year_simplified_lulc_path,
                        ]

                        hb.dask_compute(op_paths, op, crop_value_baseline_existing_ag_path)

                    # Zonal stats for difference from Baseline
                    current_output_excel_path = os.path.join(p.cur_dir, 'crop_value_baseline_existing_ag_' + policy_scenario_label + '_zonal_stats.xlsx')
                    if not hb.path_exists(current_output_excel_path):
                        df = hb.zonal_statistics_flex(crop_value_baseline_existing_ag_path,
                                                      p.gtap37_aez18_path,
                                                      zone_ids_raster_path=p.zone_ids_raster_path,
                                                      id_column_label='pyramid_id',
                                                      zones_raster_data_type=5,
                                                      values_raster_data_type=6,
                                                      zones_ndv=-9999,
                                                      values_ndv=-9999,
                                                      all_touched=None,
                                                      stats_to_retrieve='sums_counts',
                                                      assert_projections_same=False, )
                        generated_scenario_label = 'crop_value_baseline_existing_ag_compared_to_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label
                        df.rename(columns={'sums': generated_scenario_label + '_sum', 'counts': generated_scenario_label + '_count'}, inplace=True)
                        df.to_excel(current_output_excel_path)
                    else:
                        df = pd.read_excel(current_output_excel_path, index_col=0)
                    merged_df = pd.merge(merged_df, df, how='outer', left_index=True, right_index=True)


                    if policy_scenario_label != 'BAU':
                        pass

                        # Calc difference between scenario and BAU for crop_value adjusted
                        # IMPORTANT NOTE: This is really just for plotting and visualization. The shockfiles themselves are all defined relative to the baseline, not relative to bau.
                        crop_value_difference_from_bau_path = os.path.join(p.cur_dir, 'crop_value_difference_from_bau_' + policy_scenario_label + '.tif')
                        bau_lulc_path = os.path.join(p.stitched_lulc_esa_scenarios_dir, 'lulc_seals7_gtap1_' + luh_scenario_label + '_' + str(scenario_year) + '_BAU.tif')
                        if not hb.path_exists(crop_value_difference_from_bau_path):
                            hb.dask_compute([bau_crop_value_pollinator_adjusted_path, current_crop_value_pollinator_adjusted_path], lambda x, y: y - x, crop_value_difference_from_bau_path)

                        # Zonal stats for difference from BAU
                        current_output_excel_path = os.path.join(p.cur_dir, 'crop_value_difference_from_bau_' + policy_scenario_label + '_zonal_stats.xlsx')
                        if not hb.path_exists(current_output_excel_path):
                            df = hb.zonal_statistics_flex(crop_value_difference_from_bau_path,
                                                          p.gtap37_aez18_path,
                                                          zone_ids_raster_path=p.zone_ids_raster_path,
                                                          id_column_label='pyramid_id',
                                                          zones_raster_data_type=5,
                                                          values_raster_data_type=6,
                                                          zones_ndv=-9999,
                                                          values_ndv=-9999,
                                                          all_touched=None,
                                                          stats_to_retrieve='sums_counts',
                                                          assert_projections_same=False, )
                            generated_scenario_label = 'gtap2_crop_value_difference_from_bau_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label
                            df.rename(columns={'sums': generated_scenario_label + '_sum', 'counts': generated_scenario_label + '_count'}, inplace=True)
                            df.to_excel(current_output_excel_path)
                        else:
                            df = pd.read_excel(current_output_excel_path, index_col=0)
                        merged_df = pd.merge(merged_df, df, how='outer', left_index=True, right_index=True)

                        # Calc difference between scenario and BAU for crop_value adjusted ON EXISTING AG
                        # ie for crop_value on grid-cells that were agri in both lulc maps.
                        lulc_path = os.path.join(p.stitched_lulc_esa_scenarios_dir, 'lulc_seals7_gtap1_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label + '.tif')
                        crop_value_difference_from_bau_existing_ag_path = os.path.join(p.cur_dir, 'crop_value_difference_from_bau_existing_ag_' + policy_scenario_label + '.tif')

                        def op(x, y, w, z):
                            r = dask.array.where((w == 2) & (z == 2), y - x, 0.)
                            rr = (z * 0.0) + r  # HACK. Dask.array.where was returning a standard xarray rather than a rioxarray. This dumb hack makes it inherit the rioxarray parameters of z
                            return rr

                        if not hb.path_exists(crop_value_difference_from_bau_existing_ag_path):
                            op_paths = [
                                bau_crop_value_pollinator_adjusted_path,
                                crop_value_pollinator_adjusted_path,
                                lulc_path,
                                bau_lulc_path,
                            ]

                            hb.dask_compute(op_paths, op, crop_value_difference_from_bau_existing_ag_path)

                        current_output_excel_path = os.path.join(p.cur_dir, 'crop_value_difference_from_bau_existing_ag_' + policy_scenario_label + '_zonal_stats.xlsx')
                        if not hb.path_exists(current_output_excel_path):
                            df = hb.zonal_statistics_flex(crop_value_difference_from_bau_existing_ag_path,
                                                          p.gtap37_aez18_path,
                                                          zone_ids_raster_path=p.zone_ids_raster_path,
                                                          id_column_label='pyramid_id',
                                                          zones_raster_data_type=5,
                                                          values_raster_data_type=6,
                                                          zones_ndv=-9999,
                                                          values_ndv=-9999,
                                                          all_touched=None,
                                                          stats_to_retrieve='sums_counts',
                                                          assert_projections_same=False, )


                            generated_scenario_label = 'gtap2_crop_value_difference_from_bau_existing_ag_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label
                            df.rename(columns={'sums': generated_scenario_label + '_sum', 'counts': generated_scenario_label + '_count'}, inplace=True)
                            df.to_excel(current_output_excel_path)
                        else:
                            df = pd.read_excel(current_output_excel_path, index_col=0)
                        merged_df = pd.merge(merged_df, df, how='outer', left_index=True, right_index=True)

        ###########################################
        ###### Calculate the actual shock as the mean change.
        ###########################################

        scenario_shock_column_names_to_keep = []
        # scenario_shock_column_names_to_keep = ['pyramid_id', 'pyramid_ids_concatenated', 'pyramid_ids_multiplied', 'gtap37v10_pyramid_id', 'aez_pyramid_id', 'gtap37v10_pyramid_name', 'ISO3', 'AZREG', 'AEZ_COMM']
        if not hb.path_exists(p.pollination_shock_csv_path):
            baseline_generated_scenario_label = 'gtap1_baseline_2014'
            baseline_generated_scenario_label_existing_ag = 'gtap1_baseline_2014_existing_ag'

            scenario_shock_column_names_to_keep.append(baseline_generated_scenario_label + '_sum')
            # scenario_shock_column_names_to_keep.append(baseline_generated_scenario_label_existing_ag + '_sum')

            for luh_scenario_label in p.luh_scenario_labels:
                for scenario_year in p.scenario_years:
                    for policy_scenario_label in p.policy_scenario_labels:

                        generated_label = 'gtap2_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label + '_pollination_shock'
                        merged_df[generated_label] = merged_df[generated_scenario_label + '_sum'] / merged_df[baseline_generated_scenario_label + '_sum']
                        scenario_shock_column_names_to_keep.append(generated_label)

                        # # NOTE: When calculating the value only on existing cells, cannot use the sum / sum method above. Need to use the new rasters created ad calculate their mean.
                        # generated_scenario_label = 'gtap2_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label + '_existing_ag'
                        # merged_df[generated_scenario_label + '_mean'] = merged_df[generated_scenario_label_existing_ag + '_sum'] / merged_df[baseline_generated_scenario_label_existing_ag + '_count']

                        # TODOOO: ALMOST got the full sim ready to run on the new pollination method but didn't finish getting the averages here calculated.

                        # generated_scenario_label_existing_ag = 'gtap2_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label + '_existing_ag'
                        # generated_label = 'gtap2_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label + '_pollination_shock'
                        # merged_df[generated_label] = merged_df[generated_scenario_label_existing_ag + '_sum'] / merged_df[baseline_generated_scenario_label_existing_ag + '_sum']
                        # scenario_shock_column_names_to_keep.append(generated_label)

                        # generated_scenario_existing_ag_label = 'gtap2_crop_value_difference_from_baseline_existing_ag_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label
                        # generated_label = 'gtap2_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label + '_pollination_shock_existing_ag'
                        #
                        # merged_df[generated_label] = merged_df[generated_scenario_existing_ag_label + '_sum'] / merged_df[baseline_generated_scenario_label + '_sum']

                        # # Also subtract the difference with BAU for each other policy
                        # if policy_scenario_label != 'BAU':
                        #     bau_label = 'gtap2_' + luh_scenario_label + '_' + str(scenario_year) + '_BAU_pollination_shock'
                        #     scenario_label = 'gtap2_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label + '_pollination_shock'
                        #     new_label = luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label + '_shock_minus_bau'
                        #     merged_df[new_label] = merged_df[scenario_label] - merged_df[bau_label]
                        #
                        #     # generated_scenario_existing_ag_label = 'gtap2_crop_value_difference_from_bau_existing_ag_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label
                        #     # merged_df[generated_scenario_existing_ag_label + '_mean'] = merged_df[generated_scenario_existing_ag_label + '_sum'] / merged_df[generated_scenario_existing_ag_label + '_count']
                        #     # merged_df[generated_scenario_existing_ag_label + '_mean_minus_baseline'] = merged_df[generated_scenario_existing_ag_label + '_mean'] - merged_df[baseline_generated_scenario_label + '_mean']
                        #
                        #
                        #     generated_bau_label  = 'gtap2_' + luh_scenario_label + '_' + str(scenario_year) + '_BAU'
                        #     generated_scenario_label = 'gtap2_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label
                        #
                        #     generated_bau_existing_ag_label  = 'gtap2_' + luh_scenario_label + '_' + str(scenario_year) + '_BAU'
                        #     merged_df[generated_scenario_label + '_sum_div_bau'] = merged_df[generated_scenario_label + '_sum'] / merged_df[generated_bau_label + '_sum']
                        #
                        #     generated_scenario_existing_ag_label = 'gtap2_crop_value_difference_from_bau_existing_ag_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label
                        #     # merged_df[generated_scenario_existing_ag_label + '_mean_minus_bau'] = merged_df[generated_scenario_existing_ag_label + '_mean'] - merged_df[generated_scenario_existing_ag_label + '_mean']
                        #
                        #     # merged_df[generated_scenario_existing_ag_label + '_sum_div_bau'] = merged_df[generated_scenario_existing_ag_label + '_sum'] / merged_df[generated_bau_existing_ag_label + '_sum']
                        #
                        #     generated_scenario_label= 'gtap2_crop_value_difference_from_bau_' + luh_scenario_label + '_' + str(scenario_year) + '_' + policy_scenario_label


            # write to csv and gpkg
            merged_df.to_csv(hb.suri(p.pollination_shock_csv_path, 'comprehensive'))
            gdf = gpd.read_file(p.gtap37_aez18_path)
            gdf = gdf.merge(merged_df, left_on='pyramid_id', right_index=True, how='outer')
            gdf.to_file(hb.suri(p.pollination_shock_change_per_region_path, 'comprehensive'), driver='GPKG')

            merged_df = merged_df[scenario_shock_column_names_to_keep]
            merged_df.to_csv(p.pollination_shock_csv_path)
            gdf = gpd.read_file(p.gtap37_aez18_path)
            gdf = gdf.merge(merged_df, left_on='pyramid_id', right_index=True, how='outer')
            gdf.to_file(p.pollination_shock_change_per_region_path, driver='GPKG')

