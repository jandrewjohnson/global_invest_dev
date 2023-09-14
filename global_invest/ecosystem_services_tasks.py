import os, sys
import hazelbean as hb
import seals.seals_utils as seals_utils
import global_invest
from global_invest import ecosystem_services_functions
# from seals import seals_utils

def ecosystem_services(p):
    pass # Just to generate a folder

def carbon_storage(p):
    """Calculate carbon storage presentfrom LULC maps."""

    p.exhaustive_carbon_table_path = os.path.join(p.base_data_dir, "global_invest", "carbon", "exhaustive_carbon_table.csv")

    p.carbon_zones_path = os.path.join(p.base_data_dir, "global_invest", "carbon", 'carbon_zones_rasterized.tif')

    p.joined_carbon_table_stacked_with_seals_simplified_classes_path = os.path.join(p.cur_dir, 'joined_carbon_table_stacked_with_seals_simplified_classes.csv')

    # TODOO I may want to redo the run_subset_tiles functionality to instead just require a small aoi 
    # polygon (which would in principle only cover 1 tile. This is because the definition of BB gets confused
    # but also i would need to specify the difference between bb_pyramid vs bb_vector. Think about this.
    # Also note that the bb_pyramid would be determined by the processing_resolution

    if p.run_this:      
        for index, row in p.scenarios_df.iterrows():

            # NOTE THAT WE CALL THE SEALS_UTILS version of this function and not the hazelbean function
            # because we are assuming the structure of a seals scenarios.csv
            seals_utils.assign_df_row_to_object_attributes(p, row)
            p.L.info('Calculating carbon storage on ' + str(index) + ' of ' + str(len(p.scenarios_df)))

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


                

def pollination(p):
    """Calculate carbon storage presentfrom LULC maps."""

    p.exhaustive_carbon_table_path = os.path.join(p.base_data_dir, "global_invest", "carbon", "exhaustive_carbon_table.csv")

    p.carbon_zones_path = os.path.join(p.base_data_dir, "global_invest", "carbon", 'carbon_zones_rasterized.tif')

    p.joined_carbon_table_stacked_with_seals_simplified_classes_path = os.path.join(p.cur_dir, 'joined_carbon_table_stacked_with_seals_simplified_classes.csv')

    # TODOO I may want to redo the run_subset_tiles functionality to instead just require a small aoi 
    # polygon (which would in principle only cover 1 tile. This is because the definition of BB gets confused
    # but also i would need to specify the difference between bb_pyramid vs bb_vector. Think about this.
    # Also note that the bb_pyramid would be determined by the processing_resolution


    if p.run_this:
        # clip raster_zones with current bb


      
        for index, row in p.scenarios_df.iterrows():
            seals_utils.assign_df_row_to_object_attributes(p, row)
            p.L.info('Calculating carbon storage on ' + str(index) + ' of ' + str(len(p.scenarios_df)))


            if p.scenario_type == 'baseline':
                for year in p.base_years:

                    current_lulc_path = os.path.join(p.stitched_lulc_simplified_scenarios_dir  , 'lulc_' + p.lulc_src_label + '_' + p.lulc_simplification_label + '_' + p.model_label + '_' + str(year) + '.tif')
                    current_lulc_bb = hb.get_bounding_box(current_lulc_path)
                    pollination_sufficiency_output_path = os.path.join(p.cur_dir, 'pollination_sufficiency_' + p.model_label + '_' + str(year) + '.tif')
                    if not hb.path_exists(pollination_sufficiency_output_path):
                        p.L.info('Running global_invest_main.make_poll_suff on LULC: ' + str(current_lulc_path) + ' and saving results to ' + str(pollination_sufficiency_output_path))
                        ecosystem_services_functions.pollination_sufficiency(current_lulc_path, pollination_sufficiency_output_path)
            
            elif p.scenario_type != 'baseline':
                for year in p.years:
                    current_lulc_path = os.path.join(p.stitched_lulc_simplified_scenarios_dir  , 'lulc_' + p.lulc_src_label + '_' + p.lulc_simplification_label + '_' + p.exogenous_label + '_' + p.climate_label + '_' + p.model_label + '_' + p.counterfactual_label + '_' + str(year) + '.tif')
                    current_lulc_bb = hb.get_bounding_box(current_lulc_path)
                    pollination_sufficiency_output_path = os.path.join(p.cur_dir, 'pollination_sufficiency_' + p.exogenous_label + '_' + p.climate_label + '_' + p.model_label + '_' + p.counterfactual_label + '_' + str(year) + '.tif')
                    if not hb.path_exists(pollination_sufficiency_output_path):
                        p.L.info('Running global_invest_main.make_poll_suff on LULC: ' + str(current_lulc_path) + ' and saving results to ' + str(pollination_sufficiency_output_path))
                        ecosystem_services_functions.pollination_sufficiency(current_lulc_path, pollination_sufficiency_output_path)

