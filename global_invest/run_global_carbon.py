import os, sys
import hazelbean as hb
import pandas as pd

main = ''
if __name__ == '__main__':

    ### ------- ENVIRONMENT SETTINGS -------------------------------
    # Users should only need to edit lines in this ENVIRONMENT SETTINGS section

    # A ProjectFlow object is created from the Hazelbean library to organize directories and enable parallel processing.
    # project-level variables are assigned as attributes to the p object (such as in p.base_data_dir = ... below)
    # The only agrument for a project flow object is where the project directory is relative to the current_working_directory.
    # This organization, defined with extra dirs relative to the user_dir, is the EE-spec.
    user_dir = os.path.expanduser('~')        
    extra_dirs = ['Files', 'global_invest', 'projects']

    # The project_name is used to name the project directory below. Also note that
    # ProjectFlow only calculates tasks that haven't been done yet, so adding 
    # a new project_name will give a fresh directory and ensure all parts
    # are run.
    project_name = 'test_global_invest'
    project_name = project_name + '_' + hb.pretty_time()
    
    # The project-dir is where everything will be stored, in particular in an input, intermediate, or output dir
    # IMPORTANT NOTE: This should not be in a cloud-synced directory (e.g. dropbox, google drive, etc.), which
    # will either make the run fail or cause it to be very slow. The recommended place is (as coded above)
    # somewhere in the users's home directory.
    project_dir = os.path.join(user_dir, os.sep.join(extra_dirs), project_name)

    # Create the ProjectFlow Object
    p = hb.ProjectFlow(project_dir)
    
    p.run_in_parallel = 1

    # NEXT: finish making this the standard test by making it faster (just rwa), noting that the test_manual_shift file is still needed to show different permutations of regional, coarse projections.
    # Build the task tree via a building function and assign it to p
    # IF YOU WANT TO LOOK AT THE MODEL LOGIC, INSPECT THIS FUNCTION
    
    START HERE: Figure out what way we chose to do task trees in separate projects. Was it a different initialize function?
    seals_initialize_project.build_standard_run_task_tree(p)

    # Set the base data dir. The model will check here to see if it has everything it needs to run.
    # If anything is missing, it will download it. You can use the same base_data dir across multiple projects.
    # Additionally, if you're clever, you can move files generated in your tasks to the right base_data_dir
    # directory so that they are available for future projects and avoids redundant processing.
    # NOTE THAT the final directory has to be named base_data to match the naming convention on the google cloud bucket.
    p.base_data_dir = os.path.join(user_dir, 'Files/base_data')

    # In order for SEALS to download using the google_cloud_api service, you need to have a valid credentials JSON file.
    # Identify its location here. If you don't have one, email jajohns@umn.edu. The data are freely available but are very, very large
    # (and thus expensive to host), so I limit access via credentials.
    p.data_credentials_path = None
    # p.data_credentials_path = '..\\api_key_credentials.json'

    # There are different versions of the base_data in gcloud, but best open-source one is 'gtap_invest_seals_2023_04_21'
    p.input_bucket_name = None
    
    ## Set defaults and generate the scenario_definitions.csv if it doesn't exist.
    # SEALS will run based on the scenarios defined in a scenario_definitions.csv
    # If you have not run SEALS before, SEALS will generate it in your project's input_dir.
    # A useful way to get started is to to run SEALS on the test data without modification
    # and then edit the scenario_definitions.csv to your project needs.
    # Some of the other test files use different scenario definition csvs 
    # to illustrate the technique. If you point to one of these 
    # (or any one CSV that already exists), SEALS will not generate a new one.   
    p.scenario_definitions_filename = 'test_standard_scenarios.csv' 
    p.scenario_definitions_path = os.path.join(p.input_dir, p.scenario_definitions_filename)

    # If the scenarios csv doesn't exist, generate it and put it in the input_dir 
    if not hb.path_exists(p.scenario_definitions_path):
        
        # Before generating a new scenarios file, check if there's not one in the base data with the matching name.
        possible_path = p.get_path('seals', 'default_inputs', p.scenario_definitions_filename)
        if hb.path_exists(possible_path):
            hb.path_copy(possible_path, p.scenario_definitions_path)
            
        # If theres truly nothing in the base data, generate it for the default.
        else:        
            # There are multiple scenario_csv generator functions. Here we use the default.
            seals_utils.set_attributes_to_dynamic_default(p) # Default option

            # Once the attributes are set, generate the scenarios csv and put it in the input_dir.
            seals_utils.generate_scenarios_csv_and_put_in_input_dir(p)
        
        # After writing, read it it back in, cause this is how other attributes might be modified
        p.scenarios_df = pd.read_csv(p.scenario_definitions_path)
    else:
        # Read in the scenarios csv and assign the first row to the attributes of this object (in order to setup additional 
        # project attributes like the resolutions of the fine scale and coarse scale data)
        p.scenarios_df = pd.read_csv(p.scenario_definitions_path)

    # Set p attributes from df (but only the first row, cause its for initialization)
    for index, row in p.scenarios_df.iterrows():
        
        # NOTE! This also downloads any files references in the csv
        seals_utils.assign_df_row_to_object_attributes(p, row)
        break # Just get first for initialization.
        
        
    ### ------- DERIVED ATTRIBUTES ---------------------------------  
    # From here on, simple runs should not require editing any of these lines
        
    # If you want to run SEALS with the run.py file in a different directory (ie in the project dir)
    # then you need to add the path to the seals directory to the system path.
    custom_seals_path = None
    if custom_seals_path is not None: # G:/My Drive/Files/Research/seals/seals_dev/seals
        sys.path.insert(0, custom_seals_path)

    # To support loading an extent by AOI, define a country-iso3 vector here.
    aoi_creation_vector_ref_path = os.path.join('pyramids', 'countries_iso3.gpkg')
    p.countries_iso3_path = p.get_path(aoi_creation_vector_ref_path)

    # Some variables need further processing into attributes, like parsing a correspondence csv into a dict.
    seals_utils.set_derived_attributes(p)

    # Set processing resolution: determines how large of a chunk should be processed at a time
    # 4 deg is about max for 64gb memory systems
    # SEALS has two other implied resolutions: fine and coarse. In most applications, fine is 10 arcseconds (~300m at equator, based on ESACCI)
    # and coarse is based on IAM results that are 900 arcseconds (LUH2) or 1800 arcseconds (MAgPIE). Note that there is a coarser-yet
    # scale possible from e.g. GTAP-determined endogenous LUC. This is excluded in the base SEALS config.
    # These are set based on the input data, so here we will analyze two files in the input csv
    # and set the resolutions based on that, but first we need to download the files.
    p.processing_resolution = 1.0 # In degrees. Must be in pyramid_compatible_resolutions
    p.processing_resolution_arcseconds = p.processing_resolution * 3600.0 # MUST BE FLOAT

    # calibration_parameters_override_dict can be used in specific scenarios to e.g., not allow expansion of cropland into forest by overwriting the default calibration. Note that I haven't set exactly how this
    # will work if it is specific to certain zones or a single-zone that overrides all. The latter would probably be easier.
    # If the DF has a value, override. If it is None or "", keep from parameters source.
    p.calibration_parameters_override_dict = {}
    # p.calibration_parameters_override_dict['rcp45_ssp2'][2030]['BAU'] = os.path.join(p.input_dir, 'calibration_overrides', 'prevent_cropland_expansion_into_forest.xlsx')

    # TODO Switch to hb.log()
    p.L = hb.get_logger('test_run_seals')

    ### ------- RUNTIME OPTIONS ------------------------------------------------------

    p.build_overviews_and_stats = 0  # For later fast-viewing, this can be enabled to write ovr files and geotiff stats files. NYI anywhere.
    p.force_to_global_bb = 0
    p.plotting_level = 0 

    p.cython_reporting_level = 0
    p.calibration_cython_reporting_level = 0
    p.output_writing_level = 0  # >=2 writes chunk-baseline lulc
    p.write_projected_coarse_change_chunks = 0  # in the SEALS allocation, for troubleshooting, it can be useful to see what was the coarse allocation input.
    p.write_calibration_generation_arrays = 0  # in the SEALS allocation, for troubleshooting, it can be useful to see what was the coarse allocation input.

    p.years_to_convolve_override = None # Because convolving a lulc map is so computationally expensive, you may set this option to manually set which year of convolutions to use. If this is None (the default), it will convolve the base year.

    p.num_workers = None  # None sets it to max available. Otherwise, set to an integer.
    
    # Assign useful locals to project flow level
    p.user_dir = user_dir
    p.project_name = project_name
    p.project_dir = project_dir
    
    seals_initialize_project.initialize_paths(p)
    seals_initialize_project.run(p)

    result = 'Done!'


