import os
import pygeoprocessing
from osgeo import ogr, gdal
from pathlib import Path
import pandas as pd
from rapidfuzz import process
import hazelbean as hb

from global_invest.pollination.pollination_funcs import calculate_pollination_dependent_production, calculate_zonal_statistics


def pollination_gep(p):

    p.pollination_dependency_refpath = 'global_invest/Pollination/poll_dep_greater_0.xlsx'
    p.base_crop_data_refdir = 'crops/earthstat/crop_production'
    p.pollination_sufficiency_refpath = 'global_invest/Pollination/pollination_sufficiency_2017_300sec.tif'
    p.price_data_refpath = "global_invest/pollination/Price_data/Prices_E_All_Data.csv"
    p.earthstat_fao_crop_names_correspondence_refpath = 'global_invest/Pollination/earthstat_fao_crop_names_correspondence.csv'
    
    if p.run_this:
        # Define output locations
        p.pollination_output_raster_dir = os.path.join(p.cur_dir, 'pollination_dependent_crops')
        p.pollination_output_table_dir = os.path.join(p.cur_dir, 'crop_production_by_country')
        
        # Determine path from refpath.
        p.pollination_dependency_path = p.get_path(p.pollination_dependency_refpath, verbose=True)
        p.base_crop_data_dir = p.get_path(p.base_crop_data_refdir)
        p.pollination_sufficiency_path = p.get_path(p.pollination_sufficiency_refpath)
        p.price_data_path = p.get_path(p.price_data_refpath)
    
        # Initialize lists to store results as we iterate through crops.
        best_matches = []   
        total_values = []
        
        df = pd.read_excel(p.pollination_dependency_path)
        price_df = pd.read_csv(p.price_data_path)
        print(price_df.head())


        # Process each row in the DataFrame
        for _, row in df.iterrows():
            crop_filenm = row['filenm']
            dependence_ratio = row['poll.dep']
            commodity_name = row['Commodity_1']

            # Convert crop_filenm and commodity_name to lowercase for matching
            crop_filenm_lower = crop_filenm.lower()
            commodity_name_lower = commodity_name.lower()

            # Locate the crop production raster file
            crop_production_raster = os.path.join(p.base_crop_data_dir, crop_filenm, f"{crop_filenm}_production_Mg.tif")

            # Check if the raster file exists and set the status accordingly
            if os.path.exists(crop_production_raster):
                raster_status = "Found"
                # Define output raster path
                output_raster = os.path.join(p.pollination_output_raster_dir, f"pollination_dependent_{crop_filenm}.tif")
                hb.create_directories(output_raster)
                # Execute the raster operation to calculate pollination-dependent production
                pygeoprocessing.raster_calculator(
                    [(p.pollination_sufficiency_path, 1), (crop_production_raster, 1), (dependence_ratio, 'raw')],
                    calculate_pollination_dependent_production,
                    output_raster,
                    gdal.GDT_Float32,
                    None
                )

                print(f"Pollination-dependent production raster generated for {crop_filenm} and saved to {output_raster}")

                # Define the output CSV file path based on the raster file name
                output_table = os.path.join(p.pollination_output_table_dir, Path(output_raster).stem + "_crop_production_by_country.csv")
                hb.create_directories(output_table)
                # Perform zonal statistics and save the results
                zone_input_path = p.countries_vector_path
                calculate_zonal_statistics(output_raster, zone_input_path, output_table)
                
                production_df = pd.read_csv(output_table)

                # Convert the 'Item' column to lowercase for fuzzy matching
                price_df['Item_lower'] = price_df['Item'].str.lower()

                # Find the best match for the crop using either crop_filenm or commodity_name
                best_match = process.extractOne(crop_filenm_lower, price_df['Item_lower'], score_cutoff=80)
                if not best_match:
                    best_match = process.extractOne(commodity_name_lower, price_df['Item_lower'], score_cutoff=80)

                if best_match:
                    matched_item = best_match[0]
                    print(f"Best match for {crop_filenm} or {commodity_name}: {matched_item}")

                    # Store the best match in the list for later saving to CSV
                    best_matches.append({
                        'crop_filenm': crop_filenm,
                        'commodity_name': commodity_name,
                        'matched_item': matched_item,
                        'production_raster': crop_production_raster,
                        'raster_status': raster_status
                    })

                    # Filter the price data using the best match
                    crop_price_df = price_df[(price_df['Item_lower'] == matched_item) & (price_df['Element'] == 'Producer Price (USD/tonne)')]

                    # Print the filtered price data for debugging
                    print(crop_price_df)

                    # If price data is found, save it to a CSV file with both the price name and crop name
                    if not crop_price_df.empty:
                        price_name = crop_price_df['Item'].iloc[0]
                        output_price_table = os.path.join(p.pollination_output_table_dir, f"{price_name}_{crop_filenm}_price_data.csv")
                        crop_price_df.to_csv(output_price_table, index=False)
                        print(f"Saved price data for {price_name} and {crop_filenm} to: {output_price_table}")

                    # Keep only relevant columns: 'Area' (country) and the relevant year (e.g., 'Y2000')
                    crop_price_df_filtered = crop_price_df[['Area', 'Y2000']].copy()
                    print(crop_price_df_filtered)

                    # Rename the 'Y2000' column to something more descriptive
                    crop_price_df_filtered.rename(columns={'Y2000': 'Price_2000_USD_per_tonne'}, inplace=True)

                    # Standardize the formatting of the 'Zone_ID' and 'Area' columns for merging
                    production_df['Zone_ID'] = production_df['Zone_ID'].astype(str).str.strip().str.lower()
                    crop_price_df_filtered['Area'] = crop_price_df_filtered['Area'].str.strip().str.lower()

                    # Merge the crop price data with the production data based on the country
                    merged_df = pd.merge(production_df, crop_price_df_filtered, left_on='Zone_ID', right_on='Area', how='left')
                    print(merged_df)

                    # Drop the 'Area' column as it's no longer needed after the merge
                    merged_df.drop(columns=['Area'], inplace=True)

                    # Check for missing price data
                    missing_prices = merged_df[merged_df['Price_2000_USD_per_tonne'].isna()]
                    if not missing_prices.empty:
                        print(f"Missing price data for the following countries for {crop_filenm}:")
                        print(missing_prices[['Zone_ID']])

                    # Calculate the total value of crop production (Sum * Price)
                    merged_df['Total_Value_USD'] = merged_df['Sum'] * merged_df['Price_2000_USD_per_tonne']

                    # Count how many countries have price data in USD
                    price_data_count = merged_df['Price_2000_USD_per_tonne'].notna().sum()

                    # Display the first few rows of the merged DataFrame with the new column
                    print(merged_df.head())

                    # Save the merged DataFrame to a new CSV file with the total value column included
                    output_path = os.path.join(p.pollination_output_table_dir, f"{crop_filenm}_pollination_production_with_price_data.csv")
                    merged_df.to_csv(output_path, index=False)

                    print(f"The merged DataFrame with total value for {crop_filenm} has been saved to: {output_path}")

                    # Calculate and store the overall total value for the crop
                    overall_total_value = merged_df['Total_Value_USD'].sum()
                    print(f"Overall total value of {crop_filenm} production: ${overall_total_value:,.2f}")

                    # Store the total value and price data count for the crop in the list for later saving to CSV
                    total_values.append({
                        'crop_filenm': crop_filenm,
                        'commodity_name': commodity_name,
                        'overall_total_value': overall_total_value,
                        'price_data_count': price_data_count
                    })

                else:
                    print(f"No good match found for {crop_filenm} or {commodity_name} in price data.")
                    best_matches.append({
                        'crop_filenm': crop_filenm,
                        'commodity_name': commodity_name,
                        'matched_item': 'No match found',
                        'production_raster': crop_production_raster,
                        'raster_status': raster_status
                    })

            else:
                raster_status = "Not Found"
                print(f"Production raster not found for {crop_filenm}. Skipping.")
                best_matches.append({
                    'crop_filenm': crop_filenm,
                    'commodity_name': commodity_name,
                    'matched_item': 'No match found',
                    'production_raster': crop_production_raster,
                    'raster_status': raster_status
                })

        # Save all best matches to a CSV file
        best_matches_df = pd.DataFrame(best_matches)
        best_matches_csv_path = os.path.join(p.pollination_output_table_dir, "best_matches.csv")
        best_matches_df.to_csv(best_matches_csv_path, index=False)
        print(f"Best matches for all crops saved to: {best_matches_csv_path}")

        # Save all total values to a CSV file
        total_values_df = pd.DataFrame(total_values)
        total_values_csv_path = os.path.join(p.pollination_output_table_dir, "overall_total_values.csv")
        total_values_df.to_csv(total_values_csv_path, index=False)
        print(f"Overall total values for all crops saved to: {total_values_csv_path}")

        print("Processing complete for all crops.")



