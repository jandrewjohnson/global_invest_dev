import pandas as pd

# Load the two uploaded CSV files with raw strings for paths
production_data_path = r"D:\Shared drives\NatCapTEEMs\Projects\Global GEP\Ecosystem Services SubFolders\Pollination\crop_production_by_country\pear_pollination_production.csv"
price_data_path = r"D:\Shared drives\NatCapTEEMs\Projects\Global GEP\Ecosystem Services SubFolders\Pollination\Price_data\Prices_E_All_Data_ori.csv"

# Load the pollination-dependent apple production data with encoding error handling
encodings = ['utf-8', 'ISO-8859-1', 'latin1']
for encoding in encodings:
    try:
        production_df = pd.read_csv(production_data_path, encoding=encoding)
        print(f"Successfully read production data with encoding: {encoding}")
        break
    except UnicodeDecodeError as e:
        print(f"Failed to read production data with encoding {encoding}. Error: {e}")

# If none of the encodings work, raise an error
if 'production_df' not in locals():
    raise ValueError("Unable to decode the production data file with tried encodings.")

# Load the price data with a different encoding due to special characters
price_df = pd.read_csv(price_data_path, encoding='ISO-8859-1')

# Filter the price data for "Apples" and "Producer Price (USD/tonne)" in the year 2000
apple_price_df = price_df[(price_df['Item'] == 'Pears') & (price_df['Element'] == 'Producer Price (USD/tonne)')]

# Keep only relevant columns: 'Area' (country), 'Y2000' (price in USD/tonne)
apple_price_2000_df = apple_price_df[['Area', 'Y2000']].copy()

# Rename the 'Y2000' column to something more descriptive
apple_price_2000_df.rename(columns={'Y2000': 'Price_2000_USD_per_tonne'}, inplace=True)

# Standardize the formatting of the 'Zone_ID' and 'Area' columns for merging
production_df['Zone_ID'] = production_df['Zone_ID'].str.strip().str.lower()
apple_price_2000_df['Area'] = apple_price_2000_df['Area'].str.strip().str.lower()

# Merge the apple price data with the production data based on the country
merged_df = pd.merge(production_df, apple_price_2000_df, left_on='Zone_ID', right_on='Area', how='left')

# Drop the 'Area' column as it's no longer needed after the merge
merged_df.drop(columns=['Area'], inplace=True)

# Check for missing price data
missing_prices = merged_df[merged_df['Price_2000_USD_per_tonne'].isna()]
if not missing_prices.empty:
    print("Missing price data for the following countries:")
    print(missing_prices[['Zone_ID']])

# Calculate the total value of crop production (Sum * Price)
merged_df['Total_Value_USD'] = merged_df['Sum'] * merged_df['Price_2000_USD_per_tonne']

# Display the first few rows of the merged DataFrame with the new column
print(merged_df.head())

# Save the merged DataFrame to a new CSV file with the total value column included
output_path = r"D:\Shared drives\NatCapTEEMs\Projects\Global GEP\Ecosystem Services SubFolders\Pollination\crop_production_by_country\pear_pollination_production_with_price_data.csv"
merged_df.to_csv(output_path, index=False)

print(f"The merged DataFrame with total value has been saved to: {output_path}")

# Calculate and display the overall total value for all zones
overall_total_value = merged_df['Total_Value_USD'].sum()
print(f"Overall total value of apple production: ${overall_total_value:,.2f}")