import os
import pandas as pd
import functions as fn


'''This code reads the csv files, crops them to a window where the piezo voltage is just increasing,
removes NaNs and fits the piezo data linearly. It then saves a csv file with  4 columns: timestamps,
volt_laser, volt_piezo, piezo_fitted.'''

# Define the folder and file paths
folder_name = 'data6'
# Example: process files scope_24.csv to scope_25.csv
file_paths = [f"{folder_name}/scope_{i}.csv" for i in range(15, 26)]
os.makedirs(f"{folder_name}/clean_data", exist_ok=True)

# Loop through each file in the file_paths
for file_path in file_paths:
    file_name = os.path.basename(file_path).replace(
        '.csv', '')  # Extract file name without extension

    # Read the CSV file, skip the first 2 rows, and specify the data types
    data = pd.read_csv(file_path, skiprows=2, names=['timestamp', 'volt_laser', 'volt_piezo'], dtype={
                       'timestamp': float, 'volt_laser': float, 'volt_piezo': float})

    # Remove any rows with NaN values
    data_cleaned = data.dropna()

    # Extract the columns
    timestamps = data_cleaned['timestamp'].to_numpy()  # Convert to NumPy array
    # Convert to NumPy array
    volt_laser = data_cleaned['volt_laser'].to_numpy()
    # Convert to NumPy array
    volt_piezo = data_cleaned['volt_piezo'].to_numpy()

    # Crop data to one piezo cycle
    result = fn.crop_to_min_max(timestamps, volt_laser, volt_piezo)

    # Redefine the arrays with the cropped versions
    timestamps = result[0]
    volt_laser = result[1]
    volt_piezo = result[2]

    # Fit the piezo data
    piezo_fitted = fn.fit_piezo_line(timestamps, volt_piezo)

    # Adjust the file path as needed
    output_file = figure_name = f'{folder_name}/clean_data/{file_name}_cropped.csv'
    fn.save_to_csv(timestamps, volt_laser, volt_piezo,
                   piezo_fitted, output_file)
