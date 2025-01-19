import os
import pandas as pd
import matplotlib.pyplot as plt
import functions as fn


'''This code reads the csv files, crops them to a window where the piezo voltage is just increasing,
removes NaNs and fits the piezo data linearly. It then saves a csv file with  4 columns: timestamps,
transmission, volt_piezo, piezo_fitted.'''

#####################################################################################################
# section with file specific values

# Define the folder and file paths
folder_name = 'data_temp'
title = 'mod_LD00001'
df = pd.DataFrame()
# Import single channel files
channel_files = [f'{folder_name}/raw/C{i}{title}.csv' for i in range(1, 4)]
# Define channel names
column_titles = ['volt_laser', 'volt_piezo', 'volt_ld']
os.makedirs(f"{folder_name}/clean_data", exist_ok=True)
os.makedirs(f"{folder_name}/figures/simple_plots_time", exist_ok=True)
# Boolean to activate extra cropping, if necessary:
cropping = False
# You have to set the mask variable to what you need in the code
# If you are using more than 3 channels you have to modify the rest of the code

#####################################################################################################

# Loop through each file in the file_paths
for i, channel in enumerate(channel_files):
    file_name = os.path.basename(channel).replace(
        '.csv', '')  # Extract file name without extension

    # Read the CSV file, skip the first 5 rows, and specify the data types
    # Encoding is different from usual utf-8, so we need to specify it.
    data = pd.read_csv(channel, skiprows=5, names=['timestamp', 'ch'], encoding='cp1252', dtype={
                       'timestamp': float, 'ch': float})

    # Write timestamp column just once
    if df.empty:
        df['timestamp'] = data.iloc[:, 0]

    # Write channel reading to dataframe
    df[column_titles[i]] = data.iloc[:, 1]


df = df.dropna()

# Producing single numpy arrays for manipulation with functions
# This may not be the most efficient way to handle this but it works at least
timestamps = df['timestamp'].to_numpy()
volt_laser = df['volt_laser'].to_numpy()
volt_piezo = df['volt_piezo'].to_numpy()
volt_ld = df['volt_ld'].to_numpy()

# Extra cropping, use it in case of mode hopping
if cropping:
    mask = (timestamps >= -0.01)
    timestamps = timestamps[mask]
    volt_laser = volt_laser[mask]
    volt_piezo = volt_piezo[mask]
    volt_ld = volt_ld[mask]

# Cropping data to a single sweep
# NOTICE, VOLT PIEZO HAS TO BE THE FIRST ONE
volt_piezo, timestamps, volt_laser, volt_ld = fn.crop_to_min_max(
    volt_piezo, timestamps, volt_laser, volt_ld)

piezo_fitted = fn.fit_piezo_line(timestamps, volt_piezo)
ld_fitted = fn.fit_piezo_line(timestamps, volt_ld)

# Plotting, just for fun
plt.figure()
plt.scatter(timestamps, volt_laser, label='Laser intensity',
            color='blue', s=5, marker='.')
plt.scatter(timestamps, piezo_fitted/10, label='Piezo voltage/10',
            color='red', s=5, marker='.')
plt.scatter(timestamps, volt_ld, label='Diode modulation',
            color='green', s=5, marker='.')
plt.xlabel('Timestamp [s]')
plt.ylabel('Channel Value [V]')
plt.title('Timestamp vs Channel Data')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig(f'{folder_name}/figures/simple_plots_time/{title}.png')
plt.show()

# Saving data in clean_data folder
output_file = f'{folder_name}/clean_data/{title}.csv'
df = df.iloc[:len(timestamps)].reset_index(drop=True)
df['timestamp'] = timestamps
df['volt_laser'] = volt_laser
df['volt_piezo'] = piezo_fitted
df['volt_ld'] = ld_fitted

df.to_csv(output_file, index=False)
print(f"Data saved to {output_file}")
