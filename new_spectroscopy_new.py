import numpy as np
import spectroscopy_functions as sfn
import functions as fn
import pandas as pd
import os

# file settings
folder_name = 'data_temp'
figure_folder = "figures/calib"
figure_folder = os.path.join(folder_name, figure_folder)
os.makedirs(figure_folder, exist_ok=True)
save = True

# CALIBRATION
# read and plot calibration data

# load peaks
peaks_filename = "data_temp/clean_data/mod_LD00001_peaks.csv"
peaks = pd.read_csv(peaks_filename)

peaks_time = np.array(peaks['timestamp'])
peaks_mean = np.array(peaks['lor_mean'])
peaks_gamma = np.array(peaks['lor_gamma'])
peaks_freqs = np.array(peaks['freq'])

# calculate calibration coefficients
coeff1, coeff2, d_coeff1, d_coeff2 = fn.plot_fits(peaks_mean, peaks_freqs, x_label=r'lormean($\text{V}_{\text{piezo}}$)',
                                                  y_label='Frequency [Hz]', title='Fixed diode current, modulated piezo position: time calib',
                                                  file_name=f"{figure_folder}/calibration_2.png", save=True)

def linear(x, a, b):
    return a * x + b


print('Displacement of calibrated peaks from theoretical value (in MHz) for modld00000')
print((linear(peaks_mean, *coeff1) - peaks_freqs)/1e6)

filename = "data_temp/clean_data/mod_LD00001.csv"
data = pd.read_csv(filename)

frequencies = np.array(linear(data['volt_piezo'], *coeff1))

# double check frequencies same length as data and then add frequencies to data pandas array. saves to csv
if len(data) != len(frequencies):
    raise ValueError("Data and frequencies must have the same length.")

data['frequencies'] = frequencies
base, ext = os.path.splitext(filename)
new_filename_calib = f"{base}_calibrated_2{ext}"
new_filename_crop = f"{base}_calibrated_cropped_2{ext}"
data.to_csv(new_filename_calib, index=False)

gamma_factor = 2.5
cropped_data = sfn.remove_peaks2(peaks_mean, peaks_gamma, data, gamma_factor)
cropped_data.to_csv(new_filename_crop, index=False)

# plot the calibrated and cropped data
sfn.scattering(data['frequencies'], data['volt_laser'], 'frequencies', 'signal at photodiode',
               'spectroscopy of Rb, \n modulated piezo position', f"{figure_folder}/data_calibrated_2.png", save=True, fit=False)
sfn.scattering(cropped_data['frequencies'], cropped_data['volt_laser'], 'frequencies', 'signal at photodiode',
               'spectroscopy of Rb, \n modulated piezo position', f"{figure_folder}/data_calib_crop_2.png", save=True, fit=False)
