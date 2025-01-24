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
peaks_filename = "data_temp/clean_data/mod_LD00000_peaks.csv"
peaks = pd.read_csv(peaks_filename)

peaks_time = np.array(peaks['timestamp'])
peaks_mean = np.array(peaks['lor_mean'])
peaks_gamma = np.array(peaks['lor_gamma'])
peaks_freqs = np.array(peaks['freq'])

# calculate calibration coefficients
coeff1, coeff2, d_coeff1, d_coeff2 = fn.plot_fits(peaks_mean, peaks_freqs, x_label=r'lormean($\text{V}_{\text{piezo}}$)',
                                                  y_label='Frequency [Hz]', title='Fixed diode current, modulated piezo position: time calib',
                                                  file_name=f"{figure_folder}/calibration.png", save=True)

# check distances
def quadratic(x, a, b, c):
    return a*(x**2) + b*x + c

print('Displacement of calibrated peaks from theoretical value (in MHz) for modld00000')
print((quadratic(peaks_mean, *coeff2) - peaks_freqs)/1e6)



'''
#THIS PART IS NOT STRICTLY NECESSARY, I JUST CALIBRATE THE CALIBRATION FILE ITSELF AND PLOT IT
filename = "data/mod_LD00000.csv"
data = pd.read_csv(filename)
frequencies = coeff2[0] * data['volt_piezo']**2 + coeff2[1] * data['volt_piezo'] + coeff2[2]
frequencies = np.array(frequencies)
#double check frequencies same length as data and then add frequencies to data pandas array. saves to csv
if len(data) != len(frequencies):
    raise ValueError("Data and frequencies must have the same length.")
data['frequencies']= frequencies
fn.plotting(data['frequencies'], data['volt_laser'],'calibrated frequency', 'signal at photodiode', 'spectroscopy of Rb, \n constant diode current, modulated piezo position', f"{figure_folder}/calibrated_calibdata.png", save=True)
# Construct new filename and save calibrated data to new file
base, ext = os.path.splitext(filename)
new_filename_calib = f"{base}_calibrated{ext}"
data.to_csv(new_filename_calib, index=False)
'''


# CALIBRATE THE DATA, SAVE IT, LOAD PEAKS, CROP THE PEAKS ZONE, SAVE

# load data and calibrate it
filename = "data_temp/clean_data/mod_LD00001.csv"
data = pd.read_csv(filename)


frequencies = coeff2[0] * data['volt_piezo']**2 + coeff2[1] * data['volt_piezo'] + coeff2[2]

frequencies = np.array(frequencies)
# double check frequencies same length as data and then add frequencies to data pandas array. saves to csv
if len(data) != len(frequencies):
    raise ValueError("Data and frequencies must have the same length.")
data['frequencies'] = frequencies

# Construct new filename and save calibrated data to new file
base, ext = os.path.splitext(filename)
new_filename_calib = f"{base}_calibrated{ext}"
new_filename_crop = f"{base}_calibrated_cropped{ext}"
data.to_csv(new_filename_calib, index=False)

# load peaks
peaks_filename = "data_temp/clean_data/mod_LD00001_peaks.csv"
peaks = pd.read_csv(peaks_filename)

peaks_mean = np.array(peaks['lor_mean'])
peaks_gamma = np.array(peaks['lor_gamma'])
peaks_freqs = np.array(peaks['freq'])

# remove 2 * gamma_factor * gamma interval around peaks, plot and generate file with new data
gamma_factor = 2.5
cropped_data = sfn.remove_peaks2(peaks_mean, peaks_gamma, data, gamma_factor)
cropped_data.to_csv(new_filename_crop, index=False)

# plot the calibrated and cropped data
sfn.scattering(data['frequencies'], data['volt_laser'], 'frequencies', 'signal at photodiode',
               'spectroscopy of Rb, \n modulated piezo position', f"{figure_folder}/data_calibrated.png", save=True, fit=False)
sfn.scattering(cropped_data['frequencies'], cropped_data['volt_laser'], 'frequencies', 'signal at photodiode',
               'spectroscopy of Rb, \n modulated piezo position', f"{figure_folder}/data_calib_crop.png", save=True, fit=False)


# calculate freq distance between peaks
print('Displacement of calibrated peaks from theoretical value (in MHz) for modld00001:')
peaks_freqs_meas = quadratic(peaks_mean, *coeff2)
displ = (peaks_freqs_meas - peaks_freqs)
print(displ/1e6)

print('Calibrated gammas (in MHz) for modld00001:')
gamma_freq = coeff2[0] * (2*peaks_gamma*peaks_mean + peaks_gamma**2) + coeff2[1] * peaks_gamma
print(gamma_freq/1e6)

print('Displacement of calibrated peaks from theoretical value (in units of gamma) for modld00001:')
print(displ/gamma_freq)

print('Distance between calibrated peaks (in MHz):')
dist_meas = np.diff(peaks_freqs_meas)
print(dist_meas/1e6)

print('Expected distance between peaks (in MHz):')
dist_theo = np.diff(peaks_freqs)
print(dist_theo/1e6)

print('Difference of distances in MHz:')
diff = dist_meas -  dist_theo
print(diff/1e6)

print('Difference of distances in units of gamma (the left one):')
diff = dist_meas -  dist_theo
print(diff/gamma_freq[0])
