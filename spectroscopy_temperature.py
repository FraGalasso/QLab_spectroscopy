import numpy as np
import functions as fn
import pandas as pd
import os
import linewidth_functions as lfn
import uncertainties.unumpy as unp

'''
#load Rb transitions
Rb = [377104390084020.94, 377104798412020.94, 377105206740020.94, 377105909878483.7, 377106090669483.7, 377106271460483.7, 377108945610922.8, 377109126401922.8, 377109307192922.8, 377111224766631.8, 377111633094631.8, 377112041422631.8]
Rb_labels = ['21', 'cross', '22', '32', 'cross', '33', '22', 'cross', '23', '11', 'cross', '12']
'''
os.makedirs('data_temp_2/figures/calib', exist_ok=True)

# load data + plot
filename = "data_temp_2/clean_data/fullspec00000_new.csv"
data = pd.read_csv(filename, skiprows=1, names=['timestamp', 'photodiode', 'volt_piezo', 'volt_ld'], dtype={
                   'timestamp': float, 'photodiode': float, 'volt_piezo': float, 'volt_ld': float})
fn.plotting3(data['timestamp'], data['photodiode'], data['volt_piezo']/10, data['volt_ld'], 'Timestamps', 'Photodiode', 'Piezo voltage / 10',
             'Amplitude modulation', '', 'data_temp_2/figures/calib/time_constdiode.png', save=True)

# load peaks
peaks_filename = "data_temp_2/clean_data/fullspec00000_peaks.csv"
peaks = pd.read_csv(peaks_filename)

peaks_mean = np.array(peaks['lor_mean'])
peaks_gamma = np.array(peaks['lor_gamma'])
peaks_freqs = np.array(peaks['freq'])
peaks_d_mean = np.array(peaks['lor_d_mean'])

piezo = unp.uarray(peaks_mean, peaks_d_mean)


def quadratic_odr(par, x):
    return par[0] * x**2 + par[1] * x + par[2]


coeffs, dcoeffs = lfn.plot_ufloat_fit(x=piezo, y=peaks_freqs, model_func=quadratic_odr, x_label='Piezo Voltages',
                                      y_label='Frequencies [Hz]', title='', beta0=[1e7, 1e8, 1e14], file_name='data_temp_2/figures/calib/calibration.png', save=True)

with open('data_temp_2/calibration_coeffs.txt', 'w') as file:
    file.write('Calibration coefficients (a*x^2+b*x+c):\n')
    file.write(f'a = {coeffs[0]} +/- {dcoeffs[0]}\n')
    file.write(f'a = {coeffs[1]} +/- {dcoeffs[1]}\n')
    file.write(f'a = {coeffs[2]} +/- {dcoeffs[2]}\n')

frequencies = quadratic_odr(coeffs, data['volt_piezo'])

fn.plotting(frequencies/1e9, data['photodiode'], 'Calibrated frequency [GHz]', 'Signal at photodiode V',
              '', "data_temp_2/figures/calib/spectrum.png", save=True)

data['frequencies'] = frequencies
data.to_csv('data_temp_2/clean_data/fullspec00000_frequencies.csv', index=False)

# remove 2.5 * gamma_factor * gamma interval around peaks, plot and generate file with new data

gamma_factor = 2.5
cropped_data = fn.remove_peaks2(peaks_mean, peaks_gamma, data, gamma_factor)
cropped_data.to_csv('data_temp_2/clean_data/fullspec00000_frequencies_cropped.csv', index=False)

fn.scattering(cropped_data['frequencies']/1e9, cropped_data['photodiode'], 'Frequencies [GHz]', 'Signal at photodiode [V]',
              '', "data_temp_2/figures/calib/spectrum_cropped.png", save=True)

