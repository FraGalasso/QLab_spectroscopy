import matplotlib.pyplot as plt
import numpy as np
import functions as fn
import pandas as pd

'''
#load Rb transitions
Rb = [377104390084020.94, 377104798412020.94, 377105206740020.94, 377105909878483.7, 377106090669483.7, 377106271460483.7, 377108945610922.8, 377109126401922.8, 377109307192922.8, 377111224766631.8, 377111633094631.8, 377112041422631.8]
Rb_labels = ['21', 'cross', '22', '32', 'cross', '33', '22', 'cross', '23', '11', 'cross', '12']
'''


#load data + plot
filename = "data/Fixed_ld00000_peaks_fit_offset.csv"
data = pd.read_csv(filename, skiprows=1, names=['timestamp', 'photodiode', 'volt_piezo', 'volt_ld', 'offset'], dtype={'timestamp': float, 'photodiode': float, 'volt_piezo': float, 'volt_ld': float, 'offset': float})
fn.plotting3(data['timestamp'], data['photodiode'], data['volt_piezo']/10, data['volt_ld'], 'time', 'pd', 'piezo/10', 'ld', 'constant diode current, modulated piezo position', 'figures_temperature/time_constdiode.pdf', save = True)

#load peaks
peaks_filename = "data/Fixed_ld00000_peaks_fit.csv"
peaks = pd.read_csv(peaks_filename, skiprows=1, names=['indices','timestamp', 'pd_peaks','piezo_peaks','freq', 'lor_A','lor_mean','lor_gamma', 'lor_off','lor_d_A','lor_d_mean','lor_d_gamma','lor_d_off'])

peaks_time = np.array(peaks['timestamp'])
peaks_mean = np.array(peaks['lor_mean'])
peaks_gamma = np.array(peaks['lor_gamma'])
peaks_freqs = np.array(peaks['freq'])

'''
#In this section I checked weather calibrating the frequencies in the virtual time calibration (finding dv/dt) or in the piezo voltage calibration (dv/dVp) made a difference.
it didn't. 

#fitting peaks data to find freq conversions: dv/dVp (at const diode current) and dv/dt (v = nu)
#if we want to use lorentzian, the one in piezovolt doesnt work anymore
coeff1_piezo, coeff2_piezo = fn.plot_fits(peaks['piezo_peaks'], peaks['freq'], 'piezo_peaks', 'freq', 'Constant diode current, modulated piezo position', "figures_temperature/fit_constdiode3.pdf", save=True)
dv_dVp = coeff1[0]
frequencies_piezo = dv_dVp * data['volt_piezo'] + coeff1_piezo[1]
fn.plotting(frequencies_piezo, data['photodiode'],'calibrated frequency', 'signal at photodiode', 'spectroscopy of Rb, \n constant diode current, modulated piezo position', "figures_temperature/spec_constdiode_volt3.pdf", save=True)

coeff1_time, coeff2_time = fn.plot_fits(peaks['timestamps'], peaks['freq'], 'time', 'freq', 'Fixed Piezo, modulated diode current: time calib', "figures_temperature/timefit_constdiode3.pdf", save=True)
dv_dt = coeff1_time[0]
frequenceis_time = dv_dt * data['timestamp'] + coeff1_time[1]
fn.plotting(frequencies_time, data['photodiode'],'calibrated frequency', 'signal at photodiode', 'spectroscopy of Rb, \n constant diode current, modulated piezo position', "figures_temperature/spec_constdiode_time3.pdf", save=True)

#comparison:
fn.plotting(data['photodiode'], frequencies_time - frequencies_piezo, 'signal at photodiode','calibrated frequency', 'spectroscopy of Rb, \n constant diode current, modulated piezo position', "figures_temperature/spec_constdiode_diff3.pdf", save=True)

#conclusion: calibrating with time or with piezo voltage doesnt make a difference and is equivalent
so we will now proceed with the virtual time calibration only
'''

#find timestamp <-> frequency conversion. generate the array of frequencies and plot
'''
coeff1_time, coeff2_time = fn.plot_fits(peaks_mean, peaks_freqs, 'lormean', 'freq', 'Fixed diode current, modulated piezo position: time calib', "figures_temperature/timefit_constdiode3.pdf", save=True)
dv_dt = coeff1_time[0]
frequencies = dv_dt * data['timestamp'] + coeff1_time[1]
frequencies = np.array(frequencies)
fn.plotting(frequencies, data['photodiode'],'calibrated frequency', 'signal at photodiode', 'spectroscopy of Rb, \n constant diode current, modulated piezo position', "figures_temperature/spec_constdiode3.pdf", save=True)
'''

#calib quadratic
coeff1_time, coeff2_time = fn.plot_fits(peaks_mean, peaks_freqs, 'lormean', 'freq', 'Fixed diode current, modulated piezo position: time calib', "figures_temperature/timefit_constdiode3.pdf", save=True)

frequencies = coeff2_time[0] * data['timestamp']**2 + coeff2_time[1] * data['timestamp'] + coeff2_time[2]
frequencies = np.array(frequencies)
fn.plotting(frequencies, data['photodiode'],'calibrated frequency', 'signal at photodiode', 'spectroscopy of Rb, \n constant diode current, modulated piezo position', "figures_temperature/spec_constdiode3.pdf", save=True)

#double check frequencies same length as data and then add frequencies to data pandas array. saves to csv
if len(data) != len(frequencies):
    raise ValueError("Data and frequencies must have the same length.")
data['frequencies']= frequencies
data.to_csv('data/Fixed_ld00000_frequencies_3.csv', index=False)

#remove 2 * gamma_factor * gamma interval around peaks, plot and generate file with new data

gamma_factor = 2.5
cropped_data = fn.remove_peaks(peaks_mean, peaks_gamma, data, gamma_factor)

fn.scattering(cropped_data['timestamp'], cropped_data['photodiode'],'timestamp', 'signal at photodiode', 'spectroscopy of Rb, \n constant diode current, modulated piezo position', "figures_temperature/spec_cropped_time3.pdf", save=True)
fn.scattering(cropped_data['frequencies'], cropped_data['photodiode'],'frequencies', 'signal at photodiode', 'spectroscopy of Rb, \n constant diode current, modulated piezo position', "figures_temperature/spec_cropped_freq3.pdf", save=True)

cropped_data.to_csv('data/Fixed_ld00000_frequencies_cropped_3.csv', index=False)


