import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import functions as fn

#constants
kb = 1.38e-23  # J/K
c = 3e8
h = 6.62607015e-34
rb_87mass = 1.67e-27 * 87   # kg
rb_85mass = 1.67e-27 * 85   # kg
gamma87 = 2 * np.pi * 5.75e6
gamma85 = 2 * np.pi * 5.75e6 #si è la stessa ho controllato

f_22 = 377108945610922.8
f_23 = 377109307192922.8
f_cross85 = (f_22 + f_23) / 2
f_12 =  377112041422631.8
f_11 = 377111224766631.8
f_cross87 = (f_12 + f_11) /2

Isat22 = ( np.pi/3 ) * (h * gamma85 / c**2) * (f_22)
Isatcross85 = ( np.pi/3 ) * (h * gamma85 / c**2) * (f_cross85)
Isat23 = ( np.pi/3 ) * (h * gamma85 / c**2) * (f_23)

Isat12 = ( np.pi/3 ) * (h * gamma87 / c**2) * (f_12)
Isatcross87 = ( np.pi/3 ) * (h * gamma87 / c**2) * (f_cross87)
Isat11 = ( np.pi/3 ) * (h * gamma87 / c**2) * (f_11)

area_beam = 1.7671458676442586e-6 #m^2

powers = [0.079, 0.07, 0.1, 0.221, 0.092, 0.104, 0.442, 1.009, 0.313, 0.7, 0.295]
powers_cross87 = [0.079, 0.07, 0.1, 0.221, 0.092, 0.104, 0.442, 0.313, 0.295]
intensities = np.array(powers)/area_beam
intensities_cross87 = np.array(powers_cross87)/area_beam

save = True
folder_name = 'data_linewidth'
file_names = [f"intensity0000{i}" for i in range(0, 10)]
file_names.append("intensity00010")
file_paths = [f'{folder_name}/clean_data/{file_name}.csv' for file_name in file_names]
file_peaks = [f'{folder_name}/clean_data/{file_name}_peaks.csv' for file_name in file_names]


peaks_22 = []
peaks_cross85 = [] 
peaks_23 = []
peaks_11 = []
peaks_cross87 = [] 
peaks_12 = []

#calculate calibration coeffs from Intensity00000

data = pd.read_csv(file_paths[0], skiprows=1, names=['timestamp','photodiode','volt_piezo','volt_ld'])
peaks = pd.read_csv(file_peaks[0], skiprows=1, names=['indices','timestamp','pd_peaks','piezo_peaks','freq', 'lor_A','lor_mean','lor_gamma','lor_off','lor_d_A','lor_d_mean','lor_d_gamma','lor_d_off'])

peaks_gamma = np.array(peaks['lor_gamma'])
peaks_mean = np.array(peaks['lor_mean'])
peaks_freqs = np.array(peaks['freq'])

coeff1, __ = fn.plot_fits(peaks_mean, peaks_freqs, 'lormean_piezovolt', 'freq', 'Fixed diode current, modulated piezo position: time calib', f"x.pdf", save= False)
fn.ensure_dir(folder_name + '/figures')

# for each file: reads data, finds calibrated frequencies, plots, and appends peaks to peak arrays
for i, (file_name, file_path, file_peak) in enumerate(zip(file_names, file_paths, file_peaks)):
    figure_path = os.path.join(folder_name, 'figures', f"{powers[i]:.3g}_{file_name}.pdf")

    # Read the CSV file, and specify the data types
    data = pd.read_csv(file_path, skiprows=1, names=['timestamp','photodiode','volt_piezo','volt_ld','offset','frequencies'])
    peaks = pd.read_csv(file_peak, skiprows=1, names=['indices','timestamp','pd_peaks','piezo_peaks','freq', 'lor_A','lor_mean','lor_gamma','lor_off','lor_d_A','lor_d_mean','lor_d_gamma','lor_d_off'])
    
    peaks_gamma = np.array(peaks['lor_gamma'])
    peaks_mean = np.array(peaks['lor_mean'])
    peaks_freqs = np.array(peaks['freq'])
    peaks_gamma = coeff1[0]*peaks_gamma
    
    frequencies = fn.calibrate(data['volt_piezo'], coeff1)

    fn.plotting(
        frequencies,
        data['photodiode'],
        'freq',
        'Photodiode Signal',
        f'Pump Laser Power = {powers[i]}, Intensity = {intensities[i]}',
        figure_path,
        save
    )

    peaks_22.append(peaks_gamma[0])
    peaks_cross85.append(peaks_gamma[1])
    peaks_23.append(peaks_gamma[2])
    peaks_11.append(peaks_gamma[3])

    if len(peaks) == 5:
        peaks_12.append(peaks_gamma[4])
    else:
        peaks_cross87.append(peaks_gamma[4])
        peaks_12.append(peaks_gamma[5])


figure_name = os.path.join(folder_name, 'figures', 'linewidths.pdf')
figure_22 = os.path.join(folder_name, 'figures', 'linewidths22.pdf')
figure_cross85 = os.path.join(folder_name, 'figures', 'linewidthscross85.pdf')
figure_23 = os.path.join(folder_name, 'figures', 'linewidths23.pdf')
figure_11 = os.path.join(folder_name, 'figures', 'linewidths11.pdf')
figure_cross87 = os.path.join(folder_name, 'figures', 'linewidthscross87.pdf')
figure_12 = os.path.join(folder_name, 'figures', 'linewidths12.pdf')

#fn.scatter3(intensities, peaks_11, peaks_cross, peaks_12, 'intensities', 'linewidth_11', 'linewidth_cross', 'linewidth_12', 'linewidth vs pump laser intensity', figure_name, save)
fn.scattering(np.sqrt(1+intensities/Isat22), peaks_22, '(1 + I/Isat)^.5', 'linewidth', 'peak22', figure_22, save)
fn.scattering(np.sqrt(1+intensities/Isatcross85), peaks_cross85, '(1 + I/Isat)^.5', 'linewidth', 'crosspeak85', figure_cross85, save)
fn.scattering(np.sqrt(1+intensities/Isat23), peaks_23, '(1 + I/Isat)^.5', 'linewidth', 'peak23', figure_23, save)
fn.scattering(np.sqrt(1+intensities/Isat11), peaks_11, '(1 + I/Isat)^.5', 'linewidth', 'peak11', figure_11, save)
fn.scattering(np.sqrt(1+intensities_cross87/Isatcross87), peaks_cross87, '(1 + I/Isat)^.5', 'linewidth', 'crosspeak87', figure_cross87, save)
fn.scattering(np.sqrt(1+intensities/Isat12), peaks_12, '(1 + I/Isat)^.5', 'linewidth', 'peak12', figure_12, save)
