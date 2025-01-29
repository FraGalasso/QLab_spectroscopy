import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.preprocessing import MinMaxScaler
import os
from functions import fit_residuals


data = pd.read_csv('data_temp/clean_data/mod_LD00001_calibrated_cropped_2.csv')
peaks = pd.read_csv('data_temp/clean_data/mod_LD00001_peaks.csv')

#####################################################################################################
# contstants and tunable by hand parameters

kb = 1.38e-23  # J/K
rb_mass = 1.67e-27 * 87   # kg
c = 3e8
term2 = kb / (c**2 * rb_mass)

scale1_guess = 1.35
scale2_guess = 4

os.makedirs('data_temp/figures/temperature', exist_ok=True)

#####################################################################################################

frequencies = data['frequencies'].to_numpy()
photodiode = data['volt_laser'].to_numpy()
f_peaks = peaks['freq']

frequencies = frequencies[1200:]
photodiode = photodiode[1200:]

# scaling frequencies in the range (0, 1)
# scaled_frequencies = (frequencies - min_frequencies) / (max_frequencies - min_frequencies)
scaler = MinMaxScaler()

frequencies_reshaped = frequencies.reshape(-1, 1)
fs = scaler.fit_transform(frequencies_reshaped)
scaled_frequencies = fs.flatten()

# relevant scaling parameters
scale_factor = (scaler.data_max_ - scaler.data_min_)[0]
x_min = (scaler.data_min_)[0]

# guesses
slope_guess = (photodiode[-1] - photodiode[0]) / \
    (scaled_frequencies[-1] - scaled_frequencies[0])
int_guess = photodiode[-1] - slope_guess * scaled_frequencies[0] + 0.25
lb_int = -0.8 * slope_guess * photodiode[-1]


def optical_density_scaled(f, scale, temp, f0):
    return scale / np.sqrt(temp) * np.exp(- ((1-f/f0*scale_factor - x_min/f0)**2) / (2 * term2 * temp))


def transmission_temp_scaled(x, slope, intercept, scale1, scale2, temp):
    return (slope * x + intercept) * (np.exp(- optical_density_scaled(f=x, scale=scale1, temp=temp, f0=f_peaks[0])
                                             - optical_density_scaled(f=x, scale=scale2, temp=temp, f0=f_peaks[1])))


lower_bounds = [2 * slope_guess, 0.8 * int_guess, 0.5, 0.5, 273]
upper_bounds = [0, 1.2 * int_guess, np.inf, np.inf, 1000]
p0 = [slope_guess, int_guess, scale1_guess, scale2_guess, 300]

T_list_low = []
dT_list_low = []

bound_list_low = (np.arange(0.3, 2.4, 0.1) + 3.7711e5) * 1e9

for lm in bound_list_low:
    lower_mask = lm
    upper_mask = (2.9 + 3.7711e5) * 1e9

    mask = (scaled_frequencies >= ((lower_mask - x_min) / scale_factor)) & (
        scaled_frequencies <= ((upper_mask - x_min) / scale_factor))
    restricted_freq_scaled = scaled_frequencies[mask]
    restricted_pd = photodiode[mask]

    popt, pcov = curve_fit(transmission_temp_scaled, xdata=restricted_freq_scaled,
                           ydata=restricted_pd, bounds=(lower_bounds, upper_bounds), p0=p0, maxfev=10000)
    
    T_list_low.append(popt[4])
    dT_list_low.append(pcov[4, 4])

T_list_up = []
dT_list_up = []

bound_list_up = (np.arange(2.6, 3.3, 0.1) + 3.7711e5) * 1e9

for ub in bound_list_up:
    lower_mask = f_peaks[1]
    upper_mask = ub

    mask = (scaled_frequencies >= ((lower_mask - x_min) / scale_factor)) & (
        scaled_frequencies <= ((upper_mask - x_min) / scale_factor))
    restricted_freq_scaled = scaled_frequencies[mask]
    restricted_pd = photodiode[mask]

    popt, pcov = curve_fit(transmission_temp_scaled, xdata=restricted_freq_scaled,
                           ydata=restricted_pd, bounds=(lower_bounds, upper_bounds), p0=p0, maxfev=10000)
    
    T_list_up.append(popt[4])
    dT_list_up.append(pcov[4, 4])

plt.figure()
plt.errorbar(x=bound_list_low, y=T_list_low, yerr=dT_list_low, fmt='.', color='blue', ls='')
plt.xlabel('Lower bound [Hz]')
plt.ylabel('T [K]')
plt.title('Estimated temperature with upper bound at (2.9 + 3.7711e5) GHz')
plt.grid()
plt.tight_layout()
plt.savefig(f'data_temp/figures/temperature/lower_bound_sweep.png')
plt.close()



plt.figure()
plt.errorbar(x=bound_list_up, y=T_list_up, yerr=dT_list_up, fmt='.', color='blue', ls='')
plt.xlabel('Upper bound [Hz]')
plt.ylabel('T [K]')
plt.title('Estimated temperature with lower bound at second transition')
plt.grid()
plt.tight_layout()
plt.savefig(f'data_temp/figures/temperature/upper_bound_sweep_2.png')
plt.close()
