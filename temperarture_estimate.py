import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.preprocessing import MinMaxScaler
import os
from functions import fit_residuals
from functools import partial


data = pd.read_csv(
    'data_temp_2/clean_data/fullspec00000_frequencies_cropped.csv')
peaks = pd.read_csv('data_temp_2/clean_data/fullspec00000_peaks.csv')

#####################################################################################################
# contstants and tunable by hand parameters

kb = 1.38e-23  # J/K
rb_mass = 1.67e-27 * 87   # kg
c = 3e8
term2 = kb / (c**2 * rb_mass)

lin_fit_lb = (7.1 + 3.7710e5) * 1e9
lin_fit_ub = (8.15 + 3.7710e5) * 1e9
region_3_lb = (0.3 + 3.7711e5) * 1e9
region_3_ub = (2.3 + 3.7711e5) * 1e9
region_2_lb = (1.8 + 3.7711e5) * 1e9
region_2_ub = (2.3 + 3.7711e5) * 1e9
region_1_lb = (0.8 + 3.7711e5) * 1e9
region_1_ub = (1.6 + 3.7711e5) * 1e9

scale1_guess = 1.35
scale2_guess = 4

os.makedirs('data_temp_2/figures/temperature', exist_ok=True)
os.makedirs('data_temp_2/figures/temperature/residuals', exist_ok=True)
# figure_name = f'new_calib_temp_fit_{lin_fit_lb/1e9-3.7711e5:.2g}_{lin_fit_ub/1e9-3.7711e5:.2g}.png'

#####################################################################################################

frequencies = data['frequencies'].to_numpy()
photodiode = data['photodiode'].to_numpy()
# theorical peaks, with correction for calibration
f_peaks = peaks['freq'].to_numpy()

# scaling frequencies in the range (0, 1)
# scaled_frequencies = (frequencies - min_frequencies) / (max_frequencies - min_frequencies)
scaler = MinMaxScaler()

frequencies_reshaped = frequencies.reshape(-1, 1)
fs = scaler.fit_transform(frequencies_reshaped)
scaled_frequencies = fs.flatten()

# relevant scaling parameters
scale_factor = (scaler.data_max_ - scaler.data_min_)[0]
x_min = (scaler.data_min_)[0]

# restricting region of fit
mask = (scaled_frequencies >= ((lin_fit_lb - x_min) / scale_factor)) & (
    scaled_frequencies <= ((lin_fit_ub - x_min) / scale_factor))
lin_region_freq_scaled = scaled_frequencies[mask]
lin_region_pd = photodiode[mask]


def linear(x, a, b):
    return a * x + b


popt, pcov = curve_fit(linear, lin_region_freq_scaled, lin_region_pd)
fixed_slope = popt[0]
fixed_intercept = popt[1]


x_lin = np.linspace(min(scaled_frequencies), max(scaled_frequencies), 100)
y_lin = linear(x_lin, *popt)
x_lin = np.linspace(min(frequencies), max(frequencies), 100)

title = f'Fitting region: ({lin_fit_lb/1e9-3.7710e5:.2g} - {lin_fit_ub/1e9-3.7710e5:.2g} GHz) + 377.11 THz'

plt.figure()
plt.scatter(frequencies, photodiode, label='Data',
            color='blue', s=5, marker='.')
plt.plot(x_lin, y_lin, label=f'Fit result', color='red', linewidth=2)
plt.axvline(lin_fit_lb, color='green')
plt.axvline(lin_fit_ub, color='green')
plt.xlabel('Frequencies [Hz]')
plt.ylabel('Photodiode readings [V]')
plt.title(title)
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('data_temp_2/figures/temperature/linear.png')
plt.close()


def optical_density_scaled(f, scale, temp, f0):
    return scale / np.sqrt(temp) * np.exp(- ((1-f/f0*scale_factor - x_min/f0)**2) / (2 * term2 * temp))


def transmission_temp_scaled(x, scale1, scale2, temp):
    return (fixed_slope * x + fixed_intercept) * (np.exp(- optical_density_scaled(f=x, scale=scale1, temp=temp, f0=f_peaks[-2])
                                                         - optical_density_scaled(f=x, scale=scale2, temp=temp, f0=f_peaks[-1])))


lower_bounds = [0.5, 0.5, 273]
upper_bounds = [np.inf, np.inf, 1000]
p0 = [scale1_guess, scale2_guess, 330]

# restricting region of fit
mask = (scaled_frequencies >= ((region_3_lb - x_min) / scale_factor)) & (
    scaled_frequencies <= ((region_3_ub - x_min) / scale_factor))
restricted_freq_scaled = scaled_frequencies[mask]
restricted_pd = photodiode[mask]

popt, pcov = curve_fit(transmission_temp_scaled, xdata=restricted_freq_scaled,
                       ydata=restricted_pd, bounds=(lower_bounds, upper_bounds), p0=p0, maxfev=10000)
fixed_scale1 = popt[0]

f = np.linspace(min(frequencies), max(frequencies), 500)
f_sc = np.linspace(min(scaled_frequencies), max(scaled_frequencies), 500)
pd_fit = transmission_temp_scaled(f_sc, *popt)

print('Fitting region 3:')
print(f'scale1:\t\t{popt[0]} +/- {np.sqrt(pcov[0, 0])} K^0.5')
print(f'scale2:\t\t{popt[1]} +/- {np.sqrt(pcov[1, 1])} K^0.5')
print(f'temperature:\t{popt[2]} +/- {np.sqrt(pcov[2, 2])} K\n')

title = f'Fitting region: ({region_3_lb/1e9-3.7711e5:.2g} - {region_3_ub/1e9-3.7711e5:.2g} GHz) + 377.11 THz'

plt.figure()
plt.scatter(frequencies, photodiode, label='Data',
            color='blue', s=5, marker='.')
plt.plot(f, pd_fit, label=f'Fit result, T$={popt[2]:.1f}\\pm{np.sqrt(pcov[2, 2]):.1f}$K',
         color='red', linewidth=2)
plt.axvline(f_peaks[-2], color='black')
plt.axvline(f_peaks[-1], color='black')
plt.axvline(region_3_lb, color='green')
plt.axvline(region_3_ub, color='green')
plt.xlabel('Frequencies [Hz]')
plt.ylabel('Photodiode readings [V]')
plt.title(title)
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('data_temp_2/figures/temperature/region_3.png')
plt.close()

fit_residuals(func=transmission_temp_scaled, x=restricted_freq_scaled, y=restricted_pd, params=popt,
              x_label='Frequencies [Hz]', y_label='Residuals [V]', title='', file_name=f'data_temp_2/figures/temperature/residuals/resid_3.png', save=True)


def transmission_temp_scaled_region_2(
    x, scale2, temp): return transmission_temp_scaled(x, fixed_scale1, scale2, temp)


lower_bounds = [0.5, 273]
upper_bounds = [np.inf, 1000]
p0 = [scale2_guess, 330]

# restricting region of fit
mask = (scaled_frequencies >= ((region_2_lb - x_min) / scale_factor)) & (
    scaled_frequencies <= ((region_2_ub - x_min) / scale_factor))
restricted_freq_scaled = scaled_frequencies[mask]
restricted_pd = photodiode[mask]

popt, pcov = curve_fit(transmission_temp_scaled_region_2, xdata=restricted_freq_scaled,
                       ydata=restricted_pd, bounds=(lower_bounds, upper_bounds), p0=p0, maxfev=10000)
fixed_scale2 = popt[0]

f = np.linspace(min(frequencies), max(frequencies), 500)
f_sc = np.linspace(min(scaled_frequencies), max(scaled_frequencies), 500)
pd_fit = transmission_temp_scaled_region_2(f_sc, *popt)

print('Fitting region 2:')
print(f'scale2:\t\t{popt[0]} +/- {np.sqrt(pcov[0, 0])} K^0.5')
print(f'temperature:\t{popt[1]} +/- {np.sqrt(pcov[1, 1])} K\n')

title = f'Fitting region: ({region_2_lb/1e9-3.7711e5:.2g} - {region_2_ub/1e9-3.7711e5:.2g} GHz) + 377.11 THz'

plt.figure()
plt.rcParams.update({'font.size': 10})
plt.scatter(frequencies, photodiode, label='Data',
            color='blue', s=5, marker='.')
plt.plot(f, pd_fit, label=f'Fit result, T$={popt[1]:.1f}\\pm{np.sqrt(pcov[1, 1]):.1f}$K',
         color='red', linewidth=2)
plt.axvline(f_peaks[-2], color='black')
plt.axvline(f_peaks[-1], color='black')
plt.axvline(region_2_lb, color='green')
plt.axvline(region_2_ub, color='green')
plt.xlabel('Frequencies [Hz]')
plt.ylabel('Photodiode readings [V]')
plt.title(title)
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('data_temp_2/figures/temperature/region_2.png')
plt.show()

fit_residuals(func=transmission_temp_scaled_region_2, x=restricted_freq_scaled, y=restricted_pd, params=popt,
              x_label='Frequencies [Hz]', y_label='Residuals [V]', title='', file_name=f'data_temp_2/figures/temperature/residuals/resid_2.png', save=True)


def transmission_temp_scaled_region_1(
    x, scale1, temp): return transmission_temp_scaled(x, scale1, fixed_scale2, temp)


lower_bounds = [0.5, 273]
upper_bounds = [np.inf, 1000]
p0 = [scale1_guess, 330]

# restricting region of fit
mask = (scaled_frequencies >= ((region_1_lb - x_min) / scale_factor)) & (
    scaled_frequencies <= ((region_1_ub - x_min) / scale_factor))
restricted_freq_scaled = scaled_frequencies[mask]
restricted_pd = photodiode[mask]

popt, pcov = curve_fit(transmission_temp_scaled_region_1, xdata=restricted_freq_scaled,
                       ydata=restricted_pd, bounds=(lower_bounds, upper_bounds), p0=p0, maxfev=10000)
fixed_scale1 = popt[0]


f = np.linspace(min(frequencies), max(frequencies), 500)
f_sc = np.linspace(min(scaled_frequencies), max(scaled_frequencies), 500)
pd_fit = transmission_temp_scaled_region_1(f_sc, *popt)

print('Fitting region 1:')
print(f'scale1:\t\t{popt[0]} +/- {np.sqrt(pcov[0, 0])} K^0.5')
print(f'temperature:\t{popt[1]} +/- {np.sqrt(pcov[1, 1])} K\n')

title = f'Fitting region: ({region_1_lb/1e9-3.7711e5:.2g} - {region_1_ub/1e9-3.7711e5:.2g} GHz) + 377.11 THz'

plt.figure()
plt.rcParams.update({'font.size': 10})
plt.scatter(frequencies, photodiode, label='Data',
            color='blue', s=5, marker='.')
plt.plot(f, pd_fit, label=f'Fit result, T$={popt[1]:.1f}\\pm{np.sqrt(pcov[1, 1]):.1f}$K',
         color='red', linewidth=2)
plt.axvline(f_peaks[-2], color='black')
plt.axvline(f_peaks[-1], color='black')
plt.axvline(region_1_lb, color='green')
plt.axvline(region_1_ub, color='green')
plt.xlabel('Frequencies [Hz]')
plt.ylabel('Photodiode readings [V]')
plt.title(title)
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('data_temp_2/figures/temperature/region_1.png')
plt.show()

fit_residuals(func=transmission_temp_scaled_region_1, x=restricted_freq_scaled, y=restricted_pd, params=popt,
              x_label='Frequencies [Hz]', y_label='Residuals [V]', title='', file_name=f'data_temp_2/figures/temperature/residuals/resid_1.png', save=True)


def transmission_temp_scaled_only_temp(
    x, temp): return transmission_temp_scaled(x, fixed_scale1, fixed_scale2, temp)


lower_bounds = [273]
upper_bounds = [1000]
p0 = [330]

# restricting region of fit
mask = (scaled_frequencies >= ((region_3_lb - x_min) / scale_factor)) & (
    scaled_frequencies <= ((region_3_ub - x_min) / scale_factor))
restricted_freq_scaled = scaled_frequencies[mask]
restricted_pd = photodiode[mask]

popt, pcov = curve_fit(transmission_temp_scaled_only_temp, xdata=restricted_freq_scaled,
                       ydata=restricted_pd, bounds=(lower_bounds, upper_bounds), p0=p0, maxfev=10000)

f = np.linspace(min(frequencies), max(frequencies), 500)
f_sc = np.linspace(min(scaled_frequencies), max(scaled_frequencies), 500)
pd_fit = transmission_temp_scaled_only_temp(f_sc, *popt)

print('Fitting region 3, only with temperature:')
print(f'temperature:\t{popt[0]} +/- {np.sqrt(pcov[0, 0])} K\n')

title = f'Fitting region: ({region_3_lb/1e9-3.7711e5:.2g} - {region_3_ub/1e9-3.7711e5:.2g} GHz) + 377.11 THz'

plt.figure()
plt.rcParams.update({'font.size': 10})
plt.scatter(frequencies, photodiode, label='Data',
            color='blue', s=5, marker='.')
plt.plot(f, pd_fit, label=f'Fit result, T$={popt[0]:.1f}\\pm{np.sqrt(pcov[0, 0]):.1f}$K',
         color='red', linewidth=2)
plt.axvline(f_peaks[-2], color='black')
plt.axvline(f_peaks[-1], color='black')
plt.axvline(region_3_lb, color='green')
plt.axvline(region_3_ub, color='green')
plt.xlabel('Frequencies [Hz]')
plt.ylabel('Photodiode readings [V]')
plt.title(title)
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('data_temp_2/figures/temperature/only_temp.png')
plt.show()

fit_residuals(func=transmission_temp_scaled_only_temp, x=restricted_freq_scaled, y=restricted_pd, params=popt,
              x_label='Frequencies [Hz]', y_label='Residuals [V]', title='', file_name=f'data_temp_2/figures/temperature/residuals/resid_temp.png', save=True)
