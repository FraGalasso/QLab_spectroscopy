import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
# from sklearn.preprocessing import MinMaxScaler


data = pd.read_csv('data_temp/clean_data/mod_LD00001_calibrated_cropped.csv')
peaks = pd.read_csv('data_temp/clean_data/mod_LD00001_peaks.csv')

frequencies = data['frequencies'].to_numpy()/1e9
photodiode = data['volt_laser'].to_numpy()
f_peaks = (peaks['freq']+np.array([75, 251])*1e6)/1e9

frequencies = frequencies[1200:]
photodiode = photodiode[1200:]

kb = 1.38e-23  # J/K
rb_mass = 1.67e-27 * 87   # kg
c = 3e8
term2 = kb / (c**2 * rb_mass)

'''scaler = MinMaxScaler()

frequencies_reshaped = frequencies.reshape(-1, 1)
fs = scaler.fit_transform(frequencies_reshaped)
scaled_frequencies = fs.flatten()

scale_factor = scaler.data_max_ - scaler.data_min_
x_min = scaler.data_min_

peaks['scaled_freq'] = (peaks['freq'] - x_min) / scale_factor

mask = (scaled_frequencies >= ((lower_mask - x_min) / scale_factor)) & (
    scaled_frequencies <= ((upper_mask - x_min) / scale_factor))
restricted_freq_scaled = scaled_frequencies[mask]'''


def optical_density(f, scale, temp, f0):
    return scale / np.sqrt(temp) * np.exp(- ((1-f/f0)**2) / (2 * term2 * temp))


def transmission_temp(x, slope, intercept, scale1, scale2, temp):
    return (slope * x + intercept) * (np.exp(- optical_density(f=x, scale=scale1, temp=temp, f0=f_peaks[0])
                                             - optical_density(f=x, scale=scale2, temp=temp, f0=f_peaks[1])))


lower_mask = 2.55 + 3.7711e5
upper_mask = 2.85 + 3.7711e5
# lower_mask = f_peaks[1]
# upper_mask = max(frequencies)

mask = (frequencies >= lower_mask) & (frequencies <= upper_mask)
restricted_freq = frequencies[mask]
restricted_pd = photodiode[mask]

slope_guess = (photodiode[-1]-photodiode[0]) / (frequencies[-1]-frequencies[0])
int_guess = photodiode[-1] - slope_guess * frequencies[0]
lb_int = -0.8 * slope_guess * photodiode[-1]

lower_bounds = [2 * slope_guess, 0.8 * int_guess, 0.5, 0.5, 273]
upper_bounds = [0, 1.2 * int_guess, np.inf, np.inf, 1000]
p0 = [slope_guess, int_guess, 1.21, 3.22, 300]

popt, pcov = curve_fit(transmission_temp, xdata=restricted_freq,
                       ydata=restricted_pd, bounds=(lower_bounds, upper_bounds), p0=p0, maxfev=10000)

f = np.linspace(min(frequencies), max(frequencies), 500)
pd_fit = transmission_temp(f, *popt)

print(f'slope:\t{popt[0]} +/- {np.sqrt(pcov[0, 0])} V/Hz')
print(f'intercept:\t{popt[1]} +/- {np.sqrt(pcov[1, 1])} V')
print(f'scale1:\t\t{popt[2]} +/- {np.sqrt(pcov[2, 2])} K^0.5')
print(f'scale2:\t\t{popt[3]} +/- {np.sqrt(pcov[3, 3])} K^0.5')
print(f'temperature:\t{popt[4]} +/- {np.sqrt(pcov[4, 4])} K')


plt.figure()
plt.scatter(frequencies, photodiode, label='Data',
            color='blue', s=5, marker='.')
plt.plot(f, pd_fit, label=f'Fit result, T$={popt[4]:.1f}\\pm{np.sqrt(pcov[4, 4]):.1f}$K',
         color='red', linewidth=2)
plt.axvline(f_peaks[0], color='black')
plt.axvline(f_peaks[1], color='black')
plt.axvline(lower_mask, color='green')
plt.axvline(upper_mask, color='green')
plt.xlabel('Frequencies [Hz]')
plt.ylabel('Photodiode readings [V]')
plt.title('')
plt.grid()
plt.legend()
plt.tight_layout()
# plt.savefig('data9/figures/temperature/temp_fit_shoulder.pdf')
plt.show()
