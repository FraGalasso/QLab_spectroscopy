import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.preprocessing import MinMaxScaler
import os


data = pd.read_csv('data_temp/clean_data/mod_LD00001_calibrated_cropped_2.csv')
peaks = pd.read_csv('data_temp/clean_data/mod_LD00001_peaks.csv')

#####################################################################################################
# contstants and tunable by hand parameters

kb = 1.38e-23  # J/K
rb_mass = 1.67e-27 * 87   # kg
c = 3e8
term2 = kb / (c**2 * rb_mass)

lower_mask = (2.2 + 3.7711e5) * 1e9
upper_mask = (2.7 + 3.7711e5) * 1e9

# calib_correction = np.array([75, 251]) * 1e6

scale1_guess = 1.35
scale2_guess = 4

os.makedirs('data_temp/figures/temperature', exist_ok=True)
figure_name = f'new_calib_temp_fit_{lower_mask/1e9-3.7711e5:.2g}_{upper_mask/1e9-3.7711e5:.2g}.png'

#####################################################################################################

frequencies = data['frequencies'].to_numpy()
photodiode = data['volt_laser'].to_numpy()
# theorical peaks, with correction for calibration
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

# restricting region of fit
mask = (scaled_frequencies >= ((lower_mask - x_min) / scale_factor)) & (
    scaled_frequencies <= ((upper_mask - x_min) / scale_factor))
restricted_freq_scaled = scaled_frequencies[mask]
restricted_pd = photodiode[mask]


def optical_density_scaled(f, scale, temp, f0):
    return scale / np.sqrt(temp) * np.exp(- ((1-f/f0*scale_factor - x_min/f0)**2) / (2 * term2 * temp))


def transmission_temp_scaled(x, slope, intercept, scale1, scale2, temp):
    return (slope * x + intercept) * (np.exp(- optical_density_scaled(f=x, scale=scale1, temp=temp, f0=f_peaks[0])
                                             - optical_density_scaled(f=x, scale=scale2, temp=temp, f0=f_peaks[1])))


lower_bounds = [2 * slope_guess, 0.8 * int_guess, 0.5, 0.5, 273]
upper_bounds = [0, 1.2 * int_guess, np.inf, np.inf, 1000]
p0 = [slope_guess, int_guess, scale1_guess, scale2_guess, 300]

popt, pcov = curve_fit(transmission_temp_scaled, xdata=restricted_freq_scaled,
                       ydata=restricted_pd, bounds=(lower_bounds, upper_bounds), p0=p0, maxfev=10000)

f = np.linspace(min(frequencies), max(frequencies), 500)
f_sc = np.linspace(min(scaled_frequencies), max(scaled_frequencies), 500)
pd_fit = transmission_temp_scaled(f_sc, *popt)
pd_expected = transmission_temp_scaled(f_sc, *p0)

intercept_estimate = popt[1] - x_min / scale_factor * popt[0]
intercept_estimate_error = np.sqrt(
    pcov[1, 1] + pcov[0, 0] * (x_min / scale_factor)**2)

print(
    f'slope:\t\t{popt[0] / scale_factor} +/- {np.sqrt(pcov[0, 0]) / scale_factor} V/Hz')
print(
    f'intercept:\t{intercept_estimate} +/- {intercept_estimate_error} V')
print(f'scale1:\t\t{popt[2]} +/- {np.sqrt(pcov[2, 2])} K^0.5')
print(f'scale2:\t\t{popt[3]} +/- {np.sqrt(pcov[3, 3])} K^0.5')
print(f'temperature:\t{popt[4]} +/- {np.sqrt(pcov[4, 4])} K')

title = f'Fitting region: ({lower_mask/1e9-3.7711e5:.2g} - {upper_mask/1e9-3.7711e5:.2g} GHz) + 377.11 THz'

plt.figure()
plt.scatter(frequencies, photodiode, label='Data',
            color='blue', s=5, marker='.')
plt.plot(f, pd_fit, label=f'Fit result, T$={popt[4]:.1f}\\pm{np.sqrt(pcov[4, 4]):.1f}$K',
         color='red', linewidth=2)
plt.plot(f, pd_expected, label='Expected behaviour',
         color='purple', linewidth=2)
plt.axvline(f_peaks[0], color='black')
plt.axvline(f_peaks[1], color='black')
plt.axvline(lower_mask, color='green')
plt.axvline(upper_mask, color='green')
plt.xlabel('Frequencies [Hz]')
plt.ylabel('Photodiode readings [V]')
plt.title(title)
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig(f'data_temp/figures/temperature/{figure_name}')
plt.show()
