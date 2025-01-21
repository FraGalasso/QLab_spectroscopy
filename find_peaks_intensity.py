import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import fit_peaks as fp
import functions as fn


def transmission(x, slope, intercept, scale1, scale2, scale3, mean1, mean2, mean3, sigma1, sigma2, sigma3):
    return (slope * x + intercept) * (np.exp(-scale1 * np.exp(-(x - mean1)**2 / (2 * (sigma1**2))) -
                                             scale2 * np.exp(-(x - mean2)**2 / (2 * (sigma2**2))) -
                                             scale3 * np.exp(-(x - mean3)**2 / (2 * (sigma3**2)))))


#####################################################################################################
# section with file specific values

# use this to find peaks in intensity
# folder = 'data_linewidth'
# title = 'intensity00010'

# use this to find peaks in modld files
folder = 'data_temp'
title = 'mod_LD00000'
data_file = f'{folder}/clean_data/{title}.csv'

min_height = 0.017
min_distance = 200

# other stuff you have to change in the rest of the file is:
# - initial guesses and bounds for fit (lines 45-58)
# - which peaks to remove (lines 88-99)

#####################################################################################################


data = pd.read_csv(data_file)

timestamp = data['timestamp'].to_numpy()
volt_laser = data['volt_laser'].to_numpy()
volt_piezo = data['volt_piezo'].to_numpy()
volt_ld = data['volt_ld'].to_numpy()

lower_bounds = [-np.inf, 0,
                0, 0, 0,
                4.3, 6.5, 7.5,
                0, 0, 0]

upper_bounds = [0, np.inf,
                np.inf, np.inf, np.inf,
                6, 8, 9,
                np.inf, np.inf, np.inf]

p0 = [-0.072, 1.796,
      0.5, 0.07, 0.18,
      5.1, 7.3, 8.2,
      1, 1, 1]

popt, pcov = curve_fit(transmission, volt_piezo[700:],
                       volt_laser[700:], p0=p0, bounds=(lower_bounds, upper_bounds), maxfev=10000)

print(popt)

vp = np.linspace(min(volt_piezo), max(volt_piezo), 500)
fit = transmission(vp, *popt)

plt.figure()
plt.plot(volt_piezo, volt_laser, label='Laser intensity',
         color='blue', markersize=5, marker='.')
plt.plot(vp, fit, label='Fit result',
         color='red', markersize=5, marker='')
plt.xlabel('Volt Piezo [V]')
plt.ylabel('Volt Laser [V]')
plt.title(f'Fitting {title}')
plt.grid()
plt.legend()
plt.savefig(f'{folder}/figures/find_peaks/first_fit{title}.png')
plt.close()


residuals = volt_laser - transmission(volt_piezo, *popt)
peaks_indices, _ = find_peaks(
    residuals, height=min_height, distance=min_distance)
lor, cov = fp.fit_peaks_spectroscopy(
    volt_piezo, residuals, height=min_height, distance=min_distance)

# Correct peaks
peaks_indices = np.delete(peaks_indices, [3, 4, 5, 6])

lor.pop(6)
lor.pop(5)
lor.pop(4)
lor.pop(3)

cov.pop(6)
cov.pop(5)
cov.pop(4)
cov.pop(3)

piezo_peaks = np.array(volt_piezo[peaks_indices])
timestamp = np.array(timestamp[peaks_indices])
laser_peaks = np.array(volt_laser[peaks_indices])
ld_peaks = np.array(volt_ld[peaks_indices])
y_peaks = np.array(residuals[peaks_indices])

x0_list = []
A_list = []
gamma_list = []
off_list = []
dx0 = []
dA = []
dgamma = []
doff = []

for popt, pcov in zip(lor, cov):
    A_list.append(popt[0])
    x0_list.append(popt[1])
    gamma_list.append(popt[2])
    off_list.append(popt[3])
    dA.append(np.sqrt(pcov[0, 0]))
    dx0.append(np.sqrt(pcov[1, 1]))
    dgamma.append(np.sqrt(pcov[2, 2]))
    doff.append(np.sqrt(pcov[3, 3]))

x0_list = np.array(x0_list)
A_list = np.array(A_list)
gamma_list = np.array(gamma_list)
off_list = np.array(off_list)
dx0 = np.array(dx0)
dA = np.array(dA)
dgamma = np.array(dgamma)
doff = np.array(doff)

print(x0_list)

fn.plot_piezo_laser_fit(volt_piezo, residuals, f'{folder}/figures/find_peaks/residuals_{title}.png',
                        A_list, x0_list, gamma_list, off_list, piezo_peaks, y_peaks, save=True)

freq = [377108945610922.8, 377109126401922.8, 377109307192922.8,
        377111224766631.8, 377112041422631.8]

# Saving data in clean_data folder
output_file = f'{folder}/clean_data/{title}_peaks.csv'
df = pd.DataFrame()
df['indices'] = peaks_indices
df['timestamp'] = timestamp
df['pd_peaks'] = laser_peaks
df['piezo_peaks'] = piezo_peaks
df['ld_peaks'] = ld_peaks
df['freq'] = freq
df['lor_A'] = A_list
df['lor_mean'] = x0_list
df['lor_gamma'] = gamma_list
df['lor_off'] = off_list
df['lor_d_A'] = dA
df['lor_d_mean'] = dx0
df['lor_d_gamma'] = dgamma
df['lor_d_off'] = doff

df.to_csv(output_file, index=False)
print(f"Data saved to {output_file}")
