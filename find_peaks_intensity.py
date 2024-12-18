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


folder = 'data10'
title = 'intensity00010'
data_file = f'{folder}/clean_data/{title}.csv'

data = pd.read_csv(data_file)

timestamp = data['timestamp'].to_numpy()
volt_laser = data['volt_laser'].to_numpy()
volt_piezo = data['volt_piezo'].to_numpy()
volt_ld = data['volt_ld'].to_numpy()

lower_bounds = [-np.inf, -np.inf,
                0, 0, 0,
                4.5, 7.5, 8.6,
                0, 0, 0]

upper_bounds = [0, np.inf,
                np.inf, np.inf, np.inf,
                5.5, 8.5, 9.6,
                np.inf, np.inf, np.inf]

p0 = [-0.029, 1,
      0.63, 0.08, 0.25,
      5, 7.9, 9.1,
      1, 1, 1]

popt, pcov = curve_fit(transmission, volt_piezo,
                       volt_laser, p0=p0, bounds=(lower_bounds, upper_bounds), maxfev=10000)

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
peaks_indices, _ = find_peaks(residuals, height=0.01, distance=2000)
lor, cov = fp.fit_peaks_spectroscopy(
    volt_piezo, residuals, height=0.01, distance=2000)

# Correct peaks
peaks_indices = np.delete(peaks_indices, [0, 4, 5, 6, 10])

lor.pop(8)
lor.pop(4)
lor.pop(0)

cov.pop(8)
cov.pop(4)
cov.pop(0)

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
        377111224766631.8, 377111633094631.8, 377112041422631.8]

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
