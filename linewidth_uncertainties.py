import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import linewidth_functions as lfn
from uncertainties import ufloat, UFloat
from uncertainties.unumpy import nominal_values, std_devs, uarray

# file settings
save = True
plotting_spec = True
fit = 'linear'  # choose either 'linear' or 'sqrt' fit

folder_name = 'data_linewidth'
figure_folder = os.path.join(folder_name, 'figures_unc')
lfn.ensure_dir(figure_folder)

# constants:
Isat85 = 4.4876  # mW/cm2
Isat87 = 4.484  # mW/cm2
gamma87 = 2 * np.pi * 5.75  # MHz
gamma85 = 2 * np.pi * 5.75  # MHz

# XDATA: I/Isat:

# Input data
powers = [0.079, 0.07, 0.1, 0.221, 0.092, 0.104,
          0.442, 1.009, 0.313, 0.7, 0.295]  # mW
powers_cross87 = [0.079, 0.07, 0.1, 0.221,
                  0.092, 0.104, 0.442, 0.313, 0.295]  # mW

# compute arrays containing both data and errors: 5% uncertainty in power readings
powers_err = [ufloat(p, 0.05 * p) for p in powers]
powers_cross87_err = [ufloat(p, 0.05 * p) for p in powers_cross87]

# Calculate the beam area with uncertainty
radius = ufloat(0.08, 0.01)  # cm with uncertainty ±0.01
area_beam = np.pi * (radius ** 2)  # cm^2 with propagated uncertainty

# Convert powers to int
int87 = np.array(powers_err) / (area_beam)
int85 = np.array(powers_err) / (area_beam)
int_cross87 = np.array(powers_cross87_err) / (area_beam)

int87_x = gamma87 * (1 + int87/Isat87)**0.5
intcross87_x = gamma87 * (1 + int_cross87/Isat87)**0.5
int85_x = gamma85 * (1 + int85/Isat85)**0.5
intcross85_x = gamma85 * (1 + int85/Isat85)**0.5

# YDATA:

# generate filepaths:
file_names = [f"intensity0000{i}" for i in range(10)] + ["intensity00010"]
file_paths = [
    f'{folder_name}/clean_data/{file_name}.csv' for file_name in file_names]
file_peaks = [
    f'{folder_name}/clean_data/{file_name}_peaks.csv' for file_name in file_names]

# calculate calibration

data = pd.read_csv(file_paths[0], skiprows=1, names=[
                   'timestamp', 'photodiode', 'volt_piezo', 'volt_ld'])
peaks = pd.read_csv(file_peaks[0], skiprows=1, names=['indices', 'timestamp', 'pd_peaks', 'piezo_peaks',
                    'freq', 'lor_A', 'lor_mean', 'lor_gamma', 'lor_off', 'lor_d_A', 'lor_d_mean', 'lor_d_gamma', 'lor_d_off'])

peaks_mean = uarray(peaks['lor_mean'].to_numpy(),
                    peaks['lor_d_mean'].to_numpy())
peaks_freqs = peaks['freq'].to_numpy()

guesses = [1, 1]

coeff1, __ = lfn.plot_ufloat_fit(peaks_mean, peaks_freqs, lfn.lin_model, 'peaks_piezovolt', 'freq',
                                 'Fixed diode current, modulated piezo position: time calib', guesses, os.path.join(figure_folder, "calibration.png"), save)

# create storage arrays
peaks_22 = []
peaks_cross85 = []
peaks_23 = []
peaks_11 = []
peaks_cross87 = []
peaks_12 = []

# loop through files
# for each file: reads data, finds calibrated frequencies, plots, and appends peaks to peak arrays
for i, (file_name, file_path, file_peak) in enumerate(zip(file_names, file_paths, file_peaks)):
    # Read the CSV file, and specify the data types
    data = pd.read_csv(file_path, skiprows=1, names=[
                       'timestamp', 'photodiode', 'volt_piezo', 'volt_ld'])
    peaks = pd.read_csv(file_peak, skiprows=1, names=['indices', 'timestamp', 'pd_peaks', 'piezo_peaks', 'freq',
                        'lor_A', 'lor_mean', 'lor_gamma', 'lor_off', 'lor_d_A', 'lor_d_mean', 'lor_d_gamma', 'lor_d_off'])

    # Create arrays with uncertainties
    peaks_gamma = uarray(peaks['lor_gamma'].to_numpy(),
                         peaks['lor_d_gamma'].to_numpy())
    peaks_mean = uarray(peaks['lor_mean'].to_numpy(),
                        peaks['lor_d_mean'].to_numpy())
    peaks_freqs = peaks['freq'].to_numpy()

    # Process peaks_gamma: go from HWHM to FWHM and calibrate
    peaks_gamma = peaks_gamma * 2  # Convert from HWHM to FWHM
    peaks_gamma = coeff1[0] * peaks_gamma * \
        1e-6  # Calibrate and transform to MHz

    # append peak values
    peaks_22.append(peaks_gamma[0])
    peaks_cross85.append(peaks_gamma[1])
    peaks_23.append(peaks_gamma[2])
    peaks_11.append(peaks_gamma[3])

    if len(peaks) == 5:
        peaks_12.append(peaks_gamma[4])
    else:
        peaks_cross87.append(peaks_gamma[4])
        peaks_12.append(peaks_gamma[5])

    if plotting_spec:
        figure_path = os.path.join(
            figure_folder, f"{powers[i]:.3g}_{file_name}.png")
        spec_data = data['photodiode'].to_numpy()
        frequencies = lfn.calibrate(data['volt_piezo'], coeff1)
        lfn.plot_ufloat(
            frequencies,
            spec_data,
            'Frequency [Hz]',
            'Photodiode Intensity [V]',
            f'Spectroscopy with Pump Laser Power = {powers[i]} mW, Intensity = {int85[i]} mW/cm^2',
            figure_path,
            save
        )

# names for the figures
figure_name = os.path.join(figure_folder, 'linewidths.png')
figure_22 = os.path.join(figure_folder, 'linewidths22.png')
figure_cross85 = os.path.join(figure_folder, 'linewidthscross85.png')
figure_23 = os.path.join(figure_folder, 'linewidths23.png')

figure_11 = os.path.join(figure_folder, 'linewidths11.png')
figure_cross87 = os.path.join(figure_folder, 'linewidthscross87.png')
figure_12 = os.path.join(figure_folder, 'linewidths12.png')

x_label = r'$\Gamma \cdot \sqrt{ 1 + \left( \frac{I}{I_{\text{sat}}} \right) }$'
y_label = r'$\text{Linewidth [MHz]}$'

coeff22, dcoeff22 = lfn.plot_ufloat_fit(int85_x, peaks_22, lfn.lin_model, x_label, y_label,
                             r'Peak 2-2 $\text{Rb}_{85}$', guesses, figure_22, save)
coeffcross85, dcoeffcross85 = lfn.plot_ufloat_fit(intcross85_x, peaks_cross85, lfn.lin_model, x_label, y_label,
                             r'Crosspeak $\text{Rb}_{85}$', guesses, figure_cross85, save)
coeff23, dcoeff23 = lfn.plot_ufloat_fit(int85_x, peaks_23, lfn.lin_model, x_label, y_label,
                             r'Peak 2-3 $\text{Rb}_{85}$', guesses, figure_23, save)

coeff11, dcoeff11 = lfn.plot_ufloat_fit(int87_x, peaks_11, lfn.lin_model, x_label, y_label,
                             r'Peak 1-1 $\text{Rb}_{87}$', guesses, figure_11, save)
coeffcross87, dcoeffcross87 = lfn.plot_ufloat_fit(intcross87_x, peaks_cross87, lfn.lin_model, x_label, y_label,
                             r'Crosspeak $\text{Rb}_{87}$', guesses, figure_cross87, save)
coeff12, dcoeff12 = lfn.plot_ufloat_fit(int87_x, peaks_12, lfn.lin_model, x_label, y_label,
                             r'Peak 1-2 $\text{Rb}_{87}$', guesses, figure_12, save)

plt.figure(figsize=(12, 10))
plt.rcParams.update({'font.size': 16})
plt.errorbar(nominal_values(int87_x), nominal_values(peaks_11)-coeff11[1], xerr=std_devs(int87_x), yerr=std_devs(peaks_11), fmt='.', color='green', label=f'F=1$\\to$F\'=1 - ({coeff11[1]:.3g} ± {dcoeff11[1]:.3g}) V', capsize=5)
plt.errorbar(nominal_values(int87_x), nominal_values(peaks_12)-coeff12[1], xerr=std_devs(int87_x), yerr=std_devs(peaks_12), fmt='.', color='blue', label=f'F=1$\\to$F\'=2 - ({coeff12[1]:.3g} ± {dcoeff12[1]:.3g}) V', capsize=5)
plt.errorbar(nominal_values(intcross87_x), nominal_values(peaks_cross87)-coeffcross87[1], xerr=std_devs(intcross87_x), yerr=std_devs(peaks_cross87), fmt='.', color='red', label=f'Crosspeak - ({coeffcross87[1]:.3g} ± {dcoeffcross85[1]:.3g}) V', capsize=5)

fitstr11 = f'y = ({coeff11[0]:.3g} ± {dcoeff11[0]:.3g}) * x'
fitstr12 = f'y = ({coeff12[0]:.3g} ± {dcoeff12[0]:.3g}) * x'
fitstrcross87 = f'y = ({coeffcross87[0]:.3g} ± {dcoeffcross87[0]:.3g}) * x'
        
x_fit = np.linspace(min(nominal_values(int87_x)), max(nominal_values(int87_x)), 100)
y_fit11 = coeff11[0] * x_fit
y_fit12 = coeff12[0] * x_fit
y_fitcross87 = coeffcross87[0] * x_fit
plt.plot(x_fit, y_fit11, label=fitstr11, color='green', linestyle='--')
plt.plot(x_fit, y_fit12, label=fitstr12, color='blue', linestyle='--')
plt.plot(x_fit, y_fitcross87, label=fitstrcross87, color='red', linestyle='--')

plt.title(r'$^{87}$Rb power')
plt.xlabel(x_label)
plt.ylabel(y_label)

handles, labels = plt.gca().get_legend_handles_labels()

plt.xticks(rotation=45)
plt.grid()
plt.tight_layout()
plt.savefig('data_linewidth/figures_unc/rb87_1.png')
plt.show()

fig_legend = plt.figure(figsize=(3, 2))
ax = fig_legend.add_subplot(111)
ax.axis('off')
legend = ax.legend(handles, labels, loc='center', fontsize=10, frameon=False)  

fig_legend.savefig('data_linewidth/figures_unc/legend_87.png', bbox_inches='tight', dpi=300)

plt.figure(figsize=(12, 10))
plt.rcParams.update({'font.size': 16})
plt.errorbar(nominal_values(int85_x), nominal_values(peaks_22)-coeff22[1], xerr=std_devs(int85_x), yerr=std_devs(peaks_22), fmt='.', color='green', label=f'F=2$\\to$F\'=2 - ({coeff22[1]:.3g} ± {dcoeff22[1]:.3g}) V', capsize=5)
plt.errorbar(nominal_values(int85_x), nominal_values(peaks_23)-coeff23[1], xerr=std_devs(int85_x), yerr=std_devs(peaks_23), fmt='.', color='blue', label=f'F=2$\\to$F\'=3 - ({coeff23[1]:.3g} ± {dcoeff23[1]:.3g}) V', capsize=5)
plt.errorbar(nominal_values(intcross85_x), nominal_values(peaks_cross85)-coeffcross85[1], xerr=std_devs(intcross85_x), yerr=std_devs(peaks_cross85), fmt='.', color='red', label=f'Crosspeak - ({coeffcross85[1]:.3g} ± {dcoeffcross85[1]:.3g}) V', capsize=5)

fitstr22 = f'y = ({coeff22[0]:.3g} ± {dcoeff22[0]:.3g}) * x'
fitstr23 = f'y = ({coeff23[0]:.3g} ± {dcoeff23[0]:.3g}) * x'
fitstrcross85 = f'y = ({coeffcross85[0]:.3g} ± {dcoeffcross85[0]:.3g}) * x'
        
x_fit = np.linspace(min(nominal_values(int85_x)), max(nominal_values(int85_x)), 100)
y_fit22 = coeff22[0] * x_fit
y_fit23 = coeff23[0] * x_fit
y_fitcross85 = coeffcross85[0] * x_fit
plt.plot(x_fit, y_fit22, label=fitstr22, color='green', linestyle='--')
plt.plot(x_fit, y_fit23, label=fitstr23, color='blue', linestyle='--')
plt.plot(x_fit, y_fitcross85, label=fitstrcross85, color='red', linestyle='--')

plt.title(r'$^{85}$Rb power')
plt.xlabel(x_label)
plt.ylabel(y_label)

handles, labels = plt.gca().get_legend_handles_labels()

plt.xticks(rotation=45)
plt.grid()
plt.tight_layout()
plt.savefig('data_linewidth/figures_unc/rb85_1.png')
plt.show()

fig_legend = plt.figure(figsize=(3, 2))
ax = fig_legend.add_subplot(111)
ax.axis('off')
legend = ax.legend(handles, labels, loc='center', fontsize=10, frameon=False)  

fig_legend.savefig('data_linewidth/figures_unc/legend_85.png', bbox_inches='tight', dpi=300)