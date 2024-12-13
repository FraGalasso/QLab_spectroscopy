import matplotlib.pyplot as plt
import functions as fn
import numpy as np
from matplotlib import cm
from fit_peaks import lorentzian_off


def plot_freq_laser_fit(piezo_fitted, volt_laser, file_name, A, x0, gamma, off, xpeaks, ypeaks, width, save=False):
    fitted_curves = []
    for A_, x0_, gamma_, peak, w, off_ in zip(A, x0, gamma, xpeaks, width, off):
        x = np.linspace(peak - w * 1.2, peak + w * 1.2, 100)
        y = lorentzian_off(x, A_, x0_, gamma_, off_)
        fitted_curves.append((x, y))

    cmap = cm.get_cmap('Oranges')
    colors = cmap(np.linspace(0.5, 0.9, len(fitted_curves)))

    plt.figure(figsize=(12, 6))
    plt.plot(piezo_fitted, volt_laser, label='Laser Intensity vs. Frequency',
             color='green', marker='.', linestyle=None)
    plt.scatter(xpeaks, ypeaks, marker='x', label='Peak Values')
    for i, (x, y) in enumerate(fitted_curves):
        plt.plot(x, y, '--', label=f'Fitted Lorentzian {i+1}', color=colors[i])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Laser Intensity (V)')
    plt.title('Frequency vs Laser Voltage')
    plt.legend()
    plt.grid()
    plt.xticks(rotation=45)
    plt.tight_layout()
    if save:
        plt.savefig(file_name)
        plt.close()
    else:
        plt.show()


def plot_piezo_laser_fit(time, volt_laser, file_name, A, x0, gamma, off, xpeaks, ypeaks, save=False):
    fitted_curves = []
    for A_, x0_, gamma_, off_, peak in zip(A, x0, gamma, off, xpeaks):
        x = np.linspace(peak - 3 * gamma_, peak +  3 * gamma_, 100)
        y = lorentzian_off(x, A_, x0_, gamma_, off_)
        fitted_curves.append((x, y))

    cmap = cm.get_cmap('Oranges')
    colors = cmap(np.linspace(0.5, 0.9, len(fitted_curves)))

    plt.figure(figsize=(12, 6))
    plt.scatter(time, volt_laser, label='Data',
             color='green', marker='.')
    plt.scatter(xpeaks, ypeaks, marker='x', label='Peak Values')
    for i, (x, y) in enumerate(fitted_curves):
        plt.plot(x, y, '--', label=f'Fitted Lorentzian {i+1}', color=colors[i])
    plt.xlabel('Voltage Piezo (V)')
    plt.ylabel('Laser Intensity (V)')
    plt.title('Piezo Voltage vs Laser Voltage')
    plt.legend()
    plt.grid()
    plt.xticks(rotation=45)
    plt.tight_layout()
    if save:
        plt.savefig(file_name)
        plt.close()
    else:
        plt.show()


def plot_calibrated_laser(xvalues_freq, volt_laser, file_name, extra_title='', save=False):
    plt.figure(figsize=(12, 6))
    plt.scatter(xvalues_freq, volt_laser, s=5,
                label='Laser Intensity vs. freq', color='green')
    plt.xlabel('Relative frequency values (GHz)')
    plt.ylabel('Laser Intensity (V)')
    plt.title(' Laser Intensity (calibrated)' + extra_title)
    plt.legend()
    plt.grid()
    plt.ticklabel_format(style='sci', axis='x',
                         scilimits=(9, 9), useOffset=False)
    plt.tight_layout()
    if save:
        plt.savefig(file_name)
        plt.close()
    else:
        plt.show()
