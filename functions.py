import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from fit_peaks import lorentzian_off


def plot_piezo_laser_fit(time, volt_laser, file_name, A, x0, gamma, off, xpeaks, ypeaks, save=False):
    fitted_curves = []
    for A_, x0_, gamma_, off_, peak in zip(A, x0, gamma, off, xpeaks):
        x = np.linspace(peak - 3 * gamma_, peak + 3 * gamma_, 100)
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
