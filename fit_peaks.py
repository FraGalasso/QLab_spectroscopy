from scipy.signal import find_peaks, peak_widths
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import cm


def lorentzian_off(x, A, x0, gamma, offset):
    return A / (1 + ((x - x0) / gamma) ** 2) + offset

def fit_peaks_spectroscopy(x, y, height, distance):
    '''Given x and y (and peaks height and distance) returns a list of lorentzian curves
    (just the 3 parameters A, x0, gamma) and a list with the covariance matrices from the fit.'''

    x_spacing = np.mean(np.diff(x))

    # Find peaks and widths
    peaks, _ = find_peaks(y, height=height, distance=distance)
    widths_full = peak_widths(y, peaks, rel_height=1)[0]

    # Loop through each peak and fit
    params = []
    covs = []
    for peak, width in zip(peaks, widths_full):
        # Determine a fitting range around the peak, i.e. width/2
        # fit_range = min(width/2, 0.001/x_spacing) # use this for data9
        fit_range = min(width/2, 0.07/x_spacing) # use this for data10
        # fit_range = width/2
        start = max(0, peak - int(fit_range))
        end = min(len(x), peak + int(fit_range))

        # Extract data around the peak
        x_fit_range = x[start:end]
        y_fit_range = y[start:end]
        width_scaled = 2 * fit_range * x_spacing

        # Initial guess: A=height at peak, x0=peak position in x_fitted, gamma=half-width at half-maximum
        # initial_guess = [y[peak], x[peak], 0.0001, -0.01] # use this for data9
        initial_guess = [y[peak], x[peak], 0.02, 0] # use this for data10
        # initial_guess = [y[peak], x[peak], 0.001, 0]

        # Define bounds for A, x0, and gamma
        bounds = (
            # Lower bounds for [A, x0, gamma, off]
            [0, x[peak] - width_scaled, 0, -np.inf],
            # Upper bounds for [A, x0, gamma, off]
            [np.inf, x[peak] + width_scaled, width_scaled, np.inf]
        )

        try:
            popt, pcov = curve_fit(lorentzian_off, x_fit_range, y_fit_range,
                                   p0=initial_guess, bounds=bounds, maxfev=10000)
            params.append(popt)
            covs.append(pcov)
        except RuntimeError as e:
            print(
                f"Failed to fit peak at voltage = {x[peak]:.4f} due to RuntimeError: {e}")
        except Exception as e:
            print(
                f"An unexpected error occurred while fitting peak at voltage = {x[peak]:.4f}: {e}")

    return params, covs

def fit_peaks_leonardi(x, y, peaks, widths, ind_range):
    '''Given x and y, a list of peaks, FWHMs and the number of points to use in the fit, returns a list of lorentzian curves
    (just the 3 parameters A, x0, gamma) and a list with the covariance matrices from the fit.'''

    indices = [np.flatnonzero(x == pk)[0] for pk in peaks if pk in x]

    # Loop through each peak and fit
    params = []
    covs = []
    for peak_index, width in zip(indices, widths):
        # Determine a fitting range around the peak, e.g., ±1.2 * width
        start = max(0, peak_index - int(ind_range/2))
        end = min(len(x), peak_index + int(ind_range/2))

        # Extract data around the peak
        x_fit_range = x[start:end]
        y_fit_range = y[start:end]

        # Initial guess: A=height at peak, x0=peak position in x_fitted, gamma=half-width at half-maximum
        initial_guess = [y[peak_index], x[peak_index], width / 2, 0]

        # Define bounds for A, x0, gamma and off
        bounds = (
            # Lower bounds for [A, x0, gamma, off]
            [0, x[peak_index] - 2.5 * width, 0, -np.inf],
            # Upper bounds for [A, x0, gamma, off]
            [np.inf, x[peak_index] + 2.5 * width, width * 5, np.inf]
        )

        try:
            popt, pcov = curve_fit(lorentzian_off, x_fit_range, y_fit_range,
                                   bounds=bounds, maxfev=10000)
            params.append(popt)
            covs.append(pcov)
        except RuntimeError as e:
            print(
                f"Failed to fit peak at piezo_fitted = {x[peak_index]:.2f} due to RuntimeError: {e}")
        except Exception as e:
            print(
                f"An unexpected error occurred while fitting peak at piezo_fitted = {x[peak_index]:.2f}: {e}")

    return params, covs


def plot_piezo_laser_fit_leonardi(piezo_fitted, volt_laser, file_name, A, x0, gamma, off, xpeaks, ypeaks, ind_range, save=False, log=False):
    fitted_curves = []
    piezo_spacing = np.mean(np.diff(piezo_fitted))
    for A_, x0_, gamma_, off_, peak in zip(A, x0, gamma, off, xpeaks):
        x = np.linspace(peak - ind_range * piezo_spacing/2,
                        peak + ind_range * piezo_spacing/2, 100)
        y = lorentzian_off(x, A_, x0_, gamma_, off_)
        fitted_curves.append((x, y))

    cmap = cm.get_cmap('Oranges')
    colors = cmap(np.linspace(0.5, 0.9, len(fitted_curves)))

    plt.figure(figsize=(12, 6))
    plt.rcParams.update({'font.size': 16})
    plt.scatter(piezo_fitted, volt_laser, label='Laser Intensity vs. Piezo volt',
                color='green', marker='.')
    plt.scatter(xpeaks, ypeaks, marker='x', label='Peak Values')
    for i, (x, y) in enumerate(fitted_curves):
        plt.plot(x, y, '--', label=f'Fitted Lorentzian {i+1}', color=colors[i])
    plt.xlabel('Voltage Piezo (V)')
    plt.ylabel('Residuals (V)')
    if log:
        plt.yscale('log')
    # plt.legend()
    plt.grid()
    plt.xticks(rotation=45)
    plt.tight_layout()
    if save:
        plt.savefig(file_name)
        plt.close()
    else:
        plt.show()