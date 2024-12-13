from scipy.signal import find_peaks, peak_widths
import numpy as np
from scipy.optimize import curve_fit


def lorentzian(x, A, x0, gamma):
    return A / (1 + ((x - x0) / gamma) ** 2)


def lorentzian_off(x, A, x0, gamma, offset):
    return A / (1 + ((x - x0) / gamma) ** 2) + offset


def fit_peaks(x, y, height, distance):
    '''Given x and y (and peaks height and distance) returns a list of lorentzian curves
    (just the 3 parameters A, x0, gamma) and a list with the covariance matrices from the fit.'''

    x_spacing = np.mean(np.diff(x))

    # Find peaks and widths
    peaks, _ = find_peaks(y, height=height, distance=distance)
    widths_half = peak_widths(y, peaks, rel_height=0.5)[0]

    # Loop through each peak and fit
    params = []
    covs = []
    for peak, width in zip(peaks, widths_half):
        # Determine a fitting range around the peak, e.g., ±1.2 * width
        fit_range = int(width * 1.2)
        start = max(0, peak - fit_range)
        end = min(len(x), peak + fit_range)

        # Extract data around the peak
        x_fit_range = x[start:end]
        y_fit_range = y[start:end]
        width_scaled = width * x_spacing

        # Initial guess: A=height at peak, x0=peak position in x_fitted, gamma=half-width at half-maximum
        initial_guess = [y[peak], x[peak], width_scaled / 2]

        # Define bounds for A, x0, and gamma
        bounds = (
            # Lower bounds for [A, x0, gamma]
            [0, x[peak] - width_scaled, 0],
            # Upper bounds for [A, x0, gamma]
            [np.inf, x[peak] + width_scaled, width_scaled * 2]
        )

        try:
            popt, pcov = curve_fit(lorentzian, x_fit_range, y_fit_range,
                                   p0=initial_guess, bounds=bounds, maxfev=10000)
            params.append(popt)
            covs.append(pcov)
        except RuntimeError as e:
            print(
                f"Failed to fit peak at piezo_fitted = {x[peak]:.2f} due to RuntimeError: {e}")
        except Exception as e:
            print(
                f"An unexpected error occurred while fitting peak at piezo_fitted = {x[peak]:.2f}: {e}")

    return params, covs


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
