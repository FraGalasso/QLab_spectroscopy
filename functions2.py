import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import functions as fn


def peaks_fsr(piezo_voltage, laser_voltage):

    # Use find_peaks with specified height and distance
    peaks_indices, _ = find_peaks(laser_voltage, height=0.2, distance=400)

    # Convert peaks_indices to an integer NumPy array
    peaks_indices = np.array(peaks_indices, dtype=int)
    # Extract peak values
    peaks = laser_voltage[peaks_indices]
    peaks_xvalues = piezo_voltage[peaks_indices]

    return peaks_xvalues, peaks, peaks_indices


def peaks_hfsr(piezo_voltage, laser_voltage):

    # Use find_peaks with specified height and distance
    peaks_indices, _ = find_peaks(laser_voltage, height=0.08, distance=400)

    # Convert peaks_indices to an integer NumPy array
    peaks_indices = np.array(peaks_indices, dtype=int)
    # Extract peak values
    peaks = laser_voltage[peaks_indices]
    peaks_xvalues = piezo_voltage[peaks_indices]

    return peaks_xvalues, peaks, peaks_indices


def plot_fits(x, y, x_label, y_label, file_name, confocal, save=False):
    coeffs_1, coeffs_2, x_fit, lin_fit, quad_fit = fn.lin_quad_fits(x, y)

    # Format the coefficients for display
    lin_coeff_label = f"Linear Fit: y = ({coeffs_1[0]:.2e})x + ({coeffs_1[1]:.2e})"
    quad_coeff_label = (
        f"Quadratic Fit: y = ({coeffs_2[0]:.2e})x² + ({coeffs_2[1]:.2e})x + ({coeffs_2[2]:.2e})"
    )

    plt.figure(figsize=(12, 6))
    plt.scatter(x, y, label='Data', color='green', marker='x', s=30)
    plt.plot(x_fit, lin_fit, label=lin_coeff_label,
             color='blue', linestyle='--')  # Plot the linear fit
    plt.plot(x_fit, quad_fit, label=quad_coeff_label, color='red',
             linestyle='--')  # Plot the quadratic fit
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if (confocal):
        plt.title('Assuming confocality')
    else:
        plt.title('Without assuming confocality')
    plt.legend()
    plt.xticks(rotation=45)
    plt.grid()
    plt.tight_layout()
    if save:
        plt.savefig(file_name)
        plt.close()
    else:
        plt.show()
    return coeffs_1, coeffs_2
