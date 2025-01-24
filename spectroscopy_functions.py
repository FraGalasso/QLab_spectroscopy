import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

def lin_quad_fits(x, y):
    # Perform linear fit
    coeffs_1 = np.polyfit(x, y, 1)  # coeffs = [a, b] for y = ax + b
    # Create a polynomial function from the coefficients
    linear_fit = np.poly1d(coeffs_1)
    # Perform quadratic fit
    coeffs_2 = np.polyfit(x, y, 2)  # coeffs = [a, b, c] for y = ax^2 + bx + c
    # Create a polynomial function from the coefficients
    quadratic_fit = np.poly1d(coeffs_2)

    # Generate x values for the fitted curve (same x range as the original data)
    x_fit = np.linspace(min(x), max(x), 100)  # Smooth line for plotting
    lin_fit = linear_fit(x_fit)  # Calculate the fitted line
    # Calculate corresponding y values for the fitted curve
    quad_fit = quadratic_fit(x_fit)
    return coeffs_1, coeffs_2, x_fit, lin_fit, quad_fit

def plot_fits(x, y, x_label, y_label, title, file_name, save):
    coeffs_1, coeffs_2, x_fit, lin_fit, quad_fit = lin_quad_fits(x, y)

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
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
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
        
def plotting(x, y, x_label, y_label, title, file_name, save):
    plt.figure(figsize=(12, 6))
    plt.plot(x, y, label='Data', color='green')
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    plt.xticks(rotation=45)
    plt.grid()
    plt.tight_layout()
    if save:
        plt.savefig(file_name)
        plt.close()
    else:
        plt.show()

def scattering(x, y, x_label, y_label, title, file_name, save, fit):
    """
    Scatter plot function with optional linear fit.

    Parameters:
    - x: array-like, x-axis data
    - y: array-like, y-axis data
    - x_label: str, label for the x-axis
    - y_label: str, label for the y-axis
    - title: str, title of the plot
    - file_name: str, file name for saving the plot
    - save: bool, whether to save the plot or display it
    - fit: bool, whether to compute and plot a linear fit
    """
    plt.figure(figsize=(12, 6))
    plt.scatter(x, y, label='Data', color='green')

    # Optional linear fit
    if fit:
        # Perform a linear fit
        coefficients = np.polyfit(x, y, 1)  # Fit to a line y = mx + c
        linear_fit = np.poly1d(coefficients)  # Create a polynomial object
        y_fit = linear_fit(x)  # Evaluate the fit at x points

        # Plot the linear fit
        plt.plot(x, y_fit, color='red', label=f'Fit: y = {coefficients[0]:.2f}x + {coefficients[1]:.2f}')

    # Plot details
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    plt.xticks(rotation=45)
    plt.grid()
    plt.tight_layout()

    # Save or show the plot
    if save:
        plt.savefig(file_name)
        plt.close()
    else:
        plt.show()


def plotting3(x, y1, y2, y3, x_label, y1_label, y2_label, y3_label, title, file_name, save):
    plt.figure(figsize=(12, 6))
    plt.plot(x, y1, label=y1_label, color='green')
    plt.plot(x, y2, label=y2_label, color='blue')
    plt.plot(x, y3, label=y3_label, color='red')
    plt.title(title)
    plt.xlabel(x_label)
    plt.legend()
    plt.xticks(rotation=45)
    plt.grid()
    plt.tight_layout()
    if save:
        plt.savefig(file_name)
        plt.close()
    else:
        plt.show()

def remove_peaks2(peaks_mean, peaks_gamma, data, gamma_factor):
    '''
    Removes regions around peaks in peaks_mean from "data",
    where the region is defined by gamma_factor * gamma around each peak.
    '''
    new_data = data.copy()

    for mean, gamma in zip(peaks_mean, peaks_gamma):
        begin_volt = mean - gamma_factor * gamma
        end_volt = mean + gamma_factor * gamma
        # Apply the filtering
        new_data = new_data[(new_data['volt_piezo'] < begin_volt) | (new_data['volt_piezo'] > end_volt)]

    return new_data
