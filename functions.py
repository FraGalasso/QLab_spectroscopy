import numpy as np
from scipy.signal import find_peaks, peak_widths
import matplotlib.pyplot as plt
import pandas as pd

# DATA PREP


def crop_to_min_max(time, laser_voltage, piezo_voltage, ld_volt=None):
    '''Crops time and voltage readings in order to include just one
     frequency sweep from the piezo.

     If necessary crops diode modulation data. Crop is perfomed only
     following piezo data, so if you have to crop data in order to include
     just one modulation sweep, you need to invert piezo_voltage and ld_volt.'''

    if len(piezo_voltage) == 0:
        return None  # Return None if piezo_voltage is empty

    # Find indices of minimum and maximum values
    # Assuming there will be just one full sweep in the dataset
    min_index = np.argmin(piezo_voltage)
    max_index = np.argmax(piezo_voltage)

    # Ensure the possibility to analyse "backwards" sweeps too
    start_index = min(min_index, max_index)
    end_index = max(min_index, max_index)

    # Extract the data between min_index and max_index
    time_window = time[start_index:end_index + 1]
    laser_window = laser_voltage[start_index:end_index + 1]
    piezo_window = piezo_voltage[start_index:end_index + 1]

    if ld_volt is not None:
        ld_window = ld_volt[start_index:end_index + 1]
        return time_window, laser_window, piezo_window, ld_window
    else:
        return time_window, laser_window, piezo_window


def crop_to_min_max_extended(time, laser_voltage, piezo_voltage, ld_volt):
    '''Crops time and voltage readings in order to include just one
     frequency sweep from the piezo and from modulation of laser current.'''

    min_index_piezo = np.argmin(piezo_voltage)
    max_index_piezo = np.argmax(piezo_voltage)
    min_index_diode = np.argmin(ld_volt)
    max_index_diode = np.argmax(ld_volt)

    start_index_piezo = min(min_index_piezo, max_index_piezo)
    end_index_piezo = max(min_index_piezo, max_index_piezo)
    start_index_diode = min(min_index_diode, max_index_diode)
    end_index_diode = max(min_index_diode, max_index_diode)

    start_index = max(start_index_diode, start_index_piezo)
    end_index = min(end_index_diode, end_index_piezo)

    time_window = time[start_index:end_index + 1]
    laser_window = laser_voltage[start_index:end_index + 1]
    piezo_window = piezo_voltage[start_index:end_index + 1]
    ld_window = ld_volt[start_index:end_index + 1]

    return time_window, laser_window, piezo_window, ld_window


def fit_piezo_line(time, piezo_voltage):
    '''Converts timestamps in voltages on piezo.

    Returns voltages from a linear interpolation of input data.'''

    if len(time) == 0 or len(piezo_voltage) == 0 or len(time) != len(piezo_voltage):
        return None  # Return None if the input arrays are empty or of different lengths

    # Fit a line (degree 1 polynomial) to the piezo voltage data
    slope, intercept = np.polyfit(time, piezo_voltage, 1)
    piezo_fit = slope * time + intercept

    print('Fitting a*x+b:')
    print(f'slope = {slope} V/s\t intercept = {intercept} V')

    return piezo_fit
# PEAKS


def peaks(piezo_voltage, laser_voltage, height, distance):
    '''
    Finds peaks in readings from the photodiode.

    Parameters:

    piezo_voltage: voltages on the piezo, should be the cleaned version (the output of fit_piezo_line)
    laser_voltage: voltages on the photodiode
    height: min height of the peaks
    height: min number of point between a peak and the following one

    Returns: 

    peaks_xvalues: voltage values on the piezo corresponding to detected peaks
    peaks: detected peaks
    scaled_widths: width of each peak in Volts
    '''
    peaks_indices, _ = find_peaks(
        laser_voltage, height=height, distance=distance)
    peaks = laser_voltage[peaks_indices]
    peaks_xvalues = piezo_voltage[peaks_indices]

    widths = peak_widths(laser_voltage, peaks_indices, rel_height=0.5)
    piezo_voltage_spacing = np.mean(np.diff(piezo_voltage))
    scaled_widths = widths[0]*piezo_voltage_spacing

    return peaks_xvalues, peaks, scaled_widths

# CALIBRATION

def calibrate(x, coeffs1):
    new_x = coeffs1[0] * x + coeffs1[1]
    return new_x
    

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

# REMOVING PEAKS

def remove_singlepeak_data(data, beginning_time, end_time):
    '''
    Removes data and corresponding frequencies within the specified time window [beginning_time, end_time].
    '''
    # Create a mask for values outside the [beginning_time, end_time] interval
    mask = (data['timestamp'] < beginning_time) | (
        data['timestamp'] > end_time)

    # Apply the mask to filter the data and frequencies
    cropped_data = data[mask]

    return cropped_data


def remove_peaks(peaks_mean, peaks_gamma, data, gamma_factor):
    '''
    Removes regions around peaks in peaks_mean from "data" and "frequencies",
    where the region is defined by gamma_factor * gamma around each peak.
    '''

    # Make copies of the original data and frequencies
    new_data = data.copy()

    # Loop through peaks and remove regions
    for mean, gamma in zip(peaks_mean, peaks_gamma):
        begin_time = mean - gamma_factor * gamma
        end_time = mean + gamma_factor * gamma
        new_data = remove_singlepeak_data(new_data, begin_time, end_time)

    return new_data


# LOGISTICS: SAVING AND PLOTTING
def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def save_to_csv(timestamps, volt_laser, volt_piezo, piezo_fitted, output_file):
    """
    Save the cropped and fitted data to a CSV file with the specified columns.

    Parameters:
    - timestamps: Array of timestamps
    - volt_laser: Array of laser voltages
    - volt_piezo: Array of piezo voltages
    - piezo_fitted: Array of fitted piezo voltage values
    - output_file: Path to the output CSV file
    """

    # Create a DataFrame with the data
    data = {
        'timestamp': timestamps,
        'volt_laser': volt_laser,
        'volt_piezo': volt_piezo,
        'piezo_fitted': piezo_fitted
    }

    df = pd.DataFrame(data)

    # Save the DataFrame to a CSV file
    df.to_csv(output_file, index=False)
    print(f"Data saved to {output_file}")


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


def scattering(x, y, x_label, y_label, title, file_name, save):
    plt.figure(figsize=(12, 6))
    plt.scatter(x, y, label='Data', color='green')
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


def scatter3(x, y1, y2, y3, x_label, y1_label, y2_label, y3_label, title, file_name, save):
    plt.figure(figsize=(12, 6))
    plt.scatter(x, y1, label=y1_label, color='green')
    plt.scatter(x, y2, label=y2_label, color='blue')
    plt.scatter(x, y3, label=y3_label, color='red')
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
