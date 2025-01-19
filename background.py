import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import fit_peaks as fp
import functions as fn

folder = 'data_temp'
titles = ['background00000', 'no_absorption00000']

means = []
std_devs = []

for title in titles:
    data_file = f'{folder}/clean_data/{title}.csv'

    data = pd.read_csv(data_file)

    volt_laser = data['volt_laser'].to_numpy()
    volt_piezo = data['volt_piezo'].to_numpy()
    
    means.append(np.mean(volt_laser))
    std_devs.append(np.std(volt_laser))

    plt.figure()
    plt.plot(volt_piezo, volt_laser,
            color='blue', markersize=5, marker='.')
    plt.xlabel('Volt Piezo [V]')
    plt.ylabel('Volt Laser [V]')
    plt.title(title)
    plt.grid()
    plt.tight_layout()
    plt.show()

with open('data_temp/background.txt', 'w') as file:
    file.write(f'Background: {means[0]} +/- {std_devs[0]} V')
    file.write(f'\nNo absorption reading: {means[1]} +/- {std_devs[1]} V')