import characterization_data as data
import matplotlib.pyplot as plt
import os
import piecewise_regression
from uncertainties import ufloat


def char_plot_fit(xx, yy, filename, label, title, figurename, save):
    pw_fit = piecewise_regression.Fit(xx, yy, n_breakpoints=1)
    pw_results = pw_fit.get_results()
    threshold = ufloat(pw_results["estimates"]['breakpoint1']
                       ['estimate'], pw_results["estimates"]['breakpoint1']['se'])

    with open(filename, 'a') as file:
        file.write(pw_fit.summary())

    plt.figure()
    pw_fit.plot_data(color="blue", s=20, label=label)
    fitlabel = r'Fit, $\text{I}_{\text{thr}} = $' + f'{threshold} mA'
    pw_fit.plot_fit(color="red", linewidth=2, label=fitlabel)
    pw_fit.plot_breakpoints()
    pw_fit.plot_breakpoint_confidence_intervals()
    plt.title(title)
    plt.xlabel(r'$\text{I}_{\text{diode}}$ mA')
    plt.ylabel(r'P [$\mu$ W]')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    if save:
        plt.savefig(figurename)
        plt.close()
    else:
        plt.show()
    return threshold


if __name__ == '__main__':
    folder_name = 'laser_characterization'
    os.makedirs(folder_name, exist_ok=True)

    '''checking threshold current'''
    aligned_file = os.path.join(folder_name, 'aligned_char.txt')

    with open(aligned_file, 'w') as file:
        file.write('Cavity aligned on 24/09:')
    thr_aligned_1 = char_plot_fit(data.x_aligned, data.y_aligned, filename=aligned_file,
                                  label='Aligned cavity 1', title='', figurename=f'{folder_name}/aligned_1.png', save=True)

    with open(aligned_file, 'a') as file:
        file.write('Cavity aligned on 01/10:')
    thr_aligned_2 = char_plot_fit(data.x, data.y, filename=aligned_file,
                                  label='Aligned cavity 2', title='', figurename=f'{folder_name}/aligned_2.png', save=True)

    with open(aligned_file, 'a') as file:
        file.write(f'Final estimate for laser threshold:{
                   min(thr_aligned_1, thr_aligned_2)} mA')

    temp_file = os.path.join(folder_name, 'temp.txt')
    with open(temp_file, 'w') as file:
        file.write('Unaligned cavity 23°C: only metal screw')
    thr_23_1 = char_plot_fit(data.x_nonal, data.y_nonal, filename=temp_file,
                             label='Unaligned cavity: 23°C', title='', figurename=f'{folder_name}/unaligned_23_1.png', save=True)

    with open(temp_file, 'a') as file:
        file.write('Unaligned cavity 23°C: metal and plastic screws')
    thr_23_2 = char_plot_fit(data.x_nonal2, data.y_nonal2, filename=temp_file,
                             label='Unaligned cavity: 23°C', title='', figurename=f'{folder_name}/unaligned_23_2.png', save=True)

    with open(temp_file, 'a') as file:
        file.write('Unaligned cavity 40°C')
    thr_40 = char_plot_fit(data.x2, data.y2, filename=temp_file,
                           label='Unaligned cavity: 40°C', title='', figurename=f'{folder_name}/unaligned_40.png', save=True)

    with open(temp_file, 'a') as file:
        file.write('Unaligned cavity 15°C')
    thr_15 = char_plot_fit(data.x3, data.y3, filename=temp_file,
                           label='Unaligned cavity: 15°C', title='', figurename=f'{folder_name}/unaligned_15.png', save=True)
