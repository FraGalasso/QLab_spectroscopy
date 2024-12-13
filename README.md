code:

data_prep:  to be used for preparing data from the first oscilloscope. need to modify filename etc
data_prep_lecroy: to be used for preparing data from the second oscilloscope
linewidth.py: 


For LINEWIDTH ANALYSIS
dependencies: functions, data_prep_lecroy, linewidths, 
- data_prep_lecroy on all files
- Fra produced the peaks list files "intensity000x_peaks.csv"
- Marta analyzed with "linewidths_calib_plots.py"

For TEMPERATURE
dependencies:
- data_prep_lecroy on all files
- Fra rough fit and find peaks. Produces files Fixed_ld00000_peaks_fit.csv and Fixed_ld00000_peaks_fit_offset.csv
- Marta calibrates and removes 2.5 gamma around peaks using spectroscopy_temperature.py
- Fra refit the data with holes to extract temperature estimate 
