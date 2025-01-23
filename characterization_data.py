# 780 nm diode laser data
#x data: current through the diode (actually a voltage, but resistance is 1 ohm so it is the same as the current)
#y data: laser intensity measured with a powermeter. 

#list of datasets present on this file:
#2 curves with aligned cavity at 23 deg (24 sept and 1 oct)
#2 curves with the unaligned cavity (24 sept)
#2 curves at different temperatures with the cavity being unaligned : 40 deg and 15 deg


#NB the y is actually a relative measurement, as the powermeter was not calibrated and it depends on factors such as distance from diode and background. However we tried to keep the powermeter at the same distance from the diode for all data sets, and selected the wavelength on the powermeter and noticed that this makes it insensitive to background light.

#sept 24:
#temperature: 23 degrees
#diffraction grating is aligned
x_aligned = [
    0, 5, 10, 15.1, 20.0, 25.0, 30.1, 35.0, 40.0, 45.1, 46.1, 46.5, 47.0, 47.3,
    47.7, 48.0, 48.5, 49.0, 49.5, 50, 51, 52, 54, 56.1, 58.1, 60.1
]
y_aligned = [
    0.27, 0.45, 1.68, 4.73, 9.08, 14.58, 21.4, 30.1, 41.9, 58.2, 66.5, 72.4,
    84.2, 106.4, 183.2, 230, 308, 803, 920, 1035, 1219, 1414, 2070, 2480, 3040,
    3370
]

#now the metal screw is unscrewed
x_nonal = [
    0, 5, 10, 15.1, 20.0, 25.0, 30.0, 35.1, 40.1, 45.0, 46.1, 47.1, 47.5, 47.8,
    48.0, 49.0, 49.5, 50.2, 51.1, 52.1
]
y_nonal = [
    0.23, 0.38, 1.65, 4.55, 8.88, 14.43, 21.1, 27.7, 40.8, 63.4, 71.2, 94.1,
    112.5, 137.8, 202, 640, 809, 943, 1150, 1354
]
#now the plastic screw is also unscrewed
x_nonal2 = [
    0, 5.1, 9.9, 15, 20.1, 25.0, 30.2, 35.1, 39.9, 45, 46, 46.5, 47.0, 47.5,
    48.0, 48.3, 48.8, 49.1, 49.5, 50.1, 52, 54, 56, 58, 60
]
y_nonal2 = [
    0.42, 0.5, 1.53, 4.30, 8.56, 13.90, 20.7, 29.0, 39.4, 55.6, 64.2, 69.4,
    77.7, 91.7, 146.2, 210, 349, 429, 560, 760, 1290, 1860, 2360, 2850, 3250
]

#oct 1
#we realigned the diffraction grating and took another curve to see if it matches the first one. still 23 degrees.
x = [
    0, 10.0, 20.0, 30.0, 35.1, 40.0, 43.0, 45.0, 45.5, 46.0, 46.6, 47, 47.5,
    48.0, 48.5, 49.0, 49.5, 50, 51, 52, 55.1, 57.9, 60
]
y = [
    0.09, 1.27, 8.44, 21.2, 30.0, 41.3, 50.9, 60.6, 62.3, 66.1, 72.0, 79.5,
    106.0, 205.5, 300, 392, 877, 1008, 1224, 1410, 2260, 2790, 3300
]

#now, change the temperature and take more curves. NB: if the temperature is changed, the alignment will also be lost. so these next curves should be compared to the 23 degrees unaligned curve

#40 degrees:
x2 = [
    0, 10.0, 20.0, 30.1, 35.1, 40, 43.0, 44.9, 45.6, 46, 46.5, 47.0, 47.6,
    48.0, 48.5, 49, 49.5, 50.0, 50.5, 51, 51.5, 52, 52.5, 54, 56.9, 60.1
]
y2 = [
    0.06, 1.64, 9.66, 23.6, 32.2, 43.3, 52.3, 59.7, 62.0, 63.8, 67.2, 71.6,
    76.8, 81.5, 89.7, 103.7, 135, 205, 346, 513, 647, 791, 920, 1363, 2160,
    2980
]
#15 degrees:
x3 = [
    0, 10, 20.1, 30.1, 35.0, 40.1, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.1,
    52.1, 54, 56.1, 58, 60.1
]
y3 = [
    0.04, 1.34, 8.45, 21.0, 29.5, 40.5, 59.0, 66.5, 82.0, 167.4, 392, 702, 980,
    1252, 1793, 2270, 2750, 3160
]
