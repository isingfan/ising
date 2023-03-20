import math
from math import *
import numpy as np
from numpy import *
from numpy.linalg import *
import random
import matplotlib.pyplot as plt
import pandas as pd
import csv
from scipy.optimize import curve_fit
#import data
with open ("D:/WPy64-31090/notebooks/Data pairings (square lattice)/Reduced csv files/file_Tmin_1.8_Tmax_3.0_dT_0.025_Dimension_55_nmc_50000000_nmceqstart_938000_nmceq_45000000_data reduced.csv", "r") as i:
    rawdata = list(csv.reader(i, delimiter = ","))

exampledata = np.array(rawdata[0:], dtype=float)
xdata = exampledata[:,0]
ydata = exampledata[:,1]
#plot the data
plt.figure(1, dpi=120)
#plt.title("")
plt.xlabel("Temperature (Kelvin)")
plt.ylabel("Magnetization")
plt.xlim(1.7, 3.0)
plt.ylim(0, 1.2)
plt.xscale("linear")
plt.yscale("linear")
plt.plot(xdata, ydata,  'o', label = "Experimental Data")

#define constants
kset = 1.0967
Cset = 2.325
bset = 0.172
Tmin=1.7
Tmax=2.7
Tdisplaymin=1.0
Tdisplaymax=2.7
#define function
def func(x, k, C, b):
    return k*(C-x)**b

#write array for display
xdisplay=np.array([])
for i in range(int(Tdisplaymin*100), int(Tdisplaymax*100)):
    xdisplay = np.append(xdisplay, i*0.01)
#print("xdisplay is: " + str(xdisplay))
#evaluate and plot function
funcdata = func(xdata, kset, Cset, bset)
funcdisplay=np.array([])
for j in np.arange (0, int(Tdisplaymax*100-Tdisplaymin*100)):
    funcdisplay = np.append(funcdisplay, func(xdisplay[j], kset, Cset, bset))
print(funcdisplay)
plt.plot(xdisplay, funcdisplay, label="Model")
plt.legend()

#curve fit data to model
popt, pcov = curve_fit(func, xdata, ydata, bounds = (0, 10))
perr = np.sqrt(np.diag(pcov))
#print(popt)
#print(pcov)
plt.show()
