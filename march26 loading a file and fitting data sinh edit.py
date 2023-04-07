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
with open ("D:/WPy64-31090/notebooks/Data pairings (square lattice)/Reduced csv files/square/file_Tmin_1.6_Tmax_3.0_dT_0.05_Dimension_90_nmc_80000000_nmceqstart_15000000_nmceq_50000000_data_reduced_crit 2,3032.csv", "r") as i:
    rawdata = list(csv.reader(i, delimiter = ","))

exampledata = np.array(rawdata[0:], dtype=float)
xdata = exampledata[:,0]
ydata = exampledata[:,1]
#plot the data
plt.figure(1, dpi=120)
#plt.title("")
plt.xlabel("Temperature (Kelvin)")
plt.ylabel("Expected magnetization per site")
plt.xlim(1.5, 2.5)
plt.ylim(0, 1.2)
plt.xscale("linear")
plt.yscale("linear")
plt.plot(xdata, ydata,  'o', label = "Experimental Data")

#define constants
cset = 2.03
gset = 0.15
Tdisplaymin=1.6
Tdisplaymax=2.5
#define function
def func(x, c, g):
    return (1.0-(1.0/(sinh(c/x))**4))**g

#write array for display
xdisplay=np.array([])
for i in range(int(Tdisplaymin*100), int(Tdisplaymax*100)):
    xdisplay = np.append(xdisplay, i*0.01)

#evaluate and plot function
funcdata = func(xdata, cset, gset)
funcdisplay=np.array([])
for j in np.arange (0, int(Tdisplaymax*100-Tdisplaymin*100)):
    funcdisplay = np.append(funcdisplay, func(xdisplay[j], cset, gset))
plt.plot(xdisplay, funcdisplay, label="Model")
plt.legend()

#curve fit data to model
popt, pcov = curve_fit(func, xdata, ydata, bounds = (1.6, 2.5))
perr = np.sqrt(np.diag(pcov))
print(popt)
#print(pcov)
plt.show()
