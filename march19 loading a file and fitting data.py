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
with open ("D:/WPy64-31090/notebooks/file_Tmin_2.3_Tmax_3.3_dT_0.05_Dimension_30_nmc_2000000_nmceqstart_153000_nmceq_600000_lattice type_rhombus_data.csv", "r") as i:
    rawdata = list(csv.reader(i, delimiter = ","))

exampledata = np.array(rawdata[0:], dtype=float)
xdata = exampledata[:,0]
ydata = exampledata[:,1]
#plot the data
plt.figure(1, dpi=120)
#plt.title("")
plt.xlabel("Temperature (Kelvin)")
plt.ylabel("Magnetization")
plt.xlim(2.3, 3.3)
plt.ylim(0, 1.2)
plt.xscale("linear")
plt.yscale("linear")
plt.plot(xdata, ydata,  'o', label = "Experimental Data")

#define constants
kset = 1.0967
Cset = 2.325
bset = 0.172
Tdisplaymin=2.3
Tdisplaymax=3.3
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
print(popt)
#print(pcov)
plt.show()
