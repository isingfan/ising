import math
from math import *
import numpy as np
from numpy import *
from numpy.linalg import *
import random
import matplotlib.pyplot as plt
import pandas as pd
J=1 #interaction energy
lattice = "square"
L=60
N=L**2 #number of cells
print("number of sites is: " + str(N))
nmc1=70000000
nmc2=10000000 #number of monte carlo steps for each temperature value
nmc=nmc1+nmc2
#nmceq=40000 #number of states that are counted for expectation value
#Tmin=2.28
#Tmax=2.32
Tmin=1.0
Tmax=1.03
dT=0.01
#k=1.380649*(10**-23) #boltzmann constant in m^2 kg s^-2 K^-1
k=1.0

def getmagnetization(x):
    M=0
    for i in range (0, N):
        M=M+latt[i]
    return M

def randspin():
    a=random.random()
    if a < 0.5:
        return -1
    else:
        return 1

A=np.random.rand(N)
B=(A<0.5)
C=B.astype(int)

def getfirstneighbor(cell):
    i = (cell - 1)%N
    return i
def getsecondneighbor(cell):
    i = (cell + 1)%N
    return i
def getthirdneighbor(cell):
    i = (cell - L)%N
    return i
def getfourthneighbor(cell):
    i = (cell + L)%N
    return i

Tlist=[]
Mlist=[]
steplist=np.linspace(1, nmc2, num=nmc2)
indMlist=[None]*nmc2
Tarray=[]
t=Tmin
while t<=Tmax:
    Tarray.append(t)
    t=t+dT
print('Side length is:', L)
print('nmc is:', nmc)
#print('nmceq is:', nmceq)
for T in Tarray:
    print('Temperature is:', T)
    if T==Tarray[0]:
        latt=2*C-1
    else:
        latt=prevlatt
    beta=1.0/(k*T)
    acc1=math.exp(-4.0*beta*J)
    acc2=math.exp(-8.0*beta*J)
    M=0.0
    Madd=0.0
    nmccount=0
    for i in range (0, N):
        M=M+latt[i]
    if True:
        #for j in range (0, round(-abs(T-2.6)+3)*nmc):
        for j in range (0, nmc1):
            nmccount=nmccount+1
            celltoflip=random.randint(0, N-1)
            dE=2.0*J*latt[celltoflip]*(latt[getfourthneighbor(celltoflip)]+latt[getthirdneighbor(celltoflip)]+latt[getsecondneighbor(celltoflip)]+latt[getfirstneighbor(celltoflip)])
            if dE<=0.0:
                latt[celltoflip] = -latt[celltoflip]
                M=M+2.0*latt[celltoflip]
            else:
                rand=random.random()
                if dE==4.0*J:
                    if rand < acc1:
                        latt[celltoflip] = -latt[celltoflip]
                        M=M+2.0*latt[celltoflip]
                if dE==8.0*J:
                    if rand<acc2:
                        latt[celltoflip] = -latt[celltoflip]
                        M=M+2.0*latt[celltoflip]
        for j in range (0, nmc2):
            nmccount=nmccount+1
            celltoflip=random.randint(0, N-1)
            dE=2.0*J*latt[celltoflip]*(latt[getfourthneighbor(celltoflip)]+latt[getthirdneighbor(celltoflip)]+latt[getsecondneighbor(celltoflip)]+latt[getfirstneighbor(celltoflip)])
            if dE<=0.0:
                latt[celltoflip] = -latt[celltoflip]
                M=M+2.0*latt[celltoflip]
            else:
                rand=random.random()
                if dE==4.0*J:
                    if rand < acc1:
                        latt[celltoflip] = -latt[celltoflip]
                        M=M+2.0*latt[celltoflip]
                if dE==8.0*J:
                    if rand<acc2:
                        latt[celltoflip] = -latt[celltoflip]
                        M=M+2.0*latt[celltoflip]
            indMlist[j]=M
        plt.plot(steplist, indMlist)
        plt.xscale('log')
        plt.xlabel("Monte Carlo Step")
        plt.ylabel("Magnetization of state")
        plt.show()
        print('choose nmceqstart:')
        nmceqstart=int(input())
        #print('choose nmceqstart for the rest of the temperatures:')
        #nmceq=nmc-int(input())
        filesuffix="D:/WPy64-31090/notebooks/file" +"_Tmin_"+ str(Tmin) +"_Tmax_"+ str(Tmax) +"_dT_" + str(dT) + "_Side Length_" + str(L) + "_Sites number_" + str(N) + "_nmc_" + str(nmc) + "_nmceqstart_" + "_lattice type_" + lattice
        filedata = filesuffix+"_data.csv"
        for p in range(nmceqstart, nmc2):
            Madd=Madd+abs(M)
        Madd=Madd/(nmc2-nmceqstart)/N
        #if T==Tarray[1]:
        indMlist=[None]*nmc2
    Mlist.append(abs(Madd))
    print('expectation value of magnetization for temperature ', T, ' is ', abs(Madd))
    prevlatt=latt
    #indarr = np.array([steplist, indMlist])
    #indfileplot = filesuffix + "Temp" + str(T) + "indmag_plot.pdf"
    Tlist.append(T)
df = pd.DataFrame(Tlist).T
arr = np.array([Tlist, Mlist])
np.savetxt(filedata, arr.T, delimiter=",")
#loading:
#arr=np.loadtxt(filedata, delimiter=",")
fileplot=filesuffix+"_plot.pdf"
plt.xlabel("Temperature (Kelvin)")
plt.ylabel("Expected magnetization per site")
plt.plot(Tlist, Mlist, 'o')
plt.savefig(fileplot)
print('Simulation complete.')
#plt.show()
