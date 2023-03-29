import math
from math import *
import numpy as np
from numpy import *
from numpy.linalg import *
import random
import matplotlib.pyplot as plt
import pandas as pd
J=1 #interaction energy
lattice = "hexagon"
L=40
N=2*L-1 #number of cells
for x in range (1, int(L)):
    N=N + 2*(2*L-1-x)
print(N)
nmc=8000000 #number of monte carlo steps for each temperature value
#nmceq=40000 #number of states that are counted for expectation value
Tmin=1.2
Tmax=3.0
dT=0.05
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
    i = (cell - (3*L-1))%N
    return i
def getfourthneighbor(cell):
    i = (cell + (3*L-1))%N
    return i
def getfifthneighbor(cell):
    i = (cell + (3*L-2))%N
    return i
def getsixthneighbor(cell):
    i = (cell + (3*L-2))%N
    return i

Tlist=[]
Mlist=[]
steplist=[]
indMlist=[]
Tarray=[]
t=Tmin
while t<=Tmax:
    Tarray.append(t)
    t=t+dT
nmctotal=0
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
    acc3=math.exp(-12.0*beta*J)
    M=0.0
    Madd=0.0
    nmccount=0
    for i in range (0, N):
        M=M+latt[i]
    if T==Tarray[0]:
        #for j in range (0, round(-abs(T-2.6)+3)*nmc):
        for j in range (0, nmc):
            nmccount=nmccount+1
            celltoflip=random.randint(0, N-1)
            oldspin=latt[celltoflip]
            dE=2.0*J*latt[celltoflip]*(latt[getsixthneighbor(celltoflip)]+latt[getfifthneighbor(celltoflip)]+latt[getfourthneighbor(celltoflip)]+latt[getthirdneighbor(celltoflip)]+latt[getsecondneighbor(celltoflip)]+latt[getfirstneighbor(celltoflip)])
            if dE<=0.0:
                latt[celltoflip] = -latt[celltoflip]
            else:
                if dE==4.0*J:
                    rand=random.random()
                    if rand < acc1:
                        latt[celltoflip] = -latt[celltoflip]
                if dE==8.0*J:
                    rand=random.random()
                    if rand<acc2:
                        latt[celltoflip] = -latt[celltoflip]
                if dE==12.0*J:
                    rand=random.random()
                    if rand<acc3:
                        latt[celltoflip] = -latt[celltoflip]
            if oldspin != latt[celltoflip]:
                M=M+2.0*latt[celltoflip]
            indMlist.append(M)
            steplist.append(nmctotal+1)
            nmctotal=nmctotal+1
        plt.plot(steplist, indMlist)
        plt.xscale('log')
        plt.show()
        print('choose nmceqstart:')
        nmceqstart=int(input())
        print('choose nmceqstart for the rest of the temperatures:')
        nmceq=nmc-int(input())
        filesuffix="D:/WPy64-31090/notebooks/file" +"_Tmin_"+ str(Tmin) +"_Tmax_"+ str(Tmax) +"_dT_" + str(dT) + "_Side Length_" + str(L) + "_Sites number_" + str(N) + "_nmc_" + str(nmc) + "_nmceqstart_" + str(nmceqstart) + "_nmceq_" + str(nmceq) + "_lattice type_" + lattice
        filedata = filesuffix+"_data.csv"
        for p in range(nmceqstart, nmc):
            Madd=Madd+abs(M)
        Madd=Madd/(nmc-nmceqstart)/N
        #if T==Tarray[1]:
    else:
        for j in range (0, nmc):
            nmccount=nmccount+1
            celltoflip=random.randint(0, N-1)
            oldspin=latt[celltoflip]
            dE=2.0*J*latt[celltoflip]*(latt[getsixthneighbor(celltoflip)]+latt[getfifthneighbor(celltoflip)]+latt[getfourthneighbor(celltoflip)]+latt[getthirdneighbor(celltoflip)]+latt[getsecondneighbor(celltoflip)]+latt[getfirstneighbor(celltoflip)])
            if dE<=0.0:
                latt[celltoflip] = -latt[celltoflip]
            else:
                if dE==4.0*J:
                    rand=random.random()
                    if rand < acc1:
                        latt[celltoflip] = -latt[celltoflip]
                if dE==8.0*J:
                    rand=random.random()
                    if rand<acc2:
                        latt[celltoflip] = -latt[celltoflip]
                if dE==12.0*J:
                    rand=random.random()
                    if rand<acc3:
                        latt[celltoflip] = -latt[celltoflip]
            if oldspin != latt[celltoflip]:
                M=M+2.0*latt[celltoflip]
            indMlist.append(M)
            steplist.append(nmctotal+1)
            nmctotal=nmctotal+1
            if j>nmc-nmceq:
                Madd=Madd+abs(M)
        Madd=Madd/(nmceq)/N
    Mlist.append(abs(Madd))
    print('expectation value of magnetization for temperature ', T, ' is ', abs(Madd))
    prevlatt=latt
    nmctotal=0
    #indarr = np.array([steplist, indMlist])
    #indfileplot = filesuffix + "Temp" + str(T) + "indmag_plot.pdf"
    steplist=[]
    indMlist=[]
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
