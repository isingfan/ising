import math
from math import *
import numpy as np
from numpy import *
from numpy.linalg import *
import random
import matplotlib.pyplot as plt
import pandas as pd
J = 1 #interaction energy
lattice = "rhombus"
D = 30
N = D**2 #number of cells
nmc = 2000000 #number of monte carlo steps for each temperature value
nmceq = 600000 #number of states that are counted for expectation value
Tmin = 2.3
Tmax = 3.3
dT = 0.05
#k=1.380649*(10**-23) #boltzmann constant in m^2 kg s^-2 K^-1
k=1

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
    index = (cell - 1)%N
    return index
def getsecondneighbor(cell):
    index = (cell + 1)%N
    return index
def getthirdneighbor(cell):
    index = (cell - D)%N
    return index
def getfourthneighbor(cell):
    index = (cell + D)%N
    return index
def getfifthneighbor(cell):
    index = (cell + (D+1))%N
    return index
def getsixthneighbor(cell):
    index = (cell + (D-1))%N
    return index

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
print('Dimension is:', D)
print('nmc is:', nmc)
print('nmceq is:', nmceq)
for T in Tarray:
    print('Temperature is:', T)
    if T==Tarray[0]:
        latt=2*C-1
    else:
        latt=prevlatt
    beta=1/(k*T)
    acc1=math.exp(-4*beta*J)
    acc2=math.exp(-8*beta*J)
    acc3=math.exp(-12*beta*J)
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
        filesuffix="D:/WPy64-31090/notebooks/file" +"_Tmin_"+ str(Tmin) +"_Tmax_"+ str(Tmax) +"_dT_" + str(dT) + "_Dimension_" + str(D) + "_nmc_" + str(nmc) + "_nmceqstart_" + str(nmceqstart) + "_nmceq_" + str(nmceq) + "_lattice type_" + lattice
        filedata = filesuffix+"_data.csv"
        for p in range(nmceqstart, nmc):
            Madd=Madd+abs(M)
        Madd=Madd/(nmc-nmceqstart)/D**2
        #if T==Tarray[1]:
            #print('choose nmceq for the rest of the temperatures:')
            #nmceq=int(input())
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
        Madd=Madd/(nmceq)/D**2
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
plt.plot(Tlist, Mlist, 'o')
plt.savefig(fileplot)
print('Simulation complete.')
#plt.show()
