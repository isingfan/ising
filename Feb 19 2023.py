import math
from math import *
import numpy as np
from numpy import *
from numpy.linalg import *
import random
import matplotlib.pyplot as plt
import pandas as pd
J=1 #interaction energy
D=14
N=D**2 #number of cells
nmc=100000 #number of monte carlo steps for each temperature value
nmceq=95000 #number of states that are counted for expectation value
#k=1.380649*(10**-23) #boltzmann constant in m^2 kg s^-2 K^-1
k=1
Tlist=[]
Mlist=[]

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

def getleftneighbor(cell):
    index = (cell - 1)%N
    return index
def getrightneighbor(cell):
    index = (cell + 1)%N
    return index
def getbottomneighbor(cell):
    index = (cell - D)%N
    return index
def gettopneighbor(cell):
    index = (cell + D)%N
    return index
def getleftneighborvalue(cell):
    value=latt[getleftneighbor(cell)]
    return value
def getrightneighborvalue(cell):
    value=latt[getrightneighbor(cell)]
    return value
def gettopneighborvalue(cell):
    value=latt[gettopneighbor(cell)]
    return value
def getbottomneighborvalue(cell):
    value=latt[getbottomneighbor(cell)]
    return value

steplist=[]
indMlist=[]
#for custom temperature values, put values in this:
Tarray=[]
#for range of temperatures, use this:
for i in range (36, 62):
    Tarray.append(i/20)
nmctotal=0
for T in Tarray:
    print('Temperature is:', T)
    if T==Tarray[0]:
        latt=2*C-1
        print('first magnetization is: ', getmagnetization(latt))
    else:
        latt=prevlatt
    counta=0
    countb=0
    countc=0
    beta=1/(k*T)
    acc1=math.exp(-4*beta*J)
    acc2=math.exp(-8*beta*J)
    M=0.0
    Madd=0.0
    nmccount=0
    for i in range (0, N):
        M=M+latt[i]
    if T==Tarray[0]:
        for j in range (0, round(-abs(T-2.6)+3)*nmc):
            nmccount=nmccount+1
            celltoflip=random.randint(0, N-1)
            oldspin=latt[celltoflip]
            dE=2.0*J*latt[celltoflip]*(latt[gettopneighbor(celltoflip)]+latt[getbottomneighbor(celltoflip)]+latt[getrightneighbor(celltoflip)]+latt[getleftneighbor(celltoflip)])
            if dE<=0.0:
                latt[celltoflip] = -latt[celltoflip]
                counta=counta+1
            else:
                #if dE>3.9*J and dE<4.1*J:
                if dE==4.0*J:
                    rand=random.random()
                    countb=countb+1
                    if rand < acc1:
                        latt[celltoflip] = -latt[celltoflip]
                #if dE<8.1*J and dE>7.9*J:
                if dE==8.0*J:
                    rand=random.random()
                    countc=countc+1
                    if rand<acc2:
                        latt[celltoflip] = -latt[celltoflip]
            if oldspin != latt[celltoflip]:
                M=M+2.0*latt[celltoflip]
            indMlist.append(M)
            steplist.append(nmctotal+1)
            nmctotal=nmctotal+1
        plt.plot(steplist, indMlist)
        plt.show()
        print('choose nmceqstart:')
        nmceqstart=int(input())
        for p in range(nmceqstart, nmc):
            Madd=Madd+M
        Madd=Madd/(nmc-nmceqstart)/D**2
    else:
        for j in range (0, nmc):
            nmccount=nmccount+1
            celltoflip=random.randint(0, N-1)
            oldspin=latt[celltoflip]
            dE=2.0*J*latt[celltoflip]*(latt[gettopneighbor(celltoflip)]+latt[getbottomneighbor(celltoflip)]+latt[getrightneighbor(celltoflip)]+latt[getleftneighbor(celltoflip)])
            if dE<=0.0:
                latt[celltoflip] = -latt[celltoflip]
                counta=counta+1
            else:
                #if dE>3.9*J and dE<4.1*J:
                if dE==4.0*J:
                    rand=random.random()
                    countb=countb+1
                    if rand < acc1:
                        latt[celltoflip] = -latt[celltoflip]
                #if dE<8.1*J and dE>7.9*J:
                if dE==8.0*J:
                    rand=random.random()
                    countc=countc+1
                    if rand<acc2:
                        latt[celltoflip] = -latt[celltoflip]
            if oldspin != latt[celltoflip]:
                M=M+2.0*latt[celltoflip]
            indMlist.append(M)
            steplist.append(nmctotal+1)
            nmctotal=nmctotal+1
            if j>nmc-nmceq:
                Madd=Madd+M
        Madd=Madd/(nmceq)/D**2
    Mlist.append(abs(Madd))
    print('monte carlo steps performed in this temp: ', nmccount)
    #print('new last lattice magnetization is:', getmagnetization(latt))
    print('expectation value of magnetization for temperature ', T, ' is ', abs(Madd))
    print('counta, countb, and countc are:', counta, countb, countc)
    counta=0
    countb=0
    countc=0
    prevlatt=latt
    nmctotal=0
    steplist=[]
    indMlist=[]
    Tlist.append(T)
df = pd.DataFrame(Tlist).T
df.to_excel(excel_writer = "C:/Users/Shane Armstrong/Downloads/WPy64-31090/notebooks/Tlist.xlsx")
df = pd.DataFrame(Mlist).T
df.to_excel(excel_writer = "C:/Users/Shane Armstrong/Downloads/WPy64-31090/notebooks/Mlist.xlsx")
plt.plot(Tlist, Mlist, 'o')
plt.show()
