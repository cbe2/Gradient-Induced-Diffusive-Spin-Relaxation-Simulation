import matplotlib.pyplot as plt
import math
import cmath
import numpy as np
from scipy import special as sp
import sys
sys.path.append('../')
from KuboRedfieldFunctions import *

'''
This code plots the surface occupation
See the plots fom this code and further explaination at
https://www.evernote.com/l/AuEmmLn7gv1CIZv05EO63PCszdrEk36WoQQ
'''


L= 0.1 # cm
A=1.#L**2 # area of free surface (cm^2)
V=1.#L**3 #Volume of free surface (cc)
T=0.4 #Kelvin

#number density of 3He
n4=2.18e22 # in cm^-3

nl=6.4e14 #one layer surface number density (cm^-2)

#bulk number densities
n_3=np.logspace(np.log10(n4*1e-6),np.log10(n4*1e-2),1000)

#total 3He atoms
N_3=n_3*V#np.logspace(np.log10(n4*1e-6*V),np.log10(n4*1e-2*V),1000)

#computes results for classical 2D gas
Sfrac=GetSfrac(A,V,T)
n_sFrac=Sfrac*n_3*V/A #in cm^-2
n_sFrac=n_sFrac/nl # #in layers

#This tests the significance of different effects
plt.loglog(n_3/n4,n_sFrac,label="Non Degenerate, Non Interacting")
plt.loglog(n_3/n4,Getn_s_(n_3,T,0),label="Degenerate, Non Interacting")
plt.loglog(n_3/n4,Getn_s(N_3,A,V,T),label="Degenerate, Interacting")


plt.legend()
plt.ylabel("Surface Number Density (in Layers)")
plt.xlabel("Bulk 3He concentration")
plt.title("1 cc cube of liquid at 100mK")
plt.grid()

plt.show()
