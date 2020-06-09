import matplotlib.pyplot as plt
import math
import cmath
import numpy as np
from scipy import special as sp
import sys
sys.path.append('../')
from KuboRedfieldFunctions import *
from HePhysicalProperties import *

'''
This code plots Diffusion Constants at the surface

See the plots fom this code at
https://www.evernote.com/l/AuEmmLn7gv1CIZv05EO63PCszdrEk36WoQQ

'''

#Temperature range
TL=np.linspace(0.001,0.0021,100)
TH=np.linspace(.21,0.5,100)


#T=np.logspace(np.log10(1e-7),np.log10(1),1000)


print(LowTDs(0.01,0.0021))
print(LowTDs(0.1,0.021))

plt.semilogy(TL,LowTDs(0.01,TL),linewidth=3,label="Questionable Low Temp Extrapolation")
plt.semilogy(TH,HighTDs(0.01,TH),linewidth=3,label="Rough High Temp Prediction")
#plt.scatter([TL[-1],TH[0]],[LowTDs(0.1,TL[-1]),HighTDs(0.1,TH[0])],c='r')


plt.legend()
plt.ylabel(r"D (cm$^2$/sec)")
plt.xlabel("Temperature (Kelvin)")
plt.title(r"Surface Diffusion for $l_3=0.01$ layers")
plt.grid()

plt.show()
