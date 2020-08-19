import matplotlib.pyplot as plt
import math
import numpy as np
from WalkerClasses import *
from DFTSuite import *

'''
For plotting the norm as a function of time
'''


# D=1.0
# Gx=1e-3/np.sqrt(2) ;Gy=Gx*1.5; Gz=Gx*0.5
# gamma=(2*np.pi)*(3240) #radian/Gauss
# BDx=np.power((4.*D)/(gamma*Gx),1./3.)
# L=1.0*BDx #Length of cell cm


fnames=getFilenames("SimulationData/PlotFolder/")
signals=[]
for f in fnames:
    signals.append(FIDdata(f))


for S in signals:
    #if S.LN=="D=0": times=np.arange(0,S.SS/S.SR,1./S.SR)[:-2]
    times=np.arange(0,S.SS/S.SR,1./S.SR)[:-1]
    plt.plot(times,S.norms,label=S.LN,linestyle='-')#label=S.name[:-4]

#T2=(120.*D)/(((gamma*Gx)**2)*(L**4)) #Fast limit T2 time
#y=np.exp(-1.25*time/T2)

plt.title(r' Bulk Signal D=2.7 cm$^2$/sec Vs D=0 cm$^2$/sec, L=4cm, G=1 G/cm')
plt.ylabel(r'$|\langle e^{i\theta}\rangle|$')
plt.xlabel("Time (sec)")
#plt.xlim([0,0.8])
#plt.xlim([0,1.2])
plt.legend()


plt.show()
