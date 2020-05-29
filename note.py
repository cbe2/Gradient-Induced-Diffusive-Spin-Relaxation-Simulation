import matplotlib.pyplot as plt
import math
import numpy as np
from WalkerClasses import *

D=1.0
Gx=1e-3 ;Gy=Gx*1.5; Gz=Gx*0.5
gamma=(2*np.pi)*(3240) #radian/Gauss
BDx=np.power((4.*D)/(gamma*Gx),1./3.)
L=1.0*BDx #Length of cell cm

time,data=LoadData("SimulationData/FreeGxSimAbs.txt")
time2,data2=LoadData("SimulationData/SemiFreeGxSimAbs.txt")
time3,data3=LoadData("SimulationData/SemiFreeGxGyGzZ=10SimAbs.txt")

T2=(120.*D)/(((gamma*Gx)**2)*(L**4)) #Fast limit T2 time
#y=np.exp(-1.25*time/T2)
#Axy=((gamma*Gx*D)**2)/3. #Slow limit for G=xy
#Ax=(D*(gamma*Gx)**2)/3.
# Ay=(3./8.-8./(9.*np.pi))*(D*(gamma*Gy)**2)/3.
# Az=(3./8.-8./(9.*np.pi))*(D*(gamma*Gz)**2)/3.
#y=np.exp(-A*(time)**4) # for G=xy
#y=np.exp(-3*Ax*(time)**3)
#
#time2,y=LoadData("L=3_Kubo.txt")

#y=np.exp(-time*np.power(gamma*G,2)*np.power(L,4)/(120.*D))


plt.plot(time,data,label=r'Curve A: $B_z=B_0+Gx$, $r_0=[\infty,\infty,\infty]$',linewidth=3)
plt.plot(time2,data2,label=r'Curve B: $B_z=B_0+Gx$, $r_0=[0,0,0]$',linewidth=3)
plt.plot(time3,data3,label=r'Curve C: $B_z=B_0+Gx+Gy+Gz$, $r_0=[0,0,\infty]$',linewidth=3)
data=np.hstack([data,np.zeros(200)])

plt.plot(time2,data*data2*data2,linestyle='dotted',label=r'$(A)\times(B)\times(B)$',linewidth=3)
plt.title(str(r'$10^5$ Simulated Spins, SemiFree Diffusion'))
plt.ylabel(r'$|\langle e^{i\theta}\rangle|$')
plt.xlabel("Time (sec)")
plt.xlim([0,0.4])
#plt.xlim([0,1.2])
plt.legend()


plt.show()
