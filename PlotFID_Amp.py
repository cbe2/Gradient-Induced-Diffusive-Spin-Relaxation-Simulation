import matplotlib.pyplot as plt
import math
import numpy as np
from WalkerClasses import *

D=1.0
Gx=1e-3 ;Gy=Gx*1.5; Gz=Gx*0.5
gamma=(2*np.pi)*(3240) #radian/Gauss
BDx=np.power((4.*D)/(gamma*Gx),1./3.)
L=1.0*BDx #Length of cell cm

time,data=LoadData("SimulationData/Lx=1Ly=1Lz=1SimAbs.txt")
time2,data2=LoadData("SimulationData/Lx=2Ly=2Lz=2SimAbs.txt")
time3,data3=LoadData("SimulationData/Lx=3Ly=3Lz=3SimAbs.txt")
time4,data4=LoadData("SimulationData/Lx=4Ly=4Lz=4SimAbs.txt")

T2=(120.*D)/(((gamma*Gx)**2)*(L**4)) #Fast limit T2 time
#y=np.exp(-1.25*time/T2)


plt.plot(time,data,label=r'$L_x=L_y=L_z=1$')
plt.plot(time2,data2,label=r'$L_x=L_y=L_z=2$')
plt.plot(time3,data3,label=r'$L_x=L_y=L_z=3$')
plt.plot(time4,data4,label=r'$L_x=L_y=L_z=4$')
plt.title(str(r'3D Gradients Diagonal to a Cube'))
plt.ylabel(r'$|\langle e^{i\theta}\rangle|$')
plt.xlabel("Time (sec)")
plt.xlim([0,0.8])
#plt.xlim([0,1.2])
plt.legend()


plt.show()
