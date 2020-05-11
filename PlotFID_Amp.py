import matplotlib.pyplot as plt
import math
import numpy as np
from WalkerClasses import *
from BoundaryDiffFunctions import *
import time

D=1.0
G=1e-3 #Gauss/cm
gamma=(2*np.pi)*(3240) #radian/Gauss
BDx=np.power((4.*D)/(gamma*G),1./3.)
L=1.0*BDx #Length of cell cm

time,data=LoadData("N1e5G3D1L10.txt")

time2,y=LoadData("L=10_f0=250absAn.txt")

#y=np.exp(-time*np.power(gamma*G,2)*np.power(L,4)/(120.*D))

plt.plot(time,data,label=r'$10^5$ Simulated Spins')
plt.plot(time2,y,label="Kubo-Redfield Theory" )
plt.title(str(r'$L=10$'))
plt.ylabel(r'$|\langle e^{i\theta}\rangle|$')
plt.xlabel("Time (sec)")
plt.legend()


plt.show()


# path=""#"/Users/cameronerickson/Desktop/Academics/NPL/DRMRI/SignalSim/FakeSig/Signals/"
# fname="L=x.txt"
#
# data=np.imag(data)
#
# f=open(path+fname,'w')
#
# f.write("Sample Rate:\t" + str(SR)+'\n')
# f.write("Sample Size:\t" + str(SS)+'\n')
# f.write("Read Time:\t" + str(float(SS)/SR)+'\n')
# f.write("Diffusion (cm^2/sec):\t" + str(D)+'\n')
# f.write("Gradient (G/cm):\t" + str(G)+'\n')
# f.write("Data Start:\r\n")
#
# for d in data: f.write(str(d)+'\n')
#
# f.close()
#
# print("Created "+fname)


# path=""#"/Users/cameronerickson/Desktop/Academics/NPL/DRMRI/SignalSim/FakeSig/Signals/"
# fname="N5G1WallAnl.txt"
#
# y=np.imag(y)
#
# f=open(path+fname,'w')
#
# f.write("Sample Rate:\t" + str(SR)+'\n')
# f.write("Sample Size:\t" + str(SS)+'\n')
# f.write("Read Time:\t" + str(float(SS)/SR)+'\n')
# f.write("Diffusion (cm^2/sec):\t" + str(D)+'\n')
# f.write("Gradient (G/cm):\t" + str(G)+'\n')
# f.write("Data Start:\r\n")
#
# for d in y: f.write(str(d)+'\n')
#
# f.close()
#
# print("Created "+fname)
