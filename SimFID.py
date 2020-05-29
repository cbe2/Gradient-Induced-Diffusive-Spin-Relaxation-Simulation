import matplotlib.pyplot as plt
import math
import numpy as np
from WalkerClasses import *
import time
start_time = time.time()

#All units cm,sec,Gauss

#np.random.seed(1)
#random.seed(1)

N=int(1e5) #number of walkers
steps=1000
print("{0:.1g}".format(N)+" spins simulating")


Gx=1e-3 ;Gy=Gx; Gz=Gx ;D=1.;gamma=(2*np.pi)*(3240) #Bfield params
BDx=np.power((4.*D)/(gamma*Gx),1./3.)
BDy=np.power((4.*D)/(gamma*Gy),1./3.)
BDz=np.power((4.*D)/(gamma*Gz),1./3.)
print("x Boundary Thickness: "+"{0:.3g}".format(BDx))
print("y Boundary Thickness: "+"{0:.3g}".format(BDy))
print("y Boundary Thickness: "+"{0:.3g}".format(BDz))

Lx=2.*BDx #x Length of cell cm
Ly=2.*BDy #y length of cell
Lz=2.*BDz

paramsDict={
'D':D, #cm^2/sec
'r0':[0,0], #to be altered below
'dt': 1e-3,
'L': [Lx,Ly]
}

f0=250 #choosen precession frequency
print("Frequency: "+str(f0))
B0=2.*np.pi*f0/gamma


#Defines magnetic field
def _Bz(r,Gx,Gy,Gz,B0):
    return B0+Gx*r[0]+Gy*r[1]#+Gz*r[2]

Bz=lambda x: _Bz(x,Gx,Gy,Gz,B0) #So Bz only takes one input


#creates distribution of walkers
walkers=[]
for i in range(N):
    paramsDict['r0']=[np.random.uniform(0,Lx),np.random.uniform(0,Ly)]
    walkers.append(BoxWalk(paramsDict))


data=[]
times=[]
#performs the simulation
for i in range(steps):
    thetaSum=0
    for W in walkers:
        thetaSum+=np.exp(1j*W.ptcl.theta)
        W.step(Bz)
    data.append(thetaSum/N)
    times.append(i*walkers[0].dt)

data=np.asarray(data)
times=np.asarray(times)

#Plots results
plt.plot(times,np.abs(data),label="Simulated")
plt.title(str(r'Simulation'))
plt.ylabel(r'$|\langle e^{i\theta}\rangle|$')
plt.xlabel("Time (sec)")
plt.legend()

#reports computation time
TimeTaken=(time.time() - start_time)/60. #in minutes
print("{0:.2g}".format(TimeTaken)+" minutes")
print("{0:.2g}".format(TimeTaken/float(steps*N)) + " min per step*particle")

plt.show()

#---------Writes data to file---------------
SaveData(times,np.abs(data),"SimulationData/Lx=2Ly=2SimAbs.txt")

path="SimulationData/"
fname="Lx=2Ly=2.txt"

data=np.imag(data)

f=open(path+fname,'w')

f.write("Spins Simulated:\t"+str(int(N))+'\n')
f.write("Sample Rate:\t" + str(1./paramsDict['dt'])+'\n')
f.write("Sample Size:\t" + str(steps)+'\n')
f.write("Read Time:\t" + str(float(steps)*paramsDict['dt'])+'\n')
f.write("Diffusion (cm^2/sec):\t" + str(paramsDict['D'])+'\n')
f.write("X-Gradient (G/cm):\t" + str(Gx)+'\n')
f.write("Y-Gradient (G/cm):\t" + str(Gy)+'\n')
f.write("Z-Gradient (G/cm):\t" + str(Gz)+'\n')
f.write("f0 (Hz):\t"+str(f0)+'\n')
f.write("Data Start:\r\n")

for d in data: f.write(str(d)+'\n')

f.close()
print("Created "+fname)
