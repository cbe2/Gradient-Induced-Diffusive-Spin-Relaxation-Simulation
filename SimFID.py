import matplotlib.pyplot as plt
import math
import numpy as np
from WalkerClasses import *
from HePhysicalProperties import *
import time
start_time = time.time()

#All units cm,sec,Gauss

"""
Run check list:
Test Run,N,steps,gradients,D,Lengths,f0,distribution,walker type,filename
"""

#np.random.seed(1)
#random.seed(1)

N=int(1e6) #number of walkers
steps=int(1.5e3)
print("{0:.1g}".format(N)+" spins simulating")

x3=1e-4
T=0.2

Gx=4e-4;Gy=0; Gz=-1 ;gamma=(2*np.pi)*(3240) #Bfield params
D=2e-2#Db(x3,T)
BDx=bThickness(D,Gx,gamma)#np.power((4.*D)/(gamma*Gx),1./3.)
BDy=bThickness(D,Gy,gamma)#np.power((4.*D)/(gamma*Gy),1./3.)
BDz=bThickness(D,Gz,gamma)#np.power((4.*D)/(gamma*Gz),1./3.)

print("D="+"{0:.3g}".format(D))
print("------------------------------------")
print("x Boundary Thickness: "+"{0:.3g}".format(BDx))
print("y Boundary Thickness: "+"{0:.3g}".format(BDy))
print("y Boundary Thickness: "+"{0:.3g}".format(BDz))
print("------------------------------------")

Lx=4.#BDx*4. #x Length of cell cm
Ly=4.#2.*Lx#10.*BDy #y length of cell
Lz=4.#Lx#10.*BDz

print("Lx: "+"{0:.3g}".format(Lx))
print("Ly: "+"{0:.3g}".format(Ly))
print("Lz: "+"{0:.3g}".format(Lz))
print("------------------------------------")



paramsDict={
'D':D, #cm^2/sec
'r0':[0,0], #to be altered below
'dt': 5*np.power(2,12)*1e-7,
'L': [Lx,Ly]
}

f0=1./(4.*paramsDict['dt']) #choosen precession frequency at position x=y=z=0
print("Frequency: "+str(f0))
B0=2.*np.pi*f0/gamma# corresponding field


#Defines magnetic field
def _Bz(r,Gx,Gy,Gz,B0,Lx):

    return B0+Gx*r[0]+Gy*r[1]#+Gz*r[2]

Bz=lambda x: _Bz(x,Gx,Gy,Gz,B0,Lx) #So Bz only takes one input


#creates distribution of walkers
walkers=[]
for i in range(N):
    paramsDict['r0']=[np.random.uniform(0,Lx),np.random.uniform(0,Ly)]#np.random.uniform(0,-Lz)
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
fname="6-17_note_surface_Dmax"


#SaveData(times,np.abs(data),"SimulationData/"+fname+ "Abs.txt")
#
path="SimulationData/"
fname=fname+".txt"

#data=np.real(data)

f=open(path+fname,'w')

f.write("Notes:\t Simulates surface signal for design in the 6/17 note with maximum diffusion limit. \
  Saves full complex data instead of just imaginary part"+'\n')
f.write("Spins Simulated:\t"+str(int(N))+'\n')
f.write("Sample Rate:\t" + str(1./paramsDict['dt'])+'\n')
f.write("Sample Size:\t" + str(steps)+'\n')
f.write("Read Time:\t" + str(float(steps)*paramsDict['dt'])+'\n')
f.write("Diffusion (cm^2/sec):\t" + str(paramsDict['D'])+'\n')
f.write("X-Gradient (G/cm):\t" + str(Gx)+'\n')
f.write("Y-Gradient (G/cm):\t" + str(Gy)+'\n')
f.write("Z-Gradient (G/cm):\t" + str(Gz)+'\n')
f.write("Lx (cm):\t" + str(Lx)+'\n')
f.write("Ly (cm):\t" + str(Ly)+'\n')
f.write("Lz (cm):\t" + str(Lz)+'\n')
f.write("f0 (Hz):\t"+str(f0)+'\n')
f.write("Data Start:\r\n")

for d in data: f.write(str(d)+'\n')

f.close()
print("Created "+fname)
