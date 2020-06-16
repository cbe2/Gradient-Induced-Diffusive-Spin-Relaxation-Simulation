import matplotlib.pyplot as plt
import math
import numpy as np
from WalkerClasses import *
from HePhysicalProperties import *


#All units cm,sec,Gauss

T=0.4 #Kelvin
x3=1e-4 #3He concentration

Gx=1e-3 ;Gy=0; Gz=1e-3;gamma=(2*np.pi)*(3240) #Bfield params
D_b=1.#Db(x3,T) #bulk diffusion constant
D_s=1.#1e-2 #surface diffusion constant
print("Bulk Diffusion Const "+"{0:.3g}".format(D_b))
print("Surface Diffusion Const "+"{0:.3g}".format(D_s))
print("------------------------------------")

SBDx=np.power((4.*D_s)/(gamma*Gx),1./3.)
BDy=-1#np.power((4.*D)/(gamma*Gy),1./3.)
BBDz=np.power((4.*D_b)/(gamma*Gz),1./3.)
print("X Surface Boundary Thickness: "+"{0:.3g}".format(SBDx))
print("Y Surface Boundary Thickness: "+"{0:.3g}".format(BDy))
print("Z Bulk Boundary Thickness: "+"{0:.3g}".format(BBDz))
print("------------------------------------")

#Box dimensions (cm)
Lx=1.*SBDx
Ly=1.#1.*BDy
Lz=1.*BBDz

print("Lx: "+"{0:.3g}".format(Lx))
print("Ly: "+"{0:.3g}".format(Ly))
print("Lz: "+"{0:.3g}".format(Lz))
print("------------------------------------")

T2_bx=(120.*D_b)/(((gamma*Gx)**2)*(Lx**4)) #Fast limit T2 time
#T2_by=(120.*D_b)/(((gamma*Gy)**2)*(Ly**4)) #Fast limit T2 time
T2_bz=(120.*D_b)/(((gamma*Gz)**2)*(Lz**4)) #Fast limit T2 time

#T2_by=(120.*D_s)/(((gamma*Gy)**2)*(Ly**4)) #Fast limit T2 time



V=Lx*Ly*Lz
S=Lx*Ly

n4=2.18e22 # number density of 4He in cm^-3
N3=V*n4*x3

n_s=Getn_s(N3,S,V,T) #gets number of layers
print("Surface Layers: "+"{0:.3g}".format(n_s))

nl=6.4e14 #one layer surface number density (cm^-2)
print("Surface Atoms: "+"{0:.3g}".format(nl*n_s*S))

Sfrac = 1e-2#(n_s*nl)*S/N3 #Fraction on the surface N_s/(N_s+N_b)
print("Surface Fraction: "+"{0:.3g}".format(Sfrac))
print("------------------------------------")


dt=1e-3
Ttotal=T2_bx*5# sec

times=np.arange(0,Ttotal,dt)


f0=1./(4*dt) #choosen precession frequency
print("Frequency: "+str(f0))
B0=2.*np.pi*f0/gamma #holding field vaule to match f0

Bs=B0+Gx*Lx/2.+Gy*Ly/2. #effective surface field
Bb=B0+Gx*Lx/2.+Gy*Ly/2.+Gz*Lz/2. #effective bulk field

surfaceSig= np.exp(1j*gamma*Bs*times - times/T2_sx)
bulkSig=np.exp(1j*gamma*Bb*times - times/T2_bx-times/T2_bz)


data=surfaceSig*Sfrac+bulkSig

#Plots results
plt.plot(times,np.abs(data),label="Simulated")
plt.title(str(r'Simulation'))
plt.ylabel(r'$|\langle e^{i\theta}\rangle|$')
plt.xlabel("Time (sec)")
plt.legend()

plt.show()

#---------Writes data to file---------------
SaveData(times,np.abs(data),"SimulationData/HighDS=1e-2GAbs.txt")

path="SimulationData/"
fname="HighDS=1e-2G.txt"

data=np.imag(data)

f=open(path+fname,'w')


f.write("Sample Rate:\t" + str(1./dt)+'\n')
f.write("Sample Size:\t" + str(len(times))+'\n')
f.write("Read Time:\t" + str(float(Ttotal))+'\n')
f.write("Surface Diffusion (cm^2/sec):\t" + str(D_s)+'\n')
f.write("Bulk Diffusion (cm^2/sec):\t" + str(D_b)+'\n')
f.write("Surface Fraction:\t" + str(Sfrac)+'\n')
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
