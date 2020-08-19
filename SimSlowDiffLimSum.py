import matplotlib.pyplot as plt
import math
import numpy as np
from WalkerClasses import *
from HePhysicalProperties import *


#All units cm,sec,Gauss

"""
This code is for producing signals in the slow diffusion limit
"""

T=0.2 #Kelvin
x3=1e-4 #3He concentration

Gx=4e-4 ;Gy=0; Gz=1;gamma=(2*np.pi)*(3240) #Bfield params
D_b=0 #bulk diffusion constant
D_s=0 #surface diffusion constant
print("------------------------------------")
print("Bulk Diffusion Const "+"{0:.3g}".format(D_b))
print("Surface Diffusion Const "+"{0:.3g}".format(D_s))
print("------------------------------------")

SBDx=bThickness(D_s,Gx,gamma)#np.power((4.*D_s)/(gamma*np.abs(Gx)),1./3.)
SBDy=bThickness(D_s,Gy,gamma)#np.power((4.*D_s)/(gamma*np.abs(Gy)),1./3.)
BBDz=bThickness(D_b,Gz,gamma)#np.power((4.*D_b)/(gamma*np.abs(Gz)),1./3.)
print("X Surface Boundary Thickness: "+"{0:.3g}".format(SBDx))
print("Y Surface Boundary Thickness: "+"{0:.3g}".format(SBDy))
print("Z Bulk Boundary Thickness: "+"{0:.3g}".format(BBDz))
print("------------------------------------")

#Box dimensions (cm)
Lx=4.#1.*SBDx
Ly=4.#1.*SBDy
Lz=4.#1.*BBDz

print("Lx: "+"{0:.3g}".format(Lx))
print("Ly: "+"{0:.3g}".format(Ly))
print("Lz: "+"{0:.3g}".format(Lz))
print("------------------------------------")

#reciprcol of T2 times (in sec^-1)
#bulk
rT2_x=gamma*Gx*Lx/2.
rT2_y=gamma*Gy*Ly/2.#(((gamma*Gy)**2)*(Ly**4))/(120.*D_b)
rT2_z=gamma*Gz*Lz/2.#(((gamma*Gz)**2)*(Lz**4))/(120.*D_b)


V=Lx*Ly*Lz
S=Lx*Ly

n4=2.18e22 # number density of 4He in cm^-3
N3=V*n4*x3

n_s=Getn_s(N3,S,V,T) #gets number of layers, really a fct of concentration
print("Surface Layers: "+"{0:.3g}".format(n_s))

nl=6.4e14 #one layer surface number density (cm^-2)
print("Surface Atoms: "+"{0:.3g}".format(nl*n_s*S))

Sfrac = (n_s*nl)*S/N3 #Fraction on the surface N_s/(N_s+N_b)
print("Surface Fraction: "+"{0:.3g}".format(Sfrac))
print("------------------------------------")

SR=2000000.0
SS=2000
dt=1./SR
Ttotal=dt*SS# sec
print("Total time: "+"{0:.3g}".format(Ttotal))

times=np.arange(0,Ttotal,dt)


f0=1./(4*dt) #choosen precession frequency at x=y=z=0
print("Frequency: "+str(f0))
B0=2.*np.pi*f0/gamma #holding field vaule to match f0

Bs=B0+Gx*Lx/2.+Gy*Ly/2. #effective surface field (at z=0)
Bb=B0+Gx*Lx/2.+Gy*Ly/2.-Gz*Lz/2. #effective bulk field, - sign because z is negative in bulk

surfaceSig= np.exp(1j*gamma*Bs*times)*np.sinc(rT2_x*times/np.pi)*np.sinc(rT2_y*times/np.pi)
bulkSig=np.exp(1j*gamma*Bb*times)*np.sinc(rT2_x*times/np.pi)*np.sinc(rT2_y*times/np.pi)*np.sinc(rT2_z*times/np.pi)

#normallize these so at t=0 their sum is equal to 1
#surfaceSig=surfaceSig*Sfrac
#bulkSig=bulkSig*(1.-Sfrac)

#the combined signal
data=bulkSig#surfaceSig

#Plots results
plt.plot(times,np.abs(data),label="combined")
#plt.plot(times,np.abs(surfaceSig),label="surface")
#plt.plot(times,np.abs(bulkSig),label="bulk")
plt.title(str(r'Simulation'))
plt.ylabel(r'$|\langle e^{i\theta}\rangle|$')
plt.xlabel("Time (sec)")
plt.legend()

plt.show()

#---------Writes data to file---------------
#SaveData(times,np.abs(data),"SimulationData/Sf=3e-3Abs.txt")

path="SimulationData/"
fname="6-17_note_bulk_D=0.txt"

#data=np.imag(data)

f=open(path+fname,'w')

f.write("Notes:\t  SBulk slow diffusion signal (no surface) for design in the 6/17 note for the D=0 case. \
  Saves full complex data instead of just imaginary part"+'\n')
f.write("Legend Name:\tD=0\n")
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
