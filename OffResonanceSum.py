import matplotlib.pyplot as plt
import math
import numpy as np
from WalkerClasses import *
from HePhysicalProperties import *


#All units cm,sec,Gauss

"""
This code is for the note XXX
producing signals in the slow diffusion limit
"""

T=0.2 #Kelvin
x3=1e-4 #3He concentration

Gx=4e-4 ;Gy=0; Gz=0.001;gamma=(2*np.pi)*(3240) #Bfield params
theta=np.radians(90) #degree of kick in radians
B1=0.01 #strength of kicking field
tp=np.abs(theta/(gamma*B1))

D_b=1e4#Db(x3,T) #bulk diffusion constant
D_s=2e-2 #surface diffusion constant
print("------------------------------------")
print("Temp (K) "+"{0:.3g}".format(T))
print("Concentration "+"{0:.3g}".format(x3))
print("Z-gradient "+"{0:.2g}".format(Gz))
print("------------------------------------")
print("X-gradient "+"{0:.1g}".format(Gx))
print("Y-gradient "+"{0:.1g}".format(Gy))
print("Z-gradient "+"{0:.1g}".format(Gz))
print("------------------------------------")
print("Bulk Diffusion Const (cm^2/sec) "+"{0:.3g}".format(D_b))
print("Surface Diffusion Const (cm^2/sec) "+"{0:.3g}".format(D_s))
print("------------------------------------")

SBDx=bThickness(D_s,Gx,gamma)
SBDy=bThickness(D_s,Gy,gamma)
BBDz=bThickness(D_b,Gz,gamma)
print("X Surface Boundary Thickness: "+"{0:.3g}".format(SBDx))
print("Y Surface Boundary Thickness: "+"{0:.3g}".format(SBDy))
print("Z Bulk Boundary Thickness: "+"{0:.3g}".format(BBDz))
print("------------------------------------")

SRx=rThickness(B1,Gx,theta)
SRy=rThickness(B1,Gy,theta)
BRz=rThickness(B1,Gz,theta)
print("X Surface Resonance Thickness: "+"{0:.3g}".format(SRx))
print("Y Surface Resonance Thickness: "+"{0:.3g}".format(SRy))
print("Z Bulk Resonance Thickness: "+"{0:.3g}".format(BRz))
print("------------------------------------")

#Box dimensions (cm)
Lx=3.#1.*SBDx
Ly=3.#1.*SBDy
Lz=2.#length

print("Lx: "+"{0:.3g}".format(Lx))
print("Ly: "+"{0:.3g}".format(Ly))
print("Lz: "+"{0:.3g}".format(Lz))
print("------------------------------------")


print("X spatial dephasing time scale (sec): "+"{0:.3g}".format(Stime(Gx,gamma,Lx)))
print("Y spatial dephasing time scale (sec): "+"{0:.3g}".format(Stime(Gy,gamma,Ly)))
print("Z spatial dephasing time scale (sec): "+"{0:.3g}".format(Stime(Gz,gamma,Lz)))
print("------------------------------------")


print("X Surface diffusive time scale (sec): "+"{0:.3g}".format(Dtime(Gx,gamma,D_s)))
print("Y Surface diffusive time scale (sec): "+"{0:.3g}".format(Dtime(Gy,gamma,D_s)))
print("Z Bulk diffusive time scale (sec): "+"{0:.3g}".format(Dtime(Gz,gamma,D_b)))
print("------------------------------------")

print("Kicking field strength (G): "+"{0:.3g}".format(B1))
print("Kicking pulse time (milli-sec): "+"{0:.3g}".format(tp*1000))
print("------------------------------------")


V=Lx*Ly*Lz #volume of bulk
S=Lx*Ly #surface area of free surface

n4=2.18e22 # number density of 4He in cm^-3
N3=V*n4*x3

n_s=Getn_s(N3,S,V,T) #gets number of layers, really a fct of concentration
print("Surface Layers: "+"{0:.3g}".format(n_s))

nl=6.4e14 #one layer surface number density (cm^-2)
Ns=nl*n_s*S
print("Surface Atoms: "+"{0:.3g}".format(Ns))
print("Bulk Atoms: "+"{0:.3g}".format(N3))

Sfrac = (n_s*nl)*S/N3 #Fraction on the surface N_s/(N_s+N_b)
print("Surface Fraction: "+"{0:.3g}".format(Sfrac))
print("------------------------------------")

SR=50000.0
SS=50000.
dt=1./SR
Ttotal=dt*SS# sec
print("Total time: "+"{0:.3g}".format(Ttotal))

times=np.arange(0,Ttotal,dt)


f0=1./(4*dt) #choosen precession frequency at x=y=z=0
B0=50#2.*np.pi*f0/gamma #holding field vaule to match f0
f0=gamma*B0/(2.*np.pi)
print("Frequency: "+"{0:.3g}".format(f0))

Bs=B0+Gx*Lx/2.+Gy*Ly/2. #effective surface field (at z=0)
Bb=Bs

#reciprcol of T2 times (in sec^-1)
#bulk
rT2_x=gamma*Gx*Lx/2.
rT2_y=gamma*Gy*Ly/2.#(((gamma*Gy)**2)*(Ly**4))/(120.*D_b)
rT2_z=gamma*Gz*Lz/2.#(((gamma*Gz)**2)*(Lz**4))/(120.*D_b)

surfaceSig=np.sin(theta)*np.exp(1j*gamma*Bs*times)*np.sinc(rT2_x*times/np.pi)*np.sinc(rT2_y*times/np.pi)
l=Lz/BRz #length in resonance thicknesss
taus=times/tp #time scaled by pulse duration
Egamma=0.57721566490153
bulkSig=np.sin(theta)*np.exp(1j*gamma*Bb*times)*(1j/l)*(1-np.sinc(l/(2.*np.pi))*np.exp(1j*(taus+0.5)*l))*np.log(1.+1./(taus+1./(l*np.exp(Egamma))))
#np.exp(1j*gamma*Bb*times)*np.sinc(rT2_x*times/np.pi)*np.sinc(rT2_y*times/np.pi)*np.sinc(rT2_z*times/np.pi)

#Add diffusive decay components:
surfaceSig=surfaceSig*np.exp(-(1./3.)*((gamma*Gx)**2)*D_s*times**3)
surfaceSig=surfaceSig*np.exp(-(1./3.)*((gamma*Gy)**2)*D_s*times**3)
bulkSig=bulkSig*np.exp(-(1./3.)*((gamma*Gz)**2)*D_b*times**3)
bulkSig=bulkSig*np.exp(-(1./3.)*((gamma*Gx)**2)*D_b*times**3)
bulkSig=bulkSig*np.exp(-(1./3.)*((gamma*Gy)**2)*D_b*times**3)


#normallize these so at t=0 their sum is equal to 1
Adv=6.02214076e23
surfaceSig=surfaceSig*Ns/Adv#*Sfrac/(Sfrac+1.)
bulkSig=bulkSig*N3/Adv#/(1.+Sfrac)

#the combined signal
data=surfaceSig+bulkSig#

#Plots results
plt.plot(times,np.abs(data),label="combined")
#plt.plot(times,np.exp(-(1./3.)*((gamma*Gx)**2)*D_b*times**3))
plt.title(str(r'Simulation'))
plt.ylabel(r'$|\langle e^{i\theta}\rangle|$')
plt.xlabel("Time (sec)")
plt.grid()
plt.figure()
plt.title(str(r'Surface'))
plt.ylabel(r'$|\langle e^{i\theta}\rangle|$')
plt.xlabel("Time (sec)")
plt.grid()
plt.plot(times,np.abs(surfaceSig),label="surface")
plt.figure()
plt.title(str(r'Bulk'))
plt.ylabel(r'$|\langle e^{i\theta}\rangle|$')
plt.xlabel("Time (sec)")
plt.grid()
plt.plot(times,np.abs(bulkSig),label="bulk")
plt.legend()

#plt.show()

#---------Writes data to file---------------
#SaveData(times,np.abs(data),"SimulationData/Sf=3e-3Abs.txt")

path="SimulationData/"
fname="IgnoreThis.txt"

#data=np.imag(data)

# f=open(path+fname,'w')
#
# f.write("Notes:\t  First calculated signal taking into account off-resonance effects \
#   Diffusion is assumed to be an independent process. B0 field raised to 50 Gauss"+'\n')
# f.write("Legend Name:\tG=0.01 G/cm\n")
# f.write("Sample Rate:\t" + str(1./dt)+'\n')
# f.write("Sample Size:\t" + str(len(times))+'\n')
# f.write("Read Time:\t" + str(float(Ttotal))+'\n')
# f.write("Concentration (x3):\t" + str(x3)+'\n')
# f.write("Temperature (K):\t" + str(T)+'\n')
# f.write("Surface Diffusion (cm^2/sec):\t" + str(D_s)+'\n')
# f.write("Bulk Diffusion (cm^2/sec):\t" + str(D_b)+'\n')
# f.write("Surface Fraction:\t" + str(Sfrac)+'\n')
# f.write("X-Gradient (G/cm):\t" + str(Gx)+'\n')
# f.write("Y-Gradient (G/cm):\t" + str(Gy)+'\n')
# f.write("Z-Gradient (G/cm):\t" + str(Gz)+'\n')
# f.write("Lx (cm):\t" + str(Lx)+'\n')
# f.write("Ly (cm):\t" + str(Ly)+'\n')
# f.write("Lz (cm):\t" + str(Lz)+'\n')
# f.write("B1 (G):\t" + str(B1)+'\n')
# f.write("kicking angle (degrees):\t" + str(np.degrees(theta))+'\n')
# f.write("kicking time (milli-sec):\t" + str(tp*1000)+'\n')
# f.write("f0 (Hz):\t"+str(f0)+'\n')
# f.write("Data Start:\r\n")
#
# for d in data: f.write(str(d)+'\n')
#
# f.close()
# print("Created "+fname)
# print("-----------------------------------------------------")
