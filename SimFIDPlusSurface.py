import matplotlib.pyplot as plt
import math
import numpy as np
from WalkerClasses import *
from HePhysicalProperties import *
import time
start_time = time.time()

#All units cm,sec,Gauss

#np.random.seed(1)
#random.seed(1)

Ns=int(5e3) #number of surface walkers
Nb=int(5e3) #number of bulk walkers
N=Ns+Nb
steps=3300
print("\n{0:.1g}".format(N)+" spins simulating")

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


V=Lx*Ly*Lz
S=Lx*Ly

n4=2.18e22 # number density of 4He in cm^-3
N3=V*n4*x3

n_s=Getn_s(N3,S,V,T) #gets number of layers
print("Surface Layers: "+"{0:.3g}".format(n_s))

nl=6.4e14 #one layer surface number density (cm^-2)
print("Surface Atoms: "+"{0:.3g}".format(nl*n_s*S))

Sfrac = 1#(n_s*nl)*S/N3 #Fraction on the surface N_s/(N_s+N_b)
print("Surface Fraction: "+"{0:.3g}".format(Sfrac))
print("------------------------------------")

#B is for bulk
paramsDictB={
'D':D_b,
'r0':[Lx,Ly,Lz],
'dt': 2e-3,
'L': [Lx,Ly,Lz]
}

#S is for surface
paramsDictS={
'D':D_s,
'r0':[Lx,Ly],
'dt': 2e-3,
'L': [Lx,Ly]
}


f0=1./(4*paramsDictB['dt']) #choosen precession frequency
print("Frequency: "+str(f0))
B0=2.*np.pi*f0/gamma #holding field vaule to match f0


#Defines magnetic field
def _Bz(r,Gs,B0):
    B=B0
    for i in range(len(r)): B+=Gs[i]*r[i]
    return B

Bz=lambda x: _Bz(x,[Gx,Gy,Gz],B0) #So Bz only takes one input


#creates distribution of bulk walkers
walkers=[]
for i in range(Nb):
    paramsDictB['r0']=[np.random.uniform(0,Lx),np.random.uniform(0,Ly),np.random.uniform(0,Lz)]
    walkers.append(BoxWalk(paramsDictB))

#creates distribution of surface walkers
for i in range(Nb):
    paramsDictS['r0']=[np.random.uniform(0,Lx),np.random.uniform(0,Ly)]
    walkers.append(BoxWalk(paramsDictS))

data=[]
times=[]
#performs the simulation
for i in range(steps):
    thetaSumB=0
    thetaSumS=0
    for W in walkers:
        if len(W.ptcl.r)==3: thetaSumB+=np.exp(1j*W.ptcl.theta)
        if len(W.ptcl.r)==2: thetaSumS+=np.exp(1j*W.ptcl.theta)
        W.step(Bz)


    #scale the surface signal to the right proportion,
    #keeps total normalization to 1
    thetaSumS=thetaSumS*Sfrac/Ns
    thetaSumB=thetaSumB/Nb#*(1.-Sfrac)/Nb

    data.append(thetaSumS+thetaSumB)
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
SaveData(times,np.abs(data),"SimulationData/SurfaceTestSf=1e-1Abs.txt")

path="SimulationData/"
fname="Sf=1.txt"

data=np.imag(data)

f=open(path+fname,'w')

f.write("Spins Simulated:\t"+str(int(N))+'\n')
f.write("Sample Rate:\t" + str(1./paramsDictB['dt'])+'\n')
f.write("Sample Size:\t" + str(steps)+'\n')
f.write("Read Time:\t" + str(float(steps)*paramsDictB['dt'])+'\n')
f.write("Surface Diffusion (cm^2/sec):\t" + str(paramsDictS['D'])+'\n')
f.write("Bulk Diffusion (cm^2/sec):\t" + str(paramsDictB['D'])+'\n')
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
