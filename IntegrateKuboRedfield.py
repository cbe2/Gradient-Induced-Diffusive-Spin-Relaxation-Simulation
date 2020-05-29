import matplotlib.pyplot as plt
import math
import cmath
import numpy as np
from scipy import special as sp
from KuboRedfieldFunctions import *
from WalkerClasses import *

dt=1e-3
SR=1./dt#50000.0 #Sample Rate
SS=1200.#Sample Size
Ttotal=SS/SR
print("Total time: "+str(Ttotal))

times=np.arange(0,Ttotal,1./SR)# Do not use linspace for this, it shifts the freqencies

gamma=(2*np.pi)*(3240) #Hz/G

f0=250.#16000.#in Hz Boundary Frequency
D=1 #cm^2/sec
G=Gx=1e-3 #G/cm
Gy=Gx #G/cm
BDx=np.power((4.*D)/(gamma*Gx),1./3.)
#BDy=np.power((4.*D)/(gamma*Gy),1./3.)
L=Lx=BDx*3.0
Ly=0

#f0=f0-(Lx*0.5)*(gamma*Gx)/(2.*np.pi)-(Ly*0.5)*(gamma*Gy)/(2.*np.pi)

print("Boundary Thickness: "+"{0:.3g}".format(BDx))

#creates probability density for integration
res=BDx*.001 #spatial integration step
xi=0
xf=Lx
yi=0
yf=Ly
pdfx_x=np.arange(xi,xf,res) #x-axis intial number density x-coordinates
pdfy_x=np.arange(yi,yf,res) #x-axis initial number deisntiy values
#equilibrium distribution:
pdfx_y=np.zeros(len(pdfx_x))+1./float(len(pdfx_x)) #creates an equilibrium distribution
#pdfy_y=np.zeros(len(pdfy_x))+1./float(len(pdfy_x))

#signal=IntegrateSigs1W(times,G,D,f0,pdf_x,pdf_y)
signal=IntegrateFullSol(times,G,D,L,f0,pdfx_x,pdfx_y)
#signal=IntegrateFullSol2D(times,Gx,Lx,Gy,Ly,D,f0,pdfx_x,pdfx_y,pdfy_x,pdfy_y)

plt.plot(times,np.abs(signal),label="L=3")
plt.legend()
plt.show()

data=np.imag(signal)

#-------Creating File------------------------
SaveData(times,np.abs(signal),"L=3KuboAbs.txt")

path="SimulationData/"
fname="L=3Kubo.txt"


f=open(path+fname,'w')

f.write("Sample Rate:\t" + str(SR)+'\n')
f.write("Sample Size:\t" + str(SS)+'\n')
f.write("Read Time:\t" + str(SS/SR)+'\n')
f.write("Diffusion (cm^2/sec):\t" + str(D)+'\n')
f.write("Gradient (G/cm):\t" + str(Gx)+'\n')
f.write("f0 (Hz):\t"+str(f0)+'\n')
f.write("Data Start:\r\n")

for d in data: f.write(str(d)+'\n')

f.close()

print("Created "+fname)
