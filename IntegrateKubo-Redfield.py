import matplotlib.pyplot as plt
import math
import cmath
import numpy as np
from scipy import special as sp
from BoundaryDiffFunctions import *



SR=50000.#50000.0 #Sample Rate
SS=50000.#50000.0 #Sample Size
Ttotal=SS/SR
print("Total time: "+str(Ttotal))

times=np.arange(0,Ttotal,1./SR)# Do not use linspace for this, as it shifts the freqencies

gamma=(2*np.pi)*(3240) #Hz/G

f0= 16000.#16000.#in Hz Boundary Frequency
D=1 #cm^2/sec
G=0.1 #G/cm
BDx=np.power((4.*D)/(gamma*G),1./3.) #Depol has 1.25 cm, total length=7.5 to 8.3 cm (6 boundary thicknesses)
#BDx*3.0*0.5# in cm. Depol's thickness is 2.8 cm
L=BDx*6
x0=L/2.
#f0=f0-x0*(gamma*G)/(2.*np.pi) #adjusts frequency to be centered on x0
tau_d=np.power(BDx,2)/(4.*D)
T2=1./(BDx*gamma*G) #single boundary solution T2 time
print("Boundary Thickness = "+str(BDx))
print("T2 time = "+str(T2))

#parameters for integration
# res=BDx*.01
# xi=0
# xf=BDx*3.
# print(str(xf/BDx)+" Boundary Thicknesses Long")
# #pdf_x=np.arange(xi,xf,res)
#
# #equilibrium distribution:#np.zeros(len(pdf_x))+1./float(len(pdf_x))
# t=tau_d*9.0
# pdf_y=np.zeros(len(pdf_x))+1.#np.exp(-np.power(pdf_x,2.)/(4.*D*t))
# pdf_y=pdf_y/np.sum(pdf_y) #normalized


#OneWall=IntegrateSigs1W(times,G,D,f0,pdf_x,pdf_y)
#TwoWall=IntegrateSigs2W(times,G,D,f0,pdf_x,pdf_y)
#TwoWall=IntegrateSigsNW(times,G,D,f0,pdf_x,pdf_y)


#here are the four integratin
#Free=FreeDiff(times,x0,G,D,f0)
#Bdry=BoundaryDiff(times,x0,G,D,f0)
#NBdry=NearBoundaryDiff(times,x0,G,D,f0)
#TwoBdry=TwoBoundaryDiff(times,x0,G,D,L,f0)
FullBox=TwoBoundaryInt(times,G,D,L,f0,N=100)

#plt.plot(times,np.abs(Bdry),label="Free Diffusion")
plt.plot(times,np.abs(FullBox),label="Boundary Diffusion")
# plt.semilogx(times/T2,np.abs(Bdry),label="1 Wall")
# plt.semilogx(times/T2,np.abs(TwoBoundaryDiff(times,x0,G,D,BDx,f0)),label="2 Wall, L=1")
# plt.semilogx(times/T2,np.abs(TwoBoundaryDiff(times,x0,G,D,BDx*1.5,f0)),label="2 Wall, L=1.5")
# plt.semilogx(times/T2,np.abs(TwoBoundaryDiff(times,x0,G,D,BDx*2,f0)),label="2 Wall, L=2.0")
plt.title("|R(t)|")
plt.xlabel("Time (unitless)")
plt.ylabel("Magnitude (unitless)")
plt.legend()
#plt.figure()
# plt.plot(pdf_x,pdf_y)
# plt.title("spatial pdf")
# plt.legend()
plt.show()


data=np.imag(FullBox)
#data=np.append(np.flipud(data),data)
#plt.plot(np.append(times,times+Ttotal),np.abs(data),label="Flip")
#plt.show()
#print(data)
#Creating File
path=""#"/Users/cameronerickson/Desktop/Academics/NPL/DRMRI/SignalSim/FakeSig/Signals/"
fname="NewSol_L=6.txt"



f=open(path+fname,'w')

f.write("Sample Rate:\t" + str(SR)+'\n')
f.write("Sample Size:\t" + str(SS)+'\n')
f.write("Read Time:\t" + str(SS/SR)+'\n')
f.write("Diffusion (cm^2/sec):\t" + str(D)+'\n')
f.write("Gradient (G/cm):\t" + str(G)+'\n')
f.write("Data Start:\r\n")

for d in data: f.write(str(d)+'\n')

f.close()

print("Created "+fname)
