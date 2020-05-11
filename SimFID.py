import matplotlib.pyplot as plt
import math
import numpy as np
from WalkerClasses import *
from BoundaryDiffFunctions import *
import time
start_time = time.time()

#np.random.seed(5)
#random.seed(5)

N=int(1e5) #number of walkers
print(str(N)+" spins simulated")
D=1.0#diffusion constant, cm^2/sec
x0=0 #starting position

#B0=0.5 #Gauss
G=1e-3 #Gauss/cm
gamma=(2*np.pi)*(3240) #radian/Gauss
#the free T2 time
FreeT2=1./(np.power(gamma*G,2./3.)*np.power(4.*D,1./3.))
print("Free T2 time: "+str(FreeT2))
BDx=np.power((4.*D)/(gamma*G),1./3.)
print("Boundary Thickness: %3.3f" % (BDx))
L=10.0*BDx #Length of cell cm


SR=1000#100./(FreeT2)# 1000 Sample Rate
n=1#number of steps between samples
dt=(1./SR)*(1./n)
print("Time Step: "+str(dt))
print("Step Distance: "+str(np.sqrt(2.*D*dt)))
SS=1000#int(FreeT2*14*SR) #sample size

f0=250#1./(4.*dt)#Hz
print("Frequency: "+str(f0))
B0=2.*np.pi*f0/gamma


x=np.linspace(-10*L,10*L,1000)
B=B0+G*x #Linear Gradient
B=BMap(x,B)
# plt.plot(x,B(x))
# plt.show()

walkers=[]
for i in range(N): walkers.append(Walk(D,dt,L,B,np.random.uniform(0,L)))


data=[]
times=[]
t=0

while len(data)<=SS:
    t+=dt
    ThetaSum=0
    for W in walkers:
        ThetaSum+=np.exp(1j*W.ptcl.theta)
        W.step(type="Gauss")

    if walkers[0].ptcl.steps%n==0:
        data.append(ThetaSum/N)
        times.append(t)


data=np.asarray(data)
times=np.asarray(times)
#print(times)

#times=np.arange(0,Ttotal,1./SR)
y=np.exp(-np.power(gamma*G,2)*np.power(L,4)*times/(120.*D))#BoundaryDiff(times,x0,G,D,f0)
sigma=np.sqrt(2.*np.power(gamma*G,2)*D*np.power(times,3)/3.)
#yerr=(1.-np.exp(-sigma**2))/np.sqrt(2*N)
plt.plot(times,np.abs(data),label="Simulated")


#plt.plot(times,np.abs(FreeDiff(times,x0,G,D,f0)),label="Analytical, No Wall" )
#plt.plot(times,np.abs(y),label="Analytical" )
#plt.fill_between(times,(y-yerr),(y+yerr),alpha=0.3,color='orange',label=str(r'$1\sigma$')+" confidence interval")
plt.title(str(r'$10^5$ spins, $G=10^{-9}$ Gauss/cm, D=1 cm$^2$/sec, $\Delta t=10^{-3}$ sec'))
plt.ylabel(r'$|\langle e^{i\theta}\rangle|$')
plt.xlabel("Time (sec)")
plt.legend()


TimeTaken=(time.time() - start_time)/60. #in minutes
print(str(TimeTaken)+" minutes")

SaveData(times,np.abs(data),"N1e5G3D1L10.txt")

plt.show()


path=""#"/Users/cameronerickson/Desktop/Academics/NPL/DRMRI/SignalSim/FakeSig/Signals/"
fname="L=10f0=250N=1e5.txt"

data=np.imag(data)

f=open(path+fname,'w')

f.write("Sample Rate:\t" + str(SR)+'\n')
f.write("Sample Size:\t" + str(SS)+'\n')
f.write("Read Time:\t" + str(float(SS)/SR)+'\n')
f.write("Diffusion (cm^2/sec):\t" + str(D)+'\n')
f.write("Gradient (G/cm):\t" + str(G)+'\n')
f.write("Data Start:\r\n")

for d in data: f.write(str(d)+'\n')

f.close()

print("Created "+fname)


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
