import matplotlib.pyplot as plt
import math
import numpy as np
from WalkerClasses import *
from BoundaryDiffFunctions import *
import time
start_time = time.time()

np.random.seed(5)
random.seed(5)

N=10000 #number of walkers
D=1.0#diffusion constant, cm^2/sec
#dt=1e-4 #sec
L=1.0 #Length of cell cm
x0=0 #starting position
B0=5. #Gauss
G=1e-3 #Gauss/cm
gamma=(2*np.pi)*(3240)
#the free T2 time
FreeT2=1./(np.power(gamma*G,2./3.)*np.power(4.*D,1./3.))
print("Free T2 time: "+str(FreeT2))



SR=1000 #Sample Rate
n=1#number of steps between samples
dt=(1./SR)*(1./n)
SS=int(FreeT2*7*SR) #sample size


x=np.linspace(-100*L,100*L,1000)
B=G*x#B0+G*x #Linear Gradient
B=BMap(x,B)
# plt.plot(x,B(x))
# plt.show()

walkers=[]
for i in range(N): walkers.append(Walk(D,dt,L,B,x0))


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


# for i in range(SS):
#     times.append(t)
#     t+=dt
#     ThetaSum=0
#     for W in walkers:
#         ThetaSum+=np.exp(np.exp(1j*W.ptcl.theta))
#         W.step(type="Gauss")
#
#     data.append(ThetaSum/N)

data=np.asarray(data)
times=np.asarray(times)
#print(times)

#times=np.arange(0,Ttotal,1./SR)
y=np.real(FreeDiff(times,x0,G,D,gamma*B(x0)))
sigma=np.sqrt(2.*np.power(gamma*G,2)*D*np.power(times,3)/3.)
yerr=(1.-np.exp(-sigma**2))/np.sqrt(2*N)
plt.plot(times,np.real(data),label="Simulated")
plt.plot(times,y,label="Analytical" )
plt.fill_between(times,(y-yerr),(y+yerr),alpha=0.3,color='orange',label=str(r'$1\sigma$')+" confidence interval")
plt.title("Re[R], 100 spins, "+str(r'$G=1\times 10^{-3}$ Gauss/cm, D=1 cm$^2$/sec, $\Delta t=1\times 10^{-3}$ sec'))
plt.ylabel("Voltage")
plt.xlabel("Time (sec)")
plt.legend()

TimeElapsed=(time.time() - start_time)/60. #in minutes
print(str(TimeElapsed)+" minutes")

plt.show()
