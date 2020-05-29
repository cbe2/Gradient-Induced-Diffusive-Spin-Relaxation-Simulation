import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import sys
sys.path.append('../')
from WalkerClasses import *
from KuboRedfieldFunctions import *
import time
start_time = time.time()

#np.random.seed(5)
#random.seed(5)

paramsDict={
'D':1.0,
'r0':[0],
'dt': 0.5e-4,
'L': [1.]
}


N=100 #number of walkers

B0=0;G=1e-3 ;D=paramsDict['D'];gamma=(2*np.pi)*(3240) #Bfield params
FreeT2=1./(np.power(gamma*G,2./3.)*np.power(4.*D,1./3.)) #the free T2 time
print("Free T2 time: "+str(FreeT2))

paramsDict['dt']=0.01*FreeT2 #sets appropiate time scale
steps=int(FreeT2*7/paramsDict['dt']) #sets appropiate sample size

#-------Allows for N-dim interpolated Bfield, but is much slower

# #for the rectangular grid
# L=1
# x=np.linspace(-50*L,50*L,100)
# y=np.linspace(-50*L,50*L,100)
# z=np.linspace(-50*L,50*L,100)
#
# #N-Dim Bfield function
# def Bfield_func(x):
#     Gx=1e-3
#     B0=0
#     return B0+Gx*x
#
# data = Bfield_func(*np.meshgrid(x, indexing='ij', sparse=True))
# Bz_ = RegularGridInterpolator(points=[x], values=data) #Bz field funtion
# Bz=lambda x:Bz_(x)[0] #function that returns a float instead of a [float]

#----------

def Bz(r):
    Gx=1e-3
    return Gx*r[0]


#creates the random walkers
walkers=[]
for i in range(N): walkers.append(FreeWalk(paramsDict))

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

x0=paramsDict['r0'][0]
y=FreeDiff(times,x0,G,D,gamma*Bz([x0])) #Kubo-Redfield Analytical solution
sigma=np.sqrt(2.*np.power(gamma*G,2)*D*np.power(times,3)/3.)
yerr=(1.-np.exp(-sigma**2))/np.sqrt(2*N) # 1 sigma range


plt.plot(times,np.real(data),label="Simulated")
plt.plot(times,np.real(y),label="Kubo-Redfield" )
plt.fill_between(times,(y-yerr),(y+yerr),alpha=0.3,color='orange',label=str(r'$1\sigma$')+" confidence interval")
plt.title("Free Walk Test")
plt.ylabel(r'Re$[\langle e^{i\theta}\rangle]$')
plt.xlabel("Time (sec)")
plt.legend()

TimeElapsed=(time.time() - start_time)/60. #in minutes
print(str(TimeElapsed)+" minutes to simulate")

plt.show()
