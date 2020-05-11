import matplotlib.pyplot as plt
import math
import numpy as np
from WalkerClasses import *

#np.random.seed(5)
#random.seed(5)

N=int(1e3) #number of walkers
print(str(N)+" spins simulated")
D=1.0#diffusion constant, cm^2/sec
x0=0 #starting position

G=1e-3 #Gauss/cm
gamma=(2*np.pi)*(3240) #radian/Gauss
#the free T2 time
FreeT2=1./(np.power(gamma*G,2./3.)*np.power(4.*D,1./3.))
print("Free T2 time: "+str(FreeT2))
BDx=np.power((4.*D)/(gamma*G),1./3.)
print("Boundary Thickness: %3.3f" % (BDx))
L=1.0*BDx #Length of cell cm


SR=1000#100./(FreeT2)# 1000 Sample Rate
n=1#number of steps between samples
dt=(1./SR)*(1./n)
print("Time Step: "+str(dt))
print("Step Distance: "+str(np.sqrt(2.*D*dt)))
SS=0#int(FreeT2*1*SR) #sample size

f0=250#1./(4.*dt)#Hz
print("Frequency: "+str(f0))
B0=2.*np.pi*f0/gamma

x=np.linspace(-10*L,10*L,1000)
B=B0+G*x #Linear Gradient
B=BMap(x,B)

walkers=[]

for i in range(N): walkers.append(Walk(D,dt,L,B,np.random.uniform(0,L)))
for i in range(SS):
    for W in walkers: W.step(type="Gauss")


samples=[]
for W in walkers: samples.append(W.ptcl.x)
print("Sample Length")
print(len(samples))

#Histogram the data
hist,bin_edges=np.histogram(samples, bins=25)
#find bin centers
bin_centers=0.5*(bin_edges[1:]+bin_edges[:-1])

#find expected sigma and mu
total_time=walkers[1].ptcl.t
sigma=np.sqrt(2.*D*total_time)
mu=0
print("Expected sigma: "+str(sigma))

Es=[]
stds=[]
for i in range(1,len(bin_edges)):

    xf=bin_edges[i]; xi=bin_edges[i-1]
    prob=(xf-xi)/L#GetProb(xi,xf,x0,L,D,total_time,N=1000)
    Es.append(prob*N)
    stds.append(np.sqrt(prob*(1.-prob)*N))



x=np.linspace(0,L,100)
widths=(bin_centers[1]-bin_centers[0])*1
#y=GetDist(x,x0,L,D,total_time,N=100)*N*widths
y=(np.zeros(len(x))+1./L)*N*widths
plt.bar(bin_centers,hist,label="Simulation Data",ls='-',color='blue',fill=False,ec="blue",width=widths)
plt.bar(bin_centers,Es,yerr=stds,label="Expected",ls='-',color='green',fill=False,ec='green',width=widths)
plt.plot(x,y,linestyle="--",color="tab:green")
plt.legend()

obs,exp,err=binomialcull(hist,Es,stds,N)

DoF=len(obs)
chi2,pval=getPvalue(obs,exp,err,DoF)

print("Degrees of Freedom: "+str(DoF))
print("Reduced Chisquare: "+str(chi2/DoF))
print("p-value: "+str(pval))

plt.title("1000 Spins Initial Distribution")
plt.xlabel("Position")
plt.ylabel("Counts")
plt.show()
