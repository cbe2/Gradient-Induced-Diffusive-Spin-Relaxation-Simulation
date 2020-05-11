import matplotlib.pyplot as plt
import math
import numpy as np
from WalkerClasses import *
from BoundaryDiffFunctions import *
import time
start_time = time.time()

#np.random.seed(5)
#random.seed(5)

N=1000 #number of walkers
D=1.0#diffusion constant, cm^2/sec
#dt=1e-4 #sec
x0=0 #starting position
B0=5. #Gauss
G=1e-3 #Gauss/cm
gamma=(2*np.pi)*(3240)
#the free T2 time
FreeT2=1./(np.power(gamma*G,2./3.)*np.power(4.*D,1./3.))
print("Free T2 time: "+str(FreeT2))
BDx=np.power((4.*D)/(gamma*G),1./3.)
print("Boundary Thickness: %3.3f" % (BDx))
L=1.0*BDx #Length of cell cm



SR=1000 #Sample Rate
n=1#number of steps between samples
dt=(1./SR)*(1./n)
SS=int(FreeT2*2*SR) #sample size


x=np.linspace(-100*L,100*L,1000)
B=B0+G*x #Linear Gradient
B=BMap(x,B)
# plt.plot(x,B(x))
# plt.show()

walkers=[]
for i in range(N): walkers.append(Walk(D,dt,L,B,np.random.uniform(0,L)))

for i in range(SS):
    for W in walkers:
        W.step(type="Gauss")

samples=[]
for W in walkers: samples.append(W.ptcl.theta-gamma*B(x0)*W.ptcl.t)


#Histogram the data
hist,bin_edges=np.histogram(samples, bins=40)
#find bin centers
bin_centers=0.5*(bin_edges[1:]+bin_edges[:-1])


#find expected sigma and mu
t=walkers[1].ptcl.t
print("Simulation time elapsed: "+str(t))
sigma=np.sqrt(2.*np.power(gamma*G,2)*D*np.power(t,3)/3.)
mu=0
print("Expected sigma: "+str(sigma))

Es=[]
stds=[]
for i in range(1,len(bin_edges)):

    prob=stats.norm.cdf(bin_edges[i],loc=0,scale=sigma)-stats.norm.cdf(bin_edges[i-1],loc=0,scale=sigma)
    Es.append(prob*N)
    stds.append(np.sqrt(prob*(1.-prob)*N))


x=np.linspace(-3*sigma,3*sigma,1000)
widths=(bin_centers[1]-bin_centers[0])*1
y=stats.norm.pdf(x,loc=0,scale=sigma)*N*widths
plt.bar(bin_centers,hist,label="t=0.168 sec",ls='-',color='blue',fill=False,ec="blue",width=widths)
plt.bar(bin_centers,Es,yerr=stds,label="Expected",ls='-',color='green',fill=False,ec='green',width=widths)
plt.plot(x,y,linestyle="--",color="tab:green")
plt.legend()

obs,exp,err=binomialcull(hist,Es,stds,N)

DoF=len(obs)
chi2,pval=getPvalue(obs,exp,err,DoF)

print("Degrees of Freedom: "+str(DoF))
print("Reduced Chisquare: "+str(chi2/DoF))
print("p-value: "+str(pval))

plt.title("10000 spins, "+str(r'$G=1\times 10^{-3}$ Gauss/cm, D=1 cm$^2$/sec, $\Delta t=10^{-3}$ sec'))
plt.xlabel(r'$\theta-\gamma B_0 t$ (radians)')
plt.ylabel("Counts")

TimeTaken=(time.time() - start_time)/60. #in minutes
print(str(TimeTaken)+" minutes")

plt.show()
