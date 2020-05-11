import matplotlib.pyplot as plt
import math
import numpy as np
from WalkerClasses import *

np.random.seed(5)
random.seed(5)

D=1#diffusion constant, cm^2/sec
dt=0.1 #sec
N=1000

walkers=[]

for i in range(N): walkers.append(Walk(1.0,1.0))
for i in range(1000):
    for W in walkers: W.step(type="Gauss")


samples=[]
for W in walkers: samples.append(W.ptcl.x)
print("Sample Length")
print(len(samples))

#Histogram the data
hist,bin_edges=np.histogram(samples, bins=40)
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

    prob=stats.norm.cdf(bin_edges[i],loc=0,scale=sigma)-stats.norm.cdf(bin_edges[i-1],loc=0,scale=sigma)
    #print("prob: "+str(prob))
    Es.append(prob*N)
    #print("Expected: "+str(prob*N))
    stds.append(np.sqrt(prob*(1.-prob)*N))

    #print("std's above zero: "+str(Es[-1]/np.sqrt(stds[-1])))

x=np.linspace(-3*sigma,3*sigma,1000)
widths=(bin_centers[1]-bin_centers[0])*1
y=stats.norm.pdf(x,loc=0,scale=sigma)*N*widths
plt.bar(bin_centers,hist,label="data",ls='-',color='blue',fill=False,ec="blue",width=widths)
plt.bar(bin_centers,Es,yerr=stds,label="Expected",ls='-',color='green',fill=False,ec='green',width=widths)
plt.plot(x,y,linestyle="--",color="tab:green")
plt.legend()

obs,exp,err=binomialcull(hist,Es,stds,N)

DoF=len(obs)
chi2,pval=getPvalue(obs,exp,err,DoF)

print("Degrees of Freedom: "+str(DoF))
print("Reduced Chisquare: "+str(chi2/DoF))
print("p-value: "+str(pval))

plt.title("1000 Gaussian Random Walks after 1000 steps")
plt.xlabel("Position")
plt.ylabel("Counts")
plt.show()
