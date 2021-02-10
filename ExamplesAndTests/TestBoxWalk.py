import matplotlib.pyplot as plt
import math
import numpy as np
import sys
sys.path.append('../')
from WalkerClasses import *

np.random.seed(5)
random.seed(5)

paramsDict={
'D':1.0,
'r0':[0.4,0.2,0.1],
'dt': 0.5e-4,
'L': [1.,1.,.5]
}

N=1000
steps=100

walkers=[]

for i in range(N):
    #paramsDict['r0']=[np.random.uniform(0,L)]
    walkers.append(BoxWalk(paramsDict))
for i in range(steps):
    for W in walkers: W.step(lambda x:0)


samples=[]
for W in walkers: samples.append(W.ptcl.r[2])


#Histogram the data
hist,bin_edges=np.histogram(samples, bins=25)
#find bin centers
bin_centers=0.5*(bin_edges[1:]+bin_edges[:-1])

#find expected sigma and mu
total_time=walkers[1].ptcl.steps*walkers[1].dt


Es=[]
stds=[]
for i in range(1,len(bin_edges)):

    xf=bin_edges[i]; xi=bin_edges[i-1]
    x0=paramsDict['r0'][2] ; L=paramsDict['L'][2] ; D=paramsDict['D']
    prob=GetProb(xi+L/2.,xf+L/2.,x0+L/2.,L,D,total_time,N=1000)# +L/2 due to coordinate shift
    Es.append(prob*N)
    stds.append(np.sqrt(prob*(1.-prob)*N))



x=np.linspace(-L/2.,L/2.,100)
widths=(bin_centers[1]-bin_centers[0])*1
y=GetDist(x+L/2.,x0+L/2.,L,D,total_time,N=100)*N*widths
#y=(np.zeros(len(x))+1./L)*N*widths
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
