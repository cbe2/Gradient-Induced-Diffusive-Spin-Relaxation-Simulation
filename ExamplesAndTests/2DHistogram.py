import matplotlib.pyplot as plt
import math
import numpy as np
import seaborn as sns
import pandas as pd
import sys
sys.path.append('../')
from WalkerClasses import *
import time
start_time = time.time()


#np.random.seed(5)
#random.seed(5)

paramsDict={
'D':1.0,
'r0':[0.1,0.2],
'dt': 0.5e-4,
'L': [1.,1.]
}

N=10000
steps=1000

walkers=[]

for i in range(N):
    #paramsDict['r0']=[np.random.uniform(0,L)]
    walkers.append(BoxWalk(paramsDict))
for i in range(steps):
    for W in walkers: W.step(lambda x:0)


Xsamples=[]
Ysamples=[]
for W in walkers:
    Xsamples.append(W.ptcl.r[0])
    Ysamples.append(W.ptcl.r[1])

#SaveData("")

counts, xedges, yedges, im=plt.hist2d(Xsamples,Ysamples,bins=10)
plt.colorbar()


total_time=walkers[1].ptcl.steps*walkers[1].dt

print("----counts---")
print(counts)

print("xedges")
print(xedges)
print("yedges")
print(yedges)

orderedCounts=[] #order list of counts for chi squared test
Es=[] #Expected values
stds=[] #expected uncertainty
#iterate through all bins
for i in range(1,len(xedges)):
    for j in range(1,len(yedges)):

        xf=xedges[i]; xi=xedges[i-1]; yf=yedges[j]; yi=yedges[j-1]
        #print("xEdge:")
        #print(xi,xf)
        #print("yEdge")
        #print(yi,yf)
        x0=paramsDict['r0'][0] ; y0=paramsDict['r0'][1]
        Lx=paramsDict['L'][0] ; Ly=paramsDict['L'][1] ;D=paramsDict['D']
        prob=GetProb(xi,xf,x0,Lx,D,total_time,N=1000)*GetProb(yi,yf,y0,Ly,D,total_time,N=1000)
        Es.append(prob*N)
        #print("measured")
        orderedCounts.append(counts[i-1][j-1])
        #print(counts[i-1][j-1])
        stds.append(np.sqrt(prob*(1.-prob)*N))



plt.title(r'$10^{4}$ Spins')
plt.xlabel("X position (cm)")
plt.ylabel("Y position (cm)")

obs,exp,err=binomialcull(orderedCounts,Es,stds,N)

DoF=len(obs)
chi2,pval=getPvalue(obs,exp,err,DoF)

print("Degrees of Freedom: "+str(DoF))
print("Reduced Chisquare: "+str(chi2/DoF))
print("p-value: "+str(pval))


TimeTaken=(time.time() - start_time)/60. #in minutes
print("{0:.2g}".format(TimeTaken)+" minutes")
print("{0:.2g}".format(TimeTaken/float(steps*N)) + " min per step*particle")

plt.show()
