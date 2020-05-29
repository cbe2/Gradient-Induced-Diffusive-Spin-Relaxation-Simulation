import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.stats import chisquare
from scipy import stats
import sys
sys.path.append('../')
from WalkerClasses import *


mu=0
sigma=1.
N=1000
samples = np.random.normal(mu, sigma, N)



hist,bin_edges=np.histogram(samples, bins=30)
bin_centers=0.5*(bin_edges[1:]+bin_edges[:-1])

Es=[] #expected values
stds=[] #standard deviations
for i in range(1,len(bin_edges)):

    prob=stats.norm.cdf(bin_edges[i],loc=0,scale=sigma)-stats.norm.cdf(bin_edges[i-1],loc=0,scale=sigma)
    Es.append(prob*N)
    stds.append(np.sqrt(prob*(1.-prob)*N))


x=np.linspace(-3*sigma,3*sigma,1000)
widths=(bin_centers[1]-bin_centers[0])*1
y=stats.norm.pdf(x,loc=0,scale=sigma)*N*widths
plt.bar(bin_centers,hist,label="data",ls='-',color='blue',fill=False,ec="blue",width=widths)
plt.bar(bin_centers,Es,yerr=stds,label="Expected",ls='-',color='green',fill=False,ec='green',width=widths)
plt.plot(x,y,linestyle="--",color="tab:green")

obs,exp,err=binomialcull(hist,Es,stds,N)

DoF=len(obs)
chi2,pval=getPvalue(obs,exp,err,DoF)

print("Degrees of Freedom: "+str(DoF))
print("Reduced Chisquare: "+str(chi2/DoF))
print("p-value: "+str(pval))


plt.title("1000 samples")
plt.xlabel("x")
plt.ylabel("counts")
plt.legend()
plt.show()
