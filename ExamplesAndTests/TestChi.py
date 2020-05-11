import matplotlib.pyplot as plt
import math
import numpy as np
from WalkerClasses import *
from scipy.stats import chisquare
from scipy import stats

def getchisquare(Observed,Expected,Weights):
    sum=0
    for i in range(len(Observed)):
        sum+=np.power((Observed[i]-Expected[i])/Weights[i],2)

    return sum


def getPvalue(Observed,Expected,Weights,dof):

    chi2=getchisquare(Observed,Expected,Weights)
    pval=stats.chi2.sf(chi2, df=dof)

    return chi2,pval

np.random.seed(5)

mu=0
sigma=1.
N=1000
samples = np.random.normal(mu, sigma, 1000)



hist,bin_edges=np.histogram(samples, bins=30)

bin_centers=0.5*(bin_edges[1:]+bin_edges[:-1])

Es=[] #expected values
stds=[]
sigma=sigma*1.0
for i in range(1,len(bin_edges)):

    prob=stats.norm.cdf(bin_edges[i],loc=0,scale=sigma)-stats.norm.cdf(bin_edges[i-1],loc=0,scale=sigma)
    #print("prob: "+str(prob))
    Es.append(prob*N)
    stds.append(np.sqrt(prob*(1.-prob)*N))

    #print("std's above zero: "+str(Es[-1]/np.sqrt(stds[-1])))


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

#     w.step()
#     # if w.ptcl.steps%10000==0:
#     #     print(w.ptcl.x)
#     #     print(len(t_inj))
#
# print(t_inj)
