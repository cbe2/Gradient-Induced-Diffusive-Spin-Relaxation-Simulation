import matplotlib.pyplot as plt
import math
import numpy as np #np.interp
import time
import random
from scipy import stats

He3Mass=5.008234e-27  #in Kg
BoltzmannK= 1.380649e-23 #in J/K
RoomTemp=294.3 #Kelvin (70 degrees farenheit)



#This class serves as the particle
class Walker:
    def __init__(self):

        self.x=0 #meters
        self.t=0 #seconds
        self.steps=0 #integer
        self.theta=0


#Contains all relevant information for a single walk
class Walk:
    #L=length of box
    #D=diffusion constant
    #dt=time step
    def __init__(self,D,dt,L,Bfield,x0=0):

        self.ptcl=Walker()
        self.ptcl.x=x0
        self.dt=dt #seconds
        self.D=D #
        self.L=L
        self.Bfield=Bfield
        self.gamma=(2*np.pi)*(3240) #Hz/G#


    #steps the particle
    def step(self,type='const'):


        l=self.GetStep(type)

        #spin precesses at current position
        self.ptcl.theta+=self.gamma*self.Bfield(self.ptcl.x)*self.dt

        self.ptcl.x+=l
        self.ptcl.x=np.abs(self.ptcl.x) #This relfects about x=0
        self.ptcl.x=self.L-np.abs(self.ptcl.x-self.L) #reflects about x=L
        self.ptcl.t+=self.dt
        self.ptcl.steps+=1

        return

    #determines the mean free path of the step and time between steps
    def GetStep(self,type="const"):

        if type=="const":
            dt=self.dt
            D=self.D
            l=np.sqrt(2.*D*dt)

            if random.random()<= 0.5:l=-l

        if type=="Gauss":
            dt=self.dt
            D=self.D
            sigma=np.sqrt(2.*D*dt)

            l=sigma*np.random.normal(loc=0,scale=1.0)


        return l


#Given a set of spatial coordinates x and corresponding temperatures T, TempMap(x_i) returns
#the interpolated temperature at x_i. TempMap(x_i) is constant if x_i>x or x_i<x for all x.
class BMap:
    def __init__(self,x,B):

        self.x=x
        self.B=B

    def __call__(self, x):
        return np.interp(x,self.x,self.B)




#Removes data points that are binomial but not sufficiently guassian for chisqure test
#N is the number of samples that are binned to create the observed array
def binomialcull(Observed,Expected,Weights,N):

    CullObs=[]
    CullExp=[]
    CullWgths=[]
    for i in range(len(Observed)):
        delta0=Expected[i]/Weights[i]
        deltaN=(N-Expected[i])/Weights[i]

        if delta0>=3.0 and deltaN>=3. and Expected[i]>10:
            CullObs.append(Observed[i])
            CullExp.append(Expected[i])
            CullWgths.append(Weights[i])


    print(str(len(Observed)-len(CullObs))+" data points culled")

    return CullObs,CullExp,CullWgths


#returns chi square value
def getchisquare(Observed,Expected,Weights):
    sum=0
    for i in range(len(Observed)):
        sum+=np.power((Observed[i]-Expected[i])/Weights[i],2)

    return sum

#returns chisquare value and its pvalue
def getPvalue(Observed,Expected,Weights,dof):

    chi2=getchisquare(Observed,Expected,Weights)
    pval=stats.chi2.sf(chi2, df=dof)

    return chi2,pval

#returns p(x|x0) for 1D box
def GetDist(x,x0,L,D,t,N=100):
    sum=0

    for n in range(1,N+1):
        term=np.cos(np.pi*n*x/L)
        term*=np.cos(np.pi*n*x0/L)
        term*=np.exp(-np.power(np.pi*n/L,2)*D*t)
        sum+=term

    sum*=2./L
    sum+=1./L

    return sum


#returns probability of x in [xi,xf] for 1D box starting at x0
def GetProb(xi,xf,x0,L,D,t,N=1000):
    sum=0

    for n in range(1,N+1):
        term=np.sin(np.pi*n*xf/L)-np.sin(np.pi*n*xi/L)
        term*=np.cos(np.pi*n*x0/L)/n
        term*=np.exp(-np.power(np.pi*n/L,2)*D*t)
        sum+=term

    sum*=2./np.pi
    sum+=(xf-xi)/L

    return sum


def SaveData(x,y,fname):

    f=open("SavedData/"+fname,'w')

    for i in range(len(x)):
        f.write(str(x[i])+'\t'+str(y[i])+'\n')

    f.close()

    return


def LoadData(fname):
    x=[]
    y=[]

    f=open("SavedData/"+fname,'r')

    for line in f:
        line=line.split("\t")
        x.append(float(line[0]))
        y.append(float(line[1]))

    return np.asarray(x),np.asarray(y)
