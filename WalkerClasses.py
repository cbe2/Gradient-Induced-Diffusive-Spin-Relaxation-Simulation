import math
import numpy as np #np.interp
import time
import random
from scipy import stats

He3Mass=5.008234e-27  #in Kg
BoltzmannK= 1.380649e-23 #in J/K
RoomTemp=294.3 #Kelvin (70 degrees farenheit)

"""
Contains all the classes used for the random walk simulation
"""

#This class serves as the particle
class Walker:
    def __init__(self,dim):

        self.r=np.zeros(dim) #position vector: [x,y,z] in cm
        self.steps=int(0) #integer number of steps the particle has taken
        self.theta=0 #angular position (radians)


#Contains all relevant information for a general single walk with no geometery specified
class Walk:

    def __init__(self,paramsDict):

        self.ptcl=Walker(len(paramsDict['r0'])) #number of spatial dimension
        self.ptcl.r=np.asarray(paramsDict['r0'],dtype='float64') #starting position (cm)
        self.dt=paramsDict['dt'] #simulation time step (sec)
        self.D=paramsDict['D'] #diffusion constant (cm^2/sec)
        self.gamma=(2*np.pi)*(3240) #gyromagentic ratio (Hz/Gauss)

    #determines the next step in the random walk
    def getStep(self):
        sigma=np.sqrt(2.*self.D*self.dt)
        return sigma*np.random.normal(loc=0,scale=1.0,size=len(self.ptcl.r))


#Random walk with no boundaries
class FreeWalk(Walk):
    def __init__(self,paramsDict):
        Walk.__init__(self,paramsDict)

    #steps the particle and updates particle preccesion according to Bfield
    def step(self,Bfield):

        #aquires next step vector
        l=self.getStep()

        #spin precesses at current position
        self.ptcl.theta+=self.gamma*Bfield(self.ptcl.r)*self.dt

        #update's partilce position and time
        self.ptcl.r+=l
        self.ptcl.steps+=int(1)

        return

#Random walk with infinite walls along x,y, and z axes
class SemiFreeWalk(Walk):
    def __init__(self,paramsDict):
        Walk.__init__(self,paramsDict)

    #steps the particle and updates particle preccesion according to Bfield
    def step(self,Bfield):

        #aquires next step vector
        l=self.getStep()

        #spin precesses at current position
        self.ptcl.theta+=self.gamma*Bfield(self.ptcl.r)*self.dt

        #update's partilce position and time
        self.ptcl.r+=l

        #reflects about x=0
        self.ptcl.r=np.abs(self.ptcl.r)

        self.ptcl.steps+=int(1)

        return

#Random walk in box with lengths specified by the tuple L
class BoxWalk(Walk):
    def __init__(self,paramsDict):
        Walk.__init__(self,paramsDict)

        self.L=paramsDict['L'] #tuple of lengths of the box, all in cm

        if len(self.L)!=len(self.ptcl.r):
            raise ValueError('Spatial dimension does not match box dimension!')

        if np.min(self.L)<= 3.*np.sqrt(self.dt*2.*self.D):
            raise ValueError('Step size within 3 sigma of box wall!')



    #steps the particle and updates particle preccesion according to Bfield
    def step(self,Bfield):

        #aquires next step vector
        l=self.getStep()

        #spin precesses at current position
        self.ptcl.theta+=self.gamma*Bfield(self.ptcl.r)*self.dt

        self.ptcl.r+=l
        self.ptcl.steps+=int(1)
        self.reflect()

        return

    #Reflects particle according specular reflection rule
    def reflect(self):

        for i in range(len(self.L)):
            if self.ptcl.r[i]<self.L[i]:
                self.ptcl.r[i]=np.abs(self.ptcl.r[i]) #relfects about x=0
            if self.ptcl.r[i]>self.L[i]:
                self.ptcl.r[i]=self.L[i]-np.abs(self.ptcl.r[i]-self.L[i]) #reflects about x=L

        return



#Removes data points that are binomial but not sufficiently guassian for chisqure test
#N is the number of samples that are binned to create the observed array
def binomialcull(Observed,Expected,Weights,N):

    CullObs=[]
    CullExp=[]
    CullWgths=[] #standard devaition
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

    f=open(fname,'w')

    for i in range(len(x)):
        f.write(str(x[i])+'\t'+str(y[i])+'\n')

    f.close()

    return


def LoadData(fname):
    x=[]
    y=[]

    f=open(fname,'r')

    for line in f:
        line=line.split("\t")
        x.append(float(line[0]))
        y.append(float(line[1]))

    return np.asarray(x),np.asarray(y)

#returns the boundary thickness
def bThickness(D,G,gamma):

    if G==0: return np.inf
    else: return np.power((4.*D)/(gamma*np.abs(G)),1./3.)


#returns the resonance thickness (B1 in Gauss, G in Gauss/cm, theta in radians)
def rThickness(B1,G,theta):

    if G==0 or theta==0: return np.inf
    else: return B1/(G*theta)

#returns the diffusion time (G in Gauss/cm, D in cm^2/sec, gamma in Hz/G)
def Dtime(G,gamma,D):

    if G==0 or D==0: return np.inf
    else: return np.power(3./(D*(gamma*G)**2),1./3.)


#returns the spatial dephasing time (G in Gauss/cm, D in cm^2/sec, gamma in Hz/G)
def Stime(G,gamma,L):

    if G==0 or L==0: return np.inf
    else: return 2.*np.pi/(gamma*G*L)
