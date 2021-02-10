import math
import numpy as np #np.interp
import time
import random
from scipy import stats

He3Mass=5.008234e-27  #in Kg
BoltzmannK= 1.380649e-23 #in J/K
RoomTemp=294.3 #Kelvin (70 degrees farenheit)

"""
Contains all the classes used for the random walk simulation. Times in seconds, lengths in cm
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

#Random walk in box of arbitary dimension and lengths specified by the tuple L
#Box is assumed to be centered on the origin.
class BoxWalk(Walk):
    def __init__(self,paramsDict):
        Walk.__init__(self,paramsDict)

        self.L=paramsDict['L'] #tuple of lengths of the box, all in cm. Length of the tuple determines the dimension of the box

        if len(self.L)!=len(self.ptcl.r):
            raise ValueError('Spatial dimension does not match box dimension!')

        if np.min(self.L)<= 3.*np.sqrt(self.dt*2.*self.D):
            raise ValueError('Step size within 3 sigma of box wall!')

        for i in range(len(self.L)):
            if np.abs(self.ptcl.r[i])> self.L[i]/2.:
                raise ValueError('particle not starting in box!')


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
            if self.ptcl.r[i]> 0.5*self.L[i]:
                self.ptcl.r[i]=self.L[i]-self.ptcl.r[i]
            if self.ptcl.r[i]< -0.5*self.L[i]:
                self.ptcl.r[i]=-self.L[i]-self.ptcl.r[i]
        return

#3D Random walk in a rectangular cell with feed through at the top along the z-axis
#Box+feed through is assumed be symmetric about z-axis. Feed through is assumed to be infinitely high
class CellWalk(Walk):
    def __init__(self,paramsDict):
        Walk.__init__(self,paramsDict)

        self.L=paramsDict['L'] #tuple of lengths of the box, all in cm. Length of the tuple determines the dimension of the box
        self.SL=paramsDict['SL'] #stem/straw lengths (x and y lengths only, z is infinite)

        if len(self.L)!=len(self.ptcl.r):
            raise ValueError('Spatial dimension does not match box dimension!')

        if len(self.ptcl.r)!=3:
            raise ValueError('Spatial dimension must be 3!')

        if np.min(self.L)<= 3.*np.sqrt(self.dt*2.*self.D):
            raise ValueError('Step size within 3 sigma of cell wall!')

        if np.min(self.SL)<= 3.*np.sqrt(self.dt*2.*self.D):
            raise ValueError('Step size within 3 sigma of stem wall!')

        for i in range(len(self.L)):
            if np.abs(self.ptcl.r[i])> self.L[i]/2. and (not self.inStem(self.ptcl.r)):
                raise ValueError('particle not starting in cell!')


    #steps the particle and updates particle preccesion according to Bfield
    def step(self,Bfield):

        #aquires next step vector
        l=self.getStep()


        #spin precesses at current position
        self.ptcl.theta+=self.gamma*Bfield(self.ptcl.r)*self.dt

        r0=self.ptcl.r.copy() #important to copy because of mutable aliasing
        self.ptcl.r+=l
        self.ptcl.steps+=int(1)
        self.reflect(r0)

        return

    #checks to see if walker is in stem, r is 3D vector
    def inStem(self,r):
        if np.abs(r[0])<0.5*self.SL[0] and np.abs(r[1])<0.5*self.SL[1]: return True
        else: return False

    #Reflects particle according specular reflection rule
    def reflect(self,r0):

        #reflect about x and y axis
        for i in range(2):
            if self.ptcl.r[i]> 0.5*self.L[i]:
                self.ptcl.r[i]=self.L[i]-self.ptcl.r[i]
            if self.ptcl.r[i]< -0.5*self.L[i]:
                self.ptcl.r[i]=-self.L[i]-self.ptcl.r[i]

        #reflect about bottom z
        if self.ptcl.r[2]< -0.5*self.L[2]:
            self.ptcl.r[2]=-self.L[2]-self.ptcl.r[2]

        #reflect about stem if both points are above cell
        if r0[2]> 0.5*self.L[2] and self.ptcl.r[2]> 0.5*self.L[2]:
            #reflect off stem walls if needed
            for i in range(2):
                if self.ptcl.r[i]> 0.5*self.SL[i]:
                    self.ptcl.r[i]=self.SL[i]-self.ptcl.r[i]
                if self.ptcl.r[i]< -0.5*self.SL[i]:
                    self.ptcl.r[i]=-self.SL[i]-self.ptcl.r[i]

        if (r0[2]< 0.5*self.L[2] and self.ptcl.r[2]> 0.5*self.L[2]) or (r0[2]> 0.5*self.L[2] and self.ptcl.r[2]< 0.5*self.L[2]):
            #if the particle remains in the stem, do nothing
            if self.inStem(r0) and self.inStem(self.ptcl.r):
                #print("Still in stem")
                return

            # if passing cell boundary out of stem, reflect off of cell
            elif (not self.inStem(r0)) and (not self.inStem(self.ptcl.r)):
                #print("not in stem")
                self.ptcl.r[2]=self.L[2]-self.ptcl.r[2]

            # if passing cell boundary but into the stem do nothing (to be fixed)
            elif (not self.inStem(r0)) and (self.inStem(self.ptcl.r)):
                #print("entering stem")
                return

            # if passing cell boundary and out the stem do nothing (to be fixed)
            elif (self.inStem(r0)) and (not self.inStem(self.ptcl.r)):
                #print("leaving stem")
                return

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


#returns p(x|x0) for 1D box from x=0 to x=L
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


#returns probability of x in [xi,xf] for 1D box (x=0 to x=L) starting at x0
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
