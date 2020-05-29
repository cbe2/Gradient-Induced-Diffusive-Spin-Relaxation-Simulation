import matplotlib.pyplot as plt
import math
import cmath
import numpy as np
from scipy import special as sp

pi2=np.pi*2.0

u=1.660539e-27#atomic mass unit in kg
eb=2.28 #binding energy in kelvin
m3=3.0160293*u # mass of 3He in Kg

m_s=1.53*m3 #effective mass of surface
m_b=2.34*m3 #effective mass of bulk
hbar=1.0545718e-34 #J*s planks constant
k_B=1.380649e-23 #J/K boltzmanns constant


#t in seconds, x0 in cm, G in G/cm, D in cm^2/sec, f0 is the frequency at the boundary in Hz
def FreeDiff(t,x0,G,D,f0):

    gamma=(2*np.pi)*(3240) #Hz/G
    #tauD=(x0**2)/(4.*D)
    #tau=t/tauD

    #alpha=gamma*G*x0*tauD #dimensionless constant measuring boundary effects

    #phi=alpha*tau
    phi=gamma*G*x0*t

    #A=-(alpha**2)*(1./12.)*tau**3
    A=-np.power(gamma*G,2)*D*np.power(t,3)/3.

    theta=2.*np.pi*f0*t+phi

    R=np.exp(1j*theta+A)

    return R



#t in seconds, x0 in cm, G in G/cm, D in cm^2/sec, f0 is the frequency at the boundary in Hz
def BoundaryDiff(t,x0,G,D,f0):

    gamma=(2*np.pi)*(3240) #Hz/G
    #tauD=(x0**2)/(4.*D)
    #tau=t/tauD

    #alpha=gamma*G*x0*tauD #dimensionless constant measuring boundary effects

    #phi=alpha*(2./(3.*np.sqrt(np.pi)))*tau**(3./2.)

    phi=(4./3.)*G*gamma*np.power(D/np.pi,.5)*np.power(t,3./2.)

    #A=-(alpha**2)*(3./32.-2./(9.*np.pi))*tau**3
    A=-np.power(gamma*G,2)*D*np.power(t,3)*(3./8.-8./(9.*np.pi))

    theta=2*np.pi*f0*t+phi

    R=np.exp(1j*theta+A)

    return R


#t in seconds, x0 in cm, G in G/cm, D in cm^2/sec, f0 is the frequency at the boundary in Hz
#x0=starting distance from the wall, G is the gradient which points away from the wall if positive
def NearBoundaryDiff(t,x0,G,D,f0):


    if x0==0: return BoundaryDiff(t,x0,G,D,f0)

    gamma=(2*np.pi)*(3240) #Hz/G
    tauD=(x0**2)/(4.*D)
    tau=t/tauD

    alpha=gamma*G*x0*tauD #dimensionless constant measuring boundary effects

    phiTerm1=(2./(3.*np.sqrt(np.pi)))*(1.+tau)*np.sqrt(tau)*np.exp(-1./tau)
    phiTerm2=(2./3.+tau)*sp.erf(1./np.sqrt(tau))
    phi=alpha*(phiTerm1+phiTerm2-2./3.) #phase shift term

    Aterm1=(tau**3)/12.+(tau**2)/2.
    A288term1=(3.*tau**3+18.*tau**2+108.*tau+40.)*sp.erfc(1./np.sqrt(tau))
    A288term2=(6.*tau**2-88.*tau-40.)*np.sqrt(tau/np.pi)*np.exp(-1./tau)
    Aterm2=(A288term1+A288term2)/288.
    Aterm3=-0.5*(phi/alpha)**2
    A=-(alpha**2)*(Aterm1+Aterm2+Aterm3) #Amplitude term

    theta=2*np.pi*f0*t+phi

    R=np.exp(1j*theta+A)

    return R

#only the theta component of 2 boundary diff function
def _Theta(t,x0,G,D,L,f0,N=100):
    gamma=(2*np.pi)*(3240) #Hz/G

    #Determination of \phi and \Theta
    sum=0
    for i in range(N):
        n=2*i+1

        term=np.cos(np.pi*n*x0/L)/np.power(n,4)#np.cos(np.pi*n*x0/L)
        term=term*(1.-np.exp(-np.power(n*np.pi/L,2)*D*t))

        sum+=term


    sum=sum*(4.*L**3)/(D*np.pi**4)

    #phi=(integral of expectation value)*gamma*G
    phi=(L/2.)*t-sum
    phi=phi*gamma*G

    #theta=integral of 1st cumulant
    theta=2*np.pi*f0*t+phi

    return theta

# the A component of the 2 boundary diff function
def _A(t,x0,G,D,L,f0,N=100):

    gamma=(2*np.pi)*(3240) #Hz/G

    #Determination of A
    sum1=0
    sum2=0
    for i in range(N):
        n=2*i+1
        sum1+=np.exp(-np.power(n*np.pi/L,2)*D*t)/np.power(n,8)
        sum2+=(1.-np.exp(-np.power(n*np.pi/L,2)*D*t))*np.cos(np.pi*n*x0/L)/np.power(n,4)


    doubleSum=0

    #double sum ensured to have N evaluations
    Nd=int(np.sqrt(N))
    for i in range(1,Nd+1):
        for j in range(Nd):
            n=2*i
            m=2*j+1
            denom=(m**2)*((m**2-n**2)**3)
            numer=n**2+m**2

            term1=(1.-np.exp(-np.power(n*np.pi/L,2)*D*t))/np.power(n,2)

            term2=(1.-np.exp(-np.power(m*np.pi/L,2)*D*t))/np.power(m,2)

            doubleSum+= (float(numer)/float(denom))*(term1-term2)*np.cos(np.pi*n*x0/L)


    A=np.power(L,4)*t/(120.*D)
    A+=-(17.*L**6)/(20160.*(D**2))
    A+=sum1*(8.*L**6)/((np.pi**8)*(D**2))
    A+=-(sum2**2)*(8.*L**6)/((np.pi**8)*(D**2))
    A+=doubleSum*(16.*L**6)/((np.pi**8)*(D**2))

    A=A*np.power(gamma*G,2)

    return A


#Boundaries at x=0 and x=L. t in seconds, x0 in cm, G in G/cm, D in cm^2/sec, f0 is the frequency at the x=0 boundary in Hz
#x0=starting position, G is the gradient which points in the x+ direction if positive
#Point by point averaging
def TwoBoundaryDiff(t,x0,G,D,L,f0,N=100):

    gamma=(2*np.pi)*(3240) #Hz/G

    #Determination of \phi and \Theta
    sum=0
    for i in range(N):
        n=2*i+1

        term=np.cos(np.pi*n*x0/L)/np.power(n,4)#np.cos(np.pi*n*x0/L)
        term=term*(1.-np.exp(-np.power(n*np.pi/L,2)*D*t))

        sum+=term


    sum=sum*(4.*L**3)/(D*np.pi**4)

    #phi=(integral of expectation value)*gamma*G
    phi=(L/2.)*t-sum
    phi=phi*gamma*G

    #theta=integral of 1st cumulant
    theta=2*np.pi*f0*t+phi


    #Determination of A
    sum1=0
    sum2=0
    for i in range(N):
        n=2*i+1

        sum1+=np.exp(-np.power(n*np.pi/L,2)*D*t)/np.power(n,8)

        #np.cos(np.pi*n*x0/L)

        sum2+=(1.-np.exp(-np.power(n*np.pi/L,2)*D*t))*np.cos(np.pi*n*x0/L)/np.power(n,4)


    doubleSum=0

    #double sum ensured to have N evaluations
    Nd=int(np.sqrt(N))
    for i in range(1,Nd+1):
        for j in range(Nd):
            n=2*i
            m=2*j+1
            denom=(m**2)*((m**2-n**2)**3)
            numer=n**2+m**2

            term1=(1.-np.exp(-np.power(n*np.pi/L,2)*D*t))/np.power(n,2)

            term2=(1.-np.exp(-np.power(m*np.pi/L,2)*D*t))/np.power(m,2)

            doubleSum+= (float(numer)/float(denom))*(term1-term2)*np.cos(np.pi*n*x0/L)


    A=np.power(L,4)*t/(120.*D)
    A+=-(17.*L**6)/(20160.*(D**2))
    A+=sum1*(8.*L**6)/((np.pi**8)*(D**2))
    A+=-(sum2**2)*(8.*L**6)/((np.pi**8)*(D**2))
    A+=doubleSum*(16.*L**6)/((np.pi**8)*(D**2))

    A=A*np.power(gamma*G,2)




    R=np.exp(1j*theta-A)

    return R


#Boundaries at x=0 and x=L. t in seconds, x0 in cm, G in G/cm, D in cm^2/sec, f0 is the frequency at the x=0 boundary in Hz
#x0=starting position, G is the gradient which points in the x+ direction if positive
#This is "full cell" averaging
def TwoBoundaryInt(t,G,D,L,f0,N=100):

    gamma=(2*np.pi)*(3240) #Hz/G

    #phi=(integral of expectation value)*gamma*G
    phi=(L/2.)*t
    phi=phi*gamma*G

    #theta=integral of 1st cumulant
    theta=2*np.pi*f0*t+phi


    #Determination of A
    sum1=0
    for i in range(N):
        n=2*i+1
        sum1+=np.exp(-np.power(n*np.pi/L,2)*D*t)/np.power(n,8)


    A=np.power(L,4)*t/(120.*D)
    A+=-(17.*L**6)/(20160.*(D**2))
    A+=sum1*(8.*L**6)/((np.pi**8)*(D**2))

    A=A*np.power(gamma*G,2)


    R=np.exp(1j*theta-A)

    return R



#integrates the signal for one wall, pdf is discrete and normed to 1
def IntegrateSigs1W(times,G,D,f0,pdf_x,pdf_y):

    gamma=(2*np.pi)*(3240) #Hz/G

    sig=np.zeros(len(times),dtype='complex128')
    for i in range(len(pdf_x)): sig+=NearBoundaryDiff(times,pdf_x[i],G,D,f0)*pdf_y[i]

    return sig

#integrates the full solution for one wall, pdf is discrete and normed to 1
def IntegrateFullSol(times,G,D,L,f0,pdf_x,pdf_y):

    gamma=(2*np.pi)*(3240) #Hz/G

    sig=np.zeros(len(times),dtype='complex128')
    for i in range(len(pdf_x)):
        print(pdf_x[i])
        sig+=TwoBoundaryDiff(times,pdf_x[i],G,D,L,f0)*pdf_y[i]

    return sig

def IntegrateFullSol2D(times,Gx,Lx,Gy,Ly,D,f0,pdfx_x,pdfx_y,pdfy_x,pdfy_y):

    gamma=(2*np.pi)*(3240) #Hz/G

    sig=np.zeros(len(times),dtype='complex128')
    for i in range(len(pdfx_x)):
        print(pdfx_x[i])
        for j in range(len(pdfy_y)):
            A=_A(times,pdfx_x[i],Gx,D,Lx,f0,N=100)+_A(times,pdfy_x[j],Gy,D,Ly,f0,N=100)
            Theta=_Theta(times,pdfx_x[i],Gx,D,Lx,f0,N=100)+_Theta(times,pdfy_x[j],Gy,D,Ly,f0,N=100)-2*np.pi*f0*times

            R=np.exp(1j*Theta-A)

            sig+=R*pdfx_y[i]*pdfy_y[j]

    return sig


def IntegrateSigs2W(times,G,D,f0,pdf_x,pdf_y):

    gamma=(2.*np.pi)*(3240) #Hz/G

    FWG=-G #far wall gradient
    FWf0=f0+gamma*G*pdf_x[-1]/(2.*np.pi) #far wall frequency
    FWD=D #far wall diffusion constant

    sig=np.zeros(len(times),dtype='complex128')
    for i in range(len(pdf_x)):
        #If closest to x=0 wall
        if pdf_x[i]<= pdf_x[-1]/2.:
            sig+=NearBoundaryDiff(times,pdf_x[i],G,D,f0)*pdf_y[i]
            print("First wall x=" +str(pdf_x[i]))
        #If closest to x=xf wall
        else:
            FWx=pdf_x[-1]-pdf_x[i]
            print("Second wall x=" +str(pdf_x[i]) + "But use "+str(FWx))
            sig+=NearBoundaryDiff(times,FWx,FWG,FWD,FWf0)*pdf_y[i]

    return sig

def IntegrateSigsNW(times,G,D,f0,pdf_x,pdf_y):

    gamma=(2.*np.pi)*(3240) #Hz/G

    FWG=-G #far wall gradient
    FWf0=f0+gamma*G*pdf_x[-1]/(2.*np.pi) #far wall frequency
    FWD=D #far wall diffusion constant

    sig=np.zeros(len(times),dtype='complex128')
    for i in range(len(pdf_x)):
        #If closest to x=0 wall
        if pdf_x[i]<= pdf_x[-1]/2.:
            sig+=FreeDiff(times,pdf_x[i],G,D,f0)*pdf_y[i]
            print("First wall x=" +str(pdf_x[i]))
        #If closest to x=xf wall
        else:
            FWx=pdf_x[-1]-pdf_x[i]
            print("Second wall x=" +str(pdf_x[i]) + "But use "+str(FWx))
            sig+=FreeDiff(times,FWx,FWG,FWD,FWf0)*pdf_y[i]

    return sig

#Returns surface number density (in layers) as long as N>>A*n_s
#n_s is in units of cm^-2 and T is in kevlin
def Getn_s(n_b,T,N=10):

    s=0.5 #spin of system
    e_b=2.28 #Surface binding energy in Kelvin
    nl=6.4e14 #one layer surface number density (cm^-2)

    n_b=n_b*1e6 #converts to m^-2

    #surface quantum number density (m^-2)
    n_Qs=m_s*k_B*T/(2.*np.pi*hbar**2)

    #bulk quantum number density (m^-3)
    n_Qb=(m_b*k_B*T/(2.*np.pi*hbar**2))**(3./2.)

    x=n_b/((1.+2.*s)*(n_Qb))
    n_s=(1.+2.*s)*n_Qs*np.log(x*np.exp(e_b/T)+1.)

    #3He 3He interaction term
    V_0=0.23e-18 #in K*m^2


    #fixed point iteration
    for i in range(N):
        n_s=(1.+2.*s)*n_Qs*np.log(x*np.exp((e_b-0.5*V_0*n_s)/T)+1.)

    return (n_s*1e-4)/nl #converts to cm^-2 then to layers


#This function is essentially the same as Getn_s, but with different inputs
#N3= total number of atoms, S=surface area (cm^2), V=volume (cm^3)
#n=number of iterations
def Getn_s2(N3,S,V,T,n):

    s=0.5 #spin of system
    e_b=2.28 #Surface binding energy in Kelvin
    nl=6.4e14 #one layer surface number density (cm^-2)

    #n_b=n_b*1e6 #converts to m^-2

    #Gets n_s for the case when N_s<<N_b
    n_s=Getn_s(N3/V,T,10)

    #Worst case
    Sfrac=np.max(n_s*nl*S/N3)
    print("Maximum fraction on surface (this should be much smaller than 1): "+str(Sfrac))

    if Sfrac>0.1: print("Warning! n_s may not be valid")

    # print(n_s[0])
    # # #Case when N_s~N_b
    # for i in range(n):
    #     n_s=Getn_s((N3/V)-(S/V)*n_s*nl,T,10)
    #     print(n_s[0])

    return n_s #in  layers



#Gets fraction of 3He on surface assuming 2D high temp limit
#S in cm, V in cm, T in kelvin
def GetSfrac(S,V,T):

    V=V*1e-6#cm^3 -> m^3
    S=S*1e-4#cm^2 -> m^2

    alpha=np.power(S/V,2.)*np.power(m_s,2)*np.power(m_b,-3)*(2.*np.pi*np.power(hbar,2.)/k_B)
    #print(alpha)

    Tb=2.28 #binding energy in Kelvin

    frac=1./(1.+np.sqrt(T/alpha)*np.exp(-Tb/T))

    return frac

#n in layers, T in kelvin
# def GetSurfaceD(n,T):
#
#     Dref=1e-2 #cm^2/sec from paper. Taken at 0.1 layers and T=0.030 K
#
#     D=Dref*((n/0.1)**2)((0.030/T)**2)*(np.log(0.1/))
#
#     return
