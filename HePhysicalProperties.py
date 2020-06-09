import matplotlib.pyplot as plt
import math
import cmath
import numpy as np
from scipy import special as sp

'''
This file contains functions used for determining surface occpancy,
bulk diffusion, and surface diffusion


Surface Occupancy: https://www.evernote.com/l/AuEmmLn7gv1CIZv05EO63PCszdrEk36WoQQ
Diffusion Functions: https://www.evernote.com/l/AuEqWuKaidtFJZmolVFtMRqjJAsD7aFz87I

'''


pi2=np.pi*2.0

u=1.660539e-27#atomic mass unit in kg
eb=2.28 #binding energy in kelvin
m3=3.0160293*u # mass of 3He in Kg

m_s=1.53*m3 #effective mass of surface
m_b=2.34*m3 #effective mass of bulk
hbar=1.0545718e-34 #J*s planks constant
k_B=1.380649e-23 #J/K boltzmanns constant

#Returns surface number density n_s (in layers) as long as N>>A*n_s
#WARNING: this has no checks for the above condition
#n_b is the bulk number density in units of cm^-2 and T is in kevlin
def Getn_s_(n_b,T,N=10):

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


#This function is essentially the same as Getn_s, but with different inputs and added checks
#N3= total number of atoms, S=surface area (cm^2), V=volume (cm^3)
#n=number of iterations for fixed point method
def Getn_s(N3,S,V,T):

    s=0.5 #spin of system
    e_b=2.28 #Surface binding energy in Kelvin
    nl=6.4e14 #one layer surface number density (cm^-2)

    #n_b=n_b*1e6 #converts to m^-2

    #Gets n_s for the case when N_s<<N_b
    n_s=Getn_s_(N3/V,T,10)

    #Worst case
    Sfrac=np.max(n_s*nl*S/N3)
    print("Maximum fraction on surface (this should be much smaller than 1): "+"{0:.2g}".format(Sfrac))

    if Sfrac>0.1: print("Warning! n_s may not be valid")

    return n_s #in  layers



#Gets fraction of 3He on surface assuming 2D high temp limit (Nondegenerate,Noninteracting)
#S in cm, V in cm, T in kelvin
def GetSfrac(S,V,T):

    V=V*1e-6#cm^3 -> m^3
    S=S*1e-4#cm^2 -> m^2

    alpha=np.power(S/V,2.)*np.power(m_s,2)*np.power(m_b,-3)*(2.*np.pi*np.power(hbar,2.)/k_B)
    #print(alpha)

    Tb=2.28 #binding energy in Kelvin

    frac=1./(1.+np.sqrt(T/alpha)*np.exp(-Tb/T))

    return frac

#Returns low temperature surface diffusion constant as extrapolated in reference above
# l is in layers, T is in kelivn, returns D in cm^2/sec
def LowTDs(l,T):

    D=(3.5e-3)*((l/T)**2)/np.log(2.1*l/T)

    return D

#Returns high temperature surface diffusion constant as explained in reference above
# l is in layers, T is in kelivn, returns D in cm^2/sec
def HighTDs(l,T):

    D=(6.6e-5)*np.sqrt(T)/l

    return D

#Emperical Bulk Diffusion Constant
#x is 3He concentration, T is in kelvin, returns D in cm^2/sec
def Db(x,T):

    D=(1.+5.6/T)*(1e-5)/x

    return D
