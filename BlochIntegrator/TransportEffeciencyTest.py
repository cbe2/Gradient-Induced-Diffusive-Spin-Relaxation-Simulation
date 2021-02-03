import numpy as np
import matplotlib.pyplot as plt

u=1.660539e-27#atomic mass unit in kg
eb=2.28 #binding energy in kelvin
m3=3.0160293*u # mass of 3He in Kg

m_eff=2.34*m3 #effective mass of bulk
hbar=1.0545718e-34 #J*s planks constant
k_B=1.380649e-23 #J/K boltzmanns constant

#returns gamma for cylinder
#Pd=probability per bounce
#T=temp in Kelvin
#r=radius in cm of pipe
def WallGamma(Pd,R,T):

    #average v in m/s
    v=np.sqrt(8.*k_B*T/(np.pi*m_eff))

    #S/V in m^-1
    SV=2./(R*1e-2)

    return Pd*SV*v/4.

#D in cm^2/sec
#B0 in Gauss
#gradients in Gauss/cm
#returns s^-1
def SlowDLimGamma(D,B0,Gx,Gy):
    Gtot2=Gx**2+Gy**2
    return 2.*D*Gtot2/(B0**2)


#D in cm^2/sec
#B0 in Gauss
#gradients in Gauss/cm
#returns s^-1
def FastDLimGammaFastT(D,B0,Gx,Gy):
    Gtot2=Gx**2+Gy**2
    return D*Gtot2/(B0**2)

#L in cm
def FastDLimGammaSlowT(D,B0,Gx,Gy,L):
    Gtot2=Gx**2+Gy**2
    gamma=(2*np.pi)*(3240)
    return Gtot2*(gamma**2)*(L**4)/(120.*(D))


#refelcting wall at x=0, spins transported to x=L
#Gamma in s^-1
#D in cm^2/sec
#L in cm
def ReflectingWall(Gamma,D,L,x0):
    arg=np.sqrt(Gamma/D)
    return np.cosh(x0*arg)/np.cosh(L*arg)

def NoWall(Gamma,D,L,x0):
    arg=np.sqrt(Gamma/D)
    return np.exp(-L*arg)



#print(WallGamma(Pd=1e-7,R=1.5,T=0.3))


G=WallGamma(Pd=1e-6,R=1.5,T=0.3)

eff=ReflectingWall(Gamma=G,D=1000,L=200,x0=0)
print("Percent lost to wall depolarization: "+ str((1.-eff)*100))

Gxs= np.logspace(-6, -3, 1000, endpoint=True)

D=1000 #cm^2/sec
B0=0.03 #Gauss
R=1.5 #cm
L=200 #cm

Gamma1=SlowDLimGamma(D,B0,Gxs,0)
Gamma2=FastDLimGammaFastT(D,B0,Gxs,0)
Gamma3=FastDLimGammaSlowT(D,B0,Gxs,0,2*R)

plt.semilogx(Gxs,ReflectingWall(Gamma1,D,L,0),label="Slow Diffusion")
plt.semilogx(Gxs,ReflectingWall(Gamma2,D,L,0),label="Fast Diffusion Fast Perioid")
plt.semilogx(Gxs,ReflectingWall(Gamma3,D,L,0),label="Fast Diffusion Slow Perioid")

plt.ylabel("Fraction of Surviving Spins")
plt.xlabel("Gradient (Gauss/cm)")

plt.legend()
plt.show()

# plt.plot(pos*G/B1,thetas,label='Simulated Tip Angles')
#
# plt.plot(pos*G/B1,phases,label='Simulation Phases')
# plt.plot(pos*G/B1,expPhases,label='Expected Phases')
# #plt.plot(pos*G/B1,-np.degrees(pos*0.5*G*np.radians(10)/B1),label='linearphases')
# plt.ylabel("Angle (degrees)")
# plt.xlabel(r"position$\;\times\; G/B_1$ (unitless)")
# plt.title("Tip angles and Phases from 10 degree kicking Pulse")
# #plt.plot(pos*G/B1,expThetas3,label='small angle')
# #plt.plot(pos*G/B1,np.abs(np.sinc(pos*G/(B1*3.875)))**2 )
# plt.legend()
# plt.grid()
# plt.show()
