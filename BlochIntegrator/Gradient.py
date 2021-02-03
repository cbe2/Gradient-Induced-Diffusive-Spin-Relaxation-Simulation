import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R

#t=time, B0=holding field value, B1=Kicking field value, omega=RF angular freq
def Bcirc(t,B0,B1,w):
    return np.asarray([B1*np.cos(w*t),B1*np.sin(w*t),B0])

def Blinear(t,B0,B1,w):
    return np.asarray([B1*np.sin(w*t),0,B0])

#lobes=number of lobes on each side
def Bsinc(t,B0,B1,w,dw, lobes):
    tp=t-2.*np.pi*lobes/dw

    if np.abs(tp)<= 2.*np.pi*lobes/dw:
        return np.asarray([2.*B1*np.sin(w*tp)*np.sinc(dw*tp/(2.*np.pi)),0,B0])
    else:
        return np.asarray([0,0,B0])



#B in guass
def IntBloch(m0,Bfield,times,gamma,return_all=False):

    m=[]

    for i in range(len(times)-1):
        t=(times[i]+times[i+1])/2.
        dt=times[i+1]-times[i]

        b=Bfield(t)
        bmag=np.linalg.norm(b,axis=0)

        rot_axis = b/bmag
        angle=gamma*bmag*dt

        rotation=R.from_rotvec(angle*rot_axis)
        m0=rotation.apply(m0)

        if return_all==True: m.append(m0)

    if return_all==True: return np.asarray(m)
    else: return m0


m0=np.asarray([0,0,1.])

B0= 50. #gauss
B1= 0.1 #gauss (same as depol)
G=1.0 #gauss/cm
gamma=(2*np.pi)*(3240)
w0=B0*gamma
f0=w0/(2.*np.pi)
#tipping deg/sec
omega1=np.degrees(gamma*B1)#use /2 for 'linear field'
# print("Frequency")
# print(f0)
N=10 #number of points per peroid
dt=1./(f0*N)
tipAngle=10. # degree
T=tipAngle/omega1 # 100./f0 #Time in terms of number of periods
lobes=4
dl=5 #cm
dw=np.abs(gamma*G*dl)
maxT=2.*2.*np.pi*lobes/dw

times=np.arange(0,maxT,dt)

# plt.title("B1 Field Profile")
# plt.plot(times*1000,2.*B1*np.sinc(dw*(times-maxT/2)/(2.*np.pi)))
# plt.grid()
# plt.ylabel("Field Magnitude (Gauss)")
# plt.xlabel("Time (milli-sec)")
# plt.show()

#times=np.arange(0,T,dt)


pos=np.linspace(-dl*1.5,dl*1.5,300)#np.linspace(-50*B1/G,50*B1/G,50) #cm
thetas=[]
expThetas=[]
phases=[]
expPhases=[]

ms=IntBloch(m0,lambda x: Bsinc(x,B0,B1,B0*gamma,dw,lobes),times,gamma)
print("Final M on resonance")
print(ms)
print("Transverse Norm")
mTnorm=np.linalg.norm(ms[0:2])
print(mTnorm)
phi0=np.degrees(np.arctan2(ms[1],ms[0]))
print("Phi")
print(phi0)
#expThetas2=[]
#expThetas3=[]
for x in pos:

    print(x)
    localB0=B0+G*x
    omega=B0*gamma
    ms=IntBloch(m0,lambda x: Bsinc(x,localB0,B1,B0*gamma,dw,lobes),times,gamma)
    #print("M")
    #print(ms)

    mTnorm=np.linalg.norm(ms[0:2])
    phase=np.degrees(np.arctan2(ms[1],ms[0]))
    #print("phase")
    #print(phase)
    phase=phase-phi0
    #phase=np.tan(np.radians(phase))
    #print(phase)
    theta=np.degrees(np.arctan2(mTnorm,ms[2]))
    thetas.append(theta)
    phases.append(phase)

    dOmega=(G*x*gamma)

    omegaMag=np.power((G*x*gamma)**2+(gamma*B1)**2,0.5)
    #sinAngle2=(2./(1.+(G*x/B1)**2))

    #expTheta=1.-sinAngle2*(np.sin(0.5*times[-1]*omegaMag))**2
    #expThetas.append(np.degrees(np.arccos(expTheta)))
    #expThetas.append(-(expTheta-1))
    sinterm=np.abs((gamma*B1*times[-1]/2)*np.sinc(omegaMag*times[-1]/(np.pi*2.)))
    expThetas.append(2*np.degrees(np.arcsin(sinterm)))
    #expThetas3.append(np.degrees(2*sinterm))
    pNum = (dOmega/(gamma*B0))*times[-1]*np.sinc(omegaMag*times[-1]/(np.pi))
    pDenom=np.cos(omegaMag*times[-1])-1.
    o=-gamma*B0
    o0=-gamma*localB0
    o1=-gamma*B1
    oD=o0-o
    oM=np.power((oD)**2+(o1)**2,0.5)

    mx=(o1*oD/(oM**2))*(1-np.cos(oM*times[-1]))
    my=-(o1/oM)*np.sin(oM*times[-1])


    expPhases.append(np.degrees(np.arctan((oD/oM)*np.tan(oM*times[-1]/2.))))



plt.plot(pos,thetas)#,label='Simulated Tip Angles')
#plt.plot(pos*G/B1,expThetas,label='Expected Tip Angles')
#plt.plot(pos*G/B1,expThetas,label='Analytical Solution')
#print(phases)
#plt.plot(pos*G/B1,phases,label='Simulation Phases')
#plt.plot(pos*G/B1,expPhases,label='Expected Phases')
#plt.plot(pos*G/B1,-np.degrees(pos*0.5*G*np.radians(10)/B1),label='linearphases')
plt.ylabel("Angle (degrees)")
plt.xlabel(r"position (cm)")#plt.xlabel(r"position$\;\times\; G/B_1$ (unitless)")
plt.title("Simulated Tip angles from Sinc pulse, 4x2 zero crossings, bandwith=5cm ")
#plt.plot(pos*G/B1,expThetas3,label='small angle')
#plt.plot(pos*G/B1,np.abs(np.sinc(pos*G/(B1*3.875)))**2 )
plt.legend()
plt.grid()
plt.show()
