import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R

#t=time, B0=holding field value, B1=Kicking field value, omega=RF angular freq
def Bcirc(t,B0,B1,w):
    return np.asarray([B1*np.cos(w*t),B1*np.sin(w*t),B0])

def Blinear(t,B0,B1,w):
    return np.asarray([B1*np.cos(w*t),0,B0]) #proper rotating frame is -w0


#B in guass
def IntBloch(m0,Bfield,times,gamma,return_all=False):

    m=[]

    for i in range(len(times)-1):
        t=(times[i]+times[i+1])/2.
        dt=times[i+1]-times[i]

        b=Bfield(t)
        bmag=np.linalg.norm(b,axis=0)

        rot_axis = -b/bmag
        angle=gamma*bmag*dt

        rotation=R.from_rotvec(angle*rot_axis)
        m0=rotation.apply(m0)

        if return_all==True: m.append(m0)

    if return_all==True: return np.asarray(m)
    else: return m0


m0=np.asarray([0,0,1.])

B0= 50. #gauss
B1= 0.1 #gauss (same as depol)
G=0.1 #gauss/cm
gamma=(2*np.pi)*(3240)
w0=B0*gamma
f0=w0/(2.*np.pi)
omega1=np.degrees(gamma*B1)#tipping deg/sec
# print("Frequency")
# print(f0)
N=300 #number of points per peroid
dt=1./(f0*N)
tipAngle=10. # degree
T=tipAngle/omega1 # 100./f0 #Time in terms of number of periods
times=np.arange(0,T,dt)

pos=np.linspace(2*50*B1/G,2*50*B1/G,150) #cm
thetas=[]
expThetas=[]
phases=[]
expPhases=[]
finalMs=[]

ms=IntBloch(m0,lambda t: Bcirc(t,B0,B1,-B0*gamma),times,gamma)
resM=ms
phi0=np.degrees(np.arctan2(ms[1],ms[0]))

ApproxMs=[]

Xpts=[]
Ypts=[]
Xpts2=[]
Ypts2=[]


for x in pos:

    print("location")
    print(x)
    localB0=B0+G*x
    print("Field")
    print(localB0)
    omega=-B0*gamma
    ms=IntBloch(m0,lambda t: Bcirc(t,localB0,B1,omega),times,gamma)
    print("vector")
    print(ms)
    finalMs.append(ms)
    #print("M")
    #print(ms)

    #mTnorm=np.linalg.norm(ms[0:2])
    #phase=np.degrees(np.arctan2(ms[1],ms[0]))
    #print("phase")
    #print(phase)
    #phase=phase-phi0
    #phase=np.tan(np.radians(phase))
    #print(phase)
    theta=np.degrees(np.arccos(ms[2]))#np.arctan2(mTnorm,ms[2]))
    thetas.append(theta)
    #phases.append(phase)

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

    aphase=oD*times[-1]/2
    rePart=np.cos(aphase)*np.sinc(aphase/np.pi)*o1*T
    imPart=np.sin(aphase)*np.sinc(aphase/np.pi)*o1*T
    Xpts.append(rePart)
    Ypts.append(imPart)

    ApproxMs.append([rePart,imPart,0])
    # Xpts2.append(my)
    # Ypts2.append(mx)


    expPhases.append(np.degrees(np.arctan((oD/oM)*np.tan(oM*times[-1]/2.))))

#transforming Vectors to rotating frame:
finalMs=np.asarray(finalMs)
# ApproxMs=np.asarray(ApproxMs)
#
rot_axis = np.asarray([0,0,-1])
angle=-gamma*B0*T -(1.-.01)*np.pi/2.#
#
rotation=R.from_rotvec(angle*rot_axis)
finalMs=rotation.apply(finalMs)
resM=rotation.apply(resM)

# ax = plt.axes(projection='3d')
# zdata = finalMs[:,2]
# xdata = finalMs[:,0]
# ydata = finalMs[:,1]
# ax.set_xlabel('x');ax.set_ylabel('y');ax.set_zlabel('z');
# ax.plot3D(xdata, ydata, zdata);
# ax.scatter3D(resM[0],resM[1],resM[2],label="On resonance")

plt.figure()
plt.plot(pos*G/B1,thetas,label='Simulated Tip Angles')
plt.title("Off Resonance Tip Angles After 10 Degree Tip")
plt.ylabel("Angle (degrees)")
plt.xlabel(r"position$\;\times\; G/B_1$ (unitless)")
plt.grid()

plt.figure()
plt.plot(finalMs[:,0],finalMs[:,1],label='Simulation')
plt.scatter(resM[0],resM[1],label='On Resonance Simulation',c='g')
#plt.plot([resM[0],0],[resM[1],0],c='g')
plt.title("Transverse Plane Spin Compenents After 10 Degree Tip")
plt.plot(Xpts,Ypts,label="Small Angle Limit")
plt.ylabel("y-component of spin")
plt.xlabel(r"x-component of spin")
#plt.plot(Xpts2,Ypts2,label="|sinc|")
plt.legend(loc='upper left')
plt.grid()
#plt.plot(pos*G/B1,expThetas,label='Expected Tip Angles')
#plt.plot(pos*G/B1,expThetas,label='Analytical Solution')
#print(phases)
#plt.plot(pos*G/B1,phases,label='Simulation Phases')
#plt.plot(pos*G/B1,expPhases,label='Expected Phases')
#plt.plot(pos*G/B1,-np.degrees(pos*0.5*G*np.radians(10)/B1),label='linearphases')
#plt.ylabel("Angle (degrees)")
#plt.xlabel(r"position$\;\times\; G/B_1$ (unitless)")
#plt.title("Tip angles and Phases from 10 degree kicking Pulse")
#plt.plot(pos*G/B1,expThetas3,label='small angle')
#plt.plot(pos*G/B1,np.abs(np.sinc(pos*G/(B1*3.875)))**2 )

# Data for three-dimensional scattered points



#plt.legend()
#plt.grid()
plt.show()
