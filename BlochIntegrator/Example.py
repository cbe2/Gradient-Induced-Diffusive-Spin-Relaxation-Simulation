import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R

def Bcirc(t,B0,B1,omega):
    return np.asarray([B1*np.cos(w0*t),B1*np.sin(w0*t),B0])

def Blinear(t,B0,B1,omega):
    return np.asarray([B1*np.sin(w0*t),0,B0])

#B in guass
def IntBloch(m0,Bfield,times,gamma,return_all=False):

    m=[]

    for i in range(len(times)):
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
B1= 0.01 #gauss (same as depol)
gamma=(2*np.pi)*(3240)
w0=B0*gamma
f0=w0/(2.*np.pi)
#tipping deg/sec
omega1=np.degrees(gamma*B1)#use /2 for 'linear field'

print("Frequency")
print(f0)
N=100 #number of points per peroid
dt=1./(f0*N)
print(dt)
tipAngle=10. # degree
T=tipAngle/omega1#100./f0 #Time in terms of number of periods

times=np.arange(0,T,dt)

ms=IntBloch(m0,lambda x: Bcirc(x,B0,B1,w0),times,gamma)

# ms.reshape(len(times),3)
# mTnorms=np.linalg.norm(ms[:,0:2],axis=1) #transverse norms
# thetas=np.degrees(np.arctan2(ms[:,2],mTnorms))
mTnorm=np.linalg.norm(ms[0:2])
theta=np.degrees(np.arctan2(mTnorm,ms[2]))

print("Final angle")
print(theta)
print("Expected Tipping Angle")
print(tipAngle)

#plt.plot(times,thetas)

#plt.show()
