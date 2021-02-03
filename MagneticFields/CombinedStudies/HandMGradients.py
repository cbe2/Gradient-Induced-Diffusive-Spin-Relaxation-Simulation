import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import magpylib as magpy

#Y = direction normal to earth's surface
#Z = Holding field direction
#X = determined by Y and Z

'''
Creates color map of T1 times or Bz gradients
'''

#returns the gradient from a collection
#c= collection generating the B-field
#r= given positions to evaluate gradient at r=[[x1,y1,z1],[x2,y2,z2],...] units in mm
#comp=index that gives the proper component (0=x, 1=y, 2=z)
def GetGrad(c,R,step,comp):

    dx=np.asarray([step,0,0])
    dy=np.asarray([0,step,0])
    dz=np.asarray([0,0,step])

    gradients=[]
    for r in R:

        xcomp=(c.getB(r+dx)[comp]-c.getB(r)[comp])/step
        ycomp=(c.getB(r+dy)[comp]-c.getB(r)[comp])/step
        zcomp=(c.getB(r+dz)[comp]-c.getB(r)[comp])/step

        gradients.append([xcomp,ycomp,zcomp])

    return np.asarray(gradients) #units of mT/mm

a=10. # radius parameter of coils (cm)
a=a*10.# change to mm
I=1000 #current in amps
IM=100


s1 = magpy.source.current.Circular( curr = I, dim =2.*a, pos=[0,0,a/2.])
s2 = magpy.source.current.Circular( curr = I, dim =2.*a, pos=[0,0,-a/2.])
#ensuring d=a for homogenous magnetic field from hemoltz coils

s3 = magpy.source.current.Circular( curr = IM, dim =2.*a, pos=[0,0,a*np.sqrt(3)/2.])
s4 = magpy.source.current.Circular( curr = -IM, dim =2.*a, pos=[0,0,-a*np.sqrt(3)/2.])
#ensuring d=a\sqrt{3} for homogenous magnetic field from hemoltz coils

c = magpy.Collection(s1,s2,s3,s4)


# # create positions
bound=a*0.7
Np=100
#xs = np.linspace(-bound,bound,100) ;
ys = np.linspace(-bound,bound,Np) ; dys=ys[1]-ys[0]
zs = np.linspace(-bound,bound,Np) ; dzs=zs[1]-zs[0]
R = [[0,y,z] for y in ys for z in zs] #increments last variable first
#posis = [[0,y,0] for y in ys] #increments last variable first
R=np.asarray(R)
# plt.show()
#
# returns gradients in mT/mm
GBx= GetGrad(c,R,0.0001*a,0).reshape([Np,Np,3]) #y,z, [Bx,By,Bz]
GBy= GetGrad(c,R,0.0001*a,1).reshape([Np,Np,3])
#GBz=GetGrad(c,R,0.0001*a,2).reshape([Np,Np,3])
B0=np.linalg.norm(c.getB([0,0,0]))*10 # to Gauss

GBxn=np.linalg.norm(GBx,axis=2)*10*10# mT/mm --> G/cm
GByn=np.linalg.norm(GBy,axis=2)*10*10# mT/mm --> G/cm
#GBzn=np.linalg.norm(GBz,axis=2)*10*10# mT/mm --> G/cm
Gtot=(GBxn**2+GByn**2)/(B0**2)
T1=1./Gtot


fig = plt.figure(figsize=(10,5)) #fig size
# ax1 = fig.add_subplot(121,projection='3d')
ax2 = fig.add_subplot(111)

#cp=ax2.pcolor((zs-0.5*dzs)/a,(ys-0.5*dys)/a,T1,cmap='jet',norm=LogNorm(vmin=1e2, vmax=1e5))
cp=ax2.contourf(zs/a, ys/a, T1, 10, cmap='jet',norm=LogNorm(vmin=1e2, vmax=1e5) ,origin='lower')
#cp=ax2.pcolor((zs-0.5*dzs)/a,(ys-0.5*dys)/a,T1,cmap='jet',vmin=T1.min(), vmax=T1.max())
cbar=fig.colorbar(cp, ax=ax2)
cbar.ax.set_ylabel(r'local $T_1$ (sec)', rotation=270,labelpad=15)


ax2.set_title( '1000 Amp Hemholtz Coil + 100 Amp Maxwell Coil (a= 10cm, D= 1cm^2/sec)')
ax2.set_xlabel('z/a')
ax2.set_ylabel('y/a')

plt.show()
