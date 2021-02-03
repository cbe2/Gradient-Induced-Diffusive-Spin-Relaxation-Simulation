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
I=1 #current in amps


# create collection of two magnets
s1 = magpy.source.current.Circular( curr = I, dim =2.*a, pos=[0,0,a*np.sqrt(3)/2.])
s2 = magpy.source.current.Circular( curr = -I, dim =2.*a, pos=[0,0,-a*np.sqrt(3)/2.])

c = magpy.Collection(s1,s2)


# # create positions
bound=a*0.6
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
#GBx= GetGrad(c,R,0.0001*a,0).reshape([Np,Np,3]) #y,z, [Bx,By,Bz]
#GBy= GetGrad(c,R,0.0001*a,1).reshape([Np,Np,3])
G0=GetGrad(c,[[0,0,0]],0.0001*a,2)[0,2]*10*10 #G/cm
print(G0)
GBz=GetGrad(c,R,0.0001*a,2).reshape([Np,Np,3])
B0=np.linalg.norm(c.getB([0,0,0]))

#GBxn=np.linalg.norm(GBx,axis=2)*10# mT/mm --> mT/cm
#GByn=np.linalg.norm(GBy,axis=2)*10# mT/mm --> mT/cm
GBzn=np.linalg.norm(GBz,axis=2)*10*10# mT/mm --> G/cm
GBzn=GBzn/G0
#Gtot=(GBxn**2+GByn**2)/(B0**2)
#Gtot=1./Gtot


fig = plt.figure(figsize=(10,5)) #fig size
# ax1 = fig.add_subplot(121,projection='3d')
ax2 = fig.add_subplot(111)

#cp=ax2.pcolor((zs-0.5*dzs)/a,(ys-0.5*dys)/a,GBzn,cmap='jet',norm=LogNorm(vmin=1e-5, vmax=GBzn.max()))
cp=ax2.pcolor((zs-0.5*dzs)/a,(ys-0.5*dys)/a,GBzn,cmap='jet',vmin=GBzn.min(), vmax=GBzn.max())
cbar=fig.colorbar(cp, ax=ax2)
cbar.ax.set_ylabel(r'$|\nabla B_z|/G_z$', rotation=270,labelpad=15)


ax2.set_title( 'Maxwell Coil Field Gradients with I=1 amp, a=10 cm, x=0 plane')
ax2.set_xlabel('z/a')
ax2.set_ylabel('y/a')

plt.show()
