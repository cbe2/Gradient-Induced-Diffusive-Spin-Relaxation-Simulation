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


a=10 # radius parameter of coils (cm)
a=a*10.# change to mm
I=np.abs(1) #current in Amps
N=10 #number of theta points

# create collection of two magnets
angles=np.linspace(-np.pi/3.,np.pi/3.,N)

#clockwise points for arc in x-y plane in x>0
arcNearXplusV=np.asarray([[a*np.cos(theta),a*np.sin(theta),0.4*a] for theta in angles])
arcFarXplusV=np.asarray([[a*np.cos(-theta),a*np.sin(-theta),1.64*a] for theta in angles])
#connecting the arcs
NearFarV=np.asarray([[a*np.cos(angles[-1]),a*np.sin(angles[-1]),0.4*a],[a*np.cos(angles[-1]),a*np.sin(angles[-1]),1.64*a]])
FarNearV=np.asarray([[a*np.cos(angles[0]),a*np.sin(angles[0]),1.64*a],[a*np.cos(angles[0]),a*np.sin(angles[0]),0.4*a]])

#collection of wire loops
LCs=[]
for i in [I,-I,-I,I]:
    arcNearXplus = magpy.source.current.Line( curr = i, vertices=arcNearXplusV)
    arcFarXplus = magpy.source.current.Line( curr = i, vertices=arcFarXplusV)
    NearFar = magpy.source.current.Line( curr = i, vertices=NearFarV)
    FarNear = magpy.source.current.Line( curr = i, vertices=FarNearV)

    LCs.append(magpy.Collection(arcNearXplus,NearFar,arcFarXplus,FarNear))
#rotating them to proper positions
LCs[1].rotate(angle=180,axis=[1,0,0],anchor=[0,0,0])
LCs[2].rotate(angle=180,axis=[0,0,1],anchor=[0,0,0])
LCs[3].rotate(angle=180,axis=[0,1,0],anchor=[0,0,0])


c = magpy.Collection(*LCs)

#rotate collection to create Gy gradient instead of Gx
c.rotate(angle=90,axis=[0,0,1],anchor=[0,0,0])

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
