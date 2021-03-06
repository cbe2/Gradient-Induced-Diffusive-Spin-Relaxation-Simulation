import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import magpylib as magpy

#Y= direction normal to earth's surface
#Z= Holding field direction
#X= determined by Y and Z

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

# create collection of two magnets
s1 = magpy.source.current.Circular( curr = I, dim =2.*a, pos=[0,0,a*np.sqrt(3)/2.])
s2 = magpy.source.current.Circular( curr = -I, dim =2.*a, pos=[0,0,-a*np.sqrt(3)/2.])
#ensuring d=a for homogenous magnetic field from hemoltz coils
#s1.move([0,0,a/2.])
#s2.move([0,0,-a/2.])

c = magpy.Collection(s1,s2)

#magpy.displaySystem(c,direc=True)

# r=[[0,0,z] for z in np.linspace(-a*2.,a*2.,100)]
#
# r=np.asarray(r)
# Bfield=c.getB(r)*10. #convert from mT to guass
# GBz=GetGrad(c,r,0.0001*a,2).reshape([100,3]) *10*10 #G/cm
# G0=GetGrad(c,[[0,0,0]],0.0001*a,2)[0,2]*10*10 #G/cm
# GBz=GBz/G0
# GBzn=np.linalg.norm(GBz,axis=1)
#
#
# #plt.plot(r[:,2]/a,Bfield[:,2]) #plot z comp of field along z axis
# plt.plot(r[:,2]/a,GBz[:,2]) #plot z comp of field along z axis
# plt.xlabel("z/a (unitless)")
# plt.ylabel(r'$\partial_zB_z/G_z$ (unitless)')
# plt.title(r"Maxwell Coil $\partial_zB_z$ along the z-axis (a=10cm, I=1 Amp) ")
# plt.grid()
# plt.show()

# # create positions
ys = np.linspace(-a*.8,a*.8,100) ; dys=ys[1]-ys[0]
zs = np.linspace(-a*.8,a*.8,100) ; dzs=zs[1]-zs[0]
posis = [[0,y,z] for y in ys for z in zs] #increments last variable first


# # calculate field and amplitude
B = [c.getB(pos) for pos in posis] #
Bs = np.array(B).reshape([100,100,3]) #reshape to [z_pos,y_pos,[Bx,By,Bz]]
Bamp = np.linalg.norm(Bs,axis=2) #comuptes norm of all vectors

Bamp=Bamp*10# to Gauss /np.linalg.norm(c.getB([0,0,0])) #normalize to center value

#Bamp=Bamp*10 #converts from mT to Gauss

# #field at center in guass
# print("Field at center in Gauss")
# print(np.linalg.norm(c.getB([0,0,0])*10))
#
# # define figure with a 2d and a 3d axis
fig = plt.figure(figsize=(10,5)) #fig size
# ax1 = fig.add_subplot(121,projection='3d')
ax2 = fig.add_subplot(111)
#
# # add displaySystem on ax1
# magpy.displaySystem(c,subplotAx=ax1,suppress=True)
# #ax1.view_init(elev=75)
#
# amplitude plot on ax2
#shift the coordinates b/c they are used as the corners.
#cp=ax2.pcolor((zs-0.5*dzs)/a,(ys-0.5*dys)/a,Bamp,cmap='jet',norm=LogNorm(vmin=Bamp.min(), vmax=Bamp.max()))
cp=ax2.pcolor((zs-0.5*dzs)/a,(ys-0.5*dys)/a,Bamp,cmap='jet',vmin=Bamp.min(), vmax=Bamp.max())
cbar=fig.colorbar(cp, ax=ax2)
cbar.ax.set_ylabel(r'Intensity $|B|$ (Gauss)', rotation=270,labelpad=15)


#
ax2.set_title( 'Maxwell Coil Field, x=0 plane (a=10 cm, I= 1 Amp)')
ax2.set_xlabel('z/a')
ax2.set_ylabel('y/a')
#
# # plot field lines on ax2
#Z,Y = np.meshgrid(zs,ys)
Bz,By = Bs[:,:,2], Bs[:,:,1]
ax2.streamplot(zs/a,ys/a,Bz,By,color='k',density=1)
# #ax2.quiver(X,Z,U,V)

plt.show()
