import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import magpylib as magpy

#Y= direction normal to earth's surface
#Z= Holding field direction
#X= determined by Y and Z

d=10 # length of square coil (cm)
d=d*10.# change to mm
I=np.abs(1)/26. #current in Amps

#points of the square centered at the orgin with current running clockwise
squarePts=[[-d/2.,-d/2.,0],[d/2.,-d/2.,0],[d/2.,d/2.,0],[-d/2.,d/2.,0],[-d/2.,-d/2.,0]]

# create collection of two magnets
#3-coil config:
# coil0 = magpy.source.current.Line( curr = I*39, vertices=squarePts)
# coil1 = magpy.source.current.Line( curr = I*20, vertices=squarePts)
# coil2 = magpy.source.current.Line( curr = I*39, vertices=squarePts)
#
# coil0.move([0,0,.4106*d])
# coil2.move([0,0,-.4106*d])

#4-coil config
coil0 = magpy.source.current.Line( curr = I*26, vertices=squarePts)
coil1 = magpy.source.current.Line( curr = I*11, vertices=squarePts)
coil2 = magpy.source.current.Line( curr = I*11, vertices=squarePts)
coil3 = magpy.source.current.Line( curr = I*26, vertices=squarePts)

coil0.move([0,0,.5055*d])
coil1.move([0,0,.1281*d])
coil2.move([0,0,-.1281*d])
coil3.move([0,0,-.5055*d])



#ensuring d=a for homogenous magnetic field from hemoltz coils
#s1.move([0,0,a/2.])
#s2.move([0,0,-a/2.])

c = magpy.Collection(coil0,coil1,coil2,coil3)

#magpy.displaySystem(c,direc=True)

r=[[0,0,z] for z in np.linspace(-d,d,100)]

r=np.asarray(r)
Bfield=c.getB(r)*10. #convert from mT to guass
Bfield=Bfield/(c.getB([0,0,0])[2]*10)

print(c.getB([0,0,0])[2]*10)

plt.plot(r[:,2]/d,Bfield[:,2]) #plot z comp of field along z axis
plt.show()

# # create positions
# xs = np.linspace(-8,8,100) ; dxs=xs[1]-xs[0]
# zs = np.linspace(-6,6,100) ; dzs=zs[1]-zs[0]
# posis = [[x,0,z] for z in zs for x in xs] #increments x first
#
#
# plt.show()
#
# # calculate field and amplitude
# B = [c.getB(pos) for pos in posis] #
# Bs = np.array(B).reshape([100,100,3]) #reshape to [z_pos,x_pos,[Bx,By,Bz]]
# Bamp = np.linalg.norm(Bs,axis=2) #comuptes norm of all vectors
# Bamp=Bamp*10 #converts from mT to Gauss
#
# #field at center in guass
# print("Field at center in Gauss")
# print(np.linalg.norm(c.getB([0,0,0])*10))
#
# # define figure with a 2d and a 3d axis
# fig = plt.figure(figsize=(10,5)) #fig size
# ax1 = fig.add_subplot(121,projection='3d')
# ax2 = fig.add_subplot(122)
#
# # add displaySystem on ax1
# magpy.displaySystem(c,subplotAx=ax1,suppress=True)
# #ax1.view_init(elev=75)
#
# # amplitude plot on ax2
# #shift the coordinates b/c they are used as the corners.
# c=ax2.pcolor(xs-0.5*dxs,zs-0.5*dzs,Bamp,cmap='jet',norm=LogNorm(vmin=Bamp.min(), vmax=Bamp.max()))
# cbar=fig.colorbar(c, ax=ax2)
# cbar.ax.set_ylabel('Gauss', rotation=270)
#
# ax2.set_title('y=0 plane')
# ax2.set_xlabel('x')
# ax2.set_ylabel('z')
#
# # plot field lines on ax2
# X,Z = np.meshgrid(xs,zs)
# U,V = Bs[:,:,0], Bs[:,:,2] #x and z components of the field
# ax2.streamplot(X,Z,U,V,color='k',density=1)
# #ax2.quiver(X,Z,U,V)
#
# #display
# plt.show()
