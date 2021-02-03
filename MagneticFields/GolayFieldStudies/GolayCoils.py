import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import magpylib as magpy

#Creates Gradient in the +y direction (Gy=+|Gy|)

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

#print(c.sources)

#magpy.displaySystem(c,direc=True)

# r=[[0,y,0] for y in np.linspace(-a,a,100)]
#
# r=np.asarray(r)
#
# #Bfield=c.getB(r)*10. #convert from mT to guass
# GBz=GetGrad(c,r,0.0001*a,1).reshape([100,3]) *10*10 #G/cm
# G0=GetGrad(c,[[0,0,0]],0.0001*a,2)[0,1]*10*10#[0,1]*10*10 #G/cm
# GBz=GBz/G0
# # GBzn=np.linalg.norm(GBz,axis=1)
# #
# #
# # #plt.plot(r[:,2]/a,Bfield[:,2]) #plot z comp of field along z axis
# plt.plot(r[:,1]/a,GBz[:,2]) #plot z comp of field along z axis
# plt.xlabel("z/a (unitless)")
# plt.ylabel(r'$\partial_yB_z/G_y$ (unitless)')
# plt.title(r"Golay Coil $\partial_yB_z$ along the y-axis (a=10cm, I=1 Amp) ")
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
ax2.set_title( 'Golay Coil Field, x=0 plane (a=10 cm, I= 1 Amp)')
ax2.set_xlabel('z/a')
ax2.set_ylabel('y/a')
#
# # plot field lines on ax2
#Z,Y = np.meshgrid(zs,ys)
Bz,By = Bs[:,:,2], Bs[:,:,1]
ax2.streamplot(zs/a,ys/a,Bz,By,color='k',density=1)
# #ax2.quiver(X,Z,U,V)

plt.show()
