import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import magpylib as magpy

# create collection of two magnets
s1 = magpy.source.current.Circular( curr = 1, dim =8)
s1.rotate(45,[0,1,0],anchor=s1.position) #tilts loop 45 degree about y-axis that passes through its center
c = magpy.Collection(s1)

# create positions
xs = np.linspace(-8,8,100) ; dxs=xs[1]-xs[0]
zs = np.linspace(-6,6,100) ; dzs=zs[1]-zs[0]
posis = [[x,0,z] for z in zs for x in xs] #increments x first



# calculate field and amplitude
B = [c.getB(pos) for pos in posis] #
Bs = np.array(B).reshape([100,100,3]) #reshape to [z_pos,x_pos,[Bx,By,Bz]]
Bamp = np.linalg.norm(Bs,axis=2) #comuptes norm of all vectors
Bamp=Bamp*10 #converts from mT to Gauss

#field at center in guass
print("Field at center in Gauss")
print(np.linalg.norm(c.getB([0,0,0])*10))

# define figure with a 2d and a 3d axis
fig = plt.figure(figsize=(10,5)) #fig size
ax1 = fig.add_subplot(121,projection='3d')
ax2 = fig.add_subplot(122)

# add displaySystem on ax1
magpy.displaySystem(c,subplotAx=ax1,suppress=True)
#ax1.view_init(elev=75)

# amplitude plot on ax2
#shift the coordinates b/c they are used as the corners.
c=ax2.pcolor(xs-0.5*dxs,zs-0.5*dzs,Bamp,cmap='jet',norm=LogNorm(vmin=Bamp.min(), vmax=Bamp.max()))
cbar=fig.colorbar(c, ax=ax2)
cbar.ax.set_ylabel('Gauss', rotation=270)

ax2.set_title('y=0 plane')
ax2.set_xlabel('x')
ax2.set_ylabel('z')

# plot field lines on ax2
X,Z = np.meshgrid(xs,zs)
U,V = Bs[:,:,0], Bs[:,:,2] #x and z components of the field
ax2.streamplot(X,Z,U,V,color='k',density=1)
#ax2.quiver(X,Z,U,V)

#display
plt.show()
