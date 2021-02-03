import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import magpylib as magpy
import math

#Y = direction normal to earth's surface
#Z = Holding field direction
#X = determined by Y and Z

'''
Plots worst case T1 time in holding field for a sphere and cube (cube commented out)
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

#generates points on a sphere of radius r
def fibonacci_sphere(samples,r):

    points = []
    phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians

    for i in range(samples):
        y = (1 - (i / float(samples - 1)) * 2)*r  # y goes from r to -r
        radius = math.sqrt(r*r - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append([x, y, z])

    return points



a=10. # radius parameter of coils (cm)
a=a*10.# change to mm
I=1 #current in amps


# create collection of two magnets
s1 = magpy.source.current.Circular( curr = I, dim =2.*a, pos=[0,0,a/2.])
s2 = magpy.source.current.Circular( curr = I, dim =2.*a, pos=[0,0,-a/2.])

c = magpy.Collection(s1,s2)

#markerPos = [(0,0,0,'cell location')]

#magpy.displaySystem(c,markers=markerPos,direc=True)
#plt.title("Hemholtz Holding Coils")
#plt.show()
print("center field in gauss")
print(np.linalg.norm(c.getB([0,0,0]))*10)

L=np.linspace(0.1,0.6,30)#[0.1,0.2,0.3,0.4,0.5]
T1s=[]
Gmaxs=[]

for l in L:
# # create positions
    bound=a*l
    radius=a*l
    Np=500
    #xs = np.linspace(-bound,bound,Np) ;
    #ys = np.linspace(-bound,bound,Np) ; #dys=ys[1]-ys[0]
    #zs = np.linspace(-bound,bound,Np) ; #dzs=zs[1]-zs[0]
    #R = [[x,y,z] for x in xs for y in ys for z in zs] #increments last variable first
    R=fibonacci_sphere(Np,radius)


    R=np.asarray(R)
    # plt.show()
    #
    # returns gradients in mT/mm
    #GBx= GetGrad(c,R,0.0001*a,0).reshape(Np,3)#reshape([Np,Np,Np,3]) #y,z, [Bx,By,Bz]
    #GBy= GetGrad(c,R,0.0001*a,1).reshape(Np,3)#reshape([Np,Np,Np,3])
    GBz= GetGrad(c,R,0.0001*a,2).reshape(Np,3)
    B0=np.linalg.norm(c.getB([0,0,0]))

    #GBxn=np.linalg.norm(GBx,axis=1)*10#np.linalg.norm(GBx,axis=3)*10# mT/mm --> mT/cm
    #GByn=np.linalg.norm(GBy,axis=1)*10#np.linalg.norm(GBy,axis=3)*10# mT/mm --> mT/cm
    GBzn=np.linalg.norm(GBz,axis=1)*10*10#np.linalg.norm(GBy,axis=3)*10# mT/mm --> G/cm
    #Gtot=(GBxn**2+GByn**2)/(B0**2)
    #Gtot=1./Gtot
    #T1s.append(Gtot.min())
    Gmaxs.append(GBzn.max())

T1s=np.asarray(T1s)
Gmaxs=np.asarray(Gmaxs)
print(T1s)
plt.semilogy(L,Gmaxs)
#plt.title(r"Minimum $T_1$ time in a sphere of radius r (D=1 cm$^2$/sec) (a=10 cm)")
plt.title(r"Maximum $|\nabla B_z|$ in a sphere of radius r (I=1 amp) (a=10 cm)")
plt.xlabel("r/a")
plt.ylabel(r"$|\nabla B_z|$ (Guass/cm)")
plt.grid()
plt.show()
