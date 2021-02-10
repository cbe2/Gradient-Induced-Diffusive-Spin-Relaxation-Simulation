import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np
import sys
sys.path.append('../')
from WalkerClasses import *

np.random.seed(2)
#random.seed(5)

paramsDict={
'D':1,
'r0':[0,0,0.6],
'dt': 1e-4,
'L': [1.,1.,1.0],
'SL': [0.1,0.1]
}

N=1
steps=10000

walkers=[]
x=[]
y=[]
z=[]

for i in range(N):
    #paramsDict['r0']=[np.random.uniform(0,L)]
    walkers.append(CellWalk(paramsDict))
for i in range(steps):
    for W in walkers:
        x.append(W.ptcl.r[0])
        y.append(W.ptcl.r[1])
        z.append(W.ptcl.r[2])
        W.step(lambda x:0)

plt.plot(x,z,label="path")
#plt.plot(y,z,label="path")
plt.scatter(x[0],z[0],c='r',label='start')
plt.xlabel('x')
plt.ylabel('z')
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.plot(x, y, z, label='parametric curve')
# ax.legend()
plt.show()
