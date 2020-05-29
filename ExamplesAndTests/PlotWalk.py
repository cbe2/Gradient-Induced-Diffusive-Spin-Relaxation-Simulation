import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import math
import numpy as np
import sys
sys.path.append('../')
from WalkerClasses import *



paramsDict={
'D':1.0,
'r0':[0.1,0.5,0.5],
'dt': 1e-2,
'L': [1.,1.,1.]
}

w=BoxWalk(paramsDict)

x_pos=[]
y_pos=[]
z_pos=[]
t_pos=[]
for i in range(100):

    x_pos.append(w.ptcl.r[0])
    y_pos.append(w.ptcl.r[1])
    z_pos.append(w.ptcl.r[2])
    t_pos.append(float(w.ptcl.steps*w.dt))

    w.step(lambda x:0)

print(x_pos[0])
#plt.plot(t_pos,x_pos)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter3D([x_pos[0],x_pos[-1]],[y_pos[0],y_pos[-1]],[z_pos[0],z_pos[-1]],label="start/end points",c='r',s=100)
ax.plot(x_pos, y_pos,z_pos)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()
#plt.savefig('foo.png', bbox_inches='tight')
plt.show()
