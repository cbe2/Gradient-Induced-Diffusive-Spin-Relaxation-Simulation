import matplotlib.pyplot as plt
import math
import numpy as np
import seaborn as sns
import pandas as pd
import sys
sys.path.append('../')
from WalkerClasses import *
import time
start_time = time.time()

'''
This code plots a 2D histogram of the spatial positions for a
combined set of bulk and surface walkers
'''

#np.random.seed(5)
#random.seed(5)

#Box dimensions
Lx=1.;Ly=1.;Lz=1.

#B is for bulk
paramsDictB={
'D':1.0,
'r0':[Lx,Ly,Lz],
'dt': 0.5e-4,
'L': [1.,1.,1.]
}

#S is for surface
paramsDictS={
'D':1.0,
'r0':[Lx,Ly],
'dt': 0.5e-4,
'L': [1.,1.]
}

# surface spans x,y while bulk spans x,y,z. z=0 is the top of the
# bulk (corresponding to the surface) and z increases as one moves
# downward

Ns=int(5e2) #number of surface walkers
Nb=int(5e2) #number of bulk walkers
N=Ns+Nb
steps=1000

walkers=[]

#Creates Bulk walkers
for i in range(Nb):
    paramsDictB['r0']=[np.random.uniform(0,Lx),np.random.uniform(0,Ly),np.random.uniform(0,Lz)]
    walkers.append(BoxWalk(paramsDictB))

#Creates Surface Walkers
for i in range(Ns):
    paramsDictS['r0']=[np.random.uniform(0,Lx),np.random.uniform(0,Ly)]
    walkers.append(BoxWalk(paramsDictS))

# steps walkers
for i in range(steps):
    for W in walkers: W.step(lambda x:0)


Xsamples=[]
Ysamples=[]
for W in walkers:
    Xsamples.append(W.ptcl.r[1])
    if len(W.ptcl.r)==3:
        Ysamples.append(W.ptcl.r[2])
    else:
        Ysamples.append(0)

counts, xedges, yedges, im=plt.hist2d(Xsamples,Ysamples,bins=20)
plt.colorbar()


plt.title(r'$10^{3}$ Spins with $N_s/N_3=0.5$')
plt.xlabel("Y position (cm)")
plt.ylabel("Z position (cm)")
plt.ylim(1,0)

TimeTaken=(time.time() - start_time)/60. #in minutes
print("{0:.2g}".format(TimeTaken)+" minutes")
#print("{0:.2g}".format(TimeTaken/float(steps*N)) + " min per step*particle")

plt.show()
