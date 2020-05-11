import matplotlib.pyplot as plt
import math
import numpy as np
from WalkerClasses import *



w=Walk(1.0,.001)

x_pos=[]
t_pos=[]
for i in range(1000):

    x_pos.append(w.ptcl.x)
    t_pos.append(w.ptcl.t)

    w.step(type="Gauss")


plt.plot(t_pos,x_pos)
plt.title("1000 Gaussian step random walk")
plt.xlabel("Time")
plt.ylabel("Position")
plt.show()
