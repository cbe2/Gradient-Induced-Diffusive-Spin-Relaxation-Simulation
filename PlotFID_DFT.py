import matplotlib.pyplot as plt
import math
import cmath
import numpy as np
from DFTSuite import *

#aquire file names for plotting
fnames=getFilenames("SimulationData/PlotFolder/")
#fnames=['SimulationData/PlotFolder/Lx=1Ly=1Lz=1.txt','SimulationData/PlotFolder/Lx=2Ly=2Lz=2.txt','SimulationData/PlotFolder/Lx=3Ly=3Lz=3.txt','SimulationData/PlotFolder/Lx=4Ly=4Lz=4.txt',]

#imports data
signals=[]
for f in fnames:
    signals.append(FIDdata(f))
    #signals[-1].show()


#performs DFTs
pad=5 #seconds of padding on data set
sDFTs=[]
for S in signals:
    sDFTs.append(FIDFFTWindow(S,PadTime=pad))

#plots DFTs
for i in range(len(sDFTs)):
    plt.plot(sDFTs[i].Freqs-250.,sDFTs[i].Amps,label=signals[i].name[:-4],linestyle='-')



plt.xlim([-5,5])
#plt.ylim([0,1])
plt.grid()
plt.legend()#fname1,fname2])
plt.title("3D Diagonal Gradient FID Comparisons")
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude')
plt.show()
