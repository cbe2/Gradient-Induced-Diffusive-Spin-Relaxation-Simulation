import matplotlib.pyplot as plt
import math
import cmath
import numpy as np
from DFTSuite import *
from WalkerClasses import *

#aquire file names for plotting
fnames=getFilenames("SimulationData/PlotFolder/")
#fnames=['SimulationData/PlotFolder/Lx=1Ly=1Lz=1.txt','SimulationData/PlotFolder/Lx=2Ly=2Lz=2.txt','SimulationData/PlotFolder/Lx=3Ly=3Lz=3.txt','SimulationData/PlotFolder/Lx=4Ly=4Lz=4.txt',]

#imports data
signals=[]
for f in fnames:
    signals.append(FIDdata(f))
    #signals[-1].show()

# time,data=LoadData("SimulationData/2DNormalL=4Abs.txt")
# signals[0].data=signals[0].data/data

#performs DFTs
pad=0 #seconds of padding on data set
sDFTs=[]
for S in signals:
    sDFTs.append(FIDFFTWindow(S,PadTime=pad))
    #sDFTs.append(FIDFFTWindow(S,Times=[0,1],PadTime=pad))

#plots DFTs
for i in range(len(sDFTs)):
    #plt.plot(sDFTs[i].Freqs-250,-2.*np.imag(sDFTs[i].data)/sDFTs[i].N0)#,label=signals[i].name[:-4],linestyle='-')
    plt.plot(sDFTs[i].Freqs-signals[i].f0,2.*np.real(sDFTs[i].data)/sDFTs[i].N0,label=signals[i].LN)#,label=signals[i].name[:-4],linestyle='-')
    #plt.plot(sDFTs[i].Freqs,sDFTs[i].Amps,label=signals[i].name[:-4],linestyle='-')
    #plt.plot(sDFTs[i].Freqs-250,sDFTs[i].Phases)#,label=signals[i].name[:-4],linestyle='-')



#plt.xlim([-1,1])
#plt.ylim([0,1])
plt.grid()
plt.legend()#
plt.title(r' Bulk Signal D=2.7 cm$^2$/sec Vs D=0 cm$^2$/sec, L=4cm, G=1 G/cm')
plt.xlabel('Frequency (Hz) (offset is arbitrary)')
plt.ylabel('Real Part of the DFT (arbitrary units)')
#plt.ylabel('Imaginary part of DFT')
plt.show()
