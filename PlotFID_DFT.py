import matplotlib.pyplot as plt
import math
import cmath
import numpy as np
from DataSuite import *
from matplotlib.widgets import Cursor
from scipy.fftpack import fft,ifft
import sys
sys.path.insert(0, "/Users/cameronerickson/Desktop/Academics/PythonLibraries/DataCursor")
from DataCursor import DataCursor


#GradientCoilScan#SignalConfirm

#path="SignalConfirm/"#"SignalConfirm/"
#path="../FakeSig/"
#fname1="2018-11-16-14-05-20-.17.txt"
#Background=FIDdata(path+fname1)
#Background.show()

path="../../MonteCarlo/"
fname1="L=1f0=250N=2e4.txt"
Signal1=FIDdata(path+fname1)
#Signal1.show()

path="../FakeSig/"
fname2="L=3_f0=250.txt"
Signal2=FIDdata(path+fname2)
#Signal2.show()
#
#path="../../MonteCarlo/"
#fname3="N1000Anl.txt"
#Signal3=FIDdata(path+fname3)
#Signal3.show()
#
# path="../FakeSig/"
# fname4="NewSol_L=4.txt"
# Signal4=FIDdata(path+fname4)
# Signal4.show()
#
# path="../FakeSig/"
# fname5="FullSol_L=6.txt"
# Signal5=FIDdata(path+fname5)
# Signal5.show()
#
# path="../FakeSig/"
# fname6="NewSol_L=6.txt"
# Signal6=FIDdata(path+fname6)
# Signal5.show()

Tstart=0
Twindow=1


pad=0
#Bfft=FIDFFTWindow(Background,Times=[Tstart,Tstart+Twindow],PadTime=pad,BackSub=False)
Sfft1=FIDFFTWindow(Signal1,TIndexes=[0,Signal1.N-1],PadTime=pad,BackSub=False,f0=250.)
#Sfft2=FIDFFTWindow(Signal2,TIndexes=[0,Signal2.N-1],PadTime=pad,BackSub=False,f0=250)
#Sfft3=FIDFFTWindow(Signal3,TIndexes=[0,Signal3.N-1],PadTime=pad,BackSub=False)
#Sfft2=FIDFFTWindow(Signal2,Times=[Tstart,Tstart+Twindow],PadTime=pad,BackSub=False)
# Sfft3=FIDFFTWindow(Signal3,Times=[Tstart,Tstart+Twindow],PadTime=pad,BackSub=False)
# Sfft4=FIDFFTWindow(Signal4,Times=[Tstart,Tstart+Twindow],PadTime=pad,BackSub=False)
# Sfft5=FIDFFTWindow(Signal5,Times=[Tstart,Tstart+Twindow],PadTime=pad,BackSub=False)
# Sfft6=FIDFFTWindow(Signal6,Times=[Tstart,Tstart+Twindow],PadTime=pad,BackSub=False)



#print(len(Bfft.Freqs))
#print(len(Bfft.Amps))
#plt.plot(Bfft.Freqs*c,Bfft.Amps,label="Background")
plt.plot(Sfft1.XFreqs/1.,Sfft1.Amps/(Signal1.BT2*2.*np.pi),label="1e5 Spin Simulation",linestyle='-')
#plt.plot(Sfft2.XFreqs/3,Sfft2.Amps/(Signal1.BT2*2.*np.pi),label="Kubo-Redfield Theory",linestyle='-')
#plt.plot(Sfft3.Freqs,Sfft3.Amps,label="Analytical",linewidth=3,alpha=.7,c='g')
# plt.plot(Sfft4.XFreqs/4.,Sfft4.Amps/(Signal2.BT2*2.*np.pi),label="L=4 Full Cell")
# plt.plot(Sfft5.XFreqs/6.,Sfft5.Amps/(Signal2.BT2*2.*np.pi),label="L=6 Point by Point")
# plt.plot(Sfft6.XFreqs/6.,Sfft6.Amps/(Signal2.BT2*2.*np.pi),label="L=6 Full Cell ")
#plt.plot(x1/BDx1,Sfft3.Amps/alpha1,label="x=Imag")
#plt.plot(x1/BDx1,Sfft4.Amps,label="x=Real")
#plt.plot(x1/BDx1,Sfft5.Amps,label="x=1")
#plt.plot(x1/BDx1,Sfft6.Amps,label="x=2.0")


plt.xlim([-2,3])
#plt.ylim([0,1])
plt.grid()
plt.legend()#fname1,fname2])
plt.title("L=3 FID Comparisons")
plt.xlabel('Frequency (unitless position)')
plt.ylabel('Imaginary Part of DFT')
plt.show()
