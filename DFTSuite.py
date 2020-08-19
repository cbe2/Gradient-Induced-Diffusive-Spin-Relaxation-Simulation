import matplotlib.pyplot as plt
import os
import math
import cmath
import numpy as np
from matplotlib.widgets import Cursor
from scipy.fftpack import fft,ifft
import datetime
from scipy.optimize import curve_fit
from scipy.special import ellipe
from scipy import stats
import uncertainties as unc
import uncertainties.umath as umath

#FIDdata EDITED for Fake Signals

#Stores the Raw data and auxiliary information
class FIDdata:
    def __init__(self,fname):

        f=open(fname,'r')
        header=False
        V=[]
        for line in f:
            if header==False:
                line=line.split("\t")
                if line[0]=="Legend Name:":
                    self.LN=line[1][:-1]
                if line[0]=="Sample Rate:":
                    self.SR=float(line[1])
                if line[0]=="Sample Size:":
                    self.SS=float(line[1])
                if line[0]=="Read Time:":
                    self.ReadTime=float(line[1])
                if line[0]=="Gradient Coil (V):":
                    self.G=line[1]#float(line[1])
                if line[0]=="Gradient (G/cm):":
                    self.G=float(line[1])
                if line[0]=="X-Gradient (G/cm):":
                    self.G=float(line[1])
                if line[0]=="Diffusion (cm^2/sec):":
                    self.D=float(line[1])
                if line[0]=="f0 (Hz):":
                    self.f0=float(line[1])
                if line[0]=="Date  and Time:":
                    info=line[1].split("-")
                    self.Tstamp=datetime.datetime(int(info[0]), int(info[1]), int(info[2]), int(info[3]), int(info[4]), int(info[5]))
                if line[0]=="Data Start:\r\n":
                    header=True
            else:
                line=line.split("\n")
                V.append(complex(line[0]))
        f.close()
        V=np.asarray(V)
        #V=V.astype(np.float)
        self.norms=np.abs(V)
        self.data=np.real(V)
        self.N=len(V)
        self.name=fname.split("/")[-1]
        self.gamma=gamma=(2*np.pi)*(3240.) #Hz/G (radians)
        #self.dX=np.power((4.*self.D)/(gamma*self.G),1./3.) #Boundary Thickness (cm)
        print("File "+str(self.name)+" loaded")



    def __str__(self):
        return str(self.name)

    #plots the raw data with FFT bars at times [start,end]. Skips every skip data points
    def show(self,Bars=[False,False],skip=1000):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        ax1.set(title="Raw Signal", xlabel="Time (sec)", ylabel="Voltage (V)")
        time=np.linspace(0,self.N/self.SR,self.N)
        ax1.plot(time,self.norms)
        ax1.grid()
        plt.show()


#Computes the DFT of the Raw Data contained in an FIDdata class
class FIDFFTWindow:
    def __init__(self,RawData,Times=False,PadTime=0,Phases=False):

        self.SR=RawData.SR

        #If not given the start and end time of the DFT, then the DFT is performed on the entire set
        if Times==False:
            self.StartTime=float(0)/RawData.SR
            self.EndTime=float(RawData.N)/RawData.SR
            self.N0=RawData.N
            self.TimeIndexes=[0,RawData.N] #index of time slice
        #If not given the time indexes to be used for the DFT, they are determined from start and end times.
        else:
            self.StartTime=Times[0]
            self.EndTime=Times[1]
            self.TimeIndexes=[int(math.floor(Times[0]*RawData.SR)),int(math.ceil(Times[1]*RawData.SR))]
            self.N0=self.TimeIndexes[1]-self.TimeIndexes[0]

        #Computes appropiate padding length such that the total data set is a power of 2
        padnumber=int(PadTime*RawData.SR)
        if padnumber!=0:
            Nfast=int(math.pow(2, math.ceil(math.log(self.N0+padnumber, 2))))
            self.N=Nfast
        else:self.N=self.N0

        zeros=self.N-self.N0
        zerosBefore=0#zeros//2
        zerosAfter=zeros-zerosBefore


        #Performs the DFT along with mean value zero padding
        self.data=np.fft.rfft(np.lib.pad(RawData.data[self.TimeIndexes[0]:self.TimeIndexes[1]], (zerosBefore,zerosAfter), 'constant', constant_values=(0, 0)))
        #rescales the DFT data to account for padding
        self.Amps=2.0*np.abs(self.data)/self.N0
        self.Freqs=self.GetFreqs(self.N,RawData.SR)

        if Phases==True:
            self.Phases=np.arcsin(np.imag(self.data)/self.Amps)*180/np.pi
            #np.angle(self.data,deg=True)

        #Scaled Frequency, note this does not scale the amplitudes!
        #self.XFreqs=(self.Freqs-RawData.f0)*2.*np.pi/(RawData.G*RawData.gamma) #frequency in units of position with x=0 at f0
        #self.XFreqs=self.XFreqs/RawData.dX #frequency to unitless postion (in boundary thicknesses)

        #to scale amps, use amps->amps*unitlessScale
        #self.unitlessScale=(RawData.G*RawData.gamma*RawData.dX)/2.*np.pi #unitless position/Hz



    def GetIndex(self,Freq):
        return int((Freq*float(self.N))/self.SR)


    #retrieves the frequencies of the DFT
    def GetFreqs(self,N,SR):
        Neff=self.N//2+1
        freqs=[]
        for i in range(Neff):
            freq=i*(SR/N)
            freqs.append(freq)
        return np.asarray(freqs)

#returns all filenames in directory dir, which can be sorted according to the key "sorted" ( Ex: lambda fname:int(fname[0]) )
def getFilenames(dir,sorted=False):
    fnames=[]
    for file in os.listdir(dir):
        if file.endswith(".txt"):
            #print(os.path.join(file))
            fnames.append(dir+os.path.join(file))

    if sorted!=False:
        return sorted(fnames, key=sorted)
    else: return fnames
