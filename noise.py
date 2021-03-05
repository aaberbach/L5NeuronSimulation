import numpy as np
import scipy.signal as ss
import scipy
import matplotlib.pyplot as plt
import h5py
from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator
import pandas as pd
from scipy.fft import fft

def make_noise(num_traces=100,num_samples=4999, mean_fr = 1, std_fr = 1):
    num_samples = num_samples+2000
    fv = np.linspace(0, 1, 40);                                # Normalised Frequencies
    #a = 1/(1+2*fv);                                      # Amplitudes Of '1/f'
    #b = ss.firls(43, fv, a);                                   # Filter Numerator Coefficients
    B = [0.049922035, -0.095993537, 0.050612699, -0.004408786];
    A = [1, -2.494956002,   2.017265875,  -0.522189400];
    invfn = np.zeros((num_traces,num_samples))
    for i in np.arange(0,num_traces):
        wn = np.random.normal(loc=mean_fr,
 			   scale=std_fr,size=num_samples)
        invfn[i,:] = ss.lfilter(B, A, wn);                             # Create '1/f' Noise
    return invfn[:,2000:]

x = make_noise(num_samples=20000-1,mean_fr=1,num_traces=100)

z = make_noise(num_samples=20000-1,mean_fr=1,num_traces=100)

plt.figure()
plt.subplot(3,1,1)
plt.plot(z[0,:])
plt.title('firing rate trace 1')

plt.subplot(3,1,2)
plt.hist(z[0,:],bins=100)

plt.subplot(3,1,3)
dt = 1/1000
from matplotlib.mlab import psd
psds = np.zeros((z.shape[0],129))
for i in np.arange(0,z.shape[0]):
    psds[i,:], psd_frequencies = psd((z[i,:]-np.mean(z[i,:]))/np.std(z[i,:]), NFFT=256, Fs=1/dt)
plt.loglog(psd_frequencies, np.mean(psds,0))


################3
numUnits=100
rateProf=x[0,:]/1000

meanRates = np.mean(rateProf)

normRateProf= rateProf/meanRates
 
selProfs = np.ceil(np.random.rand(numUnits,1))

selRates = np.ceil(np.random.rand(numUnits,1))

spk_dict=[]

df_full = pd.DataFrame({'timestamps':[],'node_ids':[]})

for i in np.arange(0,numUnits):
    rate_temp=[];simSpks_temp=[];

    #rate_temp=normRateProf[selProfs[i],:]*meanRates[selRates[i]]
    rate_temp = normRateProf*meanRates
    numbPoints = scipy.stats.poisson(rate_temp).rvs()#Poisson number of points

    simSpks=np.where(numbPoints>0)[0]

    timestamps=simSpks
    node_ids = np.tile(i,simSpks.shape[0])
    df = pd.DataFrame({'timestamps':timestamps, 'node_ids':node_ids})
    df_full = df_full.append(df)

plt.figure()
plt.plot(df_full['timestamps'], df_full['node_ids'], '.')

#spk_dict = np.array([i for i in spk_dict], dtype=object)
#plt.figure()
#x = np.histogram(df_full['timestamps'],bins=np.arange(0,5000,1))
#y =  x[0]/10
#plt.plot(x[1][1:],(y-np.mean(y))/np.std(y),label='actual')
#plt.plot((rateProf-np.mean(rateProf))/np.std(rateProf),label='input')
#plt.legend()

#plt.figure()
#plt.plot((y-np.mean(y))/np.std(y), (rateProf-np.mean(rateProf))/np.std(rateProf),'.')

##### CROSS CORRELATION #####
#def spkxcorr(spkref,spktgt,binwidth,window):
#    binwidth=1
#    window=100
#
#    h = np.zeros((1,window))
#    for i in spkref:
#        h += np.array(np.histogram(spktgt,bins=np.arange(i-int(window/2),i+int(window/2)+1,binwidth))[0])
#    return h
#
#spktgt = df_full.loc[df_full.node_ids==1,'timestamps']
#spkref = df_full.loc[df_full.node_ids==0,'timestamps']
#
#xcorrs = np.zeros((1000,100))
#for i in np.arange(1,1001):
#    xcorrs[i-1,:] = spkxcorr(spkref,spktgt,binwidth=1,window=100)
#
#
#
#
#plt.figure()
#plt.bar(np.arange(0,100),np.mean(xcorrs,0))



plt.show()
