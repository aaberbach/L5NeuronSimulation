import numpy as np
import scipy.signal as ss
import scipy
import matplotlib.pyplot as plt
import h5py
from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator
import pandas as pd
from scipy.fftpack import fft

def zscore(x):
    return (x-np.mean(x))/np.std(x)

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

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
psds = np.zeros((1,129))
for i in np.arange(0,1):
    psds[i,:], psd_frequencies = psd((z[i,:]-np.mean(z[i,:]))/np.std(z[i,:]), NFFT=256, Fs=1/dt)
plt.loglog(psd_frequencies, np.mean(psds,0))


################3
def make_spikes(numUnits=10000,rateProfMat=z):
    rateProf=z[0,:]/1000

    meanRates = np.mean(rateProf)
    normRateProf= rateProf/meanRates

    df_full = pd.DataFrame({'timestamps':[],'node_ids':[]})

    for i in np.arange(0,numUnits):
        rate_temp=[];simSpks_temp=[];

        rate_temp = normRateProf*meanRates
        numbPoints = scipy.stats.poisson(rate_temp).rvs()#Poisson number of points

        simSpks=np.where(numbPoints>0)[0]

        timestamps=simSpks
        node_ids = np.tile(i,simSpks.shape[0])
        df = pd.DataFrame({'timestamps':timestamps, 'node_ids':node_ids})
        df_full = df_full.append(df)

    return df_full

df_full = make_spikes(numUnits=10000,rateProfMat=z)
#df_full2 = make_spikes(numUnits=1000,rateProfMat=z)
df_full3 = make_spikes(numUnits=100,rateProfMat=z)

ma_win = 50

x1 = np.histogram(df_full['timestamps'],bins=np.arange(0,20000,1))
y1 =  moving_average(x1[0],ma_win)

#x2 = np.histogram(df_full2['timestamps'],bins=np.arange(0,20000,1))
#y2 =  moving_average(x2[0],ma_win)

x3 = np.histogram(df_full3['timestamps'],bins=np.arange(0,20000,1))
y3 =  moving_average(x3[0],100)

plt.figure()


plt.plot(x1[1][:-ma_win],zscore(moving_average(z[0,:],ma_win)),label='input',color='k')
plt.plot(x1[1][:-ma_win],zscore(y1*(1000/10000)),color='b',linestyle='-',label='output FR (10000 units)')
#plt.plot(x2[1][:-ma_win],zscore(y2*(1000/1000)),color='b',linestyle='-',label='output FR (1000 units)')
plt.plot(x3[1][:-100],zscore(y3*(1000/100)),color='r',linestyle='-',label='output FR (100 units)')
plt.legend()
plt.ylabel('zscore firing rate')
plt.xlabel('time (ms)')
plt.xlim(5000,7000)

plt.show()

plt.figure()
plt.scatter(y*(1000/numUnits),moving_average(rateProf*1000,100),alpha=0.01)
plt.xlabel('actual FR (Hz)')
plt.ylabel('expected FR (Hz)')

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
