import numpy as np
import scipy.signal as ss
import matplotlib.pyplot as plt
import h5py
from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator
import pandas as pd
import scipy

def make_noise(num_samples=4999, mean_fr = 1, std_fr = 1):
    fv = np.linspace(0, 1, 20);                                # Normalised Frequencies
    a = 1/(1 + fv*2);                                      # Amplitudes Of '1/f'
    b = ss.firls(43, fv, a);                                   # Filter Numerator Coefficients
    wn = np.random.uniform(low=mean_fr-std_fr,
 			   high=mean_fr+std_fr,size=num_samples)
    invfn = ss.filtfilt(b, 1, wn);                             # Create '1/f' Noise
    return invfn

x = make_noise()


#function simSpks = GenerateSimulatedEnsemble(rateProf, numUnits)

 

#  % Generates simulated spike trains by randomly recombining the rate

#  % profile of units with their mean rates and then driving a poisson

#  % process.

 

#  % rateProf - an NxT matrix, the rate profile for each of N units across T

#  % time steps.

#  % numUnits - a scalar, the number of units you want to simulate.

numUnits=10000
rateProf=x/1000

meanRates = np.mean(rateProf)

 

#  %normRateProf = rateProf./meanRates;

normRateProf= rateProf/meanRates

 

selProfs = np.ceil(np.random.rand(numUnits,1))

selRates = np.ceil(np.random.rand(numUnits,1))

spk_dict=[]

df_full = pd.DataFrame({'timestamps':[],'node_ids':[]})

plt.figure()
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


plt.plot(df_full['timestamps'], df_full['node_ids'], '.')

#spk_dict = np.array([i for i in spk_dict], dtype=object)
plt.figure()
x = np.histogram(df_full['timestamps'],bins=np.arange(0,5000,1))
y =  x[0]/10
plt.plot(x[1][1:],(y-np.mean(y))/np.std(y),label='actual')
plt.plot((rateProf-np.mean(rateProf))/np.std(rateProf),label='input')
plt.legend()

plt.figure()
plt.plot((y-np.mean(y))/np.std(y), (rateProf-np.mean(rateProf))/np.std(rateProf),'.')

plt.show()
