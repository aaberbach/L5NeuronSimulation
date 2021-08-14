import pandas as pd
import h5py
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.signal as ss

freqs = [8,16,32,64,128]
mods = [2,10,50,100]


#fr = np.zeros((5,4))
#
#for i in np.arange(0,5):
#    for j in np.arange(0,4):   
#        s = h5py.File('./output_{}_{}p/spikes.h5'.format(freqs[i],mods[j]))
#        try:
#            fr[i,j] = s['spikes']['biophysical']['timestamps'][:].shape[0]/10
#        except:
#            fr[i,j] = 0
#
#plt.figure()
#sns.heatmap(fr,annot=True, cbar_kws={'label': 'output firing rate (Hz)'})
#plt.xticks(ticks=np.arange(0,4)+0.5, labels=['2','10','50','100'])
#plt.yticks(ticks=np.arange(0,5)+0.5, labels=['8','16','32','64','128'])
#plt.ylabel('inhibitory sin frequency (Hz)')
#plt.xlabel('modulation (%)')

#plt.figure()
###############################
freq = 16
tstop = 60

tstmps = []
for i in [50]:
    f = h5py.File('./output_{}_{}p/spikes.h5'.format(freq,i))
    try:
        tstmps.append(f['spikes']['biophysical']['timestamps'][:])
    except:
        tstmps.append([0])

t = np.arange(0,tstop,1/(freq*30))

pks,_ = ss.find_peaks(-np.sin(2*np.pi*freq*t))

rel_hist = np.zeros((len(tstmps),29))
for j in np.arange(0,len(tstmps)):
    for p in np.arange(0,pks.shape[0]-1):
        t = np.arange(0,tstop*1000,tstop*1000/t.shape[0])
        hist,_ = np.histogram(tstmps[j],bins=np.linspace(t[pks[p]],t[pks[p+1]],30))
        try:
            rel_hist[j,:] += hist
        except:

            import pdb; pdb.set_trace()

pd.DataFrame(rel_hist).to_csv('rel_hist_{}Hz.csv'.format(freq),index=False)
#rel_hist = pd.read_csv('rel_hist_60Hz.csv').values
plt.figure()
plt.bar(np.arange(0,29),rel_hist[0,:]/np.sum(rel_hist[0,:]))
plt.xticks(ticks=[0,15,30], labels=[r'-$\pi$', '0', r'$\pi$'])
plt.ylabel('p(spike)')

#plt.figure()
#plt.bar(np.arange(0,29),100*rel_hist[1,:]/np.sum(rel_hist[1,:]))
#plt.xticks(ticks=[0,15,30], labels=[r'-$\pi$', '0', r'$\pi$'])
#plt.ylabel('p(spike)')
#
#plt.figure()
#plt.bar(np.arange(0,29),100*rel_hist[2,:]/np.sum(rel_hist[2,:]))
#plt.xticks(ticks=[0,15,30], labels=[r'-$\pi$', '0', r'$\pi$'])
#plt.ylabel('p(spike)')
plt.show()
#plt.title('FR = {} Hz'.format(s['spikes']['biophysical']['timestamps'][:].shape[0]/5))


