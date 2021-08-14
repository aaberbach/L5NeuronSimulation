import pandas as pd
import numpy as np
import h5py


def spkxcorr(spkref, spktgt, binwidth, window):
    spkref = np.sort(np.array(spkref).astype(int))
    spktgt = np.sort(np.array(spktgt).astype(int))

    h = np.zeros((int(window/binwidth)+1,))
    for i in spkref:
        l = i-int(window/2)-binwidth/2
        r = i+int(window/2)+binwidth/2
        h += np.array(np.histogram(spktgt[(spktgt>=1) & (spktgt<=r)], bins=np.arange(l,r+1,binwidth))[0])
    return h

f = h5py.File('exc_stim_spikes.h5','r')

xc = np.zeros((f['spikes']['exc_stim']['node_ids'][:].max(),41))

for i in np.arange(0,3):#f['spikes']['exc_stim']['node_ids'][:].max()+1):
    precell = f['spikes']['exc_stim']['timestamps'][:][f['spikes']['exc_stim']['node_ids'][:]==0]
    postcell = f['spikes']['exc_stim']['timestamps'][:][f['spikes']['exc_stim']['node_ids'][:]==i]
    print(i)
    xc[i,:] = spkxcorr(spkref=postcell, spktgt=precell, binwidth=1, window=40)

pd.DataFrame(xc).to_csv('cell0xcorr.csv')
