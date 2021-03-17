import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

from my_plotting import *
import numpy as np
import matplotlib.pyplot as plt

def zscore(x):
    return (x-np.mean(x))/np.std(x)

def plot_raster_trace(file, seconds, show=True, shift=0):
    data = load_dataset(file)
    #import pdb; pdb.set_trace()
    ts = np.array(data['timestamps'])
    nid = np.array(data['node_ids'])

    h = np.histogram(ts,bins=np.arange(0,seconds*1000,1))
    plt.plot(np.arange(len(h[0])) + shift, zscore(h[0]), label=file)
    
    if show:
        plt.show()

#plot_raster_trace("prox_inh_stim_spikes.h5", seconds=2, show=False, shift=-4)
# plot_raster_trace("dist_inh_stim_spikes.h5", seconds=2, show=False, shift=-4)
# plot_raster_trace("exc_stim_spikes.h5", seconds=2, show=False)
# plt.legend()
# plt.show()

data = load_dataset("output/spikes.h5")
ts = np.array(data['timestamps'])
nid = np.array(data['node_ids'])
print("FR:", len(ts) / 10)

plot_v("shift_output/v_report.h5")
plt.title("With Shift")
plt.figure()
plot_v("output/v_report.h5", show=True)

plot_spikes("exc_stim_spikes.h5", time_scale=1, show=True)
plot_spikes("prox_inh_stim_spikes.h5", time_scale=1, show=True)