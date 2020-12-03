import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import my_plotting
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as s

#Bins spiking by ms and filters for a smooth gamma oscillation.
def generate_spike_gamma(file, time):
    data = my_plotting.load_dataset(file)
    timestamps = np.array(data['timestamps'])

    bins = [len(np.where(np.trunc(timestamps) == i)[0]) for i in range(time)]

    b, a = s.butter(2, [50, 80], btype = "bandpass",fs = 1000)

    gamma = s.filtfilt(b, a, bins)

    return gamma

#Analyzes how spike times line up with gamma waves.
def generate_analysis(inh_file, spike_file, time):
    gamma = generate_spike_gamma(inh_file, time)

    data = my_plotting.load_dataset(spike_file)
    node_ids = np.array(data['node_ids'])
    timestamps = np.array(data['timestamps'])

    troughs = s.find_peaks(-gamma)[0]

    new_ts = np.zeros(len(timestamps))#relative time stamps
    cycle_num = np.zeros(len(new_ts))
    for i in range(len(troughs) - 1):
        start = troughs[i]
        stop = troughs[i + 1]
        length = stop - start

        spikes = np.where((timestamps >= start) & (timestamps < stop))[0]
        cycle_num[spikes] = i
        times = timestamps[spikes]
        times = times - start
        times = times / length
        new_ts[spikes] = times

    bins = np.zeros(10)
    sep = 0.1#proportion of wave per bin.
    for i in range(int(1/sep)):
        bins[i] = len(np.where((new_ts >= (i * sep)) & (new_ts < ((i+1)*sep)))[0])
    
    proportions = bins/bins.sum()#Converts each bin to a proportion of all spikes instead of a total.

    plt.plot(np.arange(10)+0.5, proportions, color="black", label="AP probability")
    plt.plot(new_ts*(10), cycle_num/max(cycle_num), ".")
    plt.xticks([0, 5, 10], labels = ["-" + r'$\pi$', 0, r'$\pi$'])
    plt.axvline(x=5, ls="--", color = "black")
    plt.legend()
    ax = plt.gca()
    ax.axes.yaxis.set_visible(False)
    plt.show()

generate_analysis('inh_stim_spikes.h5', 'results/spikes.h5', 120000)