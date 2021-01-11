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

def clamp_analysis(inh_file, clamp_file, spike_file, time):
    generate_analysis(inh_file, spike_file, time)

    clamp_wave = avg_clamp_wave(inh_file, clamp_file, time, length=15)

    clamp_wave -= np.min(clamp_wave)
    clamp_wave /= np.max(clamp_wave)

    x = np.arange(len(clamp_wave)) / 10
    x *= (10/15)

    plt.plot(x, clamp_wave, label = "clamp current", color = "green")


#Analyzes how spike times line up with gamma waves.
def generate_analysis(inh_file, spike_file, time):
    gamma = generate_spike_gamma(inh_file, time)

    data = my_plotting.load_dataset(spike_file)
    node_ids = np.array(data['node_ids'])
    timestamps = np.array(data['timestamps'])

    troughs = s.find_peaks(-gamma)[0]

    # lens = [troughs[i+1] - troughs[i] for i in range(len(troughs)-1)]

    # import pdb; pdb.set_trace()

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
    #plt.show()

#Takes a trough-to-trough wave and lines it up to the center of an array of length max_len with padding zeroes.
def lineup_wave(wave, gamma, max_len):
    result = np.zeros(max_len)
    peak = np.argmax(gamma) * 10 + 5

    #print(len(wave), len(gamma))

    if (max_len % 2 == 0):
        raise Exception("max_len must be odd in the function lineup_wave.")

    mid = (max_len + 1) // 2

    for i in range(len(wave)):
        result[mid - (peak - i)] = wave[i]

    # plt.plot(wave)
    # plt.plot(result)
    # plt.plot(gamma)
    # plt.show()

    return result

#Lines up a voltage clamp report with gamma-rhythmic inhibition and creates an average wave.
def avg_clamp_wave(inh_file, clamp_file, time, max_len = 201, length=None):
    gamma = generate_spike_gamma(inh_file, time)
    se_data = np.array(my_plotting.load_dataset(clamp_file, groups=1))

    # plt.plot(gamma)
    # plt.plot(np.arange(len(se_data)) / 10, se_data)
    # plt.show()
    # import pdb; pdb.set_trace()
    troughs = np.array(s.find_peaks(-gamma)[0])
    avg_wave = np.zeros(max_len)

    count = 0
    for i in range(len(troughs) - 1):
        if length != None:
            wave = gamma[troughs[i]: troughs[i+1]]
            if (len(wave) == length and np.argmax(wave) == length//2):
                count += 1
                avg_wave += lineup_wave(se_data[troughs[i]*10 : troughs[i+1]*10], gamma[troughs[i] : troughs[i+1]], max_len)
        else:
            count += 1
            avg_wave += lineup_wave(se_data[troughs[i]*10 : troughs[i+1]*10], gamma[troughs[i] : troughs[i+1]], max_len)
    
    if length != None:
        # mid = ((max_len + 1) // 2) - 1
        # half_len = (length / 2) * 10
        # avg_wave = avg_wave[mid - half_len + 1 : mid + half_len]
        avg_wave = avg_wave[avg_wave != 0]

    return avg_wave / count
    #return avg_wave / (len(troughs) - 1)

#Takes a spike raster and plots a histogram of ms between spikes.
def plot_spike_separation(spike_file):
    data = my_plotting.load_dataset(spike_file)
    timestamps = np.sort(np.array(data['timestamps']))

    sep_times = []
    for i in range(1, len(timestamps)):
        if timestamps[i] - timestamps[i-1] > 40 and timestamps[i] - timestamps[i-1] < 50:
            print(timestamps[i])
        sep_times.append(timestamps[i] - timestamps[i-1])

    plt.figure()
    plt.hist(sep_times, bins = 200)
    plt.xlabel("Time between spikes (ms)")
    plt.ylabel("Frequency")
    plt.show()
        
#Takes a spike raster and plots the distribution on spikes per burst. 
#   max_gap: maximum space between two spikes in a burst.
def plot_burst_counts(spike_file, max_gap):
    data = my_plotting.load_dataset(spike_file)
    timestamps = np.sort(np.array(data['timestamps']))

    dist = {}

    #Creates a dict with [1...10] spikes per burst
    for i in range(1, 10):
        dist[i] = 0

    burst_count = 1
    for i in range(1, len(timestamps)-1):
        #If the count just reset
        if burst_count == 0:
            burst_count = 1
        #If current spike not connected to the previous.
        elif timestamps[i] - timestamps[i-1] > max_gap:
            # if burst_count == 5:
            #     print(timestamps[i-1])
            dist[burst_count] += 1
            burst_count = 0
        else:
            burst_count += 1
    #Adds the final burst to the dict.
    if burst_count > 0:
        dist[burst_count] += 1

    import pdb; pdb.set_trace()
    
    plt.figure()
    plt.bar(list(dist.keys()), list((dist.values())))
    plt.xlabel("Spikes per burst")
    plt.ylabel("Frequency")
    plt.show()

    #import pdb; pdb.set_trace()

my_plotting.plot_v('non_gamma_results/v_report.h5')
#my_plotting.plot_spikes('gamma_results/spikes.h5')
plot_spike_separation('non_gamma_results/spikes.h5')
plot_burst_counts('non_gamma_results/spikes.h5', 30)

import pdb; pdb.set_trace()
avg_wave = avg_clamp_wave('gamma_results/inh_stim_spikes.h5', "non_gamma_results/clamped/se_clamp_report.h5", 120000, length=15)
#import pdb; pdb.set_trace()
plt.plot(avg_wave)
plt.show()

clamp_analysis('gamma_results/inh_stim_spikes.h5', "non_gamma_results/clamped/se_clamp_report.h5",'non_gamma_results/spikes.h5', 120000)
plt.legend()
plt.show()

clamp_analysis('gamma_results/inh_stim_spikes.h5', "gamma_results/clamped/se_clamp_report.h5",'gamma_results/spikes.h5', 120000)
plt.legend()
plt.show()
# avg_wave = avg_clamp_wave('gamma_results/inh_stim_spikes.h5', "gamma_results/clamped/se_clamp_report.h5", 120000, length=15)
# #import pdb; pdb.set_trace()
# plt.plot(avg_wave)
# plt.show()
generate_analysis('gamma_results/inh_stim_spikes.h5', 'gamma_results/spikes.h5', 120000)
generate_analysis('gamma_results/inh_stim_spikes.h5', 'non_gamma_results/spikes.h5', 120000)