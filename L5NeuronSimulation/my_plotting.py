"""My Plotting Tool

Tool used for plotting various things throughout the project.
Helpful for navigating the BMTK h5 structuring.
"""

import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy.signal as s

def get_key(group, index=0):
    """From the list of the keys of the given group, returns the key at
    the given index.

    Parameters
    ----------
    group : h5py.Group
        the h5py group to get the key from
    index : int, optional
        index of the key in the list of the group's keys, by default 0

    Returns
    -------
    str
        desired key for the h5py group
    """
    return list(group.keys())[index]

def load_dataset(fname, groups=2):
    """Gets a dataset within the given h5 file.
    Many BMTK h5 files have one dataset within some layers of group,
    and this is a useful function for getting to that dataset.
    Assumes that each group just has one key.

    Parameters
    ----------
    fname : str
        h5 file to load
    groups : int, optional
        number of groups before the dataset, by default 2

    Returns
    -------
    h5py.Dataset
        the desired dataset
    """    
    f = h5py.File(fname, 'r')
    for i in range(groups):
        f = f[get_key(f)]
    return f

def plot_spikes(file, show=False, id_scale=-1, id_shift = 0, time_scale = 10):
    data = load_dataset(file)

    scale = 1
    if id_scale > 0:
        scale = id_scale / np.max(data['node_ids'])

    plt.plot(np.array(data['timestamps'])*time_scale,np.array(data['node_ids']) * scale + id_shift,'.')
    
    if(show):
        plt.show()

def plot_v(file, show=False, ms=False):
    """Plots the membrane potential from the given BMTK v_report.h5 file.

    Parameters
    ----------
    file : str
        location of the h5py file
    show : bool, optional
        whether to call plt.show() at the end, by default False
    ms : bool, optional
        whether to scale x by 0.1 to get ms scale, by default False
    """    
    data = load_dataset(file)

    x = np.arange(0, np.array(data['data']).shape[0])
    if ms:
        x = x / 10

    plt.plot(x, data['data'][:, 0])

    if(show):
        plt.show()

def plot_all_v(file, ms=False):
    """Plots each membrane potential in the given BMTK v_report.h5 file.

    Parameters
    ----------
    file : str
        location of the h5py file
    ms : bool, optional
        whether to scale x by 0.1 to get ms scale, by default False
    """    
    data = load_dataset(file)

    x = np.arange(0, np.array(data['data']).shape[0])
    if ms:
        x = x / 10

    for i in range(data['data'].shape[1]):
        plt.plot(x, data['data'][:, i])
        plt.show()

def plot_se(file, show=False):
    """Used to plot se_clamp_reports from BMTK.

    Parameters
    ----------
    file : str
        location of the h5py file
    show : bool, optional
        whether to call plt.show() at the end, by default False
    """    
    data = load_dataset(file, groups=1)
    plt.plot(data[:, 0])

    if(show):
        plt.show()

# def generate_spike_probs(inh_file, spike_file, time):
#     gamma = generate_spike_gamma(inh_file, time)

#     data = load_dataset(spike_file)
#     timestamps = np.array(data['timestamps'])

#     troughs = s.find_peaks(-gamma)[0]

#     n_parts = 10
#     parts = np.zeros(n_parts)
#     for i in range(len(troughs) - 1):
#         start = troughs[i]
#         part_len = (troughs[i+1] - start)/n_parts
#         for j in range(n_parts):
#             parts[j] += len(np.where((timestamps >= j*part_len + start) & (timestamps < (j+1)*part_len + start))[0])

#     parts = np.array(parts) / parts.sum()
#     #t1 = gamma[troughs[100]:troughs[101]]
#     t1 = gamma[troughs[0]:troughs[1]]
#     plt.plot(parts, label="spike probability")
#     plt.plot(np.arange(len(t1)) * (n_parts/len(t1)), t1/10, label="gamma ex.")
#     plt.legend()
#     plt.show()
#     return parts

# def generate_prob_raster(inh_file, spike_file, time):
#     gamma = generate_spike_gamma(inh_file, time)

#     data = load_dataset(spike_file)
#     node_ids = np.array(data['node_ids'])
#     timestamps = np.array(data['timestamps'])

#     troughs = s.find_peaks(-gamma)[0]

#     new_ts = np.zeros(len(timestamps))
#     #ids = np.arange(len(timestamps))
#     cycle_num = np.zeros(len(new_ts))
#     for i in range(len(troughs) - 1):
#         start = troughs[i]
#         stop = troughs[i + 1]
#         length = stop - start

#         spikes = np.where((timestamps >= start) & (timestamps < stop))[0]
#         cycle_num[spikes] = i
#         times = timestamps[spikes]
#         times = times - start
#         times = times / length
#         new_ts[spikes] = times
#         #part_len = (troughs[i+1] - start)/n_parts
#         # for j in range(n_parts):
#         #     parts[j] += len(np.where((timestamps >= j*part_len + start) & (timestamps < (j+1)*part_len + start))[0])

#     parts = np.zeros(10)
#     sep = 0.1
#     for i in range(10):
#         parts[i] = len(np.where((new_ts >= (i * sep)) & (new_ts < ((i+1)*sep)))[0])
#     #parts = np.array(parts) / parts.sum()
#     #t1 = gamma[troughs[100]:troughs[101]]
#     t1 = gamma[troughs[0]:troughs[1]]
#     #import pdb; pdb.set_trace()
#     #plt.plot(parts, label="spike probability")
#     plt.plot(np.arange(10)+0.5, (parts / parts.sum()), color="black", label="spike probability")
#     plt.plot(new_ts*(10), cycle_num/max(cycle_num), ".")
#     plt.xticks([0, 5, 10], labels = ["-" + r'$\pi$', 0, r'$\pi$'])
#     plt.axvline(x=5, ls="--", color = "black")
#     #plt.plot(np.arange(len(t1)) * (len(t1)/len(t1)), t1/3, label="gamma ex.")
#     #plt.plot(t1/3, label="gamma ex.")
#     plt.legend()
#     ax = plt.gca()
#     #ax.axes.xaxis.set_visible(False)
#     ax.axes.yaxis.set_visible(False)
#     plt.show()
#     #return parts

# def plot_spike_gamma(file, time):
#     gamma = generate_spike_gamma(file, time)
#     # troughs = s.find_peaks(-smooth)[0]

#     # parts = np.zeros(10)
#     # for i in range(len(troughs) - 1):
#     #     part_len = (troughs[i+1] - troughs[i])/10
#     #     for j in range(10):
#     #         parts[j] += len(np.where((timestamps >= j*part_len) & (timestamps < (j+1)*part_len))[0])

#     #import pdb; pdb.set_trace()
#     plt.plot(np.arange(time)*10, smooth)
