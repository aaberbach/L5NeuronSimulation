import matplotlib.pyplot as plt
import numpy as np
import h5py

def get_key(group, index=0):
    return list(group.keys())[index]

def load_dataset(file, groups=2):
    f = h5py.File(file, 'r')
    f = f[get_key(f)]
    f = f[get_key(f)]
    return f

def plot_spikes(file, show=False):
    data = load_dataset(file)
    plt.plot(np.array(data['timestamps'])*10,data['node_ids'],'.')
    
    if(show):
        plt.show()

def plot_v(file, show=False):
    data = load_dataset(file)
    plt.plot(data['data'][:, 0])

    if(show):
        plt.show()
