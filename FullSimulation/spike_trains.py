import matplotlib.pyplot as plt
import h5py
import numpy as np

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

from my_plotting import *

def bin_raster(time, file):
    hist = np.full(0, time)
    data = load_dataset(file)
    ts = data['timestamps']
    #import pdb; pdb.set_trace()
    for t in ts:
        hist[np.trunc(t)] += 1

    return hist

file = "exc_stim_spikes.h5"
hist = bin_raster(2000, file)
import pdb; pdb.set_trace()