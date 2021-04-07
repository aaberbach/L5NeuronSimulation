import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import my_plotting
import numpy as np
import matplotlib.pyplot as plt
#import scipy.signal as s
import seaborn as sb
import pandas as pd

#Takes a voltage trace and returns the magnitude of an EPSP
#resulting from an input at time time.
def calc_IPSP(c, time = 4000):
    base = c[time - 2]
    peak = min(c[time - 2:])
    return peak - base

#Takes a file with some number of voltage traces, each with one IPSP.
#Returns an array of IPSP magnitudes.
def calc_IPSPs(file, time = 4000):
    data = np.array(my_plotting.load_dataset(file)['data'])
    ipsps = np.zeros(data.shape[1])

    for i in range(len(ipsps)):
        ipsps[i] = calc_IPSP(data[:, i], time=time)

    return ipsps