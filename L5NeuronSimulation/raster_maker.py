"""
Contains the functions and class (SonataWriter) necessary for generating and saving the input spike rasters.
"""

import numpy as np
import scipy.signal as ss
import scipy
import scipy.stats as st
import matplotlib.pyplot as plt
import h5py
from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator
import pandas as pd
from scipy.fft import fft
import matplotlib
import statsmodels.api as sm


class SonataWriter:
    """Class used to dynamically writing spike rasters to an h5 file.

    Attributes
    ----------
    file : h5py.File
        file object being worked on
    group : h5py.Group
        gropu where the datasets reside
    datasets : dict
        datasets that are saved to the file

    Methods
    -------
    append_ds(vals, ds)
        appends the given values to the end of the given dataset
    append_repeat(ds, val, N)
        appends the given value N times to the end of the given dataset
    close()
        close the h5py file
    """
    def __init__(self, f_name, groups, datasets, types):
        """
        Parameters
        ----------
        f_name : str
            name of file location
        groups : list
            list of group names (str) that are layered into the h5py file
            in the order given.
        datasets : list
            list of dataset names (str)
        types : list
            list of data types that corresponds to the datasets list
        """        
        self.file = h5py.File(f_name, 'w')

        self.group = self.file
        for group in groups:
            self.group = self.group.create_group(group)

        self.datasets = {}
        for i, ds in enumerate(datasets):
            self.datasets[ds] = self.group.create_dataset(ds, data=[], dtype=types[i], chunks=True, maxshape=(None,))

    def append_ds(self, vals, ds):
        """appends the given values to the end of the given dataset

        Parameters
        ----------
        vals : list
            list of values to be appended to the dataset
        ds : str
            key of the dataset to append to
        """        
        length = len(self.datasets[ds])
        self.datasets[ds].resize((length + len(vals), ))
        self.datasets[ds][length:] = vals

    def append_repeat(self, ds, val, N):
        """appends the given value N times to the end of the given dataset

        Parameters
        ----------
        ds : str
            key of the dataset to append to
        val : [type]
            value to be appended N times
        N : int
            number of vals to append to the dataset
        """        
        self.append_ds([val for i in range(N)], ds)

    def close(self):
        """Closes the h5py File
        """        
        self.file.close()

def zscore(x):
    """z scores the given array"""
    return (x-np.mean(x))/np.std(x)

def minmax(x):
    """min max normalizes the given array"""
    return (x - np.min(x))/(np.max(x)-np.min(x))

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def make_noise(num_traces=100,num_samples=4999):
    """Creates a noise trace used in generating spike rasters.

    Parameters
    ----------
    num_traces : int, optional
        number of noise traces to create (first dimension), by default 100
    num_samples : int, optional
        length of the trace (second dimension), by default 4999

    Returns
    -------
    np.array
        noise trace
    """    
    B = [0.049922035, -0.095993537, 0.050612699, -0.004408786]
    A = [1, -2.494956002,   2.017265875,  -0.522189400]
    invfn = np.zeros((num_traces,num_samples))
    for i in np.arange(0,num_traces):
        wn = np.random.normal(loc=1,
 			   scale=0.5,size=num_samples+2000)
        invfn[i,:] = minmax(ss.lfilter(B, A, wn)[2000:])+0.5                             # Create '1/f' Noise
    return invfn


def shift_exc_noise(ts, nid, seconds, time_shift=4):
    """Creates a shifted, min-max normalized average traces of the given spike raster.

    Parameters
    ----------
    ts : list
        times (float) where spikes occur
    nid : int
        node id associated with each spike
    seconds : float
        length of the raster in seconds
    time_shift : int, optional
        how many ms to shift the average trace by, by default 4

    Returns
    -------
    [type]
        [description]
    """    
    h = np.histogram(ts,bins=np.arange(0,seconds*1000,1))

    fr_prof = h[0]/(0.001*(np.max(nid)+1))
    wrap = fr_prof[-4:]
    fr_prof[4:] = fr_prof[0:-4]
    fr_prof[0:4] = wrap

    fr_prof = minmax(fr_prof)+0.5
    return fr_prof

def make_save_spikes(writer, exp, dist, numUnits=100,rateProf=None,start_id=0,start_time=0):
    """Creates and saves spikes for the given nodes using
    the provided noise trace and a random mean firing rate generated with
    the given distribution.

    Parameters
    ----------
    writer : SonataWriter
        how the spikes are saved
    exp : bool
        whether the value from dist should be fed to np.exp()
    dist : np.array()
        array of firing rates of shape (numUnits,)
    numUnits : int, optional
        number of nodes to generate spikes for, by default 100
    rateProf : np.array(), optional
        noise trace for each unit must have numUnits rows, by default None
    start_id : int, optional
        node_id that the first unit/node should be associated with, by default 0
    start_time : int, optional
        at what time the spikes should start being generated, by default 0
    """    
    
    for i in np.arange(0,numUnits):
       
        try: 
            r = rateProf[i,:]
        except:
            import pdb; pdb.set_trace()
        r[r<0] = 0#Can't have negative firing rates.
        rate_temp=[];simSpks_temp=[]

        #Multiplies the noise trace by the randomly generated firing rate.
        if exp:
            rate_temp = r*np.exp(dist[i])
        else:
            rate_temp = r*dist[i]

        numbPoints = scipy.stats.poisson(rate_temp/1000).rvs()#Poisson number of points

        simSpks=np.where(numbPoints>0)[0]

        writer.append_repeat("node_ids", i + start_id, len(simSpks))
        writer.append_ds(simSpks + start_time, "timestamps")
