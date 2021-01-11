import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import my_plotting
import numpy as np

def compare_h5(file1, file2):
    data1 = my_plotting.load_dataset(file1)
    data2 = my_plotting.load_dataset(file2)

    import pdb; pdb.set_trace()

compare_h5('non_gamma_results/inh_stim_spikes.h5', 'inh_stim_spikes.h5')