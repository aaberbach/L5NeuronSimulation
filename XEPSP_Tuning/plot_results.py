import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import my_plotting
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as s

my_plotting.plot_spikes("exc_stim_spikes.h5", id_scale=10, id_shift = -67, time_scale = 1)
my_plotting.plot_v("output/v_report.h5", show=True, ms=True)