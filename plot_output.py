import numpy as np
import matplotlib.pyplot as plt
import h5py
from bmtk.analyzer.cell_vars import _get_cell_report, plot_report
import matplotlib.pyplot as plt
import pandas as pd 
from scipy.signal import find_peaks
import pdb

# Load data
config_file = "simulation_config.json"
raster_file = './output/spikes.h5'

mem_pot_file = './output/v_report.h5'


# load 
f = h5py.File(mem_pot_file,'r')
mem_potential = f['report']['inh_stim']['data']
plt.plot(mem_potential[:,0]+70)

# f = h5py.File('exc_stim_spikes.h5','r')
# input_spikes = f['spikes']['exc_stim']['timestamps'][:]
# plt.plot(input_spikes*10,f['spikes']['exc_stim']['node_ids'][:],'r.')
# plt.show()
#pdb.set_trace()

#df = pd.DataFrame.from_csv("PN_C.csv")


try:
    f = h5py.File(raster_file,'r')
    gids = f['spikes']['stim']['node_ids']
    timestamps = f['spikes']['stim']['timestamps']

    plt.figure()
    plt.plot(timestamps,gids,'.')
except:
    print("No spikes.")

plt.figure()
plt.plot(mem_potential[:,0])


# mem = np.array(mem_potential[:, 0])
# base = mem[5000]
# mem = mem[10000:]
# plt.figure()
# plt.plot(mem)
# subed = mem - base
# plt.figure()
# plt.plot(subed)

# print(np.min(mem[5000:]))
# print(np.min(mem[5000:6000]))
# print(mem[5000])

# print(np.trapz(subed, dx=1))


plt.show()

