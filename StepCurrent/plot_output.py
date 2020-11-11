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

f = h5py.File(mem_pot_file,'r')
mem_potential = f['report']['biophysical']['data']

try:
    f = h5py.File(raster_file,'r')
    #import pdb; pdb.set_trace()
    gids = f['spikes']['biophysical']['node_ids']
    timestamps = f['spikes']['biophysical']['timestamps']

    plt.figure()
    plt.plot(timestamps,gids,'.')
    plt.show()
    print("Spikes:", len(timestamps))
except:
    print("No spikes.")

plt.figure()
plt.plot(mem_potential[:,0], label="soma")
# plt.plot(mem_potential[:, 10], label="dend")
# plt.plot(mem_potential[:, -20], label="apic")
# plt.plot(mem_potential[:, -1], label="axon")
plt.legend()

plt.show()

