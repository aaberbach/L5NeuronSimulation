import h5py
import numpy as np
import pickle
import sys


fname = 'percent2_results.pkl'

# Load data
config_file = "simulation_config.json"
raster_file = './output/spikes.h5'
inp_file = 'stim_spikes.h5'

f = h5py.File(inp_file,'r')
inp_gids = f['spikes']['stim']['node_ids']
inp_timestamps = f['spikes']['stim']['timestamps']

num_exc = len(np.where(inp_gids.value == 0)[0])

try:
    f = h5py.File(raster_file,'r')
    gids = f['spikes']['stim']['node_ids']
    timestamps = f['spikes']['stim']['timestamps']

    num_spikes = len(gids)
except:
    num_spikes = 0

print(num_spikes/num_exc)
if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        inp = sys.argv[-1]
    else:
        raise Exception("no work" + str(sys.argv[-1]))

key = float(inp)

if num_exc == 0:
    res = 0.0
else:
    res = num_spikes / num_exc

file = open(fname, 'rb')
d = pickle.load(file)
d[key] = res
file.close()

file = open(fname, 'wb')
pickle.dump(d, file)
file.close()