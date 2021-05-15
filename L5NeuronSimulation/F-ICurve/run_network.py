from bmtk.simulator import bionet
import numpy as np
from neuron import h
import pandas as pd
import sys

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import my_plotting

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        inp = sys.argv[-1]
    else:
        raise Exception("no work" + str(sys.argv[-1]))

amp = float(inp)

config_file = 'simulation_config.json'

conf = bionet.Config.from_json(config_file, validate=True)
conf.build_env()

graph = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=graph)

sim.run()

try:
    df = pd.read_csv("FI_data.csv")
except:
    print("Making new df")
    df = pd.DataFrame(columns = ["Current", "FR"])

raster_file = './output/spikes.h5'
try:
    #f = h5py.File(raster_file,'r')
    #import pdb; pdb.set_trace()
    #gids = f['spikes']['biophysical']['node_ids']
    timestamps = my_plotting.load_dataset(raster_file)['timestamps']
    #timestamps = f['spikes']['biophysical']['timestamps']
    #f.close()
    #import pdb; pdb.set_trace()
    FR = len(np.where(np.array(timestamps) >= 1000)[0]) / 2
    # plt.figure()
    # plt.plot(timestamps,gids,'.')
    # plt.show()
    # print("Spikes:", len(timestamps))
except:
    print("No spikes.")
    FR = 0

df.loc[len(df.index)] = [amp, FR]  

df.to_csv("FI_data.csv", index=False)