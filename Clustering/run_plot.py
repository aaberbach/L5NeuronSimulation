import sys, os
from bmtk.simulator import bionet
import numpy as np
import pandas as pd
#import h5py
from neuron import h
#from scipy.stats import skew
import synapses
#from bmtk.simulator.bionet.pyfunction_cache import add_weight_function
#import pickle
#import load_processor

#load_processor.load()
np.random.seed(42)
#synapses.load()
#import pdb; pdb.set_trace()
synapses.load()
syn = synapses.syn_params_dicts()

pc = h.ParallelContext()  # object to access MPI methods
MPI_size = int(pc.nhost())
MPI_rank = int(pc.id())

# if __name__ == '__main__':
#     if __file__ != sys.argv[-1]:
#         inp = sys.argv[-1]
#     else:
#         raise Exception("no work" + str(sys.argv[-1]))

# fname = str(inp)

config_file = 'simulation_config.json'

conf = bionet.Config.from_json(config_file, validate=True)
conf.build_env()

graph = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=graph)

cells = graph.get_local_cells()

cell = cells[list(cells.keys())[0]]

psoma = cell.morphology.psoma

import matplotlib.pyplot as plt

ps = h.PlotShape(False)
ax = ps.plot(plt)
#import pdb; pdb.set_trace()
for con in cell._connections:
    #print(con._connector.postseg())
    #print(p.mark)
    #import pdb; pdb.set_trace()
    ax.mark(con._connector.postseg())

import pickle
f = open('groups.pkl', 'rb')
groups = pickle.load(f)
f.close()

from clustering import *

df = pd.read_csv("Segments.csv")
dends = df[df["Type"] == "dend"]
apics = df[(df["Type"] == "apic")]

for key in groups.keys():
    segs = df[df["Type"] == key]
    for group in groups[key]:
        plot_group(segs, group, ax, psoma=psoma)

plt.show()