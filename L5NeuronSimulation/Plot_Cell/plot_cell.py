from bmtk.builder import NetworkBuilder
import numpy as np
import sys
#import synapses   
import pandas as pd                                                                                                             

np.random.seed(2129)

net = NetworkBuilder("biophysical")


##################################################################################
###################################Pyr Type C#####################################

net.add_nodes(N=1, pop_name='Pyrc',
    potental='exc',
    model_type='biophysical',
    model_template='hoc:L5PCtemplate',
    morphology = None)


# Build and save our networks
net.build()
net.save_nodes(output_dir='network')

seconds = 0.01

from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',
                network_dir='./network',
                dt = 0.1, tstop=seconds * 1000.0,
                dL = 5,
                #report_vars=['v', 'cai'],
                # current_clamp={           # Creates a step current from 500.ms to 1500.0 ms  
                #      'amp': 0.793,
                #      #'std': [0.0, 0.0],
                #      'delay': 700,
                #      'duration': 2000,
                #      'gids':"0"
                #  },
                # clamp_reports=['se'],#Records se clamp currents.
                # se_voltage_clamp={
                #         "amps":[[0, 0, 0]],
                #         "durations": [[120000, 0, 0]],
                #         'gids': [0],
                #         'rs': [0.01],
                # },
                spikes_threshold=-10,
                #spikes_inputs=[('exct_stim', 'exc_stim_spikes.h5'), ('inh_stim', 'inh_stim_spikes.h5')],
                components_dir='../biophys_components',
                compile_mechanisms=True)

import sys, os
from bmtk.simulator import bionet
import numpy as np
import pandas as pd
#import h5py
from neuron import h
#from scipy.stats import skew
#import synapses
#from bmtk.simulator.bionet.pyfunction_cache import add_weight_function
#import pickle
#import load_processor

#load_processor.load()
np.random.seed(42)
#synapses.load()
#import pdb; pdb.set_trace()

# pc = h.ParallelContext()  # object to access MPI methods
# MPI_size = int(pc.nhost())
# MPI_rank = int(pc.id())

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
hobj = cell.hobj

import matplotlib.pyplot as plt
#from clustering import *
#import pdb; pdb.set_trace()
ps = h.PlotShape(False)
ax = ps.plot(plt)

df = pd.read_csv("Connections.csv")
#import pdb; pdb.set_trace()
# dends = df[df["Type"] == "dend"]
# apics = df[(df["Type"] == "apic") & (df["Distance"] >= 600)]

def plot_row(row):
    name = row["Name"]
    sec_type = row["Type"]
    #import pdb; pdb.set_trace()
    x = float(name.split("(")[1].split(")")[0])
    sec_id = int(name.split("[")[2].split("]")[0])

    ax.mark(getattr(hobj, sec_type)[sec_id](x))

plot_row(df.iloc[23])

np.random.seed(43)

# plot_group(apics, FunctionalGroup(apics, apics.sample().iloc[0], 8, 8, "first", partial(make_seg_sphere, radius = 100), partial(make_seg_sphere, radius = 5)), ax, psoma=psoma)
# plot_group(dends, FunctionalGroup(dends, dends.sample().iloc[0], 8, 8, "first", partial(make_seg_sphere, radius = 100), partial(make_seg_sphere, radius = 5)), ax, psoma=psoma)
plt.show()