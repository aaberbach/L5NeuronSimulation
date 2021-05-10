"""Script for running the network built in build_network.py

Also saves a file called Connections.csv that consists of information about
each synapse in the simulation.
"""

from bmtk.simulator import bionet
import numpy as np
from neuron import h
import synapses 
import pandas as pd

synapses.load()
syn = synapses.syn_params_dicts()

np.random.seed(42)

config_file = 'simulation_config.json'
conf = bionet.Config.from_json(config_file, validate=True)
conf.build_env()

graph = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=graph)

cells = graph.get_local_cells()
cell = cells[list(cells.keys())[0]]

h.distance(sec=cell.hobj.soma[0])#Makes the distances correct.

sec_types = []#soma, apic, or dend
weights = []#scaled conductances (initW)
dists = []#distance from soma
node_ids = []#node id within respective node population (exc, prox_inh, dist_inh)
names = []#full NEURON str representation of postsynaptic segment
source_pops = []#node population
release_probs = []#propability of release.

for c in cell.connections():
    con = c._connector
    source = c.source_node
    syn = con.syn()
    seg = con.postseg()
    fullsecname = seg.sec.name()

    source_pops.append(source._population)
    node_ids.append(source._node_id)

    weights.append(float(syn.initW))
    release_probs.append(float(syn.P_0))
    names.append(str(seg))
    sec_types.append(fullsecname.split(".")[1][:4])
    dists.append(float(h.distance(seg)))

df = pd.DataFrame()
df["Node ID"] = node_ids
df["Distance"] = dists
df["Conductance"] = weights
df["Type"] = sec_types
df["Name"] = names
df["Source Population"] = source_pops
df["Release Probability"] = release_probs
df.to_csv("Connections.csv", index=False)

sim.run()