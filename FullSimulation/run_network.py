from bmtk.simulator import bionet
import numpy as np
from neuron import h
import synapses 
import pandas as pd

synapses.load()
syn = synapses.syn_params_dicts()

np.random.seed(42)

# pc = h.ParallelContext()  # object to access MPI methods
# MPI_size = int(pc.nhost())
# MPI_rank = int(pc.id())

config_file = 'simulation_config.json'
try:
    conf = bionet.Config.from_json("config.json", validate=True)
except:
    conf = bionet.Config.from_json(config_file, validate=True)

conf.build_env()

graph = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=graph)

cells = graph.get_local_cells()
cell = cells[list(cells.keys())[0]]

h.distance(sec=cell.hobj.soma[0])

sec_types = []
weights = []
dists = []
node_ids = []
names = []
source_pops = []

for c in cell.connections():
    con = c._connector
    source = c.source_node
    syn = con.syn()
    seg = con.postseg()
    fullsecname = seg.sec.name()

    source_pops.append(source._population)
    node_ids.append(source._node_id)

    weights.append(float(syn.initW))
    names.append(fullsecname)
    sec_types.append(fullsecname.split(".")[1][:4])
    dists.append(float(h.distance(seg)))

df = pd.DataFrame()
df["Node ID"] = node_ids
df["Distance"] = dists
df["Conductance"] = weights
df["Type"] = sec_types
df["Name"] = names
df["Source Population"] = source_pops
#df["Node ID"] = node_ids
#import pdb; pdb.set_trace()
# df.to_csv("EPSPs.csv", index=False)
# import pdb; pdb.set_trace()
df.to_csv("Connections.csv", index=False)
# def integerize(val):
#     return int(val)

# v_mod = sim._sim_mods[0]

# v_mod._transforms["v"] = integerize
# v_mod._variables = []

import pdb; pdb.set_trace()
sim.run()
# pc.barrier()