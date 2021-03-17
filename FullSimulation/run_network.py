from bmtk.simulator import bionet
import numpy as np
from neuron import h
import synapses 

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

# def integerize(val):
#     return int(val)

# v_mod = sim._sim_mods[0]

# v_mod._transforms["v"] = integerize
# v_mod._variables = []

#import pdb; pdb.set_trace()
sim.run()
# pc.barrier()