from bmtk.simulator import bionet
import numpy as np
from neuron import h
import synapses
synapses.load()
syn = synapses.syn_params_dicts()

pc = h.ParallelContext()  # object to access MPI methods
MPI_size = int(pc.nhost())
MPI_rank = int(pc.id())

config_file = 'simulation_config.json'

conf = bionet.Config.from_json(config_file, validate=True)
conf.build_env()

graph = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=graph)

cells = graph.get_local_cells()
cell = cells[list(cells.keys())[0]]

hobj = cell.hobj
conn = cell.connections()[0]

som = h.Vector()
som.record(hobj.soma[0](0.5)._ref_v)

syn = h.Vector()
syn.record(conn._connector.postseg()._ref_v)

fac = h.Vector()
try:
    fac.record(conn._syn._ref_igaba)
except:
    fac.record(conn._syn._ref_facfactor)
import pdb; pdb.set_trace()
sim.run()

import matplotlib.pyplot as plt
plt.figure()
plt.plot(np.array(syn))
plt.title("Potential at synapse")

plt.figure()
plt.plot(np.array(som))
plt.title("Potential at soma")

plt.figure()
plt.plot(np.array(fac))

plt.show()

pc.barrier()