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

sim.run()

import matplotlib.pyplot as plt

import scipy.signal as s

volt = np.array(som)[5000:]
peaks = s.find_peaks(volt)[0]
print("Num Spikes:", len(peaks) - 1)

plt.figure()
plt.plot(np.array(syn))
plt.scatter(peaks + 5000, np.full(len(peaks), -60))
plt.title("Potential at synapse")

plt.figure()
plt.plot(np.array(som))
plt.title("Potential at soma")

plt.show()

pc.barrier()