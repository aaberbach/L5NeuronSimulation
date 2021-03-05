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
syn = conn._syn

# vs = ["r_nmda", "r_ampa", "Capoolcon", "Cdur_nmda", "AlphaTmax_nmda", "Beta_nmda", "Erev_nmda", "gbar_nmda", "W_nmda", "on_nmda", "g_nmda", "Cdur_ampa", "AlphaTmax_ampa", "Beta_ampa", "Erev_ampa", "gbar_ampa", "W_ampa", "on_ampa", "g_ampa", "ECa", "ICa", "P0", "fCa", "tauCa", "iCatotal", "Cainf", "pooldiam", "z", "lambda1", "lambda2", "threshold1", "threshold2", "fmax", "fmin", "Wmax", "Wmin", "maxChange", "normW", "scaleW", "F", "f", "tauF", "D1", "d1", "tauD1", "D2", "d2", "tauD2", "facfactor", "aACH", "bACH", "aDA", "bDA", "wACH", "wDA", "calcium", "initW", "i_nmda", "i_ampa"]
# ripples = [False, False, False, False, False, False, False, False, False, True, True, False, False, False, False, False, False, True, True, False, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, True, True]
# trackers = [h.Vector() for v in vs]
# for i in range(len(trackers)):
#     trackers[i].record(getattr(syn, "_ref_"+vs[i]))

som = h.Vector()
som.record(hobj.soma[0](0.5)._ref_v)

#r = h.Vector()
#r.record(syn._ref_rp)

syn = h.Vector()
syn.record(conn._connector.postseg()._ref_v)



# fac = h.Vector()
# try:
#     fac.record(conn._syn._ref_igaba)
# except:
#     fac.record(conn._syn._ref_facfactor)
#import pdb; pdb.set_trace()
sim.run()

import matplotlib.pyplot as plt
import pdb; pdb.set_trace()
# plt.figure()
# plt.plot(np.array(syn))
# plt.plot(np.array(trackers[-1]) * 10)
# plt.show()

plt.figure()
plt.plot(np.array(syn))
plt.title("Potential at synapse")

plt.figure()
plt.plot(np.arange(0, len(som)) / 10, np.array(som))
plt.title("PN->PN")
plt.ylabel("Voltage at Soma (mV)")
plt.xlabel("time (ms)")

# plt.figure()
# plt.title("Gabba curent FSI->PN")
# plt.xlabel("time (ms)")
# plt.ylabel("Current (nA)")
# plt.plot(np.arange(0, len(fac)) / 10, np.array(fac))

# for i, tracker in enumerate(trackers):
#     plt.figure()
#     plt.title(vs[i])
#     plt.xlim([4200, 6000])
#     plt.plot(np.array(tracker))

plt.show()

pc.barrier()