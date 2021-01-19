import sys, os
from bmtk.simulator import bionet
import numpy as np
import h5py
import synapses
import pickle
from neuron import h

np.random.seed(2129)

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        inp = sys.argv[-1]
    else:
        raise Exception("no work" + str(sys.argv[-1]))

dist = float(inp)

#synapses.set_pyr_w(weight)
synapses.load()

config_file = 'simulation_config.json'

conf = bionet.Config.from_json(config_file, validate=True)
conf.build_env()

graph = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=graph)

cells = graph.get_local_cells()

cell = cells[0]

v_ref = cell.connections()[0]._connector.postseg()._ref_v

rec = h.Vector()
rec.record(v_ref)

sim.run()

syn_volt = rec.as_numpy()

mem_pot_file = 'output/v_report.h5'

# load 
f = h5py.File(mem_pot_file,'r')
mem_potential = f['report']['exc_stim']['data']
import matplotlib.pyplot as plt
plt.plot(mem_potential[:, 0], label='soma')
plt.plot(syn_volt, label='synapse')
plt.legend()
plt.show()

low = mem_potential[3998, 0]
high = max(mem_potential[3998:, 0])
soma_mag = high - low
#print(soma_mag)

low = syn_volt[3998]
high = max(syn_volt[3998:])
dend_mag = high - low
#print(dend_mag)

attenuation = soma_mag / dend_mag
print(attenuation)

# f = open('syn_att.pkl', 'rb')
# res = pickle.load(f)
# if dist in res.keys():
#     res[dist].append(attenuation)
# else:
#     res[dist] = [attenuation]
# f.close()

# #f = open('syn_epsps.pkl', 'wb')
# f = open('syn_att.pkl', 'wb')
# pickle.dump(res, f)
# f.close()

