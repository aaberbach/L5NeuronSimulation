import sys, os
from bmtk.simulator import bionet
import numpy as np
import h5py
import pickle
from neuron import h
import pandas as pd

np.random.seed(2129)

#Calculates maximum polarization of the given arr.
#Time time is base level and finds the lowest after that.
def calc_pol(arr, time=4998):
    base = arr[time]
    trough = min(arr[time:])
    return trough - base

# if __name__ == '__main__':
#     if __file__ != sys.argv[-1]:
#         inp = sys.argv[-1]
#     else:
#         raise Exception("no work" + str(sys.argv[-1]))

# dist = float(inp)

#synapses.set_pyr_w(weight)
#synapses.load()

config_file = 'simulation_config.json'

conf = bionet.Config.from_json(config_file, validate=True)
conf.build_env()

graph = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=graph)

cells = graph.get_local_cells()

cell = cells[0]
hobj = cell.hobj

soma = cell.hobj.soma[0](0.5)

soma_v = h.Vector()
soma_v.record(soma._ref_v)

segs = cell._morph.seg_prop
dic = {}
seg_names = []
i = 0
h.distance(sec=cell.hobj.soma[0])
for sec in hobj.all:
    fullsecname = sec.name()
    sec_type = fullsecname.split(".")[1][:4]
    for seg in sec:
        i+= 1
        #if (i % 100 == 0):
        rec = h.Vector()
        rec.record(seg._ref_v)
        dic[str(seg)] = [sec_type, h.distance(seg.x), rec]
        # import pdb; pdb.set_trace()
        # seg_names.append(str(seg))
#import pdb; pdb.set_trace()
# v_ref = cell.connections()[0]._connector.postseg()._ref_v

# rec = h.Vector()
# rec.record(v_ref)

sim.run()
import pdb; pdb.set_trace()

soma_pol = calc_pol(np.array(soma_v))

df = pd.DataFrame()
attenuations = []
distances = []
parts = []
for key in dic.keys():
    parts.append(dic[key][0])
    distances.append(dic[key][1])
    attenuations.append(calc_pol(np.array(dic[key][2])) / soma_pol)

df["type"] = parts
df["distance"] = distances
df["attenuation"] = attenuations

df.to_csv("results.csv", index=False)


# syn_volt = rec.as_numpy()

# mem_pot_file = 'output/v_report.h5'

# # load 
# f = h5py.File(mem_pot_file,'r')
# mem_potential = f['report']['exc_stim']['data']
# import matplotlib.pyplot as plt
# plt.plot(mem_potential[:, 0], label='soma')
# plt.plot(syn_volt, label='synapse')
# plt.legend()
# plt.show()

# low = mem_potential[3998, 0]
# high = max(mem_potential[3998:, 0])
# soma_mag = high - low
# #print(soma_mag)

# low = syn_volt[3998]
# high = max(syn_volt[3998:])
# dend_mag = high - low
# #print(dend_mag)

# attenuation = soma_mag / dend_mag
# print(attenuation)

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

