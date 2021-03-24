import sys, os
from bmtk.simulator import bionet
import numpy as np
import h5py
import synapses
import pandas as pd
from neuron import h

np.random.seed(2129)
#np.random.seed(42)

if __name__ == '__main__':
    if __file__ != sys.argv[-1] and __file__ != sys.argv[-2]:
        m = float(sys.argv[-2])
        s = float(sys.argv[-1])
    else:
        raise Exception("The inputted command argument not functional: " + str(sys.argv[-1]) + " " + str(sys.argv[-2]))

synapses.set_pyr_w(m, s)
synapses.load()

config_file = 'simulation_config.json'

conf = bionet.Config.from_json(config_file, validate=True)
conf.build_env()

graph = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=graph)
#import pdb; pdb.set_trace()

#segs = list(graph.get_local_cells().values())[0]._morph.seg_prop
#import pdb; pdb.set_trace()
cells = graph.get_local_cells()


sec_types = []
weights = []
dists = []
node_ids = []
names = []
for cell in cells.values():
    node_ids.append(cell.node_id)
    #import pdb; pdb.set_trace()
    h.distance(sec=cell.hobj.soma[0])

    ws = []
    ns = []
    ds = []

    st = None

    for c in cell.connections():
        con = c._connector
        syn = con.syn()

        ws.append(float(syn.initW))
        seg = con.postseg()

        fullsecname = seg.sec.name()
        ns.append(fullsecname)

        if st == None:
            st = fullsecname.split(".")[1][:4]

        ds.append(float(h.distance(seg)))

    weights.append(ws)
    names.append(ns)
    dists.append(ds)

    sec_types.append(st)
    # con = cell.connections()[0]._connector
    # syn = con.syn()
    # weights.append(float(syn.initW))

    # seg = con.postseg()
    # fullsecname = seg.sec.name()
    # names.append(fullsecname)
    # #import pdb; pdb.set_trace()
    # dists.append(float(h.distance(seg)))
    # sec_types.append(fullsecname.split(".")[1][:4])
#import pdb; pdb.set_trace()
df = pd.DataFrame()
df["Node ID"] = node_ids
df["Distance"] = dists
df["Conductance"] = weights
df["Type"] = sec_types
df["Name"] = names
#df["Node ID"] = node_ids
#import pdb; pdb.set_trace()
# df.to_csv("EPSPs.csv", index=False)
# import pdb; pdb.set_trace()

sim.run()

from analyze_IPSCs import calc_IPSCs
ipscs = calc_IPSCs("output/se_clamp_report.h5")
df["IPSC"] = ipscs

df.to_csv("IPSCs.csv", index=False)
# mem_pot_file = './output/v_report.h5'

# # load 
# f = h5py.File(mem_pot_file,'r')
# mem_potential = f['report']['exc_stim']['data']
# import matplotlib.pyplot as plt
# plt.plot(mem_potential[:, 0])

# plt.show()
# #import pdb; pdb.set_trace()

# low = mem_potential[3800, 0]
# high = max(mem_potential[3800:, 0])
# mag = high - low
# print(mag)

# #f = open('syn_epsps.pkl', 'rb')
# f = open('dist_test.pkl', 'rb')
# res = pickle.load(f)
# res[weight] = mag
# f.close()

# #f = open('syn_epsps.pkl', 'wb')
# f = open('dist_test.pkl', 'wb')
# pickle.dump(res, f)
# f.close()

