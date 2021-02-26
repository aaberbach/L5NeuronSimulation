import sys, os
from bmtk.simulator import bionet
import numpy as np
import h5py
#import synapses
import pandas as pd
from neuron import h
from matplotlib import cm

np.random.seed(2129)
#np.random.seed(42)

# if __name__ == '__main__':
#     if __file__ != sys.argv[-1] and __file__ != sys.argv[-2]:
#         m = float(sys.argv[-2])
#         s = float(sys.argv[-1])
#     else:
#         raise Exception("The inputted command argument not functional: " + str(sys.argv[-1]) + " " + str(sys.argv[-2]))

# synapses.set_pyr_w(m, s)
# synapses.load()

config_file = 'simulation_config.json'

conf = bionet.Config.from_json(config_file, validate=True)
conf.build_env()

graph = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=graph)


dic = {}
dic["high"] = {"dist": [], "Ih": [], "LVAst":[], "HVA": []}
dic["low"] = {"dist": [], "Ih": [], "LVAst":[], "HVA": []}

#segs = list(graph.get_local_cells().values())[0]._morph.seg_prop
#import pdb; pdb.set_trace()
cell = list(graph.get_local_cells().values())[0]
hobj = cell.hobj
import matplotlib.pyplot as plt
ps = h.PlotShape(False)  # False tells h.PlotShape not to use NEURON's gui
h.distance(sec=cell.hobj.soma[0])
for sec in hobj.all:
    fullsecname = sec.name()
    sec_type = fullsecname.split(".")[1][:4]
    sec_id = int(fullsecname.split("[")[-1].split("]")[0])
    #sec.v = np.random.random(1)[0]*1
    #print(sec.v)
    for seg in sec:
        #import pdb; pdb.set_trace()
        dist = h.distance(seg)
        if sec_type == "soma": #light blue
            seg.v = -20
        if sec_type == "axon": #neon green
            seg.v = -30
        if sec_type == "dend": #dark blue
            seg.v = 0
        if sec_type == "apic":
            if dist < 500: #magenta
                seg.v = 20
            if dist >= 500:
                if sec_id >= 60: #Red
                    # dic["high"]["dist"].append(dist)
                    # dic["high"]["Ih"].append(seg.gIhbar_Ih)
                    # dic["high"]["LVAst"].append(seg.gCa_LVAstbar_Ca_LVAst)
                    # dic["high"]["HVA"].append(seg.gCa_HVAbar_Ca_HVA)
                    seg.v = 40
                else: #teel
                    # dic["low"]["dist"].append(dist)
                    # dic["low"]["Ih"].append(seg.gIhbar_Ih)
                    # dic["low"]["LVAst"].append(seg.gCa_LVAstbar_Ca_LVAst)
                    # dic["low"]["HVA"].append(seg.gCa_HVAbar_Ca_HVA)
                    seg.v = -10

# for key in dic["low"].keys():
#     plt.figure()
#     plt.scatter(dic["low"]["dist"], dic["low"][key])
#     plt.title("low" + str(key))

#     plt.figure()
#     plt.scatter(dic["high"]["dist"], dic["high"][key])
#     plt.title("high" + str(key))

# plt.show()

import pdb; pdb.set_trace()
cmap = cm.hsv
#import pdb; pdb.set_trace()

#ps.color_all(1)
p = ps.plot(plt, cmap=cmap)
#import pdb; pdb.set_trace()

# segs = cell._morph.seg_prop
# h.distance(sec=cell.hobj.soma[0])
# for sec in hobj.all:
#     fullsecname = sec.name()
#     sec_type = fullsecname.split(".")[1][:4]
#     sec_id = int(fullsecname.split("[")[-1].split("]")[0])
#     ps.color(0, sec=sec)
#     if sec_id >= 60 and sec_type == "apic":
#         for seg in sec:
#             dist = h.distance(seg)
#             #help(p.mark())
#             p.mark(seg).v = 1000
#             #ps.mark(seg)

#for con in cell._connections:
    #print(con._connector.postseg())
    #print(p.mark)
    #import pdb; pdb.set_trace()
    #p.mark(con._connector.postseg())
#ps.plot(plt).mark(cell.soma[0](0.5))
plt.show()

#sim.run()

