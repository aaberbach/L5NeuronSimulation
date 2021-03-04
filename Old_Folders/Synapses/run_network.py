import sys, os
from bmtk.simulator import bionet
import numpy as np
import pandas as pd
import h5py
from neuron import h
from scipy.stats import skew
import synapses
from bmtk.simulator.bionet.pyfunction_cache import add_weight_function
import pickle
import matplotlib.pyplot as plt
import load_processor

#load_processor.load()

synapses.load()
#import pdb; pdb.set_trace()

pc = h.ParallelContext()  # object to access MPI methods
MPI_size = int(pc.nhost())
MPI_rank = int(pc.id())

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        inp = sys.argv[-1]
    else:
        raise Exception("no work" + str(sys.argv[-1]))

fname = str(inp)

config_file = 'simulation_config.json'

def gaussianBL(edge_props, source, target):
    w0 = edge_props["syn_weight"]
    sigma = edge_props["weight_sigma"]
    try:
        maximum = edge_props["weight_max"]
        return min(maximum, np.random.normal(w0, sigma, 1))
    except:
        return np.random.normal(w0, sigma, 1)

def lognormal(edge_props, source, target):
    m = edge_props["syn_weight"]
    s = edge_props["weight_sigma"]
    mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
    std = np.sqrt(np.log((s/m)**2 + 1))

    try:
        maximum = edge_props["weight_max"]
        return max(min(maximum, np.random.lognormal(mean, std, 1)), 0)
    except:
        return max(0, np.random.lognormal(mean, std, 1))

add_weight_function(lognormal)
add_weight_function(gaussianBL)

conf = bionet.Config.from_json(config_file, validate=True)
conf.build_env()

graph = bionet.BioNetwork.from_config(conf)
#import pdb; pdb.set_trace()
sim = bionet.BioSimulator.from_config(conf, network=graph)
#import pdb; pdb.set_trace()

# cell = graph.get_local_cells()[1]
# memb = h.Vector()
# memb.record(cell.hobj.soma[0](0.5)._ref_v)

#import pdb; pdb.set_trace()

cells = graph.get_local_cells()
#import pdb; pdb.set_trace()
#import pdb; pdb.set_trace()
#gid_min = min(cells.keys())
syn = cells[0].connections()[0]._syn
nc = cells[0].connections()[0]._connector
#import pdb; pdb.set_trace()
capool = h.Vector()
capool.record(syn._ref_Capoolcon)

thresh1 = h.Vector()
thresh1.record(syn._ref_threshold1)

i_ampa = h.Vector()
i_ampa.record(syn._ref_i_ampa)

i_nmda = h.Vector()
i_nmda.record(syn._ref_i_nmda)

g_ampa = h.Vector()
g_ampa.record(syn._ref_g_ampa)

g_nmda = h.Vector()
g_nmda.record(syn._ref_g_nmda)

gbar_ampa = h.Vector()#CONSTANT
gbar_ampa.record(syn._ref_gbar_ampa)

gbar_nmda = h.Vector()#CONSTANT
gbar_nmda.record(syn._ref_gbar_nmda)

#import pdb; pdb.set_trace()
v = h.Vector()
v.record(nc.postseg()._ref_v)

mult = h.Vector()
mult.record(syn._ref_comb)

W_ampa = h.Vector()
W_ampa.record(syn._ref_W_ampa)


# exc_strengths = {}
# inh_strengths = {}
# comb_strens = {}

# for gid, cell in cells.items():
#     exc_strens = []
#     inh_strens = []
#     fr_comb_total = 0
#     for con in cell._connections:
#         #import pdb; pdb.set_trace()
#         #if len(np.array(con.source_node._train_vec)) > 0:
#             #import pdb; pdb.set_trace()
#         if con._edge_prop.source_population == 'exc_stim':
#             fr_comb_total += len(np.array(con.source_node._train_vec)) * con._syn.initW
#             #if len(np.array(con.source_node._train_vec)) > 0:
                
#             #exc_strens.append(con.syn_weight)
#             exc_strens.append(con._syn.initW)
#         elif con._edge_prop.source_population == 'inh_stim':
#             #inh_strens.append(con.syn_weight)
#             inh_strens.append(con._syn.initW)
#         else:
#             raise Exception("Source pop is: " + str(con._edge_prop.source_population))

#     #print(gid, ":", fr_comb_total / len(cell._connections))
#     comb_strens[gid] = fr_comb_total / len(exc_strens)
    
#     #exc_strengths[gid - gid_min] = exc_strens
#     #inh_strengths[gid - gid_min] = inh_strens
#     exc_strengths[gid] = exc_strens
#     inh_strengths[gid] = inh_strens
# #import pdb; pdb.set_trace()
print(synapses.max_exc)
sim.run()
# pc.barrier()
# #import pdb; pdb.set_trace()

# raster_file = './output/spikes.h5'

# frs = {}
# local_gids = list(exc_strengths.keys())
# #local_gids = local_gids - np.min(local_gids)
# #import pdb; pdb.set_trace()

# for key in local_gids:
#     frs[key] = 0

# try:
#     f = h5py.File(raster_file,'r')
    
#     spike_keys = list(f['spikes'].keys())
#     if len(spike_keys) > 1:
#         raise Exception("Spike keys: " + str(spike_keys))

#     spike_key = list(f['spikes'].keys())[0]
#     timestamps = f['spikes'][spike_key]['timestamps'].value
#     gids = f['spikes'][spike_key]['node_ids'].value

#     for i in range(len(gids)):
#         gid = gids[i] + min(local_gids)
#         if gid in local_gids and timestamps[i] >= 200:
#             frs[gid] += 1
# except:
#     print("No spikes.")

# df = pd.DataFrame()
# dicts = [{"gid": gid, "FR": frs[gid] / 5, "num_exc": len(exc_strengths[gid]), "num_inh": len(inh_strengths[gid]),
#             "avg_exc": np.mean(exc_strengths[gid]), "avg_inh": np.mean(inh_strengths[gid]), 
#             "max_exc": np.max(exc_strengths[gid]), "max_inh": np.max(inh_strengths[gid]),
#             "std_exc": np.std(exc_strengths[gid]), "std_inh": np.std(inh_strengths[gid]),
#             "skew_exc": skew(exc_strengths[gid]), "skew_inh": skew(inh_strengths[gid]), "comb_stren": comb_strens[gid]} for gid in local_gids]
# df = pd.DataFrame(dicts)
# #df.set_index("gid")
# df.to_csv(fname+str(MPI_rank)+'.csv', index=False)

# #import pdb; pdb.set_trace()

# pc.barrier()

# if MPI_rank == 0:
#     base_df = pd.read_csv(fname+"0.csv", index_col="gid")
#     res_df = pd.concat([base_df] + [pd.read_csv(fname+str(rank)+".csv", index_col="gid") for rank in range(1, MPI_size)])
#     frs_df = pd.read_csv('frs_temp.csv', index_col="gid")
#     res_df = res_df.join(frs_df)
#     os.remove('frs_temp.csv')
#     [os.remove(fname+str(rank)+".csv") for rank in range(MPI_size)]
#     res_df.to_csv(fname+".csv")
from bmtk.analyzer.cell_vars import plot_report

_ = plot_report(config_file='simulation_config.json', node_ids=[0], report_name='v_report')

# plt.figure()
# plt.plot(capool, label="Capoolcon")
# plt.plot(thresh1)
# plt.legend()

plt.figure()
plt.plot(i_ampa,  label="i_ampa")
plt.legend()

plt.figure()
plt.plot(i_nmda,  label="i_nmda")
plt.legend()

# plt.figure()
# plt.plot(g_ampa,  label="g_ampa")
# plt.legend()

# plt.figure()
# plt.plot(g_nmda,  label="g_nmda")
# plt.legend()

# plt.figure()
# plt.plot(W_ampa,  label="W_ampa")
# plt.legend()
# plt.show()


# plt.figure()
# plt.plot(gbar_ampa,  label="gbar_ampa")
# plt.legend()

# plt.figure()
# plt.plot(gbar_nmda,  label="gbar_nmda")
# plt.legend()

plt.figure()
plt.plot(v, label="v")
plt.legend()

plt.figure()
plt.plot(mult, label="mult")
plt.legend()

sv = [syn.sfunc(i) for i in v]
plt.figure()
plt.plot(sv, label="sv")
plt.legend()

plt.show()