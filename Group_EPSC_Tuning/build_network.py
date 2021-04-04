from bmtk.builder import NetworkBuilder
import numpy as np
import sys
import synapses
import pandas as pd

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

#Add divergent inhibition
#Change synapse counts for inhibition and excitation based on 50 um
from clustering import *

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        inp = sys.argv[-1]
    else:
        raise Exception("no work" + str(sys.argv[-1]))

# df = pd.read_csv("Segments.csv")
# types = np.array(df["Type"])
# xs = np.array(df["X"])
# ids = np.array(df["Sec ID"])
# distances = np.array(df["Distance"])

# dendrites = np.where(((types == "dend") | (types == "apic")) & (distances >= 50))[0]
# N = len(dendrites)
N = int(inp)

np.random.seed(2129)
#np.random.seed(42)

#38% basal dend
#62% apical dend
all_dend = False

n_dend = 0#int(0.38 * N)
n_apic = 0#N - n_dend
#n_apic = int(0.62 * N)

if all_dend:
    n_dend = N
else:
    n_apic = N

synapses.load()
syn = synapses.syn_params_dicts()

net = NetworkBuilder("biophysical")

net.add_nodes(N=N, pop_name='Pyrc',
    potental='exc',
    model_type='biophysical',
    model_template='hoc:L5PCtemplate',
    morphology = None)

exc_stim = NetworkBuilder('exc_stim')
exc_stim.add_nodes(N=1,
                pop_name='exc_stim',
                potential='exc',
                model_type='virtual')
                
# def connection(source, target, id):
#         if target.node_id == id:
#             return 1
#         else:
#             return 0

# for i in range(N):
#     #import pdb; pdb.set_trace()
#     net.add_edges(source=exc_stim.nodes(), target=net.nodes(),
#                     connection_rule=connection,
#                     connection_params={"id":i},
#                     syn_weight = 1,
#                     sec_id = ids[dendrites][i],
#                     delay=0.1,
#                     sec_x = xs[dendrites][i],
#                     dynamics_params='PN2PN.json',
#                     model_template=syn['PN2PN.json']['level_of_detail'])

df = pd.read_csv("../Segments_1um.csv")
#dends = df[df["Type"] == "dend"]
dends = df[(df["Type"] == "dend") & (df["Distance"] >= 50)]
apics = df[(df["Type"] == "apic")]

#Creates n_groups functional groups from the given segmenst.
#Each functional group will have only 1 cell and 8 clusters.
def make_groups(segs, n_groups):
    groups = []

    for i in range(n_groups):
        new_group = FunctionalGroup(segs, segs.sample().iloc[0], 1, 8, "", 0, partial(make_seg_sphere, radius = 100), partial(make_seg_sphere, radius = 10), syn_per_cell=[2,8])
        groups.append(new_group)

    return groups

#Sets the location of synapses based on the given cell list.
def set_location(source, target, groups):
        index = target.node_id
        seg = groups[index].cells[0].get_seg()
        return seg.bmtk_id, seg.x

#Each of these "groups" really just represent one input cell.
all_groups = make_groups(dends, n_dend) + make_groups(apics, n_apic)

#Sets the number of synapses for each input cell.
def connector_func(source, target, groups):
        index = target.node_id
        cons = groups[index].cells[0].n_syns
        return cons

# Create connections between Exc --> Pyr cells
conn = net.add_edges(source=exc_stim.nodes(), target=net.nodes(),
                connection_rule=connector_func,
                syn_weight=1,
                connection_params={'groups': all_groups},
                #target_sections=['apic', 'dend'],
                delay=0.1,
                #distance_range=[149.0, 151.0], #0.348->0.31, 0.459->0.401
                #distance_range=[50, 2000],#(2013, Pouille et al.)
                #distance_range=[1250,2000],
                #distance_range=[-500, 500],
                dynamics_params='PN2PN.json',
                model_template=syn['PN2PN.json']['level_of_detail'])

conn.add_properties(['sec_id',"sec_x"], 
                        rule=set_location,
                        rule_params={'groups': all_groups},
                        dtypes=[np.int, np.float])

# Build and save our networks
net.build()
net.save_nodes(output_dir='network')
net.save_edges(output_dir='network')

exc_stim.build()
exc_stim.save_nodes(output_dir='network')

import h5py
f = h5py.File('exc_stim_spikes.h5', 'w')
f.create_group('spikes')
f['spikes'].create_group('exc_stim')
f['spikes']['exc_stim'].create_dataset("node_ids", data=[0])
f['spikes']['exc_stim'].create_dataset("timestamps", data=[400])
f.close()

from bmtk.utils.sim_setup import build_env_bionet

holding_v = -75

build_env_bionet(base_dir='./',
                network_dir='./network',
                tstop=500.0, dt = 0.1,
                report_vars=['v'],
                spikes_threshold=-10,
                clamp_reports=['se'],#Records se clamp currents.
                se_voltage_clamp={
                     "amps":[[holding_v, holding_v, holding_v]],
                     "durations": [[500, 0, 0]],
                     'gids': "all",
                     'rs': [0.01 for i in range(N)],
                },
                spikes_inputs=[('exc_stim', 'exc_stim_spikes.h5')],
                components_dir='../biophys_components',
                compile_mechanisms=True)