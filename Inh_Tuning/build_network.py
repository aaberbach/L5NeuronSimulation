from bmtk.builder import NetworkBuilder
import numpy as np
import sys
import synapses
import pandas as pd

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

synapses.load()
syn = synapses.syn_params_dicts()

net = NetworkBuilder("biophysical")

net.add_nodes(N=N, pop_name='Pyrc',
    potental='exc',
    model_type='biophysical',
    model_template='hoc:L5PCtemplate',
    morphology = None)

inh_stim = NetworkBuilder('inh_stim')
inh_stim.add_nodes(N=1,
                pop_name='inh_stim',
                #potential='exc',
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

# Create connections between Exc --> Pyr cells
net.add_edges(source=inh_stim.nodes(), target=net.nodes(),
                connection_rule=1,
                syn_weight=1,
                target_sections=['soma'],
                delay=0.1,
                distance_range=[0, 50],#(2013, Pouille et al.)
                dynamics_params='INT2PN.json',
                model_template=syn['INT2PN.json']['level_of_detail'])

# Build and save our networks
net.build()
net.save_nodes(output_dir='network')
net.save_edges(output_dir='network')

inh_stim.build()
inh_stim.save_nodes(output_dir='network')

import h5py
f = h5py.File('inh_stim_spikes.h5', 'w')
f.create_group('spikes')
f['spikes'].create_group('inh_stim')
f['spikes']['inh_stim'].create_dataset("node_ids", data=[0])
f['spikes']['inh_stim'].create_dataset("timestamps", data=[400])
f.close()

from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',
                network_dir='./network',
                tstop=500.0, dt = 0.1,
                report_vars=['v'],
                spikes_threshold=-10,
                current_clamp={           # Creates a step current from 500.ms to 1500.0 ms  
                     'amp': 0.2,
                     'delay': 300,
                     'duration': 2000,
                     'gids':"all"
                 },
                spikes_inputs=[('inh_stim', 'inh_stim_spikes.h5')],
                components_dir='../biophys_components',
                compile_mechanisms=True)