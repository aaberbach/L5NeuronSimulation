from bmtk.builder import NetworkBuilder
import numpy as np
import sys
import synapses   
import pandas as pd                                                                                                             
synapses.load()
syn = synapses.syn_params_dicts()

np.random.seed(2129)

net = NetworkBuilder("biophysical")

##################################################################################
###################################Pyr Type C#####################################

net.add_nodes(N=1, pop_name='Pyrc',
    potental='exc',
    model_type='biophysical',
    model_template='hoc:L5PCtemplate',
    morphology = None)

##################################################################################
###################################External Networks##############################

# External excitatory inputs
exc_stim = NetworkBuilder('exc_stim')
exc_stim.add_nodes(N=1,
                pop_name='exc_stim',
                potential='exc',
                model_type='virtual')

# External inhibitory inputs
inh_stim = NetworkBuilder('inh_stim')
inh_stim.add_nodes(N=1,
                pop_name='inh_stim',
                potential='exc',
                model_type='virtual')

##################################################################################
###################################Edges##########################################

#Inhibitory on soma.
net.add_edges(source=inh_stim.nodes(), target=net.nodes(),
                connection_rule=1,
                syn_weight=1,
                delay=0.1,
                dynamics_params='INT2PN.json',
                model_template=syn['INT2PN.json']['level_of_detail'],
                distance_range=[-2000, 2000.0],
                target_sections=['somatic'])

# Create connections between Exc --> Pyr cells

#Excitatory on basal dendrites.
net.add_edges(source=exc_stim.nodes(), target=net.nodes(),
                connection_rule=1,
                syn_weight=1,
                target_sections=['dend'],
                delay=0.1,
                distance_range=[50.0, 2000.0],
                dynamics_params='PN2PN.json',
                model_template=syn['PN2PN.json']['level_of_detail'])

# #Excitatory on apical dendrites.
# net.add_edges(source=exc_stim.nodes(), target=net.nodes(),
#                 connection_rule=correct_cell,
#                 connection_params={'num_per': num_apic_exc, 'start':start},
#                 syn_weight=1,
#                 target_sections=['apic'],
#                 delay=0.1,
#                 distance_range=[50.0, 2000.0],
#                 dynamics_params='PN2PN.json',
#                 model_template=syn['PN2PN.json']['level_of_detail'])

# Build and save our networks
net.build()
net.save_nodes(output_dir='network')
net.save_edges(output_dir='network')

exc_stim.build()
exc_stim.save_nodes(output_dir='network')

inh_stim.build()
inh_stim.save_nodes(output_dir='network')

node_ids = np.full(9, 0)
spacing = 10
timestamps = np.linspace(0 + 500, 8*spacing + 500, 9)
#timestamps *= 10
#import pdb; pdb.set_trace()
import h5py
key = 'inh_stim'
f = h5py.File('stim_spikes.h5', 'w')
group = f.create_group('spikes')
group = group.create_group(key)
group.create_dataset("node_ids", data = node_ids)
group.create_dataset("timestamps", data = timestamps)
f.close()

from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',
                network_dir='./network',
                dt = 0.1, tstop=1000.0,
                report_vars=['v', 'cai'],
                dL=20,
                # current_clamp={           # Creates a step current from 500.ms to 1500.0 ms  
                #      'amp': 0.793,
                #      #'std': [0.0, 0.0],
                #      'delay': 700,
                #      'duration': 2000,
                #      'gids':"0"
                #  },
                # clamp_reports=['se'],#Records se clamp currents.
                # se_voltage_clamp={
                #         "amps":[[0, 0, 0]],
                #         "durations": [[120000, 0, 0]],
                #         'gids': [0],
                #         'rs': [0.01],
                # },
                spikes_threshold=-10,
                spikes_inputs=[(key, 'stim_spikes.h5')],
                components_dir='../biophys_components',
                compile_mechanisms=True)