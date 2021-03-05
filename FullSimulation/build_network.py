from bmtk.builder import NetworkBuilder
import numpy as np
import sys
import synapses   
import pandas as pd       

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

from clustering import *

synapses.load()
syn = synapses.syn_params_dicts()

np.random.seed(2129)

net = NetworkBuilder("biophysical")

def lognormal(m, s):
        mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
        std = np.sqrt(np.log((s/m)**2 + 1))
        return max(np.random.lognormal(mean, std, 1), 0)


scale_div = 1#10

# Dend Excitatory: 7186.0
# Dend Inhibitory: 718.0
# Apic Excitatory: 10417.0
# Apic Inhibitory: 1041.0
# Soma Inhibitory: 148

avg_syn_per_cell = 5 #Average number of synapses from each input cell.

num_dend_exc = (7186 // avg_syn_per_cell) // scale_div
num_apic_exc = (10417 // avg_syn_per_cell) // scale_div

num_dend_inh = (718 // avg_syn_per_cell) // scale_div
num_apic_inh = (1041 // avg_syn_per_cell) // scale_div
num_soma_inh = (148 // avg_syn_per_cell) // scale_div

exc_fr_mean = 0.1
exc_fr_std = 0.5
inh_fr = 7 #* scale_div

clust_per_group = 8
num_dend_groups = 70 // scale_div
num_apic_groups = 100 // scale_div

dend_groups = []
apic_groups = []

##################################################################################
###################################Pyr Type C#####################################

net.add_nodes(N=1, pop_name='Pyrc',
    potental='exc',
    model_type='biophysical',
    model_template='hoc:L5PCtemplate',
    morphology = None)

##################################################################################
###################################External Networks##############################

#########################Excitatory Inputs##########################33####33
from clustering import *
#num_exc = (num_apic_exc + num_dend_exc) #* N

# External excitatory inputs
exc_stim = NetworkBuilder('exc_stim')

df = pd.read_csv("Segments.csv")
dends = df[df["Type"] == "dend"]
apics = df[(df["Type"] == "apic")]

#Sets the number of synapses for each input cell.
def connector_func(source, target, cells):
        return [cell.n_syns for cell in cells]

#Sets the location of synapses based on the given cell list.
def set_location(source, target, cells, start_id):
        index = source.node_id - start_id
        seg = cells[index].get_seg()
        return seg.bmtk_id, seg.x

def build_nodes(stim, cells_per, clust_per, group_list, segs, base_name, n_groups, start=0):
        start_id = start

        for i in range(n_groups):
                name = base_name + str(i)

                stim.add_nodes(N=cells_per,
                        pop_name=name,
                        potential="exc",
                        model_type='virtual')

                new_group = FunctionalGroup(segs, segs.sample().iloc[0], cells_per, clust_per, name, start_id, partial(make_seg_sphere, radius = 100), partial(make_seg_sphere, radius = 10))
                group_list.append(new_group)
                start_id += cells_per

        return start_id

def build_edges(stim, group_list, base_name, target_net):
        for i in range(len(group_list)):
                name = base_name + str(i)
                group = group_list[i]

                conn = net.add_edges(source=stim.nodes(pop_name=name), target=target_net.nodes(),
                        iterator="all_to_one",
                        connection_rule=connector_func,
                        connection_params={'cells': group.cells},#NEED BE SURE IN SAME ORDER ALWYAS
                        syn_weight=1,
                        delay=0.1,
                        dynamics_params='PN2PN.json',
                        model_template=syn['PN2PN.json']['level_of_detail'],)

                conn.add_properties(['sec_id',"sec_x"], 
                        rule=set_location,
                        rule_params={'cells': group.cells, 'start_id': group.start_id},#NEED BE SURE IN SAME ORDER ALWYAS
                        dtypes=[np.int, np.float])

end = build_nodes(exc_stim, num_dend_exc//num_dend_groups, clust_per_group, dend_groups, dends, "dend", num_dend_groups)
build_nodes(exc_stim, num_apic_exc//num_apic_groups, clust_per_group, apic_groups, apics, "apic", num_apic_groups, start=end)

build_edges(exc_stim, dend_groups, "dend", net)
build_edges(exc_stim, apic_groups, "apic", net)

#########################Inhibitory Inputs#####################################

# External proximal inhibitory inputs
prox_inh_stim = NetworkBuilder('prox_inh_stim')
prox_inh_stim.add_nodes(N=num_soma_inh,
                pop_name='prox_inh_stim',
                potential='exc',
                model_type='virtual')

#Edges
net.add_edges(source=prox_inh_stim.nodes(), target=net.nodes(),
                connection_rule=1,
                #connection_params={'num_per': num_soma_inh , 'start':0},
                syn_weight=1,
                delay=0.1,
                dynamics_params='INT2PN.json',
                model_template=syn['INT2PN.json']['level_of_detail'],
                distance_range=[-2000, 2000.0],
                target_sections=['somatic'])

# External distal inhibitory inputs
dist_inh_stim = NetworkBuilder('dist_inh_stim')

dist_inh_stim.add_nodes(N=num_dend_inh,
                pop_name='dend',
                potential='exc',
                model_type='virtual')

dist_inh_stim.add_nodes(N=num_apic_inh,
                pop_name='apic',
                potential='exc',
                model_type='virtual')

#Dend edges.
net.add_edges(source=dist_inh_stim.nodes(pop_name="dend"), target=net.nodes(),
                connection_rule=1,
                syn_weight=1,
                delay=0.1,
                dynamics_params='INT2PN.json',
                model_template=syn['INT2PN.json']['level_of_detail'],
                distance_range=[50, 2000.0],
                target_sections=['dend'])

#Apic edges.
net.add_edges(source=dist_inh_stim.nodes(pop_name="apic"), target=net.nodes(),
                connection_rule=1,
                syn_weight=1,
                delay=0.1,
                dynamics_params='INT2PN.json',
                model_template=syn['INT2PN.json']['level_of_detail'],
                distance_range=[50, 2000.0],
                target_sections=['apic'])


# Build and save our networks
net.build()
net.save_nodes(output_dir='network')
net.save_edges(output_dir='network')

exc_stim.build()
exc_stim.save_nodes(output_dir='network')

prox_inh_stim.build()
prox_inh_stim.save_nodes(output_dir='network')

dist_inh_stim.build()
dist_inh_stim.save_nodes(output_dir='network')

##############################External Spike Rasters#######################

from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator
from bmtk.utils.reports.spike_trains.spikes_file_writers import write_csv

# exc_dend_frs = []
# exc_apic_frs = []

# exc_means = []
# exc_stds = []
# exc_maxs = []

# #for i in range(N):
# #fr_mean = np.random.uniform(exc_fr_mean - 0.00, exc_fr_mean + 0.00, 1)
# #fr_mean = np.random.uniform(exc_fr_mean + 0.5, exc_fr_mean + 0.55, 1)
# dend_frs = [min(float(lognormal(exc_fr_mean, exc_fr_std)), exc_fr_mean + 8*exc_fr_std) for _ in range(num_dend_exc)]
# apic_frs = [min(float(lognormal(exc_fr_mean, exc_fr_std)), exc_fr_mean + 8*exc_fr_std) for _ in range(num_apic_exc)]

# frs = dend_frs + apic_frs

# exc_dend_frs += dend_frs
# exc_apic_frs += apic_frs

# exc_means.append(np.mean(frs))
# exc_stds.append(np.std(frs))
# exc_maxs.append(np.max(frs))

# exc_frs = exc_dend_frs + exc_apic_frs

# fr_df = pd.DataFrame()
# #fr_df['gid'] = [i + num_exc+num_inh for i in range(N)]
# fr_df['gid'] = [i for i in range(N)]
# fr_df['fr_mean'] = exc_means
# fr_df['fr_std'] = exc_stds
# fr_df['fr_max'] = exc_maxs

#fr_df.to_csv('frs_temp.csv', index=False)
seconds = 1

def generate_spikes(psg, group, fr, times):
        psg.add(node_ids=range(group.start_id, group.start_id + len(group.cells)),
                        firing_rate=fr,
                        times=times)

exc_psg = PoissonSpikeGenerator(population='exc_stim')
exc_psg.add(node_ids = range(num_apic_exc + num_dend_exc),
                firing_rate=0.01,
                times=(0.0*1000, seconds*1000))
exc_psg.to_sonata('exc_stim_spikes.h5')
# exc_psg = PoissonSpikeGenerator(population='exc_stim')
# for i in range(num_exc):
#         exc_psg.add(node_ids=[i],  
#                 firing_rate=exc_frs[i]/1000,    
#                 times=(0.0*1000, seconds*1000))     
# exc_psg.to_sonata('exc_stim_spikes.h5')

# inh_psg = PoissonSpikeGenerator(population='inh_stim')
# inh_psg.add(node_ids=range(num_inh), 
#         firing_rate=inh_fr/1000,  
#         times=(0.0*1000, seconds*1000))   
# inh_psg.to_sonata('inh_stim_spikes.h5')

# from crop_raster import crop_raster
# crop_raster("rhythmic_inh_spikes.h5", 'inh_stim_spikes.h5', 120000, num_inh)


from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',
                network_dir='./network',
                dt = 0.1, tstop=seconds * 1000.0,
                report_vars=['v', 'cai'],
                dL = 5,
                spikes_threshold=-10,
                spikes_inputs=[('exc_stim', 'exc_stim_spikes.h5')],#, ('inh_stim', 'inh_stim_spikes.h5')],
                components_dir='../biophys_components',
                compile_mechanisms=True)