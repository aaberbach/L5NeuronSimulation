from bmtk.builder import NetworkBuilder
import numpy as np
import sys
import synapses   
import pandas as pd                                                                                                             
synapses.load()
syn = synapses.syn_params_dicts()

np.random.seed(2129)

#############################PLANS##############################################

# There will be some number of correlated excitatory input groups
#       -Each with own pop_name (numbered indices)
# Each will consist of some number of neurons which each have 2-8 synapses
#       -Most likely uniformly distributed
# Each input group will be given a random center of concentration, and individual neurons
# will then be distiributed randomly in some way around that center.
#       -Manner of distribution should be interchangeable
# Then, each individual neuron will distribute its synapses in some random
# way around the center it is assigned.
#       -Distribution should again be interchangeable
# There should be some way to split the input groups into smaller groups for partial clustering.

#############################PLANS##############################################

#############################PLANS##############################################

# There will be some number of correlated excitatory input groups
#       -Each with own pop_name (numbered indices)
# Each will consist of some number of neurons which each have 2-8 synapses
#       -Most likely uniformly distributed
# Each functional group will be given a random center,
# 8 clusters will be placed in some manner around that center.
# Then, each individual neuron will distribute its synapses in some random
# way to the 8 clusters of the functional group.
#       -Distribution should again be interchangeable
# There should be some way to split the input groups into smaller groups for partial clustering.
#       -Could just be a variable saying how many cluster groups each
#       functional unit is broken into.

#############################PLANS##############################################

# Bins are proportions of width.
# cdvdt=i_ampa+i_nmda+active/passive currents in the cell


# if __name__ == '__main__':
#     if __file__ != sys.argv[-1]:
#         inp = sys.argv[-1]
#     else:
#         raise Exception("no work" + str(sys.argv[-1]))

#N = 1#int(inp)

net = NetworkBuilder("biophysical")

def lognormal(m, s):
        mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
        std = np.sqrt(np.log((s/m)**2 + 1))
        return max(np.random.lognormal(mean, std, 1), 0)


scale_div = 10#10

# Dend Excitatory: 7186.0
# Dend Inhibitory: 718.0
# Apic Excitatory: 10417.0
# Apic Inhibitory: 1041.0
# Soma Inhibitory: 148

avg_syn_per_cell = 5 #Average numver of synapses from each input cell.

num_dend_exc = (7186 // avg_syn_per_cell) // scale_div
num_apic_exc = (10417 // avg_syn_per_cell) // scale_div

num_dend_inh = (718 // avg_syn_per_cell) // scale_div
num_apic_inh = (1041 // avg_syn_per_cell) // scale_div
num_soma_inh = (148 // avg_syn_per_cell) // scale_div

exc_fr_mean = 0.1
exc_fr_std = 0.5
inh_fr = 7 #* scale_div

clust_per_group = 8
num_dend_groups = 5
num_apic_groups = 2

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

df = pd.read_csv("Segments_small.csv")
dends = df[df["Type"] == "dend"]
apics = df[(df["Type"] == "apic")]

#Sets the number of synapses for each input cell.
def connector_func(source, target, cells):
        #import pdb; pdb.set_trace()
        return [cell.n_syns for cell in cells]

#Sets the location of synapses based on the given cell list.
def set_location(source, target, cells, start_id):
        #import pdb; pdb.set_trace()
        index = source.node_id - start_id
        seg = cells[index].get_seg()
        return seg.bmtk_id, seg.x

start_id = 0
cells_per_group = num_dend_exc // num_dend_groups
for i in range(num_dend_groups):
        name = "dend"+str(i)
        exc_stim.add_nodes(N=num_dend_exc//num_dend_groups,
                        pop_name=name,
                        potential="exc",
                        model_type='virtual')

        new_group = FunctionalGroup(dends, dends.sample().iloc[0], cells_per_group, clust_per_group, name, start_id, partial(make_seg_sphere, radius = 100), partial(make_seg_sphere, radius = 10))
        dend_groups.append(new_group)
        start_id += cells_per_group

for i in range(num_dend_groups):
        name = "dend"+str(i)
        group = dend_groups[i]
        conn = net.add_edges(source=exc_stim.nodes(pop_name=name), target=net.nodes(),
                iterator="all_to_one",
                connection_rule=connector_func,
                connection_params={'cells': group.cells},#NEED BE SURE IN SAME ORDER ALWYAS
                syn_weight=1,
                #target_sections=['dend'],
                delay=0.1,
                #distance_range=[50.0, 2000.0],
                dynamics_params='PN2PN.json',
                model_template=syn['PN2PN.json']['level_of_detail'],)
                #sec_id=3,
                #sec_x=0.5)

        conn.add_properties(['sec_id',"sec_x"], 
                    rule=set_location,
                    rule_params={'cells': group.cells, 'start_id': group.start_id},#NEED BE SURE IN SAME ORDER ALWYAS
                    dtypes=[np.int, np.float])

        

#import pdb; pdb.set_trace()


# exc_stim.add_nodes(N=num_exc,
#                 pop_name='exc_stim',
#                 potential='exc',
#                 model_type='virtual')

#########################Inhibitory Inputs#####################################

num_inh = num_soma_inh #* N

# External inhibitory inputs
inh_stim = NetworkBuilder('inh_stim')
inh_stim.add_nodes(N=num_inh,
                pop_name='inh_stim',
                potential='exc',
                model_type='virtual')

##################################################################################
###################################Edges##########################################

# Here we specify which set of nodes to use as sources and targets. Our source/pre-synaptic cells are all thamalus cells with the property "pop_name=tON", which in this case is every thalmus cell (We could also use source=thalamus.nodes(), or source={'level_of_detail': 'filter'}). The target/post-synaptic is all cell(s) of the "cortex" network.

# def correct_cell(source, target, num_per, start):
#         #import pdb; pdb.set_trace()
#         sid = source.node_id - start
#         tid = target.node_id
#         #import pdb; pdb.set_trace()

#         lower_bound = num_per * tid

#         upper_bound = lower_bound + num_per

#         if sid < upper_bound and sid >= lower_bound:
#                 #print("connecting cell {} to {}".format(sid+start,tid))
#                 return 1
#         else:
#                 return None

#https://nbviewer.jupyter.org/github/AllenInstitute/bmtk/blob/develop/docs/tutorial/NetworkBuilder_Intro.ipynb
#Goo to bottom for add_properties

# #Inhibitory on soma.
# net.add_edges(source=inh_stim.nodes(), target=net.nodes(),
#                 connection_rule=correct_cell,
#                 connection_params={'num_per': num_soma_inh , 'start':0},
#                 syn_weight=1,
#                 delay=0.1,
#                 dynamics_params='INT2PN.json',
#                 model_template=syn['INT2PN.json']['level_of_detail'],
#                 distance_range=[-2000, 2000.0],
#                 target_sections=['somatic'])

# #start += np.sum(inh_bounds_apic)
# start = 0

# # Create connections between Exc --> Pyr cells

# #Excitatory on basal dendrites.
# net.add_edges(source=exc_stim.nodes(), target=net.nodes(),
#                 connection_rule=correct_cell,
#                 connection_params={'num_per': num_dend_exc, 'start':start},
#                 syn_weight=1,
#                 target_sections=['dend'],
#                 delay=0.1,
#                 distance_range=[50.0, 2000.0],
#                 dynamics_params='PN2PN.json',
#                 model_template=syn['PN2PN.json']['level_of_detail'])

# #start += np.sum(exc_bounds_dend)
# start += np.sum(num_dend_exc )#* N)

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
                dL = 1,
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
                spikes_inputs=[('exc_stim', 'exc_stim_spikes.h5')],#, ('inh_stim', 'inh_stim_spikes.h5')],
                components_dir='../biophys_components',
                compile_mechanisms=True)