from bmtk.builder import NetworkBuilder
import numpy as np
import sys
import synapses   
import h5py
import pandas as pd       

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

#Add divergent inhibition
#Change synapse counts for inhibition and excitation based on 50 um
from raster_maker import *
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

# Dend Excitatory: 6509.0
# Beta Dend Inhibitory: 650.0
# Gamma Dend Inhibitory: 67.0
# Apic Excitatory: 10417.0
# Apic Inhibitory: 1041.0
# Soma Inhibitory: 148

avg_syn_per_cell = 5 #Average number of synapses from each input cell.

inh_syn_per_cell = 3

num_dend_exc = (6509 // avg_syn_per_cell) // scale_div
num_apic_exc = (10417 // avg_syn_per_cell) // scale_div

num_dend_inh = (650 // inh_syn_per_cell) // scale_div
num_apic_inh = (1041 // inh_syn_per_cell) // scale_div

num_prox_dend_inh = (67 // inh_syn_per_cell) // scale_div
num_soma_inh = (148 // inh_syn_per_cell) // scale_div

# exc_fr_mean = 2#0.1
# exc_fr_std = 0.5

prox_fr_mean = 16.9
prox_fr_std = 14.3

dist_fr_mean = 3.9
dist_fr_std = 4.9

clust_per_group = 8
num_dend_groups = 65 // scale_div
num_apic_groups = 104 // scale_div

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
                pop_name='on_soma',
                potential='exc',
                model_type='virtual')

prox_inh_stim.add_nodes(N=num_prox_dend_inh,
                pop_name='on_dend',
                potential='exc',
                model_type='virtual')

#Edges

def uniform_connect(source, target, low=1, high=5):
        #return np.random.uniform(low, high)
        return np.random.randint(low=low, high=high + 1)

#On soma.
net.add_edges(source=prox_inh_stim.nodes(pop_name='on_soma'), target=net.nodes(),
                connection_rule=uniform_connect,
                #connection_params={'num_per': num_soma_inh , 'start':0},
                syn_weight=1,
                delay=0.1,
                dynamics_params='INT2PN.json',
                model_template=syn['INT2PN.json']['level_of_detail'],
                distance_range=[-2000, 2000.0],
                target_sections=['somatic'])

#On dendrites within 50 um
net.add_edges(source=prox_inh_stim.nodes(pop_name='on_dend'), target=net.nodes(),
                connection_rule=uniform_connect,
                #connection_params={'num_per': num_soma_inh , 'start':0},
                syn_weight=1,
                delay=0.1,
                dynamics_params='INT2PN.json',
                model_template=syn['INT2PN.json']['level_of_detail'],
                distance_range=[0, 50.0],
                target_sections=['dend'])

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
                connection_rule=uniform_connect,
                syn_weight=1,
                delay=0.1,
                dynamics_params='INT2PN.json',
                model_template=syn['INT2PN.json']['level_of_detail'],
                distance_range=[50, 2000.0],
                target_sections=['dend'])

#Apic edges.
net.add_edges(source=dist_inh_stim.nodes(pop_name="apic"), target=net.nodes(),
                connection_rule=uniform_connect,
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
seconds = 2
times = (0, seconds)

import scipy.stats as st
from functools import partial

levy_dist = partial(st.levy_stable.rvs, alpha=1.37, beta=-1.00, loc=0.92, scale=0.44, size=1)

norm_dist = partial(st.norm.rvs, loc=5, scale=1, size=1)
#Generates the spike raster for a given group.
#The group has the same noise.
def gen_group_spikes(group, seconds):
        z = make_noise(num_samples=(int(seconds*1000))-1,num_traces=group.n_cells)
        df = make_spikes(exp=True, dist=levy_dist, numUnits=group.n_cells,rateProf=z[0,:])
        return df

def raster_to_sonata(node_ids, timestamps, key, file):
        f = h5py.File(file, 'w')
        group = f.create_group('spikes')
        group = group.create_group(key)
        group.create_dataset("node_ids", data = node_ids.astype("u4"))
        group.create_dataset("timestamps", data = timestamps)
        f.close()

def gen_spikes(dend_groups, apic_groups, times, file):
        node_ids = []
        timestamps = []

        length = times[1] - times[0]
        buffer = times[0]

        for group in (dend_groups + apic_groups):
                df = gen_group_spikes(group, length)
                node_ids = np.concatenate((node_ids, np.array(df["node_ids"]) + group.start_id))
                timestamps = np.concatenate((timestamps, np.array(df["timestamps"])))

        timestamps += buffer * 1000
        raster_to_sonata(node_ids, timestamps, "exc_stim", file)
        

gen_spikes(dend_groups, apic_groups, times, 'exc_stim_spikes.h5')

################ INH FR PROFILES FROM EXC RASTER #########

f = h5py.File('exc_stim_spikes.h5')
ts = f['spikes']['exc_stim']['timestamps']
nid = f['spikes']['exc_stim']['node_ids']




h = np.histogram(ts,bins=np.arange(0,seconds*1000,1))

#scale_factor = 3.47
#plt.plot(h[1][1:],scale_factor*h[0]/(0.001*(np.max(nid)+1)),label='scaled, not shifted')
#plt.title('FR: {} +/- {}'.format(np.mean(scale_factor*h[0]/(0.001*(np.max(nid)+1))),
#				 np.std(scale_factor*h[0]/(0.001*(np.max(nid)+1)))))


fr_prof = h[0]/(0.001*(np.max(nid)+1))

plt.plot(h[1][1:],0.5+(fr_prof-np.min(fr_prof))/(np.max(fr_prof)-np.min(fr_prof)),label='scaled, baseline')
time_shift = 4 # ms
wrap = fr_prof[-4:]
fr_prof[4:] = fr_prof[0:-4]
fr_prof[0:4] = wrap
plt.plot(h[1][1:],0.5+(fr_prof-np.min(fr_prof))/(np.max(fr_prof)-np.min(fr_prof)),label='scaled, time shift')
plt.legend()
plt.show()

# exc_psg = PoissonSpikeGenerator(population='exc_stim')
# exc_psg.add(node_ids = range(num_apic_exc + num_dend_exc),
#                 firing_rate=exc_fr_mean,
#                 #times=(0.0*1000, seconds*1000))
#                 times=times)
# exc_psg.to_sonata('exc_stim_spikes.h5')
# exc_psg = PoissonSpikeGenerator(population='exc_stim')
# for i in range(num_exc):
#         exc_psg.add(node_ids=[i],  
#                 firing_rate=exc_frs[i]/1000,    
#                 times=(0.0*1000, seconds*1000))     
# exc_psg.to_sonata('exc_stim_spikes.h5')

def positive_normal(mean, std):
        return max(np.random.normal(loc=mean, scale=std), 0.01)

#Makes a spike raster with each cell having its own noise trace.
def gen_inh_spikes(n_cells, mean_fr, std_fr, key, file):
        node_ids = []
        timestamps = []

        length = times[1] - times[0]
        buffer = times[0]

        for i in range(n_cells):
                z = make_noise(num_samples=(int(length*1000))-1,num_traces=1)
                df = make_spikes(False, partial(positive_normal, mean=mean_fr, std=std_fr), numUnits=1,rateProf=z[0,:])
                node_ids = np.concatenate((node_ids, np.array(df["node_ids"]) + i))
                timestamps = np.concatenate((timestamps, np.array(df["timestamps"])))

        timestamps += buffer * 1000
        raster_to_sonata(node_ids, timestamps, key, file)

gen_inh_spikes(num_soma_inh + num_prox_dend_inh, prox_fr_mean, prox_fr_std, "prox_inh_stim", 'prox_inh_stim_spikes.h5')
gen_inh_spikes(num_apic_inh + num_dend_inh, dist_fr_mean, dist_fr_std, "dist_inh_stim", 'dist_inh_stim_spikes.h5')
# inh_psg = PoissonSpikeGenerator(population='prox_inh_stim')
# inh_psg.add(node_ids=range(num_soma_inh + num_prox_dend_inh), 
#         firing_rate=inh_fr,  
#         times=times)   
# inh_psg.to_sonata('prox_inh_stim_spikes.h5')

# inh_psg = PoissonSpikeGenerator(population='dist_inh_stim')
# inh_psg.add(node_ids=range(num_apic_inh + num_dend_inh), 
#         firing_rate=inh_fr,  
#         times=times)   
# inh_psg.to_sonata('dist_inh_stim_spikes.h5')

# from crop_raster import crop_raster
# crop_raster("rhythmic_inh_spikes.h5", 'inh_stim_spikes.h5', 120000, num_inh)


from bmtk.utils.sim_setup import build_env_bionet

try:
        build_env_bionet(base_dir='./',
                        network_dir='./network',
                        dt = 0.1, tstop=seconds * 1000.0,
                        report_vars=['v'],
                        dL = 5,
                        spikes_threshold=-10,
                        spikes_inputs=[('exc_stim', 'exc_stim_spikes.h5'), ('prox_inh_stim', 'prox_inh_stim_spikes.h5'), ('dist_inh_stim', 'dist_inh_stim_spikes.h5')],
                        components_dir='../biophys_components',
                        compile_mechanisms=True,
                        config_file="config.json")
except:
        build_env_bionet(base_dir='./',
                        network_dir='./network',
                        dt = 0.1, tstop=seconds * 1000.0,
                        report_vars=['v'],
                        dL = 5,
                        spikes_threshold=-10,
                        spikes_inputs=[('exc_stim', 'exc_stim_spikes.h5'), ('prox_inh_stim', 'prox_inh_stim_spikes.h5'), ('dist_inh_stim', 'dist_inh_stim_spikes.h5')],
                        components_dir='../biophys_components',
                        compile_mechanisms=True)
