from bmtk.builder import NetworkBuilder
import numpy as np
import sys
import synapses   
import pandas as pd                                                                                                             
synapses.load()
syn = synapses.syn_params_dicts()

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        inp = sys.argv[-1]
    else:
        raise Exception("no work" + str(sys.argv[-1]))

N = int(inp)

# exc_fr = float(frs[0])
# inh_fr = float(frs[1])
# print(exc_fr, inh_fr)

#N = 10

# Initialize our network

net = NetworkBuilder("biophysical")

def lognormal(m, s):
        mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
        std = np.sqrt(np.log((s/m)**2 + 1))
        return max(np.random.lognormal(mean, std, 1), 0)

# Dend Excitatory: 8998.0
# Apic Excitatory: 12380.0
# Soma Inhibitory: 2285.0

# num_inh = [231, 405, 61]#[int(lognormal(56, 7.5)) for i in range(N)]
# print(num_inh)
# inh_bounds = []
# sum = 0
# for num in num_inh:
#         sum += num
#         inh_bounds.append(sum)

# num_inh_soma = [int(lognormal(61, 9)) for i in range(N)]
# print(num_inh_soma)
# inh_bounds_soma = []
# sum = 0
# for num in num_inh_soma:
#         sum += num
#         inh_bounds_soma.append(sum)

##############################Ecitatory###############################333

# num_exc = [1590, 2785]#[int(lognormal(1000, 80)) for i in range(N)]
# print(num_exc)
# exc_bounds = []
# sum = 0
# for num in num_exc:
#         sum += num
#         exc_bounds.append(sum)

# num_exc_apic = [int(lognormal(2785, 350)) for i in range(N)]
# print(num_exc_apic)
# exc_bounds_apic = []
# sum = 0
# for num in num_exc_apic:
#         sum += num
#         exc_bounds_apic.append(sum)

# num_exc_dend = [int(lognormal(1590, 200)) for i in range(N)]
# print(num_exc_dend)
# exc_bounds_dend = []
# sum = 0
# for num in num_exc_dend:
#         sum += num
#         exc_bounds_dend.append(sum)

scale_div =10

num_dend_exc = 7186 // scale_div
num_apic_exc = 10417 // scale_div
num_soma_inh = 1908 // scale_div

exc_fr_mean = 0.1
exc_fr_std = 0.5
inh_fr = 7 #* scale_div

##################################################################################
###################################Pyr Type C#####################################

net.add_nodes(N=N, pop_name='Pyrc',
    potental='exc',
    model_type='biophysical',
    model_template='hoc:L5PCtemplate',
    morphology = None)

# net.add_nodes(N=N, pop_name='Pyrc',
#     potental='exc',
#     model_type='biophysical',
#     model_template='ctdb:Biophys1.hoc',
#     model_processing='aibs_allactive',
#     #model_processing='my_allactive',
#     dynamics_params='L5Conductances.json',
#     #dynamics_params='491766131_fit.json',
#     morphology='Rbp4-Cre_KL100_Ai14-203503.04.01.01_527109145_m.swc')
#     #morphology='Rbp4-Cre_KL100_Ai14-203503.04.01.01_527109145_m.swc')


##################################################################################
###################################External Networks##############################

#print("Internal nodes built")

num_exc = (num_apic_exc + num_dend_exc) * N

# External excitatory inputs
exc_stim = NetworkBuilder('exc_stim')
exc_stim.add_nodes(N=num_exc,
                pop_name='exc_stim',
                potential='exc',
                model_type='virtual')

num_inh = num_soma_inh * N

# External inhibitory inputs
inh_stim = NetworkBuilder('inh_stim')
inh_stim.add_nodes(N=num_inh,
                pop_name='inh_stim',
                potential='exc',
                model_type='virtual')

##################################################################################
###################################Edges##########################################

def correct_cell(source, target, num_per, start):
        #import pdb; pdb.set_trace()
        sid = source.node_id - start
        tid = target.node_id
        #import pdb; pdb.set_trace()

        lower_bound = num_per * tid

        upper_bound = lower_bound + num_per

        if sid < upper_bound and sid >= lower_bound:
                #print("connecting cell {} to {}".format(sid+start,tid))
                return 1
        else:
                return None

#Inhibitory on soma.
net.add_edges(source=inh_stim.nodes(), target=net.nodes(),
                connection_rule=correct_cell,
                connection_params={'num_per': num_soma_inh , 'start':0},
                syn_weight=1,
                delay=0.1,
                dynamics_params='INT2PN.json',
                model_template=syn['INT2PN.json']['level_of_detail'],
                distance_range=[-2000.0, 2000.0],
                target_sections=['somatic'])

#start += np.sum(inh_bounds_apic)
start = 0

# Create connections between Exc --> Pyr cells

#Excitatory on basal dendrites.
net.add_edges(source=exc_stim.nodes(), target=net.nodes(),
                connection_rule=correct_cell,
                connection_params={'num_per': num_dend_exc, 'start':start},
                syn_weight=1,
                target_sections=['dend'],
                delay=0.1,
                distance_range=[-2000.0, 2000.0],
                dynamics_params='PN2PN.json',
                model_template=syn['PN2PN.json']['level_of_detail'])

#start += np.sum(exc_bounds_dend)
start += np.sum(num_dend_exc * N)

#Excitatory on apical dendrites.
net.add_edges(source=exc_stim.nodes(), target=net.nodes(),
                connection_rule=correct_cell,
                connection_params={'num_per': num_apic_exc, 'start':start},
                syn_weight=1,
                target_sections=['apic'],
                delay=0.1,
                distance_range=[-2000.0, 2000.0],
                dynamics_params='PN2PN.json',
                model_template=syn['PN2PN.json']['level_of_detail'])

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

exc_dend_frs = []
exc_apic_frs = []

exc_means = []
exc_stds = []
exc_maxs = []

for i in range(N):
        fr_mean = np.random.uniform(exc_fr_mean - 0.015, exc_fr_mean + 0.03, 1)
        #fr_mean = np.random.uniform(exc_fr_mean + 0.5, exc_fr_mean + 0.55, 1)
        dend_frs = [min(float(lognormal(fr_mean, exc_fr_std)), fr_mean+8*exc_fr_std) for _ in range(num_dend_exc)]
        apic_frs = [min(float(lognormal(fr_mean, exc_fr_std)), fr_mean+8*exc_fr_std) for _ in range(num_apic_exc)]

        frs = dend_frs + apic_frs

        exc_dend_frs += dend_frs
        exc_apic_frs += apic_frs

        exc_means.append(np.mean(frs))
        exc_stds.append(np.std(frs))
        exc_maxs.append(np.max(frs))

exc_frs = exc_dend_frs + exc_apic_frs
fr_df = pd.DataFrame()
#fr_df['gid'] = [i + num_exc+num_inh for i in range(N)]
fr_df['gid'] = [i for i in range(N)]
fr_df['fr_mean'] = exc_means
fr_df['fr_std'] = exc_stds
fr_df['fr_max'] = exc_maxs

fr_df.to_csv('frs_temp.csv', index=False)

exc_psg = PoissonSpikeGenerator(population='exc_stim')
for i in range(num_exc):
        exc_psg.add(node_ids=[i],  
                firing_rate=float(exc_frs[i]/1000),    
                times=(0.2*1000, 5.2*1000))     
exc_psg.to_sonata('exc_stim_spikes.h5')

inh_psg = PoissonSpikeGenerator(population='inh_stim')
inh_psg.add(node_ids=range(num_inh), 
        firing_rate=inh_fr/1000,  
        times=(0.2*1000, 5.2*1000))   
inh_psg.to_sonata('inh_stim_spikes.h5')


from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',
                network_dir='./network',
                tstop=5200.0, dt = 0.1,
                report_vars=['v', 'cai'],
                # current_clamp={           # Creates a step current from 500.ms to 1500.0 ms  
                #      'amp': 0.793,
                #      #'std': [0.0, 0.0],
                #      'delay': 700,
                #      'duration': 2000,
                #      'gids':"0"
                #  },
                spikes_threshold=-10,
                spikes_inputs=[('exc_stim', 'exc_stim_spikes.h5'), ('inh_stim', 'inh_stim_spikes.h5')],
                components_dir='../biophys_components',
                compile_mechanisms=True)