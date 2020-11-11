from bmtk.builder import NetworkBuilder
import numpy as np
import sys
import synapses                                                                                                                
synapses.load()
syn = synapses.syn_params_dicts()

# if __name__ == '__main__':
#     if __file__ != sys.argv[-1]:
#         inp = sys.argv[-1]
#     else:
#         raise Exception("no work" + str(sys.argv[-1]))

# frs = str(inp).split(",")

# exc_fr = float(frs[0])
# inh_fr = float(frs[1])
# print(exc_fr, inh_fr)

N = 1

# Initialize our network

net = NetworkBuilder("biophysical")

def lognormal(m, s):
        mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
        std = np.sqrt(np.log((s/m)**2 + 1))
        return max(np.random.lognormal(mean, std, 1), 0)

# Dend Excitatory: 8998.0
# Apic Excitatory: 12380.0
# Soma Inhibitory: 2285.0

scale_div = 10

num_dend_exc = 8998 // scale_div
num_apic_exc = 12380 // scale_div
num_soma_inh = 2285 // scale_div

exc_fr_mean = 0.1
exc_fr_std = 0.5
inh_fr = 7

##################################################################################
###################################Pyr Type C#####################################

net.add_nodes(N=1, pop_name='Pyrc',
    potental='exc',
    model_type='biophysical',
    model_template='hoc:L5PCtemplate',
    morphology = None)


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
                #print("connecting cell {} to {}".format(sid,tid))
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

exc_dend_frs = [lognormal(exc_fr_mean, exc_fr_std) for _ in range(num_dend_exc*N)]
exc_apic_frs = [lognormal(exc_fr_mean, exc_fr_std) for _ in range(num_apic_exc*N)]

exc_frs = exc_dend_frs + exc_apic_frs

from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator
from bmtk.utils.reports.spike_trains.spikes_file_writers import write_csv

exc_psg = PoissonSpikeGenerator(population='exc_stim')
exc_psg.add(node_ids=range(num_exc),  
        firing_rate=exc_frs,    
        times=(0.2, 5.2))    
exc_psg.to_sonata('exc_stim_spikes.h5')

inh_psg = PoissonSpikeGenerator(population='inh_stim')
inh_psg.add(node_ids=range(num_inh), 
        firing_rate=inh_fr,  
        times=(0.2, 5.2))   
inh_psg.to_sonata('inh_stim_spikes.h5')


from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',
                network_dir='./network',
                tstop=1200.0, dt = 0.1,
                report_vars=['v'],
                spikes_inputs=[('exc_stim', 'exc_stim_spikes.h5'), ('inh_stim', 'inh_stim_spikes.h5')],
                components_dir='biophys_components',
                compile_mechanisms=True)
