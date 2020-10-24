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

#apical area = 5064 -> 2785 exc, 405 inh
#basal area = 2891 -> 1590 exc, 231 inh
#soma area = 608 -> 61 inh

num_inh = [231, 405, 61]#[int(lognormal(56, 7.5)) for i in range(N)]
print(num_inh)
inh_bounds = []
sum = 0
for num in num_inh:
        sum += num
        inh_bounds.append(sum)

num_exc = [1590, 2785]#[int(lognormal(1000, 80)) for i in range(N)]
print(num_exc)
exc_bounds = []
sum = 0
for num in num_exc:
        sum += num
        exc_bounds.append(sum)

exc_fr = 2.8
inh_fr = 10

##################################################################################
###################################Pyr Type C#####################################

net.add_nodes(N=1, pop_name='Pyrc',
    potental='exc',
    model_type='biophysical',
    model_template='ctdb:Biophys1.hoc',
    #model_processing='aibs_allactive',
    model_processing='my_allactive',
    dynamics_params='L5Conductances.json',
    #dynamics_params='491766131_fit.json',
    morphology='L5Morphology.swc')
    #morphology='Rbp4-Cre_KL100_Ai14-203503.04.01.01_527109145_m.swc')


##################################################################################
###################################External Networks##############################

#print("Internal nodes built")

# External excitatory inputs
exc_stim = NetworkBuilder('exc_stim')
exc_stim.add_nodes(N=np.sum(num_exc),
                pop_name='exc_stim',
                potential='exc',
                model_type='virtual')

# External inhibitory inputs
inh_stim = NetworkBuilder('inh_stim')
inh_stim.add_nodes(N=np.sum(num_inh),
                pop_name='inh_stim',
                potential='exc',
                model_type='virtual')

##################################################################################
###################################Edges##########################################

def correct_cell(source, target, bounds, part):
        sid = source.node_id

        lower_bound = 0
        if part > 0:
                lower_bound = bounds[part - 1]

        upper_bound = bounds[part]

        if sid < upper_bound and sid >= lower_bound:
                #print("connecting cell {} to {}".format(sid,tid))
                return 1
        else:
                return None

#Create connections between Inh --> Pyr cells
# net.add_edges(source=inh_stim.nodes(), target=net.nodes(),
#        connection_rule=correct_cell,
#        connection_params={'bounds': inh_bounds},
#        syn_weight=1e-3,
#        #weight_function='lognormal',
#        weight_sigma=1e-3,
#        weight_max=20e-3,
#        dynamics_params='GABA_InhToExc.json',
#        model_template='Exp2Syn',
#        distance_range=[-1000.0, 1000.0],
#        target_sections=['somatic'],
#        delay=2.0)
# Create connections between Exc --> Pyr cells

#Inhibitory on soma.
net.add_edges(source=inh_stim.nodes(), target=net.nodes(),
                connection_rule=correct_cell,
                connection_params={'bounds': inh_bounds, 'part':2},
                syn_weight=1,
                delay=0.1,
                dynamics_params='INT2PN.json',
                model_template=syn['INT2PN.json']['level_of_detail'],
                distance_range=[-1000.0, 1000.0],
                target_sections=['somatic'])

#Inhibitory on basal dendrites.
net.add_edges(source=inh_stim.nodes(), target=net.nodes(),
                connection_rule=correct_cell,
                connection_params={'bounds': inh_bounds, 'part':0},
                syn_weight=1,
                delay=0.1,
                dynamics_params='INT2PN.json',
                model_template=syn['INT2PN.json']['level_of_detail'],
                distance_range=[-1000.0, 1000.0],
                target_sections=['dend'])

#Inhibitory on apical dendrites.
net.add_edges(source=inh_stim.nodes(), target=net.nodes(),
                connection_rule=correct_cell,
                connection_params={'bounds': inh_bounds, 'part':1},
                syn_weight=1,
                delay=0.1,
                dynamics_params='INT2PN.json',
                model_template=syn['INT2PN.json']['level_of_detail'],
                distance_range=[-1000.0, 1000.0],
                target_sections=['apic'])

# Create connections between Exc --> Pyr cells

#Excitatory on basal dendrites.
net.add_edges(source=exc_stim.nodes(), target=net.nodes(),
                connection_rule=correct_cell,
                connection_params={'bounds': exc_bounds, 'part':0},
                syn_weight=1,
                target_sections=['dend'],
                delay=0.1,
                distance_range=[-1000.0, 1000.0],
                dynamics_params='PN2PN.json',
                model_template=syn['PN2PN.json']['level_of_detail'])

#Excitatory on apical dendrites.
net.add_edges(source=exc_stim.nodes(), target=net.nodes(),
                connection_rule=correct_cell,
                connection_params={'bounds': exc_bounds, 'part':1},
                syn_weight=1,
                target_sections=['apic'],
                delay=0.1,
                distance_range=[-1000.0, 1000.0],
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

exc_psg = PoissonSpikeGenerator(population='exc_stim')
exc_psg.add(node_ids=range(np.sum(num_exc)),  
        firing_rate=exc_fr,    
        times=(0.2, 5.2))    
exc_psg.to_sonata('exc_stim_spikes.h5')

inh_psg = PoissonSpikeGenerator(population='inh_stim')
inh_psg.add(node_ids=range(np.sum(num_inh)), 
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
