from bmtk.builder import NetworkBuilder
import numpy as np
import sys

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        inp = sys.argv[-1]
    else:
        raise Exception("no work" + str(sys.argv[-1]))

N = int(inp)

# Initialize our network

net = NetworkBuilder("biophysical")

def lognormal(m, s):
        mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
        std = np.sqrt(np.log((s/m)**2 + 1))
        return np.random.lognormal(mean, std, 1)

num_inh = [int(lognormal(50, 13)) for i in range(N)]
print(num_inh)
inh_bounds = []
sum = 0
for num in num_inh:
        sum += num
        inh_bounds.append(sum)

num_exc = [int(lognormal(35, 10)) for i in range(N)]
print(num_exc)
exc_bounds = []
sum = 0
for num in num_exc:
        sum += num
        exc_bounds.append(sum)

exc_fr = 2
inh_fr = 10

##################################################################################
###################################Pyr Type C#####################################
        
net.add_nodes(N=N, pop_name='PyrC',
              potental='exc',
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_perisomatic',
              dynamics_params='472363762_fit.json',
              morphology='Scnn1a_473845048_m.swc')


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

def correct_cell(source, target, bounds):
        sid = source.node_id
        tid = target.node_id

        lower_bound = 0
        if tid > 0:
                lower_bound = bounds[tid - 1]

        upper_bound = bounds[tid]

        if sid < upper_bound and sid >= lower_bound:
                #print("connecting cell {} to {}".format(sid,tid))
                return 1
        else:
                return None

# Create connections between Inh --> Pyr cells
net.add_edges(source=inh_stim.nodes(), target=net.nodes(),
        connection_rule=correct_cell,
        connection_params={'bounds': inh_bounds},
        syn_weight=10.0e-03,
        weight_function='lognormal',
        weight_sigma=1.0e-03,
        weight_max=40e-03,
        dynamics_params='GABA_InhToExc.json',
        model_template='Exp2Syn',
        distance_range=[0.0, 300.0],
        target_sections=['somatic'],
        delay=2.0)

# Create connections between Exc --> Pyr cells
net.add_edges(source=exc_stim.nodes(), target=net.nodes(),
                connection_rule=correct_cell,
                connection_params={'bounds': exc_bounds},
                syn_weight=10.0e-03,
                weight_function='gaussianBL',
                weight_sigma=1.0e-03,
                weight_max=30e-03,
                target_sections=['dend'],
                delay=2.0,
                distance_range=[0.0, 300.0],
                dynamics_params='AMPA_ExcToExc.json',
                model_template='Exp2Syn')


# Build and save our networks
net.build()
net.save_nodes(output_dir='network')
net.save_edges(output_dir='network')

exc_stim.build()
exc_stim.save_nodes(output_dir='network')

inh_stim.build()
inh_stim.save_nodes(output_dir='network')


from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator

exc_psg = PoissonSpikeGenerator(population='exc_stim')
exc_psg.add(node_ids=range(np.sum(num_exc)),  
        firing_rate=int(exc_fr) / 1000,    
        times=(200.0, 1200.0))    
exc_psg.to_sonata('exc_stim_spikes.h5')

inh_psg = PoissonSpikeGenerator(population='inh_stim')
inh_psg.add(node_ids=range(np.sum(num_inh)), 
        firing_rate=int(inh_fr) / 1000,  
        times=(200.0, 1200.0))   
inh_psg.to_sonata('inh_stim_spikes.h5')


from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',
                network_dir='./network',
                tstop=1200.0, dt = 0.1,
                report_vars=['v'],
                spikes_inputs=[('exc_stim', 'exc_stim_spikes.h5'), ('inh_stim', 'inh_stim_spikes.h5')],
                components_dir='biophys_components',
                include_examples=True,
                compile_mechanisms=True)
