from bmtk.builder.networks import NetworkBuilder

net = NetworkBuilder("biophysical")
net.add_nodes(N=1, pop_name='Pyrc',
    potental='exc',
    model_type='biophysical',
    model_template='hoc:L5PCtemplate',
    morphology = None)

net.build()
net.save_nodes(output_dir='network')
    
from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',      # Where to save the scripts and config files 
                 components_dir='../biophys_components',
                 network_dir='./network',    # Location of directory containing network files
                 tstop=3000.0, dt=0.1,     # Run a simulation for 2000 ms at 0.1 ms intervals
                 report_vars=['v'],
                 dL = 5,
                 #clamp_reports=["se"], # Tells simulator we want to record membrane potential and calcium traces
                 spikes_threshold=-10,
                 compile_mechanisms=True   # Will try to compile NEURON mechanisms
                )
