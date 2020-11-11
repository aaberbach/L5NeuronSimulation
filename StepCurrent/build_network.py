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
                 #clamp_reports=["se"], # Tells simulator we want to record membrane potential and calcium traces
                 current_clamp={           # Creates a step current from 500.ms to 1500.0 ms  
                     'amp': 0.793,
                     #'std': [0.0, 0.0],
                     'delay': 700,
                     'duration': 2000,
                     'gids':"all"
                 },
                 spikes_threshold=-10,
                #  file_current_clamp={
                #      "input_file": "PN_IClamp/inputs/amps.h5"
                #  },
                #  se_voltage_clamp={
                #      "amps":[[-70, -70, -70], [-70, -70, -70]],
                #      "durations": [[2000, 2000, 2000], [2000, 2000, 2000]],
                #      'gids': [0, 1],
                #      'rs': [0.001, 0.01],
                #      'name':"PN_se_clamp"
                #  },
                 compile_mechanisms=True   # Will try to compile NEURON mechanisms
                )