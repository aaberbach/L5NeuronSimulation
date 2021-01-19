from bmtk.builder import NetworkBuilder
import numpy as np
import sys

np.random.seed(2129)

# if __name__ == '__main__':
#     if __file__ != sys.argv[-1]:
#         inp = sys.argv[-1]
#     else:
#         raise Exception("no work" + str(sys.argv[-1]))

# dist = float(inp)


# synapses.load()
# syn = synapses.syn_params_dicts()

net = NetworkBuilder("biophysical")

net.add_nodes(N=1, pop_name='Pyrc',
    potental='exc',
    model_type='biophysical',
    model_template='hoc:L5PCtemplate',
    morphology = None)

# Build and save our networks
net.build()
net.save_nodes(output_dir='network')
net.save_edges(output_dir='network')

from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',
                network_dir='./network',
                tstop=1000.0, dt = 0.1,
                report_vars=['v'],
                spikes_threshold=-10,
                current_clamp={           # Creates a step current from 500.ms to 1500.0 ms  
                     'amp': -0.05,
                     #'std': [0.0, 0.0],
                     'delay': 500,
                     'duration': 400,
                     'gids':"all"
                 },
                components_dir='../biophys_components',
                compile_mechanisms=True)