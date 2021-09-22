from bmtk.builder import NetworkBuilder
from bmtk.simulator import bionet
import numpy as np
import sys
import synapses   
import pandas as pd                                                                                                             
import h5py
import matplotlib.pyplot as plt
synapses.load()
syn = synapses.syn_params_dicts()

np.random.seed(2129)

net = NetworkBuilder("biophysical")

#Whether you want exc or inh input.
excitatory = True

##################################################################################
###################################Pyr Type C#####################################

net.add_nodes(N=1, pop_name='Pyrc',
    potental='exc',
    model_type='biophysical',
    model_template='hoc:L5PCtemplate',
    morphology = None)

##################################################################################
###################################External Networks##############################

# External excitatory inputs
exc_stim = NetworkBuilder('exc_stim')
exc_stim.add_nodes(N=1,
                pop_name='exc_stim',
                potential='exc',
                model_type='virtual')

# External inhibitory inputs
inh_stim = NetworkBuilder('inh_stim')
inh_stim.add_nodes(N=1,
                pop_name='inh_stim',
                potential='exc',
                model_type='virtual')

##################################################################################
###################################Edges##########################################

#Inhibitory.
net.add_edges(source=inh_stim.nodes(), target=net.nodes(),
                connection_rule=1,
                syn_weight=1,
                delay=0.1,
                dynamics_params='PV2PN.json',
                model_template=syn['PV2PN.json']['level_of_detail'],
                distance_range=[100, 2000.0],
                target_sections=['dend'])
                #sec_id=whatever you want,
                #sec_x=whateber you want)

# Create connections between Exc --> Pyr cells

#Excitatory.
net.add_edges(source=exc_stim.nodes(), target=net.nodes(),
                connection_rule=1,
                syn_weight=1,
                delay=0.1,
                #target_sections=['dend'],
                #distance_range=[50.0, 2000.0],
                sec_id=126,
                sec_x=0.5,
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

#SPIKE TRAIN

spacing = 20#10
num_AP = 10#13
node_ids = np.full(num_AP, 0)
timestamps = np.linspace(0 + 500, (num_AP-1)*spacing + 500, num_AP)
#timestamps *= 10
#import pdb; pdb.set_trace()
import h5py

if excitatory:
    key = "exc_stim"
else:
    key = 'inh_stim'
    
f = h5py.File('stim_spikes.h5', 'w')
group = f.create_group('spikes')
group = group.create_group(key)
group.create_dataset("node_ids", data = node_ids)
group.create_dataset("timestamps", data = timestamps)
f.close()

from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',
                network_dir='./network',
                dt = 0.1, tstop=1000.0,
                report_vars=['v', 'cai'],
                dL=5,
                # current_clamp={           # Creates a step current from 500.ms to 1500.0 ms  
                #      'amp': 0.2,
                #      #'std': [0.0, 0.0],
                #      'delay': 300,
                #      'duration': 2000,
                #      'gids':"all"
                #  },
                # clamp_reports=['se'],#Records se clamp currents.
                # se_voltage_clamp={
                #         "amps":[[0, 0, 0]],
                #         "durations": [[120000, 0, 0]],
                #         'gids': [0],
                #         'rs': [0.01],
                # },
                spikes_threshold=-10,
                spikes_inputs=[(key, 'stim_spikes.h5')],
                components_dir='../biophys_components',
                compile_mechanisms=True)

################ RUN #############################
"""Script for running the network built in build_network.py

Also saves a file called Connections.csv that consists of information about
each synapse in the simulation.
"""

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 
sys.path.insert(0, currentdir)

from bmtk.simulator import bionet
import numpy as np
from neuron import h
import synapses 
import pandas as pd
import run

try:
    np.random.seed(int(sys.argv[1]))
except:
    np.random.seed(123)

def splitcell(graph, sim):
    pc = h.ParallelContext()  # object to access MPI methods
    MPI_size = int(pc.nhost())
    MPI_rank = int(pc.id())

    h.load_file("netparmpi.hoc")
    pnm = h.ParallelNetManager(1)

    cells = graph.get_local_cells()
    cell = cells[list(cells.keys())[0]]

    if MPI_rank == 0:
        pnm.splitcell(0, 1, sec=cell.hobj.apic[50])
        # cells = graph.get_local_cells()
        # cell = cells[list(cells.keys())[0]]
        # import pdb; pdb.set_trace()
    else:
        pnm.splitcell(1, 0, sec=cell.hobj.apic[50])

def reduce_reports(graph, sim, percent = 0.1):
    """Reduces the number of segments whose variables are saved. 

    Parameters
    ----------
    graph : BioNetwork
        the BMTK network
    sim : BioSimulator
        the BMTK simulation that contains the MembraneReport
    percent : float, optional
        proportion of segments to save, by default 0.1
    """    
    recorder = None

    #Finds the correct recorder.
    for mod in sim._sim_mods:
        if type(mod) == bionet.modules.record_cellvars.MembraneReport:
            recorder = mod
    
    if recorder == None:
        raise Exception("For reduce_reports to be called, there must be a MembraneReport in the simulation.")

    raise Exception("This has not been implemented yet!")

def save_connections(graph, sim):
    """Saves Connections.csv based on the given network.

    Parameters
    ----------
    graph : BioNetwork
        the network that the connections are retrieved from
    sim : BioSimulator
        the simulation about to be run (not used in this function)
    """    
    cells = graph.get_local_cells()
    cell = cells[list(cells.keys())[0]]

    h.distance(sec=cell.hobj.soma[0])#Makes the distances correct.

    sec_types = []#soma, apic, or dend
    weights = []#scaled conductances (initW)
    dists = []#distance from soma
    node_ids = []#node id within respective node population (exc, prox_inh, dist_inh)
    names = []#full NEURON str representation of postsynaptic segment
    source_pops = []#node population
    release_probs = []#propability of release.

    for c in cell.connections():
        con = c._connector
        source = c.source_node
        syn = con.syn()
        seg = con.postseg()
        fullsecname = seg.sec.name()

        source_pops.append(source._population)
        node_ids.append(source._node_id)

        weights.append(float(syn.initW))
        release_probs.append(float(syn.P_0))
        names.append(str(seg))
        sec_types.append(fullsecname.split(".")[1][:4])
        dists.append(float(h.distance(seg)))

    df = pd.DataFrame()
    df["Node ID"] = node_ids
    df["Distance"] = dists
    df["Conductance"] = weights
    df["Type"] = sec_types
    df["Name"] = names
    df["Source Population"] = source_pops
    df["Release Probability"] = release_probs
    df.to_csv("Connections.csv", index=False)


run.run_network([save_connections], v_report_all = True)#make v_report_all True to save all segments

