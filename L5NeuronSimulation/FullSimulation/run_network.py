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


if __name__ == "__main__":
    synapses.load()
    syn = synapses.syn_params_dicts()

    np.random.seed(42)

    run.run_network([save_connections], v_report_all = False)#make v_report_all True to save all segments