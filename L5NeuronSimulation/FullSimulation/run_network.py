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
from bmtk.simulator.bionet import modules as mods
from bmtk.simulator.bionet.seclamp import SEClamp
import numpy as np
from neuron import h
import synapses 
import pandas as pd
import run
from functools import partial

try:
    np.random.seed(int(sys.argv[1]))
except:
    np.random.seed(123)

def vclamp_seg(seg, durs, amps, rs=None):
    clamp = SEClamp(amps, durs, rs=rs)

    clamp._stim = h.SEClamp(seg) 
    clamp._stim.dur1 = durs[0]
    clamp._stim.dur2 = durs[1]
    clamp._stim.dur3 = durs[2]

    clamp._stim.amp1 = amps[0]
    clamp._stim.amp2 = amps[1]
    clamp._stim.amp3 = amps[2]

    if rs != None:
        clamp._stim.rs = rs

    return clamp

def vclamp_all_segs(graph, sim, durs, amps, rs=None):
    cells = graph.get_local_cells()
    cell = cells[list(cells.keys())[0]]
    hobj = cell.hobj

    vclamps = []

    for sec in hobj.all:
        for seg in sec:
            vclamps.append(vclamp_seg(seg, durs, amps, rs))

    return vclamps

def record_all_vclamps(graph, sim, durs, amps, rs=None):
    vclamps = vclamp_all_segs(graph, sim, durs, amps, rs)

    clamp_params = {}
    clamp_params["input_type"] = "voltage_clamp"
    clamp_params["module"] = "SEClamp"
    clamp_params["node_set"] = "all"

    mod = mods.ClampReport(file_name = currentdir + "/output/se_clamp_report.h5", tmp_dir = currentdir + "/output/", variable_name = "se", **clamp_params)
    sim.add_mod(mod)

    for clamp in vclamps:
        sim._seclamps.append(clamp)


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


if __name__ == "__main__":
    synapses.load()
    syn = synapses.syn_params_dicts()

    callbacks = []
    
    save_cons = True
    vclamp_all = False

    if (save_connections):
        callbacks.append(save_connections)

    if (vclamp_all):
        save_epscs = partial(record_all_vclamps, durs = [1000000, 0, 0], amps = [0, 0, 0], rs=0.01)
        callbacks.append(save_epscs)

    callback_returns = run.run_network(callbacks, v_report_all = False)#make v_report_all True to save all segments
