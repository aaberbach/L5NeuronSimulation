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
from functools import partial
import run
import matplotlib.pyplot as plt

def store_seg_voltage(graph, sim, sec_type : str, sec_id : int, sec_x : float) -> h.Vector:
    cells = graph.get_local_cells()
    cell = cells[list(cells.keys())[0]]
    #import pdb; pdb.set_trace()
    seg_volt = h.Vector()
    seg_volt.record(getattr(cell.hobj, sec_type)[sec_id](sec_x)._ref_v)

    return seg_volt

def save_inputs(graph, sim):
    cells = graph.get_local_cells()
    cell = cells[list(cells.keys())[0]]

    h.distance(sec=cell.hobj.soma[0])

    sec_types = []
    weights = []
    dists = []
    names = []

    cons = cell.connections()
    for c in cons:
        con = c._connector
        syn = con.syn()
        weights.append(float(syn.initW))

        seg = con.postseg()
        fullsecname = seg.sec.name()
        names.append(fullsecname)
        #import pdb; pdb.set_trace()
        dists.append(float(h.distance(seg)))
        sec_types.append(fullsecname.split(".")[1][:4])

    df = pd.DataFrame()
    df["Distance"] = dists
    df["Conductance"] = weights
    df["Type"] = sec_types
    df["Name"] = names

    df.to_csv("Inputs.csv")

if __name__ == "__main__":
    if __file__ != sys.argv[-1] and __file__ != sys.argv[-2]:
        m = float(sys.argv[-2])
        s = float(sys.argv[-1])
    else:
        raise Exception("The inputted command argument not functional: " + str(sys.argv[-1]) + " " + str(sys.argv[-2]))

    synapses.set_pyr_w(m, s)
    synapses.load()
    syn = synapses.syn_params_dicts()

    store_nexus_volt = partial(store_seg_voltage, sec_type = "apic", sec_id = 36, sec_x = 1)

    callback_returns = run.run_network([store_nexus_volt, save_inputs], v_report_all = False, quit_execution=False)

    seg_volt = callback_returns[0]
    plt.plot(np.array(seg_volt))
    plt.show()
    #import pdb; pdb.set_trace()