"""Contains the general function run_network used to run a bmtk simulation.
"""
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
#sys.path.insert(0, parentdir) 
sys.path.insert(0, currentdir)

from bmtk.simulator import bionet
import numpy as np
from neuron import h
import pandas as pd

def run_network(callbacks, v_report_all=False, quit_execution=True):
    """Runs the standard bmtk simulation and call the given callbacks right before running the simulation.

    Parameters
    ----------
    callbacks : list
        list of functions to be called before sim.run().
        each function will be called with (graph, sim)
    v_report_all : bool
        whether the v_report should be set to record every section
    quit_execution : bool
        whether to call bionet.nrn.quit_execution()
    return: list
        list of returns of the given callbacks
    """    
    np.random.seed(42)

    config_file = 'simulation_config.json'
    conf = bionet.Config.from_json(config_file, validate=True)

    if v_report_all:
        conf["reports"]["v_report"]["sections"] = "all"

    conf.build_env()

    graph = bionet.BioNetwork.from_config(conf)
    sim = bionet.BioSimulator.from_config(conf, network=graph)

    callback_returns = []

    for c in callbacks:
        callback_returns.append(c(graph, sim))

    sim.run()

    if (quit_execution):
        bionet.nrn.quit_execution()

    return callback_returns
