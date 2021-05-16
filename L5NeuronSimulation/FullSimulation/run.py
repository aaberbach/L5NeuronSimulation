"""Contains the general function run_network used to run a bmtk simulation.
"""

from bmtk.simulator import bionet
import numpy as np
from neuron import h
import pandas as pd

def run_network(callbacks):
    """Runs the standard bmtk simulation and call the given callbacks right before running the simulation.

    Parameters
    ----------
    callbacks : list
        list of functions to be called before sim.run()
        each function will be called with (graph, sim)
    """    
    np.random.seed(42)

    config_file = 'simulation_config.json'
    conf = bionet.Config.from_json(config_file, validate=True)
    conf.build_env()

    graph = bionet.BioNetwork.from_config(conf)
    sim = bionet.BioSimulator.from_config(conf, network=graph)

    for c in callbacks:
        c(graph, sim)

    sim.run()