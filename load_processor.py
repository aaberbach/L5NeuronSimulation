from neuron import h
from bmtk.simulator.bionet.pyfunction_cache import add_cell_processor
from bmtk.simulator.bionet.default_setters.cell_models import set_params_peri, set_params_allactive

def add_axon(hobj):
    """Replace reconstructed axon with a stub

    :param hobj: hoc object
    """

    h.execute('create axon[2]', hobj)

    for sec in hobj.axon:
        sec.L = 30
        sec.diam = 1
        hobj.axonal.append(sec=sec)
        hobj.all.append(sec=sec)  # need to remove this comment

    hobj.axon[0].connect(hobj.soma[0], 0.5, 0)
    hobj.axon[1].connect(hobj.axon[0], 1, 0)

    h.define_shape()

def my_allactive(hobj, cell, dynamics_params):
    if dynamics_params is not None:
        add_axon(hobj)
        #set_params_peri(hobj, dynamics_params)
        set_params_allactive(hobj, dynamics_params)

    return hobj

def load():
    add_cell_processor(my_allactive, overwrite=False)