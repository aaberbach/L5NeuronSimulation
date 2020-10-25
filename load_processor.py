from neuron import h
import io
from bmtk.simulator.bionet.pyfunction_cache import add_cell_processor
from bmtk.simulator.bionet.default_setters.cell_models import set_params_peri, set_params_allactive

def add_params_allactive(hobj, params_dict):
    # params_dict = json.load(open(params_file_name, 'r'))
    passive = params_dict['passive'][0]
    genome = params_dict['genome']
    conditions = params_dict['conditions'][0]

    section_map = {}
    for sec in hobj.all:
        section_name = sec.name().split(".")[1][:4]
        if section_name in section_map:
            section_map[section_name].append(sec)
        else:
            section_map[section_name] = [sec]

    for sec in hobj.all:
        sec.insert('pas')
        # sec.insert('extracellular')

    if 'e_pas' in passive:
        e_pas_val = passive['e_pas']
        for sec in hobj.all:
            for seg in sec:
                seg.pas.e = e_pas_val

    if 'g_pas' in passive:
        for g_pas_dict in passive['g_pas']:
            g_pas_val = g_pas_dict['g_pas']
            for sec in section_map.get(g_pas_dict['section'], []):
                for seg in sec:
                    seg.pas.g = g_pas_val

    if 'ra' in passive:
        ra_val = passive['ra']
        for sec in hobj.all:
            sec.Ra = ra_val

    if 'cm' in passive:
        # print('Setting cm')
        for cm_dict in passive['cm']:
            cm = cm_dict['cm']
            for sec in section_map.get(cm_dict['section'], []):
                sec.cm = cm

    for genome_dict in genome:
        g_section = genome_dict['section']
        if genome_dict['section'] == 'glob':
            io.log_warning("There is a section called glob, probably old json file")
            continue

        g_value = float(genome_dict['value'])
        g_name = genome_dict['name']
        g_mechanism = genome_dict.get("mechanism", "")
        for sec in section_map.get(g_section, []):
            if g_mechanism != "":
                sec.insert(g_mechanism)
            setattr(sec, g_name, g_value)

    for erev in conditions['erev']:
        erev_section = erev['section']
        erev_ena = erev['ena']
        erev_ek = erev['ek']

        if erev_section in section_map:
            for sec in section_map.get(erev_section, []):
                if h.ismembrane('k_ion', sec=sec) == 1:
                    setattr(sec, 'ek', erev_ek)
                if h.ismembrane('na_ion', sec=sec) == 1:
                    setattr(sec, 'ena', erev_ena)
        else:
            io.log_warning("Can't set erev for {}, section array doesn't exist".format(erev_section))

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
        #set_params_allactive(hobj, dynamics_params)
        add_params_allactive(hobj, dynamics_params)

    return hobj

def load():
    add_cell_processor(my_allactive, overwrite=False)