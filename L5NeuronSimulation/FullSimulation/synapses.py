"""synapses.py"""
import glob
import json
import os

from bmtk.simulator.bionet.pyfunction_cache import add_synapse_model
from neuron import h
import random
import numpy as np

np.random.seed(42)
generators = []

pyrWeight_m = 0.45#0.229#0.24575#0.95
pyrWeight_s = 0.345#1.3

def lognormal(m, s):
        mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
        std = np.sqrt(np.log((s/m)**2 + 1))
        #import pdb; pdb.set_trace()
        return max(np.random.lognormal(mean, std, 1), 0.00000001)

def set_pyr_w(m, s):
    global pyrWeight_m
    global pyrWeight_s
    pyrWeight_m = m
    pyrWeight_s = s

def AMPANMDA(syn_params, sec_x, sec_id):
    """Create a bg2pyr synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """

    lsyn = h.ProbAMPANMDA2(sec_x, sec=sec_id)

    if syn_params.get('tau_r_AMPA'):
        lsyn.tau_r_AMPA = float(syn_params['tau_r_AMPA'])
    if syn_params.get('tau_d_AMPA'):
        lsyn.tau_d_AMPA = float(syn_params['tau_d_AMPA'])
    if syn_params.get('tau_r_NMDA'):
        lsyn.tau_r_NMDA = float(syn_params['tau_r_NMDA'])
    if syn_params.get('tau_d_NMDA'):
        lsyn.tau_d_NMDA = float(syn_params['tau_d_NMDA'])
    if syn_params.get('Use'):
        lsyn.Use = float(syn_params['Use'])
    if syn_params.get('Dep'):
        lsyn.Dep = float(syn_params['Dep'])
    if syn_params.get('Fac'):
        lsyn.Fac = float(syn_params['Fac'])
    if syn_params.get('e'):
        lsyn.e = float(syn_params['e'])
    if syn_params.get('initW'):
        h.distance(sec=sec_id.cell().soma[0])
        dist = h.distance(sec_id(sec_x))
        fullsecname = sec_id.name()
        sec_type = fullsecname.split(".")[1][:4]
        sec_id = int(fullsecname.split("[")[-1].split("]")[0])

        # if pyrWeight_s == 0:
        #     base = float(pyrWeight_m)
        # else:
        #     base = float(np.clip(lognormal(pyrWeight_m, pyrWeight_s), 0, 5))

        ####OLD
        # dend = lambda x: 0.9278403931213186 * ( 1.0022024845737223 ** x )
        # close_apic = lambda x: 0.9131511669645764 * ( 1.0019436631560847 ** x )
        # far_apic = lambda x: 0.16857988107990907 * ( 1.0039628707324273 ** x )
        #############

        #distance based conductance scaling functions.
        #dend = lambda x: 0.9475625702815389 * ( 1.001318965242205 ** x )
        #close_apic = lambda x: 0.8522367331040966 * ( 1.0020433032052223 ** x )
        #far_apic = lambda x: 0.09043087364217033 * ( 1.004632615014859 ** x )
        
        dend = lambda x: ( 1.001 ** x )
        close_apic = lambda x: ( 1.002 ** x )
        #far_apic = lambda x: ( 1.002 ** x )
        far_apic = lambda x: 1

        if sec_type == "dend":
            base = float(np.clip(lognormal(pyrWeight_m, pyrWeight_s), 0, 5))
            lsyn.initW = base * dend(dist)
        elif sec_type == "apic":
            if dist < 750:
                base = float(np.clip(lognormal(pyrWeight_m, pyrWeight_s), 0, 5))
                lsyn.initW = base * close_apic(dist)
            else:
                base = float(np.clip(lognormal(0.17, 0.2), 0, 5))
                lsyn.initW = base * far_apic(dist)

        lsyn.initW = np.clip(float(lsyn.initW), 0, 5)
    if syn_params.get('u0'):
        lsyn.u0 = float(syn_params['u0'])
    return lsyn


def ampanmda(syn_params, xs, secs):
    """Create a list of bg2pyr synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = AMPANMDA(syn_params, x, sec)
        syns.append(syn)
    return syns

def Bg2Pyr(syn_params, sec_x, sec_id):
    """Create a bg2pyr synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """

    lsyn = h.bg2pyr(sec_x, sec=sec_id)

    if syn_params.get('initW'):
        lsyn.initW = float(syn_params['initW'])
    if syn_params.get('taun1'):
        lsyn.taun1 = float(syn_params['taun1'])
    if syn_params.get('taun2'):
        lsyn.taun2 = float(syn_params['taun2'])
    if syn_params.get('gNMDAmax'):
        lsyn.gNMDAmax = float(syn_params['gNMDAmax'])
    if syn_params.get('enmda'):
        lsyn.enmda = float(syn_params['enmda'])
    if syn_params.get('taua1'):
        lsyn.taua1 = float(syn_params['taua1'])
    if syn_params.get('taua2'):
        lsyn.taua2 = float(syn_params['taua2'])
    if syn_params.get('gAMPAmax'):
        lsyn.gAMPAmax = float(syn_params['gAMPAmax'])
    if syn_params.get('eampa'):
        lsyn.eampa = float(syn_params['eampa'])
    return lsyn


def bg2pyr(syn_params, xs, secs):
    """Create a list of bg2pyr synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = Pyr2Pyr(syn_params, x, sec)
        syns.append(syn)
    return syns

def Pyr2Int(syn_params, sec_x, sec_id):
    """Create a pyr2int synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """

    lsyn = h.pyr2int(sec_x, sec=sec_id)

    if syn_params.get('AlphaTmax_ampa'):
        lsyn.AlphaTmax_ampa = float(syn_params['AlphaTmax_ampa']) # par.x(21)
    if syn_params.get('Beta_ampa'):
        lsyn.Beta_ampa = float(syn_params['Beta_ampa']) # par.x(22)
    if syn_params.get('Cdur_ampa'):
        lsyn.Cdur_ampa = float(syn_params['Cdur_ampa']) # par.x(23)
    if syn_params.get('gbar_ampa'):
        lsyn.gbar_ampa = float(syn_params['gbar_ampa']) # par.x(24)
    if syn_params.get('Erev_ampa'):
        lsyn.Erev_ampa = float(syn_params['Erev_ampa']) # par.x(16)

    if syn_params.get('AlphaTmax_nmda'):
        lsyn.AlphaTmax_nmda = float(syn_params['AlphaTmax_nmda']) # par.x(25)
    if syn_params.get('Beta_nmda'):
        lsyn.Beta_nmda = float(syn_params['Beta_nmda']) # par.x(26)
    if syn_params.get('Cdur_nmda'):
        lsyn.Cdur_nmda = float(syn_params['Cdur_nmda']) # par.x(27)
    if syn_params.get('gbar_nmda'):
        lsyn.gbar_nmda = float(syn_params['gbar_nmda']) # par.x(28)
    if syn_params.get('Erev_nmda'):
        lsyn.Erev_nmda = float(syn_params['Erev_nmda']) # par.x(16)
    
    if syn_params.get('initW'):
        lsyn.initW = float(syn_params['initW']) * random.uniform(0.5,1.0) # par.x(0) * rC.uniform(0.5,1.0)//rand.normal(0.5,1.5) //`rand.repick() 

    if syn_params.get('Wmax'):
        lsyn.Wmax = float(syn_params['Wmax']) * lsyn.initW # par.x(1) * lsyn.initW
    if syn_params.get('Wmin'):
        lsyn.Wmin = float(syn_params['Wmin']) * lsyn.initW # par.x(2) * lsyn.initW
    #delay = float(syn_params['initW']) # par.x(3) + delayDistance
    #lcon = new NetCon(&v(0.5), lsyn, 0, delay, 1)

    if syn_params.get('lambda1'):
        lsyn.lambda1 = float(syn_params['lambda1']) # par.x(6)
    if syn_params.get('lambda2'):
        lsyn.lambda2 = float(syn_params['lambda2']) # par.x(7)
    if syn_params.get('threshold1'):
        lsyn.threshold1 = float(syn_params['threshold1']) # par.x(8)
    if syn_params.get('threshold2'):
        lsyn.threshold2 = float(syn_params['threshold2']) # par.x(9)
    if syn_params.get('tauD1'):
        lsyn.tauD1 = float(syn_params['tauD1']) # par.x(10)
    if syn_params.get('d1'):
        lsyn.d1 = float(syn_params['d1']) # par.x(11)
    if syn_params.get('tauD2'):
        lsyn.tauD2 = float(syn_params['tauD2']) # par.x(12)
    if syn_params.get('d2'):
        lsyn.d2 = float(syn_params['d2']) # par.x(13)
    if syn_params.get('tauF'):
        lsyn.tauF = float(syn_params['tauF']) # par.x(14)
    if syn_params.get('f'):
        lsyn.f = float(syn_params['f']) # par.x(15)

    if syn_params.get('bACH'):
        lsyn.bACH = float(syn_params['bACH']) # par.x(17)
    if syn_params.get('aDA'):
        lsyn.aDA = float(syn_params['aDA']) # par.x(18)
    if syn_params.get('bDA'):
        lsyn.bDA = float(syn_params['bDA']) # par.x(19)
    if syn_params.get('wACH'):
        lsyn.wACH = float(syn_params['wACH']) # par.x(20)
    
    return lsyn


def pyr2int(syn_params, xs, secs):
    """Create a list of pyr2int synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = Pyr2Int(syn_params, x, sec)
        syns.append(syn)
    return syns

def Int2Pyr(syn_params, sec_x, sec_id):
    """Create a int2pyr synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """

    lsyn = h.int2pyr(sec_x, sec=sec_id)

    h.distance(sec=sec_id.cell().soma[0])
    dist = h.distance(sec_id(sec_x))
    fullsecname = sec_id.name()
    sec_type = fullsecname.split(".")[1][:4]
    #sec_id = int(fullsecname.split("[")[-1].split("]")[0])

    #Assigns random generator of release probability.
    r = h.Random()
    r.MCellRan4()
    r.uniform(0,1)
    lsyn.setRandObjRef(r)

    generators.append(r)

    #Assigns release probabilty and conductance based on location of the synapse.
    if sec_type == "soma":
        lsyn.P_0 = 0.25#np.clip(np.random.normal(0.877, 0.052), 0, 1)
        lsyn.initW = 0.06#62.31
    if sec_type == "dend":
        if dist <= 50:
            lsyn.P_0 = 0.25#np.clip(np.random.normal(0.877, 0.052), 0, 1)
            lsyn.initW = 0.1#62.31
        else:
            lsyn.P_0 = 0.25#np.clip(np.random.normal(0.72, 0.1), 0, 1)
            lsyn.initW = 0.1#42.6#66.6
    if sec_type == "apic":
        lsyn.P_0 = 0.25#np.clip(np.random.normal(0.72, 0.1), 0, 1)
        lsyn.initW = 0.1#118.7#168.7

    #Short Term Plasticity
    #######################
    # SOM+
    #   d1: 0.96, tauD1: 40
    # PV+
    #   d1: 0.6, tauD1: 50
    #######################

    #if sec_type == "soma":
    #    #PV+
    #    lsyn.d1 = 0.6
    #    lsyn.tauD1 = 50
    #if sec_type == "dend":
    #    if dist <= 50:
    #        #PV+
    #        lsyn.d1 = 0.6
    #        lsyn.tauD1 = 50
    #    else:
    #        #SOM+
    #        lsyn.d1 = 0.96
    #        lsyn.tauD1 = 40
    #if sec_type == "apic":
    #    #SOM+
    #    lsyn.d1 = 0.96
    #    lsyn.tauD1 = 40

    if syn_params.get('AlphaTmax_ampa'):
        lsyn.AlphaTmax_ampa = float(syn_params['AlphaTmax_ampa']) # par.x(21)
    if syn_params.get('Beta_ampa'):
        lsyn.Beta_ampa = float(syn_params['Beta_ampa']) # par.x(22)
    if syn_params.get('Cdur_ampa'):
        lsyn.Cdur_ampa = float(syn_params['Cdur_ampa']) # par.x(23)
    if syn_params.get('gbar_ampa'):
        lsyn.gbar_ampa = float(syn_params['gbar_ampa']) # par.x(24)
    if syn_params.get('Erev_ampa'):
        lsyn.Erev_ampa = float(syn_params['Erev_ampa']) # par.x(16)

    if syn_params.get('AlphaTmax_nmda'):
        lsyn.AlphaTmax_nmda = float(syn_params['AlphaTmax_nmda']) # par.x(25)
    if syn_params.get('Beta_nmda'):
        lsyn.Beta_nmda = float(syn_params['Beta_nmda']) # par.x(26)
    if syn_params.get('Cdur_nmda'):
        lsyn.Cdur_nmda = float(syn_params['Cdur_nmda']) # par.x(27)
    if syn_params.get('gbar_nmda'):
        lsyn.gbar_nmda = float(syn_params['gbar_nmda']) # par.x(28)
    if syn_params.get('Erev_nmda'):
        lsyn.Erev_nmda = float(syn_params['Erev_nmda']) # par.x(16)
    
    # if syn_params.get('initW'):
    #     #lsyn.initW = float(syn_params['initW']) * random.uniform(0.5,1.0) # par.x(0) * rC.uniform(0.5,1.0)//rand.normal(0.5,1.5) //`rand.repick() 
    #     lsyn.initW = 3*float(max(np.random.normal(36, 18), 0.01))#2 * float(np.random.normal(12, np.sqrt(2)))#float(pyrWeight)

    if syn_params.get('Wmax'):
        lsyn.Wmax = float(syn_params['Wmax']) * lsyn.initW # par.x(1) * lsyn.initW
    if syn_params.get('Wmin'):
        lsyn.Wmin = float(syn_params['Wmin']) * lsyn.initW # par.x(2) * lsyn.initW
    #delay = float(syn_params['initW']) # par.x(3) + delayDistance
    #lcon = new NetCon(&v(0.5), lsyn, 0, delay, 1)

    if syn_params.get('lambda1'):
        lsyn.lambda1 = float(syn_params['lambda1']) # par.x(6)
    if syn_params.get('lambda2'):
        lsyn.lambda2 = float(syn_params['lambda2']) # par.x(7)
    if syn_params.get('threshold1'):
        lsyn.threshold1 = float(syn_params['threshold1']) # par.x(8)
    if syn_params.get('threshold2'):
        lsyn.threshold2 = float(syn_params['threshold2']) # par.x(9)
    if syn_params.get('tauD1'):
         lsyn.tauD1 = float(syn_params['tauD1']) # par.x(10)
    if syn_params.get('d1'):
         lsyn.d1 = float(syn_params['d1']) # par.x(11)
    if syn_params.get('tauD2'):
        lsyn.tauD2 = float(syn_params['tauD2']) # par.x(12)
    if syn_params.get('d2'):
        lsyn.d2 = float(syn_params['d2']) # par.x(13)
    if syn_params.get('tauF'):
        lsyn.tauF = float(syn_params['tauF']) # par.x(14)
    if syn_params.get('f'):
        lsyn.f = float(syn_params['f']) # par.x(15)

    
    return lsyn


def int2pyr(syn_params, xs, secs):
    """Create a list of int2pyr synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    syns = []
    for x, sec in zip(xs, secs):
        syn = Int2Pyr(syn_params, x, sec)
        syns.append(syn)
    return syns


def Pyr2Pyr(syn_params, sec_x, sec_id):
    """Create a pyr2pyr synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """

    lsyn = h.pyr2pyr(sec_x, sec=sec_id)

    #Assigns random generator of release probability.
    r = h.Random()
    r.MCellRan4()
    r.uniform(0,1)
    lsyn.setRandObjRef(r)

    #A list of random generators is kept so that they are not automatically garbaged.
    generators.append(r)

    lsyn.P_0 = 0.6#np.clip(np.random.normal(0.53, 0.22), 0, 1)#Release probability

    if syn_params.get('AlphaTmax_ampa'):
        lsyn.AlphaTmax_ampa = float(syn_params['AlphaTmax_ampa']) # par.x(21)
    if syn_params.get('Beta_ampa'):
        lsyn.Beta_ampa = float(syn_params['Beta_ampa']) # par.x(22)
    if syn_params.get('Cdur_ampa'):
        lsyn.Cdur_ampa = float(syn_params['Cdur_ampa']) # par.x(23)
    if syn_params.get('gbar_ampa'):
        lsyn.gbar_ampa = float(syn_params['gbar_ampa']) # par.x(24)
    if syn_params.get('Erev_ampa'):
        lsyn.Erev_ampa = float(syn_params['Erev_ampa']) # par.x(16)

    if syn_params.get('AlphaTmax_nmda'):
        lsyn.AlphaTmax_nmda = float(syn_params['AlphaTmax_nmda']) # par.x(25)
    if syn_params.get('Beta_nmda'):
        lsyn.Beta_nmda = float(syn_params['Beta_nmda']) # par.x(26)
    if syn_params.get('Cdur_nmda'):
        lsyn.Cdur_nmda = float(syn_params['Cdur_nmda']) # par.x(27)
    if syn_params.get('gbar_nmda'):
        lsyn.gbar_nmda = float(syn_params['gbar_nmda']) # par.x(28)
    if syn_params.get('Erev_nmda'):
        lsyn.Erev_nmda = float(syn_params['Erev_nmda']) # par.x(16)
    
    if syn_params.get('initW'):
        h.distance(sec=sec_id.cell().soma[0])
        dist = h.distance(sec_id(sec_x))
        fullsecname = sec_id.name()
        sec_type = fullsecname.split(".")[1][:4]
        sec_id = int(fullsecname.split("[")[-1].split("]")[0])

        # if pyrWeight_s == 0:
        #     base = float(pyrWeight_m)
        # else:
        #     base = float(np.clip(lognormal(pyrWeight_m, pyrWeight_s), 0, 5))

        ####OLD
        # dend = lambda x: 0.9278403931213186 * ( 1.0022024845737223 ** x )
        # close_apic = lambda x: 0.9131511669645764 * ( 1.0019436631560847 ** x )
        # far_apic = lambda x: 0.16857988107990907 * ( 1.0039628707324273 ** x )
        #############

        #distance based conductance scaling functions.
        #dend = lambda x: 0.9475625702815389 * ( 1.001318965242205 ** x )
        #close_apic = lambda x: 0.8522367331040966 * ( 1.0020433032052223 ** x )
        #far_apic = lambda x: 0.09043087364217033 * ( 1.004632615014859 ** x )
        
        dend = lambda x: ( 1.00 ** x )
        close_apic = lambda x: ( 1.00 ** x )
        #far_apic = lambda x: ( 1.002 ** x )
        far_apic = lambda x: 1

        if sec_type == "dend":
            base = float(np.clip(lognormal(pyrWeight_m, pyrWeight_s), 0, 5))
            lsyn.initW = base * dend(dist)
        elif sec_type == "apic":
            if dist < 750:
                base = float(np.clip(lognormal(pyrWeight_m, pyrWeight_s), 0, 5))
                lsyn.initW = base * close_apic(dist)
            else:
                base = float(np.clip(lognormal(pyrWeight_m, pyrWeight_s), 0, 5))
                lsyn.initW = base * far_apic(dist)

        lsyn.initW = np.clip(float(lsyn.initW), 0, 5)


    if syn_params.get('Wmax'):
        lsyn.Wmax = float(syn_params['Wmax']) * lsyn.initW # par.x(1) * lsyn.initW
    if syn_params.get('Wmin'):
        lsyn.Wmin = float(syn_params['Wmin']) * lsyn.initW # par.x(2) * lsyn.initW
    #delay = float(syn_params['initW']) # par.x(3) + delayDistance
    #lcon = new NetCon(&v(0.5), lsyn, 0, delay, 1)

    if syn_params.get('lambda1'):
        lsyn.lambda1 = float(syn_params['lambda1']) # par.x(6)
    if syn_params.get('lambda2'):
        lsyn.lambda2 = float(syn_params['lambda2']) # par.x(7)
    if syn_params.get('threshold1'):
        lsyn.threshold1 = float(syn_params['threshold1']) # par.x(8)
    if syn_params.get('threshold2'):
        lsyn.threshold2 = float(syn_params['threshold2']) # par.x(9)
    if syn_params.get('tauD1'):
        lsyn.tauD1 = float(syn_params['tauD1']) # par.x(10)
    if syn_params.get('d1'):
        lsyn.d1 = float(syn_params['d1']) # par.x(11)
    if syn_params.get('tauD2'):
        lsyn.tauD2 = float(syn_params['tauD2']) # par.x(12)
    if syn_params.get('d2'):
        lsyn.d2 = float(syn_params['d2']) # par.x(13)
    if syn_params.get('tauF'):
        lsyn.tauF = float(syn_params['tauF']) # par.x(14)
    if syn_params.get('f'):
        lsyn.f = float(syn_params['f']) # par.x(15)

    if syn_params.get('bACH'):
        lsyn.bACH = float(syn_params['bACH']) # par.x(17)
    if syn_params.get('aDA'):
        lsyn.aDA = float(syn_params['aDA']) # par.x(18)
    if syn_params.get('bDA'):
        lsyn.bDA = float(syn_params['bDA']) # par.x(19)
    if syn_params.get('wACH'):
        lsyn.wACH = float(syn_params['wACH']) # par.x(20)
    
    return lsyn


def pyr2pyr(syn_params, xs, secs):
    """Create a list of pyr2pyr synapses
    :param syn_params: parameters of a synapse
    :param xs: list of normalized distances along the section
    :param secs: target sections
    :return: list of NEURON synpase objects
    """
    np.random.seed(2129)
    syns = []
    for x, sec in zip(xs, secs):
        syn = Pyr2Pyr(syn_params, x, sec)
        syns.append(syn)
    return syns


def load():
    add_synapse_model(AMPANMDA, 'ampanmda', overwrite=False)
    add_synapse_model(AMPANMDA, overwrite=False)
    add_synapse_model(Bg2Pyr, 'bg2pyr', overwrite=False)
    add_synapse_model(Bg2Pyr, overwrite=False)
    add_synapse_model(Pyr2Pyr, 'pyr2pyr', overwrite=False)
    add_synapse_model(Pyr2Pyr, overwrite=False)
    add_synapse_model(Pyr2Int, 'pyr2int', overwrite=False)
    add_synapse_model(Pyr2Int, overwrite=False)
    add_synapse_model(Int2Pyr, 'int2pyr', overwrite=False)
    add_synapse_model(Int2Pyr, overwrite=False)
    return

def syn_params_dicts(syn_dir='../biophys_components/synaptic_models'):
    """
    returns: A dictionary of dictionaries containing all
    properties in the synapse json files
    """
    files = glob.glob(os.path.join(syn_dir,'*.json'))
    data = {}
    for fh in files:
        with open(fh) as f:
            data[os.path.basename(fh)] = json.load(f) #data["filename.json"] = {"prop1":"val1",...}
    return data
