import pandas as pd
import numpy as np
from neuron import h

def make_seg_df(cell):
    df = pd.DataFrame()
    i = 0
    sec_ids = []
    full_names = []
    xs = []
    parts = []
    distances = []
    h.distance(sec=cell.hobj.soma[0])
    for sec in cell.hobj.all:
        for seg in sec:
            distances.append(h.distance(seg))
            sec_ids.append(i)
            xs.append(seg.x)
            fullsecname = sec.name()
            sec_type = fullsecname.split(".")[1][:4]
            parts.append(sec_type)
            full_names.append(str(seg))
        i += 1

    df["Sec ID"] = sec_ids
    df["X"] = xs
    df["Type"] = parts
    df["Distance"] = distances
    df["Name"] = full_names

    df.to_csv("Segments.csv", index=False)
    import pdb; pdb.set_trace()

def analyze_area(prop):
    #import pdb; pdb.set_trace()

    types = prop['type']
    areas = prop['area']
    lens = prop['length']

    soma_ids = np.where(types == 1)[0]
    dend_ids = np.where(types == 3)[0]
    apic_ids = np.where(types == 4)[0]

    soma_areas = areas[soma_ids]
    dend_areas = areas[dend_ids]
    apic_areas = areas[apic_ids]

    soma_lens = lens[soma_ids]
    dend_lens = lens[dend_ids]
    apic_lens = lens[apic_ids]
    
    print("Dend Excitatory:", np.trunc(sum(dend_lens) * 1.4))
    print("Dend Inhibitory:", np.trunc(sum(dend_lens) * 0.14))
    print("Apic Excitatory:", np.trunc(sum(apic_lens) * 1.4))
    print("Apic Inhibitory:", np.trunc(sum(apic_lens) * 0.14))
    print("Soma Inhibitory:", 148)