import pandas as pd
import numpy as np
from neuron import h

def make_seg_df(cell):
    #import pdb; pdb.set_trace()
    seg_locs = cell.morphology.seg_coords['p05']
    px = seg_locs[0]
    py = seg_locs[1]
    pz = seg_locs[2]
    df = pd.DataFrame()
    i = 0
    j = 0
    bmtk_ids = []
    sec_ids = []
    full_names = []
    xs = []
    parts = []
    distances = []
    h.distance(sec=cell.hobj.soma[0])
    for sec in cell.hobj.all:
        for seg in sec:
            distances.append(h.distance(seg))
            bmtk_ids.append(i)
            xs.append(seg.x)
            fullsecname = sec.name()
            sec_ids.append(int(fullsecname.split("[")[2].split("]")[0]))
            sec_type = fullsecname.split(".")[1][:4]
            parts.append(sec_type)
            full_names.append(str(seg))
            j += 1
        i += 1

    df["BMTK ID"] = bmtk_ids
    df["X"] = xs
    df["Type"] = parts
    df["Sec ID"] = sec_ids
    df["Distance"] = distances
    df["Coord X"] = px
    df["Coord Y"] = py
    df["Coord Z"] = pz

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