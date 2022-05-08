import pandas as pd
import numpy as np
from neuron import h

def save_degrees(cell):
    degrees = {}
    calculate_degree(h.SectionRef(sec=cell.hobj.soma[0]), degrees, 0)

    df_dict = {}
    df_dict["SectionName"] = list(degrees.keys())
    df_dict["Degrees"] = list(degrees.values())
    
    df = pd.DataFrame(df_dict)
    df.to_csv("SectionDegrees.csv", index=False)

def calculate_degree(sref, degrees, deg):
    degrees[sref.sec.name()] = deg

    for c in sref.child:
        calculate_degree(h.SectionRef(sec=c), degrees, deg+1)

def make_seg_df(cell):
    seg_locs = cell.morphology.seg_coords['p05']
    px = seg_locs[0]
    py = seg_locs[1]
    pz = seg_locs[2]
    df = pd.DataFrame()
    i = 0
    j = 0
    lens = []
    diams = []
    bmtk_ids = []
    sec_ids = []
    full_names = []
    xs = []
    parts = []
    distances = []
    elec_distances = []
    h.distance(sec=cell.hobj.soma[0])
    zz = h.Impedance()
    zz.loc(cell.hobj.soma[0](0.5))
    zz.compute(25,1)
    
    for sec in cell.hobj.all:
        for seg in sec:
            lens.append(seg.sec.L)
            diams.append(seg.sec.diam)
            distances.append(h.distance(seg))
            bmtk_ids.append(i)
            xs.append(seg.x)
            fullsecname = sec.name()
            sec_ids.append(int(fullsecname.split("[")[2].split("]")[0]))
            sec_type = fullsecname.split(".")[1][:4]
            parts.append(sec_type)
            full_names.append(str(seg))
            elec_distances.append(zz.ratio(seg))
            j += 1
        i += 1

    df["BMTK ID"] = bmtk_ids
    df["X"] = xs
    df["Type"] = parts
    df["Sec ID"] = sec_ids
    df["Distance"] = distances
    df["Section_L"] = lens
    df["Section_diam"] = diams
    df["Coord X"] = px
    df["Coord Y"] = py
    df["Coord Z"] = pz
    df["Elec_distance"] = elec_distances


    df.to_csv("Segments.csv", index=False)
    #import pdb; pdb.set_trace()

def analyze_area(prop):
    #import pdb; pdb.set_trace()

    types = prop['type']
    # areas = prop['area']
    lens = prop['length']
    dists = prop['dist']

    soma_ids = np.where(types == 1)[0]
    close_dend_ids = np.where((types == 3) & (dists < 50))[0]
    far_dend_ids = np.where((types == 3) & (dists >= 50))[0]
    apic_ids = np.where(types == 4)[0]

    print("LENGTHS:")
    print("Close Dend:", np.trunc(sum(lens[close_dend_ids])))
    print("Further Dend:", np.trunc(sum(lens[far_dend_ids])))
    print("Apic:", np.trunc(sum(lens[apic_ids])))
    print()
    # soma_areas = areas[soma_ids]
    # dend_areas = areas[dend_ids]
    # apic_areas = areas[apic_ids]

    # soma_lens = lens[soma_ids]
    # dend_lens = lens[dend_ids]
    # apic_lens = lens[apic_ids]
    
    # print("Dend Excitatory:", np.trunc(sum(dend_lens) * 1.4))
    # print("Dend Inhibitory:", np.trunc(sum(dend_lens) * 0.14))
    # print("Apic Excitatory:", np.trunc(sum(apic_lens) * 1.4))
    # print("Apic Inhibitory:", np.trunc(sum(apic_lens) * 0.14))
    # print("Soma Inhibitory:", 148)

    print("Dend Excitatory:", np.trunc(sum(lens[far_dend_ids]) * 1.4))
    print("Beta Dend Inhibitory:", np.trunc(sum(lens[far_dend_ids]) * 0.14))
    print("Gamma Dend Inhibitory:", np.trunc(sum(lens[close_dend_ids]) * 0.14))
    print("Apic Excitatory:", np.trunc(sum(lens[apic_ids]) * 1.4))
    print("Apic Inhibitory:", np.trunc(sum(lens[apic_ids]) * 0.14))
    print("Soma Inhibitory:", 148)
