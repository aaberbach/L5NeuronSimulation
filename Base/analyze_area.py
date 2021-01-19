
import numpy as np

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
    import pdb; pdb.set_trace()
    print("Dend Excitatory:", np.trunc(sum(dend_lens) * 1.4))
    print("Apic Excitatory:", np.trunc(sum(apic_lens) * 1.4))
    print("Soma Inhibitory:", np.trunc(148 + sum(dend_lens) * 0.14 + sum(apic_lens) * 0.14))
    #import pdb; pdb.set_trace()