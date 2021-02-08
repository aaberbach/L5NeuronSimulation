import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b):
    return a * (b ** x)
    
def func2(x, b):
    return b**x

def lin_func(x, m, b):
    return m*x + b

def get_opts(dists, epsps, func, plot=False, title=None, syn_ids=None):
    dists = dists - 50
    epsps = 0.4 / epsps
    popt, pcov = curve_fit(func, dists, epsps)

    if plot:
        plt.figure()
        plt.scatter(dists, func(dists, *popt))

        plt.title(title)
        print(popt)

        plt.scatter(dists, epsps)
        if type(syn_ids) != type(None):
            for i, txt in enumerate(syn_ids):
                plt.annotate(str(int(txt)), (dists[i], epsps[i]))
        #plt.show()

    return popt

def comb(ids1, ids2):
    #return np.unique(np.concatenate((ids1,ids2),0))
    return np.intersect1d(ids1, ids2)

df = pd.read_csv("info_EPSPs.csv")
names = np.array(df["Name"])
syn_ids = np.zeros(len(names))
for i in range(len(names)):
    syn_ids[i] = int(names[i].split("[")[-1].split("]")[0])

#import pdb; pdb.set_trace()

df = pd.read_csv("EPSP_files/500_0.4_EPSPs.csv")
epsps = np.array(df["EPSP"])
distances = np.array(df["Distance"])
types = np.array(df["Type"])


#include = np.where(distances < 500)[0]
dends = np.where(types == "dend")[0]
apics = np.where(types == "apic")[0]
fars = np.where(distances >= 500)[0]
closes = np.where(distances < 500)[0]
#extrafars = np.where(distances >= 600)[0]
#np.unique(np.concatenate((a,b),0))
#[0.94330582 1.00199629] < 600
#[0.08563904 1.00487715] >= 600

high_ids = np.where(syn_ids >= 60)[0]
low_ids = np.where(syn_ids < 60)[0]

get_opts(distances, epsps, func, plot=True, title="All")
plt.show()

# [0.00159513 1.01897285] Dend Lin

# [1.02451963 1.00138249] Dend Exp
# [0.59768734 1.00326839] Far Apic Exp (High IDs)
# [0.62153507 1.00248601] Far Apic Exp (Low IDs)
# [0.97717162 1.00195518] Close Apic Exp


get_opts(distances[dends], epsps[dends], func, plot=True, title="Dend")
get_opts(distances[comb(comb(apics,high_ids),fars)], epsps[comb(comb(apics,high_ids),fars)], func, plot=True ,title="Far Apic (High IDs)")
get_opts(distances[comb(comb(apics,low_ids),fars)], epsps[comb(comb(apics,low_ids),fars)], func, plot=True ,title="Far Apic (Low IDs)")
get_opts(distances[comb(apics,closes)], epsps[comb(apics,closes)], func, plot=True, title="Close Apic")

plt.show()
