import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b):
    return a * (b ** x)
    
# def func2(x, b):
#     return b**x

# def lin_func(x, m, b):
#     return m*x + b

def get_opts(dists, epsps, conductances, func, plot=False, title=None, syn_ids=None):
    dists = dists
    #import pdb; pdb.set_trace()
    epsps = conductances / epsps
    popt, pcov = curve_fit(func, dists, epsps)

    if plot:
        plt.figure()
        plt.scatter(dists, func(dists, *popt))

        plt.title(title)
        print(popt[0], "* (", popt[1], "^ x )")

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

#df = pd.read_csv("EPSP_files/500_0.4_EPSPs.csv")
#df = pd.read_csv("uniform_all_EPSPs.csv")
df = pd.read_csv("all_0.5_EPSPs.csv")

names = np.array(df["Name"])
syn_ids = np.zeros(len(names))
for i in range(len(names)):
    syn_ids[i] = int(names[i].split("[")[-1].split("]")[0])

epsps = np.array(df["EPSP"])
distances = np.array(df["Distance"])
types = np.array(df["Type"])
conductances = np.array(df["Conductance"])

# non_zeros = np.where(epsps != 0)[0]
# epsps=epsps[non_zeros]
# distances=distances[non_zeros]
# types=types[non_zeros]
# conductances=conductances[non_zeros]


#include = np.where(distances < 500)[0]
dends = np.where(types == "dend")[0]
apics = np.where(types == "apic")[0]
fars = np.where(distances >= 500)[0]
closes = np.where(distances < 500)[0]

high_ids = np.where(syn_ids >= 60)[0]
low_ids = np.where(syn_ids < 60)[0]

#Same function for all: 
#0.25004893047877735 * ( 1.0036760619164031 ^ x )
get_opts(distances, epsps, conductances, func, plot=True, title="All")
plt.show()

#Split only by distance
#0.9109627181995834 * ( 1.0018838469261027 ^ x )
get_opts(distances[closes], epsps[closes], conductances[closes], func, plot=True, title="Close")
#0.15197845065030272 * ( 1.0040953610902503 ^ x )
get_opts(distances[fars], epsps[fars], conductances[fars], func, plot=True, title="Far")
plt.show()

0.9371285272684119 * ( 1.0017027151340998 ^ x )
0.8842233502464731 * ( 1.001970475980598 ^ x )
0.15197845065030272 * ( 1.0040953610902503 ^ x )

#Split by distance and apic/dend
#0.9371285272684119 * ( 1.0017027151340998 ^ x )
get_opts(distances[dends], epsps[dends], conductances[dends], func, plot=True, title="Close Dend")
#0.8842233502464731 * ( 1.001970475980598 ^ x )
get_opts(distances[comb(closes,apics)], epsps[comb(closes,apics)], conductances[comb(closes,apics)], func, plot=True, title="Close Apic")
#0.15197845065030272 * ( 1.0040953610902503 ^ x )
get_opts(distances[fars], epsps[fars], conductances[fars], func, plot=True, title="Far")
plt.show()


# 0.9371285272684119 * ( 1.0017027151340998 ^ x )
get_opts(distances[dends], epsps[dends], conductances[dends], func, plot=True, title="Dend")
# 0.47555389722370983 * ( 1.0032837561151722 ^ x )
get_opts(distances[comb(comb(apics,high_ids),fars)], epsps[comb(comb(apics,high_ids),fars)], conductances[comb(comb(apics,high_ids),fars)], func, plot=True ,title="Far Apic (High IDs)")
# 0.5309901629694934 * ( 1.0025102142106757 ^ x )
get_opts(distances[comb(comb(apics,low_ids),fars)], epsps[comb(comb(apics,low_ids),fars)], conductances[comb(comb(apics,low_ids),fars)], func, plot=True ,title="Far Apic (Low IDs)")
# 0.8842233502464731 * ( 1.001970475980598 ^ x )
get_opts(distances[comb(apics,closes)], epsps[comb(apics,closes)], conductances[comb(apics,closes)], func, plot=True, title="Close Apic")

plt.show()
