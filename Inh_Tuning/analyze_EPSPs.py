import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import my_plotting
import numpy as np
import matplotlib.pyplot as plt
#import scipy.signal as s
import seaborn as sb
import pandas as pd

#Takes a voltage trace and returns the magnitude of an EPSP
#resulting from an input at time time.
def calc_EPSP(v, time = 4000):
    base = v[time - 2]
    peak = min(v[time - 2:])
    return peak - base

#Takes a file with some number of voltage traces, each with one EPSP.
#Returns an array of EPSP magnitudes.
def calc_EPSPs(file, time = 4000):
    data = np.array(my_plotting.load_dataset(file)['data'])
    epsps = np.zeros(data.shape[1])

    for i in range(len(epsps)):
        epsps[i] = calc_EPSP(data[:, i], time=time)

    return epsps

def save_EPSPs(file, output, time = 4000):
    epsps = calc_EPSPs(file, time)

    df = pd.read_csv(output)
    df["EPSP"] = epsps
    df.to_csv(output, index=False)

def plot_EPSPs(file, time = 4000):
    epsps = calc_EPSPs(file, time)

    # plt.figure()
    # plt.hist(epsps, bins = 400)

    # plt.figure()
    # plt.hist(epsps, bins = 400)
    # plt.xscale("log")
    print("Mean: ", np.mean(epsps))
    print("Std: ", np.std(epsps))
    plt.figure()
    sb.distplot(epsps)
    plt.xscale('log')

    plt.figure()
    sb.distplot(epsps)

    plt.show()

def plot_scatter(file, label=None):
    df = pd.read_csv(file)
    plt.scatter(df["Distance"], df["EPSP"], label=label)
    #plt.show()

# plot_scatter("EPSPs.csv")
# plot_scatter("EPSP_files/0.4_EPSPs.csv")
# plt.show()
#save_EPSPs("output/v_report.h5", "synapse_info.csv")
#plot_EPSPs("output/v_report.h5")
#plot_EPSPs("results/500-0.65-0.48v_report.h5")
#import pdb; pdb.set_trace()