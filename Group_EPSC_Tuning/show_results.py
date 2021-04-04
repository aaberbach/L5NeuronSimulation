from analyze_EPSCs import *
from scipy.stats import lognorm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb

def plot_EPSCs(file, time = 4000):
    df = pd.read_csv(file)
    epscs = np.array(df["EPSC"])
    #import pdb; pdb.set_trace()

    # conductances = np.array(df["Conductance"])
    # print("Mean: ", np.mean(conductances))
    # print("Std: ", np.std(conductances))
    # plt.figure()
    # sb.distplot(conductances)
    # plt.xscale('log')

    # plt.figure()
    # sb.distplot(conductances)

    #plt.show()
    # plt.figure()
    # plt.hist(epscs, bins = 400)

    # plt.figure()
    # plt.hist(epscs, bins = 400)
    # plt.xscale("log")
    print("Mean: ", np.mean(epscs) * 1000)
    print("Std: ", np.std(epscs) * 1000)


    plt.figure()
    sb.distplot(epscs)

    plt.show()

def plot_traces(file):
    traces = get_traces(file)
    for i in range(traces.shape[1]):
        plt.plot(traces[:, i])
        plt.show()
    

#plot_traces("output/se_clamp_report.h5")
plot_EPSCs("EPSCs.csv")
#plot_EPSCs("Tunings/0.1_0.65_EPSCs.csv")