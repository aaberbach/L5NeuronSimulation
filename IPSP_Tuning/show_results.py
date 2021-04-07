# from analyze_EPSPs import *
# from scipy.stats import lognorm
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb

from my_plotting import *

def plot_IPSPs(file, time = 4000):
    df = pd.read_csv(file)
    ipsps = np.array(df["IPSP"])


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
    # plt.hist(ipscs, bins = 400)

    # plt.figure()
    # plt.hist(ipscs, bins = 400)
    # plt.xscale("log")
    print("Mean: ", np.mean(ipsps))
    print("Std: ", np.std(ipsps))
    # plt.figure()
    # sb.distplot(ipscs, fit=lognorm, bins = 100)
    # plt.xscale('log')

    plt.figure()
    sb.distplot(ipsps)
    plt.title(file)

    plt.show()

plot_IPSPs("IPSPs.csv")

plot_all_v("output/v_report.h5")
plt.title("Voltage")

plt.show()