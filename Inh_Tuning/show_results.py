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

def plot_IPSCs(file, time = 4000):
    df = pd.read_csv(file)
    ipscs = 1000 * np.array(df["IPSC"])


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
    print("Mean: ", np.mean(ipscs))
    print("Std: ", np.std(ipscs))
    # plt.figure()
    # sb.distplot(ipscs, fit=lognorm, bins = 100)
    # plt.xscale('log')

    plt.figure()
    sb.distplot(ipscs)
    plt.title(file)

    plt.show()

plot_IPSCs("dend_27.79_0.21_IPSCs.csv")
plot_IPSCs("som_IPSCs.csv")
plot_IPSCs("som_dend_IPSCs.csv")
plot_IPSCs("dend_42.76_2_IPSCs.csv")

plot_v("output/v_report.h5", show=False)
plt.title("Voltage")

plt.figure()
plt.title("Clamp Current")
plt.ylabel("nA")

data = load_dataset("output/se_clamp_report.h5", groups=1)
print("pA:", 1e3*(data[3999, 0] - min(data[4000:, 0])))

#FSI/Clos Inhibition: 208.3±58.7 pA/pair
#   (58.7)^2 = 3(s^2) ---> s = sqrt((58.7^2)/3) = 33.89
#   Avg 3 per pair ---> mean = 69.4 pA, std = 33.89 pA
#   Mean Weight: 36.84 Std: 17.99
#LTS/Far Inhibition: 26.5±1.6 pA/pair
#   Avg 3 per pair ---> mean = 8.8 pA std = sqrt((1.6^2)/3) = 0.92
#   Mean Weight: 4.67 Std: 0.49

plot_se("output/se_clamp_report.h5", show=True)

# def plot_EPSPs(file, time = 4000):
#     df = pd.read_csv(file)
#     epsps = np.array(df["EPSP"])


#     # conductances = np.array(df["Conductance"])
#     # print("Mean: ", np.mean(conductances))
#     # print("Std: ", np.std(conductances))
#     # plt.figure()
#     # sb.distplot(conductances)
#     # plt.xscale('log')

#     # plt.figure()
#     # sb.distplot(conductances)

#     #plt.show()
#     # plt.figure()
#     # plt.hist(epsps, bins = 400)

#     # plt.figure()
#     # plt.hist(epsps, bins = 400)
#     # plt.xscale("log")
#     print("Mean: ", np.mean(epsps))
#     print("Std: ", np.std(epsps))
#     plt.figure()
#     sb.distplot(epsps, fit=lognorm, bins = 100)
#     plt.xscale('log')

#     plt.figure()
#     sb.distplot(epsps)

#     #plt.show()

# # plot_scatter("scaled/scale_0.4_EPSPs.csv", label="scaled 0.4")
# # #plot_scatter("EPSP_files/0.4_EPSPs.csv", label="base 0.4")
# # plt.axhline(y = 0.4, color = 'green', linestyle = '-', label="target1") 
# # plt.xlabel("Distance from soma (um)")
# # plt.ylabel("Somatic EPSP magnitude")
# # plt.legend()

# # plt.show()

# plot_EPSPs("Tunings/0.45_0.35_EPSPs.csv")
# #plot_EPSPs("new_scaled_results/1.1_1.4_EPSPs.csv")
# plt.show()
# #plot_EPSPs("scaled_results/0.45_0.4_EPSPs.csv")

# plot_scatter("scaled/scale_0.8_EPSPs.csv", label="scaled 0.8")
# plot_scatter("all_0.8_EPSPs.csv", label="base 0.8")
# plt.axhline(y = 0.8, color = 'green', linestyle = '-', label="target")

# plt.show()

# plot_scatter("EPSP_files/scale_0.2_EPSPs.csv", label="scaled 0.2")
# plot_scatter("EPSP_files/0.2_EPSPs.csv", label="base 0.2")
# plt.axhline(y = 0.2, color = 'green', linestyle = '-', label="target3")

# plot_scatter("EPSP_files/scale_0.8_EPSPs.csv", label="scaled 0.8")
# plot_scatter("EPSP_files/0.8_EPSPs.csv", label="base 0.8")
# plt.axhline(y = 0.8, color = 'green', linestyle = '-', label="target2")

# plot_scatter("EPSP_files/scale_0.4_EPSPs.csv", label="scaled 0.4")
# plot_scatter("EPSP_files/0.4_EPSPs.csv", label="base 0.4")
# plt.axhline(y = 0.4, color = 'green', linestyle = '-', label="target1") 
# plt.xlabel("Distance from soma (um)")
# plt.ylabel("Somatic EPSP magnitude")
# plt.legend()
# # plot_scatter("EPSP_files/0.8_EPSPs.csv")
# #plot_scatter("EPSP_files/1250_0.4_EPSPs.csv")
# plt.show()