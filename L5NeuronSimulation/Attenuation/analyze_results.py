import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def attenuation_plot(file, parts = ["dend", "apic"]):
    df = pd.read_csv(file)
    types = df["type"]
    ids = np.where((types == "dend") | (types=="apic"))[0]
    plt.scatter(df.iloc[ids]["distance"], df.iloc[ids]["attenuation"])
    plt.xlabel("Distance from the soma (um)")
    plt.ylabel("delta_V_soma / delta_V_dendrite")
    plt.title("Attenuation")
    plt.show()

attenuation_plot("results.csv")