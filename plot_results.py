import numpy as np
import matplotlib.pyplot as plt
import h5py
from bmtk.analyzer.cell_vars import _get_cell_report, plot_report
import matplotlib.pyplot as plt
import pandas as pd 
from scipy.signal import find_peaks
import pickle

f = "inp_to_out.pkl"
file = open(f, 'rb')
res = pickle.load(file)
file.close()

print(res)

plt.figure()

keys = list(res.keys())
excs = [float(frs.split(',')[0]) for frs in keys]
inhs = [float(frs.split(',')[1]) for frs in keys]

both_x = np.array([excs[i] for i in range(len(keys)) if inhs[i] != 2 and excs[i] <= 5])
both_y = np.array([res[keys[i]] for i in range(len(keys)) if inhs[i] != 2 and excs[i] <= 5])

m, b = np.polyfit(both_x, both_y, 1)

#plt.plot(both_x, m*both_x + b, color="orange")
plt.scatter(both_x, both_y, label='exc + inh', marker='.', color='orange')

exc_x = np.array([excs[i] for i in range(len(keys)) if inhs[i] == 2 and excs[i] <= 5])
exc_y = np.array([res[keys[i]] for i in range(len(keys)) if inhs[i] == 2 and excs[i] <= 5])

m, b = np.polyfit(exc_x, exc_y, 1)

#plt.plot(exc_x, m*exc_x + b, color="blue")

#import pdb; pdb.set_trace()
#plt.figure()
plt.scatter(exc_x, exc_y, label='exc only', marker='.', color='blue')
#plt.scatter(both_x, both_y, label='exc + inh', marker='.')
plt.legend()

plt.show()