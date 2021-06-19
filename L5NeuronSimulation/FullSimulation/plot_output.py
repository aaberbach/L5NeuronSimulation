import pandas as pd
import h5py
import numpy as np
import matplotlib.pyplot as plt

f = h5py.File('./output/v_report.h5')

plt.plot(f['report']['biophysical']['data'][:])
plt.show()


