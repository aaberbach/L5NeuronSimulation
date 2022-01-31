import pandas as pd
import h5py
import numpy as np
import matplotlib.pyplot as plt

f = h5py.File('./output_control_test/v_report.h5','r')
g = h5py.File('./output_nablock_test/v_report.h5','r')
h = h5py.File('./output_control_test/NaTa_t.gNaTa_t_report.h5','r')


plt.figure()
plt.plot(f['report']['biophysical']['data'][:,0],label='control')
plt.plot(g['report']['biophysical']['data'][:,0],label='nablock')

plt.twinx()
plt.plot(h['report']['biophysical']['data'][:,37],color='r')
plt.legend()
plt.show()


