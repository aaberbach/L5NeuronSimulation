import pandas as pd
import h5py
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv('Connections.csv')

f = h5py.File('./output/v_report.h5','r')
na1 = h5py.File('./output/NaTa_t.gNaTa_t_report.h5','r')
na2 = h5py.File('./output_test2/NaTa_t.gNaTa_t_report.h5','r')


plt.figure()
plt.plot(na1['report']['biophysical']['data'][:,1700],color='b')
plt.plot(na2['report']['biophysical']['data'][:,1700],color='b',alpha=0.2)
plt.plot(na1['report']['biophysical']['data'][:,100],color='r')
plt.plot(na2['report']['biophysical']['data'][:,100],color='r',alpha=0.2)


plt.figure()
plt.plot(f['report']['biophysical']['data'][:,0],label='soma')
plt.plot(f['report']['biophysical']['data'][:,1700],label='apical')
plt.plot(f['report']['biophysical']['data'][:,100],label='basal')
plt.legend()

plt.show()


