from math import pi, sin, log, exp
import matplotlib.pyplot as plt
import h5py
import numpy as np

def sweep(f_start, f_end, interval, n_steps):
    b = log(f_end/f_start) / interval
    a = 2 * pi * f_start / b
    t_list = []; g_t = [];
    for i in range(n_steps):
        delta = i / float(n_steps)
        t = interval * delta
        t_list.append(t)
        g_t.append(sin(a * exp(b * t)))
    return t_list, g_t

dt = .1 # samples per millisecond

t, g_t = sweep(f_start = 0.2, f_end = 200, interval = 50, n_steps = int(50000/dt))
g_t = np.array(g_t)
g_t = g_t + 0.09
g_t = np.concatenate((np.zeros((10000,)),g_t))

# current injection
g_t = -1*np.ones((int(1000/dt),))
g_t = np.concatenate((np.zeros((10000,)),g_t,np.zeros((10000,))))


d_t = [dt for i in range(g_t.shape[0])]

print(g_t.shape)
hf = h5py.File("zap.h5", "w")
hf.create_dataset("amplitudes", data=[g_t])
hf.create_dataset("dts", data=d_t)
hf.close()

plt.figure()
plt.plot(g_t)
plt.show()
