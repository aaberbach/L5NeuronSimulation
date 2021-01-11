from bmtk.analyzer.cell_vars import plot_report
from bmtk.analyzer import spike_trains
import h5py
import matplotlib.pyplot as plt
import numpy as np

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import my_plotting
my_plotting.plot_v("gamma_results/v_report.h5", show=True)
#_ = plot_report(config_file='simulation_config.json', report_name='gamma_results/v_report.h5')
my_plotting.plot_se("gamma_results/clamped/se_clamp_report.h5", show=True)
my_plotting.plot_v("gamma_results/clamped/v_report.h5", show=True)
my_plotting.generate_prob_raster('inh_stim_spikes.h5', 'results/spikes.h5', 120000)
#parts = my_plotting.generate_spike_probs('inh_stim_spikes.h5', 'results/spikes.h5', 120000)
import pdb; pdb.set_trace()
my_plotting.plot_spike_gamma('inh_stim_spikes.h5', 120000)
#plt.show()

my_plotting.plot_spikes('inh_stim_spikes.h5')
my_plotting.plot_v('output/v_report.h5')
plt.show()
# #import pdb; pdb.set_trace()
_ = plot_report(config_file='simulation_config.json', report_name='v_report')


spike_trains.plot_raster(config_file='simulation_config.json', spikes_file = 'inh_stim_spikes.h5')
#f = h5py.File()
#_ = plot_report(config_file='simulation_config.json', node_ids=[0], report_name='cai_report')
#_ = plot_traces(config_file='model_info/simulation_config.json', node_ids=[0], report_name='cai_report')
raster_file = './output/spikes.h5'
#try:
#import pdb; pdb.set_trace()
f = h5py.File(raster_file,'r')
#import pdb; pdb.set_trace()
#import pdb; pdb.set_trace()
key = list(f['spikes'].keys())[0]
gids = np.array(f['spikes'][key]['node_ids'])
timestamps = np.array(f['spikes'][key]['timestamps'])
#import pdb; pdb.set_trace()
f.close()

spike_counts = []
for i in range(np.max(gids) + 1):
    spike_counts.append(len(np.where(gids == i)[0]))

#import pdb; pdb.set_trace()
plt.figure()
plt.hist(np.array(spike_counts)/5)

plt.figure()
plt.plot(timestamps,gids,'.')
plt.show()
print("FR:", np.mean(spike_counts) / 5)
print(np.array(spike_counts) / 5)
#except:
    #print("No spikes.")

#from bmtk.analyzer import spike_trains

#spike_trains.plot_raster(config_file='simulation_config.json', spikes_file = 'exc_stim_spikes.h5')
#spike_trains.plot_raster(config_file='simulation_config.json', spikes_file = 'inh_stim_spikes.h5')