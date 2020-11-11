from bmtk.analyzer.cell_vars import plot_report

_ = plot_report(config_file='FRPrediction/simulation_config.json', node_ids=[0, 1], report_name='v_report')
#_ = plot_traces(config_file='model_info/simulation_config.json', node_ids=[0], report_name='cai_report')

from bmtk.analyzer import spike_trains

spike_trains.plot_raster(config_file='FRPrediction/simulation_config.json', spikes_file = 'FRPrediction/exc_stim_spikes.h5')