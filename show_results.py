from bmtk.analyzer.compartment import plot_traces

_ = plot_traces(config_file='simulation_config.json', node_ids=[5072], report_name='v_report')
#_ = plot_traces(config_file='model_info/simulation_config.json', node_ids=[0], report_name='cai_report')