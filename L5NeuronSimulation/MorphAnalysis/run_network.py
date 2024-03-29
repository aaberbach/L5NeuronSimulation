from bmtk.simulator import bionet
import numpy as np
from neuron import h

pc = h.ParallelContext()  # object to access MPI methods
MPI_size = int(pc.nhost())
MPI_rank = int(pc.id())

config_file = 'simulation_config.json'

conf = bionet.Config.from_json(config_file, validate=True)
conf.build_env()

graph = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=graph)

from analyze_area import analyze_area, make_seg_df
make_seg_df(list(graph.get_local_cells().values())[0])
#analyze_area(list(graph.get_local_cells().values())[0]._morph.seg_prop)

#sim.run()
pc.barrier()
