FullSimulation
==============

This is the folder where the core simulation is built and run.

.. toctree::
  :maxdepth: 1

  build_network


Files
-----

Connections.csv
^^^^^^^^^^^^^^^
Contains information for each synapse in the simulation. Created during the run step. 

Columns:

* Node ID: the integer id of the presynaptic node within its population (exc, prox_inh, dist_inh)
* Distance: float representing the distance (um) from the synapse to the soma
* Conductance: the weight of the synapse (already scaled by distance)
* Type: what part of the cell the synapse is on (soma, dend, apic)
* Name: the full string that NEURON associates with the postsynaptic segment.
* Source Population: the node population that the presynaptic node is a member of (exc, prox_inh, dist_inh)
* Release Probability: release probability for the synapse


FunctionalGroups.csv 
^^^^^^^^^^^^^^^^^^^^
Contains the functional group id for each excitatory presynaptic node. Created during the build step by SimulationBuilder.save_groups


NetParams.json
^^^^^^^^^^^^^^
Contains parameters for the simulation that SimulationBuilder uses.

Structure:

* "lengths": total lengths of different types of segments (um)
  
  * "basal_dist": basal dendrites more than 50 um from the soma,
  * "basal_prox": basal dendrites less than 50 um from the soma,
  * "apic": apical dendrites
  
* "syn_density": number of synapses per um
  
  * "exc": excitatory synapses,
  * "inh": inhibitory synapses
  
* "n_soma_syns": number of (PV+) synapses on the soma,
* "divergence": distributions of connections per cell pairing
  
  * "exc": uniform with \{"min", "max"\} counts,
  * "peri_inh": normal with \{"m": mean, "s": standard deviation, "min":floor, "max":cap\},
  * "basal_inh": normal with \{"m": mean, "s": standard deviation, "min":floor, "max":cap\},
  * "apic_inh": normal with \{"m": mean, "s": standard deviation, "min":floor, "max":cap\},

* "groups": properties of functional groups
  
  * "cells_per_group": number of cells per functional group,
  * "cluster_radius": radius of the sphere that clusters are constrained to,
  * "group_radius": radius of the sphere that groups are constrained to
  
* "inh_frs": inhibitory firing rates
  
  * "proximal": PV+, normal with \{"m":mean, "s":standard deviation\},
  * "distal": SOM+, normal with \{"m":mean, "s":standard deviation\}
  
* "time": \{"start":when input should start (ms),"stop":simulation run time (ms)\},
* "dL": target length of each segment,
* "dt": time (ms) between each simulation step,
* "inh_shift": how many ms the average excitation trace is shifted to make the inhibition noise trace


Segments.csv
^^^^^^^^^^^^
Contains information about every segment in the morphology. Each segment is approximately 1 um in length.

Columns:

* BMTK ID: the ID that bmtk associates with the segments's section
* X: standardized (0 to 1) distance along the segment's section
* Type: whether the segment is soma, apic, dend (basal), or axon
* Sec ID: the ID of the segments's section within the morphology
* Distance: length (um) of closest path to soma
* Coord X: x coordinate (um) of the segment
* Coord Y: y coordinate (um) of the segment
* Coord Z: z coordinate (um) of the segment


:doc:`build_network.py <build_network>`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Contains the class SimulationBuilder which performs the build step of the simulation. This uses BMTK to create files that are then loaded in for the simulation step.

run_network.py
^^^^^^^^^^^^^^

.. automodule:: L5NeuronSimulation.FullSimulation.run_network
   :members:
