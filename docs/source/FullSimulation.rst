FullSimulation
==============

This is the folder where the core simulation is built and run.

.. toctree::
  :maxdepth: 1
  :hidden:

  build_network


Files
-----

bmtk_modifications.py
^^^^^^^^^^^^^^^^^^^^^
Contains run time modications to BMTK. 
Just import and run modify_bmtk().

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

* "cells": bmtk parameters for building the cell

  * "dynamic_params": a JSON file of conductances
  * "model_template": the hoc file used to build the cell
  * "morphology": the file containging the cell morphol0gy
  * "model_processing": a function BMTK will use to modify the cell
  * "segments_file": csv file containing all of the cell's segments

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
* "file_current_clamp": name of file to use as current clamp

* "record_cellvars": cellvars to record and on what part of the cell

  * "vars": variable name
  * "locs": location on cell ("all", "soma", "axon", "apic", "dend")

* "inh_frs": inhibitory firing rates
  
  * "proximal": PV+, normal with \{"m":mean, "s":standard deviation\},
  * "distal": SOM+, normal with \{"m":mean, "s":standard deviation\}
  
* "time": \{"start":when input should start (ms),"stop":simulation run time (ms)\},
* "dL": target length of each segment,
* "dt": time (ms) between each simulation step,
* "inh_shift": how many ms the average excitation trace is shifted to make the inhibition noise trace

L2-3NetParams.csv
^^^^^^^^^^^^^^^^^
Same structure as NetParams.csv above. This is the parameters file for an L2-3 cell used.
The cell is from the Allen Cell Database: https://celltypes.brain-map.org/experiment/electrophysiology/477127614.

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

batch_expanse_build.sh/batch_expanse_run.sh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Batch scripts used to build and run the simulation on Expanse.

synapses.py
^^^^^^^^^^^
Contains functionality called in the run step of the simulation. Used to set the parameters for different synapses, including things like weight and release probability.


:doc:`build_network.py <build_network>`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Contains the class SimulationBuilder which performs the build step of the simulation. This uses BMTK to create files that are then loaded in for the simulation step.

run_network.py
^^^^^^^^^^^^^^

.. automodule:: L5NeuronSimulation.FullSimulation.run_network
   :members:
