FullSimulation
==============

This is the folder where the core simulation is built and run.

Files
-----

Connections.csv
^^^^^^^^^^^^^^^
Contains information for each synapse in the simulation. Created during the run step. Columns:

* Node ID: the integer id of the presynaptic node within its population (exc, prox_inh, dist_inh)
* Distance: float representing the distance (um) from the synapse to the soma
* Conductance: the weight of the synapse (already scaled by distance)
* Type: what part of the cell the synapse is on (soma, dend, apic)
* Name: the full string that NEURON associates with the postsynaptic segment.
* Source Population: the node population that the presynaptic node is a member of (exc, prox_inh, dist_inh)
* Release Probability: release probability for the synapse
