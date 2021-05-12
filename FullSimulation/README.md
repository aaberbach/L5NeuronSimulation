# Full Simulation Folder

This is the folder where the core simulation is run.

## Files

batch_comet_build.sh/batch_comet_run.sh - scripts used to build and run the simulation remotely.

build_network.py - Contains the class SimulationBuilder which performs the build step of the simulation. This uses BMTK to create files that are then loaded infor the simulation step.

Connections.csv - Contains information for each synapse in the simulation. Created during the run step. Columns:
- Node ID: the integer id of the presynaptic node within its population (exc, prox_inh, dist_inh)
- Distance: float representing the distance (um) from the synapse to the soma
- Conductance: the weight of the synapse (already scaled by distance)
- Type: what part of the cell the synapse is on (soma, dend, apic)
- Name: the full string that NEURON associates with the postsynaptic segment. 
- Source Population: the node population that the presynaptic node is a member of (exc, prox_inh, dist_inh)
- Release Probability: release probability for the synapse

FunctionalGroups.csv - contains the functional group id for each excitatory presynaptic node

NetParams.json - contains parameters for the simulation that SimulationBuilder uses

output_FRs.txt - contains firing rate for various combinations of synaptic weight distributions **REMOVE**

run_network.py - runs the simulation and save Connections.csv

Segments.csv - contains information about every segment in the morphology

show_results - used to plot and print various aspects of the simulation **CAN REMOVE IF NOT NEEDED**

synapses.py- used to set the correct parameters for the synapses during the run step