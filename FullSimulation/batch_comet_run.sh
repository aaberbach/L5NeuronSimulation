#!/bin/bash

#SBATCH --partition shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -A TG-DBS180005
#SBATCH --job-name=run
#SBATCH --output=run.out
#SBATCH --time 0-00:05

module purge
#module load python
module load intel
module load openmpi_ib
export PATH=$HOME/neuron/nrn/x86_64/bin:$PATH
export LD_LIBRARY_PATH=$HOME/neuron/nrn/x86_64/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$HOME/neuron/nrn/lib/python:$PYTHONPATH
export PATH=$HOME/neuron/nrn/x86_64/bin:$PATH

rm -rf output


echo "Running model at $(date)"

#mpirun nrniv -mpi -quiet -python3 run_network.py simulation_config.json
#ibrun nrniv -mpi -python run_network.py
python run_network.py

echo "Done running model at $(date)"
