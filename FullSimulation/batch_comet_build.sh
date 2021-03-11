#!/bin/bash
#SBATCH --partition shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -A TG-DBS180005
#SBATCH --job-name=build
#SBATCH --output=build.out
#SBATCH --time 0-00:05


module purge
#module load python
module load intel
module load openmpi_ib
export PYTHONPATH=$HOME/neuron/nrn/lib/python:$PYTHONPATH
export LD_LIBRARY_PATH=$HOME/neuron/nrn/x86_64/lib:$LD_LIBRARY_PATH
export PATH=$HOME/neuron/nrn/x86_64/bin:$PATH

rm exc_stim_spikes.h5
rm inh_stim_spikes.h5
rm -rf network/*

echo "Building model at $(date)"

python build_network.py

echo "Done building model at $(date)"
