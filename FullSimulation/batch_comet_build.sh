#!/bin/bash
#SBATCH --partition shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --account=umc110
#SBATCH --job-name=build
#SBATCH --mem=8G
#SBATCH --output=build.out
#SBATCH --time 0-02:00


module purge
module load slurm
module load cpu
module load intel
module load intel-mpi
module load ncurses
export PYTHONPATH=$HOME/nrn/lib/python:$PYTHONPATH
export LD_LIBRARY_PATH=$HOME/nrn/x86_64/lib:$LD_LIBRARY_PATH
export PATH=$HOME/nrn/x86_64/bin:$PATH

rm exc_stim_spikes.h5
rm inh_stim_spikes.h5
rm -rf network/*

echo "Building model at $(date)"

python build_network.py

echo "Done building model at $(date)"
