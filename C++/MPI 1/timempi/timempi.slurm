#!/bin/bash

#SBATCH --nodes=4
#SBATCH --ntasks=4                  # 4 MPI processes
#SBATCH --cpus-per-task=1           # 1 CPU per MPI task 

#SBATCH --partition=vl-parcio

#SBATCH --job-name=timempi
#SBATCH --output=timempi.out
#SBATCH --error=timempi.err

. /opt/spack/pp-2025/env.sh         # Load environment if needed

module avail openmpi                 # Load MPI if needed on your cluster

srun ./timempi
