#!/bin/bash

#SBATCH --job-name=mpi-circle
#SBATCH --nodes=2
#SBATCH --ntasks=4           
#SBATCH --cpus-per-task=1     
#SBATCH --output=circle.out
#SBATCH --error=circle.err

. /opt/spack/pp-2025/env.sh        

module load openmpi                

# Run the MPI program
srun ./circle 8
