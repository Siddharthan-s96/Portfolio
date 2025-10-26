#!/bin/bash
#SBATCH -N 4                          # Request 4 nodes
#SBATCH -n 16                         # Total 16 tasks (4 per node)
#SBATCH --job-name=timescript_job
#SBATCH --output=timescript.out       # Output from all processes
#SBATCH --time=00:05:00

srun ./timescript.sh

echo "done" >> jobscript.out
