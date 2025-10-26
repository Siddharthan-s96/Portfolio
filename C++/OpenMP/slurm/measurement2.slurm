#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48

#SBATCH --partition=vl-parcio

#SBATCH --job-name=measurement2
#SBATCH --output=measurement2.out
#SBATCH --error=measurement2.err

. /opt/spack/pp-2025/env.sh
module load hyperfine

for i in 1 2 4 8 16 32 64 128 256 512 1024 2048 4096
do
    srun hyperfine -r 3 "../pde/partdiff-row 24 2 $i 2 5000"
done
