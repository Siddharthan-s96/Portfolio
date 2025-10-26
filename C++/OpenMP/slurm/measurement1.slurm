#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48

#SBATCH --partition=vl-parcio

#SBATCH --job-name=measurement1
#SBATCH --output=measurement1.out
#SBATCH --error=measurement1.err

. /opt/spack/pp-2025/env.sh
module load hyperfine

# SERIAL version (reference)
srun hyperfine -r 3 "../pde/partdiff 1 2 4096 2 5"

# ROW-wise parallel version
for i in {1..24}; do
    srun hyperfine -r 3 "OMP_NUM_THREADS=$i ../pde/partdiff-row $i 2 4096 2 5"
done

# COLUMN-wise parallel version
for i in {1..24}; do
    srun hyperfine -r 3 "OMP_NUM_THREADS=$i ../pde/partdiff-column $i 2 4096 2 5"
done

# ELEMENT-wise parallel version
for i in {1..24}; do
    srun hyperfine -r 3 "OMP_NUM_THREADS=$i ../pde/partdiff-element $i 2 4096 2 5"
done