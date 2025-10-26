#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --nodelist=ant13

#SBATCH --partition=vl-parcio

#SBATCH --job-name=measurement1
#SBATCH --output=measurement1.out
#SBATCH --error=measurement1.err

. /opt/spack/pp-2025/env.sh
module load hyperfine

echo "===== Hardware Info ====="
hostname
lscpu
free -h
echo "========================="

for i in 1 2 3 6 12 18 24
do
    srun hyperfine -r 3 "../pde/partdiff $i 2 4096 2 75"
done
 