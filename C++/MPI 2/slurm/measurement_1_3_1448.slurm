#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --exclusive
#SBATCH --partition=vl-parcio
#SBATCH --constraint=epyc3
#SBATCH --job-name=measurement_1_3_1448
#SBATCH --output=measurement_1_3_1448.out
#SBATCH --error=measurement_1_3_1448.err

. /opt/spack/pp-2025/env.sh
module load hyperfine

echo "=== Configuration (1, 3, 1448) - 3 Processes ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo ""

echo "=== Hardware Specifications ==="
echo "CPU Information:"
lscpu | grep -E "Model name|CPU\(s\)|Thread|Core|Socket|Cache"
echo ""
echo "Memory Information:"
free -h
echo ""
echo "Node Details:"
scontrol show node $(hostname) | grep -E "CPUTot|RealMemory|Features"
echo ""

echo "=== Performance Measurement ==="
cd ../pde
hyperfine -r 3 'srun --mpi=pmi2 ./partdiff 3 2 1448 2 80'
