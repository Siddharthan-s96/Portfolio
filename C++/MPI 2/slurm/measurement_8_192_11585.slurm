#!/bin/bash
#SBATCH --nodes=8
#SBATCH --ntasks=192
#SBATCH --ntasks-per-node=24
#SBATCH --exclusive
#SBATCH --partition=vl-parcio

#SBATCH --job-name=measurement_8_192_11585
#SBATCH --output=measurement_8_192_11585.out
#SBATCH --error=measurement_8_192_11585.err

. /opt/spack/pp-2025/env.sh
module load hyperfine

echo "=== Configuration (8, 192, 11585) - 192 Processes on 8 Nodes ==="
echo "Date: $(date)"
echo "Nodes: $(srun hostname | sort | uniq | tr '\n' ' ')"
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
hyperfine -r 3 'srun --mpi=pmi2 ./partdiff 192 2 11585 2 80'
