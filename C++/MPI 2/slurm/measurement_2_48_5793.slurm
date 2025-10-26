#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=24
#SBATCH --exclusive
#SBATCH --partition=vl-parcio
#SBATCH --constraint=epyc3
#SBATCH --job-name=measurement_2_48_5793
#SBATCH --output=measurement_2_48_5793.out
#SBATCH --error=measurement_2_48_5793.err

. /opt/spack/pp-2025/env.sh
module load hyperfine

echo "=== Configuration (2, 48, 5793) - 48 Processes on 2 Nodes ==="
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
hyperfine -r 3 'srun --mpi=pmi2 ./partdiff 48 2 5793 2 80'
