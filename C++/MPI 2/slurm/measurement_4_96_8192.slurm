#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=96
#SBATCH --ntasks-per-node=24
#SBATCH --exclusive
#SBATCH --partition=vl-parcio
#SBATCH --constraint=epyc3
#SBATCH --job-name=measurement_4_96_8192
#SBATCH --output=measurement_4_96_8192.out
#SBATCH --error=measurement_4_96_8192.err

. /opt/spack/pp-2025/env.sh
module load hyperfine

echo "=== Configuration (4, 96, 8192) - 96 Processes on 4 Nodes ==="
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
hyperfine -r 3 'srun --mpi=pmi2 ./partdiff 96 2 8192 2 80'
