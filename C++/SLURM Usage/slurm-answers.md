# SLURM Command Explanation

## Question
**Why are four threads allocated if you request three with the `-c` parameter? Write a few sentences to explain this behavior.**  

## Explanation
SLURM always allocates whole CPU cores, not fractional ones.  
With hyper-threading enabled, each physical core provides two hardware threads.  
If you request 3 threads (`-c 3`), SLURM cannot allocate 1.5 cores, so it rounds up to 2 full cores.  
2 cores × 2 threads per core = **4 threads** in total.

## Demonstration Command
```bash
srun -p vl-parcio -c 4 nproc
```

### Breakdown
- **`srun`**: Launches a SLURM job step, creating an implicit allocation if needed.  
- **`-p vl-parcio`**: Specifies the partition (compute nodes) named `vl-parcio`.  
- **`-c 4`**: Requests 4 CPUs (hardware threads) for the task. SLURM allocates whole cores (2 threads each), so you get 2 physical cores = 4 threads.  
- **`nproc`**: Prints the number of available processing units; under SLURM with 4 threads reserved, it outputs `4`.

### Expected Output
```
4
```
indicating that SLURM has allocated four hardware threads.
