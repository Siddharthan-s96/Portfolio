# Parallel Execution Observations

## Question 1
**What do you notice from the timestamps? Try to explain your observations!**

### Answer
- **Parallel execution per node:**  
  You see **four** nearly-identical timestamps for each node (e.g., `ant14`), separated by only nanoseconds. This proves the four tasks on each node truly ran **in parallel**.
- **Inter-node timing skew:**  
  The timestamps between different nodes differ by small offsets (micro- to nanoseconds), indicating slight scheduling delays or clock skew across nodes.
- **Repeatability:**  
  Running the script multiple times reproduces the same pattern (groups of parallel timestamps), but at different wall-clock times.

### Sample Output (`timescript.out`)
```bash
ant14: 2025-05-10T23:52:40.237212314+02:00
ant14: 2025-05-10T23:52:40.252890604+02:00
ant14: 2025-05-10T23:52:40.252908974+02:00
ant14: 2025-05-10T23:52:40.252910844+02:00
ant15: 2025-05-10T23:52:40.237213214+02:00
ant15: 2025-05-10T23:52:40.244336948+02:00
ant15: 2025-05-10T23:52:40.244337008+02:00
ant15: 2025-05-10T23:52:40.244348198+02:00
...
(total 16 lines)
```

## Question 2
**Could the `timescript.out` file also be created inside the `timescript.sh` script? If so: How? If not: Why not?**

### Answer
Yes, it is technically possible to have each invocation of `timescript.sh` append its own output directly into `timescript.out`. For example, you could modify `timescript.sh` like this:

```bash
#!/bin/bash
HOST=$(hostname --short)
TIME=$(date --iso-8601=ns)
echo "$HOST: $TIME" | tee -a timescript.out
```

or simply:

```bash
echo "$HOST: $TIME" >> timescript.out
```

However, this approach has significant drawbacks:

1. **Race Conditions:**  
   When multiple processes write to the same file simultaneously, lines can overlap or get lost. Without file locking, the output may be corrupted or interleaved.

2. **Performance and I/O Consistency:**  
   On many clusters, compute nodes use local scratch that buffers writes before syncing to a shared filesystem. Concurrent appends can lead to missing or delayed entries in the final shared file.

3. **Maintainability:**  
   Embedding file management inside every instance of the script mixes concerns (business logic vs. logging) and makes debugging harder.

**Recommended Practice:**  
Use SLURM’s `#SBATCH --output=timescript.out` directive in your job script (`jobscript.sh`) to let SLURM handle output collection safely. SLURM serializes each task’s stdout, avoiding race conditions and ensuring consistent, ordered output without modifying `timescript.sh`.  
