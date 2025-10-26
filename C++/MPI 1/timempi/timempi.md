# MPI Timempi Program Report

## Overview

This report documents the implementation and execution of an MPI (Message Passing Interface) program called `timempi`. The program demonstrates distributed computing concepts by having multiple processes collect hostname and timestamp information, coordinate their output through process 0, and synchronize their termination.

## Solution Architecture

### Program Flow

1. **Initialization**: All processes initialize MPI and gather their local hostname and timestamp
2. **Data Collection**: Processes with ranks 1 to n send their hostname\:timestamp strings to rank 0
3. **Centralized Output**: Rank 0 receives and prints all strings in rank order
4. **Reduction Operation**: All processes participate in finding the minimum microsecond value
5. **Synchronization**: All processes wait at a barrier before termination
6. **Termination**: Each process prints its termination message

### Key MPI Concepts Demonstrated

* **Point-to-point communication**: `MPI_Send()` and `MPI_Recv()`
* **Collective operations**: `MPI_Reduce()` for finding minimum values
* **Synchronization**: `MPI_Barrier()` for coordinating process termination
* **Process identification**: `MPI_Comm_rank()` and `MPI_Comm_size()`

## Implementation Details

### Source Code Structure

```c
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <mpi.h>

#define MAX_STRING_LENGTH 512
```

### Core Functionality

#### 1. Time and Hostname Collection

```c
gettimeofday(&tv, NULL);
gethostname(hostname, MAX_STRING_LENGTH);
tm_info = localtime(&tv.tv_sec);
strftime(buffer, 26, "%Y-%m-%d %H:%M:%S", tm_info);
microseconds = tv.tv_usec;
```

#### 2. Message Passing (Ranks 1 to n)

```c
snprintf(output_string, MAX_STRING_LENGTH, "%s: %s.%06ld", hostname, buffer, microseconds);
MPI_Send(output_string, strlen(output_string) + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
```

#### 3. Centralized Output (Rank 0)

```c
for (int i = 1; i < size; i++) {
    MPI_Recv(received_string, MAX_STRING_LENGTH, MPI_CHAR, i, 0, 
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("%s\n", received_string);
}
```

#### 4. Reduction Operation

```c
MPI_Reduce(&microseconds, &min_microseconds, 1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
```

#### 5. Synchronization

```c
MPI_Barrier(MPI_COMM_WORLD);
printf("Rank %d terminating now!\n", rank);
```

## Build System

### Makefile

```makefile
CC      = mpicc
CFLAGS  = -std=c11 -Wall -Wextra -Wpedantic -O3 -g

timempi: timempi.c

clean:
	$(RM) timempi
```

The Makefile uses:

* `mpicc`: MPI C compiler wrapper
* Strict compilation flags for code quality
* Simple implicit rules for building

## Execution Environment

### SLURM Job Script

```bash
#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --partition=vl-parcio
#SBATCH --job-name=timempi
#SBATCH --output=timempi.out
#SBATCH --error=timempi.err

. /opt/spack/pp-2025/env.sh
module avail openmpi
srun ./timempi
```

### Execution Parameters

* **Nodes**: 4 compute nodes
* **Tasks**: 4 MPI processes (1 per node)
* **CPUs per task**: 1
* **Launcher**: `srun` for SLURM integration

## Results Analysis

### Program Output

```
ant18: 2025-06-14 20:57:38.845815
ant19: 2025-06-14 20:57:38.845702
ant20: 2025-06-14 20:57:38.845832
Rank 1 terminating now!
Rank 3 terminating now!
Rank 2 terminating now!
845688
Rank 0 terminating now!
```

### Output Analysis

1. **Hostname and Timestamp Display**:

   * `ant18`: Process rank 1 on node ant18
   * `ant19`: Process rank 2 on node ant19
   * `ant20`: Process rank 3 on node ant20
   * All timestamps show execution at 20:57:38 on 2025-06-14

2. **Microsecond Precision**:

   * Rank 1: 845815 microseconds
   * Rank 2: 845702 microseconds
   * Rank 3: 845832 microseconds

3. **Minimum Microsecond Value**: 845688

   * This represents the minimum across all processes (including rank 0)
   * Rank 0's microsecond value (845688) was the smallest

4. **Process Termination Order**:

   * Non-deterministic order due to parallel execution
   * All processes wait at barrier before printing termination messages
   * Rank 0 terminates last, after printing the reduction result

### Key Observations

* **Synchronization Success**: All processes properly coordinated their execution
* **Ordered Output**: Hostname\:timestamp strings displayed in rank order (1, 2, 3)
* **Correct Reduction**: Minimum microsecond value correctly computed and displayed
* **Proper Termination**: All processes synchronized before termination

## Technical Implementation Notes

### Compilation Considerations

* Used `#define _GNU_SOURCE` to enable `gethostname()` function
* Replaced `sprintf` with `snprintf` for buffer safety
* Increased buffer size to 512 bytes to prevent overflow warnings

### MPI Communication Pattern

* **Fan-in pattern**: Multiple processes send to rank 0
* **Collective operation**: All processes participate in reduction
* **Barrier synchronization**: Ensures coordinated termination

## Conclusion

The `timempi` program successfully demonstrates fundamental MPI concepts including point-to-point communication, collective operations, and process synchronization. The execution on 4 nodes shows proper distributed behavior with coordinated output and correct reduction operations. The program meets all specified requirements and provides a solid foundation for understanding MPI programming patterns.