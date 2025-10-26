# Task 3: PDE Optimizations Results

This document summarizes the **numerical correctness** and **performance improvements** of the `partdiff` solver across three optimization stages. All matrix outputs remain **identical** at every level.

## Compiler Optimization Flags

* \`\`: Disables compiler optimizations, generating straightforward code that closely follows the source. This setting preserves debug information and is used for baseline performance measurements and initial profiling.
* \`\`: Enables aggressive optimizations—including function inlining, loop unrolling, vectorization, and dead-code elimination—to maximize execution speed in production builds.

### Build Command and Flags

To compile with baseline settings (profile- and debug-enabled), use:

```makefile
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic **-O0** -g -pg
```

The `-O0` ensures minimal optimization for initial profiling; `-g` enables debugging; `-pg` enables `gprof` instrumentation.
---

## 1. Initial Implementation (O0)

### 1.1 Timing
```bash
$ time ./partdiff 1 2 64 2 5120
```
**Highlights**
- Calculation time: **82.734 s**
- Real elapsed: `1m22.744s`

```bash
real    1m22.744s
user    1m22.721s
sys     0m0.008s
```

### 1.2 Profiling with gprof
```bash
$ gprof ./partdiff gmon.out
```
**Flat profile (excerpt)**
| % time | self (sec) | calls         | name           |
|-------:|-----------:|--------------:|:---------------|
|  66.56 |      24.26 | 1             | calculate      |
|  33.28 |      12.13 | 1,379,128,320 | calculate_func |
|   0.16 |       0.06 | —             | _init          |
> **Explanation:** The nested loops in `calculate` and the two `sin()` calls per cell in `calculate_func` are the dominant cost.

---

## 2. Level 1: CPU Binding Allocation

Acquire dedicated hardware threads by binding your job to exactly **4** CPU threads on a compute node:
```bash
$ srun -p vl-parcio -c 4 --pty bash
```
This requests **4 logical processors** (hardware threads), which SLURM allocates as 2 physical cores with hyperthreading enabled. By reserving these dedicated threads, you eliminate OS scheduling variability.

After allocation, rerun under O0 and profile:
```bash
$ time ./partdiff 1 2 64 2 5120
$ gprof ./partdiff gmon.out
```

**Flat profile (excerpt)**
| % time | self (sec) | calls         | name           |
|-------:|-----------:|--------------:|:---------------|
|  71.22 |      15.44 | 1             | calculate      |
|  26.01 |       5.64 | 1,379,128,320 | calculate_func |

> **Explanation:**  
> The `-c 4` flag reserved 4 threads (2 cores × 2 hyperthreads), ensuring your process ran strictly on these threads. This isolation removes OS scheduling jitter but leaves the computational hotspots (`calculate`, `calculate_func`) unchanged.

---

## 3. Manual Optimizations in `calculate` Function

The following **manual** changes were applied directly inside `calculate` to eliminate redundant work and function-call overhead:

### 3.1 Precompute constant expressions
```c
double const h     = arguments->h;
double const pi_h  = M_PI * h;
double const coeff = (2 * M_PI * M_PI * h * h) / 4.0;
double const inv4  = 0.25;
```
> Moved repeated arithmetic (π·h, π²·h², 1/4) out of inner loops.

### 3.2 Build sine lookup table
```c
double* sin_table = malloc((N+1) * sizeof *sin_table);
for (int k = 0; k <= N; ++k)
    sin_table[k] = sin(pi_h * k);
```
> Replaced two `sin()` calls per grid cell with two array loads.

### 3.3 Inline source-term computation
Old:
```c
star += calculate_func(arguments, i, j);
```
New:
```c
star = inv4 * sum4
     + coeff * sin_table[i] * sin_table[j];
```
> Eliminated function call and merged its arithmetic into the core loop.

### 3.4 Restrict-qualified pointers
```c
double (*__restrict__ src)[N+1] = Matrix[m2];
double (*__restrict__ dst)[N+1] = Matrix[m1];
```
> Guarantees no aliasing, enabling efficient memory access.

### 3.5 Method-specific loops & buffer swap
- **Jacobi** (two-buffer “ping-pong”): uses `m1 ^= 1; m2 ^= 1;`
- **Red-Black Gauss–Seidel** (in-place, two color passes)  
> Tailored loops remove unnecessary dependencies and extra swaps.

### 3.6 Simplified termination logic
```c
if (options->termination == TERM_PREC && maxres < options->term_precision) break;
else --term_iter;
```
> Replaced counter resets with direct `break` and clean decrement.

---

## 4. Level 2: Manual Inlining & Lookup Table

By combining these manual changes, the costly function calls and heavy math were removed from the nested loops, resulting in a dramatic speedup under O0.

### 4.1 Timing
```bash
$ time ./partdiff 1 2 64 2 5120
```
**Highlights**
- Calculation time: **10.246 s**
- Real elapsed: `0m10.255s`

```bash
real    0m10.255s
user    0m10.211s
sys     0m0.002s
```

---

## 5. Level 3: Full Compiler Optimization (-O3)

### 5.1 Compilation & Run
```bash
$ gcc -O3 -g -pg -o partdiff partdiff.c
$ time ./partdiff 1 2 64 2 5120
```
**Highlights**
- Calculation time: **1.385 s**
- Real elapsed: `0m1.392s`

```bash
real    0m1.392s
user    0m1.381s
sys     0m0.002s
```

### 5.2 Profiling with gprof
```bash
$ gprof ./partdiff gmon.out
```
**Flat profile**
| % time | self (sec) | name |
|-------:|-----------:|-----:|
| 100.00 |       1.38 | main |
> **Explanation:**  
> Under `-O3`, the compiler fully **inlined** both `calculate` and `calculate_func`—the two major helper functions—directly into `main`, and eliminated all dead code. This consolidation merged every loop and arithmetic operation into a single symbol, so `gprof` reports **100%** of execution time in `main` with no separate entries for other functions. As a result, `gprof` reports **100%** of the execution time in `main`, with no profiling data for other functions.

---

## 6. Key Mathematical Changes

- **Division → Multiplication**
```c
// Old: sum4 / 4
// New: inv4 * sum4    // inv4 = 0.25

// Old: (2·π²·h²) * sin()*sin() / 4
// New: coeff * sin_table[i] * sin_table[j]
// where coeff = (2 * M_PI * M_PI * h * h) / 4.0
```
- **Sin calls → Table lookup**
```c
// Old:
double term = sin(pi_h * i) * sin(pi_h * j);

// New:
double term = sin_table[i] * sin_table[j];
```
> These mathematical optimizations remove heavy arithmetic and function calls inside the loops, replacing them with simple multiplies and array accesses.

---

## 7. Performance Summary

| Stage            | Real Time | Speedup vs O0 |
|-----------------:|----------:|--------------:|
| Initial (O0)     | 82.744 s  |       1×      |
| Level 2 (manual) | 10.255 s  |     ~8.1×     |
| Level 3 (-O3)    | 1.392 s   |    ~59.4×     |

All numerical results (matrix outputs) are unchanged. The combination of manual and compiler optimizations yields dramatic performance gains without altering correctness.
