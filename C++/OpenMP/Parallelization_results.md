# OpenMP Parallelization and Data Distribution Results

## Task: Parallelization with OpenMP (180 Points)

**Command Used:**
```bash
./partdiff 1 2 1024 2 100
```

**Results:**
- Calculation Time: `86.197204 s`
- Memory Usage: `1026.251236 MiB`
- Method: Jacobi
- Interlines: 1024
- Termination: Number of Iterations
- Iterations: 100
- Error Precision: `7.339034e-08`

---

**Parallel Row-wise Version**
```bash
./partdiff-row 4 2 1024 2 100
```
- Time: `34.437046 s`

---

**Parallel Column-wise Version**
```bash
./partdiff-column 4 2 1024 2 100
```
- Time: `34.413871 s`

---

**Parallel Element-wise Version**
```bash
./partdiff-element 4 2 1024 2 100
```
- Time: `74.294054 s`

---

## Task: Data Distribution (60 Bonus Points)

### Comparison of Versions (with 1024 Interlines, 100 Iterations)

| Version        | Time (s)     | Memory Usage (MiB) | Error Precision    |
|----------------|--------------|---------------------|---------------------|
| Serial         | 86.20        | 1026.25             | 7.339034e-08        |
| Row-wise       | 34.44        | 1026.25             | 7.339034e-08        |
| Column-wise    | 34.41        | 1026.25             | 7.339034e-08        |
| Element-wise   | 74.29        | 1026.25             | 7.339034e-08        |

### Observations

- **Speedup**: Row-wise and Column-wise provide nearly identical and best performance.
- **Element-wise**: Slower than other parallel strategies, likely due to higher thread overhead or fine-grain contention.
- **Correctness**: All versions yield the same final matrix and precision.

---

_Compiled on: 2025-06-01 16:35:18_  
