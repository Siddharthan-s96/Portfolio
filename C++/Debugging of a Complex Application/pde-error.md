# Task 2: Debugging Report – PDE Solver (`partdiff.c`)

---

## Bug 1: Incorrect Argument Access

- **Location:** `askParams()` function, Line ~86
- **Issue:** Incorrectly accessed `argv[333]`, leading to out-of-bounds memory access.
- **Original Code:**
  ```c
  ret = sscanf(argv[333], "%" SCNu64, &(options->interlines));
  ```
- **Fix:**
  ```c
  ret = sscanf(argv[3], "%" SCNu64, &(options->interlines));
  ```

---

## Bug 2: Insufficient Memory Allocation for Matrix

- **Location:** `allocateMatrices()` function, Line ~164
- **Issue:** Allocated memory for `(N+1)*(N-1)` matrix but used `(N+1)*(N+1)` in access.
- **Original Code:**
  ```c
  arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N - 1) * sizeof(double));
  ```
- **Fix:**
  ```c
  arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
  ```

---

## Bug 3: Invalid Matrix Indexing in Core Calculation

- **Location:** `calculate()` function, Line ~307
- **Issue:** Swapped matrix indexing caused a segmentation fault.
- **Original Code:**
  ```c
  star = (Matrix[m2][i - 1][j] + Matrix[j - 1][m2][i] + Matrix[m2][i][j + 1] + Matrix[m2][i + 1][j]) / 4;
  ```
- **Fix:**
  ```c
  star = (Matrix[m2][i - 1][j] + Matrix[m2][i][j - 1] + Matrix[m2][i][j + 1] + Matrix[m2][i + 1][j]) / 4;
  ```

---

## Bug 4: Memory Leak (Valgrind)

- **Location:** `main()` function
- **Issue:** Allocated matrix memory was never freed.
- **Fix:**
  ```c
  free(arguments.M);
  ```

After this fix, Valgrind reports:
```
All heap blocks were freed -- no leaks are possible
```

---

## Unused Variable

1. **Unused parameter `options`**
   - In `initMatrices` and `calculate_func`
   - Fixed using:
     ```c
     (void)options;
     ```

2. **Unused variable `h`** in `initMatrices`
   - Simply removed

---

## Memory Leak

- **Location:** `main()` function
- **Issue:** Allocated matrix memory was never freed.

### Before Fix:
```bash
==2105290== HEAP SUMMARY:
==2105290==     in use at exit: 5,235,848 bytes in 1 blocks
==2105290==   total heap usage: 2 allocs, 1 frees, 5,236,872 bytes allocated

==2105290== LEAK SUMMARY:
==2105290==    definitely lost: 5,235,848 bytes in 1 blocks
==2105290==    indirectly lost: 0 bytes in 0 blocks
==2105290==    still reachable: 0 bytes in 0 blocks
```

### Fix:
```c
free(arguments.M);
```

### After Fix:
```bash
==2105378== All heap blocks were freed -- no leaks are possible
==2105378== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)
```