
# Task 1: Memory Error (Valgrind)

---

## Error 1: Stack Return in `mistake1`

### Before Fix:
```
Segmentation fault (core dumped)
```

### Explanation:
```c
int buf[] = {1, 1, 2, 3, 4, 5};
return buf;
```
Returning a pointer to a stack array (`buf`) results in accessing invalid memory.

### Fix:
```c
int* buf = malloc(6 * sizeof(int));
for (int i = 0; i < 6; i++) buf[i] = tmp[i];
return buf;
```

---

## Error 2: Invalid Malloc Size & Free in `mistake2`

### Before Fix:
```
== Invalid free()
== Address 0x... is 4 bytes inside a block of size 16 alloc'd
== definitely lost: 72 bytes in 4 blocks
```

### Explanation:
- `malloc(sizeof(char) * 4)` is too small for 4 integers.
- `free(p[1])` points to the middle of a block.

### Fix:
```c
int* buf = malloc(4 * sizeof(int));
free(p[1] - 1);
```

---

## Error 3: Use-After-Free and OOB in `mistake3`

### Before Fix:
```
Use of uninitialised value
Invalid write of size 4
```

### Explanation:
```c
buf[4] = 3;
free(buf);
return buf;
```
- Writes past allocated memory.
- Returns pointer after it has been freed.

### Fix:
```c
int* buf = malloc(4 * sizeof(int));
buf[2] = 3;
return buf;
```

---

## Error 4: Invalid Function Pointer Cast in `mistake4`

### Before Fix:
```
Segmentation fault
```

### Explanation:
```c
int* buf = (int*)&mistake2;
buf[0] = 4;
```
Function pointer cast to data pointer and modified — undefined behavior.

### Fix:
```c
int* buf = malloc(sizeof(int));
buf[0] = 4;
```

---

## Before and After Memory Fix Results

### Before Fix:
```valgrind ./simple
=2100860== Memcheck, a memory error detector
==2100860== Copyright (C) 2002-2024, and GNU GPL'd, by Julian Seward et al.
==2100860== Using Valgrind-3.23.0 and LibVEX; rerun with -h for copyright info
==2100860== Command: ./simple
==2100860==
1: 1
2: 2
3: 3
4: 4
==2100860== Invalid free() / delete / delete[] / realloc()
==2100860==    at 0x48489FC: free (vg_replace_malloc.c:989)
==2100860==    by 0x4012CA: main (simple.c:62)
==2100860==  Address 0x4a710a4 is 4 bytes inside a block of size 16 alloc'd
==2100860==    at 0x48457A3: malloc (vg_replace_malloc.c:446)
==2100860==    by 0x4011CF: mistake2 (simple.c:21)
==2100860==    by 0x401248: main (simple.c:54)
==2100860==
==2100860==
==2100860== HEAP SUMMARY:
==2100860==     in use at exit: 72 bytes in 4 blocks
==2100860==   total heap usage: 5 allocs, 2 frees, 1,096 bytes allocated
==2100860==
==2100860== LEAK SUMMARY:
==2100860==    definitely lost: 72 bytes in 4 blocks
==2100860==    indirectly lost: 0 bytes in 0 blocks
==2100860==      possibly lost: 0 bytes in 0 blocks
==2100860==    still reachable: 0 bytes in 0 blocks
==2100860==         suppressed: 0 bytes in 0 blocks
==2100860== Rerun with --leak-check=full to see details of leaked memory
==2100860==
==2100860== For lists of detected and suppressed errors, rerun with: -s
==2100860== ERROR SUMMARY: 1 errors from 1 contexts (suppressed: 0 from 0)
```

### After Fix:
```simple.c
        free(p[0] - 1); // correct base pointer for mistake1
	    free(p[1] - 1); // correct base pointer for mistake2
	    free(p[2]);     // directly returned base from mistake3
	    free(p[3]); 
```
```valgrind ./simple
==2100932== Memcheck, a memory error detector
==2100932== Copyright (C) 2002-2024, and GNU GPL'd, by Julian Seward et al.
==2100932== Using Valgrind-3.23.0 and LibVEX; rerun with -h for copyright info
==2100932== Command: ./simple
==2100932==
1: 1
2: 2
3: 3
4: 4
==2100932==
==2100932== HEAP SUMMARY:
==2100932==     in use at exit: 0 bytes in 0 blocks
==2100932==   total heap usage: 5 allocs, 5 frees, 1,096 bytes allocated
==2100932==
==2100932== All heap blocks were freed -- no leaks are possible
==2100932==
==2100932== For lists of detected and suppressed errors, rerun with: -s
==2100932== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)
```
