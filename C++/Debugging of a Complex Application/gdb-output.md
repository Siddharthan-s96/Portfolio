# Task 1: GDB Debugging 

---

## Subtask 1: Breakpoint on `mistake1`

### Commands and Output:
```gdb
(gdb) break mistake1
Breakpoint 1 at 0x40114a: file simple.c, line 12.
(gdb) run
Breakpoint 1, mistake1 () at simple.c:12
12              int buf[] = { 1, 1, 2, 3, 4, 5 };
(gdb) print buf
$1 = {23, 0, 0, 0, 0, 0}
(gdb) print buf[2]
$2 = 0
(gdb) whatis buf
type = int [6]
(gdb) next
14              return buf;
(gdb) print buf
$3 = {1, 1, 2, 3, 4, 5}
(gdb) print buf[2]
$4 = 2
```

---

## Subtask 2: Breakpoint on `mistake2`

### Commands and Output:
```gdb
(gdb) break mistake2
Breakpoint 2 at 0x401183: file simple.c, line 20.
(gdb) run
Start it from the beginning? (y or n) y
Breakpoint 1, mistake1 () at simple.c:12
(gdb) continue
Breakpoint 2, mistake2 () at simple.c:20
20              int* buf = malloc(sizeof(char) * 4);
(gdb) print buf
$5 = (int *) 0x0
(gdb) print buf[2]
Cannot access memory at address 0x8
(gdb) whatis buf
type = int *
```

---

## Subtask 3: Continue execution, show output, stack frames, switch to frame 1, print `p`

### Commands and Output:
```gdb
(gdb) continue
Program received signal SIGSEGV, Segmentation fault.
0x00000000004011eb in mistake4 () at simple.c:44
(gdb) list
...
(gdb) backtrace
...
(gdb) frame 1
...
(gdb) print p
...
```

---

## Subtask 4: Call `mistake3`

### Commands and Output:
```gdb
(gdb) call mistake3()
$7 = (int *) 0x4052c0
```
