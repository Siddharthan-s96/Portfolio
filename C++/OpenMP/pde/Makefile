CFLAGS  = -std=c11 -Wall -Wextra -Wpedantic -O3 -g -fopenmp
LDFLAGS = -fopenmp
LDLIBS  = -lm

all: partdiff partdiff-row partdiff-column partdiff-element

partdiff: partdiff.c
	$(CC) $(CFLAGS) -DDEFAULT $(LDFLAGS) $^ $(LDLIBS) -o $@

partdiff-row: partdiff.c
	$(CC) $(CFLAGS) -DROW $(LDFLAGS) $^ $(LDLIBS) -o $@

partdiff-column: partdiff.c
	$(CC) $(CFLAGS) -DCOLUMN $(LDFLAGS) $^ $(LDLIBS) -o $@

partdiff-element: partdiff.c
	$(CC) $(CFLAGS) -DELEMENT $(LDFLAGS) $^ $(LDLIBS) -o $@

clean:
	$(RM) partdiff partdiff-row partdiff-column partdiff-element

