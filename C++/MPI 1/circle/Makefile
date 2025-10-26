CC      = mpicc
CFLAGS  = -std=c11 -Wall -Wextra -Wpedantic -O3 -g

# Target for the program
circle: circle.c
	$(CC) $(CFLAGS) -o circle circle.c

# Remove generated files
clean:
	rm -f circle
