# Simple Makefile for timempi MPI program

# Compiler
CC = mpicc

# Compiler flags
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -O3

# Executable name
TARGET = timempi

# Build rule
all: $(TARGET)

$(TARGET): timempi.c
	$(CC) $(CFLAGS) timempi.c -o $(TARGET)

# Clean rule
clean:
	rm -f $(TARGET)
