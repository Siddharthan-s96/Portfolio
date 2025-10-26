/*
 * circle.c - MPI Ring Communication Assignment (Parallel Programming)
 * 
 * This program implements ring communication with MPI: an array is split into
 * local chunks, each held by a single process. At each iteration, each process
 * sends its chunk to the next process in the ring. The program stops when the
 * first value of the original array (from rank 0) appears at the start of the
 * last process's chunk. Output is shown before and after, in rank order.
 *
 * Usage:
 *   mpicc -std=c11 -Wall -O3 -o circle circle.c
 *   mpirun -np 4 ./circle 8
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

// Global variables to store array partition information
int GLOBAL_N;         // Total array length
int* local_sizes;     // Number of elements per process
int* displacements;   // Starting index for each process in the global array

/*
 * init(N)
 * -----------
 * - Only rank 0 allocates and fills the full array with random integers [0,24].
 * - All ranks allocate a local buffer for their portion.
 * - The array is distributed using MPI_Scatterv.
 * - local_sizes and displacements arrays are filled for correct data distribution.
 * Returns: Pointer to local chunk of the array.
 */
int* init(int N)
{
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    GLOBAL_N = N;

    // Calculate partition sizes and displacements
    local_sizes = malloc(sizeof(int) * nprocs);
    displacements = malloc(sizeof(int) * nprocs);

    int base_size = N / nprocs;
    int remainder = N % nprocs;
    int offset = 0;
    for (int i = 0; i < nprocs; i++) {
        local_sizes[i] = base_size + (i < remainder ? 1 : 0);
        displacements[i] = offset;
        offset += local_sizes[i];
    }

    int* full_array = NULL;

    // Only rank 0 initializes the full array
    if (rank == 0) {
        full_array = malloc(sizeof(int) * N);
        srand(time(NULL));
        for (int i = 0; i < N; i++) {
            full_array[i] = rand() % 25; // Random values [0,24]
        }
    }

    // Each rank allocates space for its local chunk
    int* local_buf = malloc(sizeof(int) * local_sizes[rank]);

    // Distribute data using MPI_Scatterv
    MPI_Scatterv(full_array, local_sizes, displacements, MPI_INT, 
                 local_buf, local_sizes[rank], MPI_INT, 0, MPI_COMM_WORLD);

    // Free full array on rank 0 after scatter
    if (rank == 0) {
        free(full_array);
    }

    return local_buf;
}

/*
 * circle(buf)
 * -----------
 * - Rotates local array chunks in a ring.
 * - Continues until the first element of the last process equals the original
 *   first element from rank 0's chunk.
 * - Uses MPI_Sendrecv to avoid deadlock and to handle different chunk sizes.
 * Returns: Pointer to the final local chunk (may be a new buffer).
 */
int* circle(int* buf)
{
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Broadcast original first element (from rank 0) to all
    int original_first_element;
    if (rank == 0) {
        original_first_element = buf[0];
    }
    MPI_Bcast(&original_first_element, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int* current_buf = buf;
    int current_size = local_sizes[rank];
    int source_rank = rank; // Track which process's chunk we have

    int iteration = 0;
    int terminate = 0;

    do {
        iteration++;

        int next_rank = (rank + 1) % nprocs;
        int prev_rank = (rank - 1 + nprocs) % nprocs;
        // The chunk we will receive originated from:
        int recv_source_rank = (source_rank - 1 + nprocs) % nprocs;
        int recv_size = local_sizes[recv_source_rank];

        int* recv_buf = malloc(sizeof(int) * recv_size);

        MPI_Sendrecv(current_buf, current_size, MPI_INT, next_rank, 0,
                     recv_buf, recv_size, MPI_INT, prev_rank, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (current_buf != buf) free(current_buf);

        current_buf = recv_buf;
        current_size = recv_size;
        source_rank = recv_source_rank;

        // Only last process checks stop condition
        if (rank == nprocs - 1 && iteration > 1 && current_buf[0] == original_first_element) {
            terminate = 1;
        }
        MPI_Bcast(&terminate, 1, MPI_INT, nprocs - 1, MPI_COMM_WORLD);

    } while (!terminate);

    return current_buf;
}


/*
 * main()
 * ------
 * - Initializes MPI and parses command-line arguments.
 * - Calls init() to scatter the initial array.
 * - Prints the local array for each rank (BEFORE).
 * - Calls circle() to perform ring communication.
 * - Prints the final state for each rank (AFTER).
 * - Frees all memory and finalizes MPI.
 */
int main(int argc, char** argv)
{
    int N;
    int rank, nprocs;
    int* buf;

    // Initialize MPI environment
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Parse array size from command-line
    if (argc < 2) {
        if (rank == 0) {
            printf("Arguments error!\nPlease specify a buffer size.\n");
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    N = atoi(argv[1]);
    buf = init(N);

    // Print BEFORE state (each rank in order, for clean output)
    if (rank == 0) printf("\nBEFORE\n");
    MPI_Barrier(MPI_COMM_WORLD);
    for (int r = 0; r < nprocs; r++) {
        if (rank == r) {
            for (int i = 0; i < local_sizes[rank]; i++) {
                printf("rank %d: %d\n", rank, buf[i]);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Perform the ring communication and get the final buffer
    int* final_buf = circle(buf);

    // Print AFTER state (each rank in order)
    if (rank == 0) printf("\nAFTER\n");
    MPI_Barrier(MPI_COMM_WORLD);
    for (int r = 0; r < nprocs; r++) {
        if (rank == r) {
            for (int i = 0; i < local_sizes[r]; i++) {
                printf("rank %d: %d\n", rank, final_buf[i]);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Free all allocated memory
    if (final_buf != buf) {
        free(final_buf);
    }
    free(buf);
    free(local_sizes);
    free(displacements);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
