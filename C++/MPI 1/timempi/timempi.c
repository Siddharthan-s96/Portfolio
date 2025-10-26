/**
 * @file timempi.c
 * @brief MPI program that displays hostname and timestamp from multiple processes
 *
 * This program demonstrates MPI communication where:
 * - Processes with ranks 1 to n create hostname:timestamp strings
 * - All strings are sent to rank 0 for ordered output
 * - The minimum microsecond value across all processes is computed and displayed
 * - All processes synchronize before terminating
 *
 */ 

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <mpi.h>

#define MAX_STRING_LENGTH 128   /**< Maximum length for all strings */
#define HOSTNAME_LENGTH    64   /**< Maximum length for hostname   */

/**
 * @brief Main function for the timempi MPI program
 *
 * This function implements the following workflow:
 * 1. Initialize MPI and get rank/size information
 * 2. Get current timestamp and hostname for all processes
 * 3. Rank 0 receives and prints strings from other processes in order
 * 4. Find minimum microsecond value across all processes using MPI_Reduce
 * 5. Synchronize all processes using MPI_Barrier
 * 6. Each process prints termination message
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line arguments
 * @return Exit status (0 for success)
 */
int main(int argc, char *argv[]) {
    int rank, size;                           /**< MPI rank and communicator size */
    char hostname[HOSTNAME_LENGTH];           /**< Buffer for hostname */
    char timestamp_buf[32];                   /**< Buffer for formatted timestamp */
    char output_string[MAX_STRING_LENGTH];    /**< Buffer for complete output string */
    struct tm* tm_info;                       /**< Structure for time information */
    struct timeval tv;                        /**< Structure for time with microseconds */
    long microseconds;                        /**< Microseconds from current time */
    long min_microseconds;                    /**< Minimum microseconds across all processes */

    /* Initialize MPI environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);    /* Get process rank */
    MPI_Comm_size(MPI_COMM_WORLD, &size);    /* Get total number of processes */

    /* Get current time and hostname for all processes
     * (All processes need this for the microsecond reduction operation) */
    gettimeofday(&tv, NULL);                  /* Get time with microsecond precision */
    gethostname(hostname, HOSTNAME_LENGTH);   /* Get system hostname */
    tm_info = localtime(&tv.tv_sec);          /* Convert to local time */
    strftime(timestamp_buf, sizeof(timestamp_buf), "%Y-%m-%d %H:%M:%S", tm_info); /* Format timestamp */
    microseconds = tv.tv_usec;                /* Extract microseconds */

    if (rank == 0) {
        /* Process with rank 0: Receive strings from other processes and handle output */
        char received_string[MAX_STRING_LENGTH]; /**< Buffer for received strings */

        /* Receive and print strings from ranks 1 to n-1 in sequential order
         * This ensures ordered output by process rank */
        for (int i = 1; i < size; i++) {
            MPI_Recv(received_string, MAX_STRING_LENGTH, MPI_CHAR, i, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("%s\n", received_string);
        }

        /* Perform reduction to find minimum microseconds across all processes
         * Rank 0 acts as the root process and receives the result */
        MPI_Reduce(&microseconds, &min_microseconds, 1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
        printf("%ld\n", min_microseconds);

    } else {
        /* Processes with ranks 1 to n: Create and send hostname:timestamp string */

        /* Format the output string: hostname: YYYY-MM-DD HH:MM:SS.microseconds */
        snprintf(output_string, MAX_STRING_LENGTH, "%s: %s.%06ld", hostname, timestamp_buf, microseconds);

        /* Send the formatted string to rank 0 */
        MPI_Send(output_string, strlen(output_string) + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);

        /* Participate in the reduction operation (but don't receive the result) */
        MPI_Reduce(&microseconds, &min_microseconds, 1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
    }

    /* Synchronize all processes before termination messages
     * This ensures no process terminates before all output is complete */
    MPI_Barrier(MPI_COMM_WORLD);

    /* Each process outputs its termination message */
    printf("Rank %d terminating now!\n", rank);

    /* Clean up MPI environment */
    MPI_Finalize();

    return 0;
}
