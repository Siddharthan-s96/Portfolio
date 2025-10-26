/*
 * partdiff - PDE solver for Gauß-Seidel and Jacobi methods
 * Copyright (C) 1997 Thomas Ludwig
 * Copyright (C) 1997 Thomas A. Zochler
 * Copyright (C) 1997 Andreas C. Schmidt
 * Copyright (C) 2007-2010 Julian M. Kunkel
 * Copyright (C) 2010-2021 Michael Kuhn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <sys/time.h>
#include <mpi.h>

/* ************* */
/* Some defines. */
/* ************* */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_INTERLINES    100000
#define MAX_ITERATION     200000
#define MAX_THREADS       1024
#define METH_GAUSS_SEIDEL 1
#define METH_JACOBI       2
#define FUNC_F0           1
#define FUNC_FPISIN       2
#define TERM_PREC         1
#define TERM_ITER         2

struct calculation_arguments
{
	uint64_t N;            /* number of spaces between lines (lines=N+1) */
	uint64_t num_matrices; /* number of matrices */
	double   h;            /* length of a space between two lines */
	double*  M;            /* two matrices with real values */
    int      my_rows;
};

struct calculation_results
{
	uint64_t m;
	uint64_t stat_iteration; /* number of current iteration */
	double   stat_precision; /* actual precision of all slaves in iteration */
};

struct options
{
	uint64_t number;         /* Number of threads */
	uint64_t method;         /* Gauss Seidel or Jacobi method of iteration */
	uint64_t interlines;     /* matrix size = interlines*8+9 */
	uint64_t termination;    /* termination condition */
	uint64_t term_iteration; /* terminate if iteration number reached */
	double   term_precision; /* terminate if precision reached */
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time; /* time when program started */
struct timeval comp_time;  /* time when calculation completed */
static void calculateMpiJacobi(struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options);

static void
usage(char* name)
{
	printf("Usage: %s [num] [method] [lines] [func] [term] [prec/iter]\n", name);
	printf("\n");
	printf("  - num:       number of threads (1 .. %d)\n", MAX_THREADS);
	printf("  - method:    calculation method (1 .. 2)\n");
	printf("                 %1d: Gauß-Seidel\n", METH_GAUSS_SEIDEL);
	printf("                 %1d: Jacobi\n", METH_JACOBI);
	printf("  - lines:     number of interlines (0 .. %d)\n", MAX_INTERLINES);
	printf("                 matrixsize = (interlines * 8) + 9\n");
	printf("  - term:      termination condition ( 1.. 2)\n");
	printf("                 %1d: sufficient precision\n", TERM_PREC);
	printf("                 %1d: number of iterations\n", TERM_ITER);
	printf("  - prec/iter: depending on term:\n");
	printf("                 precision:  1e-4 .. 1e-20\n");
	printf("                 iterations:    1 .. %d\n", MAX_ITERATION);
	printf("\n");
	printf("Example: %s 1 2 100 2 100 \n", name);
}

static void
askParams(struct options* options, int argc, char** argv)
{
	int ret;

	if (argc < 6 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-?") == 0)
	{
		usage(argv[0]);
		exit(0);
	}

	ret = sscanf(argv[1], "%" SCNu64, &(options->number));

	if (ret != 1 || !(options->number >= 1 && options->number <= MAX_THREADS))
	{
		usage(argv[0]);
		exit(1);
	}

	ret = sscanf(argv[2], "%" SCNu64, &(options->method));

	if (ret != 1 || !(options->method == METH_GAUSS_SEIDEL || options->method == METH_JACOBI))
	{
		usage(argv[0]);
		exit(1);
	}

	ret = sscanf(argv[3], "%" SCNu64, &(options->interlines));

	if (ret != 1 || !(options->interlines <= MAX_INTERLINES))
	{
		usage(argv[0]);
		exit(1);
	}


	ret = sscanf(argv[4], "%" SCNu64, &(options->termination));

	if (ret != 1 || !(options->termination == TERM_PREC || options->termination == TERM_ITER))
	{
		usage(argv[0]);
		exit(1);
	}

	if (options->termination == TERM_PREC)
	{
		ret = sscanf(argv[5], "%lf", &(options->term_precision));

		options->term_iteration = MAX_ITERATION;

		if (ret != 1 || !(options->term_precision >= 1e-20 && options->term_precision <= 1e-4))
		{
			usage(argv[0]);
			exit(1);
		}
	}
	else
	{
		ret = sscanf(argv[5], "%" SCNu64, &(options->term_iteration));

		options->term_precision = 0;

		if (ret != 1 || !(options->term_iteration >= 1 && options->term_iteration <= MAX_ITERATION))
		{
			usage(argv[0]);
			exit(1);
		}
	}
}

/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static void
initVariables(struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N            = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h            = 1.0 / arguments->N;
    arguments->my_rows      = 0;

	results->m              = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static void
freeMatrices(struct calculation_arguments* arguments)
{
	free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static void*
allocateMemory(size_t size)
{
	void* p;

	if ((p = malloc(size)) == NULL)
	{
		printf("Memory error! (%" PRIu64 " Bytes requested)\n", size);
		exit(1);
	}

	return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static void
allocateMatrices(struct calculation_arguments* arguments)
{
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    uint64_t const N = arguments->N;
    
    if (size == 1) {
        /* Serial case - same as original */
        arguments->my_rows = N - 1;  /* All interior rows (1 to N-1) */
        arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
    } else {
        /* MPI case - distribute rows among processes */
        int total_interior_rows = N - 1;  /* Rows 1 to N-1 are interior rows */
        int base_rows = total_interior_rows / size;
        int extra_rows = total_interior_rows % size;
        
        /* Calculate how many rows this process gets */
        int my_rows = base_rows + (rank < extra_rows ? 1 : 0);
        arguments->my_rows = my_rows;
        
        /* Allocate local matrix: my_rows + 2 ghost rows (one above, one below) */
        int local_rows = my_rows + 2;
        arguments->M = allocateMemory(arguments->num_matrices * local_rows * (N + 1) * sizeof(double));
    }
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static void
initMatrices(struct calculation_arguments* arguments, struct options const* options)
{
    (void)options;

    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    uint64_t g, j, i;
    uint64_t const N = arguments->N;

    if (size == 1) {
        /* Serial case - same as original */
        typedef double(*matrix)[N + 1][N + 1];
        matrix Matrix = (matrix)arguments->M;

        for (g = 0; g < arguments->num_matrices; g++) {
            for (i = 0; i <= N; i++) {
                for (j = 0; j <= N; j++) {
                    Matrix[g][i][j] = 0.0;
                }
            }
        }
    } else {
        /* MPI case - initialize local matrix */
        uint64_t local_rows = arguments->my_rows + 2;
        typedef double(*matrix)[local_rows][N + 1];
        matrix Matrix = (matrix)arguments->M;

        for (g = 0; g < arguments->num_matrices; g++) {
            for (i = 0; i < local_rows; i++) {
                for (j = 0; j <= N; j++) {
                    Matrix[g][i][j] = 0.0;
                }
            }
        }
    }
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static void
calculate(struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size == 1) {
        /* SERIAL VERSION - Keep original algorithm exactly */
        int    i, j;        
        int    m1, m2;      
        double star;        
        double residuum;    
        double maxresiduum; 

        int const    N = arguments->N;
        double const h = arguments->h;

        double pih    = 0.0;
        double fpisin = 0.0;

        int term_iteration = options->term_iteration;

        typedef double(*matrix)[N + 1][N + 1];
        matrix Matrix = (matrix)arguments->M;

        /* initialize m1 and m2 depending on algorithm */
        if (options->method == METH_JACOBI) {
            m1 = 0;
            m2 = 1;
        } else {
            m1 = 0;
            m2 = 0;
        }

        pih    = M_PI * h;
        fpisin = 0.25 * (2 * M_PI * M_PI) * h * h;

        while (term_iteration > 0) {
            maxresiduum = 0;

            /* over all rows */
            for (i = 1; i < N; i++) {
                double fpisin_i = 0.0;
                fpisin_i = fpisin * sin(pih * (double)i);

                /* over all columns */
                for (j = 1; j < N; j++) {
                    star = 0.25 * (Matrix[m2][i - 1][j] + Matrix[m2][i][j - 1] + Matrix[m2][i][j + 1] + Matrix[m2][i + 1][j]);
                    star += fpisin_i * sin(pih * (double)j);

                    if (options->termination == TERM_PREC || term_iteration == 1) {
                        residuum    = Matrix[m2][i][j] - star;
                        residuum    = fabs(residuum);
                        maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
                    }

                    Matrix[m1][i][j] = star;
                }
            }

            results->stat_iteration++;
            results->stat_precision = maxresiduum;

            /* exchange m1 and m2 */
            i  = m1;
            m1 = m2;
            m2 = i;

            /* check for stopping calculation depending on termination method */
            if (options->termination == TERM_PREC) {
                if (maxresiduum < options->term_precision) {
                    term_iteration = 0;
                }
            } else if (options->termination == TERM_ITER) {
                term_iteration--;
            }
        }

        results->m = m2;
        
    } else {
    /* MPI VERSION - Only Jacobi method supported */
        calculateMpiJacobi(arguments, results, options);
    }  
}

/* ************************************************************************ */
/* calculate: solves the equation for MPI                                          */
/* ************************************************************************ */

static void
calculateMpiJacobi(struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int    i, j;
    int    m1 = 0, m2 = 1;
    double star, residuum, maxresiduum, global_maxresiduum;

    int const    N = arguments->N;
    double const h = arguments->h;
    int const    my_rows = arguments->my_rows;
    int local_rows = my_rows + 2;

    /* Calculate global row range for this process */
    int total_interior_rows = N - 1;
    int base_rows = total_interior_rows / size;
    int extra_rows = total_interior_rows % size;
    int global_start_row = 1 + rank * base_rows + (rank < extra_rows ? rank : extra_rows);

    typedef double(*matrix)[local_rows][N + 1];
    matrix Matrix = (matrix)arguments->M;

    double pih    = M_PI * h;
    double fpisin = 0.25 * (2 * M_PI * M_PI) * h * h;
    int term_iteration = options->term_iteration;

    while (term_iteration > 0) {
        maxresiduum = 0;

        /* Exchange ghost rows with neighbors */
        MPI_Request reqs[4];
        int nreqs = 0;

        /* Send to next process, receive from next process */
        if (rank < size - 1) {
            MPI_Isend(Matrix[m2][my_rows], N + 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &reqs[nreqs++]);
            MPI_Irecv(Matrix[m2][my_rows + 1], N + 1, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &reqs[nreqs++]);
        }
        
        /* Send to previous process, receive from previous process */
        if (rank > 0) {
            MPI_Isend(Matrix[m2][1], N + 1, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &reqs[nreqs++]);
            MPI_Irecv(Matrix[m2][0], N + 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &reqs[nreqs++]);
        }
        
        if (nreqs > 0)
            MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);

        /* Calculate new values for interior points */
        for (i = 1; i <= my_rows; i++) {
            int global_i = global_start_row + (i - 1);
            double fpisin_i = fpisin * sin(pih * (double)global_i);
            
            for (j = 1; j < N; j++) {
                star = 0.25 * (Matrix[m2][i - 1][j] + Matrix[m2][i][j - 1] + 
                              Matrix[m2][i][j + 1] + Matrix[m2][i + 1][j]);
                star += fpisin_i * sin(pih * (double)j);

                if (options->termination == TERM_PREC || term_iteration == 1) {
                    residuum = Matrix[m2][i][j] - star;
                    residuum = fabs(residuum);
                    maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
                }
                Matrix[m1][i][j] = star;
            }
        }

        /* Find global maximum residuum */
        MPI_Allreduce(&maxresiduum, &global_maxresiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        results->stat_iteration++;
        results->stat_precision = global_maxresiduum;
        
        /* Swap matrices */
        int tmp = m1; m1 = m2; m2 = tmp;

        /* Check termination condition */
        if (options->termination == TERM_PREC) {
            if (global_maxresiduum < options->term_precision)
                term_iteration = 0;
        } else if (options->termination == TERM_ITER) {
            term_iteration--;
        }
    }
    results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static void
displayStatistics(struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
	int    N    = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Calculation time:   %f s\n", time);
	printf("Memory usage:       %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Calculation method: ");

	if (options->method == METH_GAUSS_SEIDEL)
	{
		printf("Gauß-Seidel");
	}
	else if (options->method == METH_JACOBI)
	{
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:         %" PRIu64 "\n", options->interlines);

	printf("Termination:        ");

	if (options->termination == TERM_PREC)
	{
		printf("Required precision");
	}
	else if (options->termination == TERM_ITER)
	{
		printf("Number of iterations");
	}

	printf("\n");
	printf("Number iterations:  %" PRIu64 "\n", results->stat_iteration);
	printf("Error precision:    %e\n", results->stat_precision);
	printf("\n");
}

/****************************************************************************/
/** Explanation of the displayMatrix function:                             **/
/**                                                                        **/
/** The function displayMatrix outputs a Matrix                            **/
/** in a humanly readable way.                                             **/
/**                                                                        **/
/** This is achieved by only printing parts of the matrix.                 **/
/** From the matrix the first and last lines/columns and seven in between  **/
/** rows/cols are printed out.                                             **/
/****************************************************************************/
static void
displayMatrix(struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
	int x, y;

	int const interlines = options->interlines;
	int const N          = arguments->N;

	typedef double(*matrix)[N + 1][N + 1];

	matrix Matrix = (matrix)arguments->M;

	printf("Matrix:\n");

	for (y = 0; y < 9; y++)
	{
		for (x = 0; x < 9; x++)
		{
			printf("%7.4f", Matrix[results->m][y * (interlines + 1)][x * (interlines + 1)]);
		}

		printf("\n");
	}

	fflush(stdout);
}

/****************************************************************************/
/** displayMatrixMpi function:                                             **/
/****************************************************************************/

static void
displayMatrixMpi(struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
{
  int const elements = 8 * options->interlines + 9;

  int x, y;

  int local_rows = arguments->my_rows + 2;
  typedef double(*matrix)[local_rows][arguments->N + 1];
  matrix Matrix = (matrix)arguments->M;
  int m = results->m;

  MPI_Status status;

  // The global first line belongs to rank 0
  if (rank == 0) {
    from--;
  }

  // The global last line belongs to rank (size - 1)
  if (rank == size - 1) {
    to++;
  }

  if (rank == 0) {
    printf("Matrix:\n");
  }

  for (y = 0; y < 9; y++)
  {
    int line = y * (options->interlines + 1);

    if (rank == 0)
    {
      // Check if the line belongs to rank 0
      if (line < from || line > to)
      {
        // The tag is used to receive lines in the correct order
        // Matrix[m][0] is overwritten because the data is no longer needed
        MPI_Recv(Matrix[m][0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
      }
    }
    else
    {
      if (line >= from && line <= to)
      {
        // Send line to rank 0 if it belongs to this process
        // (line - from + 1) is used to caluclate the local line 
        MPI_Send(Matrix[m][line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
      }
    }

    if (rank == 0)
    {
      for (x = 0; x < 9; x++)
      {
        int col = x * (options->interlines + 1);

        if (line >= from && line <= to)
        {
          // Diese Zeile gehört zu Rang 0
          printf("%7.4f", Matrix[m][line][col]);
        }
        else
        {
          // This line belong to another rank and was received above
          printf("%7.4f", Matrix[m][0][col]);
        }
      }

      printf("\n");
    }
  }

  fflush(stdout);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main(int argc, char** argv)
{
    struct options               options;
    struct calculation_arguments arguments;
    struct calculation_results   results;
    int rank = 0, size = 1;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* Parse parameters (only rank 0) and broadcast to all processes */
    if (rank == 0) {
        askParams(&options, argc, argv);
        if (size > 1 && options.method != METH_JACOBI) {
            printf("Error: MPI version only supports Jacobi method!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    
    /* Broadcast options to all processes */
    MPI_Bcast(&options, sizeof(struct options), MPI_BYTE, 0, MPI_COMM_WORLD);

    initVariables(&arguments, &results, &options);
    allocateMatrices(&arguments);
    initMatrices(&arguments, &options);

    /* Synchronize all processes before timing */
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        gettimeofday(&start_time, NULL);

    calculate(&arguments, &results, &options);

    /* Synchronize all processes after calculation */
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        gettimeofday(&comp_time, NULL);

    /* Display results (only rank 0) */
    if (rank == 0) 
        displayStatistics(&arguments, &results, &options);
        
    /* Display matrix */
    /* Display matrix */
    if (size == 1) {
        displayMatrix(&arguments, &results, &options);
    } else {
        /* Calculate from/to for this process */
        int total_interior_rows = arguments.N - 1;
        int base_rows = total_interior_rows / size;
        int extra_rows = total_interior_rows % size;
        int from = 1 + rank * base_rows + (rank < extra_rows ? rank : extra_rows);
        int to = from + arguments.my_rows - 1;
        displayMatrixMpi(&arguments, &results, &options, rank, size, from, to);
    }

    freeMatrices(&arguments);
    MPI_Finalize();
    return 0;
}
