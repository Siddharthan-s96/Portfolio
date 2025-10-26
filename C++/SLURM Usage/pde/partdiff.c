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
#define TERM_PREC         1
#define TERM_ITER         2

struct calculation_arguments
{
	uint64_t N;            /* number of spaces between lines (lines=N+1) */
	uint64_t num_matrices; /* number of matrices */
	double   h;            /* length of a space between two lines */
	double*  M;            /* two matrices with real values */
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
	uint64_t const N = arguments->N;

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static void
initMatrices(struct calculation_arguments* arguments)
{
	uint64_t g, i, j; /* local variables for loops */

	uint64_t const N = arguments->N;

	typedef double(*matrix)[N + 1][N + 1];

	matrix Matrix = (matrix)arguments->M;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= N; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

}

/* ************************************************************************ */
/* calculate_func: calculates the interference function                     */
/* ************************************************************************ */
/* ---------------------------------------------------------------------- */
/* calculate_func: optimized math for the source term                    */
/* ---------------------------------------------------------------------- */
double
calculate_func(struct calculation_arguments const* arguments, int i, int j)
{
	double const h = arguments->h;

    return 0.25 * ((2 * M_PI * M_PI) * h * h * sin(M_PI * h * (double)i) * sin(M_PI * h * (double)j));
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static void
calculate(struct calculation_arguments const* arguments,
                      struct calculation_results       *results,
                      struct options            const *options)
{
    int    const N      = arguments->N;
    double const h      = arguments->h;
    double const pi     = M_PI;
    double const pi_h   = pi * h;
    double const coeff  = (pi * pi * h * h) / 2.0;  /* (2·π²·h²)/4 */
    double const inv4   = 0.25;

    /* 1) Build a tiny sin‐lookup: sin(pi_h * k) for k=0..N */
    double *sin_table = malloc((N+1) * sizeof *sin_table);
    if (!sin_table) { perror("malloc"); exit(EXIT_FAILURE); }
    for (int k = 0; k <= N; ++k)
        sin_table[k] = sin(pi_h * (double)k);

    /* 2) Cast the 3D matrix back to double[M+1][M+1] slices */
    typedef double (*matrix)[N+1][N+1];
    matrix Matrix = (matrix)arguments->M;

    int m1 = 0, m2 = (options->method == METH_JACOBI ? 1 : 0);
    int term_iter = options->term_iteration;                   

    while (term_iter > 0)
    {
        double maxres = 0.0;

        if (options->method == METH_JACOBI)
        {
            /* ——— Jacobi (two‐matrix ping‐pong) ——— */
            double (* __restrict__ src)[N+1] = Matrix[m2];
            double (* __restrict__ dst)[N+1] = Matrix[m1];

            for (int i = 1; i < N; ++i)
            {
                double const s_i = sin_table[i];
                for (int j = 1; j < N; ++j)
                {
                    double sum4 =
                        src[i-1][j] +
                        src[i  ][j-1] +
                        src[i  ][j+1] +
                        src[i+1][j];

                    double star = inv4 * sum4
                                  + coeff * s_i * sin_table[j];

                    double d = fabs(src[i][j] - star);
                    if (d > maxres) maxres = d;

                    dst[i][j] = star;
                }
            }

            /* swap buffers */
            m1 ^= 1;  
            m2 ^= 1;
        }
        else
        {
            /* ——— Red–Black Gauss–Seidel (in‐place) ———
             * Two “colors” so that within each inner loop
             * there are no j↦j–1 dependencies.
             */
            for (int color = 0; color < 2; ++color)
            {
                for (int i = 1; i < N; ++i)
                {
                    double const s_i = sin_table[i];
                    /* pointers to three distinct rows—no aliasing */
                    double * __restrict__ row    = Matrix[0][i];
                    double * __restrict__ row_p  = Matrix[0][i-1];
                    double * __restrict__ row_n  = Matrix[0][i+1];

                    /* pick starting j so that (i+j)%2 == color */
                    int j0 = 1 + ((i + color) & 1);
                    for (int j = j0; j < N; j += 2)
                    {
                        double sum4 =
                            row_p[j] +
                            row[j-1] +
                            row[j+1] +
                            row_n[j];

                        double star = inv4 * sum4
                                      + coeff * s_i * sin_table[j];

                        double d = fabs(row[j] - star);
                        if (d > maxres) maxres = d;

                        row[j] = star;
                    }
                }
            }
            /* in‐place, so no buffer swap */
        }

        results->stat_iteration++;
        results->stat_precision = maxres;

        /* termination check */
        if (options->termination == TERM_PREC)
        {
            if (maxres < options->term_precision)
                break;
        }
        else  /* TERM_ITER */
        {
            --term_iter;
        }
    }

    results->m = (options->method == METH_JACOBI ? m2 : 0);
    free(sin_table);
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
/** Beschreibung der Funktion displayMatrix:                               **/
/**                                                                        **/
/** Die Funktion displayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
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

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main(int argc, char** argv)
{
	struct options               options;
	struct calculation_arguments arguments;
	struct calculation_results   results;

	askParams(&options, argc, argv);

	initVariables(&arguments, &results, &options);

	allocateMatrices(&arguments);
	initMatrices(&arguments);

	gettimeofday(&start_time, NULL);
	calculate(&arguments, &results, &options);
	gettimeofday(&comp_time, NULL);

	displayStatistics(&arguments, &results, &options);
	displayMatrix(&arguments, &results, &options);

	freeMatrices(&arguments);

	return 0;
}
