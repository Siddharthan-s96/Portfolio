/*
 * rank and size are the MPI rank and the size of the communicator
 * from and to stand for the global(!) range of rows managed by this process
 *
 * Example with 9 matrix rows and 4 processes
 *  - rank 0 is responsible for lines 1-2, rank 1 for 3-4, rank 2 for 5-6 and rank 3 for 7
 *  - Rows 0 amd 8 are not included because they are not being caluclated by the program
 *  - Every process stores two border lines in its local matrix
 *  - E.g.: rank 2 has 4 rows 0-3 but only calculates 1-2 because 0 and 4 are both border lines from other processes; 
 *    rank 2 therefore is responsible for the global matrix lines 5-6
 */

static void
displayMatrixMpi(struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
{
  int const elements = 8 * options->interlines + 9;

  int x, y;

  typedef double(*matrix)[to - from + 3][arguments->N + 1];
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
