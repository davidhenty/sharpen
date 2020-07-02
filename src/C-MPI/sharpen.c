/*  Program to sharpen an image by convolving with a filter function. The
 *  filter is a combination of a Gaussian (to remove noise) and a Laplacian
 *  (to detect the edges). Input and output is via Portable Grey Map (PGM)
 *  files - note that the input file must have a specific header format.
 *
 *  In this version of the program the image processing is parallelised using
 *  MPI processes and replicated data. The master process reads in the fuzzy 
 *  image and distributes it to the other processes using an MPI communication. 
 *  The  convolution computation is distributed over all processes and the result 
 *  communicated back to the master process. Finally the master process adds the 
 *  convolution result to the fuzzy image and writes the resulting sharp image to 
 *  file.
 *
 *  Actual calculation is done in a subroutine to allow for declaration
 *  of automatic arrays of correct size.
 *
 *  David Henty, EPCC, September 2009
 *  Arno Proeme, EPCC, March 2013 (minor modifications)
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "sharpen.h"

int main(void)
{
  MPI_Comm comm;
  int rank, size;
  double tstart, tstop, time;

  char *filename;
  int xpix, ypix;

  comm = MPI_COMM_WORLD;

  MPI_Init(NULL, NULL);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  MPI_Barrier(comm);

  tstart = MPI_Wtime();

  filename = "fuzzy.pgm";

  if (rank == 0)
    {
      printf("\n");
      printf("Image sharpening code running on %d processor(s)\n", size);
      printf("\n");
      printf("Input file is: %s\n", filename);

      pgmsize(filename, &xpix, &ypix);

      printf("Image size is %d x %d\n", xpix, ypix);
      printf("\n");

      fflush(stdout);
    }
      
  MPI_Bcast(&xpix, 1, MPI_INT, 0, comm);
  MPI_Bcast(&ypix, 1, MPI_INT, 0, comm);

  dosharpen(filename, xpix, ypix, comm);

  MPI_Barrier(comm);

  tstop = MPI_Wtime();
  time  = tstop - tstart;

  if (rank == 0)
    {
      printf("Overall run time was %f seconds\n", time);
    }

  MPI_Finalize();
}
