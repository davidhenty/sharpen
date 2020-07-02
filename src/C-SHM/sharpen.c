/*  Program to sharpen an image by convolving with a filter function. The
 *  filter is a combination of a Gaussian (to remove noise) and a Laplacian
 *  (to detect the edges). Input and output is via Portable Grey Map (PGM)
 *  files - note that the input file must have a specific header format.
 *
 *  In this version of the program the image processing is parallelised using
 *  OpenSHME PEs and replicated data. The master process reads in the fuzzy 
 *  image and distributes it to the other processes using a broadcast. 
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
#include <shmem.h>
#include "sharpen.h"

/*
 * Data that must be in symmetric storage - easiest to simply delcare
 * small variables like this in the data segment.
 */


long pSync[_SHMEM_BCAST_SYNC_SIZE];
int xpix, ypix;

int main(void)
{
  int rank, size;
  double tstart, tstop, time;

  int i;

  char *filename;

  for (i = 0; i < _SHMEM_BCAST_SYNC_SIZE; i++)
    {
      pSync[i] = _SHMEM_SYNC_VALUE;
    }

  shmem_init();

  rank = shmem_my_pe();
  size = shmem_n_pes();
  shmem_barrier_all();

  tstart = wtime();

  filename = "fuzzy.pgm";

  if (rank == 0)
    {
      printf("\n");
      printf("Image sharpening code running on %d PE(s)\n", size);
      printf("\n");
      printf("Input file is: %s\n", filename);

      pgmsize(filename, &xpix, &ypix);

      printf("Image size is %d x %d\n", xpix, ypix);
      printf("\n");

      fflush(stdout);
    }
      

  shmem_broadcast32(&xpix, &xpix, 1, 0, 0, 0, size, pSync);
  shmem_barrier_all();
  shmem_broadcast32(&ypix, &ypix, 1, 0, 0, 0, size, pSync);

  dosharpen(filename, xpix, ypix);

  shmem_barrier_all();

  tstop = wtime();
  time  = tstop - tstart;

  if (rank == 0)
    {
      printf("Overall run time was %f seconds\n", time);
    }

  shmem_finalize();
}

