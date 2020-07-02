/*  Subroutine to sharpen an image by convolving with a filter function. The
 *  filter is a combination of a Gaussian (to remove noise) and a Laplacian
 *  (to detect the edges). Input and output is via Portable Grey Map (PGM)
 *  files - note that the input file must have a specific header format.
 *
 *  In this version of the program the image processing is parallelised using
 *  OpenSHMEM PEs and replicated data. The master process reads in the fuzzy 
 *  image and distributes it to the other processes using a broadcast. 
 *  The  convolution computation is distributed over all processes and the result 
 *  communicated back to the master process. Finally the master process adds the 
 *  convolution result to the fuzzy image and writes the resulting sharp image to 
 *  file.
 *
 *  David Henty, EPCC, September 2009
 *  Arno Proeme, EPCC, March 2013 (minor modifications)
 *  Dominic Sloan-Murphy, EPCC, November 2013 (more minor modifications)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <shmem.h>

#include "utilities.h"
#include "sharpen.h"

/*
 * Data that must be in symmetric storage - easiest to simply delcare
 * small variables like this in the data segment.
 */

extern long pSync[_SHMEM_BCAST_SYNC_SIZE];
extern int xpix, ypix;

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

void dosharpen(char *infile, int nx, int ny)
{
  int        d = 8; 
  /* Sets the linear range of the sharpen filter as measured from any given pixel:
     only pixels within a (2d+1)*(2d+1) square centered on the pixel are used to 
     compute its new value. */
  
  double  norm = (2*d-1)*(2*d-1);  
  double scale = 2.0;
  
  int rank, size;
  int pixcount, numput;

  int i, j, k, l;
  double tstart, tstop, time;

  double fuzzyPadded[nx+2*d][ny+2*d];  /* Will store the fuzzy input image plus additional border padding                                         */
  double sharp[nx][ny];                /* Will store the sharpened image obtained by adding rescaled convolution to the fuzzy image               */
  double sharpCropped[nx-2*d][ny-2*d]; /* Will store the sharpened image cropped to remove a border layer distorted by the algorithm              */
  
  /*
   * The two arrays below must be dynamically allocated in symmetric
   * memory in order to use OpenSHMEM. Original declarations are left
   * here as comments for clarity.
   */

  //  int fuzzy[nx][ny];                   /* Will store the fuzzy input image when it is first read in from file                                     */
  // double convolution[nx][ny];          /* Will store the convolution of the filter with the full fuzzy image                                      */

  int **fuzzy;
  double **convolution;
  double *pWrk;
  int pWrksize;

  char *outfile = "sharpened.pgm";

  rank = shmem_my_pe();
  size = shmem_n_pes();

  /* Allocate arrays in symmetric memory */

  fuzzy = symallocint2d(nx, ny);
  convolution = symallocdouble2d(nx, ny);

  /* Allocated work array needed for reduction */

  pWrksize = MAX((nx*ny)/2 + 1, _SHMEM_REDUCE_MIN_WRKDATA_SIZE);

  pWrk = (double *) shmalloc(pWrksize*sizeof(double));

  /* Initialise image arrays */
  for (i=0; i < nx; i++)
    {
      for (j=0; j < ny; j++)
        {
          fuzzy[i][j] = 0;        
          sharp[i][j] = 0.0;
        }
    }

  if (rank == 0)
    {
      printf("Using a filter of size %d x %d\n", 2*d+1, 2*d+1);
      printf("\n");

      printf("Reading image file: %s\n", infile);
      fflush(stdout);
       
      pgmread(infile, &fuzzy[0][0], nx, ny, &xpix, &ypix);
      printf("... done\n\n");
      fflush(stdout);
    }

  shmem_broadcast32(&xpix, &xpix, 1, 0, 0, 0, size, pSync);
  shmem_barrier_all(); // Needed to ensure we can reuse pSync

  shmem_broadcast32(&ypix, &ypix, 1, 0, 0, 0, size, pSync);
  shmem_barrier_all(); // Needed to ensure we can reuse pSync

  if (xpix == 0 || ypix == 0 || nx != xpix || ny != ypix)
    {
      if (rank == 0) printf("Error reading %s\n", infile);
      fflush(stdout);

      shmem_finalize();
      exit(-1);
    }

  /* Broadcast the pixel image to all processes */

  shmem_broadcast32(&fuzzy[0][0], &fuzzy[0][0], nx*ny, 0, 0, 0, size, pSync);

  for (i=0; i < nx+2*d; i++)
    {
      for (j=0; j < ny+2*d; j++)
        {
          fuzzyPadded[i][j] = 0.0;
        }
    }

  /* Transfer fuzzy image into padded array */
  for (i=0; i < nx; i++)
    {
      for (j=0; j < ny; j++)
        {
          fuzzyPadded[i+d][j+d] = fuzzy[i][j];
        }
    }

  shmem_barrier_all();
  
  /* Print out current core and node location. */

  fflush(stdout);
  printlocation();
  fflush(stdout);
    
  if (rank == 0) printf("Starting calculation ...\n\n");

  shmem_barrier_all();

  tstart = wtime();

  pixcount = 0;

  for (i=0; i < nx; i++)
    {
      for (j=0; j < ny; j++)
        {
          /* Computation of convolution allocated to processes using simple cyclic distribution
             i.e. consecutively ranked processes take turns computing convolution for consecutive pixels */
          if (pixcount%size  == rank)
            {
              for (k=-d; k <= d; k++)
                {
                  for (l= -d; l <= d; l++)
                    {
                      convolution[i][j] = convolution[i][j] + filter(d,k,l)*fuzzyPadded[i+d+k][j+d+l];
                    }
                }
            }
          pixcount += 1;
        }
    }
  
  shmem_barrier_all();

  tstop = wtime();
  time = tstop - tstart;

  if (rank == 0)
    {
      printf("... finished\n");
      printf("\n");
      fflush(stdout);
    }

  /*
   * Gather the partial convolution results computed by individual processes.
   * Could use global reductions as in MPI (see immediately below) but more
   * illustrative to use strided puts of doubles. This relies on the
   * pixels being scanned in "memory" order. Also need to be careful computimg
   * how many to send.
   *
   *   shmem_double_sum_to_all(&convolution[0][0], &convolution[0][0],
   *                          nx*ny, 0, 0, size, pWrk, pSync);
   */

  shmem_barrier_all();

  if (rank != 0)
    {
      numput = (nx*ny + (size-rank-1))/size;

      shmem_double_iput(&convolution[0][rank], &convolution[0][rank],
                        size, size, numput, 0);
    }

  shmem_barrier_all();

  /* The master process applies the filter and writes the sharpened image to file */
  if (rank == 0)
    {
      /* Add rescaled convolution to fuzzy image to obtain sharp image */
      for (i=0 ; i < nx; i++)
        {
          for (j=0; j < ny; j++)
            {
              sharp[i][j] = fuzzyPadded[i+d][j+d] - scale/norm * convolution[i][j];
            }
        }

      printf("Writing output file: %s\n", outfile);
      printf("\n");

      /* Only save the core of the sharpened image to remove edge effects */
      for (i=d ; i < nx-d; i++)
        {
          for (j=d; j < ny-d; j++)
            {
              sharpCropped[i-d][j-d] = sharp[i][j];
            }
        }
      
      pgmwrite(outfile, sharpCropped, nx-2*d, ny-2*d);

      printf("... done\n");
      printf("\n");
      printf("Calculation time was %f seconds\n", time);
      fflush(stdout);
    }
}


int **symallocint2d(int nx, int ny)
{
  int i;
  int **idata;

  idata = (int **) malloc(sizeof(int *)*nx);

  idata[0] = (int *) shmalloc(nx*ny*sizeof(int));

  for(i=1; i < nx; i++)
    {
      idata[i] = idata[i-1]+ny;
    }

  return idata;
}

double **symallocdouble2d(int nx, int ny)
{
  int i;
  double **idata;

  idata = (double **) malloc(sizeof(double *)*nx);

  idata[0] = (double *) shmalloc(nx*ny*sizeof(double));

  for(i=1; i < nx; i++)
    {
      idata[i] = idata[i-1]+ny;
    }

  return idata;
}
