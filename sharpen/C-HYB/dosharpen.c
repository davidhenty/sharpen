/*  Function to sharpen an image by convolving with a filter function. The
 *  filter is a combination of a Gaussian (to remove noise) and a Laplacian
 *  (to detect the edges). Input and output is via Portable Grey Map (PGM)
 *  files - note that the input file must have a specific header format.
 *
 *  In this version of the program the image processing is
 *  parallelised using a hybrid of MPI processes and OpenMP threads,
 *  with replicated data. The master process reads in the fuzzy image
 *  and distributes it to the other processes using an MPI
 *  communication.
 *
 *  The convolution computation is distributed over all processes and
 *  threads, and the result communicated back to the master
 *  process. Finally the master process adds the convolution result to
 *  the fuzzy image and writes the resulting sharp image to file.
 *
 *  David Henty, EPCC, September 2009
 *  Arno Proeme, EPCC, March 2013 (minor modifications)
 *  Dominic Sloan-Murphy, EPCC, November 2013 (more minor modifications)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include "utilities.h"
#include "sharpen.h"

void dosharpen(char *infile, int nx, int ny, MPI_Comm comm)
{
  int        d = 8; 
  /* Sets the linear range of the sharpen filter as measured from any given pixel:
     only pixels within a (2d+1)*(2d+1) square centered on the pixel are used to 
     compute its new value. */
  
  double  norm = (2*d-1)*(2*d-1);  
  double scale = 2.0;
  
  int rank, size;
  int xpix, ypix, pixcount;

  int nthreads, threadid;
  int globalsize, globalid;

  int i, j, k, l;
  double tstart, tstop, time;

  int fuzzy[nx][ny];                   /* Will store the fuzzy input image when it is first read in from file                                     */
  double fuzzyPadded[nx+2*d][ny+2*d];  /* Will store the fuzzy input image plus additional border padding                                         */
  double convolutionPartial[nx][ny];   /* Will store the convolution of the filter with parts of the fuzzy image computed by individual processes */
  double convolution[nx][ny];          /* Will store the convolution of the filter with the full fuzzy image                                      */
  double sharp[nx][ny];                /* Will store the sharpened image obtained by adding rescaled convolution to the fuzzy image               */
  double sharpCropped[nx-2*d][ny-2*d]; /* Will store the sharpened image cropped to remove a border layer distorted by the algorithm              */
  
  char *outfile = "sharpened.pgm";

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

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
       
      pgmread(infile, fuzzy, nx, ny, &xpix, &ypix);
      printf("... done\n\n");
      fflush(stdout);
    }

  MPI_Bcast(&xpix, 1, MPI_INT, 0, comm);
  MPI_Bcast(&ypix, 1, MPI_INT, 0, comm);

  if (xpix == 0 || ypix == 0 || nx != xpix || ny != ypix)
    {
      if (rank == 0) printf("Error reading %s\n", infile);
      fflush(stdout);

      MPI_Finalize();
      exit(-1);
    }

  /* Broadcast the pixel image to all processes */

  MPI_Bcast(fuzzy, nx*ny, MPI_INT, 0, comm);

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

  if (rank == 0) printf("Starting calculation ...\n");

  MPI_Barrier(comm);
  
  /* Print out current core and node location. */
#pragma omp parallel
{
  printlocation();
}    

  tstart = MPI_Wtime();

#pragma omp parallel default(none) \
  shared(nx, ny, d, convolutionPartial, fuzzyPadded, rank, size) \
  private(i, j, k, l, pixcount, nthreads, threadid, globalsize, globalid)
{

  nthreads = omp_get_num_threads();
  threadid = omp_get_thread_num();

  globalsize = size * nthreads;
  globalid = rank*nthreads + threadid;

  pixcount = 0;

  for (i=0; i < nx; i++)
    {
      for (j=0; j < ny; j++)
        {
          /* Computation of convolution allocated to processes using simple cyclic distribution
             i.e. consecutively ranked processes take turns computing convolution for consecutive pixels */
          if (pixcount%globalsize  == globalid)
            {
              for (k=-d; k <= d; k++)
                {
                  for (l= -d; l <= d; l++)
                    {
                      convolutionPartial[i][j] = convolutionPartial[i][j] + filter(d,k,l)*fuzzyPadded[i+d+k][j+d+l];
                    }
                }
            }
          pixcount += 1;
        }
    }
}  

  MPI_Barrier(comm);

  tstop = MPI_Wtime();
  time = tstop - tstart;

  if (rank == 0)
    {
      printf("... finished\n");
      printf("\n");
      fflush(stdout);
    }

  /* Gather the partial convolution results computed by individual processes */
  MPI_Reduce(convolutionPartial, convolution, nx*ny, MPI_DOUBLE, MPI_SUM, 0, comm);
  
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
