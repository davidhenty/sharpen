/*  Function to sharpen an image by convolving with a filter function. The
 *  filter is a combination of a Gaussian (to remove noise) and a Laplacian
 *  (to detect the edges). Input and output is via Portable Grey Map (PGM)
 *  files - note that the input file must have a specific header format.
 *
 *  In this version of the program the image processing is
 *  parallelised using CUDA for a GPU device. The CPU reads in the
 *  fuzzy image and stores it in host memory which is then copied to
 *  device memory on the GPU. CUDA kernels are launched on the GPU and
 *  the convolution computation is distributed over all CUDA threads.
 *  Finally, the CPU copies the result back to host memory, adds the
 *  convolution result to the fuzzy image and writes the resulting
 *  sharp image to file.
 *  
 *  David Henty, EPCC, September 2009
 *  Arno Proeme, EPCC, March 2013 (minor modifications)
 *  Dominic Sloan-Murphy, EPCC, June 2014 (improved consistency with other versions)
 */

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include "utilities.h"
#include "sharpen.h"

__global__ void dosharpenpixel(int nx, int ny, int d,
                    double *convolution, double *fuzzyPadded);

void dosharpen(char *infile, int nx, int ny, int verbose)
{
  int d = 8;
  /* Sets the linear range of the sharpen filter as measured from any given pixel:
     only pixels within a (2d+1)*(2d+1) square centered on the pixel are used to 
     compute its new value. */

  double  norm = (2*d-1)*(2*d-1);
  double scale = 2.0;
  
  int xpix, ypix;
  
  int i, j;
  double tstart, tstop, time;
  
  int fuzzy[nx][ny];                   /* Will store the fuzzy input image when it is first read in from file                        */
  double fuzzyPadded[nx+2*d][ny+2*d];  /* Will store the fuzzy input image plus additional border padding                            */
  double convolution[nx][ny];          /* Will store the convolution of the filter with the fuzzy image                              */
  double sharp[nx][ny];                /* Will store the sharpened image obtained by adding the  convolution to the fuzzy image      */
  double sharpCropped[nx-2*d][ny-2*d]; /* Will store the sharpened image cropped to remove a border layer distorted by the algorithm */
  
  double *d_fuzzyp, *d_conv;

  char outfile[] = "sharpened.pgm";
  
  /* Initialise image arrays */
  for (i=0; i < nx; i++)
    {
      for (j=0; j < ny; j++)
        {
          fuzzy[i][j] = 0;
          sharp[i][j] = 0.0;
          convolution[i][j] = 0.0;
        }
    }

  if (verbose)
    {
      printf("Using a filter of size %d x %d\n", 2*d+1, 2*d+1);
      printf("\n");

      printf("Reading image file: %s\n", infile);
      fflush(stdout);
    }
       
    pgmread(infile, fuzzy, nx, ny, &xpix, &ypix);

  if (verbose)
    {
      printf("... done\n\n");
      fflush(stdout);
    }  

  if (xpix == 0 || ypix == 0 || nx != xpix || ny != ypix)
    {
      printf("Error reading %s\n", infile);
      fflush(stdout);
      exit(-1);
    }
  
  /* Initialise image array */
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

  // Allocate CUDA memory
  cudaMalloc((void **) &d_conv,   nx*ny*sizeof(double));
  cudaMalloc((void **) &d_fuzzyp, (nx+2*d)*(ny+2*d)*ny*sizeof(double));

  // Copy

  cudaMemcpy(d_conv, convolution, nx*ny*sizeof(double),
             cudaMemcpyHostToDevice);

  cudaMemcpy(d_fuzzyp, fuzzyPadded, (nx+2*d)*(ny+2*d)*sizeof(double),
             cudaMemcpyHostToDevice);


  dim3 nthread = {16, 16, 1}; // 256 in a 16x16 grid
  dim3 nblock  = {(nx+nthread.x-1)/nthread.x, (ny+nthread.y-1)/nthread.y, 1};

  if (verbose)
    {
      printf("thread grid = %d x %d\n", nthread.x, nthread.y);
      printf("block  grid = %d x %d\n", nblock.x, nblock.y);

      printf("\nStarting calculation ...\n");
    }

  tstart = wtime();

  /* Start of parallel region where filter is applied to fuzzy image */

  dosharpenpixel<<<nblock, nthread>>>(nx, ny, d, d_conv, d_fuzzyp);
  cudaDeviceSynchronize();
  
  /* End of parallel region and convolution computation */
  
  tstop = wtime();
  time = tstop - tstart;
  
  cudaMemcpy(convolution, d_conv, nx*ny*sizeof(double),
             cudaMemcpyDeviceToHost);

  if (verbose)
    {
      printf("... finished\n");
      printf("\n");
      fflush(stdout);
    }
  
  /* Add rescaled convolution to fuzzy image to obtain sharp image */
  for (i=0; i < nx; i++)
    {
      for (j=0; j < ny; j++)
        {
          sharp[i][j] = fuzzyPadded[i+d][j+d] - scale/norm * convolution[i][j];
        }
    }

  if (verbose)
    {
      printf("Writing output file: %s\n", outfile);
      printf("\n");
    }
  
  /* Only save the core of the sharpened image to remove edge effects */
  for (i=d ; i < nx-d; i++)
    {
      for (j=d; j < ny-d; j++)
        {
          sharpCropped[i-d][j-d] = sharp[i][j];
        }
    }
  
  pgmwrite(outfile, sharpCropped, nx-2*d, ny-2*d);

  if (verbose)
    {
      printf("... done\n");
      printf("\n");
      printf("Calculation time was %f seconds\n", time);
      fflush(stdout);
    }

  // Free memory
  cudaFree(d_conv);
  cudaFree(d_fuzzyp);
}
