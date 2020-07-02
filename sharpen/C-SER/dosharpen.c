/*  Function to sharpen an image by convolving with a filter function. The
 *  filter is a combination of a Gaussian (to remove noise) and a Laplacian
 *  (to detect the edges). Input and output is via Portable Grey Map (PGM)
 *  files - note that the input file must have a specific header format.
 *
 *  In this traditional, serial version of the program the image processing
 *  is performed by one single-threaded process running on a single core.
 *  The fuzzy image is read in, the convolution computed and then added to 
 *  the fuzzy image. The resulting sharp image is written to file.
 *  
 *  David Henty, EPCC, September 2009
 *  Arno Proeme, EPCC, March 2013 (minor modifications)
 *  Dominic Sloan-Murphy, EPCC, June 2014 (improved consistency with other versions)
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "sharpen.h"
#include "utilities.h"

void dosharpen(char *infile, int nx, int ny)
{
  int d = 8;  
  /* Sets the linear range of the sharpen filter as measured from any given pixel:
     only pixels within a (2d+1)*(2d+1) square centered on the pixel are used to 
     compute its new value. */

  double  norm = (2*d-1)*(2*d-1);
  double scale = 2.0;
  
  int xpix, ypix, pixcount;
  
  int i, j, k, l;
  double tstart, tstop, time;
  
  int **fuzzy = int2Dmalloc(nx, ny);                   /* Will store the fuzzy input image when it is first read in from file */
  double **fuzzyPadded = double2Dmalloc(nx+2*d, ny+2*d);  /* Will store the fuzzy input image plus additional border padding */
  double **convolutionPartial = double2Dmalloc(nx, ny);   /* Will store the convolution of the filter with parts of the fuzzy image computed by individual processes */
  double **convolution = double2Dmalloc(nx, ny);          /* Will store the convolution of the filter with the full fuzzy image */
  double **sharp = double2Dmalloc(nx, ny);                /* Will store the sharpened image obtained by adding rescaled convolution to the fuzzy image */
  double **sharpCropped = double2Dmalloc(nx-2*d, ny-2*d); /* Will store the sharpened image cropped to remove a border layer distorted by the algorithm */
  
  char *outfile = "sharpened.pgm";

  /* Initialise image arrays */
  for (i=0; i < nx; i++)
    {
      for (j=0; j < ny; j++)
        {
          fuzzy[i][j] = 0;
          sharp[i][j] = 0.0;
        }
    }
  
  printf("Using a filter of size %d x %d\n", 2*d+1, 2*d+1);
  printf("\n");

  printf("Reading image file: %s\n", infile);
  fflush(stdout);
       
  pgmread(infile, &fuzzy[0][0], nx, ny, &xpix, &ypix);
  printf("... done\n\n");
  fflush(stdout);
  
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
  
  printf("Starting calculation ...\n");
  
  fflush(stdout);
  printlocation();
  fflush(stdout);
  
  tstart = wtime();
  
  pixcount = 0;
  
  for (i=0; i < nx; i++)
    {
      for (j=0; j < ny; j++)
        {
          for (k=-d; k <= d; k++)
            {
              for (l= -d; l <= d; l++)
                {
                  convolution[i][j] = convolution[i][j] + filter(d,k,l)*fuzzyPadded[i+d+k][j+d+l];
                }
            }
          pixcount += 1;
        }
    }
  
  tstop = wtime();
  time = tstop - tstart;
  
  printf("... finished\n");
  printf("\n");
  fflush(stdout);
  
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
  
  pgmwrite(outfile, &sharpCropped[0][0], nx-2*d, ny-2*d);
  
  printf("... done\n");
  printf("\n");
  printf("Calculation time was %f seconds\n", time);
  fflush(stdout);

  free(fuzzy);
  free(fuzzyPadded);
  free(convolutionPartial);
  free(convolution);
  free(sharp);
  free(sharpCropped);
}

int **int2Dmalloc(int nx, int ny)
{
  int i;
  int **idata;

  idata = (int **) malloc(nx*sizeof(int *) + nx*ny*sizeof(int));

  idata[0] = (int *) (idata + nx);

  for(i=1; i < nx; i++)
    {
      idata[i] = idata[i-1] + ny;
    }

  return idata;
}

double **double2Dmalloc(int nx, int ny)
{
  int i;
  double **ddata;

  ddata = (double **) malloc(nx*sizeof(double *) + nx*ny*sizeof(double));

  ddata[0] = (double *) (ddata + nx);

  for(i=1; i < nx; i++)
    {
      ddata[i] = ddata[i-1] + ny;
    }

  return ddata;
}
