/*  Program to sharpen an image by convolving with a filter function. 
 *  The filter is a combination of a Gaussian (to remove noise) and a 
 *  Laplacian (to detect the edges). Input and output is via Portable 
 *  Grey Map (PGM) files - note that the input file must have a specific
 *  header format.
 *
 *  In this version of the program the image processing is parallelised 
 *  using OpenMP threads. The master thread reads in the fuzzy image and
 *  stores it in shared memory. Further threads are launched, the convolution
 *  computation is distributed over all threads (including the master thread),
 *  and the result stored in shared memory. Finally the master thread adds the 
 *  convolution result to the fuzzy image and writes the resulting sharp image
 *  to file.
 *  
 *  Actual calculation is done in a subroutine to allow for declaration
 *  of automatic arrays of correct size.
 *  
 *  David Henty, EPCC, September 2009
 *  Arno Proeme, EPCC, March 2013 (minor modifications)
 */


#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "sharpen.h"

int main(void)
{
  int nthread;
  double tstart, tstop, time;
  
  char *filename;
  int xpix, ypix;
  
  nthread = omp_get_max_threads();
  
  filename = "fuzzy.pgm";
  
  printf("\n");
  printf("Image sharpening code running on %d thread(s)\n", nthread);
  printf("\n");
  printf("Input file is: %s\n", filename);
  
  pgmsize(filename, &xpix, &ypix);
  
  printf("Image size is %d x %d\n", xpix, ypix);
  printf("\n");
  
  tstart  = omp_get_wtime();
  
  dosharpen(filename, xpix, ypix);
  
  
  tstop = omp_get_wtime();
  time  = tstop - tstart;
  
  printf("Overall run time was %f seconds\n", time);
}
