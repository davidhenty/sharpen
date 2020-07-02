/*  Program to sharpen an image by convolving with a filter function. 
 *  The filter is a combination of a Gaussian (to remove noise) and a 
 *  Laplacian (to detect the edges). Input and output is via Portable 
 *  Grey Map (PGM) files - note that the input file must have a specific
 *  header format.
 *
 *  This is a serial version.
 *  
 *  Actual calculation is done in a subroutine to allow for declaration
 *  of automatic arrays of correct size.
 *  
 *  David Henty, EPCC, September 2009
 *  Arno Proeme, EPCC, March 2013 (minor modifications)
 */


#include <stdio.h>
#include <stdlib.h>
#include "sharpen.h"

int main(void)
{
  double tstart, tstop, time;
  
  char *filename;
  int xpix, ypix;
  
  filename = "fuzzy.pgm";
  
  printf("\n");
  printf("Image sharpening code running in serial\n");
  printf("\n");
  printf("Input file is: %s\n", filename);
  
  pgmsize(filename, &xpix, &ypix);
  
  printf("Image size is %d x %d\n", xpix, ypix);
  printf("\n");
  
  tstart  = wtime();
  
  dosharpen(filename, xpix, ypix);
  
  tstop = wtime();
  time  = tstop - tstart;
  
  printf("Overall run time was %f seconds\n", time);
}
