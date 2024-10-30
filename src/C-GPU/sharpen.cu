/*  Program to sharpen an image by convolving with a filter function. 
 *  The filter is a combination of a Gaussian (to remove noise) and a 
 *  Laplacian (to detect the edges). Input and output is via Portable 
 *  Grey Map (PGM) files - note that the input file must have a specific
 *  header format.
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
 *  Actual calculation is done in a subroutine to allow for declaration
 *  of automatic arrays of correct size.
 *
 *  There appears to be a significant startup cost the first time the
 *  GPU is used so the calculation is performed twice with all timings
 *  taken on the second pass.
 *  
 *  David Henty, EPCC, September 2009
 *  Arno Proeme, EPCC, March 2013 (minor modifications)
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "sharpen.h"

int main(void)
{
  double tstart, tstop, time;
  
  char filename[] = "fuzzy.pgm";
  int xpix, ypix;
  
  
  printf("\n");
  printf("Image sharpening code running on GPU with CUDA\n");
  printf("\n");
  printf("Input file is: %s\n", filename);
  
  pgmsize(filename, &xpix, &ypix);
  
  printf("Image size is %d x %d\n", xpix, ypix);
  printf("\n");
  
  // Initial dummy call to warm up the GPU

  dosharpen(filename, xpix, ypix, false);

  // Now measure the time

  tstart  = wtime();
  
  dosharpen(filename, xpix, ypix, true);
  
  tstop = wtime();
  time  = tstop - tstart;
  
  printf("Overall run time was %f seconds\n", time);
}
