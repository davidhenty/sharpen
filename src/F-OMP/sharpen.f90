!  Program to sharpen an image by convolving with a filter function. 
!  The filter is a combination of a Gaussian (to remove noise) and a 
!  Laplacian (to detect the edges). Input and output is via Portable 
!  Grey Map (PGM) files - note that the input file must have a specific
!  header format.
!
!  In this version of the program the image processing is parallelised 
!  using OpenMP threads. The master thread reads in the fuzzy image and
!  stores it in shared memory. Further threads are launched, the convolution
!  computation is distributed over all threads (including the master thread),
!  and the result stored in shared memory. Finally the master thread adds the 
!  convolution result to the fuzzy image and writes the resulting sharp image
!  to file.
! 
!  Actual calculation is done in a subroutine to allow for declaration
!  of automatic arrays of correct size.
!
!  David Henty, EPCC, August 2009
!  OpenMP version: Elena Breitmoser, EPCC, August 2009
!  Minor modifications: Arno Proeme, EPCC, March 2013
!

program sharpen

  use omp_lib

  implicit none

  integer :: nx, ny, xpix, ypix

  double precision :: tstart, tstop, time
  integer :: nthreads

  integer, parameter :: maxlen = 32
  character*(maxlen) :: filename

  filename = 'fuzzy.pgm'

  nthreads = omp_get_max_threads()

  write(*,*)
  write(*, fmt='(" Image sharpening code running on", I3, " thread(s)")') nthreads
  write(*,*)
  write(*,*) 'Input file is: ', filename

  call pgmsize(filename, xpix, ypix)
  
  write(*,fmt='(" Image size is: ", I3," x ",I3," pixels")') xpix, ypix
  write(*,*)

  nx = xpix
  ny = ypix

  tstart  = omp_get_wtime()

  call dosharpen(filename, nx, ny)

  tstop = omp_get_wtime()
  time  = tstop - tstart

  write(*,fmt='(" Overall run time was ", f6.3, " seconds")') time
  write(*,*)

end program sharpen
