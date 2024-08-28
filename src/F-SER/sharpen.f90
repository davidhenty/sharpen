!  Program to sharpen an image by convolving with a filter function. 
!  The filter is a combination of a Gaussian (to remove noise) and a 
!  Laplacian (to detect the edges). Input and output is via Portable 
!  Grey Map (PGM) files - note that the input file must have a specific
!  header format.
!
!  This is a serial version.
! 
!  Actual calculation is done in a subroutine to allow for declaration
!  of automatic arrays of correct size.
!
!  David Henty, EPCC, August 2009
!  OpenMP version: Elena Breitmoser, EPCC, August 2009
!  Minor modifications: Arno Proeme, EPCC, March 2013
!

program sharpen

  implicit none

  integer :: nx, ny, xpix, ypix

  double precision :: wtime, tstart, tstop, time

  integer, parameter :: maxlen = 32
  character*(maxlen) :: filename

  filename = 'fuzzy.pgm'

  write(*,*)
  write(*, fmt='(" Image sharpening code running in serial")')
  write(*,*)
  write(*,*) 'Input file is: ', filename

  call pgmsize(filename, xpix, ypix)
  
  write(*,fmt='(" Image size is: ", I3," x ",I3," pixels")') xpix, ypix
  write(*,*)

  nx = xpix
  ny = ypix

  tstart  = wtime()

  call dosharpen(filename, nx, ny)

  tstop = wtime()
  time  = tstop - tstart

  write(*,fmt='(" Overall run time was ", f6.3, " seconds")') time
  write(*,*)

end program sharpen
