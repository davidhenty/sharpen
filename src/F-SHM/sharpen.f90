!  Program to sharpen an image by convolving with a filter function. The
!  filter is a combination of a Gaussian (to remove noise) and a Laplacian
!  (to detect the edges). Input and output is via Portable Grey Map (PGM)
!  files - note that the input file must have a specific header format.
!
!  In this version of the program the image processing is parallelised using
!  OpenSHMEM PEs and replicated data. The master process reads in the fuzzy 
!  image and distributes it to the other processes using a broadcast.
!  The  convolution computation is distributed over all processes and the result 
!  communicated back to the master process. Finally the master process adds the 
!  convolution result to the fuzzy image and writes the resulting sharp image to 
!  file.
!
!  Actual calculation is done in a subroutine to allow for declaration
! of automatic arrays of correct size.
!
!  David Henty, EPCC, September 2009
!  Arno Proeme, EPCC, March 2013 (minor modifications)
!

program sharpen

  implicit none

  include 'shmem.fh'

! Put various data in symmetric storage so it can be used in shmem calls

  integer, save :: xpix, ypix
  integer, save :: pSync(SHMEM_BCAST_SYNC_SIZE)

  data PSYNC /SHMEM_BCAST_SYNC_SIZE*SHMEM_SYNC_VALUE/
 
  integer :: size, rank, ierr
  integer :: nx, ny

  double precision :: tstart, tstop, time
  double precision :: wtime
  
  integer, parameter :: maxlen = 32
  character*(maxlen) :: filename

  call shmem_init();

  rank = shmem_my_pe()
  size = shmem_n_pes()

  call shmem_barrier_all()

  tstart = wtime()

  filename = 'fuzzy.pgm'

  if (rank == 0) then

     write(*,*)
     write(*,fmt='(" Image sharpening code running on", I3," PE(s)")') size
     write(*,*) 'Input file is: ', filename
     flush(6)
     call pgmsize(filename, xpix, ypix)
     write(*,fmt='(" Image size is: ", I3," x ",I3," pixels")') xpix, ypix
     write(*,*)
     flush(6)
  end if

  call shmem_broadcast4(xpix, xpix, 1, 0, 0, 0, size, pSync)
  call shmem_barrier_all()
  call shmem_broadcast4(ypix, ypix, 1, 0, 0, 0, size, pSync)

  nx = xpix
  ny = ypix

  call dosharpen(filename, nx, ny)

  call shmem_barrier_all()

  tstop = wtime()

  time  = tstop - tstart

  if (rank == 0) then
     write(*,fmt='(" Overall run time was ", f6.3, " seconds")') time
     write(*,*)
     flush(6)

  end if

  call shmem_finalize()

end program sharpen
