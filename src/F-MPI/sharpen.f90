!  Program to sharpen an image by convolving with a filter function. The
!  filter is a combination of a Gaussian (to remove noise) and a Laplacian
!  (to detect the edges). Input and output is via Portable Grey Map (PGM)
!  files - note that the input file must have a specific header format.
!
!  In this version of the program the image processing is parallelised using
!  MPI processes and replicated data. The master process reads in the fuzzy 
!  image and distributes it to the other processes using an MPI communication. 
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

  use mpi

  implicit none

  integer :: comm, size, rank, ierr
  integer :: nx, ny, xpix, ypix

  double precision :: tstart, tstop, time

  integer, parameter :: maxlen = 32
  character*(maxlen) :: filename

  comm = MPI_COMM_WORLD

  call MPI_Init(ierr)

  call MPI_Comm_size(comm, size, ierr)
  call MPI_Comm_rank(comm, rank, ierr)

  call MPI_Barrier(comm, ierr)

  tstart = MPI_Wtime()

  filename = 'fuzzy.pgm'

  if (rank == 0) then

     write(*,*)
     write(*,fmt='(" Image sharpening code running on", I3," process(es)")') size
     write(*,*) 'Input file is: ', filename
     flush(6)
     call pgmsize(filename, xpix, ypix)
     write(*,fmt='(" Image size is: ", I3," x ",I3," pixels")') xpix, ypix
     write(*,*)
     flush(6)
  end if

  call MPI_Bcast(xpix, 1, MPI_INTEGER, 0, comm, ierr)
  call MPI_Bcast(ypix, 1, MPI_INTEGER, 0, comm, ierr)

  nx = xpix
  ny = ypix

  call dosharpen(filename, nx, ny, comm)

  call MPI_Barrier(comm, ierr)

  tstop = MPI_Wtime()
  time  = tstop - tstart

  if (rank == 0) then
     write(*,fmt='(" Overall run time was ", f6.3, " seconds")') time
     write(*,*)
     flush(6)

  end if

  call MPI_Finalize(ierr)

end program sharpen
