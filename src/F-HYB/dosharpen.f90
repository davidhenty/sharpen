!  Subroutine to sharpen an image by convolving with a filter function. The
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
!  David Henty, EPCC, September 2009
!  Arno Proeme, EPCC, March 2013 (minor modifications)
!  Dominic Sloan-Murphy, EPCC, November 2013 (more minor modifications)
!


subroutine dosharpen(filename, nx, ny, comm)

  use mpi
  use omp_lib

  implicit none

  integer :: comm, size, rank, ierr

  integer :: nthreads, threadid
  integer :: globalsize, globalid

  integer, parameter :: d = 8  
  ! Sets the linear range of the sharpen filter as measured from any given pixel:
  ! only pixels within a (2d+1)*(2d+1) square centered on the pixel are used to 
  ! compute its new value.

  double precision, parameter :: norm = dble((2*d-1)*(2*d-1))
  double precision, parameter :: scale = 2.0

  double precision :: tstart, tstop, time

  double precision, external ::  filter
  
  integer, parameter :: maxlen = 32
  character*(maxlen) :: filename
  integer, parameter :: iounit = 10

  integer :: i, j, k, l, nx, ny, xpix, ypix, pixcount

  integer, dimension(nx, ny) :: fuzzy                                    ! Will store the fuzzy input image when it is first read in from file
  double precision, dimension(1-d:nx+d, 1-d:ny+d) :: fuzzyPadded         ! Will store the fuzzy input image plus additional border padding
  double precision, dimension(nx, ny) :: convolutionPartial              ! Will store the convolution of the filter with parts of the fuzzy image computed by individual processes
  double precision, dimension(nx, ny) :: convolution                     ! Will store the convolution of the filter with the fuzzy image
  double precision, dimension(nx, ny) :: sharp                           ! Will store the sharpened image obtained by adding the convolution to the fuzzy image
  
  call MPI_Comm_size(comm, size, ierr)
  call MPI_Comm_rank(comm, rank, ierr)

  call MPI_Barrier(comm, ierr)

  fuzzy(:,:) = 0
  convolutionPartial(:,:) = 0.0D0

  if (rank == 0) then
     write(*,fmt='(" Using a filter of size", I3," x ",I3," pixels")') 2*d+1, 2*d+1
     write(*,*) 'Reading image file: ', filename
     flush(6)
     call pgmread(filename, fuzzy, nx, ny, xpix, ypix)
     write(*,*) '... done'
     write(*,*)
     flush(6)
  end if

  call MPI_Bcast(xpix, 1, MPI_INTEGER, 0, comm, ierr)
  call MPI_Bcast(ypix, 1, MPI_INTEGER, 0, comm, ierr)

  if (xpix == 0 .or. ypix == 0 .or. nx /= xpix .or. ny /= ypix) then

     if (rank == 0) write(*,*) 'Error reading ', filename
     flush(6)

     call MPI_Finalize(ierr)
     stop

  end if

  !  Broadcast the pixel image to all processes
  call MPI_Bcast(fuzzy, nx*ny, MPI_INTEGER, 0, comm, ierr)
  
  ! Initialise image arrays
  fuzzyPadded(:,:) = 0.0
  fuzzyPadded(1:nx, 1:ny) = fuzzy(1:nx, 1:ny)
  sharp(:,:) = 0.0
  
  if (rank == 0) write(*,*) 'Starting calculation ...'
  flush(6)
  
  call MPI_Barrier(comm, ierr)
  
!$omp parallel
  call printlocation()
!$omp end parallel   

  tstart = MPI_Wtime()
  
!$omp parallel default(none) &
!$omp shared (nx, ny, convolutionPartial, fuzzyPadded, rank, size) &
!$omp private(i, j, k, l, pixcount, nthreads, threadid, globalsize, globalid)

  nthreads = omp_get_num_threads()
  threadid = omp_get_thread_num()

  globalsize = size * nthreads
  globalid   = rank * nthreads + threadid

  pixcount = 0
  
  do i = 1, nx
     do j = 1, ny

        pixcount = pixcount + 1
        ! Computation of convolution allocated to threads using simple cyclic distribution
        ! i.e. consecutively numbered threads take turns computing convolution for consecutive pixels

        if (mod(pixcount-1, globalsize) == globalid) then   

           do k = -d, d
              do l = -d, d
                 
                 convolutionPartial(i,j) = convolutionPartial(i,j) + filter(d, k, l)*fuzzyPadded(i+k, j+l)
                 
              end do
           end do

        end if

     end do
  end do

!$omp end parallel

  call MPI_Barrier(comm, ierr)

  tstop = MPI_Wtime()
  time  = tstop - tstart

  if (rank == 0) then

     write(*,*) '.. finished'
     write(*,*)
     flush(6)

  end if

  ! Gather the partial convolution results computed by individual processes
  call MPI_Reduce(convolutionPartial, convolution, nx*ny, &
       MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)

  ! The master process applies the filter and writes the sharpened image to file 
  if (rank == 0) then

     ! Add rescaled convolution to fuzzy image to obtain sharp image
     sharp(1:nx, 1:ny) = fuzzy(1:nx, 1:ny) - scale/norm * convolution(1:nx, 1:ny)

     filename = 'sharpened.pgm'
     write(*,*) 'Writing output file: ', filename
     write(*,*)
     flush(6)

     !  Only save the core of the sharpened image to remove edge effects
     call pgmwrite(filename, sharp(d+1:nx-d,d+1:ny-d), nx-2*d, ny-2*d)
     write(*,*) '... done'
     write(*,*)
     flush(6)

  end if

  call MPI_Barrier(comm, ierr)

  if (rank == 0) then
     write(*,fmt='(" Calculation time was ", f6.3, " seconds")') time
     flush(6)
     
  end if

end subroutine dosharpen
