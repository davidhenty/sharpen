!  Subroutine to sharpen an image by convolving with a filter function. The
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
!  David Henty, EPCC, September 2009
!  Arno Proeme, EPCC, March 2013 (minor modifications)
!  Dominic Sloan-Murphy, EPCC, November 2013 (more minor modifications)
!


subroutine dosharpen(filename, nx, ny)

  implicit none

  include 'shmem.fh'

  double precision :: wtime

  integer :: size, rank, ierr, pWrksize, numput, ipe

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

  integer :: i, j, k, l, nx, ny, pixcount

  double precision, dimension(1-d:nx+d, 1-d:ny+d) :: fuzzyPadded         ! Will store the fuzzy input image plus additional border padding
  double precision, dimension(nx, ny) :: sharp                           ! Will store the sharpened image obtained by adding the convolution to the fuzzy image

! Data that needs to be in shared memory

  integer, save :: xpix, ypix

  integer, save :: pSync(SHMEM_BCAST_SYNC_SIZE)
  data PSYNC /SHMEM_BCAST_SYNC_SIZE*SHMEM_SYNC_VALUE/

  double precision :: pWrk(1)

  integer, dimension(nx, ny) :: fuzzy                                    ! Will store the fuzzy input image when it is first read in from file
  double precision, dimension(nx, ny) :: convolution                     ! Will store the convolution of the filter with the fuzzy image

! Use cray pointers to allow these arrays to be dynamically allocated

  pointer(fuzzaddr, fuzzy)
  pointer(convaddr, convolution)
  pointer(pWrkaddr, pWrk)
  
! Call symmetric allocation routine and abort on error (last param = 0)
! Remember that convolution contains doubles so need twice as many words
! Same goes for for pWrk

  call shpalloc(fuzzaddr,   nx*ny, ierr, 0)
  call shpalloc(convaddr, 2*nx*ny, ierr, 0)

  pWrksize = max((nx*ny)/2 + 1, SHMEM_REDUCE_MIN_WRKDATA_SIZE)
  call shpalloc(pWrkaddr, 2*nx*ny, ierr, 0)

  rank = shmem_my_pe()
  size = shmem_n_pes()

  call shmem_barrier_all()

  fuzzy(:,:) = 0
  convolution(:,:) = 0.0D0

  if (rank == 0) then
     write(*,fmt='(" Using a filter of size", I3," x ",I3," pixels")') 2*d+1, 2*d+1
     write(*,*) 'Reading image file: ', filename
     flush(6)
     call pgmread(filename, fuzzy, nx, ny, xpix, ypix)
     write(*,*) '... done'
     write(*,*)
     flush(6)
  end if

  call shmem_broadcast4(xpix, xpix, 1, 0, 0, 0, size, pSync)
  call shmem_barrier_all()
  call shmem_broadcast4(ypix, ypix, 1, 0, 0, 0, size, pSync)
  call shmem_barrier_all()

  if (xpix == 0 .or. ypix == 0 .or. nx /= xpix .or. ny /= ypix) then

     if (rank == 0) write(*,*) 'Error reading ', filename
     flush(6)

     call shmem_finalize()
     stop

  end if

  !  Broadcast the pixel image to all processes
  !  Seems to be a problem at the moment using this routine from Fortran so
  !  just do it by hand.
  !
  !   call shmem_broadcast4(fuzzy, fuzzy, nx*ny, 0, 0, 0, size, pSync)

  if (rank == 0) then

     do ipe = 1, size-1

        call shmem_integer_put(fuzzy, fuzzy, nx*ny, ipe)

     end do

  end if

  call shmem_barrier_all()

  ! Initialise image arrays
  fuzzyPadded(:,:) = 0.0
  fuzzyPadded(1:nx, 1:ny) = fuzzy(1:nx, 1:ny)
  sharp(:,:) = 0.0
  
  call printlocation()

  if (rank == 0) write(*,*) 'Starting calculation ...'
  flush(6)
  
  call shmem_barrier_all()
  
  tstart = wtime()
  
  pixcount = 0
  
  do j = 1, ny
     do i = 1, nx

        pixcount = pixcount + 1
        ! Computation of convolution allocated to threads using simple cyclic distribution
        ! i.e. consecutively numbered threads take turns computing convolution for consecutive pixels

        if (mod(pixcount-1, size) == rank) then   

           do k = -d, d
              do l = -d, d
                 
                 convolution(i,j) = convolution(i,j) + filter(d, k, l)*fuzzyPadded(i+k, j+l)
                 
              end do
           end do

        end if

     end do
  end do

  call shmem_barrier_all()

  tstop = wtime()
  time  = tstop - tstart

  if (rank == 0) then

     write(*,*) '.. finished'
     write(*,*)
     flush(6)

  end if

  if (pixcount /= nx*ny) then

     write(*,*) 'ERROR on rank ', rank, ': pixcount = ', pixcount, &
          ', nx*ny = ', nx*ny
     flush(6)

     call shmem_finalize()
     stop

  end if

  ! Gather the partial convolution results computed by individual processes
  ! Could use global reductions as in MPI (see immediately below) but more
  ! illustrative to use strided puts of doubles. This relies on the
  ! pixels being scanned in "memory" order. Also need to be careful computimg
  ! how many to send.
  !
  !  call shmem_real8_sum_to_all(convolution, convolution, nx*ny, &
  !                           0, 0, size, pWrk, pSync)
  !

  call shmem_barrier_all()

  if (rank /= 0) then

     numput = (nx*ny + (size-rank-1))/size;

     call shmem_double_iput(convolution(rank+1, 1), convolution(rank+1, 1), &
                            size, size, numput, 0)
  end if

  call shmem_barrier_all()

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

  call shmem_barrier_all()

  if (rank == 0) then
     write(*,fmt='(" Calculation time was ", f6.3, " seconds")') time
     flush(6)
     
  end if

end subroutine dosharpen
