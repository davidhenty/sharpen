!  Subroutine to sharpen an image by convolving with a filter function. 
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
!  David Henty, EPCC, August 2009
!  OpenMP version: Elena Breitmoser, EPCC, August 2009
!  Minor modifications: Arno Proeme, EPCC, March 2013
!

subroutine dosharpen(filename, nx, ny)

  use omp_lib

  implicit none

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

  integer :: i, j, k, l, nx, ny, xpix, ypix, pixcount, nthreads, threadid

  integer, dimension(nx, ny) :: fuzzy                             ! Will store the fuzzy input image when it is first read in from file
  double precision, dimension(1-d:nx+d, 1-d:ny+d) :: fuzzyPadded  ! Will store the fuzzy input image plus additional border padding
  double precision, dimension(1-d:nx+d, 1-d:ny+d) :: convolution  ! Will store the convolution of the filter with the fuzzy image
  double precision, dimension(nx, ny) :: sharp                    ! Will store the sharpened image obtained by adding the convolution to the fuzzy image

  fuzzy(:,:) = 0
  convolution(:,:) = 0.0D0

  write(*,fmt='(" Using a filter of size", I3," x ",I3," pixels")') 2*d+1, 2*d+1
  write(*,*) 'Reading input image file: ', filename
  flush(6)
  call pgmread(filename, fuzzy, nx, ny, xpix, ypix)  ! Loads fuzzy image from file into fuzzy image array
  write(*,*) '... done'
  write(*,*)
  flush(6)

  if (xpix == 0 .or. ypix == 0 .or. nx /= xpix .or. ny /= ypix) then

     write(*,*) 'Error reading ', filename
     flush(6)
     stop

  end if

  ! Initialise image arrays
  fuzzyPadded(:,:) = 0.0
  fuzzyPadded(1:nx, 1:ny) = fuzzy(1:nx, 1:ny)
  sharp(:,:) = 0.0

  write(*,*) 'Starting calculation ...'
  flush(6)

!$omp parallel
  call printlocation()
!$omp end parallel   
  
  tstart  = omp_get_wtime()

  ! Start of parallel region where filter is applied to fuzzy image
  !$omp parallel private(i, j, pixcount, k, l, threadid)

  nthreads  = omp_get_num_threads()
  threadid = omp_get_thread_num()

  pixcount = 0

  do i = 1, nx
     do j = 1, ny
        pixcount = pixcount + 1

        ! Computation of convolution allocated to threads using simple cyclic distribution
        ! i.e. consecutively numbered threads take turns computing convolution for consecutive pixels
        if (mod(pixcount-1, nthreads) == threadid) then

           do k = -d, d
              do l = -d, d
                 
                 convolution(i,j) = convolution(i,j) + filter(d, k, l)*fuzzyPadded(i+k, j+l)

              end do
           end do
        
        end if

     end do
  end do

  !$omp end parallel 
  !End parallel region and convolution computation

  tstop = omp_get_wtime()
  time  = tstop - tstart

  write(*,*) '.. finished'
  write(*,*)

  ! Add rescaled convolution to fuzzy image to obtain sharp image
  sharp(1:nx, 1:ny) = fuzzy(1:nx, 1:ny) - scale/norm * convolution(1:nx, 1:ny)

  filename = 'sharpened.pgm'
  write(*,*) 'Writing output image file: ', filename
  write(*,*)
  flush(6)

  !  Only save the core of the sharpened image to remove edge effects
  call pgmwrite(filename, sharp(d+1:nx-d,d+1:ny-d), nx-2*d, ny-2*d)
  write(*,*) '... done'
  write(*,*)
  write(*,fmt='(" Calculation time was ", f6.3, " seconds")') time
  flush(6)

end subroutine dosharpen
