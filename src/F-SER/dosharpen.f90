!  Subroutine to sharpen an image by convolving with a filter function. 
!  The filter is a combination of a Gaussian (to remove noise) and a 
!  Laplacian (to detect the edges). Input and output is via Portable 
!  Grey Map (PGM) files - note that the input file must have a specific
!  header format.
!
!  In this traditional, serial version of the program the image processing
!  is performed by one single-threaded process running on a single core.
!  The fuzzy image is read in, the convolution computed and then added to 
!  the fuzzy image. The resulting sharp image is written to file.
!
!  David Henty, EPCC, August 2009
!  OpenMP version: Elena Breitmoser, EPCC, August 2009
!  Minor modifications: Arno Proeme, EPCC, March 2013
!  Improved consistency with other versions: Dominic Sloan-Murphy, EPCC, June 2014
!

subroutine dosharpen(filename, nx, ny)

  implicit none

  integer, parameter :: d = 8 
  ! Sets the linear range of the sharpen filter as measured from any given pixel:
  ! only pixels within a (2d+1)*(2d+1) square centered on the pixel are used to 
  ! compute its new value.

  double precision, parameter :: norm = dble((2*d-1)*(2*d-1))
  double precision, parameter :: scale = 2.0

  double precision :: wtime, tstart, tstop, time

  double precision ::  filter

  integer, parameter :: maxlen = 32
  character*(maxlen) :: filename
  integer, parameter :: iounit = 10

  integer :: i, j, k, l, nx, ny, xpix, ypix, pixcount

  integer, allocatable, dimension(:,:) :: fuzzy                       ! Will store the fuzzy input image when it is first read in from file
  double precision, allocatable, dimension(:,:) :: fuzzyPadded        ! Will store the fuzzy input image plus additional border padding
  double precision, allocatable, dimension(:,:) :: convolution        ! Will store the convolution of the filter with the fuzzy image
  double precision, allocatable, dimension(:,:) :: sharp              ! Will store the sharpened image obtained by adding the convolution to the fuzzy image

  allocate(fuzzy(nx,ny))
  allocate(fuzzyPadded(1-d:nx+d, 1-d:ny+d))
  allocate(convolution(nx,ny))
  allocate(sharp(nx,ny))

  fuzzy(:,:) = 0
  convolution(:,:) = 0.0D0

  write(*,fmt='(" Using a filter of size", I3," x ",I3," pixels")') 2*d+1, 2*d+1
  write(*,*) 'Reading image file: ', filename
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

  call printlocation()  
  
  tstart = wtime()

  pixcount = 0

  do i = 1, nx
     do j = 1, ny
        pixcount = pixcount + 1

           do k = -d, d
              do l = -d, d
                 
                 convolution(i,j) = convolution(i,j) + filter(d, k, l)*fuzzyPadded(i+k, j+l)

              end do
           end do

     end do
  end do

  tstop = wtime()
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

  deallocate(fuzzy)
  deallocate(fuzzyPadded)
  deallocate(convolution)
  deallocate(sharp)

end subroutine dosharpen
