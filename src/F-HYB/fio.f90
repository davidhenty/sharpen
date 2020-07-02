!  Routine to find the size of a PGM image file. Note that this only
!  works if there are exactly two lines before the dimensions.

subroutine pgmsize(filename, nx, ny)

  implicit none

  character*(*) :: filename
  integer :: nx, ny

  integer, parameter :: iounit = 12

  open(unit=iounit, file=filename)

  read(iounit,*)
  read(iounit,*)
  read(iounit,*) nx, ny

  close(unit=iounit)

end subroutine pgmsize


!  Routine to read in a PGM image file as integers. Like pgmsize, this
!  only works if there are exactly two lines before the dimensions.

subroutine pgmread(filename, pixmap, nxmax, nymax, nx, ny)

  implicit none

  character*(*) :: filename
  integer :: nxmax, nymax, nx, ny
  integer, dimension(nxmax, nymax) :: pixmap

  integer i, j

  integer, parameter :: iounit = 12

  open(unit=iounit, file=filename)

  read(iounit,*)
  read(iounit,*)
  read(iounit,*) nx, ny

  if (nx .gt. nxmax .or. ny .gt. nymax) then
     write(*,*) 'Not enough space in pgmread'
     nx = 0
     ny = 0
     return
  end if

  read(iounit,*)

  read(iounit,*) ((pixmap(i,ny-j+1), i=1,nx), j=1,ny)

  close(unit=iounit)

end subroutine pgmread


!  Routine to write a PGM image file from a 2D floating point array
!  x(nx,ny). Uses unit 10 for IO.

subroutine pgmwrite(filename, x, nx, ny)

  implicit none

  character*(*) :: filename
  integer :: nx, ny

  double precision, dimension(nx, ny) :: x

  double precision, dimension(nx, ny) :: tmp
  integer,          dimension(nx, ny) :: grey

  double precision :: tmin, tmax
  double precision, parameter :: thresh = 255.0

  integer, parameter :: iounit = 10

  integer :: i, j

  tmp(:,:) = x(:,:)

  !  Find the max and min absolute values of the array

  tmin = minval(abs(tmp(:,:)))
  tmax = maxval(abs(tmp(:,:)))

  !  Scale the values appropriately so they lie between 0 and thresh

  if (tmin .lt. 0 .or. tmax .gt. thresh) then

     !    write(*,*) 'Scaling!'
     tmp(:,:) = int((thresh*((abs(tmp(:,:)-tmin))/(tmax-tmin))) + 0.5)

  else

     tmp(:,:) = int(abs(tmp(:,:)) + 0.5)

  end if

  grey(:,:) = tmp(:,:)

  !  Increase the contrast by boosting the lower values?
  !
  !  grey(:,:) = thresh * sqrt(tmp(:,:)/thresh)

  open(unit=iounit, file=filename)

  write(iounit,fmt='(''P2''/''# Written by pgmwrite'')')
  write(iounit,fmt='(i5, i5)') nx, ny
  write(iounit,*) int(thresh)
  write(iounit,fmt='(16(i3,'' ''))') ((grey(i,ny-j+1), i=1,nx), j=1,ny)  

  close(unit=iounit)

end subroutine pgmwrite
