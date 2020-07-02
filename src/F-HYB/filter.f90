double precision function filter(d, i, j)

  implicit none

  integer :: d, i, j

  double precision :: rd4sq, rsq, sigmad4sq, sigmasq, x, y, delta

  integer, parameter :: d4      =     4

  double precision,    parameter :: sigmad4 =   1.4
  double precision,    parameter :: filter0 = -40.0

  rd4sq = dble(d4*d4)
  rsq   = dble(d*d)

  sigmad4sq = sigmad4**2
  sigmasq   = sigmad4sq * (rsq/rd4sq)

  x = dble(i)
  y = dble(j)

  rsq = x*x + y*y

  delta = rsq/(2.0*sigmasq)

  filter = filter0 * (1.0-delta) * exp(-delta)

end function filter
