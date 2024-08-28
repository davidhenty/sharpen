#include <math.h>

double filter(int d, int i, int j)
{
  double rd4sq, rsq, sigmad4sq, sigmasq, x, y, delta;

  int d4 = 4;

  double sigmad4 = 1.4;
  double filter0 = -40.0;

  rd4sq = d4*d4;
  rsq   = d*d;

  sigmad4sq = sigmad4*sigmad4;
  sigmasq   = sigmad4sq * (rsq/rd4sq);

  x = (double) i;
  y = (double) j;

  rsq = x*x + y*y;

  delta = rsq/(2.0*sigmasq);

  return(filter0 * (1.0-delta) * exp(-delta));
}
