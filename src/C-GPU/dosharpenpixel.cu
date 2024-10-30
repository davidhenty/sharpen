#include <stdio.h>
#include <math.h>

__device__
double filter(int d, int i, int j);

__global__
void dosharpenpixel(int nx, int ny, int d, 
                    double *convolution, double *fuzzyPadded)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;

  int idx, idxp, idxppk;

  if (i < nx && j < ny)
    {
      idx  = i*ny + j;
      idxp = (i+d)*(ny+2*d) + (j+d);

      for (int k= -d; k <= d; k++)
        {
          idxppk = idxp + k*(ny+2*d);

          for (int l= -d; l <= d; l++)
            {
                  convolution[idx] =   convolution[idx]
                                     + filter(d,k,l)*fuzzyPadded[idxppk+l];
            }
        }
    }
}

__global__
void dosharpenpixelinline(int nx, int ny, int d, 
                          double *convolution, double *fuzzyPadded)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;

  int idx, idxp, idxppk;
  double fval;

  double rd4sq, rsq, sigmad4sq, sigmasq, x, y, delta;

  int d4 = 4;

  double sigmad4 = 1.4;
  double filter0 = -40.0;

  rd4sq = d4*d4;
  rsq   = d*d;

  sigmad4sq = sigmad4*sigmad4;
  sigmasq   = sigmad4sq * (rsq/rd4sq);

  if (i < nx && j < ny)
    {
      idx  = i*ny + j;
      idxp = (i+d)*(ny+2*d) + (j+d);

      /*      if ((i < nx/2) && (j < ny/2) || (i >= nx/2) && (j >= ny/2))
        {
          convolution[idx] = 127;
        }
      else
        {
          convolution[idx] = 255;
          } */
      
      for (int k= -d; k <= d; k++)
        {
          idxppk = idxp + k*(ny+2*d);

          for (int l= -d; l <= d; l++)
            {
              x = (double) k;
              y = (double) l;

              rsq = x*x + y*y;

              delta = rsq/(2.0*sigmasq);

              fval = filter0 * (1.0-delta) * exp(-delta);

              //              printf("i, j, k, l, idx, idxppk+l = %d, %d, %d, %d, %d, %d, convolution[idx], filter(d,k,l), fuzzyPadded[idxppk+l] = %f, %f, %f\n", i, j, k, l, idx, idxppk+l, convolution[idx], fval, fuzzyPadded[idxppk+l]);

                  convolution[idx] = convolution[idx]
                               + fval*fuzzyPadded[idxppk+l];
            }
        }
    }
}
