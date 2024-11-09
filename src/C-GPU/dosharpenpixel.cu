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
