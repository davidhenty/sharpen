from math import *
import sharpenio as io
import numpy as np
import time

"""
Sets the linear range of the sharpen filter as measured from any given pixel:
only pixels within a (2d+1)*(2d+1) square centered on the pixel are used to 
compute its new value.
"""
def dosharpen(pix, infile, outfile):
    nx = pix[0]
    ny = pix[1]
    d = 8
    norm = (2*d-1)*(2*d-1)  
    scale = 2.0
  
    # arrays automatically initialised to zero
    fuzzy = np.zeros((nx, ny), dtype=np)                    # stores the fuzzy input image when it is first read in from file
    fuzzyPadded = np.zeros((nx+2*d, ny+2*d), dtype=np.double)   # stores the fuzzy input image plus additional border padding
    convolution = np.zeros((nx, ny), dtype=np.double)           # stores the convolution of the filter with the full fuzzy image
    sharp = np.zeros((nx, ny), dtype=np.double)                 # stores the sharpened image obtained by adding rescaled convolution to the fuzzy image
    sharpCropped = np.zeros((nx-2*d, ny-2*d), dtype=np.double)  # stores the sharpened image cropped to remove a border layer distorted by the algorithm

    print("Using a filter of size ", str(2*d+1), " x ", str(2*d+1), "\n")
    print("Reading image file ", infile)
    npix = np.array(io.pgmread(infile, fuzzy, nx, ny))
    print("... done\n")

    xpix = npix[0]
    ypix = npix[1]
    
    if (0 == xpix or 0 == ypix or nx != xpix or ny != ypix):
        print("Error reading ", infile)
        raise SystemExit
  
    # transfer fuzzy image into padded array
    for i in range(nx):
        for j in range(ny):
            fuzzyPadded[i+d][j+d] = fuzzy[i][j]

    print("Starting calculation ..." )

     # print(out current core and node location)

    tstart = time.time()

    pixcount = 0

#    print("start: convolution[", 1+d, "][", d, "] = ", convolution[1+d][d])

    for i in range(nx):
        for j in range(ny):
            for k in range(-d,d+1):
                for l in range(-d,d+1):
                    convolution[i][j] = convolution[i][j] + filter(d,k,l)*fuzzyPadded[i+d+k][j+d+l]
#                    if (i == d+1 and j == d):
#                        print("k = ", k, ", l = ", l, ", convolution[", i, "][", j, "] = ", convolution[1+d][d])

            pixcount += 1
      
    tstop = time.time()
    tdiff = tstop - tstart
    print("... finished\n")
      
    # apply the filter and writes the sharpened image to file
    # add rescaled convolution to fuzzy image to obtain sharp image
    for i in range(nx):
        for j in range(ny):
            sharp[i][j] = fuzzyPadded[i+d][j+d] - scale/norm * convolution[i][j]

    print("Writing output file: ", outfile)

    # only save the core of the sharpened image to remove edge effects
    for i in range(d,nx-d):
        for j in range(d,ny-d):
            sharpCropped[i-d][j-d] = sharp[i][j]
              
    io.pgmwrite(outfile, sharpCropped, nx-2*d, ny-2*d)

    print("... done\n")
    print("Calculation time was ", '{0:.2f}'.format(tdiff), " seconds")
    
       
def filter(d, i, j):
    d4 = 4
    sigmad4 = float(1.4)
    filter0 = float(-40.0)

    rd4sq = float(d4*d4)
    rsq   = float(d*d)

    sigmad4sq = sigmad4*sigmad4
    sigmasq   = sigmad4sq*(rsq/rd4sq)

    x = float(i)
    y = float(j)

    rsq = x*x + y*y
    delta = rsq/(2.0*sigmasq)

    return filter0*(1.0-delta)*exp(-delta)
