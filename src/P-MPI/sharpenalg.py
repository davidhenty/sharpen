from mpi4py import MPI
from math import *
import sharpenio as io
import numpy as np


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
    fuzzy = np.zeros((nx, ny), dtype=np.int)                    # stores the fuzzy input image when it is first read in from file
    fuzzyPadded = np.zeros((nx+2*d, ny+2*d), dtype=np.double)   # stores the fuzzy input image plus additional border padding
    convolutionPartial = np.zeros((nx, ny), dtype=np.double)    # stores the convolution of the filter with parts of the fuzzy image computed by individual processes
    convolution = np.zeros((nx, ny), dtype=np.double)           # stores the convolution of the filter with the full fuzzy image
    sharp = np.zeros((nx, ny), dtype=np.double)                 # stores the sharpened image obtained by adding rescaled convolution to the fuzzy image
    sharpCropped = np.zeros((nx-2*d, ny-2*d), dtype=np.double)  # stores the sharpened image cropped to remove a border layer distorted by the algorithm

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    if 0 == rank:
        print "Using a filter of size ", str(2*d+1), " x ", str(2*d+1), "\n"
        print "Reading image file ", infile
        npix = np.array(io.pgmread(infile, fuzzy, nx, ny))
        print "... done\n"
    else:
        npix = np.zeros(2, dtype=np.int)
        
    comm.Bcast([npix, MPI.INT], root=0)
    xpix = npix[0]
    ypix = npix[1]
    
    if (0 == xpix or 0 == ypix or nx != xpix or ny != ypix):
        if 0 == rank:
            print "Error reading ", infile
        raise SystemExit
  
    # broadcast the pixel image to all processes
    comm.Bcast([fuzzy, MPI.INT], root=0)
  
    # transfer fuzzy image into padded array
    for i in range(nx):
        for j in range(ny):
            fuzzyPadded[i+d][j+d] = fuzzy[i][j]

    if 0 == rank:
        print "Starting calculation ..." 

    comm.Barrier()

    # print out current core and node location
    node_name = MPI.Get_processor_name()
    print "Rank ", str(rank), " of ", str(size), " on node ", node_name
    
    tstart = MPI.Wtime()

    pixcount = 0

    for i in range(nx):
        for j in range(ny):
            # computation of convolution allocated to processes using simple cyclic distribution
            # i.e. consecutively ranked processes take turns computing convolution for consecutive pixels
            if rank == pixcount%size:
                for k in range(-d,d):
                    for l in range(-d,d):
                        convolutionPartial[i][j] = convolutionPartial[i][j] + filter(d,k,l)*fuzzyPadded[i+d+k][j+d+l]
                        
            pixcount += 1
      
    comm.Barrier()

    tstop = MPI.Wtime()
    time = tstop - tstart

    if 0 == rank:
        print "... finished\n"
      
    # gather the partial convolution results computed by individual processes
    comm.Reduce([convolutionPartial, MPI.DOUBLE], [convolution, MPI.DOUBLE], op=MPI.SUM, root=0)

    # the master process applies the filter and writes the sharpened image to file
    if 0 == rank:
        # add rescaled convolution to fuzzy image to obtain sharp image
        for i in range(nx):
            for j in range(ny):
                sharp[i][j] = fuzzyPadded[i+d][j+d] - scale/norm * convolution[i][j]
      
        print "Writing output file: ", outfile

        # only save the core of the sharpened image to remove edge effects
        for i in range(d,nx-d):
            for j in range(d,ny-d):
                sharpCropped[i-d][j-d] = sharp[i][j]
              
        io.pgmwrite(outfile, sharpCropped, nx-2*d, ny-2*d)

        print "... done\n"
        print "Calculation time was ", str(time), " seconds"
    
       
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
