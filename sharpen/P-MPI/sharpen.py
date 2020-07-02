"""
Program to sharpen an image by convolving with a filter function. The
filter is a combination of a Gaussian (to remove noise) and a Laplacian
(to detect the edges). Input and output is via Portable Grey Map (PGM)
files - note that the input file must have a specific header format.

In this version of the program the image processing is parallelised using
MPI processes and replicated data. The master process reads in the fuzzy 
image and distributes it to the other processes using an MPI communication. 
The convolution computation is distributed over all processes and the result 
communicated back to the master process. Finally the master process adds the 
convolution result to the fuzzy image and writes the resulting sharp image to 
file.

David Henty, EPCC, September 2009
Arno Proeme, EPCC, March 2013 (minor modifications)
Dominic Sloan-Murphy, EPCC, November 2013 (more minor modifications)
Michael Bareford, EPCC, March 2016 (converted to python)
"""
from mpi4py import MPI
import sharpenio as io
import sharpenalg as alg
import numpy as np


infile = "fuzzy.pgm"
outfile = "sharpened.pgm"

comm = MPI.COMM_WORLD

size = comm.Get_size()
rank = comm.Get_rank()

comm.Barrier()

tstart = MPI.Wtime()

if 0 == rank:
    print "\nImage sharpening code running on ", str(size), " processor(s)\n"
    print "Input file is ", infile, "\n"
    pix = np.array(io.pgmsize(infile))
    print "Image size is ", str(pix[0]), " x ", str(pix[1]), "\n"
else:
    pix = np.zeros(2, dtype=np.int)
          
comm.Bcast([pix, MPI.INT], root=0)

alg.dosharpen(pix, infile, outfile)

comm.Barrier()

tstop = MPI.Wtime()
tdiff = tstop - tstart

if 0 == rank:
    print "Overall run time was ", str(tdiff), " seconds"
