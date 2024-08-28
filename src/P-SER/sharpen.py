"""
Program to sharpen an image by convolving with a filter function. The
filter is a combination of a Gaussian (to remove noise) and a Laplacian
(to detect the edges). Input and output is via Portable Grey Map (PGM)
files - note that the input file must have a specific header format.

This is a serial version.

David Henty, EPCC, September 2009
Arno Proeme, EPCC, March 2013 (minor modifications)
Dominic Sloan-Murphy, EPCC, November 2013 (more minor modifications)
Michael Bareford, EPCC, March 2016 (converted to python)
"""
import sharpenio as io
import sharpenalg as alg
import numpy as np
import time


infile = "fuzzy.pgm"
outfile = "sharpened.pgm"

tstart = time.time()

print("\nImage sharpening code running in serial")
print("Input file is ", infile, "\n")
pix = np.array(io.pgmsize(infile))
print("Image size is ", str(pix[0]), " x ", str(pix[1]), "\n")

alg.dosharpen(pix, infile, outfile)

tstop = time.time()
tdiff = tstop - tstart

print("Overall run time was ", '{0:.2f}'.format(tdiff), " seconds")
