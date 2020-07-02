import numpy as np


def pgmsize(fname):
    f = open(fname, 'r')
    f.readline()
    f.readline()
    line = f.readline()
    cols = line.split()
    f.close()
    return (int(cols[0]),int(cols[1])) if 2 == len(cols) else (0,0)
    

def pgmread(fname, pixmap, nxmax, nymax):

    n = pgmsize(fname)
    nx = n[0]
    ny = n[1]
    if (nx > nxmax or ny > nymax):
        print "pgmread: image larger than array.\n"
        print "nxmax, nymax, nx, ny = "
        print str(nxmax), ", ", str(nymax), ", ", str(nx), ", ", str(ny), ".\n"
        raise SystemExit

    f = open(fname, 'r')

    i = 0
    j = ny-1
    fline_cnt = 1
    
    for line in f:
        if fline_cnt > 4:
            cols = line.split()
            for k in range(len(cols)):
                pixmap[i][j] = int(cols[k])
                i = i+1 if i < nx-1 else 0
                j = j if i > 0 else j-1
                
        fline_cnt = fline_cnt + 1
        
    f.close()

    return n
  

"""
Routine to write a PGM image file from a 2D floating point array
x[nx][ny]. Because of the way C handles (or fails to handle!)
multi-dimensional arrays we have to cast the pointer to void.
"""
def pgmwrite(fname, pixmap, nx, ny):

    thresh = float(255.0)
    PIXPERLINE = 16

    # find the max and min absolute values of the array
    xmin = abs(pixmap[0][0])
    xmax = xmin
    for i in range(nx):
        for j in range(ny):
            xabs = abs(pixmap[i][j])
            if xabs < xmin:
                xmin = xabs
            elif xabs > xmax:
                xmax = xabs
  

    f = open(fname,'w')
    
    f.write("P2\n")
    f.write("# Written by pgmwrite\n")
    f.write(str(nx) + " " + str(ny) + "\n")
    f.write(str(int(thresh)) + "\n")
    
    k = 0

    for j in range(ny-1,-1,-1):
        for i in range(nx):

            tmp = pixmap[i][j]
      
            # scale the value appropriately so it lies between 0 and thresh
            if (xmin < 0 or xmax > thresh):
                tmp = int((thresh*((abs(tmp-xmin))/(xmax-xmin))) + 0.5)
            else:
                tmp = int(abs(tmp) + 0.5)
      
            grey = tmp

            # increase the contrast by boosting the lower values?
            #grey = thresh * sqrt(tmp/thresh)

            f.write(str(grey) + " ")

            if (0 == (k+1) % PIXPERLINE):
                f.write("\n")

            k = k + 1
    
    if (0 != k % PIXPERLINE):
        f.write("\n")
      
    f.close()
