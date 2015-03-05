#!/usr/bin/env python
#
# Python script to generate a Mangle polygon mask containing
# the "footprint" of the survey.
# The file is generated with one polygon per brick.
# Depending on what this mask is desired for, this may need to
# be snapped and balkanized, and/or pixelized.
#

from __future__ import print_function

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mwhite@berkeley.edu"

import math  as M
from model.datarelease import DataRelease





def write_poly(capslist,arealist,fname):
    """
    write_poly(capslist,area,fname):
    Given a list of caps (which are themselves lists) and areas per cap
    writes the simple Mangle mask as an ascii ".ply" file.
    """
    ff = open(fname,"w")
    ff.write("%d polygons\n"%len(capslist))
    ipoly = 0
    for caps,area in zip(capslist,arealist):
        ff.write("polygon %d ( %d caps, 1.0 weight, %15.10f str):\n"%\
          (ipoly,len(caps),area))
        for icap in caps:
            ff.write("%15.12f %15.12f %15.12f %15.12f\n"%\
              (icap[0],icap[1],icap[2],icap[3]))
        ipoly+=1
    ff.close()
    #




def make_caps(ramin,ramax,decmin,decmax):
    """
    make_caps(ramin,ramax,decmin,decmax):
    Returns a list of caps (and an area) which bound the region.  This
    list of caps forms the body of the polygon for a single brick.
    """
    # Convert from ra,dec in decimal degrees to (theta,phi)
    ppmin = M.pi/180.*ramin
    ppmax = M.pi/180.*ramax
    ttmin = M.pi/180.*(90.0-decmax)	# Note reversal.
    ttmax = M.pi/180.*(90.0-decmin)	# Note reversal.
    # Build up the caps.
    caps = []
    caps.append([-M.sin(ppmin), M.cos(ppmin),0., 1.0])
    caps.append([-M.sin(ppmax), M.cos(ppmax),0.,-1.0])
    caps.append([0.,0., 1.,-(1-M.cos(ttmin))])
    caps.append([0.,0.,-1.,-(1+M.cos(ttmax))])
    # Work out the area, in steradians.
    area  = M.fabs((ppmax-ppmin)*(M.cos(ttmax)-M.cos(ttmin)))
    return( (caps,area) )
    #











def make_polygon_file(fn="footprint.ply"):
    """
    make_polygon_file(fn="footprint.ply"):
    Does the actual work of making the footprint and writing the file.
    Reads in the meta-data for each brick and generates a Mangle polygon
    for each one.  Once the list is constructed, calls write_poly to write
    the mask file.
    """
    capslist,arealist = [],[]
    dr = DataRelease(version='EDR')
    for brickid in dr.observed_brickids:
        bb        = dr.brickindex.get_brick(brickid)
        caps,area = make_caps(bb.ra1,bb.ra2,bb.dec1,bb.dec2)
        capslist.append(caps)
        arealist.append(area)
    write_poly(capslist,arealist,fn)
    #





if __name__=="__main__":
    print("Generating polygon file.")
    make_polygon_file()
    #
