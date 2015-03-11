"""
sample random points with in the RA DEC range of the survy

usage:
set DECALS_IMAGING and DECALS_CACHE
ensure the milky way extinction bricks are generated (build_sdfbricks.py)

On Edison, in a job, run

  aprun -n 24 python-mpi mpirandom.py 400000 random

  400000 is the number of random points 
  random.txt is the output.

on front end, run

  python mpirandom.py 400000 random

  (will be a lot slower than with MPI)

  output format is added at the top of file as a comment

 Yu Feng <yfeng1@berkeley.edu>
 Martin White <mwhite@berkeley.edu>
""" 
from model.datarelease import DataRelease
import numpy
import sys
from sys import argv

import filehandler

def main(comm):
    N = int(argv[1])
    output = argv[2]
    dr = DataRelease()
    bi = dr.brickindex

    ramin, ramax, decmin, decmax = dr.observed_range

    u1, u2 = numpy.random.random(size=(2, N))

    RA = (ramax - ramin) * u1 + ramin
    a = 0.5 * ((numpy.cos(decmax / 180. * numpy.pi)  + 1))
    b = 0.5 * ((numpy.cos(decmin / 180. * numpy.pi)  + 1))
    u2 = (a - b) * u2 + b
    
    DEC = numpy.arccos(2.0 * u2 - 1) * (180. / numpy.pi)

    (RA, DEC), junk = dr.brickindex.optimize((RA, DEC))

    if comm is not None:
        mystart = comm.rank * len(RA) // comm.size
        myend = (comm.rank + 1)* len(RA) // comm.size
    else:
        mystart = 0
        myend = len(RA)

    sl = slice(mystart, myend)
    coord = (RA[sl], DEC[sl])
    mydepth = numpy.zeros(len(RA[sl]), dtype=('f4', 3))
    myebv = numpy.zeros(len(RA[sl]), dtype='f4')
    mydepth[:, 0] = dr.readout(coord, dr.images['depth']['r'])
    mydepth[:, 1] = dr.readout(coord, dr.images['depth']['g'])
    mydepth[:, 2] = dr.readout(coord, dr.images['depth']['z'])
    mydepth[...] = 22.5 - 2.5 * numpy.log10(5 / mydepth ** 0.5)

    if comm is not None:
        depth = comm.gather(mydepth)
        ebv = comm.gather(myebv)
        if comm.rank == 0:
            depth = numpy.concatenate(depth, axis=0)
            ebv = numpy.concatenate(ebv, axis=0)
        else:
            return
    else:
        depth = mydepth
        ebv = myebv
    depth_r = depth[:, 0]
    depth_g = depth[:, 1]
    depth_z = depth[:, 2]

    # 5 sigma detection limit in mag
    filehandler.write_file(output, dict(RA=RA, 
                DEC=DEC, 
                mdepth_r=depth_r,
                mdepth_g=depth_g,
                mdepth_z=depth_z,
                ebv=ebv))
#    numpy.savetxt(output, zip(RA, DEC, depth[:, 1], depth[:, 2], depth[:, 4], ebv), fmt='%.6f', header='#ra dec r g z ebv')

if __name__ == '__main__':    
    from sys import argv
    if 'mpi' in sys.executable:
        from mpi4py import MPI
        main(MPI.COMM_WORLD)
    else:
        main(None)
