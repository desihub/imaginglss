from model.datarelease import DataRelease
import numpy
import sys
from sys import stdout
from sys import argv
def main(comm):
    N = int(argv[1])
    dr = DataRelease()
    bi = dr.brickindex

    ramin, ramax, decmin, decmax = dr.footprint

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
    mydepth = numpy.zeros(len(RA[sl]), dtype=('f4', 6))
    mydepth[:, 1] = dr.readout(coord, dr.images['depth']['r'])
    mydepth[:, 2] = dr.readout(coord, dr.images['depth']['g'])
    mydepth[:, 4] = dr.readout(coord, dr.images['depth']['z'])
    if comm is not None:
        depth = comm.gather(mydepth)
        if comm.rank == 0:
            depth = numpy.concatenate(depth, axis=0)
        else:
            return
    else:
        depth = mydepth

    numpy.savetxt(stdout, zip(RA, DEC, depth[:, 1], depth[:, 2], depth[:, 4]), header='#ra dec r g z')

if __name__ == '__main__':    
    from sys import argv
    if 'mpi' in sys.executable:
        from mpi4py import MPI
        main(MPI.COMM_WORLD)
    else:
        main(None)
