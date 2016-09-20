from __future__ import print_function
from sys import stdout
from glob import glob
import os.path
import numpy
import bigfile

from imaginglss.utils import fits

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"


def bcast_call(comm, func, *args, **kwargs):
    if comm.rank == 0:
        value = func(*args, **kwargs)
    else:
        value = None
    return comm.bcast(value)

class CacheBuilder(object):
    def __init__(self, sweepdir, cachedir, columns):
        self.sweepdir = sweepdir
        self.columns = columns
        self.destdir = os.path.join(cachedir, 'catalogue')
        print(self.destdir)
    # build cache from MPI
    def build(self, comm):
        files = bcast_call(comm, self.listfiles)

        start = comm.rank * len(files) // comm.size
        end = (comm.rank + 1) * len(files) // comm.size

        rows = 0

        bf = bigfile.BigFileMPI(comm, self.destdir, create=True)

        fulldtype = bcast_call(comm, lambda : fits.read_table(files[0]).dtype)

        dtype = [(column, fulldtype[column]) for column in self.columns if column is not 'BRICK_PRIMARY']

        dtype.append(('BRICK_PRIMARY', '?'))
        dtype = numpy.dtype(dtype)

        for column in self.columns:
            data = []
            for i, filename in zip(range(start, end), files[start:end]):
                size =fits.size_table(filename)
                onefile = fits.read_table(filename, subset=(column, 0, size))
                onedata = numpy.empty(len(onefile), dtype=dtype[column])
                if column != 'BRICK_PRIMARY':
                    onedata[...] = onefile
                else:
                    onedata[...] = True

                data.append(onedata)

            cdata = numpy.concatenate(data, axis=0)
            bf.create_from_array(column, cdata)

    def listfiles(self):
        return (list(sorted(glob(os.path.join(self.sweepdir, '*.fits'))))
             +  list(sorted(glob(os.path.join(self.sweepdir, '*.fits.gz')))))

