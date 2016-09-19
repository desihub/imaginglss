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

        print(sorted(fulldtype.names))
        dtype = [(column, fulldtype[column]) for column in self.columns if column is not 'BRICK_PRIMARY']

        dtype.append(('BRICK_PRIMARY', '?'))

        data = []
        for i, filename in zip(range(start, end), files[start:end]):
            onefile = fits.read_table(filename)
            onedata = numpy.empty(len(onefile), dtype=dtype)
            for column in self.columns:
                if column != 'BRICK_PRIMARY':
                    onedata[column][...] = onefile[column]
                    if not column in onefile.dtype.names:
                        raise KeyError("column `%s` not found in sweep files")
                else:
                    onedata[column][...] = True

            data.append(onedata)

        for column in self.columns:
            cdata = numpy.concatenate([data1[column] for data1 in data], axis=0)
            #print(cdata.dtype.itemsize, cdata.shape, cdata.dtype.str)
            bf.create_from_array(column, cdata)

    def listfiles(self):
        return list(sorted(glob(os.path.join(self.sweepdir, '*.fits'))))

