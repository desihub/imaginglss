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

    def build(self):

        files = self.listfiles()

        bf = bigfile.BigFile(self.destdir, create=True)

        fulldtype = fits.read_table(files[0]).dtype

        dtype = [(column, fulldtype[column]) for column in self.columns if column is not 'BRICK_PRIMARY']

        dtype.append(('BRICK_PRIMARY', '?'))
        dtype = numpy.dtype(dtype)

        sizes = []
        for filename in files:
            sizes.append(fits.size_table(filename))

        sizes = numpy.array(sizes)
        offsets = numpy.concatenate([[0], numpy.cumsum(sizes)])

        blocks = {}
        for column in self.columns:
            blocks[column] = bf.create(column, dtype=dtype[column], size=sizes.sum(), Nfile=1)

        for i, filename in enumerate(files):
            onefile = fits.read_table(filename)
            for column in self.columns:
                onedata = numpy.empty(len(onefile), dtype=dtype[column])

                if column != 'BRICK_PRIMARY':
                    onedata[...] = onefile[column]
                else:
                    try:
                        onedata[...] = onefile[column]
                    except:
                        onedata[...] = True
                blocks[column].write(offsets[i], onedata)
            print(filename, 'done')
    def listfiles(self):
        return (list(sorted(glob(os.path.join(self.sweepdir, '*.fits'))))
             +  list(sorted(glob(os.path.join(self.sweepdir, '*.fits.gz')))))

