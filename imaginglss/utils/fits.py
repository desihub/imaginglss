from __future__ import print_function

"""
    Simple routines that provide a uniform interface for reading
    fits hdus via either astropy.io.fits or fitsio.

    Three functions are provided:

        read_image, read_table, read_metadata

    Warning: A copy of data is returned to detach any backreference to 
    the original fits library. Use with caution with extremely 
    large files. (> half of avail memory on machine)
"""

__author__ = "Yu Feng and Martin White"
__version__ = "0.9"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"
__all__ = ['read_image', 'read_table', 'read_metadata', 'size_table']

import numpy

def use_astropy():
    """ select the astropy backend """
    global read_image
    global read_table
    global size_table
    global read_metadata
    global backend 

    backend = "astropy"

    from astropy.io import fits

    def read_image(filename, hdu=0):
        """
            Read the zeroth image HDU from a fits file

        """
        file = fits.open(filename)
        #    A copy of data is made before the file object
        #    is dereferenced. This is to ensure no back references
        #    to the fits object and file gets closed in a timely
        #    fashion.
        return numpy.array(file[hdu].data, copy=True)

    def read_table(filename, hdu=1, subset=None):
        """ 
            Read the first HDU table from a fits file
        """
        file = fits.open(filename)
        if subset is not None:
            column, start, end = subset
            return numpy.array(file[hdu][column][start:end], copy=True)
        return numpy.array(file[hdu].data, copy=True)

    def size_table(filename, hdu=1):
        """ 
            Read the first HDU table from a fits file
        """
        file = fits.open(filename)
        return file[hdu].shape[0]

    def read_metadata(filename, hdu=0):
        """ 
            Read the metadata of a HDU table
        """
        file = fits.open(filename)
        return dict(file[hdu].header)

def use_fitsio():
    """ select the fitsio backend """
    global read_image
    global read_table
    global size_table
    global read_metadata
    global backend 

    backend = "fitsio"

    from fitsio import FITS

    def read_image(filename, hdu=0):
        """
            Read the zeroth image HDU from a fits file

        """
        file = FITS(filename, upper=True)
        #    A copy of data is made before the file object
        #    is dereferenced. This is to ensure no back references
        #    to the fits object and file gets closed in a timely
        #    fashion.
        return numpy.array(file[hdu].read(), copy=True)

    def read_table(filename, hdu=1, subset=None):
        """ 
            Read the first HDU table from a fits file
        """
        file = FITS(filename, upper=True)
        if subset is not None:
            column, start, end = subset
            return numpy.array(file[hdu][column][start:end], copy=True)
        else:
            return numpy.array(file[hdu].read(), copy=True)

    def size_table(filename, hdu=1):
        """ 
            Read the first HDU table from a fits file
        """
        file = FITS(filename, upper=True)
        return file[hdu].get_nrows()

    def read_metadata(filename, hdu=0):
        """ 
            Read the metadata of a HDU table
        """
        file = FITS(filename, upper=True)
        return dict(file[hdu].read_header())

_priorities = [use_fitsio, use_astropy]

for backend in _priorities:
    try:
        backend()
        break
    except ImportError:
        continue

if __name__ == '__main__':
    use_fitsio()
    im0 = read_image('depth-370143-g.fits.gz')
    meta0 = read_metadata('depth-370143-g.fits.gz')
    tbl0 = read_table('tractor-370143.fits')
    assert backend == 'fitsio'
    use_astropy()
    im1 = read_image('depth-370143-g.fits.gz')
    meta1 = read_metadata('depth-370143-g.fits.gz')
    tbl1 = read_table('tractor-370143.fits')

    assert backend == 'astropy'

    assert (im0[...] == im1[...]).all()
    assert set(tbl0.dtype.names) == set(tbl1.dtype.names)
    assert set(meta0.keys()) == set(meta1.keys())
    for k in meta0.keys():
        if k == "COMMENT": continue
        print (k, '::::', meta0[k], '::::', meta1[k])
        assert meta0[k] == meta1[k]
