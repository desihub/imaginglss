# The base Python object corresponding to an imaging brick.
# This object contains only the meta-data, and it calls
# imagerepo to handle looking up the pixel-level information.

from __future__ import print_function

import numpy
from utils import wcs_tangent

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"




class Brick(object):
    """
    The base (immutable) object corresponding to an imaging brick.
    This object contains only the meta-data, and it calls
    imagerepo to handle looking up the pixel-level information.
    Major methods are "readout", "query" and "revert".
    """
    def __init__(self, id, name, ra, dec, ra1, ra2, dec1, dec2):
        """
        Initialize a brick object.  Brick is immutable.
        """
        self.id   = id
        self.name = name
        self.ra   = ra
        self.dec  = dec
        self.ra1  = ra1
        self.ra2  = ra2
        self.dec1 = dec1
        self.dec2 = dec2
        #print(self.query(([self.ra], [self.dec])))
        #print(self.revert(([1799.], [1799.])))

    def __hash__(self):
        return self.id

    def __repr__(self):
        return ("Brick(id=%d, name=%s, ra=%g, dec=%g, ...)"\
            % (self.id, self.name, self.ra, self.dec ))

    def readout(self, coord, repo, default=numpy.nan):
        """
        Return image values at coord=(RA, DEC) using the imagerepository repo.
        If the image in this brick does not cover this region,
        put in `default'.
        """
        meta      = repo.metadata(self)
        img       = repo.open(self)
        RA, DEC   = coord
        value     = numpy.empty(len(RA))
        value[...]= default
        coord     = numpy.array(coord)
        xy = numpy.int32(wcs_tangent.ang2pix_hdr(coord, meta, 
                zero_offset=True))
        xy2 = numpy.int32(self.query(coord))
        mask = (xy < numpy.array(img.shape) \
            .reshape(2, 1)).all(axis=0)
        mask &= (xy >= 0).all(axis=0)
        l = numpy.ravel_multi_index(xy[:, mask], img.shape, mode='raise')
        value[mask] = img.flat[l] 
        return value
 
    def query(self, coord):
        """
        Returns the xy index of pixels at coord
        coord can be:
            (RA, DEC) tuple of arrays
             array of shape (2xN) RA DEC
        In either case returns a (2xN) array
        """
        #FIXME: other types of input
        coord = numpy.array(coord)
        out = wcs_tangent.ang2pix(coord,
                CRPIX=(1800.5, 1800.5),
                CRVAL=(self.ra, self.dec),
                CD=(-7.27777777777778E-05,0,0, 7.27777777777778E-05),
               )
        return out

    def revert(self, xy):
        """
        Returns the RA, DEC index of pixels for coord
        coord can be:
            (x, y) tuple of arrays
            array of shape (2xN) x, y
        In either case returns a (2xN) array.
        """
        #FIXME: other types of input
        xy = numpy.array(xy)
        out = wcs_tangent.pix2ang(xy,
                CRPIX=(1800.5, 1800.5),
                CRVAL=(self.ra, self.dec),
                CD=(-7.27777777777778E-05,0,0,7.27777777777778E-05),
               )
        return out
