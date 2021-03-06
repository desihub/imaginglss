"""
The base Python object corresponding to an imaging brick.
This object contains only the meta-data, and it calls
imagerepo to handle looking up the pixel-level information.

"""
from __future__ import print_function

import numpy
from ..utils import wcs_tangent

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"




class Brick(object):
    """
    The base (immutable) object corresponding to an imaging brick.
    This object contains only the meta-data, and is immutable.

    It uses :py:class:`~model.imagerepo.ImageRepo` to handle 
    looking up the pixel-level information.

    Attributes
    ----------
    index   :   integer
        Internal index used by :py:class:`~imaginglss.brickindex.BrickIndex`
    id :   integer
        ID of the brick as in the Tractor catalogue (unique).
        Note that we do not use it in indexing.
    name : string_like 
        name of the brick as in the Tractor catalogue (unique)
    ra   : float
        center coordinate ra
    dec  : float
        center coordinate dec
    ra1  : float
        min coordinate ra
    dec1 : float
        min coordinate dec
    ra2  : float
        max coordinate ra
    dec2: float
        max coordinate dec

    """
    def __init__(self, index, id, name, ra, dec, ra1, ra2, dec1, dec2):
        """
        Initialize a brick object.  Brick is immutable.
        """
        self.index   = index
        self.id   = int(id)
        self.name = name
        self.ra   = ra
        self.dec  = dec
        self.ra1  = ra1
        self.ra2  = ra2
        self.dec1 = dec1
        self.dec2 = dec2
        deg = numpy.pi / 180.
        self.area = (numpy.sin(dec2 * deg) \
            - numpy.sin(dec1 * deg )) * (ra2 - ra1) * deg \
            * 129600 / numpy.pi / (4 * numpy.pi)
        #print(self.query(([self.ra], [self.dec])))
        #print(self.revert(([1799.], [1799.])))

    def __hash__(self):
        return self.id

    def __repr__(self):
        return ("Brick(id=%d, name=%s, area=%g, ra=%g, dec=%g, ...)"\
            % (self.id, self.name, self.area, self.ra, self.dec ))

    def readout(self, coord, repo, default=numpy.nan):
        """
        Return image values from an image repository.

        Parameters
        ----------
        coord  : array_like
            Coordinate of the pixels, coord=(RA, DEC)
        repo   : :py:class:`~model.imagerepo.ImageRepo`
            Image repository to read from. Refer to DataRelease.images.
        default : float
            Default value to return if the pixel is outside of the brick.

        Returns
        -------
        values : array_like
            values read out from repo.

        """
        coord     = numpy.array(coord)
        RA, DEC   = coord
        value     = numpy.empty(len(RA))
        value[...]= default

        img   = repo.open(self)

        xy = numpy.round(self.query(repo, coord)).astype('int32')
        mask = (xy < numpy.array(img.shape) \
            .reshape(2, 1)).all(axis=0)
        mask &= (xy >= 0).all(axis=0)
        x, y= xy[:, mask]
        value[mask] = img[y, x]
        return value
 
    def query(self, repo, coord):
        """
        Query the pixel coordinate in 'xy' system from RA/DEC system.

        Note that the returned array is in a x-first-varying order.
        To index the ndarray image, use [y, x].

        Parameters
        ----------
        repo   : :py:class:`~model.imagerepo.ImageRepo`
            image repository that contains the necessary transformation
            for this brick
        coord  : array_like 
            coord = (RA, DEC) tuple of arrays
        
        Returns
        -------
        xy     : array_like
            a (2xN) array, xy=(x, y).

        """
        meta      = repo.metadata(self)
        coord     = numpy.array(coord)
        xy = wcs_tangent.ang2pix_hdr(coord, meta, 
                zero_offset=True)
        return xy

    def revert(self, repo, xy):
        """
        Look up RA/DEC coord from pixels coordinate in 'xy' system.

        Parameters
        ----------
        repo   : :py:class:`~model.imagerepo.ImageRepo`
            image repository that contains the necessary transformation
            for this brick
        xy     : array_like
        
        Returns
        -------
        coord  : array_like
            coord=(RA, DEC) for each given position.
        """
        meta      = repo.metadata(self)
        xy = numpy.array(xy)
        coord = wcs_tangent.pix2ang_hdr(xy, meta, 
                zero_offset=True)
        return coord
