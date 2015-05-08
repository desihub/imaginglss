"""
Python code to look at the imaging data, containing
interfaces to deal with "bricks" and "catalogs".
This is the "highest level" interface to the
imaging data and makes use of several lower
level objects.

"""
from __future__ import print_function

import os
import os.path
import numpy
import glob
import re
from collections import namedtuple

from ..utils import fits

from . import brickindex
from . import imagerepo
from . import catalogue
from . import schema

__author__ = "Yu Feng and Martin White"
__version__ = "0.9"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"

import warnings


class Lazy(object):
    """ Lazy initialization of object attributes.

    """
    def __init__(self, calculate_function):
        self._calculate = calculate_function

    def __get__(self, obj, _=None):
        if obj is None:
            return self
        value = self._calculate(obj)
        setattr(obj, self._calculate.func_name, value)
        return value
    
def contains(haystack, needle):
    """ test if needle is in haystack. 

        Parameters
        ----------
        haystack  : array_like
            Sorted array
        needle    : array_like
            items to look for

        Returns
        -------
        mask : array_like
             mask[i] is true only if needle[i] is in haystack;
    
        Examples
        --------
        >>> contains([1, 2, 3], [2])
        [True]

    """
    haystack = numpy.asarray(haystack)
    ind = haystack.searchsorted(needle)
    if len(haystack) == 0:
        # build a False array of the right shape
        return ind != ind
    ind.clip(0, len(haystack) - 1, ind)
    return haystack[ind] == needle

class Footprint(object):
    """ footprint of a data release.
        
        Attributes
        ----------
        bricks : list of model.brick.Brick
            A list of Bricks that are covered by the footprint
        range :  tuple
            The range of RA and DEC of all bricks
            (ramin, ramax, decmin, decmax)
        area  : float
            Covered outline area in square degrees
    """
    def __init__(self, bricks, brickindex):
        self.bricks = bricks 
        self.area = 41253. * len(bricks) / len(brickindex)
        self._covered_brickids = numpy.array(
                    [b.index for b in bricks], dtype='i8')

        # range of ra dec of covered bricks
        FootPrintRange = namedtuple('FootPrintRange', ['ramin', 'ramax', 'decmin', 'decmax'])
        if len(bricks) == 0:
            self.range = FootPrintRange(ramin=0, ramax=0, decmin=0, decmax=0)
        else:
            self.range = FootPrintRange(
                ramin=min([brick.ra1 for brick in self.bricks]),
                ramax=max([brick.ra2 for brick in self.bricks]),
                decmin=min([brick.dec1 for brick in self.bricks]),
                decmax=max([brick.dec2 for brick in self.bricks]),)

        self.brickindex = brickindex

    def __repr__(self):
        return "Footprint: len(bricks)=%d , area=%g degrees, range=%s" % (
                len(self.bricks),
                self.area,
                str(self.range)
            )

    def filter(self, coord):
        """ Remove coordinates that are not covered by the footprint 

            Parameters
            ----------
            coord : array_like
                must be compatible with (RA, DEC)

            Returns
            -------
            coord_in_footprint : array_like
                items in the input coord that is in the footprint
        
        """
        coord = numpy.array(coord)
        bid = self.brickindex.query(coord)
        mask = contains(self._covered_brickids, bid)
        return coord[:, mask]

class DataRelease(object):
    """
    The highest level interface into the data for a given imaging
    data release.  Uses several "helper" classes and has methods for
    looking at pixelized data or catalogs arranged in bricks.

    Attributes
    ----------

    brickindex : :py:class:`~model.brickindex.BrickIndex`
        an index object of all of the bricks (covering the entire sky)
    bands      : dict
        a dictionary translating from band name to integer used in Tractor catalogue
    catalogue  : :py:class:`~model.catalogue.Catalogue`
        the concatenated tractor catalogue, accessed by attributes.
    extinction : array_like
        an array storing the extinction coeffcients. 
        The array matches the defination of `DECAM_FLUX` column of 
        the catalogue, u, g, r, i, z, Y.
    images     : :py:class:`~model.imagerepo.ImageRepo`
        Image repositories. These image repositories are used by
        py:meth:`readout`.
    footprint  : :py:class:`Footprint`
        the footprint of the data release.

    Examples
    --------
    >>> dr = DataRelease()
    >>> dr.images['depth']['r']  # r band depth images.
    >>> dr.images['image']['r']  # r band coadd images.
    >>> dr.images['model']['r']  # r band tractor model images.


    """
    def __init__(self, root, cache, version):
        root = os.path.normpath(root)

        cache = os.path.join(cache, version)

        self.root = root
        self.cache = cache

        if not hasattr(schema, version):
            raise KeyError("Data Release of version %s is not supported" % version)

        myschema = getattr(schema, version)

        try:
            os.makedirs(self.cache)
        except :
            pass

        brickdata = fits.read_table(os.path.join(self.root, myschema.BRICKS_FILENAME))

        self.brickindex = brickindex.BrickIndex(brickdata)

        self.bands = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'Y':5}

        # E(B-V) to ugrizY bands, SFD98; used in tractor
        self.extinction = numpy.array([3.995, 3.214, 2.165, 1.592, 1.211, 1.064], dtype='f8')\
            .view(dtype=[(band, 'f8') for band in 'ugrizY'])[0]

        self.images = {}
        image_filenames = myschema.format_image_filenames()
        for image in image_filenames:
            if isinstance(image_filenames[image], dict):
                self.images[image] = {}
                for band in image_filenames[image]:
                    self.images[image][band] = imagerepo.ImageRepo(self.root, image_filenames[image][band])
            else:
                self.images[image] = imagerepo.ImageRepo(self.root, image_filenames[image])
        def image_filename(brick, 
            PATTERN='aux/%(pre)s/%(brickname)s/decals-%(brickname)s-ebv.fits'):
            return PATTERN % dict(pre=brick.name[:3], brickname=brick.name)
        self.images['ebv'] = imagerepo.ImageRepo(self.cache, image_filename)
            
        try: 
            _covered_brickids = numpy.fromfile(
                os.path.join(self.cache, 'covered_brickids.i8'), dtype='i8')
        except IOError:
            
            _covered_brickids = [ ]
            for roots, dirnames, filenames in \
                os.walk(os.path.join(self.root, 'tractor'), followlinks=True):
                for filename in filenames:
                    try:
                        _covered_brickids.append(
                            myschema.parse_filename(filename, self.brickindex))
                    except ValueError:
                        pass 
            _covered_brickids = numpy.array(_covered_brickids, dtype='i8')
            _covered_brickids.tofile(os.path.join(self.cache, 'covered_brickids.i8'))
            
        self._covered_brickids = _covered_brickids

        # the list of covered bricks must be sorted.
        self._covered_brickids.sort()

        bricks = [self.brickindex.get_brick(bid) for bid in self._covered_brickids]
        self.footprint = Footprint(bricks, self.brickindex) # build the footprint property

        catalogue_filenames = [
                os.path.join(self.root, 
                myschema.format_catalogue_filename(brick))
                for brick in self.footprint.bricks]

        self.catalogue = catalogue.Catalogue(
            cachedir=os.path.join(self.cache, 'catalogue'),
            filenames=catalogue_filenames,
            aliases=myschema.CATALOGUE_ALIASES
            )

    def readout(self, coord, repo, default=numpy.nan, ignore_missing=False):
        """ Readout pixels from an image.
            
            Parameters
            ----------
            coord  : array_like
                coordinates of the pixels, (RA, DEC)
            repo   : ImageRepo
                the images to read from.
            default : scalar
                value to return if the pixel is not in the footprint.
            image_missing : boolean
                When ignore_missing is True, missing brick files are treated
                as not in the footprint.

            Notes
            -----
            This is here, because we want to query multiple images
            at the same time. 
            It is also convenient to have it here to make use of
            brickindex. (ImageRepo is then just a stub with no business logic)

            Otherwise it makes more sense to
            have readout in ImageRepo.
        """
        RA, DEC = coord
        images = numpy.empty(len(RA), dtype='f4')
        images[...] = default

        bid = self.brickindex.query((RA, DEC))

        mask = contains(self._covered_brickids, bid)

        ra = RA[mask]
        dec = DEC[mask]
        if len(ra) == 0:
            # do not try to work if no point is within the
            # survey 
            return images

        coord, invarg = self.brickindex.optimize((ra, dec), return_inverse=True)
        bid = self.brickindex.query(coord)

        pixels = numpy.empty(len(bid), 'f8')
        pixels[:] = default

        
        ubid = numpy.unique(bid)

        for b in ubid:
            if b not in self._covered_brickids:
                continue
            brick = self.brickindex.get_brick(b)
            first = bid.searchsorted(b, side='left')
            last = bid.searchsorted(b, side='right')
            sl = slice(first, last)

            try:
                img = brick.readout(coord[:, sl], repo, default=default)
                #print( 'readout', b, img)
                pixels[sl] = img
            except IOError as e:
                #print( 'readout', b, self.brickindex.get_brick(b), 'error', e, coord)
                if not ignore_missing:
                    raise
                else:
                    warnings.warn(str(e), stacklevel=2)
        #
        images[mask] = pixels[invarg]
            
        return images
