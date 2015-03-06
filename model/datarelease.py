# Python code to look at the imaging data, containing
# interfaces to deal with "bricks" and "catalogs".
# This is the "highest level" interface to the
# imaging data and makes use of several lower
# level objects.

from __future__ import print_function

import os
import os.path
import numpy
import glob
import re
from collections import namedtuple

from astropy.io import fits

from . import brickindex
from . import imagerepo
from . import catalogue


__author__ = "Yu Feng and Martin White"
__version__ = "0.9"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"




class Lazy(object):
    def __init__(self, calculate_function):
        self._calculate = calculate_function

    def __get__(self, obj, _=None):
        if obj is None:
            return self
        value = self._calculate(obj)
        setattr(obj, self._calculate.func_name, value)
        return value
    
def contains(haystack, needle):
    """ returns true if needle is in haystack """

    ind = haystack.searchsorted(needle)
    ind.clip(0, len(haystack) - 1, ind)
    return haystack[ind] == needle

class DataRelease(object):
    """
    The highest level interface into the data for a given imaging
    data release.  Uses several "helper" classes and has methods for
    looking at pixelized data or catalogs arranged in bricks.
    """
    BANDS = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'Y':5}

    BRICKS_FILENAME = 'bricks.fits'
    DEPTH_FILENAME = 'coadd/depth-%(brickid)d-%(band)s.fits.gz'
    IMAGE_FILENAME = 'coadd/image-%(brickid)d-%(band)s.fits' 
    MODEL_FILENAME = 'coadd/model-%(brickid)d-%(band)s.fits'

    @staticmethod
    def getimagefilename(PATTERN, band):
        def getfilename(brick):
            return PATTERN % dict(brickid=brick.id, band=band)
        return getfilename

    @staticmethod
    def gettractorfilename(brick):
        TRACTOR_FILENAME = 'tractor/tractor-%(brickid)d.fits'
        return TRACTOR_FILENAME % dict(brickid=brick.id) 

    @staticmethod
    def gettractorfilename(brick):
        TRACTOR_FILENAME = 'tractor/tractor-%(brickid)d.fits'
        return TRACTOR_FILENAME % dict(brickid=brick.id) 

    @staticmethod
    def tractorfilename_to_brick(filename, brickindex):
        return brickindex.get_brick(
            brickindex.search_by_id(
                int(re.search('-([0123456789]+)\.', 
                os.path.basename(filename)).group(1))))

    def __init__(self, root=None, cacheroot=None, version=None):
        if root is None:
            root = os.environ.get("DECALS_IMAGING", '.') 
        self.root = root

        if cacheroot is None:
            cacheroot = os.environ.get("DECALS_CACHE", '.') 

        self.cacheroot = os.path.join(cacheroot, version)

        brickdata = fits.open(os.path.join(self.root, self.BRICKS_FILENAME))[1].data

        self.brickindex = brickindex.BrickIndex(brickdata)

        self.bands = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'Y':5}

        self.images = {}
        for band in self.bands:
            self.images[band] = dict(
                DEPTH=imagerepo.ImageRepo(self.root, self.getimagefilename(self.DEPTH_FILENAME, band)),
                IMAGE=imagerepo.ImageRepo(self.root, self.getimagefilename(self.IMAGE_FILENAME, band)),
                MODEL=imagerepo.ImageRepo(self.root, self.getimagefilename(self.MODEL_FILENAME, band)),
            )

        self.observed_bricks = [ ]
        for roots, dirnames, filenames in os.walk(os.path.join(self.root, 'tractor'), followlinks=True):
            for filename in filenames:
                if not (filename.startswith('tractor') and filename.endswith('fits')): continue
                self.observed_bricks.append(self.tractorfilename_to_brick(filename, self.brickindex))

        self._observed_brickids = self.brickindex.search_by_id(
            [ brick.id for brick in self.observed_bricks ])
        arg = self._observed_brickids.argsort()
        self._observed_brickids = self._observed_brickids[arg]
        self.observed_bricks = numpy.array(self.observed_bricks)[arg]

        # approximate area in degrees. Currently a brick is 0.25 * 0.25 deg**2
        self.observed_area = 41253. * len(self.observed_bricks) / len(self.brickindex)

        self.catalogue = catalogue.Catalogue(
            os.path.join(self.cacheroot, 'catalogue'),
            [
            os.path.join(self.root, self.gettractorfilename(brick))
            for brick in self.observed_bricks])

        # footprint of the survey
        Footprint = namedtuple('Footprint', ['ramin', 'ramax', 'decmin', 'decmax'])
        self.footprint = Footprint(
            ramin=min([brick.ra1 for brick in self.observed_bricks]),
            ramax=max([brick.ra2 for brick in self.observed_bricks]),
            decmin=min([brick.dec1 for brick in self.observed_bricks]),
            decmax=max([brick.dec2 for brick in self.observed_bricks]),)
         
    def readout(self, coord, repo, default=numpy.nan):
        """ readout pixels at coord.
            
            example:
                corrd = (RA, DEC)
                repo is fetched from self.images

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
        # watch out bid + 1
        print(bid)
        mask = contains(self._observed_brickids, bid)
        print(mask)
        ra = RA[mask]
        dec = DEC[mask]
        coord, invarg = self.brickindex.optimize((ra, dec))
        bid = self.brickindex.query(coord)

        pixels = numpy.empty(len(bid), 'f8')
        pixels[:] = default

        
        ubid = numpy.unique(bid)

        for b in ubid:
            brick = self.brickindex.get_brick(b)
            first = bid.searchsorted(b, side='left')
            last = bid.searchsorted(b, side='right')
            sl = slice(first, last)

            img = brick.readout(coord[:, sl], repo, default=default)
            pixels[sl] = img
        #
        images[mask] = pixels[invarg]
            
        return images
