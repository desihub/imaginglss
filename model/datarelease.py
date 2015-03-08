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

from utils import fits

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

class EDR:
    BRICKS_FILENAME = 'bricks.fits'
    CATALOGUE_ALIASES = [('EXTINCTION', 'DECAM_EXTINCTION')]

    @staticmethod
    def format_image_filenames():
        images = {'depth': 'coadd/depth-%(brickid)d-%(band)s.fits.gz',
         'image': 'coadd/image-%(brickid)d-%(band)s.fits',
         'model': 'coadd/model-%(brickid)d-%(band)s.fits'}
        imagerepos= {}
        for image in images:
            imagerepos[image] = {}
            for band in 'rgz':
                PATTERN = images[image]
                def getfilename(brick):
                    return PATTERN % dict(brickid=brick.id, band=band)
                imagerepos[image][band] = getfilename
        return imagerepos

    @staticmethod
    def format_catalogue_filename(brick):
        TRACTOR_FILENAME = 'tractor/tractor-%(brickid)d.fits'
        return TRACTOR_FILENAME % dict(brickid=brick.id) 

    @staticmethod
    def parse_filename(filename, brickindex):
        if not filename.endswith('.fits'): raise ValueError
        return brickindex.get_brick(
            brickindex.search_by_id(
                int(re.search('-([0123456789]+)\.', 
                os.path.basename(filename)).group(1))))

class EDR3:
    BRICKS_FILENAME = 'decals-bricks.fits'
    CATALOGUE_ALIASES = []
    
    @staticmethod
    def format_image_filenames():
        images = {
        'depth': 'coadd/%(pre)s/%(brickname)s/decals-%(brickname)s-depth-%(band)s.fits',
        'model': 'coadd/%(pre)s/%(brickname)s/decals-%(brickname)s-model-%(band)s.fits',
        'image': 'coadd/%(pre)s/%(brickname)s/decals-%(brickname)s-image-%(band)s.fits',
        }
        imagerepos= {}
        for image in images:
            imagerepos[image] = {}
            for band in 'rgz':
                PATTERN = images[image]
                def getfilename(brick):
                    return PATTERN % dict(pre=brick.name[:3], brickname=brick.name, band=band)
                imagerepos[image][band] = getfilename
        return imagerepos

    @staticmethod
    def format_catalogue_filename(brick):
        TRACTOR_FILENAME = 'tractor/%(pre)s/tractor-%(brickname)s.fits'
        return TRACTOR_FILENAME % dict(pre=brick.name[:3], brickname=brick.name) 

    @staticmethod
    def parse_filename(filename, brickindex):
        if not filename.endswith('.fits'): raise ValueError
        brickname = re.search('-([p0123456789]+)\.', 
                os.path.basename(filename)).group(1)
        bid = brickindex.search_by_name(brickname)
        return brickindex.get_brick(bid)

_configurations = {
    'EDR': EDR,
    'EDR3': EDR3
}
class DataRelease(object):
    """
    The highest level interface into the data for a given imaging
    data release.  Uses several "helper" classes and has methods for
    looking at pixelized data or catalogs arranged in bricks.
    """
    def __init__(self, root=None, cacheroot=None, version=None):
        if root is None:
            root = os.environ.get("DECALS_IMAGING", '.') 
        root = os.path.normpath(root)
        self.root = root

        if cacheroot is None:
            cacheroot = os.path.normpath(os.environ.get("DECALS_CACHE", '.'))

        if version is None:
            version = os.path.basename(root).upper()

        config = _configurations[version]


        self.cacheroot = os.path.join(cacheroot, version)

        brickdata = fits.read_table(os.path.join(self.root, config.BRICKS_FILENAME))

        self.brickindex = brickindex.BrickIndex(brickdata)

        self.bands = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'Y':5}

        self.images = {}
        image_filenames = config.format_image_filenames()
        for image in image_filenames:
            self.images[image] = {}
            for band in image_filenames[image]:
                self.images[image][band] = imagerepo.ImageRepo(self.root, image_filenames[image][band])

        self.observed_bricks = [ ]
        for roots, dirnames, filenames in \
            os.walk(os.path.join(self.root, 'tractor'), followlinks=True):
            for filename in filenames:
                try:
                    self.observed_bricks.append(
                        config.parse_filename(filename, self.brickindex))
                except ValueError:
                    pass 

        self._observed_brickids = self.brickindex.search_by_id(
            [ brick.id for brick in self.observed_bricks ])
        # the list of observed bricks must be sorted.
        arg = self._observed_brickids.argsort()
        self._observed_brickids = self._observed_brickids[arg]
        self.observed_bricks = numpy.array(self.observed_bricks)[arg]

        # approximate area in degrees. Currently a brick is 0.25 * 0.25 deg**2
        self.observed_area = 41253. * len(self.observed_bricks) / len(self.brickindex)

        self.catalogue = catalogue.Catalogue(
            cachedir=os.path.join(self.cacheroot, 'catalogue'),
            filenames=[
                os.path.join(self.root, 
                config.format_catalogue_filename(brick))
                for brick in self.observed_bricks],
            aliases=config.CATALOGUE_ALIASES
            )

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

        mask = contains(self._observed_brickids, bid)

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
