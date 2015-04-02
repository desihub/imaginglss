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
    """ returns mask:

        len(mask) == len(need), 
        mask[i] == true only if needle[i] is in haystack;
    
        haystack must be sorted.
    """

    ind = haystack.searchsorted(needle)
    ind.clip(0, len(haystack) - 1, ind)
    return haystack[ind] == needle

class EDR:
    BRICKS_FILENAME = 'bricks.fits'
    CATALOGUE_ALIASES = [('EXTINCTION', 'DECAM_MW_TRANSMISSION', lambda x: 10**(x/-2.5))]

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
                def getfilename(brick, PATTERN=PATTERN, band=band):
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
        return brickindex.search_by_id(
                int(re.search('-([0123456789]+)\.', 
                os.path.basename(filename)).group(1)))

class EDR3:
    BRICKS_FILENAME = 'decals-bricks.fits'
    CATALOGUE_ALIASES = [('DECAM_EXTINCTION', 'DECAM_MW_TRANSMISSION', lambda x: 10**(x/-2.5))]
    
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
                def getfilename(brick, PATTERN=PATTERN, band=band):
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
        return bid

class EDR4(EDR3):
    CATALOGUE_ALIASES = []
    pass

class DR1(EDR4):
    pass

_configurations = {
    'EDR': EDR,
    'EDR3': EDR3,
    'EDR4': EDR4,
    'DR1': DR1,
}
class Footprint(object):
    """ footprint of a data release.
        
        Attributes
        ----------

        bricks : list, model.brick.Brick
            covered bricks
        range :  tuple
            (ramin, ramax, decmin, decmax)
        area  : float
            covered outline area in square degrees
    """
    def __init__(self, dr):
        self.bricks = [dr.brickindex.get_brick(bid) for bid in dr._covered_brickids]
        self.area = 41253. * len(dr._covered_brickids) / len(dr.brickindex)

        # range of ra dec of covered bricks
        FootPrintRange = namedtuple('FootPrintRange', ['ramin', 'ramax', 'decmin', 'decmax'])
        self.range = FootPrintRange(
            ramin=min([brick.ra1 for brick in self.bricks]),
            ramax=max([brick.ra2 for brick in self.bricks]),
            decmin=min([brick.dec1 for brick in self.bricks]),
            decmax=max([brick.dec2 for brick in self.bricks]),)

        self._covered_brickids = dr._covered_brickids
        self.brickindex = dr.brickindex

    def __repr__(self):
        return "Footprint: len(bricks)=%d , area=%g degrees, range=%s" % (
                len(self.bricks),
                self.area,
                str(self.range)
            )
    def filter(self, coord):
        """ filter coord, remove those are not covered by the current footprint 

            returns coord_in_footprint
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

    brickindex : model.brickindex.BrickIndex
        an index object of all of the bricks (covering the entire sky)
    bands      : dict
        a dictionary translating from band name to integer used in Tractor catalogue
    catalogue  : model.catalogue.Catalogue
        the concatenated tractor catalogue, accessed by attributes.
    extinction : array_like
        an array stroing the extinction coeffcients
    images :      model.imagerepo.ImageRepo
        image repositories (depth, image, model are by bands, ebv is the extinction map)
    footprint  : Footprint 
        the footprint of the data release

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

        try:
            brickdata = fits.read_table(os.path.join(self.root, config.BRICKS_FILENAME))
        except :
            brickdata = fits.read_table(os.path.join(os.path.dirname(__file__), '..', 'fallback', 'default-bricks.fits'))

        self.brickindex = brickindex.BrickIndex(brickdata)

        self.bands = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'Y':5}

        # E(B-V) to ugrizY bands, SFD98; used in tractor
        self.extinction = numpy.array([3.995, 3.214, 2.165, 1.592, 1.211, 1.064], dtype='f8')\
            .view(dtype=[(band, 'f8') for band in 'ugrizY'])[0]

        self.images = {}
        image_filenames = config.format_image_filenames()
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
        self.images['ebv'] = imagerepo.ImageRepo(self.cacheroot, image_filename)
            
        _covered_brickids = [ ]
        for roots, dirnames, filenames in \
            os.walk(os.path.join(self.root, 'tractor'), followlinks=True):
            for filename in filenames:
                try:
                    _covered_brickids.append(
                        config.parse_filename(filename, self.brickindex))
                except ValueError:
                    pass 
        
        self._covered_brickids = numpy.array(_covered_brickids, dtype='i8')
        # the list of covered bricks must be sorted.
        self._covered_brickids.sort()

        self.footprint = Footprint(self) # build the footprint property

        self.catalogue = catalogue.Catalogue(
            cachedir=os.path.join(self.cacheroot, 'catalogue'),
            filenames=[
                os.path.join(self.root, 
                config.format_catalogue_filename(brick))
                for brick in self.footprint.bricks],
            aliases=config.CATALOGUE_ALIASES
            )

    def readout(self, coord, repo, default=numpy.nan, ignore_missing=False):
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
                pixels[sl] = img
            except IOError:
                if not ignore_missing:
                    raise
        #
        images[mask] = pixels[invarg]
            
        return images
