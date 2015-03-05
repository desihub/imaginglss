import os
import os.path
from astropy.io import fits
import numpy
import glob
import re

from . import brickindex
from . import imagerepo
from . import catalogue

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
    BANDS = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'Y':5}

    BRICKS_FILENAME = 'bricks.fits'
    DEPTH_FILENAME = 'coadd/depth-%(brickid)d-%(band)s.fits.gz'
    IMAGE_FILENAME = 'coadd/image-%(brickid)d-%(band)s.fits' 
    MODEL_FILENAME = 'coadd/model-%(brickid)d-%(band)s.fits'
    TRACTOR_FILENAME = lambda brick: \
                'tractor/tractor-%(brickid)d.fits' \
                % dict(brickid=brick.id) 

    @staticmethod
    def getimagefilename(PATTERN, band):
        def getfilename(brick):
            return PATTERN % dict(brickid=brick.id, band=band)
        return getfilename

    def __init__(self, root=None, cacheroot=None, version=None):
        if root is None:
            root = os.environ.get("DECALS_IMAGING", '.') 
        self.root = root

        if cacheroot is None:
            cacheroot = os.environ.get("DECALS_CACHE", '.') 

        self.cacheroot = os.path.join(cacheroot, version)

        bricks = fits.open(os.path.join(self.root, self.BRICKS_FILENAME))[1].data

        self.bricks = numpy.array(bricks, copy=True)
        self.brickindex = brickindex.BrickIndex(self.bricks)

        self.bands = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'Y':5}

        self.images = {}
        for band in self.bands:
            self.images[band] = dict(
                DEPTH=imagerepo.ImageRepo(self.root, self.getimagefilename(self.DEPTH_FILENAME, band)),
                IMAGE=imagerepo.ImageRepo(self.root, self.getimagefilename(self.IMAGE_FILENAME, band)),
                MODEL=imagerepo.ImageRepo(self.root, self.getimagefilename(self.MODEL_FILENAME, band)),
            )

        self.build_catalogue()

    def build_catalogue(self):
        observed_brickids = numpy.unique([
            int(re.search('-([0123456789]+)\.', fn).group(1))
                for fn in glob.glob(os.path.join(self.root, 'tractor/tractor-[0-9]*.fits'))
        ])
        self.observed_brickids = observed_brickids

        # approximate area in degrees. Currently a brick is 0.25 * 0.25 deg**2
        self.observed_area = 41253. * len(self.observed_brickids) / len(self.bricks)

        self.observed_bricks = self.bricks[self.bricks['BRICKID'].searchsorted(self.observed_brickids)]
        self.catalogue = catalogue.Catalogue(
            os.path.join(self.cacheroot, 'catalogue'),
            [
            os.path.join(self.root, 'tractor/tractor-%d.fits' % brick)
            for brick in self.observed_brickids])

         
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
        mask = contains(self.observed_brickids, bid + 1)
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
