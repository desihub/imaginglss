""" The schema of DataReleases.

    DataRelease Objects manages the migrating schema of EDR, DR. 
    For example, if a field is renamed in DR from EDR, we want to modify
    the EDR schema such that the code using DataRelease does not have to
    be modified.

    The `version` parameter of `DataRelease.__init__` picks the correct
    schema. 
"""
import re
import os.path

class Schema:
    pass

class EDR(Schema):
    """ Schema for EDR.

    """
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

class EDR3(Schema):
    """ Schema for EDR3.

    """
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
        brickname = re.search('-([pm0123456789]+)\.', 
                os.path.basename(filename)).group(1)
        bid = brickindex.search_by_name(brickname)
        return bid

class EDR4(EDR3):
    """ Schema for EDR4.
        Changed to DECAM_MW_TRANSMISSION

    """
    CATALOGUE_ALIASES = []
    pass

class DR1(EDR4):
    """ Schema for DR1. Same as EDR4.

    """
    pass

class DR1J(EDR4):
    """ Schema for DR1J. Same as EDR4.

    """
    pass

class DR2P(EDR4):
    """ Schema for DR2P. Added new images.

    """
    @staticmethod
    def format_image_filenames():
        images = {
        'depth': 'coadd/%(pre)s/%(brickname)s/decals-%(brickname)s-depth-%(band)s.fits',
        'galdepth': 'coadd/%(pre)s/%(brickname)s/decals-%(brickname)s-galdepth-%(band)s.fits',
        'nexp': 'coadd/%(pre)s/%(brickname)s/decals-%(brickname)s-nexp-%(band)s.fits',
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

    pass
