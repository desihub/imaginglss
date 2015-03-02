from astropy.io import fits
import os.path

class ImageRepo(object):
    """ Image repository.
        Fetching image data from a brick

        pattern is "image-%(brick)s-r.fits.gz"
    """
    def __init__(self, root, pattern):
        self.root = root
        self.pattern = pattern

    def open(self, brick):
        return fits.open(self.get_filename(brick))

    def get_filename(self, brick):
        return os.path.join(self.root, 
            self.pattern) % brick.name)[0].data[:]

class ImageRepoEDR(ImageRepo):
    """ EDR image repo. 
        Note that naming convention has changed from 
        EDR to EDR3

        in EDR it is %(brick)d % brick.id
    """
    def get_filename(self, brick):
        return os.path.join(self.root, 
            self.pattern) % brick.id)[0].data[:]

