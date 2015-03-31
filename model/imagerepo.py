from utils import fits
import os.path

class ImageRepo(object):
    """
    Image repository.
    This class serves as an interface to the file system, the
    layout of files and IO for brick files.
    Given a Brick object, ImageRepo can return the image data
    or meta data of the image.
    Standard users should not need to modify this class.
    """
    def __init__(self, root, pattern):
        """
        pattern is a python formatting string.
        supported keywords are

        pattern can also be a callable with signature:

        def pattern(brick):
            return filename_without_root

        %(brickid)d
        %(brickname)d

        We need to extend this for DR1.
        """
        self.root = root
        self.pattern = pattern
        self.cache = {}
        self.meta_cache = {} 

    def preload(self, bricks, **kwargs):
        """ preload the images for 
            Brick objects listed in bricks list.

            **kwargs is unused 
        """
        for b in bricks:
            self.cache[b] = self.open(b, **kwargs) 
            
    def open(self, brick, **kwargs):
        """ open and read the image for 
            Brick object brick. The image content is returned
            as an array.

            **kwargs is unused 
        """
        if brick in self.cache:
            return self.cache[brick]
        fname = self.get_filename(brick, **kwargs)
        if not os.path.exists(fname):
            raise IOError('%s does not exist' % fname)
        _ = self.metadata(brick, **kwargs)
        return fits.read_image(fname)
        
    def metadata(self, brick, **kwargs):
        """ Fetch the meta data about a brick.
        """
        if brick not in self.meta_cache:
            meta = fits.read_metadata(self.get_filename(brick, **kwargs))
            self.meta_cache[brick] = meta
        return self.meta_cache[brick]

    def get_filename(self, brick, **kwargs):
        if hasattr(self.pattern, '__call__'):
            fn = os.path.join(self.root, 
                self.pattern(brick, **kwargs))
        else:
            kwargs['brickid'] = brick.id
            kwargs['brickname'] = brick.name
            fn = os.path.join(self.root, 
                self.pattern) % kwargs
        if not os.path.exists(fn):
            fngz = fn + '.gz' 
            if os.path.exists(fngz):
                return fngz
        return fn

