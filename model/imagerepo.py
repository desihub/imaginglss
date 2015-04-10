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
        Initizlize a ImageRepo.

        Parameters
        ----------
        pattern: callable
            the function is called with brick object to generate
            a file name. String is supported for backward compatibility.
        root: string
            root path that is concatenated to the pattern.

        """
        self.root = root
        self.pattern = pattern
        self.cache = {}
        self.meta_cache = {} 

    def preload(self, bricks, **kwargs):
        """ Preload images into the cache
    
            This function loads images for given bricks into the cache.
            In generate the function is not useful, and we shall
            try to remove it soon.

            Parameters
            ----------
            bricks: list of Brick
                the bricks whose images will be loaded

        """
        for b in bricks:
            self.cache[b] = self.open(b, **kwargs) 
            
    def open(self, brick, **kwargs):
        """ Open and read an image.

            The image for Brick object brick is read into memory and returned.
            If the image is already in cache, return the cache.
            
            This function does not add the image to cache.

            Parameters
            ----------
            brick: Brick
                the brick whose image will be loaded
            
            Returns
            -------
            image: array_like
                the image as an ndarray, as the ordering in FITS files.
                (index with (y, x))
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

            Parameters
            ----------
            brick: Brick
                the brick whose metadata will be loaded
            
            Returns
            -------
            metadata: dict
                the metadata as a dictionary. Currently this is the
                FITS header
        """
        if brick not in self.meta_cache:
            meta = fits.read_metadata(self.get_filename(brick, **kwargs))
            self.meta_cache[brick] = meta
        return self.meta_cache[brick]

    def get_filename(self, brick, **kwargs):
        """ Generate a filename.

            Generate a filename for a brick. 

            Notes
            -----
            We also try to generate
            a .gz filename, if the original filename does not exist.

            Parameters
            ----------
            brick: Brick
                the brick whose filename will be generated.
             
        """
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

