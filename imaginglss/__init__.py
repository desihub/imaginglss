from .model.datarelease import DataRelease
from .model.sfdmap import SFDMap
from .model.tycho import Tycho

import os.path

class DECALS(object):
    """ DECALS configuration object

        Load DECALS configuration from a configuration file.
        The configuration file describes where and how to
        access the data release.

        The main objects for access this data release are created
        as the configuration object is parsed.

        If :code:`filename` is None, the file pointed by environment
        variable :code:`DECALS_PY_CONFIG` is used.

        Attributes
        ----------
        datarelease : :py:class:`imaginglss.model.datarelease.DataRelease`
            The datarelease object described by this configuration file

        sfdmap : :py:class:`imaginglss.model.sfdmap.SFDMap`
            The dust extinction object described by this configuration file.

        Configuation File Format
        ------------------------

        A configuration file is parsed as a python script.

        Predefined variables are:

        - :py:mod:`os.path` 
            
        - :code:`__file__` : the filename of the configration script 

        - :code:`__path__` : the dirname part of :code:`__file__`

        The file must set the following variables:

        - decals_root : the root path of the DECALS data release; 
          for example :code:`/scratch1/scratchdirs/desiproc/dr1j`

        - decals_cache : the path of the DECALS cache directory. 
          This directory shall be writable unless the cache has
          been complete / frozen.

        - decals_release : the version of the DECALS data release.
          It shall be a supported string that points to a class
          in :py:mod:`imaginglss.schema`.
        
        - dust_dir  : the location to look for dust map. for example
          /project/projectdirs/desi/software/edison/dust/v0_0/

    """
    def __init__(self, filename=None):
        if filename is None:
            filename = os.environ.get('DECALS_PY_CONFIG')

        d = {}
        d['decals_root'] = os.environ.get("DECALS_IMAGING", '.') 
        d['decals_cache'] = os.environ.get("DECALS_CACHE", '.') 
        d['dust_dir'] = os.environ.get("DUST_DIR", '.') 
        d['tycho_dir'] = os.environ.get("TYCHO_DIR", '.') 
        d['decals_release'] = os.path.basename(d['decals_root']).upper()

        if filename:
            # if a configuration file is specified.

            d['__file__'] = os.path.abspath(filename)
            d['os'] = os
            d['__path__'] = os.path.abspath(os.path.dirname(filename))

            with open(filename, 'r') as ff:
                script = ff.read()
                exec(script, d)

        self.decals_root = d['decals_root']
        self.decals_release = d['decals_release']
        self.cache_dir = d['decals_cache']
        self.dust_dir = d['dust_dir']
        self.tycho_dir = d['tycho_dir']

        self.datarelease = DataRelease(root=self.decals_root, 
                cache=self.cache_dir, version=self.decals_release,
                dustdir=self.dust_dir)
        self.sfdmap = self.datarelease.sfdmap
        self.tycho = Tycho(self.tycho_dir)
