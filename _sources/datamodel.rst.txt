Data Model
==========

The DECALS imaging data are processed in "Brick"s on the sky,
each Brick being a small, contiguous regions of sky (roughly
0.5x0.5 degrees).  The data, including information on depths,
extinction, etc. are stored in pixelized images within the
FITS files.  There are also catalogs of objects, produced by
Tractor, which are associated with each brick.

The main way to access the DECALS imaging data is via DECALS
objects defined in :py:class:`imaginglss.DECALS`. 

This is an example using the data model in an interactive Python shell.

.. code-block:: bash

    from imaginglss import DECALS
    decals = DECALS('/global/project/projectdirs/m779/imaginglss/dr2.conf.py')

    dr = decals.datarelease
    cat = decals.datarelease.catalogue

The DECALS object constructs the DECALS data release and dust extinction
map objects from environment variables or a configuration file. Please refer
to the documentation of :py:class:`imaginglss.DECALS` for details.

The :py:attr:`imaginglss.DECALS.datarelease` attribute is represents the data
release. It is an instance of
:py:class:`imaginglss.model.datarelease.DataRelease`.
Important attributes are

 - :py:attr:`~imaginglss.model.datarelease.DataRelease.footprint` : the survey footprint. 
 - :py:attr:`~imaginglss.model.datarelease.DataRelease.brickindex`: the Brick geometry (BrickIndex)
 - :py:attr:`~imaginglss.model.datarelease.DataRelease.images`    : the images released in the DR
 - :py:attr:`~imaginglss.model.datarelease.DataRelease.catalogue` : the imaging catalogue reported by tractor

The member methods are

 - :py:meth:`~imaginglss.model.datarelease.DataRelease.readout`   : reading out pixel value of an image 
   for coordinates on the sky.
 - :py:meth:`~imaginglss.model.datarelease.DataRelease.create_footprint`   : Creating a footprint for
   a given extent on the sky. The footprint will be rounded to bricks. 
 - :py:meth:`~imaginglss.model.datarelease.DataRelease.create_catalogue`   : Creating a catalogue for
   a given footprint. Currently this is not cached and can be slow.

The full catalog information is handled by :py:class:`imaginglss.model.catalogue.CachedCatalogue`, which stores the
object catalogs associated with a data release.
The catalogs are contained in many small FITS files, 
but this class caches the
information for speed and only columns that are accessed are loaded
into memory. Usually at a site where DECALS data is deployed, the cache only needs to be generated once.


The images are represented by :py:class:`imaginglss.model.imagerepo.ImageRepo` objects. 
An ImageRepo object takes care of reading the image tile of a Brick. 
An image tile contains a padded region around the Brick, 
which is sometimes mentioned as a margin, a padding, or a bleeding region.

A DataRelease has several image repositories, for example,

 - DataRelease.images['depth']['u'] : the 'u' band inverse variance depth images
 - DataRelease.images['image']['u'] : the coadded 'u' band image in nanomaggies
 - DataRelease.images['model']['u'] : the model 'u' band image in nanomaggies; 

:py:class:`imaginglss.model.brick.Brick` holds metadata for each Brick, including the central 
:py:attr:`~imaginglss.model.brick.Brick.ra`, :py:attr:`~imaginglss.model.brick.Brick.dec`, 
and the primary area(
:py:attr:`~imaginglss.model.brick.Brick.ra1`, :py:attr:`~imaginglss.model.brick.Brick.dec1`, 
:py:attr:`~imaginglss.model.brick.Brick.ra2`, :py:attr:`~imaginglss.model.brick.Brick.dec2`
) covered by the brick, excluding the padding region. 
There are routines for converting from sky coordinates to image
pixels within each brick, and routines for querying the images
handled by imagerepo. 

:py:class:`imaginglss.model.brickindex.BrickIndex` holds the metadata
of of the brick decomposition scheme. 
There are routines for converting from
sky coordinates to bricks, and acts as a factory that creates 
:py:class:`~imaginglss.model.brick.Brick` objects.

The :py:attr:`imaginglss.DECALS.sfdmap` attribute is represents the dust extinction map
. It is an instance of :py:class:`imaginglss.model.sfdmap.SFDMap`.

For additional information please refer to :doc:`modules`.

