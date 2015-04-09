Description of model used to process DECALS imaging data.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The DECALS imaging data are processed in "Brick"s on the sky,
each Brick being a small, contiguous regions of sky (roughly
0.5x0.5 degrees).  The data, including information on depths,
extinction, etc. are stored in pixelized images within the
FITS files.  There are also catalogs of objects, produced by
Tractor, which are associated with each brick.

The main way to access the DECALS imaging data is via DataRelease
objects defined in "datarelease.py". Important attributes are

 - DataRelease.footprint : the survey footprint. 
 - DataRelease.brickindex: the Brick geometry (BrickIndex)
 - DataRelease.images    : the images released in the DR
 - DataRelease.catalogue : the imaging catalogue reported by tractor

The member methods are

 - DataRelease.readout   : reading out pixel value of an image 
                           for coordinates on the sky.

Catalog information is handled by catalog.Catalogue, which stores the
object catalogs associated with a data release.
The catalogs are contained in many small FITS files, 
but this class caches the
information for speed and only columns that are accessed are loaded
into memory. Use 'forget_cache.py' to clear this cache.


The images are represented by imagerepo.ImageRepo objects. 
An ImageRepo object takes care of reading the image tile of a Brick. 
An image tile contains a padded region around the Brick, 
which is sometimes mentioned as a margin, a padding, or a bleeding region.

A DataRelease has several image repositories, for example,

 - DataRelease.images['depth']['u'] : the 'u' band inverse variance depth images
 - DataRelease.images['image']['u'] : the coadded 'u' band image in nanomaggies
 - DataRelease.images['model']['u'] : the model 'u' band image in nanomaggies; 

or what tractor thinks the sky looks like.

brick.Brick holds metadata for each Brick, including the central 'ra',
'dec', and the primary area covered by the brick, excluding the padding regsions. There are routines for converting from sky coordinates to image
pixels within each brick, and routines for querying the images
handled by imagerepo. 

brickindex.BrickIndex holds the metadata
of of the brick decomposition scheme. 
There are routines for converting from
sky coordinates to bricks, and acts as a factory that creates Brick
objects.

Helper routines. In order to reduce the number of dependecies, 
we implemented a few coordinate transformations 

 - wcs_tangent impments the WCS tangent format.
 - wcs_simplezea implement the WCS ZEA format, oriented at north / south pole

