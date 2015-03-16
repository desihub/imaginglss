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

The member methods are:
 - DataRelease.readout   : reading out pixel value of an image 
                for coordinates on the sky.

It holds the imaging catalogue 
(Catalogue), a BrickIndex for the sky decomposition, and routines
to query within the full survey area image values at any give sky
coordinate. 

The images are represented by imagerepo.ImageRepo objects. A DataRelease
has several image repositories. The ImageRepo object takes care of reading
the image tile of a Brick. An image tile contains a padded region around the
Brick, which is sometimes mentioned as a margin, a padding, or a bleeding region.

 - DataRelease.images['depth']['u'] : the 'u' band inverse variance depth images
 - DataRelease.images['image']['u'] : the coadded 'u' band image in nanomaggies
 - DataRelease.images['model']['u'] : the model 'u' band image in nanomaggies; 
or what tractor thinks the sky looks like.

The Brick class, defined in "brick.py", holds metadata for each
brick, routines for converting from sky coordinates to image
pixels within each brick, and routines for querying the images
handled by imagerepo. 

The BrickIndex class, defined in "brick_index.py", holds the metadata
of of the brick decomposition scheme, rountines for converting from
sky coordinates to bricks, and acts as a factory that creates Brick
objects.

Catalog information is handled by "catalog.py", which stores the
object catalogs associated with a data release.
The catalogs are contained in FITS files, but this class caches the
information for speed and only columns that are accessed are loaded
into memory.

The images contain only coordinate transformations in the WCS
"tangent" format.  In order to speed up these transformations, and
to make the code more stand-alone, we have implemented the pixel
to sky and sky to pixel routines in "wcs_tangent.py" as helper
routines.
