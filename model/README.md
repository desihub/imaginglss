# Description of model used to process DECALS imaging data.

The DECALS imaging data are processed in "bricks" on the sky,
each brick being a small, contiguous regions of sky (roughly
0.5x0.5 degrees).  The data, including information on depths,
extinction, etc. are stored in pixelized images within the
FITS files.  There are also catalogs of objects, produced by
Tractor, which are associated with each brick.

Our data model has the management of the individual FITS files
handled by "imagerepo.py" and the imagerepo class.

The brick class, defined in "brick.py", holds metadata for each
brick, routines for converting from sky coordinates to image
pixels within each brick, and routines for querying the images
handled by imagerepo.

The brickindex class, defined in "brick_index.py", contains...

Catalog information is handled by "catalog.py", which stores the
list of files to be requested by

The images contain only coordinate transformations in the WCS
"tangent" format.  In order to speed up these transformations, and
to make the code more stand-alone, we have implemented the pixel
to sky and sky to pixel routines in "wcs_tangent.py" as helper
routines.
