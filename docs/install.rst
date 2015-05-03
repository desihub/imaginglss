Installation
============

ImagingLSS is supposed to be used without installation.

However, some basic packages need to be set up before using the scripts

Depencency
----------

ImagingLSS is light on dependency.
It depends only on fitsio or astropy:

.. code-block:: bash

    easy_install --user fitsio
    #easy_install --user astropy

On Edison,

.. code-block:: bash

    module load fitsio

File Paths
----------
 
ImagingLSS core libary and the tools scripts defaults to use :code:`DECALS_IMAGING` as the root path of the data release.
In addition, :code:`DECALS_CACHE` must be set to a directory with large free storage space, for caching the catalogue data.
The dust extinction map shall be given in :code:`DUST_DIR` variable.

If a directory is not specified at the creation of an :py:class:`imaginglss.DataRelease` or :py:class:`imaginglss.SFDMap`,
the values from the environment variables are used.

.. code-block:: bash

    export DECALS_IMAGING=/global/project/projectdirs/cosmo/work/decam/release/edr/
    export DECALS_CACHE=$GSCRATCH/desicache
    export DUST_DIR=/project/projectdirs/desi/software/edison/dust/v0_0/
 
Getting Started
---------------

.. code-block:: bash

    from model.datarelease import DataRelease
    dr = DataRelease()

Also refer to the :doc:`examples` of the documentation.

Example Dataset
---------------

We provide a small sampling data set for testing imaginglss. The small dataset is
based on the DR1J data release of DECALS. It can be downloaded at 

http://imaginglss.s3-website-us-west-1.amazonaws.com/dr1j.tar.gz

The total size is less than 200 MB after decompressing. To deploy this dataset with
the source code tree, in imainginglss toplevel directory (that is the directory 
with :code:`.git` directory)

.. code-block:: bash

    wget http://imaginglss.s3-website-us-west-1.amazonaws.com/dr1j.tar.gz
    tar -xzvf dr1j.tar.gz

    export DECALS_IMAGING=$PWD/dr1j
    export DECALS_CACHE=$PWD/dr1j/cache
    export DUST_DIR=$PWD/dr1j
    




