Installation
============

ImagingLSS can be installed as a regular Python package. 

Since the code still evolves, the way to install ImagingLSS is

.. code-block:: bash

    git clone https://github.com/desihub/imaginglss/
    cd imaginglss
    pip install --user -e .

The package only declares the dependency for the data model.

However, some basic packages need to be set up before using the data pipeline.
We also need to notify imaginglss the location of files. These post installation
steps are documented in this section.

Depencency
----------

ImagingLSS is light on dependency.

 - numpy 1.8+. 
   Do not use 1.7.x for this bug: https://github.com/numpy/numpy/issues/3469

   To check the version of numpy, use

   .. code-block:: bash
    
       python -c 'import numpy; print numpy.version.version'

 - mpi4py, for the parallel infrastructure. 

 - fitsio 0.9.8 + for accessing FITS files.
   Do not use 0.9.7 for this bug: https://github.com/esheldon/fitsio/issues/39
   which affects building catalogue cache. 
   For now, we can install rc2 prerelease of fitsio.

   .. code-block:: bash

       pip install -U --no-deps --user https://github.com/esheldon/fitsio/archive/v0.9.8rc2.tar.gz
   
   We still do support astropy.io.fits, but the use of astropy.io.fits is not
   recommended.
 

On parallel HPC systems where files are hosted in a shared file system, 
the time it takes to launch a Python application may fluctuate badly. 
This applied to ImagingLSS, too. 
At NERSC, we setup python-mpi-bcast to eliminate this issue https://github.com/rainwoodman/python-mpi-bcast .
This will be noted in the full examples.

Data Dependency
---------------

Besides the DECALS catalogue and coadd image data, ImagingLSS also needs 

 - SFD98 dust map [todo give a reference]
 - Tycho star catalogue [http://slac.stanford.edu/~erykoff/catalogs/tycho.fit]

We provide mirrors of these files at
    
 - http://imaginglss.s3-website-us-west-1.amazonaws.com/SFD98.tar.gz
 - http://imaginglss.s3-website-us-west-1.amazonaws.com/tycho.tar.gz

Location of Data Release
------------------------
 
ImagingLSS tries to initialize the DECALS data release from a configuration file, by
either passing in the file name as an argument to :py:class:`imaginglss.DECALS` 
or by setting the environment variable :code:`DECALS_PY_CONFIG`.

Here is an example configuration file (that works on Edison):

.. code-block:: python

    # dr2.conf.py
    decals_root = "/global/project/projectdirs/cosmo/data/legacysurvey/dr2"
    decals_cache = "/project/projectdirs/m779/imaginglss/dr2/cache"
    decals_release = "DR2"
    dust_dir = "/project/projectdirs/desi/software/edison/dust/v0_0/"
    tycho_dir = "/project/projectdirs/m779/imaginglss/tycho.fit"
    
DR2 at NERSC
------------

.. attention::

   FIXME: Ellie please check the path is correct.

ImagingLSS has been prepackaged for DR2 at Edison in the following locations:

.. code-block:: bash

    /global/project/projectdirs/m779/imaginglss/dr2.conf.py


Example Dataset
---------------

For those who do not work on NERSC, 
we provide a small sampling data set that covers the EDR3 foot-print.

http://imaginglss.s3-website-us-west-1.amazonaws.com/dr1j-edr3.tar.gz 

The total size is less than 45 MB after decompressing. 

The SFD98 dust map is required for target selelection and completeness masks
The SFD98 file is somewhat larger, on the order of 100 MB.

http://imaginglss.s3-website-us-west-1.amazonaws.com/SFD98.tar.gz 

The Tycho star catalogue is required for target selelection and completeness masks.

http://imaginglss.s3-website-us-west-1.amazonaws.com/tycho.tar.gz 

To deploy this dataset with the source code tree, 
see the following steps.

.. code-block:: bash

    mkdir testdata
    cd testdata

    wget http://imaginglss.s3-website-us-west-1.amazonaws.com/dr1j-mini.tar.gz
    tar -xzvf dr1j-edr3.tar.gz

    wget http://imaginglss.s3-website-us-west-1.amazonaws.com/SFD98.tar.gz
    tar -xzvf SFD98.tar.gz

    wget http://imaginglss.s3-website-us-west-1.amazonaws.com/tycho.tar.gz
    tar -xzvf tycho.tar.gz

    cd -

