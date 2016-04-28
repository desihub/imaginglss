Installation
============

ImagingLSS can be installed as a regular Python package, but it can also be
used without installation (since the code evolves rapidly with time).

To install (this is how we test ImagingLSS at home)

.. code-block:: bash

    git clone https://github.com/desihub/imaginglss/
    cd imaginglss
    pip install --user -e .

To use without installation (this is how we test ImagingLSS at NERSC)

.. code-block:: bash

    git clone https://github.com/desihub/imaginglss/
    cd imaginglss
    export PYTHONPATH=$PWD:$PYTHONPATH
    # remember to set PYTHONPATH in new shells.

The package only declares the dependency for the data model.

However, some basic packages need to be set up before using the data pipeline.
We also need to tell ImagingLSS the location of files. These post-installation
steps are documented in the following section.

Depencency
----------

ImagingLSS is light on dependency.

 - numpy 1.8+. 
   Do not use 1.7.x due to this bug: https://github.com/numpy/numpy/issues/3469

   To check the version of numpy, use

   .. code-block:: bash
    
       python -c 'import numpy; print numpy.version.version'

 - mpi4py, for the parallel infrastructure. 

 - h5py for accessing HDF5 files. HDF5 is the internal data storage format of the
   data pipeline. The final data product is plain text.

 - fitsio 0.9.8 + for accessing FITS files.
   Do not use 0.9.7 due this bug: https://github.com/esheldon/fitsio/issues/39
   which affects building catalogue cache. 
   For now, you can install the rc2 prerelease of fitsio.

   .. code-block:: bash

       pip install -U --no-deps --user https://github.com/esheldon/fitsio/archive/v0.9.8rc2.tar.gz
   
   We still do support astropy.io.fits, but the use of astropy.io.fits is not
   recommended.
 

On parallel HPC systems where files are hosted in a shared file system, 
the time it takes to launch a Python application may fluctuate badly. 
This applies to ImagingLSS, too. 
At NERSC, we set up python-mpi-bcast to eliminate this issue https://github.com/rainwoodman/python-mpi-bcast .
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
either passing in the file name as an argument to :py:class:`imaginglss.DECALS`, 
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

ImagingLSS has been prepackaged for DR2 at Edison in the following locations.

These commands will work in JupyterHub: https://jupyter.nersc.gov .

.. code-block:: python

    from imaginglss import DECALS
    decals = DECALS('/global/project/projectdirs/m779/imaginglss/dr2.conf.py')

    dr = decals.datarelease
    cat = decals.datarelease.catalogue


Example Dataset
---------------

For those who do **not** work on NERSC, 
we provide a small sampling data set that contains a few bricks from DR2. 

http://imaginglss.s3-website-us-west-1.amazonaws.com/dr2-mini.tar.gz 

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

    wget http://imaginglss.s3-website-us-west-1.amazonaws.com/dr2-mini.tar.gz
    tar -xzvf dr2-mini.tar.gz

    wget http://imaginglss.s3-website-us-west-1.amazonaws.com/SFD98.tar.gz
    tar -xzvf SFD98.tar.gz

    wget http://imaginglss.s3-website-us-west-1.amazonaws.com/tycho.tar.gz
    tar -xzvf tycho.tar.gz

    cd -

