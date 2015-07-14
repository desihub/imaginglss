Installation
============

ImagingLSS is supposed to be used without installation.

However, some basic packages need to be set up before using the scripts

Depencency
----------

ImagingLSS is light on dependency.

 - numpy 1.8+. 
   Do not use 1.7.x for this bug: https://github.com/numpy/numpy/issues/3469

   .. code-block:: bash
    
       python -c 'import numpy; print numpy.version.version'

 - mpi4py, for the parallel infrastructure. 

 - fitsio 0.9.8 + (astropy is supported too) for accessing FITS files.
   Do not use 0.9.7 for this bug: https://github.com/esheldon/fitsio/issues/39
   which affects building catalogue cache. 
   Since 0.9.8 has not yet been released, you may need to install this
   this version of fitsio which contains this crucial fix: 
   https://github.com/esheldon/fitsio/tree/150cd
   
We do recomment pointing the variable PYTHONUSERBASE to a fast filesystem before installing
the dependencies of ImagingLSS. This shall help the speed and consistency of the python startup time.
A more comprehensive solution to faster PYTHON can be found at https://github.com/rainwoodman/python-mpi-bcast.

For example, on Edison.

.. code-block:: bash

    export PYTHONUSERBASE=$SCRATCH

    easy_install --user fitsio

Data Dependency
---------------

Besides the DECAM catalogue and coadd image data, ImagingLSS also needs 

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

    # dr1.conf.py
    decals_root = "/global/project/projectdirs/cosmo/data/legacysurvey/dr1"
    decals_cache = "/project/projectdirs/m779/imaginglss/dr1/cache"
    decals_release = "DR1"
    dust_dir = "/project/projectdirs/desi/software/edison/dust/v0_0/"
    tycho_dir = "/project/projectdirs/m779/imaginglss/tycho.fit"

.. code-block:: bash

    export DECALS_PY_CONFIG=/path/to/dr1j.conf.py

Getting Started
---------------

Before start using imaginglss, it is crucial to build the catalogue cache. 

Building the cache takes a long time, and we provide a MPI script for this
:code:`scripts/build_cache.py`. On Edison for DR1, this step can be omitted
if the above example configuration file is used.

build_cache.py is an MPI capable script with good weak scaling. 
It is bound by the file system bandwidth; 

We recommend running it with a modest number of cores, 128 or 256. Here is
an example job script that works on Edison (Note that python-mpi-bcast is used):

.. code-block:: bash

    #PBS -j eo
    #PBS -l mppwidth=256
    #PBS -q debug

    cd $PBS_O_WORKDIR

    export OMP_NUM_THREADS=1
    export DECALS_PY_CONFIG=/project/projectdirs/m779/imaginglss/dr1/dr1.conf.py

    export PYTHON_MPI_CHROOT=/dev/shm
    export PYTHON_MPI_PKGROOT=/project/projectdirs/m779/python-mpi/usg
    export PYTHON_MPI_PACKAGES=matplotlib-1.4.3.tar.gz:mpi4py-1.3.1.tar.gz:numpy-1.9.2.tar.gz:python-2.7.9.tar.gz:scipy-0.15.1.tar.gz:yfeng1-local-05132015.tar.gz:fitsio-0.9.7a.tar.gz

    aprun -n 256 -d 1 /project/projectdirs/m779/python-mpi/python-mpi \
        ./build_cache.py

For DR1 this takes about 10 mins and produces 300 to 400 GB of data in the
cache directory.  Be ready for it. 

Now we are ready for some interactive fun

.. code-block:: bash

    from imaginglss import DECALS
    decals = DECALS() # or DECALS('/path/to/dr1.conf.py')

    dr = decals.datarelease
    cat = decals.datarelease.catalogue
    print cat['RA'][:100]
    print cat['DEC'][:100]

Example Dataset
---------------

We provide a small sampling data set that is covers the EDR3 foot-print.

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

    export DECALS_PY_CONFIG=$PWD/dr1j-edr3/dr1j.conf.py

    cd -
    
DR1 on Edison
-------------

ImagingLSS has been prepackaged for DR1 at Edison in the following locations:

.. code-block:: bash

    export DECALS_PY_CONFIG=/global/project/projectdirs/m779/imaginglss/dr1/dr1.conf.py

The working copy of the source code can be accessed at that location as well.

.. code-block:: bash

    tree DECALS_PY_CONFIG=/global/project/projectdirs/m779/imaginglss/source

