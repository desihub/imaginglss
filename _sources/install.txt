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

For example, on Edison.

.. code-block:: bash

    export PYTHONUSERBASE=$SCRATCH

    easy_install --user fitsio
    #easy_install --user astropy


Location of Data Release
------------------------
 
ImagingLSS tries to initialize the DECALS data release from a configuration file, by
either passing in the file name as an argument to :py:class:`imaginglss.DECALS` 
or by setting the environment variable :code:`DECALS_PY_CONFIG`.

Here is an example configuration file:

.. code-block:: python

    # dr1j.conf.py

    decals_root = "/scratch1/scratchdirs/desiproc/dr1j"
    decals_cache = "/global/scratch2/sd/yfeng1/desicache/dr1j")
    decals_release = "DR1J"
    dust_dir = "/project/projectdirs/desi/software/edison/dust/v0_0/"

.. code-block:: bash

    export DECALS_PY_CONFIG=/path/to/dr1j.conf.py


Getting Started
---------------

Before start using imaginglss, it is crucial to build the catalogue cache. 

Building the cache takes a long time, and we provide a script :code:`scripts/build_cache.py`. 

This is an MPI capable script with good weak scaling and only bound by the file system bandwidth; 
We recommend running it with a modest number of cores, 128 or 256.

Python starts slowly on Edison, depending on the condition of the meta data server of the file systems. 
We are narrowing down to a permanant solution to this.

.. code-block:: bash

    # in a job script or an interactive job.
    # replace mpirun with aprun on Cray

    mpirun -n 256 python-mpi scripts/build_cache.py

For DR1 this takes about 10 mins and produces 300 to 400 GB of data in the
cache directory. Be ready for it. 

Now we are ready for some interactive sessions

.. code-block:: bash

    from imaginglss import DECALS
    decals = DECALS() # or DECALS('/path/to/dr1.conf.py')

    dr = decals.datarelease
    cat = decals.datarelease.catalogue
    sfd = decals.sfdmap

Also refer to the :doc:`examples` of the documentation.

Example Dataset
---------------

We provide a small sampling data set that contains only 1 brick for testing imaginglss. 
The small dataset is based on the DR1J data release of DECALS. It can be downloaded at 

http://imaginglss.s3-website-us-west-1.amazonaws.com/dr1j-mini.tar.gz 

The total size is less than 30 MB after decompressing. 

Note that the SFD98 dust map is also required for most of the ImagingLSS operations.
The SFD98 file is somewhat larger, on the order of 100 MB.

http://imaginglss.s3-website-us-west-1.amazonaws.com/SFD98.tar.gz 

To deploy this dataset with the source code tree, 
see the following steps.

.. code-block:: bash

    mkdir testdata
    cd testdata

    wget http://imaginglss.s3-website-us-west-1.amazonaws.com/dr1j-mini.tar.gz
    tar -xzvf dr1j-mini.tar.gz

    wget http://imaginglss.s3-website-us-west-1.amazonaws.com/SFD98.tar.gz
    tar -xzvf SFD98.tar.gz

    export DECALS_PY_CONFIG=$PWD/dr1j-mini/dr1j-mini.conf.py

    cd -
    

Replace dr1j-mini with edr3 to get the EDR3 dataset. 
(
http://imaginglss.s3-website-us-west-1.amazonaws.com/edr3.tar.gz 
)
EDR3 is older than DR1J, but with more bricks than the mini dataset.

DR1 on Edison
-------------

ImagingLSS has been prepackaged for DR1 at Edison in the following locations:

.. code-block:: bash

    export DECALS_PY_CONFIG=/global/project/projectdirs/m779/imaginglss/dr1/dr1.conf.py

Most of code development is happening at that location as well.

.. code-block:: bash

    tree DECALS_PY_CONFIG=/global/project/projectdirs/m779/imaginglss/source

