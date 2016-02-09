Scripts
=======

ImagingLSS provides several scripts that producing complete targeting catalogues and
the corresponding sky masks. In this section we will look at the usage of these files.

Building Catalogue Cache
------------------------

build_cache.py

Caching makes interactive analysis faster.  
We convert the small catalogue (tractor) files into several large files, one for each
column in the table.

build_cache.py must be run with MPI. Here is an example script that works on Edison.

.. code-block::

    #PBS -j eo
    #PBS -l mppwidth=256
    #PBS -q debug

    cd $PBS_O_WORKDIR

    export OMP_NUM_THREADS=1
    export DECALS_PY_CONFIG=/project/projectdirs/m779/imaginglss/dr1/dr1.conf.py

    export PYTHON_MPI_CHROOT=/dev/shm
    export PYTHON_MPI_PKGROOT=/project/projectdirs/m779/python-mpi/usg
    export PYTHON_MPI_PACKAGES=matplotlib-1.4.3.tar.gz:mpi4py-1.3.1.tar.gz:numpy-1.9.2.tar.gz:python-2.7.9.tar.gz:scipy-0.15.1.tar.gz:yfeng1-local-05132015.tar.gz
    :fitsio-0.9.7a.tar.gz

    aprun -n 256 -d 1 /project/projectdirs/m779/python-mpi/python-mpi \
        ./build_cache.py


    
Complete Target Selection
-------------------------

select_objs.py selects objects of a type. The target types are defined at http://desi.lbl.gov/trac/wiki/TargetSelection

The output is the RA DEC and magnitudes of objects. There is an option to use the Tycho star catalogue to veto objects
near a known star.

Complete Random Sky Mask
------------------------

make_random.py generates the randoms for the sky mask of an target type.

The output is the RA DEC and magnitudes limit at that location on the sky. 
There is an option to use the Tycho star catalogue to veto regions near a known star.

