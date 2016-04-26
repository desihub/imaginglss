Usage: Data Pipeline
====================

The primary usage of ImagingLSS is to produce large scale structure catalogues.
In this section, we document the data pipeline. There are the steps:

- Building the cache
- Generating the object catalogue
- Generating the random sampling of the mask
- Applying veto based on proximity to stars
- Estimating completeness for randoms and objects
- Assemble and export final text files

The first 3 steps are computational/data extensive and we provide mpi4py based
scripts; which shall go through the job queue of a computing facility.

The rest of the steps are very light weight and the head nodes can easily handle them.

We will cover them in the following sections. The examples are based off our
working configuration at NERSC.

Building the Cache
------------------
Before we start using ImagingLSS, it is crucial to build the catalogue cache. 

.. attention:: 

    On Edison for DR2, this step can be omitted 
    if the following configuration file is used:

    .. code-block:: bash
       	/project/projectdirs/m779/imaginglss/dr2.conf.py

The application to build the cache is :code:`scripts/imglss-mpi-build-cache.py`. The application
scans the files in the data release described by a configuration file (provided in
:code:`--conf` argument), and converts the columns used by ImagingLSS to a binary 
format that is much easier to use into the 'cache' directory specified in the configuration
file. For the format of the configuration file, refer to `Configuration File`.

The inline help of the script describes the usage:

.. code-block:: bash

    usage: imglss-mpi-build-cache.py [-h] [--conf CONF]

    optional arguments:
      -h, --help   show this help message and exit
      --conf CONF  Path to the imaginglss config file, default is from
                   DECALS_PY_CONFIG

Here is an example job script that works on Edison (Note that python-mpi-bcast is used). 
Submit the job script with :code:`sbatch`.

.. code-block:: bash

    #!/bin/bash
    #SBATCH -J imglss-mpi-build-cache
    #SBATCH -n 512
    #SBATCH -o imglss-mpi-build-cache.%j
    #SBATCH -p debug
    #SBATCH -t 00:30:00

    export OMP_NUM_THREADS=1

    module load python/2.7-anaconda
    source /project/projectdirs/m779/python-mpi/nersc/activate.sh

    # change the following line to where your ImagingLSS is installed
    mirror ~/source/imaginglss imaginglss scripts

    # change conf to your imaginglss configuration file
    srun -n 256 python-mpi /dev/shm/local/scripts/imglss-mpi-build-cache.py --conf /project/projectdirs/m779/imaginglss/dr2.conf.py

Generating Object Catalogue
---------------------------

:code:`scripts/imglss-mpi-select-objects.py` selects objects of a type, and writes out the objects
that we will use in the later stages of the pipeline.

The target types are defined at http://desi.lbl.gov/trac/wiki/TargetSelection

The output is the RA DEC and magnitudes of objects. 

The inline help of the script describes the usage:

.. code-block:: bash

    usage: imglss-mpi-select-objects.py [-h] 
                          [--conf CONF]
                          {QSO,LRG,ELG,BGS} output

    positional arguments:
      {QSO,LRG,ELG,BGS}
      output

    optional arguments:
      -h, --help            show this help message and exit
      --conf CONF           Path to the imaginglss config file, default is from
                            DECALS_PY_CONFIG


Here is an example job script we use on Edison to generate the LRG catalogue.
Submit the job script with :code:`sbatch`. We also encourage typing in the commands
one by one from an interactive job session, obtained via :code:`salloc`. Refer to
`http://www.nersc.gov/users/computational-systems/cori/running-jobs/interactive-jobs/`_.


.. code-block:: bash

    #!/bin/bash

    #SBATCH -J imglss-mpi-select-objects
    #SBATCH -n 512
    #SBATCH -o imglss-mpi-select-objects.%j
    #SBATCH -p debug
    #SBATCH -t 00:30:00

    export OMP_NUM_THREADS=1

    module load python/2.7-anaconda
    source /project/projectdirs/m779/python-mpi/nersc/activate.sh

    # change the following line to where your imaginglss is installed
    mirror ~/source/imaginglss imaginglss scripts

    # use without installing
    export PYTHONPATH=/dev/shm/local:$PYTHONPATH

    # change conf to your imaginglss configuration file
    srun -n 256 python-mpi /dev/shm/local/scripts/imglss-mpi-select-objects.py LRG LRG.hdf5 --conf /project/projectdirs/m779/imaginglss/dr2.conf.py


Generating Complete Random Sky Mask
-----------------------------------

imglss-mpi-make-random.py generates the randoms for the sky mask. The points will be uniform within the survey footprint.

The inline help of the script describes the usage:

.. code-block:: bash

    usage: imglss-mpi-make-random.py [-h] 
                          [--conf CONF]
                          Nran output

    positional arguments:
      Nran                  Minimum number of randoms
      output

    optional arguments:
      -h, --help            show this help message and exit
      --conf CONF           Path to the imaginglss config file, default is from
                            DECALS_PY_CONFIG


Here is an example job script we use on Edison to generate a QSO random catalogue.
Submit the job script with :code:`sbatch`. We also encourage typing in the commands
one by one from an interactive job session, obtained via :code:`salloc`. Refer to
`http://www.nersc.gov/users/computational-systems/cori/running-jobs/interactive-jobs/`_.

.. code:: 

    #!/bin/bash

    #SBATCH -J imglss-mpi-make-random
    #SBATCH -n 512
    #SBATCH -o imglss-mpi-make-random.%j
    #SBATCH -p debug
    #SBATCH -t 00:30:00

    export OMP_NUM_THREADS=1

    module load python/2.7-anaconda
    source /project/projectdirs/m779/python-mpi/nersc/activate.sh

    # change the following line to where your imaginglss is installed
    mirror ~/source/imaginglss imaginglss scripts

    # use without installing
    export PYTHONPATH=/dev/shm/local:$PYTHONPATH

    # change conf to your imaginglss configuration file
    srun -n 256 python-mpi /dev/shm/local/scripts/imglss-mpi-make-random.py 6000000 QSO-random.hdf5 --conf /project/projectdirs/m779/imaginglss/dr2.conf.py

Apply Star veto mask
--------------------

imglss-query-tycho-veto.py applies the bright star veto masks to a target or random catalogue. The veto types are defined
in imaginglss/analysis/tycho_veto.py . As you can tell, we currently only support vetoing via a Tycho2 catalogue.

The star veto mask is important for correctly building the completeness estimator.

The inline help of the script describes the usage:

.. code::

    usage: imglss-query-tycho-veto.py [-h] [--conf CONF] catalogue

    Query the TYCHOVETO flags of input data. The position is taken from the NOISES
    extension of input. The result is written to the TYCHOVETO extension of
    output. Currently, only veto by proximity to tycho stars are implemented. Each
    veto in imaginglss.analysis.tycho_veto is calculated as a column in the
    TYCHOVETO extension. Unfortunately, this script is not sufficiently smart to
    decide the correct TYCHOVETO for the target type. Therefore, no combined veto
    flag is generated.

    positional arguments:
      catalogue    HDF5 catalogue file, can be either random or objects.
                   TYCHO_VETO dataset will be added

    optional arguments:
      -h, --help   show this help message and exit
      --conf CONF  Path to the imaginglss config file, default is from
                   DECALS_PY_CONFIG


Query Completeness
------------------

imglss-query-completeness.py esitmates the fractional completeness for objects / randoms based on their depth.
A threshold confidence level is used to generate a 100% complete sample based on an object catalogue. Then
this sample is taken to model the fractional completeness. The result is appended as COMPLETENESS column to the 
catalogue.

The inline help of the script describes the usage:

.. code::

    Usage: imglss-query-completeness.py [-h]
                                        [--use-tycho-veto {BOSS_DR9,DECAM_BGS,DECAM_ELG,DECAM_LRG,DECAM_QSO}]
                                        [--sigma-z SIGMA_Z] [--sigma-g SIGMA_G]
                                        [--sigma-r SIGMA_R] [--conf CONF]
                                        {MYBGS,ELG,QSOC,LRG,QSO,QSOd,BGS} objects
                                        query

    positional arguments:
      {MYBGS,ELG,QSOC,LRG,QSO,QSOd,BGS}
      objects               object catalogue for building the completeness model.
      query                 catalogue to query completeness

    optional arguments:
      -h, --help            show this help message and exit
      --use-tycho-veto {BOSS_DR9,DECAM_BGS,DECAM_ELG,DECAM_LRG,DECAM_QSO}
      --sigma-z SIGMA_Z
      --sigma-g SIGMA_G
      --sigma-r SIGMA_R
      --conf CONF           Path to the imaginglss config file, default is from
                            DECALS_PY_CONFIG


Assemble Final Product
----------------------

imglss-export-text.py assembles a final catalogue for objects or randoms. The final product is a plain text file.
fluxes (only for objects) and depths of selected bands can be included in the final product.

We need a threshold confidence level (usually identical to the one used in imglss-query-completenesss) to filter
out poorly detected objects.

Vetoing by proximity to stars is also applied at this final stage.

The inline help of the script describes the usage:

.. code::

    usage: imglss-export-text.py [-h] [--conf CONF]
                             [--use-tycho-veto {BOSS_DR9,DECAM_BGS,DECAM_ELG,DECAM_LRG,DECAM_QSO}]
                             [--bands {Y,W4,r,u,W1,g,i,W3,z,W2} [{Y,W4,r,u,W1,g,i,W3,z,W2} ...]]
                             [--sigma-z SIGMA_Z] [--sigma-g SIGMA_G]
                             [--sigma-r SIGMA_R]
                             catalogue output

    positional arguments:
      catalogue             internal catalogue of HDF5 type.
      output                text file to store the catalogue.

    optional arguments:
      -h, --help            show this help message and exit
      --conf CONF           Path to the imaginglss config file, default is from
                            DECALS_PY_CONFIG
      --use-tycho-veto {BOSS_DR9,DECAM_BGS,DECAM_ELG,DECAM_LRG,DECAM_QSO}
      --bands {Y,W4,r,u,W1,g,i,W3,z,W2} [{Y,W4,r,u,W1,g,i,W3,z,W2} ...]
      --sigma-z SIGMA_Z
      --sigma-g SIGMA_G
      --sigma-r SIGMA_R

