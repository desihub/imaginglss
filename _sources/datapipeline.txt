Usage: Data Pipeline
====================

The primary usage of ImagingLSS is to produce large scale structure catalogues.
In this section, we document the data pipeline. There are three steps:

- Building the cache
- Generating the object catalogue
- Generating the random sampling of the mask

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

The application to build the cache is :code:`scripts/build_cache.py`. The application
scans the files in the data release described by a configuration file (provided in
:code:`--conf` argument), and converts the columns used by ImagingLSS to a binary 
format that is much easier to use into the 'cache' directory specified in the configuration
file. For the format of the configuration file, refer to `Configuration File`.

The inline help of the script describes the usage:

.. code-block:: bash

    usage: build_cache.py [-h] [--conf CONF]

    optional arguments:
      -h, --help   show this help message and exit
      --conf CONF  Path to the imaginglss config file, default is from
                   DECALS_PY_CONFIG

Here is an example job script that works on Edison (Note that python-mpi-bcast is used). 
Submit the job script with :code:`sbatch`.

.. code-block:: bash

    #!/bin/bash
    #SBATCH -J build_cache
    #SBATCH -n 512
    #SBATCH -o build_cache.%j
    #SBATCH -p debug
    #SBATCH -t 00:30:00

    export OMP_NUM_THREADS=1

    module load python/2.7-anaconda
    source /project/projectdirs/m779/python-mpi/nersc/activate.sh

    # change the following line to where your ImagingLSS is installed
    mirror ~/source/imaginglss imaginglss scripts

    # change conf to your imaginglss configuration file
    srun -n 256 python-mpi /dev/shm/local/scripts/build_cache.py --conf /project/projectdirs/m779/imaginglss/dr2.conf.py
    
Generating Object Catalogue
---------------------------

:code:`scripts/select_objs.py` selects objects of a type, and writes out the objects
that we will use in the LSS catalogue.

The target types are defined at http://desi.lbl.gov/trac/wiki/TargetSelection

The output is the RA DEC and magnitudes of objects. 
There is an option to use the Tycho star catalogue to veto objects near a known star.

The inline help of the script describes the usage:

.. code-block:: bash

    usage: select_objs.py [-h] [--sigma-z SIGMA_Z] [--sigma-g SIGMA_G]
                          [--sigma-r SIGMA_R]
                          [--with-tycho {BOSS_DR9,DECAM_BGS,DECAM_ELG,DECAM_LRG,DECAM_QSO}]
                          [--conf CONF]
                          {QSO,LRG,ELG,BGS} output

    positional arguments:
      {QSO,LRG,ELG,BGS}
      output

    optional arguments:
      -h, --help            show this help message and exit
      --sigma-z SIGMA_Z
      --sigma-g SIGMA_G
      --sigma-r SIGMA_R
      --with-tycho {BOSS_DR9,DECAM_BGS,DECAM_ELG,DECAM_LRG,DECAM_QSO}
                            Type of veto.
      --conf CONF           Path to the imaginglss config file, default is from
                            DECALS_PY_CONFIG


Here is an example job script we use on Edison to generate the LRG catalogue.
Submit the job script with :code:`sbatch`. We also encourage typing in the commands
one by one from an interactive job session, obtained via :code:`salloc`. Refer to
`http://www.nersc.gov/users/computational-systems/cori/running-jobs/interactive-jobs/`_.


.. code-block:: bash

    #!/bin/bash

    #SBATCH -J select_objs
    #SBATCH -n 512
    #SBATCH -o select_objs.%j
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
    srun -n 256 python-mpi /dev/shm/local/scripts/select_objs.py LRG LRG.txt --with-tycho DECAM_LRG --conf /project/projectdirs/m779/imaginglss/dr2.conf.py


Generating Complete Random Sky Mask
-----------------------------------

make_random.py generates the randoms for the sky mask of a target type.

The output is the RA DEC and magnitudes limit at that location on the sky. 
There is an option to use the Tycho star catalogue to veto regions near a known star.

The inline help of the script describes the usage:

.. code-block:: bash

    usage: make_random.py [-h] [--sigma-z SIGMA_Z] [--sigma-g SIGMA_G]
                          [--sigma-r SIGMA_R]
                          [--with-tycho {BOSS_DR9,DECAM_BGS,DECAM_ELG,DECAM_LRG,DECAM_QSO}]
                          [--conf CONF]
                          Nran {QSO,LRG,ELG,BGS} output

    positional arguments:
      Nran                  Minimum number of randoms
      {QSO,LRG,ELG,BGS}
      output

    optional arguments:
      -h, --help            show this help message and exit
      --sigma-z SIGMA_Z
      --sigma-g SIGMA_G
      --sigma-r SIGMA_R
      --with-tycho {BOSS_DR9,DECAM_BGS,DECAM_ELG,DECAM_LRG,DECAM_QSO}
                            Type of veto.
      --conf CONF           Path to the imaginglss config file, default is from
                            DECALS_PY_CONFIG


Here is an example job script we use on Edison to generate a QSO random catalogue.
Submit the job script with :code:`sbatch`. We also encourage typing in the commands
one by one from an interactive job session, obtained via :code:`salloc`. Refer to
`http://www.nersc.gov/users/computational-systems/cori/running-jobs/interactive-jobs/`_.

.. code:: 

    #!/bin/bash

    #SBATCH -J make_random
    #SBATCH -n 512
    #SBATCH -o make_random.%j
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
    srun -n 256 python-mpi /dev/shm/local/scripts/make_random.py 6000000 QSO QSO_rand.txt --with-tycho DECAM_QSO --conf /project/projectdirs/m779/imaginglss/dr2.conf.py
