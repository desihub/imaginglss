Install
=======

To use this imaing data model on edison:

1. set environments
  
.. code-block:: bash

    export DECALS_IMAGING=/global/project/projectdirs/cosmo/work/decam/release/edr/
    export DECALS_CACHE=$GSCRATCH/desicache
 
2. install astropy or fitsio (prefered)

.. code-block:: bash

    easy_install --user fitsio

or

.. code-block:: bash  

    easy_install --user astropy

3. try if usecase.py runs.

.. code-block:: bash

    python usecase.py

4. start hacking with

.. code-block:: bash

    from model.datarelease import DataRelease
    dr = DataRelease()
