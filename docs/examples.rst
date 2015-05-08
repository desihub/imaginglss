Examples
========

.. contents:: 
    :depth: 2

Notebooks
---------

ImagingLSS is compatible with IPython Notebook. 

However, the notebook must manipulate the python import path at the begining.
This limitation is because imaginglss is intended to be used without installation. 

Also note that querying the imaging pixels (e.g. depth) typically requires a 
high-throughput operation that cannot currently be driven in a Notebook; 
For example, to apply ELG selection for the full DR1, we implement a MPI python scripts 
:py:code:`script/select_elgs.py`.

We provide a comprensive Notebook example at

http://nbviewer.ipython.org/urls/imaginglss-git.s3.amazonaws.com/NotebookDemo.ipynb

Here is an example notebook investigating selection of ELGs based on color / mag cuts
and completeness cuts.

http://nbviewer.ipython.org/urls/imaginglss-git.s3.amazonaws.com/SelectELG.ipynb

Here is an example notebook investigating the complete area mask for ELGs.

http://nbviewer.ipython.org/urls/imaginglss-git.s3.amazonaws.com/RandomMaskELG.ipynb

Example Job script on Edison
----------------------------

.. literalinclude:: ../scripts/edison_example.sh

