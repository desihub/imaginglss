Working with Jupyter Notebook 
=============================

.. attention::

   FIXME: We need to update these notebooks.

ImagingLSS is compatible with IPython Notebook (after ImagingLSS is installed).

We provide a comprensive Notebook example at

http://nbviewer.ipython.org/urls/imaginglss-git.s3.amazonaws.com/NotebookDemo.ipynb

Here is an example notebook investigating selection of ELGs based on color / mag cuts
and completeness cuts.

http://nbviewer.ipython.org/urls/imaginglss-git.s3.amazonaws.com/SelectELG.ipynb

Here is an example notebook investigating the complete area mask for ELGs.

http://nbviewer.ipython.org/urls/imaginglss-git.s3.amazonaws.com/RandomMaskELG.ipynb

.. attention:: 

    Note that the LSS generation pipe line (select_objs.py, make_randoms.py) 
    typically requires a high-throughput operation querying a large amount of data, 
    these cannot be done with notebook. The notebooks are most useful for inspecting
    smaller chunks of data, individual bricks, etc. Refer to :doc:`datapipeline`.

At NERSC, the jupyterhub service ImagingLSS can be used with https://jupyter.nersc.gov .
Of course, one need to install imaginglss in from a Jupyter Terminal session first via pip.
Refer to :doc:`install`.
