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

    Note that the part of the LSS generation pipe line (files in scripts/imglss-mpi-*.py ) 
    typically requires a high-throughput operation querying a large amount of data, 
    these cannot be done with notebook. The notebooks are most useful for inspecting
    smaller chunks of data, individual bricks, etc. Refer to :doc:`datapipeline`.

At NERSC, ImagingLSS can be used with the Jupyter Hub service at https://jupyter.nersc.gov .
Of course, one need to properly set up imaginglss for the service.
Refer to :doc:`install` and https://github.com/bccp/imaginglss-notebooks/blob/master/NERSCJupyterGuide.ipynb

