Installation
===========

This guide covers different methods to install HistoMapTx, with Conda being the recommended approach.

Conda Installation (Recommended)
-------------------------------

The recommended way to install HistoMapTx is using Conda, which helps manage dependencies efficiently:

1. First, make sure you have Conda installed (either `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ or `Anaconda <https://www.anaconda.com/products/distribution>`_)

2. Create a new environment for HistoMapTx:

   .. code-block:: bash

      conda create -n histomap python=3.10
      conda activate histomap

3. Install required dependencies:

   .. code-block:: bash

      conda install -c conda-forge geopandas pandas numpy matplotlib plotly scipy shapely
      conda install -c conda-forge jupyterlab  # for running tutorial notebooks

4. Install HistoMapTx using pip:

   .. code-block:: bash

      pip install histomaptx

Pip Installation
--------------

If you prefer using pip, you can install HistoMapTx directly:

.. code-block:: bash

   pip install histomaptx

Note that you'll need to ensure all dependencies are properly installed, which may be more challenging than using Conda.


Troubleshooting
-------------

Common installation issues and their solutions:

GeoPandas Installation Issues
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you encounter issues installing GeoPandas, try:

.. code-block:: bash

   conda install -c conda-forge geopandas

This ensures all GeoPandas dependencies are correctly installed.

Platform-Specific Notes
^^^^^^^^^^^^^^^^^^^^^

**Windows Users**:
   Shapely and other geospatial libraries might require additional steps. Using Conda is strongly recommended.

**Mac M1/M2 Users**:
   Make sure to use the arm64 version of Conda for best performance.

Verifying Installation
-------------------

To verify HistoMapTx is correctly installed, run:

.. code-block:: python

   import histomap
   print(histomap.__version__)

If this runs without errors, your installation is successful.