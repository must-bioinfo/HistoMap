Installation
===========

This guide covers how to install HistoMapTx. We recommend using a Python virtual environment (venv).

Venv Installation (Recommended)
--------------------------------

1. Create and activate a virtual environment:

   .. code-block:: bash

      python3 -m venv histomap_env
      source histomap_env/bin/activate  # On Windows: histomap_env\Scripts\activate

2. Install HistoMapTx:

   .. code-block:: bash

      pip install histomaptx

Conda
-----

.. warning::

   Installing HistoMapTx in a Conda environment may cause issues with ``pkg_resources``
   (a dependency of ``xarray_schema``, itself a transitive dependency of ``spatialdata``).
   If you encounter a ``ModuleNotFoundError: No module named 'pkg_resources'`` error,
   use a venv instead.

