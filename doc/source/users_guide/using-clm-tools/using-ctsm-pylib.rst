
.. include:: ../substitutions.rst
.. _using-ctsm-pylib:

Installing the CTSM Python environment
======================================

Many of our Python-based tools require non-standard Python modules to be installed. To facilitate this, you can install a CTSM-specific `Conda <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html>`__ environment with the exact version of everything needed. From the top level of your CTSM checkout, simply do

.. code:: shell

   ./py_env_create

and that script will install the ``ctsm_pylib`` environment for you. If ``ctsm_pylib`` already exists, it will give you options on how to handle that.
