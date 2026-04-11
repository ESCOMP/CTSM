.. _rst_mizuRoute:

River Routing Model (MizuRoute)
====================================================

.. _Overview mizuRoute:

Overview
---------

MizuRoute is a river transport model designed for applications across local, regional and global scales that can use either regular grids or more intricate unstructured grids on Hydrologic Response Units (HRU's).
When run with regular grids those grids still need to be formatted to an unstructured grid format for use by mizuRoute, as mizuRoute currently only handles data over land and not data over ocean.
The name corresponds to the Japanese word for water, "mizu", and the English word "route," reflecting its purpose in modeling water flow through river systems.

MizuRoute is a significant advancement beyond the MOSART model used in CLM50. A few notable mizuRoute features include:

#. Ability to run on HRU's allows for more sophisticated hydrology applications.

#. Hybrid (MPI Distributed Memory + OpenMP Shared Memory) parallelization.

#. A lake model which enables the simulations of lake volumes and flow from natural and managed lakes.

#. MizuRoute like MOSART can also be run connected to the Hillslope option in CTSM.

Furhter reading
-------------------------

For more information about mizuRoute, please read `mizuRoute Technical note`_.

.. _mizuRoute Technical note: https://mizuroute.readthedocs.io/en/latest/tech_note/index.html

