.. _running-with-tillage:

.. include:: ../substitutions.rst

=====================
 Running with tillage
=====================


Cropland tillage (Sect. :numref:`decomp_mgmt_modifiers`) is set to ``'low'`` by default. This can be changed to a value of ``'off'`` (no tillage) or ``'high'`` (high-intensity tillage) for the  ``tillage_mode`` namelist option.

Depth of tillage can be changed with the ``max_tillage_depth`` parameter (meters; default 0.26).

Tillage multipliers for different soil pools and time since planting are defined on the parameter file, in variables ``bgc_till_decompk_multipliers`` (for CENTURY soil) and ``mimics_till_decompk_multipliers`` (for MIMICS soil). These variables were originally added with the script at ``tools/contrib/add_tillage_to_paramsfile.py``, which can be modified as needed to change tillage multipliers.


Example: Crop simulation with no tillage
----------------------------------------
::

   > cime/scripts/create_newcase -case IHistClm60BgcCrop_notill -res f19_g17_gl4 -compset IHistClm60BgcCrop


   > cd IHistClm60BgcCrop_notill
   > ./case.setup

   # turn off tillage
   > echo "tillage_mode = 'off'" >> user_nl_clm

Reverting fixes relative to original tillage implementation
-----------------------------------------------------------

The current implementation of tillage in CTSM is based on work by :ref:`Graham et al. (2021) <Grahametal2021>`, but with fixes to how days after planting is calculated and to default tillage depth. To run without those changes:

::

   > echo "use_original_tillage_phases = .true." >> user_nl_clm
   > echo "max_tillage_depth = 0.32" >> user_nl_clm
