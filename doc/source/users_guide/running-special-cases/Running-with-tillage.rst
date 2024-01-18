.. _running-with-tillage:

.. include:: ../substitutions.rst

=====================
 Running with tillage
=====================


Cropland tillage (Sect. :numref:`decomp_mgmt_modifiers`) can be toggled by specifying a value of ``'low'`` (low intensity) or ``'high'`` (high intensity) for the  ``tillage_mode`` namelist option. By default this option is ``'off'``.

Depth of tillage can be changed with the ``max_tillage_depth`` parameter (meters; default 0.26).

Tillage multipliers for different soil pools and time since planting are defined on the parameter file, in variables ``bgc_till_decompk_multipliers`` (for CENTURY soil) and ``mimics_till_decompk_multipliers`` (for MIMICS soil). These variables were originally added with the script at ``tools/contrib/add_tillage_to_paramsfile.py``, which can be modified as needed to change tillage multipliers.


Example: Crop simulation with tillage
-------------------------------------
::

   > cime/scripts/create_newcase -case IHistClm51BgcCrop_till -res f19_g17_gl4 -compset IHistClm51BgcCrop


   > cd IHistClm51BgcCrop_till
   > ./case.setup

   # turn on tillage ('low' or 'high'; default 'off')
   > echo "tillage_mode = 'high'" >> user_nl_clm

Reverting fixes relative to original tillage implementation
-----------------------------------------------------------

The current implementation of tillage in CTSM is based on work by :ref:`Graham et al. (2021) <Grahametal2021>`, but with fixes to how days after planting is calculated and to default tillage depth. To run without those changes:

::

   > echo "use_original_tillage_phases = .true." >> user_nl_clm
   > echo "max_tillage_depth = 0.32" >> user_nl_clm
