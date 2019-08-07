.. sectnum::

.. contents::

=====================================================
 General notes about file and directory organization
=====================================================

- Various checklists (useful for developers) are scattered throughout the repository, in
  the location that is most closely related to the given process. To aid discovery of
  these checklists, these have the naming convention ``README.CHECKLIST.*``.

- Unit tests are stored in the ``test`` subdirectory of each ``src`` directory. For
  example, unit tests of biogeophysics modules reside in
  ``src/biogeophysics/test``. Within each test subdirectory, there are further
  subdirectories for each set of unit tests. These typically correspond to a single module
  (e.g., ``src/biogeophysics/test/Irrigation_test`` contains tests of
  ``src/biogeophysics/IrrigationMod.F90``), but that is not a strict requirement.
