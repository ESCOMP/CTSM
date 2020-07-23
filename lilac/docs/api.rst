LILAC API
=========

LILAC provides a high-level API (in Fortran) for coupling to CTSM.
LILAC was built using the assumption that the Atmosphere model
component in each LILAC application would opperate as the "driver".

The atmosphere component will need to call each of the following subroutines:

  * `lilac_init`
  * `lilac_run`
  * `lilac_final`

LILAC Core
----------
.. f:autosrcfile:: core.f90

ESMF Utils
----------
.. f:autosrcfile:: esmf_utils.f90
