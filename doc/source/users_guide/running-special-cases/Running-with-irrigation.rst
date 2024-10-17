.. include:: ../substitutions.rst

.. _running-with-irrigation:

===================================
 Running with irrigation
===================================

.. todo::
    Remove refs to pre-5.0 behavior?

In CLM4.0 irrigation isn't an allowed option. In CLM4.5 irrigation can ONLY be used WITH crop. With CLM5.0 irrigation can be used whether crop is on or not -- **BUT** if crop is off, your surface datasets **HAVE** to have irrigation defined appropriately. Right now *ALL* surface datasets without crop enabled have irrigation hard-wired on. In order to create datasets with irrigation off, you'd need to make changes to ``mksurfdata_esmf`` in order to have all generic crops to be non-irrigated. To turn on irrigation in |version| we simply add ``-irrig on`` to ``CLM_BLDNML_OPTS``.

Example: Irrigation Simulation
------------------------------------------
::

   # Note here we do a CLMSP simulation as that is what has been validated
   > cd cime/scripts
   > ./create_newcase -case IRRIG -res f19_g17_gl4 -compset I1850Clm50BgcCrop
   > cd IRRIG

   # Append "-irrig on" to CLM_BLDNML_OPTS in env_run.xml (you could also use an editor)
   > ./xmlchange CLM_BLDNML_OPTS="-irrig on" -append

   > ./case.setup

   # Now build and run normally
   > ./case.build
   > ./case.submit

