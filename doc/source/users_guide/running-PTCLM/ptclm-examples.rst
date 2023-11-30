.. include:: ../substitutions.rst

.. _ptclm-examples:

==============================
 Examples of using PTCLMmkdata
==============================

Now let's give a few more complex examples using some of the options we have discussed above.

Now, let's demonstrate using a different group list, doing a spinup, running with Qian global forcing data, but using tower years to set the years to run over. This uses the options: sitegroupname, useQIAN, and QIANtower_years.

Example: Running PTCLMmkdata without tower years
------------------------------------------------
::

   > cd $CTSMROOT/tools/PTCLM
   > ./PTCLMmkdata -s US-Ha1 -d $CSMDATA --sitegroupname AmeriFlux --donot_use_tower_yrs
   > cd ../../../../../US-Ha1_ICRUCLM45BGC_QIAN
   # Now build and run normally
   ```

Finally, let's demonstrate using a generic machine (which then requires the scratchroot option), using the global grid for PFT and soil types, and setting the run length to two months.

Example: Running PTCLMmkdata with global PFT and soil types dataset
-------------------------------------------------------------------
::

   > cd $CTSMROOT/tools/PTCLM
   # Note, see the the Section called Converting AmeriFlux Data for use by PTCLMmkdata with permission information
   # to use the US-UMB data.
   > ./PTCLMmkdata -s US-UMB -d $CSMDATA --pftgrid --soilgrid
   > cd ../../../../../US-UMB_ICRUCLM45BGC
