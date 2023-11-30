.. include:: ../substitutions.rst

.. _scientific-validiation:

========================
 Scientific Validation
========================

In this section we go over what has been extensively tested and scientifically validated with |version|, and maybe more importantly what has NOT been tested and may NOT be scientifically validated. You can use all features of CLM, but need to realize that some things haven't been tested extensively or validated scientifically. When you use these features you may run into trouble doing so, and will need to do your own work to make sure the science is reasonable.

--------------------------------------------------------------
Standard Configuration and Namelist Options that are Validated
--------------------------------------------------------------

See `http://www.cesm.ucar.edu/models/cesm1.2/clm/CLM_configurations_CESM1.2.pdf <http://www.cesm.ucar.edu/models/cesm1.2/clm/CLM_configurations_CESM1.2.pdf>`_ for an explanation of what configurations are scientifically validated for |version|. For CLM4.0 changes to the science of the model are minimal since CESM1.1.1 so we expect answers to be very similar to using it.

In the sections below we go through configuration and/or namelist options or modes that the user should be especially wary of using. You are of course free to use these options, and you may find that they work functionally. Although in some cases you will find issues even with functionality of using them. If so you will need to test, debug and find solutions for these issues on your own. But in every case you will need to go through more extensive work to validate these options from a scientific standpoint. Some of these options are only for |version| while others are for both CLM4.0 AND |version| we explicitly say which they apply to.

-----------------------------------------------
Configurations that should be used with caution
-----------------------------------------------

There are some options in |version| that are available but either not tested extensively, or not scientifically evaluated. These options should be used with caution. And any options that deviate from the scientifically supported configurations can have issues. The IMPORTANT_NODES file goes into more details on this.

The IMPORTANT_NOTES (which can be found in ``$CTSMROOT/doc``) is repeated here.

.. include:: ../../../IMPORTANT_NOTES
   :literal:
