.. include:: ../substitutions.rst

.. _creating-domain-files:

*****************************
 Creating CLM domain files
*****************************

.. todo::
    Delete this page? Domain files aren't needed with nuopc.

``gen_domain`` to create a domain file for datm from a mapping file. ``gen_domain`` is a tool that is a part of CIME. The domain file is then used by BOTH DATM AND CLM to define the grid and land-mask. The general data flow is shown in two figures. :numref:`Figure mkmapdata.sh` shows the general flow for a general global case (or for a regional grid that DOES include ocean). :numref:`Figure mknoocnmap.pl` shows the use of ``mknoocnmap.pl`` (see :ref:`using-mkocnmap`) to create a regional or single-point map file that is then run through ``gen_domain`` to create the domain file for it. As stated before :numref:`Figure Data_Flow_Legend` is the legend for both of these figures. See `the $CIMEROOT/tools/mapping/gen_domain_files/README <https://github.com/ESMCI/cime/blob/master/tools/mapping/gen_domain_files/README>`_ file for more help on ``gen_domain``.

Here we create domain files for a regular global domain.

Global Domain file creation
===========================

.. _Figure Global-Domain:

.. figure:: GlobalDomain.jpeg

  Global Domain file creation

Starting from SCRIP grid files for both your atmosphere and ocean, you use ``$CIMEROOT/tools/mapping/gen_mapping_files/gen_cesm_maps.sh`` to create a mapping file between the atmosphere and ocean. That mapping file is then used as input to ``gen_domain`` to create output domain files for both atmosphere and ocean. The atmosphere domain file is then used by both CLM and DATM for I compsets, while the ocean domain file is ignored. For this process you have to define your SCRIP grid files on your own. For a regional or single-point case that doesn't include ocean see :numref:`Figure mknoocnmap.pl`. (See :numref:`Figure Global-Domain` for the legend for this figure.)

Note that the SCRIP grid file used to start this process is also used in ``mkmapdata.sh`` (see :ref:`using-mkocnmap`). Next we create domain files for a single-point or regional domain.

Domain file creation using mknoocnmap.pl
========================================

.. _Figure mknoocnmap.pl:

.. figure:: mknoocnmap.jpeg

  Domain file creation using mknoocnmap.pl

For a regular latitude/longitude grid that can be used for regional or single point simulations -- you can use ``mknoocnmap.pl``. It creates a SCRIP grid file that can then be used as input to ``mkmapdata.sh`` as well as a SCRIP mapping file that is then input to ``gen_domain``. The output of ``gen_domain`` is a atmosphere domain file used by both CLM and DATM and a ocean domain file that is ignored. (See :numref:`Figure mknoocnmap.pl` for the legend for this figure.)

In this case the process creates both SCRIP grid files to be used by ``mkmapdata.sh`` as well as the domain files that will be used by both CLM and DATM.
