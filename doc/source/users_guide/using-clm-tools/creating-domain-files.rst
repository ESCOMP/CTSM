.. _creating-domain-files:

.. include:: ../substitutions.rst

*****************************
 Creating CLM domain files
*****************************

*gen_domain* to create a domain file for datm from a mapping file. The domain file is then used by BOTH DATM AND CLM to define the grid and land-mask. The general data flow is shown in two figures. `Figure 2-4 <CLM-URL>`_ shows the general flow for a general global case (or for a regional grid that DOES include ocean). `Figure 2-5 <CLM-URL>`_ shows the use of **mknoocnmap.pl** (see `the Section called Using mknocnmap.pl to create grid and maps for single-point regional grids <CLM-URL>`_) to create a regional or single-point map file that is then run through **gen_domain** to create the domain file for it. As stated before `Figure 2-2 <CLM-URL>`_ is the legend for both of these figures. See `the 
tools/mapping/gen_domain_files/README <CLM-URL>`_ file for more help on **gen_domain**.

Here we create domain files for a regular global domain.

Figure 2-4. Global Domain file creation
=======================================

Insert figure 2-4

Starting from SCRIP grid files for both your atmosphere and ocean, you use **tools/mapping/gen_mapping_files/gen_cesm_maps.sh** to create a mapping file between the atmosphere and ocean. That mapping file is then used as input to **gen_domain** to create output domain files for both atmosphere and ocean. The atmosphere domain file is then used by both CLM and DATM for I compsets, while the ocean domain file is ignored. For this process you have to define your SCRIP grid files on your own. For a regional or single-point case that doesn't include ocean see `Figure 2-5 <CLM-URL>`_. (See `Figure 2-2 <CLM-URL>`_ for the legend for this figure.)

Note, that the SCRIP grid file used to start this process, is also used in **mkmapdata.sh** (see `the Section called Creating mapping files that mksurfdata_map will use <CLM-URL>`_). Next we create domain files for a single-point or regional domain.

Figure 2-5. Domain file creation using mknoocnmap.pl
====================================================
Insert figure 2-5

For a regular latitude/longitude grid that can be used for regional or single point simulations -- you can use **mknoocnmap.pl**. It creates a SCRIP grid file that can then be used as input to **mkmapdata.sh** as well as a SCRIP mapping file that is then input to **gen_domain**. The output of **gen_domain** is a atmosphere domain file used by both CLM and DATM and a ocean domain file that is ignored. (See `Figure 2-2 <CLM-URL>`_ for the legend for this figure.)

In this case the process creates both SCRIP grid files to be used by **mkmapdata.sh** as well as the domain files that will be used by both CLM and DATM.
