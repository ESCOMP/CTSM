.. _rst_Ecosystem Demography with FATES:

Ecosystem Demography with FATES
===================================

Ecosystem Demography
^^^^^^^^^^^^^^^^^^^^^^^

An important component of representing the land surface is accurately capturing vegetation dynamics and their impact on and interaction with the Earth system. Some models utilize dynamic global vegetation models (DGVMs), which track the fractional cover of plant functional types (PFTs) on a grid cell over time, usually driven by "climate envelopes" which use climate metrics (e.g., temperature and precipitation) to determine where PFTs can grow. However, many DGVMs still do not simulate important processes like mortality, regeneration, and plant competition, which we know are important for accurately capturing vegetation and carbon dynamics (:ref:`Fisher et al. 2015<Fisheretal2015>`).

Ecosystem demography models explicitly represent the size structure and successional state of vegetation, through direct simulation of plant growth, mortality, and regeneration. Thus, important vegetation characteristics such as vegetation canopy height, succession, and even potential biome shifts become emergent properties of the model rather than being prescribed.

FATES
^^^^^^^^^^^^^^^^^^^^

FATES is the "Functionally Assembled Terrestrial Ecosystem Simulator". It is an external module which can run within a given "Host Land Model" (HLM) like CLM. FATES was derived from the CLM Ecosystem Demography model (CLM(ED)), documented in :ref:`Fisher et al. (2015)<Fisheretal2015>`.

FATES is derived from the work of :ref:`Moorcroft et al. (2001)<mc_2001>` and is a cohort model of vegetation competition and co-existence, allowing a representation of the biosphere which accounts for the division of the land surface into successional stages, and for competition for light between height structured cohorts of representative trees of various plant functional types. Individual plants within FATES are grouped into "cohorts" of the same size and PFT, and these cohorts compete for light and resources on individual "patches" that represent different disturbance histories (:numref:`fig_FATES_vs_CLM`). This type of ecosystem heterogeneity is in contrast to the default vegetation model in CLM, which uses two (sunlit & shaded) "big leaf" canopies per PFT, each on their own patch, with no representation of within-canopy structural heterogeneity or disturbance history. 

.. _fig_FATES_vs_CLM:

.. figure:: FATES_tiling.png
  
   Comparison of grid structure in CLM-FATES and CLM. (a) An example grid of FATES patches showing individual cohorts of plants of different size and PFT. (b) those same FATES patches represented as their age-since-disturbance. (c) the equivalent gridcell in default CLM, showing instead a set of specific PFTs with fixed areas.

FATES also introduces a new organizational structure to CLM which differs from its default organization. The original hierarchical organization of CLM as described in :ref:`Oleson et al. (2013)<Olesonetal2013>` are gridcells which contain land units (e.g. vegetated, lake, urban), which contain columns (e.g., naturally vegetated), which contain PFTs. FATES replaces the PFT level of this hierarchy with patches and cohorts. Thus, gridcells contain land units, which contain columns, which contain patches, which in turn contain some number of cohorts of the same or different PFT (:numref:`fig_FATES_and_CLM_hierarchy`). 

.. _fig_FATES_and_CLM_hierarchy:

.. figure:: FATES_and_CLM_hierarchy.png
   :width: 60%
  
   Comparison of the organizational hierarchy of default CLM (left) and CLM-FATES (right).

The implementation of the Ecosystem Demography concept within FATES links the surface flux and canopy physiology concepts in CLM with numerous additional developments necessary to accommodate the new model. These include a version of the SPITFIRE (Spread and InTensity of Fire) model of :ref:`Thonicke et al. (2010)<thonickeetal2010>`, and an adoption of the concept of `Perfect Plasticity Approximation` approach of :ref:`Purves et al. 2008<purves2008>`, :ref:`Lichstein et al. 2011<lichstein2011>` and :ref:`Weng et al. 2014<weng2014>`, in accounting for the spatial arrangement of crowns. Novel algorithms accounting for the fragmentation of coarse woody debris into chemical litter streams, for the physiological optimization of canopy thickness, for the accumulation of seeds in the seed bank, for multi-layer multi-PFT radiation transfer, for drought-deciduous and cold-deciduous phenology, for carbon storage allocation, and for tree mortality under carbon stress, are also included.


FATES Reduced Complexity Modes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Currently, FATES can be run in several different "reduced complexity modes", where parts of the vegetation model are driven by input data rather than simulated. These modes can be used to facilitate calibration, test features, or run simulations more quickly. These modes are:

1. **Satellite Phenology** (SP) mode: this mode is designed to run with leaf area index (LAI), stem area index (SAI), and canopy height (HTOP) as input to the model. As such, all processes that are normally used to calculate these values are turned off (e.g., mortality, allocation, etc.)

2. **No-Competition Mode**: this mode runs with full complexity in terms of processes, but places each FATES PFT on its own patch. As such, PFTs do not compete with one another. 

3. **Fixed Biogeography Mode**: this mode turns off prognostic spatial changes in the distribution of vegetation and instead, the model uses input data to determine which PFTs are present at any given gridcell. The patch area for each PFT is derived from the input CLM surface dataset. However, please note that the PFTs in the FATES parameter file do not always map one-to-one with the CLM PFTs on the surface dataset. See the FATES parameter *fates_hlm_pft_map* on the FATES parameter file for the correct mapping of FATES to CLM PFTs.

Note that there are different combinations of no-competition and fixed biogeography mode that will result in different behaviors. See the `FATES namelist documentation <https://fates-users-guide.readthedocs.io/en/latest/user/namelist-options.html>`_ for these options.

Scientifically Supported CLM-FATES Configurations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We currently scientifally support CLM-FATES run in `carbon-only mode <https://fates-users-guide.readthedocs.io/en/latest/user/PARTEH-Modes.html>`_, and with either SP or no-competition + fixed biogeography mode.


Further reading
^^^^^^^^^^^^^^^^^^^^

For more information about FATES, including a Users Guide and Technical Note, please see the `FATES documentation <https://fates-users-guide.readthedocs.io/en/latest/index.html>`_.
