.. _rst_Dynamic Global Vegetation and FATES:

Dynamic Global Vegetation and FATES
===================================

What has changed
^^^^^^^^^^^^^^^^^^^^

- Deprecation of the dynamic global vegetation model (DGVM): The CLM5.0 model contains the legacy 'CNDV' code, which runs the CLM biogeochemistry model in combination with the LPJ-derived dynamics vegetation model introduced in CLM3. While this capacity has not technically been removed from the model, the DGVM has not been tested in the development of CLM5 and is no longer scientifically supported.

- Introduction of FATES: The Functionally Assembled Terrestrial Ecosystem Simulator (FATES) is the actively developed DGVM for the CLM5.

FATES
^^^^^^^^^^^^^^^^^^^^

FATES is the "Functionally Assembled Terrestrial Ecosystem Simulator". It is an external module which can run within a given "Host Land Model" (HLM) like CLM.

FATES was derived from the CLM Ecosystem Demography model (CLM(ED)), which was documented in :ref:`Fisher et al. (2015)<Fisheretal2015>`.

The Ecosystem Demography ('ED'), concept within FATES is derived from the work of :ref:`Moorcroft et al. (2001)<mc_2001>` and is a cohort model of vegetation competition and co-existence, allowing a representation of the biosphere which accounts for the division of the land surface into successional stages, and for competition for light between height structured cohorts of representative trees of various plant functional types.

The implementation of the Ecosystem Demography concept within FATES links the surface flux and canopy physiology concepts in CLM with numerous additional developments necessary to accommodate the new model. These include a version of the SPITFIRE (Spread and InTensity of Fire) model of :ref:`Thonicke et al. (2010)<thonickeetal2010>`, and an adoption of the concept of `Perfect Plasticity Approximation` approach of :ref:`Purves et al. 2008<purves2008>`, :ref:`Lichstein et al. 2011<lichstein2011>` and :ref:`Weng et al. 2014<weng2014>`, in accounting for the spatial arrangement of crowns. Novel algorithms accounting for the fragmentation of coarse woody debris into chemical litter streams, for the physiological optimization of canopy thickness, for the accumulation of seeds in the seed bank, for multi-layer multi-PFT radiation transfer, for drought-deciduous and cold-deciduous phenology, for carbon storage allocation, and for tree mortality under carbon stress, are also included.


FATES Reduced Complexity Modes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Currently, FATES can be run in several different "reduced complexity modes", where parts of the vegetation model are driven by input data rather than simulated. These modes can be used to facilitate calibration, test features, or run simulations more quickly. These modes are:

1. **Satellite Phenology** (SP) mode: this mode is designed to run with leaf area index (LAI), stem area index (SAI), and canopy height (HTOP) as input to the model. As such, all processes that are normally used to calculate these values are turned off (e.g., mortality, allocation, etc.)

2. **No-Competition Mode**: this mode runs with full complexity in terms of processes, but places each FATES PFT on its own patch. As such, PFTs do not compete with one another. 

3. **Fixed Biogeography Mode**: this mode turns off prognostic spatial changes in the distribution of vegetation and instead, the model uses input data to determine which PFTs are present at any given gridcell. The patch area for each PFT is derived from the input CLM surface dataset. However, please note that the PFTs in the FATES parameter file do not always map one-to-one with the CLM PFTs on the surface dataset. See the FATES parameter *fates_hlm_pft_map* on the FATES parameter file for the correct mapping of FATES to CLM PFTs.

Note that there are different combinations of no-competition and fixed biogeography mode that will result in different behaviors. See the `FATES documentation<https://fates-users-guide.readthedocs.io/en/latest/user/namelist-options.html>`_ for these options.

Scientifically Supported CLM-FATES Configurations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We currently scientifally support CLM-FATES run in `carbon-only mode<https://fates-users-guide.readthedocs.io/en/latest/user/PARTEH-Modes.html>`_, and with either SP or no-competition + fixed biogeography mode.


Further reading
^^^^^^^^^^^^^^^^^^^^

For more information about FATES, including a Users Guide and Technical Note, please see the `FATES documentation`_.

.. _FATES documentation: https://fates-users-guide.readthedocs.io/en/latest/index.html
