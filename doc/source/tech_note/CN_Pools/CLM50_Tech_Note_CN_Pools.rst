.. _rst_CN Pools:

CN Pools
===================

Introduction
-----------------

CLM includes a prognostic treatment of the terrestrial carbon and
nitrogen cycles including natural vegetation, crops, and soil biogeochemistry. The model is
fully prognostic with respect to all carbon and nitrogen state variables
in the vegetation, litter, and soil organic matter. The seasonal timing
of new vegetation growth and litterfall is also prognostic, responding
to soil and air temperature, soil water availability, daylength, and
crop management practices in
varying degrees depending on a specified phenology type or management for each PFT
(Chapter 
:numref:`rst_Vegetation Phenology and Turnover`). The
prognostic LAI, SAI,
tissue stoichiometry, and vegetation heights are
utilized by the biophysical model that couples carbon, water, and
energy cycles.

Separate state variables for C and N are tracked for leaf, live stem,
dead stem, live coarse root, dead coarse root, fine root, and grain pools
(:numref:`Figure Vegetation fluxes and pools`). Each of these pools has two corresponding
storage pools representing, respectively, short-term and long-term
storage of non-structural carbohydrates and labile nitrogen. There are
two additional carbon pools, one for the storage of growth respiration
reserves, and another used to meet excess demand for maintenance
respiration during periods with low photosynthesis. One additional
nitrogen pool tracks retranslocated nitrogen, mobilized from leaf tissue
prior to abscission and litterfall. Altogether there are 23 state
variables for vegetation carbon, and 22 for vegetation nitrogen.

.. _Figure Vegetation fluxes and pools:

.. figure:: CLMCN_pool_structure_v2_lores.png
    :width: 753px
    :height: 513px

    Vegetation fluxes and pools for carbon cycle in CLM5.

In addition to the vegetation pools, CLM includes a series of
decomposing carbon and nitrogen pools as vegetation successively
breaks down to CWD, and/or litter, and subsequently to soil organic
matter. Discussion of the decomposition model, alternate
specifications of decomposition rates, and methods to rapidly
equilibrate the decomposition model, is in Chapter 
:numref:`rst_Decomposition`.

