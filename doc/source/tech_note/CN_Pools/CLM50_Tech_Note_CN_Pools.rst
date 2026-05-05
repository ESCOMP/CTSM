.. _rst_CN Pools:

CN Pools
===================

Introduction
-----------------

CLM includes a prognostic treatment of the terrestrial carbon and nitrogen cycles including natural vegetation, crops, and soil biogeochemistry. The model is fully prognostic with respect to all carbon and nitrogen state variables in the vegetation, litter, and soil organic matter. The seasonal timing of new vegetation growth and litterfall is also prognostic, responding to soil and air temperature, soil water availability, daylength, and crop management practices in varying degrees depending on a specified phenology type or management for each PFT (Chapter :numref:`rst_Vegetation Phenology and Turnover`). The prognostic LAI, SAI, tissue stoichiometry, and vegetation heights are utilized by the biophysical model that couples carbon, water, and energy cycles.

Separate state variables for C and N are tracked for leaf, live stem, dead stem, live coarse root, dead coarse root, fine root, and grain pools (:numref:`Figure Vegetation fluxes and pools`). Each of these pools has two corresponding storage pools representing, respectively, short-term and long-term storage of non-structural carbohydrates and labile nitrogen. There are two additional carbon pools, one for the storage of growth respiration reserves, and another used to meet excess demand for maintenance respiration during periods with low photosynthesis. One additional nitrogen pool tracks retranslocated nitrogen, mobilized from leaf tissue prior to abscission and litterfall. Altogether there are 23 state variables for vegetation carbon, and 22 for vegetation nitrogen.

.. _Figure Vegetation fluxes and pools:

.. figure:: CLMCN_pool_structure_v2_lores.png
    :width: 753px
    :height: 513px

    Vegetation fluxes and pools for carbon cycle in CLM5.

In addition to the vegetation pools, CLM includes a series of decomposing carbon and nitrogen pools as vegetation successively breaks down to CWD, and/or litter, and subsequently to soil organic matter. Discussion of the decomposition model, alternate specifications of decomposition rates, and methods to rapidly equilibrate the decomposition model, is in Chapter :numref:`rst_Decomposition`.

Tissue Stoichiometry
-----------------------

As of CLM5, vegetation tissues have a flexible stoichiometry, as described in :ref:`Ghimire et al. (2016) <Ghimireetal2016>`. Each tissue has a target C\:N ratio, with the target leaf C\:N, :math:`CN_{\text{target}}^\text{pft}`, varying by plant functional type (PFT) (:numref:`Table PFT target leaf CN`). Nitrogen is allocated at each timestep to allow the plant to best match the target stoichiometry. Nitrogen downregulation of productivity acts by increasing the actual C\:N ratio of leaves when insufficient nitrogen is available to meet stoichiometric demands of leaf growth, thereby reducing the N available for photosynthesis and reducing the :math:`V_{\text{c,max25}}` and :math:`J_{\text{max25}}` terms, as described in Chapter :numref:`rst_Photosynthetic Capacity`. Details of the flexible tissue stoichiometry are described in Chapter :numref:`rst_CN Allocation`.

As of CLM5.4, the target leaf C\:N may be time-evolving, :math:`CN_{\text{target}}^{\text{pft,CO2}}`, as a logarithmic function of atmospheric CO\ :sub:`2` that we denote :math:`CN_{\text{perturb}}^{\text{CO2}}`:

.. math::
  :label: time-evolv target leaf CN

  \begin{split}
  CN_{\text{perturb}}^{\text{CO2}} &= CN_{\text{slope}}^{\text{CO2}} \cdot \ln\left(\frac{\text{CO2}_{\text{atm}}}{\text{CO2}_{\text{atm}}^{\text{ref}}}\right), \text{where } CN_{\text{perturb}}^{\text{CO2}} &\ge 0 \\
  CN_{\text{target}}^{\text{pft,CO2}} &= CN_{\text{target}}^\text{pft} + CN_{\text{perturb}}^{\text{CO2}}
  \end{split}

where :math:`CN_{\text{target}}^\text{pft}` is the time-invarying target leaf C\:N at reference CO\ :sub:`2` that depends on PFT, :math:`CN_{\text{slope}}^{\text{CO2}}` (unitless) is the slope of the function, :math:`\text{CO2}_{\text{atm}}` is atmospheric CO\ :sub:`2` in parts per million by volume (ppmv), and :math:`\text{CO2}_{\text{atm}}^{\text{ref}}` is the reference CO\ :sub:`2` (ppmv) above which atmospheric CO\ :sub:`2` begins to scale the target leaf C\:N. The optional time-evolving target leaf C\:N was documented in :ref:`Hauser et al. (2023) <Hauseretal2023>`, and its current default is off by setting :math:`CN_{\text{slope}}^{\text{CO2}} = 0`.

.. _Table PFT target leaf CN:

.. table:: Plant functional type (PFT) target leaf C:N, :math:`CN_{\text{target}}^\text{pft}`.

 +--------------------------+-----------------+
 | PFT                      | target leaf C:N |
 +==========================+=================+
 | NET Temperate            |       58.00     |
 +--------------------------+-----------------+
 | NET Boreal               |       60.24     |
 +--------------------------+-----------------+
 | NDT Boreal               |       28.92     |
 +--------------------------+-----------------+
 | BET Tropical             |       36.03     |
 +--------------------------+-----------------+
 | BET temperate            |       34.59     |
 +--------------------------+-----------------+
 | BDT tropical             |       18.63     |
 +--------------------------+-----------------+
 | BDT temperate            |       21.64     |
 +--------------------------+-----------------+
 | BDT boreal               |       17.09     |
 +--------------------------+-----------------+
 | BES temperate            |       36.42     |
 +--------------------------+-----------------+
 | BDS temperate            |       23.26     |
 +--------------------------+-----------------+
 | BDS boreal               |       21.40     |
 +--------------------------+-----------------+
 | C\ :sub:`3` arctic grass |       20.70     |
 +--------------------------+-----------------+
 | C\ :sub:`3` grass        |       29.39     |
 +--------------------------+-----------------+
 | C\ :sub:`4` grass        |       35.36     |
 +--------------------------+-----------------+
 | Temperate Corn           |       25.00     |
 +--------------------------+-----------------+
 | Spring Wheat             |       20.00     |
 +--------------------------+-----------------+
 | Temperate Soybean        |       20.00     |
 +--------------------------+-----------------+
 | Cotton                   |       20.00     |
 +--------------------------+-----------------+
 | Rice                     |       20.00     |
 +--------------------------+-----------------+
 | Sugarcane                |       25.00     |
 +--------------------------+-----------------+
 | Tropical Corn            |       25.00     |
 +--------------------------+-----------------+
 | Tropical Soybean         |       20.00     |
 +--------------------------+-----------------+
 | Miscanthus               |       20.00     |
 +--------------------------+-----------------+
 | Switchgrass              |       20.00     |
 +--------------------------+-----------------+

