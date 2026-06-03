.. _rst_Plant Mortality:

Plant Mortality
===================
What has changed
----------------------------------------------

In CLM5, gap-phase mortality in non-FATES vegetation was commonly represented using a prescribed constant mortality rate (2% yr :sup:`-1`) applied uniformly across perennial plant functional types. However, this process should be better constrained given different climate zones (:ref:`Keller et al. 2004<Kelleretal2004>`; :ref:`Sollins 1982<Sollins1982>`), different species mixtures (:ref:`Gomes et al. 2003<Gomesetal2003>`), and different size and age classes (:ref:`Busing 2005<Busing2005>`; :ref:`Law et al. 2003<Lawetal2003>`). In CLM6, mortality formulations have evolved beyond a single fixed mortality parameter for all PFTs. The annual mortality (*r_mort*) has now become a PFT-dependent variable, with some woody PFTs by default assigned a value other than 0.02 (:numref:`Table Model default annual mortality rate`).

.. _Table Model default annual mortality rate:

.. table:: Model default annual mortality rate (*r_mort*) for woody PFTs

   +---------------------------------------------------+--------------------------------------+
   | Plant functional type                             | Mortality rate (yr :sup:`-1`)        | 
   +===================================================+======================================+
   | NET Temperate                                     | 0.0211945164991821                   |
   +---------------------------------------------------+--------------------------------------+
   | NET Boreal                                        | 0.0174                               |
   +---------------------------------------------------+--------------------------------------+
   | NDT Boreal                                        | 0.0198950093389492                   |
   +---------------------------------------------------+--------------------------------------+
   | BET Tropical                                      | 0.024                                |
   +---------------------------------------------------+--------------------------------------+
   | BET Temperate                                     | 0.0199981934178915                   |
   +---------------------------------------------------+--------------------------------------+
   | BDT Tropical                                      | 0.0200001818014196                   |
   +---------------------------------------------------+--------------------------------------+
   | BDT Temperate                                     | 0.0210684434513937                   |
   +---------------------------------------------------+--------------------------------------+
   | BDT Boreal                                        | 0.024                                |
   +---------------------------------------------------+--------------------------------------+
 
PFT name abbreviations: NET = Needleleaf Evergreen Tree, NDT = Needleleaf Deciduous Tree, BET = Broadleaf Evergreen Tree, BDT = Broadleaf Deciduous Tree



Conceptual representation
----------------------------------------------
This section describes plant mortality in the biogeochemistry (BGC) configuration of CLM for non-FATES vegetation. It focuses on the representation of gap-phase mortality, which accounts for the aggregate loss of plant biomass due to processes such as disturbance, competition, age-related decline, and environmental stress. 

In this framework, mortality is represented as a first-order loss process, in which all vegetation carbon and nitrogen pools experience proportional losses over time. The equations presented in this section describe the pool-level mortality fluxes and their routing within the biogeochemical framework. 

This section does not describe mechanistic mortality processes represented in the Functionally Assembled Terrestrial Ecosystem Simulator (FATES), where mortality emerges from explicit demographic, physiological, and disturbance processes (see Chapter :numref:`rst_Dynamic Global Vegetation and FATES`). Readers interested in those formulations should refer to the `FATES documentation`_. Mortality associated with fire and land-use or harvest processes is treated separately in the Fire and Land Use Change sections (see Chapters :numref:`rst_Fire` and :numref:`rst_Transient Landcover Change`, respectively). Legacy dynamic vegetation configurations (CNDV) used related mortality formulations but are no longer actively supported; mechanistic dynamic vegetation and mortality processes in CTSM are now handled through FATES.

.. _FATES documentation: https://fates-users-guide.readthedocs.io/en/latest/index.html


Implementation
----------------------------------------------

Vegetation carbon and nitrogen dynamics are implemented using a matrix-based formulation (:ref:`Lu et al. 2020<Luetal2020>`, :ref: `Liao et al. 2023`<Liaoetal2023>), in which the evolution of vegetation pools is governed by a combination of process-specific transfer and turnover operators. Within this framework, mortality is represented as one component of the overall vegetation transfer system, alongside phenological turnover and fire-related processes.

The matrix formulation provides a compact representation of vegetation dynamics and organizes the bookkeeping of transfers among pools, but does not alter the underlying conceptual treatment of mortality as a proportional loss from vegetation pools.


Mortality fluxes and routing
----------------------------------------------

Carbon and nitrogen removed from vegetation pools are routed according to tissue type:

- Leaf and fine-root pools are transferred to litter pools (labile, cellulose, and lignin components)

- Stem and coarse-root pools are transferred to coarse woody debris pools

- Storage and transfer pools are assumed to represent labile material and are therefore transferred to litter pools

This formulation provides a simplified, bulk representation of mortality consistent with the structure of non-demographic CLM configurations. A conceptual diagram of vegetation fluxes and pools (:numref:`Figure Vegetation fluxes and pools`) can be found in Chapter :numref:`rst_CN Pools`.


Mortality fluxes leaving vegetation pools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Whole-plant mortality is represented as a first-order process with a prescribed mortality rate. For an annual mortality fraction (:math:`a_m`, yr\ :sup:`-1`), the corresponding rate per unit time (second) :math:`m` is given by:

.. math::
   :label: 33.1)

   m = \frac{a_m}{365 \times 86400}

Mortality fluxes from vegetation carbon pools are computed as:

.. math::
   :label: 33.2)

   CF_{i,mort} =CS_{i} m

and similarly for nitrogen:

.. math::
   :label: 33.3)

   NF_{i,mort} =NS_{i} m

where :math:`CF_{i}` is carbon flux, :math:`CS_{i}` is carbon state variable (or pool), :math:`NF_{i}` is nitrogen flux, :math:`NS_{i}` is nitrogen state, for each vegetation pool :math:`i`, respectively.

These fluxes are applied to all vegetation pools, including:

- displayed pools: *leaf*, *froot*, *livestem*, *deadstem*, *livecroot*, *deadcroot*
- storage pools: *leaf_stor*, *froot_stor*, *livestem_stor*, *deadstem_stor*, *livecroot_stor*, *deadcroot_stor*, *gresp_stor*
- transfer pools: *leaf_xfer*, *froot_xfer*, *livestem_xfer*, *deadstem_xfer*, *livecroot_xfer*, *deadcroot_xfer*, *gresp_xfer*
- retranslocated pools: *retrans*

where *croot* refers to coarse roots, *froot* refers to fine roots, *gresp* refers to growth respiration, *retrans* refers to retranslocated, *stor* refers to storage, and *xfer* refers to transfer. Note that *gresp* only exists in carbon pools, and *retrans* only exists in nitrogen pools. 



Aggregation to column level
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Analogous to the treatment of litterfall fluxes (Chapter :numref:`rst_Vegetation Phenology and Turnover`), mortality fluxes leaving the vegetation pools are first computed at the plant functional type (PFT) level and then aggregated to the column level according to the weighted distribution of PFT :math:`p` on the column (:math:`wcol_{p}` ), and deposited in litter or coarse woody debris pools, which are defined at the column level. 

- Non-woody tissue pools

Carbon and nitrogen fluxes from mortality of displayed leaf and fine root into litter pools are calculated as

.. math::
   :label: 33.4)

   CF_{leaf\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{leaf\_ mort} f_{lab\_ leaf,p} wcol_{p}

.. math::
   :label: 33.5)

   CF_{leaf\_ mort,lit2} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{leaf\_ mort} f_{cel\_ leaf,p} wcol_{p}

.. math::
   :label: 33.6)

   CF_{leaf\_ mort,lit3} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{leaf\_ mort} f_{lig\_ leaf,p} wcol_{p}

.. math::
   :label: 33.7)

   CF_{froot\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{froot\_ mort} f_{lab\_ froot,p} wcol_{p}

.. math::
   :label: 33.8)

   CF_{froot\_ mort,lit2} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{froot\_ mort} f_{cel\_ froot,p} wcol_{p}

.. math::
   :label: 33.9)

   CF_{froot\_ mort,lit3} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{froot\_ mort} f_{lig\_ froot,p} wcol_{p}

where :math:`{f}_{lab\_leaf,p}`, :math:`{f}_{cel\_leaf,p}`, and :math:`{f}_{lig\_leaf,p}` are the labile (lit1), cellulose/hemicellulose (lit2), and lignin (lit3) fractions of leaf litter for PFT :math:`p`, 
and the same rule applies to the fine root litter fractions.
:math:`{wcol}_{p}` is the fractional contribution of PFT :math:`p` to the column, :math:`p` is an index through the plant functional types occurring on a column, and :math:`n_{\mathrm{pft}}` is the number of PFTs present in the column. 

Nitrogen fluxes to the litter pools are assumed to follow the C:N of the senescent tissue, and so are distributed using the same fractions used for carbon fluxes:

.. math::
   :label: 33.10)

   NF_{leaf\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{leaf\_ mort} f_{lab\_ leaf,p} wcol_{p}

.. math::
   :label: 33.11)

   NF_{leaf\_ mort,lit2} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{leaf\_ mort} f_{cel\_ leaf,p} wcol_{p}

.. math::
   :label: 33.12)

   NF_{leaf\_ mort,lit3} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{leaf\_ mort} f_{lig\_ leaf,p} wcol_{p}

.. math::
   :label: 33.13)

   NF_{froot\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{froot\_ mort} f_{lab\_ froot,p} wcol_{p}

.. math::
   :label: 33.14)

   NF_{froot\_ mort,lit2} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{froot\_ mort} f_{cel\_ froot,p} wcol_{p}

.. math::
   :label: 33.15)

   NF_{froot\_ mort,lit3} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{froot\_ mort} f_{lig\_ froot,p} wcol_{p}  .


- Woody tissue pools

Carbon and nitrogen mortality fluxes from displayed live and dead stem and coarse root pools are merged to the column level and deposited in the coarse woody debris (*cwd*) pools:

.. math::
   :label: 33.16)

   CF_{livestem\_ mort,cwd} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{livestem\_ mort} wcol_{p}

.. math::
   :label: 33.17)

   CF_{deadstem\_ mort,cwd} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{deadstem\_ mort} wcol_{p}

.. math::
   :label: 33.18)

   CF_{livecroot\_ mort,cwd} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{livecroot\_ mort} wcol_{p}

.. math::
   :label: 33.19)

   CF_{deadcroot\_ mort,cwd} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{deadcroot\_ mort} wcol_{p}

.. math::
   :label: 33.20)

   NF_{livestem\_ mort,cwd} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{livestem\_ mort} wcol_{p}

.. math::
   :label: 33.21)

   NF_{deadstem\_ mort,cwd} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{deadstem\_ mort} wcol_{p}

.. math::
   :label: 33.22)

   NF_{livecroot\_ mort,cwd} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{livecroot\_ mort} wcol_{p}

.. math::
   :label: 33.23)

   NF_{deadcroot\_ mort,cwd} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{deadcroot\_ mort} wcol_{p}

- Storage and transfer pools

All vegetation storage and transfer pools for carbon and nitrogen are assumed to exist as labile pools within the plant (e.g. as carbohydrate stores, in the case of carbon pools). This assumption applies to storage and transfer pools for both non-woody and woody tissues. The mortality fluxes from these pools are therefore assumed to be deposited in the labile litter pools (:math:`{CS}_{lit1}`, :math:`{NS}_{lit1}`), after being merged to the column level. 

Carbon mortality fluxes out of storage and transfer pools are:

.. math::
   :label: 33.24)

   CF_{leaf\_ stor\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{leaf\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.25)

   CF_{froot\_ stor\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{froot\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.26)

   CF_{livestem\_ stor\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{livestem\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.27)

   CF_{deadstem\_ stor\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{deadstem\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.28)

   CF_{livecroot\_ stor\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{livecroot\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.29)

   CF_{deadcroot\_ stor\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{deadcroot\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.30)

   CF_{gresp\_ stor\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{gresp\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.31)

   CF_{leaf\_ xfer\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{leaf\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.32)

   CF_{froot\_ xfer\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{froot\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.33)

   CF_{livestem\_ xfer\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{livestem\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.34)

   CF_{deadstem\_ xfer\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{deadstem\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.35)

   CF_{livecroot\_ xfer\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{livecroot\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.36)

   CF_{deadcroot\_ xfer\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{deadcroot\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.37)

   CF_{gresp\_ xfer\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}CF_{gresp\_ xfer\_ mort} wcol_{p}  .

Nitrogen mortality fluxes out of storage and transfer pools, including the storage pool for retranslocated nitrogen, are calculated as:

.. math::
   :label: 33.38)

   NF_{leaf\_ stor\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{leaf\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.39)

   NF_{froot\_ stor\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{froot\_ stor\_ mort} wcol_{p}
   
.. math::
   :label: 33.40)

   NF_{livestem\_ stor\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{livestem\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.41)

   NF_{deadstem\_ stor\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{deadstem\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.42)

   NF_{livecroot\_ stor\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{livecroot\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.43)

   NF_{deadcroot\_ stor\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{deadcroot\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.44)

   NF_{retrans\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{retrans\_ mort} wcol_{p}

.. math::
   :label: 33.45)

   NF_{leaf\_ xfer\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{leaf\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.46)

   NF_{froot\_ xfer\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{froot\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.47)

   NF_{livestem\_ xfer\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{livestem\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.48)

   NF_{deadstem\_ xfer\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{deadstem\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.49)

   NF_{livecroot\_ xfer\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{livecroot\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.50)

   NF_{deadcroot\_ xfer\_ mort,lit1} =\sum _{p=0}^{n_{\mathrm{pft}}}NF_{deadcroot\_ xfer\_ mort} wcol_{p}  .

