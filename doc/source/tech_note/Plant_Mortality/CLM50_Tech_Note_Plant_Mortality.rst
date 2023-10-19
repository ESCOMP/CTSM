.. _rst_Plant Mortality:

Plant Mortality
===================

Plant mortality as described here applies to perennial vegetation types, and is intended to represent the death of individuals from a stand of plants due to the aggregate of processes such as wind throw, insect attack, disease, extreme temperatures or drought, and age-related decline in vigor. These processes are referred to in aggregate as "gap-phase" mortality. Mortality due to fire and anthropogenic land cover change are treated separately (see Chapters :numref:`rst_Fire` and :numref:`rst_Transient Landcover Change`, respectively).

Mortality Fluxes Leaving Vegetation Pools
----------------------------------------------

Whole-plant mortality is parameterized very simply, assuming a mortality rate of 2% yr\ :sup:`-1` for all vegetation types. This is clearly a gross oversimplification of an important process, and additional work is required to better constrain this process in different climate zones (:ref:`Keller et al. 2004<Kelleretal2004>`; :ref:`Sollins 1982<Sollins1982>`), for different species mixtures (:ref:`Gomes et al. 2003<Gomesetal2003>`), and for different size and age classes (:ref:`Busing 2005<Busing2005>`; :ref:`Law et al. 2003<Lawetal2003>`). Literature values for forest mortality rates range from at least 0.7% to 3.0% yr\ :sup:`-1`. Taking the annual rate of mortality (*am*, proportion yr\ :sup:`-1`) as 0.02, a mortality rate per second (*m*) is calculated as :math:`m={am\mathord{\left/ {\vphantom {am \left(365\cdot 86400\right)}} \right.} \left(365\cdot 86400\right)}`. All vegetation carbon and nitrogen pools for display, storage, and transfer are affected at rate *m*, with mortality fluxes out of vegetation pools eventually merged to the column level and deposited in litter pools. Mortality (*mort*) fluxes out of displayed vegetation carbon and nitrogen pools are

.. math::
   :label: 33.1)

   CF_{leaf\_ mort} =CS_{leaf} m

.. math::
   :label: 33.2)

   CF_{froot\_ mort} =CS_{froot} m

.. math::
   :label: 33.3)

   CF_{livestem\_ mort} =CS_{livestem} m

.. math::
   :label: 33.4)

   CF_{deadstem\_ mort} =CS_{deadstem} m

.. math::
   :label: 33.5)

   CF_{livecroot\_ mort} =CS_{livecroot} m

.. math::
   :label: 33.6)

   CF_{deadcroot\_ mort} =CS_{deadcroot} m

.. math::
   :label: 33.7)

   NF_{leaf\_ mort} =NS_{leaf} m

.. math::
   :label: 33.8)

   NF_{froot\_ mort} =NS_{froot} m

.. math::
   :label: 33.9)

   NF_{livestem\_ mort} =NS_{livestem} m

.. math::
   :label: 33.10)

   NF_{deadstem\_ mort} =NS_{deadstem} m

.. math::
   :label: 33.11)

   NF_{livecroot\_ mort} =NS_{livecroot} m

.. math::
   :label: 33.12)

   NF_{deadcroot\_ mort} =NS_{deadcroot} m

.. math::
   :label: 33.13)

   NF_{retrans\_ mort} =NS_{retrans} m.

where CF are carbon fluxes, CS is carbon storage, NF are nitrogen fluxes, NS is nitrogen storage, *croot* refers to coarse roots, *froot* refers to fine roots, and *retrans* refers to retranslocated.

Mortality fluxes out of carbon and nitrogen storage (*stor)* pools are

.. math::
   :label: 33.14)

   CF_{leaf\_ stor\_ mort} =CS_{leaf\_ stor} m

.. math::
   :label: 33.15)

   CF_{froot\_ stor\_ mort} =CS_{froot\_ stor} m

.. math::
   :label: 33.16)

   CF_{livestem\_ stor\_ mort} =CS_{livestem\_ stor} m

.. math::
   :label: 33.17)

   CF_{deadstem\_ stor\_ mort} =CS_{deadstem\_ stor} m

.. math::
   :label: 33.18)

   CF_{livecroot\_ stor\_ mort} =CS_{livecroot\_ stor} m

.. math::
   :label: 33.19)

   CF_{deadcroot\_ stor\_ mort} =CS_{deadcroot\_ stor} m

.. math::
   :label: 33.20)

   CF_{gresp\_ stor\_ mort} =CS_{gresp\_ stor} m

.. math::
   :label: 33.21)

   NF_{leaf\_ stor\_ mort} =NS_{leaf\_ stor} m

.. math::
   :label: 33.22)

   NF_{froot\_ stor\_ mort} =NS_{froot\_ stor} m

.. math::
   :label: 33.23)

   NF_{livestem\_ stor\_ mort} =NS_{livestem\_ stor} m

.. math::
   :label: 33.24)

   NF_{deadstem\_ stor\_ mort} =NS_{deadstem\_ stor} m

.. math::
   :label: 33.25)

   NF_{livecroot\_ stor\_ mort} =NS_{livecroot\_ stor} m

.. math::
   :label: 33.26)

   NF_{deadcroot\_ stor\_ mort} =NS_{deadcroot\_ stor} m

where *gresp* refers to growth respiration.

Mortality fluxes out of carbon and nitrogen transfer (*xfer)* growth pools are

.. math::
   :label: 33.27)

   CF_{leaf\_ xfer\_ mort} =CS_{leaf\_ xfer} m

.. math::
   :label: 33.28)

   CF_{froot\_ xfer\_ mort} =CS_{froot\_ xfer} m

.. math::
   :label: 33.29)

   CF_{livestem\_ xfer\_ mort} =CS_{livestem\_ xfer} m

.. math::
   :label: 33.30)

   CF_{deadstem\_ xfer\_ mort} =CS_{deadstem\_ xfer} m

.. math::
   :label: 33.31)

   CF_{livecroot\_ xfer\_ mort} =CS_{livecroot\_ xfer} m

.. math::
   :label: 33.32)

   CF_{deadcroot\_ xfer\_ mort} =CS_{deadcroot\_ xfer} m

.. math::
   :label: 33.33)

   CF_{gresp\_ xfer\_ mort} =CS_{gresp\_ xfer} m

.. math::
   :label: 33.34)

   NF_{leaf\_ xfer\_ mort} =NS_{leaf\_ xfer} m

.. math::
   :label: 33.35)

   NF_{froot\_ xfer\_ mort} =NS_{froot\_ xfer} m

.. math::
   :label: 33.36)

   NF_{livestem\_ xfer\_ mort} =NS_{livestem\_ xfer} m

.. math::
   :label: 33.37)

   NF_{deadstem\_ xfer\_ mort} =NS_{deadstem\_ xfer} m

.. math::
   :label: 33.38)

   NF_{livecroot\_ xfer\_ mort} =NS_{livecroot\_ xfer} m

.. math::
   :label: 33.39)

   NF_{deadcroot\_ xfer\_ mort} =NS_{deadcroot\_ xfer} m

Mortality Fluxes Merged to the Column Level
------------------------------------------------

Analogous to the treatment of litterfall fluxes, mortality fluxes leaving the vegetation pools are merged to the column level according to the weighted distribution of PFTs on the column (:math:`wcol_{p}` ), and deposited in litter and coarse woody debris pools, which are defined at the column level. Carbon and nitrogen fluxes from mortality of displayed leaf and fine root into litter pools are calculated as

.. math::
   :label: 33.40)

   CF_{leaf\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{leaf\_ mort} f_{lab\_ leaf,p} wcol_{p}

.. math::
   :label: 33.41)

   CF_{leaf\_ mort,lit2} =\sum _{p=0}^{npfts}CF_{leaf\_ mort} f_{cel\_ leaf,p} wcol_{p}

.. math::
   :label: 33.42)

   CF_{leaf\_ mort,lit3} =\sum _{p=0}^{npfts}CF_{leaf\_ mort} f_{lig\_ leaf,p} wcol_{p}

.. math::
   :label: 33.43)

   CF_{froot\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{froot\_ mort} f_{lab\_ froot,p} wcol_{p}

.. math::
   :label: 33.44)

   CF_{froot\_ mort,lit2} =\sum _{p=0}^{npfts}CF_{froot\_ mort} f_{cel\_ froot,p} wcol_{p}

.. math::
   :label: 33.45)

   CF_{froot\_ mort,lit3} =\sum _{p=0}^{npfts}CF_{froot\_ mort} f_{lig\_ froot,p} wcol_{p}

.. math::
   :label: 33.46)

   NF_{leaf\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{leaf\_ mort} f_{lab\_ leaf,p} wcol_{p}

.. math::
   :label: 33.47)

   NF_{leaf\_ mort,lit2} =\sum _{p=0}^{npfts}NF_{leaf\_ mort} f_{cel\_ leaf,p} wcol_{p}

.. math::
   :label: 33.48)

   NF_{leaf\_ mort,lit3} =\sum _{p=0}^{npfts}NF_{leaf\_ mort} f_{lig\_ leaf,p} wcol_{p}

.. math::
   :label: 33.49)

   NF_{froot\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{froot\_ mort} f_{lab\_ froot,p} wcol_{p}

.. math::
   :label: 33.50)

   NF_{froot\_ mort,lit2} =\sum _{p=0}^{npfts}NF_{froot\_ mort} f_{cel\_ froot,p} wcol_{p}

.. math::
   :label: 33.51)

   NF_{froot\_ mort,lit3} =\sum _{p=0}^{npfts}NF_{froot\_ mort} f_{lig\_ froot,p} wcol_{p}  .

where *lab* refers to labile, *cel* refers to cellulose, and *lig* refers to lignin. Carbon and nitrogen mortality fluxes from displayed live and dead stem and coarse root pools are merged to the column level and deposited in the coarse woody debris (*cwd*) pools:

.. math::
   :label: 33.52)

   CF_{livestem\_ mort,cwd} =\sum _{p=0}^{npfts}CF_{livestem\_ mort} wcol_{p}

.. math::
   :label: 33.53)

   CF_{deadstem\_ mort,cwd} =\sum _{p=0}^{npfts}CF_{deadstem\_ mort} wcol_{p}

.. math::
   :label: 33.54)

   CF_{livecroot\_ mort,cwd} =\sum _{p=0}^{npfts}CF_{livecroot\_ mort} wcol_{p}

.. math::
   :label: 33.55)

   CF_{deadcroot\_ mort,cwd} =\sum _{p=0}^{npfts}CF_{deadcroot\_ mort} wcol_{p}

.. math::
   :label: 33.56)

   NF_{livestem\_ mort,cwd} =\sum _{p=0}^{npfts}NF_{livestem\_ mort} wcol_{p}

.. math::
   :label: 33.57)

   NF_{deadstem\_ mort,cwd} =\sum _{p=0}^{npfts}NF_{deadstem\_ mort} wcol_{p}

.. math::
   :label: 33.58)

   NF_{livecroot\_ mort,cwd} =\sum _{p=0}^{npfts}NF_{livecroot\_ mort} wcol_{p}

.. math::
   :label: 33.59)

   NF_{deadcroot\_ mort,cwd} =\sum _{p=0}^{npfts}NF_{deadcroot\_ mort} wcol_{p}

All vegetation storage and transfer pools for carbon and nitrogen are assumed to exist as labile pools within the plant (e.g. as carbohydrate stores, in the case of carbon pools). This assumption applies to storage and transfer pools for both non-woody and woody tissues. The mortality fluxes from these pools are therefore assumed to be deposited in the labile litter pools (:math:`{CS}_{lit1}`, :math:`{NS}_{lit1}`), after being merged to the column level. Carbon mortality fluxes out of storage and transfer pools are:

.. math::
   :label: 33.60)

   CF_{leaf\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{leaf\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.61)

   CF_{froot\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{froot\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.62)

   CF_{livestem\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{livestem\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.63)

   CF_{deadstem\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{deadstem\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.64)

   CF_{livecroot\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{livecroot\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.65)

   CF_{deadcroot\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{deadcroot\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.66)

   CF_{gresp\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{gresp\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.67)

   CF_{leaf\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{leaf\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.68)

   CF_{froot\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{froot\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.69)

   CF_{livestem\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{livestem\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.70)

   CF_{deadstem\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{deadstem\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.71)

   CF_{livecroot\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{livecroot\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.72)

   CF_{deadcroot\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{deadcroot\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.73)

   CF_{gresp\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{gresp\_ xfer\_ mort} wcol_{p}  .

Nitrogen mortality fluxes out of storage and transfer pools, including the storage pool for retranslocated nitrogen, are calculated as:

.. math::
   :label: 33.74)

   NF_{leaf\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{leaf\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.75)

   NF_{froot\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{froot\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.76)

   NF_{livestem\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{livestem\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.77)

   NF_{deadstem\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{deadstem\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.78)

   NF_{livecroot\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{livecroot\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.79)

   NF_{deadcroot\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{deadcroot\_ stor\_ mort} wcol_{p}

.. math::
   :label: 33.80)

   NF_{retrans\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{retrans\_ mort} wcol_{p}

.. math::
   :label: 33.81)

   NF_{leaf\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{leaf\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.82)

   NF_{froot\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{froot\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.83)

   NF_{livestem\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{livestem\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.84)

   NF_{deadstem\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{deadstem\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.85)

   NF_{livecroot\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{livecroot\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 33.86)

   NF_{deadcroot\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{deadcroot\_ xfer\_ mort} wcol_{p}  .

