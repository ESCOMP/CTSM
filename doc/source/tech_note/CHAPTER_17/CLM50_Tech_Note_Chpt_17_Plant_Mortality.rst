.. _rst_Plant Mortality:

Plant Mortality
===================

Plant mortality as described here applies to perennial vegetation types,
and is intended to represent the death of individuals from a stand of
plants due to the aggregate of processes such as wind throw, insect
attack, disease, extreme temperatures or drought, and age-related
decline in vigor. These processes are referred to in aggregate as
“gap-phase” mortality. Mortality due to fire and anthropogenic land
cover change are treated separately (see Chapters 18 and 21,
respectively).

Mortality Fluxes Leaving Vegetation Pools
----------------------------------------------

Whole-plant mortality is parameterized very simply, assuming a mortality
rate of 2% yr\ :sup:`-1` for all vegetation types. This is clearly
a gross oversimplification of an important process, and additional work
is required to better constrain this process in different climate zones
(Keller et al. 2004; Sollins 1982), for different species mixtures
(Gomes et al. 2003), and for different size and age classes (Busing
2005; Law et al. 2003). Literature values for forest mortality rates
range from at least 0.7% to 3.0% yr\ :sup:`-1`. Taking the annual
rate of mortality (*am*, proportion yr\ :sup:`-1`) as 0.02, a
mortality rate per second (*m*) is calculated as
:math:`m={am\mathord{\left/ {\vphantom {am \left(365\cdot 86400\right)}} \right. \kern-\nulldelimiterspace} \left(365\cdot 86400\right)}` .
All vegetation carbon and nitrogen pools for display, storage, and
transfer are affected at rate *m*, with mortality fluxes out of
vegetation pools eventually merged to the column level and deposited in
litter pools. Mortality (*mort*) fluxes out of displayed vegetation
carbon and nitrogen pools are

.. math::
   :label: 17.1) 

   CF_{leaf\_ mort} =CS_{leaf} m

.. math::
   :label: 17.2) 

   CF_{froot\_ mort} =CS_{froot} m

.. math::
   :label: 17.3) 

   CF_{livestem\_ mort} =CS_{livestem} m

.. math::
   :label: 17.4) 

   CF_{deadstem\_ mort} =CS_{deadstem} m

.. math::
   :label: 17.5) 

   CF_{livecroot\_ mort} =CS_{livecroot} m

.. math::
   :label: 17.6) 

   CF_{deadcroot\_ mort} =CS_{deadcroot} m

.. math::
   :label: 17.7) 

   NF_{leaf\_ mort} =NS_{leaf} m

.. math::
   :label: 17.8) 

   NF_{froot\_ mort} =NS_{froot} m

.. math::
   :label: 17.9) 

   NF_{livestem\_ mort} =NS_{livestem} m

.. math::
   :label: 17.10) 

   NF_{deadstem\_ mort} =NS_{deadstem} m

.. math::
   :label: 17.11) 

   NF_{livecroot\_ mort} =NS_{livecroot} m

.. math::
   :label: 17.12) 

   NF_{deadcroot\_ mort} =NS_{deadcroot} m

.. math::
   :label: 17.13) 

   NF_{retrans\_ mort} =NS_{retrans} m.

where CF are carbon fluxes, CS is carbon storage, NF are nitrogen
fluxes, NS is nitrogen storage, *croot* refers to coarse roots, *froot*
refers to fine roots, and *retrans* refers to retranslocated.

Mortality fluxes out of carbon and nitrogen storage (*stor)* pools are

.. math::
   :label: 17.14) 

   CF_{leaf\_ stor\_ mort} =CS_{leaf\_ stor} m

.. math::
   :label: 17.15) 

   CF_{froot\_ stor\_ mort} =CS_{froot\_ stor} m

.. math::
   :label: 17.16) 

   CF_{livestem\_ stor\_ mort} =CS_{livestem\_ stor} m

.. math::
   :label: 17.17) 

   CF_{deadstem\_ stor\_ mort} =CS_{deadstem\_ stor} m

.. math::
   :label: 17.18) 

   CF_{livecroot\_ stor\_ mort} =CS_{livecroot\_ stor} m

.. math::
   :label: 17.19) 

   CF_{deadcroot\_ stor\_ mort} =CS_{deadcroot\_ stor} m

.. math::
   :label: 17.20) 

   CF_{gresp\_ stor\_ mort} =CS_{gresp\_ stor} m

.. math::
   :label: 17.21) 

   NF_{leaf\_ stor\_ mort} =NS_{leaf\_ stor} m

.. math::
   :label: 17.22) 

   NF_{froot\_ stor\_ mort} =NS_{froot\_ stor} m

.. math::
   :label: 17.23) 

   NF_{livestem\_ stor\_ mort} =NS_{livestem\_ stor} m

.. math::
   :label: 17.24) 

   NF_{deadstem\_ stor\_ mort} =NS_{deadstem\_ stor} m

.. math::
   :label: 17.25) 

   NF_{livecroot\_ stor\_ mort} =NS_{livecroot\_ stor} m

.. math::
   :label: 17.26) 

   NF_{deadcroot\_ stor\_ mort} =NS_{deadcroot\_ stor} m

where *gresp* refers to growth respiration.

Mortality fluxes out of carbon and nitrogen transfer (*xfer)* growth
pools are

.. math::
   :label: 17.27) 

   CF_{leaf\_ xfer\_ mort} =CS_{leaf\_ xfer} m

.. math::
   :label: 17.28) 

   CF_{froot\_ xfer\_ mort} =CS_{froot\_ xfer} m

.. math::
   :label: 17.29) 

   CF_{livestem\_ xfer\_ mort} =CS_{livestem\_ xfer} m

.. math::
   :label: 17.30) 

   CF_{deadstem\_ xfer\_ mort} =CS_{deadstem\_ xfer} m

.. math::
   :label: 17.31) 

   CF_{livecroot\_ xfer\_ mort} =CS_{livecroot\_ xfer} m

.. math::
   :label: 17.32) 

   CF_{deadcroot\_ xfer\_ mort} =CS_{deadcroot\_ xfer} m

.. math::
   :label: 17.33) 

   CF_{gresp\_ xfer\_ mort} =CS_{gresp\_ xfer} m

.. math::
   :label: 17.34) 

   NF_{leaf\_ xfer\_ mort} =NS_{leaf\_ xfer} m

.. math::
   :label: 17.35) 

   NF_{froot\_ xfer\_ mort} =NS_{froot\_ xfer} m

.. math::
   :label: 17.36) 

   NF_{livestem\_ xfer\_ mort} =NS_{livestem\_ xfer} m

.. math::
   :label: 17.37) 

   NF_{deadstem\_ xfer\_ mort} =NS_{deadstem\_ xfer} m

.. math::
   :label: 17.38) 

   NF_{livecroot\_ xfer\_ mort} =NS_{livecroot\_ xfer} m

.. math::
   :label: 17.39) 

   NF_{deadcroot\_ xfer\_ mort} =NS_{deadcroot\_ xfer} m

Mortality Fluxes Merged to the Column Level
------------------------------------------------

Analogous to the treatment of litterfall fluxes, mortality fluxes
leaving the vegetation pools are merged to the column level according to
the weighted distribution of PFTs on the column (:math:`wcol_{p}` ), and
deposited in litter and coarse woody debris pools, which are defined at
the column level. Carbon and nitrogen fluxes from mortality of displayed
leaf and fine root into litter pools are calculated as

.. math::
   :label: 17.40) 

   CF_{leaf\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{leaf\_ mort} f_{lab\_ leaf,p} wcol_{p}

.. math::
   :label: 17.41) 

   CF_{leaf\_ mort,lit2} =\sum _{p=0}^{npfts}CF_{leaf\_ mort} f_{cel\_ leaf,p} wcol_{p}

.. math::
   :label: 17.42) 

   CF_{leaf\_ mort,lit3} =\sum _{p=0}^{npfts}CF_{leaf\_ mort} f_{lig\_ leaf,p} wcol_{p}

.. math::
   :label: 17.43) 

   CF_{froot\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{froot\_ mort} f_{lab\_ froot,p} wcol_{p}

.. math::
   :label: 17.44) 

   CF_{froot\_ mort,lit2} =\sum _{p=0}^{npfts}CF_{froot\_ mort} f_{cel\_ froot,p} wcol_{p}

.. math::
   :label: 17.45) 

   CF_{froot\_ mort,lit3} =\sum _{p=0}^{npfts}CF_{froot\_ mort} f_{lig\_ froot,p} wcol_{p}

.. math::
   :label: 17.46) 

   NF_{leaf\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{leaf\_ mort} f_{lab\_ leaf,p} wcol_{p}

.. math::
   :label: 17.47) 

   NF_{leaf\_ mort,lit2} =\sum _{p=0}^{npfts}NF_{leaf\_ mort} f_{cel\_ leaf,p} wcol_{p}

.. math::
   :label: 17.48) 

   NF_{leaf\_ mort,lit3} =\sum _{p=0}^{npfts}NF_{leaf\_ mort} f_{lig\_ leaf,p} wcol_{p}

.. math::
   :label: 17.49) 

   NF_{froot\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{froot\_ mort} f_{lab\_ froot,p} wcol_{p}

.. math::
   :label: 17.50) 

   NF_{froot\_ mort,lit2} =\sum _{p=0}^{npfts}NF_{froot\_ mort} f_{cel\_ froot,p} wcol_{p}

.. math::
   :label: 17.51) 

   NF_{froot\_ mort,lit3} =\sum _{p=0}^{npfts}NF_{froot\_ mort} f_{lig\_ froot,p} wcol_{p}  .

where *lab* refers to labile, *cel* refers to cellulose, and *lig*
refers to lignin. Carbon and nitrogen mortality fluxes from displayed
live and dead stem and coarse root pools are merged to the column level
and deposited in the coarse woody debris (*cwd*) pools:

.. math::
   :label: 17.52) 

   CF_{livestem\_ mort,cwd} =\sum _{p=0}^{npfts}CF_{livestem\_ mort} wcol_{p}

.. math::
   :label: 17.53) 

   CF_{deadstem\_ mort,cwd} =\sum _{p=0}^{npfts}CF_{deadstem\_ mort} wcol_{p}

.. math::
   :label: 17.54) 

   CF_{livecroot\_ mort,cwd} =\sum _{p=0}^{npfts}CF_{livecroot\_ mort} wcol_{p}

.. math::
   :label: 17.55) 

   CF_{deadcroot\_ mort,cwd} =\sum _{p=0}^{npfts}CF_{deadcroot\_ mort} wcol_{p}

.. math::
   :label: 17.56) 

   NF_{livestem\_ mort,cwd} =\sum _{p=0}^{npfts}NF_{livestem\_ mort} wcol_{p}

.. math::
   :label: 17.57) 

   NF_{deadstem\_ mort,cwd} =\sum _{p=0}^{npfts}NF_{deadstem\_ mort} wcol_{p}

.. math::
   :label: 17.58) 

   NF_{livecroot\_ mort,cwd} =\sum _{p=0}^{npfts}NF_{livecroot\_ mort} wcol_{p}

.. math::
   :label: 17.59) 

   NF_{deadcroot\_ mort,cwd} =\sum _{p=0}^{npfts}NF_{deadcroot\_ mort} wcol_{p}

All vegetation storage and transfer pools for carbon and nitrogen are
assumed to exist as labile pools within the plant (e.g. as carbohydrate
stores, in the case of carbon pools). This assumption applies to storage
and transfer pools for both non-woody and woody tissues. The mortality
fluxes from these pools are therefore assumed to be deposited in the
labile litter pools (:math:`{CS}_{lit1}`, :math:`{NS}_{lit1}`),
after being merged to the column level. Carbon mortality fluxes out of
storage and transfer pools are:

.. math::
   :label: 17.60) 

   CF_{leaf\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{leaf\_ stor\_ mort} wcol_{p}

.. math::
   :label: 17.61) 

   CF_{froot\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{froot\_ stor\_ mort} wcol_{p}

.. math::
   :label: 17.62) 

   CF_{livestem\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{livestem\_ stor\_ mort} wcol_{p}

.. math::
   :label: 17.63) 

   CF_{deadstem\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{deadstem\_ stor\_ mort} wcol_{p}

.. math::
   :label: 17.64) 

   CF_{livecroot\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{livecroot\_ stor\_ mort} wcol_{p}

.. math::
   :label: 17.65) 

   CF_{deadcroot\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{deadcroot\_ stor\_ mort} wcol_{p}

.. math::
   :label: 17.66) 

   CF_{gresp\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{gresp\_ stor\_ mort} wcol_{p}

.. math::
   :label: 17.67) 

   CF_{leaf\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{leaf\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 17.68) 

   CF_{froot\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{froot\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 17.69) 

   CF_{livestem\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{livestem\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 17.70) 

   CF_{deadstem\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{deadstem\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 17.71) 

   CF_{livecroot\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{livecroot\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 17.72) 

   CF_{deadcroot\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{deadcroot\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 17.73) 

   CF_{gresp\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}CF_{gresp\_ xfer\_ mort} wcol_{p}  .

Nitrogen mortality fluxes out of storage and transfer pools, including
the storage pool for retranslocated nitrogen, are calculated as:

.. math::
   :label: 17.74) 

   NF_{leaf\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{leaf\_ stor\_ mort} wcol_{p}

.. math::
   :label: 17.75) 

   NF_{froot\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{froot\_ stor\_ mort} wcol_{p}

.. math::
   :label: 17.76) 

   NF_{livestem\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{livestem\_ stor\_ mort} wcol_{p}

.. math::
   :label: 17.77) 

   NF_{deadstem\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{deadstem\_ stor\_ mort} wcol_{p}

.. math::
   :label: 17.78) 

   NF_{livecroot\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{livecroot\_ stor\_ mort} wcol_{p}

.. math::
   :label: 17.79) 

   NF_{deadcroot\_ stor\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{deadcroot\_ stor\_ mort} wcol_{p}

.. math::
   :label: 17.80) 

   NF_{retrans\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{retrans\_ mort} wcol_{p}

.. math::
   :label: 17.81) 

   NF_{leaf\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{leaf\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 17.82) 

   NF_{froot\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{froot\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 17.83) 

   NF_{livestem\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{livestem\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 17.84) 

   NF_{deadstem\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{deadstem\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 17.85) 

   NF_{livecroot\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{livecroot\_ xfer\_ mort} wcol_{p}

.. math::
   :label: 17.86) 

   NF_{deadcroot\_ xfer\_ mort,lit1} =\sum _{p=0}^{npfts}NF_{deadcroot\_ xfer\_ mort} wcol_{p}  .


