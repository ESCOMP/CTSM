.. _rst_Plant Respiration:

Plant Respiration
=================
CLM5 includes changes to plant respiration including
 - A new leaf respiration algorithm based on Atkin et al. (2016)
 - A lower growth respiration coefficient, based on Atkin et al. (2017)

Autotrophic Respiration
----------------------------

The model treats maintenance and growth respiration fluxes separately, even though it is difficult to measure them as separate fluxes (Lavigne and Ryan, 1997; Sprugel et al., 1995). Maintenance respiration is defined as the carbon cost to support the metabolic activity of existing live tissue, while growth respiration is defined as the additional carbon cost for the synthesis of new growth.

Maintenance Respiration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Atkin et al. (2016) propose a model for leaf respiration that is based on the leaf nitrogen content per unit area (:math:`NS_{narea}` (gN m :sup:`2` leaf), with an intercept parameter that is PFT dependant, and an acclimation term that depends upon the average temperature of the previous 10 day period :math:`t_{2m,10days}`, in Celsius.

.. math::
   :label: 17.46)

   CF_{mr\_ leaf} =  i_{atkin,pft} + (NS_{narea}  0.2061) - (0.0402 (t_{2m,10days}))

The temperature dependance of leaf maintenance (dark) respiration is described in Chapter :numref:`rst_Stomatal Resistance and Photosynthesis`.

.. math::
   :label: 17.47)

   CF_{mr\_ livestem} \_ =NS_{livestem} MR_{base} MR_{Q10} ^{(T_{2m} -20)/10}

.. math::
   :label: 17.48)

   CF_{mr\_ livecroot} \_ =NS_{livecroot} MR_{base} MR_{Q10} ^{(T_{2m} -20)/10}

.. math::
   :label: 17.49)

   CF_{mr\_ froot} \_ =\sum _{j=1}^{nlevsoi}NS_{froot} rootfr_{j} MR_{base} MR_{Q10} ^{(Ts_{j} -20)/10}

where :math:`MR_{q10}` (= 2.0) is the temperature sensitivity for maintenance respiration, :math:`T_{2m}` (°C) is the air temperature at 2m height, :math:`Ts_{j}`* (°C) is the soil temperature at level *j*, and :math:`rootfr_{j}` is the fraction of fine roots distributed in soil level *j*.

.. table:: Atkin leaf respiration model intercept values.

 ========================  =============
 Plant functional type     :math:`i_{atkin}`
 ========================  =============
 NET Temperate                       1.499
 NET Boreal                          1.499
 NDT Boreal                          1.499
 BET Tropical                        1.756
 BET temperate                       1.756
 BDT tropical                        1.756
 BDT temperate                       1.756
 BDT boreal                          1.756
 BES temperate                       2.075
 BDS temperate                       2.075
 BDS boreal                          2.075
 C\ :sub:`3` arctic grass            2.196
 C\ :sub:`3` grass                   2.196
 C\ :sub:`4` grass                   2.196
 ========================  =============

Note that, for woody vegetation, maintenance respiration costs are not calculated for the dead stem and dead coarse root components. These components are assumed to consist of dead xylem cells, with no metabolic function. By separating the small live component of the woody tissue (ray parenchyma, phloem, and sheathing lateral meristem cells) from the larger fraction of dead woody tissue, it is reasonable to assume a common base maintenance respiration rate for all live tissue types.

The total maintenance respiration cost is then given as:

.. math::
   :label: 17.50)

   CF_{mr} =CF_{mr\_ leaf} +CF_{mr\_ froot} +CF_{mr\_ livestem} +CF_{mr\_ livecroot} .

.. _Growth Respiration:

Growth Respiration
^^^^^^^^^^^^^^^^^^^^^^^^^

Growth respiration is calculated as a factor of 0.11 times the total carbon allocation to new growth (:math:`CF_{growth}`, after allocating carbon for N acquisition, Chapter :numref:`rst_FUN`.) on a given timestep, based on construction costs for a range of woody and non-woody tissues, with estimates of the growth respiration flux revised downswards following (Atkin et al. 2017). For new carbon and nitrogen allocation that enters storage pools for subsequent display, it is not clear what fraction of the associated growth respiration should occur at the time of initial allocation, and what fraction should occur later, at the time of display of new growth from storage. Eddy covariance estimates of carbon fluxes in forest ecosystems suggest that the growth respiration associated with transfer of allocated carbon and nitrogen from storage into displayed tissue is not significant (Churkina et al., 2003), and so it is assumed in CLM that all of the growth respiration cost is incurred at the time of initial allocation, regardless of the fraction of allocation that is displayed immediately (i.e. regardless of the value of :math:`f_{cur}`, section 13.5). This behavior is parameterized in such a way that if future research suggests that some fraction of the growth respiration cost should be incurred at the time of display from storage, a simple parameter modification will effect the change. [1]_

.. [1]
   Parameter :math:`\text{grpnow}` in routines CNGResp and CNAllocation, currently set to 1.0, could be changed to a smaller value to transfer some portion (1 - :math:`\text{grpnow}` ) of the growth respiration forward in time to occur at the time of growth display from storage.

