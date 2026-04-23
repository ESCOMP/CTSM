.. _rst_Plant Respiration:

Plant Respiration
=================

Autotrophic Respiration
----------------------------

The model treats maintenance and growth respiration fluxes separately, even though it is difficult to measure them as separate fluxes (:ref:`Lavigne and Ryan 1997 <LavigneRyan1997>`; :ref:`Sprugel et al. 1995 <Sprugeletal1995>`). Maintenance respiration is defined as the carbon cost to support the metabolic activity of existing live tissue, while growth respiration is defined as the additional carbon cost for the synthesis of new growth.

Maintenance Respiration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Maintenance respiration is calculated separately for leaves, live stems, live coarse roots, and fine roots. For leaf maintenance respiration (:math:`CF_{mr\_ leaf}`; gC m\ :sup:`-2` s\ :sup:`-1`), :ref:`Atkin et al. (2016)<Atkin2016>` propose a model that includes an intercept parameter that is PFT-dependent (:math:`i_{atkin,pft}`; see :numref:`Table Atkin leaf respiration model intercept values` for values), a term for leaf nitrogen content per unit area (:math:`NS_{narea}`; gN m\ :sup:`-2` leaf), and an acclimation term that depends on the average temperature of the previous 10 day period (:math:`t_{2m,10days}`; °C). Leaf maintenance respiration is calculated as:

.. math::
   :label: 17.46)

   CF_{mr\_ leaf} =  i_{atkin,pft} + (NS_{narea}  0.2061) - (0.0402 (t_{2m,10days}))

The temperature dependence of leaf maintenance (dark) respiration is described in Chapter :numref:`rst_Stomatal Resistance and Photosynthesis`. 

Live stem (:math:`CF_{mr\_ livestem}`; gC m\ :sup:`-2` s\ :sup:`-1`), live coarse root (:math:`CF_{mr\_ livecroot}`; gC m\ :sup:`-2` s\ :sup:`-1`), and fine root (:math:`CF_{mr\_ froot}`; gC m\ :sup:`-2` s\ :sup:`-1`) maintenance respiration are calculated as the product of a nitrogen-scaled base respiration rate and a Q10-based temperature scaling. :math:`MR_{base}` (gC gN\ :sup:`-1` s\ :sup:`-1`) represents the base respiration rate per unit nitrogen, :math:`MR_{Q10}` (= 2.0) is the temperature sensitivity for maintenance respiration, and :math:`T_{2m}` (°C) is the air temperature at 2m height. :math:`NS_{livestem}` (gN m\ :sup:`-2` live stem), :math:`NS_{livecroot}` (gN m\ :sup:`-2` live coarse root), and :math:`NS_{froot}` (gN m\ :sup:`-2` fine root) represent nitrogen content per unit area for live stem, live coarse root, and fine root, respectively. Finally, fine root respiration is calculated as the sum of contributions from each soil level *j*. :math:`rootfr_{j}` represents the fraction of fine roots distributed in soil level *j* and :math:`Ts_{j}` (°C) is the soil temperature at level *j*. Thus,

Live stem maintenance respiration is calculated as:

.. math::
   :label: 17.47)

   CF_{mr\_ livestem} =NS_{livestem} MR_{base} MR_{Q10} ^{(T_{2m} -20)/10}

Live coarse root respiration is calculated as:

.. math::
   :label: 17.48)

   CF_{mr\_ livecroot} =NS_{livecroot} MR_{base} MR_{Q10} ^{(T_{2m} -20)/10}

Fine root maintenance respiration is calculated as:

.. math::
   :label: 17.49)

   CF_{mr\_ froot} =\sum _{j=1}^{nlevsoi}NS_{froot} rootfr_{j} MR_{base} MR_{Q10} ^{(Ts_{j} -20)/10}

The total maintenance respiration cost (:math:`CF_{mr}`; gC m\ :sup:`-2` s\ :sup:`-1`) is then calculated as the sum of leaf (:math:`CF_{mr\_ leaf}`), fine root (:math:`CF_{mr\_ froot}`), live stem (:math:`CF_{mr\_ livestem}`), and live coarse root (:math:`CF_{mr\_ livecroot}`) components:

.. math::
   :label: 17.50)

   CF_{mr} =CF_{mr\_ leaf} +CF_{mr\_ froot} +CF_{mr\_ livestem} +CF_{mr\_ livecroot}


.. _Table Atkin leaf respiration model intercept values:

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

Note that, for woody vegetation, maintenance respiration costs are not calculated for dead stem and dead coarse root components. These components are assumed to consist of dead xylem cells with no metabolic function. By separating the small live component of the woody tissue (ray parenchyma, phloem, and sheathing lateral meristem cells) from the larger fraction of dead woody tissue, it is reasonable to assume a common base maintenance respiration rate for all live tissue types.


.. _Growth Respiration:

Growth Respiration
^^^^^^^^^^^^^^^^^^^^^^^^^

Growth respiration is calculated as a factor of 0.11 times the total carbon allocated to new growth (:math:`CF_{growth}`) after allocating carbon for N acquisition (see Chapter :numref:`rst_FUN`) on a given timestep, based on construction costs for a range of woody and non-woody tissues. Estimates of the growth respiration flux were revised downswards following :ref:`Atkin et al. (2017)<Atkin2017>`. For new carbon and nitrogen allocation that enters storage pools for subsequent display, it is not clear what fraction of the associated growth respiration should occur at the time of initial allocation, and what fraction should occur later, at the time of display of new growth from storage. Eddy covariance estimates of carbon fluxes in forest ecosystems suggest that the growth respiration associated with transfer of allocated carbon and nitrogen from storage into displayed tissue is not significant (:ref:`Churkina et al. 2003 <Churkinaetal2003>`), so it is assumed in CLM that all of the growth respiration cost is incurred at the time of initial allocation, regardless of the fraction of allocation that is displayed immediately (i.e. regardless of the value of :math:`f_{cur}`, section :numref:`Carbon Allocation to New Growth`). This behavior is parameterized in such a way that if future research suggests that some fraction of the growth respiration cost should be incurred at the time of display from storage, a simple parameter modification will effect the change. [1]_

.. [1]
   Parameter :math:`\text{grpnow}` in routines CNGResp and CNAllocation, currently set to 1.0, could be changed to a smaller value to transfer some portion (1 - :math:`\text{grpnow}` ) of the growth respiration forward in time to occur at the time of growth display from storage.
