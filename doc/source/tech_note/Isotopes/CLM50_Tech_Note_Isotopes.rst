.. _rst_Carbon Isotopes:

Carbon Isotopes
===================

CLM includes a fully prognostic representation of the fluxes, storage, and isotopic discrimination of the carbon isotopes :sup:`13`\ C and :sup:`14`\ C. The implementation of the C isotopes capability takes advantage of the CLM hierarchical data structures, replicating the carbon state and flux variable structures at the column and PFT level to track total carbon and both C isotopes separately (see description of data structure hierarchy in Chapter 2). For the most part, fluxes and associated updates to carbon state variables for :sup:`13`\ C are calculated directly from the corresponding total C fluxes. Separate calculations are required in a few special cases, such as where isotopic discrimination occurs, or where the necessary isotopic ratios are undefined. The general approach for :sup:`13`\ C flux and state variable calculation is described here, followed by a description of all the places where special calculations are required.

General Form for Calculating :sup:`13`\ C and :sup:`14`\ C Flux
--------------------------------------------------------------------------------

In general, the flux of :sup:`13`\ C and corresponding to a given flux of total C (:math:`{CF}_{13C}` and :math:`{CF}_{totC}`, respectively) is determined by :math:`{CF}_{totC}`, the masses of :sup:`13`\ C and total C in the upstream pools (:math:`{CS}_{13C\_up}` and :math:`{CS}_{totC\_up}`, respectively, i.e. the pools *from which* the fluxes of :sup:`13`\ C and total C originate), and a fractionation factor, :math:`{f}_{frac}`:

.. math::
   :label: ZEqnNum629812

   CF_{13C} =\left\{\begin{array}{l} {CF_{totC} \frac{CS_{13C\_ up} }{CS_{totC\_ up} } f_{frac} \qquad {\rm for\; }CS_{totC} \ne 0} \\ {0\qquad {\rm for\; }CS_{totC} =0} \end{array}\right\}

If the :math:`{f}_{frac}` = 1.0 (no fractionation), then the fluxes :math:`{CF}_{13C}` and :math:`{CF}_{totC}` will be in simple proportion to the masses :math:`{CS}_{13C\_up}` and :math:`{CS}_{totC\_up}`. Values of :math:`{f}_{frac} < 1.0` indicate a discrimination against the heavier isotope (:sup:`13`\ C) in the flux-generating process, while :math:`{f}_{frac}` :math:`>` 1.0 would indicate a preference for the heavier isotope. Currently, in all cases where Eq. is used to calculate a :sup:`13`\ C flux, :math:`{f}_{frac}` is set to 1.0.

For :sup:`14`\ C, no fractionation is used in either the initial photosynthetic step, nor in subsequent fluxes from upstream to downstream pools; as discussed below, this is because observations of :sup:`14` C are typically described in units that implicitly correct out the fractionation of :sup:`14`\ C by referencing them to :sup:`13`\ C ratios.

Isotope Symbols, Units, and Reference Standards
----------------------------------------------------

Carbon has two primary stable isotopes, :sup:`12`\ C and :sup:`13`\ C. :sup:`12`\ C is the most abundant, comprising about 99% of all carbon. The isotope ratio of a compound, :math:`{R}_{A}`, is the mass ratio of the rare isotope to the abundant isotope

.. math::
   :label: 30.2)

   R_{A} =\frac{{}^{13} C_{A} }{{}^{12} C_{A} } .

Carbon isotope ratios are often expressed using delta notation, :math:`\delta`. The :math:`\delta^{13}`\ C value of a compound A, :math:`\delta^{13}`\ C\ :sub:`A`, is the difference between the isotope ratio of the compound, :math:`{R}_{A}`, and that of the Pee Dee Belemnite standard, :math:`{R}_{PDB}`, in parts per thousand

.. math::
   :label: 30.3)

   \delta ^{13} C_{A} =\left(\frac{R_{A} }{R_{PDB} } -1\right)\times 1000

where :math:`{R}_{PDB}` = 0.0112372, and units of :math:`\delta` are per mil (‰).

Isotopic fractionation can be expressed in several ways. One expression of the fractionation factor is with alpha (:math:`\alpha`) notation. For example, the equilibrium fractionation between two reservoirs A and B can be written as:

.. math::
   :label: 30.4)

   \alpha _{A-B} =\frac{R_{A} }{R_{B} } =\frac{\delta _{A} +1000}{\delta _{B} +1000} .

This can also be expressed using epsilon notation (:math:`\epsilon`), where

.. math::
   :label: 30.5)

   \alpha _{A-B} =\frac{\varepsilon _{A-B} }{1000} +1

In other words, if :math:`{\epsilon }_{A-B} = 4.4` ‰ , then :math:`{\alpha}_{A-B} =1.0044`.

In addition to the stable isotopes :sup:`1`\ :sup:`2`\ C and :sup:`13`\ C, the unstable isotope :sup:`14`\ C is included in CLM. :sup:`14`\ C can also be described using the delta notation:

.. math::
   :label: 30.6)

   \delta ^{14} C=\left(\frac{A_{s} }{A_{abs} } -1\right)\times 1000

However, observations of :sup:`14`\ C are typically fractionation-corrected using the following notation:

.. math::
   :label: 30.7)

   \Delta {}^{14} C=1000\times \left(\left(1+\frac{\delta {}^{14} C}{1000} \right)\frac{0.975^{2} }{\left(1+\frac{\delta {}^{13} C}{1000} \right)^{2} } -1\right)

where :math:`\delta^{14}`\ C is the measured isotopic fraction and :math:`\mathrm{\Delta}^{14}`\ C corrects for mass-dependent isotopic fractionation processes (assumed to be 0.975 for fractionation of :sup:`13`\ C by photosynthesis). CLM assumes a background preindustrial atmospheric :sup:`14`\ C /C ratio of 10\ :sup:`-12`, which is used for A\ :sub::`abs`. For the reference standard A\ :math:`{}_{abs}`, which is a plant tissue and has a :math:`\delta^{13}`\ C value is :math:`\mathrm{-}`\ 25 ‰ due to photosynthetic discrimination, :math:`\delta`\ :sup:`14`\ C = :math:`\mathrm{\Delta}`\ :sup:`14`\ C. For CLM, in order to use the :sup:`14`\ C model independently of the :sup:`13`\ C model, for the :sup:`14`\ C calculations, this fractionation is set to zero, such that the 0.975 term becomes 1, the :math:`\delta^{13}`\ C term (for the calculation of :math:`\delta^{14}`\ C only) becomes 0, and thus :math:`\delta^{14}`\ C = :math:`\mathrm{\Delta}`\ :sup:`14`\ C.

Carbon Isotope Discrimination During Photosynthesis
--------------------------------------------------------

Photosynthesis is modeled in CLM as a two-step process: diffusion of CO\ :sub:`2` into the stomatal cavity, followed by enzymatic fixation (Chapter :numref:`rst_Stomatal Resistance and Photosynthesis`). Each step is associated with a kinetic isotope effect. The kinetic isotope effect during diffusion of CO\ :sub:`2` through the stomatal opening is 4.4‰. The kinetic isotope effect during fixation of CO\ :sub:`2` with Rubisco is :math:`\sim`\ 30‰; however, since about 5-10% of carbon in C3 plants reacts with phosphoenolpyruvate carboxylase (PEPC) (Melzer and O'Leary, 1987), the net kinetic isotope effect during fixation is :math:`\sim`\ 27‰ for C3 plants. In C4 plant photosynthesis, only the diffusion effect is important. The fractionation factor equations for C3 and C4 plants are given below:

For C4 plants,

.. math::
   :label: 30.8)

   \alpha _{psn} =1+\frac{4.4}{1000}

For C3 plants,

.. math::
   :label: 30.9)

   \alpha _{psn} =1+\frac{4.4+22.6\frac{c_{i}^{*} }{pCO_{2} } }{1000}

where :math:`{\alpha }_{psn}` is the fractionation factor, and :math:`c^*_i` and pCO\ :sub:`2` are the revised intracellular and atmospheric CO\ :sub:`2` partial pressure, respectively.

As can be seen from the above equation, kinetic isotope effect during fixation of CO\ :sub:`2` is dependent on the intracellular CO\ :sub:`2` concentration, which in turn depends on the net carbon assimilation. That is calculated during the photosynthesis calculation as follows:

.. math::
   :label: 30.10)

   c_{i} =pCO_{2} -a_{n} p\frac{\left(1.4g_{s} \right)+\left(1.6g_{b} \right)}{g_{b} g_{s} }

where :math:`a_n` is net carbon assimilation during photosynthesis, :math:`p` is atmospheric pressure, :math:`g_b` is leaf boundary layer conductance, and :math:`g_s` is leaf stomatal conductance.

Isotopic fractionation code is compatible with multi-layered canopy parameterization; i.e., it is possible to calculate varying discrimination rates for each layer of a multi-layered canopy. However, as with the rest of the photosynthesis model, the number of canopy layers is currently set to one by default.

:sup:`14`\ C radioactive decay and historical atmospheric :sup:`14`\ C and :sup:`13`\ C concentrations
------------------------------------------------------------------------------------------------------

In the preindustrial biosphere, radioactive decay of :sup:`14`\ C in carbon pools allows dating of long-term age since photosynthetic uptake; while over the 20\ :math:`{}^{th}` century, radiocarbon in the atmosphere was first diluted by radiocarbon-free fossil fuels and then enriched by aboveground thermonuclear testing to approximately double its long-term mean concentration. CLM includes both of these processes to allow comparison of carbon that may vary on multiple timescales with observed values.

For radioactive decay, at each timestep all :sup:`14`\ C pools are reduced at a rate of –log/:math:`\tau`, where :math:`\tau` is the half-life (Libby half-life value of 5568 years). In order to rapidly equilibrate the long-lived pools during accelerated decomposition spinup, the radioactive decay of the accelerated pools is also accelerated by the same degree as the decomposition, such that the :sup:`14`\ C value of these pools is in equilibrium when taken out of the spinup mode.

For variation of atmospheric :sup:`14`\ C and :sup:`13`\ C over the historical period, :math:`\mathrm{\Delta}`\ :sup:`14`\ C and :math:`\mathrm{\Delta}`\ :sup:`13`\ C values can be set to either fixed concentrations or time-varying concentrations read in from a file. A default file is provided that spans the historical period (:ref:`Graven et al., 2017 <Gravenetal2017>`). For :math:`\mathrm{\Delta}`\ :sup:`14`\ C, values are provided and read in for three latitude bands (30°N--90°N, 30°S--30°N, and 30°S--90°S).

