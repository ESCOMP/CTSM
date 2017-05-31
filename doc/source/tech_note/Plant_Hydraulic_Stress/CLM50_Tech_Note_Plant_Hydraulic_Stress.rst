.. _rst_Plant Hydraulic Stress:

Plant Hydraulic Stress
======================

Hello

.. _Plant Water Supply:

Plant Water Supply
-----------------------

PHS explicitly models water transport through the vegetation according to a simple hydraulic framework following Darcy's Law for porous media flow equations. The supply equations are used to solve for vegetation water potential forced by a given transpiration and the set of layer-by-layer soil water potentials.

The water supply is discretized into segments: soil-to-root, root-to-stem, and stem-to-leaf. There are typically several (1-49) soil-to-root flows operating in parallel, one per soil layer. There are two stem-to-leaf flows operating in parallel, corresponding to the sunlit and shaded "leaves".

In general the water fluxes (e.g. soil-to-root, root-to-stem, etc) are calculated as:

.. math::
   :label: 11.1) 

   q = kA\Delta\psi


:math:`q` is the flux of water (mmH\ :sub:`2`\ O)

:math:`k` is the hydraulic conductance (s\ :sup:`-1`\ )

:math:`A` is the area basis (m\ :sup:`2`\ /m\ :sup:`2`\ )

:math:`\Delta\psi` is the gradient in water potential (mmH\ :sub:`2`\ O)

.. math::
   :label: 11.2)
  
   k=k_{max}\cdot 2^{-\left(\dfrac{\psi}{p50}\right)^{c_k}}

:math:`k_{max}` is the maximum segment conductance (s-1) 

:math:`p50` is the water potential at 50% loss of conductivity (mmH2O) 

:math:`\psi` is the water potential of the lower segment terminus (mmH2O)

.. math:: 
   :label: 11.3)

   q_{1a}=k_{1a}*\mbox{LAI}_{sun}*\left(\psi_{stem}-\psi_{sunleaf} \right) 

.. math:: 
   :label: 11.4)

   q_{1b}=k_{1b}*\mbox{LAI}_{shade}*\left(\psi_{stem}-\psi_{shadeleaf} \right) 

.. math:: 
   :label: 11.5)

   k_{1a}=k_{1a,max}*2^{-\left(\dfrac{\psi_{stem}}{p50_1}\right)^{c_k}} 

.. math::
   :label: 11.6)
  
   k_{1b}=k_{1b,max}*2^{-\left(\dfrac{\psi_{stem}}{p50_1}\right)^{c_k}}

Variables:

:math:`LAI_{sun}` = sunlit leaf area index (m2/m2)

:math:`LAI_{shade}` = shaded leaf area index (m2/m2)

:math:`\psi_{stem}` = stem water potential (mmH20)

:math:`\psi_{sunleaf}` = sunlit leaf water potential (mmH20)

:math:`\psi_{shadeleaf}` = shaded leaf water potential (mmH20)

Parameters:

:math:`k_{1a,max}` = maximum leaf conductance (s-1)

:math:`k_{1b,max}` = maximum leaf conductance (s-1)

:math:`p50_{1}` = water potential at 50% loss of conductance (mmH2O)

:math:`c_{k}` = vulnerability curve shape-fitting parameter (-)

.. math::
   :label: 11.7)
  
   q_2=k_2 \cdot SAI \cdot \left( \psi_{root} - \psi_{stem} - \Delta \psi_z  \right)

.. math::
   :label: 11.8)

   k_2=\dfrac{k_{2,max}}{z_2} \cdot 2^{-\left(\dfrac{\psi_{root}}{p50_2}\right)^{c_k}}

Variables:

:math:`SAI` = stem area index (m2/m2)

:math:`\Delta\psi_z` = gravitational potential (mmH2O)

:math:`\psi_{root}` = root water potential (mmH2O)

:math:`\psi_{stem}` = stem water potential (mmH2O)

Parameters:

:math:`k_{2,max}` = maximum stem conductivity (m/s)

:math:`p50_2` = water potential at 50% loss of conductivity (mmH2O)

:math:`z_2` = vegetation height (m)

.. math::
   :label: 11.9)

   q_{3,i}=k_{3,i}*RAI*\left(\psi_{soil,i}-\psi_{root} + \Delta\psi_{z,i} \right)

.. math::
   :label: 11.10)

   RAI=\left(LAI+SAI \right)\cdot r_i \cdot f_{root-leaf}

.. math::
   :label: 11.11)

   k_{3,i}=\dfrac{k_{root}*k_{soil-to-root}}{k_{root}+k_{soil-to-root}} 

.. math::
   :label: 11.12)

   k_{root,i}=\dfrac{k_{3,max}}{z_{3,i}}*2^{-\left(\dfrac{\psi_{soil,i}}{p50_3}\right)^{c_k}}

.. math::
   :label: 11.13)

   k_{soil-to-root} = \dfrac{k_{soil,i}}{dx_{root,i}} 

.. math::
   :label: 11.14)

   dx_{root,i} = \left(\pi*\mbox{root-length-density}_i\right)^{-0.5}    

.. math::
   :label: 11.15)

   \mbox{root-length-density} = \dfrac{\mbox{total root length}}{\mbox{soil volume}} 

Variables:

:math:`\Delta\psi_{z,i}` = change in gravitational potential from soil layer :math:`i` to surface (mmH2O)

:math:`LAI` = total leaf area index (m2/m2)

:math:`SAI` = stem area index (m2/m2) 

:math:`\psi_{soil,i}` = water potential in soil layer :math:`i` (mmH2O)

:math:`\psi_{root}` = root water potential (mmH2O)

:math:`z_{3,i}` = length of root tissue conducting path = soil layer depth + root lateral length (m)

:math:`r_i` = root fraction in soil layer :math:`i` (-)

:math:`k_{soil,i}` = Brooks-Corey soil conductivity in soil layer :math:`i` (m/s)

Parameters:

:math:`f_{root-leaf}` = root-to-shoot ratio (-)

:math:`p50_3` = water potential at 50% loss of root tissue conductance (mmH2O)

:math:`ck` = shape-fitting parameter for vulnerability curve (-)

.. _Plant Water Demand:

Plant Water Demand
-----------------------

.. math::
   :label: 11.16)

   E_{sun} = E_{sun,max}*2^{-\left(\dfrac{\psi_{sunleaf}}{p50_e}\right)^{c_k}} 

.. math::
   :label: 11.17)

   E_{shade} = E_{shade,max}*2^{-\left(\dfrac{\psi_{shadeleaf}}{p50_e}\right)^{c_k}} 

.. math::
   :label: 11.18)

   B_{t,sun} = \dfrac{g_{s,sun}}{g_{s,sun,B_t=1}} 

.. math::
   :label: 11.19)

   B_{t,shade} = \dfrac{g_{s,shade}}{g_{s,shade,B_t=1}} 

:math:`E_{sun}` = sunlit leaf transpiration (mm/s)

:math:`E_{shade}` = shaded leaf transpiration (mm/s)

:math:`E_{sun,max}` = sunlit leaf transpiration absent water stress (mm/s)

:math:`E_{shade,max}` = shaded leaf transpiration absent water stress (mm/s)

:math:`\psi_{sunleaf}` = sunlit leaf water potential (mmH2O)

:math:`\psi_{shadeleaf}` = shaded leaf water potential (mmH2O) 

:math:`g_{s,sun}` = stomatal conductance of water corresponding to :math:`E_{sun}`

:math:`g_{s,shade}` = stomatal conductance of water corresponding to :math:`E_{shade}`

:math:`g_{s,sun,max}` = stomatal conductance of water corresponding to :math:`E_{sun,max}`

:math:`g_{s,shade,max}` = stomatal conductance of water corresponding to :math:`E_{shade,max}`

.. _Vegetation Water Potential:

Vegetation Water Potential
-----------------------------

PHS explicitly models root, stem, shaded leaf, and sunlit leaf water potential at each timestep. PHS iterates to find the vegetation water potential vector :math:`\psi` that satisfies continuity in the non-linear vegetation water supply and demand equations.

.. math::
   :label: 11.20)

   \psi=\left[\psi_{sunleaf},\psi_{shadeleaf},\psi_{stem},\psi_{root}\right]

.. math::
   :label: 11.21

   \begin{aligned}
   E_{sun}&=q_{1a}\\
   E_{shade}&=q_{1b}\\
   E_{sun}+E_{shade}&=q_{1a}+q_{1b}\\
   &=q_2\\
   &=\sum_{i=1}^{nlevsoi}{q_{3,i}}
   \end{aligned}

The demand terms (left-hand side) are decreasing functions of absolute leaf water potential. As absolute leaf water potential becomes larger, water stress increases, causing a decrease in transpiration demand. The supply terms (right-hand side) are increasing functions of absolute leaf water potential. As absolute leaf water potential becomes larger, the gradients in water potential increase, causing an increase in vegetation water supply. PHS takes a Newton's method approach to iteratively solve for the vegetation water potentials that satisfy :eq:`11.21`.

.. math::
   :label: 11.22)

   ff



Yo 





