.. _rst_Photosynthetic Capacity:

Photosynthetic Capacity
=======================

.. _V\ :sub:`cmax25` and Canopy scaling:

V\ :sub:`cmax25` and Canopy scaling
--------------------------------------------

The maximum rate of carboxylation at 25 :sup:`o`\ C varies with
foliage nitrogen concentration and specific leaf area and is calculated
as in Thornton and Zimmermann (2007). At 25ºC,

.. math::
      :label: ZEqnNum217783 

   V_{c\max 25} =N_{a} F_{LNR} F_{NR} a_{R25}

where :math:`N_{a}`  is the area-based leaf nitrogen concentration (g N
m\ :sup:`-2` leaf area), :math:`F_{LNR}`  is the fraction of leaf
nitrogen in Rubisco (g N in Rubisco g\ :sup:`-1` N),
:math:`F_{NR} =7.16` is the mass ratio of total Rubisco molecular mass
to nitrogen in Rubisco (g Rubisco g\ :sup:`-1` N in Rubisco), and
:math:`a_{R25} =60` is the specific activity of Rubisco (µmol
CO\ :sub:`2` g\ :sup:`-1` Rubisco s\ :sup:`-1`).
:math:`N_{a}`  is calculated from mass-based leaf N concentration and
specific leaf area

.. math::
      :label: ZEqnNum561340 

   N_{a} =\frac{1}{CN_{L} \; SLA_{0} }

where :math:`CN_{L}`  is the leaf carbon-to-nitrogen ratio (g C
g\ :sup:`-1` N) and :math:`SLA_{0}`  is specific leaf area at the
canopy top ( m\ :sup:`2` leaf area g\ :sup:`-1` C ). Table 8.1
lists values of :math:`F_{LNR}` , :math:`CN_{L}` , and :math:`SLA_{0}` 
for each plant functional type. :math:`F_{LNR}`  was chosen to give
:math:`V_{c\max 25}`  consistent with Kattge et al. (2009), as discussed
by Bonan et al. (2011, 2012). Table 8.1 lists derived values for
:math:`V_{c\max 25}`  at the top of the canopy using :math:`SLA_{0}` .
Tropical broadleaf evergreen trees are an exception, and a higher
:math:`V_{c\max 25}`  is used to alleviate model biases (Bonan et al.
2012).

:math:`V_{c\max 25}`  is calculated separately for sunlit and shaded
leaves using an exponential profile to area-based leaf nitrogen
(:math:`N_{a}` ), as in Bonan et al. (2011). :math:`V_{c\max 25}`  at
cumulative leaf area index :math:`x` from the canopy top scales directly
with :math:`N_{a}` , which decreases exponentially with greater
cumulative leaf area, so that

.. math::
      :label: ZEqnNum745439 

   V_{c\; \max 25}^{} \left(x\right)=V_{c\; \max 25}^{} \left(0\right)e^{-K_{n} x}

where :math:`V_{c\; \max 25}^{} \left(0\right)` is defined at the top of
the canopy using :math:`SLA_{0}` , and :math:`K_{n}`  is the decay
coefficient for nitrogen. The canopy integrated value for sunlit and
shaded leaves is

.. math::
      :label: 9.20) 

   \begin{array}{rcl} {V_{c\; \max 25}^{sun} } & {=} & {\int _{0}^{L}V_{c\; \max 25}^{} \left(x\right)f_{sun} \left(x\right)\,  dx} \\ {} & {=} & {V_{c\; \max 25}^{} \left(0\right)\left[1-e^{-\left(K_{n} +K\right)L} \right]\frac{1}{K_{n} +K} } \end{array}

.. math::
      :label: 9.21) 

   \begin{array}{rcl} {V_{c\; \max 25}^{sha} } & {=} & {\int _{0}^{L}V_{c\; \max 25}^{} \left(x\right)\left[1-f_{sun} \left(x\right)\right] \, dx} \\ {} & {=} & {V_{c\; \max 25}^{} \left(0\right)\left\{\left[1-e^{-K_{n} L} \right]\frac{1}{K_{n} } -\left[1-e^{-\left(K_{n} +K\right)L} \right]\frac{1}{K_{n} +K} \right\}} \end{array}

and the average value for the sunlit and shaded leaves is

.. math::
      :label: 9.22) 

   \bar{V}_{c\; \max 25}^{sun} ={V_{c\; \max 25}^{sun} \mathord{\left/ {\vphantom {V_{c\; \max 25}^{sun}  L^{sun} }} \right. \kern-\nulldelimiterspace} L^{sun} }

.. math::
      :label: 9.23) 

   \bar{V}_{c\; \max 25}^{sha} ={V_{c\; \max 25}^{sha} \mathord{\left/ {\vphantom {V_{c\; \max 25}^{sha}  L^{sha} }} \right. \kern-\nulldelimiterspace} L^{sha} } .

This integration is over all leaf area (:math:`L`) with
:math:`f_{sun} (x)=\exp \left(-Kx\right)` and :math:`K` the direct beam
extinction coefficient (equation 4.9). Photosynthetic parameters
:math:`J_{\max 25}` , :math:`T_{p25}` , :math:`k_{p25}` , and
:math:`R_{d25}`  scale similarly.

The value :math:`K_{n} = 0.11` chosen by Bonan et al. (2011) is
consistent with observationally-derived estimates for forests, mostly
tropical, and provides a gradient in V\ :sub:`cmax` similar to
the original CLM4 specific leaf area scaling. However, Bonan et al.
(2012) showed that the sunlit/shaded canopy parameterization does not
match an explicit multi-layer canopy parameterization. The discrepancy
arises from absorption of scattered radiation by shaded leaves and can
be tuned out with higher :math:`K_{n}` . The model uses
:math:`K_{n} =0.30` to match an explicit multi-layer canopy.

:math:`V_{c\max 25}`  additionally varies with daylength (:math:`DYL`)
using the function :math:`f(DYL)`, which introduces seasonal variation
to :math:`V_{c\max }` 

.. math::
      :label: 9.24) 

   f\left(DYL\right)=\frac{\left(DYL\right)^{2} }{\left(DYL_{\max } \right)^{2} }

with :math:`0.01\le f\left(DYL\right)\le 1`. Daylength (seconds) is
given by

.. math::
      :label: 9.25) 

   DYL=2\times 13750.9871\cos ^{-1} \left[\frac{-\sin \left(lat\right)\sin \left(decl\right)}{\cos \left(lat\right)\cos \left(decl\right)} \right]

where :math:`lat` (latitude) and :math:`decl` (declination angle) are
from section 3.3. Maximum daylength (:math:`DYL_{\max }` ) is calculated
similarly but using the maximum declination angle for present-day
orbital geometry (:math:`\pm`\ 23.4667º [:math:`\pm`\ 0.409571 radians],
positive for Northern Hemisphere latitudes and negative for Southern
Hemisphere).
