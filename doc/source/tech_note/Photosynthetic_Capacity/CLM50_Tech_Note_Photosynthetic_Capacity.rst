.. _rst_Photosynthetic Capacity:

Photosynthetic Capacity
=======================

The photosynthetic capacity is represented by two key parameters: 1) the maximum rate of carboxylation at 25 °C, :math:`V_{\text{c,max25}}`; and 2) the maximum rate of electron transport at 25 °C, :math:`J_{\text{max25}}`. They are predicted by a mechanistic model of leaf utilization of nitrogen for assimilation (LUNA V1.0) (:ref:`Ali et al. 2016<Alietal2016>`) based on an optimality hypothesis to nitrogen allocation among light capture, electron transport, carboxylation, respiration and storage. Specifically, the model allocates the nitrogen by maximizing the daily net photosynthetic carbon gain under following two key assumptions:

- nitrogen allocated for light capture, electron transport and carboxylation are co-limiting;

- respiratory nitrogen is allocated to maintain dark respiration determined by :math:`V_{\text{c,max}}`.

Compared to traditional photosynthetic capacity models, a key advantage of LUNA is that the model is able to predict the potential acclimation of photosynthetic capacities at different environmental conditions as determined by temperature, radiation, CO :sub:`2` concentrations, day length, and humidity.

.. _Model inputs and parameter estimations:

Model inputs and parameter estimations
-------------------------------------------------------
The LUNA model includes the following four unitless parameters:

- :math:`J_{maxb0}` , which specifies the baseline proportion of nitrogen allocated for electron transport;
-  :math:`J_{maxb1}` , which determines response of electron transport rate to light availability;
-  :math:`t_{c,j0}` , which defines the baseline ratio of Rubisco-limited rate to light-limited rate;
-  :math:`H` , which determines the response of electron transport rate to relative humidity.

The above four parameters are estimated by fitting the LUNA model to a global compilation of >800 obervations located at different biomes, canopy locations, and time of the year from 1993-2013 (Ali et al. 2015). The model inputs are area-based leaf nitrogen content, leaf mass per unit leaf area and the driving environmental conditions (average of past 10 days) including temperature, CO :sub:`2` concentrations, daily mean and maximum radiation, relative humidity and day length. The estimated values in CLM5 for the listed parameters are 0.0311, 0.17, 0.8054, and 6.0999, repectively. In LUNA V1.0, the estimated parameter values are for C3 natural vegetations. In view that potentially large differences in photosythetic capacity could exist between crops and natural vegetations due to human selection and genetic modifications, in CLM5, the LUNA model are used only for C3 natural vegetations. The photosynthetic capacity for crops and C4 plants are thus still kept the same as CLM4.5. Namely, it is estimated based on the leaf nitrogen content, fixed RUBISCO allocations for :math:`V_{c\max 25}` and an adjusting factor to account for the impact of day length. In CLM5, the model simulates both sun-lit and shaded leaves; however, because the sun-lit and shaded leaves can changes through the day based on the sun angles, we do not differentiate the photosynthetic capacity difference for sun-lit or shaded leaves.

.. _Model structure:

Model structure
----------------------------------------------------------

.. _Plant Nitrogen:

Plant Nitrogen
''''''''''''''''''''''''''

The structure of the LUNA model is adapted from :ref:`Xu et al. (2012)<Xuetal2012>`, where the plant nitrogen at the leaf level ( :math:`\text{LNC}_{a}`; gN/ m :sup:`2` leaf) is divided into four pools: structural nitrogen( :math:`N_{\text{str}}`; gN/m :sup:`2` leaf), photosynthetic nitrogen ( :math:`N_{\text{psn}}`; gN/m :sup:`2` leaf), storage nitrogen( :math:`N_{\text{store}}`; gN/m :sup:`2` leaf), and respiratory nitrogen ( :math:`N_{\text{resp}}`; gN/m :sup:`2` leaf). Namely,

.. math::
  :label: 10.1)

   \text{LNC}_{a} = N_{\text{psn}} + N_{\text{str}}+ N_{\text{store}} + N_{\text{resp}}.

The photosynthetic nitrogen, :math:`N_{\text{psn}}`, is further divided into nitrogen for light capture ( :math:`N_{\text{lc}}`; gN/m :sup:`2` leaf), nitrogen for electron transport ( :math:`N_{\text{et}}`; gN/m :sup:`2` leaf), and nitrogen for carboxylation ( :math:`N_{\text{cb}}`; gN/m :sup:`2` leaf). Namely,

.. math::
  :label: 10.2)

   N_{\text{psn}} =N_{\text{et}} + N_{\text{cb}} + N_{\text{lc}}.

The structural nitrogen, :math:`N_{\text{str}}`, is calculated as the multiplication of leaf mass per unit area (:math:`\text{LMA}`; g biomass/m :sup:`2` leaf), and the structural nitrogen content (:math:`\text{SNC}`; gN/g biomass). Namely,

.. math::
  :label: 10.3)

   N_{\text{str}} = \text{SNC} \cdot \text{LMA}

where :math:`\text{SNC}` is set to be fixed at 0.004 (gN/g biomass), based on data on C:N ratio from dead wood (White etal.,2000), and :math:`\text{LMA}` is the inverse of specific leaf area at the canopy top (:math:`SLA_{\text{0}}`), a PFT-level parameter (:numref:`Table Plant functional type (PFT) leaf N parameters`).

.. _Table Plant functional type (PFT) leaf N parameters:

.. table:: Plant functional type (PFT) leaf N parameters.

 +----------------------------------+--------------------------+
 | PFT                              |  :math:`SLA_{\text{0}}`  |
 +==================================+==========================+
 | NET Temperate                    |          0.01000         |
 +----------------------------------+--------------------------+
 | NET Boreal                       |          0.01000         |
 +----------------------------------+--------------------------+
 | NDT Boreal                       |          0.02018         |
 +----------------------------------+--------------------------+
 | BET Tropical                     |          0.01900         |
 +----------------------------------+--------------------------+
 | BET temperate                    |          0.01900         |
 +----------------------------------+--------------------------+
 | BDT tropical                     |          0.03080         |
 +----------------------------------+--------------------------+
 | BDT temperate                    |          0.03080         |
 +----------------------------------+--------------------------+
 | BDT boreal                       |          0.03080         |
 +----------------------------------+--------------------------+
 | BES temperate                    |          0.01798         |
 +----------------------------------+--------------------------+
 | BDS temperate                    |          0.03072         |
 +----------------------------------+--------------------------+
 | BDS boreal                       |          0.02800         |
 +----------------------------------+--------------------------+
 | C\ :sub:`3` arctic grass         |          0.02100         |
 +----------------------------------+--------------------------+
 | C\ :sub:`3` grass                |          0.04024         |
 +----------------------------------+--------------------------+
 | C\ :sub:`4` grass                |          0.03846         |
 +----------------------------------+--------------------------+
 | Temperate Corn                   |          0.05000         |
 +----------------------------------+--------------------------+
 | Spring Wheat                     |          0.03500         |
 +----------------------------------+--------------------------+
 | Temperate Soybean                |          0.03500         |
 +----------------------------------+--------------------------+
 | Cotton                           |          0.03500         |
 +----------------------------------+--------------------------+
 | Rice                             |          0.03500         |
 +----------------------------------+--------------------------+
 | Sugarcane                        |          0.05000         |
 +----------------------------------+--------------------------+
 | Tropical Corn                    |          0.05000         |
 +----------------------------------+--------------------------+
 | Tropical Soybean                 |          0.03500         |
 +----------------------------------+--------------------------+
 | Miscanthus                       |          0.03500         |
 +----------------------------------+--------------------------+
 | Switchgrass                      |          0.03500         |
 +----------------------------------+--------------------------+

Notes: :math:`SLA_{\text{0}}` is the specific leaf area at the canopy top (m :sup:`2` leaf/g biomass)

We assume that plants optimize their nitrogen allocations (i.e., :math:`N_{\text{store}}`, :math:`N_{\text{resp}}`, :math:`N_{\text{lc}}`, :math:`N_{\text{et}}`, :math:`N_{\text{cb}}`) to maximize the photosynthetic carbon gain, defined as the gross photosynthesis ( :math:`A` ) minus the maintenance respiration for photosynthetic enzymes ( :math:`R_{\text{psn}}` ), under specific environmental conditions and given plant's strategy of leaf nitrogen use. Namely, the solutions of nitrogen allocations \{ :math:`N_{\text{store}}`, :math:`N_{\text{resp}}`, :math:`N_{\text{lc}}`, :math:`N_{\text{et}}`, :math:`N_{\text{cb}}` \} can be estimated as follows,

.. math::
  :label: 10.4)

  \left\{\hat{N}_{\text{{store}}}, \hat{N}_{\text{{resp}}},
    \hat{\mathrm{N}}_{\text{lc}}, \hat{N}_{\text{et}}, \hat{\mathrm{N}}_{\text{cb}}
  \right\} = \underset{\mathrm{N}_{\text{store}}\,+\,\mathrm{N}_{\text{resp}}\,+\,\mathrm{N}_{\text{lc}}\,+\,\mathrm{N}_{\text{et}}\,+\,\mathrm{N}_{\text{cb}}\,<\text{FNC}_{\mathrm{a}}}{\text{argmax}} (A-R_{\text{psn}}),

where :math:`\text{FNC}_{a}` is the functional nitrogen content defined as the total leaf nitrogen content ( :math:`\text{LNC}_{a}`) minus the structural nitrogen content ( :math:`N_{\text{str}}` ).

The gross photosynthesis, :math:`A`, was calculated with a coupled leaf gas exchange model based on the :ref:`Farquhar et al. (1980)<Farquharetal1980>` model of photosynthesis and Ball--Berry-type stomatal conductance model (Ball et al. 1987). The maintenance respiration for photosynthetic enzymes, :math:`R_{\text{psn}}`, is calculated by the multiplication of total photosynthetic nitrogen ( :math:`N_{\text{psn}}` ) and the maintenance respiration cost for photosynthetic enzymes.

Maximum electron transport rate
'''''''''''''''''''''''''''''''''

In the LUNA model, the maximum electron transport rate ( :math:`J_{\text{max}}`; :math:`{\mu} mol` electron / m :sup:`2`/s) is simulated to have a baseline allocation of nitrogen and additional nitrogen allocation to change depending on the average daytime photosynthetic active radiation (PAR; :math:`{\mu} mol` electron / m :sup:`2`/s), day length (hours) and air humidity. Specifically, the LUNA model has

.. math::
  :label: 10.5)

  J_{\text{{max}}} = J_{\text{max}0} + J_{\text{max}b1}
  f\left(\text{day length} \right)f\left(\text{humidity}
  \right)\alpha \text{PAR}

The baseline electron transport rate, :math:`J_{\text{max}0}`, is calculated as follows,

.. math::
  :label: 10.6)

  J_{\text{max}0} = J_{\text{max}b0}{\text{FNC}}_{\mathrm{a}}{\text{NUE}}_{J_{\text{{max}}}}

where :math:`J_{\text{max}b0}` (unitless) is the baseline proportion of nitrogen allocated for electron transport rate. :math:`{\text{NUE}}_{J_{\text{{max}}}}` ( :math:`{\mu} mol` electron /s/g N) is the nitrogen use efficiency of :math:`J_{\text{{max}}}`. :math:`J_{\text{max}b1}` (unitless) is a coefficient determining the response of the electron transport rate to amount of absorbed light (i.e., :math:`\alpha \text{PAR}`). :math:`f\left(\text{day length} \right)` is a function specifies the impact of day length (hours) on :math:`J_{\text{max}}` in view that longer day length has been demonstrated by previous studies to alter :math:`V_{\mathrm{c}\text{max}25}` and :math:`J_{\text{max}25}` (Bauerle et al. 2012; Comstock and Ehleringer 1986) through photoperiod sensing and regulation (e.g., Song et al. 2013). Following Bauerle et al. (2012), :math:`f\left(\text{day length} \right)` is simulated as follows,

.. math::
  :label: 10.7)

  f\left(\text{day length} \right) = \left(\frac{\text{day length}}{12} \right)^{2}.

:math:`f\left(\text{humidity} \right)` represents the impact of air humidity on :math:`J_{\text{{max}}}`. We assume that higher humidity leads to higher :math:`J_{\text{{max}}}` with less water limiation on stomta opening and that low relative humidity has a stronger impact on nitrogen allocation due to greater water limitation. When relative humidity (RH; unitless) is too low, we assume that plants are physiologically unable to reallocate nitrogen. We therefore assume that there exists a critical value of relative humidity ( :math:`RH_{0} = 0.25`; unitless), below which there is no optimal nitrogen allocation. Based on the above assumptions, we have

.. math::
  :label: 10.8)

  f\left(\text{humidity}
  \right) = \left(1-\mathrm{e}^{\left(-H
        \frac{\text{max}\left(\text{RH}-{\text{RH}}_{0}, 0 \right)}{1-\text{RH}_{0}} \right)} \right),

where :math:`H` (unitless) specifies the impact of relative humidity on electron transport rate.

The efficiency of light energy absorption (unitless), :math:`\alpha`, is calculated depending on the amount of nitrogen allocated for light capture, :math:`\mathrm{N}_{\text{lc}}`. Following Niinemets and Tenhunen (1997), the LUNA model has,

.. math::
  :label: 10.9)

  \alpha =\frac{0.292}{1+\frac{0.076}{\mathrm{N}_{\text{lc}}C_{b}}}

where 0.292 is the conversion factor from photon to electron. :math:`C_{b}` is the conversion factor (1.78) from nitrogen to chlorophyll. After we estimate :math:`J_{\text{{max}}}`, the actual electron transport rate with the daily maximum radiation ( :math:`J_{x}`) can be calculated using the empirical expression of leaf (1937),

.. math::
  :label: 10.10)

  J_{x} = \frac{\alpha \text{PAR}_{\text{max}}} {\left(1 + \frac{\alpha^{2}{\text{PAR}}_{\text{{max}}}^{2}}{J_{\text{{max}}}^{2}}
    \right)^{0.5}}

where :math:`\text{PAR}_{\text{{max}}}` ( :math:`\mu mol`/m :sup:`2`/s) is the maximum photosynthetically active radiation during the day.

Maximum rate of carboxylation
''''''''''''''''''''''''''''''

The maximum rate of carboxylation at 25°C varies with foliage nitrogen concentration and specific leaf area and is calculated as in :ref:`Thornton and Zimmermann (2007)<ThorntonZimmermann2007>`. At 25°C,

.. math::
  :label: 10.11)

   V_{c\max 25} = N_{cb} NUE_{V_{c\max 25}}

where :math:`N_{cb}` is nitrogen for carboxylation (g N m\ :sup:`-2` leaf, :numref:`Table Plant functional type (PFT) leaf N parameters`), and :math:`NUE_{V_{c\max 25}}` = 47.3 x 6.25 and is the nitrogen use efficiency for :math:`V_{c\max 25}`. The constant 47.3 is the specific Rubisco activity ( :math:`\mu` mol CO\ :sub:`2` g\ :sup:`-1` Rubisco s\ :sup:`-1`) measured at 25°C, and the constant 6.25 is the nitrogen binding factor for Rubisco (g Rubisco g\ :sup:`-1` N; :ref:`Rogers 2014<Rogers2014>`).

:math:`V_{c\max 25}` additionally varies with daylength (:math:`DYL`) using the function :math:`f(DYL)`, which introduces seasonal variation to :math:`V_{c\max }`

.. math::
  :label: 10.12)

   f\left(DYL\right)=\frac{\left(DYL\right)^{2} }{\left(DYL_{\max } \right)^{2} }

with :math:`0.01\le f\left(DYL\right)\le 1`. Daylength (seconds) is given by

.. math::
  :label: 10.13)

   DYL=2\times 13750.9871\cos ^{-1} \left[\frac{-\sin \left(lat\right)\sin \left(decl\right)}{\cos \left(lat\right)\cos \left(decl\right)} \right]

where :math:`lat` (latitude) and :math:`decl` (declination angle) are from section :numref:`Solar Zenith Angle`. Maximum daylength (:math:`DYL_{\max }` ) is calculated similarly but using the maximum declination angle for present-day orbital geometry (:math:`\pm`\ 23.4667° [:math:`\pm`\ 0.409571 radians], positive for Northern Hemisphere latitudes and negative for Southern Hemisphere).

Implementation of Photosynthetic Capacity
''''''''''''''''''''''''''''''''''''''''''

Based on :ref:`Farquhar et al. (1980)<Farquharetal1980>` and Wullschleger (1993), we can calculate the electron-limited photosynthetic rate under daily maximum radiation ( :math:`W_{jx}`) and the Rubisco-limited photosynthetic rate ( :math:`W_{\mathrm{c}}`) as follows,

.. math::
  :label: 10.14)

  W_{J_{x}} = K_{j}J_{x} ,

.. math::
  :label: 10.15)

  W_{\mathrm{c}} = K_{\mathrm{c}} V_{{\mathrm{c}, \text{max}}},

where :math:`K_{j}` and :math:`K_{\mathrm{c}}` as the conversion factors for :math:`J_{x}` and :math:`V_{{\mathrm{c}, \text{max}}}` ( :math:`V_{{\mathrm{c}, \text{max}}}` to :math:`W_{\mathrm{c}}` and :math:`J_{x}` to :math:`W_{J_{x}}`), respectively. Based on :ref:`Xu et al. (2012)<Xuetal2012>`, Maire et al. (2012) and Walker et al. (2014), we assume that :math:`W_{\mathrm{c}}` is proportional to :math:`W_{J_{x}}`. Specifically, we have

.. math::
  :label: 10.16)

  W_{\mathrm{c}}=t_{\alpha}t_{\mathrm{c}, j0}W_{J_{x}}

where :math:`t_{\mathrm{c}, j0}` is the baseline ratio of :math:`W_{\mathrm{c}}` to :math:`W_{J_{x}}`. We recognize that this ratio may change depending on the nitrogen use efficiency of carboxylation and electron transport (Ainsworth and Rogers 2007), therefore the LUNA model has the modification factor, :math:`t_{\alpha}`, to adjust baseline the ratio depending on the nitrogen use efficiency for electron vs carboxylation (:ref:`Ali et al. 2016<Alietal2016>`).

Total Respiration
'''''''''''''''''''

Following :ref:`Collatz et al. (1991)<Collatzetal1991>`, the total respiration ( :math:`R_{\mathrm{t}}`) is calculated in proportion to :math:`V_{\text{c,max}}`,

.. math::
  :label: 10.17)

  R_{\mathrm{t}} = 0.015 V_{\text{c,max}}.

Accounting for the daytime and nighttime temperature, the daily respirations is calculated as follows,

.. math::
  :label: 10.18)

   R_{\text{td}}={R}_{\mathrm{t}} [D_{\text{day}} + D_{\text{night}} f_{\mathrm{r}}{(T_{\text{night}})/f_{\mathrm{r}}{(T_{\text{day}})}}],

where :math:`D_{\text{day}}` and :math:`D_{\text{night}}` are daytime and nighttime durations in seconds. :math:`f_{\mathrm{r}}(T_{\text{night}})` and :math:`f_{\mathrm{r}}(T_{\text{day}})` are the temperature response functions for respiration (see Appendix B in :ref:`Ali et al. (2016)<Alietal2016>` for details).

.. _Numerical scheme:

Numerical scheme
---------------------------------------------------------

The LUNA model searches for the "optimal" nitrogen allocations for maximum net photosynthetic carbon gain by incrementally increase the nitrogen allocated for light capture (i.e., :math:`N_{\text{lc}}`) (see :ref:`Ali et al. (2016)<Alietal2016>` for details). We assume that plants only optimize the nitrogen allocation when they can grow (i.e., GPP>0.0). If GPP become zero under stress, then the LUNA model assume a certain amount of enzyme will decay at daily rates of 0.1, in view that the half-life time for photosynthetic enzymes are short (~7 days) (Suzuki et al. 2001). To avoid unrealistic low values of photosynthetic capacity, the decay is only limited to 50 percent of the original enzyme levels.
