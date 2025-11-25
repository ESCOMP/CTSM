.. _rst_Dust Model:

Dust Model
==============

Atmospheric dust is mobilized from the land by wind in the CLM. The most important factors determining soil erodibility and dust emission include the wind friction velocity, the vegetation cover, and the soil moisture. The latest CTSM allows users to choose between two dust emission schemes: One is Leung_2023 (Leung et al. 2023) which is the default for the CLM6 physics, and the other is Zender_2003 (:ref:`Mahowald et al. 2006<Mahowaldetal2006>`) based on the DEAD (Dust Entrainment and Deposition model of :ref:`Zender et al. (2003)<Zenderetal2003>`.

One can control the use of the dust emission scheme by setting the namelist variable: 

**dust_emis_method = 'Leung_2023'**, which is the default for CLM6 physics or later, or 

**dust_emis_method = 'Zender_2003'**, which is the default for CLM5 physics or older.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Leung et al. (2023) scheme (default for CLM6 physics or later):
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Dust emission modeling is a threshold parameterization of an aeolian (wind-driven) process. For Leung_2023, in any given timestep the soil surface wind velocity :math:`u_{*s}` has to be greater than the fluid threshold friction velocity :math:`u_{*ft}` to generate saltation and dust emission:

.. math::
   :label: 30.0

   u_{*ft} = u_{*ft0}(D_{p},\rho_{a}) f_{m}(w)

where :math:`u_{*ft0}(D_{p},\rho_{a})` is the dry fluid threshold without the soil moisture effect :math:`f_{m}`, as a function of median soil diameter :math:`D_{p}` and air density :math:`\rho_{a}`. In CTSM for Leung_2023, :math:`D_{p}` is a globally uniform number of 130 :math:`\mu` m. :math:`u_{*ft0}(D_{p},\rho_{a})` is given by Shao and Lu (2000):

.. math::
   :label: 30.01

   u_{*ft0}(D_{p},\rho_{a}) = \sqrt{\frac{A(\rho_{p} g D_{p} + \gamma / D_{p}) }{\rho_{a}} } 


where :math:`g = 9.81` m s\ :sup:`-2` is gravity, :math:`\rho_{p} = 2650` kg m\ :sup:`-3` is typical soil particle density, and :math:`A = 0.0123` and :math:`\gamma = 1.65 \times 10^{-4}` kg s\ :sup:`-2` are empirical constants. 
The soil moisture effect :math:`f_{m}(w)` is a function of gravimetric soil moisture :math:`w` (kg water / kg soil) at the topmost soil layer. 

:math:`w` is converted from the CLM volumetric soil moisture :math:`\theta` and porosity (saturation moisture :math:`\theta_{sat}`):

.. math::
   :label: 30.02

   w=\theta\frac{ \rho _{water} }{\rho_{bulk} }

where :math:`\rho_{water} = 1000 kg\ m^{-3}` is typical water density, and bulk density :math:`\rho_{bulk}` is given by 

.. math::
   :label: 30.03

   \rho_{bulk} = (1 - \theta_{sat} ) \rho_{p}

Then, the soil moisture effect :math:`f_{m}(w)` on increasing the fluid threshold is given by Fecan et al. (1999):

.. math::
   :label: 30.04

   f_{m}(w) =\left\{\begin{array}{l} {1{\rm \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; for\; }w\le w_{t} } \\ {\sqrt{1+1.21\left[100\left(w-w_{t} \right)\right]^{0.68} } {\rm \; \; for\; }w>w_{t} } \end{array}\right.


where :math:`w_{t}` increases with clay fraction :math:`f_{clay}`:

.. math::
   :label: 30.05

   w_{t} =0.01a\left(17f_{clay} +14f_{clay}^{2} \right){\rm \; \; \; \; \; \; 0}\le f_{clay} =\% clay\times 0.01\le 1


where :math:`a=f_{clay}^{-1}` for tuning purposes. Note that this is different from the paper in which a = 1 was chosen. The coefficient 0.01 is used for converting :math:`w_{t}` from % to fraction  (kg water / kg soil), :math:`\theta` is the volumetric soil moisture in the topmost soil layer (m :math:`{}^{3 }`\ water \ m\ :sup:`-3` soil) (section :numref:`Soil Water`), :math:`\rho _{liq}` is the density of liquid water (kg m\ :sup:`-3`) (:numref:`Table Physical constants`), and :math:`\rho _{bulk}` is the bulk density of soil in the top soil layer (kg m\ :sup:`-3`) defined as in section :numref:`Soil and Snow Thermal Properties` rather than as in :ref:`Zender et al. (2003)<Zenderetal2003>`. :math:`f_{clay}` is the mass fraction of clay particles in the topmost soil layer and %clay comes from the surface dataset (section :numref:`Surface Data`).

Another essential dust emission threshold is the impact / dynamic threshold :math:`u_{*it}`, which is the lowest friction velocity or wind stress to matintain saltation:

.. math::
   :label: 30.06

   u_{*it} = B_{*it} u_{*ft0}

where :math:`u_{*it}` is a constant on Earth following Kok et al. (2012). The above equations imply that :math:`u_{*it} \, < \, u_{*ft0} \le \, u_{*ft}`. This means that the winds need a bigger momentum to initiate saltation and dust emission but can reduce below :math:`u_{*ft}` and still maintain a weak dust emission flux. The emission flux goes to zero when :math:`u_{*s}` drops below :math:`u_{*it}`.

The total vertical mass emission flux of dust, :math:`F_{d}` (kg m\ :sup:`-2` s\ :sup:`-1`), from the ground into a transport mode/bin :math:`j` of aerosol is given by

.. math::
   :label: 30.07

   F_{d} = C_{tune} C_{d} f_{bare} f_{clay} \frac{ \rho_{a} (u^2_{*s} - u^2_{it} ) }{ u^2_{it} }  \left( \frac{ u^2_{*s} }{u^2_{it} } \right) ^\kappa

where :math:`C_{tune} = 0.05` is a constant, and :math:`F_{d}` is the total emission flux summed across modes/bins following a revised form of Kok et al. (2014). The dust emission flux goes to zero when :math:`u_{*s} \, < \, u_{*it}`. :math:`\rho_{a}` is surface air density from CAM (the atm model). 
:math:`\kappa` is fragmentation exponent, and :math:`C_{d}` is the dust emission coefficient (or the soil erodibility coefficient):

.. math::
   :label: 30.1

   C_{d} = C_{d0} \exp{ (-C_{e} \frac{ u_{*st} - u_{st0} }{ u_{st0} } ) }

.. math::
   :label: 30.2

   \kappa = C_{\kappa} \frac{ u_{*st} - u_{st0} }{ u_{st0} }

where :math:`C_{\kappa} = 2.7`, :math:`u_{st0} = 0.16 m s^{-1}`, :math:`C_{d0} = 4.4 \times 10^{-5}`, and :math:`C_{e} = 2.0`. :math:`F_{d}` thus scales with :math:`u^{2+\kappa}_{*s}`, where :math:`\kappa \sim 1` over major deserts and :math:`\sim 3` or higher over semiarid and nonarid regions. Since Kok et al. (2014) has not measured :math:`\kappa > 3` in their measurements, we will cap :math:`\kappa` at a maximum value (currently set as 2.5). :math:`u_{*st}` is the standardized wet fluid threshold at a typical atmospheric surface air density (Kok et al., 2014):

.. math::
   :label: 30.3

   u_{*st} = u_{*ft} \sqrt{ \rho_{a} / \rho_{0a}}

where :math:`\rho_{0a} = 1.225` kg m\ :sup:`-3`. As can be seen, :math:`u_{*st}` scales with :math:`u_{*ft}` and thus soil moisture :math:`w`. Therefore, moisture :math:`w` decreases soil erodibility :math:`C_{d}` but increases dust emission sensitivity :math:`\kappa` to the winds.

Kok et al. (2014) is different from most other dust emission parameterizations in the way that the soil erodibility :math:`C_{d}` is not a time-invariant input data but is a transient function, with erodibility increasing with reducing :math:`u^2_{ft}` (and thus implicitly soil moisture). Similarly, the fragmentation exponent :math:`\kappa` is also transient but increases with enhancing soil moisture.

The grid cell fraction of exposed bare soil suitable for dust mobilization :math:`f_{bare}` is given by

.. math::
   :label: 30.4

   f_{bare} =\left(1-f_{lake} \right)\left(1-f_{sno} \right)\left(1-f_{v} \right)\frac{w_{liq,1} }{w_{liq,1} +w_{ice,1} }

where :math:`f_{lake}` and :math:`f_{sno}` are the CLM grid cell fractions of lake (section :numref:`Surface Data`) and snow cover (section :numref:`Snow Covered Area Fraction`), all ranging from zero to one. Not mentioned by :ref:`Zender et al. (2003)<Zenderetal2003>`, :math:`w_{liq,\, 1}` and :math:`{}_{w_{ice,\, 1} }` are the CLM top soil layer liquid water and ice contents (mm) entered as a ratio expressing the decreasing ability of dust to mobilize from increasingly frozen soil. The grid cell fraction of vegetation cover,\ :math:`f_{v}`, is defined as

.. math::
   :label: 30.5

   0\le f_{v} =\frac{\mathrm{VAI}}{\mathrm{VAI_{thr}} } \le 1{\rm \; \; \; \; where\; } \mathrm{VAI_{thr}} =0.6{\rm \; m}^{2} {\rm m}^{-2}

where :math:`\mathrm{VAI}=\mathrm{LAI}+\mathrm{SAI}` is the vegetation area index as a sum of the CLM leaf and stem area index values (m :sup:`2` m\ :sup:`-2`) averaged at the land unit level so as to include all the pfts and the bare ground present in a vegetated land unit. :math:`\mathrm{LAI}` and :math:`\mathrm{SAI}` may be prescribed from the CLM input data (section :numref:`Phenology and vegetation burial by snow`) or simulated by the CLM biogeochemistry model (Chapter :numref:`rst_Vegetation Phenology and Turnover`).

On top of Kok et al. (2014), Leung_2023 introduced the soil surface friction velocity :math:`u_{*s}` as the friction velocity :math:`u_{*}` from CLM corrected by the surface roughness due to the presented rocks and vegetation on the soil surface, encapsulated by the so-called drag partition factor :math:`F_{eff}`.

.. math::
   :label: 30.6

   u_{*s} = u_{*} F_{eff}

The Leung et al. (2023) paper uses an area-weighted averaging method to determine the mean drag partitioning for a grid cell: :math:`F_{eff}^3 = A_{rock} f_{rock}^3 + A_{veg} f_{veg}^3`. :math:`A_{rock}` and :math:`A_{veg}` are fractional area cover (in fraction) from the CLM-prescribed land use from the Land Use Harmonization 2 (LUH2; Hurtt et al., 2020). :math:`F_{eff}` is thus a weighted mean of the rock drag partitioning and the vegetation drag partitioning in Leung et al. (2023). However, since CTSM has the privilege of supporting sub-grid patch-level simulations of dust emissions, we simply separate the calculations of dust emissions into the areas of bare soils and areas of the short vegetation. For a bare soil patch/PFT we use:

.. math::
   :label: 30.7

    F_{eff} = f_{rock}

And for a patch/PFT with short vegetation (shrub, grass, crop) we use:

.. math::
   :label: 30.8

    F_{eff} = f_{veg} 

The rock drag partition factor scales with the surface roughness density of rocks, captured by the aeolian roughness length :math:`z_{0a}` from the satellite-derived dataset from Prigent et al. (2005). The expression was first developed by Marticorena and Bergametti (1995), which is more valid for the low-roughness surfaces:

.. math::
   :label: 30.9

   f_{rock} = 1 - \frac{\ln(\frac{z_{0a}}{z_{0s}})}{\ln(b_{1}[\frac{X}{z_{0s}}]^{b_{2}})}

where :math:`X = 10 m` is the distance downstream the point of discontinuity in surface obstacle, :math:`b_{1} = 0.7` and :math:`b_{2} = 0.8` are coefficients (Darmenova et al., 2009), :math:`z_{0s}` is the soil roughness length. :math:`z_{0a}` is from Prigent et al. (2005) and should not be confused with the aerodynamic roughness length :math:`Z_{0}` from the model. This equation only applies for gridcells with LAI smaller than the LAI threshold for dust emission.

The vegetation drag partition factor scales with vegetation density as captured by LAI following Okin (2008) and Pierre et al. (2014):

.. math::
   :label: 30.10

   f_{veg} = \frac{K + f_{0} c }{K + c}

.. math::
   :label: 30.11

   K = 2 (\frac{1}{f_{v}} - 1) = 2 (\frac{1}{\mathrm{VAI}/\mathrm{VAI_{thr}}} - 1)

where :math:`f_{0} = 0.33` is the wind dissipation right behind the obstacle, :math:`c = 4.8 m` is the average e-folding distance for wind momentum to recover after dissipated by the plant, and :math:`K` is defined as the average normalized gap length between two obstacles given the vegetation density quantified by VAI. This equation assumes that the gap length K goes to zero as VAI increases and approaches :math:`\mathrm{VAI_{thr}}`.

Then, the fraction of time within a timestep :math:`\Delta t` when dust emission is active is accounted for by Leung_2023, known as the intermittency factor :math:`\eta`. The fraction of time :math:`\eta \in [0,1]` is due to the boundary-layer turbulence that causes sub-timestep turbulent wind fluctuations that cross the dust emission thresholds multiple times within a timestep, in turn initiating and shutting off emission at time within :math:`\Delta t`. The parameterization of :math:`\eta` comes from Comola et al. (2019), in which a statistical substepping method was proposed to account for the temporary shutoff of dust emission fluxes. 

The fraction of time :math:`\eta` is parameterized using the surface winds and thresholds at the saltation height. Therefore, the friction velocities are translated using the log law of the wall to the saltation height, which was defined as :math:`z_{sal}` = 0.1 m by Comola et al. (2019):

.. math::
   :label: 30.12

   u_{s}(z_{sal}) =  \frac{u_{*s}}{k} \ln(z_{sal}/z_{0a})

.. math::
   :label: 30.13

   u_{ft}(z_{sal}) =  \frac{u_{*ft}}{k} \ln(z_{sal}/z_{0a})

.. math::
   :label: 30.14

   u_{it}(z_{sal}) =  \frac{u_{*it}}{k} \ln(z_{sal}/z_{0a})


With saltation-height variables defined, the instantaneous wind :math:`\tilde{u}_s` is assumed by Comola to follow a Gaussian distribution with a mean equal to the mean wind speed and the spread :math:`\sigma_{u_{s}}`  parameterized by the Similarity Theory (Panofsky, 1979):

.. math::
   :label: 30.15

   \tilde{u}_s \sim N(u_s, \sigma_{u_s}) 

And the fluctuation strength is parameterized by the similarity theory:

.. math::
   :label: 30.16

   \sigma_{u_s} = u_*^s \left( 12 - 0.5 \frac{z_i}{L} \right)^{1/3}
   \quad \text{for } \left( 12 - 0.5 \frac{z_i}{L} \right) \ge 0

where :math:`z_i = 1000` m is the planetary boundary-layer height set as a constant for now, and :math:`L` is the Obukhov length. This means the instantaneous wind's fluctuation comes from both a shear contribution and a buoyancy contribution. 


Then, the total fraction of time :math:`\eta` when saltation is active within a model timestep is then formulated as

.. math::
   :label: 30.17

   \eta = 1 - P_{ft} + \alpha \left( P_{ft} - P_{it} \right)

where :math:`P_{it}` is the cumulative probability that the instantaneous wind :math:`\tilde{u}_s` does not exceed the impact threshold :math:`u_{it}`, and :math:`P_{ft}` is the cumulative probability that :math:`\tilde{u}_s` does not exceed the fluid threshold :math:`u_{ft}`. The fluid threshold crossing fraction :math:`\alpha` is defined as the rate of :math:`\tilde{u}_s` sweeping across :math:`u_{ft}` divided by the rate of sweeping across :math:`u_{it}` and :math:`u_{ft}` summed up. For instance, if :math:`\tilde{u}_s` sweeps across :math:`u_{ft}` often and does not sweep across :math:`u_{it}` much, it means that :math:`\tilde{u}_s` (and thus the timestep-mean :math:`u_s`) should be closer to :math:`u_{ft}` and generally higher than :math:`u_{it}`. Then, :math:`\alpha` is close to 1, and the fraction of time :math:`\eta` (with active emission within a timestep) should also be close to 1. :math:`\alpha` can be represented as

.. math::
   :label: 30.18

   \alpha \approx \left\{ \exp\left[ 
      \frac{u_{ft}^2 - u_{it}^2 - 2u_s(u_{ft}-u_{it})}{2 \sigma^2_{u_s} }
      \right] + 1 \right\}^{-1}

Then, the fraction of time in :math:`\Delta t` when :math:`\tilde{u}_s` is above :math:`u_{ft}` is given by :math:`1 - P_{ft}`, where

.. math::
   :label: 30.19

   P_{ft} = \frac{1}{2} \left[ 1 + \operatorname{erf}
   \left( \frac{u_{ft} - u_s}{\sqrt{2} \sigma_{u_s}} \right) \right]

And, the fraction of time in :math:`\Delta t` when :math:`\tilde{u}_s` is dropping from above :math:`u_{ft}` to within the hysteresis zone between :math:`u_{it}` and :math:`u_{it}` is given by :math:`\alpha \left( P_{ft} - P_{it} \right)`, where
 
.. math::
   :label: 30.20

   P_{it} = \frac{1}{2} \left[ 1 + \operatorname{erf}
   \left( \frac{u_{it} - u_s}{\sqrt{2} \sigma_{u_s}} \right) \right]

And so the fraction of time :math:`\eta` within :math:`\Delta t` with active emission is determined. Dust emission will be inactive when :math:`\tilde{u}_s` is below the impact threshold :math:`u_{it}` or when :math:`\tilde{u}_s` is rising from below :math:`u_{it}` to the hysteresis zone between :math:`u_{it}` and :math:`u_{it}`.


The total vertical mass emission flux of dust, :math:`F_{d}` (kg m\ :sup:`-2` s\ :sup:`-1`) is then computed with all the above equations and terms. The emission flux is then passed through the coupler to the atmospheric model (CAM) to simulate dust aerosol transport and deposition. The default aerosol model supported in the CAM6 and CAM7 physics is the Modal Aerosol Model (MAM). CAM7 uses the 5-mode MAM (MAM5), in which three modes (Aitken, accumulation, and coarse modes) by default contain dust. The total mass emission flux per grid is partitioned into the three modes following the Brittle Fragmentation Theory (BFT) by Kok et al. (2014) and later modified by Meng et al. (2022). In the current model version, the fractions of dust emission flux partitioned in the three modes are 1.65E-05, 0.021, and 0.979 for the Aitken (0.01–0.1 um), accumulation (0.1–1 um), and coarse (1–10 um) modes, respectively. These values are prescribed in the MAM code inside CAM.

In the early CESM versions, CAM employed the Bulk Aerosol Model (BAM) as the default aerosol model. Thus, the emission fluxes in CTSM is by default partitioned into 4 bins before passing to the coupler:

.. math::
   :label: 30.21

   F_{j} = F_d \sum _{i=1}^{I}M_{i,j}

where :math:`F_{j}` is the mass emission flux from the :math:`j` th aerosol bin. The bottom of this page shows how :math:`M_{i,j}` is calculated.

The current way of paritioning the emission fluxes before passing to the coupler is still being used, but the partition into the different bins is only required by BAM, not the current default MAM in CAM6 and CAM7. Therefore, for MAM, the four :math:`F_{j}` are summed up to become one total flux :math:`F_d` again. It is then redistributed to different individual modes of MAM following the Brittle Fragmentation Theory. The default CAM physics that is paired with the CLM6 physics is the CAM7 physics, which uses MAM4 (for minimal atmospheric chemistry) or MAM5 (for complex chemistry).

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
**Zender et al. (2003) scheme (default for CLM5 physics and older):**
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The total vertical mass flux of dust, :math:`F_{j}` (kg m\ :sup:`-2` s\ :sup:`-1`), from the ground into transport bin :math:`j` is given by

.. math::
   :label: 30.22

   F_{j} =TSf_{m} \alpha Q_{s} \sum _{i=1}^{I}M_{i,j}

where :math:`T` is a global factor that compensates for the DEAD model's sensitivity to horizontal and temporal resolution and equals 5 x 10\ :sup:`-4` in the CLM instead of 7 x 10\ :sup:`-4` in :ref:`Zender et al. (2003)<Zenderetal2003>`. :math:`S` is the source erodibility factor set to 1 in the CLM and serves as a place holder at this time.

The grid cell fraction of exposed bare soil suitable for dust mobilization :math:`f_{m}` is given by

.. math::
   :label: 30.23

   f_{m} =\left(1-f_{lake} \right)\left(1-f_{sno} \right)\left(1-f_{v} \right)\frac{w_{liq,1} }{w_{liq,1} +w_{ice,1} }

where :math:`f_{lake}` and :math:`f_{sno}` are the CLM grid cell fractions of lake (section :numref:`Surface Data`) and snow cover (section :numref:`Snow Covered Area Fraction`), all ranging from zero to one. Not mentioned by :ref:`Zender et al. (2003)<Zenderetal2003>`, :math:`w_{liq,\, 1}` and :math:`{}_{w_{ice,\, 1} }` are the CLM top soil layer liquid water and ice contents (mm) entered as a ratio expressing the decreasing ability of dust to mobilize from increasingly frozen soil. The grid cell fraction of vegetation cover,\ :math:`{}_{f_{v} }`, is defined as

.. math::
   :label: 30.24

   0\le f_{v} =\frac{L+S}{\left(L+S\right)_{t} } \le 1{\rm \; \; \; \; where\; }\left(L+S\right)_{t} =0.3{\rm \; m}^{2} {\rm m}^{-2}

where equation :eq:`30.3` applies only for dust mobilization and is not related to the plant functional type fractions prescribed from the CLM input data or simulated by the CLM dynamic vegetation model (Chapter 22). :math:`L` and :math:`S` are the CLM leaf and stem area index values (m :sup:`2` m\ :sup:`-2`) averaged at the land unit level so as to include all the pfts and the bare ground present in a vegetated land unit. :math:`L` and :math:`S` may be prescribed from the CLM input data (section :numref:`Phenology and vegetation burial by snow`) or simulated by the CLM biogeochemistry model (Chapter :numref:`rst_Vegetation Phenology and Turnover`).

The sandblasting mass efficiency :math:`\alpha` (m :sup:`-1`) is calculated as

.. math::
   :label: 30.25

   \alpha =100e^{\left(13.4f_{clay} -6.0\right)\ln 10} {\rm \; \; }\left\{\begin{array}{l} {f_{clay} =\% clay\times 0.01{\rm \; \; \; 0}\le \% clay\le 20} \\ {f_{clay} =20\times 0.01{\rm \; \; \; \; \; \; \; \; 20<\% }clay\le 100} \end{array}\right.

where :math:`f_{clay}` is the mass fraction of clay particles in the soil and %clay is determined from the surface dataset (section :numref:`Surface Data`). :math:`f_{clay} =0` corresponds to sand and :math:`f_{clay} =0.2` to sandy loam.

:math:`Q_{s}` is the total horizontally saltating mass flux (kg m\ :sup:`-1` s\ :sup:`-1`) of "large" particles (:numref:`Table Dust Mass fraction`), also referred to as the vertically integrated streamwise mass flux

.. math::
   :label: 30.26

   Q_{s} = \left\{
   \begin{array}{lr}
   \frac{c_{s} \rho _{atm} u_{*s}^{3} }{g} \left(1-\frac{u_{*t} }{u_{*s} } \right)\left(1+\frac{u_{*t} }{u_{*s} } \right)^{2} {\rm \; } & \qquad {\rm for\; }u_{*t} <u_{*s}  \\
   0{\rm \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; } & \qquad {\rm for\; }u_{*t} \ge u_{*s}
   \end{array}\right.

where the saltation constant :math:`c_{s}` equals 2.61 and :math:`\rho _{atm}` is the atmospheric density (kg m\ :sup:`-3`) (:numref:`Table Atmospheric input to land model`), :math:`g` the acceleration of gravity (m s\ :sup:`-2`) (:numref:`Table Physical constants`). The threshold friction speed for saltation :math:`u_{*t}` (m s\ :sup:`-1`) is technically the fluid threshold :math:`u_{*ft}` referred in the Leung_2023 scheme, since Zender_2003 does not consider the impact threshold :math:`u_{*it}`. :math:`u_{*t}` is formulated the same way as in Leung_2023, but with the dry threshold :math:`u_{*ft0}(D_{p},\rho_{a})` using the Marticorena and Bergametti (1995) formulation instead of the Shao and Lu (2000) one:

.. math::
   :label: 30.27

   u_{*t} = u_{*ft0}(D_{p},\rho_{a}) f_{m}(w) = \left[Re_{*t}^{f} \rho _{osp} gD_{osp} \left(1+\frac{6\times 10^{-7} }{\rho _{osp} gD_{osp}^{2.5} } \right)\right]^{\frac{1}{2} } \rho _{atm} ^{-\frac{1}{2} } f_{w}

where :math:`\rho _{osp}` and :math:`D_{osp}` are the density (2650 kg m\ :sup:`-3`) and diameter (75 x 10\ :math:`{}^{-6}` m) of optimal saltation particles, and :math:`f_{w}` is a factor dependent on soil moisture. The following formulas are the same as the ones for Leung_2023 above:

.. math::
   :label: 30.28

   f_{w} =\left\{\begin{array}{l} {1{\rm \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; for\; }w\le w_{t} } \\ {\sqrt{1+1.21\left[100\left(w-w_{t} \right)\right]^{0.68} } {\rm \; \; for\; }w>w_{t} } \end{array}\right.

where

.. math::
   :label: 30.29

   w_{t} =a\left(0.17f_{clay} +0.14f_{clay}^{2} \right){\rm \; \; \; \; \; \; 0}\le f_{clay} =\% clay\times 0.01\le 1

and

.. math::
   :label: 30.30

   w=\frac{\theta _{1} \rho _{liq} }{\rho _{d,1} }

where :math:`a=f_{clay}^{-1}` for tuning purposes, :math:`\theta _{1}` is the volumetric soil moisture in the top soil layer (m :math:`{}^{3 }`\ m\ :sup:`-3`) (section :numref:`Soil Water`), :math:`\rho _{liq}` is the density of liquid water (kg m\ :sup:`-3`) (:numref:`Table Physical constants`), and :math:`\rho _{d,\, 1}` is the bulk density of soil in the top soil layer (kg m\ :sup:`-3`) defined as in section :numref:`Soil and Snow Thermal Properties` rather than as in :ref:`Zender et al. (2003)<Zenderetal2003>`. :math:`Re_{*t}^{f}` from equation :eq:`30.6` is the threshold friction Reynolds factor

.. math::
   :label: 30.31

   Re_{*t}^{f} =\left\{\begin{array}{l} {\frac{0.1291^{2} }{-1+1.928Re_{*t} } {\rm \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; for\; 0.03}\le Re_{*t} \le 10} \\ {0.12^{2} \left(1-0.0858e^{-0.0617(Re_{*t} -10)} \right)^{2} {\rm \; for\; }Re_{*t} >10} \end{array}\right.

and :math:`Re_{*t}` is the threshold friction Reynolds number approximation for optimally sized particles

.. math::
   :label: 30.32

   Re_{*t} =0.38+1331\left(100D_{osp} \right)^{1.56}

In :eq:`30.5`, :math:`u_{*s}` is defined as the wind friction speed (m s\ :sup:`-1`) accounting for the Owen effect (:ref:`Owen 1964<Owen1964>`)

.. math::
   :label: 30.33

   u_{*s} = \left\{
   \begin{array}{lr}
   u_{*} & \quad {\rm \; for \;} U_{10} <U_{10,t}  \\
   u_{*} +0.003\left(U_{10} -U_{10,t} \right)^{2} & \quad {\rm \; for\; }U_{10} \ge U_{10,t}
   \end{array}\right.

where :math:`u_{*}` is the CLM wind friction speed (m s\ :sup:`-1`), also known as friction velocity (section :numref:`Monin-Obukhov Similarity Theory`), :math:`U_{10}` \ is the 10-m wind speed (m s\ :sup:`-1`) calculated as the wind speed at the top of the canopy in section 4.3 of :ref:`Bonan (1996)<Bonan1996>` but here for 10 m above the ground, and :math:`U_{10,\, t}` is the threshold wind speed at 10 m (m s\ :sup:`-1`)

.. math::
   :label: 30.34

   U_{10,t} =u_{*t} \frac{U_{10} }{u_{*} }

In equation :eq:`30.1` we sum :math:`M_{i,\, j}` over :math:`I=3` source modes :math:`i` where :math:`M_{i,\, j}` is the mass fraction of each source mode :math:`i` carried in each of *:math:`J=4`* transport bins :math:`j`

.. math::
   :label: 30.35

   M_{i,j} =\frac{m_{i} }{2} \left[{\rm erf}\left(\frac{\ln {\textstyle\frac{D_{j,\max } }{\tilde{D}_{v,i} }} }{\sqrt{2} \ln \sigma _{g,i} } \right)-{\rm erf}\left(\frac{\ln {\textstyle\frac{D_{j,\min } }{\tilde{D}_{v,i} }} }{\sqrt{2} \ln \sigma _{g,i} } \right)\right]

where :math:`m_{i}`, :math:`\tilde{D}_{v,\, i}`, and :math:`\sigma _{g,\, i}` are the mass fraction, mass median diameter, and geometric standard deviation assigned to each particle source mode :math:`i` (:numref:`Table Dust Mass fraction`), while :math:`D_{j,\, \min }` and :math:`D_{j,\, \max }` are the minimum and maximum diameters (m) in each transport bin :math:`j` (:numref:`Table Dust Minimum and maximum particle diameters`).

.. _Table Dust Mass fraction:

.. table:: Mass fraction :math:`m_{i}` , mass median diameter :math:`\tilde{D}_{v,\, i}` , and geometric standard deviation :math:`\sigma _{g,\, i}` , per dust source mode :math:`i`

 +-------------+-----------------------------+-----------------------------------+-----------------------------+
 | :math:`i`   | :math:`m_{i}`  (fraction)   | :math:`\tilde{D}_{v,\, i}`  (m)   | :math:`\sigma _{g,\, i}`    |
 +=============+=============================+===================================+=============================+
 | 1           | 0.036                       | 0.832 x 10\ :math:`{}^{-6}`       | 2.1                         |
 +-------------+-----------------------------+-----------------------------------+-----------------------------+
 | 2           | 0.957                       | 4.820 x 10\ :math:`{}^{-6}`       | 1.9                         |
 +-------------+-----------------------------+-----------------------------------+-----------------------------+
 | 3           | 0.007                       | 19.38 x 10\ :math:`{}^{-6}`       | 1.6                         |
 +-------------+-----------------------------+-----------------------------------+-----------------------------+

.. _Table Dust Minimum and maximum particle diameters:

.. table:: Minimum and maximum particle diameters in each dust transport bin :math:`j`

 +-------------+-------------------------------+-------------------------------+
 | :math:`j`   | :math:`D_{j,\, \min }`  (m)   | :math:`D_{j,\, \max }`  (m)   |
 +=============+===============================+===============================+
 | 1           | 0.1 x 10\ :math:`{}^{-6}`     | 1.0 x 10\ :math:`{}^{-6}`     |
 +-------------+-------------------------------+-------------------------------+
 | 2           | 1.0 x 10\ :math:`{}^{-6}`     | 2.5 x 10\ :math:`{}^{-6}`     |
 +-------------+-------------------------------+-------------------------------+
 | 3           | 2.5 x 10\ :math:`{}^{-6}`     | 5.0 x 10\ :math:`{}^{-6}`     |
 +-------------+-------------------------------+-------------------------------+
 | 4           | 5.0 x 10\ :math:`{}^{-6}`     | 10.0 x 10\ :math:`{}^{-6}`    |
 +-------------+-------------------------------+-------------------------------+
