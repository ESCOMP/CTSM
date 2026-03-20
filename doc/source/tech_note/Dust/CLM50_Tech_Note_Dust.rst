.. _rst_Dust Model:

Dust Model
==============

Atmospheric dust is mobilized from the land by wind in the CLM. The most important factors determining soil erodibility and dust emission include the wind friction speed, the vegetation cover, and the soil moisture The CLM dust mobilization scheme (:ref:`Mahowald et al. 2006<Mahowaldetal2006>` accounts for these factors based on the DEAD (Dust Entrainment and Deposition model of :ref:`Zender et al. (2003)<Zenderetal2003>`. Please refer to the :ref:`Zender et al. (2003)<Zenderetal2003>` article for additional information regarding the equations presented in this section.

The total vertical mass flux of dust, :math:`F_{j}` (kg m\ :sup:`-2` s\ :sup:`-1`), from the ground into transport bin :math:`j` is given by

.. math::
   :label: 29.1

   F_{j} =TSf_{m} \alpha Q_{s} \sum _{i=1}^{I}M_{i,j}

where :math:`T` is a global factor that compensates for the DEAD model's sensitivity to horizontal and temporal resolution and equals 5 x 10\ :sup:`-4` in the CLM instead of 7 x 10\ :sup:`-4` in :ref:`Zender et al. (2003)<Zenderetal2003>`. :math:`S` is the source erodibility factor set to 1 in the CLM and serves as a place holder at this time.

The grid cell fraction of exposed bare soil suitable for dust mobilization :math:`f_{m}` is given by

.. math::
   :label: 29.2

   f_{m} =\left(1-f_{lake} \right)\left(1-f_{sno} \right)\left(1-f_{v} \right)\frac{w_{liq,1} }{w_{liq,1} +w_{ice,1} }

where :math:`f_{lake}` and :math:`f_{sno}` are the CLM grid cell fractions of lake (section :numref:`Surface Data`) and snow cover (section :numref:`Snow Covered Area Fraction`), all ranging from zero to one. Not mentioned by :ref:`Zender et al. (2003)<Zenderetal2003>`, :math:`w_{liq,\, 1}` and :math:`{}_{w_{ice,\, 1} }` are the CLM top soil layer liquid water and ice contents (mm) entered as a ratio expressing the decreasing ability of dust to mobilize from increasingly frozen soil. The grid cell fraction of vegetation cover,\ :math:`{}_{f_{v} }`, is defined as

.. math::
   :label: 29.3

   0\le f_{v} =\frac{L+S}{\left(L+S\right)_{t} } \le 1{\rm \; \; \; \; where\; }\left(L+S\right)_{t} =0.3{\rm \; m}^{2} {\rm m}^{-2}

where equation :eq:`29.3` applies only for dust mobilization and is not related to the plant functional type fractions prescribed from the CLM input data or simulated by the CLM dynamic vegetation model (Chapter 22). :math:`L` and :math:`S` are the CLM leaf and stem area index values (m :sup:`2` m\ :sup:`-2`) averaged at the land unit level so as to include all the pfts and the bare ground present in a vegetated land unit. :math:`L` and :math:`S` may be prescribed from the CLM input data (section :numref:`Phenology and vegetation burial by snow`) or simulated by the CLM biogeochemistry model (Chapter :numref:`rst_Vegetation Phenology and Turnover`).

The sandblasting mass efficiency :math:`\alpha` (m :sup:`-1`) is calculated as

.. math::
   :label: 29.4

   \alpha =100e^{\left(13.4M_{clay} -6.0\right)\ln 10} {\rm \; \; }\left\{\begin{array}{l} {M_{clay} =\% clay\times 0.01{\rm \; \; \; 0}\le \% clay\le 20} \\ {M_{clay} =20\times 0.01{\rm \; \; \; \; \; \; \; \; 20<\% }clay\le 100} \end{array}\right.

where :math:`M_{clay}` is the mass fraction of clay particles in the soil and %clay is determined from the surface dataset (section :numref:`Surface Data`). :math:`M_{clay} =0` corresponds to sand and :math:`M_{clay} =0.2` to sandy loam.

:math:`Q_{s}` is the total horizontally saltating mass flux (kg m\ :sup:`-1` s\ :sup:`-1`) of "large" particles (:numref:`Table Dust Mass fraction`), also referred to as the vertically integrated streamwise mass flux

.. math::
   :label: 29.5

   Q_{s} = \left\{
   \begin{array}{lr}
   \frac{c_{s} \rho _{atm} u_{*s}^{3} }{g} \left(1-\frac{u_{*t} }{u_{*s} } \right)\left(1+\frac{u_{*t} }{u_{*s} } \right)^{2} {\rm \; } & \qquad {\rm for\; }u_{*t} <u_{*s}  \\
   0{\rm \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; } & \qquad {\rm for\; }u_{*t} \ge u_{*s}
   \end{array}\right.

where the saltation constant :math:`c_{s}` equals 2.61 and :math:`\rho _{atm}` is the atmospheric density (kg m\ :sup:`-3`) (:numref:`Table Atmospheric input to land model`), :math:`g` the acceleration of gravity (m s\ :sup:`-2`) (:numref:`Table Physical constants`). The threshold wind friction speed for saltation :math:`u_{*t}` (m s\ :sup:`-1`) is

.. math::
   :label: 29.6

   u_{*t} =f_{z} \left[Re_{*t}^{f} \rho _{osp} gD_{osp} \left(1+\frac{6\times 10^{-7} }{\rho _{osp} gD_{osp}^{2.5} } \right)\right]^{\frac{1}{2} } \rho _{atm} ^{-\frac{1}{2} } f_{w}

where :math:`f_{z}` is a factor dependent on surface roughness but set to 1 as a place holder for now, :math:`\rho _{osp}` and :math:`D_{osp}` are the density (2650 kg m\ :sup:`-3`) and diameter (75 x 10\ :math:`{}^{-6}` m) of optimal saltation particles, and :math:`f_{w}` is a factor dependent on soil moisture:

.. math::
   :label: 29.7

   f_{w} =\left\{\begin{array}{l} {1{\rm \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; for\; }w\le w_{t} } \\ {\sqrt{1+1.21\left[100\left(w-w_{t} \right)\right]^{0.68} } {\rm \; \; for\; }w>w_{t} } \end{array}\right.

where

.. math::
   :label: 29.8

   w_{t} =a\left(0.17M_{clay} +0.14M_{clay}^{2} \right){\rm \; \; \; \; \; \; 0}\le M_{clay} =\% clay\times 0.01\le 1

and

.. math::
   :label: 29.9

   w=\frac{\theta _{1} \rho _{liq} }{\rho _{d,1} }

where :math:`a=M_{clay}^{-1}` for tuning purposes, :math:`\theta _{1}` is the volumetric soil moisture in the top soil layer (m :math:`{}^{3 }`\ m\ :sup:`-3`) (section :numref:`Soil Water`), :math:`\rho _{liq}` is the density of liquid water (kg m\ :sup:`-3`) (:numref:`Table Physical constants`), and :math:`\rho _{d,\, 1}` is the bulk density of soil in the top soil layer (kg m\ :sup:`-3`) defined as in section :numref:`Soil and Snow Thermal Properties` rather than as in :ref:`Zender et al. (2003)<Zenderetal2003>`. :math:`Re_{*t}^{f}` from equation :eq:`29.6` is the threshold friction Reynolds factor

.. math::
   :label: 29.10

   Re_{*t}^{f} =\left\{\begin{array}{l} {\frac{0.1291^{2} }{-1+1.928Re_{*t} } {\rm \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; for\; 0.03}\le Re_{*t} \le 10} \\ {0.12^{2} \left(1-0.0858e^{-0.0617(Re_{*t} -10)} \right)^{2} {\rm \; for\; }Re_{*t} >10} \end{array}\right.

and :math:`Re_{*t}` is the threshold friction Reynolds number approximation for optimally sized particles

.. math::
   :label: 29.11

   Re_{*t} =0.38+1331\left(100D_{osp} \right)^{1.56}

In :eq:`29.5`, :math:`u_{*s}` is defined as the wind friction speed (m s\ :sup:`-1`) accounting for the Owen effect (:ref:`Owen 1964<Owen1964>`)

.. math::
   :label: 29.12

   u_{*s} = \left\{
   \begin{array}{lr}
   u_{*} & \quad {\rm \; for \;} U_{10} <U_{10,t}  \\
   u_{*} +0.003\left(U_{10} -U_{10,t} \right)^{2} & \quad {\rm \; for\; }U_{10} \ge U_{10,t}
   \end{array}\right.

where :math:`u_{*}` is the CLM wind friction speed (m s\ :sup:`-1`), also known as friction velocity (section :numref:`Monin-Obukhov Similarity Theory`), :math:`U_{10}` \ is the 10-m wind speed (m s\ :sup:`-1`) calculated as the wind speed at the top of the canopy in section 4.3 of :ref:`Bonan (1996)<Bonan1996>` but here for 10 m above the ground, and :math:`U_{10,\, t}` is the threshold wind speed at 10 m (m s\ :sup:`-1`)

.. math::
   :label: 29.13

   U_{10,t} =u_{*t} \frac{U_{10} }{u_{*} }

In equation :eq:`29.1` we sum :math:`M_{i,\, j}` over :math:`I=3` source modes :math:`i` where :math:`M_{i,\, j}` is the mass fraction of each source mode :math:`i` carried in each of *:math:`J=4`* transport bins :math:`j`

.. math::
   :label: 29.14

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
