.. _rst_Surface Albedos:

Surface Albedos
==================

.. _Canopy Radiative Transfer:

Canopy Radiative Transfer
-----------------------------

Radiative transfer within vegetative canopies is calculated from the two-stream approximation of :ref:`Dickinson (1983) <Dickinson1983>` and :ref:`Sellers (1985) <Sellers1985>` as described by :ref:`Bonan (1996) <Bonan1996>`

.. math::
   :label: 3.1

   -\bar{\mu }\frac{dI\, \uparrow }{d\left(L+S\right)} +\left[1-\left(1-\beta \right)\omega \right]I\, \uparrow -\omega \beta I\, \downarrow =\omega \bar{\mu }K\beta _{0} e^{-K\left(L+S\right)}

.. math::
   :label: 3.2

   \bar{\mu }\frac{dI\, \downarrow }{d\left(L+S\right)} +\left[1-\left(1-\beta \right)\omega \right]I\, \downarrow -\omega \beta I\, \uparrow =\omega \bar{\mu }K\left(1-\beta _{0} \right)e^{-K\left(L+S\right)}

where :math:`I\, \uparrow` and :math:`I\, \downarrow` are the upward and downward diffuse radiative fluxes per unit incident flux, :math:`K={G\left(\mu \right)\mathord{\left/ {\vphantom {G\left(\mu \right) \mu }} \right.} \mu }` is the optical depth of direct beam per unit leaf and stem area, :math:`\mu` is the cosine of the zenith angle of the incident beam, :math:`G\left(\mu \right)` is the relative projected area of leaf and stem elements in the direction :math:`\cos ^{-1} \mu`, :math:`\bar{\mu }` is the average inverse diffuse optical depth per unit leaf and stem area, :math:`\omega` is a scattering coefficient, :math:`\beta` and :math:`\beta _{0}` are upscatter parameters for diffuse and direct beam radiation, respectively, :math:`L` is the exposed leaf area index, and :math:`S` is the exposed stem area index (section :numref:`Phenology and vegetation burial by snow`). Given the direct beam albedo :math:`\alpha _{g,\, \Lambda }^{\mu }` and diffuse albedo :math:`\alpha _{g,\, \Lambda }` of the ground (section :numref:`Ground Albedos`), these equations are solved to calculate the fluxes, per unit incident flux, absorbed by the vegetation, reflected by the vegetation, and transmitted through the vegetation for direct and diffuse radiation and for visible (:math:`<` 0.7\ :math:`\mu {\rm m}`) and near-infrared (:math:`\geq` 0.7\ :math:`\mu {\rm m}`) wavebands. The absorbed radiation is partitioned to sunlit and shaded fractions of the canopy. The optical parameters :math:`G\left(\mu \right)`, :math:`\bar{\mu }`, :math:`\omega`, :math:`\beta`, and :math:`\beta _{0}` are calculated based on work in :ref:`Sellers (1985) <Sellers1985>` as follows.

The relative projected area of leaves and stems in the direction :math:`\cos ^{-1} \mu` is

.. math::
   :label: 3.3

   G\left(\mu \right)=\phi _{1} +\phi _{2} \mu

where :math:`\phi _{1} ={\rm 0.5}-0.633\chi _{L} -0.33\chi _{L}^{2}` and :math:`\phi _{2} =0.877\left(1-2\phi _{1} \right)` for :math:`-0.4\le \chi _{L} \le 0.6`. :math:`\chi _{L}` is the departure of leaf angles from a random distribution and equals +1 for horizontal leaves, 0 for random leaves, and –1 for vertical leaves.

The average inverse diffuse optical depth per unit leaf and stem area is

.. math::
   :label: 3.4

   \bar{\mu }=\int _{0}^{1}\frac{\mu '}{G\left(\mu '\right)}  d\mu '=\frac{1}{\phi _{2} } \left[1-\frac{\phi _{1} }{\phi _{2} } \ln \left(\frac{\phi _{1} +\phi _{2} }{\phi _{1} } \right)\right]

where :math:`\mu '` is the direction of the scattered flux.

The optical parameters :math:`\omega`, :math:`\beta`, and :math:`\beta _{0}`, which vary with wavelength (:math:`\Lambda` ), are weighted combinations of values for vegetation and snow, using the canopy snow-covered fraction :math:`f_{can,\, sno}` (Chapter :numref:`rst_Hydrology`). The optical parameters are

..
   The model determines that snow is on the canopy if :math:`T_{v} \le T_{f}`, where :math:`T_{v}` is the vegetation temperature (K) (Chapter :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`) and :math:`T_{f}` is the freezing temperature of water (K) (:numref:`Table Physical Constants`). In this case, the optical parameters are

.. math::
   :label: 3.5

   \omega _{\Lambda } =\omega _{\Lambda }^{veg} \left(1-f_{can,\, sno} \right)+\omega _{\Lambda }^{sno} f_{can,\, sno}

.. math::
   :label: 3.6

   \omega _{\Lambda } \beta _{\Lambda } =\omega _{\Lambda }^{veg} \beta _{\Lambda }^{veg} \left(1-f_{can,\, sno} \right)+\omega _{\Lambda }^{sno} \beta _{\Lambda }^{sno} f_{can,\, sno}

.. math::
   :label: 3.7

   \omega _{\Lambda } \beta _{0,\, \Lambda } =\omega _{\Lambda }^{veg} \beta _{0,\, \Lambda }^{veg} \left(1-f_{can,\, sno} \right)+\omega _{\Lambda }^{sno} \beta _{0,\, \Lambda }^{sno} f_{can,\, sno}

The snow and vegetation weights are applied to the products :math:`\omega _{\Lambda } \beta _{\Lambda }` and :math:`\omega _{\Lambda } \beta _{0,\, \Lambda }` because these products are used in the two-stream equations. If there is no snow on the canopy, this reduces to

.. math::
   :label: 3.8

   \omega _{\Lambda } =\omega _{\Lambda }^{veg}

.. math::
   :label: 3.9

   \omega _{\Lambda } \beta _{\Lambda } =\omega _{\Lambda }^{veg} \beta _{\Lambda }^{veg}

.. math::
   :label: 3.10

   \omega _{\Lambda } \beta _{0,\, \Lambda } =\omega _{\Lambda }^{veg} \beta _{0,\, \Lambda }^{veg} .

For vegetation, :math:`\omega _{\Lambda }^{veg} =\alpha _{\Lambda } +\tau _{\Lambda }`. :math:`\alpha _{\Lambda }` is a weighted combination of the leaf and stem reflectances (:math:`\alpha _{\Lambda }^{leaf},\alpha _{\Lambda }^{stem}` )

.. math::
   :label: 3.11

   \alpha _{\Lambda } =\alpha _{\Lambda }^{leaf} w_{leaf} +\alpha _{\Lambda }^{stem} w_{stem}

where :math:`w_{leaf} ={L\mathord{\left/ {\vphantom {L \left(L+S\right)}} \right.} \left(L+S\right)}` and :math:`w_{stem} ={S\mathord{\left/ {\vphantom {S \left(L+S\right)}} \right.} \left(L+S\right)}`. :math:`\tau _{\Lambda }` is a weighted combination of the leaf and stem transmittances (:math:`\tau _{\Lambda }^{leaf}, \tau _{\Lambda }^{stem}`)

.. math::
   :label: 3.12

   \tau _{\Lambda } =\tau _{\Lambda }^{leaf} w_{leaf} +\tau _{\Lambda }^{stem} w_{stem} .

The upscatter for diffuse radiation is

.. math::
   :label: 3.13

   \omega _{\Lambda }^{veg} \beta _{\Lambda }^{veg} =\frac{1}{2} \left[\alpha _{\Lambda } +\tau _{\Lambda } +\left(\alpha _{\Lambda } -\tau _{\Lambda } \right)\cos ^{2} \bar{\theta }\right]

where :math:`\bar{\theta }` is the mean leaf inclination angle relative to the horizontal plane (i.e., the angle between leaf normal and local vertical) (:ref:`Sellers (1985) <Sellers1985>`). Here, :math:`\cos \bar{\theta }` is approximated by

.. math::
   :label: 3.14

   \cos \bar{\theta }=\frac{1+\chi _{L} }{2}

Using this approximation, for vertical leaves (:math:`\chi _{L} =-1`, :math:`\bar{\theta }=90^{{\rm o}}` ), :math:`\omega _{\Lambda }^{veg} \beta _{\Lambda }^{veg} =0.5\left(\alpha _{\Lambda } +\tau _{\Lambda } \right)`, and for horizontal leaves (:math:`\chi _{L} =1`, :math:`\bar{\theta }=0^{{\rm o}}` ), :math:`\omega _{\Lambda }^{veg} \beta _{\Lambda }^{veg} =\alpha _{\Lambda }`, which agree with both :ref:`Dickinson (1983) <Dickinson1983>` and :ref:`Sellers (1985) <Sellers1985>`. For random (spherically distributed) leaves (:math:`\chi _{L} =0`, :math:`\bar{\theta }=60^{{\rm o}}` ), the approximation yields :math:`\omega _{\Lambda }^{veg} \beta _{\Lambda }^{veg} ={5\mathord{\left/ {\vphantom {5 8}} \right.} 8} \alpha _{\Lambda } +{3\mathord{\left/ {\vphantom {3 8}} \right.} 8} \tau _{\Lambda }` whereas the approximate solution of :ref:`Dickinson (1983) <Dickinson1983>` is :math:`\omega _{\Lambda }^{veg} \beta _{\Lambda }^{veg} ={2\mathord{\left/ {\vphantom {2 3}} \right.} 3} \alpha _{\Lambda } +{1\mathord{\left/ {\vphantom {1 3}} \right.} 3} \tau _{\Lambda }`. This discrepancy arises from the fact that a spherical leaf angle distribution has a true mean leaf inclination :math:`\bar{\theta }\approx 57` :ref:`(Campbell and Norman 1998) <CampbellNorman1998>` in equation :eq:`3.13`, while :math:`\bar{\theta }=60` in equation :eq:`3.14`. The upscatter for direct beam radiation is

.. math::
   :label: 3.15

   \omega _{\Lambda }^{veg} \beta _{0,\, \Lambda }^{veg} =\frac{1+\bar{\mu }K}{\bar{\mu }K} a_{s} \left(\mu \right)_{\Lambda }

where the single scattering albedo is

.. math::
   :label: 3.16

   \begin{array}{rcl} {a_{s} \left(\mu \right)_{\Lambda } } & {=} & {\frac{\omega _{\Lambda }^{veg} }{2} \int _{0}^{1}\frac{\mu 'G\left(\mu \right)}{\mu G\left(\mu '\right)+\mu 'G\left(\mu \right)}  d\mu '} \\ {} & {=} & {\frac{\omega _{\Lambda }^{veg} }{2} \frac{G\left(\mu \right)}{\max (\mu \phi _{2} +G\left(\mu \right),1e-6)} \left[1-\frac{\mu \phi _{1} }{\max (\mu \phi _{2} +G\left(\mu \right),1e-6)} \ln \left(\frac{\mu \phi _{1} +\max (\mu \phi _{2} +G\left(\mu \right),1e-6)}{\mu \phi _{1} } \right)\right].} \end{array}

Note here the restriction on :math:`\mu \phi _{2} +G\left(\mu \right)`. We have seen cases where small values can cause unrealistic single scattering albedo associated with the log calculation, thereby eventually causing a negative soil albedo.

The upward diffuse fluxes per unit incident direct beam and diffuse flux (i.e., the surface albedos) are

.. math::
   :label: 3.17

   I\, \uparrow _{\Lambda }^{\mu } =\frac{h_{1} }{\sigma } +h_{2} +h_{3}

.. math::
   :label: 3.18

   I\, \uparrow _{\Lambda } =h_{7} +h_{8} .

The downward diffuse fluxes per unit incident direct beam and diffuse radiation, respectively, are

.. math::
   :label: 3.19

   I\, \downarrow _{\Lambda }^{\mu } =\frac{h_{4} }{\sigma } e^{-K\left(L+S\right)} +h_{5} s_{1} +\frac{h_{6} }{s_{1} }

.. math::
   :label: 3.20

   I\, \downarrow _{\Lambda } =h_{9} s_{1} +\frac{h_{10} }{s_{1} } .

With reference to :numref:`Figure Radiation Schematic`, the direct beam flux transmitted through the canopy, per unit incident flux, is :math:`e^{-K\left(L+S\right)}`, and the direct beam and diffuse fluxes absorbed by the vegetation, per unit incident flux, are

.. math::
   :label: 3.21

   \vec{I}_{\Lambda }^{\mu } =1-I\, \uparrow _{\Lambda }^{\mu } -\left(1-\alpha _{g,\, \Lambda } \right)I\, \downarrow _{\Lambda }^{\mu } -\left(1-\alpha _{g,\, \Lambda }^{\mu } \right)e^{-K\left(L+S\right)}

.. math::
   :label: 3.22

   \vec{I}_{\Lambda } =1-I\, \uparrow _{\Lambda } -\left(1-\alpha _{g,\, \Lambda } \right)I\, \downarrow _{\Lambda } .

These fluxes are partitioned to the sunlit and shaded canopy using an analytical solution to the two-stream approximation for sunlit and shaded leaves :ref:`(Dai et al. 2004) <Daietal2004>`, as described by :ref:`Bonan et al. (2011) <Bonanetal2011>`. The absorption of direct beam radiation by sunlit leaves is

.. math::
   :label: 3.23

   \vec{I}_{sun,\Lambda }^{\mu } =\left(1-\omega _{\Lambda } \right)\left[1-s_{2} +\frac{1}{\bar{\mu }} \left(a_{1} +a_{2} \right)\right]

and for shaded leaves is

.. math::
   :label: 3.24

   \vec{I}_{sha,\Lambda }^{\mu } =\vec{I}_{\Lambda }^{\mu } -\vec{I}_{sun,\Lambda }^{\mu }

with

.. math::
   :label: 3.25

   a_{1} =\frac{h_{1} }{\sigma } \left[\frac{1-s_{2}^{2} }{2K} \right]+h_{2} \left[\frac{1-s_{2} s_{1} }{K+h} \right]+h_{3} \left[\frac{1-{s_{2} \mathord{\left/ {\vphantom {s_{2}  s_{1} }} \right.} s_{1} } }{K-h} \right]

.. math::
   :label: 3.26

   a_{2} =\frac{h_{4} }{\sigma } \left[\frac{1-s_{2}^{2} }{2K} \right]+h_{5} \left[\frac{1-s_{2} s_{1} }{K+h} \right]+h_{6} \left[\frac{1-{s_{2} \mathord{\left/ {\vphantom {s_{2}  s_{1} }} \right.} s_{1} } }{K-h} \right].

For diffuse radiation, the absorbed radiation for sunlit leaves is

.. math::
   :label: 3.27

   \vec{I}_{sun,\Lambda }^{} =\left[\frac{1-\omega _{\Lambda } }{\bar{\mu }} \right]\left(a_{1} +a_{2} \right)

and for shaded leaves is

.. math::
   :label: 3.28

   \vec{I}_{sha,\Lambda }^{} =\vec{I}_{\Lambda }^{} -\vec{I}_{sun,\Lambda }^{}

with

.. math::
   :label: 3.29

   a_{1} =h_{7} \left[\frac{1-s_{2} s_{1} }{K+h} \right]+h_{8} \left[\frac{1-{s_{2} \mathord{\left/ {\vphantom {s_{2}  s_{1} }} \right.} s_{1} } }{K-h} \right]

.. math::
   :label: 3.30

   a_{2} =h_{9} \left[\frac{1-s_{2} s_{1} }{K+h} \right]+h_{10} \left[\frac{1-{s_{2} \mathord{\left/ {\vphantom {s_{2}  s_{1} }} \right.} s_{1} } }{K-h} \right].

The parameters :math:`h_{1}` –:math:`h_{10}`, :math:`\sigma`, :math:`h`, :math:`s_{1}`, and :math:`s_{2}` are from :ref:`Sellers (1985) <Sellers1985>` [note the error in :math:`h_{4}` in :ref:`Sellers (1985) <Sellers1985>`]:

.. math::
   :label: 3.31

   b=1-\omega _{\Lambda } +\omega _{\Lambda } \beta _{\Lambda }

.. math::
   :label: 3.32

   c=\omega _{\Lambda } \beta _{\Lambda }

.. math::
   :label: 3.33

   d=\omega _{\Lambda } \bar{\mu }K\beta _{0,\, \Lambda }

.. math::
   :label: 3.34

   f=\omega _{\Lambda } \bar{\mu }K\left(1-\beta _{0,\, \Lambda } \right)

.. math::
   :label: 3.35

   h=\frac{\sqrt{b^{2} -c^{2} } }{\bar{\mu }}

.. math::
   :label: 3.36

   \sigma =\left(\bar{\mu }K\right)^{2} +c^{2} -b^{2}

.. math::
   :label: 3.37

   u_{1} =b-{c\mathord{\left/ {\vphantom {c \alpha _{g,\, \Lambda }^{\mu } }} \right.} \alpha _{g,\, \Lambda }^{\mu } } {\rm \; or\; }u_{1} =b-{c\mathord{\left/ {\vphantom {c \alpha _{g,\, \Lambda } }} \right.} \alpha _{g,\, \Lambda } }

.. math::
   :label: 3.38

   u_{2} =b-c\alpha _{g,\, \Lambda }^{\mu } {\rm \; or\; }u_{2} =b-c\alpha _{g,\, \Lambda }

.. math::
   :label: 3.39

   u_{3} =f+c\alpha _{g,\, \Lambda }^{\mu } {\rm \; or\; }u_{3} =f+c\alpha _{g,\, \Lambda }

.. math::
   :label: 3.40

   s_{1} =\exp \left\{-\min \left[h\left(L+S\right),40\right]\right\}

.. math::
   :label: 3.41

   s_{2} =\exp \left\{-\min \left[K\left(L+S\right),40\right]\right\}

.. math::
   :label: 3.42

   p_{1} =b+\bar{\mu }h

.. math::
   :label: 3.43

   p_{2} =b-\bar{\mu }h

.. math::
   :label: 3.44

   p_{3} =b+\bar{\mu }K

.. math::
   :label: 3.45

   p_{4} =b-\bar{\mu }K

.. math::
   :label: 3.46

   d_{1} =\frac{p_{1} \left(u_{1} -\bar{\mu }h\right)}{s_{1} } -p_{2} \left(u_{1} +\bar{\mu }h\right)s_{1}

.. math::
   :label: 3.47

   d_{2} =\frac{u_{2} +\bar{\mu }h}{s_{1} } -\left(u_{2} -\bar{\mu }h\right)s_{1}

.. math::
   :label: 3.48

   h_{1} =-dp_{4} -cf

.. math::
   :label: 3.49

   h_{2} =\frac{1}{d_{1} } \left[\left(d-\frac{h_{1} }{\sigma } p_{3} \right)\frac{\left(u_{1} -\bar{\mu }h\right)}{s_{1} } -p_{2} \left(d-c-\frac{h_{1} }{\sigma } \left(u_{1} +\bar{\mu }K\right)\right)s_{2} \right]

.. math::
   :label: 3.50

   h_{3} =\frac{-1}{d_{1} } \left[\left(d-\frac{h_{1} }{\sigma } p_{3} \right)\left(u_{1} +\bar{\mu }h\right)s_{1} -p_{1} \left(d-c-\frac{h_{1} }{\sigma } \left(u_{1} +\bar{\mu }K\right)\right)s_{2} \right]

.. math::
   :label: 3.51

   h_{4} =-fp_{3} -cd

.. math::
   :label: 3.52

   h_{5} =\frac{-1}{d_{2} } \left[\left(\frac{h_{4} \left(u_{2} +\bar{\mu }h\right)}{\sigma s_{1} } \right)+\left(u_{3} -\frac{h_{4} }{\sigma } \left(u_{2} -\bar{\mu }K\right)\right)s_{2} \right]

.. math::
   :label: 3.53

   h_{6} =\frac{1}{d_{2} } \left[\frac{h_{4} }{\sigma } \left(u_{2} -\bar{\mu }h\right)s_{1} +\left(u_{3} -\frac{h_{4} }{\sigma } \left(u_{2} -\bar{\mu }K\right)\right)s_{2} \right]

.. math::
   :label: 3.54

   h_{7} =\frac{c\left(u_{1} -\bar{\mu }h\right)}{d_{1} s_{1} }

.. math::
   :label: 3.55

   h_{8} =\frac{-c\left(u_{1} +\bar{\mu }h\right)s_{1} }{d_{1} }

.. math::
   :label: 3.56

   h_{9} =\frac{u_{2} +\bar{\mu }h}{d_{2} s_{1} }

.. math::
   :label: 3.57

   h_{10} =\frac{-s_{1} \left(u_{2} -\bar{\mu }h\right)}{d_{2} } .

Plant functional type optical properties (:numref:`Table Plant functional type optical properties`) for trees and shrubs are from :ref:`Dorman and Sellers (1989) <DormanSellers1989>`. Leaf and stem optical properties (VIS and NIR reflectance and transmittance) were derived for grasslands and crops from full optical range spectra of measured optical properties (:ref:`Asner et al. 1998 <Asneretal1998>`). Optical properties for intercepted snow (:numref:`Table Intercepted snow optical properties`) are from :ref:`Sellers et al. (1986) <Sellersetal1986>`.

.. _Table Plant functional type optical properties:

.. table:: Plant functional type optical properties

 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | Plant Functional Type            | :math:`\chi _{L}`    | :math:`\alpha _{vis}^{leaf}`    | :math:`\alpha _{nir}^{leaf}`    | :math:`\alpha _{vis}^{stem}`    | :math:`\alpha _{nir}^{stem}`    | :math:`\tau _{vis}^{leaf}`    | :math:`\tau _{nir}^{leaf}`    | :math:`\tau _{vis}^{stem}`    | :math:`\tau _{nir}^{stem}`    |
 +==================================+======================+=================================+=================================+=================================+=================================+===============================+===============================+===============================+===============================+
 | NET Temperate                    | 0.01                 | 0.07                            | 0.35                            | 0.16                            | 0.39                            | 0.05                          | 0.10                          | 0.001                         | 0.001                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | NET Boreal                       | 0.01                 | 0.07                            | 0.35                            | 0.16                            | 0.39                            | 0.05                          | 0.10                          | 0.001                         | 0.001                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | NDT Boreal                       | 0.01                 | 0.07                            | 0.35                            | 0.16                            | 0.39                            | 0.05                          | 0.10                          | 0.001                         | 0.001                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | BET Tropical                     | 0.10                 | 0.10                            | 0.45                            | 0.16                            | 0.39                            | 0.05                          | 0.25                          | 0.001                         | 0.001                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | BET temperate                    | 0.10                 | 0.10                            | 0.45                            | 0.16                            | 0.39                            | 0.05                          | 0.25                          | 0.001                         | 0.001                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | BDT tropical                     | 0.01                 | 0.10                            | 0.45                            | 0.16                            | 0.39                            | 0.05                          | 0.25                          | 0.001                         | 0.001                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | BDT temperate                    | 0.25                 | 0.10                            | 0.45                            | 0.16                            | 0.39                            | 0.05                          | 0.25                          | 0.001                         | 0.001                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | BDT boreal                       | 0.25                 | 0.10                            | 0.45                            | 0.16                            | 0.39                            | 0.05                          | 0.25                          | 0.001                         | 0.001                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | BES temperate                    | 0.01                 | 0.07                            | 0.35                            | 0.16                            | 0.39                            | 0.05                          | 0.10                          | 0.001                         | 0.001                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | BDS temperate                    | 0.25                 | 0.10                            | 0.45                            | 0.16                            | 0.39                            | 0.05                          | 0.25                          | 0.001                         | 0.001                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | BDS boreal                       | 0.25                 | 0.10                            | 0.45                            | 0.16                            | 0.39                            | 0.05                          | 0.25                          | 0.001                         | 0.001                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | C\ :sub:`3` arctic grass         | -0.30                | 0.11                            | 0.35                            | 0.31                            | 0.53                            | 0.05                          | 0.34                          | 0.120                         | 0.250                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | C\ :sub:`3` grass                | -0.30                | 0.11                            | 0.35                            | 0.31                            | 0.53                            | 0.05                          | 0.34                          | 0.120                         | 0.250                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | C\ :sub:`4` grass                | -0.30                | 0.11                            | 0.35                            | 0.31                            | 0.53                            | 0.05                          | 0.34                          | 0.120                         | 0.250                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | C\ :sub:`3` Crop                 | -0.30                | 0.11                            | 0.35                            | 0.31                            | 0.53                            | 0.05                          | 0.34                          | 0.120                         | 0.250                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | Temp Corn                        | -0.50                | 0.11                            | 0.35                            | 0.31                            | 0.53                            | 0.05                          | 0.34                          | 0.120                         | 0.250                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | Spring Wheat                     | -0.50                | 0.11                            | 0.35                            | 0.31                            | 0.53                            | 0.05                          | 0.34                          | 0.120                         | 0.250                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | Temp Soybean                     | -0.50                | 0.11                            | 0.35                            | 0.31                            | 0.53                            | 0.05                          | 0.34                          | 0.120                         | 0.250                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | Cotton                           | -0.50                | 0.11                            | 0.35                            | 0.31                            | 0.53                            | 0.05                          | 0.34                          | 0.120                         | 0.250                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | Rice                             | -0.50                | 0.11                            | 0.35                            | 0.31                            | 0.53                            | 0.05                          | 0.34                          | 0.120                         | 0.250                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | Sugarcane                        | -0.50                | 0.11                            | 0.35                            | 0.31                            | 0.53                            | 0.05                          | 0.34                          | 0.120                         | 0.250                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | Tropical Corn                    | -0.50                | 0.11                            | 0.35                            | 0.31                            | 0.53                            | 0.05                          | 0.34                          | 0.120                         | 0.250                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | Tropical Soybean                 | -0.50                | 0.11                            | 0.35                            | 0.31                            | 0.53                            | 0.05                          | 0.34                          | 0.120                         | 0.250                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | Miscanthus                       | -0.50                | 0.11                            | 0.35                            | 0.31                            | 0.53                            | 0.05                          | 0.34                          | 0.120                         | 0.250                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+
 | Switchgrass                      | -0.50                | 0.11                            | 0.35                            | 0.31                            | 0.53                            | 0.05                          | 0.34                          | 0.120                         | 0.250                         |
 +----------------------------------+----------------------+---------------------------------+---------------------------------+---------------------------------+---------------------------------+-------------------------------+-------------------------------+-------------------------------+-------------------------------+

.. _Table Intercepted snow optical properties:

.. table:: Intercepted snow optical properties

 +-----------------------------+-------+-------+
 | Parameter                   | vis   | nir   |
 +=============================+=======+=======+
 | :math:`\omega ^{sno}`       | 0.8   | 0.4   |
 +-----------------------------+-------+-------+
 | :math:`\beta ^{sno}`        | 0.5   | 0.5   |
 +-----------------------------+-------+-------+
 | :math:`\beta _{0}^{sno}`    | 0.5   | 0.5   |
 +-----------------------------+-------+-------+

.. _Ground Albedos:

Ground Albedos
------------------

The overall direct beam :math:`\alpha _{g,\, \Lambda }^{\mu }` and diffuse :math:`\alpha _{g,\, \Lambda }` ground albedos are weighted combinations of "soil" and snow albedos

.. math::
   :label: 3.58

   \alpha _{g,\, \Lambda }^{\mu } =\alpha _{soi,\, \Lambda }^{\mu } \left(1-f_{sno} \right)+\alpha _{sno,\, \Lambda }^{\mu } f_{sno}

.. math::
   :label: 3.59

   \alpha _{g,\, \Lambda } =\alpha _{soi,\, \Lambda } \left(1-f_{sno} \right)+\alpha _{sno,\, \Lambda } f_{sno}

where :math:`f_{sno}` is the fraction of the ground covered with snow (section :numref:`Snow Covered Area Fraction`).

:math:`\alpha _{soi,\, \Lambda }^{\mu }` and :math:`\alpha _{soi,\, \Lambda }` vary with glacier, lake, and soil surfaces. Glacier albedos are from :ref:`Paterson (1994) <Paterson1994>`

.. math:: \alpha _{soi,\, vis}^{\mu } =\alpha _{soi,\, vis} =0.6

.. math:: \alpha _{soi,\, nir}^{\mu } =\alpha _{soi,\, nir} =0.4.

Unfrozen lake albedos depend on the cosine of the solar zenith angle :math:`\mu`

.. math::
   :label: 3.60

   \alpha _{soi,\, \Lambda }^{\mu } =\alpha _{soi,\, \Lambda } =0.05\left(\mu +0.15\right)^{-1} .

Frozen lake albedos are from NCAR LSM (:ref:`Bonan 1996 <Bonan1996>`)

.. math:: \alpha _{soi,\, vis}^{\mu } =\alpha _{soi,\, vis} =0.60

.. math:: \alpha _{soi,\, nir}^{\mu } =\alpha _{soi,\, nir} =0.40.

As in NCAR LSM (:ref:`Bonan 1996 <Bonan1996>`), soil albedos vary with color class

.. math::
   :label: 3.61

   \alpha _{soi,\, \Lambda }^{\mu } =\alpha _{soi,\, \Lambda } =\left(\alpha _{sat,\, \Lambda } +\Delta \right)\le \alpha _{dry,\, \Lambda }

where :math:`\Delta` depends on the volumetric water content of the first soil layer :math:`\theta _{1}` (section :numref:`Soil Water`) as :math:`\Delta =0.11-0.40\theta _{1} >0`, and :math:`\alpha _{sat,\, \Lambda }` and :math:`\alpha _{dry,\, \Lambda }` are albedos for saturated and dry soil color classes (:numref:`Table Dry and saturated soil albedos`).

CLM soil colors are prescribed so that they best reproduce observed MODIS local solar noon surface albedo values at the CLM grid cell following the methods of :ref:`Lawrence and Chase (2007) <LawrenceChase2007>`. The soil colors are fitted over the range of 20 soil classes shown in :numref:`Table Dry and saturated soil albedos` and compared to the MODIS monthly local solar noon all-sky surface albedo as described in :ref:`Strahler et al. (1999) <Strahleretal1999>` and :ref:`Schaaf et al. (2002) <Schaafetal2002>`. The CLM two-stream radiation model was used to calculate the model equivalent surface albedo using climatological monthly soil moisture along with the vegetation parameters of PFT fraction, LAI, and SAI. The soil color that produced the closest all-sky albedo in the two-stream radiation model was selected as the best fit for the month. The fitted monthly soil colors were averaged over all snow-free months to specify a representative soil color for the grid cell. In cases where there was no snow-free surface albedo for the year, the soil color derived from snow-affected albedo was used to give a representative soil color that included the effects of the minimum permanent snow cover.

.. _Table Dry and saturated soil albedos:

.. table:: Dry and saturated soil albedos

 +---------------+--------+--------+--------+--------+---------------+--------+--------+--------+--------+
 |               |       Dry       |    Saturated    |               |       Dry       |    Saturated    |
 +---------------+--------+--------+--------+--------+---------------+--------+--------+--------+--------+
 | Color Class   | vis    | nir    | vis    | nir    | Color Class   | vis    | nir    | vis    | nir    |
 +---------------+--------+--------+--------+--------+---------------+--------+--------+--------+--------+
 | 1             | 0.36   | 0.61   | 0.25   | 0.50   | 11            | 0.24   | 0.37   | 0.13   | 0.26   |
 +---------------+--------+--------+--------+--------+---------------+--------+--------+--------+--------+
 | 2             | 0.34   | 0.57   | 0.23   | 0.46   | 12            | 0.23   | 0.35   | 0.12   | 0.24   |
 +---------------+--------+--------+--------+--------+---------------+--------+--------+--------+--------+
 | 3             | 0.32   | 0.53   | 0.21   | 0.42   | 13            | 0.22   | 0.33   | 0.11   | 0.22   |
 +---------------+--------+--------+--------+--------+---------------+--------+--------+--------+--------+
 | 4             | 0.31   | 0.51   | 0.20   | 0.40   | 14            | 0.20   | 0.31   | 0.10   | 0.20   |
 +---------------+--------+--------+--------+--------+---------------+--------+--------+--------+--------+
 | 5             | 0.30   | 0.49   | 0.19   | 0.38   | 15            | 0.18   | 0.29   | 0.09   | 0.18   |
 +---------------+--------+--------+--------+--------+---------------+--------+--------+--------+--------+
 | 6             | 0.29   | 0.48   | 0.18   | 0.36   | 16            | 0.16   | 0.27   | 0.08   | 0.16   |
 +---------------+--------+--------+--------+--------+---------------+--------+--------+--------+--------+
 | 7             | 0.28   | 0.45   | 0.17   | 0.34   | 17            | 0.14   | 0.25   | 0.07   | 0.14   |
 +---------------+--------+--------+--------+--------+---------------+--------+--------+--------+--------+
 | 8             | 0.27   | 0.43   | 0.16   | 0.32   | 18            | 0.12   | 0.23   | 0.06   | 0.12   |
 +---------------+--------+--------+--------+--------+---------------+--------+--------+--------+--------+
 | 9             | 0.26   | 0.41   | 0.15   | 0.30   | 19            | 0.10   | 0.21   | 0.05   | 0.10   |
 +---------------+--------+--------+--------+--------+---------------+--------+--------+--------+--------+
 | 10            | 0.25   | 0.39   | 0.14   | 0.28   | 20            | 0.08   | 0.16   | 0.04   | 0.08   |
 +---------------+--------+--------+--------+--------+---------------+--------+--------+--------+--------+

.. _Snow Albedo:

Snow Albedo
^^^^^^^^^^^^^^^^^

Snow albedo and solar absorption within each snow layer are simulated with the Snow, Ice, and Aerosol Radiative Model (SNICAR), which incorporates a two-stream radiative transfer solution from :ref:`Toon et al. (1989) <Toonetal1989>`. Albedo and the vertical absorption profile depend on solar zenith angle, albedo of the substrate underlying snow, mass concentrations of atmospheric-deposited aerosols (black carbon, mineral dust, and organic carbon), and ice effective grain size (:math:`r_{e}`), which is simulated with a snow aging routine described in section :numref:`Snow Aging`. Representation of impurity mass concentrations within the snowpack is described in section :numref:`Black and organic carbon and mineral dust within snow`. Implementation of SNICAR in CLM is also described somewhat by :ref:`Flanner and Zender (2005) <FlannerZender2005>` and :ref:`Flanner et al. (2007) <Flanneretal2007>`.

The two-stream solution requires the following bulk optical properties for each snow layer and spectral band: extinction optical depth (:math:`\tau`), single-scatter albedo (:math:`\omega`), and scattering asymmetry parameter (*g*). The snow layers used for radiative calculations are identical to snow layers applied elsewhere in CLM, except for the case when snow mass is greater than zero but no snow layers exist. When this occurs, a single radiative layer is specified to have the column snow mass and an effective grain size of freshly-fallen snow (section :numref:`Snow Aging`). The bulk optical properties are weighted functions of each constituent *k*, computed for each snow layer and spectral band as

.. math::
   :label: 3.62

   \tau =\sum _{1}^{k}\tau _{k}

.. math::
   :label: 3.63

   \omega =\frac{\sum _{1}^{k}\omega _{k} \tau _{k}  }{\sum _{1}^{k}\tau _{k}  }

.. math::
   :label: 3.64

   g=\frac{\sum _{1}^{k}g_{k} \omega _{k} \tau _{k}  }{\sum _{1}^{k}\omega _{k} \tau _{k}  }

For each constituent (ice, two black carbon species, two organic carbon species, and four dust species), :math:`\omega`, *g*, and the mass extinction cross-section :math:`\psi` (m\ :sup:`2` kg\ :sub:`-1`) are computed offline with Mie Theory, e.g., applying the computational technique from :ref:`Bohren and Huffman (1983) <BohrenHuffman1983>`. The extinction optical depth for each constituent depends on its mass extinction cross-section and layer mass, :math:`w _{k}` (kg\ m\ :sup:`-1`) as

.. math::
   :label: 3.65

   \tau _{k} =\psi _{k} w_{k}

The two-stream solution (:ref:`Toon et al. (1989) <Toonetal1989>`) applies a tri-diagonal matrix solution to produce upward and downward radiative fluxes at each layer interface, from which net radiation, layer absorption, and surface albedo are easily derived. Solar fluxes are computed in five spectral bands, listed in :numref:`Table Spectral bands and weights used for snow radiative transfer`. Because snow albedo varies strongly across the solar spectrum, it was determined that four bands were needed to accurately represent the near-infrared (NIR) characteristics of snow, whereas only one band was needed for the visible spectrum. Boundaries of the NIR bands were selected to capture broad radiative features and maximize accuracy and computational efficiency. We partition NIR (0.7-5.0 :math:`\mu` m) surface downwelling flux from CLM according to the weights listed in :numref:`Table Spectral bands and weights used for snow radiative transfer`, which are unique for diffuse and direct incident flux. These fixed weights were determined with offline hyperspectral radiative transfer calculations for an atmosphere typical of mid-latitude winter (:ref:`Flanner et al. (2007) <Flanneretal2007>`). The tri-diagonal solution includes intermediate terms that allow for easy interchange of two-stream techniques. We apply the Eddington solution for the visible band (following :ref:`Wiscombe and Warren 1980 <WiscombeWarren1980>`) and the hemispheric mean solution ((:ref:`Toon et al. (1989) <Toonetal1989>`) for NIR bands. These choices were made because the Eddington scheme works well for highly scattering media, but can produce negative albedo for absorptive NIR bands with diffuse incident flux. Delta scalings are applied to :math:`\tau`, :math:`\omega`, and :math:`g` (:ref:`Wiscombe and Warren 1980 <WiscombeWarren1980>`) in all spectral bands, producing effective values (denoted with :math:`*`) that are applied in the two-stream solution

.. math::
   :label: 3.66

   \tau ^{*} =\left(1-\omega g^{2} \right)\tau

.. math::
   :label: 3.67

   \omega ^{*} =\frac{\left(1-g^{2} \right)\omega }{1-g^{2} \omega }

.. math::
   :label: 3.68

   g^{*} =\frac{g}{1+g}

.. _Table Spectral bands and weights used for snow radiative transfer:

.. table:: Spectral bands and weights used for snow radiative transfer

 +---------------------------------------------------------+----------------------+------------------+
 | Spectral band                                           | Direct-beam weight   | Diffuse weight   |
 +=========================================================+======================+==================+
 | Band 1: 0.3-0.7\ :math:`\mu`\ m (visible)               | (1.0)                | (1.0)            |
 +---------------------------------------------------------+----------------------+------------------+
 | Band 2: 0.7-1.0\ :math:`\mu`\ m (near-IR)               | 0.494                | 0.586            |
 +---------------------------------------------------------+----------------------+------------------+
 | Band 3: 1.0-1.2\ :math:`\mu`\ m (near-IR)               | 0.181                | 0.202            |
 +---------------------------------------------------------+----------------------+------------------+
 | Band 4: 1.2-1.5\ :math:`\mu`\ m (near-IR)               | 0.121                | 0.109            |
 +---------------------------------------------------------+----------------------+------------------+
 | Band 5: 1.5-5.0\ :math:`\mu`\ m (near-IR)               | 0.204                | 0.103            |
 +---------------------------------------------------------+----------------------+------------------+

Under direct-beam conditions, singularities in the radiative approximation are occasionally approached in spectral bands 4 and 5 that produce unrealistic conditions (negative energy absorption in a layer, negative albedo, or total absorbed flux greater than incident flux). When any of these three conditions occur, the Eddington approximation is attempted instead, and if both approximations fail, the cosine of the solar zenith angle is adjusted by 0.02 (conserving incident flux) and a warning message is produced. This situation occurs in only about 1 in 10 :sup:`6` computations of snow albedo. After looping over the five spectral bands, absorption fluxes and albedo are averaged back into the bulk NIR band used by the rest of CLM.

Soil albedo (or underlying substrate albedo), which is defined for visible and NIR bands, is a required boundary condition for the snow radiative transfer calculation. Currently, the bulk NIR soil albedo is applied to all four NIR snow bands. With ground albedo as a lower boundary condition, SNICAR simulates solar absorption in all snow layers as well as the underlying soil or ground. With a thin snowpack, penetrating solar radiation to the underlying soil can be quite large and heat cannot be released from the soil to the atmosphere in this situation. Thus, if the snowpack has total snow depth less than 0.1 m (:math:`z_{sno} < 0.1`) and there are no explicit snow layers, the solar radiation is absorbed by the top soil layer. If there is a single snow layer, the solar radiation is absorbed in that layer. If there is more than a single snow layer, 75% of the solar radiation is absorbed in the top snow layer, and 25% is absorbed in the next lowest snow layer. This prevents unrealistic soil warming within a single timestep.

The radiative transfer calculation is performed twice for each column containing a mass of snow greater than :math:`1 \times 10^{-30}` kg\ m\ :sup:`-2` (excluding lake and urban columns); once each for direct-beam and diffuse incident flux. Absorption in each layer :math:`i` of pure snow is initially recorded as absorbed flux per unit incident flux on the ground (:math:`S_{sno,\, i}` ), as albedos must be calculated for the next timestep with unknown incident flux. The snow absorption fluxes that are used for column temperature calculations are

.. math::
   :label: 3.69

   S_{g,\, i} =S_{sno,\, i} \left(1-\alpha _{sno} \right)

This weighting is performed for direct-beam and diffuse, visible and NIR fluxes. After the ground-incident fluxes (transmitted through the vegetation canopy) have been calculated for the current time step (sections :numref:`Canopy Radiative Transfer` and :numref:`Solar Fluxes`), the layer absorption factors (:math:`S_{g,\, i}`) are multiplied by the ground-incident fluxes to produce solar absorption (W m\ :sup:`-2`) in each snow layer and the underlying ground.

.. _Snowpack Optical Properties:

Snowpack Optical Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Ice optical properties for the five spectral bands are derived offline and stored in a namelist-defined lookup table for online retrieval (see CLM5.0 User's Guide). Mie properties are first computed at fine spectral resolution (470 bands), and are then weighted into the five bands applied by CLM according to incident solar flux, :math:`I^{\downarrow } (\lambda )`. For example, the broadband mass-extinction cross section (:math:`\bar{\psi }`) over wavelength interval :math:`\lambda _{1}` to :math:`\lambda _{2}` is

.. math::
   :label: 3.70

   \bar{\psi }=\frac{\int _{\lambda _{1} }^{\lambda _{2} }\psi \left(\lambda \right) I^{\downarrow } \left(\lambda \right){\rm d}\lambda }{\int _{\lambda _{1} }^{\lambda _{2} }I^{\downarrow } \left(\lambda \right){\rm d}\lambda  }

Broadband single-scatter albedo (:math:`\bar{\omega }`) is additionally weighted by the diffuse albedo for a semi-infinite snowpack (:math:`\alpha _{sno}`)

.. math::
   :label: 3.71

   \bar{\omega }=\frac{\int _{\lambda _{1} }^{\lambda _{2} }\omega (\lambda )I^{\downarrow } ( \lambda )\alpha _{sno} (\lambda ){\rm d}\lambda }{\int _{\lambda _{1} }^{\lambda _{2} }I^{\downarrow } ( \lambda )\alpha _{sno} (\lambda ){\rm d}\lambda }

Inclusion of this additional albedo weight was found to improve accuracy of the five-band albedo solutions (relative to 470-band solutions) because of the strong dependence of optically-thick snowpack albedo on ice grain single-scatter albedo (:ref:`Flanner et al. (2007) <Flanneretal2007>`). The lookup tables contain optical properties for lognormal distributions of ice particles over the range of effective radii: 30\ :math:`\mu`\ m :math:`< r _{e} < \text{1500} \mu \text{m}`, at 1 :math:`\mu` m resolution. Single-scatter albedos for the end-members of this size range are listed in :numref:`Table Single-scatter albedo values used for snowpack impurities and ice`.

Optical properties for black carbon are described in :ref:`Flanner et al. (2007) <Flanneretal2007>`. Single-scatter albedo, mass extinction cross-section, and asymmetry parameter values for all snowpack species, in the five spectral bands used, are listed in :numref:`Table Single-scatter albedo values used for snowpack impurities and ice`, :numref:`Table Mass extinction values`, and :numref:`Table Asymmetry scattering parameters used for snowpack impurities and ice`. These properties were also derived with Mie Theory, using various published sources of indices of refraction and assumptions about particle size distribution. Weighting into the five CLM spectral bands was determined only with incident solar flux, as in equation :eq:`3.69`.

.. _Table Single-scatter albedo values used for snowpack impurities and ice:

.. table:: Single-scatter albedo values used for snowpack impurities and ice

 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Species                                                        | Band 1   | Band 2   | Band 3   | Band 4   | Band 5   |
 +================================================================+==========+==========+==========+==========+==========+
 | Hydrophilic black carbon                                       | 0.516    | 0.434    | 0.346    | 0.276    | 0.139    |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Hydrophobic black carbon                                       | 0.288    | 0.187    | 0.123    | 0.089    | 0.040    |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Hydrophilic organic carbon                                     | 0.997    | 0.994    | 0.990    | 0.987    | 0.951    |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Hydrophobic organic carbon                                     | 0.963    | 0.921    | 0.860    | 0.814    | 0.744    |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Dust 1                                                         | 0.979    | 0.994    | 0.993    | 0.993    | 0.953    |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Dust 2                                                         | 0.944    | 0.984    | 0.989    | 0.992    | 0.983    |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Dust 3                                                         | 0.904    | 0.965    | 0.969    | 0.973    | 0.978    |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Dust 4                                                         | 0.850    | 0.940    | 0.948    | 0.953    | 0.955    |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Ice (:math:`r _{e}` = 30 :math:`\mu` m)                        | 0.9999   | 0.9999   | 0.9992   | 0.9938   | 0.9413   |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Ice (:math:`r _{e}` = 1500 :math:`\mu` m)                      | 0.9998   | 0.9960   | 0.9680   | 0.8730   | 0.5500   |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+

.. _Table Mass extinction values:

.. table:: Mass extinction values (m\ :sup:`2` kg\ :sup:`-1`) used for snowpack impurities and ice

 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Species                                                        | Band 1   | Band 2   | Band 3   | Band 4   | Band 5   |
 +================================================================+==========+==========+==========+==========+==========+
 | Hydrophilic black carbon                                       | 25369    | 12520    | 7739     | 5744     | 3527     |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Hydrophobic black carbon                                       | 11398    | 5923     | 4040     | 3262     | 2224     |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Hydrophilic organic carbon                                     | 37774    | 22112    | 14719    | 10940    | 5441     |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Hydrophobic organic carbon                                     | 3289     | 1486     | 872      | 606      | 248      |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Dust 1                                                         | 2687     | 2420     | 1628     | 1138     | 466      |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Dust 2                                                         | 841      | 987      | 1184     | 1267     | 993      |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Dust 3                                                         | 388      | 419      | 400      | 397      | 503      |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Dust 4                                                         | 197      | 203      | 208      | 205      | 229      |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Ice (:math:`r _{e}` = 30 :math:`\mu` m)                        | 55.7     | 56.1     | 56.3     | 56.6     | 57.3     |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Ice (:math:`r _{e}` = 1500 :math:`\mu` m)                      | 1.09     | 1.09     | 1.09     | 1.09     | 1.1      |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+

.. _Table Asymmetry scattering parameters used for snowpack impurities and ice:

.. table:: Asymmetry scattering parameters used for snowpack impurities and ice.

 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Species                                                        | Band 1   | Band 2   | Band 3   | Band 4   | Band 5   |
 +================================================================+==========+==========+==========+==========+==========+
 | Hydrophilic black carbon                                       | 0.52     | 0.34     | 0.24     | 0.19     | 0.10     |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Hydrophobic black carbon                                       | 0.35     | 0.21     | 0.15     | 0.11     | 0.06     |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Hydrophilic organic carbon                                     | 0.77     | 0.75     | 0.72     | 0.70     | 0.64     |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Hydrophobic organic carbon                                     | 0.62     | 0.57     | 0.54     | 0.51     | 0.44     |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Dust 1                                                         | 0.69     | 0.72     | 0.67     | 0.61     | 0.44     |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Dust 2                                                         | 0.70     | 0.65     | 0.70     | 0.72     | 0.70     |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Dust 3                                                         | 0.79     | 0.75     | 0.68     | 0.63     | 0.67     |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Dust 4                                                         | 0.83     | 0.79     | 0.77     | 0.76     | 0.73     |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Ice (:math:`r _{e}` = 30\ :math:`\mu`\ m)                      | 0.88     | 0.88     | 0.88     | 0.88     | 0.90     |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+
 | Ice (:math:`r _{e}` = 1500\ :math:`\mu`\ m)                    | 0.89     | 0.90     | 0.90     | 0.92     | 0.97     |
 +----------------------------------------------------------------+----------+----------+----------+----------+----------+

.. _Snow Aging:

Snow Aging
^^^^^^^^^^^^^^^^

Snow aging is represented as evolution of the ice effective grain size (:math:`r_{e}`). Previous studies have shown that use of spheres which conserve the surface area-to-volume ratio (or specific surface area) of ice media composed of more complex shapes produces relatively small errors in simulated hemispheric fluxes (e.g., :ref:`Grenfell and Warren 1999 <GrenfellWarren1999>`). Effective radius is the surface area-weighted mean radius of an ensemble of spherical particles and is directly related to specific surface area (*SSA*) as :math:`r_{e} ={3\mathord{\left/ {\vphantom {3 \left(\rho _{ice} SSA\right)}} \right.} \left(\rho _{ice} SSA\right)}`, where :math:`\rho_{ice}` is the density of ice. Hence, :math:`r_{e}` is a simple and practical metric for relating the snowpack microphysical state to dry snow radiative characteristics.

Wet snow processes can also drive rapid changes in albedo. The presence of liquid water induces rapid coarsening of the surrounding ice grains (e.g., :ref:`Brun 1989 <Brun1989>`), and liquid water tends to refreeze into large ice clumps that darken the bulk snowpack. The presence of small liquid drops, by itself, does not significantly darken snowpack, as ice and water have very similar indices of refraction throughout the solar spectrum. Pooled or ponded water, however, can significantly darken snowpack by greatly reducing the number of refraction events per unit mass. This influence is not currently accounted for.

The net change in effective grain size occurring each time step is represented in each snow layer as a summation of changes caused by dry snow metamorphism (:math:`dr_{e,dry}`), liquid water-induced metamorphism (:math:`dr_{e,wet}`), refreezing of liquid water, and addition of freshly-fallen snow. The mass of each snow layer is partitioned into fractions of snow carrying over from the previous time step (:math:`f_{old}`), freshly-fallen snow (:math:`f_{new}`), and refrozen liquid water (:math:`f_{rfz}`), such that snow :math:`r_{e}` is updated each time step *t* as

.. math::
   :label: 3.72

   r_{e} \left(t\right)=\left[r_{e} \left(t-1\right)+dr_{e,\, dry} +dr_{e,\, wet} \right]f_{old} +r_{e,\, 0} f_{new} +r_{e,\, rfz} f_{rfrz}

Here, the effective radius of freshly-fallen snow (:math:`r_{e,0}`) is based on a simple linear temperature-relationship. Below -30 degrees Celsius, a minimum value is enforced of 54.5 :math:`\mu` m (corresponding to a specific surface area of 60 m\ :sup:`2` kg\ :sup:`-1`). Above 0 degrees Celsius, a maximum value is enforced of 204.5 :math:`\mu` m. Between -30 and 0 a linear ramp is used.

The effective radius of refrozen liquid water (:math:`r_{e,rfz}`) is set to 1000\ :math:`\mu` m.

Dry snow aging is based on a microphysical model described by :ref:`Flanner and Zender (2006) <FlannerZender2006>`. This model simulates diffusive vapor flux amongst collections of ice crystals with various size and inter-particle spacing. Specific surface area and effective radius are prognosed for any combination of snow temperature, temperature gradient, density, and initial size distribution. The combination of warm snow, large temperature gradient, and low density produces the most rapid snow aging, whereas aging proceeds slowly in cold snow, regardless of temperature gradient and density. Because this model is currently too computationally expensive for inclusion in climate models, we fit parametric curves to model output over a wide range of snow conditions and apply these parameters in CLM. The functional form of the parametric equation is

.. math::
   :label: 3.73

   \frac{dr_{e,\, dry} }{dt} =\left(\frac{dr_{e} }{dt} \right)_{0} \left(\frac{\eta }{\left(r_{e} -r_{e,\, 0} \right)+\eta } \right)^{{1\mathord{\left/ {\vphantom {1 \kappa }} \right.} \kappa } }

The parameters :math:`{(\frac{dr_{e}}{dt}})_{0}`, :math:`\eta`, and :math:`\kappa` are retrieved interactively from a lookup table with dimensions corresponding to snow temperature, temperature gradient, and density. The domain covered by this lookup table includes temperature ranging from 223 to 273 K, temperature gradient ranging from 0 to 300 K m\ :sup:`-1`, and density ranging from 50 to 400 kg m\ :sup:`-3`. Temperature gradient is calculated at the midpoint of each snow layer *n*, using mid-layer temperatures (:math:`T_{n}`) and snow layer thicknesses (:math:`dz_{n}`), as

.. math::
   :label: 3.74

   \left(\frac{dT}{dz} \right)_{n} =\frac{1}{dz_{n} } abs\left[\frac{T_{n-1} dz_{n} +T_{n} dz_{n-1} }{dz_{n} +dz_{n-1} } +\frac{T_{n+1} dz_{n} +T_{n} dz_{n+1} }{dz_{n} +dz_{n+1} } \right]

For the bottom snow layer (:math:`n=0`), :math:`T_{n+1}` is taken as the temperature of the top soil layer, and for the top snow layer it is assumed that :math:`T_{n-1}` = :math:`T_{n}`.

The contribution of liquid water to enhanced metamorphism is based on parametric equations published by :ref:`Brun (1989) <Brun1989>`, who measured grain growth rates under different liquid water contents. This relationship, expressed in terms of :math:`r_{e} (\mu \text{m})` and subtracting an offset due to dry aging, depends on the mass liquid water fraction :math:`f_{liq}` as

.. math::
   :label: 3.75

   \frac{dr_{e} }{dt} =\frac{10^{18} C_{1} f_{liq} ^{3} }{4\pi r_{e} ^{2} }

The constant *C*\ :sub:`1` is 4.22\ :math:`\times`\ 10\ :sup:`-13`, and: :math:`f_{liq} =w_{liq} /(w_{liq} +w_{ice} )`\ (Chapter :numref:`rst_Snow Hydrology`).

In cases where snow mass is greater than zero, but a snow layer has not yet been defined, :math:`r_{e}` is set to :math:`r_{e,0}`. When snow layers are combined or divided, :math:`r_{e}` is calculated as a mass-weighted mean of the two layers, following computations of other state variables (section :numref:`Snow Layer Combination and Subdivision`). Finally, the allowable range of :math:`r_{e}`, corresponding to the range over which Mie optical properties have been defined, is 30-1500\ :math:`\mu` m.

.. _Solar Zenith Angle:

Solar Zenith Angle
----------------------

The CLM uses the same formulation for solar zenith angle as the Community Atmosphere Model. The cosine of the solar zenith angle :math:`\mu` is

.. math::
   :label: 3.76

   \mu =\sin \phi \sin \delta -\cos \phi \cos \delta \cos h

where :math:`h` is the solar hour angle (radians) (24 hour periodicity), :math:`\delta` is the solar declination angle (radians), and :math:`\phi` is latitude (radians) (positive in Northern Hemisphere). The solar hour angle :math:`h` (radians) is

.. math::
   :label: 3.77

   h=2\pi d+\theta

where :math:`d` is calendar day (:math:`d=0.0` at 0Z on January 1), and :math:`\theta` is longitude (radians) (positive east of the Greenwich meridian).

The solar declination angle :math:`\delta` is calculated as in :ref:`Berger (1978a,b) <Berger1978a>` and is valid for one million years past or hence, relative to 1950 A.D. The orbital parameters may be specified directly or the orbital parameters are calculated for the desired year. The required orbital parameters to be input by the user are the obliquity of the Earth :math:`\varepsilon` (degrees, :math:`-90^{\circ } <\varepsilon <90^{\circ }` ), Earth's eccentricity :math:`e` (:math:`0.0<e<0.1`), and the longitude of the perihelion relative to the moving vernal equinox :math:`\tilde{\omega }` (:math:`0^{\circ } <\tilde{\omega }<360^{\circ }` ) (unadjusted for the apparent orbit of the Sun around the Earth (:ref:`Berger et al. 1993 <Bergeretal1993>`)). The solar declination :math:`\delta` (radians) is

.. math::
   :label: 3.78

   \delta =\sin ^{-1} \left[\sin \left(\varepsilon \right)\sin \left(\lambda \right)\right]

where :math:`\varepsilon` is Earth's obliquity and :math:`\lambda` is the true longitude of the Earth.

The obliquity of the Earth :math:`\varepsilon` (degrees) is

.. math::
   :label: 3.79

   \varepsilon =\varepsilon *+\sum _{i=1}^{i=47}A_{i}  \cos \left(f_{i} t+\delta _{i} \right)

where :math:`\varepsilon *` is a constant of integration (:numref:`Table Orbital parameters`), :math:`A_{i}`, :math:`f_{i}`, and :math:`\delta _{i}` are amplitude, mean rate, and phase terms in the cosine series expansion (:ref:`Berger (1978a,b) <Berger1978a>`, and :math:`t=t_{0} -1950` where :math:`t_{0}` is the year. The series expansion terms are not shown here but can be found in the source code file shr\_orb\_mod.F90.

The true longitude of the Earth :math:`\lambda` (radians) is counted counterclockwise from the vernal equinox (:math:`\lambda =0` at the vernal equinox)

.. math::
   :label: 3.80

   \lambda =\lambda _{m} +\left(2e-\frac{1}{4} e^{3} \right)\sin \left(\lambda _{m} -\tilde{\omega }\right)+\frac{5}{4} e^{2} \sin 2\left(\lambda _{m} -\tilde{\omega }\right)+\frac{13}{12} e^{3} \sin 3\left(\lambda _{m} -\tilde{\omega }\right)

where :math:`\lambda _{m}` is the mean longitude of the Earth at the vernal equinox, :math:`e` is Earth's eccentricity, and :math:`\tilde{\omega }` is the longitude of the perihelion relative to the moving vernal equinox. The mean longitude :math:`\lambda _{m}` is

.. math::
   :label: 3.81

   \lambda _{m} =\lambda _{m0} +\frac{2\pi \left(d-d_{ve} \right)}{365}

where :math:`d_{ve} =80.5` is the calendar day at vernal equinox (March 21 at noon), and

.. math::
   :label: 3.82

   \lambda _{m0} =2\left[\left(\frac{1}{2} e+\frac{1}{8} e^{3} \right)\left(1+\beta \right)\sin \tilde{\omega }-\frac{1}{4} e^{2} \left(\frac{1}{2} +\beta \right)\sin 2\tilde{\omega }+\frac{1}{8} e^{3} \left(\frac{1}{3} +\beta \right)\sin 3\tilde{\omega }\right]

where :math:`\beta =\sqrt{1-e^{2} }`. Earth's eccentricity :math:`e` is

.. math::
   :label: 3.83

   e=\sqrt{\left(e^{\cos } \right)^{2} +\left(e^{\sin } \right)^{2} }

where

.. math::
   :label: 3.84

   \begin{array}{l} {e^{\cos } =\sum _{j=1}^{19}M_{j} \cos \left(g_{j} t+B_{j} \right) ,} \\ {e^{\sin } =\sum _{j=1}^{19}M_{j} \sin \left(g_{j} t+B_{j} \right) } \end{array}

are the cosine and sine series expansions for :math:`e`, and :math:`M_{j}`, :math:`g_{j}`, and :math:`B_{j}` are amplitude, mean rate, and phase terms in the series expansions (:ref:`Berger (1978a,b) <Berger1978a>`). The longitude of the perihelion relative to the moving vernal equinox :math:`\tilde{\omega }` (degrees) is

.. math::
   :label: 3.85

   \tilde{\omega }=\Pi \frac{180}{\pi } +\psi

where :math:`\Pi` is the longitude of the perihelion measured from the reference vernal equinox (i.e., the vernal equinox at 1950 A.D.) and describes the absolute motion of the perihelion relative to the fixed stars, and :math:`\psi` is the annual general precession in longitude and describes the absolute motion of the vernal equinox along Earth's orbit relative to the fixed stars. The general precession :math:`\psi` (degrees) is

.. math::
   :label: 3.86

   \psi =\frac{\tilde{\psi }t}{3600} +\zeta +\sum _{i=1}^{78}F_{i}  \sin \left(f_{i} ^{{'} } t+\delta _{i} ^{{'} } \right)

where :math:`\tilde{\psi }` (arcseconds) and :math:`\zeta` (degrees) are constants (:numref:`Table Orbital parameters`), and :math:`F_{i}`, :math:`f_{i} ^{{'} }`, and :math:`\delta _{i} ^{{'} }` are amplitude, mean rate, and phase terms in the sine series expansion (:ref:`Berger (1978a,b) <Berger1978a>`)). The longitude of the perihelion :math:`\Pi` (radians) depends on the sine and cosine series expansions for the eccentricity :math:`e`\ as follows:

.. math::
   :label: 3.87

   \Pi =\left\{\begin{array}{lr}
   0 & \qquad {\rm for\; -1}\times {\rm 10}^{{\rm -8}} \le e^{\cos } \le 1\times 10^{-8} {\rm \; and\; }e^{\sin } = 0 \\
   1.5\pi & \qquad {\rm for\; -1}\times {\rm 10}^{{\rm -8}} \le e^{\cos } \le 1\times 10^{-8} {\rm \; and\; }e^{\sin } < 0 \\
   0.5\pi & \qquad {\rm for\; -1}\times {\rm 10}^{{\rm -8}} \le e^{\cos } \le 1\times 10^{-8} {\rm \; and\; }e^{\sin } > 0 \\
   \tan ^{-1} \left[\frac{e^{\sin } }{e^{\cos } } \right]+\pi & \qquad {\rm for\; }e^{\cos } <{\rm -1}\times {\rm 10}^{{\rm -8}}  \\
   \tan ^{-1} \left[\frac{e^{\sin } }{e^{\cos } } \right]+2\pi & \qquad {\rm for\; }e^{\cos } >{\rm 1}\times {\rm 10}^{{\rm -8}} {\rm \; and\; }e^{\sin } <0 \\
   \tan ^{-1} \left[\frac{e^{\sin } }{e^{\cos } } \right] & \qquad {\rm for\; }e^{\cos } >{\rm 1}\times {\rm 10}^{{\rm -8}} {\rm \; and\; }e^{\sin } \ge 0
   \end{array}\right\}.

The numerical solution for the longitude of the perihelion :math:`\tilde{\omega }` is constrained to be between 0 and 360 degrees (measured from the autumn equinox). A constant 180 degrees is then added to :math:`\tilde{\omega }` because the Sun is considered as revolving around the Earth (geocentric coordinate system) (:ref:`Berger et al. 1993 <Bergeretal1993>`)).

.. _Table Orbital parameters:

.. table:: Orbital parameters

 +--------------------------------------+-------------+
 | Parameter                            |             |
 +======================================+=============+
 | :math:`\varepsilon *`                | 23.320556   |
 +--------------------------------------+-------------+
 | :math:`\tilde{\psi }` (arcseconds)   | 50.439273   |
 +--------------------------------------+-------------+
 | :math:`\zeta`  (degrees)             | 3.392506    |
 +--------------------------------------+-------------+
