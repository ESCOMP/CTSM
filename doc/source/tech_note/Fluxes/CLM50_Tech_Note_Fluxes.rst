.. _rst_Momentum, Sensible Heat, and Latent Heat Fluxes:

Momentum, Sensible Heat, and Latent Heat Fluxes
==================================================

The zonal :math:`\tau _{x}` and meridional :math:`\tau _{y}` momentum fluxes (kg m\ :sup:`-1` s\ :sup:`-2`), sensible heat flux :math:`H` (W m\ :sup:`-2`), and water vapor flux :math:`E` (kg m\ :sup:`-2` s\ :sup:`-1`) between the atmosphere at reference height :math:`z_{atm,\, x}` (m) [where :math:`x` is height for wind (momentum) (:math:`m`), temperature (sensible heat) (:math:`h`), and humidity (water vapor) (:math:`w`); with zonal and meridional winds :math:`u_{atm}` and :math:`v_{atm}` (m s\ :sup:`-1`), potential temperature :math:`\theta _{atm}` (K), and specific humidity :math:`q_{atm}` (kg kg\ :sup:`-1`)] and the surface [with :math:`u_{s}`, :math:`v_{s}`, :math:`\theta _{s}`, and :math:`q_{s}` ] are

.. math::
   :label: 5.1

   \tau _{x} =-\rho _{atm} \frac{\left(u_{atm} -u_{s} \right)}{r_{am} }

.. math::
   :label: 5.2

   \tau _{y} =-\rho _{atm} \frac{\left(v_{atm} -v_{s} \right)}{r_{am} }

.. math::
   :label: 5.3

   H=-\rho _{atm} C_{p} \frac{\left(\theta _{atm} -\theta _{s} \right)}{r_{ah} }

.. math::
   :label: 5.4

   E=-\rho _{atm} \frac{\left(q_{atm} -q_{s} \right)}{r_{aw} } .

These fluxes are derived in the next section from Monin-Obukhov similarity theory developed for the surface layer (i.e., the nearly constant flux layer above the surface sublayer). In this derivation, :math:`u_{s}` and :math:`v_{s}` are defined to equal zero at height :math:`z_{0m} +d` (the apparent sink for momentum) so that :math:`r_{am}` is the aerodynamic resistance (s m\ :sup:`-1`) for momentum between the atmosphere at height :math:`z_{atm,\, m}` and the surface at height :math:`z_{0m} +d`. Thus, the momentum fluxes become

.. math::
   :label: 5.5

   \tau _{x} =-\rho _{atm} \frac{u_{atm} }{r_{am} }

.. math::
   :label: 5.6

   \tau _{y} =-\rho _{atm} \frac{v_{atm} }{r_{am} } .

Likewise, :math:`\theta _{s}` and :math:`q_{s}` are defined at heights :math:`z_{0h} +d` and :math:`z_{0w} +d` (the apparent sinks for heat and water vapor, respectively :math:`r_{aw}` are the aerodynamic resistances (s m\ :sup:`-1`) to sensible heat and water vapor transfer between the atmosphere at heights :math:`z_{atm,\, h}` and :math:`z_{atm,\, w}` and the surface at heights :math:`z_{0h} +d` and :math:`z_{0w} +d`, respectively. The specific heat capacity of air :math:`C_{p}` (J kg\ :sup:`-1` K\ :sup:`-1`) is a constant (:numref:`Table Physical constants`). The atmospheric potential temperature used here is

.. math::
   :label: 5.7

   \theta _{atm} =T_{atm} +\Gamma _{d} z_{atm,\, h}

where :math:`T_{atm}` is the air temperature (K) at height :math:`z_{atm,\, h}` and :math:`\Gamma _{d} =0.0098` K m\ :sup:`-1` is the negative of the dry adiabatic lapse rate [this expression is first-order equivalent to :math:`\theta _{atm} =T_{atm} \left({P_{srf} \mathord{\left/ {\vphantom {P_{srf} P_{atm} }} \right.} P_{atm} } \right)^{{R_{da} \mathord{\left/ {\vphantom {R_{da} C_{p} }} \right.} C_{p} } }` (:ref:`Stull 1988 <Stull1988>`), where :math:`P_{srf}` is the surface pressure (Pa), :math:`P_{atm}` is the atmospheric pressure (Pa), and :math:`R_{da}` is the gas constant for dry air (J kg\ :sup:`-1` K\ :sup:`-1`) (:numref:`Table Physical constants`)]. By definition, :math:`\theta _{s} =T_{s}`. The density of moist air (kg m\ :sup:`-3`) is

.. math::
   :label: 5.8

   \rho _{atm} =\frac{P_{atm} -0.378e_{atm} }{R_{da} T_{atm} }

where the atmospheric vapor pressure :math:`e_{atm}` (Pa) is derived from the atmospheric specific humidity :math:`q_{atm}`

.. math::
   :label: 5.9

   e_{atm} =\frac{q_{atm} P_{atm} }{0.622+0.378q_{atm} } .

.. _Monin-Obukhov Similarity Theory:

Monin-Obukhov Similarity Theory
-----------------------------------

The surface vertical kinematic fluxes of momentum :math:`\overline{u'w'}` and :math:`\overline{v'w'}` (m\ :sup:`2` s\ :sub:`-2`), sensible heat :math:`\overline{\theta 'w'}` (K m s :sup:`-1`), and latent heat :math:`\overline{q'w'}` (kg kg\ :sup:`-1` m s\ :sup:`-1`), where :math:`u'`, :math:`v'`, :math:`w'`, :math:`\theta '`, and :math:`q'` are zonal horizontal wind, meridional horizontal wind, vertical velocity, potential temperature, and specific humidity turbulent fluctuations about the mean, are defined from Monin-Obukhov similarity applied to the surface layer. This theory states that when scaled appropriately, the dimensionless mean horizontal wind speed, mean potential temperature, and mean specific humidity profile gradients depend on unique functions of :math:`\zeta =\frac{z-d}{L}` (:ref:`Zeng et al. 1998<Zengetal1998>`) as

.. math::
   :label: 5.10

   \frac{k\left(z-d\right)}{u_{*} } \frac{\partial \left|{\it u}\right|}{\partial z} =\phi _{m} \left(\zeta \right)

.. math::
   :label: 5.11

   \frac{k\left(z-d\right)}{\theta _{*} } \frac{\partial \theta }{\partial z} =\phi _{h} \left(\zeta \right)

.. math::
   :label: 5.12

   \frac{k\left(z-d\right)}{q_{*} } \frac{\partial q}{\partial z} =\phi _{w} \left(\zeta \right)

where :math:`z` is height in the surface layer (m), :math:`d` is the displacement height (m), :math:`L` is the Monin-Obukhov length scale (m) that accounts for buoyancy effects resulting from vertical density gradients (i.e., the atmospheric stability), k is the von Karman constant (:numref:`Table Physical constants`), and :math:`\left|{\it u}\right|` is the atmospheric wind speed (m s\ :sup:`-1`). :math:`\phi _{m}`, :math:`\phi _{h}`, and :math:`\phi _{w}` are universal (over any surface) similarity functions of :math:`\zeta` that relate the constant fluxes of momentum, sensible heat, and latent heat to the mean profile gradients of :math:`\left|{\it u}\right|`, :math:`\theta`, and :math:`q` in the surface layer. In neutral conditions, :math:`\phi _{m} =\phi _{h} =\phi _{w} =1`. The velocity (i.e., friction velocity) :math:`u_{*}` (m s\ :sup:`-1`), temperature :math:`\theta _{*}` (K), and moisture :math:`q_{*}` (kg kg\ :sup:`-1`) scales are

.. math::
   :label: 5.13

   u_{*}^{2} =\sqrt{\left(\overline{u'w'}\right)^{2} +\left(\overline{v'w'}\right)^{2} } =\frac{\left|{\it \tau }\right|}{\rho _{atm} }

.. math::
   :label: 5.14

   \theta _{*} u_{*} =-\overline{\theta 'w'}=-\frac{H}{\rho _{atm} C_{p} }

.. math::
   :label: 5.15

   q_{*} u_{*} =-\overline{q'w'}=-\frac{E}{\rho _{atm} }

where :math:`\left|{\it \tau }\right|` is the shearing stress (kg m\ :sup:`-1` s\ :sup:`-2`), with zonal and meridional components :math:`\overline{u'w'}=-\frac{\tau _{x} }{\rho _{atm} }` and :math:`\overline{v'w'}=-\frac{\tau _{y} }{\rho _{atm} }`, respectively, :math:`H` is the sensible heat flux (W m\ :sup:`-2`) and :math:`E` is the water vapor flux (kg m\ :sup:`-2` s\ :sup:`-1`).

The length scale :math:`L` is the Monin-Obukhov length defined as

.. math::
   :label: 5.16

   L=-\frac{u_{*}^{3} }{k\left(\frac{g}{\overline{\theta _{v,\, atm} }} \right)\theta '_{v} w'} =\frac{u_{*}^{2} \overline{\theta _{v,\, atm} }}{kg\theta _{v*} }

where :math:`g` is the acceleration of gravity (m s\ :sup:`-2`) (:numref:`Table Physical constants`), and :math:`\overline{\theta _{v,\, atm} }=\overline{\theta _{atm} }\left(1+0.61q_{atm} \right)` is the reference virtual potential temperature. :math:`L>0` indicates stable conditions. :math:`L<0` indicates unstable conditions. :math:`L=\infty` for neutral conditions. The temperature scale :math:`\theta _{v*}` is defined as

.. math::
   :label: 5.17

   \theta _{v*} u_{*} =\left[\theta _{*} \left(1+0.61q_{atm} \right)+0.61\overline{\theta _{atm} }q_{*} \right]u_{*}

where :math:`\overline{\theta _{atm} }` is the atmospheric potential temperature.

Following :ref:`Panofsky and Dutton (1984)<PanofskyDutton1984>`, the differential equations for :math:`\phi _{m} \left(\zeta \right)`, :math:`\phi _{h} \left(\zeta \right)`, and :math:`\phi _{w} \left(\zeta \right)` can be integrated formally without commitment to their exact forms. Integration between two arbitrary heights in the surface layer :math:`z_{2}` and :math:`z_{1}` (:math:`z_{2} >z_{1}` ) with horizontal winds :math:`\left|{\it u}\right|_{1}` and :math:`\left|{\it u}\right|_{2}`, potential temperatures :math:`\theta _{1}` and :math:`\theta _{2}`, and specific humidities :math:`q_{1}` and :math:`q_{2}` results in

.. math::
   :label: 5.18

   \left|{\it u}\right|_{2} -\left|{\it u}\right|_{1} =\frac{u_{*} }{k} \left[\ln \left(\frac{z_{2} -d}{z_{1} -d} \right)-\psi _{m} \left(\frac{z_{2} -d}{L} \right)+\psi _{m} \left(\frac{z_{1} -d}{L} \right)\right]

.. math::
   :label: 5.19

   \theta _{2} -\theta _{1} =\frac{\theta _{*} }{k} \left[\ln \left(\frac{z_{2} -d}{z_{1} -d} \right)-\psi _{h} \left(\frac{z_{2} -d}{L} \right)+\psi _{h} \left(\frac{z_{1} -d}{L} \right)\right]

.. math::
   :label: 5.20

   q_{2} -q_{1} =\frac{q_{*} }{k} \left[\ln \left(\frac{z_{2} -d}{z_{1} -d} \right)-\psi _{w} \left(\frac{z_{2} -d}{L} \right)+\psi _{w} \left(\frac{z_{1} -d}{L} \right)\right].

The functions :math:`\psi _{m} \left(\zeta \right)`, :math:`\psi _{h} \left(\zeta \right)`, and :math:`\psi _{w} \left(\zeta \right)` are defined as

.. math::
   :label: 5.21

   \psi _{m} \left(\zeta \right)=\int _{{z_{0m} \mathord{\left/ {\vphantom {z_{0m}  L}} \right.} L} }^{\zeta }\frac{\left[1-\phi _{m} \left(x\right)\right]}{x} \, dx

.. math::
   :label: 5.22

   \psi _{h} \left(\zeta \right)=\int _{{z_{0h} \mathord{\left/ {\vphantom {z_{0h}  L}} \right.} L} }^{\zeta }\frac{\left[1-\phi _{h} \left(x\right)\right]}{x} \, dx

.. math::
   :label: 5.23

   \psi _{w} \left(\zeta \right)=\int _{{z_{0w} \mathord{\left/ {\vphantom {z_{0w}  L}} \right.} L} }^{\zeta }\frac{\left[1-\phi _{w} \left(x\right)\right]}{x} \, dx

where :math:`z_{0m}`, :math:`z_{0h}`, and :math:`z_{0w}` are the roughness lengths (m) for momentum, sensible heat, and water vapor, respectively.

Defining the surface values

.. math:: \left|{\it u}\right|_{1} =0{\rm \; at\; }z_{1} =z_{0m} +d,

.. math:: \theta _{1} =\theta _{s} {\rm \; at\; }z_{1} =z_{0h} +d,{\rm \; and}

.. math:: q_{1} =q_{s} {\rm \; at\; }z_{1} =z_{0w} +d,

and the atmospheric values at :math:`z_{2} =z_{atm,\, x}`

.. math::
   :label: 5.24

   \left|{\it u}\right|_{2} =V_{a} {\rm =\; }\sqrt{u_{atm}^{2} +v_{atm}^{2} +U_{c}^{2} } \ge 1,

.. math:: \theta _{2} =\theta _{atm} {\rm ,\; and}

.. math:: q_{2} =q_{atm} {\rm ,\; }

the integral forms of the flux-gradient relations are

.. math::
   :label: 5.25

   V_{a} =\frac{u_{*} }{k} \left[\ln \left(\frac{z_{atm,\, m} -d}{z_{0m} } \right)-\psi _{m} \left(\frac{z_{atm,\, m} -d}{L} \right)+\psi _{m} \left(\frac{z_{0m} }{L} \right)\right]

.. math::
   :label: 5.26

   \theta _{atm} -\theta _{s} =\frac{\theta _{*} }{k} \left[\ln \left(\frac{z_{atm,\, h} -d}{z_{0h} } \right)-\psi _{h} \left(\frac{z_{atm,\, h} -d}{L} \right)+\psi _{h} \left(\frac{z_{0h} }{L} \right)\right]

.. math::
   :label: 5.27

   q_{atm} -q_{s} =\frac{q_{*} }{k} \left[\ln \left(\frac{z_{atm,\, w} -d}{z_{0w} } \right)-\psi _{w} \left(\frac{z_{atm,\, w} -d}{L} \right)+\psi _{w} \left(\frac{z_{0w} }{L} \right)\right].

The constraint :math:`V_{a} \ge 1` is required simply for numerical reasons to prevent :math:`H` and :math:`E` from becoming small with small wind speeds. The convective velocity :math:`U_{c}` accounts for the contribution of large eddies in the convective boundary layer to surface fluxes as follows

.. math::
   :label: 5.28

   U_{c} = \left\{
   \begin{array}{ll}
   0 & \qquad \zeta \ge {\rm 0} \quad {\rm (stable)} \\
   \beta w_{*} & \qquad \zeta < 0 \quad {\rm (unstable)}
   \end{array} \right\}

where :math:`w_{*}` is the convective velocity scale

.. math::
   :label: 5.29

   w_{*} =\left(\frac{-gu_{*} \theta _{v*} z_{i} }{\overline{\theta _{v,\, atm} }} \right)^{{1\mathord{\left/ {\vphantom {1 3}} \right.} 3} } ,

:math:`z_{i} =1000` is the convective boundary layer height (m), and :math:`\beta =1`.

The momentum flux gradient relations are (:ref:`Zeng et al. 1998 <Zengetal1998>`)

.. math::
   :label: 5.30

   \begin{array}{llr}
   \phi _{m} \left(\zeta \right)=0.7k^{{2\mathord{\left/ {\vphantom {2 3}} \right.} 3} } \left(-\zeta \right)^{{1\mathord{\left/ {\vphantom {1 3}} \right.} 3} } & \qquad {\rm for\; }\zeta <-1.574 & \ {\rm \; (very\; unstable)} \\
   \phi _{m} \left(\zeta \right)=\left(1-16\zeta \right)^{-{1\mathord{\left/ {\vphantom {1 4}} \right.} 4} } & \qquad {\rm for\; -1.574}\le \zeta <0 & \ {\rm \; (unstable)} \\
   \phi _{m} \left(\zeta \right)=1+5\zeta & \qquad {\rm for\; }0\le \zeta \le 1& \ {\rm \; (stable)} \\
   \phi _{m} \left(\zeta \right)=5+\zeta & \qquad {\rm for\; }\zeta  >1 & \ {\rm\; (very\; stable).}
   \end{array}

The sensible and latent heat flux gradient relations are (:ref:`Zeng et al. 1998 <Zengetal1998>`)

.. math::
   :label: 5.31

   \begin{array}{llr}
   \phi _{h} \left(\zeta \right)=\phi _{w} \left(\zeta \right)=0.9k^{{4\mathord{\left/ {\vphantom {4 3}} \right.} 3} } \left(-\zeta \right)^{{-1\mathord{\left/ {\vphantom {-1 3}} \right.} 3} } & \qquad {\rm for\; }\zeta <-0.465 & \ {\rm \; (very\; unstable)} \\
   \phi _{h} \left(\zeta \right)=\phi _{w} \left(\zeta \right)=\left(1-16\zeta \right)^{-{1\mathord{\left/ {\vphantom {1 2}} \right.} 2} } & \qquad {\rm for\; -0.465}\le \zeta <0 & \ {\rm \; (unstable)} \\
   \phi _{h} \left(\zeta \right)=\phi _{w} \left(\zeta \right)=1+5\zeta & \qquad {\rm for\; }0\le \zeta \le 1 & \ {\rm \; (stable)} \\
   \phi _{h} \left(\zeta \right)=\phi _{w} \left(\zeta \right)=5+\zeta & \qquad {\rm for\; }\zeta  >1 & \ {\rm \; (very\; stable).}
   \end{array}

To ensure continuous functions of :math:`\phi _{m} \left(\zeta \right)`, :math:`\phi _{h} \left(\zeta \right)`, and :math:`\phi _{w} \left(\zeta \right)`, the simplest approach (i.e., without considering any transition regimes) is to match the relations for very unstable and unstable conditions at :math:`\zeta _{m} =-1.574` for :math:`\phi _{m} \left(\zeta \right)` and :math:`\zeta _{h} =\zeta _{w} =-0.465` for :math:`\phi _{h} \left(\zeta \right)=\phi _{w} \left(\zeta \right)` (:ref:`Zeng et al. 1998 <Zengetal1998>`). The flux gradient relations can be integrated to yield wind profiles for the following conditions:

Very unstable :math:`\left(\zeta <-1.574\right)`

.. math::
   :label: 5.32

   V_{a} =\frac{u_{*} }{k} \left\{\left[\ln \frac{\zeta _{m} L}{z_{0m} } -\psi _{m} \left(\zeta _{m} \right)\right]+1.14\left[\left(-\zeta \right)^{{1\mathord{\left/ {\vphantom {1 3}} \right.} 3} } -\left(-\zeta _{m} \right)^{{1\mathord{\left/ {\vphantom {1 3}} \right.} 3} } \right]+\psi _{m} \left(\frac{z_{0m} }{L} \right)\right\}

Unstable :math:`\left(-1.574\le \zeta <0\right)`

.. math::
   :label: 5.33

   V_{a} =\frac{u_{*} }{k} \left\{\left[\ln \frac{z_{atm,\, m} -d}{z_{0m} } -\psi _{m} \left(\zeta \right)\right]+\psi _{m} \left(\frac{z_{0m} }{L} \right)\right\}

Stable :math:`\left(0\le \zeta \le 1\right)`

.. math::
   :label: 5.34

   V_{a} =\frac{u_{*} }{k} \left\{\left[\ln \frac{z_{atm,\, m} -d}{z_{0m} } +5\zeta \right]-5\frac{z_{0m} }{L} \right\}

Very stable :math:`\left(\zeta >1\right)`

.. math::
   :label: 5.35

   V_{a} =\frac{u_{*} }{k} \left\{\left[\ln \frac{L}{z_{0m} } +5\right]+\left[5\ln \zeta +\zeta -1\right]-5\frac{z_{0m} }{L} \right\}

where

.. math::
   :label: 5.36

   \psi _{m} \left(\zeta \right)=2\ln \left(\frac{1+x}{2} \right)+\ln \left(\frac{1+x^{2} }{2} \right)-2\tan ^{-1} x+\frac{\pi }{2}

and

:math:`x=\left(1-16\zeta \right)^{{1\mathord{\left/ {\vphantom {1 4}} \right.} 4} }` .

The potential temperature profiles are:

Very unstable :math:`\left(\zeta <-0.465\right)`

.. math::
   :label: 5.37

   \theta _{atm} -\theta _{s} =\frac{\theta _{*} }{k} \left\{\left[\ln \frac{\zeta _{h} L}{z_{0h} } -\psi _{h} \left(\zeta _{h} \right)\right]+0.8\left[\left(-\zeta _{h} \right)^{{-1\mathord{\left/ {\vphantom {-1 3}} \right.} 3} } -\left(-\zeta \right)^{{-1\mathord{\left/ {\vphantom {-1 3}} \right.} 3} } \right]+\psi _{h} \left(\frac{z_{0h} }{L} \right)\right\}

Unstable :math:`\left(-0.465\le \zeta <0\right)`

.. math::
   :label: 5.38

   \theta _{atm} -\theta _{s} =\frac{\theta _{*} }{k} \left\{\left[\ln \frac{z_{atm,\, h} -d}{z_{0h} } -\psi _{h} \left(\zeta \right)\right]+\psi _{h} \left(\frac{z_{0h} }{L} \right)\right\}

Stable :math:`\left(0\le \zeta \le 1\right)`

.. math::
   :label: 5.39

   \theta _{atm} -\theta _{s} =\frac{\theta _{*} }{k} \left\{\left[\ln \frac{z_{atm,\, h} -d}{z_{0h} } +5\zeta \right]-5\frac{z_{0h} }{L} \right\}

Very stable :math:`\left(\zeta >1\right)`

.. math::
   :label: 5.40

   \theta _{atm} -\theta _{s} =\frac{\theta _{*} }{k} \left\{\left[\ln \frac{L}{z_{0h} } +5\right]+\left[5\ln \zeta +\zeta -1\right]-5\frac{z_{0h} }{L} \right\}.

The specific humidity profiles are:

Very unstable :math:`\left(\zeta <-0.465\right)`

.. math::
   :label: 5.41

   q_{atm} -q_{s} =\frac{q_{*} }{k} \left\{\left[\ln \frac{\zeta _{w} L}{z_{0w} } -\psi _{w} \left(\zeta _{w} \right)\right]+0.8\left[\left(-\zeta _{w} \right)^{{-1\mathord{\left/ {\vphantom {-1 3}} \right.} 3} } -\left(-\zeta \right)^{{-1\mathord{\left/ {\vphantom {-1 3}} \right.} 3} } \right]+\psi _{w} \left(\frac{z_{0w} }{L} \right)\right\}

Unstable :math:`\left(-0.465\le \zeta <0\right)`

.. math::
   :label: 5.42

   q_{atm} -q_{s} =\frac{q_{*} }{k} \left\{\left[\ln \frac{z_{atm,\, w} -d}{z_{0w} } -\psi _{w} \left(\zeta \right)\right]+\psi _{w} \left(\frac{z_{0w} }{L} \right)\right\}

Stable :math:`\left(0\le \zeta \le 1\right)`

.. math::
   :label: 5.43

   q_{atm} -q_{s} =\frac{q_{*} }{k} \left\{\left[\ln \frac{z_{atm,\, w} -d}{z_{0w} } +5\zeta \right]-5\frac{z_{0w} }{L} \right\}

Very stable :math:`\left(\zeta >1\right)`

.. math::
   :label: 5.44

   q_{atm} -q_{s} =\frac{q_{*} }{k} \left\{\left[\ln \frac{L}{z_{0w} } +5\right]+\left[5\ln \zeta +\zeta -1\right]-5\frac{z_{0w} }{L} \right\}

where

.. math::
   :label: 5.45

   \psi _{h} \left(\zeta \right)=\psi _{w} \left(\zeta \right)=2\ln \left(\frac{1+x^{2} }{2} \right).

Using the definitions of :math:`u_{*}`, :math:`\theta _{*}`, and :math:`q_{*}`, an iterative solution of these equations can be used to calculate the surface momentum, sensible heat, and water vapor flux using atmospheric and surface values for :math:`\left|{\it u}\right|`, :math:`\theta`, and :math:`q` except that :math:`L` depends on :math:`u_{*}`, :math:`\theta _{*}`, and :math:`q_{*}`. However, the bulk Richardson number

.. math::
   :label: 5.46

   R_{iB} =\frac{\theta _{v,\, atm} -\theta _{v,\, s} }{\overline{\theta _{v,\, atm} }} \frac{g\left(z_{atm,\, m} -d\right)}{V_{a}^{2} }

is related to :math:`\zeta` (:ref:`Arya 2001 <Arya2001>`) as

.. math::
   :label: 5.47

   R_{iB} =\zeta \left[\ln \left(\frac{z_{atm,\, h} -d}{z_{0h} } \right)-\psi _{h} \left(\zeta \right)\right]\left[\ln \left(\frac{z_{atm,\, m} -d}{z_{0m} } \right)-\psi _{m} \left(\zeta \right)\right]^{-2} .

Using :math:`\phi _{h} =\phi _{m}^{2} =\left(1-16\zeta \right)^{-{1\mathord{\left/ {\vphantom {1 2}} \right.} 2} }` for unstable conditions and :math:`\phi _{h} =\phi _{m} =1+5\zeta` for stable conditions to determine :math:`\psi _{m} \left(\zeta \right)` and :math:`\psi _{h} \left(\zeta \right)`, the inverse relationship :math:`\zeta =f\left(R_{iB} \right)` can be solved to obtain a first guess for :math:`\zeta` and thus :math:`L` from

.. math::
   :label: 5.48

   \begin{array}{lcr}
   \zeta =\frac{R_{iB} \ln \left(\frac{z_{atm,\, m} -d}{z_{0m} } \right)}{1-5\min \left(R_{iB} ,0.19\right)} & \qquad 0.01\le \zeta \le 2 & \qquad {\rm for\; }R_{iB} \ge 0 {\rm \; (neutral\; or\; stable)} \\
   \zeta =R_{iB} \ln \left(\frac{z_{atm,\, m} -d}{z_{0m} } \right) & \qquad -100\le \zeta \le -0.01 & \qquad {\rm for\; }R_{iB} <0 \ {\rm \; (unstable)}
   \end{array}.

Upon iteration (section :numref:`Numerical Implementation`), the following is used to determine :math:`\zeta` and thus :math:`L`

.. math::
   :label: 5.49

   \zeta =\frac{\left(z_{atm,\, m} -d\right)kg\theta _{v*} }{u_{*}^{2} \overline{\theta _{v,\, atm} }}

where

.. math::

   \begin{array}{cr}
   0.01\le \zeta \le 2 & \qquad {\rm for\; }\zeta \ge 0{\rm \; (neutral\; or\; stable)} \\
   {\rm -100}\le \zeta \le {\rm -0.01} & \qquad {\rm for\; }\zeta <0{\rm \; (unstable)}
   \end{array}.

The difference in virtual potential air temperature between the reference height and the surface is

.. math::
   :label: 5.50

   \theta _{v,\, atm} -\theta _{v,\, s} =\left(\theta _{atm} -\theta _{s} \right)\left(1+0.61q_{atm} \right)+0.61\overline{\theta _{atm} }\left(q_{atm} -q_{s} \right).

The momentum, sensible heat, and water vapor fluxes between the surface and the atmosphere can also be written in the form

.. math::
   :label: 5.51

   \tau _{x} =-\rho _{atm} \frac{\left(u_{atm} -u_{s} \right)}{r_{am} }

.. math::
   :label: 5.52

   \tau _{y} =-\rho _{atm} \frac{\left(v_{atm} -v_{s} \right)}{r_{am} }

.. math::
   :label: 5.53

   H=-\rho _{atm} C_{p} \frac{\left(\theta _{atm} -\theta _{s} \right)}{r_{ah} }

.. math::
   :label: 5.54

   E=-\rho _{atm} \frac{\left(q_{atm} -q_{s} \right)}{r_{aw} }

where the aerodynamic resistances (s m\ :sup:`-1`) are

.. math::
   :label: 5.55

   r_{am} =\frac{V_{a} }{u_{*}^{2} } =\frac{1}{k^{2} V_{a} } \left[\ln \left(\frac{z_{atm,\, m} -d}{z_{0m} } \right)-\psi _{m} \left(\frac{z_{atm,\, m} -d}{L} \right)+\psi _{m} \left(\frac{z_{0m} }{L} \right)\right]^{2}

.. math::
   :label: 5.56

   \begin{array}{l} {r_{ah} =\frac{\theta _{atm} -\theta _{s} }{\theta _{*} u_{*} } =\frac{1}{k^{2} V_{a} } \left[\ln \left(\frac{z_{atm,\, m} -d}{z_{0m} } \right)-\psi _{m} \left(\frac{z_{atm,\, m} -d}{L} \right)+\psi _{m} \left(\frac{z_{0m} }{L} \right)\right]} \\ {\qquad \left[\ln \left(\frac{z_{atm,\, h} -d}{z_{0h} } \right)-\psi _{h} \left(\frac{z_{atm,\, h} -d}{L} \right)+\psi _{h} \left(\frac{z_{0h} }{L} \right)\right]} \end{array}

.. math::
   :label: 5.57

   \begin{array}{l} {r_{aw} =\frac{q_{atm} -q_{s} }{q_{*} u_{*} } =\frac{1}{k^{2} V_{a} } \left[\ln \left(\frac{z_{atm,\, m} -d}{z_{0m} } \right)-\psi _{m} \left(\frac{z_{atm,\, m} -d}{L} \right)+\psi _{m} \left(\frac{z_{0m} }{L} \right)\right]} \\ {\qquad \left[\ln \left(\frac{z_{atm,\, {\it w}} -d}{z_{0w} } \right)-\psi _{w} \left(\frac{z_{atm,\, w} -d}{L} \right)+\psi _{w} \left(\frac{z_{0w} }{L} \right)\right]} \end{array}.

A 2-m height "screen" temperature is useful for comparison with observations

.. math::
   :label: 5.58

   T_{2m} =\theta _{s} +\frac{\theta _{*} }{k} \left[\ln \left(\frac{2+z_{0h} }{z_{0h} } \right)-\psi _{h} \left(\frac{2+z_{0h} }{L} \right)+\psi _{h} \left(\frac{z_{0h} }{L} \right)\right]

where for convenience, "2-m" is defined as 2 m above the apparent sink for sensible heat (:math:`z_{0h} +d`). Similarly, a 2-m height specific humidity is defined as

.. math::
   :label: 5.59

   q_{2m} =q_{s} +\frac{q_{*} }{k} \left[\ln \left(\frac{2+z_{0w} }{z_{0w} } \right)-\psi _{w} \left(\frac{2+z_{0w} }{L} \right)+\psi _{w} \left(\frac{z_{0w} }{L} \right)\right].

Relative humidity is

.. math::
   :label: 5.60

   RH_{2m} =\min \left(100,\, \frac{q_{2m} }{q_{sat}^{T_{2m} } } \times 100\right)

where :math:`q_{sat}^{T_{2m} }` is the saturated specific humidity at the 2-m temperature :math:`T_{2m}` (section :numref:`Saturation Vapor Pressure`).

A 10-m wind speed is calculated as (note that this is not consistent with the 10-m wind speed calculated for the dust model as described in Chapter :numref:`rst_Dust Model`)

.. math::
   :label: 5.61

   u_{10m} =\left\{\begin{array}{l} {V_{a} \qquad z_{atm,\, m} \le 10} \\ {V_{a} -\frac{u_{*} }{k} \left[\ln \left(\frac{z_{atm,\, m} -d}{10+z_{0m} } \right)-\psi _{m} \left(\frac{z_{atm,\, m} -d}{L} \right)+\psi _{m} \left(\frac{10+z_{0m} }{L} \right)\right]\qquad z_{atm,\, m} >10} \end{array}\right\}

.. _Sensible and Latent Heat Fluxes for Non-Vegetated Surfaces:

Sensible and Latent Heat Fluxes for Non-Vegetated Surfaces
--------------------------------------------------------------

Surfaces are considered non-vegetated for the surface flux calculations if leaf plus stem area index :math:`L+S<0.05` (section :numref:`Phenology and vegetation burial by snow`). By definition, this includes bare soil and glaciers. The solution for lakes is described in Chapter :numref:`rst_Lake Model`. For these surfaces, the surface may be exposed to the atmosphere, snow covered, and/or surface water covered, so that the sensible heat flux :math:`H_{g}` (W m\ :sup:`-2`) is, with reference to :numref:`Figure Schematic diagram of sensible heat fluxes`,

.. math::
   :label: 5.62

   H_{g} =\left(1-f_{sno} -f_{h2osfc} \right)H_{soil} +f_{sno} H_{snow} +f_{h2osfc} H_{h2osfc}

where :math:`\left(1-f_{sno} -f_{h2osfc} \right)`, :math:`f_{sno}`, and :math:`f_{h2osfc}` are the exposed, snow covered, and surface water covered fractions of the grid cell. The individual fluxes based on the temperatures of the soil :math:`T_{1}`, snow :math:`T_{snl+1}`, and surface water :math:`T_{h2osfc}` are

.. math::
   :label: 5.63

   H_{soil} =-\rho _{atm} C_{p} \frac{\left(\theta _{atm} -T_{1} \right)}{r_{ah} }

.. math::
   :label: 5.64

   H_{sno} =-\rho _{atm} C_{p} \frac{\left(\theta _{atm} -T_{snl+1} \right)}{r_{ah} }

.. math::
   :label: 5.65

   H_{h2osfc} =-\rho _{atm} C_{p} \frac{\left(\theta _{atm} -T_{h2osfc} \right)}{r_{ah} }

where :math:`\rho _{atm}` is the density of atmospheric air (kg m\ :sup:`-3`), :math:`C_{p}` is the specific heat capacity of air (J kg\ :sup:`-1` K\ :sup:`-1`) (:numref:`Table Physical constants`), :math:`\theta _{atm}` is the atmospheric potential temperature (K), and :math:`r_{ah}` is the aerodynamic resistance to sensible heat transfer (s m\ :sup:`-1`).

The water vapor flux :math:`E_{g}` (kg m\ :sup:`-2` s\ :sup:`-1`) is, with reference to :numref:`Figure Schematic diagram of latent heat fluxes`,

.. math::
   :label: 5.66

   E_{g} =\left(1-f_{sno} -f_{h2osfc} \right)E_{soil} +f_{sno} E_{snow} +f_{h2osfc} E_{h2osfc}

.. math::
   :label: 5.67

   E_{soil} =-\frac{\rho _{atm} \left(q_{atm} -q_{soil} \right)}{r_{aw} + r_{soil}}

.. math::
   :label: 5.68

   E_{sno} =-\frac{\rho _{atm} \left(q_{atm} -q_{sno} \right)}{r_{aw} }

.. math::
   :label: 5.69

   E_{h2osfc} =-\frac{\rho _{atm} \left(q_{atm} -q_{h2osfc} \right)}{r_{aw} }

where :math:`q_{atm}` is the atmospheric specific humidity (kg kg\ :sup:`-1`), :math:`q_{soil}`, :math:`q_{sno}`, and :math:`q_{h2osfc}` are the specific humidities (kg kg\ :sup:`-1`) of the soil, snow, and surface water, respectively, :math:`r_{aw}` is the aerodynamic resistance to water vapor transfer (s m\ :sup:`-1`), and :math:`r _{soi}` is the soil resistance to water vapor transfer (s m\ :sup:`-1`). The specific humidities of the snow :math:`q_{sno}` and surface water :math:`q_{h2osfc}` are assumed to be at the saturation specific humidity of their respective temperatures

.. math::
   :label: 5.70

   q_{sno} =q_{sat}^{T_{snl+1} }

.. math::
   :label: 5.71

   q_{h2osfc} =q_{sat}^{T_{h2osfc} }

The specific humidity of the soil surface :math:`q_{soil}` is assumed to be proportional to the saturation specific humidity

.. math::
   :label: 5.72

   q_{soil} =\alpha _{soil} q_{sat}^{T_{1} }

where :math:`q_{sat}^{T_{1} }` is the saturated specific humidity at the soil surface temperature :math:`T_{1}` (section :numref:`Saturation Vapor Pressure`). The factor :math:`\alpha _{soil}` is a function of the surface soil water matric potential :math:`\psi` as in :ref:`Philip (1957)<Philip1957>`

.. math::
   :label: 5.73

   \alpha _{soil} =\exp \left(\frac{\psi _{1} g}{1\times 10^{3} R_{wv} T_{1} } \right)

where :math:`R_{wv}` is the gas constant for water vapor (J kg\ :sup:`-1` K\ :sup:`-1`) (:numref:`Table Physical constants`), :math:`g` is the gravitational acceleration (m s\ :sup:`-2`) (:numref:`Table Physical constants`), and :math:`\psi _{1}` is the soil water matric potential of the top soil layer (mm). The soil water matric potential :math:`\psi _{1}` is

.. math::
   :label: 5.74

   \psi _{1} =\psi _{sat,\, 1} s_{1}^{-B_{1} } \ge -1\times 10^{8}

where :math:`\psi _{sat,\, 1}` is the saturated matric potential (mm) (section :numref:`Hydraulic Properties`), :math:`B_{1}` is the :ref:`Clapp and Hornberger (1978) <ClappHornberger1978>` parameter (section :numref:`Hydraulic Properties`), and :math:`s_{1}` is the wetness of the top soil layer with respect to saturation. The surface wetness :math:`s_{1}` is a function of the liquid water and ice content

.. math::
   :label: 5.75

   s_{1} =\frac{1}{\Delta z_{1} \theta _{sat,\, 1} } \left[\frac{w_{liq,\, 1} }{\rho _{liq} } +\frac{w_{ice,\, 1} }{\rho _{ice} } \right]\qquad 0.01\le s_{1} \le 1.0

where :math:`\Delta z_{1}` is the thickness of the top soil layer (m), :math:`\rho _{liq}` and :math:`\rho _{ice}` are the density of liquid water and ice (kg m\ :sup:`-3`) (:numref:`Table Physical constants`), :math:`w_{liq,\, 1}` and :math:`w_{ice,\, 1}` are the mass of liquid water and ice of the top soil layer (kg m\ :sup:`-2`) (Chapter :numref:`rst_Hydrology`), and :math:`\theta _{sat,\, 1}` is the saturated volumetric water content (i.e., porosity) of the top soil layer (mm\ :sup:`3` mm\ :sup:`-3`) (section :numref:`Hydraulic Properties`). If :math:`q_{sat}^{T_{1} } >q_{atm}` and :math:`q_{atm} >q_{soil}`, then :math:`q_{soil} =q_{atm}` and :math:`\frac{dq_{soil} }{dT} =0`. This prevents large increases (decreases) in :math:`q_{soil}` for small increases (decreases) in soil moisture in very dry soils.

The resistance to water vapor transfer occurring within the soil matrix :math:`r_{soil}` (s m\ :sup:`-1`) is

.. math::
   :label: 5.76

   r_{soil} = \frac{DSL}{D_{v} \tau}

where :math:`DSL` is the thickness of the dry surface layer (m), :math:`D_{v}` is the molecular diffusivity of water vapor in air (m\ :sup:`2` s\ :sup:`-2`) and :math:`\tau` (*unitless*) describes the tortuosity of the vapor flow paths through the soil matrix (:ref:`Swenson and Lawrence 2014 <SwensonLawrence2014>`).

The thickness of the dry surface layer is given by

.. math::
   :label: 5.77

   DSL =
   \begin{array}{lr}
   D_{max} \ \frac{\left( \theta_{init} - \theta_{1}\right)}
   {\left(\theta_{init} - \theta_{air}\right)} & \qquad \theta_{1} < \theta_{init} \\
   0 &  \qquad \theta_{1} \ge \theta_{init}
   \end{array}

where :math:`D_{max}` is a parameter specifying the length scale of the maximum DSL thickness (default value = 15 mm), :math:`\theta_{init}` (mm\ :sup:`3` mm\ :sup:`-3`) is the moisture value at which the DSL initiates, :math:`\theta_{1}` (mm\ :sup:`3` mm\ :sup:`-3`) is the moisture value of the top model soil layer, and :math:`\theta_{air}` (mm\ :sup:`3` mm\ :sup:`-3`) is the 'air dry' soil moisture value (:ref:`Dingman 2002 <Dingman2002>`):

.. math::
   :label: 5.78

   \theta_{air} = \Phi \left( \frac{\Psi_{sat}}{\Psi_{air}} \right)^{\frac{1}{B_{1}}} \ .

where :math:`\Phi` is the porosity (mm\ :sup:`3` mm\ :sup:`-3`), :math:`\Psi_{sat}` is the saturated soil matric potential (mm), :math:`\Psi_{air} = 10^{7}` mm is the air dry matric potential, and :math:`B_{1}` is a function of soil texture (section :numref:`Hydraulic Properties`).

The soil tortuosity is

.. math::
   :label: 5.79

   \tau = \Phi^{2}_{air}\left(\frac{\Phi_{air}}{\Phi}\right)^{\frac{3}{B_{1}}}

where :math:`\Phi_{air}` (mm\ :sup:`3` mm\ :sup:`-3`) is the air filled pore space

.. math::
   :label: 5.80

   \Phi_{air} = \Phi - \theta_{air} \ .

:math:`D_{v}` depends on temperature

.. math::
   :label: 5.81

   D_{v} = 2.12 \times 10^{-5} \left(\frac{T_{1}}{T_{f}}\right)^{1.75} \ .

where :math:`T_{1}` (K) is the temperature of the top soil layer and :math:`T_{f}` (K) is the freezing temperature of water (:numref:`Table Physical Constants`).

The roughness lengths used to calculate :math:`r_{am}`, :math:`r_{ah}`, and :math:`r_{aw}` are :math:`z_{0m} =z_{0m,\, g}`, :math:`z_{0h} =z_{0h,\, g}`, and :math:`z_{0w} =z_{0w,\, g}`. The displacement height :math:`d=0`. The momentum roughness length is :math:`z_{0m,\, g} =0.0023` for glaciers without snow (:math:`f_{sno} =0) {\rm }`, and :math:`z_{0m,\, g} =0.00085` for bare soil surfaces without snow (:math:`f_{sno} =0) {\rm }` (:ref:`Meier et al. (2022) <Meieretal2022>`).

For bare soil and glaciers with snow ( :math:`f_{sno} > 0` ), the momentum roughness length is evaluated based on accumulated snow melt :math:`M_{a} {\rm }` (:ref:`Meier et al. (2022) <Meieretal2022>`). For :math:`M_{a} >=1\times 10^{-5}`

.. math::
   :label: 5.81a

   z_{0m,\, g} =\exp (b_{1} \tan ^{-1} \left[\frac{log_{10} (M_{a}) + 0.23)} {0.08}\right] + b_{4})\times 10^{-3}

where :math:`M_{a}` is accumulated snow melt (meters water equivalent), :math:`b_{1} =1.4` and :math:`b_{4} =-0.31`. For :math:`M_{a} <1\times 10^{-5}`

.. math::
   :label: 5.81b

   z_{0m,\, g} =\exp (-b_{1} 0.5 \pi + b_{4})\times 10^{-3}

Accumulated snow melt :math:`M_{a}` at the current time step :math:`t` is defined as

.. math::
   :label: 5.81c

   M ^{t}_{a} = M ^{t-1}_{a} - (q ^{t}_{sno} \Delta t + q ^{t}_{snowmelt} \Delta t)\times 10^{-3}

where :math:`M ^{t}_{a}` and :math:`M ^{t-1}_{a}` are the accumulated snowmelt at the current time step and previous time step, respectively (m), :math:`q ^{t}_{sno} \Delta t` is the freshly fallen snow (mm), and :math:`q ^{t}_{snowmelt} \Delta t` is the melted snow (mm).

The scalar roughness lengths (:math:`z_{0q,\, g}` for latent heat and :math:`z_{0h,\ g}` for sensible heat) are calculated as (:ref:`Meier et al. (2022) <Meieretal2022>`)

.. math::
   :label: 5.82

   z_{0h,\, g}=z_{0q,\, g}=\frac{70 \nu}{u_{*}} \exp (-\beta {u_{*}} ^{0.5} |{\theta_{*}}| ^{0.25} )

where :math:`\beta` = 7.2, and :math:`\theta_{*}` is the potential temperature scale.

The numerical solution for the fluxes of momentum, sensible heat, and water vapor flux from non-vegetated surfaces proceeds as follows:

#. An initial guess for the wind speed :math:`V_{a}`  is obtained from
   :eq:`5.24` assuming an initial convective velocity :math:`U_{c} =0` m
   s\ :sup:`-1` for stable conditions
   (:math:`\theta _{v,\, atm} -\theta _{v,\, s} \ge 0` as evaluated from
   :eq:`5.50` ) and :math:`U_{c} =0.5` for unstable conditions
   (:math:`\theta _{v,\, atm} -\theta _{v,\, s} <0`).

#. An initial guess for the Monin-Obukhov length :math:`L` is obtained
   from the bulk Richardson number using :eq:`5.46` and :eq:`5.48`.

#. The following system of equations is iterated three times:

#. Friction velocity :math:`u_{*}`  (:eq:`5.32`, :eq:`5.33`, :eq:`5.34`, :eq:`5.35`)

#. Potential temperature scale :math:`\theta _{*}`  (:eq:`5.37` , :eq:`5.38`, :eq:`5.39`, :eq:`5.40`)

#. Humidity scale :math:`q_{*}`  (:eq:`5.41`, :eq:`5.42`, :eq:`5.43`, :eq:`5.44`)

#. Roughness lengths for sensible :math:`z_{0h,\, g}`  and latent heat
   :math:`z_{0w,\, g}`  (:eq:`5.81a` , :eq:`5.81b` , :eq:`5.82`)

#. Virtual potential temperature scale :math:`\theta _{v*}`  ( :eq:`5.17`)

#. Wind speed including the convective velocity, :math:`V_{a}`  ( :eq:`5.24`)

#. Monin-Obukhov length :math:`L` (:eq:`5.49`)

#. Aerodynamic resistances :math:`r_{am}` , :math:`r_{ah}` , and
   :math:`r_{aw}`  (:eq:`5.55`, :eq:`5.56`, :eq:`5.57`)

#. Momentum fluxes :math:`\tau _{x}` , :math:`\tau _{y}`  (:eq:`5.5`, :eq:`5.6`)

#. Sensible heat flux :math:`H_{g}`  (:eq:`5.62`)

#. Water vapor flux :math:`E_{g}`  (:eq:`5.66`)

#. 2-m height air temperature :math:`T_{2m}`  and specific humidity
   :math:`q_{2m}`  (:eq:`5.58` , :eq:`5.59`)

The partial derivatives of the soil surface fluxes with respect to ground temperature, which are needed for the soil temperature calculations (section :numref:`Numerical Solution Temperature`) and to update the soil surface fluxes (section :numref:`Update of Ground Sensible and Latent Heat Fluxes`), are

.. math::
   :label: 5.83

   \frac{\partial H_{g} }{\partial T_{g} } =\frac{\rho _{atm} C_{p} }{r_{ah} }

.. math::
   :label: 5.84

   \frac{\partial E_{g} }{\partial T_{g} } =\frac{\beta _{soi} \rho _{atm} }{r_{aw} } \frac{dq_{g} }{dT_{g} }

where

.. math::
   :label: 5.85

   \frac{dq_{g} }{dT_{g} } =\left(1-f_{sno} -f_{h2osfc} \right)\alpha _{soil} \frac{dq_{sat}^{T_{soil} } }{dT_{soil} } +f_{sno} \frac{dq_{sat}^{T_{sno} } }{dT_{sno} } +f_{h2osfc} \frac{dq_{sat}^{T_{h2osfc} } }{dT_{h2osfc} } .

The partial derivatives :math:`\frac{\partial r_{ah} }{\partial T_{g} }` and :math:`\frac{\partial r_{aw} }{\partial T_{g} }`, which cannot be determined analytically, are ignored for :math:`\frac{\partial H_{g} }{\partial T_{g} }` and :math:`\frac{\partial E_{g} }{\partial T_{g} }`.

.. _Sensible and Latent Heat Fluxes and Temperature for Vegetated Surfaces:

Sensible and Latent Heat Fluxes and Temperature for Vegetated Surfaces
--------------------------------------------------------------------------

In the case of a vegetated surface, the sensible heat :math:`H` and water vapor flux :math:`E` are partitioned into vegetation and ground fluxes that depend on vegetation :math:`T_{v}` and ground :math:`T_{g}` temperatures in addition to surface temperature :math:`T_{s}` and specific humidity :math:`q_{s}`. Because of the coupling between vegetation temperature and fluxes, Newton-Raphson iteration is used to solve for the vegetation temperature and the sensible heat and water vapor fluxes from vegetation simultaneously using the ground temperature from the previous time step. In section :numref:`Theory`, the equations used in the iteration scheme are derived. Details on the numerical scheme are provided in section :numref:`Numerical Implementation`.

.. _Theory:

Theory
^^^^^^^^^^^^

The air within the canopy is assumed to have negligible capacity to store heat so that the sensible heat flux :math:`H` between the surface at height :math:`z_{0h} +d` and the atmosphere at height :math:`z_{atm,\, h}` must be balanced by the sum of the sensible heat from the vegetation :math:`H_{v}` and the ground :math:`H_{g}`

.. math::
   :label: 5.86

   H=H_{v} +H_{g}

where, with reference to :numref:`Figure Schematic diagram of sensible heat fluxes`,

.. math::
   :label: 5.87

   H=-\rho _{atm} C_{p} \frac{\left(\theta _{atm} -T_{s} \right)}{r_{ah} }

.. math::
   :label: 5.88

   H_{v} =-\rho _{atm} C_{p} \left(T_{s} -T_{v} \right)\frac{\left(L+S\right)}{r_{b} }

.. math::
   :label: 5.89

   H_{g} =\left(1-f_{sno} -f_{h2osfc} \right)H_{soil} +f_{sno} H_{snow} +f_{h2osfc} H_{h2osfc} \ ,

where

.. math::
   :label: 5.90

   H_{soil} =-\rho _{atm} C_{p} \frac{\left(T_{s} -T_{1} \right)}{r_{ah} ^{{'} } }

.. math::
   :label: 5.91

   H_{sno} =-\rho _{atm} C_{p} \frac{\left(T_{s} -T_{snl+1} \right)}{r_{ah} ^{{'} } }

.. math::
   :label: 5.92

   H_{h2osfc} =-\rho _{atm} C_{p} \frac{\left(T_{s} -T_{h2osfc} \right)}{r_{ah} ^{{'} } }

where :math:`\rho _{atm}` is the density of atmospheric air (kg m\ :sup:`-3`), :math:`C_{p}` is the specific heat capacity of air (J kg\ :sup:`-1` K\ :sup:`-1`) (:numref:`Table Physical constants`), :math:`\theta _{atm}` is the atmospheric potential temperature (K), and :math:`r_{ah}` is the aerodynamic resistance to sensible heat transfer (s m\ :sup:`-1`).

Here, :math:`T_{s}` is the surface temperature at height :math:`z_{0h} +d`, also referred to as the canopy air temperature. :math:`L` and :math:`S` are the exposed leaf and stem area indices (section :numref:`Phenology and vegetation burial by snow`), :math:`r_{b}` is the leaf boundary layer resistance (s m\ :sup:`-1`), and :math:`r_{ah} ^{{'} }` is the aerodynamic resistance (s m\ :sup:`-1`) to heat transfer between the ground at height :math:`z_{0h} ^{{'} }` and the canopy air at height :math:`z_{0h} +d`.

.. _Figure Schematic diagram of sensible heat fluxes:

.. figure:: image1.png

 Figure Schematic diagram of sensible heat fluxes for (a)
 non-vegetated surfaces and (b) vegetated surfaces.

.. _Figure Schematic diagram of latent heat fluxes:

.. figure:: image2.png

 Figure Schematic diagram of water vapor fluxes for (a)
 non-vegetated surfaces and (b) vegetated surfaces.

Equations :eq:`5.86` - :eq:`5.89` can be solved for the canopy air temperature :math:`T_{s}`

.. math::
   :label: 5.93

   T_{s} =\frac{c_{a}^{h} \theta _{atm} +c_{g}^{h} T_{g} +c_{v}^{h} T_{v} }{c_{a}^{h} +c_{g}^{h} +c_{v}^{h} }

where

.. math::
   :label: 5.94

   c_{a}^{h} =\frac{1}{r_{ah} }

.. math::
   :label: 5.95

   c_{g}^{h} =\frac{1}{r_{ah} ^{{'} } }

.. math::
   :label: 5.96

   c_{v}^{h} =\frac{\left(L+S\right)}{r_{b} }

are the sensible heat conductances from the canopy air to the atmosphere, the ground to canopy air, and leaf surface to canopy air, respectively (m s\ :sup:`-1`).

When the expression for :math:`T_{s}` is substituted into equation :eq:`5.88`, the sensible heat flux from vegetation :math:`H_{v}` is a function of :math:`\theta _{atm}`, :math:`T_{g}`, and :math:`T_{v}`

.. math::
   :label: 5.97

   H_{v} = -\rho _{atm} C_{p} \left[c_{a}^{h} \theta _{atm} +c_{g}^{h} T_{g} -\left(c_{a}^{h} +c_{g}^{h} \right)T_{v} \right]\frac{c_{v}^{h} }{c_{a}^{h} +c_{v}^{h} +c_{g}^{h} } .

Similarly, the expression for :math:`T_{s}` can be substituted into equations :eq:`5.89`, :eq:`5.90`, :eq:`5.91`, and :eq:`5.92` to obtain the sensible heat flux from ground :math:`H_{g}`

.. math::
   :label: 5.98

   H_{g} = -\rho _{atm} C_{p} \left[c_{a}^{h} \theta _{atm} +c_{v}^{h} T_{v} -\left(c_{a}^{h} +c_{v}^{h} \right)T_{g} \right]\frac{c_{g}^{h} }{c_{a}^{h} +c_{v}^{h} +c_{g}^{h} } .

The air within the canopy is assumed to have negligible capacity to store water vapor so that the water vapor flux :math:`E` between the surface at height :math:`z_{0w} +d` and the atmosphere at height :math:`z_{atm,\, w}` must be balanced by the sum of the water vapor flux from the vegetation :math:`E_{v}` and the ground :math:`E_{g}`

.. math::
   :label: 5.99

   E = E_{v} +E_{g}

where, with reference to :numref:`Figure Schematic diagram of latent heat fluxes`,

.. math::
   :label: 5.100

   E = -\rho _{atm} \frac{\left(q_{atm} -q_{s} \right)}{r_{aw} }

.. math::
   :label: 5.101

   E_{v} = -\rho _{atm} \frac{\left(q_{s} -q_{sat}^{T_{v} } \right)}{r_{total} }

.. math::
   :label: 5.102

   E_{g} = \left(1-f_{sno} -f_{h2osfc} \right)E_{soil} +f_{sno} E_{snow} +f_{h2osfc} E_{h2osfc} \ ,

where

.. math::
   :label: 5.103

   E_{soil} = -\rho _{atm} \frac{\left(q_{s} -q_{soil} \right)}{r_{aw} ^{{'} } +r_{soil} }

.. math::
   :label: 5.104

   E_{sno} = -\rho _{atm} \frac{\left(q_{s} -q_{sno} \right)}{r_{aw} ^{{'} } +r_{soil} }

.. math::
   :label: 5.105

   E_{h2osfc} = -\rho _{atm} \frac{\left(q_{s} -q_{h2osfc} \right)}{r_{aw} ^{{'} } +r_{soil} }

where :math:`q_{atm}` is the atmospheric specific humidity (kg kg\ :sup:`-1`), :math:`r_{aw}` is the aerodynamic resistance to water vapor transfer (s m\ :sup:`-1`), :math:`q_{sat}^{T_{v} }` (kg kg\ :sup:`-1`) is the saturation water vapor specific humidity at the vegetation temperature (section :numref:`Saturation Vapor Pressure`), :math:`q_{g}`, :math:`q_{sno}`, and :math:`q_{h2osfc}` are the specific humidities of the soil, snow, and surface water (section :numref:`Sensible and Latent Heat Fluxes for Non-Vegetated Surfaces`), :math:`r_{aw} ^{{'} }` is the aerodynamic resistance (s m\ :sup:`-1`) to water vapor transfer between the ground at height :math:`z_{0w} ^{{'} }` and the canopy air at height :math:`z_{0w} +d`, and :math:`r_{soil}` (:eq:`5.76`) is a resistance to diffusion through the soil (s m\ :sup:`-1`). :math:`r_{total}` is the total resistance to water vapor transfer from the canopy to the canopy air and includes contributions from leaf boundary layer and sunlit and shaded stomatal resistances :math:`r_{b}`, :math:`r_{s}^{sun}`, and :math:`r_{s}^{sha}` (:numref:`Figure Schematic diagram of latent heat fluxes`). The water vapor flux from vegetation is the sum of water vapor flux from wetted leaf and stem area :math:`E_{v}^{w}` (evaporation of water intercepted by the canopy) and transpiration from dry leaf surfaces :math:`E_{v}^{t}`

.. math::
   :label: 5.106

   E_{v} =E_{v}^{w} +E_{v}^{t} .

Equations :eq:`5.99` - :eq:`5.102` can be solved for the canopy specific humidity :math:`q_{s}`

.. math::
   :label: 5.107

   q_{s} =\frac{c_{a}^{w} q_{atm} +c_{g}^{w} q_{g} +c_{v}^{w} q_{sat}^{T_{v} } }{c_{a}^{w} +c_{v}^{w} +c_{g}^{w} }

where

.. math::
   :label: 5.108

   c_{a}^{w} =\frac{1}{r_{aw} }

.. math::
   :label: 5.109

   c_{v}^{w} =\frac{\left(L+S\right)}{r_{b} } r''

.. math::
   :label: 5.110

   c_{g}^{w} =\frac{1}{r_{aw} ^{{'} } +r_{soil} }

are the water vapor conductances from the canopy air to the atmosphere, the leaf to canopy air, and ground to canopy air, respectively. The term :math:`r''` is determined from contributions by wet leaves and transpiration and limited by available water and potential evaporation as

.. math::
   :label: 5.111

   r'' = \left\{
   \begin{array}{lr}
   \min \left(f_{wet} +r_{dry} ^{{'} {'} } ,\, \frac{E_{v}^{w,\, pot} r_{dry} ^{{'} {'} } +\frac{W_{can} }{\Delta t} }{E_{v}^{w,\, pot} } \right) & \qquad E_{v}^{w,\, pot} >0,\, \beta _{t} >0 \\
   \min \left(f_{wet} ,\, \frac{E_{v}^{w,\, pot} r_{dry} ^{{'} {'} } +\frac{W_{can} }{\Delta t} }{E_{v}^{w,\, pot} } \right) & \qquad E_{v}^{w,\, pot} >0,\, \beta _{t} \le 0 \\
   1 & \qquad E_{v}^{w,\, pot} \le 0
   \end{array}\right\}

where :math:`f_{wet}` is the fraction of leaves and stems that are wet (section :numref:`Canopy Water`), :math:`W_{can}` is canopy water (kg m\ :sup:`-2`) (section :numref:`Canopy Water`), :math:`\Delta t` is the time step (s), and :math:`\beta _{t}` is a soil moisture function limiting transpiration (Chapter :numref:`rst_Stomatal Resistance and Photosynthesis`). The potential evaporation from wet foliage per unit wetted area is

.. math::
   :label: 5.112

   E_{v}^{w,\, pot} =-\frac{\rho _{atm} \left(q_{s} -q_{sat}^{T_{v} } \right)}{r_{b} } .

The term :math:`r_{dry} ^{{'} {'} }`  is

.. math::
   :label: 5.113

   r_{dry} ^{{'} {'} } =\frac{f_{dry} r_{b} }{L} \left(\frac{L^{sun} }{r_{b} +r_{s}^{sun} } +\frac{L^{sha} }{r_{b} +r_{s}^{sha} } \right)

where :math:`f_{dry}` is the fraction of leaves that are dry (section :numref:`Canopy Water`), :math:`L^{sun}` and :math:`L^{sha}` are the sunlit and shaded leaf area indices (section :numref:`Solar Fluxes`), and :math:`r_{s}^{sun}` and :math:`r_{s}^{sha}` are the sunlit and shaded stomatal resistances (s m\ :sup:`-1`) (Chapter :numref:`rst_Stomatal Resistance and Photosynthesis`).

When the expression for :math:`q_{s}` is substituted into equation :eq:`5.101`, the water vapor flux from vegetation :math:`E_{v}` is a function of :math:`q_{atm}`, :math:`q_{g}`, and :math:`q_{sat}^{T_{v} }`

.. math::
   :label: 5.114

   E_{v} =-\rho _{atm} \left[c_{a}^{w} q_{atm} +c_{g}^{w} q_{g} -\left(c_{a}^{w} +c_{g}^{w} \right)q_{sat}^{T_{v} } \right]\frac{c_{v}^{w} }{c_{a}^{w} +c_{v}^{w} +c_{g}^{w} } .

Similarly, the expression for :math:`q_{s}` can be substituted into :eq:`5.84` to obtain the water vapor flux from the ground beneath the canopy :math:`E_{g}`

.. math::
   :label: 5.115

   E_{g} =-\rho _{atm} \left[c_{a}^{w} q_{atm} +c_{v}^{w} q_{sat}^{T_{v} } -\left(c_{a}^{w} +c_{v}^{w} \right)q_{g} \right]\frac{c_{g}^{w} }{c_{a}^{w} +c_{v}^{w} +c_{g}^{w} } .

The aerodynamic resistances to heat (moisture) transfer between the ground at height :math:`z_{0h} ^{{'} }` (:math:`z_{0w} ^{{'} }` ) and the canopy air at height :math:`z_{0h} +d` (:math:`z_{0w} +d`) are

.. math::
   :label: 5.116

   r_{ah} ^{{'} } =r_{aw} ^{{'} } =\frac{1}{C_{s} U_{av} }

where

.. math::
   :label: 5.117

   U_{av} =V_{a} \sqrt{\frac{1}{r_{am} V_{a} } } =u_{*}

is the magnitude of the wind velocity incident on the leaves (equivalent here to friction velocity) (m s\ :sup:`-1`) and :math:`C_{s}` is the turbulent transfer coefficient between the underlying soil and the canopy air. :math:`C_{s}` is obtained by interpolation between values for dense canopy and bare soil (:ref:`Zeng et al. 2005 <Zengetal2005>`)

.. math::
   :label: 5.118

   C_{s} =C_{s,\, bare} W+C_{s,\, dense} (1-W)

where the weight :math:`W` is

.. math::
   :label: 5.119

   W=e^{-\left(L+S\right)} .

The dense canopy turbulent transfer coefficient (:ref:`Dickinson et al. 1993 <Dickinsonetal1993>`) is

.. math::
   :label: 5.120)

   C_{s,\, dense} =0.004 \ .

The bare soil turbulent transfer coefficient is

.. math::
   :label: 5.121

   C_{s,\, bare} =\frac{k}{a} \left(\frac{z_{0m,\, g} U_{av} }{\upsilon } \right)^{-0.45}

where the kinematic viscosity of air :math:`\upsilon =1.5\times 10^{-5}` m\ :sup:`2` s\ :sup:`-1` and :math:`a=0.13`.

The leaf boundary layer resistance :math:`r_{b}`  is

.. math::
   :label: 5.122

   r_{b} =\frac{1}{C_{v} } \left({U_{av} \mathord{\left/ {\vphantom {U_{av}  d_{leaf} }} \right.} d_{leaf} } \right)^{{-1\mathord{\left/ {\vphantom {-1 2}} \right.} 2} }

where :math:`C_{v} =0.01` m\ s\ :sup:`-1/2` is the turbulent transfer coefficient between the canopy surface and canopy air, and :math:`d_{leaf}` is the characteristic dimension of the leaves in the direction of wind flow (:numref:`Table Plant functional type aerodynamic parameters`).

The partial derivatives of the fluxes from the soil beneath the canopy with respect to ground temperature, which are needed for the soil temperature calculations (section :numref:`Numerical Solution Temperature`) and to update the soil surface fluxes (section :numref:`Update of Ground Sensible and Latent Heat Fluxes`), are

.. math::
   :label: 5.123

   \frac{\partial H_{g} }{\partial T_{g} } = \frac{\rho _{atm} C_{p} }{r'_{ah} } \frac{c_{a}^{h} +c_{v}^{h} }{c_{a}^{h} +c_{v}^{h} +c_{g}^{h} }

.. math::
   :label: 5.124

   \frac{\partial E_{g} }{\partial T_{g} } = \frac{\rho _{atm} }{r'_{aw} +r_{soil} } \frac{c_{a}^{w} +c_{v}^{w} }{c_{a}^{w} +c_{v}^{w} +c_{g}^{w} } \frac{dq_{g} }{dT_{g} } .

The partial derivatives :math:`\frac{\partial r'_{ah} }{\partial T_{g} }` and :math:`\frac{\partial r'_{aw} }{\partial T_{g} }`, which cannot be determined analytically, are ignored for :math:`\frac{\partial H_{g} }{\partial T_{g} }` and :math:`\frac{\partial E_{g} }{\partial T_{g} }`.

The roughness lengths used to calculate :math:`r_{am}`, :math:`r_{ah}`, and :math:`r_{aw}` from :eq:`5.55`, :eq:`5.56`, and :eq:`5.57` are :math:`z_{0m} =z_{0m,\, v}`, :math:`z_{0h} =z_{0h,\, v}`, and :math:`z_{0w} =z_{0w,\, v}`. 

The vegetation roughness lengths and displacement height :math:`d` are from :ref:`Meier et al. (2022) <Meieretal2022>`

.. math::
   :label: 5.125

   z_{0m,\, v} = z_{0h,\, v} =z_{0w,\, v} = z_{top} (1 - \frac{d} {z_{top} } ) \exp (\psi_{h} - \frac{k U_{h}} {u_{*} } )

where :math:`z_{top}` is canopy top height (m) (:numref:`Table Plant functional type canopy top and bottom heights`), :math:`k` is the von Karman constant (:numref:`Table Physical constants`), and :math:`\psi_{h}` is the roughness sublayer influence function

.. math::
   :label: 5.125a

   \psi_{h} = \ln(c_{w}) - 1 + c_{w}^{-1}

where :math:`c_{w}` is a pft-dependent constant (:numref:`Table Plant functional type aerodynamic parameters`).

The ratio of wind speed at canopy height to friction velocity, :math:`\frac{U_{h}} {u_{*}}` is derived from an implicit function of the roughness density :math:`\lambda`

.. math::
   :label: 5.125b

   \frac{U_{h}} {u_{*} } =(C_{S} + \lambda C_{R})^{0.5} \exp(\frac{\min \left(\lambda, \lambda_{\max}\right) c U_{h}} {2 u_{*}})

where :math:`C_{S}` represents the drag coefficient of the ground in the absence of vegetation, :math:`C_{R}` is the drag coefficient of an isolated roughness element (plant), and :math:`c` is an empirical constant. These three are pft-dependent parameters (:numref:`Table Plant functional type aerodynamic parameters`). :math:`\lambda_{max}` is the maximum :math:`\lambda` above which :math:`\frac{U_{h}} {u_{*}}` becomes constant. :math:`\lambda_{max}` is set to the value of :math:`\lambda` for which :eq:`5.125b`, in the absence of :math:`\lambda_{max}`, would have its minimum. :math:`\lambda_{max}` is also a pft-dependent parameter (:numref:`Table Plant functional type aerodynamic parameters`). :eq:`5.125b` can be written as

.. math::
   :label: 5.125c

   X \exp(-X) =(C_{S} + \lambda C_{R})^{0.5} c \frac{\lambda} {2 }

where 

.. math::
   :label: 5.125d

   X =\frac{c \lambda U_{h}} {2 u_{*} }.

:math:`X` and therefore :math:`\frac{U_{h}} {u_{*}}` can be solved for iteratively where the initial value of :math:`X` is 

.. math::
   :label: 5.125e

   X_{i=0} =(C_{S} + \lambda C_{R})^{0.5} c \frac{\lambda} {2 }

and the next value of :math:`X` at :math:`i+1` is

.. math::
   :label: 5.125f

   X_{i+1} =(C_{S} + \lambda C_{R})^{0.5} c \frac{\lambda} {2 } \exp(X_{i}).

:math:`X` is updated until :math:`\frac{U_{h}} {u_{*}}` converges to within :math:`1 \times 10^{-4}` between iterations.

:math:`\lambda` is set to half the total single-sided area of all canopy elements, here defined as the vegetation area index (VAI) defined as the sum of leaf (:math:`L`) and stem area index (:math:`S`), subject to a maximum of :math:`\lambda_{max}` and a minimum limit applied for numerical stability

.. math::
   :label: 5.126
   
   \lambda = \frac{\min(\max(1 \times 10^{-5}, VAI), \lambda_{max})} {2 }

The displacement height :math:`d` is

.. math::
   :label: 5.127

   d = z_{top}\left[1- \frac{1-\exp(-(c_{d1} 2 \lambda)^{0.5}} {(c_{d1} 2 \lambda)^{0.5} }\right]

where :math:`c_{d1} =7.5`.

.. _Table Plant functional type aerodynamic parameters:

.. table:: Plant functional type aerodynamic parameters

 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | Plant functional type            | :math:`d_{leaf}`  (m) | :math:`c_{w}`           | :math:`C_{S}`           | :math:`C_{R}`           | :math:`c`               | :math:`\lambda_{max}`   |
 +==================================+=======================+=========================+=========================+=========================+=========================+=========================+
 | NET Temperate                    | 0.04                  | 9                       | 0.003                   | 0.05                    | 0.09                    | 4.55                    |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | NET Boreal                       | 0.04                  | 9                       | 0.003                   | 0.05                    | 0.09                    | 4.55                    |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | NDT Boreal                       | 0.04                  | 9                       | 0.003                   | 0.05                    | 0.09                    | 4.55                    |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | BET Tropical                     | 0.04                  | 3                       | 0.01                    | 0.14                    | 0.01                    | 7.87                    |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | BET temperate                    | 0.04                  | 3                       | 0.01                    | 0.14                    | 0.01                    | 7.87                    |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | BDT tropical                     | 0.04                  | 1                       | 0.013                   | 0.13                    | 0.06                    | 8.88                    |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | BDT temperate                    | 0.04                  | 1                       | 0.013                   | 0.13                    | 0.06                    | 8.88                    |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | BDT boreal                       | 0.04                  | 1                       | 0.013                   | 0.13                    | 0.06                    | 8.88                    |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | BES temperate                    | 0.04                  | 20                      | 0.001                   | 0.05                    | 0.12                    | 3.07                    |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | BDS temperate                    | 0.04                  | 20                      | 0.001                   | 0.05                    | 0.12                    | 3.07                    |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | BDS boreal                       | 0.04                  | 20                      | 0.001                   | 0.05                    | 0.12                    | 3.07                    |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | C\ :sub:`3` arctic grass         | 0.04                  | 19                      | 0.001                   | 0.05                    | 0.08                    | 4.61                    |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | C\ :sub:`3` grass                | 0.04                  | 19                      | 0.001                   | 0.05                    | 0.08                    | 4.61                    |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | C\ :sub:`4` grass                | 0.04                  | 19                      | 0.001                   | 0.05                    | 0.08                    | 4.61                    |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | Crop R                           | 0.04                  | 3.5                     | 0.001                   | 0.05                    | 0.04                    | 5.3                     |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | Crop I                           | 0.04                  | 3.5                     | 0.001                   | 0.05                    | 0.04                    | 5.3                     |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | Corn R                           | 0.04                  | 3.5                     | 0.001                   | 0.05                    | 0.04                    | 5.3                     |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | Corn I                           | 0.04                  | 3.5                     | 0.001                   | 0.05                    | 0.04                    | 5.3                     |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | Temp Cereal R                    | 0.04                  | 3.5                     | 0.001                   | 0.05                    | 0.04                    | 5.3                     |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | Temp Cereal I                    | 0.04                  | 3.5                     | 0.001                   | 0.05                    | 0.04                    | 5.3                     |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | Winter Cereal R                  | 0.04                  | 3.5                     | 0.001                   | 0.05                    | 0.04                    | 5.3                     |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | Winter Cereal I                  | 0.04                  | 3.5                     | 0.001                   | 0.05                    | 0.04                    | 5.3                     |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | Soybean R                        | 0.04                  | 3.5                     | 0.001                   | 0.05                    | 0.04                    | 5.3                     |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | Soybean I                        | 0.04                  | 3.5                     | 0.001                   | 0.05                    | 0.04                    | 5.3                     |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | Miscanthus R                     | 0.04                  | 3.5                     | 0.001                   | 0.05                    | 0.04                    | 5.3                     |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | Miscanthus I                     | 0.04                  | 3.5                     | 0.001                   | 0.05                    | 0.04                    | 5.3                     |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | Switchgrass R                    | 0.04                  | 3.5                     | 0.001                   | 0.05                    | 0.04                    | 5.3                     |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 | Switchgrass I                    | 0.04                  | 3.5                     | 0.001                   | 0.05                    | 0.04                    | 5.3                     |
 +----------------------------------+-----------------------+-------------------------+-------------------------+-------------------------+-------------------------+-------------------------+
 
.. _Numerical Implementation:

Numerical Implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Canopy energy conservation gives

.. math::
   :label: 5.128

   -\overrightarrow{S}_{v} +\overrightarrow{L}_{v} \left(T_{v} \right)+H_{v} \left(T_{v} \right)+\lambda E_{v} \left(T_{v} \right)=0

where :math:`\overrightarrow{S}_{v}` is the solar radiation absorbed by the vegetation (section :numref:`Solar Fluxes`), :math:`\overrightarrow{L}_{v}` is the net longwave radiation absorbed by vegetation (section :numref:`Longwave Fluxes`), and :math:`H_{v}` and :math:`\lambda E_{v}` are the sensible and latent heat fluxes from vegetation, respectively. The term :math:`\lambda` is taken to be the latent heat of vaporization :math:`\lambda _{vap}` (:numref:`Table Physical constants`).

:math:`\overrightarrow{L}_{v}`, :math:`H_{v}`, and :math:`\lambda E_{v}` depend on the vegetation temperature :math:`T_{v}`. The Newton-Raphson method for finding roots of non-linear systems of equations can be applied to iteratively solve for :math:`T_{v}` as

.. math::
   :label: 5.129

   \Delta T_{v} =\frac{\overrightarrow{S}_{v} -\overrightarrow{L}_{v} -H_{v} -\lambda E_{v} }{\frac{\partial \overrightarrow{L}_{v} }{\partial T_{v} } +\frac{\partial H_{v} }{\partial T_{v} } +\frac{\partial \lambda E_{v} }{\partial T_{v} } }

where :math:`\Delta T_{v} =T_{v}^{n+1} -T_{v}^{n}` and the subscript "n" indicates the iteration.

The partial derivatives are

.. math::
   :label: 5.130

   \frac{\partial \overrightarrow{L}_{v} }{\partial T_{v} } =4\varepsilon _{v} \sigma \left[2-\varepsilon _{v} \left(1-\varepsilon _{g} \right)\right]T_{v}^{3}

.. math::
   :label: 5.131

   \frac{\partial H_{v} }{\partial T_{v} } =\rho _{atm} C_{p} \left(c_{a}^{h} +c_{g}^{h} \right)\frac{c_{v}^{h} }{c_{a}^{h} +c_{v}^{h} +c_{g}^{h} }

.. math::
   :label: 5.132

   \frac{\partial \lambda E_{v} }{\partial T_{v} } =\lambda \rho _{atm} \left(c_{a}^{w} +c_{g}^{w} \right)\frac{c_{v}^{w} }{c_{a}^{w} +c_{v}^{w} +c_{g}^{w} } \frac{dq_{sat}^{T_{v} } }{dT_{v} } .

The partial derivatives :math:`\frac{\partial r_{ah} }{\partial T_{v} }` and :math:`\frac{\partial r_{aw} }{\partial T_{v} }`, which cannot be determined analytically, are ignored for :math:`\frac{\partial H_{v} }{\partial T_{v} }` and :math:`\frac{\partial \lambda E_{v} }{\partial T_{v} }`. However, if :math:`\zeta` changes sign more than four times during the temperature iteration, :math:`\zeta =-0.01`. This helps prevent "flip-flopping" between stable and unstable conditions. The total water vapor flux :math:`E_{v}`, transpiration flux :math:`E_{v}^{t}`, and sensible heat flux :math:`H_{v}` are updated for changes in leaf temperature as

.. math::
   :label: 5.133

   E_{v} =-\rho _{atm} \left[c_{a}^{w} q_{atm} +c_{g}^{w} q_{g} -\left(c_{a}^{w} +c_{g}^{w} \right)\left(q_{sat}^{T_{v} } +\frac{dq_{sat}^{T_{v} } }{dT_{v} } \Delta T_{v} \right)\right]\frac{c_{v}^{w} }{c_{a}^{w} +c_{v}^{w} +c_{g}^{w} }

.. math::
   :label: 5.134

   E_{v}^{t} =-r_{dry} ^{{'} {'} } \rho _{atm} \left[c_{a}^{w} q_{atm} +c_{g}^{w} q_{g} -\left(c_{a}^{w} +c_{g}^{w} \right)\left(q_{sat}^{T_{v} } +\frac{dq_{sat}^{T_{v} } }{dT_{v} } \Delta T_{v} \right)\right]\frac{c_{v}^{h} }{c_{a}^{w} +c_{v}^{w} +c_{g}^{w} }

.. math::
   :label: 5.135

   H_{v} =-\rho _{atm} C_{p} \left[c_{a}^{h} \theta _{atm} +c_{g}^{h} T_{g} -\left(c_{a}^{h} +c_{g}^{h} \right)\left(T_{v} +\Delta T_{v} \right)\right]\frac{c_{v}^{h} }{c_{a}^{h} +c_{v}^{h} +c_{g}^{h} } .

The numerical solution for vegetation temperature and the fluxes of momentum, sensible heat, and water vapor flux from vegetated surfaces proceeds as follows:

#. Initial values for canopy air temperature and specific humidity are obtained from

   .. math::
      :label: 5.136

      T_{s} =\frac{T_{g} +\theta _{atm} }{2}

   .. math::
      :label: 5.137

      q_{s} =\frac{q_{g} +q_{atm} }{2} .

#. An initial guess for the wind speed :math:`V_{a}` is obtained from :eq:`5.24` assuming an initial convective velocity :math:`U_{c} =0` m s\ :sup:`-1` for stable conditions (:math:`\theta _{v,\, atm} -\theta _{v,\, s} \ge 0` as evaluated from :eq:`5.50` ) and :math:`U_{c} =0.5` for unstable conditions (:math:`\theta _{v,\, atm} -\theta _{v,\, s} <0`).

#. An initial guess for the Monin-Obukhov length :math:`L` is obtained from the bulk Richardson number using equations :eq:`5.46` and :eq:`5.48`.

#. Iteration proceeds on the following system of equations:

#. Friction velocity :math:`u_{*}` (:eq:`5.32`, :eq:`5.33`, :eq:`5.34`, :eq:`5.35`)

#. Ratio :math:`\frac{\theta _{*} }{\theta _{atm} -\theta _{s} }` (:eq:`5.37`, :eq:`5.38`, :eq:`5.39`, :eq:`5.40`)

#. Ratio :math:`\frac{q_{*} }{q_{atm} -q_{s} }` (:eq:`5.41`, :eq:`5.42`, :eq:`5.43`, :eq:`5.44`)

#. Aerodynamic resistances :math:`r_{am}`, :math:`r_{ah}`, and :math:`r_{aw}` (:eq:`5.55`, :eq:`5.56`, :eq:`5.57`)

#. Magnitude of the wind velocity incident on the leaves :math:`U_{av}` (:eq:`5.117` )

#. Leaf boundary layer resistance :math:`r_{b}` (:eq:`5.122` )

#. Aerodynamic resistances :math:`r_{ah} ^{{'} }` and :math:`r_{aw} ^{{'} }`(:eq:`5.116` )

#. Sunlit and shaded stomatal resistances :math:`r_{s}^{sun}` and :math:`r_{s}^{sha}` (Chapter :numref:`rst_Stomatal Resistance and Photosynthesis`)

#. Sensible heat conductances :math:`c_{a}^{h}`, :math:`c_{g}^{h}`, and :math:`c_{v}^{h}` (:eq:`5.94`, :eq:`5.95`, :eq:`5.96`)

#. Latent heat conductances :math:`c_{a}^{w}`, :math:`c_{v}^{w}`, and :math:`c_{g}^{w}` (:eq:`5.108`, :eq:`5.109`, :eq:`5.110`)

#. Sensible heat flux from vegetation :math:`H_{v}` (:eq:`5.97` )

#. Latent heat flux from vegetation :math:`\lambda E_{v}` (:eq:`5.101` )

#. If the latent heat flux has changed sign from the latent heat flux computed at the previous iteration (:math:`\lambda E_{v} ^{n+1} \times \lambda E_{v} ^{n} <0`), the latent heat flux is constrained to be 10% of the computed value. The difference between the constrained and computed value (:math:`\Delta _{1} =0.1\lambda E_{v} ^{n+1} -\lambda E_{v} ^{n+1}` ) is added to the sensible heat flux later.

#. Change in vegetation temperature :math:`\Delta T_{v}` (:eq:`5.129` ) and update the vegetation temperature as :math:`T_{v}^{n+1} =T_{v}^{n} +\Delta T_{v}`. :math:`T_{v}` is constrained to change by no more than 1°K in one iteration. If this limit is exceeded, the energy error is

   .. math::
      :label: 5.138

      \Delta _{2} =\overrightarrow{S}_{v} -\overrightarrow{L}_{v} -\frac{\partial \overrightarrow{L}_{v} }{\partial T_{v} } \Delta T_{v} -H_{v} -\frac{\partial H_{v} }{\partial T_{v} } \Delta T_{v} -\lambda E_{v} -\frac{\partial \lambda E_{v} }{\partial T_{v} } \Delta T_{v}

where :math:`\Delta T_{v} =1{\rm \; or\; }-1`. The error :math:`\Delta _{2}` is added to the sensible heat flux later.

#. Water vapor flux :math:`E_{v}` (:eq:`5.133` )

#. Transpiration :math:`E_{v}^{t}` (:eq:`5.134` if :math:`\beta_{t} >0`, otherwise :math:`E_{v}^{t} =0`)

#. The water vapor flux :math:`E_{v}` is constrained to be less than or equal to the sum of transpiration :math:`E_{v}^{t}` and the water available from wetted leaves and stems :math:`{W_{can} \mathord{\left/ {\vphantom {W_{can} \Delta t}} \right.} \Delta t}`. The energy error due to this constraint is

   .. math::
      :label: 5.139

      \Delta _{3} =\max \left(0,\, E_{v} -E_{v}^{t} -\frac{W_{can} }{\Delta t} \right).

The error :math:`\lambda \Delta _{3}` is added to the sensible heat flux later.

#. Sensible heat flux :math:`H_{v}` (:eq:`5.135` ). The three energy error terms, :math:`\Delta _{1}`, :math:`\Delta _{2}`, and :math:`\lambda \Delta _{3}` are also added to the sensible heat flux.

#. The saturated vapor pressure :math:`e_{i}` (Chapter :numref:`rst_Stomatal Resistance and Photosynthesis`), saturated specific humidity :math:`q_{sat}^{T_{v} }` and its derivative :math:`\frac{dq_{sat}^{T_{v} } }{dT_{v} }` at the leaf surface (section :numref:`Saturation Vapor Pressure`), are re-evaluated based on the new :math:`T_{v}`.

#. Canopy air temperature :math:`T_{s}` (:eq:`5.93` )

#. Canopy air specific humidity :math:`q_{s}` (:eq:`5.107` )

#. Temperature difference :math:`\theta _{atm} -\theta _{s}`

#. Specific humidity difference :math:`q_{atm} -q_{s}`

#. Potential temperature scale :math:`\theta _{*} =\frac{\theta _{*} }{\theta _{atm} -\theta _{s} } \left(\theta _{atm} -\theta _{s} \right)` where :math:`\frac{\theta _{*} }{\theta _{atm} -\theta _{s} }` was calculated earlier in the iteration #. Humidity scale :math:`q_{*} =\frac{q_{*} }{q_{atm} -q_{s} } \left(q_{atm} -q_{s} \right)` where :math:`\frac{q_{*} }{q_{atm} -q_{s} }` was calculated earlier in the iteration #. Virtual potential temperature scale :math:`\theta _{v*}` (:eq:`5.17` )

#. Wind speed including the convective velocity, :math:`V_{a}` (:eq:`5.24` )

#. Monin-Obukhov length :math:`L` (:eq:`5.49` )

#. The iteration is stopped after two or more steps if :math:`\tilde{\Delta }T_{v} <0.01` and :math:`\left|\lambda E_{v}^{n+1} -\lambda E_{v}^{n} \right|<0.1` where :math:`\tilde{\Delta }T_{v} =\max \left(\left|T_{v}^{n+1} -T_{v}^{n} \right|,\, \left|T_{v}^{n} -T_{v}^{n-1} \right|\right)`, or after forty iterations have been carried out.

#. Momentum fluxes :math:`\tau _{x}`, :math:`\tau _{y}` (:eq:`5.5`, :eq:`5.6`)

#. Sensible heat flux from ground :math:`H_{g}` (:eq:`5.89` )

#. Water vapor flux from ground :math:`E_{g}` (:eq:`5.102` )

#. 2-m height air temperature :math:`T_{2m}`, specific humidity :math:`q_{2m}`, relative humidity :math:`RH_{2m}` \ (:eq:`5.58`, :eq:`5.59`, :eq:`5.60`)

.. _Update of Ground Sensible and Latent Heat Fluxes:

Update of Ground Sensible and Latent Heat Fluxes
----------------------------------------------------

The sensible and water vapor heat fluxes derived above for bare soil and soil beneath canopy are based on the ground surface temperature from the previous time step :math:`T_{g}^{n}` and are used as the surface forcing for the solution of the soil temperature equations (section :numref:`Numerical Solution Temperature`). This solution yields a new ground surface temperature :math:`T_{g}^{n+1}`. The ground sensible and water vapor fluxes are then updated for :math:`T_{g}^{n+1}` as

.. math::
   :label: 5.140

   H'_{g} =H_{g} +\left(T_{g}^{n+1} -T_{g}^{n} \right)\frac{\partial H_{g} }{\partial T_{g} }

.. math::
   :label: 5.141

   E'_{g} =E_{g} +\left(T_{g}^{n+1} -T_{g}^{n} \right)\frac{\partial E_{g} }{\partial T_{g} }

where :math:`H_{g}`, :math:`E_{g}`, :math:`\frac{\partial H_{g} }{\partial T_{g} }`, and :math:`\frac{\partial E_{g} }{\partial T_{g} }` are the sensible heat and water vapor fluxes and their partial derivatives derived from equations :eq:`5.62`, :eq:`5.66`, :eq:`5.83`, and :eq:`5.84` for non-vegetated surfaces and equations :eq:`5.89`, :eq:`5.102`, :eq:`5.123`, and :eq:`5.124` for vegetated surfaces using :math:`T_{g}^{n}`. One further adjustment is made to :math:`H'_{g}` and :math:`E'_{g}`. If the soil moisture in the top snow/soil layer is not sufficient to support the updated ground evaporation, i.e., if :math:`E'_{g} > 0` and :math:`f_{evap} < 1` where

.. math::
   :label: 5.142

   f_{evap} =\frac{{\left(w_{ice,\; snl+1} +w_{liq,\, snl+1} \right)\mathord{\left/ {\vphantom {\left(w_{ice,\; snl+1} +w_{liq,\, snl+1} \right) \Delta t}} \right.} \Delta t} }{\sum _{j=1}^{npft}\left(E'_{g} \right)_{j} \left(wt\right)_{j}  } \le 1,

an adjustment is made to reduce the ground evaporation accordingly as

.. math::
   :label: 5.143

   E''_{g} =f_{evap} E'_{g} .

The term :math:`\sum _{j=1}^{npft}\left(E'_{g} \right)_{j} \left(wt\right)_{j}` is the sum of :math:`E'_{g}` over all evaporating PFTs where :math:`\left(E'_{g} \right)_{j}` is the ground evaporation from the :math:`j^{th}` PFT on the column, :math:`\left(wt\right)_{j}` is the relative area of the :math:`j^{th}` PFT with respect to the column, and :math:`npft` is the number of PFTs on the column. :math:`w_{ice,\, snl+1}` and :math:`w_{liq,\, snl+1}` are the ice and liquid water contents (kg m\ :sup:`-2`) of the top snow/soil layer (Chapter :numref:`rst_Hydrology`). Any resulting energy deficit is assigned to sensible heat as

.. math::
   :label: 5.144

   H''_{g} =H_{g} +\lambda \left(E'_{g} -E''_{g} \right).

The ground water vapor flux :math:`E''_{g}` is partitioned into evaporation of liquid water from snow/soil :math:`q_{seva}` (kg\ m\ :sup:`-2` s\ :sup:`-1`), sublimation from snow/soil ice :math:`q_{subl}` (kg m\ :sup:`-2` s\ :sup:`-1`), liquid dew on snow/soil :math:`q_{sdew}` (kg m\ :sup:`-2` s\ :sup:`-1`), or frost on snow/soil :math:`q_{frost}` (kg m\ :sup:`-2` s\ :sup:`-1`) as

.. math::
   :label: 5.145

   q_{seva} =\max \left(E''_{sno} \frac{w_{liq,\, snl+1} }{w_{ice,\; snl+1} +w_{liq,\, snl+1} } ,0\right)\qquad E''_{sno} \ge 0,\, w_{ice,\; snl+1} +w_{liq,\, snl+1} >0

.. math::
   :label: 5.146

   q_{subl} =E''_{sno} -q_{seva} \qquad E''_{sno} \ge 0

.. math::
   :label: 5.147

   q_{sdew} =\left|E''_{sno} \right|\qquad E''_{sno} <0{\rm \; and\; }T_{g} \ge T_{f}

.. math::
   :label: 5.148

   q_{frost} =\left|E''_{sno} \right|\qquad E''_{sno} <0{\rm \; and\; }T_{g} <T_{f} .

The loss or gain in snow mass due to :math:`q_{seva}`, :math:`q_{subl}`, :math:`q_{sdew}`, and :math:`q_{frost}` on a snow surface are accounted for during the snow hydrology calculations (Chapter :numref:`rst_Snow Hydrology`). The loss of soil and surface water due to :math:`q_{seva}` is accounted for in the calculation of infiltration (section :numref:`Infiltration`), while losses or gains due to :math:`q_{subl}`, :math:`q_{sdew}`, and :math:`q_{frost}` on a soil surface are accounted for following the sub-surface drainage calculations (section :numref:`Lateral Sub-surface Runoff`).

The ground heat flux :math:`G` is calculated as

.. math::
   :label: 5.149

   G=\overrightarrow{S}_{g} -\overrightarrow{L}_{g} -H_{g} -\lambda E_{g}

where :math:`\overrightarrow{S}_{g}` is the solar radiation absorbed by the ground (section :numref:`Solar Fluxes`), :math:`\overrightarrow{L}_{g}` is the net longwave radiation absorbed by the ground (section :numref:`Longwave Fluxes`)

.. math::
   :label: 5.150

   \vec{L}_{g} =L_{g} \uparrow -\delta _{veg} \varepsilon _{g} L_{v} \, \downarrow -\left(1-\delta _{veg} \right)\varepsilon _{g} L_{atm} \, \downarrow +4\varepsilon _{g} \sigma \left(T_{g}^{n} \right)^{3} \left(T_{g}^{n+1} -T_{g}^{n} \right),

where

.. math::
   :label: 5.151

   L_{g} \uparrow =\varepsilon _{g} \sigma \left[\left(1-f_{sno} -f_{h2osfc} \right)\left(T_{1}^{n} \right)^{4} +f_{sno} \left(T_{sno}^{n} \right)^{4} +f_{h2osfc} \left(T_{h2osfc}^{n} \right)^{4} \right]

and :math:`H_{g}` and :math:`\lambda E_{g}` are the sensible and latent heat fluxes after the adjustments described above.

When converting ground water vapor flux to an energy flux, the term :math:`\lambda` is arbitrarily assumed to be

.. math::
   :label: 5.152

   \lambda =\left\{\begin{array}{l} {\lambda _{sub} \qquad {\rm if\; }w_{liq,\, snl+1} =0{\rm \; and\; }w_{ice,\, snl+1} >0} \\ {\lambda _{vap} \qquad {\rm otherwise}} \end{array}\right\}

where :math:`\lambda _{sub}` and :math:`\lambda _{vap}` are the latent heat of sublimation and vaporization, respectively (J (kg\ :sup:`-1`) (:numref:`Table Physical constants`). When converting vegetation water vapor flux to an energy flux, :math:`\lambda _{vap}` is used.

The system balances energy as

.. math::
   :label: 5.153

   \overrightarrow{S}_{g} +\overrightarrow{S}_{v} +L_{atm} \, \downarrow -L\, \uparrow -H_{v} -H_{g} -\lambda _{vap} E_{v} -\lambda E_{g} -G=0.

.. _Saturation Vapor Pressure:

Saturation Vapor Pressure
-----------------------------

Saturation vapor pressure :math:`e_{sat}^{T}` (Pa) and its derivative :math:`\frac{de_{sat}^{T} }{dT}`, as a function of temperature :math:`T` (°C), are calculated from the eighth-order polynomial fits of :ref:`Flatau et al. (1992) <Flatauetal1992>`

.. math::
   :label: 5.154

   e_{sat}^{T} =100\left[a_{0} +a_{1} T+\cdots +a_{n} T^{n} \right]

.. math::
   :label: 5.155

   \frac{de_{sat}^{T} }{dT} =100\left[b_{0} +b_{1} T+\cdots +b_{n} T^{n} \right]

where the coefficients for ice are valid for :math:`-75\, ^{\circ } {\rm C}\le T<0\, ^{\circ } {\rm C}` and the coefficients for water are valid for :math:`0\, ^{\circ } {\rm C}\le T\le 100\, ^{\circ } {\rm C}` (:numref:`Table Coefficients for saturation vapor pressure` and :numref:`Table Coefficients for derivative of esat`). The saturated water vapor specific humidity :math:`q_{sat}^{T}` and its derivative :math:`\frac{dq_{sat}^{T} }{dT}` are

.. math::
   :label: 5.156

   q_{sat}^{T} =\frac{0.622e_{sat}^{T} }{P_{atm} -0.378e_{sat}^{T} }

.. math::
   :label: 5.157

   \frac{dq_{sat}^{T} }{dT} =\frac{0.622P_{atm} }{\left(P_{atm} -0.378e_{sat}^{T} \right)^{2} } \frac{de_{sat}^{T} }{dT} .

.. _Table Coefficients for saturation vapor pressure:

.. table:: Coefficients for :math:`e_{sat}^{T}`

 +------------------+------------------------------------------+----------------------------------------+
 |                  | water                                    | ice                                    |
 +==================+==========================================+========================================+
 | :math:`a_{0}`    | 6.11213476                               | 6.11123516                             |
 +------------------+------------------------------------------+----------------------------------------+
 | :math:`a_{1}`    | 4.44007856 :math:`\times 10^{-1}`        | 5.03109514\ :math:`\times 10^{-1}`     |
 +------------------+-------------------------------------------+---------------------------------------+
 | :math:`a_{2}`    | 1.43064234 :math:`\times 10^{-2}`        | 1.88369801\ :math:`\times 10^{-2}`     |
 +------------------+-------------------------------------------+---------------------------------------+
 | :math:`a_{3}`    | 2.64461437 :math:`\times 10^{-4}`        | 4.20547422\ :math:`\times 10^{-4}`     |
 +------------------+-------------------------------------------+---------------------------------------+
 | :math:`a_{4}`    | 3.05903558 :math:`\times 10^{-6}`        | 6.14396778\ :math:`\times 10^{-6}`     |
 +------------------+-------------------------------------------+---------------------------------------+
 | :math:`a_{5}`    | 1.96237241 :math:`\times 10^{-8}`        | 6.02780717\ :math:`\times 10^{-8}`     |
 +------------------+-------------------------------------------+---------------------------------------+
 | :math:`a_{6}`    | 8.92344772 :math:`\times 10^{-11}`       | 3.87940929\ :math:`\times 10^{-10}`    |
 +------------------+-------------------------------------------+---------------------------------------+
 | :math:`a_{7}`    | -3.73208410 :math:`\times 10^{-13}`      | 1.49436277\ :math:`\times 10^{-12}`    |
 +------------------+-------------------------------------------+---------------------------------------+
 | :math:`a_{8}`    | 2.09339997 :math:`\times 10^{-16}`       | 2.62655803\ :math:`\times 10^{-15}`    |
 +------------------+------------------------------------------+----------------------------------------+

.. _Table Coefficients for derivative of esat:

.. table:: Coefficients for :math:`\frac{de_{sat}^{T} }{dT}`

 +------------------+----------------------------------------+----------------------------------------+
 |                  | water                                  | ice                                    |
 +==================+========================================+========================================+
 | :math:`b_{0}`    | 4.44017302\ :math:`\times 10^{-1}`     | 5.03277922\ :math:`\times 10^{-1}`     |
 +------------------+----------------------------------------+----------------------------------------+
 | :math:`b_{1}`    | 2.86064092\ :math:`\times 10^{-2}`     | 3.77289173\ :math:`\times 10^{-2}`     |
 +------------------+----------------------------------------+----------------------------------------+
 | :math:`b_{2}`    | 7.94683137\ :math:`\times 10^{-4}`     | 1.26801703\ :math:`\times 10^{-3}`     |
 +------------------+----------------------------------------+----------------------------------------+
 | :math:`b_{3}`    | 1.21211669\ :math:`\times 10^{-5}`     | 2.49468427\ :math:`\times 10^{-5}`     |
 +------------------+----------------------------------------+----------------------------------------+
 | :math:`b_{4}`    | 1.03354611\ :math:`\times 10^{-7}`     | 3.13703411\ :math:`\times 10^{-7}`     |
 +------------------+----------------------------------------+----------------------------------------+
 | :math:`b_{5}`    | 4.04125005\ :math:`\times 10^{-10}`    | 2.57180651\ :math:`\times 10^{-9}`     |
 +------------------+----------------------------------------+----------------------------------------+
 | :math:`b_{6}`    | -7.88037859 :math:`\times 10^{-13}`    | 1.33268878\ :math:`\times 10^{-11}`    |
 +------------------+----------------------------------------+----------------------------------------+
 | :math:`b_{7}`    | -1.14596802 :math:`\times 10^{-14}`    | 3.94116744\ :math:`\times 10^{-14}`    |
 +------------------+----------------------------------------+----------------------------------------+
 | :math:`b_{8}`    | 3.81294516\ :math:`\times 10^{-17}`    | 4.98070196\ :math:`\times 10^{-17}`    |
 +------------------+----------------------------------------+----------------------------------------+

