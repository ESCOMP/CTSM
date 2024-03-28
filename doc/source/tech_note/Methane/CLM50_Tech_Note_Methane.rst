.. _rst_Methane Model:

Methane Model
=================

The representation of processes in the methane biogeochemical model integrated in CLM [CLM4Me; (:ref:`Riley et al. 2011a<Rileyetal2011a>`)] is based on several previously published models (:ref:`Cao et al. 1996<Caoetal1996>`; :ref:`Petrescu et al. 2010<Petrescuetal2010>`; :ref:`Tianet al. 2010<Tianetal2010>`; :ref:`Walter et al. 2001<Walteretal2001>`; :ref:`Wania et al. 2010<Waniaetal2010>`; :ref:`Zhang et al. 2002<Zhangetal2002>`; :ref:`Zhuang et al. 2004<Zhuangetal2004>`). Although the model has similarities with these precursor models, a number of new process representations and parameterization have been integrated into CLM.

Mechanistically modeling net surface CH\ :sub:`4` emissions requires representing a complex and interacting series of processes. We first (section :numref:`Methane Model Structure and Flow`) describe the overall model structure and flow of information in the CH\ :sub:`4` model, then describe the methods used to represent: CH\ :sub:`4` mass balance; CH\ :sub:`4` production; ebullition; aerenchyma transport; CH\ :sub:`4` oxidation; reactive transport solution, including boundary conditions, numerical solution, water table interface, etc.; seasonal inundation effects; and impact of seasonal inundation on CH\ :sub:`4` production.

.. _Methane Model Structure and Flow:

Methane Model Structure and Flow
-------------------------------------

The driver routine for the methane biogeochemistry calculations (ch4, in ch4Mod.F) controls the initialization of boundary conditions, inundation, and impact of redox conditions; calls to routines to calculate CH\ :sub:`4` production, oxidation, transport through aerenchyma, ebullition, and the overall mass balance (for unsaturated and saturated soils and, if desired, lakes); resolves changes to CH\ :sub:`4` calculations associated with a changing inundated fraction; performs a mass balance check; and calculates the average gridcell CH\ :sub:`4` production, oxidation, and exchanges with the atmosphere.

.. _Governing Mass-Balance Relationship:

Governing Mass-Balance Relationship
----------------------------------------

The model (:numref:`Figure Methane Schematic`) accounts for CH\ :sub:`4` production in the anaerobic fraction of soil (*P*, mol m\ :sup:`-3` s\ :sup:`-1`), ebullition (*E*, mol m\ :sup:`-3` s\ :sup:`-1`), aerenchyma transport (*A*, mol m\ :sup:`-3` s\ :sup:`-1`), aqueous and gaseous diffusion (:math:`{F}_{D}`, mol m\ :sup:`-2` s\ :sup:`-1`), and oxidation (*O*, mol m\ :sup:`-3` s\ :sup:`-1`) via a transient reaction diffusion equation:

.. math::
   :label: 24.1

   \frac{\partial \left(RC\right)}{\partial t} =\frac{\partial F_{D} }{\partial z} +P\left(z,t\right)-E\left(z,t\right)-A\left(z,t\right)-O\left(z,t\right)

Here *z* (m) represents the vertical dimension, *t* (s) is time, and *R* accounts for gas in both the aqueous and gaseous phases:\ :math:`R = \epsilon _{a} +K_{H} \epsilon _{w}`, with :math:`\epsilon _{a}`, :math:`\epsilon _{w}`, and :math:`K_{H}` (-) the air-filled porosity, water-filled porosity, and partitioning coefficient for the species of interest, respectively, and :math:`C` represents CH\ :sub:`4` or O\ :sub:`2` concentration with respect to water volume (mol m\ :sup:`-3`).

An analogous version of equation :eq:`24.1` is concurrently solved for O\ :sub:`2`, but with the following differences relative to CH\ :sub:`4`: *P* = *E* = 0 (i.e., no production or ebullition), and the oxidation sink includes the O\ :sub:`2` demanded by methanotrophs, heterotroph decomposers, nitrifiers, and autotrophic root respiration.

As currently implemented, each gridcell contains an inundated and a non-inundated fraction. Therefore, equation :eq:`24.1` is solved four times for each gridcell and time step: in the inundated and non-inundated fractions, and for CH\ :sub:`4` and O\ :sub:`2`. If desired, the CH\ :sub:`4` and O\ :sub:`2` mass balance equation is solved again for lakes (Chapter 9). For non-inundated areas, the water table interface is defined at the deepest transition from greater than 95% saturated to less than 95% saturated that occurs above frozen soil layers. The inundated fraction is allowed to change at each time step, and the total soil CH\ :sub:`4` quantity is conserved by evolving CH\ :sub:`4` to the atmosphere when the inundated fraction decreases, and averaging a portion of the non-inundated concentration into the inundated concentration when the inundated fraction increases.

.. _Figure Methane Schematic:

.. figure:: image1.png

 Schematic representation of biological and physical processes integrated in CLM that affect the net CH\ :sub:`4`
 surface flux (:ref:`Riley et al. 2011a<Rileyetal2011a>`). (left)
 Fully inundated portion of a CLM gridcell and (right) variably saturated portion of a gridcell.

.. _CH4 Production:

CH\ :sub:`4` Production
----------------------------------

Because CLM does not currently specifically represent wetland plant functional types or soil biogeochemical processes, we used gridcell-averaged decomposition rates as proxies. Thus, the upland (default) heterotrophic respiration is used to estimate the wetland decomposition rate after first dividing off the O\ :sub:`2` limitation. The O\ :sub:`2` consumption associated with anaerobic decomposition is then set to the unlimited version so that it will be reduced appropriately during O\ :sub:`2` competition. CH\ :sub:`4` production at each soil level in the anaerobic portion (i.e., below the water table) of the column is related to the gridcell estimate of heterotrophic respiration from soil and litter (R\ :sub:`H`; mol C m\ :sup:`-2` s\ :sub:`-1`) corrected for its soil temperature (:math:`{T}_{s}`) dependence, soil temperature through a :math:`{A}_{10}` factor (:math:`f_{T}`), pH (:math:`f_{pH}`), redox potential (:math:`f_{pE}`), and a factor accounting for the seasonal inundation fraction (*S*, described below):

.. math::
   :label: 24.2

   P=R_{H} f_{CH_{4} } f_{T} f_{pH} f_{pE} S.

Here, :math:`f_{CH_{4} }` is the baseline ratio between CO\ :sub:`2` and CH\ :sub:`4` production (all parameters values are given in :numref:`Table Methane Parameter descriptions`). Currently, :math:`f_{CH_{4} }` is modified to account for our assumptions that methanogens may have a higher Q\ :math:`{}_{10}` than aerobic decomposers; are not N limited; and do not have a low-moisture limitation.

When the single BGC soil level is used in CLM (Chapter :numref:`rst_Decomposition`), the temperature factor, :math:`f_{T}`, is set to 0 for temperatures equal to or below freezing, even though CLM allows heterotrophic respiration below freezing. However, if the vertically resolved BGC soil column is used, CH\ :sub:`4` production continues below freezing because liquid water stress limits decomposition. The base temperature for the :math:`{Q}_{10}` factor, :math:`{T}_{B}`, is 22°C and effectively modified the base :math:`f_{CH_{4}}` value.

For the single-layer BGC version, :math:`{R}_{H}` is distributed among soil levels by assuming that 50% is associated with the roots (using the CLM PFT-specific rooting distribution) and the rest is evenly divided among the top 0.28 m of soil (to be consistent with CLM's soil decomposition algorithm). For the vertically resolved BGC version, the prognosed distribution of :math:`{R}_{H}` is used to estimate CH\ :sub:`4` production.

The factor :math:`f_{pH}` is nominally set to 1, although a static spatial map of *pH* can be used to determine this factor (:ref:`Dunfield et al. 1993<Dunfieldetal1993>`) by applying:

.. math::
   :label: 24.3

   f_{pH} =10^{-0.2235pH^{2} +2.7727pH-8.6} .

The :math:`f_{pE}` factor assumes that alternative electron acceptors are reduced with an e-folding time of 30 days after inundation. The default version of the model applies this factor to horizontal changes in inundated area but not to vertical changes in the water table depth in the upland fraction of the gridcell. We consider both :math:`f_{pH}` and :math:`f_{pE}` to be poorly constrained in the model and identify these controllers as important areas for model improvement.

As a non-default option to account for CH\ :sub:`4` production in anoxic microsites above the water table, we apply the Arah and Stephen (1998) estimate of anaerobic fraction:

.. math::
   :label: 24.4

   \varphi =\frac{1}{1+\eta C_{O_{2} } } .

Here, :math:`\varphi` is the factor by which production is inhibited above the water table (compared to production as calculated in equation :eq:`24.2`, :math:`C_{O_{2}}` (mol m\ :sup:`-3`) is the bulk soil oxygen concentration, and :math:`\eta` = 400 mol m\ :sup:`-3`.

The O\ :sub:`2` required to facilitate the vertically resolved heterotrophic decomposition and root respiration is estimated assuming 1 mol O\ :sub:`2` is required per mol CO\ :sub:`2` produced. The model also calculates the O\ :sub:`2` required during nitrification, and the total O\ :sub:`2` demand is used in the O\ :sub:`2` mass balance solution.

.. _Table Methane Parameter descriptions:

.. table:: Parameter descriptions and sensitivity analysis ranges applied in the methane model

 +--------------+----------------------------+----------------------------------------------+--------------------------------------------------------------------------------------------------+---------------------------------------------+--------------------------------------------------------------------------------------------+
 | Mechanism    | Parameter                  | Baseline Value                               | Range for Sensitivity Analysis                                                                   | Units                                       | Description                                                                                |
 +==============+============================+==============================================+==================================================================================================+=============================================+============================================================================================+
 | Production   | :math:`{Q}_{10}`           | 2                                            | 1.5 – 4                                                                                          | -                                           | CH\ :sub:`4` production :math:`{Q}_{10}`                                                   |
 +--------------+----------------------------+----------------------------------------------+--------------------------------------------------------------------------------------------------+---------------------------------------------+--------------------------------------------------------------------------------------------+
 |              | :math:`f_{pH}`             | 1                                            | On, off                                                                                          | -                                           | Impact of pH on CH\ :sub:`4` production                                                    |
 +--------------+----------------------------+----------------------------------------------+--------------------------------------------------------------------------------------------------+---------------------------------------------+--------------------------------------------------------------------------------------------+
 |              | :math:`f_{pE}`             | 1                                            | On, off                                                                                          | -                                           | Impact of redox potential on CH\ :sub:`4` production                                       |
 +--------------+----------------------------+----------------------------------------------+--------------------------------------------------------------------------------------------------+---------------------------------------------+--------------------------------------------------------------------------------------------+
 |              | *S*                        | Varies                                       | NA                                                                                               | -                                           | Seasonal inundation factor                                                                 |
 +--------------+----------------------------+----------------------------------------------+--------------------------------------------------------------------------------------------------+---------------------------------------------+--------------------------------------------------------------------------------------------+
 |              | :math:`\beta`              | 0.2                                          | NA                                                                                               | -                                           | Effect of anoxia on decomposition rate (used to calculate *S* only)                        |
 +--------------+----------------------------+----------------------------------------------+--------------------------------------------------------------------------------------------------+---------------------------------------------+--------------------------------------------------------------------------------------------+
 |              | :math:`f_{CH_{4} }`        | 0.2                                          | NA                                                                                               | -                                           | Ratio between CH\ :sub:`4` and CO\ :sub:`2` production below the water table               |
 +--------------+----------------------------+----------------------------------------------+--------------------------------------------------------------------------------------------------+---------------------------------------------+--------------------------------------------------------------------------------------------+
 | Ebullition   | :math:`{C}_{e,max}`        | 0.15                                         | NA                                                                                               | mol m\ :sup:`-3`                            | CH\ :sub:`4` concentration to start ebullition                                             |
 +--------------+----------------------------+----------------------------------------------+--------------------------------------------------------------------------------------------------+---------------------------------------------+--------------------------------------------------------------------------------------------+
 |              | :math:`{C}_{e,min}`        | 0.15                                         | NA                                                                                               | -                                           | CH\ :sub:`4` concentration to end ebullition                                               |
 +--------------+----------------------------+----------------------------------------------+--------------------------------------------------------------------------------------------------+---------------------------------------------+--------------------------------------------------------------------------------------------+
 | Diffusion    | :math:`f_{D_{0} }`         | 1                                            | 1, 10                                                                                            | m\ :sup:`2` s\ :sup:`-1`                    | Diffusion coefficient multiplier (Table 24.2)                                              |
 +--------------+----------------------------+----------------------------------------------+--------------------------------------------------------------------------------------------------+---------------------------------------------+--------------------------------------------------------------------------------------------+
 | Aerenchyma   | *p*                        | 0.3                                          | NA                                                                                               | -                                           | Grass aerenchyma porosity                                                                  |
 +--------------+----------------------------+----------------------------------------------+--------------------------------------------------------------------------------------------------+---------------------------------------------+--------------------------------------------------------------------------------------------+
 |              | *R*                        | 2.9\ :math:`\times`\ 10\ :sup:`-3` m         | NA                                                                                               | m                                           | Aerenchyma radius                                                                          |
 +--------------+----------------------------+----------------------------------------------+--------------------------------------------------------------------------------------------------+---------------------------------------------+--------------------------------------------------------------------------------------------+
 |              | :math:`{r}_{L}`            | 3                                            | NA                                                                                               | -                                           | Root length to depth ratio                                                                 |
 +--------------+----------------------------+----------------------------------------------+--------------------------------------------------------------------------------------------------+---------------------------------------------+--------------------------------------------------------------------------------------------+
 |              | :math:`{F}_{a}`            | 1                                            | 0.5 – 1.5                                                                                        | -                                           | Aerenchyma conductance multiplier                                                          |
 +--------------+----------------------------+----------------------------------------------+--------------------------------------------------------------------------------------------------+---------------------------------------------+--------------------------------------------------------------------------------------------+
 | Oxidation    | :math:`K_{CH_{4} }`        | 5 x 10\ :sup:`-3`                            | 5\ :math:`\times`\ 10\ :math:`{}^{-4}`\ :math:`{}_{ }`- 5\ :math:`\times`\ 10\ :sup:`-2`         | mol m\ :sup:`-3`                            | CH\ :sub:`4` half-saturation oxidation coefficient (wetlands)                              |
 +--------------+----------------------------+----------------------------------------------+--------------------------------------------------------------------------------------------------+---------------------------------------------+--------------------------------------------------------------------------------------------+
 |              | :math:`K_{O_{2} }`         | 2 x 10\ :sup:`-2`                            | 2\ :math:`\times`\ 10\ :sup:`-3` - 2\ :math:`\times`\ 10\ :sup:`-1`                              | mol m\ :sup:`-3`                            | O\ :sub:`2` half-saturation oxidation coefficient                                          |
 +--------------+----------------------------+----------------------------------------------+--------------------------------------------------------------------------------------------------+---------------------------------------------+--------------------------------------------------------------------------------------------+
 |              | :math:`R_{o,\max }`        | 1.25 x 10\ :math:`{}^{-5}`                   | 1.25\ :math:`\times`\ 10\ :math:`{}^{-6}` - 1.25\ :math:`\times`\ 10\ :math:`{}^{-4}`            | mol m\ :sup:`-3` s\ :sup:`-1`               | Maximum oxidation rate (wetlands)                                                          |
 +--------------+----------------------------+----------------------------------------------+--------------------------------------------------------------------------------------------------+---------------------------------------------+--------------------------------------------------------------------------------------------+

Ebullition
---------------

Briefly, the simulated aqueous CH\ :sub:`4` concentration in each soil level is used to estimate the expected equilibrium gaseous partial pressure (:math:`C_{e}` ), as a function of temperature and depth below the water table, by first estimating the Henry's law partitioning coefficient (:math:`k_{h}^{C}` ) by the method described in :ref:`Wania et al. (2010)<Waniaetal2010>`:

.. math::
   :label: 24.5

   \log \left(\frac{1}{k_{H} } \right)=\log k_{H}^{s} -\frac{1}{C_{H} } \left(\frac{1}{T} -\frac{1}{T^{s} } \right)

.. math::
   :label: 24.6

   k_{h}^{C} =Tk_{H} R_{g}

.. math::
   :label: 24.7

   C_{e} =\frac{C_{w} R_{g} T}{\theta _{s} k_{H}^{C} p}

where :math:`C_{H}` \ is a constant, :math:`R_{g}` is the universal gas constant, :math:`k_{H}^{s}` is Henry's law partitioning coefficient at standard temperature (:math:`T^{s}` ),\ :math:`C_{w}` \ is local aqueous CH\ :sub:`4` concentration, and *p* is pressure.

The local pressure is calculated as the sum of the ambient pressure, water pressure down to the local depth, and pressure from surface ponding (if applicable). When the CH\ :sub:`4` partial pressure exceeds 15% of the local pressure (:ref:`Baird et al. 2004<Bairdetal2004>`; :ref:`Strack et al. 2006<Stracketal2006>`; :ref:`Wania et al. 2010<Waniaetal2010>`), bubbling occurs to remove CH\ :sub:`4` to below this value, modified by the fraction of CH\ :sub:`4` in the bubbles [taken as 57%; (:ref:`Kellner et al. 2006<Kellneretal2006>`; :ref:`Wania et al. 2010<Waniaetal2010>`)]. Bubbles are immediately added to the surface flux for saturated columns and are placed immediately above the water table interface in unsaturated columns.

.. _Aerenchyma Transport:

Aerenchyma Transport
-------------------------

Aerenchyma transport is modeled in CLM as gaseous diffusion driven by a concentration gradient between the specific soil layer and the atmosphere and, if specified, by vertical advection with the transpiration stream. There is evidence that pressure driven flow can also occur, but we did not include that mechanism in the current model.

The diffusive transport through aerenchyma (*A*, mol m\ :sup:`-2` s\ :sup:`-1`) from each soil layer is represented in the model as:

.. math::
   :label: 24.8

   A=\frac{C\left(z\right)-C_{a} }{{\raise0.7ex\hbox{$ r_{L} z $}\!\mathord{\left/ {\vphantom {r_{L} z D}} \right.}\!\lower0.7ex\hbox{$ D $}} +r_{a} } pT\rho _{r} ,

where *D* is the free-air gas diffusion coefficient (m\ :sup:`2` s\ :sup:`-1`); *C(z)* (mol m\ :sup:`-3`) is the gaseous concentration at depth *z* (m); :math:`r_{L}` is the ratio of root length to depth; *p* is the porosity (-); *T* is specific aerenchyma area (m\ :sup:`2` m\ :sup:`-2`); :math:`{r}_{a}` is the aerodynamic resistance between the surface and the atmospheric reference height (s m\ :sup:`-1`); and :math:`\rho _{r}` is the rooting density as a function of depth (-). The gaseous concentration is calculated with Henry's law as described in equation :eq:`24.7`.

Based on the ranges reported in :ref:`Colmer (2003)<Colmer2003>`, we have chosen baseline aerenchyma porosity values of 0.3 for grass and crop PFTs and 0.1 for tree and shrub PFTs:

.. math::
   :label: 24.9

   T=\frac{4 f_{N} N_{a}}{0.22} \pi R^{2} .

Here :math:`N_{a}` is annual net primary production (NPP, mol m\ :sup:`-2` s\ :sup:`-1`); *R* is the aerenchyma radius (2.9 :math:`\times`\ 10\ :sup:`-3` m); :math:`{f}_{N}` is the belowground fraction of annual NPP; and the 0.22 factor represents the amount of C per tiller. O\ :sub:`2` can also diffuse in from the atmosphere to the soil layer via the reverse of the same pathway, with the same representation as Equation :eq:`24.8` but with the gas diffusivity of oxygen.

CLM also simulates the direct emission of CH\ :sub:`4` from leaves to the atmosphere via transpiration of dissolved methane. We calculate this flux (:math:`F_{CH_{4} -T}`; mol m\ :math:`{}^{-}`\ :sup:`2` s\ :sup:`-1`) using the simulated soil water methane concentration (:math:`C_{CH_{4},j}` (mol m\ :sup:`-3`)) in each soil layer *j* and the CLM predicted transpiration (:math:`F_{T}` ) for each PFT, assuming that no methane was oxidized inside the plant tissue:

.. math::
   :label: 24.10

   F_{CH_{4} -T} =\sum _{j}\rho _{r,j} F_{T} C_{CH_{4} ,j}  .

.. _CH4 Oxidation:

CH\ :sub:`4` Oxidation
---------------------------------

CLM represents CH\ :sub:`4` oxidation with double Michaelis-Menten kinetics (:ref:`Arah and Stephen 1998<ArahStephen1998>`; :ref:`Segers 1998<Segers1998>`), dependent on both the gaseous CH\ :sub:`4` and O\ :sub:`2` concentrations:

.. math::
   :label: 24.11

   R_{oxic} =R_{o,\max } \left[\frac{C_{CH_{4} } }{K_{CH_{4} } +C_{CH_{4} } } \right]\left[\frac{C_{O_{2} } }{K_{O_{2} } +C_{O_{2} } } \right]Q_{10} F_{\vartheta }

where :math:`K_{CH_{4} }` and :math:`K_{O_{2} }` \ are the half saturation coefficients (mol m\ :sup:`-3`) with respect to CH\ :sub:`4` and O\ :sub:`2` concentrations, respectively; :math:`R_{o,\max }` is the maximum oxidation rate (mol m\ :sup:`-3` s\ :sup:`-1`); and :math:`{Q}_{10}` specifies the temperature dependence of the reaction with a base temperature set to 12 °C. The soil moisture limitation factor :math:`F_{\theta }` is applied above the water table to represent water stress for methanotrophs. Based on the data in :ref:`Schnell and King (1996)<SchnellKing1996>`, we take :math:`F_{\theta } = {e}^{-P/{P}_{c}}`, where *P* is the soil moisture potential and :math:`{P}_{c} = -2.4 \times {10}^{5}` mm.

.. _Reactive Transport Solution:

Reactive Transport Solution
--------------------------------

The solution to equation :eq:`24.11` is solved in several sequential steps: resolve competition for CH\ :sub:`4` and O\ :sub:`2` (section :numref:`Competition for CH4and O2`); add the ebullition flux into the layer directly above the water table or into the atmosphere; calculate the overall CH\ :sub:`4` or O\ :sub:`2` source term based on production, aerenchyma transport, ebullition, and oxidation; establish boundary conditions, including surface conductance to account for snow, ponding, and turbulent conductances and bottom flux condition (section :numref:`CH4 and O2 Source Terms`); calculate diffusivity (section :numref:`Aqueous and Gaseous Diffusion`); and solve the resulting mass balance using a tridiagonal solver (section :numref:`Crank-Nicholson Solution Methane`).

.. _Competition for CH4and O2:

Competition for CH\ :sub:`4` and O\ :sub:`2`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For each time step, the unlimited CH\ :sub:`4` and O\ :sub:`2` demands in each model depth interval are computed. If the total demand over a time step for one of the species exceeds the amount available in a particular control volume, the demand from each process associated with the sink is scaled by the fraction required to ensure non-negative concentrations. Since the methanotrophs are limited by both CH\ :sub:`4` and O\ :sub:`2`, the stricter limitation is applied to methanotroph oxidation, and then the limitations are scaled back for the other processes. The competition is designed so that the sinks must not exceed the available concentration over the time step, and if any limitation exists, the sinks must sum to this value. Because the sinks are calculated explicitly while the transport is semi-implicit, negative concentrations can occur after the tridiagonal solution. When this condition occurs for O\ :sub:`2`, the concentrations are reset to zero; if it occurs for CH\ :sub:`4`, the surface flux is adjusted and the concentration is set to zero if the adjustment is not too large.

.. _CH4 and O2 Source Terms:

CH\ :sub:`4` and O\ :sub:`2` Source Terms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The overall CH\ :sub:`4` net source term consists of production, oxidation at the base of aerenchyma, transport through aerenchyma, methanotrophic oxidation, and ebullition (either to the control volume above the water table if unsaturated or directly to the atmosphere if saturated). For O\ :sub:`2` below the top control volume, the net source term consists of O\ :sub:`2` losses from methanotrophy, SOM decomposition, and autotrophic respiration, and an O\ :sub:`2` source through aerenchyma.

.. _Aqueous and Gaseous Diffusion:

Aqueous and Gaseous Diffusion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For gaseous diffusion, we adopted the temperature dependence of molecular free-air diffusion coefficients (:math:`{D}_{0}` (m\ :sup:`2` s\ :sup:`-1`)) as described by :ref:`Lerman (1979) <Lerman1979>` and applied by :ref:`Wania et al. (2010)<Waniaetal2010>` (:numref:`Table Temperature dependence of aqueous and gaseous diffusion`).

.. _Table Temperature dependence of aqueous and gaseous diffusion:

.. table:: Temperature dependence of aqueous and gaseous diffusion coefficients for CH\ :sub:`4` and O\ :sub:`2`

 +----------------------------------------------------------+----------------------------------------------------------+--------------------------------------------------------+
 | :math:`{D}_{0}` (cm\ :sup:`2` s\ :sup:`-1`)              | CH\ :sub:`4`                                             | O\ :sub:`2`                                            |
 +==========================================================+==========================================================+========================================================+
 | Aqueous                                                  | 0.9798 + 0.02986\ *T* + 0.0004381\ *T*\ :sup:`2`         | 1.172+ 0.03443\ *T* + 0.0005048\ *T*\ :sup:`2`         |
 +----------------------------------------------------------+----------------------------------------------------------+--------------------------------------------------------+
 | Gaseous                                                  | 0.1875 + 0.0013\ *T*                                     | 0.1759 + 0.00117\ *T*                                  |
 +----------------------------------------------------------+----------------------------------------------------------+--------------------------------------------------------+

Gaseous diffusivity in soils also depends on the molecular diffusivity, soil structure, porosity, and organic matter content. :ref:`Moldrup et al. (2003)<Moldrupetal2003>`, using observations across a range of unsaturated mineral soils, showed that the relationship between effective diffusivity (:math:`D_{e}` (m\ :sup:`2` s\ :sup:`-1`)) and soil properties can be represented as:

.. math::
   :label: 24.12

   D_{e} =D_{0} \theta _{a}^{2} \left(\frac{\theta _{a} }{\theta _{s} } \right)^{{\raise0.7ex\hbox{$ 3 $}\!\mathord{\left/ {\vphantom {3 b}} \right.}\!\lower0.7ex\hbox{$ b $}} } ,

where :math:`\theta _{a}` and :math:`\theta _{s}` are the air-filled and total (saturated water-filled) porosities (-), respectively, and *b* is the slope of the water retention curve (-). However, :ref:`Iiyama and Hasegawa (2005)<IiyamaHasegawa2005>` have shown that the original Millington-Quirk (:ref:`Millington and Quirk 1961<MillingtonQuirk1961>`) relationship matched measurements more closely in unsaturated peat soils:

.. math::
   :label: 24.13

   D_{e} =D_{0} \frac{\theta _{a} ^{{\raise0.7ex\hbox{$ 10 $}\!\mathord{\left/ {\vphantom {10 3}} \right.}\!\lower0.7ex\hbox{$ 3 $}} } }{\theta _{s} ^{2} }

In CLM, we applied equation :eq:`24.12` for soils with zero organic matter content and equation :eq:`24.13` for soils with more than 130 kg m\ :sup:`-3` organic matter content. A linear interpolation between these two limits is applied for soils with SOM content below 130 kg m\ :sup:`-3`. For aqueous diffusion in the saturated part of the soil column, we applied (:ref:`Moldrup et al. (2003)<Moldrupetal2003>`):

.. math::
   :label: 24.14

   D_{e} =D_{0} \theta _{s} ^{2} .

To simplify the solution, we assumed that gaseous diffusion dominates above the water table interface and aqueous diffusion below the water table interface. Descriptions, baseline values, and dimensions for parameters specific to the CH\ :sub:`4` model are given in :numref:`Table Methane Parameter descriptions`. For freezing or frozen soils below the water table, diffusion is limited to the remaining liquid (CLM allows for some freezing point depression), and the diffusion coefficients are scaled by the volume-fraction of liquid. For unsaturated soils, Henry's law equilibrium is assumed at the interface with the water table.

.. _Boundary Conditions:

Boundary Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^

We assume the CH\ :sub:`4` and O\ :sub:`2` surface fluxes can be calculated from an effective conductance and a gaseous concentration gradient between the atmospheric concentration and either the gaseous concentration in the first soil layer (unsaturated soils) or in equilibrium with the water (saturated soil\ :math:`w\left(C_{1}^{n} -C_{a} \right)` and :math:`w\left(C_{1}^{n+1} -C_{a} \right)` for the fully explicit and fully implicit cases, respectively (however, see :ref:`Tang and Riley (2013)<TangRiley2013>` for a more complete representation of this process). Here, *w* is the surface boundary layer conductance as calculated in the existing CLM surface latent heat calculations. If the top layer is not fully saturated, the :math:`\frac{D_{m1} }{\Delta x_{m1} }` term is replaced with a series combination: :math:`\left[\frac{1}{w} +\frac{\Delta x_{1} }{D_{1} } \right]^{-1}`, and if the top layer is saturated, this term is replaced with :math:`\left[\frac{K_{H} }{w} +\frac{\frac{1}{2} \Delta x_{1} }{D_{1} } \right]^{-1}`, where :math:`{K}_{H}` is the Henry's law equilibrium constant.

When snow is present, a resistance is added to account for diffusion through the snow based on the Millington-Quirk expression :eq:`24.13` and CLM's prediction of the liquid water, ice, and air fractions of each snow layer. When the soil is ponded, the diffusivity is assumed to be that of methane in pure water, and the resistance as the ratio of the ponding depth to diffusivity. The overall conductance is taken as the series combination of surface, snow, and ponding resistances. We assume a zero flux gradient at the bottom of the soil column.

.. _Crank-Nicholson Solution Methane:

Crank-Nicholson Solution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Equation :eq:`24.1` is solved using a Crank-Nicholson solution (:ref:`Press et al., 1992<Pressetal1992>`), which combines fully explicit and implicit representations of the mass balance. The fully explicit decomposition of equation :eq:`24.1` can be written as

.. math::
   :label: 24.15

   \frac{R_{j}^{n+1} C_{j}^{n+1} -R_{j}^{n} C_{j}^{n} }{\Delta t} =\frac{1}{\Delta x_{j} } \left[\frac{D_{p1}^{n} }{\Delta x_{p1}^{} } \left(C_{j+1}^{n} -C_{j}^{n} \right)-\frac{D_{m1}^{n} }{\Delta x_{m1}^{} } \left(C_{j}^{n} -C_{j-1}^{n} \right)\right]+S_{j}^{n} ,

where *j* refers to the cell in the vertically discretized soil column (increasing downward), *n* refers to the current time step, :math:`\Delta`\ *t* is the time step (s), *p1* is *j+½*, *m1* is *j-½*, and :math:`S_{j}^{n}` is the net source at time step *n* and position *j*, i.e., :math:`S_{j}^{n} =P\left(j,n\right)-E\left(j,n\right)-A\left(j,n\right)-O\left(j,n\right)`. The diffusivity coefficients are calculated as harmonic means of values from the adjacent cells. Equation :eq:`24.15` is solved for gaseous and aqueous concentrations above and below the water table, respectively. The *R* term ensure the total mass balance in both phases is properly accounted for. An analogous relationship can be generated for the fully implicit case by replacing *n* by *n+1* on the *C* and *S* terms of equation :eq:`24.15`. Using an average of the fully implicit and fully explicit relationships gives:

.. math::
   :label: 24.16

   \begin{array}{l} {-\frac{1}{2\Delta x_{j} } \frac{D_{m1}^{} }{\Delta x_{m1}^{} } C_{j-1}^{n+1} +\left[\frac{R_{j}^{n+1} }{\Delta t} +\frac{1}{2\Delta x_{j} } \left(\frac{D_{p1}^{} }{\Delta x_{p1}^{} } +\frac{D_{m1}^{} }{\Delta x_{m1}^{} } \right)\right]C_{j}^{n+1} -\frac{1}{2\Delta x_{j} } \frac{D_{p1}^{} }{\Delta x_{p1}^{} } C_{j+1}^{n+1} =} \\ {\frac{R_{j}^{n} }{\Delta t} +\frac{1}{2\Delta x_{j} } \left[\frac{D_{p1}^{} }{\Delta x_{p1}^{} } \left(C_{j+1}^{n} -C_{j}^{n} \right)-\frac{D_{m1}^{} }{\Delta x_{m1}^{} } \left(C_{j}^{n} -C_{j-1}^{n} \right)\right]+\frac{1}{2} \left[S_{j}^{n} +S_{j}^{n+1} \right]} \end{array},

Equation :eq:`24.16` is solved with a standard tridiagonal solver, i.e.:

.. math::
   :label: 24.17

   aC_{j-1}^{n+1} +bC_{j}^{n+1} +cC_{j+1}^{n+1} =r,

with coefficients specified in equation :eq:`24.16`.

Two methane balance checks are performed at each timestep to insure that the diffusion solution and the time-varying aggregation over inundated and non-inundated areas strictly conserves methane molecules (except for production minus consumption) and carbon atoms.

.. _Interface between water table and unsaturated zone:

Interface between water table and unsaturated zone
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We assume Henry's Law equilibrium at the interface between the saturated and unsaturated zone and constant flux from the soil element below the interface to the center of the soil element above the interface. In this case, the coefficients are the same as described above, except for the soil element above the interface:

.. math:: \frac{D_{p1} }{\Delta x_{p1} } =\left[K_{H} \frac{\Delta x_{j} }{2D_{j} } +\frac{\Delta x_{j+1} }{2D_{j+1} } \right]^{-1}

.. math:: b=\left[\frac{R_{j}^{n+1} }{\Delta t} +\frac{1}{2\Delta x_{j} } \left(K_{H} \frac{D_{p1}^{} }{\Delta x_{p1} } +\frac{D_{m1}^{} }{\Delta x_{m1} } \right)\right]

.. math::
   :label: 24.18

   r=\frac{R_{j}^{n} }{\Delta t} C_{j}^{n} +\frac{1}{2\Delta x_{j} } \left[\frac{D_{p1}^{} }{\Delta x_{p1} } \left(C_{j+1}^{n} -K_{H} C_{j}^{n} \right)-\frac{D_{m1}^{} }{\Delta x_{m1} } \left(C_{j}^{n} -C_{j-1}^{n} \right)\right]+\frac{1}{2} \left[S_{j}^{n} +S_{j}^{n+1} \right]

and the soil element below the interface:

.. math:: \frac{D_{m1} }{\Delta x_{m1} } =\left[K_{H} \frac{\Delta x_{j-1} }{2D_{j-1} } +\frac{\Delta x_{j} }{2D_{j} } \right]^{-1}

.. math:: a=-K_{H} \frac{1}{2\Delta x_{j} } \frac{D_{m1}^{} }{\Delta x_{m1} }

.. math::
   :label: 24.19

   r=\frac{R_{j}^{n} }{\Delta t} +C_{j}^{n} +\frac{1}{2\Delta x_{j} } \left[\frac{D_{p1}^{} }{\Delta x_{p1} } \left(C_{j+1}^{n} -C_{j}^{n} \right)-\frac{D_{m1}^{} }{\Delta x_{m1} } \left(C_{j}^{n} -K_{H} C_{j-1}^{n} \right)\right]+\frac{1}{2} \left[S_{j}^{n} +S_{j}^{n+1} \right]

.. _Inundated Fraction Prediction:

Inundated Fraction Prediction
----------------------------------

A simplified dynamic representation of spatial inundation based on recent work by :ref:`Prigent et al. (2007)<Prigentetal2007>` is used. :ref:`Prigent et al. (2007)<Prigentetal2007>` described a multi-satellite approach to estimate the global monthly inundated fraction (:math:`{F}_{i}`) over an equal area grid (0.25 :math:`\circ` \ :math:`\times`\ 0.25\ :math:`\circ` at the equator) from 1993 - 2000. They suggested that the IGBP estimate for inundation could be used as a measure of sensitivity of their detection approach at low inundation. We therefore used the sum of their satellite-derived :math:`{F}_{i}` and the constant IGBP estimate when it was less than 10% to perform a simple inversion for the inundated fraction for methane production (:math:`{f}_{s}`). The method optimized two parameters (:math:`{fws}_{slope}` and :math:`{fws}_{intercept}`) for each grid cell in a simple model based on simulated total water storage (:math:`{TWS}`):

.. math::
   :label: 24.20

   f_{s} =fws_{slope} TWS  + fws_{intercept} .

These parameters were evaluated at the 0.5° resolution, and aggregated for coarser simulations. Ongoing work in the hydrology submodel of CLM may alleviate the need for this crude simplification of inundated fraction in future model versions.

.. _Seasonal Inundation:

Seasonal Inundation
------------------------

A simple scaling factor is used to mimic the impact of seasonal inundation on CH\ :sub:`4` production (see appendix B in :ref:`Riley et al. (2011a)<Rileyetal2011a>` for a discussion of this simplified expression):

.. math::
   :label: 24.21

   S=\frac{\beta \left(f-\bar{f}\right)+\bar{f}}{f} ,S\le 1.

Here, *f* is the instantaneous inundated fraction, :math:`\bar{f}` is the annual average inundated fraction (evaluated for the previous calendar year) weighted by heterotrophic respiration, and :math:`\beta` is the anoxia factor that relates the fully anoxic decomposition rate to the fully oxygen-unlimited decomposition rate, all other conditions being equal.

