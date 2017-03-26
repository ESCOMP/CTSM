.. math:: 188

**NCAR/TN-503+STR**

**NCAR Technical Note**

**\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_**

****

**July 2013**

****

**Technical Description of version 4.5 of the Community Land Model
(CLM)**

***Coordinating Lead Authors***

**Keith W. Oleson, David M. Lawrence**

****

***Lead Authors***

**Gordon B. Bonan, Beth Drewniak, Maoyi Huang, Charles D. Koven, Samuel
Levis, Fang Li, William J. Riley, Zachary M. Subin, Sean C. Swenson,
Peter E. Thornton**

****

***Contributing Authors***

**Anil Bozbiyik, Rosie Fisher, Colette L. Heald, Erik Kluzek,
Jean-Francois Lamarque, Peter J. Lawrence, L. Ruby Leung, William
Lipscomb, Stefan Muszala, Daniel M. Ricciuto, William Sacks, Ying Sun,
Jinyun Tang, Zong-Liang Yang**

****

****

****

**NCAR Earth System Laboratory**

**Climate and Global Dynamics
Division\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_**

**NATIONAL CENTER FOR ATMOSPHERIC RESEARCH**

**P. O. Box 3000**

**BOULDER, COLORADO 80307-3000**

**ISSN Print Edition 2153-2397**

**ISSN Electronic Edition 2153-2400**

****

**NCAR TECHNICAL NOTES**

http://library.ucar.edu/research/publish-technote\ ****

****

The Technical Notes series provides an outlet for a variety of NCAR
Manuscripts that contribute in specialized ways to the body of
scientific knowledge but that are not yet at a point of a formal
journal, monograph or book publication. Reports in this series are
issued by the NCAR scientific divisions, serviced by OpenSky and
operated through the NCAR Library. Designation symbols for the series
include:

**EDD – Engineering, Design, or Development Reports**

Equipment descriptions, test results, instrumentation,

and operating and maintenance manuals.

**IA – Instructional Aids**

** Instruction manuals, bibliographies, film supplements,

and other research or instructional aids.

**PPR – Program Progress Reports**

** Field program reports, interim and working reports,

survey reports, and plans for experiments.

**PROC – Proceedings**

** Documentation or symposia, colloquia, conferences,

workshops, and lectures. (Distribution maybe limited to

attendees).

**STR – Scientific and Technical Reports**

Data compilations, theoretical and numerical

investigations, and experimental results.

The National Center for Atmospheric Research (NCAR) is operated by the
nonprofit University Corporation for Atmospheric Research (UCAR) under
the sponsorship of the National Science Foundation. Any opinions,
findings, conclusions, or recommendations expressed in this publication
are those of the author(s) and do not necessarily reflect the views of
the National Science Foundation.

National Center for Atmospheric Research

P. O. Box 3000

Boulder, Colorado 80307-300

****

ii

**NCAR/TN-503+STR**

**NCAR Technical Note**

**\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_**

****

**July 2013**

****

**Technical Description of version 4.5 of the Community Land Model
(CLM)**

***Coordinating Lead Authors***

**Keith W. Oleson, David M. Lawrence**

****

***Lead Authors***

**Gordon B. Bonan, Beth Drewniak, Maoyi Huang, Charles D. Koven, Samuel
Levis, Fang Li, William J. Riley, Zachary M. Subin, Sean C. Swenson,
Peter E. Thornton**

****

***Contributing Authors***

**Anil Bozbiyik, Rosie Fisher, Colette L. Heald, Erik Kluzek,
Jean-Francois Lamarque, Peter J. Lawrence, L. Ruby Leung, William
Lipscomb, Stefan Muszala, Daniel M. Ricciuto, William Sacks, Ying Sun,
Jinyun Tang, Zong-Liang Yang**

****

****

****

**NCAR Earth System Laboratory**

**Climate and Global Dynamics
Division\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_**

**NATIONAL CENTER FOR ATMOSPHERIC RESEARCH**

**P. O. Box 3000**

**BOULDER, COLORADO 80307-3000**

**ISSN Print Edition 2153-2397**

**ISSN Electronic Edition 2153-2400**

xi

**TABLE OF CONTENTS**

****

**1. Introduction 11.1 Model History 11.1.1 Inception of CLM 11.1.2 CLM2
31.1.3 CLM3 51.1.4 CLM3.5 61.1.5 CLM4 71.1.6 CLM4.5 81.2 Biogeophysical
and Biogeochemical Processes 112. Surface Characterization and Model
Input Requirements 142.1 Surface Characterization 142.1.1 Surface
Heterogeneity and Data Structure 142.1.2 Vegetation Composition 172.1.3
Vegetation Structure 192.1.4 Phenology and vegetation burial by snow
212.2 Model Input Requirements 212.2.1 Atmospheric Coupling 212.2.2
Initialization 272.2.3 Surface Data 282.2.4 Adjustable Parameters and
Physical Constants 353. Surface Albedos 373.1 Canopy Radiative Transfer
373.2 Ground Albedos 463.2.1 Snow Albedo 483.2.2 Snowpack Optical
Properties 523.2.3 Snow Aging 563.3 Solar Zenith Angle 594. Radiative
Fluxes 634.1 Solar Fluxes 634.2 Longwave Fluxes 675. Momentum, Sensible
Heat, and Latent Heat Fluxes 715.1 Monin-Obukhov Similarity Theory 735.2
Sensible and Latent Heat Fluxes for Non-Vegetated Surfaces 825.3
Sensible and Latent Heat Fluxes and Temperature for Vegetated Surfaces
885.3.1 Theory 885.3.2 Numerical Implementation 1025.4 Update of Ground
Sensible and Latent Heat Fluxes 1075.5 Saturation Vapor Pressure 1106.
Soil and Snow Temperatures 1136.1 Numerical Solution 1146.2 Phase Change
1256.2.1 Soil and Snow Layers 1256.2.2 Surface Water 1296.3 Soil and
Snow Thermal Properties 1307. Hydrology 1347.1 Canopy Water 1357.2 Snow
1377.2.1 Snow Covered Area Fraction 1397.2.2 Ice Content 1407.2.3 Water
Content 1427.2.4 Black and organic carbon and mineral dust within snow
1437.2.5 Initialization of snow layer 1467.2.6 Snow Compaction 1467.2.7
Snow Layer Combination and Subdivision 1497.2.7.1 Combination 1497.2.7.2
Subdivision 1527.3 Surface Runoff, Surface Water Storage, and
Infiltration 1537.3.1 Surface Runoff 1537.3.2 Surface Water Storage
1557.3.3 Infiltration 1567.4 Soil Water 1587.4.1 Hydraulic Properties
1607.4.2 Numerical Solution 1637.4.2.1 Equilibrium soil matric potential
and volumetric moisture 1697.4.2.2 Equation set for layer i=1 1717.4.2.3
Equation set for layers i=2,…,N\ :sub:`levsoi` -1 1717.4.2.4 Equation
set for layers i=N\ :sub:`levsoi` ,…N\ :sub:`levsoi` +1 1727.5 Frozen
Soils and Perched Water Table 1747.6 Groundwater-Soil Water Interactions
1757.7 Runoff from glaciers and snow-capped surfaces 1787.8 The Variable
Infiltration Capacity parameterizations as a hydrologic option 1798.
Stomatal Resistance and Photosynthesis 1838.1 Stomatal resistance 1838.2
Photosynthesis 1868.3 Vcmax25 and canopy scaling 1918.4 Soil water
stress 1938.5 Numerical implementation 1979. Lake Model 2009.1
Discretization 2019.2 External Data 2029.3 Surface Albedo 2029.4 Surface
Fluxes and Surface Temperature 2039.4.1 Overview of Changes from CLM4
2039.4.2 Surface Properties 2039.4.3 Surface Flux Solution 2059.5 Lake
Temperature 2119.5.1 Introduction 2119.5.2 Overview of Changes from CLM4
2129.5.3 Boundary Conditions 2139.5.4 Eddy Diffusivity and Thermal
Conductivities 2139.5.5 Radiation Penetration 2169.5.6 Heat Capacities
2179.5.7 Crank-Nicholson Solution 2179.5.8 Phase Change 2199.5.9
Convection 2209.5.10 Energy Conservation 2239.6 Lake Hydrology 2239.6.1
Overview 2239.6.2 Water Balance 2249.6.3 Precipitation, Evaporation, and
Runoff 2259.6.4 Soil Hydrology 2269.6.5 Modifications to Snow Layer
Logic 22710. Glaciers 22910.1 Overview 22910.2 Multiple elevation class
scheme 23110.3 Computation of the surface mass balance 23211. River
Transport Model (RTM) 23512. Urban Model (CLMU) 23913. Carbon and
Nitrogen Pools, Allocation, and Respiration 24413.1 Introduction 24413.2
Carbon Allocation for Maintenance Respiration Costs 24613.3 Carbon and
Nitrogen Stoichiometry of New Growth 24813.4 Deployment of
retranslocated nitrogen 25213.5 Plant nitrogen uptake from soil mineral
nitrogen pool 25313.6 Final carbon and nitrogen allocation 25313.7
Autotrophic Respiration 25613.7.1 Maintenance Respiration 25613.7.2
Growth Respiration 25714. Vegetation Phenology 25914.1 General Phenology
Flux Parameterization 25914.1.1 Onset Periods 26014.1.2 Offset Periods
26214.1.3 Background Onset Growth 26414.1.4 Background Litterfall
26514.1.5 Livewood Turnover 26614.2 Evergreen Phenology 26714.3
Seasonal-Deciduous Phenology 26814.3.1 Seasonal-Deciduous Onset Trigger
26814.3.2 Seasonal-Deciduous Offset Trigger 27114.4 Stress-Deciduous
Phenology 27114.4.1 Stress-Deciduous Onset Triggers 27114.4.2
Stress-Deciduous Offset Triggers 27314.4.3 Stress-Deciduous: Long
Growing Season 27414.5 Litterfall Fluxes Merged to the Column Level
27615. Decomposition 27815.1 CLM-CN Pool Structure, Rate Constants and
Parameters 28115.2 Century-based Pool Structure, Rate Constants and
Parameters 28515.3 Environmental modifiers on decomposition rate 28615.4
N-limitation of Decomposition Fluxes 28915.5 N Competition between plant
uptake and soil immobilization fluxes 29215.6 Final Decomposition Fluxes
29315.7 Vertical Distribution and Transport of Decomposing C and N pools
29515.8 Model Equilibration 29616. External Nitrogen Cycle 29816.1
Atmospheric Nitrogen Deposition 29816.2 Biological Nitrogen Fixation
29916.3 Nitrification and Denitrification Losses of Nitrogen 30116.3.1
CLM-CN formulation 30116.3.2 Century-based formulation 30416.4 Leaching
Losses of Nitrogen 30516.5 Losses of Nitrogen Due to Fire 30717. Plant
Mortality 30817.1 Mortality Fluxes Leaving Vegetation Pools 30817.2
Mortality Fluxes Merged to the Column Level 31118. Fire 31618.1 Non-peat
fires outside cropland and tropical closed forest 31618.1.1 Fire counts
31618.1.2 Average spread area of a fire 32018.1.3 Fire impact 32318.2
Agricultural fires 32518.3 Deforestation fires 32618.4 Peat fires 32919.
Methane Model 33219.1 Methane Model Structure and Flow 33219.2 Governing
Mass-Balance Relationship 33319.3 CH4 Production 33419.4 Ebullition
33819.5 Aerenchyma Transport 33819.6 CH4 Oxidation 34019.7 Reactive
Transport Solution 34019.7.1 Competition for CH4 and O2 34119.7.2 CH4
and O2 Source Terms 34119.7.3 Aqueous and Gaseous Diffusion 34219.7.4
Boundary Conditions 34319.7.5 Crank-Nicholson Solution 34419.7.6
Interface between water table and unsaturated zone 34519.8 Inundated
Fraction Prediction 34619.9 Seasonal Inundation 34720. Crops and
Irrigation 34820.1 Summary of CLM4.5 updates relative to the CLM4.0
34820.2 The crop model 34820.2.1 Introduction 34820.2.2 Crop plant
functional types 34920.2.3 Phenology 35020.2.3.1 Planting 35120.2.3.2
Leaf emergence 35220.2.3.3 Grain fill 35320.2.3.4 Harvest 35320.2.4
Allocation 35320.2.4.1 Leaf emergence to grain fill 35420.2.4.2 Grain
fill to harvest 35420.2.5 General comments 35520.3 The irrigation model
36020.4 The details about what is new in CLM4.5 36120.4.1 Interactive
irrigation for corn, temperate cereals, and soybean 36120.4.2
Interactive fertilization 36320.4.3 Biological nitrogen fixation for
soybeans 36420.4.4 Modified C:N ratios for crops 36520.4.5 Nitrogen
retranslocation for crops 36520.4.6 Separate reproductive pool 36721.
Transient Landcover Change 36921.1 Annual Transient Land Cover Data and
Time Interpolation 36921.2 Mass and Energy Conservation 37121.3 Annual
Transient Land Cover Dataset Development 37221.3.1 UNH Transient Land
Use and Land Cover Change Dataset 37221.3.2 Representing Land Use and
Land Cover Change in CLM 37421.3.3 Present Day PFT Dataset 37521.3.4
Potential PFT Distribution 37621.3.5 Transient Land Cover Change Dataset
37721.3.6 Forest Harvest Dataset Changes 37822. Dynamic Global
Vegetation Model 38122.1 Establishment and survival 38222.2 Light
competition 38322.3 CN processes modified for the CNDV coupling 38323.
Biogenic Volatile Organic Compounds (BVOCs) 38624. Dust Model 38825.
Carbon Isotopes 39325.1 General Form for Calculating 13C and 14C Flux
39325.2 Isotope Symbols, Units, and Reference Standards 39425.3 Carbon
Isotope Discrimination During Photosynthesis 39625.4 14C radioactive
decay and historical atmospheric 14C concentrations 39826. Offline Mode
40027. References 405**

**LIST OF FIGURES**

Figure 1.1. Land biogeophysical, biogeochemical, and landscape processes
simulated by CLM (adapted from Lawrence et al. (2011) for CLM4.5).
13Figure 2.1. Configuration of the CLM subgrid hierarchy. 15Figure 4.1.
Schematic diagram of (a) direct beam radiation, (b) diffuse solar
radiation, and (c) longwave radiation absorbed, transmitted, and
reflected by vegetation and ground. 64Figure 5.1. Schematic diagram of
sensible heat fluxes for (a) non-vegetated surfaces and (b) vegetated
surfaces. 91Figure 5.2. Schematic diagram of water vapor fluxes for (a)
non-vegetated surfaces and (b) vegetated surfaces. 92Figure 6.1.
Schematic diagram of numerical scheme used to solve for soil
temperature. 119Figure 7.1. Hydrologic processes represented in CLM.
135Figure 7.2. Example of three layer snow pack (snl=-3). 137Figure 7.3.
Schematic diagram of numerical scheme used to solve for soil water
fluxes. 165Figure 12.1. Schematic representation of the urban land unit.
242Figure 12.2. Schematic of urban and atmospheric model coupling.
243Figure 13.1. Vegetation fluxes and pools. 245Figure 13.2: Carbon and
nitrogen pools. 246Figure 14.1. Example of annual phenology cycle for
seasonal deciduous. 260Figure 14.2. Example fluxes and pools sizes for
an onset growth period of 15 days, with initial transfer pool size of
100 gC m-2 and a timestep of one hour. a) Flux leaving transfer pool
(e.g. CFleaf\_xfer,leaf). b) Carbon content of transfer pool and its
associated display pool (e.g. CSleaf\_xfer and CSleaf, respectively).
262Figure 14.3. Example fluxes and pool sizes for an offset (litterfall)
period of 15 days, with initial display pool size of 100 gC m-2 and a
timestep of one hour. a) Litterfall flux (e.g CFleaf,litter). b) Carbon
content of display pool and litter pool through the litterfall period,
ignoring the losses from litter pool due to decomposition during this
period. 264Figure 15.1. Schematic of decomposition model in CLM.
279Figure 15.2. Pool structure, transitions, respired fractions (numbers
at end of arrows), and turnover times (numbers in boxes) for the 2
alternate soil decomposition models included in CLM. 281Figure 16.1.
Biological nitrogen fixation as a function of annual net primary
production. 300Figure 19.1. Schematic representation of biological and
physical processes integrated in CLM that affect the net CH4 surface
flux. (left) Fully inundated portion of a CLM gridcell and (right)
variably saturated portion of a gridcell. 334Figure 21.1. Schematic of
land cover change impacts on CLM carbon pools and fluxes. 379Figure
21.2. Schematic of translation of annual UNH land units to CLM4 plant
functional types. 380Figure 25.1. Atmospheric
:math:`\mathrm{\Delta}`\ 14C used to drive 14C model over the historical
period. 399\ ****

**LIST OF TABLES**

Table 2.1. Plant functional types 18Table 2.2. Prescribed plant
functional type heights 20Table 2.3. Atmospheric input to land model
23Table 2.4. Land model output to atmospheric model 26Table 2.5. Surface
data required for CLM4.5 and their base spatial resolution 29Table 2.6.
Physical constants 36Table 3.1. Plant functional type optical properties
44Table 3.2. Intercepted snow optical properties 46Table 3.3. Dry and
saturated soil albedos 48Table 3.4. Spectral bands and weights used for
snow radiative transfer 51Table 3.5. Single-scatter albedo values used
for snowpack impurities and ice 54Table 3.6. Mass extinction values (m2
kg-1) used for snowpack impurities and ice. 55Table 3.7. Asymmetry
scattering parameters used for snowpack impurities and ice. 56Table 3.8.
Orbital parameters 62Table 5.1. Plant functional type aerodynamic
parameters 101Table 5.2. Coefficients for e\ :sub:`sat`\ :sup:`T`
111Table 5.3. Coefficients for 112Table 6.1. Soil layer structure.
115Table 7.1. Meltwater scavenging efficiency for particles within snow
146Table 7.2. Minimum and maximum thickness of snow layers (m) 152Table
8.1. Plant functional type (PFT) photosynthetic parameters. 185Table
8.2. Temperature dependence parameters for C3 photosynthesis. 190Table
8.3. Plant functional type root distribution parameters. 196Table 13.1.
Allocation and carbon:nitrogen ratio parameters 250Table 15.1.
Decomposition rate constants for litter and SOM pools, C:N ratios, and
acceleration parameters (see section 15.8 for explanation) for the
CLM-CN decomposition pool structure. 283Table 15.2. Respiration
fractions for litter and SOM pools 284Table 15.3. Respiration fractions
for litter and SOM pools for Century-based structure 285Table 15.4.
Turnover times, C:N ratios, and acceleration parameters (see section
15.8 for explanation) for the Century-based decomposition cascade.
286Table 18.1. PFT-specific combustion completeness and fire mortality
factors. 331Table 19.1. Parameter descriptions and sensitivity analysis
ranges applied in the methane model. 337Table 19.2. Temperature
dependence of aqueous and gaseous diffusion coefficients for CH4 and O2.
342Table 20.1. Crop plant functional types (pfts) in CLM4.5CNcrop and
their parameters relating to phenology and morphology. Numbers in the
first column correspond to the list of pfts in Table 2.1. 357Table 20.2.
Crop pfts in CLM4.5CNcrop and their parameters relating to allocation.
Numbers in the first column correspond to the list of pfts in Table 2.1.
359Table 20.3. Pre- and post-grain fill C:N ratios for crop leaf, stem,
fine root, and reproductive pools. 367Table 22.1. Plant functional type
(PFT) biogeography rules with respect to climate. 385Table 24.1. Mass
fraction m\ :sub:`i` , mass median diameter :sub:`v, i` , and geometric
standard deviation :sub:`g, i` , per dust source mode i 392Table 24.2.
Minimum and maximum particle diameters in each dust transport bin j
392\ ****

ACKNOWLEDGEMENTS

The authors would like to acknowledge the substantial contributions of
the following members of the Land Model and Biogeochemistry Working
Groups to the development of the Community Land Model since its
inception in 1996: Benjamin Andre, Ian Baker, Michael Barlage, Mike
Bosilovich, Marcia Branstetter, Tony Craig, Aiguo Dai, Yongjiu Dai, Mark
Decker, Scott Denning, Robert Dickinson, Paul Dirmeyer, Jared Entin, Jay
Famiglietti, Johannes Feddema, Mark Flanner, Jon Foley, Andrew Fox, Inez
Fung, David Gochis, Alex Guenther, Tim Hoar, Forrest Hoffman, Paul
Houser, Trish Jackson, Brian Kauffman, Silvia Kloster, Natalie Mahowald,
Jiafu Mao, Lei Meng, Sheri Michelson, Guo-Yue Niu, Adam Phillips, Taotao
Qian, Jon Radakovich, James Randerson, Nan Rosenbloom, Steve Running,
Koichi Sakaguchi, Adam Schlosser, Andrew Slater, Reto Stöckli, Quinn
Thomas, Mariana Vertenstein, Nicholas Viovy, Aihui Wang, Guiling Wang,
Charlie Zender, Xiaodong Zeng, and Xubin Zeng.

The authors also thank the following people for their review of this
document: Jonathan Buzan, Kyla Dahlin, Sanjiv Kumar, Hanna Lee, Danica
Lombardozzi, Quinn Thomas, and Will Wieder.

Current affiliations for the authors are as follows:

K.W. Oleson, D.M. Lawrence, G.B. Bonan, S. Levis, S.C. Swenson, R.
Fisher, E. Kluzek, J.-F. Lamarque, P.J. Lawrence, S. Muszala, and W.
Sacks (National Center for Atmospheric Research); B. Drewniak (Argonne
National Laboratory); M. Huang, L.R. Leung (Pacific Northwest National
Laboratory); C.D. Koven, W.J. Riley, and J. Tang (Lawrence Berkeley
National Laboratory); F. Li (Chinese Academy of Sciences); Z.M. Subin
(Princeton University); P.E. Thornton and D.M. Ricciuto (Oak Ridge
National Laboratory); A. Bozbiyik (Bern University); C. Heald
(Massachusetts Institute of Technology), W. Lipscomb (Los Alamos
National Laboratory); Ying Sun and Z.-L. Yang (University of Texas at
Austin)

.. math:: 188

.. math:: 7

Introduction
===============

The purpose of this technical note is to describe the biogeophysical and
biogeochemical parameterizations and numerical implementation of version
4.5 of the Community Land Model (CLM4.5). Scientific justification and
evaluation of these parameterizations can be found in the referenced
scientific papers (Chapter 27). This technical note and the CLM4.5
User’s Guide together provide the user with the scientific description
and operating instructions for CLM.

Model History 
---------------

Inception of CLM
^^^^^^^^^^^^^^^^^^^^^^

The early development of the Community Land Model can be described as
the merging of a community-developed land model focusing on
biogeophysics and a concurrent effort at NCAR to expand the NCAR Land
Surface Model (NCAR LSM, Bonan 1996) to include the carbon cycle,
vegetation dynamics, and river routing. The concept of a
community-developed land component of the Community Climate System Model
(CCSM) was initially proposed at the CCSM Land Model Working Group
(LMWG) meeting in February 1996. Initial software specifications and
development focused on evaluating the best features of three existing
land models: the NCAR LSM (Bonan 1996, 1998) used in the Community
Climate Model (CCM3) and the initial version of CCSM; the Institute of
Atmospheric Physics, Chinese Academy of Sciences land model (IAP94) (Dai
and Zeng 1997); and the Biosphere-Atmosphere Transfer Scheme (BATS)
(Dickinson et al. 1993) used with CCM2. A scientific steering committee
was formed to review the initial specifications of the design provided
by Robert Dickinson, Gordon Bonan, Xubin Zeng, and Yongjiu Dai and to
facilitate further development. Steering committee members were selected
so as to provide guidance and expertise in disciplines not generally
well-represented in land surface models (e.g., carbon cycling,
ecological modeling, hydrology, and river routing) and included
scientists from NCAR, the university community, and government
laboratories (R. Dickinson, G. Bonan, X. Zeng, Paul Dirmeyer, Jay
Famiglietti, Jon Foley, and Paul Houser).

The specifications for the new model, designated the Common Land Model,
were discussed and agreed upon at the June 1998 CCSM Workshop LMWG
meeting. An initial code was developed by Y. Dai and was examined in
March 1999 by Mike Bosilovich, P. Dirmeyer, and P. Houser. At this point
an extensive period of code testing was initiated. Keith Oleson, Y. Dai,
Adam Schlosser, and P. Houser presented preliminary results of offline
1-dimensional testing at the June 1999 CCSM Workshop LMWG meeting.
Results from more extensive offline testing at plot, catchment, and
large scale (up to global) were presented by Y. Dai, A. Schlosser, K.
Oleson, M. Bosilovich, Zong-Liang Yang, Ian Baker, P. Houser, and P.
Dirmeyer at the LMWG meeting hosted by COLA (Center for
Ocean-Land-Atmosphere Studies) in November 1999. Field data used for
validation included sites adopted by the Project for Intercomparison of
Land-surface Parameterization Schemes (Henderson-Sellers et al. 1993)
(Cabauw, Valdai, Red-Arkansas river basin) and others [FIFE (Sellers et
al. 1988), BOREAS (Sellers et al. 1995), HAPEX-MOBILHY (André et al.
1986), ABRACOS (Gash et al. 1996), Sonoran Desert (Unland et al. 1996),
GSWP (Dirmeyer et al. 1999)]. Y. Dai also presented results from a
preliminary coupling of the Common Land Model to CCM3, indicating that
the land model could be successfully coupled to a climate model.

Results of coupled simulations using CCM3 and the Common Land Model were
presented by X. Zeng at the June 2000 CCSM Workshop LMWG meeting.
Comparisons with the NCAR LSM and observations indicated major
improvements to the seasonality of runoff, substantial reduction of a
summer cold bias, and snow depth. Some deficiencies related to runoff
and albedo were noted, however, that were subsequently addressed. Z.-L.
Yang and I. Baker demonstrated improvements in the simulation of snow
and soil temperatures. Sam Levis reported on efforts to incorporate a
river routing model to deliver runoff to the ocean model in CCSM. Soon
after the workshop, the code was delivered to NCAR for implementation
into the CCSM framework. Documentation for the Common Land Model is
provided by Dai et al. (2001) while the coupling with CCM3 is described
in Zeng et al. (2002). The model was introduced to the modeling
community in Dai et al. (2003).

CLM2
^^^^^^^^^^

Concurrent with the development of the Common Land Model, the NCAR LSM
was undergoing further development at NCAR in the areas of carbon
cycling, vegetation dynamics, and river routing. The preservation of
these advancements necessitated several modifications to the Common Land
Model. The biome-type land cover classification scheme was replaced with
a plant functional type (PFT) representation with the specification of
PFTs and leaf area index from satellite data (Oleson and Bonan 2000;
Bonan et al. 2002a, b). This also required modifications to
parameterizations for vegetation albedo and vertical burying of
vegetation by snow. Changes were made to canopy scaling, leaf
physiology, and soil water limitations on photosynthesis to resolve
deficiencies indicated by the coupling to a dynamic vegetation model.
Vertical heterogeneity in soil texture was implemented to improve
coupling with a dust emission model. A river routing model was
incorporated to improve the fresh water balance over oceans. Numerous
modest changes were made to the parameterizations to conform to the
strict energy and water balance requirements of CCSM. Further
substantial software development was also required to meet coding
standards. The resulting model was adopted in May 2002 as the Community
Land Model (CLM2) for use with the Community Atmosphere Model (CAM2, the
successor to CCM3) and version 2 of the Community Climate System Model
(CCSM2).

K. Oleson reported on initial results from a coupling of CCM3 with CLM2
at the June 2001 CCSM Workshop LMWG meeting. Generally, the CLM2
preserved most of the improvements seen in the Common Land Model,
particularly with respect to surface air temperature, runoff, and snow.
These simulations are documented in Bonan et al. (2002a). Further small
improvements to the biogeophysical parameterizations, ongoing software
development, and extensive analysis and validation within CAM2 and CCSM2
culminated in the release of CLM2 to the community in May 2002.

Following this release, Peter Thornton implemented changes to the model
structure required to represent carbon and nitrogen cycling in the
model. This involved changing data structures from a single vector of
spatially independent sub-grid patches to one that recognizes three
hierarchical scales within a model grid cell: land unit, snow/soil
column, and PFT. Furthermore, as an option, the model can be configured
so that PFTs can share a single soil column and thus “compete” for
water. This version of the model (CLM2.1) was released to the community
in February 2003. CLM2.1, without the compete option turned on, produced
only round off level changes when compared to CLM2.

CLM3
^^^^^^^^^^

CLM3 implemented further software improvements related to performance
and model output, a re-writing of the code to support vector-based
computational platforms, and improvements in biogeophysical
parameterizations to correct deficiencies in the coupled model climate.
Of these parameterization improvements, two were shown to have a
noticeable impact on simulated climate. A variable aerodynamic
resistance for heat/moisture transfer from ground to canopy air that
depends on canopy density was implemented. This reduced unrealistically
high surface temperatures in semi-arid regions. The second improvement
added stability corrections to the diagnostic 2-m air temperature
calculation which reduced biases in this temperature. Competition
between PFTs for water, in which PFTs share a single soil column, is the
default mode of operation in this model version. CLM3 was released to
the community in June 2004. Dickinson et al. (2006) describe the climate
statistics of CLM3 when coupled to CCSM3.0. Hack et al. (2006) provide
an analysis of selected features of the land hydrological cycle.
Lawrence et al. (2007) examine the impact of changes in CLM3
hydrological parameterizations on partitioning of evapotranspiration
(ET) and its effect on the timescales of ET response to precipitation
events, interseasonal soil moisture storage, soil moisture memory, and
land-atmosphere coupling. Qian et al. (2006) evaluate CLM3’s performance
in simulating soil moisture content, runoff, and river discharge when
forced by observed precipitation, temperature and other atmospheric
data.

CLM3.5
^^^^^^^^^^^^

Although the simulation of land surface climate by CLM3 was in many ways
adequate, most of the unsatisfactory aspects of the simulated climate
noted by the above studies could be traced directly to deficiencies in
simulation of the hydrological cycle. In 2004, a project was initiated
to improve the hydrology in CLM3 as part of the development of CLM
version 3.5. A selected set of promising approaches to alleviating the
hydrologic biases in CLM3 were tested and implemented. These included
new surface datasets based on Moderate Resolution Imaging
Spectroradiometer (MODIS) products, new parameterizations for canopy
integration, canopy interception, frozen soil, soil water availability,
and soil evaporation, a TOPMODEL-based model for surface and subsurface
runoff, a groundwater model for determining water table depth, and the
introduction of a factor to simulate nitrogen limitation on plant
productivity. Oleson et al. (2008a) show that CLM3.5 exhibits
significant improvements over CLM3 in its partitioning of global ET
which result in wetter soils, less plant water stress, increased
transpiration and photosynthesis, and an improved annual cycle of total
water storage. Phase and amplitude of the runoff annual cycle is
generally improved. Dramatic improvements in vegetation biogeography
result when CLM3.5 is coupled to a dynamic global vegetation model.
Stöckli et al. (2008) examine the performance of CLM3.5 at local scales
by making use of a network of long-term ground-based ecosystem
observations [FLUXNET (Baldocchi et al. 2001)]. Data from 15 FLUXNET
sites were used to demonstrate significantly improved soil hydrology and
energy partitioning in CLM3.5. CLM3.5 was released to the community in
May, 2007.

CLM4
^^^^^^^^^^

The motivation for the next version of the model, CLM4, was to
incorporate several recent scientific advances in the understanding and
representation of land surface processes, expand model capabilities, and
improve surface and atmospheric forcing datasets (Lawrence et al. 2011).
Included in the first category are more sophisticated representations of
soil hydrology and snow processes. In particular, new treatments of soil
column-groundwater interactions, soil evaporation, aerodynamic
parameters for sparse/dense canopies, vertical burial of vegetation by
snow, snow cover fraction and aging, black carbon and dust deposition,
and vertical distribution of solar energy for snow were implemented.
Major new capabilities in the model include a representation of the
carbon-nitrogen cycle (CLM4CN, see next paragraph for additional
information), the ability to model land cover change in a transient
mode, inclusion of organic soil and deep soil into the existing mineral
soil treatment to enable more realistic modeling of permafrost, an urban
canyon model to contrast rural and urban energy balance and climate
(CLMU), and an updated biogenic volatile organic compounds (BVOC) model.
Other modifications of note include refinement of the global PFT,
wetland, and lake distributions, more realistic optical properties for
grasslands and croplands, and an improved diurnal cycle and spectral
distribution of incoming solar radiation to force the model in offline
mode.

Many of the ideas incorporated into the carbon and nitrogen cycle
component of CLM4 derive from the earlier development of the offline
ecosystem process model Biome-BGC (Biome BioGeochemical Cycles),
originating at the Numerical Terradynamic Simulation Group (NTSG) at the
University of Montana, under the guidance of Prof. Steven Running.
Biome-BGC itself is an extension of an earlier model, Forest-BGC
**(Running and Coughlan, 1988; Running and Gower, 1991)**, which
simulates water, carbon, and, to a limited extent, nitrogen fluxes for
forest ecosystems. Forest-BGC was designed to be driven by remote
sensing inputs of vegetation structure, and so used a diagnostic
(prescribed) leaf area index, or, in the case of the dynamic allocation
version of the model **(Running and Gower, 1991)**, prescribed maximum
leaf area index.

Biome-BGC expanded on the Forest-BGC logic by introducing a more
mechanistic calculation of leaf and canopy scale photosynthesis **(Hunt
and Running, 1992)**, and extending the physiological parameterizations
to include multiple woody and non-woody vegetation types **(Hunt et al.
1996; Running and Hunt, 1993)**. Later versions of Biome-BGC introduced
more mechanistic descriptions of belowground carbon and nitrogen cycles,
nitrogen controls on photosynthesis and decomposition, sunlit and shaded
canopies, vertical gradient in leaf morphology, and explicit treatment
of fire and harvest disturbance and regrowth dynamics **(Kimball et al.
1997; Thornton, 1998; Thornton et al. 2002; White et al. 2000)**.
Biome-BGC version 4.1.2 **(Thornton et al. 2002)** provided a point of
departure for integrating new biogeochemistry components into CLM4.

CLM4 was released to the community in June, 2010 along with the
Community Climate System Model version 4 (CCSM4). CLM4 is used in CCSM4,
CESM1, CESM1.1, and remains available as the default land component
model option for coupled simulations in CESM1.2.

CLM4.5
^^^^^^^^^^^^

The motivations for the development of CLM4.5 (the model version
described in this Technical Description) were similar to those for CLM4:
incorporate several recent scientific advances in the understanding and
representation of land surface processes, expand model capabilities, and
improve surface and atmospheric forcing datasets.

Specifically, several parameterizations were revised to reflect new
scientific understanding and in an attempt to reduce biases identified
in CLM4 simulations including low soil carbon stocks especially in the
Arctic, excessive tropical GPP and unrealistically low Arctic GPP, a dry
soil bias in Arctic soils, unrealistically high LAI in the tropics, a
transient 20\ :math:`{}^{th}` century carbon response that was
inconsistent with observational estimates, and several other more minor
problems or biases.

The main modifications include updates to canopy processes including a
revised canopy radiation scheme and canopy scaling of leaf processes,
co-limitations on photosynthesis, revisions to photosynthetic parameters
(Bonan et al. 2011), **** temperature acclimation of photosynthesis, and
improved stability of the iterative solution in the photosynthesis and
stomatal conductance model (Sun et al. 2012). Hydrology updates include
modifications such that hydraulic properties of frozen soils are
determined by liquid water content only rather than total water content
and the introduction of an ice impedance function, and other corrections
that increase the consistency between soil water state and water table
position and allow for a perched water table above icy permafrost ground
(Swenson et al. 2012). A new snow cover fraction parameterization is
incorporated that reflects the hysteresis in fractional snow cover for a
given snow depth between accumulation and melt phases (Swenson and
Lawrence, 2012). The lake model in CLM4 is replaced with a completely
revised and more realistic lake model (Subin et al. 2012a). A surface
water store is introduced, replacing the wetland land unit and
permitting prognostic wetland distribution modeling, and the surface
energy fluxes are calculated separately (Swenson and Lawrence, 2012) for
snow-covered, water-covered, and snow/water-free portions of vegetated
and crop land units, and snow-covered and snow-free portions of glacier
land units. Globally constant river flow velocity is replaced with
variable flow velocity based on mean grid cell slope. A vertically
resolved soil biogeochemistry scheme is introduced with base
decomposition rates modified by soil temperature, water, and oxygen
limitations and also including vertical mixing of soil carbon and
nitrogen due to bioturbation, cryoturbation, and diffusion (Koven et al.
2013). The litter and soil carbon and nitrogen pool structure as well as
nitrification and denitrification are modified based on the Century
model and biological fixation is revised to distribute fixation more
realistically over the year (Koven et al. 2013). The fire model is
replaced with a model that includes representations of natural and
anthropogenic triggers and suppression as well as agricultural,
deforestation, and peat fires (Li et al. 2012a,b; Li et al. 2013a). The
biogenic volatile organic compounds model is updated to MEGAN2.1
(Guenther et al. 2012).

Additions to the model include a methane production, oxidation, and
emissions model (Riley et al. 2011a) and an extension of the crop model
to include interactive fertilization, organ pools (Drewniak et al.
2013), and irrigation (Sacks et al. 2009). Elements of the Variable
Infiltration Capacity (VIC) model are included as an alternative
optional runoff generation scheme (Li et al. 2011). There is also an
option to run with a multilayer canopy (Bonan et al. 2012). Multiple
urban density classes, rather than the single dominant urban density
class used in CLM4, are modeled in the urban land unit. Carbon
(:math:`{}^{13}`\ C and :math:`{}^{14}`\ C) isotopes are enabled (Koven
et al. 2013). Minor changes include a switch of the C3 Arctic grass and
shrub phenology from stress deciduous to seasonal deciduous and a change
in the glacier bare ice albedo to better reflect recent estimates.
Finally, the carbon and nitrogen cycle spinup is accelerated and
streamlined with a revised spinup method, though the spinup timescale
remains long.

Finally, the predominantly low resolution input data for provided with
CLM4 to create CLM4 surface datasets is replaced with newer and higher
resolution input datasets where possible (see section 2.2.3 for
details). The default meteorological forcing dataset provided with CLM4
(Qian et al. 2006) is replaced with the 1901-2010 CRUNCEP forcing
dataset (see Chapter 26) for CLM4.5, though users can also still use the
Qian et al. (2006) dataset or other alternative forcing datasets.

CLM4.5 was released to the community in June 2013 along with the
Community Earth System Model version 1.2 (CESM1.2).

Biogeophysical and Biogeochemical Processes
-----------------------------------------------

Biogeophysical and biogeochemical processes are simulated for each
subgrid land unit, column, and plant functional type (PFT) independently
and each subgrid unit maintains its own prognostic variables (see
section 2.1.1 for definitions of subgrid units). The same atmospheric
forcing is used to force all subgrid units within a grid cell. The
surface variables and fluxes required by the atmosphere are obtained by
averaging the subgrid quantities weighted by their fractional areas. The
processes simulated include (Figure 1.1):

#. Surface characterization including land type heterogeneity and
   ecosystem structure (Chapter 2)

#. Absorption, reflection, and transmittance of solar radiation (Chapter
   3, 4)

#. Absorption and emission of longwave radiation (Chapter 4)

#. Momentum, sensible heat (ground and canopy), and latent heat (ground
   evaporation, canopy evaporation, transpiration) fluxes (Chapter 5)

#. Heat transfer in soil and snow including phase change (Chapter 6)

#. Canopy hydrology (interception, throughfall, and drip) (Chapter 7)

#. Snow hydrology (snow accumulation and melt, compaction, water
   transfer between snow layers) (Chapter 7)

#. Soil hydrology (surface runoff, infiltration, redistribution of water
   within the column, sub-surface drainage, groundwater) (Chapter 7)

#. Stomatal physiology and photosynthesis (Chapter 8)

#. place Lake temperatures and fluxes (Chapter 9)

#. Glacier processes (Chapter 10)

#. Routing of runoff from rivers to ocean (Chapter 11)

#. Urban energy balance and climate (Chapter 12)

#. Vegetation carbon and nitrogen allocation and respiration (Chapter
   13)

#. Vegetation phenology (Chapter 14)

#. Soil and litter carbon decomposition (Chapter 15)

#. Nitrogen cycling including deposition, biological fixation,
   denitrification, leaching, and losses due to fire (Chapter 16)

#. Plant mortality (Chapter 17)

#. Fire ignition and suppression, including natural, deforestation, and
   agricultural fire (Chapter 18)

#. Methane production, oxidation, and emissions (Chapter 19)

#. Crop dynamics and irrigation (Chapter 20)

#. Land cover and land use change including wood harvest (Chapter 21)

#. Dynamic global vegetation distribution (Chapter 22)

#. Biogenic volatile organic compound emissions (Chapter 23)

#. Dust mobilization and deposition (Chapter 24)

#. Carbon isotope fractionation (Chapter 25)

Figure 1.1. Land biogeophysical, biogeochemical, and landscape processes
simulated by CLM (adapted from Lawrence et al. (2011) for CLM4.5).

**|image|**

.. |image| image:: image1
