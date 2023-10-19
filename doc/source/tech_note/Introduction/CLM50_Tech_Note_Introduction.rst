**February 2018**

**Technical Description of version 5.0 of the Community Land Model
(CLM)**

***Coordinating Lead Authors***

**David Lawrence, Rosie Fisher, Charles Koven, Keith Oleson, Sean Swenson, Mariana Vertenstein**

***Lead Authors***

**Ben Andre, Gordon Bonan, Bardan Ghimire, Leo van Kampenhout, Daniel Kennedy, Erik Kluzek, Ryan Knox, Peter Lawrence, Fang Li, Hongyi Li, Danica Lombardozzi, Yaqiong Lu, Justin Perket, William Riley, William Sacks, Mingjie Shi, Will Wieder, Chonggang Xu**

***Contributing Authors***

**Ashehad Ali, Andrew Badger, Gautam Bisht, Patrick Broxton, Michael Brunke, Jonathan Buzan, Yanyan Cheng, Martyn Clark, Tony Craig, Kyla Dahlin, Beth Drewniak, Louisa Emmons, Josh Fisher, Mark Flanner, Pierre Gentine, Maoyi Huang, Jan Lenaerts, Sam Levis,
L. Ruby Leung, William Lipscomb, Jon Pelletier, Daniel M. Ricciuto, Ben Sanderson, Jacquelyn Shuman, Andrew Slater, Zachary Subin, Jinyun Tang, Ahmed Tawfik, Quinn Thomas, Simone Tilmes, Francis Vitt, Xubin Zeng**

The National Center for Atmospheric Research (NCAR) is operated by the nonprofit University Corporation for Atmospheric Research (UCAR) under the sponsorship of the National Science Foundation. Any opinions, findings, conclusions, or recommendations expressed in this publication are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.

National Center for Atmospheric Research
P. O. Box 3000, Boulder, Colorado 80307-300

**LIST OF FIGURES**

- :numref:`Figure Land processes` Land biogeophysical, biogeochemical, and landscape processes simulated by CLM (adapted from :ref:`Lawrence et al. (2011)<Lawrenceetal2011>` for CLM5.0).

- :numref:`Figure CLM subgrid hierarchy` Configuration of the CLM subgrid hierarchy.

- :numref:`Figure Radiation Schematic` Schematic diagram of (a) direct beam radiation, (b) diffuse solar radiation, and (c) longwave radiation absorbed, transmitted, and reflected by vegetation and ground.

- :numref:`Figure Schematic diagram of sensible heat fluxes` Schematic diagram of sensible heat fluxes for (a) non-vegetated surfaces and (b) vegetated surfaces.

- :numref:`Figure Schematic diagram of latent heat fluxes` Schematic diagram of water vapor fluxes for (a) non-vegetated surfaces and (b) vegetated surfaces.

- :numref:`Figure Soil Temperature Schematic`. Schematic diagram of numerical scheme used to solve for soil temperature.

- :numref:`Figure Hydrologic processes` Hydrologic processes represented in CLM.

- :numref:`Figure Water flux schematic` Schematic diagram of numerical scheme used to solve for soil water fluxes.

- :numref:`Figure three layer snow pack` Example of three layer snow pack (snl=-3).

- :numref:`Figure MOSART conceptual diagram` MOSART conceptual diagram.

- :numref:`Figure Schematic representation of the urban landunit` Schematic representation of the urban land unit.

- :numref:`Figure Schematic of urban and atmospheric model coupling` Schematic of urban and atmospheric model coupling.

- :numref:`Figure Schematic of THESIS urban properties tool` Schematic of THESIS urban properties tool.

- :numref:`Figure Vegetation fluxes and pools` Vegetation fluxes and pools.

- :numref:`Figure annual phenology cycle` Example of annual phenology cycle for seasonal deciduous.

- :numref:`Figure Schematic of decomposition model in CLM` Schematic of decomposition model in CLM.

- :numref:`Figure Pool structure` Pool structure, transitions, respired fractions, and turnover times for the 2 alternate soil decomposition models included in CLM.

- :numref:`Figure Biological nitrogen fixation` Biological nitrogen fixation as a function of annual net primary production.

- :numref:`Figure Methane Schematic` Schematic representation of biological and physical processes integrated in CLM that affect the net CH4 surface flux.

- :numref:`Figure Schematic of land cover change` Schematic of land cover change impacts on CLM carbon pools and fluxes.

- :numref:`Figure Schematic of translation of annual LUH2 land units` Schematic of translation of annual UNH land units to CLM plant functional types.

**LIST OF TABLES**

- :numref:`Table Plant functional types` Plant functional types

- :numref:`Table Plant functional type canopy top and bottom heights` Plant functional type canopy top and bottom heights

- :numref:`Table Soil layer structure` Soil layer structure

- :numref:`Table Atmospheric input to land model` Atmospheric input to land model

- :numref:`Table Land model output to atmospheric model` Land model output to atmospheric model

- :numref:`Table Surface data required for CLM and their base spatial resolution` Surface data required for CLM and their base spatial resolution

- :numref:`Table Physical constants` Physical constants

- :numref:`Table Plant functional type optical properties` Plant functional type optical properties

- :numref:`Table Intercepted snow optical properties` Intercepted snow optical properties

- :numref:`Table Dry and saturated soil albedos` Dry and saturated soil albedos

- :numref:`Table Spectral bands and weights used for snow radiative transfer` Spectral bands and weights used for snow radiative transfer

- :numref:`Table Single-scatter albedo values used for snowpack impurities and ice` Single-scatter albedo values used for snowpack impurities and ice

- :numref:`Table Mass extinction values` Mass extinction values (m2 kg-1) used for snowpack impurities and ice.

- :numref:`Table Asymmetry scattering parameters used for snowpack impurities and ice` Asymmetry scattering parameters used for snowpack impurities and ice.

- :numref:`Table Orbital parameters` Orbital parameters

- :numref:`Table Plant functional type aerodynamic parameters` Plant functional type aerodynamic parameters

- :numref:`Table Coefficients for saturation vapor pressure` Coefficients for e\ :sub:`sat`\ :sup:`T`

- :numref:`Table Coefficients for derivative of esat` Coefficients for the derivative of e\ :sub:`sat`\ :sup:`T`

- :numref:`Table Meltwater scavenging` Meltwater scavenging efficiency for particles within snow

- :numref:`Table snow layer thickness` Minimum and maximum thickness of snow layers (m)

- :numref:`Table Plant functional type (PFT) stomatal conductance parameters` Plant functional type (PFT) stomatal conductance parameters.

- :numref:`Table Temperature dependence parameters for C3 photosynthesis` Temperature dependence parameters for C3 photosynthesis.

- :numref:`Table Plant functional type root distribution parameters` Plant functional type root distribution parameters.

- :numref:`Table MOSART Parameters` List of parameters in the global hydrography dataset.

- :numref:`Table Allocation and CN ratio parameters` Allocation and carbon:nitrogen ratio parameters

- :numref:`Table Decomposition rate constants` Decomposition rate constants for litter and SOM pools, C:N ratios, and acceleration parameters for the CLM-CN decomposition pool structure.

- :numref:`Table Respiration fractions for litter and SOM pools` Respiration fractions for litter and SOM pools

- :numref:`Table Turnover times` Turnover times, C:N ratios, and acceleration parameters for the Century-based decomposition cascade.

- :numref:`Table Respiration fractions for Century-based structure` Respiration fractions for litter and SOM pools for Century-based structure

- :numref:`Table PFT-specific combustion completeness and fire mortality factors` PFT-specific combustion completeness and fire mortality factors.

- :numref:`Table Methane Parameter descriptions`  Parameter descriptions and sensitivity analysis ranges applied in the methane model.

- :numref:`Table Temperature dependence of aqueous and gaseous diffusion` Temperature dependence of aqueous and gaseous diffusion coefficients for CH4 and O2.

- :numref:`Table Crop plant functional types` Crop plant functional types (PFTs).

- :numref:`Table Crop phenology parameters`  Crop phenology and morphology parameters.

- :numref:`Table Crop allocation parameters` Crop allocation parameters.

- :numref:`Table Dust Mass fraction` Mass fraction m\ :sub:`i` , mass median diameter :sub:`v, i` , and geometric standard deviation :sub:`g, i` , per dust source mode i

- :numref:`Table Dust Minimum and maximum particle diameters` Minimum and maximum particle diameters in each dust transport bin j

**ACKNOWLEDGEMENTS**

The authors would like to acknowledge the substantial contributions of the following members of the Land Model and Biogeochemistry Working Groups to the development of the Community Land Model since its inception in 1996: Benjamin Andre, Ian Baker, Michael Barlage, Mike Bosilovich, Marcia Branstetter, Tony Craig, Aiguo Dai, Yongjiu Dai, Mark Decker, Scott Denning, Robert Dickinson, Paul Dirmeyer, Jared Entin, Jay Famiglietti, Johannes Feddema, Mark Flanner, Jon Foley, Andrew Fox, Inez Fung, David Gochis, Alex Guenther, Tim Hoar, Forrest Hoffman, Paul Houser, Trish Jackson, Brian Kauffman, Silvia Kloster, Natalie Mahowald, Jiafu Mao, Lei Meng, Sheri Michelson, Guo-Yue Niu, Adam Phillips, Taotao Qian, Jon Radakovich, James Randerson, Nan Rosenbloom, Steve Running, Koichi Sakaguchi, Adam Schlosser, Andrew Slater, Reto Stöckli, Ying Sun, Quinn Thomas, Peter Thornton, Mariana Vertenstein, Nicholas Viovy, Aihui Wang, Guiling Wang, Zong-Liang Yang, Charlie Zender, Xiaodong Zeng, and Xubin Zeng.

.. _rst_Introduction:

Introduction
=================

The purpose of this document is to fully describe the biogeophysical and biogeochemical parameterizations and numerical implementation of version 5.0 of the Community Land Model (CLM5.0). Scientific justification and evaluation of these parameterizations can be found in the referenced scientific papers (:ref:`rst_References`). This document and the CLM5.0 User's Guide together provide the user with the scientific description and operating instructions for CLM.

Model History
---------------

Inception of CLM
^^^^^^^^^^^^^^^^^^^^^^

The early development of the Community Land Model can be described as the merging of a community-developed land model focusing on biogeophysics and a concurrent effort at NCAR to expand the NCAR Land Surface Model (NCAR LSM, :ref:`Bonan 1996<Bonan1996>`) to include the carbon cycle, vegetation dynamics, and river routing. The concept of a community-developed land component of the Community Climate System Model (CCSM) was initially proposed at the CCSM Land Model Working Group (LMWG) meeting in February 1996. Initial software specifications and development focused on evaluating the best features of three existing land models: the NCAR LSM (:ref:`Bonan 1996, 1998<Bonan1996>`) used in the Community Climate Model (CCM3) and the initial version of CCSM; the Institute of Atmospheric Physics, Chinese Academy of Sciences land model (IAP94) (:ref:`Dai and Zeng 1997<DaiZeng1997>`); and the Biosphere-Atmosphere Transfer Scheme (BATS) (:ref:`Dickinson et al. 1993<Dickinsonetal1993>`) used with CCM2. A scientific steering committee was formed to review the initial specifications of the design provided by Robert Dickinson, Gordon Bonan, Xubin Zeng, and Yongjiu Dai and to facilitate further development. Steering committee members were selected so as to provide guidance and expertise in disciplines not generally well-represented in land surface models (e.g., carbon cycling, ecological modeling, hydrology, and river routing) and included scientists from NCAR, the university community, and government laboratories (R. Dickinson, G. Bonan, X. Zeng, Paul Dirmeyer, Jay Famiglietti, Jon Foley, and Paul Houser).

The specifications for the new model, designated the Common Land Model, were discussed and agreed upon at the June 1998 CCSM Workshop LMWG meeting. An initial code was developed by Y. Dai and was examined in March 1999 by Mike Bosilovich, P. Dirmeyer, and P. Houser. At this point an extensive period of code testing was initiated. Keith Oleson, Y. Dai, Adam Schlosser, and P. Houser presented preliminary results of offline 1-dimensional testing at the June 1999 CCSM Workshop LMWG meeting. Results from more extensive offline testing at plot, catchment, and large scale (up to global) were presented by Y. Dai, A. Schlosser, K. Oleson, M. Bosilovich, Zong-Liang Yang, Ian Baker, P. Houser, and P. Dirmeyer at the LMWG meeting hosted by COLA (Center for Ocean-Land-Atmosphere Studies) in November 1999. Field data used for validation included sites adopted by the Project for Intercomparison of Land-surface Parameterization Schemes (:ref:`Henderson-Sellers et al. 1993<Henderson-Sellersetal1993>`) (Cabauw, Valdai, Red-Arkansas river basin) and others [FIFE (:ref:`Sellers et al. 1988<Sellersetal1988>`), BOREAS :ref:`(Sellers et al. 1995<Sellersetal1995>`), HAPEX-MOBILHY (:ref:`André et al. 1986<Andréetal1986>`), ABRACOS (:ref:`Gash et al. 1996<Gashetal1996>`), Sonoran Desert (:ref:`Unland et al. 1996<Unlandetal1996>`), GSWP (:ref:`Dirmeyer et al. 1999<Dirmeyeretal1999>`)]. Y. Dai also presented results from a preliminary coupling of the Common Land Model to CCM3, indicating that the land model could be successfully coupled to a climate model.

Results of coupled simulations using CCM3 and the Common Land Model were presented by X. Zeng at the June 2000 CCSM Workshop LMWG meeting. Comparisons with the NCAR LSM and observations indicated major improvements to the seasonality of runoff, substantial reduction of a summer cold bias, and snow depth. Some deficiencies related to runoff and albedo were noted, however, that were subsequently addressed. Z.-L. Yang and I. Baker demonstrated improvements in the simulation of snow and soil temperatures. Sam Levis reported on efforts to incorporate a river routing model to deliver runoff to the ocean model in CCSM. Soon after the workshop, the code was delivered to NCAR for implementation into the CCSM framework. Documentation for the Common Land Model is provided by :ref:`Dai et al. (2001)<Daietal2001>` while the coupling with CCM3 is described in :ref:`Zeng et al. (2002)<Zengetal2002>`. The model was introduced to the modeling community in :ref:`Dai et al. (2003)<Daietal2003>`.

CLM2
^^^^^^^^^^

Concurrent with the development of the Common Land Model, the NCAR LSM was undergoing further development at NCAR in the areas of carbon cycling, vegetation dynamics, and river routing. The preservation of these advancements necessitated several modifications to the Common Land Model. The biome-type land cover classification scheme was replaced with a plant functional type (PFT) representation with the specification of PFTs and leaf area index from satellite data (:ref:`Oleson and Bonan 2000<OlesonBonan2000>`; :ref:`Bonan et al. 2002a, b<Bonanetal2002a>`). This also required modifications to parameterizations for vegetation albedo and vertical burying of vegetation by snow. Changes were made to canopy scaling, leaf physiology, and soil water limitations on photosynthesis to resolve deficiencies indicated by the coupling to a dynamic vegetation model. Vertical heterogeneity in soil texture was implemented to improve coupling with a dust emission model. A river routing model was incorporated to improve the fresh water balance over oceans. Numerous modest changes were made to the parameterizations to conform to the strict energy and water balance requirements of CCSM. Further substantial software development was also required to meet coding standards. The resulting model was adopted in May 2002 as the Community Land Model (CLM2) for use with the Community Atmosphere Model (CAM2, the successor to CCM3) and version 2 of the Community Climate System Model (CCSM2).

K. Oleson reported on initial results from a coupling of CCM3 with CLM2 at the June 2001 CCSM Workshop LMWG meeting. Generally, the CLM2 preserved most of the improvements seen in the Common Land Model, particularly with respect to surface air temperature, runoff, and snow. These simulations are documented in :ref:`Bonan et al. (2002a)<Bonanetal2002a>`. Further small improvements to the biogeophysical parameterizations, ongoing software development, and extensive analysis and validation within CAM2 and CCSM2 culminated in the release of CLM2 to the community in May 2002.

Following this release, Peter Thornton implemented changes to the model structure required to represent carbon and nitrogen cycling in the model. This involved changing data structures from a single vector of spatially independent sub-grid patches to one that recognizes three hierarchical scales within a model grid cell: land unit, snow/soil column, and PFT. Furthermore, as an option, the model can be configured so that PFTs can share a single soil column and thus "compete" for water. This version of the model (CLM2.1) was released to the community in February 2003. CLM2.1, without the compete option turned on, produced only round off level changes when compared to CLM2.

CLM3
^^^^^^^^^^

CLM3 implemented further software improvements related to performance and model output, a re-writing of the code to support vector-based computational platforms, and improvements in biogeophysical parameterizations to correct deficiencies in the coupled model climate. Of these parameterization improvements, two were shown to have a noticeable impact on simulated climate. A variable aerodynamic resistance for heat/moisture transfer from ground to canopy air that depends on canopy density was implemented. This reduced unrealistically high surface temperatures in semi-arid regions. The second improvement added stability corrections to the diagnostic 2-m air temperature calculation which reduced biases in this temperature. Competition between PFTs for water, in which PFTs share a single soil column, is the default mode of operation in this model version. CLM3 was released to the community in June 2004. :ref:`Dickinson et al. (2006)<Dickinsonetal2006>` describe the climate statistics of CLM3 when coupled to CCSM3.0. :ref:`Hack et al. (2006)<Hacketal2006>` provide an analysis of selected features of the land hydrological cycle. :ref:`Lawrence et al. (2007)<Lawrenceetal2007>` examine the impact of changes in CLM3 hydrological parameterizations on partitioning of evapotranspiration (ET) and its effect on the timescales of ET response to precipitation events, interseasonal soil moisture storage, soil moisture memory, and land-atmosphere coupling. :ref:`Qian et al. (2006)<Qianetal2006>` evaluate CLM3's performance in simulating soil moisture content, runoff, and river discharge when forced by observed precipitation, temperature and other atmospheric data.

CLM3.5
^^^^^^^^^^^^

Although the simulation of land surface climate by CLM3 was in many ways adequate, most of the unsatisfactory aspects of the simulated climate noted by the above studies could be traced directly to deficiencies in simulation of the hydrological cycle. In 2004, a project was initiated to improve the hydrology in CLM3 as part of the development of CLM version 3.5. A selected set of promising approaches to alleviating the hydrologic biases in CLM3 were tested and implemented. These included new surface datasets based on Moderate Resolution Imaging Spectroradiometer (MODIS) products, new parameterizations for canopy integration, canopy interception, frozen soil, soil water availability, and soil evaporation, a TOPMODEL-based model for surface and subsurface runoff, a groundwater model for determining water table depth, and the introduction of a factor to simulate nitrogen limitation on plant productivity. :ref:`Oleson et al. (2008a)<Olesonetal2008a>` show that CLM3.5 exhibits significant improvements over CLM3 in its partitioning of global ET which result in wetter soils, less plant water stress, increased transpiration and photosynthesis, and an improved annual cycle of total water storage. Phase and amplitude of the runoff annual cycle is generally improved. Dramatic improvements in vegetation biogeography result when CLM3.5 is coupled to a dynamic global vegetation model. :ref:`Stöckli et al. (2008)<Stocklietal2008>` examine the performance of CLM3.5 at local scales by making use of a network of long-term ground-based ecosystem observations [FLUXNET (:ref:`Baldocchi et al. 2001<Baldocchietal2001>`)]. Data from 15 FLUXNET sites were used to demonstrate significantly improved soil hydrology and energy partitioning in CLM3.5. CLM3.5 was released to the community in May, 2007.

CLM4
^^^^^^^^^^

The motivation for the next version of the model, CLM4, was to incorporate several recent scientific advances in the understanding and representation of land surface processes, expand model capabilities, and improve surface and atmospheric forcing datasets (:ref:`Lawrence et al. 2011<Lawrenceetal2011>`). Included in the first category are more sophisticated representations of soil hydrology and snow processes. In particular, new treatments of soil column-groundwater interactions, soil evaporation, aerodynamic parameters for sparse/dense canopies, vertical burial of vegetation by snow, snow cover fraction and aging, black carbon and dust deposition, and vertical distribution of solar energy for snow were implemented. Major new capabilities in the model include a representation of the carbon-nitrogen cycle (CLM4CN, see next paragraph for additional information), the ability to model land cover change in a transient mode, inclusion of organic soil and deep soil into the existing mineral soil treatment to enable more realistic modeling of permafrost, an urban canyon model to contrast rural and urban energy balance and climate (CLMU), and an updated biogenic volatile organic compounds (BVOC) model. Other modifications of note include refinement of the global PFT, wetland, and lake distributions, more realistic optical properties for grasslands and croplands, and an improved diurnal cycle and spectral distribution of incoming solar radiation to force the model in land-only mode.

Many of the ideas incorporated into the carbon and nitrogen cycle component of CLM4 derive from the earlier development of the land-only ecosystem process model Biome-BGC (Biome BioGeochemical Cycles), originating at the Numerical Terradynamic Simulation Group (NTSG) at the University of Montana, under the guidance of Prof. Steven Running. Biome-BGC itself is an extension of an earlier model, Forest-BGC (:ref:`Running and Coughlan, 1988<RunningCoughlan1988>`; :ref:`Running and Gower, 1991<RunningGower1991>`), which simulates water, carbon, and, to a limited extent, nitrogen fluxes for forest ecosystems. Forest-BGC was designed to be driven by remote sensing inputs of vegetation structure, and so used a diagnostic (prescribed) leaf area index, or, in the case of the dynamic allocation version of the model (:ref:`Running and Gower, 1991<RunningGower1991>`), prescribed maximum leaf area index.

Biome-BGC expanded on the Forest-BGC logic by introducing a more mechanistic calculation of leaf and canopy scale photosynthesis (:ref:`Hunt and Running, 1992<Huntrunning1992>`), and extending the physiological parameterizations to include multiple woody and non-woody vegetation types (:ref:`Hunt et al. 1996<Huntetal1996>`; :ref:`Running and Hunt, 1993<RunningHunt1993>`). Later versions of Biome-BGC introduced more mechanistic descriptions of belowground carbon and nitrogen cycles, nitrogen controls on photosynthesis and decomposition, sunlit and shaded canopies, vertical gradient in leaf morphology, and explicit treatment of fire and harvest disturbance and regrowth dynamics (:ref:`Kimball et al. 1997<Kimballetal1997>`; :ref:`Thornton, 1998<Thornton1998>`; :ref:`Thornton et al. 2002<Thorntonetal2002>`; :ref:`White et al. 2000<Whiteetal2000>`). Biome-BGC version 4.1.2 (:ref:`Thornton et al. 2002<Thorntonetal2002>`) provided a point of departure for integrating new biogeochemistry components into CLM4.

CLM4 was released to the community in June, 2010 along with the Community Climate System Model version 4 (CCSM4). CLM4 is used in CCSM4, CESM1, CESM1.1, and remains available as the default land component model option for coupled simulations in CESM1.2.

CLM4.5
^^^^^^^^^^^^

The motivations for the development of CLM4.5 were similar to those for CLM4: incorporate several recent scientific advances in the understanding and representation of land surface processes, expand model capabilities, and improve surface and atmospheric forcing datasets.

Specifically, several parameterizations were revised to reflect new scientific understanding and in an attempt to reduce biases identified in CLM4 simulations including low soil carbon stocks especially in the Arctic, excessive tropical GPP and unrealistically low Arctic GPP, a dry soil bias in Arctic soils, unrealistically high LAI in the tropics, a transient 20\ :math:`{}^{th}` century carbon response that was inconsistent with observational estimates, and several other more minor problems or biases.

The main modifications include updates to canopy processes including a revised canopy radiation scheme and canopy scaling of leaf processes, co-limitations on photosynthesis, revisions to photosynthetic parameters (:ref:`Bonan et al. 2011<Bonanetal2011>`), temperature acclimation of photosynthesis, and improved stability of the iterative solution in the photosynthesis and stomatal conductance model (:ref:`Sun et al. 2012<Sunetal2012>`). Hydrology updates included modifications such that hydraulic properties of frozen soils are determined by liquid water content only rather than total water content and the introduction of an ice impedance function, and other corrections that increase the consistency between soil water state and water table position and allow for a perched water table above icy permafrost ground (:ref:`Swenson et al. 2012<Swensonetal2012>`). A new snow cover fraction parameterization is incorporated that reflects the hysteresis in fractional snow cover for a given snow depth between accumulation and melt phases (:ref:`Swenson and Lawrence, 2012<SwensonLawrence2012>`). The lake model in CLM4 was replaced with a completely revised and more realistic lake model (:ref:`Subin et al. 2012a<Subinetal2012a>`). A surface water store was introduced, replacing the wetland land unit and permitting prognostic wetland distribution modeling. The surface energy fluxes are calculated separately (:ref:`Swenson and Lawrence, 2012<SwensonLawrence2012>`) for snow-covered, water-covered, and snow/water-free portions of vegetated and crop land units, and snow-covered and snow-free portions of glacier land units. Globally constant river flow velocity is replaced with variable flow velocity based on mean grid cell slope. A vertically resolved soil biogeochemistry scheme is introduced with base decomposition rates modified by soil temperature, water, and oxygen limitations and also including vertical mixing of soil carbon and nitrogen due to bioturbation, cryoturbation, and diffusion (:ref:`Koven et al. 2013<Kovenetal2013>`). The litter and soil carbon and nitrogen pool structure as well as nitrification and denitrification that were modified based on the Century model. Biological fixation was revised to distribute fixation more realistically over the year (:ref:`Koven et al. 2013<Kovenetal2013>`). The fire model was replaced with a model that includes representations of natural and anthropogenic triggers and suppression as well as agricultural, deforestation, and peat fires (:ref:`Li et al. 2012a,b<Lietal2012a>`; :ref:`Li et al. 2013a<Lietal2013a>`). The biogenic volatile organic compounds model is updated to MEGAN2.1 (:ref:`Guenther et al. 2012<Guentheretal2012>`).

Additions to the model include a methane production, oxidation, and emissions model (:ref:`Riley et al. 2011a<Rileyetal2011a>`) and an extension of the crop model to include interactive fertilization, organ pools (:ref:`Drewniak et al. 2013<Drewniaketal2013>`), and irrigation (:ref:`Sacks et al. 2009<Sacksetal2009>`). Elements of the Variable Infiltration Capacity (VIC) model are included as an alternative optional runoff generation scheme (:ref:`Li et al. 2011<Lietal2011>`). There is also an option to run with a multilayer canopy (:ref:`Bonan et al. 2012<Bonanetal2012>`). Multiple urban density classes, rather than the single dominant urban density class used in CLM4, are modeled in the urban land unit. Carbon (:math:`{}^{13}`\ C and :math:`{}^{14}`\ C) isotopes are enabled (:ref:`Koven et al. 2013<Kovenetal2013>`). Minor changes include a switch of the C3 Arctic grass and shrub phenology from stress deciduous to seasonal deciduous and a change in the glacier bare ice albedo to better reflect recent estimates. Finally, the carbon and nitrogen cycle spinup is accelerated and streamlined with a revised spinup method, though the spinup timescale remains long.

Finally, the predominantly low resolution input data for provided with CLM4 to create CLM4 surface datasets is replaced with newer and higher resolution input datasets where possible (see section :numref:`Surface Data` for details). The default meteorological forcing dataset provided with CLM4 (:ref:`Qian et al. 2006)<Qianetal2006>` is replaced with the 1901-2010 CRUNCEP forcing dataset (see Chapter :numref:`rst_Land-Only Mode`) for CLM4.5, though users can also still use the :ref:`Qian et al. (2006)<Qianetal2006>` dataset or other alternative forcing datasets.

CLM4.5 was released to the community in June 2013 along with the Community Earth System Model version 1.2 (CESM1.2).

CLM5.0
^^^^^^^^^^^^

Developments for CLM5.0 build on the progress made in CLM4.5. Most major components of the model have been updated with particularly notable changes made to soil and plant hydrology, snow density, river modeling, carbon and nitrogen cycling and coupling, and crop modeling. Much of the focus of development centered on a push towards more mechanistic treatment of key processes, in addition to more comprehensive and explicit representation of land use and land-cover change. Prior versions of CLM included relatively few options for physics parameterizations or structure. In CLM5, where new parameterizations or model decisions were made, in most cases, the CLM4.5 parameterization was maintained so that users could switch back and forth between different parameterizations via namelist control where appropriate or desirable. Throughout the CLM5 Technical Descpription, in general only the default parameterization for any given process is described. Readers are referred to the CLM4.5 or CLM4 Technical Descriptions for detailed descriptions of non-default parameterizations.

The hydrology updates include the introduction of a dry surface layer-based soil evaporation resistance parameterization :ref:`(Swenson and Lawrence, 2014)<SwensonLawrence2014>` and a revised canopy interception parameterization. Canopy interception is now divided into liquid and solid phases, with the intercepted snow subject to unloading events due to wind or above-freezing temperatures. The snow-covered fraction of the canopy is used within the canopy radiation and surface albedo calculation. Instead of applying a spatially uniform soil thickness, soil thickness can vary in space :ref:`(Brunke et al. 2016<Brunkeetal2016>` and :ref:`Swenson and Lawrence, 2015)<SwensonLawrence2015>` and is set to values within a range of 0.4m to 8.5m depth, derived from a spatially explicit soil thickness data product :ref:`(Pelletier et al., 2016)<Pelletieretal2016>`. The explicit treatment of soil thickness allows for the deprecation of the unconfined aquifer parameterization used in CLM4.5, which is replaced with a zero flux boundary condition and explicit modeling of both the saturated and unsaturated zones. The default model soil layer resolution is increased, especially within the top 3m, to more explicitly represent active layer thickness within the permafrost zone. Rooting profiles were used inconsistently in CLM4.5 with :ref:`Zeng (2001)<Zeng2001>` profiles used for water and :ref:`Jackson et al. (1996)<Jacksonetal1996>` profiles used for carbon inputs. For CLM5, the Jackson et al. (1996) rooting profiles are used for both water and carbon. Roots are deepened for the broadleaf evergreen tropical tree and broadleaf deciduous tropical tree types. Finally, an adaptive time-stepping solution to the Richard's equation is introduced, which improves the accuracy and stability of the numerical soil water solution. The River Transport Model (RTM) is replaced with the Model for Scale Adaptive River Transport (MOSART, :ref:`Li et al., 2013b)<Lietal2013b>` in which surface runoff is routed across hillslopes and then discharged along with subsurface runoff into a tributary subnetwork before entering the main channel.

Several changes are included that are mainly targeted at improving the simulation of surface mass balance over ice sheets. The fresh snow density parameterization is updated to more realistically capture temperature effects and to additionally account for wind effects on new snow density :ref:`(van Kampenhout et al., 2017)<vanKampenhoutetal2017>`. The maximum number of snow layers and snow amount is increased from 5 layers and 1m snow water equivalent to 12 layers and 10m snow water equivalent to allow for the formation of firn in regions of persistent snow-cover (e.g., glaciers and ice sheets) :ref:`(van Kampenhout et al., 2017)<vanKampenhoutetal2017>`. The CISM2 ice sheet model is included for Greenland by default. The ice sheet does not evolve for typical configurations, but ice sheet evolution can be turned on by choosing an appropriate compset. The introduction in CLM5 of the capability to dynamically adjust landunit weights means that a glacier can initiate, grow, shrink, or disappear during a simulation when ice evolution is active. That is, there are two-way feedbacks between CLM and CISM. Multiple elevation classes (10 elevation classes by default) and associated temperature, rain/snow partitioning, and downwelling longwave downscaling are used for glacier landunits to account for the strong topographic elevation heterogeneity over glaciers and ice sheets.

A plant hydraulic stress routine is introduced which explicitly models water transport through the vegetation according to a simple hydraulic framework (Kennedy et al., to be submitted). The water supply equations are used to solve for vegetation water potential forced by transpiration demand and a set of layer-by-layer soil water potentials. Stomatal conductance, therefore, is a function of prognostic leaf water potential. Water stress is calculated as the ratio of attenuated stomatal conductance to maximum stomatal conductance. An emergent feature of the plant hydraulics is soil hydraulic redistribution. In CLM5, maximum stomatal conductance is obtained from the Medlyn conductance model :ref:`(Medlyn et al., 2011)<Medlynetal2011>`, rather than the Ball-Berry stomatal conductance model that was utilized in CLM4.5 and prior versions of the model. The Medlyn stomatal conductance model is preferred mainly for it's more realistic behavior at low humidity levels :ref:`(Rogers et al., 2017)<Rogersetal2017>`. The stress deciduous vegetation phenology trigger is augmented with a antecedent precipitation requirement :ref:`(Dahlin et al. 2015)<Dahlinetal2015>`.

Plant nutrient dynamics are substantially updated to resolve several deficiencies with the CLM4 and CLM4.5 nutrient cycling representation. The Fixation and Update of Nitrogen (FUN) model based on the work of :ref:`Fisher et al. (2010)<Fisheretal2010>`, :ref:`Brzostek et al. (2014)<Brzosteketal2014>`, and :ref:`Shi et al. (2016)<Shietal2016>` is incorporated. The concept of FUN is that in most cases, N uptake requires the expenditure of energy in the form of carbon, and further, that there are numerous potential sources of N in the environment which a plant may exchange for carbon. The ratio of carbon expended to N acquired is therefore the cost, or exchange rate, of N acquisition. FUN calculates the rate of symbiotic N fixation, with this N passed straight to the plant, not the mineral N pool. Separately, CLM5 also calculates rates of symbiotic (or free living) N fixation as a function of evapotranspiration (:ref:`Cleveland et al. 1999 <Clevelandetal1999>`), which is added to the soil inorganic ammonium (NH\ :sub:`4`\ :sup:`+`) pool. The static plant carbon:nitrogen (C:N) ratios utilized in CLM4 and CLM4.5 are replaced with variable plant C:N ratios which allows plants to adjust their C:N ratio, and therefore their leaf nitrogen content, with the cost of N uptake :ref:`(Ghimire et al. 2016)<Ghimireetal2016>`. The implementation of a flexible C:N ratio means that the model no longer relies on instantaneous downregulation of potential photosynthesis rates based on soil mineral nitrogen availability to represent nutrient limitation. Furthermore, stomatal conductance is now based on the N-limited photosynthesis rather than on potential photosynthesis. Finally, the Leaf Use of Nitrogen for Assimilation (LUNA, :ref:`Xu et al., 2012<Xuetal2012>` and :ref:`Ali et al., 2016)<Alietal2016>` model is incorporated. The LUNA model calculates photosynthetic capacity based on optimization of the use of leaf nitrogen under different environmental conditions such that light capture, carboxylation, and respiration are co-limiting.

CLM5 applies a fixed allocation scheme for woody vegetation. The decision to use a fixed allocation scheme in CLM5, rather than a dynamic NPP-based allocation scheme, as was used in CLM4 and CLM4.5, was driven by the fact that observations indicate that biomass saturates with increasing productivity, in contrast to the behavior in CLM4 and CLM4.5 where biomass continuously increases with increasing productivity (:ref:`Negron-Juarez et al., 2015<NegronJuarezetal2015>`). Soil carbon decomposition processes are unchanged in CLM5, but a new metric for apparent soil carbon turnover times (:ref:`Koven et al., 2017 <Kovenetal2017>`) suggested parameter changes that produce a weak intrinsic depth limitation on soil carbon turnover rates (rather than the strong depth limitaiton in CLM4.5) and that the thresholds for soil moisture limitation on soil carbon turnover rates in dry soils should be set at a wetter soil moisture level than that used in CLM4.5.

Representation of human management of the land (agriculture, wood harvest) is augmented in several ways. The CLM4.5 crop model is extended to operate globally through the addition of rice, sugarcane, tropical varieties of corn and soybean :ref:`(Badger and Dirmeyer, 2015<BadgerandDirmeyer2015>` and :ref:`Levis et al., 2016)<Levisetal2016>`, and perennial bioenergy crops :ref:`(Cheng et al., 2019)<Chengetal2019>`. These crop types are added to the existing temperate corn, temperate soybean, spring wheat, and cotton crop types. Fertilization rates and irrigation equipped area updated annually based on crop type and geographic region through an input dataset. The irrigation trigger is updated. Additional minor changes include crop phenological triggers that vary by latitude for selected crop types, grain C and N is now removed at harvest to a 1-year product pool with the carbon for the next season's crop seed removed from the grain carbon at harvest. A fraction of leaf/livestem C and N from bioenergy crops is removed at harvest to the biofuel feedstock pools and added to the 1-year product pool. Through the introduction of the capability to dynamically adjust landunit weights during a simulation, the crop model can now be run coincidentally with prescribed land use, which significantly expands the capabilities of the model. Mass-based rather than area-based wood harvest is applied. Several heat stress indices for both urban and rural areas are calculated and output by default :ref:`(Buzan et al., 2015)<Buzanetal2015>`. A more sophisticated and realistic building space heating and air conditioning submodel that prognoses interior building air temperature and includes more realistic space heating and air conditioning wasteheat factors is incorporated.

The fire model is the same as utilized in CLM4.5 except that a modified scheme is used to estimate the dependence of fire occurrence and spread on fuel wetness for non-peat fires outside cropland and tropical closed forests :ref:`(Li and Lawrence, 2017)<LiLawrence2017>` and the dependence of agricultural fires on fuel load is removed.

Included with the release of CLM5.0 is a functionally supported version of the Functionally-Assembled Terrestrial Ecosystem Simulator (FATES, :ref:`Fisher et al., 2015)<Fisheretal2015>`. A major motivation of FATES is to allow the prediction of biome boundaries directly from plant physiological traits via their competitive interactions. FATES is a cohort model of vegetation competition and co-existence, allowing a representation of the biosphere which accounts for the division of the land surface into successional stages, and for competition for light between height structured cohorts of representative trees of various plant functional types. FATES is not active by default in CLM5.0.

Note that the classical dynamic global vegetation model (CLM-DGVM) that has been available within CLM4 and CLM4.5 remains available, though it is largely untested. The technical description of the CLM-DGVM can be found within the CLM4.5 Technical Description (:ref:`Oleson et al. 2013)<Olesonetal2013>`.

During the course of the development of CLM5.0, it became clear that the increasing complexity of the model combined with the increasing number and range of model development projects required updates to the underlying CLM infrastructure. Many such software improvements are included in CLM5 including a partial transition to an object-oriented modular software structure. Many hard coded model parameters have been extracted into either the parameter file or the CLM namelist, which allows users to more readily calibrate the model for use at specific locations or to conduct parameter sensitivity studies. As part of the effort to increase the scientific utility of the code, in most instances older generation parameterizations (i.e., the parameterizations available in CLM4 or CLM4.5) are retained under namelist switches, allowing the user to revert to CLM4.5 from the same code base or to revert individual parameterizations where the old parameterizations are compatible with the new code. Finally, multiple vertical soil layer structures are defined and it is relatively easy to add additional structures.

Biogeophysical and Biogeochemical Processes
-----------------------------------------------

Biogeophysical and biogeochemical processes are simulated for each subgrid land unit, column, and plant functional type (PFT) independently and each subgrid unit maintains its own prognostic variables (see section :numref:`Surface Heterogeneity and Data Structure` for definitions of subgrid units). The same atmospheric forcing is used to force all subgrid units within a grid cell. The surface variables and fluxes required by the atmosphere are obtained by averaging the subgrid quantities weighted by their fractional areas. The processes simulated include (:numref:`Figure Land processes`):

#. Surface characterization including land type heterogeneity and ecosystem structure (Chapter :numref:`rst_Surface Characterization, Vertical Discretization, and Model Input Requirements`)

#. Absorption, reflection, and transmittance of solar radiation (Chapter :numref:`rst_Surface Albedos`, :numref:`rst_Radiative Fluxes`)

#. Absorption and emission of longwave radiation (Chapter :numref:`rst_Radiative Fluxes`)

#. Momentum, sensible heat (ground and canopy), and latent heat (ground evaporation, canopy evaporation, transpiration) fluxes (Chapter :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`)

#. Heat transfer in soil and snow including phase change (Chapter :numref:`rst_Soil and Snow Temperatures`)

#. Canopy hydrology (interception, throughfall, and drip) (Chapter :numref:`rst_Hydrology`)

#. Soil hydrology (surface runoff, infiltration, redistribution of water within the column, sub-surface drainage, groundwater) (Chapter :numref:`rst_Hydrology`)

#. Snow hydrology (snow accumulation and melt, compaction, water transfer between snow layers) (Chapter :numref:`rst_Snow Hydrology`)

#. Stomatal physiology, photosythetic capacity, and photosynthesis (Chapters :numref:`rst_Stomatal Resistance and Photosynthesis` and :numref:`rst_Photosynthetic Capacity`)

#. Plant hydraulics (Chapter :numref:`rst_Plant Hydraulics`)

#. Lake temperatures and fluxes (Chapter :numref:`rst_Lake Model`)

#. Glacier processes (Chapter :numref:`rst_Glaciers`)

#. River routing and river flow (Chapter :numref:`rst_River Transport Model (RTM)`)

#. Urban energy balance and climate (Chapter :numref:`rst_Urban Model (CLMU)`)

#. Vegetation carbon and nitrogen allocation (Chapter :numref:`rst_CN Allocation`)

#. Vegetation phenology (Chapter :numref:`rst_Vegetation Phenology and Turnover`)

#. Plant respiration (Chapter :numref:`rst_Plant Respiration`)

#. Soil and litter carbon decomposition (Chapter :numref:`rst_Decomposition`)

#. Fixation and uptake of nitrogen (Chapter :numref:`rst_FUN`)

#. External nitrogen cycling including deposition, denitrification, leaching, and losses due to fire (Chapter :numref:`rst_External Nitrogen Cycle`)

#. Plant mortality (Chapter :numref:`rst_Plant Mortality`)

#. Fire ignition, suppression, spread, and emissions, including natural, deforestation, and agricultural fire (Chapter :numref:`rst_Fire`)

#. Methane production, oxidation, and emissions (Chapter :numref:`rst_Methane Model`)

#. Crop dynamics, irrigation, and fertilization (Chapter :numref:`rst_Crops and Irrigation`)

#. Land cover and land use change including wood harvest (Chapter :numref:`rst_Transient Landcover Change`)

#. Biogenic volatile organic compound emissions (Chapter :numref:`rst_Biogenic Volatile Organic Compounds (BVOCs)`)

#. Dust mobilization and deposition (Chapter :numref:`rst_Dust Model`)

#. Carbon isotope fractionation (Chapter :numref:`rst_Carbon Isotopes`)

.. _Figure Land processes:

.. figure:: image1.png

 Land biogeophysical, biogeochemical, and landscape processes simulated by CLM (adapted from :ref:`Lawrence et al. (2011)<Lawrenceetal2011>` for CLM5.0).
