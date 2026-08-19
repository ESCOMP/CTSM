.. _rst_Urban Model (CLMU converted from docx):

**************************************
Urban Model (CLMU converted from docx)
**************************************

Introduction
===============

This technical note describes the physical parameterizations and numerical implementation of a Community Land Model Urban (CLMU) parameterization as coupled to version 4 of the Community Land Model (CLM4). CLM4 serves as the land surface model component of the Community Atmosphere Model (CAM) and the Community Climate System Model (CCSM). This note documents the global implementation of the urban model. Other model versions may exist for specific applications.

Chapters 1-5 constitute the description of the urban parameterization when coupled to or CCSM, while Chapter 6 describes processes that pertain specifically to the operation of the urban parameterization in offline mode (uncoupled to an atmospheric model). Chapter 7 describes efforts to evaluate the urban model. The model formulation and some quantitative and qualitative evaluation are also documented in Oleson et al. (2008a, 2008b). A heat island mitigation study using the model is presented in Oleson et al. (2010a). Note that CLMU and CLM4 have some parameterizations in common (e.g., snow and sub-surface hydrology). This technical note contains material duplicated from the CLM4 technical note (Oleson et al. 2010b) where appropriate. This is done so that users interested in just the urban model do not have to refer to the CLM4 technical note.

Model Overview
------------------

Motivation
~~~~~~~~~~~~~~~~

Land use and land cover change is increasingly being recognized as an important yet poorly quantified component of global climate change (Houghton et al. 2001). Land use/cover change mechanisms include both the transformation of natural land surfaces to those serving human needs (i.e., direct anthropogenic change) (e.g., the conversion of tropical forest to agriculture) as well as changes in land cover on longer time-scales due to biogeophysical feedbacks between the atmosphere and the land (i.e., indirect change) (Cramer et al. 2001, Foley et al. 2005). Global and regional models have been used extensively to investigate the effects of direct and indirect land use/cover change mechanisms on climate (Copeland et al. 1996, Stohlgren et al. 1998, Betts 2001, Eastman et al. 2001, Bounoua et al. 2002, Pielke et al. 2002, Fu 2003, Myhre and Myhre 2003, Narisma and Pitman 2003, Wang et al. 2003, Brovkin et al. 2004, Mathews et al. 2004, Feddema et al. 2005). However, all of these studies have focused on land use/cover related to changes in vegetation types. Urbanization, or the expansion of built-up areas, is an important yet less studied aspect of anthropogenic land use/cover change in climate science.

Although currently only about 1-3% of the global land surface is urbanized, the spatial extent and intensity of urban development is expected to increase dramatically in the future (Shepherd 2005, CIESIN et al. 2004). More than one-half of the world’s population currently lives in urban areas and in Europe, North America, and at least 80% of the population resides in urban areas (Elvidge et al. 2004). Policymakers and the public are most interested in the effects of climate change on people where they live. Because urban and non-urban areas may have different sensitivities to climate change, it is possible that the true climate change signal within urban areas may only be estimated if urban areas are explicitly modeled in climate change simulations (Best 2006). Indeed, the “footprint” of urbanization on climate can be detected from surface observations and satellite data (Changnon 1992, Kalnay and Cai 2003, Zhou et al. 2004, Jin et al. 2005). Changnon (1992) points out that the average urban warming over the last 100 years in certain regions is comparable to the increase in global surface temperature predicted by climate models for the next 100 years. Thus, it is important for developers of land surface models to begin to consider the parameterization of urban surfaces.

Urbanization now appropriates significant proportions of land in certain regions. For example, the expansion of service-based industries and conversion of farmland for housing in the Chicago area has increased the amount of developed land from about 800 square miles in 1973 to 1000 square miles in 1992 (Auch et al. 2004). A T85 resolution climate model grid cell (the resolution of the CCSM3 climate change simulations submitted for the IPCC AR4) encompassing the Chicago region represents about 7100 square miles, which suggests that this grid cell should be modeled as about 14% urbanized land. For mesoscale or regional models, where grid cells are on the order of a few kilometers, an urban area this size will occupy a significant number of grid cells that would otherwise be modeled as natural surfaces. The now common use of multiple “tiles” in models enables the co-existence of multiple surface types within a single gridcell. Thus, urban areas should and can be included in a global climate model (Best 2006).

Numerical modeling of the urban energy budget was first attempted nearly 40 years ago [see Brown (2000) for a comprehensive historical overview of modeling efforts]. However, until recently, most modern land surface models [i.e., second- or third-generation models (Sellers et al. 1997)] have not formally included urban parameterizations. Masson (2006) classifies urban parameterizations in three general categories: 1) empirical models, 2) vegetation models, with and without drag terms, adapted to include an urban canopy, 3) single-layer and multi-layer models that include a three-dimensional representation of the urban canopy. Empirical models (e.g., Oke and Cleugh 1987) rely on statistical relationships determined from observed data. As such, they are generally limited to the range of conditions experienced during the observation campaign. Vegetation models adapted for the urban canopy generally focus on modifying important surface parameters to better represent urban surfaces [e.g., surface albedo, roughness length, displacement height, surface emissivity, heat capacity, thermal conductivity (Taha 1999, Atkinson 2003, Best 2005].

These relatively simple approaches (i.e., categories 1 and 2 above) may arguably be justified based on the fact that detail in complex models may be lost when averaged to a coarse grid (Taha 1999). However, they may not have sufficient functionality to be suitable for inclusion in global climate models and may require the global derivation of parameters that are difficult to interpret physically [e.g., the surface type-dependent empirical coefficients for storage heat flux in the Objective Hysteresis Model (Grimmond et al. 1991)]. Furthermore, such approaches may not fully describe the fundamental processes that determine urban effects on climate (Piringer et al. 2002). For example, cities are known to have unique characteristics that cause them to be warmer than surrounding rural areas, an effect known as the urban heat island (Oke 1987). In the absence of anthropogenic heat flux, the urban heat island is thought to be greatest on clear, calm nights when local conditions generally dominate over synoptic. Candidate causes for this phenomenon include decreased surface longwave radiation loss and increased absorption of solar radiation because of canyon geometry, anthropogenic emissions of heat, reduction of evapotranspiration due to the replacement of vegetation with impervious surfaces, increased downwelling longwave radiation from the atmosphere due to pollution and warmer atmospheric temperatures, increased storage of sensible heat within urban materials, and reduced transfer of heat due to sheltering from buildings (Oke 1982, Oke 1987, Oke et al. 1991). Single-layer or multi-layer urban canopy models are likely needed to investigate the relative contribution of these factors to the heat island effect (Piringer et al. 2002). For example, specification of an urban albedo may provide no insight into the effects of the individual albedo of roofs, walls, and roads and the interaction of shortwave radiation between these surfaces that yields urban albedos that are typically lower than those of most rural sites. Similarly, assessments of the effectiveness of techniques proposed to ameliorate heat islands, such as “green roofs” or tree planting, require more detailed models.

On the other hand, the level of complexity in a model is limited by the availability of data that the model requires, the computational burden imposed, and difficulty in understanding the complex behavior of the model. Here, following recent developments in detailed urban parameterizations designed for mesoscale models (Masson 2000, Martilli et al. 2002, Grimmond and Oke 2002, Kusaka and Kimura 2004, Otte et al. 2004, Dandou et al. 2005), we describe a model that is simple enough to be compatible with structural, computational, and data constraints of a land surface model coupled to a global climate model, yet complex enough to enable exploration of physically-based processes known to be important in determining urban climatology. Several of the parameterizations are based on the Town Energy Balance (TEB) Model (Masson 2000, Masson et al. 2002, Lemonsu et al. 2004).

Urban Ecosystems and Climate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Characteristics of urban ecosystems and their effects on climate are summarized in Landsberg (1981), Oke (1987), Bonan (2002), and Arnfield (2003). Urban ecosystems can significantly alter the radiative, thermal, moisture, and aerodynamic characteristics of a region. The three-dimensional structure and geometrical arrangement of building walls and horizontal surfaces such as roads, sidewalks, parking lots, etc. combine to reduce the albedo of urban surfaces due to radiation trapping. Unlike solar radiation reflected from a horizontal surface, solar radiation impinging on urban surfaces such as walls and roads can experience multiple reflections and absorptions, resulting in increased absorption of radiation. Similarly, longwave radiation emitted by urban surfaces can be re-absorbed by these surfaces resulting in less longwave radiation loss to the atmosphere. The ratio of building height to canyon floor width is important in determining the degree of radiation trapping (Oke 1981, Oke et al. 1991).

The materials used for the construction of buildings and roads (e.g., dense concrete and asphalt) generally have higher heat capacity and thermal conductivity than some natural surfaces such as dry soils (Oke 1987). This results in higher thermal admittance and contributes to the ability of urban surfaces to store sensible heat during the day and release it at night. The importance of thermal properties in contributing to differences between urban and rural sites depends on the types of materials used in urban construction, the contrast in thermal admittance between the urban region and surrounding rural environs, and the building geometry which establishes the relative surface area and importance of roof, walls, and canyon floor (Oke et al. 1991).

Energy consumption due to building heating and cooling, manufacturing, transportation, and human metabolism releases waste heat to the urban environment. Such anthropogenic sources of heat can be substantial in some cases and should be accounted for in studies of the urban energy budget. As an extreme example, Ichinose et al. (1999) found that the total anthropogenic heat flux in central exceeded 400 W m-2 in daytime and a maximum value of 1590 W m-2 in winter. The contribution of waste heat sources from building heating and cooling may depend on population density, external climate, and socio-economic factors such as human adaptability and comfort levels, and economic status. The presence of insulation, characterized by low thermal admittance, may reduce the contribution of waste heat from heating and cooling. Waste heat fluxes from transportation have a distinct diurnal cycle due to morning/evening rush hours (Sailor and Lu, 2004). Generally, human metabolism contributes less than 5% of total anthropogenic flux in the (Sailor and Lu, 2004).

The urban surface is characterized by a preponderance of impervious surfaces, which reduce water storage capacity and surface moisture availability (Oke 1982). The evapotranspiration flux in urban regions is thus generally lower compared to vegetated surfaces, which may increase surface and air temperatures. On the other hand, vegetated surfaces within urban areas are frequently irrigated (e.g., lawns and parks) resulting in more water availability and higher latent heat fluxes than might be expected from natural vegetation. The presence and amount of vegetated or pervious surfaces can influence the magnitude of the heat island effect (Sailor 1995, Upmanis et al. 1998, Avissar 1996). Impervious surfaces also affect the hydrological cycle by reducing infiltration compared to rural areas, thereby converting more precipitation into surface runoff (Oke 1987, Bonan 2002).

The arrangement of large roughness elements (e.g., buildings, trees) in an urban region generally increases the frictional drag of the surface on the atmospheric winds and thus reduces the mean wind speed and turbulent mixing within the urban canopy compared to more open rural areas (Oke 1987). A notable exception to this may occur during periods of weak regional winds when warm urban air creates low-level rural-urban breezes. Lower within-canopy winds can reduce total turbulent heat transport from urban surfaces and increase their surface temperature. The synoptic wind speed is an important control on the urban heat island (Landsberg 1981). Higher winds may effectively remove heat faster than the urban fabric generates it.

The geographic location of urban areas and the characteristics of the surrounding rural area influence the urban climate. For instance, many tropical heat islands are smaller than expected based on population size. Where cities are surrounded by wet rural surfaces, slower cooling by these rural surfaces due to higher thermal admittance may reduce heat island magnitudes, especially in tropical climates (Oke et al. 1991). Local wind systems may impact urban climates as well. For example, coastal cities may experience cooling of urban temperatures when ocean surface temperatures are cooler than the land and winds blow onshore. Cold-air drainage from surrounding mountainous areas may reduce urban warming as well at certain times (Comrie 2000).

Urban regions have increased downward longwave radiation from the overlying atmosphere due to trapping and re-emission from polluted layers and/or from vertical advection of warm surface air above the city. Reduced incoming solar radiation due to reflection from atmospheric aerosols may compensate for this increase in longwave forcing. Note that in order to model these particular urban effects, the land model must also deliver biogeochemical fluxes (e.g., particulates, sulphur compounds, hydrocarbons, etc.) to the atmospheric model in addition to heat and moisture fluxes. The atmospheric model must then be able to diffuse or transport these trace species and determine their interaction with radiation and clouds. It has also been established that urban regions have effects on clouds and precipitation although the underlying mechanisms are still being debated. Climate modeling systems with detailed urban parameterizations may help to understand these mechanisms (Shepherd 2005).

As mentioned briefly in the previous section, many of the characteristics of the urban ecosystem discussed above contribute to one of the most striking effects of the urban environment on climate, the heat island effect. The present model is designed to represent the urban energy balance and provide insight into issues such as the urban heat island, its causes and potential mitigation strategies, as well as the effects of climate change on urban areas. When coupled to an atmospheric model, interactions between the urban surface and the atmosphere can be investigated.

Atmospheric Coupling and Model Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The atmospheric model within CCSM requires fluxes of sensible and latent heat and momentum between the surface and lowest atmospheric model level as well as emitted longwave and reflected shortwave radiation (Figure 1.1). These must be provided at a time step that resolves the diurnal cycle. Over other types of land surfaces, the fluxes are determined by current parameterizations in CLM. An objective of this technical note is to describe a set of parameterizations that determines the fluxes from an urban surface. The vertical spatial domain of the urban model extends from the top of the urban canopy layer (UCL) down to the depth of zero vertical heat flux in the ground (Oke 1987). The current state of the atmosphere and downwelling fluxes (Table 1.1) at a given time step is used to force the urban model. The urban model provides fluxes that are area-averaged with other land cover (e.g., forests, cropland) if present within the grid cell. The area-averaged fluxes (Table 1.2) are used as lower boundary conditions by the atmospheric model.

Land surface heterogeneity in the Community Land Model (CLM) is represented as a nested subgrid hierarchy (Figure 1.2) in which grid cells are composed of multiple landunits, snow/soil columns, and plant functional types (PFTs). Each grid cell can have a different number of landunits, each landunit can have a different number of columns, and each column can have multiple PFTs. The first subgrid level, the landunit, is intended to capture the broadest spatial patterns of subgrid heterogeneity. The model described here is designed to represent urban landunits. Further division of the urban surface into urban landuse classes such as, for example, city core, industrial/commercial, and suburban is possible by specifying these classes as individual landunits.

The representation of the urban landunit is based on the canyon concept of Oke (1987). In this approach, the considerable complexity of the urban surface is reduced to a single urban canyon consisting of a canyon floor of width :math:`W` bordered by two facing buildings of height :math:`H` (Figure 1.3). Although the canyon floor is intended to represent various surfaces such as roads, parking lots, sidewalks, and residential lawns, etc., for convenience we henceforth refer to the canyon floor as a road. The urban canyon consists of roof, sunlit and shaded wall, and pervious and impervious road, each of which are treated as columns within the landunit (Figure 1.2). The impervious road is intended to represent surfaces that are impervious to water infiltration (e.g., roads, parking lots, sidewalks) while the pervious road is intended to represent surfaces such as residential lawns and parks which may have active hydrology.

The approach used here to represent pervious surfaces is different than many urban schemes designed for use within mesoscale and global models. Most urban schemes use a separate land surface model scheme to represent the effects of pervious surfaces on urban climate. For example, the urban surface in the mesoscale model Meso-NH is modeled using the TEB and ISBA (Interactions between Soil, Biosphere, and Atmosphere) schemes for urban and pervious (e.g., vegetated) surfaces, respectively (Lemonsu and Masson 2002). Fluxes from each scheme are combined according to their relative areas. A comparable approach could be implemented using the CLM scheme for vegetated surfaces; however, this presents several disadvantages for our application. First, the pervious surface would need to be assigned to an additional landunit and specially identified to distinguish it from the other vegetated landunit within the gridcell. Second, the pervious and urban landunits would then need to be aggregated according to their relative areas in a post-processing sense to estimate the composite urban effects. Third, in the Meso-NH approach, the pervious surface only interacts indirectly with the canyon air through its influence on the atmospheric model. Here, including the pervious surface within the urban canyon solves these difficulties. Thus, the pervious surface is an integral part of the urban system and interacts directly with UCL air properties such as temperature and specific humidity. Yet, implementation of a sophisticated scheme for the pervious surface, such as the vegetation scheme in CLM, within the urban canyon is problematic because of computational and data requirements. Here, we choose a simplified bulk parameterization scheme to represent latent heat flux from pervious urban surfaces (Chapter 3).

Note that the urban columns interact radiatively with one another through multiple exchanges of longwave and shortwave radiation (chapter 2). The heat and moisture fluxes from each surface interact with each other through a bulk air mass that represents air in the UCL for which specific humidity and temperature are predicted (chapter 3). We model the UCL plus the air above the roof (Figure 1.1). This allows for mixing of above-roof air with canyon air.

Figure 1.1. Schematic of urban and atmospheric coupling. The urban model is forced by the atmospheric wind ( :math:`{u_{atm}}` ), temperature ( :math:`{T_{atm}}` ), specific humidity ( :math:`{q_{atm}}` ), precipitation ( :math:`{P_{atm}}` ), solar ( :math:`{S_{atm}}\, \downarrow` ) and longwave ( :math:`{L_{atm}}\, \downarrow` ) radiation at reference height :math:`{z'_{atm}}` . Fluxes from the urban landunit to the atmosphere are turbulent sensible ( :math:`H` ) and latent heat ( :math:`\lambda E` ), momentum ( :math:`\tau` ), albedo ( :math:`I \uparrow` ), emitted longwave ( :math:`L \uparrow` ), and absorbed shortwave ( :math:`\overrightarrow S` ) radiation. Air temperature ( :math:`{T_{ac}}` ), specific humidity ( :math:`{q_{ac}}` ), and wind speed ( :math:`{u_c}` ) within the urban canopy layer are diagnosed by the urban model. :math:`H` is the average building height.

Figure 1.2. CLM subgrid hierarchy emphasizing the structure of urban landunits.

Figure 1.3. The urban canyon.

Table 1.1. Atmospheric input to urban model

1The atmospheric reference height received from the atmospheric model :math:`{z'_{atm}}` is assumed to be the height above the surface defined as the roughness length :math:`{z_0}` plus displacement height :math:`{z_d}` . Thus, the reference height used for flux computations (chapter 3) is :math:`{z_{atm}} = {z'_{atm}} + {z_0} + {z_d}` . The reference heights for temperature, wind, and specific humidity ( :math:`{z_{atm,\,h}}` , :math:`{z_{atm,\,m}}` , :math:`{z_{atm,\,w}}` ) are required. These are set equal to :math:`{z_{atm}}` .

2The provides convective and large-scale liquid and solid precipitation, which are added to yield total liquid precipitation :math:`{q_{rain}}` and solid precipitation :math:`{q_{sno}}` .

3These are provided by the atmospheric model but not used by the urban model.

Density of air ( :math:`\rho _{atm}` ) (kg m-3) is also required but is calculated directly from

.. math::

 \rho _{atm} = \frac{P_{atm} - 0.378e_{atm}}{R_{da}T_{atm}}

where :math:`{P_{atm}}` is atmospheric pressure (Pa), :math:`{e_{atm}}` is atmospheric vapor pressure (Pa), :math:`{R_{da}}` is the gas constant for dry air (J kg-1 K-1) (Table 1.4), and :math:`{T_{atm}}` is the atmospheric temperature (K). The atmospheric vapor pressure :math:`{e_{atm}}` is derived from atmospheric specific humidity :math:`{q_{atm}}` (kg kg-1) as :math:`{e_{atm}} = \frac{{{q_{atm}}{P_{atm}}}}{{0.622 + 0.378{q_{atm}}}}` .

Table 1.2. Urban model output to atmospheric model

1 :math:`\lambda` is either the latent heat of vaporization :math:`{\lambda _{vap}}` or latent heat of sublimation :math:`{\lambda _{sub}}` (J kg-1) (Table 1.4) depending on the thermal state of surface water on the roof, pervious and impervious road.

2These are set to zero for urban areas.

Biogeophysical Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Biogeophysical processes are simulated for each of the five urban columns and each column maintains its own prognostic variables (e.g., surface temperature). The processes simulated include:

Absorption and reflection of solar radiation (chapter 2)

Absorption, reflection, and emission of longwave radiation (chapter 2)

Momentum, sensible heat, and latent heat fluxes (chapter 3)

Anthropogenic heat fluxes to the canyon air due to waste heat from building heating/air conditioning (chapter 3). An example of parameterizing traffic heat fluxes is given in Oleson et al. (2008b), however, traffic heat fluxes are not currently included in the global implementation of the model.

Heat transfer in roofs, building walls, and the road including phase change (chapter 4)

Hydrology [roofs - storage of liquid and solid precipitation (ponding and dew), surface runoff; walls – hydrologically inactive; impervious road – storage of liquid and solid precipitation (ponding and dew), surface runoff; pervious road - infiltration, surface runoff, sub-surface drainage, redistribution of water within the column] (chapter 5).

Model Requirements
----------------------

Initialization
~~~~~~~~~~~~~~~~~~~~

Initialization of the urban model (i.e., providing the model with initial temperature and moisture states) depends on the type of run (startup or restart) (see the CLM4 User’s Guide). An initial run starts the model from either initial conditions that are set internally in the Fortran code (referred to as arbitrary initial conditions) or from an initial conditions dataset that enables the model to start from a spun up state (i.e., where the urban landunit is in equilibrium with the simulated climate). In restart runs, the model is continued from a previous simulation and initialized from a restart file that ensures that the output is bit-for-bit the same as if the previous simulation had not stopped. The fields that are required from the restart or initial conditions files can be obtained by examining the code. Arbitrary initial conditions are specified as follows.

All urban columns consist of fifteen layers to be consistent with CLM4. Generally, temperature calculations are done over all layers, :math:`{N_{levgrnd}} = 15` , while hydrology calculations for the pervious road are done over the top ten layers, :math:`{N_{levsoi}} = 10` , the bottom five layers being specified as bedrock. Pervious and impervious road are initialized with temperatures (surface :math:`{T_g}` , and layers :math:`{T_i}` , for layers :math:`i = 1, \ldots ,{N_{levgrnd}}` ) of 274 K. Roof, sunwall, and shadewall are initialized to 292K. This relatively high temperature is to avoid initialization shock from large space heating/air conditioning and waste heat fluxes. All surfaces are initialized with no snow ( :math:`{W_{sno}} = 0` ). Roof and impervious road are initialized with no ponded water, while the pervious road soil layers :math:`i = 1, \ldots ,{N_{levsoi}}` are initialized with volumetric soil water content :math:`{\theta _i} = 0.3` mm3 mm-3 and layers :math:`i = {N_{levsoi}} + 1, \ldots ,{N_{levgrnd}}` are initialized :math:`{\theta _i} = 0.0` mm3 mm-3. The soil liquid water and ice contents are initialized as :math:`{w_{liq,\,i}} = \Delta {z_i}{\rho _{liq}}{\theta _i}` and :math:`{w_{ice,\,i}} = 0` , where :math:`{\rho _{liq}}` is the density liquid water (kg m-3) (Table 1.4). The pervious road is initialized with water stored in the unconfined aquifer and unsaturated soil 

.. math::

 {W_a} = {W_t} = 4800

 mm and water table depth 

.. math::

 {z_\nabla } = 4.8

 m.

Surface Data
~~~~~~~~~~~~~~~~~~

Required input data for urban landunits are listed in Table 1.3. This data is provided by the surface dataset at the required spatial resolution (see the CLM4 User’s Guide). Present day global urban extent and urban properties were developed by Jackson et al. (2010). Urban extent, defined for four classes [tall building district (TBD), and high, medium, and low density (HD, MD, LD)] was derived from LandScan 2004, a population density dataset derived from census data, nighttime lights satellite observations, road proximity, and slope (Dobson et al., 2000). The urban extent data is aggregated from the original 1 km resolution to a 0.5° by 0.5° global grid. For this particular implementation, only the sum of the TBD, HD, and MD classes are used to define urban extent as the LD class is highly rural and likely better modeled as a vegetated surface.

For each of 33 distinct regions across the globe, thermal (e.g., heat capacity and thermal conductivity), radiative (e.g., albedo and emissivity) and morphological (e.g., height to width ratio, roof fraction, average building height, and pervious fraction of the canyon floor) properties of roof/wall/road are provided by Jackson et al. (2010) for each of the four density classes. Building interior minimum and maximum temperatures are prescribed based on climate and socioeconomic considerations. Urban parameters are determined for the 0.5° by 0.5° global grid based on the dominant density class by area. This prevents potentially unrealistic parameter values that may result if the density classes are averaged. As a result, the current global representation of urban is almost exclusively medium density. Future implementations of the model could represent each of the density classes as a separate landunit. The surface dataset creation routines (see CLM4 User’s Guide) aggregate the data to the desired resolution. It is surmised that the MODIS-based vegetation dataset used in CLM4 classifies built areas as bare soil, thus the urban extent preferentially replaces bare soil when it exists within the grid cell. A very small minimum threshold of 0.1% of the grid cell by area is used to resolve urban areas. An elevation threshold of 2200 m is used to eliminate urban areas where the grid cell surface elevation is significantly higher than the elevation the cities are actually at because of the coarse spatial resolution of the model. This prevents overestimates of anthropogenic heating in winter due to unrealistically cold temperatures.

Table 1.3. Input data required for the urban model

1This fraction is relative to the canyon floor.

2Required for layers 

.. math::

 i = 1,{N_{imprvrd}}

, derived from grid cell soil texture for other layers (section 4.3).

3Derived from grid cell soil texture ( :math:`\% sand,\,\% clay` ) (section 4.3).

4Obtained from grid cell soil texture ( :math:`\% sand,\,\% clay` ).

Physical Constants
~~~~~~~~~~~~~~~~~~~~~~~~

Physical constants, shared by all of the components in the CCSM, are presented in Table 1.4. Not all constants are necessarily used by the urban model.

Table 1.4. Physical constants

1Not shared by other components of the coupled modeling system.

Albedos and Radiative Fluxes
============================

The effects of geometry on the radiation balance of urban surfaces are a key driver of urban-rural energy balance differences (Oke et al. 1991). Shadowing of urban surfaces affects the incident radiation and thus temperature. Similar to vegetated surfaces, multiple reflections of radiation between urban surfaces must be accounted for (Harman et al. 2004). The net solar radiation and net longwave radiation, the net of which is the net radiation, are needed for each urban surface to drive turbulent and ground heat fluxes. The atmospheric model also requires radiative fluxes and albedo from the urban landunit, which are appropriately averaged with other landunits within the gridcell. The urban canyon unit is used to represent these radiative processes. Several simplifying assumptions are made. The effects of absorption, emission, and scattering of radiation by the canyon air are neglected and surfaces are assumed to be isotropic.

Albedo
----------

The albedo of each urban surface is a weighted combination of snow-free “ground” albedo and snow albedo. Only roof and road surfaces are affected by snow. The direct beam :math:`\alpha _{u,\,\Lambda }^\mu` and diffuse :math:`{\alpha _{u,\,\Lambda }}` albedos (where :math:`u` denotes roof, impervious or pervious road) are



.. math::

 \alpha _{u,\,\Lambda }^\mu = \alpha _{g,\,\Lambda }^\mu \left( {1 - {f_{u,\,sno}}} \right) + \alpha _{sno,\,\Lambda }^\mu {f_{u,\,sno}}





.. math::

 {\alpha _{u,\,\Lambda }} = {\alpha _{g,\,\Lambda }}\left( {1 - {f_{u,\,sno}}} \right) + {\alpha _{sno,\,\Lambda }}{f_{u,\,sno}}



where :math:`{f_{u,\,sno}}` is the fraction of the urban surface covered with snow which is calculated from (Bonan 1996)



.. math::

 {f_{u,\,sno}} = \frac{{{z_{u,\,sno}}}}{{0.05}} \leqslant 1

	.

The direct and diffuse “ground” albedos, :math:`\alpha _{g,\,\Lambda }^\mu` and :math:`{\alpha _{g,\,\Lambda }}` , where :math:`\Lambda` denotes either the visible (VIS) or near-infrared (NIR) waveband, are provided by the surface dataset (Table 1.3), and :math:`{z_{u,\,sno}}` is the depth of snow (m) (section 5.1). An estimate of snow albedo is made based on the parameterization of (1989) in which albedo depends on solar zenith angle, grain size, and soot content (e.g., as adopted by the Land Surface Model (LSM) (Bonan 1996)). Here, however, several simplifying assumptions are made due to uncertainties in how to apply such a parameterization to urban surfaces. A snow grain radius of 100 :math:`\mu m` (new powder snow, aged a few days) and a soot mass fraction of 1.5 :math:`\times {10^{ - 5}}` (arrived at by noting that the LSM global soot mass fraction is 5 :math:`\times {10^{ - 6}}` and Chylek et al. (1987) observed that soot concentrations in urban snowpacks averaged three times the concentration in rural snowpacks) are assumed. Direct and diffuse albedos are assumed to be equal. This yields :math:`\alpha _{sno,\,VIS}^\mu = {\alpha _{sno,\,VIS}} = 0.66` and :math:`\alpha _{sno,\,NIR}^\mu = {\alpha _{sno,\,NIR}} = 0.56` which fall about in the middle of the range given by Oke (1987).

Incident direct solar radiation
-----------------------------------

Unlike the horizontal roof surface, the direct beam solar radiation received by the walls and the road must be adjusted for orientation and shadowing. The analytical solution given below follows Masson (2000). First, let :math:`\theta` be the angle between the sun direction and the along-canyon axis and consider the case where the along-canyon axis is perpendicular to the sun direction ($\theta = {\pi \mathord{\left/

{\vphantom {\pi 2}} \right.

\kern-\nulldelimiterspace} 2} :math:`). In this case, as shown in Figure 2.1, if the solar zenith angle` \mu :math:`is greater than the critical solar zenith angle` {\mu _0} :math:`(` {\mu _0} = {\tan ^{ - 1}}\left( {{W \mathord{\left/

{\vphantom {W H}} \right.

\kern-\nulldelimiterspace} H}} \right) :math:`), the road is in full shade, and the sunlit wall is in partial sun. Conversely, if` \mu $ is less than 

.. math::

 {\vec L_{uc}} - \left( {{L_{uc}} \uparrow - {L_{atm}} \downarrow } \right) = 0

, the road is in partial sun and the sunlit wall is in full sun. Note that, radiatively, the pervious and impervious road are treated the same, although their albedos are specified separately and may differ (Table 1.3).

Figure 2.1. Elevation (side) view of direct beam solar radiation incident on urban canyon surfaces for solar zenith angle :math:`\mu > {\mu _0}` (top) and :math:`\mu \leqslant {\mu _0}` (bottom). :math:`{S_{atm}}\, \downarrow _\Lambda ^\mu` is the direct beam incident solar radiation incident on a horizontal surface from the atmosphere. The along-canyon axis is assumed to be perpendicular to the sun direction.

If the direct beam solar radiation received by a horizontal surface (i.e., as received by the roof) is :math:`{S_{atm}}\, \downarrow _\Lambda ^\mu` , then the solar radiation on the wall in full illumination ( :math:`\mu \leqslant {\mu _0}` ) is ${{\left( {{S_{atm}}\, \downarrow _\Lambda ^\mu \cos i} \right)} \mathord{\left/

{\vphantom {{\left( {{S_{atm}}\, \downarrow _\Lambda ^\mu \cos i} \right)} {\cos \mu }}} \right.

\kern-\nulldelimiterspace} {\cos \mu }} :math:`where` i :math:`is the incidence angle (Figure 2.1). Since` \cos i = \cos (90 - \mu ) = \sin \mu $, the solar radiation on the sunlit wall is

${S_{sunwall}}\, \downarrow _\Lambda ^\mu \left( {\theta = {\pi \mathord{\left/

{\vphantom {\pi 2}} \right.

\kern-\nulldelimiterspace} 2}} \right) = \tan \left( \mu \right){S_{atm}}\, \downarrow _\Lambda ^\mu & \mu \leqslant {\mu _0}$.

Note that this is twice the radiation received by the wall in Masson (2000) because here we force the other (shaded) wall to receive no solar radiation ( :math:`{S_{shdwall}}\, \downarrow _\Lambda ^\mu = 0` ). In the case of :math:`\mu > {\mu _0}` , the illuminated fraction is ${{\left( {H - y} \right)} \mathord{\left/

{\vphantom {{\left( {H - y} \right)} H}} \right.

\kern-\nulldelimiterspace} H} :math:`and` {S_{sunwall}}\, \downarrow _\Lambda ^\mu = \left[ {{{\left( {H - y} \right)} \mathord{\left/

{\vphantom {{\left( {H - y} \right)} H}} \right.

\kern-\nulldelimiterspace} H}} \right]\tan \mu \,{S_{atm}}\, \downarrow _\Lambda ^\mu :math:`. Since` \tan \mu = {W \mathord{\left/

{\vphantom {W {\left( {H - y} \right)}}} \right.

\kern-\nulldelimiterspace} {\left( {H - y} \right)}}$ this simplifies to

${S_{sunwall}}\, \downarrow _\Lambda ^\mu \left( {\theta = {\pi \mathord{\left/

{\vphantom {\pi 2}} \right.

\kern-\nulldelimiterspace} 2}} \right) = \frac{W}{H}{S_{atm}}\, \downarrow _\Lambda ^\mu & \mu > {\mu _0}$.

Since the road is a horizontal surface, ${S_{road}}\, \downarrow _\Lambda ^\mu = \left[ {{{\left( {W - x} \right)} \mathord{\left/

{\vphantom {{\left( {W - x} \right)} W}} \right.

\kern-\nulldelimiterspace} W}} \right]{S_{atm}}\, \downarrow _\Lambda ^\mu :math:`for` \mu \leqslant {\mu _0} :math:`. Since` x = H\tan \mu $, the direct solar radiation incident on the road (pervious and impervious) is

${S_{road}}\, \downarrow _\Lambda ^\mu \left( {\theta = {\pi \mathord{\left/

{\vphantom {\pi 2}} \right.

\kern-\nulldelimiterspace} 2}} \right) = \left\{ \begin{gathered}

0 & \mu > {\mu _0} \hfill \\

\left( {1 - \frac{H}{W}\tan \mu } \right){S_{atm}}\, \downarrow _\Lambda ^\mu & \mu \leqslant {\mu _0} \hfill \\

\end{gathered} \right\}$.

Equations and for the walls and equation for the road can now be expanded to account for any canyon orientation ($0 \leqslant \theta \leqslant {\pi \mathord{\left/

{\vphantom {\pi 2}} \right.

\kern-\nulldelimiterspace} 2} :math:`). If` \theta :math:`is the angle between the sun direction and the along-canyon axis (Figure 2.2), then the expression for the incidence angle is now` \cos i = \sin \mu \sin \theta $ and equation becomes

 :math:`{S_{sunwall}}\, \downarrow _\Lambda ^\mu \left( \theta \right) = \sin \theta \tan \mu \,{S_{atm}}\, \downarrow _\Lambda ^\mu & \mu \leqslant {\mu _0}` .

Figure 2.2. Plan view of direct beam solar radiation incident on urban canyon surfaces. :math:`{S_{atm}}\, \downarrow _\Lambda ^\mu` is the direct beam incident solar radiation incident on a horizontal surface from the atmosphere. :math:`\theta` is the angle between the along-canyon axis and the sun direction.

For the case of :math:`\mu > {\mu _0}` , ${S_{sunwall}}\, \downarrow _\Lambda ^\mu \left( \theta \right) = \left[ {{{\left( {H - y} \right)} \mathord{\left/

{\vphantom {{\left( {H - y} \right)} H}} \right.

\kern-\nulldelimiterspace} H}} \right]\sin \theta \tan \mu \,{S_{atm}}\, \downarrow _\Lambda ^\mu :math:`. However, now` \tan \mu = {{\left( {{W \mathord{\left/

{\vphantom {W {\sin \theta }}} \right.

\kern-\nulldelimiterspace} {\sin \theta }}} \right)} \mathord{\left/

{\vphantom {{\left( {{W \mathord{\left/

{\vphantom {W {\sin \theta }}} \right.

\kern-\nulldelimiterspace} {\sin \theta }}} \right)} {\left( {H - y} \right)}}} \right.

\kern-\nulldelimiterspace} {\left( {H - y} \right)}}$ and thus

 :math:`{S_{sunwall}}\, \downarrow _\Lambda ^\mu \left( \theta \right) = \frac{W}{H}{S_{atm}}\, \downarrow _\Lambda ^\mu & \mu > {\mu _0}` .

Similarly, for the road ( :math:`\mu \leqslant {\mu _0}` ), ${S_{road}}\, \downarrow _\Lambda ^\mu \left( \theta \right) = \left[ {{{\left( {{W \mathord{\left/

{\vphantom {W {\sin \theta }}} \right.

\kern-\nulldelimiterspace} {\sin \theta }} - x} \right)} \mathord{\left/

{\vphantom {{\left( {{W \mathord{\left/

{\vphantom {W {\sin \theta }}} \right.

\kern-\nulldelimiterspace} {\sin \theta }} - x} \right)} {\left( {{W \mathord{\left/

{\vphantom {W {\sin \theta }}} \right.

\kern-\nulldelimiterspace} {\sin \theta }}} \right)}}} \right.

\kern-\nulldelimiterspace} {\left( {{W \mathord{\left/

{\vphantom {W {\sin \theta }}} \right.

\kern-\nulldelimiterspace} {\sin \theta }}} \right)}}} \right]{S_{atm}}\, \downarrow _\Lambda ^\mu :math:`with` x = H\tan \mu $ simplifies to

${S_{road}}\, \downarrow _\Lambda ^\mu \left( \theta \right) = \left\{ \begin{gathered}

0 & \mu > {\mu _0} \hfill \\

\left( {1 - \frac{H}{W}\sin \theta \tan \mu } \right){S_{atm}}\, \downarrow _\Lambda ^\mu & \mu \leqslant {\mu _0} \hfill \\

\end{gathered} \right\}$.

Note that the critical solar zenith angle is now

${\mu _0} = {\tan ^{ - 1}}\left( {\frac{{{W \mathord{\left/

{\vphantom {W {\sin \theta }}} \right.

\kern-\nulldelimiterspace} {\sin \theta }}}}{H}} \right)$.

Equations , , and are integrated over all canyon orientations ($0 \leqslant \theta \leqslant {\pi \mathord{\left/

{\vphantom {\pi 2}} \right.

\kern-\nulldelimiterspace} 2} :math:`). The integration is done in two parts, first from` \theta = 0 :math:`to` \theta = {\theta _0} :math:`, and second from` \theta = {\theta _0} :math:`to` \theta = {\pi \mathord{\left/

{\vphantom {\pi 2}} \right.

\kern-\nulldelimiterspace} 2} :math:`, where` {\theta _0}$ is the critical canyon orientation for which the road is no longer illuminated. This can be derived from Equation and is



.. math::

 {\theta _0} = {\sin ^{ - 1}}\left[ {\min \left( {\frac{W}{{H\tan \mu }},\,1} \right)} \right]

.

The integrations thus are

 :math:`{S_{sunwall}}\, \downarrow _\Lambda ^\mu = \frac{4}{{2\pi }}\int\limits_0^{{\theta _0}} {\sin \theta \tan \mu \,{S_{atm}}\, \downarrow _\Lambda ^\mu d\theta } + \frac{4}{{2\pi }}\int\limits_{{\theta _0}}^{\frac{\pi }{2}} {\frac{W}{H}\,{S_{atm}}\, \downarrow _\Lambda ^\mu d\theta }` 

and

 :math:`{S_{road}}\, \downarrow _\Lambda ^\mu = \frac{4}{{2\pi }}\int\limits_0^{{\theta _0}} {\left( {1 - \frac{H}{W}\sin \theta \tan \mu } \right)\,{S_{atm}}\, \downarrow _\Lambda ^\mu d\theta }` .

The direct beam solar radiation incident on the roof, walls and road is therefore

 :math:`{S_{roof}}\, \downarrow _\Lambda ^\mu = {S_{atm}}\, \downarrow _\Lambda ^\mu` ,

 :math:`{S_{shdwall}}\, \downarrow _\Lambda ^\mu = 0` ,

 :math:`{S_{sunwall}}\, \downarrow _\Lambda ^\mu = 2{S_{atm}}\, \downarrow _\Lambda ^\mu \left[ {\frac{W}{H}\left( {\frac{1}{2} - \frac{{{\theta _0}}}{\pi }} \right) + \frac{1}{\pi }\tan \mu \left( {1 - \cos {\theta _0}} \right)} \right]` ,

 :math:`{S_{road}}\, \downarrow _\Lambda ^\mu = {S_{imprvrd}}\, \downarrow _\Lambda ^\mu = {S_{prvrd}}\, \downarrow _\Lambda ^\mu = {S_{atm}}\, \downarrow _\Lambda ^\mu \left[ {\frac{{2{\theta _0}}}{\pi } - \frac{2}{\pi }\frac{H}{W}\tan \mu \left( {1 - \cos {\theta _0}} \right)} \right]` .

The direct incident solar radiation conserves energy as

$\begin{gathered}

{S_{atm}}\, \downarrow _\Lambda ^\mu = {f_{roof}}{S_{roof}}\, \downarrow _\Lambda ^\mu + \left( {1 - {f_{roof}}} \right) \hfill \\

& \left[ {{S_{imprvrd}}\, \downarrow _\Lambda ^\mu \left( {1 - {f_{prvrd}}} \right) + {S_{prvrd}}\, \downarrow _\Lambda ^\mu {f_{prvrd}} + \frac{H}{W}\left( {{S_{sunwall}}\, \downarrow _\Lambda ^\mu + {S_{shdwall}}\, \downarrow _\Lambda ^\mu } \right)} \right] \hfill \\

\end{gathered} $.

Note that the factor ${H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}$ for the sunlit wall and shaded wall converts the flux from watts per meter squared of wall area to watts per meter squared of ground area.

View factors
----------------

The interaction of diffuse radiation (i.e., longwave and scattered solar radiation) between urban surfaces depends on angle (view) factors, i.e., the fraction of diffusely distributed energy leaving one “surface” (e.g., sky) that arrives at another surface (e.g., wall) (Sparrow and Cess 1978). If :math:`{E_{ij}}` is the diffuse radiative flux density on surface :math:`j` that originated from surface :math:`i` and :math:`{E_i}` is the radiative flux from surface :math:`i` , then



.. math::

 {E_{ij}} = {F_{ij}}{E_i}



where :math:`{F_{ij}}` is the view factor. The view factors depend only on the geometrical configurations of the involved surfaces. A table of view factors for various configurations is provided in Appendix A of Sparrow and Cess (1978). For instance, the view factor for the radiation from the wall to the sky can be derived from configuration nine of Appendix A. If :math:`d{A_1}` is an infinitesimal element on surface 1 (i.e., wall) and :math:`{A_2}` is a finite surface (i.e., sky) (Figure 2.3), then the angle factor :math:`{F_{d{A_1} - {A_2}}}` for diffuse radiation leaving element :math:`d{A_1}` and arriving at :math:`{A_2}` is



.. math::

 {F_{d{A_1} - {A_2}}} = \frac{1}{{2\pi }}\left( {{{\tan }^{ - 1}}\frac{1}{Y} - AY{{\tan }^{ - 1}}A} \right)



where $A = {1 \mathord{\left/

{\vphantom {1 {\sqrt {{X^2} + {Y^2}} }}} \right.

\kern-\nulldelimiterspace} {\sqrt {{X^2} + {Y^2}} }} :math:`,` X = {a \mathord{\left/

{\vphantom {a b}} \right.

\kern-\nulldelimiterspace} b} :math:`, and` Y = c/b :math:`. Following Sakakibara (1996) and Kusaka et al. (2001), for an infinitely long canyon,` b = \infty :math:`,` a = W :math:`, and so the wall-sky view factor at distance` c$ from a point on the wall to the canyon top is



.. math::

 {\Psi _{wall - sky\,|\,c}} = \frac{1}{2}\left( {1 - \frac{c}{{\sqrt {{c^2} + {W^2}} }}} \right)

.

The total wall-sky view factor can be found by integrating the above equation over the height of the wall as



.. math::

 {\Psi _{wall - sky}} = \frac{1}{H}\int\limits_{c = 0}^{c = H} {\frac{1}{2}\left( {1 - \frac{c}{{\sqrt {{c^2} + {W^2}} }}} \right)dc} = \frac{{\frac{1}{2}\left( {\frac{H}{W} + 1 - \sqrt {1 + {{\left( {\frac{H}{W}} \right)}^2}} } \right)}}{{\frac{H}{W}}}

.

By the reciprocity rule ( :math:`{A_1}{F_{{A_1} - {A_2}}} = {A_2}{F_{{A_2} - {A_1}}}` ) (Sparrow and Cess 1978), the sky-wall view factor is



.. math::

 {\Psi _{sky - wall}} = \frac{H}{W}{\Psi _{wall - sky}}

.

When applied to equation , :math:`{\Psi _{sky - wall}}` will yield a flux density to the wall in terms of per unit sky area. In the radiation computations detailed below, the diffuse fluxes for the walls are solved in terms of per unit wall area. Dividing equation by the height to width ratio converts the view factor to per unit wall area. Thus,



.. math::

 {\Psi _{sky - wall}} = \frac{{\frac{1}{2}\left( {\frac{H}{W} + 1 - \sqrt {1 + {{\left( {\frac{H}{W}} \right)}^2}} } \right)}}{{\frac{H}{W}}}



Similarly, the view factor for radiation from the sky to the road and from road to sky can be solved and is



.. math::

 {\Psi _{sky - road}} = \frac{W}{W}{\Psi _{road - sky}} = {\Psi _{road - sky}} = \sqrt {1 + {{\left( {\frac{H}{W}} \right)}^2}} - \frac{H}{W}

.

By symmetry,



.. math::

 {\Psi _{wall - road}} = {\Psi _{wall - sky}}

,

and the other view factors can be deduced from conservation of energy as



.. math::

 {\Psi _{road - wall}} = \frac{1}{2}\left( {1 - {\Psi _{road - sky}}} \right)

,



.. math::

 {\Psi _{wall - wall}} = 1 - {\Psi _{wall - sky}} - {\Psi _{wall - road}}

.

The view factors are presented graphically in Figure 2.4. Note that the view factors for radiation from the walls to the other surfaces sum to one ( :math:`{\Psi _{wall - wall}} + {\Psi _{wall - road}} + {\Psi _{wall - sky}} = 1` ). Similarly, the view factors for radiation from the road to the other surfaces also sum to one ( :math:`{\Psi _{road - wall}} + {\Psi _{road - wall}} + {\Psi _{road - sky}} = 1` ). As Harman et al. (2004) notes, at low height to width ratios, the road-sky view factor is close to one, the wall-wall view factor is close to zero, and the wall sky view factor is close to one half. However, at these low height to width ratios, the wall area is small compared to the road or sky area, indicating that most of the radiative exchange occurs between the road and sky, as it would for a flat surface. At height to width ratios greater than one, most of the radiative interactions take place between the two walls and the wall and the road. These view factors are consistent with those given by both Masson (2000) and Harman et al. (2004).

Figure 2.3. Schematic representation of angle (view) factor between infinitesimal element :math:`d{A_1}` (e.g., a point on the wall) and finite surface :math:`{A_2}` (e.g., the sky) (after Sparrow and Cess (1978)).

Figure 2.4. View factors as a function of canyon height to width ratio. :math:`{\Psi _{road - sky}}` is the fraction of radiation reaching the sky from the road, :math:`{\Psi _{road - wall}}` is the fraction of radiation reaching the wall from the road, :math:`{\Psi _{wall - sky}}` is the fraction of radiation reaching the sky from the wall, :math:`{\Psi _{wall - road}}` is the fraction of radiation reaching the road from the wall, and :math:`{\Psi _{wall - wall}}` is the fraction of radiation reaching the wall from the opposite wall.

Incident diffuse solar radiation
------------------------------------

The two view factors needed to compute the incident diffuse solar radiation are :math:`{\Psi _{sky - road}}` (equation ) and :math:`{\Psi _{sky - wall}}` (equation ). The diffuse solar radiation incident on roof, walls and road is then

 :math:`{S_{roof}}\, \downarrow _\Lambda ^{} = {S_{atm}}\, \downarrow _\Lambda ^{}` ,

 :math:`{S_{imprvrd}}\, \downarrow _\Lambda ^{} = {S_{prvrd}}\, \downarrow _\Lambda ^{} = {S_{atm}}\, \downarrow _\Lambda ^{}{\Psi _{sky - road}}` ,

 :math:`{S_{shdwall}}\, \downarrow _\Lambda ^{} = {S_{atm}}\, \downarrow _\Lambda ^{}{\Psi _{sky - wall}}` ,

 :math:`{S_{sunwall}}\, \downarrow _\Lambda ^{} = {S_{atm}}\, \downarrow _\Lambda ^{}{\Psi _{sky - wall}}` .

The diffuse incident solar radiation conserves energy as

$\begin{gathered}

{S_{atm}}\, \downarrow _\Lambda ^{} = {f_{roof}}{S_{roof}}\, \downarrow _\Lambda ^{} + \left( {1 - {f_{roof}}} \right) \hfill \\

& \left[ {{S_{imprvrd}}\, \downarrow _\Lambda ^{}\left( {1 - {f_{prvrd}}} \right) + {S_{prvrd}}\, \downarrow _\Lambda ^{}{f_{prvrd}} + \frac{H}{W}\left( {{S_{sunwall}}\, \downarrow _\Lambda ^{} + {S_{shdwall}}\, \downarrow _\Lambda ^{}} \right)} \right] \hfill \\

\end{gathered} $.

Absorbed and reflected solar radiation
------------------------------------------

The direct and diffuse net (absorbed) and reflected solar radiation for the roof is



.. math::

 \vec S_{roof,\,\Lambda }^\mu = {S_{roof}}\, \downarrow _\Lambda ^\mu \left( {1 - \alpha _{roof,\,\Lambda }^\mu } \right)





.. math::

 \vec S_{roof,\,\Lambda }^{} = {S_{roof}}\, \downarrow _\Lambda ^{}\left( {1 - \alpha _{roof,\,\Lambda }^{}} \right)





.. math::

 {S_{roof}} \uparrow _\Lambda ^\mu = {S_{roof}}\, \downarrow _\Lambda ^\mu \left( {\alpha _{roof,\,\Lambda }^\mu } \right)





.. math::

 {S_{roof}} \uparrow _\Lambda ^{} = {S_{roof}}\, \downarrow _\Lambda ^{}\left( {\alpha _{roof,\,\Lambda }^{}} \right)

.

The net (absorbed) and reflected solar radiation for walls and road and the reflected solar radiation to the sky are determined numerically by allowing for multiple reflections until a convergence criteria is met to ensure radiation is conserved. The reflected radiation from each urban surface is absorbed and re-reflected by the other urban surfaces. For example, the radiation scattered from the sunlit wall to the road, the shaded wall, and the sky depends on the view factors :math:`{\Psi _{wall - road}}` , :math:`{\Psi _{wall - wall}}` , and :math:`{\Psi _{wall - sky}}` , respectively (Figure 2.4). The multiple reflections are accounted for in five steps:

Determine the initial absorption and reflection by each urban surface and distribute this radiation to the sky, road, and walls according to view factors.

Determine the amount of radiation absorbed and reflected by each urban surface after the initial reflection. The solar radiation reflected from the walls to the road is projected to road area by multiplying by the height to width ratio and the solar radiation reflected from the road to the walls is projected to wall area by dividing by the height to width ratio.

The absorbed radiation for the :math:`{i^{th}}` reflection is added to the total absorbed by each urban surface.

The reflected solar radiation for the :math:`{i^{th}}` reflection is distributed to the sky, road, and walls according to view factors.

The reflected solar radiation to the sky for the :math:`{i^{th}}` reflection is added to the total reflected solar radiation.

Steps 2-5 are repeated until a convergence criterion (absorbed radiation per unit incoming solar radiation for a given reflection is less than :math:`1 \times {10^{ - 5}}` ) is met to ensure radiation is conserved. Direct beam and diffuse radiation are solved independently but follow the same solution steps. The solution below is for the direct beam component.

The initial direct beam absorption ( :math:`i = 0` ) (step 1) by each urban surface is



.. math::

 \vec S_{imprvrd,\,\Lambda ,\,i = 0}^\mu = {S_{imprvrd}}\, \downarrow _\Lambda ^\mu \left( {1 - \alpha _{imprvrd,\,\Lambda }^\mu } \right)

,



.. math::

 \vec S_{prvrd,\,\Lambda ,i = 0}^\mu = {S_{prvrd}}\, \downarrow _\Lambda ^\mu \left( {1 - \alpha _{prvrd,\,\Lambda }^\mu } \right)

,



.. math::

 \vec S_{sunwall,\,\Lambda ,\,i = 0}^\mu = {S_{sunwall}}\, \downarrow _\Lambda ^\mu \left( {1 - \alpha _{sunwall,\,\Lambda }^\mu } \right)

,



.. math::

 \vec S_{shdwall,\,\Lambda ,\,i = 0}^\mu = {S_{shdwall}}\, \downarrow _\Lambda ^\mu \left( {1 - \alpha _{shdwall,\,\Lambda }^\mu } \right)

,



.. math::

 \vec S_{road,\,\Lambda ,\,i = 0}^\mu = \vec S_{imprvrd,\,\Lambda ,\,i = 0}^\mu \left( {1 - {f_{prvrd}}} \right) + \vec S_{prvrd,\,\Lambda ,\,i = 0}^\mu {f_{prvrd}}



where, for example, 

.. math::

 {S_{imprvrd}}\, \downarrow _\Lambda ^\mu

 is the incident direct solar radiation for the impervious road (equation ) and 

.. math::

 \alpha _{imprvrd,\,\Lambda }^\mu

 is the direct albedo for the impervious road after adjustment for snow (section 2.1). Similarly, the initial reflections from each urban surface are



.. math::

 S_{imprvrd}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = {S_{imprvrd}}\, \downarrow _\Lambda ^\mu \left( {\alpha _{imprvrd,\,\Lambda }^\mu } \right)

,



.. math::

 S_{prvrd}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = {S_{prvrd}}\, \downarrow _\Lambda ^\mu \left( {\alpha _{prvrd,\,\Lambda }^\mu } \right)

,



.. math::

 S_{road}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = {S_{imprvrd}}\, \downarrow _\Lambda ^\mu \left( {1 - {f_{prvrd}}} \right) + {S_{prvrd}}\, \downarrow _\Lambda ^\mu {f_{prvrd}}





.. math::

 S_{sunwall}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = {S_{sunwall}}\, \downarrow _\Lambda ^\mu \left( {\alpha _{sunwall,\,\Lambda }^\mu } \right)

,



.. math::

 S_{shdwall}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = {S_{shdwall}}\, \downarrow _\Lambda ^\mu \left( {\alpha _{shdwall,\,\Lambda }^\mu } \right)

,

The initial reflected solar radiation is distributed to sky, walls, and road according to view factors as



.. math::

 S_{imprvrd - sky}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = S_{imprvrd}^{} \uparrow _{\Lambda ,\,i = 0}^\mu {\Psi _{road - sky}}





.. math::

 S_{imprvrd - sunwall}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = S_{imprvrd}^{} \uparrow _{\Lambda ,\,i = 0}^\mu {\Psi _{road - wall}}





.. math::

 S_{imprvrd - shdwall}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = S_{imprvrd}^{} \uparrow _{\Lambda ,\,i = 0}^\mu {\Psi _{road - wall}}





.. math::

 S_{prvrd - sky}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = S_{prvrd}^{} \uparrow _{\Lambda ,\,i = 0}^\mu {\Psi _{road - sky}}





.. math::

 S_{prvrd - sunwall}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = S_{prvrd}^{} \uparrow _{\Lambda ,\,i = 0}^\mu {\Psi _{road - wall}}





.. math::

 S_{prvrd - shdwall}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = S_{prvrd}^{} \uparrow _{\Lambda ,\,i = 0}^\mu {\Psi _{road - wall}}





.. math::

 S_{road - sky}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = S_{road}^{} \uparrow _{\Lambda ,\,i = 0}^\mu {\Psi _{road - sky}}





.. math::

 S_{road - sunwall}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = S_{road}^{} \uparrow _{\Lambda ,\,i = 0}^\mu {\Psi _{road - wall}}





.. math::

 S_{road - shdwall}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = S_{road}^{} \uparrow _{\Lambda ,\,i = 0}^\mu {\Psi _{road - wall}}





.. math::

 S_{sunwall - sky}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = S_{sunwall}^{} \uparrow _{\Lambda ,\,i = 0}^\mu {\Psi _{wall - sky}}





.. math::

 S_{sunwall - road}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = S_{sunwall}^{} \uparrow _{\Lambda ,\,i = 0}^\mu {\Psi _{wall - road}}





.. math::

 S_{sunwall - shdwall}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = S_{sunwall}^{} \uparrow _{\Lambda ,\,i = 0}^\mu {\Psi _{wall - wall}}





.. math::

 S_{shdwall - sky}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = S_{shdwall}^{} \uparrow _{\Lambda ,\,i = 0}^\mu {\Psi _{wall - sky}}





.. math::

 S_{shdwall - road}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = S_{shdwall}^{} \uparrow _{\Lambda ,\,i = 0}^\mu {\Psi _{wall - road}}





.. math::

 S_{shdwall - sunwall}^{} \uparrow _{\Lambda ,\,i = 0}^\mu = S_{shdwall}^{} \uparrow _{\Lambda ,\,i = 0}^\mu {\Psi _{wall - wall}}



The direct beam solar radiation absorbed by each urban surface after the :math:`{i^{th}}` reflection (steps 2 and 3) is

\[\begin{gathered}

\vec S_{imprvrd,\,\Lambda ,\,i}^\mu = \vec S_{imprvrd,\,\Lambda ,\,i - 1}^\mu + \hfill \\

& \left( {1 - \alpha _{imprvrd,\,\Lambda }^\mu } \right)\left( {{S_{sunwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu + {S_{shdwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right)\frac{H}{W} \hfill \\

\end{gathered} \]

\[\begin{gathered}

\vec S_{prvrd,\,\Lambda ,\,i}^\mu = \vec S_{prvrd,\,\Lambda ,\,i - 1}^\mu + \hfill \\

& \left( {1 - \alpha _{prvrd,\,\Lambda }^\mu } \right)\left( {{S_{sunwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu + {S_{shdwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right)\frac{H}{W} \hfill \\

\end{gathered} \]

\[\begin{gathered}

\vec S_{sunwall,\,\Lambda ,\,i}^\mu = \vec S_{sunwall,\,\Lambda ,\,i - 1}^\mu + \hfill \\

& \left( {1 - \alpha _{sunwall,\,\Lambda }^\mu } \right)\left( {\frac{{{S_{road - sunwall}}\, \uparrow _{\Lambda ,\,i - 1}^\mu }}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}} + {S_{shdwall - sunwall}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right) \hfill \\

\end{gathered} \]

\[\begin{gathered}

\vec S_{shdwall,\,\Lambda ,\,i}^\mu = \vec S_{shdwall,\,\Lambda ,\,i - 1}^\mu + \hfill \\

& \left( {1 - \alpha _{shdwall,\,\Lambda }^\mu } \right)\left( {\frac{{{S_{road - shdwall}}\, \uparrow _{\Lambda ,\,i - 1}^\mu }}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}} + {S_{sunwall - shdwall}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right) \hfill \\

\end{gathered} \]

The radiation from the walls to the road (

.. math::

 {S_{sunwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu

,

.. math::

 {S_{shdwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu

) is in W m-2 of wall area and must be converted to W m-2 of road area by multiplying by the height to width ratio. Similarly, the radiation from the road to the walls must be converted from W m-2 of road area to W m-2 of wall area by dividing by the height to width ratio. The direct beam solar radiation reflected by each urban surface after the :math:`{i^{th}}` reflection is distributed to sky, road, and walls (step 4) according to



.. math::

 S_{imprvrd - sky}^{} \uparrow _{\Lambda ,\,i}^\mu = \alpha _{imprvrd,\,\Lambda }^\mu \left( {{S_{sunwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu + {S_{shdwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right)\frac{H}{W}{\Psi _{road - sky}}





.. math::

 S_{imprvrd - sunwall}^{} \uparrow _{\Lambda ,\,i}^\mu = \alpha _{imprvrd,\,\Lambda }^\mu \left( {{S_{sunwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu + {S_{shdwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right)\frac{H}{W}{\Psi _{road - wall}}





.. math::

 S_{imprvrd - shdwall}^{} \uparrow _{\Lambda ,\,i}^\mu = \alpha _{imprvrd,\,\Lambda }^\mu \left( {{S_{sunwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu + {S_{shdwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right)\frac{H}{W}{\Psi _{road - wall}}





.. math::

 S_{prvrd - sky}^{} \uparrow _{\Lambda ,\,i}^\mu = \alpha _{prvrd,\,\Lambda }^\mu \left( {{S_{sunwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu + {S_{shdwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right)\frac{H}{W}{\Psi _{road - sky}}





.. math::

 S_{prvrd - sunwall}^{} \uparrow _{\Lambda ,\,i}^\mu = \alpha _{prvrd,\,\Lambda }^\mu \left( {{S_{sunwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu + {S_{shdwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right)\frac{H}{W}{\Psi _{road - wall}}





.. math::

 S_{prvrd - shdwall}^{} \uparrow _{\Lambda ,\,i}^\mu = \alpha _{prvrd,\,\Lambda }^\mu \left( {{S_{sunwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu + {S_{shdwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right)\frac{H}{W}{\Psi _{road - wall}}



\[S_{road - sky}^{} \uparrow _{\Lambda ,\,i}^\mu = \left[ \begin{gathered}

\alpha _{imprvrd,\,\Lambda }^\mu \left( {{S_{sunwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu + {S_{shdwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right) \hfill \\

\left( {1 - {f_{prvrd}}} \right) + \alpha _{prvrd,\,\Lambda }^\mu \hfill \\

\left( {{S_{sunwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu + {S_{shdwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right){f_{prvrd}} \hfill \\

\end{gathered} \right]\frac{H}{W}{\Psi _{road - sky}}\]

\[S_{road - sunwall}^{} \uparrow _{\Lambda ,\,i}^\mu = \left[ \begin{gathered}

\alpha _{imprvrd,\,\Lambda }^\mu \left( {{S_{sunwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu + {S_{shdwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right) \hfill \\

\left( {1 - {f_{prvrd}}} \right) + \alpha _{prvrd,\,\Lambda }^\mu \hfill \\

\left( {{S_{sunwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu + {S_{shdwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right){f_{prvrd}} \hfill \\

\end{gathered} \right]\frac{H}{W}{\Psi _{road - wall}}\]

\[S_{road - shdwall}^{} \uparrow _{\Lambda ,\,i}^\mu = \left[ \begin{gathered}

\alpha _{imprvrd,\,\Lambda }^\mu \left( {{S_{sunwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu + {S_{shdwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right) \hfill \\

\left( {1 - {f_{prvrd}}} \right) + \alpha _{prvrd,\,\Lambda }^\mu \hfill \\

\left( {{S_{sunwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu + {S_{shdwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right){f_{prvrd}} \hfill \\

\end{gathered} \right]\frac{H}{W}{\Psi _{road - wall}}\]

\[S_{sunwall - sky}^{} \uparrow _{\Lambda ,\,i}^\mu = \alpha _{sunwall,\,\Lambda }^\mu \left( {\frac{{{S_{road - sunwall}}\, \uparrow _{\Lambda ,\,i - 1}^\mu }}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}} + {S_{shdwall - sunwall}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right){\Psi _{wall - sky}}\]

\[S_{sunwall - road}^{} \uparrow _{\Lambda ,\,i}^\mu = \alpha _{sunwall,\,\Lambda }^\mu \left( {\frac{{{S_{road - sunwall}}\, \uparrow _{\Lambda ,\,i - 1}^\mu }}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}} + {S_{shdwall - sunwall}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right){\Psi _{wall - road}}\]

\[S_{sunwall - shdwall}^{} \uparrow _{\Lambda ,\,i}^\mu = \alpha _{sunwall,\,\Lambda }^\mu \left( {\frac{{{S_{road - sunwall}}\, \uparrow _{\Lambda ,\,i - 1}^\mu }}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}} + {S_{shdwall - sunwall}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right){\Psi _{wall - wall}}\]

\[S_{shdwall - sky}^{} \uparrow _{\Lambda ,\,i}^\mu = \alpha _{shdwall,\,\Lambda }^\mu \left( {\frac{{{S_{road - shdwall}}\, \uparrow _{\Lambda ,\,i - 1}^\mu }}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}} + {S_{sunwall - shdwall}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right){\Psi _{wall - sky}}\]

\[S_{shdwall - road}^{} \uparrow _{\Lambda ,\,i}^\mu = \alpha _{shdwall,\,\Lambda }^\mu \left( {\frac{{{S_{road - shdwall}}\, \uparrow _{\Lambda ,\,i - 1}^\mu }}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}} + {S_{sunwall - shdwall}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right){\Psi _{wall - road}}\]

\[S_{shdwall - sunwall}^{} \uparrow _{\Lambda ,\,i}^\mu = \alpha _{shdwall,\,\Lambda }^\mu \left( {\frac{{{S_{road - shdwall}}\, \uparrow _{\Lambda ,\,i - 1}^\mu }}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}} + {S_{sunwall - shdwall}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right){\Psi _{wall - wall}}\].

The reflected solar radiation to the sky is added to the total reflected solar radiation (step 5) for each urban surface as



.. math::

 S_{imprvrd}^{} \uparrow _{\Lambda ,\,i + 1}^\mu = S_{imprvrd}^{} \uparrow _{\Lambda ,\,i - 1}^\mu + S_{imprvrd - sky}^{} \uparrow _{\Lambda ,\,i}^\mu





.. math::

 S_{prvrd}^{} \uparrow _{\Lambda ,\,i + 1}^\mu = S_{prvrd}^{} \uparrow _{\Lambda ,\,i - 1}^\mu + S_{prvrd - sky}^{} \uparrow _{\Lambda ,\,i}^\mu





.. math::

 S_{sunwall}^{} \uparrow _{\Lambda ,\,i + 1}^\mu = S_{sunwall}^{} \uparrow _{\Lambda ,\,i - 1}^\mu + S_{sunwall - sky}^{} \uparrow _{\Lambda ,\,i}^\mu





.. math::

 S_{shdwall}^{} \uparrow _{\Lambda ,\,i + 1}^\mu = S_{shdwall}^{} \uparrow _{\Lambda ,\,i - 1}^\mu + S_{shdwall - sky}^{} \uparrow _{\Lambda ,\,i}^\mu

.

The system of equations (Equations -) is iterated for :math:`i = 50` reflections or until the absorption for the :math:`{i^{th}}` reflection is less than a nominal amount



.. math::

 \max \left( {\frac{{\vec S_{road,\,\Lambda ,\,i}^\mu }}{{{S_{atm}}\, \downarrow _\Lambda ^\mu }},\frac{{\vec S_{sunwall,\,\Lambda ,\,i}^\mu }}{{{S_{atm}}\, \downarrow _\Lambda ^\mu }},\frac{{\vec S_{shdwall,\,\Lambda ,\,i}^\mu }}{{{S_{atm}}\, \downarrow _\Lambda ^\mu }}} \right) < 1 \times {10^{ - 5}}



where 

.. math::

 \vec S_{sunwall,\,\Lambda ,\,i}^\mu

 (equation ) and 

.. math::

 \vec S_{shdwall,\,\Lambda ,\,i}^\mu

 (equation ) are the direct beam solar radiation absorbed by the sunlit wall and shaded wall on the :math:`{i^{th}}` reflection, and

\[\begin{gathered}

\vec S_{road,\,\Lambda ,\,i}^\mu = \left( {1 - \alpha _{imprvrd,\,\Lambda }^\mu } \right)\left( {{S_{sunwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu + {S_{shdwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right)\frac{H}{W}\left( {1 - {f_{prvrd}}} \right) \hfill \\

& + \left( {1 - \alpha _{prvrd,\,\Lambda }^\mu } \right)\left( {{S_{sunwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu + {S_{shdwall - road}}\, \uparrow _{\Lambda ,\,i - 1}^\mu } \right)\frac{H}{W}{f_{prvrd}} \hfill \\

\end{gathered} \]

is the direct beam solar radiation absorbed by the road on the :math:`{i^{th}}` reflection.

The total direct beam and diffuse solar radiation reflected by the urban canyon (walls and road) is

\[\begin{gathered}

S_{uc}^{} \uparrow _\Lambda ^\mu = S_{imprvrd}^{} \uparrow _{\Lambda ,\,i = n + 1}^\mu \left( {1 - {f_{prvrd}}} \right) + S_{prvrd}^{} \uparrow _{\Lambda ,\,i = n + 1}^\mu {f_{prvrd}} \hfill \\

& + \left( {S_{sunwall}^{} \uparrow _{\Lambda ,\,i = n + 1}^\mu + S_{shdwall}^{} \uparrow _{\Lambda ,\,i = n + 1}^\mu } \right)\frac{H}{W} \hfill \\

\end{gathered} \]

\[\begin{gathered}

S_{uc}^{} \uparrow _\Lambda ^{} = S_{imprvrd}^{} \uparrow _{\Lambda ,\,i = n + 1}^{}\left( {1 - {f_{prvrd}}} \right) + S_{prvrd}^{} \uparrow _{\Lambda ,\,i = n + 1}^{}{f_{prvrd}} \hfill \\

& + \left( {S_{sunwall}^{} \uparrow _{\Lambda ,\,i = n + 1}^{} + S_{shdwall}^{} \uparrow _{\Lambda ,\,i = n + 1}^{}} \right)\frac{H}{W} \hfill \\

\end{gathered} \]

while the total absorbed is

\[\begin{gathered}

\vec S_{uc,\,\Lambda }^\mu = \vec S_{imprvrd,\,\Lambda ,\,i = n}^\mu \left( {1 - {f_{prvrd}}} \right) + \vec S_{prvrd,\,\Lambda ,\,i = n}^\mu \,{f_{prvrd}} \hfill \\

& + \left( {\vec S_{sunwall,\,\Lambda ,\,i = n}^\mu + \vec S_{shdwall,\,\Lambda ,\,i = n}^\mu } \right)\frac{H}{W} \hfill \\

\end{gathered} \]

\[\begin{gathered}

\vec S_{uc,\,\Lambda }^{} = \vec S_{imprvrd,\,\Lambda ,\,i = n}^{}\left( {1 - {f_{prvrd}}} \right) + \vec S_{prvrd,\,\Lambda ,\,i = n}^{}\,{f_{prvrd}} \hfill \\

& + \left( {\vec S_{sunwall,\,\Lambda ,\,i = n}^{} + \vec S_{shdwall,\,\Lambda ,\,i = n}^{}} \right)\frac{H}{W} \hfill \\

\end{gathered} \].

Solar radiation in the urban canyon is conserved as

$\begin{gathered}

{S_{road}}\, \downarrow _\Lambda ^\mu + \left( {{S_{sunwall}}\, \downarrow _\Lambda ^\mu + {S_{shdwall}}\, \downarrow _\Lambda ^\mu } \right)\frac{H}{W} + {S_{road}}\, \downarrow _\Lambda ^{} + \left( {{S_{sunwall}}\, \downarrow _\Lambda ^{} + {S_{shdwall}}\, \downarrow _\Lambda ^{}} \right)\frac{H}{W} \hfill \\
* \left( {\vec S_{uc,\,\Lambda }^\mu + \vec S_{uc,\,\Lambda }^{} + S_{uc}^{} \uparrow _\Lambda ^\mu + S_{uc}^{} \uparrow _\Lambda ^{}} \right) = 0 \hfill \\

\end{gathered} $.

The direct beam and diffuse urban canyon albedos are



.. math::

 \alpha _{uc,\,\Lambda }^\mu = \frac{{S_{uc}^{} \uparrow _\Lambda ^\mu }}{{{S_{road}}\, \downarrow _\Lambda ^\mu + \left( {{S_{sunwall}}\, \downarrow _\Lambda ^\mu + {S_{shdwall}}\, \downarrow _\Lambda ^\mu } \right)\frac{H}{W}}}

,



.. math::

 \alpha _{uc,\,\Lambda }^{} = \frac{{S_{uc}^{} \uparrow _\Lambda ^{}}}{{{S_{road}}\, \downarrow _\Lambda ^{} + \left( {{S_{sunwall}}\, \downarrow _\Lambda ^{} + {S_{shdwall}}\, \downarrow _\Lambda ^{}} \right)\frac{H}{W}}}

.

The total absorbed solar radiation for the urban canopy (road, walls, and roof) is



.. math::

 \vec S = \sum\limits_\Lambda {\left[ {{W_{roof}}\left( {\vec S_{roof,\,\Lambda }^\mu + \vec S_{roof,\,\Lambda }^{}} \right) + \left( {1 - {W_{roof}}} \right)\left( {\vec S_{uc,\,\Lambda }^\mu + \vec S_{uc,\,\Lambda }^{}} \right)} \right]}



Figure 2.5 shows the solar radiation absorbed by urban surfaces for a range of height to width ratios and two solar zenith angles. The absorbed solar radiation for the roof is independent of height to width ratio and solar zenith angle. At both solar zenith angles, the absorbed solar radiation for the road decreases rapidly with increasing height to width ratio as the buildings shade more of the road. The shaded wall absorbs less solar radiation than the sunlit wall because it receives only diffuse radiation from the sun and reflected radiation from the walls and road. The sunlit wall absorbs more solar radiation at larger solar zenith angles for height to width ratios less than about three because the incidence angle of the radiation is closer to zero (Figure 2.1). The sum of the absorbed solar radiation for road, sunlit wall, and shaded wall, after converting the wall fluxes to per unit ground area, is the canyon absorbed solar radiation. The absorbed solar radiation for the canyon increases slowly with increasing height to width ratio.

Figure 2.5. Solar radiation absorbed by urban surfaces for solar zenith angles of 30º (top) and 60º (bottom). The atmospheric solar radiation is :math:`{S_{atm}}\, \downarrow _\Lambda ^\mu = 400` and :math:`{S_{atm}}\, \downarrow _\Lambda ^{} = 200` W m-2. Note that the sunlit and shaded wall fluxes are per unit wall area. The solar radiation absorbed by the canyon is the sum of road and wall fluxes after converting the walls fluxes to per unit ground area using the height to width ratio.

The canyon albedo (excluding the roof albedo) shown in Figure 2.6 has the same functional relationships with solar zenith angle and height to width ratio as TEB (Masson 2000). In general, the direct and diffuse canyon albedo decreases with height to width ratio as more solar radiation is trapped and absorbed within the canyon. The trapping of solar radiation is less effective at larger solar zenith angles. At these large solar zenith angles and small height to width ratio, the albedo increases because the higher albedo walls dominate the radiative exchange.

Figure 2.6. Direct beam and diffuse albedo of the urban canyon (walls and road) as a function of height to width ratio from 0.1 to 3.0 in increments of 0.1 and solar zenith angles from 0º to 85º in increments of 5º. The atmospheric solar radiation is :math:`{S_{atm}}\, \downarrow _\Lambda ^\mu = 400` and :math:`{S_{atm}}\, \downarrow _\Lambda ^{} = 200` W m-2.

Incident longwave radiation
-------------------------------

Similar to incident diffuse solar radiation, the longwave radiation incident on walls and roads depends on view factors. The longwave radiation incident on roof, walls and road is

 :math:`{L_{roof}}\, \downarrow = {L_{atm}} \downarrow` ,

 :math:`{L_{imprvrd}}\, \downarrow = {L_{prvrd}}\, \downarrow = {L_{atm}}\, \downarrow {\Psi _{sky - road}}` ,

 :math:`{L_{shdwall}}\, \downarrow = {L_{atm}}\, \downarrow {\Psi _{sky - wall}}` ,

 :math:`{L_{sunwall}}\, \downarrow = {L_{atm}}\, \downarrow {\Psi _{sky - wall}}` 

where :math:`{L_{atm}}\, \downarrow` is the longwave radiation from the atmosphere. The incident longwave radiation conserves energy as

$\begin{gathered}

{L_{atm}}\, \downarrow = {f_{roof}}{L_{roof}}\, \downarrow + \left( {1 - {f_{roof}}} \right) \hfill \\

& \left[ {{L_{imprvrd}}\, \downarrow \left( {1 - {f_{prvrd}}} \right) + {L_{prvrd}}\, \downarrow {f_{prvrd}} + \frac{H}{W}\left( {{L_{sunwall}}\, \downarrow + {L_{shdwall}}\, \downarrow } \right)} \right] \hfill \\

\end{gathered} $.

Absorbed, reflected, and emitted longwave radiation
-------------------------------------------------------

Emitted longwave radiation, a function of surface temperature and emissivity, must also be considered in addition to reflection and absorption when determining the longwave interactions within the canyon. The net longwave radiation (W m-2) (positive toward the atmosphere) for the roof is simply



.. math::

 \vec L_{roof}^{} = {L_{roof}}\, \uparrow - {L_{atm}} \downarrow



where



.. math::

 {L_{roof}}\, \uparrow = {\varepsilon _{roof}}\sigma {\left( {{T_{roof}}} \right)^4} + \left( {1 - {\varepsilon _{roof}}} \right){L_{atm}} \downarrow



is the emitted plus reflected longwave radiation from the roof, :math:`{\varepsilon _{roof}}` is the emissivity of the roof, :math:`\sigma` is the Stefan-Boltzmann constant (W m-2 K-4) (Table 1.4), and :math:`{T_{roof}}` is the temperature of the roof (section 4).

Similar to albedo, the emissivity of each urban surface is a weighted combination of snow-free surface and snow emissivity. Only roof and road surfaces are affected by snow as



.. math::

 {\varepsilon _u} = {\varepsilon _u}\left( {1 - {f_{u,\,sno}}} \right) + {\varepsilon _{sno}}{f_{u,\,sno}}



where :math:`{\varepsilon _u}` is the emissivity of :math:`u =` roof, pervious and impervious road (Table 1.3), :math:`{\varepsilon _{sno}} = 0.97` is the emissivity of snow (Oleson et al. 2004), and :math:`{f_{u,\,sno}}` is the fraction of the urban surface covered with snow (equation ).

As with solar radiation, the longwave interactions within the urban canyon are determined numerically by allowing for multiple reflections until a convergence criteria is met (the absorbed longwave radiation for a given reflection is less than :math:`1 \times {10^{ - 3}}` ). The following equations assume that absorptivity equals emissivity.

The initial reflected ( :math:`r` ) longwave radiation from each urban surface is



.. math::

 L_{imprvrd,\,i = 0}^{}\mathop \uparrow \limits^r = \left( {1 - {\varepsilon _{imprvrd}}} \right){L_{imprvrd}} \downarrow

,



.. math::

 L_{prvrd,\,i = 0}^{}\mathop \uparrow \limits^r = \left( {1 - {\varepsilon _{prvrd}}} \right){L_{prvrd}} \downarrow

,



.. math::

 L_{road,\,i = 0}^{}\mathop \uparrow \limits^r = L_{imprvrd,\,i = 0}^{}\mathop \uparrow \limits^r \left( {1 - {f_{prvrd}}} \right) + L_{prvrd,\,i = 0}^{}\mathop \uparrow \limits^r {f_{prvrd}}





.. math::

 L_{sunwall,\,i = 0}^{}\mathop \uparrow \limits^r = \left( {1 - {\varepsilon _{wall}}} \right){L_{sunwall}} \downarrow

,



.. math::

 L_{shdwall,\,i = 0}^{}\mathop \uparrow \limits^r = \left( {1 - {\varepsilon _{wall}}} \right){L_{shdwall}} \downarrow

.

The emitted ( :math:`e` ) longwave radiation from each surface is



.. math::

 L_{imprvrd}^{}\mathop \uparrow \limits^e = {\varepsilon _{imprvrd}}\sigma {\left( {{T_{imprvrd}}} \right)^4}

,



.. math::

 L_{prvrd}^{}\mathop \uparrow \limits^e = {\varepsilon _{prvrd}}\sigma {\left( {{T_{prvrd}}} \right)^4}

,



.. math::

 L_{road}^{}\mathop \uparrow \limits^e = {\varepsilon _{imprvrd}}\sigma {\left( {{T_{imprvrd}}} \right)^4}\left( {1 - {f_{prvrd}}} \right) + {\varepsilon _{prvrd}}\sigma {\left( {{T_{prvrd}}} \right)^4}{f_{prvrd}}

,



.. math::

 L_{sunwall}^{}\mathop \uparrow \limits^e = {\varepsilon _{wall}}\sigma {\left( {{T_{sunwall}}} \right)^4}

,



.. math::

 L_{shdwall}^{}\mathop \uparrow \limits^e = {\varepsilon _{wall}}\sigma {\left( {{T_{shdwall}}} \right)^4}

.

The initial reflected longwave radiation is distributed to sky, walls, and road according to view factors as



.. math::

 L_{imprvrd - sky,\,i = 0}^{}\mathop \uparrow \limits^r = L_{imprvrd,\,i = 0}^{}\mathop \uparrow \limits^r {\Psi _{road - sky}}

,



.. math::

 L_{prvrd - sky,\,i = 0}^{}\mathop \uparrow \limits^r = L_{prvrd,\,i = 0}^{}\mathop \uparrow \limits^r {\Psi _{road - sky}}

,



.. math::

 L_{road - sunwall,\,i = 0}^{}\mathop \uparrow \limits^r = L_{road,\,i = 0}^{}\mathop \uparrow \limits^r {\Psi _{road - wall}}

,



.. math::

 L_{road - shdwall,\,i = 0}^{}\mathop \uparrow \limits^r = L_{road,\,i = 0}^{}\mathop \uparrow \limits^r {\Psi _{road - wall}}

,



.. math::

 L_{sunwall - sky,\,i = 0}^{}\mathop \uparrow \limits^r = L_{sunwall,\,i = 0}^{}\mathop \uparrow \limits^r {\Psi _{wall - sky}}

,



.. math::

 L_{sunwall - road,\,i = 0}^{}\mathop \uparrow \limits^r = L_{sunwall,\,i = 0}^{}\mathop \uparrow \limits^r {\Psi _{wall - road}}

,



.. math::

 L_{sunwall - shdwall,\,i = 0}^{}\mathop \uparrow \limits^r = L_{sunwall,\,i = 0}^{}\mathop \uparrow \limits^r {\Psi _{wall - wall}}

,



.. math::

 L_{shdwall - sky,\,i = 0}^{}\mathop \uparrow \limits^r = L_{shdwall,\,i = 0}^{}\mathop \uparrow \limits^r {\Psi _{wall - sky}}

,



.. math::

 L_{shdwall - road,\,i = 0}^{}\mathop \uparrow \limits^r = L_{shdwall,\,i = 0}^{}\mathop \uparrow \limits^r {\Psi _{wall - road}}

,



.. math::

 L_{shdwall - sunwall,\,i = 0}^{}\mathop \uparrow \limits^r = L_{shdwall,\,i = 0}^{}\mathop \uparrow \limits^r {\Psi _{wall - wall}}

.

The emitted longwave radiation is distributed to sky, walls, and road according to view factors as



.. math::

 L_{imprvrd - sky}^{}\mathop \uparrow \limits^e = L_{imprvrd}^{}\mathop \uparrow \limits^e {\Psi _{road - sky}}

,



.. math::

 L_{prvrd - sky}^{}\mathop \uparrow \limits^e = L_{prvrd}^{}\mathop \uparrow \limits^e {\Psi _{road - sky}}

,



.. math::

 L_{road - sunwall}^{}\mathop \uparrow \limits^e = L_{road}^{}\mathop \uparrow \limits^e {\Psi _{road - wall}}

,



.. math::

 L_{road - shdwall}^{}\mathop \uparrow \limits^e = L_{road}^{}\mathop \uparrow \limits^e {\Psi _{road - wall}}

,



.. math::

 L_{sunwall - sky}^{}\mathop \uparrow \limits^e = L_{sunwall}^{}\mathop \uparrow \limits^e {\Psi _{wall - sky}}

,



.. math::

 L_{sunwall - road}^{}\mathop \uparrow \limits^e = L_{sunwall}^{}\mathop \uparrow \limits^e {\Psi _{wall - road}}

,



.. math::

 L_{sunwall - shdwall}^{}\mathop \uparrow \limits^e = L_{sunwall}^{}\mathop \uparrow \limits^e {\Psi _{wall - wall}}

,

\[\]	

.. math::

 L_{shdwall - sky}^{}\mathop \uparrow \limits^e = L_{shdwall}^{}\mathop \uparrow \limits^e {\Psi _{wall - sky}}

,



.. math::

 L_{shdwall - road}^{}\mathop \uparrow \limits^e = L_{shdwall}^{}\mathop \uparrow \limits^e {\Psi _{wall - road}}

,



.. math::

 L_{shdwall - sunwall}^{}\mathop \uparrow \limits^e = L_{shdwall}^{}\mathop \uparrow \limits^e {\Psi _{wall - wall}}

.

The initial absorption (net longwave) ( :math:`i = 0` ) by each urban surface is



.. math::

 \vec L_{imprvrd,\,i = 0}^{} = {L_{imprvrd}}\mathop \uparrow \limits^e - {\varepsilon _{imprvrd}}{L_{imprvrd}} \downarrow

,



.. math::

 \vec L_{prvrd,\,i = 0}^{} = {L_{prvrd}}\mathop \uparrow \limits^e - {\varepsilon _{prvrd}}{L_{prvrd}} \downarrow

,



.. math::

 \vec L_{sunwall,\,i = 0}^{} = {L_{sunwall}}\mathop \uparrow \limits^e - {\varepsilon _{wall}}{L_{sunwall}} \downarrow

,



.. math::

 \vec L_{shdwall,\,i = 0}^{} = {L_{shdwall}}\mathop \uparrow \limits^e - {\varepsilon _{wall}}{L_{shdwall}} \downarrow

.

The initial emitted plus reflected longwave radiation to the sky is



.. math::

 L_{imprvrd - sky,\,i = 0}^{}\mathop \uparrow \limits^{} = L_{imprvrd - sky}^{}\mathop \uparrow \limits^e + L_{imprvrd - sky,\,i = 0}^{}\mathop \uparrow \limits^r

,



.. math::

 L_{prvrd - sky,\,i = 0}^{}\mathop \uparrow \limits^{} = L_{prvrd - sky}^{}\mathop \uparrow \limits^e + L_{prvrd - sky,\,i = 0}^{}\mathop \uparrow \limits^r

,



.. math::

 L_{sunwall - sky,\,i = 0}^{}\mathop \uparrow \limits^{} = L_{sunwall - sky}^{}\mathop \uparrow \limits^e + L_{sunwall - sky,\,i = 0}^{}\mathop \uparrow \limits^r

,



.. math::

 L_{shdwall - sky,\,i = 0}^{}\mathop \uparrow \limits^{} = L_{shdwall - sky}^{}\mathop \uparrow \limits^e + L_{shdwall - sky,\,i = 0}^{}\mathop \uparrow \limits^r

.

The net longwave radiation absorbed by each urban surface after the :math:`{i^{th}}` reflection is

\[\vec L_{imprvrd,\,i}^{} = {\varepsilon _{imprvrd}}\left( \begin{gathered}

{L_{sunwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - road}}\mathop \uparrow \limits^e \hfill \\

+ {L_{shdwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - road}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right)\frac{H}{W}\],

\[\vec L_{prvrd,\,i}^{} = {\varepsilon _{prvrd}}\left( \begin{gathered}

{L_{sunwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - road}}\mathop \uparrow \limits^e \hfill \\

+ {L_{shdwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - road}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right)\frac{H}{W}\],



.. math::

 \vec L_{road,\,i}^{} = \vec L_{imprvrd,\,i}^{}\left( {1 - {f_{prvrd}}} \right) + \vec L_{prvrd,\,i}^{}{f_{prvrd}}

,

\[\vec L_{sunwall,\,i}^{} = {\varepsilon _{wall}}\left( \begin{gathered}

\frac{{{L_{road - sunwall,\,i - 1}}\mathop \uparrow \limits^r + {L_{road - sunwall}}\mathop \uparrow \limits^e }}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}} \hfill \\

+ {L_{shdwall - sunwall,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - sunwall}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right)\],

\[\vec L_{shdwall,\,i}^{} = {\varepsilon _{wall}}\left( \begin{gathered}

\frac{{{L_{road - shdwall,\,i - 1}}\mathop \uparrow \limits^r + {L_{road - shdwall}}\mathop \uparrow \limits^e }}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}} \hfill \\

+ {L_{sunwall - shdwall,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - shdwall}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right)\].

The longwave radiation from each urban surface after the :math:`{i^{th}}` reflection is distributed to sky, road, and walls according to

\[L_{imprvrd - sky,\,i}^{} \uparrow = \left( {1 - {\varepsilon _{imprvrd}}} \right)\left( \begin{gathered}

{L_{sunwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - road}}\mathop \uparrow \limits^e \hfill \\

+ {L_{shdwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - road}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right)\frac{H}{W}{\Psi _{road - sky}}\],

\[L_{imprvrd - sunwall,\,i}^{} \uparrow = \left( {1 - {\varepsilon _{imprvrd}}} \right)\left( \begin{gathered}

{L_{sunwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - road}}\mathop \uparrow \limits^e \hfill \\

+ {L_{shdwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - road}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right)\frac{H}{W}{\Psi _{road - wall}}\],

\[L_{imprvrd - shdwall,\,i}^{} \uparrow = \left( {1 - {\varepsilon _{imprvrd}}} \right)\left( \begin{gathered}

{L_{sunwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - road}}\mathop \uparrow \limits^e \hfill \\

+ {L_{shdwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - road}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right)\frac{H}{W}{\Psi _{road - wall}}\],

\[L_{prvrd - sky,\,i}^{} \uparrow = \left( {1 - {\varepsilon _{prvrd}}} \right)\left( \begin{gathered}

{L_{sunwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - road}}\mathop \uparrow \limits^e \hfill \\

+ {L_{shdwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - road}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right)\frac{H}{W}{\Psi _{road - sky}}\],

\[L_{prvrd - sunwall,\,i}^{} \uparrow = \left( {1 - {\varepsilon _{prvrd}}} \right)\left( \begin{gathered}

{L_{sunwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - road}}\mathop \uparrow \limits^e \hfill \\

+ {L_{shdwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - road}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right)\frac{H}{W}{\Psi _{road - wall}}\],

\[L_{prvrd - shdwall,\,i}^{} \uparrow = \left( {1 - {\varepsilon _{prvrd}}} \right)\left( \begin{gathered}

{L_{sunwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - road}}\mathop \uparrow \limits^e \hfill \\

+ {L_{shdwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - road}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right)\frac{H}{W}{\Psi _{road - wall}}\],

\[L_{road - sky,\,i}^{} \uparrow = \left[ \begin{gathered}

\left( {1 - {\varepsilon _{imprvrd}}} \right)\left( \begin{gathered}

{L_{sunwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - road}}\mathop \uparrow \limits^e \hfill \\

+ {L_{shdwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - road}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right) \hfill \\

\times \frac{H}{W}\left( {1 - {f_{prvrd}}} \right) + \left( {1 - {\varepsilon _{prvrd}}} \right) \hfill \\

\times \left( \begin{gathered}

{L_{sunwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - road}}\mathop \uparrow \limits^e \hfill \\

+ {L_{shdwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - road}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right)\frac{H}{W}{f_{prvrd}} \hfill \\

\end{gathered} \right]{\Psi _{road - sky}}\],

\[L_{road - sunwall,\,i}^{} \uparrow = \left[ \begin{gathered}

\left( {1 - {\varepsilon _{imprvrd}}} \right)\left( \begin{gathered}

{L_{sunwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - road}}\mathop \uparrow \limits^e \hfill \\

+ {L_{shdwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - road}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right) \hfill \\

\times \frac{H}{W}\left( {1 - {f_{prvrd}}} \right) + \left( {1 - {\varepsilon _{prvrd}}} \right) \hfill \\

\times \left( \begin{gathered}

{L_{sunwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - road}}\mathop \uparrow \limits^e \hfill \\

+ {L_{shdwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - road}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right)\frac{H}{W}{f_{prvrd}} \hfill \\

\end{gathered} \right]{\Psi _{road - wall}}\],

\[L_{road - shdwall,\,i}^{} \uparrow = \left[ \begin{gathered}

\left( {1 - {\varepsilon _{imprvrd}}} \right)\left( \begin{gathered}

{L_{sunwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - road}}\mathop \uparrow \limits^e \hfill \\

+ {L_{shdwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - road}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right) \hfill \\

\times \frac{H}{W}\left( {1 - {f_{prvrd}}} \right) + \left( {1 - {\varepsilon _{prvrd}}} \right) \hfill \\

\times \left( \begin{gathered}

{L_{sunwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - road}}\mathop \uparrow \limits^e \hfill \\

+ {L_{shdwall - road,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - road}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right)\frac{H}{W}{f_{prvrd}} \hfill \\

\end{gathered} \right]{\Psi _{road - wall}}\],

\[L_{sunwall - sky,\,i}^{} \uparrow = \left( {1 - {\varepsilon _{wall}}} \right)\left( \begin{gathered}

\frac{{{L_{road - sunwall,\,i - 1}}\mathop \uparrow \limits^r + {L_{road - sunwall}}\mathop \uparrow \limits^e }}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}} \hfill \\

+ {L_{shdwall - sunwall,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - sunwall}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right){\Psi _{wall - sky}}\],

\[L_{sunwall - road,\,i}^{} \uparrow = \left( {1 - {\varepsilon _{wall}}} \right)\left( \begin{gathered}

\frac{{{L_{road - sunwall,\,i - 1}}\mathop \uparrow \limits^r + {L_{road - sunwall}}\mathop \uparrow \limits^e }}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}} \hfill \\

+ {L_{shdwall - sunwall,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - sunwall}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right){\Psi _{wall - road}}\],

\[L_{sunwall - shdwall,\,i}^{} \uparrow = \left( {1 - {\varepsilon _{wall}}} \right)\left( \begin{gathered}

\frac{{{L_{road - sunwall,\,i - 1}}\mathop \uparrow \limits^r + {L_{road - sunwall}}\mathop \uparrow \limits^e }}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}} \hfill \\

+ {L_{shdwall - sunwall,\,i - 1}}\mathop \uparrow \limits^r + {L_{shdwall - sunwall}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right){\Psi _{wall - wall}}\],

\[L_{shdwall - sky,\,i}^{} \uparrow = \left( {1 - {\varepsilon _{wall}}} \right)\left( \begin{gathered}

\frac{{{L_{road - shdwall,\,i - 1}}\mathop \uparrow \limits^r + {L_{road - shdwall}}\mathop \uparrow \limits^e }}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}} \hfill \\

+ {L_{sunwall - shdwall,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - shdwall}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right){\Psi _{wall - sky}}\],

\[L_{shdwall - road,\,i}^{} \uparrow = \left( {1 - {\varepsilon _{wall}}} \right)\left( \begin{gathered}

\frac{{{L_{road - shdwall,\,i - 1}}\mathop \uparrow \limits^r + {L_{road - shdwall}}\mathop \uparrow \limits^e }}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}} \hfill \\

+ {L_{sunwall - shdwall,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - shdwall}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right){\Psi _{wall - road}}\],

\[L_{shdwall - sunwall,\,i}^{} \uparrow = \left( {1 - {\varepsilon _{wall}}} \right)\left( \begin{gathered}

\frac{{{L_{road - shdwall,\,i - 1}}\mathop \uparrow \limits^r + {L_{road - shdwall}}\mathop \uparrow \limits^e }}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}} \hfill \\

+ {L_{sunwall - shdwall,\,i - 1}}\mathop \uparrow \limits^r + {L_{sunwall - shdwall}}\mathop \uparrow \limits^e \hfill \\

\end{gathered} \right){\Psi _{wall - wall}}\].

Note that the emitted longwave term in equations - only applies to the first iteration. Subsequent iterations do not include this term, i.e.,

\[\begin{gathered}

{L_{road - sunwall}}\mathop \uparrow \limits^e = {L_{road - shdwall}}\mathop \uparrow \limits^e = {L_{sunwall - road}}\mathop \uparrow \limits^e = {L_{shdwall - road}}\mathop \uparrow \limits^e \hfill \\

= {L_{shdwall - sunwall}}\mathop \uparrow \limits^e = {L_{sunwall - shdwall}}\mathop \uparrow \limits^e = 0 \hfill \\

\end{gathered} \].

The reflected longwave radiation to the sky is added to the total upward longwave radiation for each urban surface as



.. math::

 {L_{imprvrd,\,i + 1}} \uparrow = {L_{imprvrd,\,i - 1}} \uparrow + {L_{imprvrd - sky,\,i}} \uparrow

,



.. math::

 {L_{prvrd,\,i + 1}} \uparrow = {L_{prvrd,\,i - 1}} \uparrow + {L_{prvrd - sky,\,i}} \uparrow

,



.. math::

 {L_{sunwall,\,i + 1}} \uparrow = {L_{sunwall,\,i - 1}} \uparrow + {L_{sunwall - sky,\,i}} \uparrow

,



.. math::

 {L_{shdwall,\,i + 1}} \uparrow = {L_{shdwall,\,i - 1}} \uparrow + {L_{shdwall - sky,\,i}} \uparrow

.

The net longwave at each iteration is added to the total net longwave for each urban surface as



.. math::

 {\vec L_{imprvrd,\,i + 1}} = {\vec L_{imprvrd,\,i - 1}} + {\vec L_{imprvrd,\,i}}

,



.. math::

 {\vec L_{prvrd,\,i + 1}} = {\vec L_{prvrd,\,i - 1}} + {\vec L_{prvrd,\,i}}

,



.. math::

 {\vec L_{sunwall,\,i + 1}} = {\vec L_{sunwall,\,i - 1}} + {\vec L_{sunwall,\,i}}

,



.. math::

 {\vec L_{shdwall,\,i + 1}} = {\vec L_{shdwall,\,i - 1}} + {\vec L_{shdwall,\,i}}

.

The system of equations (equations -) is iterated for :math:`i = 50` reflections or until the absorption for the :math:`{i^{th}}` reflection is less than a nominal amount



.. math::

 \max \left( {\vec L_{road,\,i}^{},\vec L_{sunwall,\,i}^{},\vec L_{shdwall,\,i}^{}} \right) < 1 \times {10^{ - 3}}

.

The net longwave radiation for the urban canyon (walls and road) is



.. math::

 \vec L_{uc}^{} = \vec L_{imprvrd,\,n + 1}^{}\left( {1 - {f_{prvrd}}} \right) + \vec L_{prvrd,\,n + 1}^{}{f_{prvrd}} + \left( {\vec L_{sunwall,\,n + 1}^{} + \vec L_{shdwall,\,n + 1}^{}} \right)\frac{H}{W}

.

while the total reflected plus emitted longwave radiation is

\[\begin{gathered}

{L_{uc}} \uparrow = {L_{imprvrd,\,n + 1}} \uparrow \left( {1 - {f_{prvrd}}} \right) + {L_{prvrd,\,n + 1}} \uparrow {f_{prvrd}} \hfill \\

& + \left( {{L_{sunwall,\,n + 1}} \uparrow + {L_{shdwall,\,n + 1}} \uparrow } \right)\frac{H}{W} \hfill \\

\end{gathered} \].

Longwave radiation in the urban canyon is conserved as



.. math::

 {\vec L_{uc}} - \left( {{L_{uc}} \uparrow - {L_{atm}} \downarrow } \right) = 0

.

The total net longwave radiation for the urban canopy (road, walls, and roof) is



.. math::

 \vec L = {W_{roof}}\vec L_{roof}^{} + \left( {1 - {W_{roof}}} \right)\vec L_{uc}^{}



Figure 2.7 shows the net longwave radiation for urban surfaces for two different emissivity configurations. A positive net longwave means that the outgoing longwave exceeds the incoming longwave from the atmosphere. The net longwave radiation for the roof is independent of height to width ratio and increases with higher emissivity. The net longwave radiation for the road and walls decreases rapidly with increasing height to width ratio as more longwave radiation is trapped within the canyon. The walls have lower net longwave radiation than the road because their sky view factors are smaller. The two walls behave identically with respect to net longwave radiation as long as temperatures are the same. The sum of the net longwave radiation for road, sunlit wall, and shaded wall, after converting the wall fluxes to per unit ground area, is the canyon net longwave radiation. The net longwave radiation for the canyon increases slowly with increasing height to width ratio because of the larger surface area of the walls.

Figure 2.7. Net longwave radiation (positive to the atmosphere) for urban surfaces for two different emissivity configurations. The atmospheric longwave radiation is :math:`{L_{atm}}\, \downarrow = 340` W m-2 and the temperature of each surface is 292.16 K. Note that the wall fluxes (shaded and sunlit) are per unit wall area. The net longwave radiation for the canyon is the sum of road and wall fluxes after converting the walls fluxes to per unit ground area using the height to width ratio.

Solar Zenith Angle
----------------------

The formulation for solar zenith angle is thoroughly documented in Oleson et al. (2010b) (see section 3.3) and does not differ for urban surfaces.

Heat and Momentum Fluxes
===========================

The net radiation for the urban canopy ( :math:`\vec S - \vec L` , where :math:`\vec S` is the net solar radiation absorbed by the urban canopy (section 2.5) and :math:`\vec L` is the net longwave radiation (section 2.7)) must be balanced by the sum of the turbulent and ground (storage) heat fluxes as



.. math::

 \overrightarrow S - \vec L = H + \lambda E + G



where :math:`H` is the sensible heat flux (W m-2), :math:`E` is the water vapor flux (kg m-2 s-1), :math:`G` is the ground heat flux, and :math:`\lambda` is the latent heat of vaporization (or sublimation). The urban surfaces have unique radiative, thermal and hydrologic properties and environments. Thus, their sensible and latent heat fluxes are likely to be very different from each other. For example, the pervious road may have significant latent heat flux compared to the walls, which are assumed to be hydrologically inactive. Thus, the fluxes from individual urban surfaces must be modeled separately. However, CLM directly interacts with the atmospheric model at only the lowest atmospheric layer, which is well above the roof level of the urban model at the horizontal scales to be modeled. As a consequence, fluxes from individual urban surfaces must be combined to obtain the total sensible and latent heat flux to be provided to the atmospheric model. Allowing the urban surface fluxes to interact with each other through a bulk urban air mass is an acceptable approach analogous to the simulation of vegetated canopy fluxes (Figure 3.1). This also allows for the solution of UCL air temperature and humidity, which are of interest in many applications. The approach shown in Figure 3.1 is slightly different from that of Masson (2000) in that here, fluxes from the roof interact directly with the UCL air whereas in Masson (2000) the roof and urban canyon are modeled as two independent sources of heat and moisture fluxes to the atmosphere. Here, we assume that the actual roofs are at various heights in the UCL and hence interact directly with the well-mixed UCL air.

Figure 3.1. Schematic diagram of sensible and latent heat fluxes for the urban canopy.

In general, the zonal :math:`{\tau _x}` and meridional :math:`{\tau _y}` momentum fluxes (kg m-1 s-2), sensible heat flux :math:`H` , and water vapor flux :math:`E` between the atmosphere at reference height :math:`{z_{atm,\,x}}` (m) [where :math:`x` is height for wind (momentum) ( :math:`m` ), temperature (sensible heat) ( :math:`h` ), and humidity (water vapor) ( :math:`w` ); with zonal and meridional winds :math:`{u_{atm}}` and :math:`{v_{atm}}` (m s-1), potential temperature :math:`{\theta _{atm}}` (K), and specific humidity :math:`{q_{atm}}` (kg kg-1)] and a surface [with :math:`{u_s}` , :math:`{v_s}` , :math:`{\theta _s}` , and :math:`{q_s}` ] are



.. math::

 {\tau _x} = - {\rho _{atm}}\frac{{\left( {{u_{atm}} - {u_s}} \right)}}{{{r_{am}}}}





.. math::

 {\tau _y} = - {\rho _{atm}}\frac{{\left( {{v_{atm}} - {v_s}} \right)}}{{{r_{am}}}}





.. math::

 H = - {\rho _{atm}}{C_p}\frac{{\left( {{\theta _{atm}} - {\theta _s}} \right)}}{{{r_{ah}}}}





.. math::

 E = - {\rho _{atm}}\frac{{\left( {{q_{atm}} - {q_s}} \right)}}{{{r_{aw}}}}

.

These fluxes are derived from Monin-Obukhov similarity theory developed for the inertial sub-layer (i.e., the nearly constant flux layer above the roughness sub-layer). In this derivation, :math:`{u_s}` and :math:`{v_s}` are defined to equal zero at height :math:`{z_{0m}} + d` (the apparent sink for momentum) so that :math:`{r_{am}}` is the aerodynamic resistance (s m-1) for momentum between the atmosphere at height :math:`{z_{atm,\,m}}` and the surface at height :math:`{z_{0m}} + d` . Thus, the momentum fluxes become



.. math::

 {\tau _x} = - {\rho _{atm}}\frac{{{u_{atm}}}}{{{r_{am}}}}





.. math::

 {\tau _y} = - {\rho _{atm}}\frac{{{v_{atm}}}}{{{r_{am}}}}

.

Likewise, :math:`{\theta _s}` and :math:`{q_s}` are defined at heights :math:`{z_{0h}} + d` and :math:`{z_{0w}} + d` (the apparent sinks for heat and water vapor, respectively). Consequently, :math:`{r_{ah}}` and :math:`{r_{aw}}` are the aerodynamic resistances (s m-1) to sensible heat and water vapor transfer between the atmosphere at heights :math:`{z_{atm,\,h}}` and :math:`{z_{atm,\,w}}` and the surface at heights :math:`{z_{0h}} + d` and :math:`{z_{0w}} + d` , respectively. The specific heat capacity of air :math:`{C_p}` (J kg-1 K-1) is a constant (Table 1.4). The atmospheric potential temperature used here is



.. math::

 {\theta _{atm}} = {T_{atm}} + {\Gamma _d}{z_{atm,\,h}}



where :math:`{T_{atm}}` is the air temperature (K) at height :math:`{z_{atm,\,h}}` and :math:`{\Gamma _d} = 0.0098` K m-1 is the negative of the dry adiabatic lapse rate [this expression is first-order equivalent to ${\theta _{atm}} = {T_{atm}}{\left( {{{{P_{srf}}} \mathord{\left/

{\vphantom {{{P_{srf}}} {{P_{atm}}}}} \right.

\kern-\nulldelimiterspace} {{P_{atm}}}}} \right)^{{{{R_{da}}} \mathord{\left/

{\vphantom {{{R_{da}}} {{C_p}}}} \right.

\kern-\nulldelimiterspace} {{C_p}}}}} :math:`(Stull 1988), where` {P_{srf}} :math:`is the surface pressure (Pa),` {P_{atm}} :math:`is the atmospheric pressure (Pa), and` {R_{da}} :math:`is the gas constant for dry air (J kg-1 K-1) (Table 1.4)]. By definition,` {\theta _s} = {T_s}$. The density of moist air (kg m-3) is



.. math::

 {\rho _{atm}} = \frac{{{P_{atm}} - 0.378{e_{atm}}}}{{{R_{da}}{T_{atm}}}}



where the atmospheric vapor pressure :math:`{e_{atm}}` (Pa) is derived from the atmospheric specific humidity :math:`{q_{atm}}` 



.. math::

 {e_{atm}} = \frac{{{q_{atm}}{P_{atm}}}}{{0.622 + 0.378{q_{atm}}}}

.

Monin-Obukhov Similarity Theory
-----------------------------------

The surface vertical kinematic fluxes of momentum 

.. math::

 \overline {u'w'}

 and 

.. math::

 \overline {v'w'}

 (m2 s-2), sensible heat 

.. math::

 \overline {\theta 'w'}

 (K m s-1), and latent heat 

.. math::

 \overline {q'w'}

 (kg kg-1 m s-1), where :math:`u'` , :math:`v'` , :math:`w'` , :math:`\theta '` , and :math:`q'` are zonal horizontal wind, meridional horizontal wind, vertical velocity, potential temperature, and specific humidity turbulent fluctuations about the mean, are defined from Monin-Obukhov similarity applied to the surface layer. This theory states that when scaled appropriately, the dimensionless mean horizontal wind speed, mean potential temperature, and mean specific humidity profile gradients depend on unique functions of :math:`\zeta = \frac{{z - d}}{L}` (Zeng et al. 1998) as



.. math::

 \frac{{k\left( {z - d} \right)}}{{{u_ * }}}\frac{{\partial \left| u \right|}}{{\partial z}} = {\phi _m}\left( \zeta \right)





.. math::

 \frac{{k\left( {z - d} \right)}}{{{\theta _ * }}}\frac{{\partial \theta }}{{\partial z}} = {\phi _h}\left( \zeta \right)





.. math::

 \frac{{k\left( {z - d} \right)}}{{{q_ * }}}\frac{{\partial q}}{{\partial z}} = {\phi _w}\left( \zeta \right)



where :math:`z` is height in the surface layer (m), :math:`d` is the displacement height (m), :math:`L` is the Monin-Obukhov length scale (m) that accounts for buoyancy effects resulting from vertical density gradients (i.e., the atmospheric stability), k is the von Karman constant (Table 1.4), and :math:`\left| u \right|` is the atmospheric wind speed (m s-1). :math:`{\phi _m}` , :math:`{\phi _h}` , and :math:`{\phi _w}` are universal (over any surface) similarity functions of :math:`\zeta` that relate the constant fluxes of momentum, sensible heat, and latent heat to the mean profile gradients of :math:`\left| u \right|` , :math:`\theta` , and :math:`q` in the surface layer. In neutral conditions, :math:`{\phi _m} = {\phi _h} = {\phi _w} = 1` . The velocity (i.e., friction velocity) :math:`{u_ * }` (m s-1), temperature :math:`{\theta _ * }` (K), and moisture :math:`{q_ * }` (kg kg-1) scales are



.. math::

 u_ * ^2 = \sqrt {{{\left( {\overline {u'w'} } \right)}^2} + {{\left( {\overline {v'w'} } \right)}^2}} = \frac{{\left| \tau \right|}}{{{\rho _{atm}}}}





.. math::

 {\theta _ * }{u_ * } = - \overline {\theta 'w'} = - \frac{H}{{{\rho _{atm}}{C_p}}}



 :math:`{q_ * }{u_ * } = - \overline {q'w'} = - \frac{E}{{{\rho _{atm}}}}` 

where :math:`\left| \tau \right|` is the shearing stress (kg m-1 s-2), with zonal and meridional components :math:`\overline {u'w'} = - \frac{{{\tau _x}}}{{{\rho _{atm}}}}` and :math:`\overline {v'w'} = - \frac{{{\tau _y}}}{{{\rho _{atm}}}}` , respectively, :math:`H` is the sensible heat flux (W m-2) and :math:`E` is the water vapor flux (kg m-2 s-1).

The dimensionless length scale :math:`L` is the Monin-Obukhov length defined as



.. math::

 L = - \frac{{u_ * ^3}}{{k\left( {\frac{g}{{\overline {{\theta _{v,\,atm}}} }}} \right){{\theta '}_v}w'}} = \frac{{u_ * ^2\overline {{\theta _{v,\,atm}}} }}{{kg{\theta _{v * }}}}



where :math:`g` is the acceleration of gravity (m s-2) (Table 1.4), and :math:`\overline {{\theta _{v,\,atm}}} = \overline {{\theta _{atm}}} \left( {1 + 0.61{q_{atm}}} \right)` is the reference virtual potential temperature. :math:`L > 0` indicates stable conditions. :math:`L < 0` indicates unstable conditions. :math:`L = \infty` for neutral conditions. The temperature scale :math:`{\theta _{v * }}` is defined as



.. math::

 {\theta _{v * }}{u_ * } = \left[ {{\theta _ * }\left( {1 + 0.61{q_{atm}}} \right) + 0.61\overline {{\theta _{atm}}} {q_ * }} \right]{u_ * }



where :math:`\overline {{\theta _{atm}}}` is the atmospheric potential temperature.

Following Panofsky and Dutton (1984), the differential equations for :math:`{\phi _m}\left( \zeta \right)` , :math:`{\phi _h}\left( \zeta \right)` , and :math:`{\phi _w}\left( \zeta \right)` can be integrated formally without commitment to their exact forms. Integration between two arbitrary heights in the surface layer :math:`{z_2}` and :math:`{z_1}` ( :math:`{z_2} > {z_1}` ) with horizontal winds :math:`{\left| u \right|_1}` and :math:`{\left| u \right|_2}` , potential temperatures :math:`{\theta _1}` and :math:`{\theta _2}` , and specific humidities :math:`{q_1}` and :math:`{q_2}` results in



.. math::

 {\left| u \right|_2} - {\left| u \right|_1} = \frac{{{u_ * }}}{k}\left[ {\ln \left( {\frac{{{z_2} - d}}{{{z_1} - d}}} \right) - {\psi _m}\left( {\frac{{{z_2} - d}}{L}} \right) + {\psi _m}\left( {\frac{{{z_1} - d}}{L}} \right)} \right]





.. math::

 {\theta _2} - {\theta _1} = \frac{{{\theta _ * }}}{k}\left[ {\ln \left( {\frac{{{z_2} - d}}{{{z_1} - d}}} \right) - {\psi _h}\left( {\frac{{{z_2} - d}}{L}} \right) + {\psi _h}\left( {\frac{{{z_1} - d}}{L}} \right)} \right]





.. math::

 {q_2} - {q_1} = \frac{{{q_ * }}}{k}\left[ {\ln \left( {\frac{{{z_2} - d}}{{{z_1} - d}}} \right) - {\psi _w}\left( {\frac{{{z_2} - d}}{L}} \right) + {\psi _w}\left( {\frac{{{z_1} - d}}{L}} \right)} \right]

.

The functions :math:`{\psi _m}\left( \zeta \right)` , :math:`{\psi _h}\left( \zeta \right)` , and :math:`{\psi _w}\left( \zeta \right)` are defined as

\[{\psi _m}\left( \zeta \right) = \int_{{{{z_{0m}}} \mathord{\left/

{\vphantom {{{z_{0m}}} L}} \right.

\kern-\nulldelimiterspace} L}}^\zeta {\frac{{\left[ {1 - {\phi _m}\left( x \right)} \right]}}{x}\,dx} \]

\[{\psi _h}\left( \zeta \right) = \int_{{{{z_{0h}}} \mathord{\left/

{\vphantom {{{z_{0h}}} L}} \right.

\kern-\nulldelimiterspace} L}}^\zeta {\frac{{\left[ {1 - {\phi _h}\left( x \right)} \right]}}{x}\,dx} \]

\[{\psi _w}\left( \zeta \right) = \int_{{{{z_{0w}}} \mathord{\left/

{\vphantom {{{z_{0w}}} L}} \right.

\kern-\nulldelimiterspace} L}}^\zeta {\frac{{\left[ {1 - {\phi _w}\left( x \right)} \right]}}{x}\,dx} \]

where :math:`{z_{0m}}` , :math:`{z_{0h}}` , and :math:`{z_{0w}}` are the roughness lengths (m) for momentum, sensible heat, and water vapor, respectively.

Defining the surface values



.. math::

 {\left| u \right|_1} = 0{\text{ at }}{z_1} = {z_{0m}} + d,





.. math::

 {\theta _1} = {\theta _s}{\text{ at }}{z_1} = {z_{0h}} + d,{\text{ and}}





.. math::

 {q_1} = {q_s}{\text{ at }}{z_1} = {z_{0w}} + d,



and the atmospheric values at :math:`{z_2} = {z_{atm,\,x}}` 



.. math::

 {\left| u \right|_2} = {V_a}{\text{ = }}\sqrt {u_{atm}^2 + v_{atm}^2 + U_c^2} \geqslant 1,





.. math::

 {\theta _2} = {\theta _{atm}}{\text{, and}}





.. math::

 {q_2} = {q_{atm}}{\text{, }}



the integral forms of the flux-gradient relations are



.. math::

 {V_a} = \frac{{{u_ * }}}{k}\left[ {\ln \left( {\frac{{{z_{atm,\,m}} - d}}{{{z_{0m}}}}} \right) - {\psi _m}\left( {\frac{{{z_{atm,\,m}} - d}}{L}} \right) + {\psi _m}\left( {\frac{{{z_{0m}}}}{L}} \right)} \right]





.. math::

 {\theta _{atm}} - {\theta _s} = \frac{{{\theta _ * }}}{k}\left[ {\ln \left( {\frac{{{z_{atm,\,h}} - d}}{{{z_{0h}}}}} \right) - {\psi _h}\left( {\frac{{{z_{atm,\,h}} - d}}{L}} \right) + {\psi _h}\left( {\frac{{{z_{0h}}}}{L}} \right)} \right]





.. math::

 {q_{atm}} - {q_s} = \frac{{{q_ * }}}{k}\left[ {\ln \left( {\frac{{{z_{atm,\,w}} - d}}{{{z_{0w}}}}} \right) - {\psi _w}\left( {\frac{{{z_{atm,\,w}} - d}}{L}} \right) + {\psi _w}\left( {\frac{{{z_{0w}}}}{L}} \right)} \right]

.

The constraint 

.. math::

 {V_a} \geqslant 1

 is required simply for numerical reasons to prevent :math:`H` and :math:`E` from becoming small with small wind speeds. The convective velocity :math:`{U_c}` accounts for the contribution of large eddies in the convective boundary layer to surface fluxes as follows

\[\begin{gathered}

{U_c} = 0 & \zeta \geqslant {\text{0}} & {\text{(stable)}} \hfill \\

{U_c} = \beta {w_*} & \zeta < 0 & {\text{(unstable}}) \hfill \\

\end{gathered} \]

where :math:`{w_*}` is the convective velocity scale

\[{w_*} = {\left( {\frac{{ - g{u_*}{\theta_{v*}}{z_i}}}{{\overline {{\theta _{v,\,atm}}} }}} \right)^{{1 \mathord{\left/

{\vphantom {1 3}} \right.

\kern-\nulldelimiterspace} 3}}}\],

 :math:`{z_i} = 1000` is the convective boundary layer height (m), and :math:`\beta = 1` .

The momentum flux gradient relations are (Zeng et al. 1998)

\[\begin{gathered}

{\phi _m}\left( \zeta \right) = 0.7{k^{{2 \mathord{\left/

{\vphantom {2 3}} \right.

\kern-\nulldelimiterspace} 3}}}{\left( { - \zeta } \right)^{{1 \mathord{\left/

{\vphantom {1 3}} \right.

\kern-\nulldelimiterspace} 3}}} & {\text{for }}\zeta < - 1.574{\text{ (very unstable)}} \hfill \\

{\phi _m}\left( \zeta \right) = {\left( {1 - 16\zeta } \right)^{ - {1 \mathord{\left/

{\vphantom {1 4}} \right.

\kern-\nulldelimiterspace} 4}}} & {\text{for - 1}}{\text{.574}} \leqslant \zeta < 0{\text{ (unstable)}} \hfill \\

{\phi _m}\left( \zeta \right) = 1 + 5\zeta & {\text{for }}0 \leqslant \zeta \leqslant 1{\text{ (stable)}} \hfill \\

{\phi _m}\left( \zeta \right) = 5 + \zeta & {\text{for }}\zeta {\text{ > 1 (very stable)}}{\text{.}} \hfill \\

\end{gathered} \]

The sensible and latent heat flux gradient relations are (Zeng et al. 1998)

\[\begin{gathered}

{\phi _h}\left( \zeta \right) = {\phi _w}\left( \zeta \right) = 0.9{k^{{4 \mathord{\left/

{\vphantom {4 3}} \right.

\kern-\nulldelimiterspace} 3}}}{\left( { - \zeta } \right)^{{{ - 1} \mathord{\left/

{\vphantom {{ - 1} 3}} \right.

\kern-\nulldelimiterspace} 3}}} & {\text{for }}\zeta < - 0.465{\text{ (very unstable)}} \hfill \\

{\phi _h}\left( \zeta \right) = {\phi _w}\left( \zeta \right) = {\left( {1 - 16\zeta } \right)^{ - {1 \mathord{\left/

{\vphantom {1 2}} \right.

\kern-\nulldelimiterspace} 2}}} & {\text{for - 0}}{\text{.465}} \leqslant \zeta < 0{\text{ (unstable)}} \hfill \\

{\phi _h}\left( \zeta \right) = {\phi _w}\left( \zeta \right) = 1 + 5\zeta & {\text{for }}0 \leqslant \zeta \leqslant 1{\text{ (stable)}} \hfill \\

{\phi _h}\left( \zeta \right) = {\phi _w}\left( \zeta \right) = 5 + \zeta & {\text{for }}\zeta {\text{ > 1 (very stable)}}{\text{.}} \hfill \\

\end{gathered} \]

To ensure continuous functions of :math:`{\phi _m}\left( \zeta \right)` , :math:`{\phi _h}\left( \zeta \right)` , and :math:`{\phi _w}\left( \zeta \right)` , the simplest approach (i.e., without considering any transition regimes) is to match the relations for very unstable and unstable conditions at :math:`{\zeta _m} = - 1.574` for :math:`{\phi _m}\left( \zeta \right)` and :math:`{\zeta _h} = {\zeta _w} = - 0.465` for :math:`{\phi _h}\left( \zeta \right) = {\phi _w}\left( \zeta \right)` (Zeng et al. 1998). The flux gradient relations can be integrated to yield wind profiles for the following conditions:

Very unstable :math:`\left( {\zeta < - 1.574} \right)` 

\[{V_a} = \frac{{{u_*}}}{k}\left\{ {\left[ {\ln \frac{{{\zeta _m}L}}{{{z_{0m}}}} - {\psi _m}\left( {{\zeta _m}} \right)} \right] + 1.14\left[ {{{\left( { - \zeta } \right)}^{{1 \mathord{\left/

{\vphantom {1 3}} \right.

\kern-\nulldelimiterspace} 3}}} - {{\left( { - {\zeta _m}} \right)}^{{1 \mathord{\left/

{\vphantom {1 3}} \right.

\kern-\nulldelimiterspace} 3}}}} \right] + {\psi _m}\left( {\frac{{{z_{0m}}}}{L}} \right)} \right\}\]

Unstable :math:`\left( { - 1.574 \leqslant \zeta < 0} \right)` 



.. math::

 {V_a} = \frac{{{u_*}}}{k}\left\{ {\left[ {\ln \frac{{{z_{atm,\,m}} - d}}{{{z_{0m}}}} - {\psi _m}\left( \zeta \right)} \right] + {\psi _m}\left( {\frac{{{z_{0m}}}}{L}} \right)} \right\}



Stable :math:`\left( {0 \leqslant \zeta \leqslant 1} \right)` 



.. math::

 {V_a} = \frac{{{u_*}}}{k}\left\{ {\left[ {\ln \frac{{{z_{atm,\,m}} - d}}{{{z_{0m}}}} + 5\zeta } \right] - 5\frac{{{z_{0m}}}}{L}} \right\}



Very stable :math:`\left( {\zeta > 1} \right)` 



.. math::

 {V_a} = \frac{{{u_*}}}{k}\left\{ {\left[ {\ln \frac{L}{{{z_{0m}}}} + 5} \right] + \left[ {5\ln \zeta + \zeta - 1} \right] - 5\frac{{{z_{0m}}}}{L}} \right\}



where



.. math::

 {\psi _m}\left( \zeta \right) = 2\ln \left( {\frac{{1 + x}}{2}} \right) + \ln \left( {\frac{{1 + {x^2}}}{2}} \right) - 2{\tan ^{ - 1}}x + \frac{\pi }{2}



and $x = {\left( {1 - 16\zeta } \right)^{{1 \mathord{\left/

{\vphantom {1 4}} \right.

\kern-\nulldelimiterspace} 4}}}$.

The potential temperature profiles are:

Very unstable :math:`\left( {\zeta < - 0.465} \right)` 

\[{\theta _{atm}} - {\theta _s} = \frac{{{\theta _*}}}{k}\left\{ {\left[ {\ln \frac{{{\zeta _h}L}}{{{z_{0h}}}} - {\psi _h}\left( {{\zeta _h}} \right)} \right] + 0.8\left[ {{{\left( { - {\zeta _h}} \right)}^{{{ - 1} \mathord{\left/

{\vphantom {{ - 1} 3}} \right.

\kern-\nulldelimiterspace} 3}}} - {{\left( { - \zeta } \right)}^{{{ - 1} \mathord{\left/

{\vphantom {{ - 1} 3}} \right.

\kern-\nulldelimiterspace} 3}}}} \right] + {\psi _h}\left( {\frac{{{z_{0h}}}}{L}} \right)} \right\}\]

Unstable :math:`\left( { - 0.465 \leqslant \zeta < 0} \right)` 



.. math::

 {\theta _{atm}} - {\theta _s} = \frac{{{\theta _*}}}{k}\left\{ {\left[ {\ln \frac{{{z_{atm,\,h}} - d}}{{{z_{0h}}}} - {\psi _h}\left( \zeta \right)} \right] + {\psi _h}\left( {\frac{{{z_{0h}}}}{L}} \right)} \right\}



Stable :math:`\left( {0 \leqslant \zeta \leqslant 1} \right)` 



.. math::

 {\theta _{atm}} - {\theta _s} = \frac{{{\theta _*}}}{k}\left\{ {\left[ {\ln \frac{{{z_{atm,\,h}} - d}}{{{z_{0h}}}} + 5\zeta } \right] - 5\frac{{{z_{0h}}}}{L}} \right\}



Very stable :math:`\left( {\zeta > 1} \right)` 



.. math::

 {\theta _{atm}} - {\theta _s} = \frac{{{\theta _*}}}{k}\left\{ {\left[ {\ln \frac{L}{{{z_{0h}}}} + 5} \right] + \left[ {5\ln \zeta + \zeta - 1} \right] - 5\frac{{{z_{0h}}}}{L}} \right\}

.

The specific humidity profiles are:

Very unstable :math:`\left( {\zeta < - 0.465} \right)` 

\[{q_{atm}} - {q_s} = \frac{{{q_*}}}{k}\left\{ {\left[ {\ln \frac{{{\zeta _w}L}}{{{z_{0w}}}} - {\psi _w}\left( {{\zeta _w}} \right)} \right] + 0.8\left[ {{{\left( { - {\zeta _w}} \right)}^{{{ - 1} \mathord{\left/

{\vphantom {{ - 1} 3}} \right.

\kern-\nulldelimiterspace} 3}}} - {{\left( { - \zeta } \right)}^{{{ - 1} \mathord{\left/

{\vphantom {{ - 1} 3}} \right.

\kern-\nulldelimiterspace} 3}}}} \right] + {\psi _w}\left( {\frac{{{z_{0w}}}}{L}} \right)} \right\}\]

Unstable :math:`\left( { - 0.465 \leqslant \zeta < 0} \right)` 



.. math::

 {q_{atm}} - {q_s} = \frac{{{q_*}}}{k}\left\{ {\left[ {\ln \frac{{{z_{atm,\,w}} - d}}{{{z_{0w}}}} - {\psi _w}\left( \zeta \right)} \right] + {\psi _w}\left( {\frac{{{z_{0w}}}}{L}} \right)} \right\}



Stable :math:`\left( {0 \leqslant \zeta \leqslant 1} \right)` 



.. math::

 {q_{atm}} - {q_s} = \frac{{{q_*}}}{k}\left\{ {\left[ {\ln \frac{{{z_{atm,\,w}} - d}}{{{z_{0w}}}} + 5\zeta } \right] - 5\frac{{{z_{0w}}}}{L}} \right\}



Very stable :math:`\left( {\zeta > 1} \right)` 



.. math::

 {q_{atm}} - {q_s} = \frac{{{q_*}}}{k}\left\{ {\left[ {\ln \frac{L}{{{z_{0w}}}} + 5} \right] + \left[ {5\ln \zeta + \zeta - 1} \right] - 5\frac{{{z_{0w}}}}{L}} \right\}



where



.. math::

 {\psi _h}\left( \zeta \right) = {\psi _w}\left( \zeta \right) = 2\ln \left( {\frac{{1 + {x^2}}}{2}} \right)

.

Using the definitions of :math:`{u_ * }` , :math:`{\theta _ * }` , and :math:`{q_ * }` , an iterative solution of these equations can be used to calculate the surface momentum, sensible heat, and water vapor flux using atmospheric and surface values for :math:`\left| u \right|` , :math:`\theta` , and :math:`q` except that :math:`L` depends on :math:`{u_ * }` , :math:`{\theta _ * }` , and :math:`{q_ * }` . However, the bulk number



.. math::

 {R_{iB}} = \frac{{{\theta _{v,\,atm}} - {\theta _{v,\,s}}}}{{\overline {{\theta _{v,\,atm}}} }}\frac{{g\left( {{z_{atm,\,m}} - d} \right)}}{{V_a^2}}



is related to :math:`\zeta` (Arya 2001) as



.. math::

 {R_{iB}} = \zeta \left[ {\ln \left( {\frac{{{z_{atm,\,h}} - d}}{{{z_{0h}}}}} \right) - {\psi _h}\left( \zeta \right)} \right]{\left[ {\ln \left( {\frac{{{z_{atm,\,m}} - d}}{{{z_{0m}}}}} \right) - {\psi _m}\left( \zeta \right)} \right]^{ - 2}}

.

Using ${\phi _h} = \phi _m^2 = {\left( {1 - 16\zeta } \right)^{ - {1 \mathord{\left/

{\vphantom {1 2}} \right.

\kern-\nulldelimiterspace} 2}}} :math:`for unstable conditions and` {\phi _h} = {\phi _m} = 1 + 5\zeta :math:`for stable conditions to determine` {\psi _m}\left( \zeta \right) :math:`and` {\psi _h}\left( \zeta \right) :math:`, the inverse relationship` \zeta = f\left( {{R_{iB}}} \right) :math:`can be solved to obtain a first guess for` \zeta :math:`and thus` L$ from

\[\begin{gathered}

\zeta = \frac{{{R_{iB}}\ln \left( {\frac{{{z_{atm,\,m}} - d}}{{{z_{0m}}}}} \right)}}{{1 - 5\min \left( {{R_{iB}},0.19} \right)}} & 0.01 \leqslant \zeta \leqslant 2 & {\text{for }}{R_{iB}} \geqslant 0{\text{ (neutral or stable)}} \hfill \\

\zeta = {R_{iB}}\ln \left( {\frac{{{z_{atm,\,m}} - d}}{{{z_{0m}}}}} \right) & - 100 \leqslant \zeta \leqslant - 0.01 & {\text{for }}{R_{iB}} < 0{\text{ (unstable)}} \hfill \\

\end{gathered} \].

Upon iteration, the following is used to determine :math:`\zeta` and thus :math:`L` 



.. math::

 \zeta = \frac{{\left( {{z_{atm,\,m}} - d} \right)kg{\theta _{v * }}}}{{u_ * ^2\overline {{\theta _{v,\,atm}}} }}



where

\[\begin{gathered}
#. 01 \leqslant \zeta \leqslant 2 & {\text{for }}\zeta \geqslant 0{\text{ (neutral or stable)}} \hfill \\

{\text{ - 100}} \leqslant \zeta \leqslant {\text{ - 0}}{\text{.01}} & {\text{for }}\zeta < 0{\text{ (unstable)}} \hfill \\

\end{gathered} \].

The momentum, sensible heat, and water vapor fluxes between the surface and the atmosphere can also be written in the form



.. math::

 {\tau _x} = - {\rho _{atm}}\frac{{\left( {{u_{atm}} - {u_s}} \right)}}{{{r_{am}}}}





.. math::

 {\tau _y} = - {\rho _{atm}}\frac{{\left( {{v_{atm}} - {v_s}} \right)}}{{{r_{am}}}}





.. math::

 H = - {\rho _{atm}}{C_p}\frac{{\left( {{\theta _{atm}} - {\theta _s}} \right)}}{{{r_{ah}}}}





.. math::

 E = - {\rho _{atm}}\frac{{\left( {{q_{atm}} - {q_s}} \right)}}{{{r_{aw}}}}



where :math:`{r_{am}}` , :math:`{r_{ah}}` , and :math:`{r_{aw}}` are the aerodynamic resistances for momentum, sensible heat and latent heat, respectively (s m-1).

Sensible and Latent Heat and Momentum Fluxes
------------------------------------------------

The solution for the heat and momentum fluxes is presented in roughly the order in which the equations are solved in the Fortran code.

Roughness Length and Displacement Height
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The roughness length and displacement height for the urban canopy are needed. Grimmond and Oke (1999) review approaches to calculate these parameters from morphometric methods. Here, we use the Macdonald et al. (1998) approach, which appears to be a reasonable compromise between minimizing input requirements and yielding acceptable results. The subscript “canopy” is used to distinguish between an aerodynamic parameter for the urban canopy versus a parameter for an individual urban surface (e.g., roof).

The canopy displacement height :math:`{d_{canopy}}` (m) is

 :math:`{d_{canopy}} = H\left[ {1 + {\alpha ^{ - {\lambda _P}}}({\lambda _P} - 1)} \right]` 

where :math:`H` is the canyon (roof) height (m) (Table 1.3), :math:`\alpha = {\text{4}}{\text{.43}}` is an empirical coefficient, and :math:`{\lambda _P}` is the plan area index. The plan area index :math:`{\lambda _p}` is

${\lambda _p} = \frac{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W} + 1}}$

where ${H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}$ is the height to width ratio of the urban canyon (Table 1.3).

The canopy roughness length 

.. math::

 {z_{0m,\,canopy}}

 (m) for momentum is

 :math:`{z_{0m,\,canopy}} = H\left( {1 - \frac{{{d_{canopy}}}}{H}} \right)\exp \left\{ { - {{\left[ {0.5B\frac{{{C_D}}}{{{k^2}}}\left( {1 - \frac{{{d_{canopy}}}}{H}} \right){\lambda _F}} \right]}^{ - 0.5}}} \right\}` 

where :math:`B = 1` is a correction to the drag coefficient to account for variable obstacle shapes and flow conditions, :math:`{C_D} = 1.2` is the depth-integrated mean drag coefficient for surface-mounted cubes in a shear flow, :math:`k` is the von Karman constant, and :math:`{\lambda _F}` is the frontal area index. The frontal area index :math:`{\lambda _F}` is

\[{\lambda _F} = \left( {1 - {\lambda _P}} \right)\left( {{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}} \right)\sqrt {\frac{{{B_L}{\lambda _P}}}{{{B_S}}}} \]

where ${{{B_S}} \mathord{\left/

{\vphantom {{{B_S}} {{B_L}}}} \right.

\kern-\nulldelimiterspace} {{B_L}}} :math:`is the building shortside to longside ratio (here set equal to` {\lambda _P}$).

Several checks are made to ensure that the derived aerodynamic parameters are consistent with the canyon structure and atmospheric forcing. First, the canyon height :math:`H` minus the canopy displacement height :math:`{d_{canopy}}` must be greater than the canopy roughness length :math:`{z_{0m,\,canopy}}` . Second, the atmospheric wind forcing height :math:`{z_{atm,\,m}}` (Table 1.1) minus the canopy displacement height :math:`{d_{canopy}}` must be greater than the canopy roughness length :math:`{z_{0m,\,canopy}}` . Note that :math:`{z_{0m,\,canopy}} = {z_{0h,\,canopy}} = {z_{0w,\,canopy}}` and :math:`{z_{atm}} = {z'_{atm}} + {z_{0,\,canopy}} + {z_{d,\,canopy}}` (Table 1.1) where :math:`{z'_{atm}}` is the reference height from the atmospheric model.

Wind Speed in the Urban Canyon
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Following Masson (2000) and Lemonsu et al. (2004), the wind speed in the canyon is the combination of the mean horizontal canyon wind :math:`{U_{can}}` (m s-1) and the turbulent (vertical) wind :math:`{W_{can}}` (m s-1)

 :math:`{U_{ac}} = \sqrt {{U_{can}}^2 + {W_{can}}^2}` .

To calculate the horizontal wind speed in the canyon :math:`{U_{can}}` (m s-1), a horizontal wind speed at the top of the canyon is derived by assuming a logarithmic wind profile from the atmospheric reference height to the canyon top. The wind is then extrapolated to a height inside the canyon using an exponential profile. For skimming flow (${H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W} \geqslant 1 :math:`) (Oke 1987), a zero` {U_{can}}$ is assumed when the mean flow is perpendicular to the canyon orientation. After integration over 360º (to account for all street orientations),

\[{U_{can}} = {V_r}\frac{2}{\pi }\frac{{\ln \left( {\frac{{H - {d_{canopy}}}}{{{z_{0m,\,canopy}}}}} \right)}}{{\ln \left( {\frac{{{z_{atm,\,m}} - {d_{canopy}}}}{{{z_{0m,\,canopy}}}}} \right)}}\exp \left[ { - 0.5\left( {{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}} \right)\left( {1 - \frac{{{H_w}}}{H}} \right)} \right]\]

where :math:`{H_w}` is the height at which the wind speed is estimated (Table 1.3). For isolated roughness flow (${H \mathord{\left/

{\vphantom {H {W < {\text{0}}{\text{.5}}}}} \right.

\kern-\nulldelimiterspace} {W < {\text{0}}{\text{.5}}}}$), the wind speed in the canyon is assumed to be independent of the orientation of the mean atmospheric flow above the canyon level,

\[{U_{can}} = {V_r}\frac{{\ln \left( {\frac{{H - {d_{canopy}}}}{{{z_{0m,\,canopy}}}}} \right)}}{{\ln \left( {\frac{{{z_{atm,\,m}} - {d_{canopy}}}}{{{z_{0m,\,canopy}}}}} \right)}}\exp \left( { - 0.5\left( {{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}} \right)\left( {1 - \frac{{{H_w}}}{H}} \right)} \right)\].

For wake interference flow (${{0.5 \leqslant H} \mathord{\left/

{\vphantom {{0.5 \leqslant H} {W < {\text{1}}{\text{.0}}}}} \right.

\kern-\nulldelimiterspace} {W < {\text{1}}{\text{.0}}}}$),

\[\begin{gathered}

{U_{can}} = {V_r}\left[ {1 + 2\left( {\frac{2}{\pi } - 1} \right)\left( {\frac{H}{W} - \frac{1}{2}} \right)} \right]\frac{{\ln \left( {\frac{{H - {d_{canopy}}}}{{{z_{0m,\,canopy}}}}} \right)}}{{\ln \left( {\frac{{{z_{atm,\,m}} - {d_{canopy}}}}{{{z_{0m,\,canopy}}}}} \right)}} \hfill \\

& \times \exp \left( { - 0.5\left( {{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}} \right)\left( {1 - \frac{{{H_w}}}{H}} \right)} \right) \hfill \\

\end{gathered} \].

The magnitude of the reference level atmospheric wind is

 :math:`{V_r} = \sqrt {u_{atm}^2 + v_{atm}^2} \geqslant 1` 

where zonal and meridional winds :math:`{u_{atm}}` and :math:`{v_{atm}}` (m s-1) are at height :math:`{z_{atm,\,m}}` . The turbulent (vertical) wind :math:`{W_{can}}` (m s-1) is assumed to be equal to the friction velocity (Masson 2000), which is determined from the solution for turbulent fluxes (section 3.2.3).

Iterative Solution for Urban Canopy Air Temperature and Humidity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because of the interdependence between fluxes, aerodynamic resistances, and canyon air temperature and humidity, an iterative solution for the UCL air is devised.

An initial guess for the wind speed :math:`{V_a}` (equation is obtained assuming an initial convective velocity :math:`{U_c} = 0` m s-1 for stable conditions and :math:`{U_c} = 0.5` for unstable conditions. Stable conditions ( :math:`{\theta _{v,\,atm}} - {\theta _{v,\,s}} \geqslant 0` ) and unstable conditions ( :math:`{\theta _{v,\,atm}} - {\theta _{v,\,s}} < 0` ) are evaluated from the difference in virtual potential air temperature between the reference height and the surface where

 :math:`{\theta _{v,\,atm}} - {\theta _{v,\,s}} = \left( {{\theta _{atm}} - {\theta _s}} \right)\left( {1 + 0.61{q_{atm}}} \right) + 0.61\overline {{\theta _{atm}}} \left( {{q_{atm}} - {q_s}} \right)` .

Here, :math:`{\theta _s} = {T_{ac}}` and :math:`{q_s} = {q_{ac}}` where :math:`{T_{ac}}` is the air temperature in the UCL (K) and :math:`{q_{ac}}` is the specific humidity in the UCL (kg kg-1) (Figure 3.1). The air temperature and specific humidity from the previous time step are used. The temperature :math:`{\theta _{atm}}` is defined by equation , :math:`\overline {{\theta _{atm}}}` is the atmospheric potential temperature (Table 1.1), and :math:`{q_{atm}}` is the atmospheric specific humidity (kg kg-1) (Table 1.1). An initial guess for the Monin-Obukhov length :math:`L` is obtained from the bulk Richardson number using equations and .

The iterative solution begins with the friction velocity :math:`{u_ * }` , potential temperature scale :math:`{\theta _ * }` , and humidity scale :math:`{q_ * }` being calculated from equations -. Now that the friction velocity has been determined, the wind in the urban canopy, :math:`{U_{ac}}` , is calculated from equation . The aerodynamic resistances (s m-1) to momentum, sensible heat, and latent heat transfer between the UCL air and the atmosphere are



.. math::

 {r_{am}} = \frac{{{V_a}}}{{u_ * ^2}} = \frac{1}{{{k^2}{V_a}}}{\left[ {\ln \left( {\frac{{{z_{atm,\,m}} - {d_{canopy}}}}{{{z_{0m,\,canopy}}}}} \right) - {\psi _m}\left( {\frac{{{z_{atm,\,m}} - {d_{canopy}}}}{L}} \right) + {\psi _m}\left( {\frac{{{z_{0m,\,canopy}}}}{L}} \right)} \right]^2}

,

\[\begin{gathered}

{r_{ah}} = \frac{{{\theta_{atm}} - {\theta_s}}}{{{\theta_* }{u_*}}} = \frac{1}{{{k^2}{V_a}}}\left[ {\ln \left( {\frac{{{z_{atm,\,m}} - {d_{canopy}}}}{{{z_{0m,\,canopy}}}}} \right) - {\psi _m}\left( {\frac{{{z_{atm,\,m}} - {d_{canopy}}}}{L}} \right) + {\psi_m}\left( {\frac{{{z_{0m,\,canopy}}}}{L}} \right)} \right] \hfill \\

& \left[ {\ln \left( {\frac{{{z_{atm,\,h}} - {d_{canopy}}}}{{{z_{0h,\;canopy}}}}} \right) - {\psi_h}\left( {\frac{{{z_{atm,\,h}} - {d_{canopy}}}}{L}} \right) + {\psi_h}\left( {\frac{{{z_{0h,\,canopy}}}}{L}} \right)} \right] \hfill \\

\end{gathered} \],

\[\begin{gathered}

{r_{aw}} = \frac{{{q_{atm}} - {q_s}}}{{{q_*}{u_*}}} = \frac{1}{{{k^2}{V_a}}}\left[ {\ln \left( {\frac{{{z_{atm,\,m}} - {d_{canopy}}}}{{{z_{0m,\,canopy}}}}} \right) - {\psi_m}\left( {\frac{{{z_{atm,\,m}} - {d_{canopy}}}}{L}} \right) + {\psi_m}\left( {\frac{{{z_{0m,\,canopy}}}}{L}} \right)} \right] \hfill \\

& \left[ {\ln \left( {\frac{{{z_{atm,\,w}} - {d_{canopy}}}}{{{z_{0w,\,canopy}}}}} \right) - {\psi_w}\left( {\frac{{{z_{atm,\,w}} - {d_{canopy}}}}{L}} \right) + {\psi_w}\left( {\frac{{{z_{0w,\,canopy}}}}{L}} \right)} \right] \hfill \\

\end{gathered} \].

The resistances to sensible heat and latent heat transfer between canyon surfaces (roof, sunlit and shaded wall, pervious and impervious road) and the UCL depend only on canyon wind speed following Masson (2000). Thus, the surface resistances, :math:`{r_{s,\,roof}}` , :math:`{r_{s,\,sunwall}}` , :math:`{r_{s,\,shdwall}}` , :math:`{r_{s,\,prvrd}}` , :math:`{r_{s,\,imprvrd}}` , (s m-1) are identical and are determined from (Rowley et al. 1930)

 :math:`{r_{s,\,u}} = \frac{{{\rho _{atm}}{C_p}}}{{11.8 + 4.2{U_{ac}}}}` .

The UCL air temperature and specific humidity are determined by solving the following systems of equations. For :math:`{T_{ac}}` 



.. math::

 {H_{roof}} = - {\rho _{atm}}{C_p}\frac{{{T_{ac}} - {T_{g,\,roof}}}}{{{r_{s,\,roof}}}}

,



.. math::

 {H_{prvrd}} = - {\rho _{atm}}{C_p}\frac{{{T_{ac}} - {T_{g,\,prvrd}}}}{{{r_{s,\,prvrd}}}}

,



.. math::

 {H_{imprvrd}} = - {\rho _{atm}}{C_p}\frac{{{T_{ac}} - {T_{g,\,imprvrd}}}}{{{r_{s,\,imprvrd}}}}

,



.. math::

 {H_{sunwall}} = - {\rho _{atm}}{C_p}\frac{{{T_{ac}} - {T_{g,\,sunwall}}}}{{{r_{s,\,sunwall}}}}

,



.. math::

 {H_{shdwall}} = - {\rho _{atm}}{C_p}\frac{{{T_{ac}} - {T_{g,\,shdwall}}}}{{{r_{s,\,shdwall}}}}

,

\[\begin{gathered}

H = - {\rho _{atm}}{C_p}\frac{{{\theta _{atm}} - {T_{ac}}}}{{{r_{ah}}}} \hfill \\

& = {W_{roof}}{H_{roof}} + \left( {1 - {W_{roof}}} \right) \times \hfill \\

& & \left[ {{f_{prvrd}}{H_{prvrd}} + \left( {1 - {f_{prvrd}}} \right){H_{imprvrd}} + \frac{H}{W}{H_{sunwall}} + \frac{H}{W}{H_{shdwall}}} \right] \hfill \\

\end{gathered} \]

where :math:`H` is sensible heat flux (W m-2) and :math:`{T_g}` is the surface temperature of each urban surface. The term :math:`{W_{roof}}` is the relative contribution of roof fluxes to the total urban landunit flux (Table 1.3). The term :math:`1 - {W_{roof}}` is then the relative contribution of the canyon to the total urban landunit flux. The term :math:`{f_{prvrd}}` is the fraction of road that is pervious (Table 1.3) and the term :math:`1 - {f_{prvrd}}` is the fraction of the road that is impervious. Note that the factor ${H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}$ for the sunwall and shadewall converts the flux from watts per meter squared of surface area to watts per meter squared of ground area.

In Oleson et al. (2008a), an additional heat flux :math:`{H_{wasteheat}}` , the sensible heat flux from waste heat generated by space heating and air conditioning, was included in equation . However, if this flux is large enough, the numerical solution may become unstable because of the canopy air has no heat capacity and the heat capacities of the roofs and walls are relatively small. Instead, this heat flux is added to the net heat flux for the canyon floor (section 4.1).

Equations - can be solved for the UCL air temperature as

\[{T_{ac}} = \frac{{\left( \begin{gathered}

c_a^h{\theta _{atm}} + {c_{roof}}{T_{g,\,roof}} + {c_{prvrd}}{T_{g,\,prvrd}} + \hfill \\

{c_{imprvrd}}{T_{g,\,imprvrd}} + {c_{sunwall}}{T_{g,\,sunwall}} + {c_{shdwall}}{T_{g,\,shdwall}} \hfill \\

\end{gathered} \right)}}{{c_a^h + {c_{roof}} + {c_{prvrd}} + {c_{imprvrd}} + {c_{sunwall}} + {c_{shdwall}}}}\]

where :math:`c_a^h` is the sensible heat conductance from the UCL to the atmosphere (${1 \mathord{\left/

{\vphantom {1 {{r_{ah}}}}} \right.

\kern-\nulldelimiterspace} {{r_{ah}}}} :math:`), and` {c_{roof}} :math:`,` {c_{prvrd}} :math:`,` {c_{imprvrd}} :math:`,` {c_{sunwall}} :math:`, and` {c_{shdwall}} :math:`are the weighted heat conductances from urban surfaces to UCL air [` {{{W_{roof}}} \mathord{\left/

{\vphantom {{{W_{roof}}} {{r_{s,\,roof}}}}} \right.

\kern-\nulldelimiterspace} {{r_{s,\,roof}}}} :math:`,` {{{W_{prvrd}}} \mathord{\left/

{\vphantom {{{W_{prvrd}}} {{r_{s,\,prvrd}}}}} \right.

\kern-\nulldelimiterspace} {{r_{s,\,prvrd}}}} :math:`,` {{{W_{imprvrd}}} \mathord{\left/

{\vphantom {{{W_{imprvrd}}} {{r_{s,\,imprvrd}}}}} \right.

\kern-\nulldelimiterspace} {{r_{s,\,imprvrd}}}} :math:`,` {{{W_{sunwall}}} \mathord{\left/

{\vphantom {{{W_{sunwall}}} {{r_{s,\,sunwall}}}}} \right.

\kern-\nulldelimiterspace} {{r_{s,\,sunwall}}}}$, \[{{{W_{shawall}}} \mathord{\left/

{\vphantom {{{W_{shawall}}} {{r_{s,\,shawall}}}}} \right.

\kern-\nulldelimiterspace} {{r_{s,\,shawall}}}}\], respectively, where :math:`{W_{prvrd}} = \left( {1 - {W_{roof}}} \right){f_{prvrd}}` , :math:`{W_{imprvrd}} = \left( {1 - {W_{roof}}} \right)\left( {1 - {f_{prvrd}}} \right)` , ${W_{sunwall}} = \left( {1 - {W_{roof}}} \right)\left( {{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}} \right) :math:`, and` {W_{shdwall}} = \left( {1 - {W_{roof}}} \right)\left( {{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}} \right)$].

Similarly, the system of equations for the UCL air specific humidity, :math:`{q_{ac}}` , is



.. math::

 {E_{roof}} = - {\rho _{atm}}\frac{{{f_{wet,\,roof}}\left( {{q_{ac}} - {q_{g,\,roof}}} \right)}}{{{r_{s,\,roof}}}}

,



.. math::

 {E_{prvrd}} = - {\rho _{atm}}\frac{{{q_{ac}} - {q_{g,\,prvrd}}}}{{{r_{s,\,prvrd}}}}

,



.. math::

 {E_{imprvrd}} = - {\rho _{atm}}\frac{{{f_{wet,\,imprvrd}}\left( {{q_{ac}} - {q_{g,\,imprvrd}}} \right)}}{{{r_{s,\,imprvrd}}}}

,



.. math::

 {E_{sunwall}} = 0

,



.. math::

 {E_{shdwall}} = 0

,

\[\begin{gathered}

E = - {\rho _{atm}}\frac{{{q_{atm}} - {q_{ac}}}}{{{r_{aw}}}} \hfill \\

& = {W_{roof}}{E_{roof}} + \left( {1 - {W_{roof}}} \right)\left[ {{f_{prvrd}}{E_{prvrd}} + \left( {1 - {f_{prvrd}}} \right){E_{imprvrd}}} \right] \hfill \\

\end{gathered} \]

where :math:`E` is water vapor flux (kg m-2 s-1) and :math:`{q_g}` is the specific humidity at each urban surface (kg kg-1). Note that the latent heat flux from the sunlit and shaded walls is zero. The term :math:`{f_{wet}}` is the fraction of the roof or impervious road surface that is wet. If there is dew formation ( :math:`{q_{ac}} - {q_g} > 0` ), then :math:`{f_{wet}} = 1` . If there is snow on the surface ( :math:`{z_{sno}} > 0` ), :math:`{f_{wet}}` is determined from the snow depth :math:`{z_{sno}}` as



.. math::

 {f_{wet}} = \frac{{{z_{sno}}}}{{0.05}} \leqslant 1

.

In the absence of snow,

\[{f_{wet}} = {\left( {\frac{{{w_{liq,\,snl + 1}} + {w_{ice,\,snl + 1}}}}{{{w_{pond,\max }}}}} \right)^{{2 \mathord{\left/

{\vphantom {2 3}} \right.

\kern-\nulldelimiterspace} 3}}} \leqslant 1\]

where :math:`{w_{liq,\,snl + 1}}` and :math:`{w_{ice,\,snl + 1}}` are the mass of ice and liquid water (kg m-2) stored on top of the urban surface and :math:`{w_{pond,\max }}` is the maximum amount of water that the surface can hold (Chapter 5). This latter formulation is analogous to the treatment of the wetted fraction of the vegetated canopy in CLM (Oleson et al. 2004).

In equations and , the specific humidity of the roof and the impervious road surfaces, :math:`{q_{g,\,roof}}` and :math:`{q_{g,\,imprvrd}}` , is set to the saturated specific humidity evaluated at their respective surface temperatures, :math:`q_{sat}^{{T_{g,\,roof}}}` and :math:`q_{sat}^{{T_{g,\,imprvrd}}}` (section 3.3).

As noted in section 1.1.3, a simplified bulk parameterization approach is used to represent evaporation from the pervious surface. The pervious road specific humidity, :math:`{q_{g,\,prvrd}}` , is evaluated as a function of the wetness of the soil column. This allows all of the soil moisture to potentially be available for evaporation. The specific humidity is



.. math::

 {q_{g,\,prvrd}} = \alpha q_{sat}^{{T_g}}

,

where :math:`q_{sat}^{{T_g}}` is the saturated specific humidity at the surface temperature :math:`{T_g}` (section 4.1). The factor :math:`\alpha` is a weighted combination of values for the soil column and snow



.. math::

 \alpha = {\alpha _{soi}}\left( {1 - {f_{sno}}} \right) + {\alpha _{sno}}{f_{sno}}



where :math:`{f_{sno}}` is the fraction of ground covered by snow (equation ), and :math:`{\alpha _{sno}} = 1.0` . The term :math:`{\alpha _{soi}}` is a function ranging from one when the soil column is wet to zero when the soil is dry

 :math:`{\alpha _{soi}} = \sum\limits_{i = 1}^{{N_{levsoi}}} {{w_i}{r_i}}` 

where :math:`{w_i}` is a soil wetness factor for layer :math:`i` , and :math:`{r_i}` is the relative contribution of each layer. The wetness factor :math:`{w_i}` is

\[{w_i} = \left\{ \begin{gathered}

\frac{{{\theta _{liq,\,i}} - {\theta _{dry,\,i}}}}{{{\theta _{opt,\,i}} - {\theta _{dry,\,i}}}} & {\text{for }}{T_i} \geqslant {T_f} \hfill \\

0 & {\text{for }}{T_i} < {T_f} \hfill \\

\end{gathered} \right\}\]

where :math:`{\theta _{liq,\,i}} - {\theta _{dry,\,i}} \geqslant 0` and



.. math::

 {r_i} = 0.1 & {\text{for }}i = 1, \ldots ,{N_{levsoi}}

.

The term :math:`{\theta _{dry}}` is the volumetric water content at which evapotranspiration ceases and :math:`{\theta _{opt}}` is the optimal water content



.. math::

 {\theta _{dry,\,i}} = {\theta _{sat,\;i}}{\left( {\frac{{ - 316230}}{{{\psi _{sat,\,i}}}}} \right)^{ - \frac{1}{{{B_i}}}}}





.. math::

 {\theta _{opt,\,i}} = {\theta _{sat,\,i}}{\left( {\frac{{ - 158490}}{{{\psi _{sat,\,i}}}}} \right)^{ - \frac{1}{{{B_i}}}}}



where :math:`{\theta _{sat,\,i}}` is the water content at saturation (i.e., porosity), :math:`{\psi _{sat,\,i}}` is the saturated soil matric potential (mm), and :math:`{B_i}` is the Clapp-Hornberger exponent (section 5.3.1). The soil volumetric liquid water content :math:`{\theta _{liq,\,i}}` is



.. math::

 {\theta _{liq,\,i}} = \frac{{{w_{liq,\,i}}}}{{\Delta {z_i}{\rho _{liq}}}} \leqslant {\theta _{sat,\,i}} - {\theta _{ice,\,i}}



where :math:`{w_{liq,\,i}}` is the mass of liquid water (kg m-2), :math:`\Delta {z_i}` is the layer thickness, :math:`{\rho _{liq}}` is the density of liquid water (kg m-3) (Table 1.4), and :math:`{\theta _{ice,\,i}}` is the volumetric ice content



.. math::

 {\theta _{ice,\,i}} = \frac{{{w_{ice,\,i}}}}{{\Delta {z_i}{\rho _{ice}}}} \leqslant {\theta _{sat,\,i}}



where :math:`{w_{ice,\,i}}` is the mass of ice (kg m-2) and :math:`{\rho _{ice}}` is the density of ice (kg m-3) (Table 1.4). If :math:`q_{sat}^{{T_g}} > {q_{atm}}` and :math:`{q_{atm}} > {q_{g,\,prvrd}}` , then :math:`{q_{g,\,prvrd}} = {q_{atm}}` and :math:`\frac{{d{q_{g,\,prvrd}}}}{{d{T_g}}} = 0` .

The UCL specific humidity is then



.. math::

 {q_{ac}} = \frac{{c_a^w{q_{atm}} + {c_{roof}}{f_{wet,\,roof}}{q_{g,\,roof}} + {c_{prvrd}}{q_{g,\,prvrd}} + {c_{imprvrd}}{f_{wet,\,imprvrd}}{q_{g,\,imprvrd}}}}{{c_a^w + {f_{wet,\,roof}}{c_{roof}} + {c_{prvrd}} + {f_{wet,\,imprvrd}}{c_{imprvrd}}}}



where :math:`c_a^w` is the latent heat conductance from the UCL air to the atmosphere (${1 \mathord{\left/

{\vphantom {1 {{r_{aw}}}}} \right.

\kern-\nulldelimiterspace} {{r_{aw}}}} :math:`), and` {c_{roof}} :math:`,` {c_{prvrd}} :math:`, and` {c_{imprvrd}} :math:`are the weighted heat conductances from urban surfaces to UCL air [` {{{W_{roof}}} \mathord{\left/

{\vphantom {{{W_{roof}}} {{r_{s,\,roof}}}}} \right.

\kern-\nulldelimiterspace} {{r_{s,\,roof}}}} :math:`,` {{{W_{prvrd}}} \mathord{\left/

{\vphantom {{{W_{prvrd}}} {{r_{s,\,prvrd}}}}} \right.

\kern-\nulldelimiterspace} {{r_{s,\,prvrd}}}} :math:`,` {{{W_{imprvrd}}} \mathord{\left/

{\vphantom {{{W_{imprvrd}}} {{r_{s,\,imprvrd}}}}} \right.

\kern-\nulldelimiterspace} {{r_{s,\,imprvrd}}}} :math:`, respectively, where` {W_{prvrd}} = \left( {1 - {W_{roof}}} \right){f_{prvrd}} :math:`,` {W_{imprvrd}} = \left( {1 - {W_{roof}}} \right)\left( {1 - {f_{prvrd}}} \right)$].

The stability is then updated using the new UCL air temperature and specific humidity as follows. The potential temperature, specific humidity, and virtual potential temperature scales, :math:`{\theta _ * }` , :math:`{q_ * }` , and 

.. math::

 {\theta _{v * }}

, are reevaluated using equations - and . The wind speed including the convective velocity is reevaluated using equations and -. The Monin-Obukhov length is updated from equation . This sequence of calculations is repeated for a total of three times beginning with the calculation of the friction velocity :math:`{u_ * }` (equations -).

Final Fluxes and Adjustments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The sensible and latent heat fluxes and momentum flux from urban surfaces are then calculated from equations -, -, and - using the updated UCL air temperature and specific humidity. The water vapor flux from the pervious road, :math:`{E_{prvrd}}` , is assigned to ground evaporation, :math:`{E_{g,\,prvrd}}` , or a evapotranspiration term, :math:`E_{prvrd}^{et}` , as follows

\[\begin{gathered}

{E_{g,\,prvrd}} = {E_{prvrd}} & {\text{for }}{q_s} - {q_{g,\;prvrd}} > 0{\text{ or }}{f_{sno}} > 0{\text{ or }}{\alpha _{soi}} = 0 \hfill \\

E_{prvrd}^{et} = {E_{prvrd}} & {\text{otherwise}} \hfill \\

\end{gathered} \]

This ensures that dew can form on snow or soil surfaces and that snow can sublimate. Otherwise, the evaporation is assigned to an evapotranspiration term in which the water for evaporation is removed from all soil layers which have sufficient liquid water (section 5.3).

The partial derivatives of the urban surface fluxes with respect to surface temperatures, which are needed for the soil temperature calculation and to update the urban surface fluxes, are



.. math::

 \frac{{\partial {H_{roof}}}}{{\partial {T_{g,\,roof}}}} = \frac{{{\rho _{atm}}{C_p}\left( {c_a^h + \frac{{{c_{prvrd}}}}{{{W_{prvrd}}}} + \frac{{{c_{imprvrd}}}}{{{W_{imprvrd}}}} + \frac{{{c_{sunwall}}}}{{{W_{sunwall}}}} + \frac{{{c_{shdwall}}}}{{{W_{shdwall}}}}} \right)\frac{{{c_{roof}}}}{{{W_{roof}}}}}}{{c_a^h + \frac{{{c_{roof}}}}{{{W_{roof}}}} + \frac{{{c_{prvrd}}}}{{{W_{prvrd}}}} + \frac{{{c_{imprvrd}}}}{{{W_{imprvrd}}}} + \frac{{{c_{sunwall}}}}{{{W_{sunwall}}}} + \frac{{{c_{shdwall}}}}{{{W_{shdwall}}}}}}





.. math::

 \frac{{\partial {H_{prvrd}}}}{{\partial {T_{g,\,prvrd}}}} = \frac{{{\rho _{atm}}{C_p}\left( {c_a^h + \frac{{{c_{roof}}}}{{{W_{roof}}}} + \frac{{{c_{imprvrd}}}}{{{W_{imprvrd}}}} + \frac{{{c_{sunwall}}}}{{{W_{sunwall}}}} + \frac{{{c_{shdwall}}}}{{{W_{shdwall}}}}} \right)\frac{{{c_{prvrd}}}}{{{W_{prvrd}}}}}}{{c_a^h + \frac{{{c_{roof}}}}{{{W_{roof}}}} + \frac{{{c_{prvrd}}}}{{{W_{prvrd}}}} + \frac{{{c_{imprvrd}}}}{{{W_{imprvrd}}}} + \frac{{{c_{sunwall}}}}{{{W_{sunwall}}}} + \frac{{{c_{shdwall}}}}{{{W_{shdwall}}}}}}





.. math::

 \frac{{\partial {H_{imprvrd}}}}{{\partial {T_{g,\,imprvrd}}}} = \frac{{{\rho _{atm}}{C_p}\left( {c_a^h + \frac{{{c_{roof}}}}{{{W_{roof}}}} + \frac{{{c_{prvrd}}}}{{{W_{prvrd}}}} + \frac{{{c_{sunwall}}}}{{{W_{sunwall}}}} + \frac{{{c_{shdwall}}}}{{{W_{shdwall}}}}} \right)\frac{{{c_{imprvrd}}}}{{{W_{imprvrd}}}}}}{{c_a^h + \frac{{{c_{roof}}}}{{{W_{roof}}}} + \frac{{{c_{prvrd}}}}{{{W_{prvrd}}}} + \frac{{{c_{imprvrd}}}}{{{W_{imprvrd}}}} + \frac{{{c_{sunwall}}}}{{{W_{sunwall}}}} + \frac{{{c_{shdwall}}}}{{{W_{shdwall}}}}}}





.. math::

 \frac{{\partial {H_{sunwall}}}}{{\partial {T_{g,\,sunwall}}}} = \frac{{{\rho _{atm}}{C_p}\left( {c_a^h + \frac{{{c_{roof}}}}{{{W_{roof}}}} + \frac{{{c_{prvrd}}}}{{{W_{prvrd}}}} + \frac{{{c_{imprvrd}}}}{{{W_{imprvrd}}}} + \frac{{{c_{shdwall}}}}{{{W_{shdwall}}}}} \right)\frac{{{c_{sunwall}}}}{{{W_{sunwall}}}}}}{{c_a^h + \frac{{{c_{roof}}}}{{{W_{roof}}}} + \frac{{{c_{prvrd}}}}{{{W_{prvrd}}}} + \frac{{{c_{imprvrd}}}}{{{W_{imprvrd}}}} + \frac{{{c_{sunwall}}}}{{{W_{sunwall}}}} + \frac{{{c_{shdwall}}}}{{{W_{shdwall}}}}}}





.. math::

 \frac{{\partial {H_{shdwall}}}}{{\partial {T_{g,\,shdwall}}}} = \frac{{{\rho _{atm}}{C_p}\left( {c_a^h + \frac{{{c_{roof}}}}{{{W_{roof}}}} + \frac{{{c_{prvrd}}}}{{{W_{prvrd}}}} + \frac{{{c_{imprvrd}}}}{{{W_{imprvrd}}}} + \frac{{{c_{sunwall}}}}{{{W_{sunwall}}}}} \right)\frac{{{c_{shdwall}}}}{{{W_{shdwall}}}}}}{{c_a^h + \frac{{{c_{roof}}}}{{{W_{roof}}}} + \frac{{{c_{prvrd}}}}{{{W_{prvrd}}}} + \frac{{{c_{imprvrd}}}}{{{W_{imprvrd}}}} + \frac{{{c_{sunwall}}}}{{{W_{sunwall}}}} + \frac{{{c_{shdwall}}}}{{{W_{shdwall}}}}}}

,



.. math::

 \frac{{\partial {E_{roof}}}}{{\partial {T_{g,\,roof}}}} = \frac{{{\rho _{atm}}\left( {c_a^w + \frac{{{c_{prvrd}}}}{{{W_{prvrd}}}} + \frac{{{f_{wet,\,imprvrd}}{c_{imprvrd}}}}{{{W_{imprvrd}}}}} \right)\frac{{{f_{wet,\,roof}}{c_{roof}}}}{{{W_{roof}}}}}}{{c_a^w + \frac{{{f_{wet,\,roof}}{c_{roof}}}}{{{W_{roof}}}} + \frac{{{c_{prvrd}}}}{{{W_{prvrd}}}} + \frac{{{f_{wet,\,imprvrd}}{c_{imprvrd}}}}{{{W_{imprvrd}}}}}}\frac{{d{q_{g,\,roof}}}}{{d{T_{g,\,roof}}}}

,



.. math::

 \frac{{\partial {E_{prvrd}}}}{{\partial {T_{g,\,prvrd}}}} = \frac{{{\rho _{atm}}\left( {c_a^w + \frac{{{f_{wet,\,roof}}{c_{roof}}}}{{{W_{roof}}}} + \frac{{{f_{wet,\,imprvrd}}{c_{imprvrd}}}}{{{W_{imprvrd}}}}} \right)\frac{{{c_{prvrd}}}}{{{W_{prvrd}}}}}}{{c_a^w + \frac{{{f_{wet,\,roof}}{c_{roof}}}}{{{W_{roof}}}} + \frac{{{c_{prvrd}}}}{{{W_{prvrd}}}} + \frac{{{f_{wet,\,imprvrd}}{c_{imprvrd}}}}{{{W_{imprvrd}}}}}}\frac{{d{q_{g,\,prvrd}}}}{{d{T_{g,\,prvrd}}}}

,



.. math::

 \frac{{\partial {E_{imprvrd}}}}{{\partial {T_{g,\,imprvrd}}}} = \frac{{{\rho _{atm}}\left( {c_a^w + \frac{{{f_{wet,\;roof}}{c_{roof}}}}{{{W_{roof}}}} + \frac{{{c_{prvrd}}}}{{{W_{prvrd}}}}} \right)\frac{{{f_{wet,\,imprvrd}}{c_{imprvrd}}}}{{{W_{imprvrd}}}}}}{{c_a^w + \frac{{{f_{wet,\,roof}}{c_{roof}}}}{{{W_{roof}}}} + \frac{{{c_{prvrd}}}}{{{W_{prvrd}}}} + \frac{{{f_{wet,\,imprvrd}}{c_{imprvrd}}}}{{{W_{imprvrd}}}}}}\frac{{d{q_{g,\,imprvrd}}}}{{d{T_{g,\,imprvrd}}}}

,



.. math::

 \frac{{\partial {E_{sunwall}}}}{{\partial {T_{g,\,sunwall}}}} = 0

,



.. math::

 \frac{{\partial {E_{shdwall}}}}{{\partial {T_{g,\,shdwall}}}} = 0

.

The 2-m air temperature diagnostic is set equal to the UCL air temperature :math:`{T_{ac}}` and the 2-m specific humidity diagnostic is set equal to the UCL specific humidity :math:`{q_{ac}}` . Relative humidity of the UCL air is



.. math::

 R{H_{ac}} = \min \left( {100,\,\frac{{{q_{ac}}}}{{q_{sat}^{{T_{ac}}}}} \times 100} \right)



where :math:`q_{sat}^{{T_{ac}}}` is the saturated specific humidity at UCL air temperature :math:`{T_{ac}}` (section 3.3).

The sensible heat and water vapor fluxes are based on the urban surface temperature from the previous time step, :math:`T_g^n` , and are used as the surface forcing for the solution of the soil temperature equations (section 4). This solution yields a new surface temperature :math:`T_g^{n + 1}` . The sensible heat and water vapor fluxes are updated to :math:`{H'_g}` and :math:`{E'_g}` for the new temperature as



.. math::

 {H'_g} = {H_g} + \left( {T_g^{n + 1} - T_g^n} \right)\frac{{\partial {H_g}}}{{\partial {T_g}}}

,



.. math::

 {E'_g} = {E_g} + \left( {T_g^{n + 1} - T_g^n} \right)\frac{{\partial {E_g}}}{{\partial {T_g}}}



where :math:`{H_g}` and :math:`{E_g}` are the sensible heat and water vapor fluxes derived above, and :math:`g` denotes each of the five urban surfaces. One further adjustment is made to the fluxes for the roof, pervious and impervious road. If the surface moisture (i.e., the ponded water in the case of the roof and impervious road, and top layer moisture for the pervious road) is not sufficient to supported the updated evaporation, i.e., if :math:`{E'_g} > 0` and :math:`{f_{evap}} < 1` where

\[{f_{evap}} = \frac{{{{\left( {{w_{ice,\;snl + 1}} + {w_{liq,\,snl + 1}}} \right)} \mathord{\left/

{\vphantom {{\left( {{w_{ice,\;snl + 1}} + {w_{liq,\,snl + 1}}} \right)} {\Delta t}}} \right.

\kern-\nulldelimiterspace} {\Delta t}}}}{{{{E'}_g}}} \leqslant 1\],

an adjustment is made to reduce the ground evaporation accordingly as



.. math::

 {E''_g} = {f_{evap}}{E'_g}

.

 :math:`{w_{ice,\,snl + 1}}` and :math:`{w_{liq,\,snl + 1}}` are the ice and liquid water contents (kg m-2) of the top layer. Any resulting energy deficit is assigned to sensible heat as



.. math::

 {H''_g} = {H_g} + \lambda \left( {{{E'}_g} - {{E''}_g}} \right)

.

The water vapor flux :math:`{E''_g}` is partitioned into evaporation of liquid water :math:`{q_{seva}}` , sublimation from ice :math:`{q_{subl}}` , liquid dew :math:`{q_{sdew}}` , or frost :math:`{q_{frost}}` (all in kg m-2 s-1) as



.. math::

 {q_{seva}} = \max \left( {{{E''}_g}\frac{{{w_{liq,\,snl + 1}}}}{{{w_{ice,\;snl + 1}} + {w_{liq,\,snl + 1}}}},\,0} \right) & {E''_g} \geqslant 0,\,{w_{ice,\,snl + 1}} + {w_{liq,\,snl + 1}} > 0





.. math::

 {q_{subl}} = {E''_g} - {q_{seva}} & {E''_g} \geqslant 0





.. math::

 {q_{sdew}} = \left| {{{E''}_g}} \right| & {E''_g} < 0{\text{ and }}{T_g} \geqslant {T_f}





.. math::

 {q_{frost}} = \left| {{{E''}_g}} \right| & {E''_g} < 0{\text{ and }}{T_g} < {T_f}

.

The loss or gain in snow mass due to :math:`{q_{seva}}` , :math:`{q_{subl}}` , :math:`{q_{sdew}}` , and :math:`{q_{frost}}` on a snow surface are accounted for during the snow hydrology calculations (section 5.1). The loss of surface water from non-snow surfaces due to :math:`{q_{seva}}` is accounted for in the calculation of infiltration (section 5.2), while losses or gains due to :math:`{q_{subl}}` , :math:`{q_{sdew}}` , and :math:`{q_{frost}}` from non-snow surfaces are accounted for following sub-surface drainage calculations (section 5.4).

The ground or storage heat flux :math:`G` for each urban surface is calculated as



.. math::

 G = {\overrightarrow S _g} - {\overrightarrow L _g} - {H_g} - \lambda {E_g} + {H_{wasteheat,\,g}} + {H_{aircond,\,g}}



where :math:`{\overrightarrow S _g}` is the absorbed solar radiation (section 2.5), :math:`{H_g}` and :math:`\lambda {E_g}` are the sensible and latent heat fluxes after the adjustments described above, and :math:`{H_{wasteheat,\,g}}` and :math:`{H_{aircond,\,g}}` are the waste heat and heat removed by air conditioning (section 4.1). The net longwave radiation :math:`{\overrightarrow L _g}` is updated for the change in surface temperature as



.. math::

 \vec L_g^{n + 1} = \vec L_g^n + 4{\varepsilon _g}\sigma {\left( {T_g^n} \right)^3}\left( {T_g^{n + 1} - T_g^n} \right)

.

When converting water vapor flux to an energy flux, the term :math:`\lambda` is arbitrarily assumed to be

$\lambda = \left\{ \begin{gathered}

{\lambda _{sub}} & {\text{if }}{w_{liq,\,snl + 1}} = 0{\text{ and }}{w_{ice,\,snl + 1}} > 0 \hfill \\

{\lambda _{vap}} & {\text{otherwise}} \hfill \\

\end{gathered} \right\}$

where :math:`{\lambda _{sub}}` and :math:`{\lambda _{vap}}` are the latent heat of sublimation and vaporization, respectively (J kg-1) (Table 1.4).

Saturation Specific Humidity
--------------------------------

Saturation vapor pressure :math:`e_{sat}^T` (Pa) and its derivative :math:`\frac{{de_{sat}^T}}{{dT}}` , as a function of temperature :math:`T` (ºC), are calculated from the eighth-order polynomial fits of Flatau et al. (1992)



.. math::

 e_{sat}^T = 100\left[ {{a_0} + {a_1}T + \cdots + {a_n}{T^n}} \right]





.. math::

 \frac{{de_{sat}^T}}{{dT}} = 100\left[ {{b_0} + {b_1}T + \cdots + {b_n}{T^n}} \right]



where the coefficients for ice are valid for :math:`- 75{\,^ \circ }{\text{C}} \leqslant T < 0{\,^ \circ }{\text{C}}` and the coefficients for water are valid for :math:`0{\,^ \circ }{\text{C}} \leqslant T \leqslant 100{\,^ \circ }{\text{C}}` (Table 3.1 and 3.2). The saturated water vapor specific humidity :math:`q_{sat}^T` and its derivative :math:`\frac{{dq_{sat}^T}}{{dT}}` are



.. math::

 q_{sat}^T = \frac{{0.622e_{sat}^T}}{{{P_{atm}} - 0.378e_{sat}^T}}





.. math::

 \frac{{dq_{sat}^T}}{{dT}} = \frac{{0.622{P_{atm}}}}{{{{\left( {{P_{atm}} - 0.378e_{sat}^T} \right)}^2}}}\frac{{de_{sat}^T}}{{dT}}

.

Table 3.1. Coefficients for :math:`e_{sat}^T` 

Table 3.2. Coefficients for :math:`\frac{{de_{sat}^T}}{{dT}}` 

Roof, Wall, Road, and Snow Temperatures
==========================================

The first law of heat conduction is



.. math::

 F = - \lambda \nabla T



where :math:`F` is the amount of heat conducted across a unit cross-sectional area in unit time (W m-2), :math:`\lambda` is thermal conductivity (W m-1 K-1), and :math:`\nabla T` is the spatial gradient of temperature (K m-1). In one-dimensional form



.. math::

 {F_z} = - \lambda \frac{{\partial T}}{{\partial z}}



where :math:`z` is in the vertical direction (m) and is positive downward and :math:`{F_z}` is positive upward. To account for non-steady or transient conditions, the principle of energy conservation in the form of the continuity equation is invoked as



.. math::

 c\frac{{\partial T}}{{\partial t}} = - \frac{{\partial {F_z}}}{{\partial z}}



where :math:`c` is the volumetric snow/soil heat capacity (J m-3 K-1) and :math:`t` is time (s). Combining equations and yields the second law of heat conduction in one-dimensional form



.. math::

 c\frac{{\partial T}}{{\partial t}} = \frac{\partial }{{\partial z}}\left[ {\lambda \frac{{\partial T}}{{\partial z}}} \right]

.

The nature of the solution of this equation depends on the type of urban surface. The solution for pervious and impervious roads follows the solution for CLM soils where the equation is solved numerically for a fifteen-layer column with up to five overlying layers of snow with the boundary conditions of :math:`h` as the heat flux into the surface layer from the overlying atmosphere and zero heat flux at the bottom of the soil column. In the case of pervious roads, the temperature profile is calculated first without phase change and then readjusted for phase change (section 4.2). For impervious roads, however, the moisture content of all layers is zero. Phase change then only takes place in the ponded surface water. The roof consists of a fifteen-layer column with potential ponded surface water including up to a five layer snow pack, however, the bottom boundary condition is a non-zero flux governed by prescribed controls on the internal building temperature. The walls are modeled similarly to roofs except for the absence of ponded water or snow.

Numerical Solution
----------------------

Roofs and walls are discretized into fifteen layers where the depth of layer :math:`i` , or node depth, :math:`{z_i}` (m), is



.. math::

 {z_i} = \left( {i - 0.5} \right)\left( {\frac{{\Delta z}}{{{N_{levgrnd}}}}} \right)



where :math:`\Delta z` is the total thickness of the roof or wall (Table 1.3) and :math:`{N_{levgrnd}} = 15` is the number of layers. The thickness of each layer :math:`\Delta {z_i}` (m) is

\[\Delta {z_i} = \left\{ \begin{gathered}
#. 5\left( {{z_1} + {z_2}} \right) & i = 1 \hfill \\
#. 5\left( {{z_{i + 1}} - {z_{i - 1}}} \right) & i = 2, \ldots ,{N_{levgrnd}} - 1 \hfill \\

{z_{{N_{levgrnd}}}} - {z_{{N_{levgrnd}} - 1}} & i = {N_{levgrnd}} \hfill \\

\end{gathered} \right\}\].

The depths at the layer interfaces :math:`{z_{h,\,i}}` (m) are

\[{z_{h,\,i}} = \left\{ \begin{gathered}

0 & i = 0 \hfill \\
#. 5\left( {{z_i} + {z_{i + 1}}} \right) & i = 1, \ldots ,{N_{levgrnd}} - 1 \hfill \\

{z_{{N_{levgrnd}}}} + 0.5\Delta {z_{{N_{levgrnd}}}} & i = {N_{levgrnd}} \hfill \\

\end{gathered} \right\}\].

Pervious and impervious road are discretized into fifteen layers as well with node depth



.. math::

 {z_i} = {f_s}\left\{ {\exp \left[ {0.5\left( {i - 0.5} \right)} \right] - 1} \right\}



where :math:`{f_s} = 0.025` is a scaling factor. Layer thicknesses and interface depths are calculated from equations and .

The overlying snow pack for the roof and road is modeled with up to five layers depending on the total snow depth. The layers from top to bottom are indexed in the Fortran code as :math:`i = - 4, - 3, - 2, - 1,0` , which permits the accumulation or ablation of snow at the top of the snow pack without renumbering the layers. Layer :math:`i = 0` is the snow layer next to the urban surface and layer :math:`i = snl + 1` is the top layer, where the variable :math:`snl` is the negative of the number of snow layers. The number of snow layers and the thickness of each layer is a function of snow depth :math:`{z_{sno}}` (m) as follows.

\[\left\{ \begin{gathered}

snl = - 1 \hfill \\

\Delta {z_0} = {z_{sno}} & {\text{for 0}}{\text{.01}} \leqslant {{\text{z}}_{{\text{sno}}}} \leqslant 0.03 \hfill \\

\end{gathered} \right\}\],

\[\left\{ \begin{gathered}

snl = - 2 \hfill \\

\Delta {z_{ - 1}} = {{{z_{sno}}} \mathord{\left/

{\vphantom {{{z_{sno}}} 2}} \right.

\kern-\nulldelimiterspace} 2} & {\text{for 0}}{\text{.03}}\,{\text{ < }}\,{{\text{z}}_{{\text{sno}}}} \leqslant 0.04 \hfill \\

\Delta {z_0} = \Delta {z_{ - 1}} \hfill \\

\end{gathered} \right\}\],

\[\left\{ \begin{gathered}

snl = - 2 \hfill \\

\Delta {z_{ - 1}} = 0.02 & {\text{for 0}}{\text{.04}}\,{\text{ < }}\,{{\text{z}}_{{\text{sno}}}} \leqslant 0.07 \hfill \\

\Delta {z_0} = {z_{sno}} - \Delta {z_{ - 1}} \hfill \\

\end{gathered} \right\}\],

\[\left\{ \begin{gathered}

snl = - 3 \hfill \\

\Delta {z_{ - 2}} = 0.02 \hfill \\

\Delta {z_{ - 1}} = {{\left( {{z_{sno}} - 0.02} \right)} \mathord{\left/

{\vphantom {{\left( {{z_{sno}} - 0.02} \right)} 2}} \right.

\kern-\nulldelimiterspace} 2} & {\text{for 0}}{\text{.07}}\,{\text{ < }}\,{{\text{z}}_{{\text{sno}}}} \leqslant 0.12 \hfill \\

\Delta {z_0} = \Delta {z_{ - 1}} \hfill \\

\end{gathered} \right\}\],

\[\left\{ \begin{gathered}

snl = - 3 \hfill \\

\Delta {z_{ - 2}} = 0.02 \hfill \\

\Delta {z_{ - 1}} = 0.05 & {\text{for 0}}{\text{.12}}\,{\text{ < }}\,{{\text{z}}_{{\text{sno}}}} \leqslant 0.18 \hfill \\

\Delta {z_0} = {z_{sno}} - \Delta {z_{ - 2}} - \Delta {z_{ - 1}} \hfill \\

\end{gathered} \right\}\],

\[\left\{ \begin{gathered}

snl = - 4 \hfill \\

\Delta {z_{ - 3}} = 0.02 \hfill \\

\Delta {z_{ - 2}} = 0.05 & {\text{for 0}}{\text{.18}}\,{\text{ < }}\,{{\text{z}}_{{\text{sno}}}} \leqslant 0.29 \hfill \\

\Delta {z_{ - 1}} = {{\left( {{z_{sno}} - \Delta {z_{ - 3}} - \Delta {z_{ - 2}}} \right)} \mathord{\left/

{\vphantom {{\left( {{z_{sno}} - \Delta {z_{ - 3}} - \Delta {z_{ - 2}}} \right)} 2}} \right.

\kern-\nulldelimiterspace} 2} \hfill \\

\Delta {z_0} = \Delta {z_{ - 1}} \hfill \\

\end{gathered} \right\}\],

\[\left\{ \begin{gathered}

snl = - 4 \hfill \\

\Delta {z_{ - 3}} = 0.02 \hfill \\

\Delta {z_{ - 2}} = 0.05 & {\text{for 0}}{\text{.29}}\,{\text{ < }}\,{{\text{z}}_{{\text{sno}}}} \leqslant 0.41 \hfill \\

\Delta {z_{ - 1}} = 0.11 \hfill \\

\Delta {z_0} = {z_{sno}} - \Delta {z_{ - 3}} - \Delta {z_{ - 2}} - \Delta {z_{ - 1}} \hfill \\

\end{gathered} \right\}\],

\[\left\{ \begin{gathered}

snl = - 5 \hfill \\

\Delta {z_{ - 4}} = 0.02 \hfill \\

\Delta {z_{ - 3}} = 0.05 & {\text{for 0}}{\text{.41}}\,{\text{ < }}\,{{\text{z}}_{{\text{sno}}}} \leqslant 0.64 \hfill \\

\Delta {z_{ - 2}} = 0.11 \hfill \\

\Delta {z_{ - 1}} = {{\left( {{z_{sno}} - \Delta {z_{ - 4}} - \Delta {z_{ - 3}} - \Delta {z_{ - 2}}} \right)} \mathord{\left/

{\vphantom {{\left( {{z_{sno}} - \Delta {z_{ - 4}} - \Delta {z_{ - 3}} - \Delta {z_{ - 2}}} \right)} 2}} \right.

\kern-\nulldelimiterspace} 2} \hfill \\

\Delta {z_0} = \Delta {z_{ - 1}} \hfill \\

\end{gathered} \right\}\],

\[\left\{ \begin{gathered}

snl = - 5 \hfill \\

\Delta {z_{ - 4}} = 0.02 \hfill \\

\Delta {z_{ - 3}} = 0.05 & {\text{for 0}}{\text{.64}}\,{\text{ < }}\,{{\text{z}}_{{\text{sno}}}} \hfill \\

\Delta {z_{ - 2}} = 0.11 \hfill \\

\Delta {z_{ - 1}} = 0.23 \hfill \\

\Delta {z_0} = {z_{sno}} - \Delta {z_{ - 4}} - \Delta {z_{ - 3}} - \Delta {z_{ - 2}} - \Delta {z_{ - 1}} \hfill \\

\end{gathered} \right\}\].

The node depths, which are located at the midpoint of the snow layers, and the layer interfaces are both referenced from the urban surface and are defined as negative values



.. math::

 {z_i} = {z_{h,\,i}} - 0.5\Delta {z_i} & i = snl + 1, \ldots ,0





.. math::

 {z_{h,\,i}} = {z_{h,\,i + 1}} - \Delta {z_{i + 1}} & i = snl, \ldots , - 1

.

Note that :math:`{z_{h,\,0}}` , the interface between the bottom snow layer and the top urban layer, is zero. Thermal properties (i.e., temperature :math:`{T_i}` [K]; thermal conductivity :math:`{\lambda _i}` [W m-1 K-1]; volumetric heat capacity :math:`{c_i}` [J m-3 K-1]) are defined for layers at the node depths (Figure 4.1) and for snow layers at the layer midpoints.

In general, for a zero-flux bottom boundary condition, the heat flux :math:`{F_i}` (W m-2) from layer :math:`i` to layer :math:`i + 1` is



.. math::

 {F_i} = - \lambda \left[ {{z_{h,\,i}}} \right]\left( {\frac{{{T_i} - {T_{i + 1}}}}{{{z_{i + 1}} - {z_i}}}} \right)



where the thermal conductivity at the interface :math:`\lambda \left[ {{z_{h,\,i}}} \right]` is

\[\lambda \left[ {{z_{h,\,i}}} \right] = \left\{ \begin{gathered}

\frac{{{\lambda _i}{\lambda _{i + 1}}\left( {{z_{i + 1}} - {z_i}} \right)}}{{{\lambda _i}\left( {{z_{i + 1}} - {z_{h,\,i}}} \right) + {\lambda _{i + 1}}\left( {{z_{h,\,i}} - {z_i}} \right)}} & i = snl + 1, \ldots ,{N_{levgrnd}} - 1 \hfill \\

0 & i = {N_{levgrnd}} \hfill \\

\end{gathered} \right\}\].

For a non-zero flux bottom boundary condition, 

.. math::

 \lambda \left[ {{z_{h,\,i = {N_{levgrnd}}}}} \right] = {\lambda _{i = {N_{levgrnd}}}}

. These equations are derived, with reference to Figure 4.1, assuming that the heat flux from :math:`i` (depth :math:`{z_i}` ) to the interface between :math:`i` and :math:`i + 1` (depth :math:`{z_{h,\,i}}` ) equals the heat flux from the interface to :math:`i + 1` (depth :math:`{z_{i + 1}}` ), i.e.,



.. math::

 - {\lambda _i}\frac{{{T_i} - {T_m}}}{{{z_{h,\,i}} - {z_i}}} = - {\lambda _{i + 1}}\frac{{{T_m} - {T_{i + 1}}}}{{{z_{i + 1}} - {z_{h,\,i}}}}



where :math:`{T_m}` is the temperature at the interface of layers :math:`i` and :math:`i + 1` . Solving equation for :math:`{T_m}` and substituting :math:`{T_m}` back into the left side of equation yields equations and .

Figure 4.1. Schematic diagram of numerical scheme used to solve for layer temperatures. Shown are three layers, :math:`i - 1` , :math:`i` , and :math:`i + 1` . The thermal conductivity :math:`\lambda` , specific heat capacity :math:`c` , and temperature :math:`T` are defined at the layer node depth :math:`z` . :math:`{T_m}` is the interface temperature. The thermal conductivity :math:`\lambda \left[ {{z_h}} \right]` is defined at the interface of two layers :math:`{z_h}` . The layer thickness is :math:`\Delta z` . The heat fluxes :math:`{F_{i - 1}}` and :math:`{F_i}` are defined as positive upwards.

The energy balance for the :math:`{i^{th}}` layer is



.. math::

 \frac{{{c_i}\Delta {z_i}}}{{\Delta t}}\left( {T_i^{n + 1} - T_i^n} \right) = - {F_{i - 1}} + {F_i}



where the superscripts :math:`n` and :math:`n + 1` indicate values at the beginning and end of the time step, respectively, and :math:`\Delta t` is the time step (s). This equation is solved using the Crank-Nicholson method, which combines the explicit method with fluxes evaluated at :math:`n` ( :math:`F_{i - 1}^n,F_i^n` ) and the implicit method with fluxes evaluated at :math:`n + 1` ( :math:`F_{i - 1}^{n + 1},F_i^{n + 1}` )



.. math::

 \frac{{{c_i}\Delta {z_i}}}{{\Delta t}}\left( {T_i^{n + 1} - T_i^n} \right) = \alpha \left( { - F_{i - 1}^n + F_i^n} \right) + \left( {1 - \alpha } \right)\left( { - F_{i - 1}^{n + 1} + F_i^{n + 1}} \right)



where :math:`\alpha = 0.5` , resulting in a tridiagonal system of equations



.. math::

 {r_i} = {a_i}T_{i - 1}^{n + 1} + {b_i}T_i^{n + 1} + {c_i}T_{i + 1}^{n + 1}



where :math:`{a_i}` , :math:`{b_i}` , and :math:`{c_i}` are the subdiagonal, diagonal, and superdiagonal elements in the tridiagonal matrix and :math:`{r_i}` is a column vector of constants.

For the top layer :math:`i = snl + 1` , the heat flux from the overlying atmosphere into the surface layer 

.. math::

 h

 (W m-2, defined as positive into the surface) is



.. math::

 {h^{n + 1}} = - \alpha F_{i - 1}^n - \left( {1 - \alpha } \right)F_{i - 1}^{n + 1}

.

The energy balance for layer :math:`i = snl + 1` is then



.. math::

 \frac{{{c_i}\Delta {z_i}}}{{\Delta t}}\left( {T_i^{n + 1} - T_i^n} \right) = {h^{n + 1}} + \alpha F_i^n + \left( {1 - \alpha } \right)F_i^{n + 1}

.

The heat flux 

.. math::

 h

 at :math:`n + 1` may be approximated as follows



.. math::

 {h^{n + 1}} = {h^n} + \frac{{\partial h}}{{\partial {T_i}}}\left( {T_i^{n + 1} - T_i^n} \right)

.

The resulting equations are

\[\begin{gathered}

\frac{{{c_i}\Delta {z_i}}}{{\Delta t}}\left( {T_i^{n + 1} - T_i^n} \right) = {h^n} + \frac{{\partial h}}{{\partial {T_i}}}\left( {T_i^{n + 1} - {T_i}} \right) \\
* \alpha \frac{{\lambda \left[ {{z_{h,\,i}}} \right]\left( {T_i^n - T_{i + 1}^n} \right)}}{{{z_{i + 1}} - {z_i}}} - \left( {1 - \alpha } \right)\frac{{\lambda \left[ {{z_{h,\,i}}} \right]\left( {T_i^{n + 1} - T_{i + 1}^{n + 1}} \right)}}{{{z_{i + 1}} - {z_i}}} \\

\end{gathered} \]



.. math::

 {a_i} = 0





.. math::

 {b_i} = 1 + \frac{{\Delta t}}{{{c_i}\Delta {z_i}}}\left[ {\left( {1 - \alpha } \right)\frac{{\lambda \left[ {{z_{h,\,i}}} \right]}}{{{z_{i + 1}} - {z_i}}} - \frac{{\partial h}}{{\partial {T_i}}}} \right]





.. math::

 {c_i} = - \left( {1 - \alpha } \right)\frac{{\Delta t}}{{{c_i}\Delta {z_i}}}\frac{{\lambda \left[ {{z_{h,\,i}}} \right]}}{{{z_{i + 1}} - {z_i}}}





.. math::

 {r_i} = T_i^n + \frac{{\Delta t}}{{{c_i}\Delta {z_i}}}\left[ {{h^n} - \frac{{\partial h}}{{\partial {T_i}}}T_i^n + \alpha {F_i}} \right]



where



.. math::

 {F_i} = - \lambda \left[ {{z_{h,\,i}}} \right]\left( {\frac{{T_i^n - T_{i + 1}^n}}{{{z_{i + 1}} - {z_i}}}} \right)

.

The heat flux into each urban surface 

.. math::

 h

 is



.. math::

 h = {\overrightarrow S _g} - {\overrightarrow L _g} - {H_g} - \lambda {E_g} + {H_{wasteheat,\,g}} + {H_{aircond,\,g}}



where :math:`{\overrightarrow S _g}` is the absorbed solar radiation (section 2.5), :math:`{\overrightarrow L _g}` is the net longwave radiation (section 2.7), and :math:`{H_g}` and :math:`\lambda {E_g}` are the sensible and latent heat fluxes (section 3.2). The terms :math:`{H_{wasteheat,\,g}}` and :math:`{H_{aircond,\,g}}` are the waste heat from space heating/air conditioning and heat removed by air conditioning applied only to the pervious and impervious road

\[\begin{gathered}

{H_{wasteheat,\,prvrd}} = {H_{wasteheat,\,imprvrd}} = \frac{{{H_{wasteheat}}}}{{1 - {W_{roof}}}} \hfill \\

{H_{wasteheat,\,sunwall}} = {H_{wasteheat,\,shdwall}} = {H_{wasteheat,\,roof}} = 0 \hfill \\

{H_{aircond,\,prvrd}} = {H_{aircond,\,imprvrd}} = \frac{{{H_{aircond}}}}{{1 - {W_{roof}}}} \hfill \\

{H_{aircond,\,sunwall}} = {H_{aircond,\,shdwall}} = {H_{aircond,\,roof}} = 0 \hfill \\

\end{gathered} \].

where :math:`{H_{wasteheat}}` and :math:`{H_{aircond}}` are the total waste heat and heat removed by air conditioning from equations and . Note that for the pervious road, the latent heat is always the total latent heat regardless of its partitioning into ground evaporation or transpiration (section 3.2.4). The partial derivative of the heat flux 

.. math::

 h

 with respect to surface temperature is



.. math::

 \frac{{\partial h}}{{\partial {T_g}}} = - \frac{{\partial {{\overrightarrow L }_g}}}{{\partial {T_g}}} - \frac{{\partial {H_g}}}{{\partial {T_g}}} - \frac{{\partial \lambda {E_g}}}{{\partial {T_g}}}



where the partial derivative of the net longwave radiation is



.. math::

 \frac{{\partial {{\overrightarrow L }_g}}}{{\partial {T_g}}} = 4{\varepsilon _g}\sigma {\left( {T_g^n} \right)^3}



and the partial derivatives of the sensible and latent heat fluxes are given by equations -. :math:`\sigma` is the Stefan-Boltzmann constant (W m-2 K-4) (Table 1.4) and :math:`{\varepsilon _g}` is the surface emissivity (section 2.7).

The top layer for roofs and walls is thin enough such that the layer-averaged temperature calculated above is considered to be equivalent to the surface temperature :math:`T_g^{n + 1}` . For pervious and impervious road, the top layer temperature has somewhat reduced diurnal amplitude compared with surface temperature. An accurate surface temperature is provided that compensates for this effect and numerical error by tuning the heat capacity of the top layer (through adjustment of the layer thickness) to give an exact match to the analytic solution for diurnal heating. The layer thickness for :math:`i = snl + 1` is given by



.. math::

 \Delta {z_{i*}} = 0.5\left[ {{z_i} - {z_{h,\,i - 1}} + {c_a}\left( {{z_{i + 1}} - {z_{h,\,i - 1}}} \right)} \right]



where :math:`{c_a}` is a tunable parameter, varying from 0 to 1, and is taken as 0.34 by comparing the numerical solution with the analytic solution (Z.-L. Yang 1998, unpublished manuscript). For pervious and impervious road, :math:`\Delta {z_{i*}}` is used in place of :math:`\Delta {z_i}` for :math:`i = snl + 1` in equations -.

For the pervious and impervious road, the boundary condition at the bottom is zero heat flux, :math:`{F_i} = 0` , resulting in, for :math:`i = {N_{levgrnd}}` ,



.. math::

 \frac{{{c_i}\Delta {z_i}}}{{\Delta t}}\left( {T_i^{n + 1} - T_i^n} \right) = \alpha \frac{{\lambda \left[ {{z_{h,\,i - 1}}} \right]\left( {T_{i - 1}^n - T_i^n} \right)}}{{{z_i} - {z_{i - 1}}}} + \left( {1 - \alpha } \right)\frac{{\lambda \left[ {{z_{h,\,i - 1}}} \right]\left( {T_{i - 1}^{n + 1} - T_i^{n + 1}} \right)}}{{{z_i} - {z_{i - 1}}}}





.. math::

 {a_i} = - \left( {1 - \alpha } \right)\frac{{\Delta t}}{{{c_i}\Delta {z_i}}}\frac{{\lambda \left[ {{z_{h,\,i - 1}}} \right]}}{{{z_i} - {z_{i - 1}}}}





.. math::

 {b_i} = 1 + \left( {1 - \alpha } \right)\frac{{\Delta t}}{{{c_i}\Delta {z_i}}}\frac{{\lambda \left[ {{z_{h,\,i - 1}}} \right]}}{{{z_i} - {z_{i - 1}}}}





.. math::

 {c_i} = 0





.. math::

 {r_i} = T_i^n - \alpha \frac{{\Delta t}}{{{c_i}\Delta {z_i}}}{F_{i - 1}}



where



.. math::

 {F_{i - 1}} = - \frac{{\lambda \left[ {{z_{h,\,i - 1}}} \right]}}{{{z_i} - {z_{i - 1}}}}\left( {T_{i - 1}^n - T_i^n} \right)

.

For the roof and walls, the boundary condition at the bottom is the internal building temperature :math:`{T_{iB}}` , constrained as :math:`{T_{iB,\,\max }} \geqslant {T_{iB}} \geqslant {T_{iB,\,\min }}` , where :math:`{T_{iB,\,\max }}` and :math:`{T_{iB,\,\min }}` are prescribed maximum and minimum internal building temperatures (Table 1.3). The internal building temperature :math:`T_{iB}^{}` is determined from a weighted combination of the inner layer wall and roof temperatures as



.. math::

 {T_{iB}} = \frac{{H\left( {T_{i = {N_{levgrnd}},\,shdwall}^n + T_{i = {N_{levgrnd}},\,sunwall}^n} \right) + {L_{roof}}T_{i = {N_{levgrnd}},\,roof}^n}}{{2H + {L_{roof}}}}



where :math:`H` is the building height and :math:`{L_{roof}}` is the length of the roof in an infinite canyon configuration

\[L = \left( {\frac{H}{{{H \mathord{\left/

{\vphantom {H W}} \right.

\kern-\nulldelimiterspace} W}}}} \right)\left( {\frac{{{W_{roof}}}}{{1 - {W_{roof}}}}} \right)\].

This boundary condition yields, for :math:`i = {N_{levgrnd}}` ,

\[\begin{gathered}

\frac{{{c_i}\Delta {z_i}}}{{\Delta t}}\left( {T_i^{n + 1} - T_i^n} \right) = - \alpha \frac{{\lambda \left[ {{z_{h,\,i}}} \right]\left( {T_i^n - T_{i + 1}^n} \right)}}{{{z_{h,\,i}} - {z_i}}} + \alpha \frac{{\lambda \left[ {{z_{h,\,i - 1}}} \right]\left( {T_{i - 1}^n - T_i^n} \right)}}{{{z_i} - {z_{i - 1}}}} \\
* \left( {1 - \alpha } \right)\frac{{\lambda \left[ {{z_{h,\,i}}} \right]\left( {T_i^{n + 1} - T_{i + 1}^{n + 1}} \right)}}{{{z_{h,\,i}} - {z_i}}} + \left( {1 - \alpha } \right)\frac{{\lambda \left[ {{z_{h,\,i - 1}}} \right]\left( {T_{i - 1}^{n + 1} - T_i^{n + 1}} \right)}}{{{z_i} - {z_{i - 1}}}} \\

\end{gathered} \]



.. math::

 {a_i} = - \left( {1 - \alpha } \right)\frac{{\Delta t}}{{{c_i}\Delta {z_i}}}\frac{{\lambda \left[ {{z_{h,\,i - 1}}} \right]}}{{{z_i} - {z_{i - 1}}}}





.. math::

 {b_i} = 1 + \left( {1 - \alpha } \right)\frac{{\Delta t}}{{{c_i}\Delta {z_i}}}\left[ {\frac{{\lambda \left[ {{z_{h,\,i - 1}}} \right]}}{{{z_i} - {z_{i - 1}}}} + \frac{{\lambda \left[ {{z_{h,\,i}}} \right]}}{{{z_{h,\,i}} - {z_i}}}} \right]





.. math::

 {c_i} = 0





.. math::

 {r_i} = T_i^n + \alpha \frac{{\Delta t}}{{{c_i}\Delta {z_i}}}\left( {{F_i} - \alpha {F_{i - 1}}} \right)



where



.. math::

 {F_i} = - \lambda \left[ {{z_{h,\,i}}} \right]\left( {\frac{{\alpha T_i^n - T_{iB}^n}}{{{z_{h,\,i}} - {z_i}}}} \right)

,



.. math::

 {F_{i - 1}} = - \frac{{\lambda \left[ {{z_{h,\,i - 1}}} \right]}}{{{z_i} - {z_{i - 1}}}}\left( {T_{i - 1}^n - T_i^n} \right)

.

For the interior snow/soil layers of all surfaces, :math:`snl + 1 < i < {N_{nlevgrnd}}` ,

\[\begin{gathered}

\frac{{{c_i}\Delta {z_i}}}{{\Delta t}}\left( {T_i^{n + 1} - T_i^n} \right) = - \alpha \frac{{\lambda \left[ {{z_{h,\,i}}} \right]\left( {T_i^n - T_{i + 1}^n} \right)}}{{{z_{i + 1}} - {z_i}}} + \alpha \frac{{\lambda \left[ {{z_{h,\,i - 1}}} \right]\left( {T_{i - 1}^n - T_i^n} \right)}}{{{z_i} - {z_{i - 1}}}} \\
* \left( {1 - \alpha } \right)\frac{{\lambda \left[ {{z_{h,\,i}}} \right]\left( {T_i^{n + 1} - T_{i + 1}^{n + 1}} \right)}}{{{z_{i + 1}} - {z_i}}} + \left( {1 - \alpha } \right)\frac{{\lambda \left[ {{z_{h,\,i - 1}}} \right]\left( {T_{i - 1}^{n + 1} - T_i^{n + 1}} \right)}}{{{z_i} - {z_{i - 1}}}} \\

\end{gathered} \]



.. math::

 {a_i} = - \left( {1 - \alpha } \right)\frac{{\Delta t}}{{{c_i}\Delta {z_i}}}\frac{{\lambda \left[ {{z_{h,\,i - 1}}} \right]}}{{{z_i} - {z_{i - 1}}}}





.. math::

 {b_i} = 1 + \left( {1 - \alpha } \right)\frac{{\Delta t}}{{{c_i}\Delta {z_i}}}\left[ {\frac{{\lambda \left[ {{z_{h,\,i - 1}}} \right]}}{{{z_i} - {z_{i - 1}}}} + \frac{{\lambda \left[ {{z_{h,\,i}}} \right]}}{{{z_{i + 1}} - {z_i}}}} \right]





.. math::

 {c_i} = - \left( {1 - \alpha } \right)\frac{{\Delta t}}{{{c_i}\Delta {z_i}}}\frac{{\lambda \left[ {{z_{h,\,i}}} \right]}}{{{z_{i + 1}} - {z_i}}}





.. math::

 {r_i} = T_i^n + \alpha \frac{{\Delta t}}{{{c_i}\Delta {z_i}}}\left( {{F_i} - {F_{i - 1}}} \right)

.

The heating or cooling flux applied to the roof, and sunlit and shaded wall is

\[\begin{gathered}

{F_{heat}} = \left\{ \begin{gathered}

\left| {\alpha F_{i = {N_{levgrnd}}}^n + \left( {1 - \alpha } \right)F_{i = {N_{levgrnd}}}^{n + 1}} \right| & {T_{iB}} < {T_{\min }} \hfill \\

0 & & {T_{iB}} \geqslant {T_{\min }} \hfill \\

\end{gathered} \right\} \hfill \\

\hfill \\

\end{gathered} \]

\[{F_{cool}} = \left\{ \begin{gathered}

\left| {\alpha F_{i = {N_{levgrnd}}}^n + \left( {1 - \alpha } \right)F_{i = {N_{levgrnd}}}^{n + 1}} \right| & {T_{iB}} > {T_{\min }} \hfill \\

0 & {T_{iB}} \leqslant {T_{\min }} \hfill \\

\end{gathered} \right\}\]

where



.. math::

 F_{i = {N_{levgrnd}}}^n = - \frac{{\lambda \left[ {{z_{h,\,i = {N_{levgrnd}}}}} \right]}}{{{z_{h,\,i = {N_{levgrnd}}}} - {z_{i = {N_{levgrnd}}}}}}\left( {T_{i = {N_{levgrnd}}}^n - {T_{iB}}} \right)





.. math::

 F_{i = {N_{levgrnd}}}^{n + 1} = - \frac{{\lambda \left[ {{z_{h,\,i = {N_{levgrnd}}}}} \right]}}{{{z_{h,\,i = {N_{levgrnd}}}} - {z_{i = {N_{levgrnd}}}}}}\left( {T_{i = {N_{levgrnd}}}^{n + 1} - {T_{iB}}} \right)

.

The total waste heat from space heating/air conditioning is

\[\begin{gathered}

{H_{wasteheat}} = {W_{roof}}\left( {{f_{\,heat}}{F_{heat,\,roof}} + {f_{cool}}{F_{cool,\,roof}}} \right) + \hfill \\

& \left( {1 - {W_{roof}}} \right)\frac{H}{W}\left( \begin{gathered}

{f_{\,heat}}{F_{heat,\,sunwall}} + {f_{cool}}{F_{cool,\,sunwall}} + \hfill \\

{f_{heat}}{F_{heat,\,shdwall}} + {f_{cool}}{F_{cool,\,shdwall}} \hfill \\

\end{gathered} \right) \leqslant {H_{wasteheat,\,\max }} \hfill \\

& \hfill \\

\end{gathered} \]

where ${f_{\,heat}} = {1 \mathord{\left/

{\vphantom {1 {0.75}}} \right.

\kern-\nulldelimiterspace} {0.75}} :math:`and` {f_{cool}} = {1 \mathord{\left/

{\vphantom {1 {0.25}}} \right.

\kern-\nulldelimiterspace} {0.25}} :math:`are factors describing the efficiency of space heating/air conditioning systems and` {H_{wasteheat,\,\max }} = 100$ W m-2 is a maximum limit on waste heat at any given time step. The heat removed by air conditioning is



.. math::

 {H_{aircond}} = {F_{cool}}

.

Phase Change
----------------

Phase change may take place in any snow/soil layers of the pervious road and in the ponded water on roofs and impervious road. Note that the ponded water is treated as part of the top layer. Upon solution of the tridiagonal equation set (Press et al. 1992), the temperatures are evaluated to determine if phase change will take place as

$\begin{gathered}

T_i^{n + 1} > {T_f}{\text{ and }}{w_{ice,\,i}} > 0 & i = snl + 1, \ldots ,{N_{levgrnd}} & {\text{melting}} \hfill \\

T_i^{n + 1} < {T_f}{\text{ and }}{w_{liq,\,i}} > 0 & i = snl + 1, \ldots ,0 & {\text{freezing}} \hfill \\

T_i^{n + 1} < {T_f}{\text{ and }}{w_{liq,\,i}} > {w_{liq,\,\max ,\,i}} & i = 1, \ldots ,{N_{levgrnd}} & {\text{freezing}} \hfill \\

\end{gathered} $

where :math:`T_i^{n + 1}` is the layer temperature after solution of the tridiagonal equation set, :math:`{w_{ice,\,i}}` and :math:`{w_{liq,\,i}}` are the mass of ice and liquid water (kg m-2) in each layer, respectively, and :math:`{T_f}` is the freezing temperature of water (K) (Table 1.4). For the freezing process in the layers of the pervious road, the concept of supercooled soil water from Niu and Yang (2006) is adopted. The supercooled soil water is the liquid water that coexists with ice over a wide range of temperatures below freezing and is implemented through a freezing point depression equation

\[{w_{liq,\,\max ,\,i}} = \Delta {z_i}{\theta _{sat,\,i}}{\left[ {\frac{{{{10}^3}{L_f}\left( {{T_f} - {T_i}} \right)}}{{g{T_i}{\psi _{sat,\,i}}}}} \right]^{{{ - 1} \mathord{\left/

{\vphantom {{ - 1} {{B_i}}}} \right.

\kern-\nulldelimiterspace} {{B_i}}}}} & {T_i} < {T_f}\]

where :math:`{w_{liq,\,\max ,\,i}}` is the maximum liquid water in layer :math:`i` (kg m-2) when the soil temperature :math:`{T_i}` is below the freezing temperature :math:`{T_f}` , :math:`{L_f}` is the latent heat of fusion (J kg-1) (Table 1.4), :math:`g` is the gravitational acceleration (m s-2) (Table 1.4), and :math:`{\psi _{sat,\,i}}` and :math:`{B_i}` are the soil texture-dependent saturated matric potential (mm) and Clapp and Hornberger (1978) exponent (section 5.3.1). Equation applies to pervious road only, for roof and impervious road :math:`{w_{liq,\,\max ,\,i}} = 0` .

For the special case when snow is present (snow mass :math:`{W_{sno}} > 0` ) but there are no explicit snow layers ( :math:`snl = 0` ) (i.e., there is not enough snow present to meet the minimum snow depth requirement of 0.01 m), snow melt will take place for soil layer :math:`i = 1` if the soil layer temperature is greater than the freezing temperature ( :math:`T_1^{n + 1} > {T_f}` ).

The rate of phase change is assessed from the energy excess (or deficit) needed to change :math:`{T_i}` to freezing temperature, :math:`{T_f}` . The excess or deficit of energy :math:`{H_i}` (W m-2) is determined as follows

\[{H_i} = \left\{ \begin{gathered}

h + \frac{{\partial h}}{{\partial T}}\left( {{T_f} - T_i^n} \right) + \alpha F_i^n + \left( {1 - \alpha } \right)F_i^{n + 1} \hfill \\
* \frac{{{c_i}\Delta {z_i}}}{{\Delta t}}\left( {{T_f} - T_i^n} \right) & i = snl + 1 \hfill \\

\alpha \left( {F_i^n - F_{i - 1}^n} \right) + \left( {1 - \alpha } \right)\left( {F_i^{n + 1} - F_{i - 1}^{n + 1}} \right) \hfill \\
* \frac{{{c_i}\Delta {z_i}}}{{\Delta t}}\left( {{T_f} - T_i^n} \right) & i = snl + 2, \ldots ,{N_{levgrnd}} \hfill \\

\end{gathered} \right\}\]

where :math:`F_i^{n + 1}` and :math:`F_{i - 1}^{n + 1}` are calculated from equations and using :math:`T_i^{n + 1}` . For roof and walls, :math:`F_{i = {N_{levgrnd}}}^{n + 1}` is calculated from equation . If the melting criteria is met (equation ) and :math:`{H_m} = \frac{{{H_i}\Delta t}}{{{L_f}}} > 0` , then the ice mass is readjusted as



.. math::

 w_{ice,\,i}^{n + 1} = w_{ice,\,i}^n - {H_m} \geqslant 0 & i = snl + 1, \ldots ,{N_{levgrnd}}

.

If the freezing criteria is met (equation ) and :math:`{H_m} < 0` , then the ice mass is readjusted for 

.. math::

 i = snl + 1, \ldots ,0

 as



.. math::

 w_{ice,\,i}^{n + 1} = \min \left( {w_{liq,\,i}^n + w_{ice,\,i}^n,w_{ice,\,i}^n - {H_m}} \right)



and for 

.. math::

 i = 1, \ldots ,{N_{levgrnd}}

 as

\[w_{ice,\,i}^{n + 1} = \left\{ \begin{gathered}

\min \left( {w_{liq,\,i}^n + w_{ice,\,i}^n - w_{liq,\,\max ,\,i}^n,\,w_{ice,\,i}^n - {H_m}} \right) & w_{liq,\,i}^n + w_{ice,\,i}^n \geqslant w_{liq,\,\max ,\,i}^n{\text{ }} \hfill \\

{\text{0}} & w_{liq,\,i}^n + w_{ice,\,i}^n < w_{liq,\,\max ,\,i}^n{\text{ }}\, \hfill \\

\end{gathered} \right\}\]

Liquid water mass is readjusted as



.. math::

 w_{liq,\,i}^{n + 1} = w_{liq,\,i}^n + w_{ice,\,i}^n - w_{ice,\,i}^{n + 1} \geqslant 0

.

Because part of the energy :math:`{H_i}` may not be consumed in melting or released in freezing, the energy is recalculated as



.. math::

 {H_{i * }} = {H_i} - \frac{{{L_f}\left( {w_{ice,\,i}^n - w_{ice,\,i}^{n + 1}} \right)}}{{\Delta t}}



and this energy is used to cool or warm the layer (if :math:`\left| {{H_{i * }}} \right| > 0` ) as

\[T_i^{n + 1} = \left\{ \begin{gathered}

{T_f} + {{\frac{{\Delta t}}{{{c_i}\Delta {z_i}}}{H_{i * }}} \mathord{\left/

{\vphantom {{\frac{{\Delta t}}{{{c_i}\Delta {z_i}}}{H_{i * }}} {\left( {1 - \frac{{\Delta t}}{{{c_i}\Delta {z_i}}}\frac{{\partial h}}{{\partial T}}} \right)}}} \right.

\kern-\nulldelimiterspace} {\left( {1 - \frac{{\Delta t}}{{{c_i}\Delta {z_i}}}\frac{{\partial h}}{{\partial T}}} \right)}} & i = snl + 1 \hfill \\

{T_f} + \frac{{\Delta t}}{{{c_i}\Delta {z_i}}}{H_{i * }} & i = snl + 2, \ldots {N_{levgrnd}} \hfill \\

\end{gathered} \right\}\].

For the special case when snow is present ( :math:`{W_{sno}} > 0` ), there are no explicit snow layers ( :math:`snl = 0` ), and :math:`\frac{{{H_1}\Delta t}}{{{L_f}}} > 0` (melting), the snow mass :math:`{W_{sno}}` (kg m-2) is reduced according to



.. math::

 W_{sno}^{n + 1} = W_{sno}^n - \frac{{{H_1}\Delta t}}{{{L_f}}} \geqslant 0

.

The snow depth is reduced proportionally



.. math::

 z_{sno}^{n + 1} = \frac{{W_{sno}^{n + 1}}}{{W_{sno}^n}}z_{sno}^n

.

Again, because part of the energy may not be consumed in melting, the energy for the surface layer :math:`i = 1` is recalculated as



.. math::

 {H_{1 * }} = {H_1} - \frac{{{L_f}\left( {W_{sno}^n - W_{sno}^{n + 1}} \right)}}{{\Delta t}}

.

If there is excess energy ( :math:`{H_{1 * }} > 0` ), this energy becomes available to the top layer as



.. math::

 {H_1} = {H_{1 * }}

.

The ice mass, liquid water content, and temperature of the top layer are then determined from equations , , and using the recalculated energy from equation . Snow melt :math:`{M_{1S}}` (kg m-2 s-1) and phase change energy :math:`{E_{p,\,1S}}` (W m-2) for this special case are



.. math::

 {M_{1S}} = \frac{{W_{sno}^n - W_{sno}^{n + 1}}}{{\Delta t}} \geqslant 0





.. math::

 {E_{p,\,1S}} = {L_f}{M_{1S}}

.

The total energy of phase change :math:`{E_p}` (W m-2) for the column is



.. math::

 {E_p} = {E_{p,\,1S}} + \sum\limits_{i = snl + 1}^{i = {N_{levgrnd}}} {{E_{p,i}}}



where



.. math::

 {E_{p,\,i}} = {L_f}\frac{{\left( {w_{ice,\,i}^n - w_{ice,\,i}^{n + 1}} \right)}}{{\Delta t}}

.

The total snow melt :math:`M` (kg m-2 s-1) is



.. math::

 M = {M_{1S}} + \sum\limits_{i = snl + 1}^{i = 0} {{M_i}}



where



.. math::

 {M_i} = \frac{{\left( {w_{ice,\,i}^n - w_{ice,\,i}^{n + 1}} \right)}}{{\Delta t}} \geqslant 0

.

The solution for temperatures conserves energy as



.. math::

 G - {E_p} - \sum\limits_{i = snl + 1}^{i = {N_{levgrnd}}} {\frac{{{c_i}\Delta {z_i}}}{{\Delta t}}} \left( {T_i^{n + 1} - T_i^n} \right) + \left[ {\alpha F_{i = {N_{levgrnd}}}^n + \left( {1 - \alpha } \right)F_{i = {N_{levgrnd}}}^{n + 1}} \right] = 0



where :math:`G` is the ground heat flux (section 3.2.4) and the last term is the non-zero flux bottom boundary condition (roofs and walls only).

Thermal Properties
----------------------

The thermal conductivities and heat capacities for roofs, walls, and :math:`i = 1, \ldots ,{N_{imprvrd}}` layers of the impervious road are specified by the surface dataset as described in section 1.2.2 and Table 1.3. The :math:`i = {N_{imprvrd}} + 1, \ldots ,{N_{levgrnd}}` layers of impervious road and the pervious road layers consist of soil or bedrock whose thermal properties are described below. In CLM4, organic matter modifies soil properties according to Lawrence and Slater (2008). Urban soils are assumed to have no organic matter so the equations below are shown in their reduced form. Note that the moisture content of the impervious road soil layers is maintained at zero.

Soil thermal conductivity :math:`{\lambda _i}` (W m-1 K-1) is from Farouki (1981)

\[\begin{gathered}

{\lambda _i} = \left\{ \begin{gathered}

{K_{e,\,i}}{\lambda _{sat,\,i}} + \left( {1 - {K_{e,\,i}}} \right){\lambda _{dry,\,i}} & {S_{r,\,i}} > 1 \times {10^{ - 7}} \hfill \\

{\lambda _{dry,\,i}} & {S_{r,\,i}} \leqslant 1 \times {10^{ - 7}} \hfill \\

\end{gathered} \right\} & i = 1, \ldots ,{N_{levsoi}} \hfill \\

{\lambda _i} = {\lambda _{bedrock}} & i = {N_{levsoi}} + 1, \ldots {N_{levgrnd}} \hfill \\

\end{gathered} \]

where :math:`{\lambda _{sat,\,i}}` is the saturated thermal conductivity, :math:`{\lambda _{dry,\,i}}` is the dry thermal conductivity, :math:`{K_{e,\,i}}` is the Kersten number, :math:`{S_{r,\,i}}` is the wetness of the soil with respect to saturation, and :math:`{\lambda _{bedrock}} = 3` W m-1 K-1 is the thermal conductivity assumed for the deep ground layers (typical of saturated granitic rock; Clauser and Huenges, 1995). The saturated thermal conductivity :math:`{\lambda _{sat,\,i}}` (W m-1 K-1) depends on the thermal conductivities of the soil solid, liquid water, and ice constituents

\[{\lambda _{sat,\,i}} = \left\{ \begin{gathered}

\lambda _{s,\,i}^{1 - {\theta _{sat,\,i}}}\lambda _{liq}^{{\theta _{sat,\,i}}} & {T_i} \geqslant {T_f} \hfill \\

\lambda _{s,\,i}^{1 - {\theta _{sat,\,i}}}\lambda _{liq}^{{\theta _{sat,\,i}}}\lambda _{ice}^{{\theta _{sat,\,i}} - {\theta _{liq,\,i}}} & {T_i} < {T_f} \hfill \\

\end{gathered} \right\}\]

where the thermal conductivity of soil solids :math:`{\lambda _{s,\,i}}` varies with the sand and clay content



.. math::

 {\lambda _{s,\,i}} = \frac{{8.80{\text{ }}{{\left( {\% sand} \right)}_i} + {\text{2}}{\text{.92 }}{{\left( {\% clay} \right)}_i}}}{{{{\left( {\% sand} \right)}_i} + {{\left( {\% clay} \right)}_i}}}

,

and :math:`{\theta _{sat,\,i}}` is the volumetric water content at saturation (porosity) (section 5.3.1). The thermal conductivity of dry natural soil :math:`{\lambda _{dry,\,i}}` (W m-1 K-1) depends on the bulk density :math:`{\rho _{d,\,i}} = 2700\left( {1 - {\theta _{sat,\,i}}} \right)` (kg m-3) as



.. math::

 {\lambda _{dry,\,i}} = \frac{{0.135{\rho _{d,\,i}} + 64.7}}{{2700 - 0.947{\rho _{d,\,i}}}}

.

The Kersten number :math:`{K_{e,\,i}}` is a function of the degree of saturation :math:`{S_r}` and phase of water

\[{K_{e,\,i}} = \left\{ \begin{gathered}

\log \left( {{S_{r,\,i}}} \right) + 1 \geqslant 0 & {T_i} \geqslant {T_f} \hfill \\

{S_{r,\,i}} & {T_i} < {T_f} \hfill \\

\end{gathered} \right\}\]

where

 :math:`{S_{r,\,i}} = \left( {\frac{{{w_{liq,\,i}}}}{{{\rho _{liq}}\Delta {z_i}}} + \frac{{{w_{ice,\,i}}}}{{{\rho _{ice}}\Delta {z_i}}}} \right)\frac{1}{{{\theta _{sat,\,i}}}} = \frac{{{\theta _{liq,\,i}} + {\theta _{ice,\,i}}}}{{{\theta _{sat,\,i}}}} \leqslant 1` .

Thermal conductivity :math:`{\lambda _i}` (W m-1 K-1) for snow is from (1991)



.. math::

 {\lambda _i} = {\lambda _{air}} + \left( {7.75 \times {{10}^{ - 5}}{\rho _{sno,\,i}} + 1.105 \times {{10}^{ - 6}}\rho _{sno,\,i}^2} \right)\left( {{\lambda _{ice}} - {\lambda _{air}}} \right)



where :math:`{\lambda _{air}}` and :math:`{\lambda _{ice}}` are the thermal conductivities of air and ice (Table 1.4) and :math:`{\rho _{sno,\,i}}` is the bulk density of snow (kg m-3)



.. math::

 {\rho _{sno,\,i}} = \frac{{{w_{ice,\,i}} + {w_{liq,\,i}}}}{{\Delta {z_i}}}

.

The volumetric heat capacity :math:`{c_i}` (J m-3 K-1) for soil is from de Vries (1963) and depends on the heat capacities of the soil solid, liquid water, and ice constituents



.. math::

 {c_i} = {c_{s,\,i}}\left( {1 - {\theta _{sat,\,i}}} \right) + \frac{{{w_{ice,\,i}}}}{{\Delta {z_i}}}{C_{ice}} + \frac{{{w_{liq,\,i}}}}{{\Delta {z_i}}}{C_{liq}}



where the heat capacity of soil solids :math:`{c_{s,\,i}}` (J m-3 K-1) is

\[\begin{gathered}

{c_{s,\,i}} = \left( {\frac{{2.128{\text{ }}{{\left( {\% sand} \right)}_i} + {\text{2}}{\text{.385 }}{{\left( {\% clay} \right)}_i}}}{{{{\left( {\% sand} \right)}_i} + {{\left( {\% clay} \right)}_i}}}} \right) \times {10^6} & i = 1, \ldots ,{N_{levsoi}} \hfill \\

{c_{s,\,i}} = {c_{s,\,bedrock}} & i = {N_{levsoi}} + 1, \ldots ,{N_{levgrnd}} \hfill \\

\end{gathered} \]

and :math:`{C_{liq}}` and :math:`{C_{ice}}` are the specific heat capacities (J kg-1 K-1) of liquid water and ice, respectively (Table 1.4) and :math:`{c_{s,\,bedrock}} = 2 \times {10^6}` J m-3 K-1 is the heat capacity of bedrock. For snow



.. math::

 {c_i} = \frac{{{w_{ice,\,i}}}}{{\Delta {z_i}}}{C_{ice}} + \frac{{{w_{liq,\,i}}}}{{\Delta {z_i}}}{C_{liq}}

.

For the special case when snow is present ( :math:`{W_{sno}} > 0` ) but there are no explicit snow layers ( :math:`snl = 0` ), the heat capacity of the top layer is a blend of ice and soil heat capacity



.. math::

 {c_1} = c_1^ * + \frac{{{C_{ice}}{W_{sno}}}}{{\Delta {z_1}}}



where :math:`c_1^ *` is calculated from equation .

Hydrology
============

The hydrology for the pervious road generally follows that of CLM4 for bare soil surfaces and includes snow accumulation and melt, water transfer between snow layers, infiltration, evaporation, surface runoff, sub-surface drainage, redistribution within the soil column, and groundwater discharge and recharge to simulate changes in snow water :math:`\Delta {W_{sno}}` , soil water :math:`\Delta {w_{liq,\,i}}` , soil ice :math:`\Delta {w_{ice,\,i}}` , and water in the unconfined aquifer :math:`\Delta {W_a}` (all in kg m-2 or mm of H2O) (Figure 5.1). The water balance of the pervious road is

\[\Delta {W_{sno}} + \sum\limits_{i = 1}^{{N_{levsoi}}} {\left( {\Delta {w_{liq,\,i}} + \Delta {w_{ice,\,i}}} \right) + \Delta {W_a} = \left( \begin{gathered}

{q_{rain}} + {q_{sno}} - {E_{prvrd}} - {q_{over}} - {q_{drai}} \hfill \\
* {q_{rgwl}} - {q_{snwcp,\,ice}} \hfill \\

\end{gathered} \right)} \Delta t\]

where :math:`{q_{rain}}` is liquid part of precipitation, :math:`{q_{sno}}` is solid part of precipitation, :math:`{E_{prvrd}}` is the total evaporation (chapter 3), :math:`{q_{over}}` is surface runoff (section 5.2), :math:`{q_{drai}}` is sub-surface drainage (section 5.4), :math:`{q_{rgwl}}` and :math:`{q_{snwcp,\,ice}}` are liquid and solid runoff due to snow capping (section 5.5) (all in kg m-2 s-1), :math:`{N_{levsoi}}` is the number of soil layers, and :math:`\Delta t` is the time step (s). In general, snow capping will not be invoked for urban areas, but is described here for completeness.

Figure 5.1. Hydrologic processes simulated for the pervious road. Evaporation is supplied by all soil layers. An unconfined aquifer is added to the bottom of the soil column. The depth to the water table is :math:`{z_\nabla }` (m). Changes in aquifer water content :math:`{W_a}` (mm) are controlled by the balance between drainage from the aquifer water :math:`{q_{drai}}` and the aquifer recharge rate :math:`{q_{recharge}}` (kg m-2 s-1) (defined as positive from soil to aquifer).

The roof and the impervious road are hydrologically inactive except for their capacity to intercept, store, and evaporate a limited amount of liquid precipitation (1 kg m-2), and snow. Logistically, the storage of liquid precipitation is accounted for in the top layer :math:`i = 1` . The water in excess of this storage capacity is routed to surface runoff. These surfaces are also allowed to intercept solid precipitation (snow) and store this until the snowpack is melted or sublimated. No sub-surface drainage is allowed. The water balance of the roof and impervious road is



.. math::

 \Delta {W_{sno}} + \Delta {w_{liq,\,1}} + \Delta {w_{ice,\,1}} = \left( {{q_{rain}} + {q_{sno}} - {E_{roof}} - {q_{over}} - {q_{rgwl}} - {q_{snwcp,\,ice}}} \right)\Delta t





.. math::

 \Delta {W_{sno}} + \Delta {w_{liq,\,1}} + \Delta {w_{ice,\,1}} = \left( {{q_{rain}} + {q_{sno}} - {E_{imprvrd}} - {q_{over}} - {q_{rgwl}} - {q_{snwcp,\,ice}}} \right)\Delta t



where :math:`\Delta {w_{liq,\,1}}` and :math:`\Delta {w_{ice,\,1}}` are the liquid water and ice stored on the top of the urban surface. The sunlit and shaded walls are hydrologically inactive.

The rate of liquid and solid precipitation reaching the urban surface (kg m-2 s-1) is



.. math::

 {q_{grnd,liq}} = {q_{rain}}





.. math::

 {q_{grnd,ice}} = {q_{sno}}

.

Solid precipitation reaching the surface, :math:`{q_{grnd,\,ice}}\Delta t` , is added immediately to the snow pack (section 5.1). The liquid part, :math:`{q_{grnd,\,liq}}\Delta t` is added after surface fluxes, temperatures, soil water, and runoff have been determined.

Snow
--------

The parameterizations for snow are based primarily on Anderson (1976), (1991), and Dai and Zeng (1997). Snow can have up to five layers. These layers are indexed in the Fortran code as :math:`i = - 4, - 3, - 2, - 1,0` where layer :math:`i = 0` is the snow layer next to the top soil layer and layer :math:`i = - 4` is the top layer of a five-layer snow pack. Since the number of snow layers varies according to the snow depth, we use the notation :math:`snl + 1` to describe the top layer of snow for the variable layer snow pack, where :math:`snl` is the negative of the number of snow layers. Refer to Figure 5.2 for an example of the snow layer structure for a three layer snow pack.

Figure 5.2. Example of three layer snow pack ( :math:`snl = - 3` ). Shown are three snow layers, :math:`i = - 2` , :math:`i = - 1` , and :math:`i = 0` . The layer node depth is :math:`z` , the layer interface is :math:`{z_h}` , and the layer thickness is :math:`\Delta z` .

The state variables for snow are the mass of water :math:`{w_{liq,i}}` (kg m-2), mass of ice :math:`{w_{ice,i}}` (kg m-2), layer thickness :math:`\Delta {z_i}` (m), and temperature :math:`{T_i}` (chapter 4). The water vapor phase is neglected. Snow can also exist in the model without being represented by explicit snow layers. This occurs when the snowpack is less than a specified minimum snow depth ( :math:`{z_{sno}} < 0.01` m). In this case, the state variable is the mass of snow :math:`{W_{sno}}` (kg m-2).

The next two sections (5.1.1 and 5.1.2) describe the ice and water content of the snow pack assuming that at least one snow layer exists. See section 5.1.3 for a description of how a snow layer is initialized. Snow compaction is described in section 5.1.4 and snow layer combination and subdivision in section 5.1.5.

Ice Content
~~~~~~~~~~~~~~~~~

The conservation equation for mass of ice in snow layers is

\[\frac{{\partial {w_{ice,\,i}}}}{{\partial t}} = \left\{ \begin{gathered}

{q_{ice,\,i - 1}} - \frac{{{{\left( {\Delta {w_{ice,\,i}}} \right)}_p}}}{{\Delta t}} & i = snl + 1 \hfill \\
* \frac{{{{\left( {\Delta {w_{ice,\,i}}} \right)}_p}}}{{\Delta t}} & i = snl + 2, \ldots ,0 \hfill \\

\end{gathered} \right\}\]

where :math:`{q_{ice,\,i - 1}}` is the rate of ice accumulation from precipitation or frost or the rate of ice loss from sublimation (kg m-2 s-1) in the top layer and \[{{{{\left( {\Delta {w_{ice,\,i}}} \right)}_p}} \mathord{\left/

{\vphantom {{{{\left( {\Delta {w_{ice,\,i}}} \right)}_p}} {\Delta t}}} \right.

\kern-\nulldelimiterspace} {\Delta t}}\] is the change in ice due to phase change (melting rate) (section 4.2). The term :math:`{q_{ice,\,i - 1}}` is calculated in two steps as

 :math:`{q_{ice,\,i - 1}} = {q_{grnd,\,ice}} + \left( {{q_{frost}} - {q_{subl}}} \right)` 

where :math:`{q_{grnd,\,ice}}` is the rate of solid precipitation reaching the surface and :math:`{q_{frost}}` and :math:`{q_{subl}}` are gains due to frost and losses due to sublimation, respectively (section 3.2.4). In the first step, a new snow depth :math:`{z_{sno}}` (m) is calculated from



.. math::

 z_{sno}^{n + 1} = z_{sno}^n + \Delta {z_{sno}}



where



.. math::

 \Delta {z_{sno}} = \frac{{{q_{grnd,\,ice}}\Delta t}}{{{\rho _{sno}}}}



and :math:`{\rho _{sno}}` is the bulk density of newly fallen snow (kg m-3) ( 1976)

\[{\rho _{sno}} = \left\{ \begin{gathered}

50 + 1.7{\left( {17} \right)^{1.5}} & {T_{atm}} > {T_f} + 2 \hfill \\

50 + 1.7{\left( {{T_{atm}} - {T_f} + 15} \right)^{1.5}} & {T_f} - 15 < {T_{atm}} \leqslant {T_f} + 2 \hfill \\

50 & {T_{atm}} \leqslant {T_f} - 15 \hfill \\

\end{gathered} \right\}\]

where :math:`{T_{atm}}` is the atmospheric temperature (K), and :math:`{T_f}` is the freezing temperature of water (K) (Table 1.4). The mass of snow :math:`{W_{sno}}` is



.. math::

 W_{sno}^{n + 1} = W_{sno}^n + {q_{grnd,\,ice}}\Delta t

.

The ice content of the top layer and the layer thickness are updated as



.. math::

 w_{ice,\,snl + 1}^{n + 1} = w_{ice,\,snl + 1}^n + {q_{grnd,\,ice}}\Delta t





.. math::

 \Delta z_{snl + 1}^{n + 1} = \Delta z_{snl + 1}^n + \Delta {z_{sno}}

.

In the second step, after surface fluxes and temperatures have been determined (chapters 3 and 4), :math:`{w_{ice,\,snl + 1}}` is updated for frost or sublimation as



.. math::

 w_{ice,\,snl + 1}^{n + 1} = w_{ice,\,snl + 1}^n + \left( {{q_{frost}} - {q_{subl}}} \right)\Delta t

.

If :math:`w_{ice,\,snl + 1}^{n + 1} < 0` upon solution of equation , the ice content is reset to zero and the liquid water content :math:`{w_{liq,\,snl + 1}}` is reduced by the amount required to bring :math:`w_{ice,\,snl + 1}^{n + 1}` up to zero. The snow water equivalent :math:`{W_{sno}}` is capped to not exceed 1000 kg m-2. If the addition of :math:`{q_{frost}}` were to result in :math:`{W_{sno}} > 1000` kg m-2, the frost term :math:`{q_{frost}}` is instead added to the ice runoff term :math:`{q_{snwcp,\,ice}}` (section 5.5).

Water Content
~~~~~~~~~~~~~~~~~~~

The conservation equation for mass of water in snow layers is



.. math::

 \frac{{\partial {w_{liq,\,i}}}}{{\partial t}} = \left( {{q_{liq,\,i - 1}} - {q_{liq,\,i}}} \right) + \frac{{{{\left( {\Delta {w_{liq,\,i}}} \right)}_p}}}{{\Delta t}}



where :math:`{q_{liq,\,i - 1}}` is the flow of liquid water into layer :math:`i` from the layer above, :math:`{q_{liq,\,i}}` is the flow of water out of layer :math:`i` to the layer below, \[{{{{\left( {\Delta {w_{liq,\,i}}} \right)}_p}} \mathord{\left/

{\vphantom {{{{\left( {\Delta {w_{liq,\,i}}} \right)}_p}} {\Delta t}}} \right.

\kern-\nulldelimiterspace} {\Delta t}}\] is the change in liquid water due to phase change (melting rate) (section 4.2). For the top snow layer only,



.. math::

 {q_{liq,\,i - 1}} = {q_{grnd,\,liq}} + \left( {{q_{sdew}} - {q_{seva}}} \right)



where :math:`{q_{grnd,\,liq}}` is the rate of liquid precipitation reaching the snow, :math:`{q_{seva}}` is the evaporation of liquid water and :math:`{q_{sdew}}` is the liquid dew (section 3.2.4). After surface fluxes and temperatures have been determined (chapters 3 and 4), :math:`{w_{liq,\,snl + 1}}` is updated for the liquid precipitation reaching the ground and dew or evaporation as



.. math::

 w_{liq,\,snl + 1}^{n + 1} = w_{liq,\,snl + 1}^n + \left( {{q_{grnd,\,liq}} + {q_{sdew}} - {q_{seva}}} \right)\Delta t

.

When the liquid water within a snow layer exceeds the layer’s holding capacity, the excess water is added to the underlying layer, limited by the effective porosity ( :math:`1 - {\theta _{ice}}` ) of the layer. The flow of water is assumed to be zero ( :math:`{q_{liq,\,i}} = 0` ) if the effective porosity of either of the two layers (

.. math::

 1 - {\theta _{ice,\,i}}{\text{ and }}1 - {\theta _{ice,\,i + 1}}

) is less than :math:`{\theta _{imp}} = 0.05` , the water impermeable volumetric water content. Thus, water flow between layers, :math:`{q_{liq,\,i}}` , for 

.. math::

 i = snl + 1, \ldots ,0

 is initially calculated as



.. math::

 {q_{liq,\,i}} = \frac{{{\rho _{liq}}\left[ {{\theta _{liq,\,i}} - {S_r}\left( {1 - {\theta _{ice,\,i}}} \right)} \right]\Delta {z_i}}}{{\Delta t}} \geqslant 0



where the volumetric liquid water :math:`{\theta _{liq,\,i}}` and ice :math:`{\theta _{ice,\,i}}` contents are



.. math::

 {\theta _{ice,\,i}} = \frac{{{w_{ice,\,i}}}}{{\Delta {z_i}{\rho _{ice}}}} \leqslant 1





.. math::

 {\theta _{liq,\,i}} = \frac{{{w_{liq,\,i}}}}{{\Delta {z_i}{\rho _{liq}}}} \leqslant 1 - {\theta _{ice,\,i}}

,

and :math:`{S_r} = 0.033` is the irreducible water saturation (snow holds a certain amount of liquid water due to capillary retention after drainage has ceased ( 1976)). The water holding capacity of the underlying layer limits the flow of water :math:`{q_{liq,\,i}}` calculated in equation , unless the underlying layer is the surface layer, as



.. math::

 {q_{liq,\,i}} \leqslant \frac{{{\rho _{liq}}\left[ {1 - {\theta _{ice,\,i + 1}} - {\theta _{liq,\,i + 1}}} \right]\Delta {z_{i + 1}}}}{{\Delta t}} & i = snl + 1, \ldots , - 1

.

The volumetric liquid water content :math:`{\theta _{liq,\,i}}` is updated as



.. math::

 \theta _{liq,\,i}^{n + 1} = \theta _{liq,\,i}^n + \left( {{q_{i - 1}} - {q_i}} \right)\Delta t

.

Equations - are solved sequentially from top (

.. math::

 i = snl + 1

) to bottom (

.. math::

 i = 0

) snow layer in each time step. The total flow of liquid water reaching the urban surface is then :math:`{q_{liq,\,0}}` .

Initialization of snow layer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If there are no existing snow layers ( :math:`snl + 1 = 1` ) but :math:`{z_{sno}} \geqslant 0.01` m after accounting for solid precipitation :math:`{q_{sno}}` , then a snow layer is initialized ( :math:`snl = - 1` ) as follows

\[\begin{gathered}

\Delta {z_0} = {z_{sno}} \hfill \\

{z_o} = - 0.5\Delta {z_0} \hfill \\

{z_{h,\, - 1}} = - \Delta {z_0} \hfill \\

{T_0} = \min \left( {{T_f},{T_{atm}}} \right) \hfill \\

{w_{ice,\,0}} = {W_{sno}} \hfill \\

{w_{liq,\,0}} = 0 \hfill \\

\end{gathered} \].

Snow Compaction
~~~~~~~~~~~~~~~~~~~~~

Snow compaction is initiated after the hydrology calculations [surface runoff (section 5.2), infiltration (section 5.2), soil water (section 5.3), groundwater-soilwater interactions (section 5.4)] are complete. Compaction of snow includes three types of processes: destructive metamorphism of new snow (crystal breakdown due to wind or thermodynamic stress); snow load or overburden (pressure); and melting (changes in snow structure due to melt-freeze cycles plus changes in crystals due to liquid water). The total fractional compaction rate for each snow layer :math:`{C_{R,\,i}}` (s-1) is the sum of the three compaction processes



.. math::

 {C_{R,\,i}} = \frac{1}{{\Delta {z_i}}}\frac{{\partial \Delta {z_i}}}{{\partial t}} = {C_{R1,\,i}} + {C_{R2,\,i}} + {C_{R3,\,i}}

.

Compaction is not allowed if the layer is saturated



.. math::

 1 - \left( {\frac{{{w_{ice,\,i}}}}{{\Delta {z_i}{\rho _{ice}}}} + \frac{{{w_{liq,\,i}}}}{{\Delta {z_i}{\rho _{liq}}}}} \right) \leqslant 0.001



or if the ice content is below a minimum value ( :math:`{w_{ice,\,i}} \leqslant 0.1` ).

Compaction as a result of destructive metamorphism :math:`{C_{R1,\;i}}` (s-1) is temperature dependent ( 1976)



.. math::

 {C_{R1,\,i}} = {\left[ {\frac{1}{{\Delta {z_i}}}\frac{{\partial \Delta {z_i}}}{{\partial t}}} \right]_{metamorphism}} = - {c_3}{c_1}{c_2}\exp \left[ { - {c_4}\left( {{T_f} - {T_i}} \right)} \right]



where :math:`{c_3} = 2.777 \times {10^{ - 6}}` (s-1) is the fractional compaction rate for :math:`{T_i} = {T_f}` , :math:`{c_4} = 0.04` K-1, and

\[\begin{gathered}

{c_1} = 1 & \frac{{{w_{ice,\,i}}}}{{\Delta {z_i}}} \leqslant 100{\text{ kg }}{{\text{m}}^{{\text{ - 3}}}} \hfill \\

{c_1} = \exp \left[ { - 0.046\left( {\frac{{{w_{ice,\,i}}}}{{\Delta {z_i}}} - 100} \right)} \right] & \frac{{{w_{ice,\,i}}}}{{\Delta {z_i}}} > 100{\text{ kg }}{{\text{m}}^{{\text{ - 3}}}} \hfill \\

{c_2} = 2 & \frac{{{w_{liq,\,i}}}}{{\Delta {z_i}}} > 0.01 \hfill \\

{c_2} = 1 & \frac{{{w_{liq,\,i}}}}{{\Delta {z_i}}} \leqslant 0.01 \hfill \\

\end{gathered} \]

where ${{{w_{ice,\,i}}} \mathord{\left/

{\vphantom {{{w_{ice,\,i}}} {\Delta {z_i}}}} \right.

\kern-\nulldelimiterspace} {\Delta {z_i}}} :math:`and` {{{w_{liq,\,i}}} \mathord{\left/

{\vphantom {{{w_{liq,\,i}}} {\Delta {z_i}}}} \right.

\kern-\nulldelimiterspace} {\Delta {z_i}}}$ are the bulk densities of liquid water and ice (kg m-3).

The compaction rate as a result of overburden :math:`{C_{R2,\;i}}` (s-1) is a linear function of the snow load pressure :math:`{P_{s,\,i}}` (kg m-2) ( 1976)



.. math::

 {C_{R2,\,i}} = {\left[ {\frac{1}{{\Delta {z_i}}}\frac{{\partial \Delta {z_i}}}{{\partial t}}} \right]_{overburden}} = - \frac{{{P_{s,\,i}}}}{\eta }



where :math:`\eta` is a viscosity coefficient (kg s m-2) that varies with density and temperature as



.. math::

 \eta = {\eta _0}\exp \left[ {{c_5}\left( {{T_f} - {T_i}} \right) + {c_6}\frac{{{w_{ice,\,i}}}}{{\Delta {z_i}}}} \right]



where :math:`{\eta _0} = 9 \times {10^5}` kg s m-2, and :math:`{c_5} = 0.08` K-1, :math:`{c_6} = 0.023` m3 kg-1 are constants. The snow load pressure :math:`{P_{s,\,i}}` is calculated for each layer as the sum of the ice :math:`{w_{ice,\,i}}` and liquid water contents :math:`{w_{liq,\,i}}` of the layers above plus half the ice and liquid water contents of the layer being compacted



.. math::

 {P_{s,\,i}} = \frac{{\left( {{w_{ice,\,i}} + {w_{liq,\,i}}} \right)}}{2} + \sum\limits_{j = snl + 1}^{j = i - 1} {\left( {{w_{ice,\,j}} + {w_{liq,\,j}}} \right)}

.

The compaction rate due to melting :math:`{C_{R3,\;i}}` (s-1) is taken to be the ratio of the change in snow ice fraction after the melting to the fraction before melting



.. math::

 {C_{R3,\,i}} = {\left[ {\frac{1}{{\Delta {z_i}}}\frac{{\partial \Delta {z_i}}}{{\partial t}}} \right]_{melt}} = - \frac{1}{{\Delta t}}\max \left( {0,\frac{{f_{ice,\,i}^n - f_{ice,\,i}^{n + 1}}}{{f_{ice,\,i}^n}}} \right)



where the fraction of ice :math:`{f_{ice,\,i}}` is



.. math::

 {f_{ice,\,i}} = \frac{{{w_{ice,\,i}}}}{{{w_{ice,\,i}} + {w_{liq,\,i}}}}



and melting is identified during the phase change calculations (section 4.2).

The snow layer thickness after compaction is then



.. math::

 \Delta z_i^{n + 1} = \Delta z_i^n\left( {1 + {C_{R,\,i}}\Delta t} \right)

.

Snow Layer Combination and Subdivision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After the determination of snow temperature including phase change (chapter 4), snow hydrology (sections 5.2.1, 5.2.2, and 5.2.3), and the compaction calculations (5.2.4), the number of snow layers is adjusted by either combining or subdividing layers. The combination and subdivision of snow layers is based on Jordan (1991).
#. 1.5.1 Combination

If a snow layer has nearly melted or if its thickness :math:`\Delta {z_i}` is less than the prescribed minimum thickness :math:`\Delta {z_{\min }}` (Table 5.1), the layer is combined with a neighboring layer. The overlying or underlying layer is selected as the neighboring layer according to the following rules

If the top layer is being removed, it is combined with the underlying layer

If the underlying layer is not snow, the layer is combined with the overlying layer

If the layer is nearly completely melted, the layer is combined with the underlying layer

If none of the above rules apply, the layer is combined with the thinnest neighboring layer.

A first pass is made through all snow layers to determine if any layer is nearly melted ( :math:`{w_{ice,\,i}} \leqslant 0.1` ). If so, the remaining liquid water and ice content of layer :math:`i` is combined with the underlying neighbor :math:`i + 1` as



.. math::

 {w_{liq,\,i + 1}} = {w_{liq,\,i + 1}} + {w_{liq,\,i}}





.. math::

 {w_{ice,\,i + 1}} = {w_{ice,\,i + 1}} + {w_{ice,\,i}}

.

This includes the snow layer directly above the urban surface. In this case, the liquid water and ice content of the melted snow layer is added as ponded water/ice on the urban surface layer. The layer properties, :math:`{T_i}` , :math:`{w_{ice,\,i}}` , :math:`{w_{liq,\,i}}` , :math:`\Delta {z_i}` , are then re-indexed so that the layers above the eliminated layer are shifted down by one and the number of snow layers is decremented accordingly.

At this point, if there are no explicit snow layers remaining ( :math:`snl = 0` ), the snow water equivalent :math:`{W_{sno}}` and snow depth :math:`{z_{sno}}` are set to zero, otherwise, :math:`{W_{sno}}` and :math:`{z_{sno}}` are re-calculated as



.. math::

 {W_{sno}} = \sum\limits_{i = snl + 1}^{i = 0} {\left( {{w_{ice,\,i}} + {w_{liq,\,i}}} \right)}





.. math::

 {z_{sno}} = \sum\limits_{i = snl + 1}^{i = 0} {\Delta {z_i}}

.

If the snow depth :math:`0 < {z_{sno}} < 0.01` m, the number of snow layers is set to zero, the total ice content of the snowpack :math:`\sum\limits_{i = snl + 1}^{i = 0} {{w_{ice,\;i}}}` is assigned to :math:`{W_{sno}}` , and the total liquid water :math:`\sum\limits_{i = snl + 1}^{i = 0} {{w_{liq,\;i}}}` is assigned to the urban surface layer. Otherwise, the layers are combined according to the rules above.

When two snow layers are combined (denoted here as :math:`1` and :math:`2` ), their thickness combination ( :math:`c` ) is



.. math::

 \Delta {z_c} = \Delta {z_1} + \Delta {z_2}

,

their mass combination is



.. math::

 {w_{liq,\,c}} = {w_{liq,\,1}} + {w_{liq,\,2}}





.. math::

 {w_{ice,\,c}} = {w_{ice,\,1}} + {w_{ice,\,2}}

,

and their temperatures are combined as



.. math::

 {T_c} = {T_f} + \frac{{{h_c} - {L_f}{w_{liq,\,c}}}}{{{C_{ice}}{w_{ice,\,c}} + {C_{liq}}{w_{liq,\,c}}}}



where :math:`{h_c} = {h_1} + {h_2}` is the combined enthalpy :math:`{h_i}` of the two layers where



.. math::

 {h_i} = \left( {{C_{ice}}{w_{ice,\,i}} + {C_{liq}}{w_{liq,\,i}}} \right)\left( {{T_i} - {T_f}} \right) + {L_f}{w_{liq,\,i}}

.

In these equations, :math:`{L_f}` is the latent heat of fusion (J kg-1) and :math:`{C_{liq}}` and :math:`{C_{ice}}` are the specific heat capacities (J kg-1 K-1) of liquid water and ice, respectively (Table 1.4). After layer combination, the node depths and layer interfaces (Figure 5.2) are recalculated from



.. math::

 {z_i} = {z_{h,\,i}} - 0.5\Delta {z_i} & i = 0, \ldots ,snl + 1





.. math::

 {z_{h,\,i - 1}} = {z_{h,\,i}} - \Delta {z_i} & i = 0, \ldots ,snl + 1



where :math:`\Delta {z_i}` is the layer thickness.

Table 5.1. Minimum and maximum thickness of snow layers (m)

The maximum snow layer thickness, :math:`\Delta {z_{\max }}` , depends on the number of layers, :math:`{N_l}` and :math:`{N_u}` .
#. 1.5.2 Subdivision

The snow layers are subdivided when the layer thickness exceeds a prescribed maximum thickness :math:`\Delta {z_{\max }}` with lower and upper bounds that depend on the number of snow layers (Table 5.1). For example, if there is only one layer, then the maximum thickness of that layer is 0.03 m, however, if there is more than one layer, then the maximum thickness of the top layer is 0.02 m. Layers are checked sequentially from top to bottom for this limit. If there is only one snow layer and its thickness is greater than 0.03 m (Table 5.1), the layer is subdivided into two layers of equal thickness, liquid water and ice contents, and temperature. If there is an existing layer below the layer to be subdivided, the thickness :math:`\Delta {z_i}` , liquid water and ice contents, :math:`{w_{liq,\;i}}` and :math:`{w_{ice,\;i}}` , and temperature :math:`{T_i}` of the excess snow are combined with the underlying layer according to equations -. If there is no underlying layer after adjusting the layer for the excess snow, the layer is subdivided into two layers of equal thickness, liquid water and ice contents. The vertical snow temperature profile is maintained by calculating the slope between the layer above the splitting layer ( :math:`{T_1}` ) and the splitting layer ( :math:`{T_2}` ) and constraining the new temperatures ( :math:`T_2^{n + 1}` , :math:`T_3^{n + 1}` ) to lie along this slope. The temperature of the lower layer is first evaluated from

\[{T'_3} = T_2^n - \left( {\frac{{T_1^n - T_2^n}}{{{{\left( {\Delta z_1^n + \Delta z_2^n} \right)} \mathord{\left/

{\vphantom {{\left( {\Delta z_1^n + \Delta z_2^n} \right)} 2}} \right.

\kern-\nulldelimiterspace} 2}}}} \right)\left( {\frac{{\Delta z_2^{n + 1}}}{2}} \right)\],

then adjusted as,

\[\begin{gathered}

T_3^{n + 1} = T_2^n & {{T'}_3} \geqslant {T_f} \hfill \\

T_2^{n + 1} = T_2^n + \left( {\frac{{T_1^n - T_2^n}}{{{{\left( {\Delta {z_1} + \Delta z_2^n} \right)} \mathord{\left/

{\vphantom {{\left( {\Delta {z_1} + \Delta z_2^n} \right)} 2}} \right.

\kern-\nulldelimiterspace} 2}}}} \right)\left( {\frac{{\Delta z_2^{n + 1}}}{2}} \right) & {{T'}_3} < {T_f} \hfill \\

\end{gathered} \]

where here the subscripts 1, 2, and 3 denote three layers numbered from top to bottom. After layer subdivision, the node depths and layer interfaces are recalculated from equations and .

Surface Runoff and Infiltration
-----------------------------------

For the roof and impervious road, water on these surfaces in excess of a maximum ponding limit :math:`{w_{pond,\max }} = 1` (kg m-2) is routed to surface runoff as

\[\begin{gathered}

{q_{over}} = \frac{{{w_{liq,\,1}}}}{{\Delta t}} + {q_{liq,\,0}} - {q_{seva}} - \frac{{{w_{pond,\max }}}}{{\Delta t}} \geqslant 0 & snl = 0 \hfill \\

{q_{over}} = {q_{liq,\,0}} & snl < 0 \hfill \\

\end{gathered} \].

where :math:`{q_{liq,\,0}}` is the rate of liquid water reaching the surface from rain (section 5.1) and/or snowmelt (section 5.1.2) and :math:`{q_{seva}}` is the evaporation of liquid water from the top layer (section 3.4). The liquid water content of the top layer is adjusted to

\[\begin{gathered}

{w_{liq,\,1}} = {w_{pond,\max }} & {q_{over}} > 0 \hfill \\

{w_{liq,\,1}} = {w_{liq,\,1}} + \left( {{q_{liq,\,0}} - {q_{seva}}} \right)\Delta t \geqslant 0 & {q_{over}} = 0 \hfill \\

\end{gathered} \].

For the pervious road, the simple TOPMODEL-based (Beven and Kirkby 1979) runoff model (SIMTOP) described by Niu et al. (2005) is implemented. A key concept underlying this approach is that of fractional saturated/impermeable area :math:`{f_{sat}}` , which is determined by the topographic characteristics and soil moisture state of a grid cell. The surface runoff consists of overland flow due to saturation excess (Dunne runoff) and infiltration excess (Hortonian runoff) mechanisms



.. math::

 {q_{over}} = {f_{sat}}{q_{liq,\,0}} + \left( {1 - {f_{sat}}} \right)\max \left( {0,\,{q_{liq,\,0}} - {q_{infl,\,\max }}} \right)



where :math:`{q_{liq,\,0}}` is liquid precipitation reaching the ground plus any melt water from snow (kg m-2 s-1) and :math:`{q_{infl,\,\max }}` is a maximum soil infiltration capacity (kg m-2 s-1). In Niu et al. (2005), :math:`{f_{sat}}` was a function of soil moisture whose potential or maximum value, :math:`{f_{\max }}` , was solely determined by topographic characteristics. Niu and Yang (2006) modified the expression for :math:`{f_{sat}}` to include a dependence on impermeable area fraction in frozen soil, :math:`{f_{frz,\,1}}` , of the top :math:`i = 1` soil layer as



.. math::

 {f_{sat}} = \left( {1 - {f_{frz,\,1}}} \right){f_{\max }}\exp \left( { - 0.5{f_{over}}{z_\nabla }} \right) + {f_{frz,\,1}}



where :math:`{f_{\max }}` is the maximum saturated fraction, :math:`{f_{over}}` is a decay factor (m-1), and :math:`{z_\nabla }` is the water table depth (m) (section 5.4). The maximum saturated fraction, :math:`{f_{\max }}` , is defined as the discrete cumulative distribution function (CDF) of the topographic index when the grid cell mean water table depth is zero. Thus, :math:`{f_{\max }}` is the percent of pixels in a grid cell whose topographic index is larger than or equal to the grid cell mean topographic index. It should be calculated explicitly from the CDF at each grid cell at the resolution that the model is run. However, because this is a computationally intensive task for global applications, :math:`{f_{\max }}` is calculated once from the CDF at a spatial resolution of 0.5° by 0.5° following Niu et al. (2005) and then area-averaged to the desired resolution. The 0.5° resolution is compatible with the resolution of other CLM input surface datasets (e.g., plant functional types, leaf area index). The decay factor :math:`{f_{over}}` for global simulations was determined through sensitivity analysis and comparison with observed runoff to be 0.5 m-1.

The impermeable fraction :math:`{f_{frz,\,i}}` is parameterized as a function of soil ice content (Niu and Yang 2006)



.. math::

 {f_{frz,\,i}} = \frac{{\exp \left[ { - \alpha \left( {1 - \frac{{{w_{ice,\,i}}}}{{{w_{ice,\,i}} + {w_{liq,\,i}}}}} \right)} \right] - \exp \left( { - \alpha } \right)}}{{1 - \exp \left( { - \alpha } \right)}}



where :math:`\alpha = 3` is an adjustable scale-dependent parameter, and :math:`{w_{ice,\,i}}` and :math:`{w_{liq,\,i}}` are the ice and liquid water contents of soil layer :math:`i` (kg m-2).

The maximum infiltration capacity :math:`{q_{infl,\,\max }}` in equation is determined from soil texture and soil moisture (Entekhabi and Eagleson 1989) as



.. math::

 {q_{infl,\,\max }} = {k_{sat,\,1}}\left[ {1 + v\left( {s - 1} \right)} \right]

.

The liquid water content of the top soil layer relative to effective porosity and adjusted for saturated fraction is determined from

\[\begin{gathered}

s = \frac{{\frac{{{\theta _{liq,\,1}}}}{{\max \left( {{\theta _{imp}},\,{\theta _{sat,\,1}} - {\theta _{ice,\,1}}} \right)}} - {f_{sat}}}}{{1 - {f_{sat}}}} \geqslant 0 & \frac{{{\theta _{liq,\,1}}}}{{\max \left( {{\theta _{imp}},\,{\theta _{sat,\,1}} - {\theta _{ice,\,1}}} \right)}} \geqslant 0.01 \hfill \\

& \,1 - {f_{sat}} \geqslant 0.01 \hfill \\

\end{gathered} \]

where :math:`{\theta _{liq,\,1}}` and :math:`{\theta _{ice,\,1}}` are the volumetric liquid water and ice contents of the top soil layer, and :math:`{\theta _{imp}} = 0.05` is a minimum effective porosity. The variable :math:`v` is



.. math::

 v = - {\left( {\frac{{d\psi }}{{ds}}} \right)_{s = 1}}\frac{1}{{0.5\Delta {z_1}}}



where :math:`\Delta {z_1}` is the thickness of the top soil layer (mm) and



.. math::

 {\left( {\frac{{d\psi }}{{ds}}} \right)_{s = 1}} = - {B_1}{\psi _{sat,\,1}}

.

The saturated hydraulic conductivity :math:`{k_{sat,\,1}}` (kg m-2 s-1), volumetric water content at saturation (i.e., porosity) :math:`{\theta _{sat,\,1}}` , Clapp and Hornberger (1978) exponent :math:`{B_1}` , and saturated soil matric potential :math:`{\psi _{sat,\,1}}` (mm) are determined from soil texture (section 5.3.1).

Infiltration into the surface soil layer of the pervious road is defined as the residual of the surface water balance



.. math::

 {q_{infl}} = {q_{liq,\,0}} - {q_{over}} - {q_{seva}}



when no snow layers exist, and



.. math::

 {q_{infl}} = {q_{liq,\,0}} - {q_{over}}



when at least one snow layer is present.

The infiltration for urban surfaces other than pervious road is



.. math::

 {q_{infl}} = 0

.

Soil Water for the Pervious Road
------------------------------------

Soil water for the pervious road is predicted from a multi-layer model, in which the vertical soil moisture transport is governed by infiltration, surface and sub-surface runoff, gradient diffusion, gravity, evapotranspiration through root extraction, and interactions with groundwater (Figure 5.1). Vegetation is not represented explicitly, however, the total evaporation calculated in section 3.2.4, if not assigned to surface evaporation, is removed from each soil layer through an evapotranspiration loss ( :math:`s` in the equation below). The following derivation generally follows that of Z.-L. Yang (1998, unpublished manuscript) with modifications by Zeng and Decker (2009).

For one-dimensional vertical water flow in soils, the conservation of mass is stated as



.. math::

 \frac{{\partial \theta }}{{\partial t}} = - \frac{{\partial q}}{{\partial z}} - Q



where :math:`\theta` is the volumetric soil water content (mm3 of water mm-3 of soil), :math:`t` is time (s), :math:`z` is height above some datum in the soil column (mm) (positive upwards), :math:`q` is soil water flux (kg m-2 s-1 or mm s-1) (positive upwards), and :math:`Q` is a soil moisture sink term (mm of water mm-1 of soil s-1) (ET loss). This equation is solved numerically by dividing the soil column into multiple layers in the vertical and integrating downward over each layer with an upper boundary condition of the infiltration flux into the top soil layer :math:`{q_{infl}}` and a lower boundary condition specified as zero flux.

The soil water flux :math:`q` in equation can be described by Darcy’s law

 :math:`q = - k\frac{{\partial {\psi _h}}}{{\partial z}}` 

where :math:`k` is the hydraulic conductivity (mm s-1), and :math:`{\psi _h}` is the hydraulic potential (mm). The hydraulic potential is

 :math:`{\psi _h} = {\psi _m} + {\psi _z}` 

where :math:`{\psi _m}` is the soil matric potential (mm) (which is related to the adsorptive and capillary forces within the soil matrix), and :math:`{\psi _z}` is the gravitational potential (mm) (the vertical distance from an arbitrary reference elevation to a point in the soil). If the reference elevation is the soil surface, then :math:`{\psi _z} = z` . Letting :math:`\psi = {\psi _m}` , Darcy’s law becomes



.. math::

 q = - k\left[ {\frac{{\partial \left( {\psi + z} \right)}}{{\partial z}}} \right]

.

Darcy’s equation can be further manipulated to yield



.. math::

 q = - k\left[ {\frac{{\partial \left( {\psi + z} \right)}}{{\partial z}}} \right] = - k\left( {\frac{{\partial \psi }}{{\partial z}} + 1} \right) = - k\left( {\frac{{\partial \theta }}{{\partial z}}\frac{{\partial \psi }}{{\partial \theta }} + 1} \right)

.

Substitution of this equation into equation with :math:`Q = 0` , yields the Richards equation



.. math::

 \frac{{\partial \theta }}{{\partial t}} = \frac{\partial }{{\partial z}}\left[ {k\left( {\frac{{\partial \theta }}{{\partial z}}\frac{{\partial \psi }}{{\partial \theta }}} \right) + 1} \right]

.

Zeng and Decker (2009) note that this :math:`\theta` -based form of the Richards equation cannot maintain the hydrostatic equilibrium soil moisture distribution because of the truncation errors of the finite-difference numerical scheme. They show that this deficiency can be overcome by subtracting the equilibrium state from equation as



.. math::

 q = - k\left[ {\frac{{\partial \left( {\psi + z - C} \right)}}{{\partial z}}} \right]



where :math:`C` is a constant hydraulic potential above the water table :math:`{z_\nabla }` 



.. math::

 C = {\psi _E} + z = {\psi _{sat}}{\left[ {\frac{{{\theta _E}\left( z \right)}}{{{\theta _{sat}}}}} \right]^{ - B}} + z = {\psi _{sat}} + {z_\nabla }



so that



.. math::

 q = - k\left[ {\frac{{\partial \left( {\psi - {\psi _E}} \right)}}{{\partial z}}} \right]

.

where :math:`{\psi _E}` is the equilibrium soil matric potential (mm). Substitution of equations and into equation yields Zeng and Decker’s (2009) modified Richards equation



.. math::

 \frac{{\partial \theta }}{{\partial t}} = \frac{\partial }{{\partial z}}\left[ {k\left( {\frac{{\partial \left( {\psi - {\psi _E}} \right)}}{{\partial z}}} \right)} \right] - Q



where the soil moisture source/sink term :math:`Q` is now included.

Hydraulic Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~

The hydraulic conductivity :math:`{k_i}` (mm s-1) and the soil matric potential :math:`{\psi _i}` (mm) for layer :math:`i` vary with volumetric soil water :math:`{\theta _i}` and soil texture (

.. math::

 \% san{d_i}

 and 

.. math::

 \% cla{y_i}

, section 1.2.2) based on the work of Clapp and Hornberger (1978) and Cosby et al. (1984). In CLM4, organic matter modifies soil properties according to Lawrence and Slater (2008). Urban soils are assumed to have no organic matter so the equations below are shown in their reduced form.

The hydraulic conductivity is defined at the depth of the interface of two adjacent layers :math:`{z_{h,\,i}}` (Figure 5.3) and is a function of the saturated hydraulic conductivity 

.. math::

 {k_{sat}}\left[ {{z_{h,\,i}}} \right]

, the total (ice plus liquid) volumetric soil moisture of the two layers :math:`{\theta _i}` and :math:`{\theta _{i + 1}}` and the impermeable fraction :math:`{f_{frz,\,i}}` 

\[k\left[ {{z_{h,\,i}}} \right] = \left\{ \begin{gathered}

\left( {1 - \frac{{{f_{frz,\,i}} + {f_{frz,\,i + 1}}}}{2}} \right){k_{sat}}\left[ {{z_{h,\,i}}} \right]{\left[ {\frac{{0.5\left( {{\theta _{\,i}} + {\theta _{\,i + 1}}} \right)}}{{0.5\left( {{\theta _{sat,\,i}} + {\theta _{sat,\,i + 1}}} \right)}}} \right]^{2{B_i} + 3}} & 1 \leqslant i \leqslant {N_{levsoi}} - 1 \hfill \\

\left( {1 - {f_{frz,\,i}}} \right){k_{sat}}\left[ {{z_{h,\,i}}} \right]{\left( {\frac{{{\theta _{\,i}}}}{{{\theta _{sat,\,i}}}}} \right)^{2{B_i} + 3}} & i = {N_{levsoi}} \hfill \\

\end{gathered} \right\}\]

where :math:`{f_{frz,\,i}}` is defined in equation . The saturated hydraulic conductivity 

.. math::

 {k_{sat}}\left[ {{z_{h,\,i}}} \right]

 (mm s-1) depends on soil texture (Cosby et al. 1984) as



.. math::

 {k_{sat}}\left[ {{z_{h,\,i}}} \right] = 0.0070556 \times {10^{ - 0.884 + 0.0153{{\left( {\% sand} \right)}_i}}}

.

The water content at saturation (i.e., porosity) is



.. math::

 {\theta _{sat,\,i}} = 0.489 - 0.00126{\left( {\% sand} \right)_i}



and the exponent “ :math:`B` ” is



.. math::

 {B_i} = 2.91 + 0.159{\left( {\% clay} \right)_i}

.

The soil matric potential (mm) is defined at the node depth :math:`{z_i}` of each layer :math:`i` (Figure 5.3)



.. math::

 {\psi _i} = {\psi _{sat,\,i}}{\left( {\frac{{{\theta _{\,i}}}}{{{\theta _{sat,\,i}}}}} \right)^{ - {B_i}}} \geqslant - 1 \times {10^8} & 0.01 \leqslant \frac{{{\theta _i}}}{{{\theta _{sat,\,i}}}} \leqslant 1



where the saturated soil matric potential (mm) is



.. math::

 {\psi _{sat,\,i}} = - 10.0 \times {10^{1.88 - 0.0131{{\left( {\% sand} \right)}_i}}}

.

Numerical Solution
~~~~~~~~~~~~~~~~~~~~~~~~

With reference to Figure 5.3, the equation for conservation of mass (equation can be integrated over each layer as



.. math::

 \int\limits_{ - {z_{h,\,i}}}^{ - {z_{h,\,i - 1}}} {\frac{{\partial \theta }}{{\partial t}}\,} dz = - \int\limits_{ - {z_{h,\,i}}}^{ - {z_{h,\,i - 1}}} {\frac{{\partial q}}{{\partial z}}} \,dz - \int\limits_{ - {z_{h,\,i}}}^{ - {z_{h,\,i - 1}}} {Q\,dz}

.

Note that the integration limits are negative since :math:`z` is defined as positive upward from the soil surface. This equation can be written as



.. math::

 \Delta {z_i}\frac{{\partial {\theta _{liq,\,i}}}}{{\partial t}} = - {q_{i - 1}} + {q_i} - {e_i}



where :math:`{q_i}` is the flux of water across interface :math:`{z_{h,\,i}}` , :math:`{q_{i - 1}}` is the flux of water across interface :math:`{z_{h,\,i - 1}}` , and :math:`{e_i}` is a layer-averaged soil moisture sink term (ET loss) defined as positive for flow out of the layer (mm s-1). Taking the finite difference with time and evaluating the fluxes implicitly at time :math:`n + 1` yields



.. math::

 \frac{{\Delta {z_i}\Delta {\theta _{liq,\,i}}}}{{\Delta t}} = - q_{i - 1}^{n + 1} + q_i^{n + 1} - {e_i}



where :math:`\Delta {\theta _{liq,\,i}} = \theta _{liq,\,i}^{n + 1} - \theta _{liq,\,i}^n` is the change in volumetric soil liquid water of layer :math:`i` in time :math:`\Delta t` and :math:`\Delta {z_i}` is the thickness of layer :math:`i` (mm).

The water removed by evapotranspiration in each layer :math:`{e_i}` is a function of the total evapotranspiration :math:`E_{prvrd}^{et}` (section 3.2.4) and the effective root fraction :math:`{r_{e,\,i}}` 



.. math::

 {e_i} = {r_{e,\,i}}E_{prvrd}^{et}

.

The effective root fraction :math:`{r_{e,\,i}}` is

\[{r_{e,\,i}} = \left\{ \begin{gathered}

\frac{{{r_i}{w_i}}}{{{\alpha _{soi}}}} & {\alpha _{soi}} > 0 \hfill \\

0 & {\alpha _{soi}} = 0 \hfill \\

\end{gathered} \right\}\]

where :math:`{r_i}` is the fraction of roots in layer :math:`i` (equation ), :math:`{w_i}` is a soil wetness factor for layer :math:`i` (equation ), and :math:`{\alpha _{soi}}` is a wetness factor for the total soil column (equation (section 3.2.3)).

Figure 5.3. Schematic diagram of numerical scheme used to solve for soil water fluxes. Shown are three soil layers, :math:`i - 1` , :math:`i` , and :math:`i + 1` . The soil matric potential :math:`\psi` and volumetric soil water :math:`{\theta _{liq}}` are defined at the layer node depth :math:`z` . The hydraulic conductivity :math:`k\left[ {{z_h}} \right]` is defined at the interface of two layers :math:`{z_h}` . The layer thickness is :math:`\Delta z` . The soil water fluxes :math:`{q_{i - 1}}` and :math:`{q_i}` are defined as positive upwards. The soil moisture sink term :math:`e` (ET loss) is defined as positive for flow out of the layer.

The soil water fluxes in equation , which are a function of :math:`{\theta _{liq,\,i}}` and :math:`{\theta _{liq,\,i + 1}}` because of their dependence on hydraulic conductivity and soil matric potential, can be linearized about :math:`\partial \theta` using a series expansion as



.. math::

 q_i^{n + 1} = q_i^n + \frac{{\partial {q_i}}}{{\partial {\theta _{liq,\,i}}}}\Delta {\theta _{liq,\,i}} + \frac{{\partial {q_i}}}{{\partial {\theta _{liq,\,i + 1}}}}\Delta {\theta _{liq,\,i + 1}}





.. math::

 q_{i - 1}^{n + 1} = q_{i - 1}^n + \frac{{\partial {q_{i - 1}}}}{{\partial {\theta _{liq,\,i - 1}}}}\Delta {\theta _{liq,\,i - 1}} + \frac{{\partial {q_{i - 1}}}}{{\partial {\theta _{liq,\,i}}}}\Delta {\theta _{liq,\,i}}

.

Substitution of these expressions for :math:`q_i^{n + 1}` and :math:`q_{i - 1}^{n + 1}` into equation results in a general tridiagonal equation set of the form



.. math::

 {r_i} = {a_i}\Delta {\theta _{liq,\,i - 1}} + {b_i}\Delta {\theta _{liq,\,i}} + {c_i}\Delta {\theta _{liq,\,i + 1}}



where



.. math::

 {a_i} = - \frac{{\partial {q_{i - 1}}}}{{\partial {\theta _{liq,\,i - 1}}}}





.. math::

 {b_i} = \frac{{\partial {q_i}}}{{\partial {\theta _{liq,\,i}}}} - \frac{{\partial {q_{i - 1}}}}{{\partial {\theta _{liq,\,i}}}} - \frac{{\Delta {z_i}}}{{\Delta t}}





.. math::

 {c_i} = \frac{{\partial {q_i}}}{{\partial {\theta _{liq,\,i + 1}}}}





.. math::

 {r_i} = q_{i - 1}^n - q_i^n + {e_i}

.

The tridiagonal equation set is solved over :math:`i = 1, \ldots ,{N_{levsoi}} + 1` where the layer :math:`i = {N_{levsoi}} + 1` is a virtual layer representing the aquifer.

The finite-difference forms of the fluxes and partial derivatives in equations - can be obtained from equation as



.. math::

 q_{i - 1}^n = - k\left[ {{z_{h,\,i - 1}}} \right]\left[ {\frac{{\left( {{\psi _{i - 1}} - {\psi _i}} \right) + \left( {{\psi _{E,\,i}} - {\psi _{E,\,i - 1}}} \right)}}{{{z_i} - {z_{i - 1}}}}} \right]





.. math::

 q_i^n = - k\left[ {{z_{h,\,i}}} \right]\left[ {\frac{{\left( {{\psi _i} - {\psi _{i + 1}}} \right) + \left( {{\psi _{E,\,i + 1}} - {\psi _{E,\,i}}} \right)}}{{{z_{i + 1}} - {z_i}}}} \right]





.. math::

 \frac{{\partial {q_{i - 1}}}}{{\partial {\theta _{liq,\,i - 1}}}} = - \left[ {\frac{{k\left[ {{z_{h,\,i - 1}}} \right]}}{{{z_i} - {z_{i - 1}}}}\frac{{\partial {\psi _{i - 1}}}}{{\partial {\theta _{liq,\,i - 1}}}}} \right] - \frac{{\partial k\left[ {{z_{h,\,i - 1}}} \right]}}{{\partial {\theta _{liq,\,i - 1}}}}\left[ {\frac{{\left( {{\psi _{i - 1}} - {\psi _i}} \right) + \left( {{\psi _{E,\,i}} - {\psi _{E,\,i - 1}}} \right)}}{{{z_i} - {z_{i - 1}}}}} \right]





.. math::

 \frac{{\partial {q_{i - 1}}}}{{\partial {\theta _{liq,\,i}}}} = \left[ {\frac{{k\left[ {{z_{h,\,i - 1}}} \right]}}{{{z_i} - {z_{i - 1}}}}\frac{{\partial {\psi _i}}}{{\partial {\theta _{liq,\,i}}}}} \right] - \frac{{\partial k\left[ {{z_{h,\,i - 1}}} \right]}}{{\partial {\theta _{liq,\,i}}}}\left[ {\frac{{\left( {{\psi _{i - 1}} - {\psi _i}} \right) + \left( {{\psi _{E,\,i}} - {\psi _{E,\,i - 1}}} \right)}}{{{z_i} - {z_{i - 1}}}}} \right]





.. math::

 \frac{{\partial {q_i}}}{{\partial {\theta _{liq,\,i}}}} = - \left[ {\frac{{k\left[ {{z_{h,\,i}}} \right]}}{{{z_{i + 1}} - {z_i}}}\frac{{\partial {\psi _i}}}{{\partial {\theta _{liq,\,i}}}}} \right] - \frac{{\partial k\left[ {{z_{h,\,i}}} \right]}}{{\partial {\theta _{liq,\,i}}}}\left[ {\frac{{\left( {{\psi _i} - {\psi _{i + 1}}} \right) + \left( {{\psi _{E,\,i + 1}} - {\psi _{E,\,i}}} \right)}}{{{z_{i + 1}} - {z_i}}}} \right]





.. math::

 \frac{{\partial {q_i}}}{{\partial {\theta _{liq,\,i + 1}}}} = \left[ {\frac{{k\left[ {{z_{h,\,i}}} \right]}}{{{z_{i + 1}} - {z_i}}}\frac{{\partial {\psi _{i + 1}}}}{{\partial {\theta _{liq,\,i + 1}}}}} \right] - \frac{{\partial k\left[ {{z_{h,\,i}}} \right]}}{{\partial {\theta _{liq,\,i + 1}}}}\left[ {\frac{{\left( {{\psi _i} - {\psi _{i + 1}}} \right) + \left( {{\psi _{E,\,i + 1}} - {\psi _{E,\,i}}} \right)}}{{{z_{i + 1}} - {z_i}}}} \right]

.

The derivatives of the soil matric potential at the node depth are derived from equation



.. math::

 \frac{{\partial {\psi _{i - 1}}}}{{\partial {\theta _{liq,\,\,i - 1}}}} = - {B_{i - 1}}\frac{{{\psi _{i - 1}}}}{{{\theta _{\,\,i - 1}}}}





.. math::

 \frac{{\partial {\psi _i}}}{{\partial {\theta _{\,liq,\,i}}}} = - {B_i}\frac{{{\psi _i}}}{{{\theta _i}}}





.. math::

 \frac{{\partial {\psi _{i + 1}}}}{{\partial {\theta _{liq,\,i + 1}}}} = - {B_{i + 1}}\frac{{{\psi _{i + 1}}}}{{{\theta _{\,i + 1}}}}



with the constraint :math:`0.01\,{\theta _{sat,\,i}} \leqslant {\theta _{\,i}} \leqslant {\theta _{sat,\,i}}` .

The derivatives of the hydraulic conductivity at the layer interface are derived from equation

\[\begin{gathered}

\frac{{\partial k\left[ {{z_{h,\,i - 1}}} \right]}}{{\partial {\theta _{liq,\,i - 1}}}} = \frac{{\partial k\left[ {{z_{h,\,i - 1}}} \right]}}{{\partial {\theta _{liq,\,i}}}} = \left( {1 - \frac{{{f_{frz,\,i - 1}} + {f_{frz,\,i}}}}{2}} \right)\left( {2{B_{i - 1}} + 3} \right){k_{sat}}\left[ {{z_{h,\,i - 1}}} \right] \times \hfill \\

& {\left[ {\frac{{0.5\left( {{\theta _{\,i - 1}} + {\theta _{\,i}}} \right)}}{{0.5\left( {{\theta _{sat,\,i - 1}} + {\theta _{sat,\,i}}} \right)}}} \right]^{2{B_{i - 1}} + 2}}\left( {\frac{{0.5}}{{{\theta _{sat,\,i - 1}}}}} \right) \hfill \\

\end{gathered} \]

\[\begin{gathered}

\frac{{\partial k\left[ {{z_{h,\,i}}} \right]}}{{\partial {\theta _{liq,\,i}}}} = \frac{{\partial k\left[ {{z_{h,\,i}}} \right]}}{{\partial {\theta _{liq,\,i + 1}}}} = \left( {1 - \frac{{{f_{frz,\,i}} + {f_{frz,\,i + 1}}}}{2}} \right)\left( {2{B_i} + 3} \right){k_{sat}}\left[ {{z_{h,\,i}}} \right] \times \hfill \\

& {\left[ {\frac{{0.5\left( {{\theta _{\,i}} + {\theta _{\,i + 1}}} \right)}}{{0.5\left( {{\theta _{sat,\,i}} + {\theta _{sat,\,i + 1}}} \right)}}} \right]^{2{B_i} + 2}}\left( {\frac{{0.5}}{{{\theta _{sat,\,i}}}}} \right) \hfill \\

\end{gathered} \].
#. 3.2.1 Equilibrium soil matric potential and volumetric moisture

The equilibrium soil matric potential :math:`{\psi _E}` can be derived from equation as



.. math::

 {\psi _E} = {\psi _{sat}}{\left( {\frac{{{\theta _E}\left( z \right)}}{{{\theta _{sat}}}}} \right)^{ - B}}



and the equilibrium volumetric water content :math:`{\theta _E}\left( z \right)` at depth :math:`z` can also be derived as



.. math::

 {\theta _E}\left( z \right) = {\theta _{sat}}{\left( {\frac{{{\psi _{sat}} + {z_\nabla } - z}}{{{\psi _{sat}}}}} \right)^{ - \frac{1}{B}}}

.

Here, the soil matric potentials, the water table depth :math:`{z_\nabla }` and the soil depths have units of mm. For the finite-difference scheme, a layer-average equilibrium volumetric water content is used in equation and can be obtained from



.. math::

 \overline {{\theta _{E,\,i}}} = \int\limits_{{z_{h,\,i - 1}}}^{{z_{h,\,i}}} {\frac{{{\theta _E}\left( z \right)}}{{{z_{h,\,i}} - {z_{h,\,i - 1}}}}} \,dz



which when integrated yields



.. math::

 \overline {{\theta _{E,\,i}}} = \frac{{{\theta _{sat,\,i}}{\psi _{sat,\,i}}}}{{\left( {{z_{h,\,i}} - {z_{h,\,i - 1}}} \right)\left( {1 - \frac{1}{{{B_i}}}} \right)}}\left[ {{{\left( {\frac{{{\psi _{sat,\,i}} - {z_\nabla } + {z_{h,\,i}}}}{{{\psi _{sat,\,i}}}}} \right)}^{1 - \frac{1}{{{B_i}}}}} - {{\left( {\frac{{{\psi _{sat,\,i}} - {z_\nabla } + {z_{h,\,i - 1}}}}{{{\psi _{sat,\,i}}}}} \right)}^{1 - \frac{1}{{{B_i}}}}}} \right]

.

Equation is valid when the water table :math:`{z_\nabla }` is deeper than both interface depths :math:`{z_{h,\,i - 1}}` and :math:`{z_{h,\,i}}` . Since the water table can be within the soil column, the equation is modified if the water table is within soil layer :math:`i` ( :math:`{z_{h,\,i - 1}} < {z_\nabla } < {z_{h,\,i}}` ) as a weighted average of the saturated part and the unsaturated part



.. math::

 \overline {{\theta _{E,\,i}}} = \overline {{\theta _{E,\,sat,\,i}}} \left( {\frac{{{z_{h,\,i}} - {z_\nabla }}}{{{z_{h,\,i}} - {z_{h,\,i - 1}}}}} \right) + \overline {{\theta _{E,\,unsat,\,i}}} \left( {\frac{{{z_\nabla } - {z_{h,\,i - 1}}}}{{{z_{h,\,i}} - {z_{h,\,i - 1}}}}} \right)



where 

.. math::

 \overline {{\theta _{E,\,sat,\,i}}} = {\theta _{sat,\,i}}

 and the unsaturated part 

.. math::

 \overline {{\theta _{E,\,unsat,\,i}}}

 is



.. math::

 \overline {{\theta _{E,\,unsat,\,i}}} = \frac{{{\theta _{sat,\,i}}{\psi _{sat,\,i}}}}{{\left( {{z_\nabla } - {z_{h,\,i - 1}}} \right)\left( {1 - \frac{1}{{{B_i}}}} \right)}}\left[ {1 - {{\left( {\frac{{{\psi _{sat,\,i}} - {z_\nabla } + {z_{h,\,i - 1}}}}{{{\psi _{sat,\,i}}}}} \right)}^{1 - \frac{1}{{{B_i}}}}}} \right]

.

If :math:`{z_\nabla } < {z_{h,\,i - 1}}` , then 

.. math::

 \overline {{\theta _{E,\,i}}} = \overline {{\theta _{E,\,sat,\,i}}} = {\theta _{sat,\,i}}

. If the water table is below the soil column ( :math:`{z_\nabla } > {z_{h,\,{N_{levsoi}}}}` ), an equilibrium volumetric soil moisture is calculated for a virtual layer :math:`i = {N_{levsoi}} + 1` as



.. math::

 \overline {{\theta _{E,\,i = {N_{levsoi + 1}}}}} = \frac{{{\theta _{sat,i - 1}}{\psi _{sat,\,i - 1}}}}{{\left( {{z_\nabla } - {z_{h,\,i - 1}}} \right)\left( {1 - \frac{1}{{{B_{i - 1}}}}} \right)}}\left[ {1 - {{\left( {\frac{{{\psi _{sat,\,i - 1}} - {z_\nabla } + {z_{h,\,i - 1}}}}{{{\psi _{sat,\,i - 1}}}}} \right)}^{1 - \frac{1}{{{B_{i - 1}}}}}}} \right]



The equilibrium volumetric soil moisture is constrained by



.. math::

 0 \leqslant \overline {{\theta _{E,\,i}}} \leqslant {\theta _{sat,\,i}}



The equilibrium soil matric potential is then



.. math::

 {\psi _{E,\,i}} = {\psi _{sat,\,i}}{\left( {\frac{{\overline {{\theta _{E,\,i}}} }}{{{\theta _{sat,\,i}}}}} \right)^{ - {B_i}}} \geqslant - 1 \times {10^8} & \frac{{\overline {{\theta _{E,\,i}}} }}{{{\theta _{sat,\,i}}}} \geqslant 0.01


#. 3.2.2 Equation set for layer :math:`i = 1`

For the top soil layer ( :math:`i = 1` ), the boundary condition is the infiltration rate (section 5.2), :math:`q_{i - 1}^{n + 1} = - q_{infl}^{n + 1}` , and the water balance equation is



.. math::

 \frac{{\Delta {z_i}\Delta {\theta _{liq,\,i}}}}{{\Delta t}} = q_{infl}^{n + 1} + q_i^{n + 1} - {e_i}

.

After grouping like terms, the coefficients of the tridiagonal set of equations for :math:`i = 1` are



.. math::

 {a_i} = 0





.. math::

 {b_i} = \frac{{\partial {q_i}}}{{\partial {\theta _{liq,\,i}}}} - \frac{{\Delta {z_i}}}{{\Delta t}}





.. math::

 {c_i} = \frac{{\partial {q_i}}}{{\partial {\theta _{liq,\,i + 1}}}}





.. math::

 {r_i} = q_{infl}^{n + 1} - q_i^n + {e_i}

.
#. 3.2.3 Equation set for layers :math:`i = 2, \ldots ,{N_{levsoi}} - 1`

The coefficients of the tridiagonal set of equations for :math:`i = 2, \ldots ,{N_{levsoi}} - 1` are



.. math::

 {a_i} = - \frac{{\partial {q_{i - 1}}}}{{\partial {\theta _{liq,\,i - 1}}}}





.. math::

 {b_i} = \frac{{\partial {q_i}}}{{\partial {\theta _{liq,\,i}}}} - \frac{{\partial {q_{i - 1}}}}{{\partial {\theta _{liq,\,i}}}} - \frac{{\Delta {z_i}}}{{\Delta t}}





.. math::

 {c_i} = \frac{{\partial {q_i}}}{{\partial {\theta _{liq,\,i + 1}}}}





.. math::

 {r_i} = q_{i - 1}^n - q_i^n + {e_i}

.
#. 3.2.4 Equation set for layers :math:`i = {N_{levsoi}}, \ldots {N_{levsoi}} + 1`

For the lowest soil layer ( :math:`i = {N_{levsoi}}` ), the bottom boundary condition depends on the depth of the water table. If the water table is within the soil column ( :math:`{z_\nabla } \leqslant {z_{h,\,{N_{levsoi}}}}` ), a zero-flux bottom boundary condition is applied ( :math:`q_i^n = 0` ) and the coefficients of the tridiagonal set of equations for :math:`i = {N_{levsoi}}` are



.. math::

 {a_i} = - \frac{{\partial {q_{i - 1}}}}{{\partial {\theta _{liq,\,i - 1}}}}





.. math::

 {b_i} = - \frac{{\partial {q_{i - 1}}}}{{\partial {\theta _{liq,\,i}}}} - \frac{{\Delta {z_i}}}{{\Delta t}}





.. math::

 {c_i} = 0





.. math::

 {r_i} = q_{i - 1}^n + {e_i}

.

The coefficients for the aquifer layer :math:`i = {N_{levsoi}} + 1` are then



.. math::

 {a_i} = 0





.. math::

 {b_i} = - \frac{{\Delta {z_i}}}{{\Delta t}}





.. math::

 {c_i} = 0





.. math::

 {r_i} = 0

.

If the water table is below the soil column ( :math:`{z_\nabla } > {z_{h,\,{N_{levsoi}}}}` ), the coefficients for :math:`i = {N_{levsoi}}` are



.. math::

 {a_i} = - \frac{{\partial {q_{i - 1}}}}{{\partial {\theta _{liq,\,i - 1}}}}





.. math::

 {b_i} = \frac{{\partial {q_i}}}{{\partial {\theta _{liq,\,i}}}} - \frac{{\partial {q_{i - 1}}}}{{\partial {\theta _{liq,\,i}}}} - \frac{{\Delta {z_i}}}{{\Delta t}}





.. math::

 {c_i} = \frac{{\partial {q_i}}}{{\partial {\theta _{liq,\,i + 1}}}}





.. math::

 {r_i} = q_{i - 1}^n - q_i^n + {e_i}

.

The :math:`i = {N_{levsoi}} + 1` terms are evaluated using



.. math::

 {\psi _{{N_{levsoi}} + 1}} = {\psi _{sat,\,{N_{levsoi}}}}{\left[ {{s_{{N_{levsoi}} + 1}}} \right]^{ - {B_{{N_{levsoi}}}}}} \geqslant - 1 \times {10^8}



 :math:`{z_{{N_{levsoi}} + 1}} = 0.5\left( {{z_\nabla } + {z_{{N_{levsoi}}}}} \right)` 

where



.. math::

 {s_{{N_{levsoi}} + 1}} = 0.5\left( {\frac{{{\theta _{sat,\,{N_{levsoi}}}} + {\theta _{{N_{levsoi}}}}}}{{{\theta _{sat,\,{N_{levsoi}}}}}}} \right) & 0.01 \leqslant {s_{{N_{levsoi}} + 1}} \leqslant 1

,

 :math:`{\psi _{E,\,{N_{levsoi}} + 1}}` is evaluated from equations and , and



.. math::

 \frac{{\partial {\psi _{{N_{levsoi}} + 1}}}}{{\partial {\theta _{liq,\,{N_{levsoi}} + 1}}}} = - {B_{{N_{levsoi}}}}\frac{{{\psi _{{N_{levsoi}} + 1}}}}{{{s_{\,{N_{levsoi}}}}{\theta _{sat,\,{N_{levsoi}}}}}}

.

The coefficients for the aquifer layer :math:`i = {N_{levsoi}} + 1` are then



.. math::

 {a_i} = - \frac{{\partial {q_{i - 1}}}}{{\partial {\theta _{liq,\,i - 1}}}}





.. math::

 {b_i} = - \frac{{\partial {q_{i - 1}}}}{{\partial {\theta _{liq,\,i}}}} - \frac{{\Delta {z_i}}}{{\Delta t}}





.. math::

 {c_i} = 0





.. math::

 {r_i} = q_{i - 1}^n

.

Upon solution of the tridiagonal equation set (Press et al. 1992), the liquid water contents are updated as follows



.. math::

 w_{liq,\,i}^{n + 1} = w_{liq,\,i}^n + \Delta {\theta _{liq,\,i}}\Delta {z_i} & i = 1, \ldots ,{N_{levsoi}}

.

The volumetric water content is



.. math::

 {\theta _i} = \frac{{{w_{liq,\,i}}}}{{\Delta {z_i}{\rho _{liq}}}} + \frac{{{w_{ice,\,i}}}}{{\Delta {z_i}{\rho _{ice}}}}

.

Groundwater-Soil Water Interactions for the Pervious Road
-------------------------------------------------------------

Drainage or sub-surface runoff for the pervious road is based on the SIMTOP scheme (Niu et al. 2005) with a modification to account for reduced drainage in frozen soils. In the work of Niu et al. (2005), the drainage :math:`{q_{drai}}` (kg m-2 s-1) was formulated as



.. math::

 {q_{drai}} = {q_{drai,\,\max }}\exp \left( { - {f_{drai}}{z_\nabla }} \right)

.	
Here, the water table depth :math:`{z_\nabla }` has units of meters. To restrict drainage in frozen soils, Niu et al. (2005) added the following condition

 :math:`{q_{drai}} = 0 & {\text{for }}{w_{ice,\,{N_{levsoi}}}} > {w_{liq,\,{N_{levsoi}}}}` .

In preliminary testing it was found that a more gradual restriction of drainage was required so that the water table depth remained dynamic under partially frozen conditions. The following modification is made to equation



.. math::

 {q_{drai}} = \left( {1 - {f_{imp}}} \right){q_{drai,\,\max }}\exp \left( { - {f_{drai}}{z_\nabla }} \right)



where :math:`{f_{imp}}` is the fraction of impermeable area determined from the ice content of the soil layers interacting with the water table



.. math::

 {f_{imp}} = \frac{{\exp \left[ { - \alpha \left( {1 - \frac{{\sum\limits_{i = jwt}^{i = {N_{levsoi}}} {\frac{{{w_{ice,\,i}}}}{{{w_{ice,\,i}} + {w_{liq,\,i}}}}\Delta {z_i}} }}{{\sum\limits_{i = jwt}^{i = {N_{levsoi}}} {\Delta {z_i}} }}} \right)} \right] - \exp \left( { - \alpha } \right)}}{{1 - \exp \left( { - \alpha } \right)}} \geqslant 0



where :math:`\alpha = 3` is an adjustable scale-dependent parameter, :math:`jwt` is the index of the layer directly above the water table, :math:`{w_{ice,\,i}}` and :math:`{w_{liq,\,i}}` are the ice and liquid water contents of soil layer :math:`i` (kg m-2), and :math:`\Delta {z_i}` is the layer thickness (m). This expression is functionally the same as that used to determine the impermeable fraction (equation ). In equation , the decay factor :math:`{f_{drai}} = 2.5` m-1 and the maximum drainage when the water table depth is at the surface :math:`{q_{drai,\,\max }} = 5.5 \times {10^{ - 3}}` kg m-2 s-1 were determined for global simulations through sensitivity analysis and comparison with observed runoff.

Determination of water table depth :math:`{z_\nabla }` is based on work by Niu et al. (2007). In this approach, a groundwater component is added in the form of an unconfined aquifer lying below the soil column (Figure 5.1). The groundwater solution is dependent on whether the water table is within or below the soil column. Two water stores are used to account for these solutions. The first, :math:`{W_a}` , is the water stored in the unconfined aquifer (mm) and is proportional to the change in water table depth when the water table is below the lower boundary of the hydrologically-active soil column. The second, :math:`{W_t}` , is the actual groundwater which can include water within the soil column. When the water table is below the soil column :math:`{W_t} = {W_a}` . When the water table is within the soil column, :math:`{W_a}` is constant because there is no water exchange between the soil column and the underlying aquifer, while :math:`{W_t}` varies with soil moisture conditions.

In either case, :math:`{W_t}` is first updated as



.. math::

 W_t^{n + 1} = W_t^n + \left( {{q_{recharge}} - {q_{drai}}} \right)\Delta t



where :math:`\Delta t` is the model time step (s), :math:`{q_{recharge}}` is the recharge to the aquifer (kg m-2 s-1), and the drainage :math:`{q_{drai}}` calculated from equation is equivalent to the groundwater discharge.

For the case when the water table is below the soil column, the water stored in the unconfined aquifer :math:`{W_a}` (mm) is updated as



.. math::

 W_a^{n + 1} = W_a^n + \left( {{q_{recharge}} - {q_{drai}}} \right)\Delta t



and :math:`W_t^{n + 1}` is reset as :math:`W_t^{n + 1} = W_a^{n + 1}` . The recharge rate is defined as positive when water enters the aquifer



.. math::

 {q_{recharge}} = \frac{{\Delta {\theta _{liq,\,{N_{levsoi}} + 1}}\Delta {z_{{N_{levsoi}} + 1}}}}{{\Delta t}}



where 

.. math::

 \Delta {\theta _{liq,\,{N_{levsoi}} + 1}} = \theta _{liq,\,{N_{levsoi}} + 1}^{n + 1} - \theta _{liq,\,{N_{levsoi}} + 1}^n

 is the change in liquid water content for layer :math:`i = {N_{levsoi}} + 1` calculated from the solution of the soil water equations (section 5.3), and :math:`\Delta {z_{{N_{levsoi}} + 1}}` (mm) is



.. math::

 \Delta {z_{{N_{levsoi}} + 1}} = z_\nabla ^n - {z_{h,\,{N_{levsoi}}}}

.

The water table depth is calculated from the aquifer water storage scaled by the average specific yield :math:`{S_y} = 0.2` [the fraction of water volume that can be drained by gravity in an unconfined aquifer (Niu et al. 2007)]

 :math:`{z_\nabla } = {z_{h,\,{N_{levsoi}}}} + 25 - \frac{{{W_a}}}{{{{10}^3}{S_y}}}` .

The form of equation originates from the assumption that the initial amount of water in the aquifer is 4800 mm and the corresponding water table depth is one meter below the bottom of the soil column. The water table depth is at the bottom of the soil column ( :math:`{z_\nabla } = {z_{h,\,{N_{levsoi}}}}` ) when the aquifer water is at its prescribed maximum value (5000 mm). The bottom soil layer liquid water content is updated for excess aquifer water as



.. math::

 w_{liq,\,{N_{levsoi}}}^{n + 1} = w_{liq,\,{N_{levsoi}}}^n + \max \left( {0,\,{W_a} - 5000} \right)



and aquifer water is reset to :math:`{W_a} \leqslant 5000` .

For the case when the water table is within the soil column, there is no water exchange between the soil column and the underlying aquifer. However, variations of the water table depth are still calculated as

\[{z_{nabla}} = \left\{ \begin{gathered}

{z_{h,\,jwt + 1}} - \left[ {\frac{{{W_t} - {{10}^3} \times 25{S_y} - \sum\limits_{i = jwt + 2}^{{N_{levsoi}}} {\Delta {z_i}\left( {{\theta _{sat,\,i}} - {\theta _{ice,\,i}}} \right)} }}{{{{10}^3}\left( {{\theta _{sat,\,jwt + 1}} - {\theta _{ice,\,jwt + 1}}} \right)}}} \right] & jwt = 1, \ldots {N_{levsoi}} - 2 \hfill \\

{z_{h,\,jwt + 1}} - \left[ {\frac{{{W_t} - {{10}^3} \times 25{S_y}}}{{{{10}^3}\left( {{\theta _{sat,\,jwt + 1}} - {\theta _{ice,\,jwt + 1}}} \right)}}} \right] & jwt = {N_{levsoi}} - 1 \hfill \\

\end{gathered} \right\}\]

where :math:`jwt` is the index of the layer directly above the water table, and limits are placed on the water table depth as :math:`0.05 \leqslant {z_\nabla } \leqslant 80` . In the work of Niu et al. (2007), the water table depth in this case was calculated with the specific yield determined by the volume of air pores (the pore space not filled with water) within the soil to convert :math:`{W_t}` to a water table depth. However, this was found to result in unstable water table calculations for a significant proportion of grid cells in global simulations. More specifically, when repeatedly forcing the model with a single year of atmospheric data, the temporal evolution of water table depth was significantly different from year to year for some grid cells, with occasional rapid (within a few days) movement of the water table to the soil surface in some cases. This occurred in grid cells with soil water contents near saturation because of the small amount of available pore space. This had deleterious implications for stability of surface fluxes and temperature. In equation , the calculation is based on effective porosity (

.. math::

 {\theta _{sat,\,i}} - {\theta _{ice,\,i}} \geqslant 0.01

) only. Although less defensible from a physical viewpoint, the approach stabilizes the water table calculation for these grid cells and eliminates unrealistic oscillations in surface fluxes and temperature.

In this case, the drainage :math:`{q_{drai}}` is extracted from the soil liquid water in layers within the water table. The partitioning of drainage from these layers is proportional to the layer thickness-weighted hydraulic conductivity as



.. math::

 w_{liq,\,i}^{n + 1} = w_{liq,\,i}^n - \frac{{{q_{drai}}k\left[ {{z_{h,\,i}}} \right]\Delta t\Delta {z_i}}}{{\sum\limits_{i = jwt + 1}^{i = {N_{levsoi}}} {k\left[ {{z_{h,\,i}}} \right]\Delta {z_i}} }} & i = jwt + 1, \ldots ,{N_{levsoi}}



where :math:`\Delta t` is the time step (s).

After the above calculations, two numerical adjustments are implemented to keep the liquid water content of each soil layer ( :math:`{w_{liq,\,i}}` ) within physical constraints of :math:`w_{liq}^{\min } \leqslant {w_{liq,\,i}} \leqslant \left( {{\theta _{sat,\,i}} - {\theta _{ice,\,i}}} \right)\Delta {z_i}` where :math:`w_{liq}^{\min } = 0.01` (mm). First, beginning with the bottom soil layer :math:`i = {N_{levsoi}}` , any excess liquid water in each soil layer ( :math:`w_{liq,\,i}^{excess} = {w_{liq,\,i}} - \left( {{\theta _{sat,\,i}} - {\theta _{ice,\,i}}} \right)\Delta {z_i} \geqslant 0` ) is successively added to the layer above. Any excess liquid water that remains after saturating the entire soil column (plus a maximum surface ponding depth :math:`w_{liq}^{pond} = 10` kg m-2 s-1), is added to drainage :math:`{q_{drai}}` . Second, to prevent negative :math:`{w_{liq,\,i}}` , each layer is successively brought up to :math:`{w_{liq,\,i}} = w_{liq,\,\min }^{}` by taking the required amount of water from the layer below. If this results in :math:`{w_{liq,\,{N_{levsoi}}}} < w_{liq}^{\min }` , then the layers above are searched in succession for the required amount of water ( :math:`w_{liq}^{\min } - {w_{liq,\,{N_{levsoi}}}}` ) and removed from those layers subject to the constraint :math:`{w_{liq,\,i}} \geqslant w_{liq}^{\min }` . If sufficient water is not found, then the water is removed from :math:`{W_t}` and :math:`{q_{drai}}` .

The surface layer liquid water and ice contents for roof, pervious and impervious road are then updated for dew :math:`{q_{sdew}}` , frost :math:`{q_{frost}}` , or sublimation :math:`{q_{subl}}` (section 3.4) as



.. math::

 w_{liq,\,1}^{n + 1} = w_{liq,\,1}^n + {q_{sdew}}\Delta t





.. math::

 w_{ice,\,1}^{n + 1} = w_{ice,\,1}^n + {q_{frost}}\Delta t





.. math::

 w_{ice,\,1}^{n + 1} = w_{ice,\,1}^n - {q_{subl}}\Delta t

.

Sublimation of ice is limited to the amount of ice available.

Runoff from snow-capping
----------------------------

As with other surfaces, urban surfaces are constrained to have a snow water equivalent :math:`{W_{sno}} \leqslant 1000` kg m-2. For snow-capped surfaces, the solid and liquid precipitation reaching the snow surface and dew in solid or liquid form, is separated into solid :math:`{q_{snwcp,ice}}` and liquid :math:`{q_{snwcp,liq}}` runoff terms



.. math::

 {q_{snwcp,ice}} = {q_{grnd,ice}} + {q_{frost}}





.. math::

 {q_{snwcp,liq}} = {q_{grnd,liq}} + {q_{dew}}



and snow pack properties are unchanged. The :math:`{q_{snwcp,ice}}` runoff is sent to the River Transport Model (RTM) where it is routed to the ocean as an ice stream and, if applicable, the ice is melted there. The 

.. math::

 {q_{snwcp,liq}}

 runoff is assigned to the runoff term :math:`{q_{rgwl}}` (e.g. :math:`{q_{rgwl}} = {q_{snwcp,liq}}` ) and included in the liquid water runoff sent to RTM.

Offline Mode
===============

In offline mode (uncoupled to an atmospheric model), the atmospheric forcing required by CLM (Table 1.1) is supplied by observed datasets. The standard forcing provided with the model is a 57-year (1948-2004) dataset that is described in Qian et al. (2006) though alternative observed forcing datasets could also be used. The forcing data is ingested into a data atmosphere model in three “streams”; precipitation ( :math:`P` ) (mm s-1), solar radiation ( :math:`{S_{atm}}` ) (W m-2), and four other fields [atmospheric pressure :math:`{P_{atm}}` (Pa), atmospheric specific humidity :math:`{q_{atm}}` (kg kg-1), atmospheric temperature :math:`{T_{atm}}` (K), and atmospheric wind :math:`{W_{atm}}` (m s-1)]. These are separate streams because they are handled differently according to the type of field and the temporal resolution at which they are provided. In the Qian et al. (2006) dataset, the precipitation stream is provided at six hour intervals and the data atmosphere model prescribes the same precipitation rate for each model time step within the six hour period. The four fields that are grouped together in another stream (pressure, humidity, temperature, and wind) are provided at three hour intervals and the data atmosphere model linearly interpolates these fields to the time step of the model.

The total solar radiation is provided at six hour intervals. The data is fit to the model time step using a diurnal function that depends on the cosine of the solar zenith angle :math:`\mu` to provide a smoother diurnal cycle of solar radiation and to ensure that all of the solar radiation supplied by the six-hourly forcing data is actually used. The solar radiation at model time step :math:`{t_M}` is

\[\begin{gathered}

{S_{atm}}\left( {{t_M}} \right) = \frac{{\frac{{\Delta {t_{FD}}}}{{\Delta {t_M}}}{S_{atm}}\left( {{t_{FD}}} \right)\mu \left( {{t_M}} \right)}}{{\sum\limits_{i = 1}^{\frac{{\Delta {t_{FD}}}}{{\Delta {t_M}}}} {\mu \left( {{t_{{M_i}}}} \right)} }} & {\text{for }}\mu \left( {{t_M}} \right) > 0.001 \hfill \\

{S_{atm}}\left( {{t_M}} \right) = 0 & {\text{for }}\mu \left( {{t_M}} \right) \leqslant 0.001 \hfill \\

\end{gathered} \]

where :math:`\Delta {t_{FD}}` is the time step of the forcing data (6 hours :math:`\times` 3600 seconds hour-1 = 21600 seconds), :math:`\Delta {t_M}` is the model time step (seconds), :math:`{S_{atm}}\left( {{t_{FD}}} \right)` is the six-hourly solar radiation from the forcing data (W m-2), and :math:`\mu \left( {{t_M}} \right)` is the cosine of the solar zenith angle at model time step :math:`{t_M}` (section 2.8). The term in the denominator of equation is the sum of the cosine of the solar zenith angle for each model time step falling within the six hour period. For numerical purposes, 

.. math::

 \mu \left( {{t_{{M_i}}}} \right) \geqslant 0.001

.

The total incident solar radiation :math:`{S_{atm}}` at the model time step :math:`{t_M}` is then split into near-infrared and visible radiation and partitioned into direct and diffuse according to factors derived from one year’s worth of hourly CAM output from CAM version cam3_5_55 as

 :math:`{S_{atm}}\, \downarrow _{vis}^\mu = {R_{vis}}\left( {\alpha {S_{atm}}} \right)` 

 :math:`{S_{atm}}\, \downarrow _{nir}^\mu = {R_{nir}}\left[ {\left( {1 - \alpha } \right){S_{atm}}} \right]` 



.. math::

 {S_{atm}}\,{ \downarrow _{vis}} = \left( {1 - {R_{vis}}} \right)\left( {\alpha {S_{atm}}} \right)





.. math::

 {S_{atm}}\,{ \downarrow _{nir}} = \left( {1 - {R_{nir}}} \right)\left[ {\left( {1 - \alpha } \right){S_{atm}}} \right]

.

where :math:`\alpha` , the ratio of visible to total incident solar radiation, is assumed to be



.. math::

 \alpha = \frac{{{S_{atm}}\, \downarrow _{vis}^\mu + {S_{atm}}\, \downarrow _{vis}^{}}}{{{S_{atm}}}} = 0.5

.

The ratio of direct to total incident radiation in the visible :math:`{R_{vis}}` is



.. math::

 {R_{vis}} = {a_0} + {a_1} \times \alpha {S_{atm}} + {a_2} \times {\left( {\alpha {S_{atm}}} \right)^2} + {a_3} \times {\left( {\alpha {S_{atm}}} \right)^3} & 0.01 \leqslant {R_{vis}} \leqslant 0.99



and in the near-infrared :math:`{R_{nir}}` is



.. math::

 {R_{nir}} = {b_0} + {b_1} \times \left( {1 - \alpha } \right){S_{atm}} + {b_2} \times {\left[ {\left( {1 - \alpha } \right){S_{atm}}} \right]^2} + {b_3} \times {\left[ {\left( {1 - \alpha } \right){S_{atm}}} \right]^3} & 0.01 \leqslant {R_{nir}} \leqslant 0.99



where :math:`{a_0} = 0.17639,\,{a_1} = 0.00380,\,{a_2} = - 9.0039 \times {10^{ - 6}},\,{a_3} = 8.1351 \times {10^{ - 9}}` and :math:`{b_0} = 0.29548,{b_1} = 0.00504,{b_2} = - 1.4957 \times {10^{ - 5}},{b_3} = 1.4881 \times {10^{ - 8}}` are coefficients from polynomial fits to the data.

The additional atmospheric forcing variables required by Table 1.1 are derived as follows. The atmospheric reference height :math:`{z'_{atm}}` (m) is set to 30 m. The directional wind components are derived as \[{u_{atm}} = {v_{atm}} = {{{W_{atm}}} \mathord{\left/

{\vphantom {{{W_{atm}}} {\sqrt 2 }}} \right.

\kern-\nulldelimiterspace} {\sqrt 2 }}\]. The potential temperature :math:`\overline {{\theta _{atm}}}` (K) is set to the atmospheric temperature :math:`{T_{atm}}` . The atmospheric longwave radiation :math:`{L_{atm}}\, \downarrow` (W m-2) is derived from the atmospheric vapor pressure :math:`{e_{atm}}` and temperature :math:`{T_{atm}}` (Idso 1981) as

 :math:`{L_{atm}}\, \downarrow = 0.70 + 5.95 \times {10^{ - 5}} \times 0.01{e_{atm}}\exp \left( {\frac{{1500}}{{{T_{atm}}}}} \right)\sigma T_{atm}^4` 

where



.. math::

 {e_{atm}} = \frac{{{P_{atm}}{q_{atm}}}}{{0.622 + 0.378{q_{atm}}}}



and :math:`\sigma` is the Stefan-Boltzmann constant (W m-2 K-4) (Table 1.4). The fraction of precipitation :math:`P` (mm s-1) falling as rain and/or snow is



.. math::

 {q_{rain}} = P\left( {{f_P}} \right)

,



.. math::

 {q_{snow}} = P\left( {1 - {f_P}} \right)



where



.. math::

 {f_P} = 0 < 0.5\left( {{T_{atm}} - {T_f}} \right) < 1

.

If the user wishes to provide atmospheric forcing data from another source, the data format outlined above will need to be followed with the following exceptions. The data atmosphere model will accept a user-supplied relative humidity :math:`RH` (%) and derive specific humidity :math:`{q_{atm}}` (kg kg-1) from

 :math:`{q_{atm}} = \frac{{0.622{e_{atm}}}}{{{P_{atm}} - 0.378{e_{atm}}}}` 

where the atmospheric vapor pressure :math:`{e_{atm}}` (Pa) is derived from the water ( :math:`{T_{atm}} > {T_f}` ) or ice ( :math:`{T_{atm}} \leqslant {T_f}` ) saturation vapor pressure :math:`e_{sat}^{{T_{atm}}}` as :math:`{e_{atm}} = \frac{{RH}}{{100}}e_{sat}^{{T_{atm}}}` where :math:`{T_f}` is the freezing temperature of water (K) (Table 1.4), and :math:`{P_{atm}}` is the pressure at height :math:`{z_{atm}}` (Pa). The data atmosphere model will also accept a user-supplied dew point temperature :math:`{T_{dew}}` (K) and derive specific humidity :math:`{q_{atm}}` from

 :math:`{q_{atm}} = \frac{{0.622e_{sat}^{{T_{dew}}}}}{{{P_{atm}} - 0.378e_{sat}^{{T_{dew}}}}}` .

Here, :math:`e_{sat}^T` , the saturation vapor pressure as a function of temperature, is derived from Lowe’s (1977) polynomials (section 3.3). If not provided by the user, the atmospheric pressure :math:`{P_{atm}}` (Pa) is set equal to the standard atmospheric pressure :math:`{P_{std}} = 101325` Pa, and surface pressure :math:`{P_{srf}}` (Pa) is set equal to :math:`{P_{atm}}` .

The user may provide the total direct and diffuse solar radiation, :math:`{S_{atm}}\,{ \downarrow ^\mu }` and :math:`{S_{atm}}\, \downarrow` . These will be time-interpolated using the procedure described above and then each term equally apportioned into the visible and near-infrared wavebands (e.g., :math:`{S_{atm}}\, \downarrow _{vis}^\mu = 0.5{S_{atm}}\,{ \downarrow ^\mu }` , :math:`{S_{atm}}\, \downarrow _{nir}^\mu = 0.5{S_{atm}}\,{ \downarrow ^\mu }` ).

Evaluation
=============

Oleson et al. (2008a, b) describe efforts to evaluate the urban model. This includes a quantitative evaluation of model performance at two specific urban sites, an examination of the robustness of the model through sensitivity studies, and a qualitative evaluation of the urban climate produced by the model, with a focus on the characteristics of the simulated heat island. An additional evaluation component not appearing in these two papers is presented below.

Nighttime longwave radiation and surface temperature
--------------------------------------------------------

Nighttime net longwave radiation and air temperature data for an urban canyon in the district of (49ºN, 123ºW) (Nunez and Oke, 1976, 1977) are used to examine the longwave radiation budget and surface temperatures simulated by the model. The canyon is oriented north-south and is located in a mixed light industrial and residential district. The canyon is 79m long, 7.54m wide, and the east and west walls are 7.31m and 5.59m in height, respectively. Walls are concrete, painted flat white with no windows. The canyon floor consists of a 3-5 cm layer of gravel and clay. Weather conditions on the night of September 9-10, 1973 were clear and calm. Air temperature and net longwave radiation measured at about 0.3m above the midpoint of the canyon floor and from the mid-height of each wall are compared with simulated canyon floor and wall surface temperature and net longwave radiation.

The observation site has been used to validate other urban models such as SHIM (Surface Heat Island Model) (Johnson et al. 1991), the Town Energy Budget (TEB) scheme (Masson 2000), NSLUCM (Noah land surface model/Single-layer Urban Canopy Model) (Kusaka et al. 2001), and VUCM (Vegetated Urban Canopy Model) (Lee and Park 2007). Published data from Lee and Park (2007) were used to determine input parameters for the urban model as these data appeared to produce the best simulations compared to observations (Table 7.1). The canyon floor was modeled as a sandy clay soil with no moisture content. No anthropogenic fluxes were prescribed. Atmospheric wind speed at 10m height was set to 2 m s-1 and specific humidity to 0.01 kg kg-1 throughout the simulation (Lee and Park 2007). Atmospheric air temperature was initialized at 19 ºC and set to the calculated canyon air temperature on subsequent time steps to maintain a neutral temperature profile (no thermal turbulent fluxes between the canyon and the atmosphere) (Masson 2000). Specific humidity of canyon air is set to the atmospheric specific humidity. Downward longwave radiation was initialized to 339 W m-2 and decreased linearly with the atmospheric air temperature (Masson 2000). Initial wall and canyon floor temperatures were set to 18.35ºC and 18.5ºC, respectively per Johnson et al. (1991).

Table 7.1. Urban model parameters for the site

Figure 7.1 shows the simulated surface temperatures and net longwave radiation for the walls and canyon floor compared to observations. The urban model does a good job simulating the nighttime cooling of canyon surfaces (note that the simulated west and east wall surface temperatures are the same). Temperature differences from observations are less than 1ºC at all times. Net longwave radiation is also well simulated, differences from observations are less than about 3 W m-2 for the west wall and canyon floor. The simulated net longwave radiation for the east wall is biased high by up to 7 W m-2. These results are quite similar to those from VUCM and generally slightly better than the models of Masson (2000), Johnson et al. (1991), and Kusaka et al. (2001) which generally have warmer surface temperatures as noted by Lee and Park (2007). However, one important difference between Lee and Park (2007) and the other studies is that the thermal admittance prescribed for the canyon floor is substantially lower in VUCM. When higher thermal admittance is prescribed in the urban model, warmer surface temperatures are simulated consistent with the other studies.

Figure 7.1. Simulated surface temperatures (solid lines) and net longwave radiation (dashed lines) compared to observations (circles) for a) west (east-facing) wall, b) east wall, and c) canyon floor for the night of September 9-10, 1973 in an urban canyon in the Grandview district of Vancouver, British Columbia. Observed data were digitized from Figure 5 in Johnson et al. (1991).

References
=============

, E.A. 1976. A point energy and mass balance model of a snow cover. NOAA Technical Report NWS 19, Office of Hydrology, National Weather Service,

Arnfield, A.J. 2003. Two decades of urban climate research: a review of turbulence, exchanges of energy and water, and the urban heat island. Int. J. Climatol. 23:1-26.

Arya, S.P. 2001. Introduction to Meteorology. Academic Press, .

Atkinson, B.W. 2003. Numerical modeling of urban heat-island intensity. Bound.-Layer Meteor. 109:285-310.

Auch, R., Taylor, J., and Acevedo, W. 2004. Urban growth in American cities: Glimpses of U.S. urbanization. Circular 1252, Geological Department of the Interior, 52 pp.

Avissar, R. 1996. Potential effects of vegetation on the urban thermal environment. Atmos. Environ. 30:437-448.

Best, M.J. 2005. Representing urban areas within operational numerical weather prediction models. Bound.-Layer Meteor. 114:91-109.

Best, M.J. 2006. Progress towards better weather forecasts for city dwellers: from short range to climate change. Theor. Appl. Climatol. 84:47-55. DOI:10.1007/s00704-005-0143-2.

Betts, R.A. 2001. Biogeophysical impacts of land use on present-day climate: near-surface temperature change and radiative forcing. Atmos. Sci. Lett. 2:39-51. DOI:10.1006/asle.2001.0023.

Beven, K.J., and Kirkby, M.J. 1979. A physically based variable contributing area model of basin hydrology. Hydrol. Sci. Bull. 24: 43-69.

Bonan, G.B. 1996. A land surface model (LSM version 1.0) for ecological, hydrological, and atmospheric studies: technical description and user’s guide. NCAR Technical Note NCAR/TN-417+STR, for Atmospheric Research, Boulder, CO, 150 pp.

Bonan, G.B. 2002. Ecological climatology: concepts and applications. Press, 678 pp.

Bounoua, L., DeFries, R., Collatz, G.J., Sellers, P., and Khan, H. 2002. Effects of land cover conversion on surface climate. Clim. Change 52:29-64.

Brovkin, V., Sitch, S., von Bloh, W., Claussen, M., Bauer, E., and Cramer, W. 2004. Role of land cover changes for atmospheric CO2 increase and climate change during the last 150 years. Global Change Biol. 10:1253-1266. DOI:10.1111/j.1365-2486.2004.00812.x.

Brown, M. 2000. Urban parameterizations for mesoscale meteorological models. pp. 193-255. In: Z. Boybeyi (editor) Mesoscale Atmospheric Dispersion. WIT Press, Southampton, .
#. Inadvertent weather modification in urban areas: Lessons for global climate change. Bull. Amer. Meteor. Soc. 73:619-627.

Chylek, P., Srivastava, V., Cahenzli, L., Pinnick, R.G., Dod, R.L., Novakov, T., Cook, T.L., and Hinds, B.D. 1987. Aerosol and graphitic carbon content of snow. J. Geophys. Res. 92:9801-9809.

CIESIN (Center for International Earth Science Information Network), ; International Food Policy Research Institute (IPFRI), the World Bank; and Centro Internacional de Agricultura Tropical (CIAT), 2004. Global Rural-Urban Mapping Project (GRUMP): Urban Extents. , : CIESIN, . Available at http://beta.sedac.ciesin.columbia.edu/gpw.

Clapp, R.B., and Hornberger, G.M. 1978. Empirical equations for some soil hydraulic properties. Water Resour. Res. 14:601-604.

Clauser, C., and Huenges, E. 1995. Thermal conductivity of rocks and minerals. pp. 105-126. In: T. J. Ahrens (editor) Rock Physics and Phase Relations: A Handbook of Physical Constants. Washington, D.C.

Comrie, A.C. 2000. Mapping a wind-modified urban heat island in (with comments on integrating research and undergraduate learning). Bull. Amer. Meteor. Soc. 81:2417-2431.

Copeland, J.H., Pielke, R.A., and Kittel, T.G.F. 1996. Potential climatic impacts of vegetation change: a regional modeling study. J. Geophys. Res. 101:7409-7418.

Cosby, B.J., Hornberger, G.M., Clapp, R.B., and Ginn, T.R. 1984. A statistical exploration of the relationships of soil moisture characteristics to the physical properties of soils. Water Resour. Res. 20:682-690.

Cramer, W., et al. 2001. Global response of terrestrial ecosystem structure and function to CO2 and climate change: results from six dynamic global vegetation models. Global Change Biol. 7:357-373.

Dai, Y., and Zeng, Q. 1997. A land surface model (IAP94) for climate studies. Part I: formulation and validation in off-line experiments. Adv. Atmos. Sci. 14:433-460.

Dandou, A., Tombrou, M., Akylas, E., Soulakellis, N., and Bossioli, E. 2005. Development and evaluation of an urban parameterization scheme in the Penn State/NCAR Mesoscale Model (MM5). J. Geophys. Res. 110:D10102. DOI: 10.1029/2004JD005192.

Dobson, J.E., Bright, E.A., Coleman, P.R., Durfee, R.C., and Worley, B.A. 2000. LandScan: A global population database for estimating populations at risk. Photogrammetric Engineering and Remote Sensing 66(7):849-857.

de Vries, D.A. 1963. Thermal Properties of Soils. In: W.R. van Wijk (editor) Physics of the Plant Environment. North-Holland, .

Eastman, J.L., Coughenour, M.B., and Pielke Sr., R.A. 2001. The regional effects of CO2 and landscape change using a coupled plant and meteorological model. Global Change Biol. 7:797-815.

Elvidge, C.D., Sutton, P.C., Wagner, T.W., Ryzner, R., Vogelmann, J.E., Goetz, S.J., Smith, A.J., Jantz, C., Seto, K.C., Imhoff, M.L., Wang, Y.Q., Milesi, C. and Nemani, R. 2004. Urbanization. pp. 315-328. In: G. Gutman, A. C. Janetos, C.O. Justice, E.F. Moran, J.F. Mustard, R.R. Rindfuss, D. Skole, B.L. Turner II, and M.A Cochrane (editors), Land Change Science: Observing, Monitoring and Understanding Trajectories of Change on the Earth’s Surface. Kluwer Academic Publishers, The Netherlands.

Entekhabi, D., and Eagleson, P.S. 1989. Land surface hydrology parameterization for atmospheric general circulation models including subgrid scale spatial variability. J. Climate 2:816-831.

Farouki, O.T. 1981. The thermal properties of soils in cold regions. Cold Regions Sci. and Tech. 5:67-75.

Feddema, J.J., Oleson, K.W., Bonan, G., Mearns, L.O., Buja, L.E., Meehl, G.A., and Washington, W.M. 2005: The importance of land cover change in simulating future climates. Science 310:1674-1678.

Flatau, P.J., Walko, R.L., and Cotton, W.R. 1992. Polynomial fits to saturation vapor pressure. J. Appl. Meteor. 31:1507-1513.

Foley, J.A., et al. 2005. Global consequences of land use. Science 309:570-574.

Fu, C. 2003. Potential impacts of human-induced land cover change in monsoon. Glob. Planet. Change 37:219-229. DOI:10.1016/S0921-8181(02)00207-2.

Grimmond, C.S.B., Cleugh, H.A., and Oke, T.R. 1991. An objective urban heat storage model and its comparison with other schemes. Atmos. Environ. 25B:311-326.

Grimmond, C.S.B., and Oke, T.R. 1999. Aerodynamic properties of urban areas derived from analysis of surface form. J. Appl. Meteor. 38:1262-1292.

Grimmond, C.S.B., and Oke, T.R. 2002. Turbulent heat fluxes in urban areas: observations and a local-scale urban meteorological parameterization scheme (LUMPS). J. Appl. Meteor. 41:792-810.

Harman, I.N., Best, M.J., and Belcher, S.E. 2004. Radiative exchange in an urban street canyon. Bound.-Layer Meteor. 110:301-316.

Houghton, J.T., Ding, Y., Griggs, D.J., Noguer, M., van der Linden, P.J., Dai, X., Maskell, K., and Johnson, C.A. (editors) 2001. Climate Change 2001: The Scientific Basis. Press, 881 pp.

Ichinose, T., Shimodozono, K., Hanaki, K. 1999. Impact of anthropogenic heat on urban climate in Tokyo. Atmos. Environ. 33:3897-3909.

Idso, S.B. 1981. A set of equations for full spectrum and 8- to 14- :math:`\mu` m and 10.5- to 12.5- :math:`\mu` m thermal radiation from cloudless skies. Water Resour. Res. 17:295-304.

Jackson, T.L., Feddema, J.J., Oleson, K.W., Bonan, G.B., and Bauer, J.T. 2010. Parameterization of urban characteristics for global climate modeling. A. Assoc. Am. Geog., in press.

Jin, M., Dickinson, R.E., and Zhang, D.-L. 2005. The footprint of urban areas on global climate as characterized by MODIS. J. Climate 18:1551-1565.

Johnson, G.T., Oke, T.R., Lyons, T.J., Steyn, D.G., Watson, I.D., and Voogt, J.A. 1991. Simulation of surface urban heat islands under ‘ideal’ conditions at night. Part 1: Theory and tests against field data. Bound.-Layer Meteor. 56:275-294.

Jordan, R. 1991. A One-dimensional Temperature Model for a Snow Cover: Technical Documentation for SNTHERM.89. Army Cold Regions Research and Engineering Laboratory, Special Report 91-16.

Kalnay, E., and Cai, M. 2003. Impact of urbanization and land-use change on climate. Nature 423:528-531.

Kusaka, H., Kondo, H., Kikegawa, Y., and Kimura, F. 2001. A simple single-layer urban canopy model for atmospheric models: comparison with multi-layer and slab models. Bound.-Layer Meteor. 101:329-358.

Kusaka, H., and Kimura, F. 2004. Thermal effects of urban canyon structure on the nocturnal heat island: numerical experiment using a mesoscale model coupled with an urban canopy model. J. Appl. Meteor. 43:1899-1910.

Landsberg, H.E. 1981. The Urban Climate. , Academic Press, 275 pp.

Lawrence, D.M., and Slater, A.G. 2008. Incorporating organic soil into a global climate model. Clim. Dyn. 30. DOI:10.1007/s00382-007-0278-1.

Lee, S.-H., and Park, S.-U. 2007. A vegetated urban canopy model for meteorological and environmental modeling. Bound.-Layer Meteor. DOI:10:1007/s10546-007-9221-6.

Lemonsu, A., and Masson, V. 2002. Simulation of a summer urban breeze over . Bound.-Layer Meteor. 104:463-490.

Lemonsu, A., Grimmond, C.S.B., and Masson, V. 2004. Modeling the surface energy balance of the core of an old Mediterranean city: Marseille. J. Appl. Meteor. 43:312-327.

Lowe, P.R. 1977. An approximating polynomial for the computation of saturation vapor pressure. J. Appl. Meteor. 16:100-103.

Macdonald, R.W., Griffiths, R.F., and Hall, D.J. 1998. An improved method for the estimation of surface roughness of obstacle arrays. Atmos. Environ. 32:1857-1864.

Marshall, S.E. 1989. A physical parameterization of snow albedo for use in climate models. NCAR Cooperative Thesis NCAR/CT-123, for Atmospheric Research, Boulder, CO.

Martilli, A., Clappier, A., Rotach, M.W. 2002. An urban surface exchange parameterization for mesoscale models. Bound.-Layer Meteor. 104:261-304.

Masson, V. 2000. A physically-based scheme for the urban energy budget in atmospheric models. Bound.-Layer Meteor. 94:357-397.

Masson, V., Grimmond, C.S.B., and Oke, T.R. 2002. Evaluation of the Town Energy Balance (TEB) scheme with direct measurements from dry districts in two cities. J. Appl. Meteor. 41:1011-1026.

Masson, V. 2006. Urban surface modeling and the meso-scale impact of cities. Theor. Appl. Climatol. 84:35-45.

Mathews, H.D., Weaver, A.J., Meissner, K.J., Gillett, N.P., and Eby, M. 2004. Natural and anthropogenic climate change: incorporating historical land cover change, vegetation dynamics and the global carbon cycle. Clim. Dyn. 22:461-479. DOI:10.1007/s00382-004-0392-2.

Myhre, G., and Myhre, A. 2003. Uncertainties in radiative forcing due to surface albedo changes caused by land-use changes. J. Climate 16:1511-1524.

Narisma, G.T., and Pitman, A.J. 2003. The impact of 2000 years of land cover change on the Australian near-surface climate. J. Hydrometeor. 4:424-436.

Niu, G.-Y., Yang, Z.-L., Dickinson, R.E., and Gulden, L.E. 2005. A simple TOPMODEL-based runoff parameterization (SIMTOP) for use in global climate models. J. Geophys. Res. 110:D21106. DOI: 10.1029/2005JD006111.

Niu, G.-Y., and Yang, Z.-L. 2006. Effects of frozen soil on snowmelt runoff and soil water storage at a continental scale. J. Hydrometeor. 7:937-952.

Niu, G.-Y., Yang, Z.-L., Dickinson, R.E., Gulden, L.E., and Su, H. 2007. Development of a simple groundwater model for use in climate models and evaluation with Gravity Recovery and Climate Experiment data. J. Geophys. Res. 112:D07103. DOI:10.1029/2006JD007522.

Nunez, M., and Oke, T.R. 1976. Long-wave radiative flux divergence and nocturnal cooling of the urban atmosphere. II: Within and urban canyon. Bound.-Layer Meteor. 10:121-135.

Nunez, M., and Oke, T.R. 1977. The energy balance of an urban canyon. J. Appl. Meterorol. 16:11-19.

Oke, T.R. 1981. Canyon geometry and the nocturnal urban heat island: comparison of scale model and field observations. J. Climatol. 1:237-254.

Oke, T.R. 1982. The energetic basis of the urban heat island. Quart. J. Royal Meteor. Soc. 108:1-24.

Oke, T.R. 1987. Boundary Layer Climates (2nd edition). , Routledge, 435 pp.

Oke, T.R., and Cleugh, H.A. 1987. Urban heat storage derived as energy balance residuals. Bound.-Layer Meteor. 39:233-245.

Oke, T.R., Johnson, G.T., Steyn, D.G., and Watson, I.D. 1991. Simulation of surface urban heat islands under “ideal” conditions at night, part 2: diagnosis of causation. Bound.-Layer Meteor. 56:339-358.

Oleson, K.W., et al. 2004. Technical description of the Community Land Model (CLM). NCAR Technical Note NCAR/TN-461+STR, National Center for Atmospheric Research, Boulder, CO, 173 pp.

Oleson, K.W., Bonan, G.B., Feddema, J., Vertenstein, M., and Grimmond, C.S.B. 2008a. An urban parameterization for a global climate model. 1. Formulation and evaluation for two cities. J. Appl. Meteor. Clim. 47:1038-1060.

Oleson, K.W., Bonan, G.B., Feddema, J., and Vertenstein, M. 2008b. An urban parameterization for a global climate model. 2. Sensitivity to input parameters and the simulated urban heat island in offline simulations. J. Appl. Meteor. Clim. 47:1061-1076.

Oleson, K.W., Bonan, G.B., and Feddema, J. 2010a. The effects of white roofs on urban temperature in a global climate model. Geophys. Res. Lett. 37:L03701. DOI:10.1029/2009GL042194.

Oleson, K.W., et al. 2010b. Technical description of version 4.0 of the Community Land Model (CLM). NCAR Technical Note NCAR/TN-478+STR, National Center for Atmospheric Research, Boulder, CO, 257 pp.

Otte, T.L., Lacser, A., Dupont, S., and Ching, J.K.S. 2004. Implementation of an urban canopy parameterization in a mesoscale meteorological model. J. Appl. Meteor. 43:1648-1665.

Panofsky, H.A., and Dutton, J.A. 1984. Atmospheric Turbulence: Models and Methods for Engineering Applications. John Wiley and Sons, .

Pielke Sr., R.A., Marland, G., Betts, R.A., Chase, T.N., Eastman, J.L., Niles, J.O., Niyogi, D.D.S., and Running, S.W. 2002. The influence of land-use change and landscape dynamics on the climate system: relevance to climate-change policy beyond the radiative effect of greenhouse gases. Phil. Trans. R. Soc. Lond. A, 360:1705-1719.

Piringer, M., Grimmond, C.S.B., Joffre, S.M., Mestayer, P., Middleton, D.R., Rotach, A. Baklanov, M.W., De Ridder, K., Ferreira, J., Guilloteau, E., Karppinen, A., Martilli, A., Masson, V., and Tombrou, M. 2002. Investigating the surface energy balance in urban areas - recent advances and future needs. Water, Air, Soil Pollution: Focus, 2:1-16.

Press, W.H., Teukolsky, S.A., Vetterling, W.T., and Flannery, B.P. 1992. Numerical Recipes in FORTRAN: The Art of Scientific Computing. Press, .

Qian, T., Dai, A., Trenberth, K.E., and Oleson, K.W. 2006. Simulation of global land surface conditions from 1948 to 2004: Part I: Forcing data and evaluations. J. Hydrometeor. 7:953-975.

Rowley, F.B., Algren, A.B., and Blackshaw, J.L. 1930. Surface conductances as affected by air velocity, temperature, and character of surface. ASHRAE Trans. 36:429-446.

Sailor, D.J. 1995. Simulated urban climate response to modifications in surface albedo and vegetative cover. J. App. Meteor. 34:1694-1704.

Sailor, D.J., and Lu, L. 2004. A top-down methodology for developing diurnal and seasonal anthropogenic heating profiles for urban areas. Atmos. Environ. 38:2737-2748.

Sakakibara, Y. 1996. A numerical study of the effect of urban geometry upon the surface energy budget. Atmos. Environ. 30:487-496.

Sellers, P.J., Dickinson, R.E., Randall, D.A., Betts, A.K., Hall, F.G., Berry, J.A., Collatz, G.J., Denning, A.S., Mooney, H.A., Nobre, C.A., Sato, N., Field, C.B., and Henderson-Sellers, A. 1997. Modeling the exchanges of energy, water, and carbon between continents and the atmosphere. Science 275:502-509.

Shepherd, J.M. 2005. A review of current investigations of urban-induced rainfall and recommendations for the future. Earth Interact. 9:1-27.

Sparrow, E.M., and Cess, R.D. 1978. Radiation Heat Transfer. Hemisphere Publishing Corporation, 366 pp.

Stohlgren, T.J., Chase, T.N., Pielke Sr., R.A., Kittel, T.G.F., and Baron, J.S. 1998. Evidence that local land use practices influence regional climate, vegetation, and stream flow patterns in adjacent natural areas. Global Change Biol. 4:495-504.

Stull, R.B. 1988. An Introduction to Boundary Layer Meteorology. Kluwer Academic Publishers, .

Taha, H. 1999. Modifying a mesoscale meteorological model to better incorporate urban heat storage: a bulk-parameterization approach. J. Appl. Meteor. 38:466-473.

Upmanis, H., Eliasson, I., and Lindqvist, S. 1998. The influence of green areas on nocturnal temperatures in a high latitude city (). Int. J. Climatol. 18:681-700.

Wang, H., Pitman, A.J., Zhao, M., and Leemans, R. 2003. The impact of land-cover modification on the June meteorology of since 1700, simulated using a regional climate model. Int. J. Climatol. 23:511-527.

Yang, Z.-L. 1998. Technical note of a 10-layer soil moisture and temperature model. Unpublished manuscript.

Zeng, X., Zhao, M., and R.E. 1998. Intercomparison of bulk aerodynamic algorithms for the computation of sea surface fluxes using the TOGA COARE and TAO data. J. Climate 11:2628-2644.

Zeng, X., and Decker, M. 2009. Improving the numerical solution of soil moisture-based Richards equation for land models with a deep or shallow water table. J. Hydrometeor. 10:308-319.

Zhou, L., Dickinson, R.E., Tian, Y., Fang, J., Li, Q., Kaufmann, R.K., Tucker, C.J., and Myneni, R.B. 2004. Evidence for a significant urbanization effect on climate in China. Proc. Natl. Acad. Sci. U.S.A. 101:9540-9544.
