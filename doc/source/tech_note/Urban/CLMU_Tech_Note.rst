.. _rst_Urban Model (converted from docx):

*********************************
Urban Model (converted from docx)
*********************************

Introduction
============

This technical note describes the physical parameterizations and numerical implementation of a Community Land Model Urban (CLMU) parameterization as coupled to version 4 of the Community Land Model (CLM4). CLM4 serves as the land surface model component of the Community Atmosphere Model (CAM) and the Community Climate System Model (CCSM). This note documents the global implementation of the urban model. Other model versions may exist for specific applications.

Chapters 1-5 constitute the description of the urban parameterization when coupled to CAM or CCSM, while Chapter 6 describes processes that pertain specifically to the operation of the urban parameterization in offline mode (uncoupled to an atmospheric model). Chapter 7 describes efforts to evaluate the urban model. The model formulation and some quantitative and qualitative evaluation are also documented in Oleson et al. (2008a, 2008b). A heat island mitigation study using the model is presented in Oleson et al. (2010a). Note that CLMU and CLM4 have some parameterizations in common (e.g., snow and sub-surface hydrology). This technical note contains material duplicated from the CLM4 technical note (Oleson et al. 2010b) where appropriate. This is done so that users interested in just the urban model do not have to refer to the CLM4 technical note.

Model Overview
--------------

Motivation
~~~~~~~~~~

Land use and land cover change is increasingly being recognized as an important yet poorly quantified component of global climate change (Houghton et al. 2001). Land use/cover change mechanisms include both the transformation of natural land surfaces to those serving human needs (i.e., direct anthropogenic change) (e.g., the conversion of tropical forest to agriculture) as well as changes in land cover on longer time-scales due to biogeophysical feedbacks between the atmosphere and the land (i.e., indirect change) (Cramer et al. 2001, Foley et al. 2005). Global and regional models have been used extensively to investigate the effects of direct and indirect land use/cover change mechanisms on climate (Copeland et al. 1996, Stohlgren et al. 1998, Betts 2001, Eastman et al. 2001, Bounoua et al. 2002, Pielke et al. 2002, Fu 2003, Myhre and Myhre 2003, Narisma and Pitman 2003, Wang et al. 2003, Brovkin et al. 2004, Mathews et al. 2004, Feddema et al. 2005). However, all of these studies have focused on land use/cover related to changes in vegetation types. Urbanization, or the expansion of built-up areas, is an important yet less studied aspect of anthropogenic land use/cover change in climate science.

Although currently only about 1-3% of the global land surface is urbanized, the spatial extent and intensity of urban development is expected to increase dramatically in the future (Shepherd 2005, CIESIN et al. 2004). More than one-half of the world's population currently lives in urban areas and in Europe, North America, and Japan at least 80% of the population resides in urban areas (Elvidge et al. 2004). Policymakers and the public are most interested in the effects of climate change on people where they live. Because urban and non-urban areas may have different sensitivities to climate change, it is possible that the true climate change signal within urban areas may only be estimated if urban areas are explicitly modeled in climate change simulations (Best 2006). Indeed, the "footprint" of urbanization on climate can be detected from surface observations and satellite data (Changnon 1992, Kalnay and Cai 2003, Zhou et al. 2004, Jin et al. 2005). Changnon (1992) points out that the average urban warming over the last 100 years in certain regions is comparable to the increase in global surface temperature predicted by climate models for the next 100 years. Thus, it is important for developers of land surface models to begin to consider the parameterization of urban surfaces.

Urbanization now appropriates significant proportions of land in certain regions. For example, the expansion of service-based industries and conversion of farmland for housing in the Chicago area has increased the amount of developed land from about 800 square miles in 1973 to 1000 square miles in 1992 (Auch et al. 2004). A T85 resolution climate model grid cell (the resolution of the CCSM3 climate change simulations submitted for the IPCC AR4) encompassing the Chicago region represents about 7100 square miles, which suggests that this grid cell should be modeled as about 14% urbanized land. For mesoscale or regional models, where grid cells are on the order of a few kilometers, an urban area this size will occupy a significant number of grid cells that would otherwise be modeled as natural surfaces. The now common use of multiple "tiles" in models enables the co-existence of multiple surface types within a single gridcell. Thus, urban areas should and can be included in a global climate model (Best 2006).

Numerical modeling of the urban energy budget was first attempted nearly 40 years ago [see Brown (2000) for a comprehensive historical overview of modeling efforts]. However, until recently, most modern land surface models [i.e., second- or third-generation models (Sellers et al. 1997)] have not formally included urban parameterizations. Masson (2006) classifies urban parameterizations in three general categories: 1) empirical models, 2) vegetation models, with and without drag terms, adapted to include an urban canopy, 3) single-layer and multi-layer models that include a three-dimensional representation of the urban canopy. Empirical models (e.g., Oke and Cleugh 1987) rely on statistical relationships determined from observed data. As such, they are generally limited to the range of conditions experienced during the observation campaign. Vegetation models adapted for the urban canopy generally focus on modifying important surface parameters to better represent urban surfaces [e.g., surface albedo, roughness length, displacement height, surface emissivity, heat capacity, thermal conductivity (Taha 1999, Atkinson 2003, Best 2005].

These relatively simple approaches (i.e., categories 1 and 2 above) may arguably be justified based on the fact that detail in complex models may be lost when averaged to a coarse grid (Taha 1999). However, they may not have sufficient functionality to be suitable for inclusion in global climate models and may require the global derivation of parameters that are difficult to interpret physically [e.g., the surface type-dependent empirical coefficients for storage heat flux in the Objective Hysteresis Model (Grimmond et al. 1991)]. Furthermore, such approaches may not fully describe the fundamental processes that determine urban effects on climate (Piringer et al. 2002). For example, cities are known to have unique characteristics that cause them to be warmer than surrounding rural areas, an effect known as the urban heat island (Oke 1987). In the absence of anthropogenic heat flux, the urban heat island is thought to be greatest on clear, calm nights when local conditions generally dominate over synoptic. Candidate causes for this phenomenon include decreased surface longwave radiation loss and increased absorption of solar radiation because of canyon geometry, anthropogenic emissions of heat, reduction of evapotranspiration due to the replacement of vegetation with impervious surfaces, increased downwelling longwave radiation from the atmosphere due to pollution and warmer atmospheric temperatures, increased storage of sensible heat within urban materials, and reduced transfer of heat due to sheltering from buildings (Oke 1982, Oke 1987, Oke et al. 1991). Single-layer or multi-layer urban canopy models are likely needed to investigate the relative contribution of these factors to the heat island effect (Piringer et al. 2002). For example, specification of an urban albedo may provide no insight into the effects of the individual albedo of roofs, walls, and roads and the interaction of shortwave radiation between these surfaces that yields urban albedos that are typically lower than those of most rural sites. Similarly, assessments of the effectiveness of techniques proposed to ameliorate heat islands, such as "green roofs" or tree planting, require more detailed models.

On the other hand, the level of complexity in a model is limited by the availability of data that the model requires, the computational burden imposed, and difficulty in understanding the complex behavior of the model. Here, following recent developments in detailed urban parameterizations designed for mesoscale models (Masson 2000, Martilli et al. 2002, Grimmond and Oke 2002, Kusaka and Kimura 2004, Otte et al. 2004, Dandou et al. 2005), we describe a model that is simple enough to be compatible with structural, computational, and data constraints of a land surface model coupled to a global climate model, yet complex enough to enable exploration of physically-based processes known to be important in determining urban climatology. Several of the parameterizations are based on the Town Energy Balance (TEB) Model (Masson 2000, Masson et al. 2002, Lemonsu et al. 2004).

Urban Ecosystems and Climate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Characteristics of urban ecosystems and their effects on climate are summarized in Landsberg (1981), Oke (1987), Bonan (2002), and Arnfield (2003). Urban ecosystems can significantly alter the radiative, thermal, moisture, and aerodynamic characteristics of a region. The three-dimensional structure and geometrical arrangement of building walls and horizontal surfaces such as roads, sidewalks, parking lots, etc. combine to reduce the albedo of urban surfaces due to radiation trapping. Unlike solar radiation reflected from a horizontal surface, solar radiation impinging on urban surfaces such as walls and roads can experience multiple reflections and absorptions, resulting in increased absorption of radiation. Similarly, longwave radiation emitted by urban surfaces can be re-absorbed by these surfaces resulting in less longwave radiation loss to the atmosphere. The ratio of building height to canyon floor width is important in determining the degree of radiation trapping (Oke 1981, Oke et al. 1991).

The materials used for the construction of buildings and roads (e.g., dense concrete and asphalt) generally have higher heat capacity and thermal conductivity than some natural surfaces such as dry soils (Oke 1987). This results in higher thermal admittance and contributes to the ability of urban surfaces to store sensible heat during the day and release it at night. The importance of thermal properties in contributing to differences between urban and rural sites depends on the types of materials used in urban construction, the contrast in thermal admittance between the urban region and surrounding rural environs, and the building geometry which establishes the relative surface area and importance of roof, walls, and canyon floor (Oke et al. 1991).

Energy consumption due to building heating and cooling, manufacturing, transportation, and human metabolism releases waste heat to the urban environment. Such anthropogenic sources of heat can be substantial in some cases and should be accounted for in studies of the urban energy budget. As an extreme example, Ichinose et al. (1999) found that the total anthropogenic heat flux in central Tokyo exceeded 400 W m\ :sup:`-2` in daytime and a maximum value of 1590 W m\ :sup:`-2` in winter. The contribution of waste heat sources from building heating and cooling may depend on population density, external climate, and socio-economic factors such as human adaptability and comfort levels, and economic status. The presence of insulation, characterized by low thermal admittance, may reduce the contribution of waste heat from heating and cooling. Waste heat fluxes from transportation have a distinct diurnal cycle due to morning/evening rush hours (Sailor and Lu, 2004). Generally, human metabolism contributes less than 5% of total anthropogenic flux in the U.S. (Sailor and Lu, 2004).

The urban surface is characterized by a preponderance of impervious surfaces, which reduce water storage capacity and surface moisture availability (Oke 1982). The evapotranspiration flux in urban regions is thus generally lower compared to vegetated surfaces, which may increase surface and air temperatures. On the other hand, vegetated surfaces within urban areas are frequently irrigated (e.g., lawns and parks) resulting in more water availability and higher latent heat fluxes than might be expected from natural vegetation. The presence and amount of vegetated or pervious surfaces can influence the magnitude of the heat island effect (Sailor 1995, Upmanis et al. 1998, Avissar 1996). Impervious surfaces also affect the hydrological cycle by reducing infiltration compared to rural areas, thereby converting more precipitation into surface runoff (Oke 1987, Bonan 2002).

The arrangement of large roughness elements (e.g., buildings, trees) in an urban region generally increases the frictional drag of the surface on the atmospheric winds and thus reduces the mean wind speed and turbulent mixing within the urban canopy compared to more open rural areas (Oke 1987). A notable exception to this may occur during periods of weak regional winds when warm urban air creates low-level rural-urban breezes. Lower within-canopy winds can reduce total turbulent heat transport from urban surfaces and increase their surface temperature. The synoptic wind speed is an important control on the urban heat island (Landsberg 1981). Higher winds may effectively remove heat faster than the urban fabric generates it.

The geographic location of urban areas and the characteristics of the surrounding rural area influence the urban climate. For instance, many tropical heat islands are smaller than expected based on population size. Where cities are surrounded by wet rural surfaces, slower cooling by these rural surfaces due to higher thermal admittance may reduce heat island magnitudes, especially in tropical climates (Oke et al. 1991). Local wind systems may impact urban climates as well. For example, coastal cities may experience cooling of urban temperatures when ocean surface temperatures are cooler than the land and winds blow onshore. Cold-air drainage from surrounding mountainous areas may reduce urban warming as well at certain times (Comrie 2000).

Urban regions have increased downward longwave radiation from the overlying atmosphere due to trapping and re-emission from polluted layers and/or from vertical advection of warm surface air above the city. Reduced incoming solar radiation due to reflection from atmospheric aerosols may compensate for this increase in longwave forcing. Note that in order to model these particular urban effects, the land model must also deliver biogeochemical fluxes (e.g., particulates, sulphur compounds, hydrocarbons, etc.) to the atmospheric model in addition to heat and moisture fluxes. The atmospheric model must then be able to diffuse or transport these trace species and determine their interaction with radiation and clouds. It has also been established that urban regions have effects on clouds and precipitation although the underlying mechanisms are still being debated. Climate modeling systems with detailed urban parameterizations may help to understand these mechanisms (Shepherd 2005).

As mentioned briefly in the previous section, many of the characteristics of the urban ecosystem discussed above contribute to one of the most striking effects of the urban environment on climate, the heat island effect. The present model is designed to represent the urban energy balance and provide insight into issues such as the urban heat island, its causes and potential mitigation strategies, as well as the effects of climate change on urban areas. When coupled to an atmospheric model, interactions between the urban surface and the atmosphere can be investigated.

Atmospheric Coupling and Model Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The atmospheric model within CCSM requires fluxes of sensible and latent heat and momentum between the surface and lowest atmospheric model level as well as emitted longwave and reflected shortwave radiation (:numref:`fig-atm-urban-coupling`). These must be provided at a time step that resolves the diurnal cycle. Over other types of land surfaces, the fluxes are determined by current parameterizations in CLM. An objective of this technical note is to describe a set of parameterizations that determines the fluxes from an urban surface. The vertical spatial domain of the urban model extends from the top of the urban canopy layer (UCL) down to the depth of zero vertical heat flux in the ground (Oke 1987). The current state of the atmosphere and downwelling fluxes (:numref:`table-atm-input`) at a given time step is used to force the urban model. The urban model provides fluxes that are area-averaged with other land cover (e.g., forests, cropland) if present within the grid cell. The area-averaged fluxes (:numref:`table-atm-output`) are used as lower boundary conditions by the atmospheric model.

Land surface heterogeneity in the Community Land Model (CLM) is represented as a nested subgrid hierarchy (:numref:`fig-clm-subgrid-hierarchy`) in which grid cells are composed of multiple landunits, snow/soil columns, and plant functional types (PFTs). Each grid cell can have a different number of landunits, each landunit can have a different number of columns, and each column can have multiple PFTs. The first subgrid level, the landunit, is intended to capture the broadest spatial patterns of subgrid heterogeneity. The model described here is designed to represent urban landunits. Further division of the urban surface into urban landuse classes such as, for example, city core, industrial/commercial, and suburban is possible by specifying these classes as individual landunits.

The representation of the urban landunit is based on the canyon concept of Oke (1987). In this approach, the considerable complexity of the urban surface is reduced to a single urban canyon consisting of a canyon floor of width :math:`W` bordered by two facing buildings of height :math:`H` (:numref:`fig-urban-canyon`). Although the canyon floor is intended to represent various surfaces such as roads, parking lots, sidewalks, and residential lawns, etc., for convenience we henceforth refer to the canyon floor as a road. The urban canyon consists of roof, sunlit and shaded wall, and pervious and impervious road, each of which are treated as columns within the landunit (:numref:`fig-clm-subgrid-hierarchy`). The impervious road is intended to represent surfaces that are impervious to water infiltration (e.g., roads, parking lots, sidewalks) while the pervious road is intended to represent surfaces such as residential lawns and parks which may have active hydrology.

The approach used here to represent pervious surfaces is different than many urban schemes designed for use within mesoscale and global models. Most urban schemes use a separate land surface model scheme to represent the effects of pervious surfaces on urban climate. For example, the urban surface in the mesoscale model Meso-NH is modeled using the TEB and ISBA (Interactions between Soil, Biosphere, and Atmosphere) schemes for urban and pervious (e.g., vegetated) surfaces, respectively (Lemonsu and Masson 2002). Fluxes from each scheme are combined according to their relative areas. A comparable approach could be implemented using the CLM scheme for vegetated surfaces; however, this presents several disadvantages for our application. First, the pervious surface would need to be assigned to an additional landunit and specially identified to distinguish it from the other vegetated landunit within the gridcell. Second, the pervious and urban landunits would then need to be aggregated according to their relative areas in a post-processing sense to estimate the composite urban effects. Third, in the Meso-NH approach, the pervious surface only interacts indirectly with the canyon air through its influence on the atmospheric model. Here, including the pervious surface within the urban canyon solves these difficulties. Thus, the pervious surface is an integral part of the urban system and interacts directly with UCL air properties such as temperature and specific humidity. Yet, implementation of a sophisticated scheme for the pervious surface, such as the vegetation scheme in CLM, within the urban canyon is problematic because of computational and data requirements. Here, we choose a simplified bulk parameterization scheme to represent latent heat flux from pervious urban surfaces (Chapter 3).

Note that the urban columns interact radiatively with one another through multiple exchanges of longwave and shortwave radiation (chapter 2). The heat and moisture fluxes from each surface interact with each other through a bulk air mass that represents air in the UCL for which specific humidity and temperature are predicted (chapter 3). We model the UCL plus the air above the roof (:numref:`fig-atm-urban-coupling`). This allows for mixing of above-roof air with canyon air.

.. figure:: image1.jpeg
   :width: 5.85417in
   :height: 4.67708in
   :name: fig-atm-urban-coupling

   Schematic of urban and atmospheric coupling. The urban model is forced by the atmospheric wind (:math:`u_{atm}`), temperature (:math:`T_{atm}`), specific humidity (:math:`q_{atm}`), precipitation (:math:`P_{atm}`), solar (:math:`S_{atm} \downarrow`) and longwave (:math:`L_{atm} \downarrow`) radiation at reference height :math:`z_{atm}^{'}`. Fluxes from the urban landunit to the atmosphere are turbulent sensible (:math:`H`) and latent heat (:math:`\lambda E`), momentum (:math:`\tau`), albedo (:math:`I \uparrow`), emitted longwave (:math:`L \uparrow`), and absorbed shortwave (:math:`\overrightarrow{S}`) radiation. Air temperature (:math:`T_{ac}`), specific humidity (:math:`q_{ac}`), and wind speed (:math:`u_{c}`) within the urban canopy layer are diagnosed by the urban model. :math:`H` is the average building height.

.. figure:: image2.jpeg
   :width: 5.98958in
   :height: 4.48958in
   :name: fig-clm-subgrid-hierarchy

   CLM subgrid hierarchy emphasizing the structure of urban landunits.

.. figure:: image3.jpeg
   :width: 5.86458in
   :height: 3.9375in
   :name: fig-urban-canyon

   The urban canyon.

.. list-table:: Atmospheric input to urban model
   :widths: auto
   :header-rows: 0
   :name: table-atm-input


   * - :sup:`1`\ Reference height
     - :math:`z_{atm}^{'}`
     - m
   * - Zonal wind at :math:`z_{atm}`
     - :math:`u_{atm}`
     - m s\ :sup:`-1`
   * - Meridional wind at :math:`z_{atm}`
     - :math:`v_{atm}`
     - m s\ :sup:`-1`
   * - Potential temperature
     - :math:`\overline{\theta_{atm}}`
     - K
   * - Specific humidity at :math:`z_{atm}`
     - :math:`q_{atm}`
     - kg kg\ :sup:`-1`
   * - Pressure at :math:`z_{atm}`
     - :math:`P_{atm}`
     - Pa
   * - Temperature at :math:`z_{atm}`
     - :math:`T_{atm}`
     - K
   * - Incident longwave radiation
     - :math:`L_{atm} \downarrow`
     - W m\ :sup:`-2`
   * - :sup:`2`\ Liquid precipitation
     - :math:`q_{rain}`
     - mm s\ :sup:`-1`
   * - :sup:`2`\ Solid precipitation
     - :math:`q_{sno}`
     - mm s\ :sup:`-1`
   * - Incident direct beam visible solar radiation
     - :math:`S_{atm} \downarrow_{vis}^{\mu}`
     - W m\ :sup:`-2`
   * - Incident direct beam near-infrared solar radiation
     - :math:`S_{atm} \downarrow_{nir}^{\mu}`
     - W m\ :sup:`-2`
   * - Incident diffuse visible solar radiation
     - :math:`S_{atm} \downarrow_{vis}`
     - W m\ :sup:`-2`
   * - Incident diffuse near-infrared solar radiation
     - :math:`S_{atm} \downarrow_{nir}`
     - W m\ :sup:`-2`
   * - :sup:`3`\ Carbon dioxide (CO\ :sub:`2`) concentration
     - :math:`c_{a}`
     - ppmv
   * - :sup:`3`\ Aerosol deposition rate
     - :math:`D_{sp}`
     - kg m\ :sup:`-2` s\ :sup:`-1`
   * - :sup:`3`\ Nitrogen deposition rate
     - :math:`NF_{ndep\_ sminn}`
     - g (N) m\ :sup:`-2` yr\ :sup:`-1`


..

   :sup:`1`\ The atmospheric reference height received from the atmospheric model :math:`z_{atm}^{'}` is assumed to be the height above the surface defined as the roughness length :math:`z_{0}` plus displacement height :math:`z_{d}`. Thus, the reference height used for flux computations (chapter 3) is :math:`z_{atm} = z_{atm}^{'} + z_{0} + z_{d}`. The reference heights for temperature, wind, and specific humidity (:math:`z_{atm,h}`, :math:`z_{atm,m}`, :math:`z_{atm,w}`) are required. These are set equal to\ :math:`z_{atm}`.

   :sup:`2`\ The CAM provides convective and large-scale liquid and solid precipitation, which are added to yield total liquid precipitation :math:`q_{rain}` and solid precipitation :math:`q_{sno}`.

   :sup:`3`\ These are provided by the atmospheric model but not used by the urban model.

Density of air (:math:`\rho_{atm}`) (kg m\ :sup:`-3`) is also required but is calculated directly from :math:`\rho_{atm} = \frac{P_{atm} - 0.378e_{atm}}{R_{da}T_{atm}}` where :math:`P_{atm}` is atmospheric pressure (Pa), :math:`e_{atm}` is atmospheric vapor pressure (Pa), :math:`R_{da}` is the gas constant for dry air (J kg\ :sup:`-1` K\ :sup:`-1`) (:numref:`table-physical-constants`), and :math:`T_{atm}` is the atmospheric temperature (K). The atmospheric vapor pressure :math:`e_{atm}` is derived from atmospheric specific humidity :math:`q_{atm}` (kg kg\ :sup:`-1`) as :math:`e_{atm} = \frac{q_{atm}P_{atm}}{0.622 + 0.378q_{atm}}`.

.. list-table:: Urban model output to atmospheric model
   :widths: auto
   :header-rows: 0
   :name: table-atm-output


   * - :sup:`1`\ Latent heat flux
     - :math:`\lambda E`
     - W m\ :sup:`-2`
   * - Sensible heat flux
     - :math:`H`
     - W m\ :sup:`-2`
   * - Water vapor flux
     - :math:`E`
     - mm s\ :sup:`-1`
   * - Zonal momentum flux
     - :math:`\tau_{x}`
     - kg m\ :sup:`-1` s\ :sup:`-2`
   * - Meridional momentum flux
     - :math:`\tau_{y}`
     - kg m\ :sup:`-1` s\ :sup:`-2`
   * - Emitted longwave radiation
     - :math:`L \uparrow`
     - W m\ :sup:`-2`
   * - Direct beam visible albedo
     - :math:`I \uparrow_{vis}^{\mu}`
     - -
   * - Direct beam near-infrared albedo
     - :math:`I \uparrow_{nir}^{\mu}`
     - -
   * - Diffuse visible albedo
     - :math:`I \uparrow_{vis}`
     - -
   * - Diffuse near-infrared albedo
     - :math:`I \uparrow_{nir}`
     - -
   * - Absorbed solar radiation
     - :math:`\overrightarrow{S}`
     - W m\ :sup:`-2`
   * - Radiative temperature
     - :math:`T_{rad}`
     - K
   * - Temperature at 2 meter height
     - :math:`T_{2m}`
     - K
   * - Specific humidity at 2 meter height
     - :math:`q_{2m}`
     - kg kg\ :sup:`-1`
   * - Snow water equivalent
     - :math:`W_{sno}`
     - m
   * - Aerodynamic resistance
     - :math:`r_{am}`
     - s m\ :sup:`-1`
   * - Friction velocity
     - :math:`u_{*}`
     - m s\ :sup:`-1`
   * - :sup:`2`\ Dust flux
     - :math:`F_{j}`
     - kg m\ :sup:`-2` s\ :sup:`-1`
   * - :sup:`2`\ Net ecosystem exchange
     - NEE
     - kgCO\ :sub:`2` m\ :sup:`-2` s\ :sup:`-1`


:sup:`1`\ :math:`\lambda` is either the latent heat of vaporization :math:`\lambda_{vap}` or latent heat of sublimation :math:`\lambda_{sub}` (J kg\ :sup:`-1`) (:numref:`table-physical-constants`) depending on the thermal state of surface water on the roof, pervious and impervious road.

:sup:`2`\ These are set to zero for urban areas.

Biogeophysical Processes
~~~~~~~~~~~~~~~~~~~~~~~~

Biogeophysical processes are simulated for each of the five urban columns and each column maintains its own prognostic variables (e.g., surface temperature). The processes simulated include:

-  Absorption and reflection of solar radiation (chapter 2)

-  Absorption, reflection, and emission of longwave radiation (chapter 2)

-  Momentum, sensible heat, and latent heat fluxes (chapter 3)

-  Anthropogenic heat fluxes to the canyon air due to waste heat from building heating/air conditioning (chapter 3). An example of parameterizing traffic heat fluxes is given in Oleson et al. (2008b), however, traffic heat fluxes are not currently included in the global implementation of the model.

-  Heat transfer in roofs, building walls, and the road including phase change (chapter 4)

-  Hydrology [roofs - storage of liquid and solid precipitation (ponding and dew), surface runoff; walls – hydrologically inactive; impervious road – storage of liquid and solid precipitation (ponding and dew), surface runoff; pervious road - infiltration, surface runoff, sub-surface drainage, redistribution of water within the column] (chapter 5).

Model Requirements
------------------

Initialization
~~~~~~~~~~~~~~

Initialization of the urban model (i.e., providing the model with initial temperature and moisture states) depends on the type of run (startup or restart) (see the CLM4 User's Guide). An initial run starts the model from either initial conditions that are set internally in the Fortran code (referred to as arbitrary initial conditions) or from an initial conditions dataset that enables the model to start from a spun up state (i.e., where the urban landunit is in equilibrium with the simulated climate). In restart runs, the model is continued from a previous simulation and initialized from a restart file that ensures that the output is bit-for-bit the same as if the previous simulation had not stopped. The fields that are required from the restart or initial conditions files can be obtained by examining the code. Arbitrary initial conditions are specified as follows.

All urban columns consist of fifteen layers to be consistent with CLM4. Generally, temperature calculations are done over all layers, :math:`N_{levgrnd} = 15`, while hydrology calculations for the pervious road are done over the top ten layers, :math:`N_{levsoi} = 10`, the bottom five layers being specified as bedrock. Pervious and impervious road are initialized with temperatures (surface :math:`T_{g}`, and layers :math:`T_{i}`, for layers :math:`i = 1,\ldots,N_{levgrnd}`) of 274 K. Roof, sunwall, and shadewall are initialized to 292K. This relatively high temperature is to avoid initialization shock from large space heating/air conditioning and waste heat fluxes. All surfaces are initialized with no snow (:math:`W_{sno} = 0`). Roof and impervious road are initialized with no ponded water, while the pervious road soil layers :math:`i = 1,\ldots,N_{levsoi}` are initialized with volumetric soil water content :math:`\theta_{i} = 0.3` mm\ :sup:`3` mm\ :sup:`-3` and layers :math:`i = N_{levsoi} + 1,\ldots,N_{levgrnd}` are initialized :math:`\theta_{i} = 0.0` mm\ :sup:`3` mm\ :sup:`-3`. The soil liquid water and ice contents are initialized as :math:`w_{liq,i} = \Delta z_{i}\rho_{liq}\theta_{i}` and :math:`w_{ice,i} = 0`, where :math:`\rho_{liq}` is the density liquid water (kg m\ :sup:`-3`) (:numref:`table-physical-constants`). The pervious road is initialized with water stored in the unconfined aquifer and unsaturated soil :math:`W_{a} = W_{t} = 4800` mm and water table depth :math:`z_{\nabla} = 4.8` m.

Surface Data
~~~~~~~~~~~~

Required input data for urban landunits are listed in :numref:`table-input-data` This data is provided by the surface dataset at the required spatial resolution (see the CLM4 User's Guide). Present day global urban extent and urban properties were developed by Jackson et al. (2010). Urban extent, defined for four classes [tall building district (TBD), and high, medium, and low density (HD, MD, LD)] was derived from LandScan 2004, a population density dataset derived from census data, nighttime lights satellite observations, road proximity, and slope (Dobson et al., 2000). The urban extent data is aggregated from the original 1 km resolution to a 0.5° by 0.5° global grid. For this particular implementation, only the sum of the TBD, HD, and MD classes are used to define urban extent as the LD class is highly rural and likely better modeled as a vegetated surface.

For each of 33 distinct regions across the globe, thermal (e.g., heat capacity and thermal conductivity), radiative (e.g., albedo and emissivity) and morphological (e.g., height to width ratio, roof fraction, average building height, and pervious fraction of the canyon floor) properties of roof/wall/road are provided by Jackson et al. (2010) for each of the four density classes. Building interior minimum and maximum temperatures are prescribed based on climate and socioeconomic considerations. Urban parameters are determined for the 0.5° by 0.5° global grid based on the dominant density class by area. This prevents potentially unrealistic parameter values that may result if the density classes are averaged. As a result, the current global representation of urban is almost exclusively medium density. Future implementations of the model could represent each of the density classes as a separate landunit. The surface dataset creation routines (see CLM4 User's Guide) aggregate the data to the desired resolution. It is surmised that the MODIS-based vegetation dataset used in CLM4 classifies built areas as bare soil, thus the urban extent preferentially replaces bare soil when it exists within the grid cell. A very small minimum threshold of 0.1% of the grid cell by area is used to resolve urban areas. An elevation threshold of 2200 m is used to eliminate urban areas where the grid cell surface elevation is significantly higher than the elevation the cities are actually at because of the coarse spatial resolution of the model. This prevents overestimates of anthropogenic heating in winter due to unrealistically cold temperatures.

.. list-table:: Input data required for the urban model
   :widths: auto
   :header-rows: 0
   :name: table-input-data


   * - Data
     - Symbol
     - Units
   * - Percent urban
     - -
     - %
   * - Canyon height to width ratio
     - :math:`\frac{H}{W}`
     - -
   * - Roof fraction
     - :math:`W_{roof}`
     - -
   * - :sup:`1`\ Pervious road fraction
     - :math:`f_{prvrd}`
     - -
   * - Emissivity of roof
     - :math:`\varepsilon_{roof}`
     - -
   * - Emissivity of impervious road
     - :math:`\varepsilon_{imprvrd}`
     - -
   * - Emissivity of pervious road
     - :math:`\varepsilon_{prvrd}`
     - -
   * - Emissivity of sunlit and shaded walls
     - :math:`\varepsilon_{wall}`
     - -
   * - Building height
     - :math:`H`
     - m
   * - Roof albedo – visible direct
     - :math:`\alpha_{roof,vis}^{\mu}`
     - -
   * - Roof albedo – visible diffuse
     - :math:`\alpha_{roof,vis}`
     - -
   * - Roof albedo – near-infrared direct
     - :math:`\alpha_{roof,nir}^{\mu}`
     - -
   * - Roof albedo – near-infrared diffuse
     - :math:`\alpha_{roof,nir}`
     - -
   * - Wall albedo – visible direct
     - :math:`\alpha_{walls,vis}^{\mu}`
     - -
   * - Wall albedo – visible diffuse
     - :math:`\alpha_{walls,vis}`
     - -
   * - Wall albedo – near-infrared direct
     - :math:`\alpha_{walls,nir}^{\mu}`
     - -
   * - Wall albedo – near-infrared diffuse
     - :math:`\alpha_{walls,nir}`
     - -
   * - Impervious road albedo – visible direct
     - :math:`\alpha_{imprvrd,vis}^{\mu}`
     - -
   * - Impervious road albedo – visible diffuse
     - :math:`\alpha_{imprvrd,vis}`
     - -
   * - Impervious road albedo – near-infrared direct
     - :math:`\alpha_{imprvrd,nir}^{\mu}`
     - -
   * - Impervious road albedo – near-infrared diffuse
     - :math:`\alpha_{imprvrd,nir}`
     - -
   * - Pervious road albedo – visible direct
     - :math:`\alpha_{prvrd,vis}^{\mu}`
     - -
   * - Pervious road albedo – visible diffuse
     - :math:`\alpha_{prvrd,vis}`
     - -
   * - Pervious road albedo – near-infrared direct
     - :math:`\alpha_{prvrd,nir}^{\mu}`
     - -
   * - Pervious road albedo – near-infrared diffuse
     - :math:`\alpha_{prvrd,nir}`
     - -
   * - Roof thermal conductivity
     - :math:`\lambda_{roof,i}`
     - W m\ :sup:`-1` K\ :sup:`-1`
   * - Wall thermal conductivity
     - :math:`\lambda_{wall,i}`
     - W m\ :sup:`-1` K\ :sup:`-1`
   * - :sup:`2`\ Impervious road thermal conductivity
     - :math:`\lambda_{imprvrd,i}`
     - W m\ :sup:`-1` K\ :sup:`-1`
   * - :sup:`3`\ Pervious road thermal conductivity
     - :math:`\lambda_{prvrd,i}`
     - W m\ :sup:`-1` K\ :sup:`-1`
   * - Roof volumetric heat capacity
     - :math:`c_{roof,i}`
     - J m\ :sup:`-3` K\ :sup:`-1`
   * - Wall volumetric heat capacity
     - :math:`c_{wall,i}`
     - J m\ :sup:`-3` K\ :sup:`-1`
   * - :sup:`2`\ Impervious road volumetric heat capacity
     - :math:`c_{imprvrd,i}`
     - J m\ :sup:`-3` K\ :sup:`-1`
   * - :sup:`3`\ Pervious road volumetric heat capacity
     - :math:`c_{imprvrd,i}`
     - J m\ :sup:`-3` K\ :sup:`-1`
   * - Maximum interior building temperature
     - :math:`T_{iB,\max}`
     - K
   * - Minimum interior building temperature
     - :math:`T_{iB,\min}`
     - K
   * - Height of wind source in canyon
     - :math:`H_{w}`
     - m
   * - Number of impervious road layers
     - :math:`N_{imprvrd}`
     - -
   * - Wall thickness
     - :math:`\Delta z_{wall}`
     - m
   * - Roof thickness
     - :math:`\Delta z_{roof}`
     - m
   * - :sup:`4`\ Percent sand, percent clay of pervious road (soil)
     - :math:`\% sand,\% clay`
     - %
   * - Grid cell latitude and longitude
     - :math:`\varphi,\theta`
     - degrees


:sup:`1`\ This fraction is relative to the canyon floor.

:sup:`2`\ Required for layers :math:`i = 1,N_{imprvrd}`, derived from grid cell soil texture for other layers (section 4.3).

:sup:`3`\ Derived from grid cell soil texture (:math:`\% sand,\% clay`) (section 4.3).

:sup:`4`\ Obtained from grid cell soil texture (:math:`\% sand,\% clay`).

Physical Constants
~~~~~~~~~~~~~~~~~~

Physical constants, shared by all of the components in the CCSM, are presented in :numref:`table-physical-constants` Not all constants are necessarily used by the urban model.

.. list-table:: Physical constants
   :widths: auto
   :header-rows: 0
   :name: table-physical-constants


   * - Pi
     - :math:`\pi`
     - 3.14159265358979323846
     - -
   * - Acceleration of gravity
     - :math:`g`
     - 9.80616
     - m s\ :sup:`-2`
   * - Standard pressure
     - :math:`P_{std}`
     - 101325
     - Pa
   * - Stefan-Boltzmann constant
     - :math:`\sigma`
     - 5.67\ :math:`\times 10^{- 8}`
     - W m\ :sup:`-2` K\ :sup:`-4`
   * - Boltzmann constant
     - :math:`\kappa`
     - 1.38065\ :math:`\times 10^{- 23}`
     - J K\ :sup:`-1` molecule\ :sup:`-1`
   * - Avogadro's number
     - :math:`N_{A}`
     - 6.02214\ :math:`\times 10^{26}`
     - molecule kmol\ :sup:`-1`
   * - Universal gas constant
     - :math:`R_{gas}`
     - :math:`N_{A}\kappa`
     - J K\ :sup:`-1` kmol\ :sup:`-1`
   * - Molecular weight of dry air
     - :math:`MW_{da}`
     - 28.966
     - kg kmol\ :sup:`-1`
   * - Dry air gas constant
     - :math:`R_{da}`
     - :math:`\frac{R_{gas}}{MW_{da}}`
     - J K\ :sup:`-1` kg\ :sup:`-1`
   * - Molecular weight of water vapor
     - :math:`MW_{wv}`
     - 18.016
     - kg kmol\ :sup:`-1`
   * - Water vapor gas constant
     - :math:`R_{wv}`
     - :math:`\frac{R_{gas}}{MW_{wv}}`
     - J K\ :sup:`-1` kg\ :sup:`-1`
   * - Von Karman constant
     - :math:`k`
     - 0.4
     - -
   * - Freezing temperature of fresh water
     - :math:`T_{f}`
     - 273.15
     - K
   * - Density of liquid water
     - :math:`\rho_{liq}`
     - 1000
     - kg m\ :sup:`-3`
   * - Density of ice
     - :math:`\rho_{ice}`
     - 917
     - kg m\ :sup:`-3`
   * - Specific heat capacity of dry air
     - :math:`C_{p}`
     - 1.00464\ :math:`\times 10^{3}`
     - J kg\ :sup:`-1` K\ :sup:`-1`
   * - Specific heat capacity of water
     - :math:`C_{liq}`
     - 4.188\ :math:`\times 10^{3}`
     - J kg\ :sup:`-1` K\ :sup:`-1`
   * - Specific heat capacity of ice
     - :math:`C_{ice}`
     - 2.11727\ :math:`\times 10^{3}`
     - J kg\ :sup:`-1` K\ :sup:`-1`
   * - Latent heat of vaporization
     - :math:`\lambda_{vap}`
     - 2.501\ :math:`\times 10^{6}`
     - J kg\ :sup:`-1`
   * - Latent heat of fusion
     - :math:`L_{f}`
     - 3.337\ :math:`\times 10^{5}`
     - J kg\ :sup:`-1`
   * - Latent heat of sublimation
     - :math:`\lambda_{sub}`
     - :math:`\lambda_{vap} + L_{f}`
     - J kg\ :sup:`-1`
   * - :sup:`1`\ Thermal conductivity of water
     - :math:`\lambda_{liq}`
     - 0.6
     - W m\ :sup:`-1` K\ :sup:`-1`
   * - :sup:`1`\ Thermal conductivity of ice
     - :math:`\lambda_{ice}`
     - 2.29
     - W m\ :sup:`-1` K\ :sup:`-1`
   * - :sup:`1`\ Thermal conductivity of air
     - :math:`\lambda_{air}`
     - 0.023
     - W m\ :sup:`-1` K\ :sup:`-1`
   * - Radius of the earth
     - :math:`R_{e}`
     - 6.37122\ :math:`\times 10^{6}`
     - m


:sup:`1`\ Not shared by other components of the coupled modeling system.

Albedos and Radiative Fluxes
============================

The effects of geometry on the radiation balance of urban surfaces are a key driver of urban-rural energy balance differences (Oke et al. 1991). Shadowing of urban surfaces affects the incident radiation and thus temperature. Similar to vegetated surfaces, multiple reflections of radiation between urban surfaces must be accounted for (Harman et al. 2004). The net solar radiation and net longwave radiation, the net of which is the net radiation, are needed for each urban surface to drive turbulent and ground heat fluxes. The atmospheric model also requires radiative fluxes and albedo from the urban landunit, which are appropriately averaged with other landunits within the gridcell. The urban canyon unit is used to represent these radiative processes. Several simplifying assumptions are made. The effects of absorption, emission, and scattering of radiation by the canyon air are neglected and surfaces are assumed to be isotropic.

Albedo
------

The albedo of each urban surface is a weighted combination of snow-free "ground" albedo and snow albedo. Only roof and road surfaces are affected by snow. The direct beam :math:`\alpha_{u,\Lambda}^{\mu}` and diffuse :math:`\alpha_{u,\Lambda}` albedos (where :math:`u` denotes roof, impervious or pervious road) are

.. math::
   :label: eq-0103

   \alpha_{u,\Lambda}^{\mu} = \alpha_{g,\Lambda}^{\mu}\left( 1 - f_{u,sno} \right) + \alpha_{sno,\Lambda}^{\mu}f_{u,sno}

.. math::
   :label: eq-0104

   \alpha_{u,\Lambda} = \alpha_{g,\Lambda}\left( 1 - f_{u,sno} \right) + \alpha_{sno,\Lambda}f_{u,sno}

where :math:`f_{u,sno}` is the fraction of the urban surface covered with snow which is calculated from (Bonan 1996)

:math:`f_{u,sno} = \frac{z_{u,sno}}{0.05} \leq 1` .

The direct and diffuse "ground" albedos, :math:`\alpha_{g,\Lambda}^{\mu}` and :math:`\alpha_{g,\Lambda}`, where :math:`\Lambda` denotes either the visible (VIS) or near-infrared (NIR) waveband, are provided by the surface dataset (:numref:`table-input-data`), and :math:`z_{u,sno}` is the depth of snow (m) (section 5.1). An estimate of snow albedo is made based on the parameterization of Marshall (1989) in which albedo depends on solar zenith angle, grain size, and soot content (e.g., as adopted by the Land Surface Model (LSM) (Bonan 1996)). Here, however, several simplifying assumptions are made due to uncertainties in how to apply such a parameterization to urban surfaces. A snow grain radius of 100 :math:`\mu m` (new powder snow, aged a few days) and a soot mass fraction of 1.5\ :math:`\times 10^{- 5}` (arrived at by noting that the LSM global soot mass fraction is 5\ :math:`\times 10^{- 6}` and Chylek et al. (1987) observed that soot concentrations in urban snowpacks averaged three times the concentration in rural snowpacks) are assumed. Direct and diffuse albedos are assumed to be equal. This yields :math:`\alpha_{sno,VIS}^{\mu} = \alpha_{sno,VIS} = 0.66` and :math:`\alpha_{sno,NIR}^{\mu} = \alpha_{sno,NIR} = 0.56` which fall about in the middle of the range given by Oke (1987).

Incident direct solar radiation
-------------------------------

Unlike the horizontal roof surface, the direct beam solar radiation received by the walls and the road must be adjusted for orientation and shadowing. The analytical solution given below follows Masson (2000). First, let :math:`\theta` be the angle between the sun direction and the along-canyon axis and consider the case where the along-canyon axis is perpendicular to the sun direction (:math:`\theta = \frac{\pi}{2}`). In this case, as shown in :numref:`fig-solar-elevation-view`, if the solar zenith angle :math:`\mu` is greater than the critical solar zenith angle :math:`\mu_{0}` (:math:`\mu_{0} = \tan^{- 1}\left( \frac{W}{H} \right)`), the road is in full shade, and the sunlit wall is in partial sun. Conversely, if :math:`\mu` is less than :math:`{\overrightarrow{L}}_{uc} - \left( L_{uc} \uparrow - L_{atm} \downarrow \right) = 0`, the road is in partial sun and the sunlit wall is in full sun. Note that, radiatively, the pervious and impervious road are treated the same, although their albedos are specified separately and may differ (:numref:`table-input-data`).

.. figure:: image4.jpeg
   :width: 5.40625in
   :height: 6.78125in
   :name: fig-solar-elevation-view

   Elevation (side) view of direct beam solar radiation incident on urban canyon surfaces for solar zenith angle :math:`\mu > \mu_{0}` (top) and :math:`\mu \leq \mu_{0}` (bottom). :math:`S_{atm} \downarrow_{\Lambda}^{\mu}` is the direct beam incident solar radiation incident on a horizontal surface from the atmosphere. The along-canyon axis is assumed to be perpendicular to the sun direction.

If the direct beam solar radiation received by a horizontal surface (i.e., as received by the roof) is :math:`S_{atm} \downarrow_{\Lambda}^{\mu}`, then the solar radiation on the wall in full illumination (:math:`\mu \leq \mu_{0}`) is :math:`\frac{\left( S_{atm} \downarrow_{\Lambda}^{\mu}\cos i \right)}{\cos\mu}` where :math:`i` is the incidence angle (:numref:`fig-solar-elevation-view`). Since :math:`\cos i = \cos(90 - \mu) = \sin\mu`, the solar radiation on the sunlit wall is

.. math::
   :label: eq-0105

   S_{sunwall} \downarrow_{\Lambda}^{\mu}\left( \theta = \frac{\pi}{2} \right) = \tan(\mu)S_{atm} \downarrow_{\Lambda}^{\mu} \mu \leq \mu_{0}

Note that this is twice the radiation received by the wall in Masson (2000) because here we force the other (shaded) wall to receive no solar radiation (:math:`S_{shdwall} \downarrow_{\Lambda}^{\mu} = 0`). In the case of :math:`\mu > \mu_{0}`, the illuminated fraction is :math:`\frac{(H - y)}{H}` and :math:`S_{sunwall} \downarrow_{\Lambda}^{\mu} = \left\lbrack \frac{(H - y)}{H} \right\rbrack\tan\mu S_{atm} \downarrow_{\Lambda}^{\mu}`. Since :math:`\tan\mu = \frac{W}{(H - y)}` this simplifies to

.. math::
   :label: eq-0106

   S_{sunwall} \downarrow_{\Lambda}^{\mu}\left( \theta = \frac{\pi}{2} \right) = \frac{W}{H}S_{atm} \downarrow_{\Lambda}^{\mu} \mu > \mu_{0}

Since the road is a horizontal surface, :math:`S_{road} \downarrow_{\Lambda}^{\mu} = \left\lbrack \frac{(W - x)}{W} \right\rbrack S_{atm} \downarrow_{\Lambda}^{\mu}` for :math:`\mu \leq \mu_{0}`. Since :math:`x = H\tan\mu`, the direct solar radiation incident on the road (pervious and impervious) is



.. math::
   :label: eq-0001

   S_{road} \downarrow_{\Lambda}^{\mu}\left( \theta = \frac{\pi}{2} \right) = \left\{ \begin{aligned}
   & 0 \mu > \mu_{0} \\
   & \left( 1 - \frac{H}{W}\tan\mu \right)S_{atm} \downarrow_{\Lambda}^{\mu} \mu \leq \mu_{0}
   \end{aligned} \right\}


Equations :eq:`eq-0105` and :eq:`eq-0106` for the walls and equation :eq:`eq-0001` for the road can now be expanded to account for any canyon orientation (:math:`0 \leq \theta \leq \frac{\pi}{2}`). If :math:`\theta` is the angle between the sun direction and the along-canyon axis (:numref:`fig-solar-plan-view`), then the expression for the incidence angle is now :math:`\cos i = \sin\mu\sin\theta` and equation :eq:`eq-0105` becomes

.. math::
   :label: eq-0107

   S_{sunwall} \downarrow_{\Lambda}^{\mu}(\theta) = \sin\theta\tan\mu S_{atm} \downarrow_{\Lambda}^{\mu} \mu \leq \mu_{0}

.. figure:: image5.jpeg
   :width: 6in
   :height: 3.98958in
   :name: fig-solar-plan-view

   Plan view of direct beam solar radiation incident on urban canyon surfaces. :math:`S_{atm} \downarrow_{\Lambda}^{\mu}` is the direct beam incident solar radiation incident on a horizontal surface from the atmosphere. :math:`\theta` is the angle between the along-canyon axis and the sun direction.

For the case of :math:`\mu > \mu_{0}`, :math:`S_{sunwall} \downarrow_{\Lambda}^{\mu}(\theta) = \left\lbrack \frac{(H - y)}{H} \right\rbrack\sin\theta\tan\mu S_{atm} \downarrow_{\Lambda}^{\mu}`. However, now :math:`\tan\mu = \frac{\left( \frac{W}{\sin\theta} \right)}{(H - y)}` and thus

.. math::
   :label: eq-0108

   S_{sunwall} \downarrow_{\Lambda}^{\mu}(\theta) = \frac{W}{H}S_{atm} \downarrow_{\Lambda}^{\mu} \mu > \mu_{0}

Similarly, for the road (:math:`\mu \leq \mu_{0}`), :math:`S_{road} \downarrow_{\Lambda}^{\mu}(\theta) = \left\lbrack \frac{\left( \frac{W}{\sin\theta} - x \right)}{\left( \frac{W}{\sin\theta} \right)} \right\rbrack S_{atm} \downarrow_{\Lambda}^{\mu}` with :math:`x = H\tan\mu` simplifies to



.. math::
   :label: eq-0002

   S_{road} \downarrow_{\Lambda}^{\mu}(\theta) = \left\{ \begin{aligned}
   & 0 \mu > \mu_{0} \\
   & \left( 1 - \frac{H}{W}\sin\theta\tan\mu \right)S_{atm} \downarrow_{\Lambda}^{\mu} \mu \leq \mu_{0}
   \end{aligned} \right\}


Note that the critical solar zenith angle is now

.. math::
   :label: eq-0109

   \mu_{0} = \tan^{- 1}\left( \frac{\frac{W}{\sin\theta}}{H} \right)

Equations :eq:`eq-0107`, :eq:`eq-0108`, and :eq:`eq-0002` are integrated over all canyon orientations (:math:`0 \leq \theta \leq \frac{\pi}{2}`). The integration is done in two parts, first from :math:`\theta = 0` to :math:`\theta = \theta_{0}`, and second from :math:`\theta = \theta_{0}` to :math:`\theta = \frac{\pi}{2}`, where :math:`\theta_{0}` is the critical canyon orientation for which the road is no longer illuminated. This can be derived from Equation :eq:`eq-0109` and is

.. math::
   :label: eq-0110

   \theta_{0} = \sin^{- 1}\left\lbrack \min\left( \frac{W}{H\tan\mu},1 \right) \right\rbrack

The integrations thus are

.. math::
   :label: eq-0111

   S_{sunwall} \downarrow_{\Lambda}^{\mu} = \frac{4}{2\pi}\int_{0}^{\theta_{0}}{\sin\theta\tan\mu S_{atm} \downarrow_{\Lambda}^{\mu}d\theta} + \frac{4}{2\pi}\int_{\theta_{0}}^{\frac{\pi}{2}}{\frac{W}{H}S_{atm} \downarrow_{\Lambda}^{\mu}d\theta}

and

.. math::
   :label: eq-0112

   S_{road} \downarrow_{\Lambda}^{\mu} = \frac{4}{2\pi}\int_{0}^{\theta_{0}}{\left( 1 - \frac{H}{W}\sin\theta\tan\mu \right)S_{atm} \downarrow_{\Lambda}^{\mu}d\theta}

The direct beam solar radiation incident on the roof, walls and road is therefore

:math:`S_{roof} \downarrow_{\Lambda}^{\mu} = S_{atm} \downarrow_{\Lambda}^{\mu}`,

:math:`S_{shdwall} \downarrow_{\Lambda}^{\mu} = 0`,

:math:`S_{sunwall} \downarrow_{\Lambda}^{\mu} = 2S_{atm} \downarrow_{\Lambda}^{\mu}\left\lbrack \frac{W}{H}\left( \frac{1}{2} - \frac{\theta_{0}}{\pi} \right) + \frac{1}{\pi}\tan\mu\left( 1 - \cos\theta_{0} \right) \right\rbrack`,

.. math::
   :label: eq-0113

   S_{road} \downarrow_{\Lambda}^{\mu} = S_{imprvrd} \downarrow_{\Lambda}^{\mu} = S_{prvrd} \downarrow_{\Lambda}^{\mu} = S_{atm} \downarrow_{\Lambda}^{\mu}\left\lbrack \frac{2\theta_{0}}{\pi} - \frac{2}{\pi}\frac{H}{W}\tan\mu\left( 1 - \cos\theta_{0} \right) \right\rbrack

The direct incident solar radiation conserves energy as



.. math::
   :label: eq-0003

   S_{atm} \downarrow_{\Lambda}^{\mu} = f_{roof}S_{roof} \downarrow_{\Lambda}^{\mu} + \left( 1 - f_{roof} \right)
   \left\lbrack S_{imprvrd} \downarrow_{\Lambda}^{\mu}\left( 1 - f_{prvrd} \right) + S_{prvrd} \downarrow_{\Lambda}^{\mu}f_{prvrd} + \frac{H}{W}\left( S_{sunwall} \downarrow_{\Lambda}^{\mu} + S_{shdwall} \downarrow_{\Lambda}^{\mu} \right) \right\rbrack


Note that the factor :math:`\frac{H}{W}` for the sunlit wall and shaded wall converts the flux from watts per meter squared of wall area to watts per meter squared of ground area.

View factors
------------

The interaction of diffuse radiation (i.e., longwave and scattered solar radiation) between urban surfaces depends on angle (view) factors, i.e., the fraction of diffusely distributed energy leaving one "surface" (e.g., sky) that arrives at another surface (e.g., wall) (Sparrow and Cess 1978). If :math:`E_{ij}` is the diffuse radiative flux density on surface :math:`j` that originated from surface :math:`i` and :math:`E_{i}` is the radiative flux from surface :math:`i`, then

.. math::
   :label: eq-0114

   E_{ij} = F_{ij}E_{i}

where :math:`F_{ij}` is the view factor. The view factors depend only on the geometrical configurations of the involved surfaces. A table of view factors for various configurations is provided in Appendix A of Sparrow and Cess (1978). For instance, the view factor for the radiation from the wall to the sky can be derived from configuration nine of Appendix A. If :math:`dA_{1}` is an infinitesimal element on surface 1 (i.e., wall) and :math:`A_{2}` is a finite surface (i.e., sky) (:numref:`fig-view-factor-schematic`), then the angle factor :math:`F_{dA_{1} - A_{2}}` for diffuse radiation leaving element :math:`dA_{1}` and arriving at :math:`A_{2}` is

.. math::
   :label: eq-0115

   F_{dA_{1} - A_{2}} = \frac{1}{2\pi}\left( \tan^{- 1}\frac{1}{Y} - AY\tan^{- 1}A \right)

where :math:`A = \frac{1}{\sqrt{X^{2} + Y^{2}}}`, :math:`X = \frac{a}{b}`, and :math:`Y = c/b`. Following Sakakibara (1996) and Kusaka et al. (2001), for an infinitely long canyon, :math:`b = \infty`, :math:`a = W`, and so the wall-sky view factor at distance :math:`c` from a point on the wall to the canyon top is

.. math::
   :label: eq-0116

   \Psi_{wall - sky|c} = \frac{1}{2}\left( 1 - \frac{c}{\sqrt{c^{2} + W^{2}}} \right)

The total wall-sky view factor can be found by integrating the above equation over the height of the wall as

.. math::
   :label: eq-0117

   \Psi_{wall - sky} = \frac{1}{H}\int_{c = 0}^{c = H}{\frac{1}{2}\left( 1 - \frac{c}{\sqrt{c^{2} + W^{2}}} \right)dc} = \frac{\frac{1}{2}\left( \frac{H}{W} + 1 - \sqrt{1 + \left( \frac{H}{W} \right)^{2}} \right)}{\frac{H}{W}}

By the reciprocity rule (:math:`A_{1}F_{A_{1} - A_{2}} = A_{2}F_{A_{2} - A_{1}}`) (Sparrow and Cess 1978), the sky-wall view factor is

.. math::
   :label: eq-0118

   \Psi_{sky - wall} = \frac{H}{W}\Psi_{wall - sky}

When applied to equation :eq:`eq-0114` , :math:`\Psi_{sky - wall}` will yield a flux density to the wall in terms of per unit sky area. In the radiation computations detailed below, the diffuse fluxes for the walls are solved in terms of per unit wall area. Dividing equation :eq:`eq-0118` by the height to width ratio converts the view factor to per unit wall area. Thus,

.. math::
   :label: eq-0119

   \Psi_{sky - wall} = \frac{\frac{1}{2}\left( \frac{H}{W} + 1 - \sqrt{1 + \left( \frac{H}{W} \right)^{2}} \right)}{\frac{H}{W}}

Similarly, the view factor for radiation from the sky to the road and from road to sky can be solved and is

.. math::
   :label: eq-0120

   \Psi_{sky - road} = \frac{W}{W}\Psi_{road - sky} = \Psi_{road - sky} = \sqrt{1 + \left( \frac{H}{W} \right)^{2}} - \frac{H}{W}

By symmetry,

:math:`\Psi_{wall - road} = \Psi_{wall - sky}`,

and the other view factors can be deduced from conservation of energy as

:math:`\Psi_{road - wall} = \frac{1}{2}\left( 1 - \Psi_{road - sky} \right)`,

.. math::
   :label: eq-0121

   \Psi_{wall - wall} = 1 - \Psi_{wall - sky} - \Psi_{wall - road}

The view factors are presented graphically in :numref:`fig-view-factors-hw-ratio` Note that the view factors for radiation from the walls to the other surfaces sum to one (:math:`\Psi_{wall - wall} + \Psi_{wall - road} + \Psi_{wall - sky} = 1`). Similarly, the view factors for radiation from the road to the other surfaces also sum to one (:math:`\Psi_{road - wall} + \Psi_{road - wall} + \Psi_{road - sky} = 1`). As Harman et al. (2004) notes, at low height to width ratios, the road-sky view factor is close to one, the wall-wall view factor is close to zero, and the wall sky view factor is close to one half. However, at these low height to width ratios, the wall area is small compared to the road or sky area, indicating that most of the radiative exchange occurs between the road and sky, as it would for a flat surface. At height to width ratios greater than one, most of the radiative interactions take place between the two walls and the wall and the road. These view factors are consistent with those given by both Masson (2000) and Harman et al. (2004).

.. figure:: image6.jpeg
   :width: 6in
   :height: 4.79167in
   :name: fig-view-factor-schematic

   Schematic representation of angle (view) factor between infinitesimal element :math:`dA_{1}` (e.g., a point on the wall) and finite surface :math:`A_{2}` (e.g., the sky) (after Sparrow and Cess (1978)).

.. figure:: image7.png
   :width: 5.98958in
   :height: 5.625in
   :name: fig-view-factors-hw-ratio

   View factors as a function of canyon height to width ratio. :math:`\Psi_{road - sky}` is the fraction of radiation reaching the sky from the road, :math:`\Psi_{road - wall}` is the fraction of radiation reaching the wall from the road, :math:`\Psi_{wall - sky}` is the fraction of radiation reaching the sky from the wall, :math:`\Psi_{wall - road}` is the fraction of radiation reaching the road from the wall, and :math:`\Psi_{wall - wall}` is the fraction of radiation reaching the wall from the opposite wall.

Incident diffuse solar radiation
--------------------------------

The two view factors needed to compute the incident diffuse solar radiation are :math:`\Psi_{sky - road}` (equation :eq:`eq-0119` ) and :math:`\Psi_{sky - wall}` (equation :eq:`eq-0120` ). The diffuse solar radiation incident on roof, walls and road is then

:math:`S_{roof} \downarrow_{\Lambda} = S_{atm} \downarrow_{\Lambda}`,

:math:`S_{imprvrd} \downarrow_{\Lambda} = S_{prvrd} \downarrow_{\Lambda} = S_{atm} \downarrow_{\Lambda}\Psi_{sky - road}`,

:math:`S_{shdwall} \downarrow_{\Lambda} = S_{atm} \downarrow_{\Lambda}\Psi_{sky - wall}`,

.. math::
   :label: eq-0122

   S_{sunwall} \downarrow_{\Lambda} = S_{atm} \downarrow_{\Lambda}\Psi_{sky - wall}

The diffuse incident solar radiation conserves energy as



.. math::
   :label: eq-0004

   S_{atm} \downarrow_{\Lambda} = f_{roof}S_{roof} \downarrow_{\Lambda} + \left( 1 - f_{roof} \right)
   \left\lbrack S_{imprvrd} \downarrow_{\Lambda}\left( 1 - f_{prvrd} \right) + S_{prvrd} \downarrow_{\Lambda}f_{prvrd} + \frac{H}{W}\left( S_{sunwall} \downarrow_{\Lambda} + S_{shdwall} \downarrow_{\Lambda} \right) \right\rbrack


2.5 Absorbed and reflected solar radiation
------------------------------------------

The direct and diffuse net (absorbed) and reflected solar radiation for the roof is

.. math::
   :label: eq-0123

   {\overrightarrow{S}}_{roof,\Lambda}^{\mu} = S_{roof} \downarrow_{\Lambda}^{\mu}\left( 1 - \alpha_{roof,\Lambda}^{\mu} \right)

.. math::
   :label: eq-0124

   {\overrightarrow{S}}_{roof,\Lambda} = S_{roof} \downarrow_{\Lambda}\left( 1 - \alpha_{roof,\Lambda} \right)

.. math::
   :label: eq-0125

   S_{roof} \uparrow_{\Lambda}^{\mu} = S_{roof} \downarrow_{\Lambda}^{\mu}\left( \alpha_{roof,\Lambda}^{\mu} \right)

.. math::
   :label: eq-0126

   S_{roof} \uparrow_{\Lambda} = S_{roof} \downarrow_{\Lambda}\left( \alpha_{roof,\Lambda} \right)

The net (absorbed) and reflected solar radiation for walls and road and the reflected solar radiation to the sky are determined numerically by allowing for multiple reflections until a convergence criteria is met to ensure radiation is conserved. The reflected radiation from each urban surface is absorbed and re-reflected by the other urban surfaces. For example, the radiation scattered from the sunlit wall to the road, the shaded wall, and the sky depends on the view factors :math:`\Psi_{wall - road}`, :math:`\Psi_{wall - wall}`, and :math:`\Psi_{wall - sky}`, respectively (:numref:`fig-view-factors-hw-ratio`). The multiple reflections are accounted for in five steps:

1. Determine the initial absorption and reflection by each urban surface and distribute this radiation to the sky, road, and walls according to view factors.

2. Determine the amount of radiation absorbed and reflected by each urban surface after the initial reflection. The solar radiation reflected from the walls to the road is projected to road area by multiplying by the height to width ratio and the solar radiation reflected from the road to the walls is projected to wall area by dividing by the height to width ratio.

3. The absorbed radiation for the :math:`i^{th}` reflection is added to the total absorbed by each urban surface.

4. The reflected solar radiation for the :math:`i^{th}` reflection is distributed to the sky, road, and walls according to view factors.

5. The reflected solar radiation to the sky for the :math:`i^{th}` reflection is added to the total reflected solar radiation.

Steps 2-5 are repeated until a convergence criterion (absorbed radiation per unit incoming solar radiation for a given reflection is less than :math:`1 \times 10^{- 5}`) is met to ensure radiation is conserved. Direct beam and diffuse radiation are solved independently but follow the same solution steps. The solution below is for the direct beam component.

The initial direct beam absorption (:math:`i = 0`) (step 1) by each urban surface is

:math:`{\overrightarrow{S}}_{imprvrd,\Lambda,i = 0}^{\mu} = S_{imprvrd} \downarrow_{\Lambda}^{\mu}\left( 1 - \alpha_{imprvrd,\Lambda}^{\mu} \right)`,

:math:`{\overrightarrow{S}}_{prvrd,\Lambda,i = 0}^{\mu} = S_{prvrd} \downarrow_{\Lambda}^{\mu}\left( 1 - \alpha_{prvrd,\Lambda}^{\mu} \right)`,

:math:`{\overrightarrow{S}}_{sunwall,\Lambda,i = 0}^{\mu} = S_{sunwall} \downarrow_{\Lambda}^{\mu}\left( 1 - \alpha_{sunwall,\Lambda}^{\mu} \right)`,

:math:`{\overrightarrow{S}}_{shdwall,\Lambda,i = 0}^{\mu} = S_{shdwall} \downarrow_{\Lambda}^{\mu}\left( 1 - \alpha_{shdwall,\Lambda}^{\mu} \right)`,

.. math::
   :label: eq-0127

   {\overrightarrow{S}}_{road,\Lambda,i = 0}^{\mu} = {\overrightarrow{S}}_{imprvrd,\Lambda,i = 0}^{\mu}\left( 1 - f_{prvrd} \right) + {\overrightarrow{S}}_{prvrd,\Lambda,i = 0}^{\mu}f_{prvrd}

where, for example, :math:`S_{imprvrd} \downarrow_{\Lambda}^{\mu}` is the incident direct solar radiation for the impervious road (equation :eq:`eq-0113` ) and :math:`\alpha_{imprvrd,\Lambda}^{\mu}` is the direct albedo for the impervious road after adjustment for snow (section 2.1). Similarly, the initial reflections from each urban surface are

:math:`S_{imprvrd} \uparrow_{\Lambda,i = 0}^{\mu} = S_{imprvrd} \downarrow_{\Lambda}^{\mu}\left( \alpha_{imprvrd,\Lambda}^{\mu} \right)`,

:math:`S_{prvrd} \uparrow_{\Lambda,i = 0}^{\mu} = S_{prvrd} \downarrow_{\Lambda}^{\mu}\left( \alpha_{prvrd,\Lambda}^{\mu} \right)`,

.. math::
   :label: eq-0128

   S_{road} \uparrow_{\Lambda,i = 0}^{\mu} = S_{imprvrd} \downarrow_{\Lambda}^{\mu}\left( 1 - f_{prvrd} \right) + S_{prvrd} \downarrow_{\Lambda}^{\mu}f_{prvrd}

:math:`S_{sunwall} \uparrow_{\Lambda,i = 0}^{\mu} = S_{sunwall} \downarrow_{\Lambda}^{\mu}\left( \alpha_{sunwall,\Lambda}^{\mu} \right)`,

:math:`S_{shdwall} \uparrow_{\Lambda,i = 0}^{\mu} = S_{shdwall} \downarrow_{\Lambda}^{\mu}\left( \alpha_{shdwall,\Lambda}^{\mu} \right)`,

The initial reflected solar radiation is distributed to sky, walls, and road according to view factors as

.. math::
   :label: eq-0129

   S_{imprvrd - sky} \uparrow_{\Lambda,i = 0}^{\mu} = S_{imprvrd} \uparrow_{\Lambda,i = 0}^{\mu}\Psi_{road - sky}

.. math::
   :label: eq-0130

   S_{imprvrd - sunwall} \uparrow_{\Lambda,i = 0}^{\mu} = S_{imprvrd} \uparrow_{\Lambda,i = 0}^{\mu}\Psi_{road - wall}

.. math::
   :label: eq-0131

   S_{imprvrd - shdwall} \uparrow_{\Lambda,i = 0}^{\mu} = S_{imprvrd} \uparrow_{\Lambda,i = 0}^{\mu}\Psi_{road - wall}

.. math::
   :label: eq-0132

   S_{prvrd - sky} \uparrow_{\Lambda,i = 0}^{\mu} = S_{prvrd} \uparrow_{\Lambda,i = 0}^{\mu}\Psi_{road - sky}

.. math::
   :label: eq-0133

   S_{prvrd - sunwall} \uparrow_{\Lambda,i = 0}^{\mu} = S_{prvrd} \uparrow_{\Lambda,i = 0}^{\mu}\Psi_{road - wall}

.. math::
   :label: eq-0134

   S_{prvrd - shdwall} \uparrow_{\Lambda,i = 0}^{\mu} = S_{prvrd} \uparrow_{\Lambda,i = 0}^{\mu}\Psi_{road - wall}

.. math::
   :label: eq-0135

   S_{road - sky} \uparrow_{\Lambda,i = 0}^{\mu} = S_{road} \uparrow_{\Lambda,i = 0}^{\mu}\Psi_{road - sky}

.. math::
   :label: eq-0136

   S_{road - sunwall} \uparrow_{\Lambda,i = 0}^{\mu} = S_{road} \uparrow_{\Lambda,i = 0}^{\mu}\Psi_{road - wall}

.. math::
   :label: eq-0137

   S_{road - shdwall} \uparrow_{\Lambda,i = 0}^{\mu} = S_{road} \uparrow_{\Lambda,i = 0}^{\mu}\Psi_{road - wall}

.. math::
   :label: eq-0138

   S_{sunwall - sky} \uparrow_{\Lambda,i = 0}^{\mu} = S_{sunwall} \uparrow_{\Lambda,i = 0}^{\mu}\Psi_{wall - sky}

.. math::
   :label: eq-0139

   S_{sunwall - road} \uparrow_{\Lambda,i = 0}^{\mu} = S_{sunwall} \uparrow_{\Lambda,i = 0}^{\mu}\Psi_{wall - road}

.. math::
   :label: eq-0140

   S_{sunwall - shdwall} \uparrow_{\Lambda,i = 0}^{\mu} = S_{sunwall} \uparrow_{\Lambda,i = 0}^{\mu}\Psi_{wall - wall}

.. math::
   :label: eq-0141

   S_{shdwall - sky} \uparrow_{\Lambda,i = 0}^{\mu} = S_{shdwall} \uparrow_{\Lambda,i = 0}^{\mu}\Psi_{wall - sky}

.. math::
   :label: eq-0142

   S_{shdwall - road} \uparrow_{\Lambda,i = 0}^{\mu} = S_{shdwall} \uparrow_{\Lambda,i = 0}^{\mu}\Psi_{wall - road}

.. math::
   :label: eq-0143

   S_{shdwall - sunwall} \uparrow_{\Lambda,i = 0}^{\mu} = S_{shdwall} \uparrow_{\Lambda,i = 0}^{\mu}\Psi_{wall - wall}

The direct beam solar radiation absorbed by each urban surface after the :math:`i^{th}` reflection (steps 2 and 3) is



.. math::
   :label: eq-0005

   {\overrightarrow{S}}_{imprvrd,\Lambda,i}^{\mu} = {\overrightarrow{S}}_{imprvrd,\Lambda,i - 1}^{\mu} +
   \left( 1 - \alpha_{imprvrd,\Lambda}^{\mu} \right)\left( S_{sunwall - road} \uparrow_{\Lambda,i - 1}^{\mu} + S_{shdwall - road} \uparrow_{\Lambda,i - 1}^{\mu} \right)\frac{H}{W}





.. math::
   :label: eq-0006

   {\overrightarrow{S}}_{prvrd,\Lambda,i}^{\mu} = {\overrightarrow{S}}_{prvrd,\Lambda,i - 1}^{\mu} +
   \left( 1 - \alpha_{prvrd,\Lambda}^{\mu} \right)\left( S_{sunwall - road} \uparrow_{\Lambda,i - 1}^{\mu} + S_{shdwall - road} \uparrow_{\Lambda,i - 1}^{\mu} \right)\frac{H}{W}





.. math::
   :label: eq-0007

   {\overrightarrow{S}}_{sunwall,\Lambda,i}^{\mu} = {\overrightarrow{S}}_{sunwall,\Lambda,i - 1}^{\mu} +
   \left( 1 - \alpha_{sunwall,\Lambda}^{\mu} \right)\left( \frac{S_{road - sunwall} \uparrow_{\Lambda,i - 1}^{\mu}}{\frac{H}{W}} + S_{shdwall - sunwall} \uparrow_{\Lambda,i - 1}^{\mu} \right)





.. math::
   :label: eq-0008

   {\overrightarrow{S}}_{shdwall,\Lambda,i}^{\mu} = {\overrightarrow{S}}_{shdwall,\Lambda,i - 1}^{\mu} +
   \left( 1 - \alpha_{shdwall,\Lambda}^{\mu} \right)\left( \frac{S_{road - shdwall} \uparrow_{\Lambda,i - 1}^{\mu}}{\frac{H}{W}} + S_{sunwall - shdwall} \uparrow_{\Lambda,i - 1}^{\mu} \right)



The radiation from the walls to the road (:math:`S_{sunwall - road} \uparrow_{\Lambda,i - 1}^{\mu}`,\ :math:`S_{shdwall - road} \uparrow_{\Lambda,i - 1}^{\mu}`) is in W m\ :sup:`-2` of wall area and must be converted to W m\ :sup:`-2` of road area by multiplying by the height to width ratio. Similarly, the radiation from the road to the walls must be converted from W m\ :sup:`-2` of road area to W m\ :sup:`-2` of wall area by dividing by the height to width ratio. The direct beam solar radiation reflected by each urban surface after the :math:`i^{th}` reflection is distributed to sky, road, and walls (step 4) according to

.. math::
   :label: eq-0144

   S_{imprvrd - sky} \uparrow_{\Lambda,i}^{\mu} = \alpha_{imprvrd,\Lambda}^{\mu}\left( S_{sunwall - road} \uparrow_{\Lambda,i - 1}^{\mu} + S_{shdwall - road} \uparrow_{\Lambda,i - 1}^{\mu} \right)\frac{H}{W}\Psi_{road - sky}

.. math::
   :label: eq-0145

   S_{imprvrd - sunwall} \uparrow_{\Lambda,i}^{\mu} = \alpha_{imprvrd,\Lambda}^{\mu}\left( S_{sunwall - road} \uparrow_{\Lambda,i - 1}^{\mu} + S_{shdwall - road} \uparrow_{\Lambda,i - 1}^{\mu} \right)\frac{H}{W}\Psi_{road - wall}

.. math::
   :label: eq-0146

   S_{imprvrd - shdwall} \uparrow_{\Lambda,i}^{\mu} = \alpha_{imprvrd,\Lambda}^{\mu}\left( S_{sunwall - road} \uparrow_{\Lambda,i - 1}^{\mu} + S_{shdwall - road} \uparrow_{\Lambda,i - 1}^{\mu} \right)\frac{H}{W}\Psi_{road - wall}

.. math::
   :label: eq-0147

   S_{prvrd - sky} \uparrow_{\Lambda,i}^{\mu} = \alpha_{prvrd,\Lambda}^{\mu}\left( S_{sunwall - road} \uparrow_{\Lambda,i - 1}^{\mu} + S_{shdwall - road} \uparrow_{\Lambda,i - 1}^{\mu} \right)\frac{H}{W}\Psi_{road - sky}

.. math::
   :label: eq-0148

   S_{prvrd - sunwall} \uparrow_{\Lambda,i}^{\mu} = \alpha_{prvrd,\Lambda}^{\mu}\left( S_{sunwall - road} \uparrow_{\Lambda,i - 1}^{\mu} + S_{shdwall - road} \uparrow_{\Lambda,i - 1}^{\mu} \right)\frac{H}{W}\Psi_{road - wall}

.. math::
   :label: eq-0149

   S_{prvrd - shdwall} \uparrow_{\Lambda,i}^{\mu} = \alpha_{prvrd,\Lambda}^{\mu}\left( S_{sunwall - road} \uparrow_{\Lambda,i - 1}^{\mu} + S_{shdwall - road} \uparrow_{\Lambda,i - 1}^{\mu} \right)\frac{H}{W}\Psi_{road - wall}



.. math::
   :label: eq-0009

   S_{road - sky} \uparrow_{\Lambda,i}^{\mu} = \left\lbrack \begin{aligned}
   & \alpha_{imprvrd,\Lambda}^{\mu}\left( S_{sunwall - road} \uparrow_{\Lambda,i - 1}^{\mu} + S_{shdwall - road} \uparrow_{\Lambda,i - 1}^{\mu} \right) \\
   & \left( 1 - f_{prvrd} \right) + \alpha_{prvrd,\Lambda}^{\mu} \\
   & \left( S_{sunwall - road} \uparrow_{\Lambda,i - 1}^{\mu} + S_{shdwall - road} \uparrow_{\Lambda,i - 1}^{\mu} \right)f_{prvrd}
   \end{aligned} \right\rbrack\frac{H}{W}\Psi_{road - sky}





.. math::
   :label: eq-0010

   S_{road - sunwall} \uparrow_{\Lambda,i}^{\mu} = \left\lbrack \begin{aligned}
   & \alpha_{imprvrd,\Lambda}^{\mu}\left( S_{sunwall - road} \uparrow_{\Lambda,i - 1}^{\mu} + S_{shdwall - road} \uparrow_{\Lambda,i - 1}^{\mu} \right) \\
   & \left( 1 - f_{prvrd} \right) + \alpha_{prvrd,\Lambda}^{\mu} \\
   & \left( S_{sunwall - road} \uparrow_{\Lambda,i - 1}^{\mu} + S_{shdwall - road} \uparrow_{\Lambda,i - 1}^{\mu} \right)f_{prvrd}
   \end{aligned} \right\rbrack\frac{H}{W}\Psi_{road - wall}





.. math::
   :label: eq-0011

   S_{road - shdwall} \uparrow_{\Lambda,i}^{\mu} = \left\lbrack \begin{aligned}
   & \alpha_{imprvrd,\Lambda}^{\mu}\left( S_{sunwall - road} \uparrow_{\Lambda,i - 1}^{\mu} + S_{shdwall - road} \uparrow_{\Lambda,i - 1}^{\mu} \right) \\
   & \left( 1 - f_{prvrd} \right) + \alpha_{prvrd,\Lambda}^{\mu} \\
   & \left( S_{sunwall - road} \uparrow_{\Lambda,i - 1}^{\mu} + S_{shdwall - road} \uparrow_{\Lambda,i - 1}^{\mu} \right)f_{prvrd}
   \end{aligned} \right\rbrack\frac{H}{W}\Psi_{road - wall}



.. math::
   :label: eq-0150

   S_{sunwall - sky} \uparrow_{\Lambda,i}^{\mu} = \alpha_{sunwall,\Lambda}^{\mu}\left( \frac{S_{road - sunwall} \uparrow_{\Lambda,i - 1}^{\mu}}{\frac{H}{W}} + S_{shdwall - sunwall} \uparrow_{\Lambda,i - 1}^{\mu} \right)\Psi_{wall - sky}

.. math::
   :label: eq-0151

   S_{sunwall - road} \uparrow_{\Lambda,i}^{\mu} = \alpha_{sunwall,\Lambda}^{\mu}\left( \frac{S_{road - sunwall} \uparrow_{\Lambda,i - 1}^{\mu}}{\frac{H}{W}} + S_{shdwall - sunwall} \uparrow_{\Lambda,i - 1}^{\mu} \right)\Psi_{wall - road}

.. math::
   :label: eq-0152

   S_{sunwall - shdwall} \uparrow_{\Lambda,i}^{\mu} = \alpha_{sunwall,\Lambda}^{\mu}\left( \frac{S_{road - sunwall} \uparrow_{\Lambda,i - 1}^{\mu}}{\frac{H}{W}} + S_{shdwall - sunwall} \uparrow_{\Lambda,i - 1}^{\mu} \right)\Psi_{wall - wall}

.. math::
   :label: eq-0153

   S_{shdwall - sky} \uparrow_{\Lambda,i}^{\mu} = \alpha_{shdwall,\Lambda}^{\mu}\left( \frac{S_{road - shdwall} \uparrow_{\Lambda,i - 1}^{\mu}}{\frac{H}{W}} + S_{sunwall - shdwall} \uparrow_{\Lambda,i - 1}^{\mu} \right)\Psi_{wall - sky}

.. math::
   :label: eq-0154

   S_{shdwall - road} \uparrow_{\Lambda,i}^{\mu} = \alpha_{shdwall,\Lambda}^{\mu}\left( \frac{S_{road - shdwall} \uparrow_{\Lambda,i - 1}^{\mu}}{\frac{H}{W}} + S_{sunwall - shdwall} \uparrow_{\Lambda,i - 1}^{\mu} \right)\Psi_{wall - road}

.. math::
   :label: eq-0155

   S_{shdwall - sunwall} \uparrow_{\Lambda,i}^{\mu} = \alpha_{shdwall,\Lambda}^{\mu}\left( \frac{S_{road - shdwall} \uparrow_{\Lambda,i - 1}^{\mu}}{\frac{H}{W}} + S_{sunwall - shdwall} \uparrow_{\Lambda,i - 1}^{\mu} \right)\Psi_{wall - wall}

The reflected solar radiation to the sky is added to the total reflected solar radiation (step 5) for each urban surface as

.. math::
   :label: eq-0156

   S_{imprvrd} \uparrow_{\Lambda,i + 1}^{\mu} = S_{imprvrd} \uparrow_{\Lambda,i - 1}^{\mu} + S_{imprvrd - sky} \uparrow_{\Lambda,i}^{\mu}

.. math::
   :label: eq-0157

   S_{prvrd} \uparrow_{\Lambda,i + 1}^{\mu} = S_{prvrd} \uparrow_{\Lambda,i - 1}^{\mu} + S_{prvrd - sky} \uparrow_{\Lambda,i}^{\mu}

.. math::
   :label: eq-0158

   S_{sunwall} \uparrow_{\Lambda,i + 1}^{\mu} = S_{sunwall} \uparrow_{\Lambda,i - 1}^{\mu} + S_{sunwall - sky} \uparrow_{\Lambda,i}^{\mu}

.. math::
   :label: eq-0159

   S_{shdwall} \uparrow_{\Lambda,i + 1}^{\mu} = S_{shdwall} \uparrow_{\Lambda,i - 1}^{\mu} + S_{shdwall - sky} \uparrow_{\Lambda,i}^{\mu}

The system of equations (Equations :eq:`eq-0005`-(2.85)) is iterated for :math:`i = 50` reflections or until the absorption for the :math:`i^{th}` reflection is less than a nominal amount

.. math::
   :label: eq-0160

   \max\left( \frac{{\overrightarrow{S}}_{road,\Lambda,i}^{\mu}}{S_{atm} \downarrow_{\Lambda}^{\mu}},\frac{{\overrightarrow{S}}_{sunwall,\Lambda,i}^{\mu}}{S_{atm} \downarrow_{\Lambda}^{\mu}},\frac{{\overrightarrow{S}}_{shdwall,\Lambda,i}^{\mu}}{S_{atm} \downarrow_{\Lambda}^{\mu}} \right) < 1 \times 10^{- 5}

where :math:`{\overrightarrow{S}}_{sunwall,\Lambda,i}^{\mu}` (equation :eq:`eq-0007` ) and :math:`{\overrightarrow{S}}_{shdwall,\Lambda,i}^{\mu}` (equation :eq:`eq-0008` ) are the direct beam solar radiation absorbed by the sunlit wall and shaded wall on the :math:`i^{th}` reflection, and



.. math::
   :label: eq-0012

   {\overrightarrow{S}}_{road,\Lambda,i}^{\mu} = \left( 1 - \alpha_{imprvrd,\Lambda}^{\mu} \right)\left( S_{sunwall - road} \uparrow_{\Lambda,i - 1}^{\mu} + S_{shdwall - road} \uparrow_{\Lambda,i - 1}^{\mu} \right)\frac{H}{W}\left( 1 - f_{prvrd} \right)
   + \left( 1 - \alpha_{prvrd,\Lambda}^{\mu} \right)\left( S_{sunwall - road} \uparrow_{\Lambda,i - 1}^{\mu} + S_{shdwall - road} \uparrow_{\Lambda,i - 1}^{\mu} \right)\frac{H}{W}f_{prvrd}



is the direct beam solar radiation absorbed by the road on the :math:`i^{th}` reflection.

The total direct beam and diffuse solar radiation reflected by the urban canyon (walls and road) is



.. math::
   :label: eq-0013

   S_{uc} \uparrow_{\Lambda}^{\mu} = S_{imprvrd} \uparrow_{\Lambda,i = n + 1}^{\mu}\left( 1 - f_{prvrd} \right) + S_{prvrd} \uparrow_{\Lambda,i = n + 1}^{\mu}f_{prvrd}
   + \left( S_{sunwall} \uparrow_{\Lambda,i = n + 1}^{\mu} + S_{shdwall} \uparrow_{\Lambda,i = n + 1}^{\mu} \right)\frac{H}{W}





.. math::
   :label: eq-0014

   S_{uc} \uparrow_{\Lambda} = S_{imprvrd} \uparrow_{\Lambda,i = n + 1}\left( 1 - f_{prvrd} \right) + S_{prvrd} \uparrow_{\Lambda,i = n + 1}f_{prvrd}
   + \left( S_{sunwall} \uparrow_{\Lambda,i = n + 1} + S_{shdwall} \uparrow_{\Lambda,i = n + 1} \right)\frac{H}{W}



while the total absorbed is



.. math::
   :label: eq-0015

   {\overrightarrow{S}}_{uc,\Lambda}^{\mu} = {\overrightarrow{S}}_{imprvrd,\Lambda,i = n}^{\mu}\left( 1 - f_{prvrd} \right) + {\overrightarrow{S}}_{prvrd,\Lambda,i = n}^{\mu}f_{prvrd}
   + \left( {\overrightarrow{S}}_{sunwall,\Lambda,i = n}^{\mu} + {\overrightarrow{S}}_{shdwall,\Lambda,i = n}^{\mu} \right)\frac{H}{W}




.. math::
   :label: eq-0016

   {\overrightarrow{S}}_{uc,\Lambda} = {\overrightarrow{S}}_{imprvrd,\Lambda,i = n}\left( 1 - f_{prvrd} \right) + {\overrightarrow{S}}_{prvrd,\Lambda,i = n}f_{prvrd}
   + \left( {\overrightarrow{S}}_{sunwall,\Lambda,i = n} + {\overrightarrow{S}}_{shdwall,\Lambda,i = n} \right)\frac{H}{W}


Solar radiation in the urban canyon is conserved as



.. math::
   :label: eq-0017

   S_{road} \downarrow_{\Lambda}^{\mu} + \left( S_{sunwall} \downarrow_{\Lambda}^{\mu} + S_{shdwall} \downarrow_{\Lambda}^{\mu} \right)\frac{H}{W} + S_{road} \downarrow_{\Lambda} + \left( S_{sunwall} \downarrow_{\Lambda} + S_{shdwall} \downarrow_{\Lambda} \right)\frac{H}{W}
   - \left( {\overrightarrow{S}}_{uc,\Lambda}^{\mu} + {\overrightarrow{S}}_{uc,\Lambda} + S_{uc} \uparrow_{\Lambda}^{\mu} + S_{uc} \uparrow_{\Lambda} \right) = 0


The direct beam and diffuse urban canyon albedos are

:math:`\alpha_{uc,\Lambda}^{\mu} = \frac{S_{uc} \uparrow_{\Lambda}^{\mu}}{S_{road} \downarrow_{\Lambda}^{\mu} + \left( S_{sunwall} \downarrow_{\Lambda}^{\mu} + S_{shdwall} \downarrow_{\Lambda}^{\mu} \right)\frac{H}{W}}`,

.. math::
   :label: eq-0161

   \alpha_{uc,\Lambda} = \frac{S_{uc} \uparrow_{\Lambda}}{S_{road} \downarrow_{\Lambda} + \left( S_{sunwall} \downarrow_{\Lambda} + S_{shdwall} \downarrow_{\Lambda} \right)\frac{H}{W}}

The total absorbed solar radiation for the urban canopy (road, walls, and roof) is

.. math::
   :label: eq-0162

   \overrightarrow{S} = \sum_{\Lambda}^{}\left\lbrack W_{roof}\left( {\overrightarrow{S}}_{roof,\Lambda}^{\mu} + {\overrightarrow{S}}_{roof,\Lambda} \right) + \left( 1 - W_{roof} \right)\left( {\overrightarrow{S}}_{uc,\Lambda}^{\mu} + {\overrightarrow{S}}_{uc,\Lambda} \right) \right\rbrack

:numref:`fig-absorbed-solar-radiation` shows the solar radiation absorbed by urban surfaces for a range of height to width ratios and two solar zenith angles. The absorbed solar radiation for the roof is independent of height to width ratio and solar zenith angle. At both solar zenith angles, the absorbed solar radiation for the road decreases rapidly with increasing height to width ratio as the buildings shade more of the road. The shaded wall absorbs less solar radiation than the sunlit wall because it receives only diffuse radiation from the sun and reflected radiation from the walls and road. The sunlit wall absorbs more solar radiation at larger solar zenith angles for height to width ratios less than about three because the incidence angle of the radiation is closer to zero (:numref:`fig-solar-elevation-view`). The sum of the absorbed solar radiation for road, sunlit wall, and shaded wall, after converting the wall fluxes to per unit ground area, is the canyon absorbed solar radiation. The absorbed solar radiation for the canyon increases slowly with increasing height to width ratio.

.. figure:: image8.png
   :width: 4.21875in
   :height: 6.40625in
   :name: fig-absorbed-solar-radiation

   Solar radiation absorbed by urban surfaces for solar zenith angles of 30° (top) and 60° (bottom). The atmospheric solar radiation is :math:`S_{atm} \downarrow_{\Lambda}^{\mu} = 400` and :math:`S_{atm} \downarrow_{\Lambda} = 200` W m\ :sup:`-2`. Note that the sunlit and shaded wall fluxes are per unit wall area. The solar radiation absorbed by the canyon is the sum of road and wall fluxes after converting the walls fluxes to per unit ground area using the height to width ratio.

The canyon albedo (excluding the roof albedo) shown in :numref:`fig-canyon-albedo` has the same functional relationships with solar zenith angle and height to width ratio as TEB (Masson 2000). In general, the direct and diffuse canyon albedo decreases with height to width ratio as more solar radiation is trapped and absorbed within the canyon. The trapping of solar radiation is less effective at larger solar zenith angles. At these large solar zenith angles and small height to width ratio, the albedo increases because the higher albedo walls dominate the radiative exchange.

.. figure:: image9.png
   :width: 5.98958in
   :height: 4.14583in
   :name: fig-canyon-albedo

   Direct beam and diffuse albedo of the urban canyon (walls and road) as a function of height to width ratio from 0.1 to 3.0 in increments of 0.1 and solar zenith angles from 0° to 85° in increments of 5°. The atmospheric solar radiation is :math:`S_{atm} \downarrow_{\Lambda}^{\mu} = 400` and :math:`S_{atm} \downarrow_{\Lambda} = 200` W m\ :sup:`-2`.

2.6 Incident longwave radiation
-------------------------------

Similar to incident diffuse solar radiation, the longwave radiation incident on walls and roads depends on view factors. The longwave radiation incident on roof, walls and road is

:math:`L_{roof} \downarrow = L_{atm} \downarrow`,

:math:`L_{imprvrd} \downarrow = L_{prvrd} \downarrow = L_{atm} \downarrow \Psi_{sky - road}`,

:math:`L_{shdwall} \downarrow = L_{atm} \downarrow \Psi_{sky - wall}`,

.. math::
   :label: eq-0163

   L_{sunwall} \downarrow = L_{atm} \downarrow \Psi_{sky - wall}

where :math:`L_{atm} \downarrow` is the longwave radiation from the atmosphere. The incident longwave radiation conserves energy as



.. math::
   :label: eq-0018

   L_{atm} \downarrow = f_{roof}L_{roof} \downarrow + \left( 1 - f_{roof} \right)
   \left\lbrack L_{imprvrd} \downarrow \left( 1 - f_{prvrd} \right) + L_{prvrd} \downarrow f_{prvrd} + \frac{H}{W}\left( L_{sunwall} \downarrow + L_{shdwall} \downarrow \right) \right\rbrack


Absorbed, reflected, and emitted longwave radiation
---------------------------------------------------

Emitted longwave radiation, a function of surface temperature and emissivity, must also be considered in addition to reflection and absorption when determining the longwave interactions within the canyon. The net longwave radiation (W m\ :sup:`-2`) (positive toward the atmosphere) for the roof is simply

.. math::
   :label: eq-0164

   {\overrightarrow{L}}_{roof} = L_{roof} \uparrow - L_{atm} \downarrow

where

.. math::
   :label: eq-0165

   L_{roof} \uparrow = \varepsilon_{roof}\sigma\left( T_{roof} \right)^{4} + \left( 1 - \varepsilon_{roof} \right)L_{atm} \downarrow

is the emitted plus reflected longwave radiation from the roof, :math:`\varepsilon_{roof}` is the emissivity of the roof, :math:`\sigma` is the Stefan-Boltzmann constant (W m\ :sup:`-2` K\ :sup:`-4`) (:numref:`table-physical-constants`), and :math:`T_{roof}` is the temperature of the roof (section 4).

Similar to albedo, the emissivity of each urban surface is a weighted combination of snow-free surface and snow emissivity. Only roof and road surfaces are affected by snow as

.. math::
   :label: eq-0166

   \varepsilon_{u} = \varepsilon_{u}\left( 1 - f_{u,sno} \right) + \varepsilon_{sno}f_{u,sno}

where :math:`\varepsilon_{u}` is the emissivity of :math:`u =` roof, pervious and impervious road (:numref:`table-input-data`), :math:`\varepsilon_{sno} = 0.97` is the emissivity of snow (Oleson et al. 2004), and :math:`f_{u,sno}` is the fraction of the urban surface covered with snow (equation (2.3)).

As with solar radiation, the longwave interactions within the urban canyon are determined numerically by allowing for multiple reflections until a convergence criteria is met (the absorbed longwave radiation for a given reflection is less than :math:`1 \times 10^{- 3}`). The following equation (2.3) assume that absorptivity equals emissivity.

The initial reflected (:math:`r`) longwave radiation from each urban surface is

:math:`L_{imprvrd,i = 0}\overset{r}{\uparrow} = \left( 1 - \varepsilon_{imprvrd} \right)L_{imprvrd} \downarrow`,

:math:`L_{prvrd,i = 0}\overset{r}{\uparrow} = \left( 1 - \varepsilon_{prvrd} \right)L_{prvrd} \downarrow`,

.. math::
   :label: eq-0167

   L_{road,i = 0}\overset{r}{\uparrow} = L_{imprvrd,i = 0}\overset{r}{\uparrow}\left( 1 - f_{prvrd} \right) + L_{prvrd,i = 0}\overset{r}{\uparrow}f_{prvrd}

:math:`L_{sunwall,i = 0}\overset{r}{\uparrow} = \left( 1 - \varepsilon_{wall} \right)L_{sunwall} \downarrow`,

.. math::
   :label: eq-0168

   L_{shdwall,i = 0}\overset{r}{\uparrow} = \left( 1 - \varepsilon_{wall} \right)L_{shdwall} \downarrow

The emitted (:math:`e`) longwave radiation from each surface is

:math:`L_{imprvrd}\overset{e}{\uparrow} = \varepsilon_{imprvrd}\sigma\left( T_{imprvrd} \right)^{4}`,

:math:`L_{prvrd}\overset{e}{\uparrow} = \varepsilon_{prvrd}\sigma\left( T_{prvrd} \right)^{4}`,

:math:`L_{road}\overset{e}{\uparrow} = \varepsilon_{imprvrd}\sigma\left( T_{imprvrd} \right)^{4}\left( 1 - f_{prvrd} \right) + \varepsilon_{prvrd}\sigma\left( T_{prvrd} \right)^{4}f_{prvrd}`,

:math:`L_{sunwall}\overset{e}{\uparrow} = \varepsilon_{wall}\sigma\left( T_{sunwall} \right)^{4}`,

.. math::
   :label: eq-0169

   L_{shdwall}\overset{e}{\uparrow} = \varepsilon_{wall}\sigma\left( T_{shdwall} \right)^{4}

The initial reflected longwave radiation is distributed to sky, walls, and road according to view factors as

:math:`L_{imprvrd - sky,i = 0}\overset{r}{\uparrow} = L_{imprvrd,i = 0}\overset{r}{\uparrow}\Psi_{road - sky}`,

:math:`L_{prvrd - sky,i = 0}\overset{r}{\uparrow} = L_{prvrd,i = 0}\overset{r}{\uparrow}\Psi_{road - sky}`,

:math:`L_{road - sunwall,i = 0}\overset{r}{\uparrow} = L_{road,i = 0}\overset{r}{\uparrow}\Psi_{road - wall}`,

:math:`L_{road - shdwall,i = 0}\overset{r}{\uparrow} = L_{road,i = 0}\overset{r}{\uparrow}\Psi_{road - wall}`,

:math:`L_{sunwall - sky,i = 0}\overset{r}{\uparrow} = L_{sunwall,i = 0}\overset{r}{\uparrow}\Psi_{wall - sky}`,

:math:`L_{sunwall - road,i = 0}\overset{r}{\uparrow} = L_{sunwall,i = 0}\overset{r}{\uparrow}\Psi_{wall - road}`,

:math:`L_{sunwall - shdwall,i = 0}\overset{r}{\uparrow} = L_{sunwall,i = 0}\overset{r}{\uparrow}\Psi_{wall - wall}`,

:math:`L_{shdwall - sky,i = 0}\overset{r}{\uparrow} = L_{shdwall,i = 0}\overset{r}{\uparrow}\Psi_{wall - sky}`,

:math:`L_{shdwall - road,i = 0}\overset{r}{\uparrow} = L_{shdwall,i = 0}\overset{r}{\uparrow}\Psi_{wall - road}`,

.. math::
   :label: eq-0170

   L_{shdwall - sunwall,i = 0}\overset{r}{\uparrow} = L_{shdwall,i = 0}\overset{r}{\uparrow}\Psi_{wall - wall}

The emitted longwave radiation is distributed to sky, walls, and road according to view factors as

:math:`L_{imprvrd - sky}\overset{e}{\uparrow} = L_{imprvrd}\overset{e}{\uparrow}\Psi_{road - sky}`,

:math:`L_{prvrd - sky}\overset{e}{\uparrow} = L_{prvrd}\overset{e}{\uparrow}\Psi_{road - sky}`,

:math:`L_{road - sunwall}\overset{e}{\uparrow} = L_{road}\overset{e}{\uparrow}\Psi_{road - wall}`,

:math:`L_{road - shdwall}\overset{e}{\uparrow} = L_{road}\overset{e}{\uparrow}\Psi_{road - wall}`,

:math:`L_{sunwall - sky}\overset{e}{\uparrow} = L_{sunwall}\overset{e}{\uparrow}\Psi_{wall - sky}`,

:math:`L_{sunwall - road}\overset{e}{\uparrow} = L_{sunwall}\overset{e}{\uparrow}\Psi_{wall - road}`,

:math:`L_{sunwall - shdwall}\overset{e}{\uparrow} = L_{sunwall}\overset{e}{\uparrow}\Psi_{wall - wall}`,

:math:`L_{shdwall - sky}\overset{e}{\uparrow} = L_{shdwall}\overset{e}{\uparrow}\Psi_{wall - sky}`,

:math:`L_{shdwall - road}\overset{e}{\uparrow} = L_{shdwall}\overset{e}{\uparrow}\Psi_{wall - road}`,

.. math::
   :label: eq-0171

   L_{shdwall - sunwall}\overset{e}{\uparrow} = L_{shdwall}\overset{e}{\uparrow}\Psi_{wall - wall}

The initial absorption (net longwave) (:math:`i = 0`) by each urban surface is

:math:`{\overrightarrow{L}}_{imprvrd,i = 0} = L_{imprvrd}\overset{e}{\uparrow} - \varepsilon_{imprvrd}L_{imprvrd} \downarrow`,

:math:`{\overrightarrow{L}}_{prvrd,i = 0} = L_{prvrd}\overset{e}{\uparrow} - \varepsilon_{prvrd}L_{prvrd} \downarrow`,

:math:`{\overrightarrow{L}}_{sunwall,i = 0} = L_{sunwall}\overset{e}{\uparrow} - \varepsilon_{wall}L_{sunwall} \downarrow`,

.. math::
   :label: eq-0172

   {\overrightarrow{L}}_{shdwall,i = 0} = L_{shdwall}\overset{e}{\uparrow} - \varepsilon_{wall}L_{shdwall} \downarrow

The initial emitted plus reflected longwave radiation to the sky is

:math:`L_{imprvrd - sky,i = 0}\overset{}{\uparrow} = L_{imprvrd - sky}\overset{e}{\uparrow} + L_{imprvrd - sky,i = 0}\overset{r}{\uparrow}`,

:math:`L_{prvrd - sky,i = 0}\overset{}{\uparrow} = L_{prvrd - sky}\overset{e}{\uparrow} + L_{prvrd - sky,i = 0}\overset{r}{\uparrow}`,

:math:`L_{sunwall - sky,i = 0}\overset{}{\uparrow} = L_{sunwall - sky}\overset{e}{\uparrow} + L_{sunwall - sky,i = 0}\overset{r}{\uparrow}`,

.. math::
   :label: eq-0173

   L_{shdwall - sky,i = 0}\overset{}{\uparrow} = L_{shdwall - sky}\overset{e}{\uparrow} + L_{shdwall - sky,i = 0}\overset{r}{\uparrow}

The net longwave radiation absorbed by each urban surface after the :math:`i^{th}` reflection is



.. math::
   :label: eq-0019

   {\overrightarrow{L}}_{imprvrd,i} = \varepsilon_{imprvrd}\left( \begin{aligned}
   & L_{sunwall - road,i - 1}\overset{r}{\uparrow} + L_{sunwall - road}\overset{e}{\uparrow} \\
   & + L_{shdwall - road,i - 1}\overset{r}{\uparrow} + L_{shdwall - road}\overset{e}{\uparrow}
   \end{aligned} \right)\frac{H}{W}


.. math::
   :label: eq-0020

   {\overrightarrow{L}}_{prvrd,i} = \varepsilon_{prvrd}\left( \begin{aligned}
   & L_{sunwall - road,i - 1}\overset{r}{\uparrow} + L_{sunwall - road}\overset{e}{\uparrow} \\
   & + L_{shdwall - road,i - 1}\overset{r}{\uparrow} + L_{shdwall - road}\overset{e}{\uparrow}
   \end{aligned} \right)\frac{H}{W}


.. math::
   :label: eq-0021

   {\overrightarrow{L}}_{road,i} = {\overrightarrow{L}}_{imprvrd,i}\left( 1 - f_{prvrd} \right) + {\overrightarrow{L}}_{prvrd,i}f_{prvrd}`


.. math::
   :label: eq-0022

   {\overrightarrow{L}}_{sunwall,i} = \varepsilon_{wall}\left( \begin{aligned}
   & \frac{L_{road - sunwall,i - 1}\overset{r}{\uparrow} + L_{road - sunwall}\overset{e}{\uparrow}}{\frac{H}{W}} \\
   & + L_{shdwall - sunwall,i - 1}\overset{r}{\uparrow} + L_{shdwall - sunwall}\overset{e}{\uparrow}
   \end{aligned} \right)


.. math::
   :label: eq-0023

   {\overrightarrow{L}}_{shdwall,i} = \varepsilon_{wall}\left( \begin{aligned}
   & \frac{L_{road - shdwall,i - 1}\overset{r}{\uparrow} + L_{road - shdwall}\overset{e}{\uparrow}}{\frac{H}{W}} \\
   & + L_{sunwall - shdwall,i - 1}\overset{r}{\uparrow} + L_{sunwall - shdwall}\overset{e}{\uparrow}
   \end{aligned} \right)


The longwave radiation from each urban surface after the :math:`i^{th}` reflection is distributed to sky, road, and walls according to



.. math::
   :label: eq-0024

   L_{imprvrd - sky,i} \uparrow = \left( 1 - \varepsilon_{imprvrd} \right)\left( \begin{aligned}
   & L_{sunwall - road,i - 1}\overset{r}{\uparrow} + L_{sunwall - road}\overset{e}{\uparrow} \\
   & + L_{shdwall - road,i - 1}\overset{r}{\uparrow} + L_{shdwall - road}\overset{e}{\uparrow}
   \end{aligned} \right)\frac{H}{W}\Psi_{road - sky}

,



.. math::
   :label: eq-0025

   L_{imprvrd - sunwall,i} \uparrow = \left( 1 - \varepsilon_{imprvrd} \right)\left( \begin{aligned}
   & L_{sunwall - road,i - 1}\overset{r}{\uparrow} + L_{sunwall - road}\overset{e}{\uparrow} \\
   & + L_{shdwall - road,i - 1}\overset{r}{\uparrow} + L_{shdwall - road}\overset{e}{\uparrow}
   \end{aligned} \right)\frac{H}{W}\Psi_{road - wall}

,



.. math::
   :label: eq-0026

   L_{imprvrd - shdwall,i} \uparrow = \left( 1 - \varepsilon_{imprvrd} \right)\left( \begin{aligned}
   & L_{sunwall - road,i - 1}\overset{r}{\uparrow} + L_{sunwall - road}\overset{e}{\uparrow} \\
   & + L_{shdwall - road,i - 1}\overset{r}{\uparrow} + L_{shdwall - road}\overset{e}{\uparrow}
   \end{aligned} \right)\frac{H}{W}\Psi_{road - wall}

,



.. math::
   :label: eq-0027

   L_{prvrd - sky,i} \uparrow = \left( 1 - \varepsilon_{prvrd} \right)\left( \begin{aligned}
   & L_{sunwall - road,i - 1}\overset{r}{\uparrow} + L_{sunwall - road}\overset{e}{\uparrow} \\
   & + L_{shdwall - road,i - 1}\overset{r}{\uparrow} + L_{shdwall - road}\overset{e}{\uparrow}
   \end{aligned} \right)\frac{H}{W}\Psi_{road - sky}

,



.. math::
   :label: eq-0028

   L_{prvrd - sunwall,i} \uparrow = \left( 1 - \varepsilon_{prvrd} \right)\left( \begin{aligned}
   & L_{sunwall - road,i - 1}\overset{r}{\uparrow} + L_{sunwall - road}\overset{e}{\uparrow} \\
   & + L_{shdwall - road,i - 1}\overset{r}{\uparrow} + L_{shdwall - road}\overset{e}{\uparrow}
   \end{aligned} \right)\frac{H}{W}\Psi_{road - wall}

,



.. math::
   :label: eq-0029

   L_{prvrd - shdwall,i} \uparrow = \left( 1 - \varepsilon_{prvrd} \right)\left( \begin{aligned}
   & L_{sunwall - road,i - 1}\overset{r}{\uparrow} + L_{sunwall - road}\overset{e}{\uparrow} \\
   & + L_{shdwall - road,i - 1}\overset{r}{\uparrow} + L_{shdwall - road}\overset{e}{\uparrow}
   \end{aligned} \right)\frac{H}{W}\Psi_{road - wall}

,



.. math::
   :label: eq-0030

   L_{road - sky,i} \uparrow = \left\lbrack \begin{aligned}
   & \left( 1 - \varepsilon_{imprvrd} \right)\left( \begin{aligned}
   & L_{sunwall - road,i - 1}\overset{r}{\uparrow} + L_{sunwall - road}\overset{e}{\uparrow} \\
   & + L_{shdwall - road,i - 1}\overset{r}{\uparrow} + L_{shdwall - road}\overset{e}{\uparrow}
   \end{aligned} \right) \\
   & \times \frac{H}{W}\left( 1 - f_{prvrd} \right) + \left( 1 - \varepsilon_{prvrd} \right) \\
   & \times \left( \begin{aligned}
   & L_{sunwall - road,i - 1}\overset{r}{\uparrow} + L_{sunwall - road}\overset{e}{\uparrow} \\
   & + L_{shdwall - road,i - 1}\overset{r}{\uparrow} + L_{shdwall - road}\overset{e}{\uparrow}
   \end{aligned} \right)\frac{H}{W}f_{prvrd}
   \end{aligned} \right\rbrack\Psi_{road - sky}

,



.. math::
   :label: eq-0031

   L_{road - sunwall,i} \uparrow = \left\lbrack \begin{aligned}
   & \left( 1 - \varepsilon_{imprvrd} \right)\left( \begin{aligned}
   & L_{sunwall - road,i - 1}\overset{r}{\uparrow} + L_{sunwall - road}\overset{e}{\uparrow} \\
   & + L_{shdwall - road,i - 1}\overset{r}{\uparrow} + L_{shdwall - road}\overset{e}{\uparrow}
   \end{aligned} \right) \\
   & \times \frac{H}{W}\left( 1 - f_{prvrd} \right) + \left( 1 - \varepsilon_{prvrd} \right) \\
   & \times \left( \begin{aligned}
   & L_{sunwall - road,i - 1}\overset{r}{\uparrow} + L_{sunwall - road}\overset{e}{\uparrow} \\
   & + L_{shdwall - road,i - 1}\overset{r}{\uparrow} + L_{shdwall - road}\overset{e}{\uparrow}
   \end{aligned} \right)\frac{H}{W}f_{prvrd}
   \end{aligned} \right\rbrack\Psi_{road - wall}

,



.. math::
   :label: eq-0032

   L_{road - shdwall,i} \uparrow = \left\lbrack \begin{aligned}
   & \left( 1 - \varepsilon_{imprvrd} \right)\left( \begin{aligned}
   & L_{sunwall - road,i - 1}\overset{r}{\uparrow} + L_{sunwall - road}\overset{e}{\uparrow} \\
   & + L_{shdwall - road,i - 1}\overset{r}{\uparrow} + L_{shdwall - road}\overset{e}{\uparrow}
   \end{aligned} \right) \\
   & \times \frac{H}{W}\left( 1 - f_{prvrd} \right) + \left( 1 - \varepsilon_{prvrd} \right) \\
   & \times \left( \begin{aligned}
   & L_{sunwall - road,i - 1}\overset{r}{\uparrow} + L_{sunwall - road}\overset{e}{\uparrow} \\
   & + L_{shdwall - road,i - 1}\overset{r}{\uparrow} + L_{shdwall - road}\overset{e}{\uparrow}
   \end{aligned} \right)\frac{H}{W}f_{prvrd}
   \end{aligned} \right\rbrack\Psi_{road - wall}

,



.. math::
   :label: eq-0033

   L_{sunwall - sky,i} \uparrow = \left( 1 - \varepsilon_{wall} \right)\left( \begin{aligned}
   & \frac{L_{road - sunwall,i - 1}\overset{r}{\uparrow} + L_{road - sunwall}\overset{e}{\uparrow}}{\frac{H}{W}} \\
   & + L_{shdwall - sunwall,i - 1}\overset{r}{\uparrow} + L_{shdwall - sunwall}\overset{e}{\uparrow}
   \end{aligned} \right)\Psi_{wall - sky}

,



.. math::
   :label: eq-0034

   L_{sunwall - road,i} \uparrow = \left( 1 - \varepsilon_{wall} \right)\left( \begin{aligned}
   & \frac{L_{road - sunwall,i - 1}\overset{r}{\uparrow} + L_{road - sunwall}\overset{e}{\uparrow}}{\frac{H}{W}} \\
   & + L_{shdwall - sunwall,i - 1}\overset{r}{\uparrow} + L_{shdwall - sunwall}\overset{e}{\uparrow}
   \end{aligned} \right)\Psi_{wall - road}

,



.. math::
   :label: eq-0035

   L_{sunwall - shdwall,i} \uparrow = \left( 1 - \varepsilon_{wall} \right)\left( \begin{aligned}
   & \frac{L_{road - sunwall,i - 1}\overset{r}{\uparrow} + L_{road - sunwall}\overset{e}{\uparrow}}{\frac{H}{W}} \\
   & + L_{shdwall - sunwall,i - 1}\overset{r}{\uparrow} + L_{shdwall - sunwall}\overset{e}{\uparrow}
   \end{aligned} \right)\Psi_{wall - wall}

,



.. math::
   :label: eq-0036

   L_{shdwall - sky,i} \uparrow = \left( 1 - \varepsilon_{wall} \right)\left( \begin{aligned}
   & \frac{L_{road - shdwall,i - 1}\overset{r}{\uparrow} + L_{road - shdwall}\overset{e}{\uparrow}}{\frac{H}{W}} \\
   & + L_{sunwall - shdwall,i - 1}\overset{r}{\uparrow} + L_{sunwall - shdwall}\overset{e}{\uparrow}
   \end{aligned} \right)\Psi_{wall - sky}

,



.. math::
   :label: eq-0037

   L_{shdwall - road,i} \uparrow = \left( 1 - \varepsilon_{wall} \right)\left( \begin{aligned}
   & \frac{L_{road - shdwall,i - 1}\overset{r}{\uparrow} + L_{road - shdwall}\overset{e}{\uparrow}}{\frac{H}{W}} \\
   & + L_{sunwall - shdwall,i - 1}\overset{r}{\uparrow} + L_{sunwall - shdwall}\overset{e}{\uparrow}
   \end{aligned} \right)\Psi_{wall - road}

,



.. math::
   :label: eq-0038

   L_{shdwall - sunwall,i} \uparrow = \left( 1 - \varepsilon_{wall} \right)\left( \begin{aligned}
   & \frac{L_{road - shdwall,i - 1}\overset{r}{\uparrow} + L_{road - shdwall}\overset{e}{\uparrow}}{\frac{H}{W}} \\
   & + L_{sunwall - shdwall,i - 1}\overset{r}{\uparrow} + L_{sunwall - shdwall}\overset{e}{\uparrow}
   \end{aligned} \right)\Psi_{wall - wall}


Note that the emitted longwave term in equations :eq:`eq-0019`-:eq:`eq-0038` only applies to the first iteration. Subsequent iterations do not include this term, i.e.,



.. math::
   :label: eq-0039

   L_{road - sunwall}\overset{e}{\uparrow} = L_{road - shdwall}\overset{e}{\uparrow} = L_{sunwall - road}\overset{e}{\uparrow} = L_{shdwall - road}\overset{e}{\uparrow}
   = L_{shdwall - sunwall}\overset{e}{\uparrow} = L_{sunwall - shdwall}\overset{e}{\uparrow} = 0


The reflected longwave radiation to the sky is added to the total upward longwave radiation for each urban surface as

:math:`L_{imprvrd,i + 1} \uparrow = L_{imprvrd,i - 1} \uparrow + L_{imprvrd - sky,i} \uparrow`,

:math:`L_{prvrd,i + 1} \uparrow = L_{prvrd,i - 1} \uparrow + L_{prvrd - sky,i} \uparrow`,

:math:`L_{sunwall,i + 1} \uparrow = L_{sunwall,i - 1} \uparrow + L_{sunwall - sky,i} \uparrow`,

.. math::
   :label: eq-0174

   L_{shdwall,i + 1} \uparrow = L_{shdwall,i - 1} \uparrow + L_{shdwall - sky,i} \uparrow

The net longwave at each iteration is added to the total net longwave for each urban surface as

:math:`{\overrightarrow{L}}_{imprvrd,i + 1} = {\overrightarrow{L}}_{imprvrd,i - 1} + {\overrightarrow{L}}_{imprvrd,i}`,

:math:`{\overrightarrow{L}}_{prvrd,i + 1} = {\overrightarrow{L}}_{prvrd,i - 1} + {\overrightarrow{L}}_{prvrd,i}`,

:math:`{\overrightarrow{L}}_{sunwall,i + 1} = {\overrightarrow{L}}_{sunwall,i - 1} + {\overrightarrow{L}}_{sunwall,i}`,

.. math::
   :label: eq-0175

   {\overrightarrow{L}}_{shdwall,i + 1} = {\overrightarrow{L}}_{shdwall,i - 1} + {\overrightarrow{L}}_{shdwall,i}

The system of equations (equations :eq:`eq-0019`-:eq:`eq-0175`) is iterated for :math:`i = 50` reflections or until the absorption for the :math:`i^{th}` reflection is less than a nominal amount

.. math::
   :label: eq-0176

   \max\left( {\overrightarrow{L}}_{road,i},{\overrightarrow{L}}_{sunwall,i},{\overrightarrow{L}}_{shdwall,i} \right) < 1 \times 10^{- 3}

The net longwave radiation for the urban canyon (walls and road) is

.. math::
   :label: eq-0177

   {\overrightarrow{L}}_{uc} = {\overrightarrow{L}}_{imprvrd,n + 1}\left( 1 - f_{prvrd} \right) + {\overrightarrow{L}}_{prvrd,n + 1}f_{prvrd} + \left( {\overrightarrow{L}}_{sunwall,n + 1} + {\overrightarrow{L}}_{shdwall,n + 1} \right)\frac{H}{W}

while the total reflected plus emitted longwave radiation is



.. math::
   :label: eq-0040

   L_{uc} \uparrow = L_{imprvrd,n + 1} \uparrow \left( 1 - f_{prvrd} \right) + L_{prvrd,n + 1} \uparrow f_{prvrd}
   + \left( L_{sunwall,n + 1} \uparrow + L_{shdwall,n + 1} \uparrow \right)\frac{H}{W}


Longwave radiation in the urban canyon is conserved as

.. math::
   :label: eq-0178

   {\overrightarrow{L}}_{uc} - \left( L_{uc} \uparrow - L_{atm} \downarrow \right) = 0

The total net longwave radiation for the urban canopy (road, walls, and roof) is

.. math::
   :label: eq-0179

   \overrightarrow{L} = W_{roof}{\overrightarrow{L}}_{roof} + \left( 1 - W_{roof} \right){\overrightarrow{L}}_{uc}

:numref:`fig-net-longwave-radiation` shows the net longwave radiation for urban surfaces for two different emissivity configurations. A positive net longwave means that the outgoing longwave exceeds the incoming longwave from the atmosphere. The net longwave radiation for the roof is independent of height to width ratio and increases with higher emissivity. The net longwave radiation for the road and walls decreases rapidly with increasing height to width ratio as more longwave radiation is trapped within the canyon. The walls have lower net longwave radiation than the road because their sky view factors are smaller. The two walls behave identically with respect to net longwave radiation as long as temperatures are the same. The sum of the net longwave radiation for road, sunlit wall, and shaded wall, after converting the wall fluxes to per unit ground area, is the canyon net longwave radiation. The net longwave radiation for the canyon increases slowly with increasing height to width ratio because of the larger surface area of the walls.

.. figure:: image10.png
   :width: 4.42708in
   :height: 6.08333in
   :name: fig-net-longwave-radiation

   Net longwave radiation (positive to the atmosphere) for urban surfaces for two different emissivity configurations. The atmospheric longwave radiation is :math:`L_{atm} \downarrow = 340` W m\ :sup:`-2` and the temperature of each surface is 292.16 K. Note that the wall fluxes (shaded and sunlit) are per unit wall area. The net longwave radiation for the canyon is the sum of road and wall fluxes after converting the walls fluxes to per unit ground area using the height to width ratio.

Solar Zenith Angle
------------------

The formulation for solar zenith angle is thoroughly documented in Oleson et al. (2010b) (see section 3.3) and does not differ for urban surfaces.

.. _section-2:

Heat and Momentum Fluxes
========================

The net radiation for the urban canopy (:math:`\overrightarrow{S} - \overrightarrow{L}`, where :math:`\overrightarrow{S}` is the net solar radiation absorbed by the urban canopy (section 2.5) and :math:`\overrightarrow{L}` is the net longwave radiation (section 2.7)) must be balanced by the sum of the turbulent and ground (storage) heat fluxes as

.. math::
   :label: eq-0180

   \overrightarrow{S} - \overrightarrow{L} = H + \lambda E + G

where :math:`H` is the sensible heat flux (W m\ :sup:`-2`), :math:`E` is the water vapor flux (kg m\ :sup:`-2` s\ :sup:`-1`), :math:`G` is the ground heat flux, and :math:`\lambda` is the latent heat of vaporization (or sublimation). The urban surfaces have unique radiative, thermal and hydrologic properties and environments. Thus, their sensible and latent heat fluxes are likely to be very different from each other. For example, the pervious road may have significant latent heat flux compared to the walls, which are assumed to be hydrologically inactive. Thus, the fluxes from individual urban surfaces must be modeled separately. However, CLM directly interacts with the atmospheric model at only the lowest atmospheric layer, which is well above the roof level of the urban model at the horizontal scales to be modeled. As a consequence, fluxes from individual urban surfaces must be combined to obtain the total sensible and latent heat flux to be provided to the atmospheric model. Allowing the urban surface fluxes to interact with each other through a bulk urban air mass is an acceptable approach analogous to the simulation of vegetated canopy fluxes (:numref:`fig-sensible-latent-heat-schematic`). This also allows for the solution of UCL air temperature and humidity, which are of interest in many applications. The approach shown in :numref:`fig-sensible-latent-heat-schematic` is slightly different from that of Masson (2000) in that here, fluxes from the roof interact directly with the UCL air whereas in Masson (2000) the roof and urban canyon are modeled as two independent sources of heat and moisture fluxes to the atmosphere. Here, we assume that the actual roofs are at various heights in the UCL and hence interact directly with the well-mixed UCL air.

.. figure:: image11.jpeg
   :width: 5.78125in
   :height: 5.76042in
   :name: fig-sensible-latent-heat-schematic

   Schematic diagram of sensible and latent heat fluxes for the urban canopy.

In general, the zonal :math:`\tau_{x}` and meridional :math:`\tau_{y}` momentum fluxes (kg m\ :sup:`-1` s\ :sup:`-2`), sensible heat flux :math:`H`, and water vapor flux :math:`E` between the atmosphere at reference height :math:`z_{atm,x}` (m) [where :math:`x` is height for wind (momentum) (:math:`m`), temperature (sensible heat) (:math:`h`), and humidity (water vapor) (:math:`w`); with zonal and meridional winds :math:`u_{atm}` and :math:`v_{atm}` (m s\ :sup:`-1`), potential temperature :math:`\theta_{atm}` (K), and specific humidity :math:`q_{atm}` (kg kg\ :sup:`-1`)] and a surface [with :math:`u_{s}`, :math:`v_{s}`, :math:`\theta_{s}`, and :math:`q_{s}`] are

.. math::
   :label: eq-0181

   \tau_{x} = - \rho_{atm}\frac{\left( u_{atm} - u_{s} \right)}{r_{am}}

.. math::
   :label: eq-0182

   \tau_{y} = - \rho_{atm}\frac{\left( v_{atm} - v_{s} \right)}{r_{am}}

.. math::
   :label: eq-0183

   H = - \rho_{atm}C_{p}\frac{\left( \theta_{atm} - \theta_{s} \right)}{r_{ah}}

.. math::
   :label: eq-0184

   E = - \rho_{atm}\frac{\left( q_{atm} - q_{s} \right)}{r_{aw}}

These fluxes are derived from Monin-Obukhov similarity theory developed for the inertial sub-layer (i.e., the nearly constant flux layer above the roughness sub-layer). In this derivation, :math:`u_{s}` and :math:`v_{s}` are defined to equal zero at height :math:`z_{0m} + d` (the apparent sink for momentum) so that :math:`r_{am}` is the aerodynamic resistance (s m\ :sup:`-1`) for momentum between the atmosphere at height :math:`z_{atm,m}` and the surface at height :math:`z_{0m} + d`. Thus, the momentum fluxes become

.. math::
   :label: eq-0185

   \tau_{x} = - \rho_{atm}\frac{u_{atm}}{r_{am}}

.. math::
   :label: eq-0186

   \tau_{y} = - \rho_{atm}\frac{v_{atm}}{r_{am}}

Likewise, :math:`\theta_{s}` and :math:`q_{s}` are defined at heights :math:`z_{0h} + d` and :math:`z_{0w} + d` (the apparent sinks for heat and water vapor, respectively). Consequently, :math:`r_{ah}` and :math:`r_{aw}` are the aerodynamic resistances (s m\ :sup:`-1`) to sensible heat and water vapor transfer between the atmosphere at heights :math:`z_{atm,h}` and :math:`z_{atm,w}` and the surface at heights :math:`z_{0h} + d` and :math:`z_{0w} + d`, respectively. The specific heat capacity of air :math:`C_{p}` (J kg\ :sup:`-1` K\ :sup:`-1`) is a constant (:numref:`table-physical-constants`). The atmospheric potential temperature used here is

.. math::
   :label: eq-0187

   \theta_{atm} = T_{atm} + \Gamma_{d}z_{atm,h}

where :math:`T_{atm}` is the air temperature (K) at height :math:`z_{atm,h}` and :math:`\Gamma_{d} = 0.0098` K m\ :sup:`-1` is the negative of the dry adiabatic lapse rate [this expression is first-order equivalent to :math:`\theta_{atm} = T_{atm}\left( \frac{P_{srf}}{P_{atm}} \right)^{\frac{R_{da}}{C_{p}}}` (Stull 1988), where :math:`P_{srf}` is the surface pressure (Pa), :math:`P_{atm}` is the atmospheric pressure (Pa), and :math:`R_{da}` is the gas constant for dry air (J kg\ :sup:`-1` K\ :sup:`-1`) (:numref:`table-physical-constants`)]. By definition, :math:`\theta_{s} = T_{s}`. The density of moist air (kg m\ :sup:`-3`) is

.. math::
   :label: eq-0188

   \rho_{atm} = \frac{P_{atm} - 0.378e_{atm}}{R_{da}T_{atm}}

where the atmospheric vapor pressure :math:`e_{atm}` (Pa) is derived from the atmospheric specific humidity :math:`q_{atm}`

.. math::
   :label: eq-0189

   e_{atm} = \frac{q_{atm}P_{atm}}{0.622 + 0.378q_{atm}}

Monin-Obukhov Similarity Theory
-------------------------------

The surface vertical kinematic fluxes of momentum :math:`\overline{u^{'}w^{'}}` and :math:`\overline{v^{'}w^{'}}` (m\ :sup:`2` s\ :sup:`-2`), sensible heat :math:`\overline{\theta^{'}w^{'}}` (K m s\ :sup:`-1`), and latent heat :math:`\overline{q^{'}w^{'}}` (kg kg\ :sup:`-1` m s\ :sup:`-1`), where :math:`u^{'}`, :math:`v^{'}`, :math:`w^{'}`, :math:`\theta^{'}`, and :math:`q^{'}` are zonal horizontal wind, meridional horizontal wind, vertical velocity, potential temperature, and specific humidity turbulent fluctuations about the mean, are defined from Monin-Obukhov similarity applied to the surface layer. This theory states that when scaled appropriately, the dimensionless mean horizontal wind speed, mean potential temperature, and mean specific humidity profile gradients depend on unique functions of :math:`\zeta = \frac{z - d}{L}` (Zeng et al. 1998) as

.. math::
   :label: eq-0190

   \frac{k(z - d)}{u_{*}}\frac{\partial|u|}{\partial z} = \varphi_{m}(\zeta)

.. math::
   :label: eq-0191

   \frac{k(z - d)}{\theta_{*}}\frac{\partial\theta}{\partial z} = \varphi_{h}(\zeta)

.. math::
   :label: eq-0192

   \frac{k(z - d)}{q_{*}}\frac{\partial q}{\partial z} = \varphi_{w}(\zeta)

where :math:`z` is height in the surface layer (m), :math:`d` is the displacement height (m), :math:`L` is the Monin-Obukhov length scale (m) that accounts for buoyancy effects resulting from vertical density gradients (i.e., the atmospheric stability), k is the von Karman constant (:numref:`table-physical-constants`), and :math:`|u|` is the atmospheric wind speed (m s\ :sup:`-1`). :math:`\varphi_{m}`, :math:`\varphi_{h}`, and :math:`\varphi_{w}` are universal (over any surface) similarity functions of :math:`\zeta` that relate the constant fluxes of momentum, sensible heat, and latent heat to the mean profile gradients of :math:`|u|`, :math:`\theta`, and :math:`q` in the surface layer. In neutral conditions, :math:`\varphi_{m} = \varphi_{h} = \varphi_{w} = 1`. The velocity (i.e., friction velocity) :math:`u_{*}` (m s\ :sup:`-1`), temperature :math:`\theta_{*}` (K), and moisture :math:`q_{*}` (kg kg\ :sup:`-1`) scales are

.. math::
   :label: eq-0193

   u_{*}^{2} = \sqrt{\left( \overline{u^{'}w^{'}} \right)^{2} + \left( \overline{v^{'}w^{'}} \right)^{2}} = \frac{|\tau|}{\rho_{atm}}

.. math::
   :label: eq-0194

   \theta_{*}u_{*} = - \overline{\theta^{'}w^{'}} = - \frac{H}{\rho_{atm}C_{p}}

.. math::
   :label: eq-0195

   q_{*}u_{*} = - \overline{q^{'}w^{'}} = - \frac{E}{\rho_{atm}}

where :math:`|\tau|` is the shearing stress (kg m\ :sup:`-1` s\ :sup:`-2`), with zonal and meridional components :math:`\overline{u^{'}w^{'}} = - \frac{\tau_{x}}{\rho_{atm}}` and :math:`\overline{v^{'}w^{'}} = - \frac{\tau_{y}}{\rho_{atm}}`, respectively, :math:`H` is the sensible heat flux (W m\ :sup:`-2`) and :math:`E` is the water vapor flux (kg m\ :sup:`-2` s\ :sup:`-1`).

The dimensionless length scale :math:`L` is the Monin-Obukhov length defined as

.. math::
   :label: eq-0196

   L = - \frac{u_{*}^{3}}{k\left( \frac{g}{\overline{\theta_{v,atm}}} \right)\theta_{v}^{'}w^{'}} = \frac{u_{*}^{2}\overline{\theta_{v,atm}}}{kg\theta_{v*}}

where :math:`g` is the acceleration of gravity (m s\ :sup:`-2`) (:numref:`table-physical-constants`), and :math:`\overline{\theta_{v,atm}} = \overline{\theta_{atm}}\left( 1 + 0.61q_{atm} \right)` is the reference virtual potential temperature. :math:`L > 0` indicates stable conditions. :math:`L < 0` indicates unstable conditions. :math:`L = \infty` for neutral conditions. The temperature scale :math:`\theta_{v*}` is defined as

.. math::
   :label: eq-0197

   \theta_{v*}u_{*} = \left\lbrack \theta_{*}\left( 1 + 0.61q_{atm} \right) + 0.61\overline{\theta_{atm}}q_{*} \right\rbrack u_{*}

where :math:`\overline{\theta_{atm}}` is the atmospheric potential temperature.

Following Panofsky and Dutton (1984), the differential equations for :math:`\varphi_{m}(\zeta)`, :math:`\varphi_{h}(\zeta)`, and :math:`\varphi_{w}(\zeta)` can be integrated formally without commitment to their exact forms. Integration between two arbitrary heights in the surface layer :math:`z_{2}` and :math:`z_{1}` (:math:`z_{2} > z_{1}`) with horizontal winds :math:`|u|_{1}` and :math:`|u|_{2}`, potential temperatures :math:`\theta_{1}` and :math:`\theta_{2}`, and specific humidities :math:`q_{1}` and :math:`q_{2}` results in

.. math::
   :label: eq-0198

   |u|_{2} - |u|_{1} = \frac{u_{*}}{k}\left\lbrack \ln\left( \frac{z_{2} - d}{z_{1} - d} \right) - \psi_{m}\left( \frac{z_{2} - d}{L} \right) + \psi_{m}\left( \frac{z_{1} - d}{L} \right) \right\rbrack

.. math::
   :label: eq-0199

   \theta_{2} - \theta_{1} = \frac{\theta_{*}}{k}\left\lbrack \ln\left( \frac{z_{2} - d}{z_{1} - d} \right) - \psi_{h}\left( \frac{z_{2} - d}{L} \right) + \psi_{h}\left( \frac{z_{1} - d}{L} \right) \right\rbrack

.. math::
   :label: eq-0200

   q_{2} - q_{1} = \frac{q_{*}}{k}\left\lbrack \ln\left( \frac{z_{2} - d}{z_{1} - d} \right) - \psi_{w}\left( \frac{z_{2} - d}{L} \right) + \psi_{w}\left( \frac{z_{1} - d}{L} \right) \right\rbrack

The functions :math:`\psi_{m}(\zeta)`, :math:`\psi_{h}(\zeta)`, and :math:`\psi_{w}(\zeta)` are defined as

.. math::
   :label: eq-0201

   \psi_{m}(\zeta) = \int_{\frac{z_{0m}}{L}}^{\zeta}{\frac{\left\lbrack 1 - \varphi_{m}(x) \right\rbrack}{x}dx}

.. math::
   :label: eq-0202

   \psi_{h}(\zeta) = \int_{\frac{z_{0h}}{L}}^{\zeta}{\frac{\left\lbrack 1 - \varphi_{h}(x) \right\rbrack}{x}dx}

.. math::
   :label: eq-0203

   \psi_{w}(\zeta) = \int_{\frac{z_{0w}}{L}}^{\zeta}{\frac{\left\lbrack 1 - \varphi_{w}(x) \right\rbrack}{x}dx}

where :math:`z_{0m}`, :math:`z_{0h}`, and :math:`z_{0w}` are the roughness lengths (m) for momentum, sensible heat, and water vapor, respectively.

Defining the surface values

:math:`|u|_{1} = 0\text{ at }z_{1} = z_{0m} + d,`

:math:`\theta_{1} = \theta_{s}\text{ at }z_{1} = z_{0h} + d,\text{ and}`

:math:`q_{1} = q_{s}\text{ at }z_{1} = z_{0w} + d,`

and the atmospheric values at :math:`z_{2} = z_{atm,x}`

.. math::
   :label: eq-0204

   |u|_{2} = V_{a}\text{= }\sqrt{u_{atm}^{2} + v_{atm}^{2} + U_{c}^{2}} \geq 1,

:math:`\theta_{2} = \theta_{atm}\text{, and}`

:math:`q_{2} = q_{atm}\text{, }`

the integral forms of the flux-gradient relations are

.. math::
   :label: eq-0205

   V_{a} = \frac{u_{*}}{k}\left\lbrack \ln\left( \frac{z_{atm,m} - d}{z_{0m}} \right) - \psi_{m}\left( \frac{z_{atm,m} - d}{L} \right) + \psi_{m}\left( \frac{z_{0m}}{L} \right) \right\rbrack

.. math::
   :label: eq-0206

   \theta_{atm} - \theta_{s} = \frac{\theta_{*}}{k}\left\lbrack \ln\left( \frac{z_{atm,h} - d}{z_{0h}} \right) - \psi_{h}\left( \frac{z_{atm,h} - d}{L} \right) + \psi_{h}\left( \frac{z_{0h}}{L} \right) \right\rbrack

.. math::
   :label: eq-0207

   q_{atm} - q_{s} = \frac{q_{*}}{k}\left\lbrack \ln\left( \frac{z_{atm,w} - d}{z_{0w}} \right) - \psi_{w}\left( \frac{z_{atm,w} - d}{L} \right) + \psi_{w}\left( \frac{z_{0w}}{L} \right) \right\rbrack

The constraint :math:`V_{a} \geq 1` is required simply for numerical reasons to prevent :math:`H` and :math:`E` from becoming small with small wind speeds. The convective velocity :math:`U_{c}` accounts for the contribution of large eddies in the convective boundary layer to surface fluxes as follows



.. math::
   :label: eq-0041

   U_{c} = 0 \zeta \geq 0 (\text{stable})
   U_{c} = \beta w_{*} \zeta < 0 (\text{unstable})



where :math:`w_{*}` is the convective velocity scale

:math:`w_{*} = \left( \frac{- gu_{*}\theta_{v*}z_{i}}{\overline{\theta_{v,atm}}} \right)^{\frac{1}{3}}`,

:math:`z_{i} = 1000` is the convective boundary layer height (m), and :math:`\beta = 1`.

The momentum flux gradient relations are (Zeng et al. 1998)



.. math::
   :label: eq-0042

   \varphi_{m}(\zeta) = 0.7k^{\frac{2}{3}}( - \zeta)^{\frac{1}{3}} \text{for }\zeta < - 1.574\ (\text{very unstable})
   

\ 

.. math::
   :label: eq-0043

   \begin{cases}
   \varphi_{m}(\zeta) = (1 - 16\zeta)^{- \frac{1}{4}} & \text{for -1.574} \leq \zeta < 0\ (\text{unstable}) \\
   \varphi_{m}(\zeta) = 1 + 5\zeta & \text{for }0 \leq \zeta \leq 1\ (\text{stable})
   \end{cases}

\ :math:`\varphi_{m}(\zeta) = 5 + \zeta \text{for }\zeta\text{>1 }(\text{very stable}).`

The sensible and latent heat flux gradient relations are (Zeng et al. 1998)



.. math::
   :label: eq-0044

   \varphi_{h}(\zeta) = \varphi_{w}(\zeta) = 0.9k^{\frac{4}{3}}( - \zeta)^{\frac{- 1}{3}} \text{for }\zeta < - 0.465\ (\text{very unstable})
   

\ 

.. math::
   :label: eq-0045

   \begin{cases}
   \varphi_{h}(\zeta) = \varphi_{w}(\zeta) = (1 - 16\zeta)^{- \frac{1}{2}} & \text{for -0.465} \leq \zeta < 0\ (\text{unstable}) \\
   \varphi_{h}(\zeta) = \varphi_{w}(\zeta) = 1 + 5\zeta & \text{for }0 \leq \zeta \leq 1\ (\text{stable})
   \end{cases}

\ :math:`\varphi_{h}(\zeta) = \varphi_{w}(\zeta) = 5 + \zeta \text{for }\zeta\text{>1 }(\text{very stable}).`

To ensure continuous functions of :math:`\varphi_{m}(\zeta)`, :math:`\varphi_{h}(\zeta)`, and :math:`\varphi_{w}(\zeta)`, the simplest approach (i.e., without considering any transition regimes) is to match the relations for very unstable and unstable conditions at :math:`\zeta_{m} = - 1.574` for :math:`\varphi_{m}(\zeta)` and :math:`\zeta_{h} = \zeta_{w} = - 0.465` for :math:`\varphi_{h}(\zeta) = \varphi_{w}(\zeta)` (Zeng et al. 1998). The flux gradient relations can be integrated to yield wind profiles for the following conditions:

Very unstable :math:`(\zeta < - 1.574)`

.. math::
   :label: eq-0208

   V_{a} = \frac{u_{*}}{k}\left\{ \left\lbrack \ln\frac{\zeta_{m}L}{z_{0m}} - \psi_{m}\left( \zeta_{m} \right) \right\rbrack + 1.14\left\lbrack ( - \zeta)^{\frac{1}{3}} - \left( - \zeta_{m} \right)^{\frac{1}{3}} \right\rbrack + \psi_{m}\left( \frac{z_{0m}}{L} \right) \right\}

Unstable :math:`( - 1.574 \leq \zeta < 0)`

.. math::
   :label: eq-0209

   V_{a} = \frac{u_{*}}{k}\left\{ \left\lbrack \ln\frac{z_{atm,m} - d}{z_{0m}} - \psi_{m}(\zeta) \right\rbrack + \psi_{m}\left( \frac{z_{0m}}{L} \right) \right\}

Stable :math:`(0 \leq \zeta \leq 1)`

.. math::
   :label: eq-0210

   V_{a} = \frac{u_{*}}{k}\left\{ \left\lbrack \ln\frac{z_{atm,m} - d}{z_{0m}} + 5\zeta \right\rbrack - 5\frac{z_{0m}}{L} \right\}

Very stable :math:`(\zeta > 1)`

.. math::
   :label: eq-0211

   V_{a} = \frac{u_{*}}{k}\left\{ \left\lbrack \ln\frac{L}{z_{0m}} + 5 \right\rbrack + \left\lbrack 5\ln\zeta + \zeta - 1 \right\rbrack - 5\frac{z_{0m}}{L} \right\}

where

.. math::
   :label: eq-0212

   \psi_{m}(\zeta) = 2\ln\left( \frac{1 + x}{2} \right) + \ln\left( \frac{1 + x^{2}}{2} \right) - 2\tan^{- 1}x + \frac{\pi}{2}

and :math:`x = (1 - 16\zeta)^{\frac{1}{4}}`.

The potential temperature profiles are:

Very unstable :math:`(\zeta < - 0.465)`

.. math::
   :label: eq-0213

   \theta_{atm} - \theta_{s} = \frac{\theta_{*}}{k}\left\{ \left\lbrack \ln\frac{\zeta_{h}L}{z_{0h}} - \psi_{h}\left( \zeta_{h} \right) \right\rbrack + 0.8\left\lbrack \left( - \zeta_{h} \right)^{\frac{- 1}{3}} - ( - \zeta)^{\frac{- 1}{3}} \right\rbrack + \psi_{h}\left( \frac{z_{0h}}{L} \right) \right\}

Unstable :math:`( - 0.465 \leq \zeta < 0)`

.. math::
   :label: eq-0214

   \theta_{atm} - \theta_{s} = \frac{\theta_{*}}{k}\left\{ \left\lbrack \ln\frac{z_{atm,h} - d}{z_{0h}} - \psi_{h}(\zeta) \right\rbrack + \psi_{h}\left( \frac{z_{0h}}{L} \right) \right\}

Stable :math:`(0 \leq \zeta \leq 1)`

.. math::
   :label: eq-0215

   \theta_{atm} - \theta_{s} = \frac{\theta_{*}}{k}\left\{ \left\lbrack \ln\frac{z_{atm,h} - d}{z_{0h}} + 5\zeta \right\rbrack - 5\frac{z_{0h}}{L} \right\}

Very stable :math:`(\zeta > 1)`

.. math::
   :label: eq-0216

   \theta_{atm} - \theta_{s} = \frac{\theta_{*}}{k}\left\{ \left\lbrack \ln\frac{L}{z_{0h}} + 5 \right\rbrack + \left\lbrack 5\ln\zeta + \zeta - 1 \right\rbrack - 5\frac{z_{0h}}{L} \right\}

The specific humidity profiles are:

Very unstable :math:`(\zeta < - 0.465)`

.. math::
   :label: eq-0217

   q_{atm} - q_{s} = \frac{q_{*}}{k}\left\{ \left\lbrack \ln\frac{\zeta_{w}L}{z_{0w}} - \psi_{w}\left( \zeta_{w} \right) \right\rbrack + 0.8\left\lbrack \left( - \zeta_{w} \right)^{\frac{- 1}{3}} - ( - \zeta)^{\frac{- 1}{3}} \right\rbrack + \psi_{w}\left( \frac{z_{0w}}{L} \right) \right\}

Unstable :math:`( - 0.465 \leq \zeta < 0)`

.. math::
   :label: eq-0218

   q_{atm} - q_{s} = \frac{q_{*}}{k}\left\{ \left\lbrack \ln\frac{z_{atm,w} - d}{z_{0w}} - \psi_{w}(\zeta) \right\rbrack + \psi_{w}\left( \frac{z_{0w}}{L} \right) \right\}

Stable :math:`(0 \leq \zeta \leq 1)`

.. math::
   :label: eq-0219

   q_{atm} - q_{s} = \frac{q_{*}}{k}\left\{ \left\lbrack \ln\frac{z_{atm,w} - d}{z_{0w}} + 5\zeta \right\rbrack - 5\frac{z_{0w}}{L} \right\}

Very stable :math:`(\zeta > 1)`

.. math::
   :label: eq-0220

   q_{atm} - q_{s} = \frac{q_{*}}{k}\left\{ \left\lbrack \ln\frac{L}{z_{0w}} + 5 \right\rbrack + \left\lbrack 5\ln\zeta + \zeta - 1 \right\rbrack - 5\frac{z_{0w}}{L} \right\}

where

.. math::
   :label: eq-0221

   \psi_{h}(\zeta) = \psi_{w}(\zeta) = 2\ln\left( \frac{1 + x^{2}}{2} \right)

Using the definitions of :math:`u_{*}`, :math:`\theta_{*}`, and :math:`q_{*}`, an iterative solution of these equations can be used to calculate the surface momentum, sensible heat, and water vapor flux using atmospheric and surface values for :math:`|u|`, :math:`\theta`, and :math:`q` except that :math:`L` depends on :math:`u_{*}`, :math:`\theta_{*}`, and :math:`q_{*}`. However, the bulk Richardson number

.. math::
   :label: eq-0222

   R_{iB} = \frac{\theta_{v,atm} - \theta_{v,s}}{\overline{\theta_{v,atm}}}\frac{g\left( z_{atm,m} - d \right)}{V_{a}^{2}}

is related to :math:`\zeta` (Arya 2001) as

.. math::
   :label: eq-0223

   R_{iB} = \zeta\left\lbrack \ln\left( \frac{z_{atm,h} - d}{z_{0h}} \right) - \psi_{h}(\zeta) \right\rbrack\left\lbrack \ln\left( \frac{z_{atm,m} - d}{z_{0m}} \right) - \psi_{m}(\zeta) \right\rbrack^{- 2}

Using :math:`\varphi_{h} = \varphi_{m}^{2} = (1 - 16\zeta)^{- \frac{1}{2}}` for unstable conditions and :math:`\varphi_{h} = \varphi_{m} = 1 + 5\zeta` for stable conditions to determine :math:`\psi_{m}(\zeta)` and :math:`\psi_{h}(\zeta)`, the inverse relationship :math:`\zeta = f\left( R_{iB} \right)` can be solved to obtain a first guess for :math:`\zeta` and thus :math:`L` from



.. math::
   :label: eq-0046

   \zeta = \frac{R_{iB}\ln\left( \frac{z_{atm,m} - d}{z_{0m}} \right)}{1 - 5\min\left( R_{iB},0.19 \right)} 0.01 \leq \zeta \leq 2 \text{for }R_{iB} \geq 0\ (\text{neutral or stable})
   \zeta = R_{iB}\ln\left( \frac{z_{atm,m} - d}{z_{0m}} \right)  - 100 \leq \zeta \leq - 0.01 \text{for }R_{iB} < 0\ (\text{unstable})


Upon iteration, the following is used to determine :math:`\zeta` and thus :math:`L`

.. math::
   :label: eq-0224

   \zeta = \frac{\left( z_{atm,m} - d \right)kg\theta_{v*}}{u_{*}^{2}\overline{\theta_{v,atm}}}

where



.. math::
   :label: eq-0047

   0.01 \leq \zeta \leq 2 \text{for }\zeta \geq 0\ (\text{neutral or stable})
   \text{-100} \leq \zeta \leq \text{-0.01} \text{for }\zeta < 0\ (\text{unstable})


The momentum, sensible heat, and water vapor fluxes between the surface and the atmosphere can also be written in the form

.. math::
   :label: eq-0225

   \tau_{x} = - \rho_{atm}\frac{\left( u_{atm} - u_{s} \right)}{r_{am}}

.. math::
   :label: eq-0226

   \tau_{y} = - \rho_{atm}\frac{\left( v_{atm} - v_{s} \right)}{r_{am}}

.. math::
   :label: eq-0227

   H = - \rho_{atm}C_{p}\frac{\left( \theta_{atm} - \theta_{s} \right)}{r_{ah}}

.. math::
   :label: eq-0228

   E = - \rho_{atm}\frac{\left( q_{atm} - q_{s} \right)}{r_{aw}}

where :math:`r_{am}`, :math:`r_{ah}`, and :math:`r_{aw}` are the aerodynamic resistances for momentum, sensible heat and latent heat, respectively (s m\ :sup:`-1`).

Sensible and Latent Heat and Momentum Fluxes
--------------------------------------------

The solution for the heat and momentum fluxes is presented in roughly the order in which the equations are solved in the Fortran code.

Roughness Length and Displacement Height
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The roughness length and displacement height for the urban canopy are needed. Grimmond and Oke (1999) review approaches to calculate these parameters from morphometric methods. Here, we use the Macdonald et al. (1998) approach, which appears to be a reasonable compromise between minimizing input requirements and yielding acceptable results. The subscript "canopy" is used to distinguish between an aerodynamic parameter for the urban canopy versus a parameter for an individual urban surface (e.g., roof).

The canopy displacement height :math:`d_{canopy}` (m) is

.. math::
   :label: eq-0229

   d_{canopy} = H\left\lbrack 1 + \alpha^{- \lambda_{P}}(\lambda_{P} - 1) \right\rbrack

where :math:`H` is the canyon (roof) height (m) (:numref:`table-input-data`), :math:`\alpha = \text{4.43}` is an empirical coefficient, and :math:`\lambda_{P}` is the plan area index. The plan area index :math:`\lambda_{p}` is

.. math::
   :label: eq-0230

   \lambda_{p} = \frac{\frac{H}{W}}{\frac{H}{W} + 1}

where :math:`\frac{H}{W}` is the height to width ratio of the urban canyon (:numref:`table-input-data`).

The canopy roughness length :math:`z_{0m,canopy}` (m) for momentum is

.. math::
   :label: eq-0231

   z_{0m,canopy} = H\left( 1 - \frac{d_{canopy}}{H} \right)\exp\left\{ - \left\lbrack 0.5B\frac{C_{D}}{k^{2}}\left( 1 - \frac{d_{canopy}}{H} \right)\lambda_{F} \right\rbrack^{- 0.5} \right\}

where :math:`B = 1` is a correction to the drag coefficient to account for variable obstacle shapes and flow conditions, :math:`C_{D} = 1.2` is the depth-integrated mean drag coefficient for surface-mounted cubes in a shear flow, :math:`k` is the von Karman constant, and :math:`\lambda_{F}` is the frontal area index. The frontal area index :math:`\lambda_{F}` is

.. math::
   :label: eq-0232

   \lambda_{F} = \left( 1 - \lambda_{P} \right)\left( \frac{H}{W} \right)\sqrt{\frac{B_{L}\lambda_{P}}{B_{S}}}

where :math:`\frac{B_{S}}{B_{L}}` is the building shortside to longside ratio (here set equal to :math:`\lambda_{P}`).

Several checks are made to ensure that the derived aerodynamic parameters are consistent with the canyon structure and atmospheric forcing. First, the canyon height :math:`H` minus the canopy displacement height :math:`d_{canopy}` must be greater than the canopy roughness length :math:`z_{0m,canopy}`. Second, the atmospheric wind forcing height :math:`z_{atm,m}` (:numref:`table-atm-input`) minus the canopy displacement height :math:`d_{canopy}` must be greater than the canopy roughness length\ :math:`z_{0m,canopy}`. Note that :math:`z_{0m,canopy} = z_{0h,canopy} = z_{0w,canopy}` and :math:`z_{atm} = z_{atm}^{'} + z_{0,canopy} + z_{d,canopy}` (:numref:`table-atm-input`) where :math:`z_{atm}^{'}` is the reference height from the atmospheric model.

Wind Speed in the Urban Canyon
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Following Masson (2000) and Lemonsu et al. (2004), the wind speed in the canyon is the combination of the mean horizontal canyon wind :math:`U_{can}` (m s\ :sup:`-1`) and the turbulent (vertical) wind :math:`W_{can}` (m s\ :sup:`-1`)

.. math::
   :label: eq-0233

   U_{ac} = \sqrt{{U_{can}}^{2} + {W_{can}}^{2}}

To calculate the horizontal wind speed in the canyon :math:`U_{can}` (m s\ :sup:`-1`), a horizontal wind speed at the top of the canyon is derived by assuming a logarithmic wind profile from the atmospheric reference height to the canyon top. The wind is then extrapolated to a height inside the canyon using an exponential profile. For skimming flow (:math:`\frac{H}{W} \geq 1`) (Oke 1987), a zero :math:`U_{can}` is assumed when the mean flow is perpendicular to the canyon orientation. After integration over 360° (to account for all street orientations),

.. math::
   :label: eq-0234

   U_{can} = V_{r}\frac{2}{\pi}\frac{\ln\left( \frac{H - d_{canopy}}{z_{0m,canopy}} \right)}{\ln\left( \frac{z_{atm,m} - d_{canopy}}{z_{0m,canopy}} \right)}\exp\left\lbrack - 0.5\left( \frac{H}{W} \right)\left( 1 - \frac{H_{w}}{H} \right) \right\rbrack

where :math:`H_{w}` is the height at which the wind speed is estimated (:numref:`table-input-data`). For isolated roughness flow (:math:`\frac{H}{W < \text{0.5}}`), the wind speed in the canyon is assumed to be independent of the orientation of the mean atmospheric flow above the canyon level,

.. math::
   :label: eq-0235

   U_{can} = V_{r}\frac{\ln\left( \frac{H - d_{canopy}}{z_{0m,canopy}} \right)}{\ln\left( \frac{z_{atm,m} - d_{canopy}}{z_{0m,canopy}} \right)}\exp\left( - 0.5\left( \frac{H}{W} \right)\left( 1 - \frac{H_{w}}{H} \right) \right)

For wake interference flow (:math:`\frac{0.5 \leq H}{W < \text{1.0}}`),



.. math::
   :label: eq-0048

   U_{can} = V_{r}\left\lbrack 1 + 2\left( \frac{2}{\pi} - 1 \right)\left( \frac{H}{W} - \frac{1}{2} \right) \right\rbrack\frac{\ln\left( \frac{H - d_{canopy}}{z_{0m,canopy}} \right)}{\ln\left( \frac{z_{atm,m} - d_{canopy}}{z_{0m,canopy}} \right)}
   \times \exp\left( - 0.5\left( \frac{H}{W} \right)\left( 1 - \frac{H_{w}}{H} \right) \right)


The magnitude of the reference level atmospheric wind is

.. math::
   :label: eq-0236

   V_{r} = \sqrt{u_{atm}^{2} + v_{atm}^{2}} \geq 1

where zonal and meridional winds :math:`u_{atm}` and :math:`v_{atm}` (m s\ :sup:`-1`) are at height :math:`z_{atm,m}`. The turbulent (vertical) wind :math:`W_{can}` (m s\ :sup:`-1`) is assumed to be equal to the friction velocity (Masson 2000), which is determined from the solution for turbulent fluxes (section 3.2.3).

Iterative Solution for Urban Canopy Air Temperature and Humidity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because of the interdependence between fluxes, aerodynamic resistances, and canyon air temperature and humidity, an iterative solution for the UCL air is devised.

An initial guess for the wind speed :math:`V_{a}` (equation :eq:`eq-0204` is obtained assuming an initial convective velocity :math:`U_{c} = 0` m s\ :sup:`-1` for stable conditions and :math:`U_{c} = 0.5` for unstable conditions. Stable conditions (:math:`\theta_{v,atm} - \theta_{v,s} \geq 0`) and unstable conditions (:math:`\theta_{v,atm} - \theta_{v,s} < 0`) are evaluated from the difference in virtual potential air temperature between the reference height and the surface where

.. math::
   :label: eq-0237

   \theta_{v,atm} - \theta_{v,s} = \left( \theta_{atm} - \theta_{s} \right)\left( 1 + 0.61q_{atm} \right) + 0.61\overline{\theta_{atm}}\left( q_{atm} - q_{s} \right)

Here, :math:`\theta_{s} = T_{ac}` and :math:`q_{s} = q_{ac}` where :math:`T_{ac}` is the air temperature in the UCL (K) and :math:`q_{ac}` is the specific humidity in the UCL (kg kg\ :sup:`-1`) (:numref:`fig-sensible-latent-heat-schematic`). The air temperature and specific humidity from the previous time step are used. The temperature :math:`\theta_{atm}` is defined by equation :eq:`eq-0187` , :math:`\overline{\theta_{atm}}` is the atmospheric potential temperature (:numref:`table-atm-input`), and :math:`q_{atm}` is the atmospheric specific humidity (kg kg\ :sup:`-1`) (:numref:`table-atm-input`). An initial guess for the Monin-Obukhov length :math:`L` is obtained from the bulk Richardson number using equations :eq:`eq-0222` and :eq:`eq-0046`.

The iterative solution begins with the friction velocity :math:`u_{*}`, potential temperature scale :math:`\theta_{*}`, and humidity scale :math:`q_{*}` being calculated from equations :eq:`eq-0208`-:eq:`eq-0221`. Now that the friction velocity has been determined, the wind in the urban canopy, :math:`U_{ac}`, is calculated from equation :eq:`eq-0233`. The aerodynamic resistances (s m\ :sup:`-1`) to momentum, sensible heat, and latent heat transfer between the UCL air and the atmosphere are

.. math::
   :label: eq-0238

   r_{am} = \frac{V_{a}}{u_{*}^{2}} = \frac{1}{k^{2}V_{a}}\left\lbrack \ln\left( \frac{z_{atm,m} - d_{canopy}}{z_{0m,canopy}} \right) - \psi_{m}\left( \frac{z_{atm,m} - d_{canopy}}{L} \right) + \psi_{m}\left( \frac{z_{0m,canopy}}{L} \right) \right\rbrack^{2}



.. math::
   :label: eq-0049

   r_{ah} = \frac{\theta_{atm} - \theta_{s}}{\theta_{*}u_{*}} = \frac{1}{k^{2}V_{a}}\left\lbrack \ln\left( \frac{z_{atm,m} - d_{canopy}}{z_{0m,canopy}} \right) - \psi_{m}\left( \frac{z_{atm,m} - d_{canopy}}{L} \right) + \psi_{m}\left( \frac{z_{0m,canopy}}{L} \right) \right\rbrack
   \left\lbrack \ln\left( \frac{z_{atm,h} - d_{canopy}}{z_{0h,\mspace{6mu} canopy}} \right) - \psi_{h}\left( \frac{z_{atm,h} - d_{canopy}}{L} \right) + \psi_{h}\left( \frac{z_{0h,canopy}}{L} \right) \right\rbrack



.. math::
   :label: eq-0050

   r_{aw} = \frac{q_{atm} - q_{s}}{q_{*}u_{*}} = \frac{1}{k^{2}V_{a}}\left\lbrack \ln\left( \frac{z_{atm,m} - d_{canopy}}{z_{0m,canopy}} \right) - \psi_{m}\left( \frac{z_{atm,m} - d_{canopy}}{L} \right) + \psi_{m}\left( \frac{z_{0m,canopy}}{L} \right) \right\rbrack
   \left\lbrack \ln\left( \frac{z_{atm,w} - d_{canopy}}{z_{0w,canopy}} \right) - \psi_{w}\left( \frac{z_{atm,w} - d_{canopy}}{L} \right) + \psi_{w}\left( \frac{z_{0w,canopy}}{L} \right) \right\rbrack


The resistances to sensible heat and latent heat transfer between canyon surfaces (roof, sunlit and shaded wall, pervious and impervious road) and the UCL depend only on canyon wind speed following Masson (2000). Thus, the surface resistances, :math:`r_{s,roof}`, :math:`r_{s,sunwall}`, :math:`r_{s,shdwall}`, :math:`r_{s,prvrd}`, :math:`r_{s,imprvrd}`, (s m\ :sup:`-1`) are identical and are determined from (Rowley et al. 1930)

.. math::
   :label: eq-0239

   r_{s,u} = \frac{\rho_{atm}C_{p}}{11.8 + 4.2U_{ac}}

The UCL air temperature and specific humidity are determined by solving the following systems of equations. For :math:`T_{ac}`

:math:`H_{roof} = - \rho_{atm}C_{p}\frac{T_{ac} - T_{g,roof}}{r_{s,roof}}`,

:math:`H_{prvrd} = - \rho_{atm}C_{p}\frac{T_{ac} - T_{g,prvrd}}{r_{s,prvrd}}`,

:math:`H_{imprvrd} = - \rho_{atm}C_{p}\frac{T_{ac} - T_{g,imprvrd}}{r_{s,imprvrd}}`,

:math:`H_{sunwall} = - \rho_{atm}C_{p}\frac{T_{ac} - T_{g,sunwall}}{r_{s,sunwall}}`,

:math:`H_{shdwall} = - \rho_{atm}C_{p}\frac{T_{ac} - T_{g,shdwall}}{r_{s,shdwall}}`,

.. math::
   :label: eq-0240

   H = - \rho_{atm}C_{p}\frac{\theta_{atm} - T_{ac}}{r_{ah}} = W_{roof}H_{roof} + \left( 1 - W_{roof} \right) \times \left\lbrack f_{prvrd}H_{prvrd} + \left( 1 - f_{prvrd} \right)H_{imprvrd} + \frac{H}{W}H_{sunwall} + \frac{H}{W}H_{shdwall} \right\rbrack

where :math:`H` is sensible heat flux (W m\ :sup:`-2`) and :math:`T_{g}` is the surface temperature of each urban surface. The term :math:`W_{roof}` is the relative contribution of roof fluxes to the total urban landunit flux (:numref:`table-input-data`). The term :math:`1 - W_{roof}` is then the relative contribution of the canyon to the total urban landunit flux. The term :math:`f_{prvrd}` is the fraction of road that is pervious (:numref:`table-input-data`) and the term :math:`1 - f_{prvrd}` is the fraction of the road that is impervious. Note that the factor :math:`\frac{H}{W}` for the sunwall and shadewall converts the flux from watts per meter squared of surface area to watts per meter squared of ground area.

In Oleson et al. (2008a), an additional heat flux :math:`H_{wasteheat}`, the sensible heat flux from waste heat generated by space heating and air conditioning, was included in Equations (3.69)-:eq:`eq-0240` *.* However, if this flux is large enough, the numerical solution may become unstable because of the canopy air has no heat capacity and the heat capacities of the roofs and walls are relatively small. Instead, this heat flux is added to the net heat flux for the canyon floor (section 4.1).

equation :eq:`eq-0240` can be solved for the UCL air temperature as



.. math::
   :label: eq-0051

   T_{ac} = \frac{\left( \begin{aligned}
   & c_{a}^{h}\theta_{atm} + c_{roof}T_{g,roof} + c_{prvrd}T_{g,prvrd} + \\
   & c_{imprvrd}T_{g,imprvrd} + c_{sunwall}T_{g,sunwall} + c_{shdwall}T_{g,shdwall}
   \end{aligned} \right)}{c_{a}^{h} + c_{roof} + c_{prvrd} + c_{imprvrd} + c_{sunwall} + c_{shdwall}}



where :math:`c_{a}^{h}` is the sensible heat conductance from the UCL to the atmosphere (:math:`\frac{1}{r_{ah}}`), and :math:`c_{roof}`, :math:`c_{prvrd}`, :math:`c_{imprvrd}`, :math:`c_{sunwall}`, and :math:`c_{shdwall}` are the weighted heat conductances from urban surfaces to UCL air [:math:`\frac{W_{roof}}{r_{s,roof}}`, :math:`\frac{W_{prvrd}}{r_{s,prvrd}}`, :math:`\frac{W_{imprvrd}}{r_{s,imprvrd}}`, :math:`\frac{W_{sunwall}}{r_{s,sunwall}}`, :math:`\frac{W_{shawall}}{r_{s,shawall}}`, respectively, where :math:`W_{prvrd} = \left( 1 - W_{roof} \right)f_{prvrd}`, :math:`W_{imprvrd} = \left( 1 - W_{roof} \right)\left( 1 - f_{prvrd} \right)`, :math:`W_{sunwall} = \left( 1 - W_{roof} \right)\left( \frac{H}{W} \right)`, and :math:`W_{shdwall} = \left( 1 - W_{roof} \right)\left( \frac{H}{W} \right)`].

Similarly, the system of equations for the UCL air specific humidity, :math:`q_{ac}`, is

:math:`E_{roof} = - \rho_{atm}\frac{f_{wet,roof}\left( q_{ac} - q_{g,roof} \right)}{r_{s,roof}}`,

:math:`E_{prvrd} = - \rho_{atm}\frac{q_{ac} - q_{g,prvrd}}{r_{s,prvrd}}`,

:math:`E_{imprvrd} = - \rho_{atm}\frac{f_{wet,imprvrd}\left( q_{ac} - q_{g,imprvrd} \right)}{r_{s,imprvrd}}`,

:math:`E_{sunwall} = 0`,

:math:`E_{shdwall} = 0`,



.. math::
   :label: eq-0052

   E = - \rho_{atm}\frac{q_{atm} - q_{ac}}{r_{aw}}
   = W_{roof}E_{roof} + \left( 1 - W_{roof} \right)\left\lbrack f_{prvrd}E_{prvrd} + \left( 1 - f_{prvrd} \right)E_{imprvrd} \right\rbrack



where :math:`E` is water vapor flux (kg m\ :sup:`-2` s\ :sup:`-1`) and :math:`q_{g}` is the specific humidity at each urban surface (kg kg\ :sup:`-1`). Note that the latent heat flux from the sunlit and shaded walls is zero. The term :math:`f_{wet}` is the fraction of the roof or impervious road surface that is wet. If there is dew formation (:math:`q_{ac} - q_{g} > 0`), then :math:`f_{wet} = 1`. If there is snow on the surface (:math:`z_{sno} > 0`), :math:`f_{wet}` is determined from the snow depth :math:`z_{sno}` as

.. math::
   :label: eq-0241

   f_{wet} = \frac{z_{sno}}{0.05} \leq 1

In the absence of snow,

.. math::
   :label: eq-0242

   f_{wet} = \left( \frac{w_{liq,snl + 1} + w_{ice,snl + 1}}{w_{pond,\max}}()^{\frac{2}{3}} \right)

where :math:`w_{liq,snl + 1}` and :math:`w_{ice,snl + 1}` are the mass of ice and liquid water (kg m\ :sup:`-2`) stored on top of the urban surface and :math:`w_{pond,\max}` is the maximum amount of water that the surface can hold (Chapter 5). This latter formulation is analogous to the treatment of the wetted fraction of the vegetated canopy in CLM (Oleson et al. 2004).

In equations (3.76) and (3.78), the specific humidity of the roof and the impervious road surfaces, :math:`q_{g,roof}` and :math:`q_{g,imprvrd}`, is set to the saturated specific humidity evaluated at their respective surface temperatures, :math:`q_{sat}^{T_{g,roof}}` and :math:`q_{sat}^{T_{g,imprvrd}}` (section 3.3).

As noted in section 1.1.3, a simplified bulk parameterization approach is used to represent evaporation from the pervious surface. The pervious road specific humidity, :math:`q_{g,prvrd}`, is evaluated as a function of the wetness of the soil column. This allows all of the soil moisture to potentially be available for evaporation. The specific humidity is

:math:`q_{g,prvrd} = \alpha q_{sat}^{T_{g}}`,

where :math:`q_{sat}^{T_{g}}` is the saturated specific humidity at the surface temperature :math:`T_{g}` (section 4.1). The factor :math:`\alpha` is a weighted combination of values for the soil column and snow

.. math::
   :label: eq-0243

   \alpha = \alpha_{soi}\left( 1 - f_{sno} \right) + \alpha_{sno}f_{sno}

where :math:`f_{sno}` is the fraction of ground covered by snow (equation (2.3) ), and :math:`\alpha_{sno} = 1.0`. The term :math:`\alpha_{soi}` is a function ranging from one when the soil column is wet to zero when the soil is dry

.. math::
   :label: eq-0244

   \alpha_{soi} = \sum_{i = 1}^{N_{levsoi}}{w_{i}r_{i}}

where :math:`w_{i}` is a soil wetness factor for layer :math:`i`, and :math:`r_{i}` is the relative contribution of each layer. The wetness factor :math:`w_{i}` is



.. math::
   :label: eq-0053

   w_{i} = \left\{ \begin{aligned}
   & \frac{\theta_{liq,i} - \theta_{dry,i}}{\theta_{opt,i} - \theta_{dry,i}} \text{for }T_{i} \geq T_{f} \\
   & 0 \text{for }T_{i} < T_{f}
   \end{aligned} \right\}



where :math:`\theta_{liq,i} - \theta_{dry,i} \geq 0` and

.. math::
   :label: eq-0245

   r_{i} = 0.1 \text{for }i = 1,\ldots,N_{levsoi}

The term :math:`\theta_{dry}` is the volumetric water content at which evapotranspiration ceases and :math:`\theta_{opt}` is the optimal water content

.. math::
   :label: eq-0246

   \theta_{dry,i} = \theta_{sat,\mspace{6mu} i}\left( \frac{- 316230}{\psi_{sat,i}} \right)^{- \frac{1}{B_{i}}}

.. math::
   :label: eq-0247

   \theta_{opt,i} = \theta_{sat,i}\left( \frac{- 158490}{\psi_{sat,i}} \right)^{- \frac{1}{B_{i}}}

where :math:`\theta_{sat,i}` is the water content at saturation (i.e., porosity), :math:`\psi_{sat,i}` is the saturated soil matric potential (mm), and :math:`B_{i}` is the Clapp-Hornberger exponent (section 5.3.1). The soil volumetric liquid water content :math:`\theta_{liq,i}` is

.. math::
   :label: eq-0248

   \theta_{liq,i} = \frac{w_{liq,i}}{\Delta z_{i}\rho_{liq}} \leq \theta_{sat,i} - \theta_{ice,i}

where :math:`w_{liq,i}` is the mass of liquid water (kg m\ :sup:`-2`), :math:`\Delta z_{i}` is the layer thickness, :math:`\rho_{liq}` is the density of liquid water (kg m\ :sup:`-3`) (:numref:`table-physical-constants`), and :math:`\theta_{ice,i}` is the volumetric ice content

.. math::
   :label: eq-0249

   \theta_{ice,i} = \frac{w_{ice,i}}{\Delta z_{i}\rho_{ice}} \leq \theta_{sat,i}

where :math:`w_{ice,i}` is the mass of ice (kg m\ :sup:`-2`) and :math:`\rho_{ice}` is the density of ice (kg m\ :sup:`-3`) (:numref:`table-physical-constants`). If :math:`q_{sat}^{T_{g}} > q_{atm}` and :math:`q_{atm} > q_{g,prvrd}`, then :math:`q_{g,prvrd} = q_{atm}` and :math:`\frac{dq_{g,prvrd}}{dT_{g}} = 0`.

The UCL specific humidity is then

.. math::
   :label: eq-0250

   q_{ac} = \frac{c_{a}^{w}q_{atm} + c_{roof}f_{wet,roof}q_{g,roof} + c_{prvrd}q_{g,prvrd} + c_{imprvrd}f_{wet,imprvrd}q_{g,imprvrd}}{c_{a}^{w} + f_{wet,roof}c_{roof} + c_{prvrd} + f_{wet,imprvrd}c_{imprvrd}}

where :math:`c_{a}^{w}` is the latent heat conductance from the UCL air to the atmosphere (:math:`\frac{1}{r_{aw}}`), and :math:`c_{roof}`, :math:`c_{prvrd}`, and :math:`c_{imprvrd}` are the weighted heat conductances from urban surfaces to UCL air [:math:`\frac{W_{roof}}{r_{s,roof}}`, :math:`\frac{W_{prvrd}}{r_{s,prvrd}}`, :math:`\frac{W_{imprvrd}}{r_{s,imprvrd}}`, respectively, where :math:`W_{prvrd} = \left( 1 - W_{roof} \right)f_{prvrd}`, :math:`W_{imprvrd} = \left( 1 - W_{roof} \right)\left( 1 - f_{prvrd} \right)`].

The stability is then updated using the new UCL air temperature and specific humidity as follows. The potential temperature, specific humidity, and virtual potential temperature scales, :math:`\theta_{*}`, :math:`q_{*}`, and :math:`\theta_{v*}`, are reevaluated using equations :eq:`eq-0213`-:eq:`eq-0220` and :eq:`eq-0197`. The wind speed including the convective velocity is reevaluated using equations :eq:`eq-0204` and :eq:`eq-0041`-(3.30). The Monin-Obukhov length is updated from equation :eq:`eq-0224`. This sequence of calculations is repeated for a total of three times beginning with the calculation of the friction velocity :math:`u_{*}` (equations :eq:`eq-0208`-:eq:`eq-0211`).

Final Fluxes and Adjustments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The sensible and latent heat fluxes and momentum flux from urban surfaces are then calculated from equations (3.69)-(3.73), (3.76)-(3.80), and :eq:`eq-0185`-:eq:`eq-0186` using the updated UCL air temperature and specific humidity. The water vapor flux from the pervious road, :math:`E_{prvrd}`, is assigned to ground evaporation, :math:`E_{g,prvrd}`, or a evapotranspiration term, :math:`E_{prvrd}^{et}`, as follows



.. math::
   :label: eq-0054

   E_{g,prvrd} = E_{prvrd} \text{for }q_{s} - q_{g,\mspace{6mu} prvrd} > 0\text{ or }f_{sno} > 0\text{ or }\alpha_{soi} = 0
   E_{prvrd}^{et} = E_{prvrd} \text{otherwise}



This ensures that dew can form on snow or soil surfaces and that snow can sublimate. Otherwise, the evaporation is assigned to an evapotranspiration term in which the water for evaporation is removed from all soil layers which have sufficient liquid water (section 5.3).

The partial derivatives of the urban surface fluxes with respect to surface temperatures, which are needed for the soil temperature calculation and to update the urban surface fluxes, are

.. math::
   :label: eq-0251

   \frac{\partial H_{roof}}{\partial T_{g,roof}} = \frac{\rho_{atm}C_{p}\left( c_{a}^{h} + \frac{c_{prvrd}}{W_{prvrd}} + \frac{c_{imprvrd}}{W_{imprvrd}} + \frac{c_{sunwall}}{W_{sunwall}} + \frac{c_{shdwall}}{W_{shdwall}} \right)\frac{c_{roof}}{W_{roof}}}{c_{a}^{h} + \frac{c_{roof}}{W_{roof}} + \frac{c_{prvrd}}{W_{prvrd}} + \frac{c_{imprvrd}}{W_{imprvrd}} + \frac{c_{sunwall}}{W_{sunwall}} + \frac{c_{shdwall}}{W_{shdwall}}}

.. math::
   :label: eq-0252

   \frac{\partial H_{prvrd}}{\partial T_{g,prvrd}} = \frac{\rho_{atm}C_{p}\left( c_{a}^{h} + \frac{c_{roof}}{W_{roof}} + \frac{c_{imprvrd}}{W_{imprvrd}} + \frac{c_{sunwall}}{W_{sunwall}} + \frac{c_{shdwall}}{W_{shdwall}} \right)\frac{c_{prvrd}}{W_{prvrd}}}{c_{a}^{h} + \frac{c_{roof}}{W_{roof}} + \frac{c_{prvrd}}{W_{prvrd}} + \frac{c_{imprvrd}}{W_{imprvrd}} + \frac{c_{sunwall}}{W_{sunwall}} + \frac{c_{shdwall}}{W_{shdwall}}}

.. math::
   :label: eq-0253

   \frac{\partial H_{imprvrd}}{\partial T_{g,imprvrd}} = \frac{\rho_{atm}C_{p}\left( c_{a}^{h} + \frac{c_{roof}}{W_{roof}} + \frac{c_{prvrd}}{W_{prvrd}} + \frac{c_{sunwall}}{W_{sunwall}} + \frac{c_{shdwall}}{W_{shdwall}} \right)\frac{c_{imprvrd}}{W_{imprvrd}}}{c_{a}^{h} + \frac{c_{roof}}{W_{roof}} + \frac{c_{prvrd}}{W_{prvrd}} + \frac{c_{imprvrd}}{W_{imprvrd}} + \frac{c_{sunwall}}{W_{sunwall}} + \frac{c_{shdwall}}{W_{shdwall}}}

.. math::
   :label: eq-0254

   \frac{\partial H_{sunwall}}{\partial T_{g,sunwall}} = \frac{\rho_{atm}C_{p}\left( c_{a}^{h} + \frac{c_{roof}}{W_{roof}} + \frac{c_{prvrd}}{W_{prvrd}} + \frac{c_{imprvrd}}{W_{imprvrd}} + \frac{c_{shdwall}}{W_{shdwall}} \right)\frac{c_{sunwall}}{W_{sunwall}}}{c_{a}^{h} + \frac{c_{roof}}{W_{roof}} + \frac{c_{prvrd}}{W_{prvrd}} + \frac{c_{imprvrd}}{W_{imprvrd}} + \frac{c_{sunwall}}{W_{sunwall}} + \frac{c_{shdwall}}{W_{shdwall}}}

:math:`\frac{\partial H_{shdwall}}{\partial T_{g,shdwall}} = \frac{\rho_{atm}C_{p}\left( c_{a}^{h} + \frac{c_{roof}}{W_{roof}} + \frac{c_{prvrd}}{W_{prvrd}} + \frac{c_{imprvrd}}{W_{imprvrd}} + \frac{c_{sunwall}}{W_{sunwall}} \right)\frac{c_{shdwall}}{W_{shdwall}}}{c_{a}^{h} + \frac{c_{roof}}{W_{roof}} + \frac{c_{prvrd}}{W_{prvrd}} + \frac{c_{imprvrd}}{W_{imprvrd}} + \frac{c_{sunwall}}{W_{sunwall}} + \frac{c_{shdwall}}{W_{shdwall}}}`,

:math:`\frac{\partial E_{roof}}{\partial T_{g,roof}} = \frac{\rho_{atm}\left( c_{a}^{w} + \frac{c_{prvrd}}{W_{prvrd}} + \frac{f_{wet,imprvrd}c_{imprvrd}}{W_{imprvrd}} \right)\frac{f_{wet,roof}c_{roof}}{W_{roof}}}{c_{a}^{w} + \frac{f_{wet,roof}c_{roof}}{W_{roof}} + \frac{c_{prvrd}}{W_{prvrd}} + \frac{f_{wet,imprvrd}c_{imprvrd}}{W_{imprvrd}}}\frac{dq_{g,roof}}{dT_{g,roof}}`,

:math:`\frac{\partial E_{prvrd}}{\partial T_{g,prvrd}} = \frac{\rho_{atm}\left( c_{a}^{w} + \frac{f_{wet,roof}c_{roof}}{W_{roof}} + \frac{f_{wet,imprvrd}c_{imprvrd}}{W_{imprvrd}} \right)\frac{c_{prvrd}}{W_{prvrd}}}{c_{a}^{w} + \frac{f_{wet,roof}c_{roof}}{W_{roof}} + \frac{c_{prvrd}}{W_{prvrd}} + \frac{f_{wet,imprvrd}c_{imprvrd}}{W_{imprvrd}}}\frac{dq_{g,prvrd}}{dT_{g,prvrd}}`,

:math:`\frac{\partial E_{imprvrd}}{\partial T_{g,imprvrd}} = \frac{\rho_{atm}\left( c_{a}^{w} + \frac{f_{wet,\mspace{6mu} roof}c_{roof}}{W_{roof}} + \frac{c_{prvrd}}{W_{prvrd}} \right)\frac{f_{wet,imprvrd}c_{imprvrd}}{W_{imprvrd}}}{c_{a}^{w} + \frac{f_{wet,roof}c_{roof}}{W_{roof}} + \frac{c_{prvrd}}{W_{prvrd}} + \frac{f_{wet,imprvrd}c_{imprvrd}}{W_{imprvrd}}}\frac{dq_{g,imprvrd}}{dT_{g,imprvrd}}`,

:math:`\frac{\partial E_{sunwall}}{\partial T_{g,sunwall}} = 0`,

.. math::
   :label: eq-0255

   \frac{\partial E_{shdwall}}{\partial T_{g,shdwall}} = 0

The 2-m air temperature diagnostic is set equal to the UCL air temperature :math:`T_{ac}`\ and the 2-m specific humidity diagnostic is set equal to the UCL specific humidity :math:`q_{ac}`. Relative humidity of the UCL air is

.. math::
   :label: eq-0256

   RH_{ac} = \min\left( 100,\frac{q_{ac}}{q_{sat}^{T_{ac}}} \times 100 \right)

where :math:`q_{sat}^{T_{ac}}` is the saturated specific humidity at UCL air temperature :math:`T_{ac}` (section 3.3).

The sensible heat and water vapor fluxes are based on the urban surface temperature from the previous time step, :math:`T_{g}^{n}`, and are used as the surface forcing for the solution of the soil temperature equations (section 4). This solution yields a new surface temperature :math:`T_{g}^{n + 1}`. The sensible heat and water vapor fluxes are updated to :math:`H_{g}^{'}` and :math:`E_{g}^{'}`\ for the new temperature as

:math:`H_{g}^{'} = H_{g} + \left( T_{g}^{n + 1} - T_{g}^{n} \right)\frac{\partial H_{g}}{\partial T_{g}}`,

.. math::
   :label: eq-0257

   E_{g}^{'} = E_{g} + \left( T_{g}^{n + 1} - T_{g}^{n} \right)\frac{\partial E_{g}}{\partial T_{g}}

where :math:`H_{g}` and :math:`E_{g}` are the sensible heat and water vapor fluxes derived above, and :math:`g` denotes each of the five urban surfaces. One further adjustment is made to the fluxes for the roof, pervious and impervious road. If the surface moisture (i.e., the ponded water in the case of the roof and impervious road, and top layer moisture for the pervious road) is not sufficient to supported the updated evaporation, i.e., if :math:`E_{g}^{'} > 0` and :math:`f_{evap} < 1` where

:math:`f_{evap} = \frac{\frac{\left( w_{ice,\mspace{6mu} snl + 1} + w_{liq,snl + 1} \right)}{\Delta t}}{E_{g}^{'}} \leq 1`,

an adjustment is made to reduce the ground evaporation accordingly as

.. math::
   :label: eq-0258

   E_{g}^{''} = f_{evap}E_{g}^{'}

:math:`w_{ice,snl + 1}` and :math:`w_{liq,snl + 1}` are the ice and liquid water contents (kg m\ :sup:`-2`) of the top layer. Any resulting energy deficit is assigned to sensible heat as

.. math::
   :label: eq-0259

   H_{g}^{''} = H_{g} + \lambda\left( E_{g}^{'} - E_{g}^{''} \right)

The water vapor flux :math:`E_{g}^{''}` is partitioned into evaporation of liquid water :math:`q_{seva}`, sublimation from ice :math:`q_{subl}`, liquid dew :math:`q_{sdew}`, or frost :math:`q_{frost}` (all in kg m\ :sup:`-2` s\ :sup:`-1`) as

.. math::
   :label: eq-0260

   q_{seva} = \max\left( E_{g}^{''}\frac{w_{liq,snl + 1}}{w_{ice,\mspace{6mu} snl + 1} + w_{liq,snl + 1}},0 \right) E_{g}^{''} \geq 0,w_{ice,snl + 1} + w_{liq,snl + 1} > 0

.. math::
   :label: eq-0261

   q_{subl} = E_{g}^{''} - q_{seva} E_{g}^{''} \geq 0

.. math::
   :label: eq-0262

   q_{sdew} = \left| E_{g}^{''} \right| E_{g}^{''} < 0\text{ and }T_{g} \geq T_{f}

.. math::
   :label: eq-0263

   q_{frost} = \left| E_{g}^{''} \right| E_{g}^{''} < 0\text{ and }T_{g} < T_{f}

The loss or gain in snow mass due to :math:`q_{seva}`, :math:`q_{subl}`, :math:`q_{sdew}`, and :math:`q_{frost}` on a snow surface are accounted for during the snow hydrology calculations (section 5.1). The loss of surface water from non-snow surfaces due to :math:`q_{seva}` is accounted for in the calculation of infiltration (section 5.2), while losses or gains due to :math:`q_{subl}`, :math:`q_{sdew}`, and :math:`q_{frost}` from non-snow surfaces are accounted for following sub-surface drainage calculations (section 5.4).

The ground or storage heat flux :math:`G` for each urban surface is calculated as

.. math::
   :label: eq-0264

   G = {\overrightarrow{S}}_{g} - {\overrightarrow{L}}_{g} - H_{g} - \lambda E_{g} + H_{wasteheat,g} + H_{aircond,g}

where :math:`{\overrightarrow{S}}_{g}` is the absorbed solar radiation (section 2.5), :math:`H_{g}` and :math:`\lambda E_{g}` are the sensible and latent heat fluxes after the adjustments described above, and :math:`H_{wasteheat,g}` and :math:`H_{aircond,g}` are the waste heat and heat removed by air conditioning (section 4.1). The net longwave radiation :math:`{\overrightarrow{L}}_{g}` is updated for the change in surface temperature as

.. math::
   :label: eq-0265

   {\overrightarrow{L}}_{g}^{n + 1} = {\overrightarrow{L}}_{g}^{n} + 4\varepsilon_{g}\sigma\left( T_{g}^{n} \right)^{3}\left( T_{g}^{n + 1} - T_{g}^{n} \right)

When converting water vapor flux to an energy flux, the term :math:`\lambda` is arbitrarily assumed to be



.. math::
   :label: eq-0055

   \lambda = \left\{ \begin{aligned}
   & \lambda_{sub} \text{if }w_{liq,snl + 1} = 0\text{ and }w_{ice,snl + 1} > 0 \\
   & \lambda_{vap} \text{otherwise}
   \end{aligned} \right\}



where :math:`\lambda_{sub}` and :math:`\lambda_{vap}` are the latent heat of sublimation and vaporization, respectively (J kg\ :sup:`-1`) (:numref:`table-physical-constants`).

Saturation Specific Humidity
----------------------------

Saturation vapor pressure :math:`e_{sat}^{T}` (Pa) and its derivative :math:`\frac{de_{sat}^{T}}{dT}`, as a function of temperature :math:`T` (°C), are calculated from the eighth-order polynomial fits of Flatau et al. (1992)

.. math::
   :label: eq-0266

   e_{sat}^{T} = 100\left\lbrack a_{0} + a_{1}T + \cdots + a_{n}T^{n} \right\rbrack

.. math::
   :label: eq-0267

   \frac{de_{sat}^{T}}{dT} = 100\left\lbrack b_{0} + b_{1}T + \cdots + b_{n}T^{n} \right\rbrack

where the coefficients for ice are valid for :math:`- 75^{\circ}C \leq T < 0^{\circ}C` and the coefficients for water are valid for :math:`0^{\circ}C \leq T \leq 100^{\circ}C` (:numref:`table-esat-coefficients` and :numref:`table-desat-coefficients`). The saturated water vapor specific humidity :math:`q_{sat}^{T}` and its derivative :math:`\frac{dq_{sat}^{T}}{dT}` are

.. math::
   :label: eq-0268

   q_{sat}^{T} = \frac{0.622e_{sat}^{T}}{P_{atm} - 0.378e_{sat}^{T}}

.. math::
   :label: eq-0269

   \frac{dq_{sat}^{T}}{dT} = \frac{0.622P_{atm}}{\left( P_{atm} - 0.378e_{sat}^{T} \right)^{2}}\frac{de_{sat}^{T}}{dT}

.. list-table:: Coefficients for :math:`e_{sat}^{T}`
   :widths: auto
   :header-rows: 0
   :name: table-esat-coefficients


   * - 
     - water
     - ice
   * - :math:`a_{0}`
     - 6.11213476
     - 6.11123516
   * - :math:`a_{1}`
     - 4.44007856\ :math:`\times 10^{- 1}`
     - 5.03109514\ :math:`\times 10^{- 1}`
   * - :math:`a_{2}`
     - 1.43064234\ :math:`\times 10^{- 2}`
     - 1.88369801\ :math:`\times 10^{- 2}`
   * - :math:`a_{3}`
     - 2.64461437\ :math:`\times 10^{- 4}`
     - 4.20547422\ :math:`\times 10^{- 4}`
   * - :math:`a_{4}`
     - 3.05903558\ :math:`\times 10^{- 6}`
     - 6.14396778\ :math:`\times 10^{- 6}`
   * - :math:`a_{5}`
     - 1.96237241\ :math:`\times 1{0^{-}}^{8}`
     - 6.02780717\ :math:`\times 10^{- 8}`
   * - :math:`a_{6}`
     - 8.92344772\ :math:`\times 10^{- 11}`
     - 3.87940929\ :math:`\times 10^{- 10}`
   * - :math:`a_{7}`
     - -3.73208410\ :math:`\times 10^{- 13}`
     - 1.49436277\ :math:`\times 10^{- 12}`
   * - :math:`a_{8}`
     - 2.09339997\ :math:`\times 10^{- 16}`
     - 2.62655803\ :math:`\times 10^{- 15}`


.. list-table:: Coefficients for :math:`\frac{de_{sat}^{T}}{dT}`
   :widths: auto
   :header-rows: 0
   :name: table-desat-coefficients


   * - 
     - water
     - ice
   * - :math:`b_{0}`
     - 4.44017302\ :math:`\times 10^{- 1}`
     - 5.03277922\ :math:`\times 10^{- 1}`
   * - :math:`b_{1}`
     - 2.86064092\ :math:`\times 10^{- 2}`
     - 3.77289173\ :math:`\times 10^{- 2}`
   * - :math:`b_{2}`
     - 7.94683137\ :math:`\times 10^{- 4}`
     - 1.26801703\ :math:`\times 10^{- 3}`
   * - :math:`b_{3}`
     - 1.21211669\ :math:`\times 10^{- 5}`
     - 2.49468427\ :math:`\times 10^{- 5}`
   * - :math:`b_{4}`
     - 1.03354611\ :math:`\times 10^{- 7}`
     - 3.13703411\ :math:`\times 10^{- 7}`
   * - :math:`b_{5}`
     - 4.04125005\ :math:`\times 10^{- 10}`
     - 2.57180651\ :math:`\times 10^{- 9}`
   * - :math:`b_{6}`
     - -7.88037859\ :math:`\times 10^{- 13}`
     - 1.33268878\ :math:`\times 10^{- 11}`
   * - :math:`b_{7}`
     - -1.14596802\ :math:`\times 10^{- 14}`
     - 3.94116744\ :math:`\times 10^{- 14}`
   * - :math:`b_{8}`
     - 3.81294516\ :math:`\times 10^{- 17}`
     - 4.98070196\ :math:`\times 10^{- 17}`


Roof, Wall, Road, and Snow Temperatures
=======================================

The first law of heat conduction is

.. math::
   :label: eq-0270

   F = - \lambda\nabla T

where :math:`F` is the amount of heat conducted across a unit cross-sectional area in unit time (W m\ :sup:`-2`), :math:`\lambda` is thermal conductivity (W m\ :sup:`-1` K\ :sup:`-1`), and :math:`\nabla T` is the spatial gradient of temperature (K m\ :sup:`-1`). In one-dimensional form

.. math::
   :label: eq-0271

   F_{z} = - \lambda\frac{\partial T}{\partial z}

where :math:`z` is in the vertical direction (m) and is positive downward and :math:`F_{z}` is positive upward. To account for non-steady or transient conditions, the principle of energy conservation in the form of the continuity equation is invoked as

.. math::
   :label: eq-0272

   c\frac{\partial T}{\partial t} = - \frac{\partial F_{z}}{\partial z}

where :math:`c` is the volumetric snow/soil heat capacity (J m\ :sup:`-3` K\ :sup:`-1`) and :math:`t` is time (s). Combining equations :eq:`eq-0271` and :eq:`eq-0272` yields the second law of heat conduction in one-dimensional form

.. math::
   :label: eq-0273

   c\frac{\partial T}{\partial t} = \frac{\partial}{\partial z}\left\lbrack \lambda\frac{\partial T}{\partial z} \right\rbrack

The nature of the solution of this equations :eq:`eq-0271` and :eq:`eq-0272` depends on the type of urban surface. The solution for pervious and impervious roads follows the solution for CLM soils where the equation is solved numerically for a fifteen-layer column with up to five overlying layers of snow with the boundary conditions of :math:`h` as the heat flux into the surface layer from the overlying atmosphere and zero heat flux at the bottom of the soil column. In the case of pervious roads, the temperature profile is calculated first without phase change and then readjusted for phase change (section 4.2). For impervious roads, however, the moisture content of all layers is zero. Phase change then only takes place in the ponded surface water. The roof consists of a fifteen-layer column with potential ponded surface water including up to a five layer snow pack, however, the bottom boundary condition is a non-zero flux governed by prescribed controls on the internal building temperature. The walls are modeled similarly to roofs except for the absence of ponded water or snow.

Numerical Solution
------------------

Roofs and walls are discretized into fifteen layers where the depth of layer :math:`i`, or node depth, :math:`z_{i}` (m), is

.. math::
   :label: eq-0274

   z_{i} = (i - 0.5)\left( \frac{\Delta z}{N_{levgrnd}} \right)

where :math:`\Delta z` is the total thickness of the roof or wall (:numref:`table-input-data`) and :math:`N_{levgrnd} = 15` is the number of layers. The thickness of each layer :math:`\Delta z_{i}` (m) is



.. math::
   :label: eq-0056

   \Delta z_{i} = \left\{ \begin{aligned}
   & 0.5\left( z_{1} + z_{2} \right) i = 1 \\
   & 0.5\left( z_{i + 1} - z_{i - 1} \right) i = 2,\ldots,N_{levgrnd} - 1 \\
   & z_{N_{levgrnd}} - z_{N_{levgrnd} - 1} i = N_{levgrnd}
   \end{aligned} \right\}


The depths at the layer interfaces :math:`z_{h,i}` (m) are



.. math::
   :label: eq-0057

   z_{h,i} = \left\{ \begin{aligned}
   & 0 i = 0 \\
   & 0.5\left( z_{i} + z_{i + 1} \right) i = 1,\ldots,N_{levgrnd} - 1 \\
   & z_{N_{levgrnd}} + 0.5\Delta z_{N_{levgrnd}} i = N_{levgrnd}
   \end{aligned} \right\}


Pervious and impervious road are discretized into fifteen layers as well with node depth

.. math::
   :label: eq-0275

   z_{i} = f_{s}\left\{ \exp\left\lbrack 0.5(i - 0.5) \right\rbrack - 1 \right\}

where :math:`f_{s} = 0.025` is a scaling factor. Layer thicknesses and interface depths are calculated from equations :eq:`eq-0056` and :eq:`eq-0057`.

The overlying snow pack for the roof and road is modeled with up to five layers depending on the total snow depth. The layers from top to bottom are indexed in the Fortran code as :math:`i = - 4, - 3, - 2, - 1,0`, which permits the accumulation or ablation of snow at the top of the snow pack without renumbering the layers. Layer :math:`i = 0` is the snow layer next to the urban surface and layer :math:`i = snl + 1` is the top layer, where the variable :math:`snl` is the negative of the number of snow layers. The number of snow layers and the thickness of each layer is a function of snow depth :math:`z_{sno}` (m) as follows.



.. math::
   :label: eq-0058

   \left\{ \begin{aligned}
   & snl = - 1 \\
   & \Delta z_{0} = z_{sno} \text{for 0.01} \leq z_{\text{sno}} \leq 0.03
   \end{aligned} \right\}



.. math::
   :label: eq-0059

   \left\{ \begin{aligned}
   & snl = - 2 \\
   & \Delta z_{- 1} = \frac{z_{sno}}{2} \text{for 0.03} < z_{\text{sno}} \leq 0.04 \\
   & \Delta z_{0} = \Delta z_{- 1}
   \end{aligned} \right\}



.. math::
   :label: eq-0060

   \left\{ \begin{aligned}
   & snl = - 2 \\
   & \Delta z_{- 1} = 0.02 \text{for 0.04} < z_{\text{sno}} \leq 0.07 \\
   & \Delta z_{0} = z_{sno} - \Delta z_{- 1}
   \end{aligned} \right\}



.. math::
   :label: eq-0061

   \left\{ \begin{aligned}
   & snl = - 3 \\
   & \Delta z_{- 2} = 0.02 \\
   & \Delta z_{- 1} = \frac{\left( z_{sno} - 0.02 \right)}{2} \text{for 0.07} < z_{\text{sno}} \leq 0.12 \\
   & \Delta z_{0} = \Delta z_{- 1}
   \end{aligned} \right\}



.. math::
   :label: eq-0062

   \left\{ \begin{aligned}
   & snl = - 3 \\
   & \Delta z_{- 2} = 0.02 \\
   & \Delta z_{- 1} = 0.05 \text{for 0.12} < z_{\text{sno}} \leq 0.18 \\
   & \Delta z_{0} = z_{sno} - \Delta z_{- 2} - \Delta z_{- 1}
   \end{aligned} \right\}



.. math::
   :label: eq-0063

   \left\{ \begin{aligned}
   & snl = - 4 \\
   & \Delta z_{- 3} = 0.02 \\
   & \Delta z_{- 2} = 0.05 \text{for 0.18} < z_{\text{sno}} \leq 0.29 \\
   & \Delta z_{- 1} = \frac{\left( z_{sno} - \Delta z_{- 3} - \Delta z_{- 2} \right)}{2} \\
   & \Delta z_{0} = \Delta z_{- 1}
   \end{aligned} \right\}



.. math::
   :label: eq-0064

   \left\{ \begin{aligned}
   & snl = - 4 \\
   & \Delta z_{- 3} = 0.02 \\
   & \Delta z_{- 2} = 0.05 \text{for 0.29} < z_{\text{sno}} \leq 0.41 \\
   & \Delta z_{- 1} = 0.11 \\
   & \Delta z_{0} = z_{sno} - \Delta z_{- 3} - \Delta z_{- 2} - \Delta z_{- 1}
   \end{aligned} \right\}



.. math::
   :label: eq-0065

   \left\{ \begin{aligned}
   & snl = - 5 \\
   & \Delta z_{- 4} = 0.02 \\
   & \Delta z_{- 3} = 0.05 \text{for 0.41} < z_{\text{sno}} \leq 0.64 \\
   & \Delta z_{- 2} = 0.11 \\
   & \Delta z_{- 1} = \frac{\left( z_{sno} - \Delta z_{- 4} - \Delta z_{- 3} - \Delta z_{- 2} \right)}{2} \\
   & \Delta z_{0} = \Delta z_{- 1}
   \end{aligned} \right\}



.. math::
   :label: eq-0066

   \left\{ \begin{aligned}
   & snl = - 5 \\
   & \Delta z_{- 4} = 0.02 \\
   & \Delta z_{- 3} = 0.05 \text{for 0.64} < z_{\text{sno}} \\
   & \Delta z_{- 2} = 0.11 \\
   & \Delta z_{- 1} = 0.23 \\
   & \Delta z_{0} = z_{sno} - \Delta z_{- 4} - \Delta z_{- 3} - \Delta z_{- 2} - \Delta z_{- 1}
   \end{aligned} \right\}


The node depths, which are located at the midpoint of the snow layers, and the layer interfaces are both referenced from the urban surface and are defined as negative values

.. math::
   :label: eq-0276

   z_{i} = z_{h,i} - 0.5\Delta z_{i} i = snl + 1,\ldots,0

.. math::
   :label: eq-0277

   z_{h,i} = z_{h,i + 1} - \Delta z_{i + 1} i = snl,\ldots, - 1

Note that :math:`z_{h,0}`, the interface between the bottom snow layer and the top urban layer, is zero. Thermal properties (i.e., temperature :math:`T_{i}` [K]; thermal conductivity :math:`\lambda_{i}` [W m\ :sup:`-1` K\ :sup:`-1`]; volumetric heat capacity :math:`c_{i}` [J m\ :sup:`-3` K\ :sup:`-1`]) are defined for layers at the node depths (:numref:`fig-numerical-scheme-layer-temps`) and for snow layers at the layer midpoints.

In general, for a zero-flux bottom boundary condition, the heat flux :math:`F_{i}` (W m\ :sup:`-2`) from layer :math:`i` to layer :math:`i + 1` is

.. math::
   :label: eq-0278

   F_{i} = - \lambda\left\lbrack z_{h,i} \right\rbrack\left( \frac{T_{i} - T_{i + 1}}{z_{i + 1} - z_{i}} \right)

where the thermal conductivity at the interface :math:`\lambda\left\lbrack z_{h,i} \right\rbrack` is



.. math::
   :label: eq-0067

   \lambda\left\lbrack z_{h,i} \right\rbrack = \left\{ \begin{aligned}
   & \frac{\lambda_{i}\lambda_{i + 1}\left( z_{i + 1} - z_{i} \right)}{\lambda_{i}\left( z_{i + 1} - z_{h,i} \right) + \lambda_{i + 1}\left( z_{h,i} - z_{i} \right)} i = snl + 1,\ldots,N_{levgrnd} - 1 \\
   & 0 i = N_{levgrnd}
   \end{aligned} \right\}


For a non-zero flux bottom boundary condition, :math:`\lambda\left\lbrack z_{h,i = N_{levgrnd}} \right\rbrack = \lambda_{i = N_{levgrnd}}`. These equations are derived, with reference to :numref:`fig-numerical-scheme-layer-temps`, assuming that the heat flux from :math:`i` (depth :math:`z_{i}`) to the interface between :math:`i` and :math:`i + 1` (depth :math:`z_{h,i}`) equals the heat flux from the interface to :math:`i + 1` (depth :math:`z_{i + 1}`), i.e.,

.. math::
   :label: eq-0279

   - \lambda_{i}\frac{T_{i} - T_{m}}{z_{h,i} - z_{i}} = - \lambda_{i + 1}\frac{T_{m} - T_{i + 1}}{z_{i + 1} - z_{h,i}}

where :math:`T_{m}` is the temperature at the interface of layers :math:`i` and :math:`i + 1`. Solving equation :eq:`eq-0279` for :math:`T_{m}` and substituting :math:`T_{m}` back into the left side of equation :eq:`eq-0279` yields equations :eq:`eq-0278` and :eq:`eq-0067`.

.. figure:: image12.jpeg
   :width: 6in
   :height: 4.47917in
   :name: fig-numerical-scheme-layer-temps

   Schematic diagram of numerical scheme used to solve for layer temperatures. Shown are three layers, :math:`i - 1`, :math:`i`, and :math:`i + 1`. The thermal conductivity :math:`\lambda`, specific heat capacity :math:`c`, and temperature :math:`T` are defined at the layer node depth :math:`z`. :math:`T_{m}` is the interface temperature. The thermal conductivity :math:`\lambda\left\lbrack z_{h} \right\rbrack` is defined at the interface of two layers :math:`z_{h}`. The layer thickness is :math:`\Delta z`. The heat fluxes :math:`F_{i - 1}` and :math:`F_{i}` are defined as positive upwards.

The energy balance for the :math:`i^{th}` layer is

.. math::
   :label: eq-0280

   \frac{c_{i}\Delta z_{i}}{\Delta t}\left( T_{i}^{n + 1} - T_{i}^{n} \right) = - F_{i - 1} + F_{i}

where the superscripts :math:`n` and :math:`n + 1` indicate values at the beginning and end of the time step, respectively, and :math:`\Delta t` is the time step (s). This equation is solved using the Crank-Nicholson method, which combines the explicit method with fluxes evaluated at :math:`n` (:math:`F_{i - 1}^{n},F_{i}^{n}`) and the implicit method with fluxes evaluated at :math:`n + 1` (:math:`F_{i - 1}^{n + 1},F_{i}^{n + 1}`)

.. math::
   :label: eq-0281

   \frac{c_{i}\Delta z_{i}}{\Delta t}\left( T_{i}^{n + 1} - T_{i}^{n} \right) = \alpha\left( - F_{i - 1}^{n} + F_{i}^{n} \right) + (1 - \alpha)\left( - F_{i - 1}^{n + 1} + F_{i}^{n + 1} \right)

where :math:`\alpha = 0.5`, resulting in a tridiagonal system of equations

.. math::
   :label: eq-0282

   r_{i} = a_{i}T_{i - 1}^{n + 1} + b_{i}T_{i}^{n + 1} + c_{i}T_{i + 1}^{n + 1}

where :math:`a_{i}`, :math:`b_{i}`, and :math:`c_{i}` are the subdiagonal, diagonal, and superdiagonal elements in the tridiagonal matrix and :math:`r_{i}` is a column vector of constants.

For the top layer :math:`i = snl + 1`, the heat flux from the overlying atmosphere into the surface layer :math:`h` (W m\ :sup:`-2`, defined as positive into the surface) is

.. math::
   :label: eq-0283

   h^{n + 1} = - \alpha F_{i - 1}^{n} - (1 - \alpha)F_{i - 1}^{n + 1}

The energy balance for layer :math:`i = snl + 1` is then

.. math::
   :label: eq-0284

   \frac{c_{i}\Delta z_{i}}{\Delta t}\left( T_{i}^{n + 1} - T_{i}^{n} \right) = h^{n + 1} + \alpha F_{i}^{n} + (1 - \alpha)F_{i}^{n + 1}

The heat flux :math:`h` at :math:`n + 1` may be approximated as follows

.. math::
   :label: eq-0285

   h^{n + 1} = h^{n} + \frac{\partial h}{\partial T_{i}}\left( T_{i}^{n + 1} - T_{i}^{n} \right)

The resulting equations are

.. figure:: image13.png

.. math::
   :label: eq-0286

   a_{i} = 0

.. math::
   :label: eq-0287

   b_{i} = 1 + \frac{\Delta t}{c_{i}\Delta z_{i}}\left\lbrack (1 - \alpha)\frac{\lambda\left\lbrack z_{h,i} \right\rbrack}{z_{i + 1} - z_{i}} - \frac{\partial h}{\partial T_{i}} \right\rbrack

.. math::
   :label: eq-0288

   c_{i} = - (1 - \alpha)\frac{\Delta t}{c_{i}\Delta z_{i}}\frac{\lambda\left\lbrack z_{h,i} \right\rbrack}{z_{i + 1} - z_{i}}

.. math::
   :label: eq-0289

   r_{i} = T_{i}^{n} + \frac{\Delta t}{c_{i}\Delta z_{i}}\left\lbrack h^{n} - \frac{\partial h}{\partial T_{i}}T_{i}^{n} + \alpha F_{i} \right\rbrack

where

.. math::
   :label: eq-0290

   F_{i} = - \lambda\left\lbrack z_{h,i} \right\rbrack\left( \frac{T_{i}^{n} - T_{i + 1}^{n}}{z_{i + 1} - z_{i}} \right)

The heat flux into each urban surface :math:`h` is

.. math::
   :label: eq-0291

   h = {\overrightarrow{S}}_{g} - {\overrightarrow{L}}_{g} - H_{g} - \lambda E_{g} + H_{wasteheat,g} + H_{aircond,g}

where :math:`{\overrightarrow{S}}_{g}` is the absorbed solar radiation (section 2.5), :math:`{\overrightarrow{L}}_{g}` is the net longwave radiation (section 2.7), and :math:`H_{g}` and :math:`\lambda E_{g}` are the sensible and latent heat fluxes (section 3.2). The terms :math:`H_{wasteheat,g}` and :math:`H_{aircond,g}` are the waste heat from space heating/air conditioning and heat removed by air conditioning applied only to the pervious and impervious road



.. math::
   :label: eq-0068

   H_{wasteheat,prvrd} = H_{wasteheat,imprvrd} = \frac{H_{wasteheat}}{1 - W_{roof}}
   

\ 

.. math::
   :label: eq-0069

   \begin{cases}
   H_{wasteheat,sunwall} = H_{wasteheat,shdwall} = H_{wasteheat,roof} = 0 \\
   H_{aircond,prvrd} = H_{aircond,imprvrd} = \frac{H_{aircond}}{1 - W_{roof}}
   \end{cases}


\

.. math::
   :label: eq-0070

   H_{aircond,sunwall} = H_{aircond,shdwall} = H_{aircond,roof} = 0

where :math:`H_{wasteheat}` and :math:`H_{aircond}` are the total waste heat and heat removed by air conditioning from equations :eq:`eq-0074` and :eq:`eq-0314`. Note that for the pervious road, the latent heat is always the total latent heat regardless of its partitioning into ground evaporation or transpiration (section 3.2.4). The partial derivative of the heat flux :math:`h` with respect to surface temperature is

.. math::
   :label: eq-0292

   \frac{\partial h}{\partial T_{g}} = - \frac{\partial{\overrightarrow{L}}_{g}}{\partial T_{g}} - \frac{\partial H_{g}}{\partial T_{g}} - \frac{\partial\lambda E_{g}}{\partial T_{g}}

where the partial derivative of the net longwave radiation is

.. math::
   :label: eq-0293

   \frac{\partial{\overrightarrow{L}}_{g}}{\partial T_{g}} = 4\varepsilon_{g}\sigma\left( T_{g}^{n} \right)^{3}

and the partial derivatives of the sensible and latent heat fluxes are given by equations :eq:`eq-0251`-:eq:`eq-0255`. :math:`\sigma` is the Stefan-Boltzmann constant (W m\ :sup:`-2` K\ :sup:`-4`) (:numref:`table-physical-constants`) and :math:`\varepsilon_{g}` is the surface emissivity (section 2.7).

The top layer for roofs and walls is thin enough such that the layer-averaged temperature calculated above is considered to be equivalent to the surface temperature :math:`T_{g}^{n + 1}`. For pervious and impervious road, the top layer temperature has somewhat reduced diurnal amplitude compared with surface temperature. An accurate surface temperature is provided that compensates for this effect and numerical error by tuning the heat capacity of the top layer (through adjustment of the layer thickness) to give an exact match to the analytic solution for diurnal heating. The layer thickness for :math:`i = snl + 1` is given by

.. math::
   :label: eq-0294

   \Delta z_{i*} = 0.5\left\lbrack z_{i} - z_{h,i - 1} + c_{a}\left( z_{i + 1} - z_{h,i - 1} \right) \right\rbrack

where :math:`c_{a}` is a tunable parameter, varying from 0 to 1, and is taken as 0.34 by comparing the numerical solution with the analytic solution (Z.-L. Yang 1998, unpublished manuscript). For pervious and impervious road, :math:`\Delta z_{i*}` is used in place of :math:`\Delta z_{i}` for :math:`i = snl + 1` in equations (4.20)-:eq:`eq-0289`.

For the pervious and impervious road, the boundary condition at the bottom is zero heat flux, :math:`F_{i} = 0`, resulting in, for :math:`i = N_{levgrnd}`,

.. math::
   :label: eq-0295

   \frac{c_{i}\Delta z_{i}}{\Delta t}\left( T_{i}^{n + 1} - T_{i}^{n} \right) = \alpha\frac{\lambda\left\lbrack z_{h,i - 1} \right\rbrack\left( T_{i - 1}^{n} - T_{i}^{n} \right)}{z_{i} - z_{i - 1}} + (1 - \alpha)\frac{\lambda\left\lbrack z_{h,i - 1} \right\rbrack\left( T_{i - 1}^{n + 1} - T_{i}^{n + 1} \right)}{z_{i} - z_{i - 1}}

.. math::
   :label: eq-0296

   a_{i} = - (1 - \alpha)\frac{\Delta t}{c_{i}\Delta z_{i}}\frac{\lambda\left\lbrack z_{h,i - 1} \right\rbrack}{z_{i} - z_{i - 1}}

.. math::
   :label: eq-0297

   b_{i} = 1 + (1 - \alpha)\frac{\Delta t}{c_{i}\Delta z_{i}}\frac{\lambda\left\lbrack z_{h,i - 1} \right\rbrack}{z_{i} - z_{i - 1}}

.. math::
   :label: eq-0298

   c_{i} = 0

.. math::
   :label: eq-0299

   r_{i} = T_{i}^{n} - \alpha\frac{\Delta t}{c_{i}\Delta z_{i}}F_{i - 1}

where

.. math::
   :label: eq-0300

   F_{i - 1} = - \frac{\lambda\left\lbrack z_{h,i - 1} \right\rbrack}{z_{i} - z_{i - 1}}\left( T_{i - 1}^{n} - T_{i}^{n} \right)

For the roof and walls, the boundary condition at the bottom is the internal building temperature :math:`T_{iB}`, constrained as :math:`{TiB_{iB,\min}}_{iB,\max}`, where :math:`T_{iB,\max}` and :math:`T_{iB,\min}` are prescribed maximum and minimum internal building temperatures (:numref:`table-input-data`). The internal building temperature :math:`T_{iB}` is determined from a weighted combination of the inner layer wall and roof temperatures as

.. math::
   :label: eq-0301

   T_{iB} = \frac{H\left( T_{i = N_{levgrnd},shdwall}^{n} + T_{i = N_{levgrnd},sunwall}^{n} \right) + L_{roof}T_{i = N_{levgrnd},roof}^{n}}{2H + L_{roof}}

where :math:`H` is the building height and :math:`L_{roof}` is the length of the roof in an infinite canyon configuration

.. math::
   :label: eq-0302

   L = \left( \frac{H}{\frac{H}{W}} \right)\left( \frac{W_{roof}}{1 - W_{roof}} \right)

This boundary condition yields, for :math:`i = N_{levgrnd}`,

.. figure:: image14.png

.. math::
   :label: eq-0303

   a_{i} = - (1 - \alpha)\frac{\Delta t}{c_{i}\Delta z_{i}}\frac{\lambda\left\lbrack z_{h,i - 1} \right\rbrack}{z_{i} - z_{i - 1}}

.. math::
   :label: eq-0304

   b_{i} = 1 + (1 - \alpha)\frac{\Delta t}{c_{i}\Delta z_{i}}\left\lbrack \frac{\lambda\left\lbrack z_{h,i - 1} \right\rbrack}{z_{i} - z_{i - 1}} + \frac{\lambda\left\lbrack z_{h,i} \right\rbrack}{z_{h,i} - z_{i}} \right\rbrack

.. math::
   :label: eq-0305

   c_{i} = 0

.. math::
   :label: eq-0306

   r_{i} = T_{i}^{n} + \alpha\frac{\Delta t}{c_{i}\Delta z_{i}}\left( F_{i} - \alpha F_{i - 1} \right)

where

:math:`F_{i} = - \lambda\left\lbrack z_{h,i} \right\rbrack\left( \frac{\alpha T_{i}^{n} - T_{iB}^{n}}{z_{h,i} - z_{i}} \right)`,

.. math::
   :label: eq-0307

   F_{i - 1} = - \frac{\lambda\left\lbrack z_{h,i - 1} \right\rbrack}{z_{i} - z_{i - 1}}\left( T_{i - 1}^{n} - T_{i}^{n} \right)

For the interior snow/soil layers of all surfaces, :math:`snl + 1 < i < N_{nlevgrnd}`,

.. figure:: image15.png

.. math::
   :label: eq-0308

   a_{i} = - (1 - \alpha)\frac{\Delta t}{c_{i}\Delta z_{i}}\frac{\lambda\left\lbrack z_{h,i - 1} \right\rbrack}{z_{i} - z_{i - 1}}

.. math::
   :label: eq-0309

   b_{i} = 1 + (1 - \alpha)\frac{\Delta t}{c_{i}\Delta z_{i}}\left\lbrack \frac{\lambda\left\lbrack z_{h,i - 1} \right\rbrack}{z_{i} - z_{i - 1}} + \frac{\lambda\left\lbrack z_{h,i} \right\rbrack}{z_{i + 1} - z_{i}} \right\rbrack

.. math::
   :label: eq-0310

   c_{i} = - (1 - \alpha)\frac{\Delta t}{c_{i}\Delta z_{i}}\frac{\lambda\left\lbrack z_{h,i} \right\rbrack}{z_{i + 1} - z_{i}}

.. math::
   :label: eq-0311

   r_{i} = T_{i}^{n} + \alpha\frac{\Delta t}{c_{i}\Delta z_{i}}\left( F_{i} - F_{i - 1} \right)

The heating or cooling flux applied to the roof, and sunlit and shaded wall is



.. math::
   :label: eq-0071

   F_{heat} = \left\{ \begin{aligned}
   & \left| \alpha F_{i = N_{levgrnd}}^{n} + (1 - \alpha)F_{i = N_{levgrnd}}^{n + 1} \right| T_{iB} < T_{\min} \\
   & 0  T_{iB} \geq T_{\min}\{\}
   \end{aligned} \middle| \right\}





.. math::
   :label: eq-0072

   F_{cool} = \left\{ \begin{aligned}
   & \left| \alpha F_{i = N_{levgrnd}}^{n} + (1 - \alpha)F_{i = N_{levgrnd}}^{n + 1} \right| T_{iB} > T_{\min} \\
   & 0 T_{iB} \leq T_{\min}\{\}
   \end{aligned} \right\}



where

.. math::
   :label: eq-0312

   F_{i = N_{levgrnd}}^{n} = - \frac{\lambda\left\lbrack z_{h,i = N_{levgrnd}} \right\rbrack}{z_{h,i = N_{levgrnd}} - z_{i = N_{levgrnd}}}\left( T_{i = N_{levgrnd}}^{n} - T_{iB} \right)

.. math::
   :label: eq-0313

   F_{i = N_{levgrnd}}^{n + 1} = - \frac{\lambda\left\lbrack z_{h,i = N_{levgrnd}} \right\rbrack}{z_{h,i = N_{levgrnd}} - z_{i = N_{levgrnd}}}\left( T_{i = N_{levgrnd}}^{n + 1} - T_{iB} \right)

The total waste heat from space heating/air conditioning is



.. math::
   :label: eq-0073

   H_{wasteheat} = W_{roof}\left( f_{heat}F_{heat,roof} + f_{cool}F_{cool,roof} \right) +
   

\ 

.. math::
   :label: eq-0074

    \left( 1 - W_{roof} \right)\frac{H}{W}\left( \begin{aligned}
    & f_{heat}F_{heat,sunwall} + f_{cool}F_{cool,sunwall} + \\
    & f_{heat}F_{heat,shdwall} + f_{cool}F_{cool,shdwall}
   \end{aligned} \right) \leq H_{wasteheat,\max}


where :math:`f_{heat} = \frac{1}{0.75}` and :math:`f_{cool} = \frac{1}{0.25}` are factors describing the efficiency of space heating/air conditioning systems and :math:`H_{wasteheat,\max}` W m\ :sup:`-2` is a maximum limit on waste heat at any given time step. The heat removed by air conditioning is

.. math::
   :label: eq-0314

   H_{aircond} = F_{cool}

Phase Change
------------

Phase change may take place in any snow/soil layers of the pervious road and in the ponded water on roofs and impervious road. Note that the ponded water is treated as part of the top layer. Upon solution of the tridiagonal equation set (Press et al. 1992), the temperatures are evaluated to determine if phase change will take place as



.. math::
   :label: eq-0075

   T_{i}^{n + 1} > T_{f}\text{ and }w_{ice,i} > 0 i = snl + 1,\ldots,N_{levgrnd} \text{melting}
   

\ 

.. math::
   :label: eq-0076

   T_{i}^{n + 1} < T_{f}\text{ and }w_{liq,i} > 0 i = snl + 1,\ldots,0 \text{freezing}


\

.. math::
   :label: eq-0077

    T_{i}^{n + 1} < T_{f}\text{ and }w_{liq,i} > w_{liq,{max,} i} i = 1,\ldots,N_{levgrnd} \text{freezing}


where :math:`T_{i}^{n + 1}` is the layer temperature after solution of the tridiagonal equation set, :math:`w_{ice,i}` and :math:`w_{liq,i}` are the mass of ice and liquid water (kg m\ :sup:`-2`) in each layer, respectively, and :math:`T_{f}` is the freezing temperature of water (K) (:numref:`table-physical-constants`). For the freezing process in the layers of the pervious road, the concept of supercooled soil water from Niu and Yang (2006) is adopted. The supercooled soil water is the liquid water that coexists with ice over a wide range of temperatures below freezing and is implemented through a freezing point depression equation

.. math::
   :label: eq-0315

   w_{liq,{max,}i} = \Delta z_{i}\theta_{sat,i}\left\lbrack \frac{10^{3}L_{f}\left( T_{f} - T_{i} \right)}{gT_{i}\psi_{sat,i}} \right\rbrack^{\frac{- 1}{B_{i}}} T_{i} < T_{f}

where :math:`w_{liq,{max,}i}` is the maximum liquid water in layer :math:`i` (kg m\ :sup:`-2`) when the soil temperature :math:`T_{i}` is below the freezing temperature :math:`T_{f}`, :math:`L_{f}` is the latent heat of fusion (J kg\ :sup:`-1`) (:numref:`table-physical-constants`), :math:`g` is the gravitational acceleration (m s\ :sup:`-2`) (:numref:`table-physical-constants`), and :math:`\psi_{sat,i}` and :math:`B_{i}` are the soil texture-dependent saturated matric potential (mm) and Clapp and Hornberger (1978) exponent (section 5.3.1). Equation :eq:`eq-0315` applies to pervious road only, for roof and impervious road :math:`w_{liq,{max,}i} = 0`.

For the special case when snow is present (snow mass :math:`W_{sno} > 0`) but there are no explicit snow layers (:math:`snl = 0`) (i.e., there is not enough snow present to meet the minimum snow depth requirement of 0.01 m), snow melt will take place for soil layer :math:`i = 1` if the soil layer temperature is greater than the freezing temperature (:math:`T_{1}^{n + 1} > T_{f}`).

The rate of phase change is assessed from the energy excess (or deficit) needed to change :math:`T_{i}` to freezing temperature, :math:`T_{f}`. The excess or deficit of energy :math:`H_{i}` (W m\ :sup:`-2`) is determined as follows



.. math::
   :label: eq-0078

   H_{i} = \left\{ \begin{aligned}
   & h + \frac{\partial h}{\partial T}\left( T_{f} - T_{i}^{n} \right) + \alpha F_{i}^{n} + (1 - \alpha)F_{i}^{n + 1} \\
   & - \frac{c_{i}\Delta z_{i}}{\Delta t}\left( T_{f} - T_{i}^{n} \right) i = snl + 1 \\
   & \alpha\left( F_{i}^{n} - F_{i - 1}^{n} \right) + (1 - \alpha)\left( F_{i}^{n + 1} - F_{i - 1}^{n + 1} \right) \\
   & - \frac{c_{i}\Delta z_{i}}{\Delta t}\left( T_{f} - T_{i}^{n} \right) i = snl + 2,\ldots,N_{levgrnd}
   \end{aligned} \right\}



where :math:`F_{i}^{n + 1}` and :math:`F_{i - 1}^{n + 1}` are calculated from equations :eq:`eq-0290` and :eq:`eq-0300` using :math:`T_{i}^{n + 1}`. For roof and walls, :math:`F_{i = N_{levgrnd}}^{n + 1}` is calculated from equation :eq:`eq-0313`. If the melting criteria is met (equation (4.57)) and :math:`H_{m} = \frac{H_{i}\Delta t}{L_{f}} > 0`, then the ice mass is readjusted as

.. math::
   :label: eq-0316

   w_{ice,i}^{n + 1} = w_{ice,i}^{n} - H_{m} \geq 0 i = snl + 1,\ldots,N_{levgrnd}

If the freezing criteria is met (equation (4.57)) and :math:`H_{m} < 0`, then the ice mass is readjusted for :math:`i = snl + 1,\ldots,0` as

.. math::
   :label: eq-0317

   w_{ice,i}^{n + 1} = \min\left( w_{liq,i}^{n} + w_{ice,i}^{n},w_{ice,i}^{n} - H_{m} \right)

and for :math:`i = 1,\ldots,N_{levgrnd}` as



.. math::
   :label: eq-0079

   w_{ice,i}^{n + 1} = \left\{ \begin{aligned}
   & \min\left( w_{liq,i}^{n} + w_{ice,i}^{n} - w_{liq,{max,}i}^{n},w_{ice,i}^{n} - H_{m} \right) w_{liq,i}^{n} + w_{ice,i}^{n} \geq w_{liq,{max,}i}^{n}\  \\
   & 0 w_{liq,i}^{n} + w_{ice,i}^{n} < w_{liq,{max,}i}^{n}\text{  }
   \end{aligned} \right\}



Liquid water mass is readjusted as

.. math::
   :label: eq-0318

   w_{liq,i}^{n + 1} = w_{liq,i}^{n} + w_{ice,i}^{n} - w_{ice,i}^{n + 1} \geq 0

Because part of the energy :math:`H_{i}` may not be consumed in melting or released in freezing, the energy is recalculated as

.. math::
   :label: eq-0319

   H_{i*} = H_{i} - \frac{L_{f}\left( w_{ice,i}^{n} - w_{ice,i}^{n + 1} \right)}{\Delta t}

and this energy is used to cool or warm the layer (if :math:`\left| H_{i*} \right| > 0`) as



.. math::
   :label: eq-0080

   T_{i}^{n + 1} = \left\{ \begin{aligned}
   & T_{f} + \frac{\frac{\Delta t}{c_{i}\Delta z_{i}}H_{i*}}{\left( 1 - \frac{\Delta t}{c_{i}\Delta z_{i}}\frac{\partial h}{\partial T} \right)} i = snl + 1 \\
   & T_{f} + \frac{\Delta t}{c_{i}\Delta z_{i}}H_{i*} i = snl + 2,\ldots N_{levgrnd}
   \end{aligned} \right\}

.

For the special case when snow is present (:math:`W_{sno} > 0`), there are no explicit snow layers (:math:`snl = 0`), and :math:`\frac{H_{1}\Delta t}{L_{f}} > 0` (melting), the snow mass :math:`W_{sno}` (kg m\ :sup:`-2`) is reduced according to

.. math::
   :label: eq-0320

   W_{sno}^{n + 1} = W_{sno}^{n} - \frac{H_{1}\Delta t}{L_{f}} \geq 0

The snow depth is reduced proportionally

.. math::
   :label: eq-0321

   z_{sno}^{n + 1} = \frac{W_{sno}^{n + 1}}{W_{sno}^{n}}z_{sno}^{n}

Again, because part of the energy may not be consumed in melting, the energy for the surface layer :math:`i = 1` is recalculated as

.. math::
   :label: eq-0322

   H_{1*} = H_{1} - \frac{L_{f}\left( W_{sno}^{n} - W_{sno}^{n + 1} \right)}{\Delta t}

If there is excess energy (:math:`H_{1*} > 0`), this energy becomes available to the top layer as

.. math::
   :label: eq-0323

   H_{1} = H_{1*}

The ice mass, liquid water content, and temperature of the top layer are then determined from equations :eq:`eq-0316`, :eq:`eq-0318`, and :eq:`eq-0080` using the recalculated energy from equation :eq:`eq-0323`. Snow melt :math:`M_{1S}` (kg m\ :sup:`-2` s\ :sup:`-1`) and phase change energy :math:`E_{p,1S}` (W m\ :sup:`-2`) for this special case are

.. math::
   :label: eq-0324

   M_{1S} = \frac{W_{sno}^{n} - W_{sno}^{n + 1}}{\Delta t} \geq 0

.. math::
   :label: eq-0325

   E_{p,1S} = L_{f}M_{1S}

The total energy of phase change :math:`E_{p}` (W m\ :sup:`-2`) for the column is

.. math::
   :label: eq-0326

   E_{p} = E_{p,1S} + \sum_{i = snl + 1}^{i = N_{levgrnd}}E_{p,i}

where

.. math::
   :label: eq-0327

   E_{p,i} = L_{f}\frac{\left( w_{ice,i}^{n} - w_{ice,i}^{n + 1} \right)}{\Delta t}

The total snow melt :math:`M` (kg m\ :sup:`-2` s\ :sup:`-1`) is

.. math::
   :label: eq-0328

   M = M_{1S} + \sum_{i = snl + 1}^{i = 0}M_{i}

where

.. math::
   :label: eq-0329

   M_{i} = \frac{\left( w_{ice,i}^{n} - w_{ice,i}^{n + 1} \right)}{\Delta t} \geq 0

The solution for temperatures conserves energy as

.. math::
   :label: eq-0330

   G - E_{p} - \sum_{i = snl + 1}^{i = N_{levgrnd}}\frac{c_{i}\Delta z_{i}}{\Delta t}\left( T_{i}^{n + 1} - T_{i}^{n} \right) + \left\lbrack \alpha F_{i = N_{levgrnd}}^{n} + (1 - \alpha)F_{i = N_{levgrnd}}^{n + 1} \right\rbrack = 0

where :math:`G` is the ground heat flux (section 3.2.4) and the last term is the non-zero flux bottom boundary condition (roofs and walls only).

Thermal Properties
------------------

The thermal conductivities and heat capacities for roofs, walls, and :math:`i = 1,\ldots,N_{imprvrd}` layers of the impervious road are specified by the surface dataset as described in section 1.2.2 and :numref:`table-input-data` The :math:`i = N_{imprvrd} + 1,\ldots,N_{levgrnd}` layers of impervious road and the pervious road layers consist of soil or bedrock whose thermal properties are described below. In CLM4, organic matter modifies soil properties according to Lawrence and Slater (2008). Urban soils are assumed to have no organic matter so the equations below are shown in their reduced form. Note that the moisture content of the impervious road soil layers is maintained at zero.

Soil thermal conductivity :math:`\lambda_{i}` (W m\ :sup:`-1` K\ :sup:`-1`) is from Farouki (1981)



.. math::
   :label: eq-0081

   \lambda_{i} = \left\{ \begin{aligned}
   & K_{e,i}\lambda_{sat,i} + \left( 1 - K_{e,i} \right)\lambda_{dry,i} S_{r,i} > 1 \times 10^{- 7} \\
   & \lambda_{dry,i} S_{r,i} \leq 1 \times 10^{- 7}
   \end{aligned} \right\} i = 1,\ldots,N_{levsoi}
   \lambda_{i} = \lambda_{bedrock} i = N_{levsoi} + 1,\ldots N_{levgrnd}



where :math:`\lambda_{sat,i}` is the saturated thermal conductivity, :math:`\lambda_{dry,i}` is the dry thermal conductivity, :math:`K_{e,i}` is the Kersten number, :math:`S_{r,i}` is the wetness of the soil with respect to saturation, and :math:`\lambda_{bedrock} = 3` W m\ :sup:`-1` K\ :sup:`-1` is the thermal conductivity assumed for the deep ground layers (typical of saturated granitic rock; Clauser and Huenges, 1995). The saturated thermal conductivity :math:`\lambda_{sat,i}` (W m\ :sup:`-1` K\ :sup:`-1`) depends on the thermal conductivities of the soil solid, liquid water, and ice constituents



.. math::
   :label: eq-0082

   \lambda_{sat,i} = \left\{ \begin{aligned}
   & \lambda_{s,i}^{1 - \theta_{sat,i}}\lambda_{liq}^{\theta_{sat,i}} T_{i} \geq T_{f} \\
   & \lambda_{s,i}^{1 - \theta_{sat,i}}\lambda_{liq}^{\theta_{sat,i}}\lambda_{ice}^{\theta_{sat,i} - \theta_{liq,i}} T_{i} < T_{f}
   \end{aligned} \right\}



where the thermal conductivity of soil solids :math:`\lambda_{s,i}` varies with the sand and clay content

:math:`\lambda_{s,i} = \frac{8.80\ (\% sand)_{i} + \text{2.92 }(\% clay)_{i}}{(\% sand)_{i} + (\% clay)_{i}}`,

and :math:`\theta_{sat,i}` is the volumetric water content at saturation (porosity) (section 5.3.1). The thermal conductivity of dry natural soil :math:`\lambda_{dry,i}` (W m\ :sup:`-1` K\ :sup:`-1`) depends on the bulk density :math:`\rho_{d,i} = 2700\left( 1 - \theta_{sat,i} \right)` (kg m\ :sup:`-3`) as

.. math::
   :label: eq-0331

   \lambda_{dry,i} = \frac{0.135\rho_{d,i} + 64.7}{2700 - 0.947\rho_{d,i}}

The Kersten number :math:`K_{e,i}` is a function of the degree of saturation :math:`S_{r}` and phase of water



.. math::
   :label: eq-0083

   K_{e,i} = \left\{ \begin{aligned}
   & \log\left( S_{r,i} \right) + 1 \geq 0 T_{i} \geq T_{f} \\
   & S_{r,i} T_{i} < T_{f}
   \end{aligned} \right\}



where

.. math::
   :label: eq-0332

   S_{r,i} = \left( \frac{w_{liq,i}}{\rho_{liq}\Delta z_{i}} + \frac{w_{ice,i}}{\rho_{ice}\Delta z_{i}} \right)\frac{1}{\theta_{sat,i}} = \frac{\theta_{liq,i} + \theta_{ice,i}}{\theta_{sat,i}} \leq 1

Thermal conductivity :math:`\lambda_{i}` (W m\ :sup:`-1` K\ :sup:`-1`) for snow is from Jordan (1991)

.. math::
   :label: eq-0333

   \lambda_{i} = \lambda_{air} + \left( 7.75 \times 10^{- 5}\rho_{sno,i} + 1.105 \times 10^{- 6}\rho_{sno,i}^{2} \right)\left( \lambda_{ice} - \lambda_{air} \right)

where :math:`\lambda_{air}` and :math:`\lambda_{ice}` are the thermal conductivities of air and ice (:numref:`table-physical-constants`) and :math:`\rho_{sno,i}` is the bulk density of snow (kg m\ :sup:`-3`)

.. math::
   :label: eq-0334

   \rho_{sno,i} = \frac{w_{ice,i} + w_{liq,i}}{\Delta z_{i}}

The volumetric heat capacity :math:`c_{i}` (J m\ :sup:`-3` K\ :sup:`-1`) for soil is from de Vries (1963) and depends on the heat capacities of the soil solid, liquid water, and ice constituents

.. math::
   :label: eq-0335

   c_{i} = c_{s,i}\left( 1 - \theta_{sat,i} \right) + \frac{w_{ice,i}}{\Delta z_{i}}C_{ice} + \frac{w_{liq,i}}{\Delta z_{i}}C_{liq}

where the heat capacity of soil solids :math:`c_{s,i}` (J m\ :sup:`-3` K\ :sup:`-1`) is



.. math::
   :label: eq-0084

   c_{s,i} = \left( \frac{2.128\ (\% sand)_{i} + \text{2.385 }(\% clay)_{i}}{(\% sand)_{i} + (\% clay)_{i}} \right) \times 10^{6} i = 1,\ldots,N_{levsoi}
   c_{s,i} = c_{s,bedrock} i = N_{levsoi} + 1,\ldots,N_{levgrnd}



and :math:`C_{liq}` and :math:`C_{ice}` are the specific heat capacities (J kg\ :sup:`-1` K\ :sup:`-1`) of liquid water and ice, respectively (:numref:`table-physical-constants`) and :math:`c_{s,bedrock} = 2 \times 10^{6}` J m\ :sup:`-3` K\ :sup:`-1` is the heat capacity of bedrock. For snow

.. math::
   :label: eq-0336

   c_{i} = \frac{w_{ice,i}}{\Delta z_{i}}C_{ice} + \frac{w_{liq,i}}{\Delta z_{i}}C_{liq}

For the special case when snow is present (:math:`W_{sno} > 0`) but there are no explicit snow layers (:math:`snl = 0`), the heat capacity of the top layer is a blend of ice and soil heat capacity

.. math::
   :label: eq-0337

   c_{1} = c_{1}^{*} + \frac{C_{ice}W_{sno}}{\Delta z_{1}}

where :math:`c_{1}^{*}` is calculated from equation :eq:`eq-0335`.

Hydrology
=========

The hydrology for the pervious road generally follows that of CLM4 for bare soil surfaces and includes snow accumulation and melt, water transfer between snow layers, infiltration, evaporation, surface runoff, sub-surface drainage, redistribution within the soil column, and groundwater discharge and recharge to simulate changes in snow water :math:`\Delta W_{sno}`, soil water :math:`\Delta w_{liq,i}`, soil ice :math:`\Delta w_{ice,i}`, and water in the unconfined aquifer :math:`\Delta W_{a}` (all in kg m\ :sup:`-2` or mm of H\ :sub:`2`\ O) (:numref:`fig-hydrologic-processes-pervious-road`). The water balance of the pervious road is



.. math::
   :label: eq-0085

   \Delta W_{sno} + \sum_{i = 1}^{N_{levsoi}}{\left( \Delta w_{liq,i} + \Delta w_{ice,i} \right) + \Delta W_{a} = \left( \begin{aligned}
   & q_{rain} + q_{sno} - E_{prvrd} - q_{over} - q_{drai} \\
   & - q_{rgwl} - q_{snwcp,ice}
   \end{aligned} \right)}\Delta t



where :math:`q_{rain}` is liquid part of precipitation, :math:`q_{sno}` is solid part of precipitation, :math:`E_{prvrd}` is the total evaporation (chapter 3), :math:`q_{over}` is surface runoff (section 5.2), :math:`q_{drai}` is sub-surface drainage (section 5.4), :math:`q_{rgwl}` and :math:`q_{snwcp,ice}` are liquid and solid runoff due to snow capping (section 5.5) (all in kg m\ :sup:`-2` s\ :sup:`-1`), :math:`N_{levsoi}` is the number of soil layers, and :math:`\Delta t` is the time step (s). In general, snow capping will not be invoked for urban areas, but is described here for completeness.

.. figure:: image16.jpeg
   :width: 5.60417in
   :height: 4.82292in
   :name: fig-hydrologic-processes-pervious-road

   Hydrologic processes simulated for the pervious road. Evaporation is supplied by all soil layers. An unconfined aquifer is added to the bottom of the soil column. The depth to the water table is :math:`z_{\nabla}` (m). Changes in aquifer water content :math:`W_{a}` (mm) are controlled by the balance between drainage from the aquifer water :math:`q_{drai}` and the aquifer recharge rate :math:`q_{recharge}` (kg m\ :sup:`-2` s\ :sup:`-1`) (defined as positive from soil to aquifer).

The roof and the impervious road are hydrologically inactive except for their capacity to intercept, store, and evaporate a limited amount of liquid precipitation (1 kg m\ :sup:`-2`), and snow. Logistically, the storage of liquid precipitation is accounted for in the top layer :math:`i = 1`. The water in excess of this storage capacity is routed to surface runoff. These surfaces are also allowed to intercept solid precipitation (snow) and store this until the snowpack is melted or sublimated. No sub-surface drainage is allowed. The water balance of the roof and impervious road is

.. math::
   :label: eq-0338

   \Delta W_{sno} + \Delta w_{liq,1} + \Delta w_{ice,1} = \left( q_{rain} + q_{sno} - E_{roof} - q_{over} - q_{rgwl} - q_{snwcp,ice} \right)\Delta t

.. math::
   :label: eq-0339

   \Delta W_{sno} + \Delta w_{liq,1} + \Delta w_{ice,1} = \left( q_{rain} + q_{sno} - E_{imprvrd} - q_{over} - q_{rgwl} - q_{snwcp,ice} \right)\Delta t

where :math:`\Delta w_{liq,1}` and :math:`\Delta w_{ice,1}` are the liquid water and ice stored on the top of the urban surface. The sunlit and shaded walls are hydrologically inactive.

The rate of liquid and solid precipitation reaching the urban surface (kg m\ :sup:`-2` s\ :sup:`-1`) is

.. math::
   :label: eq-0340

   q_{grnd,liq} = q_{rain}

.. math::
   :label: eq-0341

   q_{grnd,ice} = q_{sno}

Solid precipitation reaching the surface, :math:`q_{grnd,ice}\Delta t`, is added immediately to the snow pack (section 5.1). The liquid part, :math:`q_{grnd,liq}\Delta t` is added after surface fluxes, temperatures, soil water, and runoff have been determined.

Snow
----

The parameterizations for snow are based primarily on Anderson (1976), Jordan (1991), and Dai and Zeng (1997). Snow can have up to five layers. These layers are indexed in the Fortran code as :math:`i = - 4, - 3, - 2, - 1,0` where layer :math:`i = 0` is the snow layer next to the top soil layer and layer :math:`i = - 4` is the top layer of a five-layer snow pack. Since the number of snow layers varies according to the snow depth, we use the notation :math:`snl + 1` to describe the top layer of snow for the variable layer snow pack, where :math:`snl` is the negative of the number of snow layers. Refer to :numref:`fig-snow-pack-example` for an example of the snow layer structure for a three layer snow pack.

.. figure:: image17.jpeg
   :width: 6in
   :height: 4.80208in
   :name: fig-snow-pack-example

   Example of three layer snow pack (:math:`snl = - 3`). Shown are three snow layers, :math:`i = - 2`, :math:`i = - 1`, and :math:`i = 0`. The layer node depth is :math:`z`, the layer interface is :math:`z_{h}`, and the layer thickness is :math:`\Delta z`.

The state variables for snow are the mass of water :math:`w_{liq,i}` (kg m\ :sup:`-2`), mass of ice :math:`w_{ice,i}` (kg m\ :sup:`-2`), layer thickness :math:`\Delta z_{i}` (m), and temperature :math:`T_{i}` (chapter 4). The water vapor phase is neglected. Snow can also exist in the model without being represented by explicit snow layers. This occurs when the snowpack is less than a specified minimum snow depth (:math:`z_{sno} < 0.01` m). In this case, the state variable is the mass of snow :math:`W_{sno}` (kg m\ :sup:`-2`).

The next two sections (5.1.1 and 5.1.2) describe the ice and water content of the snow pack assuming that at least one snow layer exists. See section 5.1.3 for a description of how a snow layer is initialized. Snow compaction is described in section 5.1.4 and snow layer combination and subdivision in section 5.1.5.

Ice Content
~~~~~~~~~~~

The conservation equation for mass of ice in snow layers is



.. math::
   :label: eq-0086

   \frac{\partial w_{ice,i}}{\partial t} = \left\{ \begin{aligned}
   & q_{ice,i - 1} - \frac{\left( \Delta w_{ice,i} \right)_{p}}{\Delta t} i = snl + 1 \\
   & - \frac{\left( \Delta w_{ice,i} \right)_{p}}{\Delta t} i = snl + 2,\ldots,0
   \end{aligned} \right\}



where :math:`q_{ice,i - 1}` is the rate of ice accumulation from precipitation or frost or the rate of ice loss from sublimation (kg m\ :sup:`-2` s\ :sup:`-1`) in the top layer and :math:`\frac{\left( \Delta w_{ice,i} \right)_{p}}{\Delta t}` is the change in ice due to phase change (melting rate) (section 4.2). The term :math:`q_{ice,i - 1}` is calculated in two steps as

.. math::
   :label: eq-0342

   q_{ice,i - 1} = q_{grnd,ice} + \left( q_{frost} - q_{subl} \right)

where :math:`q_{grnd,ice}` is the rate of solid precipitation reaching the surface and :math:`q_{frost}` and :math:`q_{subl}` are gains due to frost and losses due to sublimation, respectively (section 3.2.4). In the first step, a new snow depth :math:`z_{sno}` (m) is calculated from

.. math::
   :label: eq-0343

   z_{sno}^{n + 1} = z_{sno}^{n} + \Delta z_{sno}

where

.. math::
   :label: eq-0344

   \Delta z_{sno} = \frac{q_{grnd,ice}\Delta t}{\rho_{sno}}

and :math:`\rho_{sno}` is the bulk density of newly fallen snow (kg m\ :sup:`-3`) (Anderson 1976)



.. math::
   :label: eq-0087

   \rho_{sno} = \left\{ \begin{aligned}
   & 50 + 1.7(17)^{1.5} T_{atm} > T_{f} + 2 \\
   & 50 + 1.7\left( T_{atm} - T_{f} + 15 \right)^{1.5} T_{f} - 15 < T_{atm} \leq T_{f} + 2 \\
   & 50 T_{atm} \leq T_{f} - 15
   \end{aligned} \right\}



where :math:`T_{atm}` is the atmospheric temperature (K), and :math:`T_{f}` is the freezing temperature of water (K) (:numref:`table-physical-constants`). The mass of snow :math:`W_{sno}` is

.. math::
   :label: eq-0345

   W_{sno}^{n + 1} = W_{sno}^{n} + q_{grnd,ice}\Delta t

The ice content of the top layer and the layer thickness are updated as

.. math::
   :label: eq-0346

   w_{ice,snl + 1}^{n + 1} = w_{ice,snl + 1}^{n} + q_{grnd,ice}\Delta t

.. math::
   :label: eq-0347

   \Delta z_{snl + 1}^{n + 1} = \Delta z_{snl + 1}^{n} + \Delta z_{sno}

In the second step, after surface fluxes and temperatures have been determined (chapters 3 and 4), :math:`w_{ice,snl + 1}` is updated for frost or sublimation as

.. math::
   :label: eq-0348

   w_{ice,snl + 1}^{n + 1} = w_{ice,snl + 1}^{n} + \left( q_{frost} - q_{subl} \right)\Delta t

If :math:`w_{ice,snl + 1}^{n + 1} < 0` upon solution of equation :eq:`eq-0348`, the ice content is reset to zero and the liquid water content :math:`w_{liq,snl + 1}` is reduced by the amount required to bring :math:`w_{ice,snl + 1}^{n + 1}` up to zero. The snow water equivalent :math:`W_{sno}` is capped to not exceed 1000 kg m\ :sup:`-2`. If the addition of :math:`q_{frost}` were to result in :math:`W_{sno} > 1000` kg m\ :sup:`-2`, the frost term :math:`q_{frost}` is instead added to the ice runoff term :math:`q_{snwcp,ice}` (section 5.5).

Water Content
~~~~~~~~~~~~~

The conservation equation :eq:`eq-0348` for mass of water in snow layers is

.. math::
   :label: eq-0349

   \frac{\partial w_{liq,i}}{\partial t} = \left( q_{liq,i - 1} - q_{liq,i} \right) + \frac{\left( \Delta w_{liq,i} \right)_{p}}{\Delta t}

where :math:`q_{liq,i - 1}` is the flow of liquid water into layer :math:`i` from the layer above, :math:`q_{liq,i}` is the flow of water out of layer :math:`i` to the layer below, :math:`\frac{\left( \Delta w_{liq,i} \right)_{p}}{\Delta t}` is the change in liquid water due to phase change (melting rate) (section 4.2). For the top snow layer only,

.. math::
   :label: eq-0350

   q_{liq,i - 1} = q_{grnd,liq} + \left( q_{sdew} - q_{seva} \right)

where :math:`q_{grnd,liq}` is the rate of liquid precipitation reaching the snow, :math:`q_{seva}` is the evaporation of liquid water and :math:`q_{sdew}` is the liquid dew (section 3.2.4). After surface fluxes and temperatures have been determined (chapters 3 and 4), :math:`w_{liq,snl + 1}` is updated for the liquid precipitation reaching the ground and dew or evaporation as

.. math::
   :label: eq-0351

   w_{liq,snl + 1}^{n + 1} = w_{liq,snl + 1}^{n} + \left( q_{grnd,liq} + q_{sdew} - q_{seva} \right)\Delta t

When the liquid water within a snow layer exceeds the layer's holding capacity, the excess water is added to the underlying layer, limited by the effective porosity (:math:`1 - \theta_{ice}`) of the layer. The flow of water is assumed to be zero (:math:`q_{liq,i} = 0`) if the effective porosity of either of the two layers (:math:`1 - \theta_{ice,i}\text{ and }1 - \theta_{ice,i + 1}`) is less than :math:`\theta_{imp} = 0.05`, the water impermeable volumetric water content. Thus, water flow between layers, :math:`q_{liq,i}`, for :math:`i = snl + 1,\ldots,0` is initially calculated as

.. math::
   :label: eq-0352

   q_{liq,i} = \frac{\rho_{liq}\left\lbrack \theta_{liq,i} - S_{r}\left( 1 - \theta_{ice,i} \right) \right\rbrack\Delta z_{i}}{\Delta t} \geq 0

where the volumetric liquid water :math:`\theta_{liq,i}` and ice :math:`\theta_{ice,i}` contents are

.. math::
   :label: eq-0353

   \theta_{ice,i} = \frac{w_{ice,i}}{\Delta z_{i}\rho_{ice}} \leq 1

:math:`\theta_{liq,i} = \frac{w_{liq,i}}{\Delta z_{i}\rho_{liq}} \leq 1 - \theta_{ice,i}`,

and :math:`S_{r} = 0.033` is the irreducible water saturation (snow holds a certain amount of liquid water due to capillary retention after drainage has ceased (Anderson 1976)). The water holding capacity of the underlying layer limits the flow of water :math:`q_{liq,i}` calculated in Equations :eq:`eq-0352`-:eq:`eq-0355` , unless the underlying layer is the surface layer, as

.. math::
   :label: eq-0354

   q_{liq,i} \leq \frac{\rho_{liq}\left\lbrack 1 - \theta_{ice,i + 1} - \theta_{liq,i + 1} \right\rbrack\Delta z_{i + 1}}{\Delta t} i = snl + 1,\ldots, - 1

The volumetric liquid water content :math:`\theta_{liq,i}` is updated as

.. math::
   :label: eq-0355

   \theta_{liq,i}^{n + 1} = \theta_{liq,i}^{n} + \left( q_{i - 1} - q_{i} \right)\Delta t

equation :eq:`eq-0352` are solved sequentially from top (:math:`i = snl + 1`) to bottom (:math:`i = 0`) snow layer in each time step. The total flow of liquid water reaching the urban surface is then :math:`q_{liq,0}`.

Initialization of snow layer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If there are no existing snow layers (:math:`snl + 1 = 1`) but :math:`z_{sno} \geq 0.01` m after accounting for solid precipitation :math:`q_{sno}`, then a snow layer is initialized (:math:`snl = - 1`) as follows



.. math::
   :label: eq-0088

   \Delta z_{0} = z_{sno}
   

\ 

.. math::
   :label: eq-0089

   \begin{cases}
   z_{o} = - 0.5\Delta z_{0} \\
   z_{h, - 1} = - \Delta z_{0} \\
   T_{0} = \min\left( T_{f},T_{atm} \right) \\
   w_{ice,0} = W_{sno}
   \end{cases}


\

.. math::
   :label: eq-0090

   w_{liq,0} = 0

Snow Compaction
~~~~~~~~~~~~~~~

Snow compaction is initiated after the hydrology calculations [surface runoff (section 5.2), infiltration (section 5.2), soil water (section 5.3), groundwater-soilwater interactions (section 5.4)] are complete. Compaction of snow includes three types of processes: destructive metamorphism of new snow (crystal breakdown due to wind or thermodynamic stress); snow load or overburden (pressure); and melting (changes in snow structure due to melt-freeze cycles plus changes in crystals due to liquid water). The total fractional compaction rate for each snow layer :math:`C_{R,i}` (s\ :sup:`-1`) is the sum of the three compaction processes

.. math::
   :label: eq-0356

   C_{R,i} = \frac{1}{\Delta z_{i}}\frac{\partial\Delta z_{i}}{\partial t} = C_{R1,i} + C_{R2,i} + C_{R3,i}

Compaction is not allowed if the layer is saturated

.. math::
   :label: eq-0357

   1 - \left( \frac{w_{ice,i}}{\Delta z_{i}\rho_{ice}} + \frac{w_{liq,i}}{\Delta z_{i}\rho_{liq}} \right) \leq 0.001

or if the ice content is below a minimum value (:math:`w_{ice,i} \leq 0.1`).

Compaction as a result of destructive metamorphism :math:`C_{R1,\mspace{6mu} i}` (s\ :sup:`-1`) is temperature dependent (Anderson 1976)

.. math::
   :label: eq-0358

   C_{R1,i} = \left\lbrack \frac{1}{\Delta z_{i}}\frac{\partial\Delta z_{i}}{\partial t} \right\rbrack_{metamorphism} = - c_{3}c_{1}c_{2}\exp\left\lbrack - c_{4}\left( T_{f} - T_{i} \right) \right\rbrack

where :math:`c_{3} = 2.777 \times 10^{- 6}` (s\ :sup:`-1`) is the fractional compaction rate for :math:`T_{i} = T_{f}`, :math:`c_{4} = 0.04` K\ :sup:`-1`, and



.. math::
   :label: eq-0091

   c_{1} = 1 \frac{w_{ice,i}}{\Delta z_{i}} \leq 100\text{ kg }\text{m}^{\text{-3}}
   

\ 

.. math::
   :label: eq-0092

   \begin{cases}
   c_{1} = \exp\left\lbrack - 0.046\left( \frac{w_{ice,i}}{\Delta z_{i}} - 100 \right) \right\rbrack & \frac{w_{ice,i}}{\Delta z_{i}} > 100\text{ kg }\text{m}^{\text{-3}} \\
   c_{2} = 2 & \frac{w_{liq,i}}{\Delta z_{i}} > 0.01
   \end{cases}

\ :math:`c_{2} = 1 \frac{w_{liq,i}}{\Delta z_{i}} \leq 0.01`

where :math:`\frac{w_{ice,i}}{\Delta z_{i}}` and :math:`\frac{w_{liq,i}}{\Delta z_{i}}` are the bulk densities of liquid water and ice (kg m\ :sup:`-3`).

The compaction rate as a result of overburden :math:`C_{R2,\mspace{6mu} i}` (s\ :sup:`-1`) is a linear function of the snow load pressure :math:`P_{s,i}` (kg m\ :sup:`-2`) (Anderson 1976)

.. math::
   :label: eq-0359

   C_{R2,i} = \left\lbrack \frac{1}{\Delta z_{i}}\frac{\partial\Delta z_{i}}{\partial t} \right\rbrack_{overburden} = - \frac{P_{s,i}}{\eta}

where :math:`\eta` is a viscosity coefficient (kg s m\ :sup:`-2`) that varies with density and temperature as

.. math::
   :label: eq-0360

   \eta = \eta_{0}\exp\left\lbrack c_{5}\left( T_{f} - T_{i} \right) + c_{6}\frac{w_{ice,i}}{\Delta z_{i}} \right\rbrack

where :math:`\eta_{0} = 9 \times 10^{5}` kg s m\ :sup:`-2`, and :math:`c_{5} = 0.08` K\ :sup:`-1`, :math:`c_{6} = 0.023` m\ :sup:`3` kg\ :sup:`-1` are constants. The snow load pressure :math:`P_{s,i}` is calculated for each layer as the sum of the ice :math:`w_{ice,i}` and liquid water contents :math:`w_{liq,i}` of the layers above plus half the ice and liquid water contents of the layer being compacted

.. math::
   :label: eq-0361

   P_{s,i} = \frac{\left( w_{ice,i} + w_{liq,i} \right)}{2} + \sum_{j = snl + 1}^{j = i - 1}\left( w_{ice,j} + w_{liq,j} \right)

The compaction rate due to melting :math:`C_{R3,\mspace{6mu} i}` (s\ :sup:`-1`) is taken to be the ratio of the change in snow ice fraction after the melting to the fraction before melting

.. math::
   :label: eq-0362

   C_{R3,i} = \left\lbrack \frac{1}{\Delta z_{i}}\frac{\partial\Delta z_{i}}{\partial t} \right\rbrack_{melt} = - \frac{1}{\Delta t}\max\left( 0,\frac{f_{ice,i}^{n} - f_{ice,i}^{n + 1}}{f_{ice,i}^{n}} \right)

where the fraction of ice :math:`f_{ice,i}` is

.. math::
   :label: eq-0363

   f_{ice,i} = \frac{w_{ice,i}}{w_{ice,i} + w_{liq,i}}

and melting is identified during the phase change calculations (section 4.2).

The snow layer thickness after compaction is then

.. math::
   :label: eq-0364

   \Delta z_{i}^{n + 1} = \Delta z_{i}^{n}\left( 1 + C_{R,i}\Delta t \right)

Snow Layer Combination and Subdivision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After the determination of snow temperature including phase change (chapter 4), snow hydrology (sections 5.2.1, 5.2.2, and 5.2.3), and the compaction calculations (5.2.4), the number of snow layers is adjusted by either combining or subdividing layers. The combination and subdivision of snow layers is based on Jordan (1991).

Combination
^^^^^^^^^^^

If a snow layer has nearly melted or if its thickness :math:`\Delta z_{i}` is less than the prescribed minimum thickness :math:`\Delta z_{\min}` (:numref:`table-snow-layer-thickness`), the layer is combined with a neighboring layer. The overlying or underlying layer is selected as the neighboring layer according to the following rules

-  If the top layer is being removed, it is combined with the underlying layer

-  If the underlying layer is not snow, the layer is combined with the overlying layer

-  If the layer is nearly completely melted, the layer is combined with the underlying layer

-  If none of the above rules apply, the layer is combined with the thinnest neighboring layer.

A first pass is made through all snow layers to determine if any layer is nearly melted (:math:`w_{ice,i} \leq 0.1`). If so, the remaining liquid water and ice content of layer :math:`i` is combined with the underlying neighbor :math:`i + 1` as

.. math::
   :label: eq-0365

   w_{liq,i + 1} = w_{liq,i + 1} + w_{liq,i}

.. math::
   :label: eq-0366

   w_{ice,i + 1} = w_{ice,i + 1} + w_{ice,i}

This includes the snow layer directly above the urban surface. In this case, the liquid water and ice content of the melted snow layer is added as ponded water/ice on the urban surface layer. The layer properties, :math:`T_{i}`, :math:`w_{ice,i}`, :math:`w_{liq,i}`, :math:`\Delta z_{i}`, are then re-indexed so that the layers above the eliminated layer are shifted down by one and the number of snow layers is decremented accordingly.

At this point, if there are no explicit snow layers remaining (:math:`snl = 0`), the snow water equivalent :math:`W_{sno}` and snow depth :math:`z_{sno}` are set to zero, otherwise, :math:`W_{sno}` and :math:`z_{sno}` are re-calculated as

.. math::
   :label: eq-0367

   W_{sno} = \sum_{i = snl + 1}^{i = 0}\left( w_{ice,i} + w_{liq,i} \right)

.. math::
   :label: eq-0368

   z_{sno} = \sum_{i = snl + 1}^{i = 0}{\Delta z_{i}}

If the snow depth :math:`0 < z_{sno} < 0.01` m, the number of snow layers is set to zero, the total ice content of the snowpack :math:`\sum_{i = snl + 1}^{i = 0}w_{ice,\mspace{6mu} i}` is assigned to :math:`W_{sno}`, and the total liquid water :math:`\sum_{i = snl + 1}^{i = 0}w_{liq,\mspace{6mu} i}` is assigned to the urban surface layer. Otherwise, the layers are combined according to the rules above.

When two snow layers are combined (denoted here as :math:`1` and :math:`2`), their thickness combination (:math:`c`) is

:math:`\Delta z_{c} = \Delta z_{1} + \Delta z_{2}`,

their mass combination is

.. math::
   :label: eq-0369

   w_{liq,c} = w_{liq,1} + w_{liq,2}

:math:`w_{ice,c} = w_{ice,1} + w_{ice,2}`,

and their temperatures are combined as

.. math::
   :label: eq-0370

   T_{c} = T_{f} + \frac{h_{c} - L_{f}w_{liq,c}}{C_{ice}w_{ice,c} + C_{liq}w_{liq,c}}

where :math:`h_{c} = h_{1} + h_{2}` is the combined enthalpy :math:`h_{i}` of the two layers where

.. math::
   :label: eq-0371

   h_{i} = \left( C_{ice}w_{ice,i} + C_{liq}w_{liq,i} \right)\left( T_{i} - T_{f} \right) + L_{f}w_{liq,i}

In these equations, :math:`L_{f}` is the latent heat of fusion (J kg\ :sup:`-1`) and :math:`C_{liq}` and :math:`C_{ice}` are the specific heat capacities (J kg\ :sup:`-1` K\ :sup:`-1`) of liquid water and ice, respectively (:numref:`table-physical-constants`). After layer combination, the node depths and layer interfaces (:numref:`fig-snow-pack-example`) are recalculated from

.. math::
   :label: eq-0372

   z_{i} = z_{h,i} - 0.5\Delta z_{i} i = 0,\ldots,snl + 1

.. math::
   :label: eq-0373

   z_{h,i - 1} = z_{h,i} - \Delta z_{i} i = 0,\ldots,snl + 1

where :math:`\Delta z_{i}` is the layer thickness.

.. list-table:: Minimum and maximum thickness of snow layers (m)
   :widths: auto
   :header-rows: 1
   :name: table-snow-layer-thickness


   * - Layer
     - :math:`\Delta z_{\min}`
     - :math:`N_{l}`
     - :math:`N_{u}`
     - :math:`\left( \Delta z_{\max}()_{l} \right)`
     - :math:`\left( \Delta z_{\max}()_{u} \right)`
   * - 1 (top)
     - 0.010
     - 1
     - >1
     - 0.03
     - 0.02
   * - 2
     - 0.015
     - 2
     - >2
     - 0.07
     - 0.05
   * - 3
     - 0.025
     - 3
     - >3
     - 0.18
     - 0.11
   * - 4
     - 0.055
     - 4
     - >4
     - 0.41
     - 0.23
   * - 5 (bottom)
     - 0.115
     - 5
     - -
     - -
     - -


The maximum snow layer thickness, :math:`\Delta z_{\max}`, depends on the number of layers, :math:`N_{l}` and :math:`N_{u}`.

Subdivision
^^^^^^^^^^^

The snow layers are subdivided when the layer thickness exceeds a prescribed maximum thickness :math:`\Delta z_{\max}` with lower and upper bounds that depend on the number of snow layers (:numref:`table-snow-layer-thickness`). For example, if there is only one layer, then the maximum thickness of that layer is 0.03 m, however, if there is more than one layer, then the maximum thickness of the top layer is 0.02 m. Layers are checked sequentially from top to bottom for this limit. If there is only one snow layer and its thickness is greater than 0.03 m (:numref:`table-snow-layer-thickness`), the layer is subdivided into two layers of equal thickness, liquid water and ice contents, and temperature. If there is an existing layer below the layer to be subdivided, the thickness :math:`\Delta z_{i}`, liquid water and ice contents, :math:`w_{liq,\mspace{6mu} i}` and :math:`w_{ice,\mspace{6mu} i}`, and temperature :math:`T_{i}` of the excess snow are combined with the underlying layer according to equations (5.38)-:eq:`eq-0370`. If there is no underlying layer after adjusting the layer for the excess snow, the layer is subdivided into two layers of equal thickness, liquid water and ice contents. The vertical snow temperature profile is maintained by calculating the slope between the layer above the splitting layer (:math:`T_{1}`) and the splitting layer (:math:`T_{2}`) and constraining the new temperatures (:math:`T_{2}^{n + 1}`, :math:`T_{3}^{n + 1}`) to lie along this slope. The temperature of the lower layer is first evaluated from

:math:`T_{3}^{'} = T_{2}^{n} - \left( \frac{T_{1}^{n} - T_{2}^{n}}{\frac{\left( \Delta z_{1}^{n} + \Delta z_{2}^{n} \right)}{2}} \right)\left( \frac{\Delta z_{2}^{n + 1}}{2} \right)`,

then adjusted as,



.. math::
   :label: eq-0093

   T_{3}^{n + 1} = T_{2}^{n} T_{3}^{'} \geq T_{f}
   T_{2}^{n + 1} = T_{2}^{n} + \left( \frac{T_{1}^{n} - T_{2}^{n}}{\frac{\left( \Delta z_{1} + \Delta z_{2}^{n} \right)}{2}} \right)\left( \frac{\Delta z_{2}^{n + 1}}{2} \right) T_{3}^{'} < T_{f}



where here the subscripts 1, 2, and 3 denote three layers numbered from top to bottom. After layer subdivision, the node depths and layer interfaces are recalculated from equations :eq:`eq-0372` and :eq:`eq-0373`.

Surface Runoff and Infiltration
-------------------------------

For the roof and impervious road, water on these surfaces in excess of a maximum ponding limit :math:`w_{pond,\max}` (kg m\ :sup:`-2`) is routed to surface runoff as



.. math::
   :label: eq-0094

   q_{over} = \frac{w_{liq,1}}{\Delta t} + q_{liq,0} - q_{seva} - \frac{w_{pond,\max}}{\Delta t}
   q_{over} = q_{liq,0} snl < 0


where :math:`q_{liq,0}` is the rate of liquid water reaching the surface from rain (section 5.1) and/or snowmelt (section 5.1.2) and :math:`q_{seva}` is the evaporation of liquid water from the top layer (section 3.4). The liquid water content of the top layer is adjusted to



.. math::
   :label: eq-0095

   w_{liq,1} = {wover}_{pond,\max}
   w_{liq,1} = w_{liq,1} + \left( q_{liq,0} - q_{seva} \right)\Delta t \geq 0 q_{over} = 0


For the pervious road, the simple TOPMODEL-based (Beven and Kirkby 1979) runoff model (SIMTOP) described by Niu et al. (2005) is implemented. A key concept underlying this approach is that of fractional saturated/impermeable area :math:`f_{sat}`, which is determined by the topographic characteristics and soil moisture state of a grid cell. The surface runoff consists of overland flow due to saturation excess (Dunne runoff) and infiltration excess (Hortonian runoff) mechanisms

.. math::
   :label: eq-0374

   q_{over} = f_{sat}q_{liq,0} + \left( 1 - f_{sat} \right)\max\left( 0,q_{liq,0} - q_{infl,\max}() \right)

where :math:`q_{liq,0}` is liquid precipitation reaching the ground plus any melt water from snow (kg m\ :sup:`-2` s\ :sup:`-1`) and :math:`q_{infl,\max}` is a maximum soil infiltration capacity (kg m\ :sup:`-2` s\ :sup:`-1`). In Niu et al. (2005), :math:`f_{sat}` was a function of soil moisture whose potential or maximum value,\ :math:`f_{\max}`, was solely determined by topographic characteristics. Niu and Yang (2006) modified the expression for :math:`f_{sat}` to include a dependence on impermeable area fraction in frozen soil, :math:`f_{frz,1}`, of the top :math:`i = 1` soil layer as

.. math::
   :label: eq-0375

   f_{sat} = \left( 1 - f_{frz,1} \right){f{\exp\left( - 0.5f_{over}z_{\nabla} \right)}_{frz,1}}_{\max}

where :math:`f_{\max}` is the maximum saturated fraction, :math:`f_{over}` is a decay factor (m\ :sup:`-1`), and :math:`z_{\nabla}` is the water table depth (m) (section 5.4). The maximum saturated fraction, :math:`f_{\max}`, is defined as the discrete cumulative distribution function (CDF) of the topographic index when the grid cell mean water table depth is zero. Thus, :math:`f_{\max}` is the percent of pixels in a grid cell whose topographic index is larger than or equal to the grid cell mean topographic index. It should be calculated explicitly from the CDF at each grid cell at the resolution that the model is run. However, because this is a computationally intensive task for global applications, :math:`f_{\max}` is calculated once from the CDF at a spatial resolution of 0.5° by 0.5° following Niu et al. (2005) and then area-averaged to the desired resolution. The 0.5° resolution is compatible with the resolution of other CLM input surface datasets (e.g., plant functional types, leaf area index). The decay factor :math:`f_{over}` for global simulations was determined through sensitivity analysis and comparison with observed runoff to be 0.5 m\ :sup:`-1`.

The impermeable fraction :math:`f_{frz,i}` is parameterized as a function of soil ice content (Niu and Yang 2006)

.. math::
   :label: eq-0376

   f_{frz,i} = \frac{\exp\left\lbrack - \alpha\left( 1 - \frac{w_{ice,i}}{w_{ice,i} + w_{liq,i}} \right) \right\rbrack - \exp( - \alpha)}{1 - \exp( - \alpha)}

where :math:`\alpha = 3` is an adjustable scale-dependent parameter, and :math:`w_{ice,i}` and :math:`w_{liq,i}` are the ice and liquid water contents of soil layer :math:`i` (kg m\ :sup:`-2`).

The maximum infiltration capacity :math:`q_{infl,\max}` in equation :eq:`eq-0374` is determined from soil texture and soil moisture (Entekhabi and Eagleson 1989) as

.. math::
   :label: eq-0377

   {qsat,1\left\lbrack 1 + v(s - 1) \right\rbrack}_{infl,\max}

The liquid water content of the top soil layer relative to effective porosity and adjusted for saturated fraction is determined from



.. math::
   :label: eq-0096

   s = \frac{\frac{\theta_{liq,1}}{\max\left( \theta_{imp},\theta_{sat,1} - \theta_{ice,1} \right)} - f_{sat}}{1 - f_{sat}} \geq 0 \frac{\theta_{liq,1}}{\max\left( \theta_{imp},\theta_{sat,1} - \theta_{ice,1} \right)} \geq 0.01
   1 - f_{sat} \geq 0.01



where :math:`\theta_{liq,1}` and :math:`\theta_{ice,1}` are the volumetric liquid water and ice contents of the top soil layer, and :math:`\theta_{imp} = 0.05` is a minimum effective porosity. The variable :math:`v` is

.. math::
   :label: eq-0378

   v = - \left( \frac{d\psi}{ds} \right)_{s = 1}\frac{1}{0.5\Delta z_{1}}

where :math:`\Delta z_{1}` is the thickness of the top soil layer (mm) and

.. math::
   :label: eq-0379

   \left( \frac{d\psi}{ds} \right)_{s = 1} = - B_{1}\psi_{sat,1}

The saturated hydraulic conductivity :math:`k_{sat,1}` (kg m\ :sup:`-2` s\ :sup:`-1`), volumetric water content at saturation (i.e., porosity) :math:`\theta_{sat,1}`, Clapp and Hornberger (1978) exponent :math:`B_{1}`, and saturated soil matric potential :math:`\psi_{sat,1}` (mm) are determined from soil texture (section 5.3.1).

Infiltration into the surface soil layer of the pervious road is defined as the residual of the surface water balance

.. math::
   :label: eq-0380

   q_{infl} = q_{liq,0} - q_{over} - q_{seva}

when no snow layers exist, and

.. math::
   :label: eq-0381

   q_{infl} = q_{liq,0} - q_{over}

when at least one snow layer is present.

The infiltration for urban surfaces other than pervious road is

.. math::
   :label: eq-0382

   q_{infl} = 0

Soil Water for the Pervious Road
--------------------------------

Soil water for the pervious road is predicted from a multi-layer model, in which the vertical soil moisture transport is governed by infiltration, surface and sub-surface runoff, gradient diffusion, gravity, evapotranspiration through root extraction, and interactions with groundwater (:numref:`fig-hydrologic-processes-pervious-road`). Vegetation is not represented explicitly, however, the total evaporation calculated in section 3.2.4, if not assigned to surface evaporation, is removed from each soil layer through an evapotranspiration loss (:math:`s` in the equation below). The following derivation generally follows that of Z.-L. Yang (1998, unpublished manuscript) with modifications by Zeng and Decker (2009).

For one-dimensional vertical water flow in soils, the conservation of mass is stated as

.. math::
   :label: eq-0383

   \frac{\partial\theta}{\partial t} = - \frac{\partial q}{\partial z} - Q

where :math:`\theta` is the volumetric soil water content (mm\ :sup:`3` of water mm\ :sup:`-3` of soil), :math:`t` is time (s), :math:`z` is height above some datum in the soil column (mm) (positive upwards), :math:`q` is soil water flux (kg m\ :sup:`-2` s\ :sup:`-1` or mm s\ :sup:`-1`) (positive upwards), and :math:`Q` is a soil moisture sink term (mm of water mm\ :sup:`-1` of soil s\ :sup:`-1`) (ET loss). This equation is solved numerically by dividing the soil column into multiple layers in the vertical and integrating downward over each layer with an upper boundary condition of the infiltration flux into the top soil layer :math:`q_{infl}` and a lower boundary condition specified as zero flux.

The soil water flux :math:`q` in equation :eq:`eq-0383` can be described by Darcy's law

.. math::
   :label: eq-0384

   q = - k\frac{\partial\psi_{h}}{\partial z}

where :math:`k` is the hydraulic conductivity (mm s\ :sup:`-1`), and :math:`\psi_{h}` is the hydraulic potential (mm). The hydraulic potential is

.. math::
   :label: eq-0385

   \psi_{h} = \psi_{m} + \psi_{z}

where :math:`\psi_{m}` is the soil matric potential (mm) (which is related to the adsorptive and capillary forces within the soil matrix), and :math:`\psi_{z}` is the gravitational potential (mm) (the vertical distance from an arbitrary reference elevation to a point in the soil). If the reference elevation is the soil surface, then :math:`\psi_{z} = z`. Letting :math:`\psi = \psi_{m}`, Darcy's law becomes

.. math::
   :label: eq-0386

   q = - k\left\lbrack \frac{\partial(\psi + z)}{\partial z} \right\rbrack

Darcy's equation :eq:`eq-0383` can be further manipulated to yield

.. math::
   :label: eq-0387

   q = - k\left\lbrack \frac{\partial(\psi + z)}{\partial z} \right\rbrack = - k\left( \frac{\partial\psi}{\partial z} + 1 \right) = - k\left( \frac{\partial\theta}{\partial z}\frac{\partial\psi}{\partial\theta} + 1 \right)

Substitution of this equation into equation :eq:`eq-0383` with :math:`Q = 0`, yields the Richards equation

.. math::
   :label: eq-0388

   \frac{\partial\theta}{\partial t} = \frac{\partial}{\partial z}\left\lbrack k\left( \frac{\partial\theta}{\partial z}\frac{\partial\psi}{\partial\theta} \right) + 1 \right\rbrack

Zeng and Decker (2009) note that this :math:`\theta` -based form of the Richards equation cannot maintain the hydrostatic equilibrium soil moisture distribution because of the truncation errors of the finite-difference numerical scheme. They show that this deficiency can be overcome by subtracting the equilibrium state from equation :eq:`eq-0386` as

.. math::
   :label: eq-0389

   q = - k\left\lbrack \frac{\partial(\psi + z - C)}{\partial z} \right\rbrack

where :math:`C` is a constant hydraulic potential above the water table :math:`z_{\nabla}`

.. math::
   :label: eq-0390

   C = \psi_{E} + z = \psi_{sat}\left\lbrack \frac{\theta_{E}(z)}{\theta_{sat}} \right\rbrack^{- B} + z = \psi_{sat} + z_{\nabla}

so that

.. math::
   :label: eq-0391

   q = - k\left\lbrack \frac{\partial\left( \psi - \psi_{E} \right)}{\partial z} \right\rbrack

where :math:`\psi_{E}`\ is the equilibrium soil matric potential (mm). Substitution of equations :eq:`eq-0390` and :eq:`eq-0389` into equation :eq:`eq-0388` yields Zeng and Decker's (2009) modified Richards equation

.. math::
   :label: eq-0392

   \frac{\partial\theta}{\partial t} = \frac{\partial}{\partial z}\left\lbrack k\left( \frac{\partial\left( \psi - \psi_{E} \right)}{\partial z} \right) \right\rbrack - Q

where the soil moisture source/sink term :math:`Q` is now included.

Hydraulic Properties
~~~~~~~~~~~~~~~~~~~~

The hydraulic conductivity :math:`k_{i}` (mm s\ :sup:`-1`) and the soil matric potential :math:`\psi_{i}` (mm) for layer :math:`i` vary with volumetric soil water :math:`\theta_{i}` and soil texture (:math:`\% sand_{i}` and :math:`\% clay_{i}`, section 1.2.2) based on the work of Clapp and Hornberger (1978) and Cosby et al. (1984). In CLM4, organic matter modifies soil properties according to Lawrence and Slater (2008). Urban soils are assumed to have no organic matter so the equations :eq:`eq-0390` and :eq:`eq-0389` below are shown in their reduced form.

The hydraulic conductivity is defined at the depth of the interface of two adjacent layers :math:`z_{h,i}` (:numref:`fig-soil-water-flux-scheme`) and is a function of the saturated hydraulic conductivity :math:`k_{sat}\left\lbrack z_{h,i} \right\rbrack`, the total (ice plus liquid) volumetric soil moisture of the two layers :math:`\theta_{i}` and :math:`\theta_{i + 1}` and the impermeable fraction :math:`f_{frz,i}`



.. math::
   :label: eq-0097

   k\left\lbrack z_{h,i} \right\rbrack = \left\{ \begin{aligned}
   & \left( 1 - \frac{f_{frz,i} + f_{frz,i + 1}}{2} \right)k_{sat}\left\lbrack z_{h,i} \right\rbrack\left\lbrack \frac{0.5\left( \theta_{i} + \theta_{i + 1} \right)}{0.5\left( \theta_{sat,i} + \theta_{sat,i + 1} \right)} \right\rbrack^{2B_{i} + 3} 1 \leq i \leq N_{levsoi} - 1 \\
   & \left( 1 - f_{frz,i} \right)k_{sat}\left\lbrack z_{h,i} \right\rbrack\left( \frac{\theta_{i}}{\theta_{sat,i}} \right)^{2B_{i} + 3} i = N_{levsoi}
   \end{aligned} \right\}



where :math:`f_{frz,i}` is defined in equation :eq:`eq-0376`. The saturated hydraulic conductivity :math:`k_{sat}\left\lbrack z_{h,i} \right\rbrack` (mm s\ :sup:`-1`) depends on soil texture (Cosby et al. 1984) as

.. math::
   :label: eq-0393

   k_{sat}\left\lbrack z_{h,i} \right\rbrack = 0.0070556 \times 10^{- 0.884 + 0.0153(\% sand)_{i}}

The water content at saturation (i.e., porosity) is

.. math::
   :label: eq-0394

   \theta_{sat,i} = 0.489 - 0.00126(\% sand)_{i}

and the exponent "\ :math:`B`\ " is

.. math::
   :label: eq-0395

   B_{i} = 2.91 + 0.159(\% clay)_{i}

The soil matric potential (mm) is defined at the node depth :math:`z_{i}` of each layer :math:`i` (:numref:`fig-soil-water-flux-scheme`)

.. math::
   :label: eq-0396

   \psi_{i} = \psi_{sat,i}\left( \frac{\theta_{i}}{\theta_{sat,i}} \right)^{- B_{i}} \geq - 1 \times 10^{8} 0.01 \leq \frac{\theta_{i}}{\theta_{sat,i}} \leq 1

where the saturated soil matric potential (mm) is

.. math::
   :label: eq-0397

   \psi_{sat,i} = - 10.0 \times 10^{1.88 - 0.0131(\% sand)_{i}}

.. _numerical-solution-1:

Numerical Solution
~~~~~~~~~~~~~~~~~~

With reference to :numref:`fig-soil-water-flux-scheme`, the equation for conservation of mass (equation :eq:`eq-0383`) can be integrated over each layer as

.. math::
   :label: eq-0398

   \int_{- z_{h,i}}^{- z_{h,i - 1}}\frac{\partial\theta}{\partial t}dz = - \int_{- z_{h,i}}^{- z_{h,i - 1}}\frac{\partial q}{\partial z}dz - \int_{- z_{h,i}}^{- z_{h,i - 1}}{Qdz}

Note that the integration limits are negative since :math:`z` is defined as positive upward from the soil surface. This equation :eq:`eq-0383` can be written as

.. math::
   :label: eq-0399

   \Delta z_{i}\frac{\partial\theta_{liq,i}}{\partial t} = - q_{i - 1} + q_{i} - e_{i}

where :math:`q_{i}` is the flux of water across interface :math:`z_{h,i}`, :math:`q_{i - 1}` is the flux of water across interface :math:`z_{h,i - 1}`, and :math:`e_{i}` is a layer-averaged soil moisture sink term (ET loss) defined as positive for flow out of the layer (mm s\ :sup:`-1`). Taking the finite difference with time and evaluating the fluxes implicitly at time :math:`n + 1` yields

.. math::
   :label: eq-0400

   \frac{\Delta z_{i}\Delta\theta_{liq,i}}{\Delta t} = - q_{i - 1}^{n + 1} + q_{i}^{n + 1} - e_{i}

where :math:`\Delta\theta_{liq,i} = \theta_{liq,i}^{n + 1} - \theta_{liq,i}^{n}` is the change in volumetric soil liquid water of layer :math:`i` in time :math:`\Delta t`\ and :math:`\Delta z_{i}` is the thickness of layer :math:`i` (mm).

The water removed by evapotranspiration in each layer :math:`e_{i}` is a function of the total evapotranspiration :math:`E_{prvrd}^{et}` (section 3.2.4) and the effective root fraction :math:`r_{e,i}`

.. math::
   :label: eq-0401

   e_{i} = r_{e,i}E_{prvrd}^{et}

The effective root fraction :math:`r_{e,i}` is



.. math::
   :label: eq-0098

   r_{e,i} = \left\{ \begin{aligned}
   & \frac{r_{i}w_{i}}{\alpha_{soi}} \alpha_{soi} > 0 \\
   & 0 \alpha_{soi} = 0
   \end{aligned} \right\}



where :math:`r_{i}` is the fraction of roots in layer :math:`i` (equation :eq:`eq-0245`), :math:`w_{i}` is a soil wetness factor for layer :math:`i` (equation (3.87)), and :math:`\alpha_{soi}` is a wetness factor for the total soil column (equation :eq:`eq-0244` (section 3.2.3)).

.. figure:: image18.jpeg
   :width: 6in
   :height: 5.17708in
   :name: fig-soil-water-flux-scheme

   Schematic diagram of numerical scheme used to solve for soil water fluxes. Shown are three soil layers, :math:`i - 1`, :math:`i`, and :math:`i + 1`. The soil matric potential :math:`\psi` and volumetric soil water :math:`\theta_{liq}` are defined at the layer node depth :math:`z`. The hydraulic conductivity :math:`k\left\lbrack z_{h} \right\rbrack` is defined at the interface of two layers :math:`z_{h}`. The layer thickness is :math:`\Delta z`. The soil water fluxes :math:`q_{i - 1}` and :math:`q_{i}` are defined as positive upwards. The soil moisture sink term :math:`e` (ET loss) is defined as positive for flow out of the layer.

The soil water fluxes in equation :eq:`eq-0400` , which are a function of :math:`\theta_{liq,i}` and :math:`\theta_{liq,i + 1}` because of their dependence on hydraulic conductivity and soil matric potential, can be linearized about :math:`\partial\theta` using a Taylor series expansion as

.. math::
   :label: eq-0402

   q_{i}^{n + 1} = q_{i}^{n} + \frac{\partial q_{i}}{\partial\theta_{liq,i}}\Delta\theta_{liq,i} + \frac{\partial q_{i}}{\partial\theta_{liq,i + 1}}\Delta\theta_{liq,i + 1}

.. math::
   :label: eq-0403

   q_{i - 1}^{n + 1} = q_{i - 1}^{n} + \frac{\partial q_{i - 1}}{\partial\theta_{liq,i - 1}}\Delta\theta_{liq,i - 1} + \frac{\partial q_{i - 1}}{\partial\theta_{liq,i}}\Delta\theta_{liq,i}

Substitution of these expressions for :math:`q_{i}^{n + 1}` and :math:`q_{i - 1}^{n + 1}` into equation :eq:`eq-0400` results in a general tridiagonal equation set of the form

.. math::
   :label: eq-0404

   r_{i} = a_{i}\Delta\theta_{liq,i - 1} + b_{i}\Delta\theta_{liq,i} + c_{i}\Delta\theta_{liq,i + 1}

where

.. math::
   :label: eq-0405

   a_{i} = - \frac{\partial q_{i - 1}}{\partial\theta_{liq,i - 1}}

.. math::
   :label: eq-0406

   b_{i} = \frac{\partial q_{i}}{\partial\theta_{liq,i}} - \frac{\partial q_{i - 1}}{\partial\theta_{liq,i}} - \frac{\Delta z_{i}}{\Delta t}

.. math::
   :label: eq-0407

   c_{i} = \frac{\partial q_{i}}{\partial\theta_{liq,i + 1}}

.. math::
   :label: eq-0408

   r_{i} = q_{i - 1}^{n} - q_{i}^{n} + e_{i}

The tridiagonal equation set is solved over :math:`i = 1,\ldots,N_{levsoi} + 1` where the layer :math:`i = N_{levsoi} + 1` is a virtual layer representing the aquifer.

The finite-difference forms of the fluxes and partial derivatives in equations :eq:`eq-0451`-:eq:`eq-0447` can be obtained from equation :eq:`eq-0391` as

.. math::
   :label: eq-0409

   q_{i - 1}^{n} = - k\left\lbrack z_{h,i - 1} \right\rbrack\left\lbrack \frac{\left( \psi_{i - 1} - \psi_{i} \right) + \left( \psi_{E,i} - \psi_{E,i - 1} \right)}{z_{i} - z_{i - 1}} \right\rbrack

.. math::
   :label: eq-0410

   q_{i}^{n} = - k\left\lbrack z_{h,i} \right\rbrack\left\lbrack \frac{\left( \psi_{i} - \psi_{i + 1} \right) + \left( \psi_{E,i + 1} - \psi_{E,i} \right)}{z_{i + 1} - z_{i}} \right\rbrack

.. math::
   :label: eq-0411

   \frac{\partial q_{i - 1}}{\partial\theta_{liq,i - 1}} = - \left\lbrack \frac{k\left\lbrack z_{h,i - 1} \right\rbrack}{z_{i} - z_{i - 1}}\frac{\partial\psi_{i - 1}}{\partial\theta_{liq,i - 1}} \right\rbrack - \frac{\partial k\left\lbrack z_{h,i - 1} \right\rbrack}{\partial\theta_{liq,i - 1}}\left\lbrack \frac{\left( \psi_{i - 1} - \psi_{i} \right) + \left( \psi_{E,i} - \psi_{E,i - 1} \right)}{z_{i} - z_{i - 1}} \right\rbrack

.. math::
   :label: eq-0412

   \frac{\partial q_{i - 1}}{\partial\theta_{liq,i}} = \left\lbrack \frac{k\left\lbrack z_{h,i - 1} \right\rbrack}{z_{i} - z_{i - 1}}\frac{\partial\psi_{i}}{\partial\theta_{liq,i}} \right\rbrack - \frac{\partial k\left\lbrack z_{h,i - 1} \right\rbrack}{\partial\theta_{liq,i}}\left\lbrack \frac{\left( \psi_{i - 1} - \psi_{i} \right) + \left( \psi_{E,i} - \psi_{E,i - 1} \right)}{z_{i} - z_{i - 1}} \right\rbrack

.. math::
   :label: eq-0413

   \frac{\partial q_{i}}{\partial\theta_{liq,i}} = - \left\lbrack \frac{k\left\lbrack z_{h,i} \right\rbrack}{z_{i + 1} - z_{i}}\frac{\partial\psi_{i}}{\partial\theta_{liq,i}} \right\rbrack - \frac{\partial k\left\lbrack z_{h,i} \right\rbrack}{\partial\theta_{liq,i}}\left\lbrack \frac{\left( \psi_{i} - \psi_{i + 1} \right) + \left( \psi_{E,i + 1} - \psi_{E,i} \right)}{z_{i + 1} - z_{i}} \right\rbrack

.. math::
   :label: eq-0414

   \frac{\partial q_{i}}{\partial\theta_{liq,i + 1}} = \left\lbrack \frac{k\left\lbrack z_{h,i} \right\rbrack}{z_{i + 1} - z_{i}}\frac{\partial\psi_{i + 1}}{\partial\theta_{liq,i + 1}} \right\rbrack - \frac{\partial k\left\lbrack z_{h,i} \right\rbrack}{\partial\theta_{liq,i + 1}}\left\lbrack \frac{\left( \psi_{i} - \psi_{i + 1} \right) + \left( \psi_{E,i + 1} - \psi_{E,i} \right)}{z_{i + 1} - z_{i}} \right\rbrack

The derivatives of the soil matric potential at the node depth are derived from equation :eq:`eq-0396`

.. math::
   :label: eq-0415

   \frac{\partial\psi_{i - 1}}{\partial\theta_{liq,i - 1}} = - B_{i - 1}\frac{\psi_{i - 1}}{\theta_{i - 1}}

.. math::
   :label: eq-0416

   \frac{\partial\psi_{i}}{\partial\theta_{liq,i}} = - B_{i}\frac{\psi_{i}}{\theta_{i}}

.. math::
   :label: eq-0417

   \frac{\partial\psi_{i + 1}}{\partial\theta_{liq,i + 1}} = - B_{i + 1}\frac{\psi_{i + 1}}{\theta_{i + 1}}

with the constraint :math:`0.01\theta_{sat,i} \leq \theta_{i} \leq \theta_{sat,i}`.

The derivatives of the hydraulic conductivity at the layer interface are derived from equation :eq:`eq-0097`



.. math::
   :label: eq-0099

   \frac{\partial k\left\lbrack z_{h,i - 1} \right\rbrack}{\partial\theta_{liq,i - 1}} = \frac{\partial k\left\lbrack z_{h,i - 1} \right\rbrack}{\partial\theta_{liq,i}} = \left( 1 - \frac{f_{frz,i - 1} + f_{frz,i}}{2} \right)\left( 2B_{i - 1} + 3 \right)k_{sat}\left\lbrack z_{h,i - 1} \right\rbrack \times
   \left\lbrack \frac{0.5\left( \theta_{i - 1} + \theta_{i} \right)}{0.5\left( \theta_{sat,i - 1} + \theta_{sat,i} \right)} \right\rbrack^{2B_{i - 1} + 2}\left( \frac{0.5}{\theta_{sat,i - 1}} \right)





.. math::
   :label: eq-0100

   \frac{\partial k\left\lbrack z_{h,i} \right\rbrack}{\partial\theta_{liq,i}} = \frac{\partial k\left\lbrack z_{h,i} \right\rbrack}{\partial\theta_{liq,i + 1}} = \left( 1 - \frac{f_{frz,i} + f_{frz,i + 1}}{2} \right)\left( 2B_{i} + 3 \right)k_{sat}\left\lbrack z_{h,i} \right\rbrack \times
   \left\lbrack \frac{0.5\left( \theta_{i} + \theta_{i + 1} \right)}{0.5\left( \theta_{sat,i} + \theta_{sat,i + 1} \right)} \right\rbrack^{2B_{i} + 2}\left( \frac{0.5}{\theta_{sat,i}} \right)


Equilibrium soil matric potential and volumetric moisture
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The equilibrium soil matric potential :math:`\psi_{E}` can be derived from equation :eq:`eq-0390` as

.. math::
   :label: eq-0418

   \psi_{E} = \psi_{sat}\left( \frac{\theta_{E}(z)}{\theta_{sat}} \right)^{- B}

and the equilibrium volumetric water content :math:`\theta_{E}(z)` at depth :math:`z` can also be derived as

.. math::
   :label: eq-0419

   \theta_{E}(z) = \theta_{sat}\left( \frac{\psi_{sat} + z_{\nabla} - z}{\psi_{sat}} \right)^{- \frac{1}{B}}

Here, the soil matric potentials, the water table depth :math:`z_{\nabla}` and the soil depths have units of mm. For the finite-difference scheme, a layer-average equilibrium volumetric water content is used in equation :eq:`eq-0418` and can be obtained from

.. math::
   :label: eq-0420

   \overline{\theta_{E,i}} = \int_{z_{h,i - 1}}^{z_{h,i}}\frac{\theta_{E}(z)}{z_{h,i} - z_{h,i - 1}}dz

which when integrated yields

.. math::
   :label: eq-0421

   \overline{\theta_{E,i}} = \frac{\theta_{sat, i}\psi_{sat,i}}{\left( z_{h,i} - z_{h,i - 1} \right)\left( 1 - \frac{1}{B_{i}} \right)}\left\lbrack \left( \frac{\psi_{sat,i} - z_{\nabla} + z_{h,i}}{\psi_{sat,i}} \right)^{1 - \frac{1}{B_{i}}} - \left( \frac{\psi_{sat,i} - z_{\nabla} + z_{h,i - 1}}{\psi_{sat,i}} \right)^{1 - \frac{1}{B_{i}}} \right\rbrack

equation :eq:`eq-0421` is valid when the water table :math:`z_{\nabla}` is deeper than both interface depths :math:`z_{h,i - 1}` and :math:`z_{h,i}`. Since the water table can be within the soil column, the equation is modified if the water table is within soil layer :math:`i` (:math:`z_{h,i - 1} < z_{\nabla} < z_{h,i}`) as a weighted average of the saturated part and the unsaturated part

.. math::
   :label: eq-0422

   \overline{\theta_{E,i}} = \overline{\theta_{E,sat,i}}\left( \frac{z_{h,i} - z_{\nabla}}{z_{h,i} - z_{h,i - 1}} \right) + \overline{\theta_{E,unsat,i}}\left( \frac{z_{\nabla} - z_{h,i - 1}}{z_{h,i} - z_{h,i - 1}} \right)

where :math:`\overline{\theta_{E,sat,i}} = \theta_{sat,i}` and the unsaturated part :math:`\overline{\theta_{E,unsat,i}}` is

.. math::
   :label: eq-0423

   \overline{\theta_{E,unsat,i}} = \frac{\theta_{sat, i}\psi_{sat,i}}{\left( z_{\nabla} - z_{h,i - 1} \right)\left( 1 - \frac{1}{B_{i}} \right)}\left\lbrack 1 - \left( \frac{\psi_{sat,i} - z_{\nabla} + z_{h,i - 1}}{\psi_{sat,i}} \right)^{1 - \frac{1}{B_{i}}} \right\rbrack

If :math:`z_{\nabla} < z_{h,i - 1}`, then :math:`\overline{\theta_{E,i}} = \overline{\theta_{E,sat,i}} = \theta_{sat,i}`. If the water table is below the soil column (:math:`z_{\nabla} > z_{h,N_{levsoi}}`), an equilibrium volumetric soil moisture is calculated for a virtual layer :math:`i = N_{levsoi} + 1` as

.. math::
   :label: eq-0424

   \overline{\theta_{E,i = N_{levsoi + 1}}} = \frac{\theta_{sat,i - 1}\psi_{sat,i - 1}}{\left( z_{\nabla} - z_{h,i - 1} \right)\left( 1 - \frac{1}{B_{i - 1}} \right)}\left\lbrack 1 - \left( \frac{\psi_{sat,i - 1} - z_{\nabla} + z_{h,i - 1}}{\psi_{sat,i - 1}} \right)^{1 - \frac{1}{B_{i - 1}}} \right\rbrack

The equilibrium volumetric soil moisture is constrained by

.. math::
   :label: eq-0425

   0 \leq \overline{\theta_{E,i}} \leq \theta_{sat,i}

The equilibrium soil matric potential is then

.. math::
   :label: eq-0426

   \psi_{E,i} = \psi_{sat,i}\left( \frac{\overline{\theta_{E,i}}}{\theta_{sat,i}} \right)^{- B_{i}} \geq - 1 \times 10^{8} \frac{\overline{\theta_{E,i}}}{\theta_{sat,i}} \geq 0.01

Equation set for layer :math:`\mathbf{i = 1}`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the top soil layer (:math:`i = 1`), the boundary condition is the infiltration rate (section 5.2), :math:`q_{i - 1}^{n + 1} = - q_{infl}^{n + 1}`, and the water balance equation is

.. math::
   :label: eq-0427

   \frac{\Delta z_{i}\Delta\theta_{liq,i}}{\Delta t} = q_{infl}^{n + 1} + q_{i}^{n + 1} - e_{i}

After grouping like terms, the coefficients of the tridiagonal set of equations for :math:`i = 1` are

.. math::
   :label: eq-0428

   a_{i} = 0

.. math::
   :label: eq-0429

   b_{i} = \frac{\partial q_{i}}{\partial\theta_{liq,i}} - \frac{\Delta z_{i}}{\Delta t}

.. math::
   :label: eq-0430

   c_{i} = \frac{\partial q_{i}}{\partial\theta_{liq,i + 1}}

.. math::
   :label: eq-0431

   r_{i} = q_{infl}^{n + 1} - q_{i}^{n} + e_{i}

Equation set for layers :math:`\mathbf{i = 2,}\mathbf{\ldots}\mathbf{,}\mathbf{N}_{\mathbf{levsoi}}\mathbf{-}\mathbf{1}`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The coefficients of the tridiagonal set of equations for :math:`i = 2,\ldots,N_{levsoi} - 1` are

.. math::
   :label: eq-0432

   a_{i} = - \frac{\partial q_{i - 1}}{\partial\theta_{liq,i - 1}}

.. math::
   :label: eq-0433

   b_{i} = \frac{\partial q_{i}}{\partial\theta_{liq,i}} - \frac{\partial q_{i - 1}}{\partial\theta_{liq,i}} - \frac{\Delta z_{i}}{\Delta t}

.. math::
   :label: eq-0434

   c_{i} = \frac{\partial q_{i}}{\partial\theta_{liq,i + 1}}

.. math::
   :label: eq-0435

   r_{i} = q_{i - 1}^{n} - q_{i}^{n} + e_{i}

Equation set for layers :math:`\mathbf{i =}\mathbf{N}_{\mathbf{levsoi}}\mathbf{,}\mathbf{\ldots}\mathbf{N}_{\mathbf{levsoi}}\mathbf{+ 1}`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the lowest soil layer (:math:`i = N_{levsoi}`), the bottom boundary condition depends on the depth of the water table. If the water table is within the soil column (:math:`z_{\nabla} \leq z_{h,N_{levsoi}}`), a zero-flux bottom boundary condition is applied (:math:`q_{i}^{n} = 0`) and the coefficients of the tridiagonal set of equations for :math:`i = N_{levsoi}` are

.. math::
   :label: eq-0436

   a_{i} = - \frac{\partial q_{i - 1}}{\partial\theta_{liq,i - 1}}

.. math::
   :label: eq-0437

   b_{i} = - \frac{\partial q_{i - 1}}{\partial\theta_{liq,i}} - \frac{\Delta z_{i}}{\Delta t}

.. math::
   :label: eq-0438

   c_{i} = 0

.. math::
   :label: eq-0439

   r_{i} = q_{i - 1}^{n} + e_{i}

The coefficients for the aquifer layer :math:`i = N_{levsoi} + 1` are then

.. math::
   :label: eq-0440

   a_{i} = 0

.. math::
   :label: eq-0441

   b_{i} = - \frac{\Delta z_{i}}{\Delta t}

.. math::
   :label: eq-0442

   c_{i} = 0

.. math::
   :label: eq-0443

   r_{i} = 0

If the water table is below the soil column (:math:`z_{\nabla} > z_{h,N_{levsoi}}`), the coefficients for :math:`i = N_{levsoi}` are

.. math::
   :label: eq-0444

   a_{i} = - \frac{\partial q_{i - 1}}{\partial\theta_{liq,i - 1}}

.. math::
   :label: eq-0445

   b_{i} = \frac{\partial q_{i}}{\partial\theta_{liq,i}} - \frac{\partial q_{i - 1}}{\partial\theta_{liq,i}} - \frac{\Delta z_{i}}{\Delta t}

.. math::
   :label: eq-0446

   c_{i} = \frac{\partial q_{i}}{\partial\theta_{liq,i + 1}}

.. math::
   :label: eq-0447

   r_{i} = q_{i - 1}^{n} - q_{i}^{n} + e_{i}

The :math:`i = N_{levsoi} + 1` terms are evaluated using

.. math::
   :label: eq-0448

   \psi_{N_{levsoi} + 1} = \psi_{sat,N_{levsoi}}\left\lbrack s_{N_{levsoi} + 1} \right\rbrack^{- B_{N_{levsoi}}} \geq - 1 \times 10^{8}

.. math::
   :label: eq-0449

   z_{N_{levsoi} + 1} = 0.5\left( z_{\nabla} + z_{N_{levsoi}} \right)

where

:math:`s_{N_{levsoi} + 1} = 0.5\left( \frac{\theta_{sat,N_{levsoi}} + \theta_{N_{levsoi}}}{\theta_{sat,N_{levsoi}}} \right) 0.01 \leq s_{N_{levsoi} + 1} \leq 1`,

:math:`\psi_{E,N_{levsoi} + 1}` is evaluated from equations :eq:`eq-0424` and :eq:`eq-0426`, and

.. math::
   :label: eq-0450

   \frac{\partial\psi_{N_{levsoi} + 1}}{\partial\theta_{liq,N_{levsoi} + 1}} = - B_{N_{levsoi}}\frac{\psi_{N_{levsoi} + 1}}{s_{N_{levsoi}}\theta_{sat,N_{levsoi}}}

The coefficients for the aquifer layer :math:`i = N_{levsoi} + 1` are then

.. math::
   :label: eq-0451

   a_{i} = - \frac{\partial q_{i - 1}}{\partial\theta_{liq,i - 1}}

.. math::
   :label: eq-0452

   b_{i} = - \frac{\partial q_{i - 1}}{\partial\theta_{liq,i}} - \frac{\Delta z_{i}}{\Delta t}

.. math::
   :label: eq-0453

   c_{i} = 0

.. math::
   :label: eq-0454

   r_{i} = q_{i - 1}^{n}

Upon solution of the tridiagonal equations :eq:`eq-0424` and :eq:`eq-0426` set (Press et al. 1992), the liquid water contents are updated as follows

.. math::
   :label: eq-0455

   w_{liq,i}^{n + 1} = w_{liq,i}^{n} + \Delta\theta_{liq,i}\Delta z_{i} i = 1,\ldots,N_{levsoi}

The volumetric water content is

.. math::
   :label: eq-0456

   \theta_{i} = \frac{w_{liq,i}}{\Delta z_{i}\rho_{liq}} + \frac{w_{ice,i}}{\Delta z_{i}\rho_{ice}}

Groundwater-Soil Water Interactions for the Pervious Road
---------------------------------------------------------

Drainage or sub-surface runoff for the pervious road is based on the SIMTOP scheme (Niu et al. 2005) with a modification to account for reduced drainage in frozen soils. In the work of Niu et al. (2005), the drainage :math:`q_{drai}` (kg m\ :sup:`-2` s\ :sup:`-1`) was formulated as

| :math:`q_{drai} = {q\exp\left( - f_{drai}z_{\nabla} \right)}_{drai,\max}`.
| Here, the water table depth :math:`z_{\nabla}` has units of meters. To restrict drainage in frozen soils, Niu et al. (2005) added the following condition

.. math::
   :label: eq-0457

   q_{drai} = 0 \text{for  }w_{ice,N_{levsoi}} > w_{liq,N_{levsoi}}

In preliminary testing it was found that a more gradual restriction of drainage was required so that the water table depth remained dynamic under partially frozen conditions. The following modification is made to equation :eq:`eq-0458`

.. math::
   :label: eq-0458

   q_{drai} = \left( 1 - f_{imp} \right){q\exp\left( - f_{drai}z_{\nabla} \right)}_{drai,\max}

where :math:`f_{imp}` is the fraction of impermeable area determined from the ice content of the soil layers interacting with the water table

.. math::
   :label: eq-0459

   f_{imp} = \frac{\exp\left\lbrack - \alpha\left( 1 - \frac{\sum_{i = jwt}^{i = N_{levsoi}}{\frac{w_{ice,i}}{w_{ice,i} + w_{liq,i}}\Delta z_{i}}}{\sum_{i = jwt}^{i = N_{levsoi}}{\Delta z_{i}}} \right) \right\rbrack - \exp( - \alpha)}{1 - \exp( - \alpha)} \geq 0

where :math:`\alpha = 3` is an adjustable scale-dependent parameter, :math:`jwt` is the index of the layer directly above the water table, :math:`w_{ice,i}` and :math:`w_{liq,i}` are the ice and liquid water contents of soil layer :math:`i` (kg m\ :sup:`-2`), and :math:`\Delta z_{i}` is the layer thickness (m). This expression is functionally the same as that used to determine the impermeable fraction (equation :eq:`eq-0376`). In equation (5.140), the decay factor :math:`f_{drai} = 2.5` m\ :sup:`-1` and the maximum drainage when the water table depth is at the surface :math:`{q - 3}_{drai,\max}` kg m\ :sup:`-2` s\ :sup:`-1` were determined for global simulations through sensitivity analysis and comparison with observed runoff.

Determination of water table depth :math:`z_{\nabla}` is based on work by Niu et al. (2007). In this approach, a groundwater component is added in the form of an unconfined aquifer lying below the soil column (:numref:`fig-hydrologic-processes-pervious-road`). The groundwater solution is dependent on whether the water table is within or below the soil column. Two water stores are used to account for these solutions. The first, :math:`W_{a}`, is the water stored in the unconfined aquifer (mm) and is proportional to the change in water table depth when the water table is below the lower boundary of the hydrologically-active soil column. The second, :math:`W_{t}`, is the actual groundwater which can include water within the soil column. When the water table is below the soil column :math:`W_{t} = W_{a}`. When the water table is within the soil column, :math:`W_{a}` is constant because there is no water exchange between the soil column and the underlying aquifer, while :math:`W_{t}` varies with soil moisture conditions.

In either case, :math:`W_{t}` is first updated as

.. math::
   :label: eq-0460

   W_{t}^{n + 1} = W_{t}^{n} + \left( q_{recharge} - q_{drai} \right)\Delta t

where :math:`\Delta t` is the model time step (s), :math:`q_{recharge}` is the recharge to the aquifer (kg m\ :sup:`-2` s\ :sup:`-1`), and the drainage :math:`q_{drai}` calculated from equation (5.140) is equivalent to the groundwater discharge.

For the case when the water table is below the soil column, the water stored in the unconfined aquifer :math:`W_{a}` (mm) is updated as

.. math::
   :label: eq-0461

   W_{a}^{n + 1} = W_{a}^{n} + \left( q_{recharge} - q_{drai} \right)\Delta t

and :math:`W_{t}^{n + 1}` is reset as :math:`W_{t}^{n + 1} = W_{a}^{n + 1}`. The recharge rate is defined as positive when water enters the aquifer

.. math::
   :label: eq-0462

   q_{recharge} = \frac{\Delta\theta_{liq,N_{levsoi} + 1}\Delta z_{N_{levsoi} + 1}}{\Delta t}

where :math:`\Delta\theta_{liq,N_{levsoi} + 1} = \theta_{liq,N_{levsoi} + 1}^{n + 1} - \theta_{liq,N_{levsoi} + 1}^{n}` is the change in liquid water content for layer :math:`i = N_{levsoi} + 1` calculated from the solution of the soil water equation (5.140) (section 5.3), and :math:`\Delta z_{N_{levsoi} + 1}` (mm) is

.. math::
   :label: eq-0463

   \Delta z_{N_{levsoi} + 1} = z_{\nabla}^{n} - z_{h,N_{levsoi}}

The water table depth is calculated from the aquifer water storage scaled by the average specific yield :math:`S_{y} = 0.2` [the fraction of water volume that can be drained by gravity in an unconfined aquifer (Niu et al. 2007)]

.. math::
   :label: eq-0464

   z_{\nabla} = z_{h,N_{levsoi}} + 25 - \frac{W_{a}}{10^{3}S_{y}}

The form of equation :eq:`eq-0464` originates from the assumption that the initial amount of water in the aquifer is 4800 mm and the corresponding water table depth is one meter below the bottom of the soil column. The water table depth is at the bottom of the soil column (:math:`z_{\nabla} = z_{h,N_{levsoi}}`) when the aquifer water is at its prescribed maximum value (5000 mm). The bottom soil layer liquid water content is updated for excess aquifer water as

.. math::
   :label: eq-0465

   w_{liq,N_{levsoi}}^{n + 1} = w_{liq,N_{levsoi}}^{n} + \max\left( 0,W_{a} - 5000 \right)

and aquifer water is reset to :math:`W_{a} \leq 5000`.

For the case when the water table is within the soil column, there is no water exchange between the soil column and the underlying aquifer. However, variations of the water table depth are still calculated as



.. math::
   :label: eq-0101

   z_{\nabla} = \left\{ \begin{aligned}
   & z_{h,jwt + 1} - \left\lbrack \frac{W_{t} - 10^{3} \times 25S_{y} - \sum_{i = jwt + 2}^{N_{levsoi}}{\Delta z_{i}\left( \theta_{sat,i} - \theta_{ice,i} \right)}}{10^{3}\left( \theta_{sat,jwt + 1} - \theta_{ice,jwt + 1} \right)} \right\rbrack jwt = 1,\ldots N_{levsoi} - 2 \\
   & z_{h,jwt + 1} - \left\lbrack \frac{W_{t} - 10^{3} \times 25S_{y}}{10^{3}\left( \theta_{sat,jwt + 1} - \theta_{ice,jwt + 1} \right)} \right\rbrack jwt = N_{levsoi} - 1
   \end{aligned} \right\}



where :math:`jwt` is the index of the layer directly above the water table, and limits are placed on the water table depth as :math:`0.05 \leq z_{\nabla} \leq 80`. In the work of Niu et al. (2007), the water table depth in this case was calculated with the specific yield determined by the volume of air pores (the pore space not filled with water) within the soil to convert :math:`W_{t}` to a water table depth. However, this was found to result in unstable water table calculations for a significant proportion of grid cells in global simulations. More specifically, when repeatedly forcing the model with a single year of atmospheric data, the temporal evolution of water table depth was significantly different from year to year for some grid cells, with occasional rapid (within a few days) movement of the water table to the soil surface in some cases. This occurred in grid cells with soil water contents near saturation because of the small amount of available pore space. This had deleterious implications for stability of surface fluxes and temperature. In equation :eq:`eq-0101` , the calculation is based on effective porosity (:math:`\theta_{sat,i} - \theta_{ice,i} \geq 0.01`) only. Although less defensible from a physical viewpoint, the approach stabilizes the water table calculation for these grid cells and eliminates unrealistic oscillations in surface fluxes and temperature.

In this case, the drainage :math:`q_{drai}` is extracted from the soil liquid water in layers within the water table. The partitioning of drainage from these layers is proportional to the layer thickness-weighted hydraulic conductivity as

.. math::
   :label: eq-0466

   w_{liq,i}^{n + 1} = w_{liq,i}^{n} - \frac{q_{drai}k\left\lbrack z_{h,i} \right\rbrack\Delta t\Delta z_{i}}{\sum_{i = jwt + 1}^{i = N_{levsoi}}{k\left\lbrack z_{h,i} \right\rbrack\Delta z_{i}}} i = jwt + 1,\ldots,N_{levsoi}

where :math:`\Delta t` is the time step (s).

After the above calculations, two numerical adjustments are implemented to keep the liquid water content of each soil layer (:math:`w_{liq,i}`) within physical constraints of :math:`w_{liq}^{min_{liq,i}\left( \theta_{sat,i} - \theta_{ice,i} \right)_{i}}` where :math:`w_{liq}^{\min}` (mm). First, beginning with the bottom soil layer :math:`i = N_{levsoi}`, any excess liquid water in each soil layer (:math:`w_{liq,i}^{excess} = w_{liq,i} - \left( \theta_{sat,i} - \theta_{ice,i} \right)\Delta z_{i} \geq 0`) is successively added to the layer above. Any excess liquid water that remains after saturating the entire soil column (plus a maximum surface ponding depth :math:`w_{liq}^{pond} = 10` kg m\ :sup:`-2` s\ :sup:`-1`), is added to drainage :math:`q_{drai}`. Second, to prevent negative :math:`w_{liq,i}`, each layer is successively brought up to :math:`w_{liq,i} = w_{liq,\min}` by taking the required amount of water from the layer below. If this results in :math:`w_{liq,N_{levsoi}} < w_{liq}^{\min}`, then the layers above are searched in succession for the required amount of water (:math:`w_{liq}^{min_{liq,N_{levsoi}}}`) and removed from those layers subject to the constraint :math:`w_{liq,i} \geq w_{liq}^{\min}`. If sufficient water is not found, then the water is removed from :math:`W_{t}` and :math:`q_{drai}`.

The surface layer liquid water and ice contents for roof, pervious and impervious road are then updated for dew :math:`q_{sdew}`, frost :math:`q_{frost}`, or sublimation :math:`q_{subl}` (section 3.4) as

.. math::
   :label: eq-0467

   w_{liq,1}^{n + 1} = w_{liq,1}^{n} + q_{sdew}\Delta t

.. math::
   :label: eq-0468

   w_{ice,1}^{n + 1} = w_{ice,1}^{n} + q_{frost}\Delta t

.. math::
   :label: eq-0469

   w_{ice,1}^{n + 1} = w_{ice,1}^{n} - q_{subl}\Delta t

Sublimation of ice is limited to the amount of ice available.

Runoff from snow-capping
------------------------

As with other surfaces, urban surfaces are constrained to have a snow water equivalent :math:`W_{sno} \leq 1000` kg m\ :sup:`-2`. For snow-capped surfaces, the solid and liquid precipitation reaching the snow surface and dew in solid or liquid form, is separated into solid :math:`q_{snwcp,ice}`\ and liquid :math:`q_{snwcp,liq}` runoff terms

.. math::
   :label: eq-0470

   q_{snwcp,ice} = q_{grnd,ice} + q_{frost}

.. math::
   :label: eq-0471

   q_{snwcp,liq} = q_{grnd,liq} + q_{dew}

and snow pack properties are unchanged. The :math:`q_{snwcp,ice}` runoff is sent to the River Transport Model (RTM) where it is routed to the ocean as an ice stream and, if applicable, the ice is melted there. The :math:`q_{snwcp,liq}` runoff is assigned to the runoff term :math:`q_{rgwl}` (e.g. :math:`q_{rgwl} = q_{snwcp,liq}`) and included in the liquid water runoff sent to RTM.

Offline Mode
============

In offline mode (uncoupled to an atmospheric model), the atmospheric forcing required by CLM (:numref:`table-atm-input`) is supplied by observed datasets. The standard forcing provided with the model is a 57-year (1948-2004) dataset that is described in Qian et al. (2006) though alternative observed forcing datasets could also be used. The forcing data is ingested into a data atmosphere model in three "streams"; precipitation (:math:`P`) (mm s\ :sup:`-1`), solar radiation (:math:`S_{atm}`) (W m\ :sup:`-2`), and four other fields [atmospheric pressure :math:`P_{atm}` (Pa), atmospheric specific humidity :math:`q_{atm}` (kg kg\ :sup:`-1`), atmospheric temperature :math:`T_{atm}` (K), and atmospheric wind :math:`W_{atm}` (m s\ :sup:`-1`)]. These are separate streams because they are handled differently according to the type of field and the temporal resolution at which they are provided. In the Qian et al. (2006) dataset, the precipitation stream is provided at six hour intervals and the data atmosphere model prescribes the same precipitation rate for each model time step within the six hour period. The four fields that are grouped together in another stream (pressure, humidity, temperature, and wind) are provided at three hour intervals and the data atmosphere model linearly interpolates these fields to the time step of the model.

The total solar radiation is provided at six hour intervals. The data is fit to the model time step using a diurnal function that depends on the cosine of the solar zenith angle :math:`\mu` to provide a smoother diurnal cycle of solar radiation and to ensure that all of the solar radiation supplied by the six-hourly forcing data is actually used. The solar radiation at model time step :math:`t_{M}` is



.. math::
   :label: eq-0102

   S_{atm}\left( t_{M} \right) = \frac{\frac{\Delta t_{FD}}{\Delta t_{M}}S_{atm}\left( t_{FD} \right)\mu\left( t_{M} \right)}{\sum_{i = 1}^{\frac{\Delta t_{FD}}{\Delta t_{M}}}{\mu\left( t_{M_{i}} \right)}} \text{for }\mu\left( t_{M} \right) > 0.001
   S_{atm}\left( t_{M} \right) = 0 \text{for }\mu\left( t_{M} \right) \leq 0.001



where :math:`\Delta t_{FD}` is the time step of the forcing data (6 hours :math:`\times` 3600 seconds hour\ :sup:`-1` = 21600 seconds), :math:`\Delta t_{M}` is the model time step (seconds), :math:`S_{atm}\left( t_{FD} \right)` is the six-hourly solar radiation from the forcing data (W m\ :sup:`-2`), and :math:`\mu\left( t_{M} \right)` is the cosine of the solar zenith angle at model time step :math:`t_{M}` (section 2.8). The term in the denominator of equation :eq:`eq-0102` is the sum of the cosine of the solar zenith angle for each model time step falling within the six hour period. For numerical purposes, :math:`\mu\left( t_{M_{i}} \right) \geq 0.001`.

The total incident solar radiation :math:`S_{atm}` at the model time step :math:`t_{M}` is then split into near-infrared and visible radiation and partitioned into direct and diffuse according to factors derived from one year's worth of hourly CAM output from CAM version cam3_5_55 as

.. math::
   :label: eq-0472

   S_{atm} \downarrow_{vis}^{\mu} = R_{vis}\left( \alpha S_{atm} \right)

.. math::
   :label: eq-0473

   S_{atm} \downarrow_{nir}^{\mu} = R_{nir}\left\lbrack (1 - \alpha)S_{atm} \right\rbrack

.. math::
   :label: eq-0474

   S_{atm} \downarrow_{vis} = \left( 1 - R_{vis} \right)\left( \alpha S_{atm} \right)

.. math::
   :label: eq-0475

   S_{atm} \downarrow_{nir} = \left( 1 - R_{nir} \right)\left\lbrack (1 - \alpha)S_{atm} \right\rbrack

where :math:`\alpha`, the ratio of visible to total incident solar radiation, is assumed to be

.. math::
   :label: eq-0476

   \alpha = \frac{S_{atm} \downarrow_{vis}^{\mu} + S_{atm} \downarrow_{vis}}{S_{atm}} = 0.5

The ratio of direct to total incident radiation in the visible :math:`R_{vis}` is

.. math::
   :label: eq-0477

   R_{vis} = a_{0} + a_{1} \times \alpha S_{atm} + a_{2} \times \left( \alpha S_{atm} \right)^{2} + a_{3} \times \left( \alpha S_{atm} \right)^{3} 0.01 \leq R_{vis} \leq 0.99

and in the near-infrared :math:`R_{nir}` is

.. math::
   :label: eq-0478

   R_{nir} = b_{0} + b_{1} \times (1 - \alpha)S_{atm} + b_{2} \times \left\lbrack (1 - \alpha)S_{atm} \right\rbrack^{2} + b_{3} \times \left\lbrack (1 - \alpha)S_{atm} \right\rbrack^{3} 0.01 \leq R_{nir} \leq 0.99

where :math:`a_{0} = 0.17639,a_{1} = 0.00380,a_{2} = - 9.0039 \times 10^{- 6},a_{3} = 8.1351 \times 10^{- 9}` and :math:`b_{0} = 0.29548,b_{1} = 0.00504,b_{2} = - 1.4957 \times 10^{- 5},b_{3} = 1.4881 \times 10^{- 8}` are coefficients from polynomial fits to the CAM data.

The additional atmospheric forcing variables required by :numref:`table-atm-input` are derived as follows. The atmospheric reference height :math:`z_{atm}^{'}` (m) is set to 30 m. The directional wind components are derived as :math:`u_{atm} = v_{atm} = \frac{W_{atm}}{\sqrt{2}}`. The potential temperature :math:`\overline{\theta_{atm}}` (K) is set to the atmospheric temperature :math:`T_{atm}`. The atmospheric longwave radiation :math:`L_{atm} \downarrow` (W m\ :sup:`-2`) is derived from the atmospheric vapor pressure :math:`e_{atm}` and temperature :math:`T_{atm}` (Idso 1981) as

.. math::
   :label: eq-0479

   L_{atm} \downarrow = 0.70 + 5.95 \times 10^{- 5} \times 0.01e_{atm}\exp\left( \frac{1500}{T_{atm}} \right)\sigma T_{atm}^{4}

where

.. math::
   :label: eq-0480

   e_{atm} = \frac{P_{atm}q_{atm}}{0.622 + 0.378q_{atm}}

and :math:`\sigma` is the Stefan-Boltzmann constant (W m\ :sup:`-2` K\ :sup:`-4`) (:numref:`table-physical-constants`). The fraction of precipitation :math:`P` (mm s\ :sup:`-1`) falling as rain and/or snow is

:math:`q_{rain} = P\left( f_{P} \right)`,

.. math::
   :label: eq-0481

   q_{snow} = P\left( 1 - f_{P} \right)

where

.. math::
   :label: eq-0482

   f_{P} = 0 < 0.5\left( T_{atm} - T_{f} \right) < 1

If the user wishes to provide atmospheric forcing data from another source, the data format outlined above will need to be followed with the following exceptions. The data atmosphere model will accept a user-supplied relative humidity :math:`RH` (%) and derive specific humidity :math:`q_{atm}` (kg kg\ :sup:`-1`) from

.. math::
   :label: eq-0483

   q_{atm} = \frac{0.622e_{atm}}{P_{atm} - 0.378e_{atm}}

where the atmospheric vapor pressure :math:`e_{atm}` (Pa) is derived from the water (:math:`T_{atm} > T_{f}`) or ice (:math:`T_{atm} \leq T_{f}`) saturation vapor pressure :math:`e_{sat}^{T_{atm}}` as :math:`e_{atm} = \frac{RH}{100}e_{sat}^{T_{atm}}` where :math:`T_{f}` is the freezing temperature of water (K) (:numref:`table-physical-constants`), and :math:`P_{atm}` is the pressure at height :math:`z_{atm}` (Pa). The data atmosphere model will also accept a user-supplied dew point temperature :math:`T_{dew}` (K) and derive specific humidity :math:`q_{atm}` from

.. math::
   :label: eq-0484

   q_{atm} = \frac{0.622e_{sat}^{T_{dew}}}{P_{atm} - 0.378e_{sat}^{T_{dew}}}

Here, :math:`e_{sat}^{T}`, the saturation vapor pressure as a function of temperature, is derived from Lowe's (1977) polynomials (section 3.3). If not provided by the user, the atmospheric pressure :math:`P_{atm}` (Pa) is set equal to the standard atmospheric pressure :math:`P_{std} = 101325` Pa, and surface pressure :math:`P_{srf}` (Pa) is set equal to\ :math:`P_{atm}`.

The user may provide the total direct and diffuse solar radiation, :math:`S_{atm} \downarrow^{\mu}` and :math:`S_{atm} \downarrow`. These will be time-interpolated using the procedure described above and then each term equally apportioned into the visible and near-infrared wavebands (e.g., :math:`S_{atm} \downarrow_{vis}^{\mu} = 0.5S_{atm} \downarrow^{\mu}`, :math:`S_{atm} \downarrow_{nir}^{\mu} = 0.5S_{atm} \downarrow^{\mu}`).

Evaluation
==========

Oleson et al. (2008a, b) describe efforts to evaluate the urban model. This includes a quantitative evaluation of model performance at two specific urban sites, an examination of the robustness of the model through sensitivity studies, and a qualitative evaluation of the urban climate produced by the model, with a focus on the characteristics of the simulated heat island. An additional evaluation component not appearing in these two papers is presented below.

Nighttime longwave radiation and surface temperature
----------------------------------------------------

Nighttime net longwave radiation and air temperature data for an urban canyon in the Grandview district of Vancouver, British Columbia (49°N, 123°W) (Nunez and Oke, 1976, 1977) are used to examine the longwave radiation budget and surface temperatures simulated by the model. The canyon is oriented north-south and is located in a mixed light industrial and residential district. The canyon is 79m long, 7.54m wide, and the east and west walls are 7.31m and 5.59m in height, respectively. Walls are concrete, painted flat white with no windows. The canyon floor consists of a 3-5 cm layer of gravel and clay. Weather conditions on the night of September 9-10, 1973 were clear and calm. Air temperature and net longwave radiation measured at about 0.3m above the midpoint of the canyon floor and from the mid-height of each wall are compared with simulated canyon floor and wall surface temperature and net longwave radiation.

The observation site has been used to validate other urban models such as SHIM (Surface Heat Island Model) (Johnson et al. 1991), the Town Energy Budget (TEB) scheme (Masson 2000), NSLUCM (Noah land surface model/Single-layer Urban Canopy Model) (Kusaka et al. 2001), and VUCM (Vegetated Urban Canopy Model) (Lee and Park 2007). Published data from Lee and Park (2007) were used to determine input parameters for the urban model as these data appeared to produce the best simulations compared to observations (:numref:`table-grandview-parameters`). The canyon floor was modeled as a sandy clay soil with no moisture content. No anthropogenic fluxes were prescribed. Atmospheric wind speed at 10m height was set to 2 m s\ :sup:`-1` and specific humidity to 0.01 kg kg\ :sup:`-1` throughout the simulation (Lee and Park 2007). Atmospheric air temperature was initialized at 19°C and set to the calculated canyon air temperature on subsequent time steps to maintain a neutral temperature profile (no thermal turbulent fluxes between the canyon and the atmosphere) (Masson 2000). Specific humidity of canyon air is set to the atmospheric specific humidity. Downward longwave radiation was initialized to 339 W m\ :sup:`-2` and decreased linearly with the atmospheric air temperature (Masson 2000). Initial wall and canyon floor temperatures were set to 18.35°C and 18.5°C, respectively per Johnson et al. (1991).

.. list-table:: Urban model parameters for the Grandview site
   :widths: auto
   :header-rows: 0
   :name: table-grandview-parameters


   * - Data
     - Symbol
     - Default Value
     - Units
   * - Percent urban
     - -
     - 100
     - %
   * - Canyon height to width ratio
     - :math:`\frac{H}{W}`
     - 0.85
     - -
   * - Roof fraction
     - :math:`W_{roof}`
     - 0.00
     - -
   * - Pervious road fraction
     - :math:`f_{prvrd}`
     - 1.00
     - -
   * - Emissivity of pervious road
     - :math:`\varepsilon_{imprvrd}`
     - 0.98
     - -
   * - Emissivity of sunlit and shaded walls
     - :math:`\varepsilon_{wall}`
     - 0.94
     - -
   * - Building height
     - :math:`H`
     - 6.45
     - m
   * - Wall thermal conductivity
     - :math:`\lambda_{wall,i = 1,10}`
     - 0.81
     - W m\ :sup:`-1` K\ :sup:`-1`
   * - Pervious road thermal conductivity
     - :math:`\lambda_{imprvrd,i = 1,10}`
     - Soil texture (Oleson et al. 2004)
     - W m\ :sup:`-1` K\ :sup:`-1`
   * - Wall volumetric heat capacity
     - :math:`c_{wall,i = 1,10}`
     - 1.0
     - MJ m\ :sup:`-3` K\ :sup:`-1`
   * - Pervious road volumetric heat capacity
     - :math:`c_{imprvrd,i = 1,10}`
     - Soil texture (Oleson et al. 2004)
     - MJ m\ :sup:`-3` K\ :sup:`-1`
   * - Percent sand, percent clay of pervious road (soil)
     - :math:`\% sand,\% clay`
     - 52% sand, 48% clay
     - %
   * - Wall thickness
     - :math:`\Delta z_{wall}`
     - 0.3
     - m


..

:numref:`fig-grandview-validation` shows the simulated surface temperatures and net longwave radiation for the walls and canyon floor compared to observations. The urban model does a good job simulating the nighttime cooling of canyon surfaces (note that the simulated west and east wall surface temperatures are the same). Temperature differences from observations are less than 1°C at all times. Net longwave radiation is also well simulated, differences from observations are less than about 3 W m\ :sup:`-2` for the west wall and canyon floor. The simulated net longwave radiation for the east wall is biased high by up to 7 W m\ :sup:`-2`. These results are quite similar to those from VUCM and generally slightly better than the models of Masson (2000), Johnson et al. (1991), and Kusaka et al. (2001) which generally have warmer surface temperatures as noted by Lee and Park (2007). However, one important difference between Lee and Park (2007) and the other studies is that the thermal admittance prescribed for the canyon floor is substantially lower in VUCM. When higher thermal admittance is prescribed in the urban model, warmer surface temperatures are simulated consistent with the other studies.

   .. figure:: image19.png
      :width: 5.29167in
      :height: 6.83333in
      :name: fig-grandview-validation

      Simulated surface temperatures (solid lines) and net longwave radiation (dashed lines) compared to observations (circles) for a) west (east-facing) wall, b) east wall, and c) canyon floor for the night of September 9-10, 1973 in an urban canyon in the Grandview district of Vancouver, British Columbia. Observed data were digitized from Figure 5 in Johnson et al. (1991).

References
==========

Anderson, E.A. 1976. A point energy and mass balance model of a snow cover. NOAA Technical Report NWS 19, Office of Hydrology, National Weather Service, Silver Spring, MD.

Arnfield, A.J. 2003. Two decades of urban climate research: a review of turbulence, exchanges of energy and water, and the urban heat island. Int. J. Climatol. 23:1-26.

Arya, S.P. 2001. Introduction to Meteorology. Academic Press, San Diego, CA.

Atkinson, B.W. 2003. Numerical modeling of urban heat-island intensity. Bound.-Layer Meteor. 109:285-310.

Auch, R., Taylor, J., and Acevedo, W. 2004. Urban growth in American cities: Glimpses of U.S. urbanization. Circular 1252, U.S. Geological Survey, U.S. Department of the Interior, 52 pp.

Avissar, R. 1996. Potential effects of vegetation on the urban thermal environment. Atmos. Environ. 30:437-448.

Best, M.J. 2005. Representing urban areas within operational numerical weather prediction models. Bound.-Layer Meteor. 114:91-109.

Best, M.J. 2006. Progress towards better weather forecasts for city dwellers: from short range to climate change. Theor. Appl. Climatol. 84:47-55. DOI:10.1007/s00704-005-0143-2.

Betts, R.A. 2001. Biogeophysical impacts of land use on present-day climate: near-surface temperature change and radiative forcing. Atmos. Sci. Lett. 2:39-51. DOI:10.1006/asle.2001.0023.

Beven, K.J., and Kirkby, M.J. 1979. A physically based variable contributing area model of basin hydrology. Hydrol. Sci. Bull. 24: 43-69.

Bonan, G.B. 1996. A land surface model (LSM version 1.0) for ecological, hydrological, and atmospheric studies: technical description and user's guide. NCAR Technical Note NCAR/TN-417+STR, National Center for Atmospheric Research, Boulder, CO, 150 pp.

Bonan, G.B. 2002. Ecological climatology: concepts and applications. Cambridge University Press, 678 pp.

Bounoua, L., DeFries, R., Collatz, G.J., Sellers, P., and Khan, H. 2002. Effects of land cover conversion on surface climate. Clim. Change 52:29-64.

Brovkin, V., Sitch, S., von Bloh, W., Claussen, M., Bauer, E., and Cramer, W. 2004. Role of land cover changes for atmospheric CO\ :sub:`2` increase and climate change during the last 150 years. Global Change Biol. 10:1253-1266. DOI:10.1111/j.1365-2486.2004.00812.x.

Brown, M. 2000. Urban parameterizations for mesoscale meteorological models. pp. 193-255. In: Z. Boybeyi (editor) Mesoscale Atmospheric Dispersion. WIT Press, Southampton, Boston.

Changnon, S.A. 1992. Inadvertent weather modification in urban areas: Lessons for global climate change. Bull. Amer. Meteor. Soc. 73:619-627.

Chylek, P., Srivastava, V., Cahenzli, L., Pinnick, R.G., Dod, R.L., Novakov, T., Cook, T.L., and Hinds, B.D. 1987. Aerosol and graphitic carbon content of snow. J. Geophys. Res. 92:9801-9809.

CIESIN (Center for International Earth Science Information Network), Columbia University; International Food Policy Research Institute (IPFRI), the World Bank; and Centro Internacional de Agricultura Tropical (CIAT), 2004. Global Rural-Urban Mapping Project (GRUMP): Urban Extents. Palisades, NY: CIESIN, Columbia University. Available at http://beta.sedac.ciesin.columbia.edu/gpw.

Clapp, R.B., and Hornberger, G.M. 1978. Empirical equations for some soil hydraulic properties. Water Resour. Res. 14:601-604.

Clauser, C., and Huenges, E. 1995. Thermal conductivity of rocks and minerals. pp. 105-126. In: T. J. Ahrens (editor) Rock Physics and Phase Relations: A Handbook of Physical Constants. Washington, D.C.

Comrie, A.C. 2000. Mapping a wind-modified urban heat island in Tucson, Arizona (with comments on integrating research and undergraduate learning). Bull. Amer. Meteor. Soc. 81:2417-2431.

Copeland, J.H., Pielke, R.A., and Kittel, T.G.F. 1996. Potential climatic impacts of vegetation change: a regional modeling study. J. Geophys. Res. 101:7409-7418.

Cosby, B.J., Hornberger, G.M., Clapp, R.B., and Ginn, T.R. 1984. A statistical exploration of the relationships of soil moisture characteristics to the physical properties of soils. Water Resour. Res. 20:682-690.

Cramer, W., et al. 2001. Global response of terrestrial ecosystem structure and function to CO\ :sub:`2` and climate change: results from six dynamic global vegetation models. Global Change Biol. 7:357-373.

Dai, Y., and Zeng, Q. 1997. A land surface model (IAP94) for climate studies. Part I: formulation and validation in off-line experiments. Adv. Atmos. Sci. 14:433-460.

Dandou, A., Tombrou, M., Akylas, E., Soulakellis, N., and Bossioli, E. 2005. Development and evaluation of an urban parameterization scheme in the Penn State/NCAR Mesoscale Model (MM5). J. Geophys. Res. 110:D10102. DOI: 10.1029/2004JD005192.

Dobson, J.E., Bright, E.A., Coleman, P.R., Durfee, R.C., and Worley, B.A. 2000. LandScan: A global population database for estimating populations at risk. Photo­grammetric Engineering and Remote Sensing 66(7):849-857.

de Vries, D.A. 1963. Thermal Properties of Soils. In: W.R. van Wijk (editor) Physics of the Plant Environment. North-Holland, Amsterdam.

Eastman, J.L., Coughenour, M.B., and Pielke Sr., R.A. 2001. The regional effects of CO\ :sub:`2` and landscape change using a coupled plant and meteorological model. Global Change Biol. 7:797-815.

Elvidge, C.D., Sutton, P.C., Wagner, T.W., Ryzner, R., Vogelmann, J.E., Goetz, S.J., Smith, A.J., Jantz, C., Seto, K.C., Imhoff, M.L., Wang, Y.Q., Milesi, C. and Nemani, R. 2004. Urbanization. pp. 315-328. In: G. Gutman, A. C. Janetos, C.O. Justice, E.F. Moran, J.F. Mustard, R.R. Rindfuss, D. Skole, B.L. Turner II, and M.A Cochrane (editors), Land Change Science: Observing, Monitoring and Understanding Trajectories of Change on the Earth's Surface. Kluwer Academic Publishers, The Netherlands.

Entekhabi, D., and Eagleson, P.S. 1989. Land surface hydrology parameterization for atmospheric general circulation models including subgrid scale spatial variability. J. Climate 2:816-831.

Farouki, O.T. 1981. The thermal properties of soils in cold regions. Cold Regions Sci. and Tech. 5:67-75.

Feddema, J.J., Oleson, K.W., Bonan, G., Mearns, L.O., Buja, L.E., Meehl, G.A., and Washington, W.M. 2005: The importance of land cover change in simulating future climates. Science 310:1674-1678.

Flatau, P.J., Walko, R.L., and Cotton, W.R. 1992. Polynomial fits to saturation vapor pressure. J. Appl. Meteor. 31:1507-1513.

Foley, J.A., et al. 2005. Global consequences of land use. Science 309:570-574.

Fu, C. 2003. Potential impacts of human-induced land cover change in East Asia monsoon. Glob. Planet. Change 37:219-229. DOI:10.1016/S0921-8181(02)00207-2.

Grimmond, C.S.B., Cleugh, H.A., and Oke, T.R. 1991. An objective urban heat storage model and its comparison with other schemes. Atmos. Environ. 25B:311-326.

Grimmond, C.S.B., and Oke, T.R. 1999. Aerodynamic properties of urban areas derived from analysis of surface form. J. Appl. Meteor. 38:1262-1292.

Grimmond, C.S.B., and Oke, T.R. 2002. Turbulent heat fluxes in urban areas: observations and a local-scale urban meteorological parameterization scheme (LUMPS). J. Appl. Meteor. 41:792-810.

Harman, I.N., Best, M.J., and Belcher, S.E. 2004. Radiative exchange in an urban street canyon. Bound.-Layer Meteor. 110:301-316.

Houghton, J.T., Ding, Y., Griggs, D.J., Noguer, M., van der Linden, P.J., Dai, X., Maskell, K., and Johnson, C.A. (editors) 2001. Climate Change 2001: The Scientific Basis. Cambridge University Press, 881 pp.

Ichinose, T., Shimodozono, K., Hanaki, K. 1999. Impact of anthropogenic heat on urban climate in Tokyo. Atmos. Environ. 33:3897-3909.

Idso, S.B. 1981. A set of equations for full spectrum and 8- to 14-:math:`\mu`\ m and 10.5- to 12.5-:math:`\mu`\ m thermal radiation from cloudless skies. Water Resour. Res. 17:295-304.

Jackson, T.L., Feddema, J.J., Oleson, K.W., Bonan, G.B., and Bauer, J.T. 2010. Parameterization of urban characteristics for global climate modeling. A. Assoc. Am. Geog., in press.

Jin, M., Dickinson, R.E., and Zhang, D.-L. 2005. The footprint of urban areas on global climate as characterized by MODIS. J. Climate 18:1551-1565.

Johnson, G.T., Oke, T.R., Lyons, T.J., Steyn, D.G., Watson, I.D., and Voogt, J.A. 1991. Simulation of surface urban heat islands under 'ideal' conditions at night. Part 1: Theory and tests against field data. Bound.-Layer Meteor. 56:275-294.

Jordan, R. 1991. A One-dimensional Temperature Model for a Snow Cover: Technical Documentation for SNTHERM.89. U.S. Army Cold Regions Research and Engineering Laboratory, Special Report 91-16.

Kalnay, E., and Cai, M. 2003. Impact of urbanization and land-use change on climate. Nature 423:528-531.

Kusaka, H., Kondo, H., Kikegawa, Y., and Kimura, F. 2001. A simple single-layer urban canopy model for atmospheric models: comparison with multi-layer and slab models. Bound.-Layer Meteor. 101:329-358.

Kusaka, H., and Kimura, F. 2004. Thermal effects of urban canyon structure on the nocturnal heat island: numerical experiment using a mesoscale model coupled with an urban canopy model. J. Appl. Meteor. 43:1899-1910.

Landsberg, H.E. 1981. The Urban Climate. New York, Academic Press, 275 pp.

Lawrence, D.M., and Slater, A.G. 2008. Incorporating organic soil into a global climate model. Clim. Dyn. 30. DOI:10.1007/s00382-007-0278-1.

Lee, S.-H., and Park, S.-U. 2007. A vegetated urban canopy model for meteorological and environmental modeling. Bound.-Layer Meteor. DOI:10:1007/s10546-007-9221-6.

Lemonsu, A., and Masson, V. 2002. Simulation of a summer urban breeze over Paris. Bound.-Layer Meteor. 104:463-490.

Lemonsu, A., Grimmond, C.S.B., and Masson, V. 2004. Modeling the surface energy balance of the core of an old Mediterranean city: Marseille. J. Appl. Meteor. 43:312-327.

Lowe, P.R. 1977. An approximating polynomial for the computation of saturation vapor pressure. J. Appl. Meteor. 16:100-103.

Macdonald, R.W., Griffiths, R.F., and Hall, D.J. 1998. An improved method for the estimation of surface roughness of obstacle arrays. Atmos. Environ. 32:1857-1864.

Marshall, S.E. 1989. A physical parameterization of snow albedo for use in climate models. NCAR Cooperative Thesis NCAR/CT-123, National Center for Atmospheric Research, Boulder, CO.

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

Oke, T.R. 1987. Boundary Layer Climates (2\ :sup:`nd` edition). London, Routledge, 435 pp.

Oke, T.R., and Cleugh, H.A. 1987. Urban heat storage derived as energy balance residuals. Bound.-Layer Meteor. 39:233-245.

Oke, T.R., Johnson, G.T., Steyn, D.G., and Watson, I.D. 1991. Simulation of surface urban heat islands under "ideal" conditions at night, part 2: diagnosis of causation. Bound.-Layer Meteor. 56:339-358.

Oleson, K.W., et al. 2004. Technical description of the Community Land Model (CLM). NCAR Technical Note NCAR/TN-461+STR, National Center for Atmospheric Research, Boulder, CO, 173 pp.

Oleson, K.W., Bonan, G.B., Feddema, J., Vertenstein, M., and Grimmond, C.S.B. 2008a. An urban parameterization for a global climate model. 1. Formulation and evaluation for two cities. J. Appl. Meteor. Clim. 47:1038-1060.

Oleson, K.W., Bonan, G.B., Feddema, J., and Vertenstein, M. 2008b. An urban parameterization for a global climate model. 2. Sensitivity to input parameters and the simulated urban heat island in offline simulations. J. Appl. Meteor. Clim. 47:1061-1076.

Oleson, K.W., Bonan, G.B., and Feddema, J. 2010a. The effects of white roofs on urban temperature in a global climate model. Geophys. Res. Lett. 37:L03701. DOI:10.1029/2009GL042194.

Oleson, K.W., et al. 2010b. Technical description of version 4.0 of the Community Land Model (CLM). NCAR Technical Note NCAR/TN-478+STR, National Center for Atmospheric Research, Boulder, CO, 257 pp.

Otte, T.L., Lacser, A., Dupont, S., and Ching, J.K.S. 2004. Implementation of an urban canopy parameterization in a mesoscale meteorological model. J. Appl. Meteor. 43:1648-1665.

Panofsky, H.A., and Dutton, J.A. 1984. Atmospheric Turbulence: Models and Methods for Engineering Applications. John Wiley and Sons, New York.

Pielke Sr., R.A., Marland, G., Betts, R.A., Chase, T.N., Eastman, J.L., Niles, J.O., Niyogi, D.D.S., and Running, S.W. 2002. The influence of land-use change and landscape dynamics on the climate system: relevance to climate-change policy beyond the radiative effect of greenhouse gases. Phil. Trans. R. Soc. Lond. A, 360:1705-1719.

Piringer, M., Grimmond, C.S.B., Joffre, S.M., Mestayer, P., Middleton, D.R., Rotach, A. Baklanov, M.W., De Ridder, K., Ferreira, J., Guilloteau, E., Karppinen, A., Martilli, A., Masson, V., and Tombrou, M. 2002. Investigating the surface energy balance in urban areas - recent advances and future needs. Water, Air, Soil Pollution: Focus, 2:1-16.

Press, W.H., Teukolsky, S.A., Vetterling, W.T., and Flannery, B.P. 1992. Numerical Recipes in FORTRAN: The Art of Scientific Computing. Cambridge University Press, New York.

Qian, T., Dai, A., Trenberth, K.E., and Oleson, K.W. 2006. Simulation of global land surface conditions from 1948 to 2004: Part I: Forcing data and evaluations. J. Hydrometeor. 7:953-975.

Rowley, F.B., Algren, A.B., and Blackshaw, J.L. 1930. Surface conductances as affected by air velocity, temperature, and character of surface. ASHRAE Trans. 36:429-446.

Sailor, D.J. 1995. Simulated urban climate response to modifications in surface albedo and vegetative cover. J. App. Meteor. 34:1694-1704.

Sailor, D.J., and Lu, L. 2004. A top-down methodology for developing diurnal and seasonal anthropogenic heating profiles for urban areas. Atmos. Environ. 38:2737-2748.

Sakakibara, Y. 1996. A numerical study of the effect of urban geometry upon the surface energy budget. Atmos. Environ. 30:487-496.

Sellers, P.J., Dickinson, R.E., Randall, D.A., Betts, A.K., Hall, F.G., Berry, J.A., Collatz, G.J., Denning, A.S., Mooney, H.A., Nobre, C.A., Sato, N., Field, C.B., and Henderson-Sellers, A. 1997. Modeling the exchanges of energy, water, and carbon between continents and the atmosphere. Science 275:502-509.

Shepherd, J.M. 2005. A review of current investigations of urban-induced rainfall and recommendations for the future. Earth Interact. 9:1-27.

Sparrow, E.M., and Cess, R.D. 1978. Radiation Heat Transfer. Hemisphere Publishing Corporation, 366 pp.

Stohlgren, T.J., Chase, T.N., Pielke Sr., R.A., Kittel, T.G.F., and Baron, J.S. 1998. Evidence that local land use practices influence regional climate, vegetation, and stream flow patterns in adjacent natural areas. Global Change Biol. 4:495-504.

Stull, R.B. 1988. An Introduction to Boundary Layer Meteorology. Kluwer Academic Publishers, Dordrecht.

Taha, H. 1999. Modifying a mesoscale meteorological model to better incorporate urban heat storage: a bulk-parameterization approach. J. Appl. Meteor. 38:466-473.

Upmanis, H., Eliasson, I., and Lindqvist, S. 1998. The influence of green areas on nocturnal temperatures in a high latitude city (Göteborg, Sweden). Int. J. Climatol. 18:681-700.

Wang, H., Pitman, A.J., Zhao, M., and Leemans, R. 2003. The impact of land-cover modification on the June meteorology of China since 1700, simulated using a regional climate model. Int. J. Climatol. 23:511-527.

Yang, Z.-L. 1998. Technical note of a 10-layer soil moisture and temperature model. Unpublished manuscript.

Zeng, X., Zhao, M., and Dickinson R.E. 1998. Intercomparison of bulk aerodynamic algorithms for the computation of sea surface fluxes using the TOGA COARE and TAO data. J. Climate 11:2628-2644.

Zeng, X., and Decker, M. 2009. Improving the numerical solution of soil moisture-based Richards equation for land models with a deep or shallow water table. J. Hydrometeor. 10:308-319.

Zhou, L., Dickinson, R.E., Tian, Y., Fang, J., Li, Q., Kaufmann, R.K., Tucker, C.J., and Myneni, R.B. 2004. Evidence for a significant urbanization effect on climate in China. Proc. Natl. Acad. Sci. U.S.A. 101:9540-9544.

