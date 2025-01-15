.. _rst_Surface Characterization, Vertical Discretization, and Model Input Requirements:

Surface Characterization, Vertical Discretization, and Model Input Requirements
===================================================================================

.. _Surface Characterization:

Surface Characterization
-----------------------------

.. _Surface Heterogeneity and Data Structure:

Surface Heterogeneity and Data Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Spatial land surface heterogeneity in CLM is represented as a nested subgrid hierarchy in which grid cells are composed of multiple land units, snow/soil columns, and PFTs (:numref:`Figure CLM subgrid hierarchy`). Each grid cell can have a different number of land units, each land unit can have a different number of columns, and each column can have multiple PFTs. The first subgrid level, the land unit, is intended to capture the broadest spatial patterns of subgrid heterogeneity. The current land units are glacier, lake, urban, vegetated, and crop (when the crop model option is turned on). The land unit level can be used to further delineate these patterns. For example, the urban land unit is divided into density classes representing the tall building district, high density, and medium density urban areas.

The second subgrid level, the column, is intended to capture potential variability in the soil and snow state variables within a single land unit. For example, the vegetated land unit could contain several columns with independently evolving vertical profiles of soil water and temperature. Similarly, the managed vegetation land unit can be divided into two columns, irrigated and non-irrigated. The default snow/soil column is represented by 25 layers for ground (with up to 20 of these layers classified as soil layers and the remaining layers classified as bedrock layers) and up to 10 layers for snow, depending on snow depth. The central characteristic of the column subgrid level is that this is where the state variables for water and energy in the soil and snow are defined, as well as the fluxes of these components within the soil and snow. Regardless of the number and type of PFTs occupying space on the column, the column physics operates with a single set of upper boundary fluxes, as well as a single set of transpiration fluxes from multiple soil levels. These boundary fluxes are weighted averages over all PFTs. Currently, for lake and vegetated land units, a single column is assigned to each land unit. The crop land unit is split into irrigated and unirrigated columns with a single crop occupying each column. The urban land units have five columns (roof, sunlit walls and shaded walls, and pervious and impervious canyon floor) (Oleson et al. 2010b). The glacier land unit is separated into up to 10 elevation classes.

.. _Figure CLM subgrid hierarchy:

.. Figure:: image1.png

  Configuration of the CLM subgrid hierarchy.  Box in upper right shows hypothetical subgrid distribution for a single grid cell.  Note that the Crop land unit is only used when the model is run with the crop model active. Abbreviations: TBD – Tall Building District; HD – High Density; MD – Medium Density, G – Glacier, L – Lake, U – Urban, C – Crop, V – Vegetated, PFT – Plant Functional Type, Irr – Irrigated, UIrr – Unirrigated.  Red arrows indicate allowed land unit transitions.  Purple arrows indicate allowed patch-level transitions.

The third subgrid level is referred to as the patch level. Patches can be PFTs or bare ground on the vegetated land unit and crop functional types (CFTs) on the crop land unit. The patch level is intended to capture the biogeophysical and biogeochemical differences between broad categories of plants in terms of their functional characteristics. On the vegetated land unit, up to 16 possible PFTs that differ in physiology and structure may coexist on a single column. All fluxes to and from the surface are defined at the PFT level, as are the vegetation state variables (e.g. vegetation temperature and canopy water storage). On the crop land unit, typically, different crop types can be represented on each crop land unit column (see Chapter :numref:`rst_Crops and Irrigation` for details).

In addition to state and flux variable data structures for conserved components at each subgrid level (e.g., energy, water, carbon), each subgrid level also has a physical state data structure for handling quantities that are not involved in conservation checks (diagnostic variables). For example, the urban canopy air temperature and humidity are defined through physical state variables at the land unit level, the number of snow layers and the soil roughness lengths are defined as physical state variables at the column level, and the leaf area index and the fraction of canopy that is wet are defined as physical state variables at the PFT level.

The standard configuration of the model subgrid hierarchy is illustrated in :numref:`Figure CLM subgrid hierarchy`. Here, only four PFTs are shown associated with the single column beneath the vegetated land unit but up to sixteen are possible. The crop land unit is present only when the crop model is active.

Note that the biogeophysical processes related to soil and snow require PFT level properties to be aggregated to the column level. For example, the net heat flux into the ground is required as a boundary condition for the solution of snow/soil temperatures (Chapter :numref:`rst_Soil and Snow Temperatures`). This column level property must be determined by aggregating the net heat flux from all PFTs sharing the column. This is generally accomplished in the model by computing a weighted sum of the desired quantity over all PFTs whose weighting depends on the PFT area relative to all PFTs, unless otherwise noted in the text.

.. _Vegetation Composition:

Vegetation Composition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Vegetated surfaces are comprised of up to 15 possible plant functional types (PFTs) plus bare ground (:numref:`Table Plant functional types`). An additional PFT is added if the irrigation model is active and six additional PFTs are added if the crop model is active (Chapter :numref:`rst_Crops and Irrigation`). In :numref:`Table Plant functional types`, IVT = 0,14 refers to the index of PCT_NAT_VEG on the surface dataset while IVT = 15,18 refers to the index of PCT_CFT on the surface dataset. These plant types differ in leaf and stem optical properties that determine reflection, transmittance, and absorption of solar radiation (:numref:`Table Plant functional type optical properties`), root distribution parameters that control the uptake of water from the soil (:numref:`Table Plant functional type root distribution parameters`), aerodynamic parameters that determine resistance to heat, moisture, and momentum transfer (:numref:`Table Plant functional type aerodynamic parameters`), and photosynthetic parameters that determine stomatal resistance, photosynthesis, and transpiration (:numref:`Table Plant functional type (PFT) stomatal conductance parameters`, :numref:`Table Temperature dependence parameters for C3 photosynthesis`). The composition and abundance of PFTs within a grid cell can either be prescribed as time-invariant fields (e.g., using the present day dataset described in section 21.3.3) or can evolve with time if the model is run in transient landcover mode (Chapter :numref:`rst_Transient Landcover Change`).

.. _Table Plant functional types:

.. table:: Plant functional types

 +-----+--------------------------------------------------------------+-------------------+
 | IVT | Plant functional type                                        | Acronym           |
 +=====+==============================================================+===================+
 | 0   | Bare Ground                                                  | NET Temperate     |
 +-----+--------------------------------------------------------------+-------------------+
 | 1   | Needleleaf evergreen tree – temperate                        | NET Temperate     |
 +-----+--------------------------------------------------------------+-------------------+
 | 2   | Needleleaf evergreen tree - boreal                           | NET Boreal        |
 +-----+--------------------------------------------------------------+-------------------+
 | 3   | Needleleaf deciduous tree – boreal                           | NDT Boreal        |
 +-----+--------------------------------------------------------------+-------------------+
 | 4   | Broadleaf evergreen tree – tropical                          | BET Tropical      |
 +-----+--------------------------------------------------------------+-------------------+
 | 5   | Broadleaf evergreen tree – temperate                         | BET Temperate     |
 +-----+--------------------------------------------------------------+-------------------+
 | 6   | Broadleaf deciduous tree – tropical                          | BDT Tropical      |
 +-----+--------------------------------------------------------------+-------------------+
 | 7   | Broadleaf deciduous tree – temperate                         | BDT Temperate     |
 +-----+--------------------------------------------------------------+-------------------+
 | 8   | Broadleaf deciduous tree – boreal                            | BDT Boreal        |
 +-----+--------------------------------------------------------------+-------------------+
 | 9   | Broadleaf evergreen shrub - temperate                        | BES Temperate     |
 +-----+--------------------------------------------------------------+-------------------+
 | 10  | Broadleaf deciduous shrub – temperate                        | BDS Temperate     |
 +-----+--------------------------------------------------------------+-------------------+
 | 11  | Broadleaf deciduous shrub – boreal                           | BDS Boreal        |
 +-----+--------------------------------------------------------------+-------------------+
 | 12  | C\ :sub:`3` arctic grass                                     | -                 |
 +-----+--------------------------------------------------------------+-------------------+
 | 13  | C\ :sub:`3` grass                                            | -                 |
 +-----+--------------------------------------------------------------+-------------------+
 | 14  | C\ :sub:`4` grass                                            | -                 |
 +-----+--------------------------------------------------------------+-------------------+
 | 15  | C\ :sub:`3` Unmanaged Rainfed Crop                           | UCrop UIrr        |
 +-----+--------------------------------------------------------------+-------------------+
 | 16  | :sup:`1`\ C\ :sub:`3` Unmanaged Irrigated Crop               | UCrop Irr         |
 +-----+--------------------------------------------------------------+-------------------+
 | 17  | :sup:`2`\ Managed Rainfed Crop                               | Crop UIrr         |
 +-----+--------------------------------------------------------------+-------------------+
 | 18  | :sup:`2`\ Managed Irrigated Crop                             | Crop Irr          |
 +-----+--------------------------------------------------------------+-------------------+

:sup:`1`\ Only used if irrigation is active (Chapter :numref:`rst_Crops and Irrigation`).
:sup:`2`\ Only used if crop model is active (see Chapter :numref:`rst_Crops and Irrigation` for list of represented crops).

.. _Vegetation Structure:

Vegetation Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^

Vegetation structure is defined by leaf and stem area indices (:math:`L,\, S`) and canopy top and bottom heights (:math:`z_{top}`,\ :math:`z_{bot}` ). Separate leaf and stem area indices and canopy heights are prescribed or calculated for each PFT. Daily leaf and stem area indices are obtained from griddeddatasets of monthly values (section :numref:`Surface Data`). Canopy top and bottom heights for trees are from ICESat (:ref:`Simard et al. (2011) <Simardetal2011>`). Canopy top and bottom heights for short vegetation are obtained from gridded datasets but are invariant in space and time and were obtained from PFT-specific values (:ref:`Bonan et al. (2002a) <Bonanetal2002a>`) (:numref:`Table Plant functional type canopy top and bottom heights`). When the biogeochemistry model is active, vegetation state (LAI, SAI, canopy top and bottom heights) are calculated prognostically (see Chapter :numref:`rst_Vegetation Phenology and Turnover`).

.. _Table Plant functional type canopy top and bottom heights:

.. table:: Plant functional type canopy top and bottom heights

 +--------------------------------------------------------------+-------------------+-------------------+
 | Plant functional type                                        | :math:`z_{top}`   | :math:`z_{bot}`   |
 +==============================================================+===================+===================+
 | BES Temperate                                                | 0.5               | 0.1               |
 +--------------------------------------------------------------+-------------------+-------------------+
 | BDS Temperate                                                | 0.5               | 0.1               |
 +--------------------------------------------------------------+-------------------+-------------------+
 | BDS Boreal                                                   | 0.5               | 0.1               |
 +--------------------------------------------------------------+-------------------+-------------------+
 | C\ :sub:`3` arctic grass                                     | 0.5               | 0.01              |
 +--------------------------------------------------------------+-------------------+-------------------+
 | C\ :sub:`3` grass                                            | 0.5               | 0.01              |
 +--------------------------------------------------------------+-------------------+-------------------+
 | C\ :sub:`4` grass                                            | 0.5               | 0.01              |
 +--------------------------------------------------------------+-------------------+-------------------+
 | UCrop UIrr                                                   | 0.5               | 0.01              |
 +--------------------------------------------------------------+-------------------+-------------------+
 | UCrop Irr                                                    | 0.5               | 0.01              |
 +--------------------------------------------------------------+-------------------+-------------------+
 | Crop UIrr                                                    | 0.5               | 0.01              |
 +--------------------------------------------------------------+-------------------+-------------------+
 | Crop Irr                                                     | 0.5               | 0.01              |
 +--------------------------------------------------------------+-------------------+-------------------+

.. _Phenology and vegetation burial by snow:

Phenology and vegetation burial by snow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When the biogeochemistry model is inactive, leaf and stem area indices (m\ :sup:`2` leaf area m\ :sup:`-2` ground area) are updated daily by linearly interpolating between monthly values. Monthly PFT leaf area index values are developed from the 1-km MODIS-derived monthly grid cell average leaf area index of :ref:`Myneni et al. (2002) <Mynenietal2002>`, as described in :ref:`Lawrence and Chase (2007) <LawrenceChase2007>`. Stem area ndex is calculated from the monthly PFT leaf area index using the methods of :ref:`Zeng et al. (2002) <Zengetal2002>`. The leaf and stem area indices are adjusted for vertical burying by snow (:ref:`Wang and Zeng 2009 <WangZeng2009>`) as

.. math::
   :label: 2.1

   A=A^{*} ( 1-f_{veg}^{sno} )

where :math:`A^{*}` is the leaf or stem area before adjustment for snow, :math:`A` is the remaining exposed leaf or stem area, :math:`f_{veg}^{sno}` is the vertical fraction of vegetation covered by snow

.. math::
   :label: 2.2

   {f_{veg}^{sno} = \frac{z_{sno} -z_{bot} }{z_{top} -z_{bot} }         \qquad {\rm for\; tree\; and\; shrub}} \\
   {f_{veg}^{sno} = \frac{\min \left(z_{sno} ,\, z_{c} \right)}{z_{c} } \qquad {\rm for\; grass\; and\; crop}}

where :math:`z_{sno} -z_{bot} \ge 0,{\rm \; }0\le f_{veg}^{sno} \le 1`, :math:`z_{sno}` is the depth of snow (m) (Chapter :numref:`rst_Snow Hydrology`), and :math:`z_{c} = 0.2` is the snow depth when short vegetation is assumed to be completely buried by snow (m). For numerical reasons, exposed leaf and stem area are set to zero if less than 0.05. If the sum of exposed leaf and stem area is zero, then the surface is treated as snow-covered ground.

.. _Vertical Discretization:

Vertical Discretization
----------------------------
..
 (this was taken from Initialization; is it still needed?
 Vegetated and glacier land units have fifteen vertical layers, while lakes have ten. For soil points, temperature calculations are done over all layers, :math:`N_{levgrnd} =15`, while hydrology calculations are done over the top ten layers, :math:`N_{levsoi} =10`, the bottom five layers being specified as bedrock.

.. _Soil Layers:

Soil Layers
^^^^^^^^^^^^^^^^^^^^^^^^^^

The soil column can be discretized into an arbitrary number of layers. The default vertical discretization (:numref:`Table Soil layer structure`) uses :math:`N_{levgrnd} = 25` layers, of which :math:`N_{levsoi} = 20` are hydrologically and biogeochemically active. The deepest 5 layers are only included in the thermodynamical calculations (:ref:`Lawrence et al. 2008 <Lawrenceetal2008>`) described in Chapter :numref:`rst_Soil and Snow Temperatures`.

The layer structure of the soil is described by the node depth, :math:`z_{i}` (m), the thickness of each layer, :math:`\Delta z_{i}` (m), and the depths at the layer interfaces :math:`z_{h,\, i}` (m).

.. _Table Soil layer structure:

.. table:: Soil layer structure

 +---------------+------------------+------------------------+------------------------+
 | Layer         | :math:`z_{i}`    | :math:`\Delta z_{i}`   | :math:`z_{h,\, i}`     |
 +===============+==================+========================+========================+
 |    1          |   0.010          |   0.020                |   0.020                |
 +---------------+------------------+------------------------+------------------------+
 |    2          |   0.040          |   0.040                |   0.060                |
 +---------------+------------------+------------------------+------------------------+
 |    3          |   0.090          |   0.060                |   0.120                |
 +---------------+------------------+------------------------+------------------------+
 |    4          |   0.160          |   0.080                |   0.200                |
 +---------------+------------------+------------------------+------------------------+
 |    5          |   0.260          |   0.120                |   0.320                |
 +---------------+------------------+------------------------+------------------------+
 |    6          |   0.400          |   0.160                |   0.480                |
 +---------------+------------------+------------------------+------------------------+
 |    7          |   0.580          |   0.200                |   0.680                |
 +---------------+------------------+------------------------+------------------------+
 |    8          |   0.800          |   0.240                |   0.920                |
 +---------------+------------------+------------------------+------------------------+
 |    9          |   1.060          |   0.280                |   1.200                |
 +---------------+------------------+------------------------+------------------------+
 |   10          |   1.360          |   0.320                |   1.520                |
 +---------------+------------------+------------------------+------------------------+
 |   11          |   1.700          |   0.360                |   1.880                |
 +---------------+------------------+------------------------+------------------------+
 |   12          |   2.080          |   0.400                |   2.280                |
 +---------------+------------------+------------------------+------------------------+
 |   13          |   2.500          |   0.440                |   2.720                |
 +---------------+------------------+------------------------+------------------------+
 |   14          |   2.990          |   0.540                |   3.260                |
 +---------------+------------------+------------------------+------------------------+
 |   15          |   3.580          |   0.640                |   3.900                |
 +---------------+------------------+------------------------+------------------------+
 |   16          |   4.270          |   0.740                |   4.640                |
 +---------------+------------------+------------------------+------------------------+
 |   17          |   5.060          |   0.840                |   5.480                |
 +---------------+------------------+------------------------+------------------------+
 |   18          |   5.950          |   0.940                |   6.420                |
 +---------------+------------------+------------------------+------------------------+
 |   19          |   6.940          |   1.040                |   7.460                |
 +---------------+------------------+------------------------+------------------------+
 |   20          |   8.030          |   1.140                |   8.600                |
 +---------------+------------------+------------------------+------------------------+
 |   21          |   9.795          |   2.390                |  10.990                |
 +---------------+------------------+------------------------+------------------------+
 |   22          |  13.328          |   4.676                |  15.666                |
 +---------------+------------------+------------------------+------------------------+
 |   23          |  19.483          |   7.635                |  23.301                |
 +---------------+------------------+------------------------+------------------------+
 |   24          |  28.871          |  11.140                |  34.441                |
 +---------------+------------------+------------------------+------------------------+
 |   25          |  41.998          |  15.115                |  49.556                |
 +---------------+------------------+------------------------+------------------------+

Layer node depth (:math:`z_{i}` ), thickness (:math:`\Delta z_{i}` ), and depth at layer interface (:math:`z_{h,\, i}` ) for default soil column. All in meters.

.. _Depth to Bedrock:

Depth to Bedrock
^^^^^^^^^^^^^^^^^^^^^^^^^^

The hydrologically and biogeochemically active portion of the soil column can be restricted to a thickness less than that of the maximum soil depth. By providing a depth-to-bedrock dataset, which may vary spatially, the number of layers used in the hydrologic and biogeochemical calculations, :math:`N_{bedrock}`, may be specified, subject to the constraint :math:`N_{bedrock} \le N_{levsoi}`. The default depth-to-bedrock values are from :ref:`Pelletier et al. [2016]<Pelletieretal2016>`.

.. _Model Input Requirements:

Model Input Requirements
----------------------------

.. _Atmospheric Coupling:

Atmospheric Coupling
^^^^^^^^^^^^^^^^^^^^^^^^^^

The current state of the atmosphere (:numref:`Table Atmospheric input to land model`) at a given time step is used to force the land model. This atmospheric state is provided by an atmospheric model in coupled mode or from an observed dataset in land-only mode (Chapter :numref:`rst_Land-Only Mode`). The land model then initiates a full set of calculations for surface energy, constituent, momentum, and radiative fluxes. The land model calculations are implemented in two steps. The land model proceeds with the calculation of surface energy, constituent, momentum, and radiative fluxes using the snow and soil hydrologic states from the previous time step. The land model then updates the soil and snow hydrology calculations based on these fluxes. These fields are passed to the atmosphere (:numref:`Table Land model output to atmospheric model`). The albedos sent to the atmosphere are for the solar zenith angle at the next time step but with surface conditions from the current time step.

.. _Table Atmospheric input to land model:

.. table:: Atmospheric input to land model

 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | Field                                                | variable name                                  | units                                           |
 +======================================================+================================================+=================================================+
 | :sup:`1`\ Reference height                           | :math:`z'_{atm}`                               | m                                               |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | Atmosphere model's surface height                    | :math:`z_{surf,atm}`                           | m                                               |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | Zonal wind at :math:`z_{atm}`                        | :math:`u_{atm}`                                | m s\ :sup:`-1`                                  |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | Meridional wind at :math:`z_{atm}`                   | :math:`v_{atm}`                                | m s\ :sup:`-1`                                  |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | Potential temperature                                | :math:`\overline{\theta _{atm} }`              | K                                               |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | Specific humidity at :math:`z_{atm}`                 | :math:`q_{atm}`                                | kg kg\ :sup:`-1`                                |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | Pressure at :math:`z_{atm}`                          | :math:`P_{atm}`                                | Pa                                              |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | Temperature at :math:`z_{atm}`                       | :math:`T_{atm}`                                | K                                               |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | Incident longwave radiation                          | :math:`L_{atm} \, \downarrow`                  | W m\ :sup:`-2`                                  |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | :sup:`2`\ Liquid precipitation                       | :math:`q_{rain}`                               | mm s\ :sup:`-1`                                 |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | :sup:`2`\ Solid precipitation                        | :math:`q_{sno}`                                | mm s\ :sup:`-1`                                 |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | Incident direct beam visible solar radiation         | :math:`S_{atm} \, \downarrow _{vis}^{\mu }`    | W m\ :sup:`-2`                                  |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | Incident direct beam near-infrared solar radiation   | :math:`S_{atm} \, \downarrow _{nir}^{\mu }`    | W m\ :sup:`-2`                                  |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | Incident diffuse visible solar radiation             | :math:`S_{atm} \, \downarrow _{vis}`           | W m\ :sup:`-2`                                  |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | Incident diffuse near-infrared solar radiation       | :math:`S_{atm} \, \downarrow _{nir}`           | W m\ :sup:`-2`                                  |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | Carbon dioxide (CO\ :sub:`2`) concentration          | :math:`c_{a}`                                  | ppmv                                            |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | :sup:`3`\ Aerosol deposition rate                    | :math:`D_{sp}`                                 | kg m\ :sup:`-2` s\ :sup:`-1`                    |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | :sup:`4`\ Nitrogen deposition rate                   | :math:`NF_{ndep\_ s{\it min}n}`                | g (N) m\ :sup:`-2` yr\ :sup:`-1`                |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+
 | :sup:`5`\ Lightning frequency                        | :math:`I_{l}`                                  | flash km\ :sup:`-2` hr\ :sup:`-1`               |
 +------------------------------------------------------+------------------------------------------------+-------------------------------------------------+

:sup:`1`\ The atmospheric reference height received from the atmospheric model :math:`z'_{atm}` is assumed to be the height above the surface as defined by the roughness length :math:`z_{0}` plus displacement height :math:`d`. Thus, the reference height used for flux computations (Chapter :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`) is :math:`z_{atm} =z'_{atm} +z_{0} +d`. The reference heights for temperature, wind, and specific humidity (:math:`z_{atm,\, h}`, :math:`z_{atm,\, {\it m}}`, :math:`z_{atm,\, w}` ) are required. These are set equal to\ :math:`z_{atm}`.

:sup:`2`\ CAM provides convective and large-scale liquid and solid precipitation, which are added to yield total liquid precipitation :math:`q_{rain}` and solid precipitation :math:`q_{sno}`. However, in CLM5, the atmosphere's partitioning into liquid and solid precipitation is ignored. Instead, CLM repartitions total precipitation using a linear ramp. For most landunits, this ramp generates all snow below :math:`0 ^{\circ} C`, all rain above :math:`2 ^{\circ} C`, and a mix of rain and snow for intermediate temperatures. For glaciers, the end points of the ramp are :math:`-2 ^{\circ} C` and :math:`0 ^{\circ} C`, respectively. Changes to the phase of precipitation are accompanied by a sensible heat flux (positive or negative) to conserve energy.

:sup:`3`\ There are 14 aerosol deposition rates required depending on species and affinity for bonding with water; 8 of these are dust deposition rates (dry and wet rates for 4 dust size bins, :math:`D_{dst,\, dry1},\, D_{dst,\, dry2},\, D_{dst,\, dry3},\, D_{dst,\, dry4}`, :math:`D_{dst,\, \, wet1},D_{dst,\, wet2},\, D_{dst,wet3},\, D_{dst,\, wet4}` ), 3 are black carbon deposition rates (dry and wet hydrophilic and dry hydrophobic rates, :math:`D_{bc,\, dryhphil},\, D_{bc,\, wethphil},\, D_{bc,\, dryhphob}` ), and 3 are organic carbon deposition rates (dry and wet hydrophilic and dry hydrophobic rates, :math:`D_{oc,\, dryhphil},\, D_{oc,\, wethphil},\, D_{oc,\, dryhphob}` ). These fluxes are computed interactively by the atmospheric model (when prognostic aerosol representation is active) or are prescribed from a time-varying (annual cycle or transient), globally-gridded deposition file defined in the namelist (see the CLM4.5 User's Guide). Aerosol deposition rates were calculated in a transient 1850-2009 CAM simulation (at a resolution of 1.9x2.5x26L) with interactive chemistry (troposphere and stratosphere) driven by CCSM3 20\ :sup:`th` century sea-surface temperatures and emissions (:ref:`Lamarque et al. 2010<Lamarqueetal2010>`) for short-lived gases and aerosols; observed concentrations were specified for methane, N\ :sub:`2`\ O, the ozone-depleting substances (CFCs),and CO\ :sub:`2`. The fluxes are used by the snow-related parameterizations (Chapters :numref:`rst_Surface Albedos` and :numref:`rst_Snow Hydrology`).

:sup:`4`\ The nitrogen deposition rate is required by the biogeochemistry model when active and represents the total deposition of mineral nitrogen onto the land surface, combining deposition of NO\ :sub:`y` and NH\ :sub:`x`. The rate is supplied either as a time-invariant spatially-varying annual mean rate or time-varying for a transient simulation. Nitrogen deposition rates were calculated from the same CAM chemistry simulation that generated the aerosol deposition rates.

:sup:`5`\ Climatological 3-hourly lightning frequency at :math:`\sim`\ 1.8° resolution is provided, which was calculated via bilinear interpolation from 1995-2011 NASA LIS/OTD grid product v2.2 (http://ghrc.msfc.nasa.gov) 2-hourly, 2.5° lightning frequency data. In future versions of the model, lightning data may be obtained directly from the atmosphere model.

Density of air (:math:`\rho _{atm}` ) (kg m\ :sup:`-3`) is also required but is calculated directly from :math:`\rho _{atm} =\frac{P_{atm} -0.378e_{atm} }{R_{da} T_{atm} }` where :math:`P_{atm}` is atmospheric pressure (Pa), :math:`e_{atm}` is atmospheric vapor pressure (Pa), :math:`R_{da}` is the gas constant for dry air (J kg\ :sup:`-1` K\ :sup:`-1`) (:numref:`Table Physical constants`), and :math:`T_{atm}` is the atmospheric temperature (K). The atmospheric vapor pressure :math:`e_{atm}` is derived from atmospheric specific humidity :math:`q_{atm}` (kg kg\ :sup:`-1`) as :math:`e_{atm} =\frac{q_{atm} P_{atm} }{0.622+0.378q_{atm} }`.

The O\ :sub:`2` partial pressure (Pa) is required but is calculated from molar ratio and the atmospheric pressure :math:`P_{atm}` as :math:`o_{i} =0.209P_{atm}`.

.. _Table Land model output to atmospheric model:

.. table:: Land model output to atmospheric model

 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Field                                 | Variable name                                  | units                                                        |
 +=======================================+================================================+==============================================================+
 | :sup:`1`\ Latent heat flux            | :math:`\lambda _{vap} E_{v} +\lambda E_{g}`    | W m\ :sup:`-2`                                               |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Sensible heat flux                    | :math:`H_{v} +H_{g}`                           | W m\ :sup:`-2`                                               |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Water vapor flux                      | :math:`E_{v} +E_{g}`                           | mm s\ :sup:`-1`                                              |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Zonal momentum flux                   | :math:`\tau _{x}`                              | kg m\ :sup:`-1` s\ :sup:`-2`                                 |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Meridional momentum flux              | :math:`\tau _{y}`                              | kg m\ :sup:`-1` s\ :sup:`-2`                                 |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Emitted longwave radiation            | :math:`L\, \uparrow`                           | W m\ :sup:`-2`                                               |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Direct beam visible albedo            | :math:`I\, \uparrow _{vis}^{\mu }`             | -                                                            |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Direct beam near-infrared albedo      | :math:`I\, \uparrow _{nir}^{\mu }`             | -                                                            |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Diffuse visible albedo                | :math:`I\, \uparrow _{vis}`                    | -                                                            |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Diffuse near-infrared albedo          | :math:`I\, \uparrow _{nir}`                    | -                                                            |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Absorbed solar radiation              | :math:`\vec{S}`                                | W m\ :sup:`-2`                                               |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Radiative temperature                 | :math:`T_{rad}`                                | K                                                            |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Temperature at 2 meter height         | :math:`T_{2m}`                                 | K                                                            |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Specific humidity at 2 meter height   | :math:`q_{2m}`                                 | kg kg\ :sup:`-1`                                             |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Wind speed at 10 meter height         | :math:`u_{10m}`                                | m s\ :sup:`-1`                                               |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Snow water equivalent                 | :math:`W_{sno}`                                | m                                                            |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Aerodynamic resistance                | :math:`r_{am}`                                 | s m\ :sup:`-1`                                               |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Friction velocity                     | :math:`u_{*}`                                  | m s\ :sup:`-1`                                               |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | :sup:`2`\ Dust flux                   | :math:`F_{j}`                                  | kg m\ :sup:`-2` s\ :sup:`-1`                                 |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+
 | Net ecosystem exchange                | NEE                                            | kgCO\ :sub:`2` m\ :sup:`-2` s\ :sup:`-1`                     |
 +---------------------------------------+------------------------------------------------+--------------------------------------------------------------+

:sup:`1`\ :math:`\lambda _{vap}` is the latent heat of vaporization (J kg\ :sup:`-1`) (:numref:`Table Physical constants`) and :math:`\lambda` is either the latent heat of vaporization :math:`\lambda _{vap}` or latent heat of sublimation :math:`\lambda _{sub}` (J kg\ :sup:`-1`) (:numref:`Table Physical constants`) depending on the liquid water and ice content of the top snow/soil layer (section 5.4).

:sup:`2`\ There are :math:`j=1,\ldots,4` dust transport bins.

.. _Initialization:

Initialization
^^^^^^^^^^^^^^^^^^^^

Initialization of the land model (i.e., providing the model with initial temperature and moisture states) depends on the type of run (startup or restart) (see the CLM4.5 User's Guide). A startup run starts the model from either initial conditions that are set internally in the Fortran code (referred to as arbitrary initial conditions) or from an initial conditions dataset that enables the model to start from a spun up state (i.e., where the land is in equilibrium with the simulated climate). In restart runs, the model is continued from a previous simulation and initialized from a restart file that ensures that the output is bit-for-bit the same as if the previous simulation had not stopped. The fields that are required from the restart or initial conditions files can be obtained by examining the code. Arbitrary initial conditions are specified as follows.

Soil points are initialized with surface ground temperature :math:`T_{g}` and soil layer temperature :math:`T_{i}`, for :math:`i=1,\ldots,N_{levgrnd}`, of 274 K, vegetation temperature :math:`T_{v}` of 283 K, no snow or canopy water (:math:`W_{sno} =0`, :math:`W_{can} =0`), and volumetric soil water content :math:`\theta _{i} =0.15` mm\ :sup:`3` mm\ :sup:`-3` for layers :math:`i=1,\ldots,N_{levsoi}` and :math:`\theta _{i} =0.0` mm\ :sup:`3` mm\ :sup:`-3` for layers :math:`i=N_{levsoi} +1,\ldots,N_{levgrnd}`. placeLake temperatures (:math:`T_{g}` and :math:`T_{i}` ) are initialized at 277 K and :math:`W_{sno} =0`.

Glacier temperatures (:math:`T_{g} =T_{snl+1}` and :math:`T_{i}` for :math:`i=snl+1,\ldots,N_{levgrnd}` where :math:`snl` is the negative of the number of snow layers, i.e., :math:`snl` ranges from –5 to 0) are initialized to 250 K with a snow water equivalent :math:`W_{sno} =1000` mm, snow depth :math:`z_{sno} =\frac{W_{sno} }{\rho _{sno} }` (m) where :math:`\rho _{sno} =250` kg m\ :sup:`-3` is an initial estimate for the bulk density of snow, and :math:`\theta _{i}` \ =1.0 for :math:`i=1,\ldots,N_{levgrnd}`. The snow layer structure (e.g., number of snow layers :math:`snl` and layer thickness) is initialized based on the snow depth (section 6.1). The snow liquid water and ice contents (kg m\ :sup:`-2`) are initialized as :math:`w_{liq,\, i} =0` and :math:`w_{ice,\, i} =\Delta z_{i} \rho _{sno}`, respectively, where :math:`i=snl+1,\ldots,0` are the snow layers, and :math:`\Delta z_{i}` is the thickness of snow layer :math:`i` (m). The soil liquid water and ice contents are initialized as :math:`w_{liq,\, i} =0` and :math:`w_{ice,\, i} =\Delta z_{i} \rho _{ice} \theta _{i}` for :math:`T_{i} \le T_{f}`, and :math:`w_{liq,\, i} =\Delta z_{i} \rho _{liq} \theta _{i}` and :math:`w_{ice,\, i} =0` for :math:`T_{i} >T_{f}`, where :math:`\rho _{ice}` and :math:`\rho _{liq}` are the densities of ice and liquid water (kg m\ :sup:`-3`) (:numref:`Table Physical constants`), and :math:`T_{f}` is the freezing temperature of water (K) (:numref:`Table Physical constants`). All vegetated and glacier land units are initialized with water stored in the unconfined aquifer and unsaturated soil :math:`W_{a} =4000` mm and water table depth :math:`z_{\nabla }` at five meters below the soil column.

.. _Surface Data:

Surface Data
^^^^^^^^^^^^^^^^^^

Required surface data for each land grid cell are listed in :numref:`Table Surface data required for CLM and their base spatial resolution` and include the glacier, lake, and urban fractions of the grid cell (vegetated and crop occupy the remainder), the fractional cover of each plant functional type (PFT), monthly leaf and stem area index and canopy top and bottom heights for each PFT, soil color, soil texture, soil organic matter density, maximum fractional saturated area, slope, elevation, biogenic volatile organic compounds (BVOCs) emissions factors, population density, gross domestic production, peat area fraction, and peak month of agricultural burning. Optional surface data include crop irrigation and managed crops. All fields are aggregated to the model's grid from high-resolution input datasets ( :numref:`Table Surface data required for CLM and their base spatial resolution`) that are obtained from a variety of sources described below.

.. _Table Surface data required for CLM and their base spatial resolution:

.. table:: Surface data required for CLM and their base spatial resolution

 +--------------------------------------------+---------------------------+
 | Surface Field                              | Resolution                |
 +============================================+===========================+
 | Percent glacier                            | 0.05°                     |
 +--------------------------------------------+---------------------------+
 | Percent lake and lake depth                | 0.05°                     |
 +--------------------------------------------+---------------------------+
 | Percent urban                              | 0.05°                     |
 +--------------------------------------------+---------------------------+
 | Percent plant functional types (PFTs)      | 0.05°                     |
 +--------------------------------------------+---------------------------+
 | Monthly leaf and stem area index           | 0.5°                      |
 +--------------------------------------------+---------------------------+
 | Canopy height (top, bottom)                | 0.5°                      |
 +--------------------------------------------+---------------------------+
 | Soil color                                 | 0.5°                      |
 +--------------------------------------------+---------------------------+
 | Percent sand, percent clay                 | 0.083°                    |
 +--------------------------------------------+---------------------------+
 | Soil organic matter density                | 0.083°                    |
 +--------------------------------------------+---------------------------+
 | Maximum fractional saturated area          | 0.125°                    |
 +--------------------------------------------+---------------------------+
 | Elevation                                  | 1km                       |
 +--------------------------------------------+---------------------------+
 | Slope                                      | 1km                       |
 +--------------------------------------------+---------------------------+
 | Biogenic Volatile Organic Compounds        | 0.5°                      |
 +--------------------------------------------+---------------------------+
 | Crop Irrigation                            | 0.083°                    |
 +--------------------------------------------+---------------------------+
 | Managed crops                              | 0.5°                      |
 +--------------------------------------------+---------------------------+
 | Population density                         | 0.5°                      |
 +--------------------------------------------+---------------------------+
 | Gross domestic production                  | 0.5°                      |
 +--------------------------------------------+---------------------------+
 | Peat area fraction                         | 0.5°                      |
 +--------------------------------------------+---------------------------+
 | Peak month of agricultural waste burning   | 0.5°                      |
 +--------------------------------------------+---------------------------+

At the base spatial resolution of 0.05°, the percentage of each PFT is defined with respect to the vegetated portion of the grid cell and the sum of the PFTs is 100%. The percent lake, glacier, and urban at their base resolution are specified with respect to the entire grid cell. The surface dataset creation routines re-adjust the PFT percentages to ensure that the sum of all land cover types in the grid cell sum to 100%. A minimum threshold of 0.1% of the grid cell by area is required for urban areas.

The percentage glacier mask was derived from vector data of global glacier and ice sheet spatial coverage. Vector data for glaciers (ice caps, icefields and mountain glaciers) were taken from the first globally complete glacier inventory, the Randolph Glacier Inventory version 1.0 (RGIv1.0: :ref:`Arendt et al. 2012 <Arendtetal2012>`). Vector data for the Greenland Ice Sheet were provided by Frank Paul and Tobias Bolch (University of Zurich: :ref:`Rastner et al. 2012 <Rastneretal2012>`). Antarctic Ice Sheet data were provided by Andrew Bliss (University of Alaska) and were extracted from the Scientific Committee on Antarctic Research (SCAR) Antarctic Digital Database version 5.0. Floating ice is only provided for the Antarctic and does not include the small area of Arctic ice shelves. High spatial resolution vector data were then processed to determine the area of glacier, ice sheet and floating ice within 30-second grid cells globally. The 30-second glacier, ice sheet and Antarctic ice shelf masks were subsequently draped over equivalent-resolution GLOBE topography (Global Land One-km Base Elevation Project, Hastings et al. 1999) to extract approximate ice-covered elevations of ice-covered regions. Grid cells flagged as land-ice in the mask but ocean in GLOBE (typically, around ice sheets at high latitudes) were designated land-ice with an elevation of 0 meters. Finally, the high-resolution mask/topography datasets were aggregated and processed into three 3-minute datasets: 3-minute fractional areal land ice coverage (including both glaciers and ice sheets); 3-minute distributions of areal glacier fractional coverage by elevation and areal ice sheet fractional coverage by elevation. Ice fractions were binned at 100 meter intervals, with bin edges defined from 0 to 6000 meters (plus one top bin encompassing all remaining high-elevation ice, primarily in the Himalaya). These distributions by elevation are used to divide each glacier land unit into columns based on elevation class.

When running with the CISM ice sheet model, CISM dictates glacier areas and elevations in its domain, overriding the values specified by CLM's datasets. In typical CLM5 configurations, this means that CISM dictates glacier areas and elevations over Greenland.

Percent lake and lake depth are area-averaged from the 90-second resolution data of :ref:`Kourzeneva (2009, 2010) <Kourzeneva2009>` to the 0.05° resolution using the MODIS land-mask. Percent urban is derived from LandScan 2004, a population density dataset derived from census data, nighttime lights satellite observations, road proximity and slope (:ref:`Dobson et al. 2000 <Dobsonetal2000>`) as described by :ref:`Jackson et al. (2010) <Jacksonetal2010>` at 1km resolution and aggregated to 0.05°. A number of urban radiative, thermal, and morphological fields are also required and are obtained from :ref:`Jackson et al. (2010) <Jacksonetal2010>`. Their description can be found in Table 3 of the Community Land Model Urban (CLMU) technical note (:ref:`Oleson et al. 2010b <Olesonetal2010b>`).

Percent PFTs are derived from MODIS satellite data as described in :ref:`Lawrence and Chase (2007) <LawrenceChase2007>` (section 21.3.3). Prescribed PFT leaf area index is derived from the MODIS satellite data of :ref:`Myneni et al. (2002) <Mynenietal2002>` using the de-aggregation methods described in :ref:`Lawrence and Chase (2007) <LawrenceChase2007>` (section 2.2.3). Prescribed PFT stem area index is derived from PFT leaf area index phenology combined with the methods of :ref:`Zeng et al. (2002) <Zengetal2002>`. Prescribed canopy top and bottom heights are from :ref:`Bonan (1996) <Bonan1996>` as described in :ref:`Bonan et al. (2002b) <Bonanetal2002b>`. If the biogeochemistry model is active, it supplies the leaf and stem area index and canopy top and bottom heights dynamically, and the prescribed values are ignored.

Soil color determines dry and saturated soil albedo (section :numref:`Ground Albedos`). Soil colors are from :ref:`Lawrence and Chase (2007) <LawrenceChase2007>`.

The soil texture and organic matter content determine soil thermal and hydrologic properties (sections 6.3 and 7.4.1). The International Geosphere-Biosphere Programme (IGBP) soil dataset (Global Soil Data Task 2000) of 4931 soil mapping units and their sand and clay content for each soil layer were used to create a mineral soil texture dataset :ref:`(Bonan et al. 2002b) <Bonanetal2002b>`. Soil organic matter data is merged from two sources. The majority of the globe is from ISRIC-WISE (:ref:`Batjes, 2006 <Batjes2006>`). The high latitudes come from the 0.25° version of the Northern Circumpolar Soil Carbon Database (:ref:`Hugelius et al. 2012 <Hugeliusetal2012>`). Both datasets report carbon down to 1m depth. Carbon is partitioned across the top seven CLM4 layers (:math:`\sim`\ 1m depth) as in :ref:`Lawrence and Slater (2008) <LawrenceSlater2008>`.

The maximum fractional saturated area (:math:`f_{\max }` ) is used in determining surface runoff and infiltration (section 7.3). Maximum fractional saturated area at 0.125° resolution is calculated from 1-km compound topographic indices (CTIs) based on the USGS HYDRO1K dataset (:ref:`Verdin and Greenlee 1996 <VerdinGreenlee1996>`) following the algorithm in :ref:`Niu et al. (2005) <Niuetal2005>`. :math:`f_{\max }` is the ratio between the number of 1-km pixels with CTIs equal to or larger than the mean CTI and the total number of pixels in a 0.125° grid cell. See section 7.3.1 and :ref:`Li et al. (2013b) <Lietal2013b>` for further details. Slope and elevation are also obtained from the USGS HYDRO1K 1-km dataset (:ref:`Verdin and Greenlee 1996 <VerdinGreenlee1996>`). Slope is used in the surface water parameterization (section :numref:`Surface Water Storage`), and elevation is used to calculate the grid cell standard deviation of topography for the snow cover fraction parameterization (section :numref:`Snow Covered Area Fraction`).

Biogenic Volatile Organic Compounds emissions factors are from the Model of Emissions of Gases and Aerosols from Nature version 2.1 (MEGAN2.1; :ref:`Guenther et al. 2012 <Guentheretal2012>`).

The default list of PFTs includes an unmanaged crop treated as a second C3 grass (:numref:`Table Plant functional types`). The unmanaged crop has grid cell fractional cover assigned from MODIS satellite data (:ref:`Lawrence and Chase (2007) <LawrenceChase2007>`). A managed crop option uses grid cell fractional cover from the present-day crop dataset of :ref:`Ramankutty and Foley (1998) <RamankuttyFoley1998>` (CLM4CNcrop). Managed crops are assigned in the proportions given by :ref:`Ramankutty and Foley (1998) <RamankuttyFoley1998>` without exceeding the area previously assigned to the unmanaged crop. The unmanaged crop continues to occupy any of its original area that remains and continues to be handled just by the CN part of CLM4CNcrop. The managed crop types (corn, soybean, and temperate cereals) were chosen based on the availability of corresponding algorithms in AgroIBIS (:ref:`Kucharik et al. 2000 <Kuchariketal2000>`; :ref:`Kucharik and Brye 2003 <KucharikBrye2003>`). Temperate cereals include wheat, barley, and rye here. All temperate cereals are treated as summer crops (like spring wheat, for example) at this time. Winter cereals (such as winter wheat) may be introduced in a future version of the model.

To allow crops to coexist with natural vegetation in a grid cell and be treated by separate models (i.e., CLM4.5BGCcrop versus the Dynamic Vegetation version (CLM4.5BGCDV)), we separate the vegetated land unit into a naturally vegetated land unit and a human managed land unit. PFTs in the naturally vegetated land unit share one soil column and compete for water (default CLM setting). PFTs in the human managed land unit do not share soil columns and thus permit for differences in land management between crops.

CLM includes the option to irrigate cropland areas that are equipped for irrigation. The application of irrigation responds dynamically to climate (see Chapter :numref:`rst_Crops and Irrigation`). In CLM, irrigation is implemented for the C3 generic crop only. When irrigation is enabled, the cropland area of each grid cell is divided into an irrigated and unirrigated fraction according to a dataset of areas equipped for irrigation (:ref:`Siebert et al. (2005) <Siebertetal2005>`). The area of irrigated cropland in each grid cell is given by the smaller of the grid cell's total cropland area, according to the default CLM4 dataset, and the grid cell's area equipped for irrigation. The remainder of the grid cell's cropland area (if any) is then assigned to unirrigated cropland. Irrigated and unirrigated crops are placed on separate soil columns, so that irrigation is only applied to the soil beneath irrigated crops.

Several input datasets are required for the fire model (:ref:`Li et al. 2013a <Lietal2013a>`) including population density, gross domestic production, peat area fraction, and peak month of agricultural waste burning. Population density at 0.5° resolution for 1850-2100 combines 5-min resolution decadal population density data for 1850–1980 from the Database of the Global Environment version 3.1 (HYDEv3.1) with 0.5° resolution population density data for 1990, 1995, 2000, and 2005 from the Gridded Population of the World version 3 dataset (GPWv3) (CIESIN, 2005). Gross Domestic Production (GDP) per capita in 2000 at 0.5° is from :ref:`Van Vuuren et al. (2006) <VanVuurenetal2006>`, which is the base-year GDP data for IPCC-SRES and derived from country-level World Bank's World Development Indicators (WDI) measured in constant 1995 US$ (:ref:`World Bank, 2004 <WorldBank2004>`) and the UN Statistics Database (:ref:`UNSTAT, 2005 <UNSTAT2005>`). The peatland area fraction at 0.5° resolution is derived from three vector datasets: peatland data in Indonesia and Malaysian Borneo (:ref:`Olson et al. 2001 <Olsonetal2001>`); peatland data in Canada (:ref:`Tarnocai et al. 2011 <Tarnocaietal2011>`); and bog, fen and mire data in boreal regions (north of 45°N) outside Canada provided by the Global Lakes and Wetlands Database (GLWD) (:ref:`Lehner and Döll, 2004 <LehnerDoll2004>`). The climatological peak month for agricultural waste burning is from :ref:`van der Werf et al. (2010) <vanderWerfetal2010>`.

.. _Adjustable Parameters and Physical Constants:

Adjustable Parameters and Physical Constants
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Values of certain adjustable parameters inherent in the biogeophysical or biogeochemical parameterizations have either been obtained from the literature or calibrated based on comparisons with observations. These are described in the text. Physical constants, generally shared by all of the components in the coupled modeling system, are presented in :numref:`Table Physical constants`.

.. _Table Physical constants:

.. csv-table:: Physical constants
   :header: "description", "name", "value", "units"
   :widths: 40, 20, 20, 20

   "Pi", :math:`\pi`, 3.14159265358979323846, "\-"
   "Acceleration of gravity", :math:`g`, 9.80616,  m s\ :sup:`-2`
   "Standard pressure", :math:`P_{std}`, 101325, "Pa"
   "Stefan-Boltzmann constant", :math:`\sigma`, 5.67 :math:`\times 10^{-8}`, W m :sup:`-2` K :math:`{}^{-4}`
   "Boltzmann constant", :math:`\kappa`, 1.38065 :math:`\times 10^{-23}`, J K :sup:`-1` molecule :sup:`-1`
   "Avogadro's number", :math:`N_{A}`, 6.02214 :math:`\times 10^{26}`, molecule kmol\ :sup:`-1`
   "Universal gas constant", :math:`R_{gas}`, :math:`N_{A} \kappa`, J K :sup:`-1` kmol :sup:`-1`
   "Molecular weight of dry air", :math:`MW_{da}`, 28.966, kg kmol :sup:`-1`
   "Dry air gas constant", :math:`R_{da}`, :math:`{R_{gas} \mathord{\left/ {\vphantom {R_{gas}  MW_{da} }} \right.} MW_{da} }`, J K :sup:`-1` kg :sup:`-1`
   "Molecular weight of water vapor", :math:`MW_{wv}`, 18.016, kg kmol :sup:`-1`
   "Water vapor gas constant", :math:`R_{wv}`, :math:`{R_{gas} \mathord{\left/ {\vphantom {R_{gas}  MW_{wv} }} \right.} MW_{wv} }`, J K :sup:`-1` kg :sup:`-1`
   "Von Karman constant", :math:`k`, 0.4, "\-"
   "Freezing temperature of fresh water", :math:`T_{f}`, 273.15, K
   "Density of liquid water", :math:`\rho _{liq}`, 1000, kg m :sup:`-3`
   "Density of ice", :math:`\rho _{ice}`, 917, kg m :sup:`-3`
   "Specific heat capacity of dry air", :math:`C_{p}`, 1.00464 :math:`\times 10^{3}`, J kg :sup:`-1` K :sup:`-1`
   "Specific heat capacity of water", :math:`C_{liq}`, 4.188 :math:`\times 10^{3}`, J kg :sup:`-1` K :sup:`-1`
   "Specific heat capacity of ice", :math:`C_{ice}`, 2.11727 :math:`\times 10^{3}`, J kg :sup:`-1` K :sup:`-1`
   "Latent heat of vaporization", :math:`\lambda _{vap}`, 2.501 :math:`\times 10^{6}`, J kg :sup:`-1`
   "Latent heat of fusion", :math:`L_{f}`, 3.337 :math:`\times 10^{5}`, J kg :sup:`-1`
   "Latent heat of sublimation", :math:`\lambda _{sub}`, :math:`\lambda _{vap} +L_{f}`, J kg :sup:`-1`
   :sup:`1` "Thermal conductivity of water", :math:`\lambda _{liq}`, 0.57, W m :sup:`-1` K :sup:`-1`
   :sup:`1` "Thermal conductivity of ice", :math:`\lambda _{ice}`, 2.29, W m :sup:`-1` K :sup:`-1`
   :sup:`1` "Thermal conductivity of air", :math:`\lambda _{air}`, 0.023 W m :sup:`-1` K :sup:`-1`
   "Radius of the earth", :math:`R_{e}`, 6.37122, :math:`\times 10^{6}` m

:sup:`1`\ Not shared by other components of the coupled modeling system.

