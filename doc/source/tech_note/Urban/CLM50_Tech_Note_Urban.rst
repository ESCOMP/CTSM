.. _rst_Urban Model (CLMU):

Urban Model (CLMU)
======================

At the global scale, and at the coarse spatial resolution of current
climate models, urbanization has negligible impact on climate. However,
the urban parameterization (CLMU; :ref:`Oleson et al. (2008b) <Olesonetal2008b>`;
:ref:`Oleson et al. (2008c) <Olesonetal2008c>`) allows
simulation of the urban environment within a climate model, and
particularly the temperature where people live. As such, the urban model
allows scientific study of how climate change affects the urban heat
island and possible urban planning and design strategies to mitigate
warming (e.g., white roofs).

Urban areas in CLM are represented by up to three urban landunits per
gridcell according to density class. The urban landunit is based on the
“urban canyon” concept of :ref:`Oke (1987) <Oke1987>` in which 
the canyon geometry is
described by building height (:math:`H`) and street width (:math:`W`)
(:numref:`Figure schematic representation of the urban landunit`). The canyon system 
consists of roofs, walls, and canyon
floor. Walls are further divided into shaded and sunlit components. The
canyon floor is divided into pervious (e.g., to represent residential
lawns, parks) and impervious (e.g., to represent roads, parking lots,
sidewalks) fractions. Vegetation is not explicitly modeled for the
pervious fraction; instead evaporation is parameterized by a simplified
bulk scheme.

Each of the five urban surfaces is treated as a column within the
landunit (:numref:`Figure schematic representation of the urban landunit`). 
Radiation parameterizations account for trapping
of solar and longwave radiation inside the canyon. Momentum fluxes are
determined for the urban landunit using a roughness length and
displacement height appropriate for the urban canyon and stability
formulations from CLM. A one-dimensional heat conduction equation is
solved numerically for a multiple-layer (:math:`N_{levurb} =10`) column
to determine conduction fluxes into and out of canyon surfaces. 

A new building energy model has been developed for CLM5.0.  It accounts
for the conduction of heat through interior surfaces (roof, sunlit and
shaded walls, and floors), convection (sensible heat exchange) between 
interior surfaces and building air, longwave radiation exchange between
interior surfaces, and ventilation (natural infiltration and exfiltration).
Idealized HAC systems are assumed where the system capacity is infinite and
the system supplies the amount of energy needed to keep the indoor air 
temperature (:math:`T_{iB}`) within maximum and minimum emperatures
(:math:`T_{iB,\, \max } ,\, T_{iB,\, \min }` ), thus explicitly
resolving space heating and air conditioning fluxes. Anthropogenic sources
of waste heat (:math:`Q_{H,\, waste}` ) from HAC that account for inefficiencies
in the heating and air conditioning equipment and from energy lost in the 
conversion of primary energy sources to end use energy are derived from 
:ref:`Sivak (2013) <Sivak2013>`.  These sources of waste heat are incorporated 
as modifications to the canyon energy budget.

Turbulent [sensible heat (:math:`Q_{H,\, u}` ) and
latent heat (:math:`Q_{E,\, u}` )] and storage (:math:`Q_{S,\, u}` )
heat fluxes and surface (:math:`T_{u,\, s}` ) and internal
(:math:`T_{u,\, i=1,\, N_{levgrnd} }` ) temperatures are determined for
each urban surface :math:`u`. Hydrology on the roof and canyon floor is
simulated and walls are hydrologically inactive. A snowpack can form on
the active surfaces. A certain amount of liquid water is allowed to pond
on these surfaces which supports evaporation. Water in excess of the
maximum ponding depth runs off
(:math:`R_{roof} ,\, R_{imprvrd} ,\, R_{prvrd}` ).

The heat and moisture fluxes from each surface interact with each other
through a bulk air mass that represents air in the urban canopy layer
for which specific humidity (:math:`q_{ac}` ) and temperature
(:math:`T_{ac}` ) are prognosed (:numref:`Figure schematic of urban and atmospheric model coupling`).
The air temperature can
be compared with that from surrounding vegetated/soil (rural) surfaces
in the model to ascertain heat island characteristics. As with other
landunits, the CLMU is forced either with output from a host atmospheric
model (e.g., the Community Atmosphere Model (CAM)) or
observed forcing (e.g., reanalysis or field observations). The urban
model produces sensible, latent heat, and momentum fluxes, emitted
longwave, and reflected solar radiation, which are area-averaged with
fluxes from non-urban “landunits” (e.g., vegetation, lakes) to supply
grid cell averaged fluxes to the atmospheric model.

Present day global urban extent and urban properties were developed by
:ref:`Jackson et al. (2010) <Jacksonetal2010>`. Urban extent, defined for four classes [tall
building district (TBD), and high, medium, and low density (HD, MD,
LD)], was derived from LandScan 2004, a population density dataset
derived from census data, nighttime lights satellite observations, road
proximity, and slope (:ref:`Dobson et al. 2000 <Dobsonetal2000>`). The urban extent data for
TBD, HD, and MD classes are aggregated from the original 1 km resolution
to both a 0.05\ :sup:`o` by 0.05\ :sup:`o` global grid
for high-resolution studies or a 0.5\ :sup:`o` by
0.5\ :sup:`o` grid. For the current implementation, the LD class
is not used because it is highly rural and better modeled as a
vegetated/soil surface. Although the TBD, HD, and MD classes are
represented as individual urban landunits, urban model history output is
currently a weighted average of the output for individual classes.

For each of 33 distinct regions across the globe, thermal (e.g., heat
capacity and thermal conductivity), radiative (e.g., albedo and
emissivity) and morphological (e.g., height to width ratio, roof
fraction, average building height, and pervious fraction of the canyon
floor) properties are provided for each of the density classes. Building
interior minimum and maximum temperatures are prescribed based on
climate and socioeconomic considerations. The surface dataset creation
routines (see CLM5.0 User’s Guide) aggregate the data to the desired
resolution.

An optional urban properties dataset, including a tool that allows for generating future
urban development scenarios is also available (:ref:`Oleson and Feddema (2018) <OlesonFeddema2018>`).
This will become the default dataset in future model versions.
As described in :ref:`Oleson and Feddema (2018) <OlesonFeddema2018>` the urban properties dataset
in :ref:`Jackson et al. (2010) <Jacksonetal2010>` was modified with respect to wall and roof thermal
properties to correct for biases in heat transfer due to layer and building type averaging. 
Further changes to the dataset reflect the need for scenario development, thus allowing for 
the creation of hypothetical wall types, and the easier interchange of wall facets.
The new urban properties tool is available as part of the Toolbox for Human-Earth System 
Integration & Scaling (THESIS) tool set 
(http://www.cgd.ucar.edu/iam/projects/thesis/thesis-urbanproperties-tool.html; 
:ref:`Feddema and Kauffman (2016) <FeddemaKauffman2016>`). The driver script (urban_prop.csh) 
specifies three input csv files (by default, mat_prop.csv, 
lam_spec.csv, and city_spec.csv; (:numref:`Figure schematic of THESIS urban properties tool`)) 
that describe the morphological, radiative, and thermal properties of urban areas, and 
generates a global dataset at 0.05° latitude by longitude in NetCDF format (urban_properties_data.05deg.nc).
A standalone NCL routine (gen_data_clm.ncl) can be run separately after the mksurfdata_map tool creates 
the CLM surface dataset.  This creates a supplementary streams file of setpoints for the maximum 
interior building temperature at yearly time resolution.

.. Figure 12.1. Schematic representation of the urban land unit

.. _Figure schematic representation of the urban landunit:

.. figure:: image1.png

 Schematic representation of the urban land unit. See the text for description of notation. Incident, reflected, and net solar and longwave radiation are calculated for each individual surface but are not shown for clarity.

.. Figure 12.2. Schematic of urban and atmospheric model coupling

.. _Figure schematic of urban and atmospheric model coupling:

.. Figure:: image2.png

 Schematic of urban and atmospheric model coupling.  The urban model is forced by the atmospheric model wind (:math:`u_{atm}` ), temperature (:math:`T_{atm}` ), specific humidity (:math:`q_{atm}` ), precipitation (:math:`P_{atm}` ), solar (:math:`S_{atm} \, \downarrow` ) and longwave (:math:`L_{atm} \, \downarrow` ) radiation at reference height :math:`z'_{atm}`  (section :numref:`Atmospheric Coupling`). Fluxes from the urban landunit to the atmosphere are turbulent sensible (:math:`H`) and latent heat (:math:`\lambda E`), momentum (:math:`\tau` ), albedo (:math:`I\uparrow` ), emitted longwave (:math:`L\uparrow` ), and absorbed shortwave (:math:`\vec{S}`) radiation. Air temperature (:math:`T_{ac}` ), specific humidity (:math:`q_{ac}` ), and wind speed (:math:`u_{c}` ) within the urban canopy layer are diagnosed by the urban model. :math:`H` is the average building height.

.. Figure 12.3. Schematic of THESIS urban properties tool

.. _Figure schematic of THESIS urban properties tool:

.. Figure:: image3.png

 Schematic of THESIS urban properties tool.  Executable scripts are in orange, input files are blue, and output files are green.  Items within the black box outline are either read in as input, executed, or output by the driver script (urban_prop.csh).  


The urban model that was first released as a component of CLM4.0 is separately
described in the urban technical note (:ref:`Oleson et al. (2010b) <Olesonetal2010b>`).
The main changes in the urban model from CLM4.0 to CLM4.5 were 1)
an expansion of the single urban landunit to up to three landunits per
grid cell stratified by urban density types, 2) the number of urban
layers for roofs and walls was no longer constrained to be equal to the
number of ground layers, 3) space heating and air conditioning wasteheat
factors were set to zero by default so that the user could customize
these factors for their own application, 4) the elevation threshold used
to eliminate urban areas in the surface dataset creation routines was
increased from 2200 meters to 2600 meters, 5) hydrologic and thermal
calculations for the pervious road followed CLM4.5 parameterizations.

The main changes in the urban model from CLM4.5 to CLM5.0 are 1) a more 
sophisticated and realistic building space heating and air conditioning 
submodel that prognoses interior building air temperature and includes more
realistic space heating and air conditioning wasteheat factors (see above), 2) the maximum
building temperature (which determines air conditioning demand) is now read in
from a namelist-defined file which allows for dynamic control of this input 
variable.  The maximum building temperatures that are defined in 
:ref:`Jackson et al. (2010) <Jacksonetal2010>` are implemented in year 1950 (thus
air conditioning is off in prior years) and air conditioning is turned off in year
2100 (because the buildings are not suitable for air conditioning in some extreme
global warming scenarios), 3) an optional updated urban properties dataset and new 
scenario tool.  These features are described in more detail in :ref:`Oleson and Feddema (2018) <OlesonFeddema2018>`. 
In addition, a module of heat stress indices calculated online
in the model that can be used to assess human thermal comfort for rural and urban
areas has been added.  This last development is described and evaluated by 
:ref:`Buzan et al. (2015) <Buzanetal2015>`.
