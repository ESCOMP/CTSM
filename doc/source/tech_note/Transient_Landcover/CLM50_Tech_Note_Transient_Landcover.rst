.. _rst_Transient Landcover Change:

Transient Land Use and Land Cover Change
========================================

CLM includes a treatment of mass and energy fluxes associated with prescribed temporal land use and land cover change (LULCC). The model uses an annual time series of the spatial distribution of the natural and crop land units of each grid cell, in combination with the distribution of PFTs and CFTs that exist in those land units. Additional land use is prescribed through annual crop-specific management of nitrogen fertilizer and irrigation (described further in :numref:`rst_Crops and Irrigation`), and through wood harvest on tree PFTs. For changes in the distributions of natural and crop vegetation, CLM diagnoses the change in area of the PFTs and CFTs on January 1 of each model year and then performs mass and energy balance accounting necessary to represent the expansion and contraction of the PFT and CFT areas. The biogeophysical impacts of LULCC are simulated through changes in surface properties which in turn impact the surface albedo, hydrology, and roughness which then impact fluxes of energy, moisture and momentum to the atmosphere under the altered properties. Additionally, changes in energy and moisture associated with changes in the natural and crop vegetation distribution are accounted for through small fluxes to the river and atmosphere. The biogeochemical impacts of LULCC are simulated through changes in CLM carbon pools and fluxes (see also Chapter :numref:`rst_CN Pools`).

CLM can also respond to changes in ice sheet areas and elevations when it is coupled to an evolving ice sheet model (in the CESM context, this is the Community Ice Sheet Model, CISM; see also Chapter :numref:`rst_Glaciers`). Conservation of water, energy, carbon and nitrogen is handled similarly for glacier-vegetation transitions as for natural vegetation-crop transitions.

.. _Transient land use and land cover data:

Annual Transient Land Use and Land Cover Data
---------------------------------------------

The changes in area over time associated with changes in natural and crop vegetation and the land use on that vegetation are prescribed through a forcing dataset, referred to here as the *landuse.timeseries* dataset. The *landuse.timeseries* dataset consists of an annual time series of global grids, where each annual time slice describes the fractional area occupied by all PFTs and CFTs along with the nitrogen fertilizer and irrigation fraction of each crop CFT, and the annual wood harvest applied to tree PFTs. Changes in area of PFTs and CFTs are performed annually on the first time step of January 1 of the year. Wood harvest for each PFT is also performed on the first time step of the year. Fertilizer application and irrigation for each CFT are performed at each model time step depending on rules from the crop model. Fertilizer application rates are set annually. The irrigation fraction is also set annually; irrigated crops are placed on separate columns from their unirrigated counterparts, so changes in irrigated fraction triggers the changes in subgrid areas discussed below (sections :numref:`Transient landcover reconciling changes in area` and :numref:`Transient landcover mass and energy conservation`).

As a special case, when the time dimension of the *landuse.timeseries* dataset starts at a later year than the current model time step, the first time slice from the *landuse.timeseries* dataset is used to represent the current time step PFT and CFT fractional area distributions. Similarly, when the time dimension of the *landuse.timeseries* dataset stops at an earlier year than the current model time step, the last time slice of the *landuse.timeseries* dataset is used. Thus, the simulation will have invariant representations of PFT and CFT distributions through time for the periods prior to and following the time duration of the *landuse.timeseries* dataset, with transient PFT and CFT distributions during the period covered by the *landuse.timeseries* dataset.

.. _Transient landcover reconciling changes in area:

Reconciling Changes in Area
---------------------------

In the first time step of January 1, changes in land unit weights can potentially come from two sources: Changes in the area of the crop land unit come from the *landuse.timeseries* dataset (section :numref:`Transient land use and land cover data`), and changes in the area of the glacier land unit come from the ice sheet model. The areas of other land units are then adjusted so that the total land unit area remains 100%.

If the total land unit area of glaciers and crops has decreased, then the natural vegetated landunit is increased to fill in the abandoned land. If the total land unit area of glaciers and crops has increased, then other land unit areas are decreased in a specified order until the total is once again 100%. The order of decrease is: natural vegetation, crop, urban medium density, urban high density, urban tall building district, wetland, lake.

These rules have two important implications:

1. We always match CISM's glacier areas exactly, even if that means a disagreement with prescribed crop areas. This is needed for conservation when CISM is evolving in two-way-coupled mode.

2. For land units other than crop, glacier and natural vegetation, their areas can decrease (due to encroaching crops or glaciers), but can never increase. So, for example, if a grid cell starts as 5% lake, crops expand to fill the entire grid cell, then later crop area decreases, the lake area will not return: instead, the abandoned cropland will become entirely natural vegetation.

For all levels of the subgrid hierarchy (land unit, column and patch), we only track net changes in area, not gross transitions. So, for example, if part of a gridcell experiences an increase in glacier area while another part of that gridcell experiences an equal decrease in glacier area (in the same glacier elevation class), CLM acts as if there were no changes. As another example, consider a gridcell containing natural vegetation, crop and glacier. If there is a decrease in glacier area and an equal increase in crop area, CLM will assume that the crop expands into the old glacier area, and nothing happened to the natural vegetation area. A more realistic alternative would be that the crop expanded into natural vegetation, and natural vegetation expanded into glacier. The final areas will be correct in these cases, but the adjustments of carbon and nitrogen states (section :numref:`Transient landcover carbon and nitrogen conservation`) will be less accurate than what would be obtained with a full tracking of gross transitions.

.. _Transient landcover mass and energy conservation:

Mass and Energy Conservation
----------------------------

.. _Transient landcover water and energy conservation:

Water and Energy Conservation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When subgrid areas change, the water and energy states remain unchanged on a per-area basis. This can lead to changes in the total gridcell water and energy content.

For example, consider a gridcell with two columns: column 1 has a water mass of 1 kg m\ :sup:`-2` and column 2 has a water mass of 2 kg m\ :sup:`-2` for a given water state variable, where these are expressed per unit column area. If column 1 increases in area at the expense of column 2, then column 1 will still have a water mass of 1 kg m\ :sup:`-2`, but now expressed over the new column area. This results in a decrease in the total gridcell water content.

Water and energy are conserved by summing up the total water and energy content of each gridcell before and after a change in area. Differences in liquid and ice water content are balanced by liquid and ice runoff terms, which can be either positive or negative. (Negative runoff is effectively a withdrawal of water from the ocean.) Differences in energy content are balanced by a sensible heat flux term, which again can be either positive or negative. These balancing fluxes are spread evenly throughout the following year.

There is a special case when a given crop column type newly comes into existence - for example, when temperate corn first comes into existence in a gridcell. In this case, the column's below-ground temperature and water states are copied from the natural vegetated column in its gridcell, so that these state variables begin in a close-to-spun-up state. Other state variables (most of which spin up relatively quickly) begin at their cold start initialization values. This initialization is not necessary for the two other land unit types that currently can grow - natural vegetation and glacier: Those land unit types are always active, even when they have zero area on the gridcell, so their state variables will be spun up immediately when they come into existence. After this initialization, the conservation code described above takes effect.

.. _Transient landcover carbon and nitrogen conservation:

Carbon and Nitrogen Conservation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Because of the long timescales involved with below-ground carbon and nitrogen dynamics, it is more important that these state variables be adjusted properly when subgrid areas change. Carbon and nitrogen variables are adjusted with the following three-step process:

(1) Patch-level (i.e., vegetation) state variables are adjusted for any changes in patch areas; this may lead to fluxes into column-level (i.e., soil) state variables (2) Column-level (i.e., soil) state variables are updated based on the fluxes generated in (1)

(3) Column-level (i.e., soil) state variables are adjusted for any changes in column areas First, patch-level (i.e., vegetation) state variables are adjusted for any changes in patch areas. This includes changes in column or land unit areas, even if the relative proportions of each patch remain constant: the relevant quantities are the patch weights relative to the gridcell.

For a patch that decreases in area, the carbon and nitrogen density on the remaining patch area remains the same as before (i.e., expressed as g per m\ :sup:`2` patch area). Because the area has decreased, this represents a decrease in total carbon or nitrogen mass (i.e., expressed as g per m\ :sup:`2` gridcell area). The lost mass meets a variety of fates: some is immediately lost to the atmosphere, some is sent to product pools (which are lost to the atmosphere over longer time scales), and some is sent to litter pools.

For a patch that increases in area, the carbon and nitrogen density on the new patch area is decreased in order to conserve mass. This decrease is basically proportional to the relative increase in patch area. However, a small amount of seed carbon and nitrogen is added to the leaf and dead stem pools in the new patch area.

Next, column-level (i.e., soil) state variables are updated based on any fluxes to soil pools due to decreases in patch areas. This step is needed so that any lost vegetation carbon and nitrogen is conserved when column areas are changing.

Finally, column-level state variables are adjusted for any changes in column areas. Similarly to patches, for a column that decreases in area, the carbon and nitrogen density on the remaining column area remains the same as before (i.e., expressed as g per m\ :sup:`2` column area). This represents a decrease in total carbon or nitrogen mass on the gridcell, and this lost mass is tracked for each gridcell. After these mass losses are summed for all shrinking columns, they are distributed amongst the growing columns in order to conserve mass. Thus, a growing column's new carbon density will be a weighted sum of its original carbon density and the carbon densities of all shrinking columns in its gridcell.

This operation makes some simplifying assumptions. First, as described in section :numref:`Transient landcover reconciling changes in area`, we only track net area changes, not gross changes. Second, we assume that growing columns all grow proportionally into each of the shrinking columns.

Non-vegetated land units (e.g., glacier) do not typically track soil carbon and nitrogen. When columns from these land units initially shrink, they are assumed to contribute zero carbon and nitrogen. However, when they grow into previously-vegetated areas, they store any pre-existing soil carbon and nitrogen from the shrinking columns. This stored carbon and nitrogen will remain unchanged until the column later shrinks, at which point it will contribute to the carbon and nitrogen in the growing columns (exactly as would happen for a vegetated column).

In contrast to water and energy (section :numref:`Transient landcover water and energy conservation`), no special treatment is needed for carbon and nitrogen states in columns that newly come into existence. The state of a new column is derived from a weighted average of the states of shrinking columns. This behavior falls out from the above general rules.

Annual Transient Land Cover Dataset Development
----------------------------------------------------

This section describes the development of the *landuse.timeseries* dataset. Development of this dataset involves the translation of harmonized datasets of LULCC for the historical period and for the different Shared Socioeconomic Pathway (SSP) - Representative Concentration Pathway (RCP) scenarios. Additionally, LULCC time series are to be generated for the Last Millennium and the extension beyond 2100 experiments of CMIP6.

LUH2 Transient Land Use and Land Cover Change Dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To coordinate the processing and consistency of LULCC data between the historical period (1850-2015) and the six SSP-RCP (2016-2100) scenarios derived from Integrated Assessment Models (IAM), the University of Maryland and the University of New Hampshire research groups (Louise Chini, George Hurtt, Steve Frolking and Ritvik Sahajpal; luh.umd.edu) produced a new version of the Land Use Harmonized version 2 (LUH2) transient datasets for use with Earth System Model simulations. The new data sets are the product of the Land Use Model Intercomparison Project (LUMIP; https://cmip.ucar.edu/lumip) as part of the Coupled Model Intercomparison Project 6 (CMIP6). The historical component of the transient LULCC dataset has agriculture and urban land use based on HYDE 3.2 with wood harvest based on FAO, Landsat and other sources, for the period 850-2015. The SSP-RCP transient LULCC components (2015-2100) are referred to as the LUH2 Future Scenario datasets. The LULCC information is provided at 0.25 degree grid resolution and includes fractional grid cell coverage by the 12 land units of:

Primary Forest, Secondary Forest, Primary Non-Forest, Secondary Non-Forest,

Pasture, Rangeland, Urban,

C3 Annual Crop, C4 Annual Crop, C3 Perennial Crop, C4 Perennial Crop, and C3 Nitrogen Fixing Crop.

The new land unit format is an improvement on the CMIP5 LULCC datasets as they: provide Forest and Non Forest information in combination with Primary and Secondary land; differentiate between Pasture and Rangelands for grazing livestock; and specify annual details on the types of Crops grown and management practices applied in each grid cell. Like the CMIP5 LULCC datasets Primary vegetation represents the fractional area of a grid cell with vegetation undisturbed by human activities. Secondary vegetation represents vegetated areas that have recovered from some human disturbance; this could include re-vegetation of pasture and crop areas as well as primary vegetation areas that have been logged. In this manner the land units can change through deforestation from Forested to Non Forested land and in the opposite direction from Non Forested to Forested land through reforestation or afforestation without going through the Crop, Pasture or Rangeland states.

The LUH2 dataset provides a time series of land cover states as well as a transition matrices that describes the annual fraction of land that is transformed from one land unit category to another (e.g. Primary Forest to C3 Annual Crop, Pasture to C3 Perrenial Crop, etc.; Lawrence et al. 2016). Included in these transition matrices is the total conversion of one land cover type to another referred to as Gross LULCC. This value can be larger than the sum of the changes in the state of a land unit from one time period to the next known as the Net LULCC. This difference is possible as land unit changes can occur both from the land unit and to the land unit at the same time. An example of this difference occurs with shifting cultivation where Secondary Forest can be converted to C3 Annual Crop at the same time as C3 Annual Crop is abandoned to Secondary Forest.

The transition matrices also provide harmonized prescriptions of wood harvest both in area of the grid cell harvested and in the amount of biomass carbon harvested. The wood harvest biomass amount includes a 30% slash component inline with the CMIP5 LULCC data described in (Hurtt et al. 2011). The harvest area and carbon amounts are prescribed for the five classes of: Primary Forest, Primary Non-Forest, Secondary Mature Forest, Secondary Young Forest, and Secondary Non-Forest.

Additional land use management is prescribed on the Crop land units for nitrogen fertilization and irrigation equipped land. The fertilizer application and the the irrigation fraction is prescribed for each Crop land unit in a grid cell individually for each year of the time series. The wood harvest and crop management are both prescribed spatially on the same 0.25 degree grid as the land use class transitions.

Representing LUH2 Land Use and Land Cover Change in CLM5
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To represent the LUH2 transient LULCC dataset in CLM5, the annual fractional composition of the twelve land units specified in the dataset needs to be faithfully represented with a corresponding PFT and CFT mosaics of CLM. CLM5 represents the land surface as a hierarchy of sub-grid types: glacier; lake; urban; vegetated land; and crop land. The vegetated land is further divided into a mosaic of Plant Functional Types (PFTs), while the crop land is divided into a mosaic of Crop Functional Types (CFTs).

To support this translation task the CLM5 Land Use Data tool has been built that extends the methods described in Lawrence et al (2012) to include all the new functionality of CMIP6 and CLM5 LULCC. The tool translates each of the LUH2 land units for a given year into fractional PFT and CFT values based on the current day CLM5 data for the land unit in that grid cell. The current day land unit descriptions are generated from from 1km resolution MODIS, MIRCA2000, ICESAT, AVHRR, SRTM, and CRU climate data products combined with reference year LUH2 land unit data, usually set to 2005. Where the land unit does not exist in a grid cell for the current day, the land unit description is generated from nearest neighbors with an inverse distance weighted search algorithm.

The Land Use Data tool produces raw vegetation, crop, and management data files which are combined with other raw land surface data to produce the CLM5 initial surface dataset and the dynamic *landuse.timeseries* dataset with the CLM5 mksurfdata_esmf tool. The schematic of this entire process from LUH2 time series and high resolution current day data to the output of CLM5 surface datasets from the mksurfdata_esmf tool is shown in Figure 21.2.

The methodology for creating the CLM5 transient PFT and CFT dataset is based on four steps which are applied across all of the historical and future time series. The first step involves generating the current day descriptions of natural and managed vegetation PFTs at 1km resolution from the global source datasets, and the current day description of crop CFTs at the 10km resolution from the MIRCA 2000 datasets. The second step combines the current day (2005) LUH2 land units with the current day CLM5 PFT and CFT distributions to get CLM5 land unit descriptions in either PFTs or CFTs at the LUH2 resolution of 0.25 degrees. The third step involves combining the LUH2 land unit time series with the CLM5 PFT and CFT descriptions for that land unit to generate the CLM5 raw PFT and CFT time series in the *landuse.timeseries* file. At this point in the process management information in terms of fertilizer, irrigation and wood harvest are added to the CLM5 PFT and CFT data to complete the CLM5 raw PFT and CFT files. The final step is to combine these files with the other raw CLM5 surface data files in the mksurfdata_esmf tool.

.. _Figure Schematic of land cover change:

.. figure:: image1.png

 Schematic of land cover change impacts on CLM carbon pools and fluxes.

.. _Figure Schematic of translation of annual LUH2 land units:

.. figure:: image2.png

 Schematic of translation of annual LUH2 land units to CLM5 plant and crop functional types.

.. _Figure Workflow of CLM5 Land Use Data Tool and mksurfdata_esmf Tool:

.. figure:: image3.png

 Workflow of CLM5 Land Use Data Tool and mksurfdata_esmf Tool
