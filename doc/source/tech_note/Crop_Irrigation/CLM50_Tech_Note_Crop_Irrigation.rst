.. _rst_Crops and Irrigation:

Crops and Irrigation
====================

.. _Summary of CLM5.0 updates relative to the CLM4.5:

Summary of CLM5.0 updates relative to the CLM4.5
------------------------------------------------

We describe here the complete crop and irrigation parameterizations that appear in CLM5.0. Corresponding information for CLM4.5 appeared in the CLM4.5 Technical Note (:ref:`Oleson et al. 2013 <Olesonetal2013>`).

CLM5.0 includes the following new updates to the CROP option, where CROP refers to the interactive crop management model and is included as an option with the BGC configuration:

- New crop functional types

- All crop areas are actively managed

- Fertilization rates updated based on crop type and geographic region

- New Irrigation triggers

- Phenological triggers vary by latitude for some crop types

- Ability to simulate transient crop management

- Adjustments to allocation and phenological parameters

- Crops reaching their maximum LAI triggers the grain fill phase

- Grain C and N pools are included in a 1-year product pool

- C for annual crop seeding comes from the grain C pool

- Initial seed C for planting is increased from 1 to 3 g C/m^2

These updates appear in detail in the sections below. Many also appear in :ref:`Levis et al. (2016) <Levisetal2016>`.

Available new features since the CLM5 release
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Addition of bioenergy crops
- Ability to customize crop calendars (sowing windows/dates, maturity requirements) using stream files
- Cropland soil tillage
- Crop residue removal

.. _The crop model:

The crop model: cash and bioenergy crops
----------------------------------------

Introduction
^^^^^^^^^^^^

Groups developing Earth System Models generally account for the human footprint on the landscape in simulations of historical and future climates. Traditionally we have represented this footprint with natural vegetation types and particularly grasses because they resemble many common crops. Most modeling efforts have not incorporated more explicit representations of land management such as crop type, planting, harvesting, tillage, fertilization, and irrigation, because global scale datasets of these factors have lagged behind vegetation mapping. As this begins to change, we increasingly find models that will simulate the biogeophysical and biogeochemical effects not only of natural but also human-managed land cover.

AgroIBIS is a state-of-the-art land surface model with options to simulate dynamic vegetation (:ref:`Kucharik et al. 2000 <Kuchariketal2000>`) and interactive crop management (:ref:`Kucharik and Brye 2003 <KucharikBrye2003>`). The interactive crop management parameterizations from AgroIBIS (March 2003 version) were coupled as a proof-of-concept to the Community Land Model version 3 [CLM3.0, :ref:`Oleson et al. (2004) <Olesonetal2004>` ] (not published), then coupled to the CLM3.5 (:ref:`Levis et al. 2009 <Levisetal2009>`) and later released to the community with CLM4CN (:ref:`Levis et al. 2012 <Levisetal2012>`), and CLM4.5BGC. Additional updates after the release of CLM4.5 were available by request (:ref:`Levis et al. 2016 <Levisetal2016>`), and those are now incorporated into CLM5.

With interactive crop management and, therefore, a more accurate representation of agricultural landscapes, we hope to improve the CLM's simulated biogeophysics and biogeochemistry. These advances may improve fully coupled simulations with the Community Earth System Model (CESM), while helping human societies answer questions about changing food, energy, and water resources in response to climate, environmental, land use, and land management change (e.g., :ref:`Kucharik and Brye 2003 <KucharikBrye2003>`; :ref:`Lobell et al. 2006 <Lobelletal2006>`). As implemented here, the crop model uses the same physiology as the natural vegetation but with uses different crop-specific parameter values, phenology, and allocation, as well as fertilizer and irrigation management.

.. _Crop plant functional types:

Crop plant functional types
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To allow crops to coexist with natural vegetation in a grid cell, the vegetated land unit is separated into a naturally vegetated land unit and a managed crop land unit. Unlike the plant functional types (PFTs) in the naturally vegetated land unit, the managed crop PFTs in the managed crop land unit do not share soil columns and thus permit for differences in the land management between crops. Each crop type has a rainfed and an irrigated PFT that are on independent soil columns. Crop grid cell coverage is assigned from satellite data (similar to all natural PFTs), and the managed crop type proportions within the crop area is based on the dataset created by :ref:`Portmann et al. (2010)<Portmannetal2010>` for present day. New in CLM5, crop area is extrapolated through time using the dataset provided by Land Use Model Intercomparison Project (LUMIP), which is part of CMIP6 Land use timeseries (:ref:`Lawrence et al. 2016 <Lawrenceetal2016>`). For more details about how crop distributions are determined, see Chapter :numref:`rst_Transient Landcover Change`.

CLM5 includes ten actively managed crop types (temperate soybean, tropical soybean, temperate corn, tropical corn, spring wheat, cotton, rice, sugarcane, miscanthus, and switchgrass) that are chosen based on the availability of corresponding algorithms in AgroIBIS and as developed by :ref:`Badger and Dirmeyer (2015)<BadgerandDirmeyer2015>` and described by :ref:`Levis et al. (2016)<Levisetal2016>`, or from available observations as described by :ref:`Cheng et al. (2019)<Chengetal2019>`. The representations of sugarcane, rice, cotton, tropical corn, and tropical soy were new in CLM5; miscanthus and switchgrass were added after the CLM5 release. Sugarcane and tropical corn are both C4 plants and are therefore represented using the temperate corn functional form. Tropical soybean uses the temperate soybean functional form, while rice and cotton use the wheat functional form. In tropical regions, parameter values were developed for the Amazon Basin, and planting date window is shifted by six months relative to the Northern Hemisphere. Plantation areas of bioenergy crops are projected to expand throughout the 21st century as a major energy source to replace fossil fuels and mitigate climate change. Miscanthus and switchgrass are perennial bioenergy crops and have quite different physiological traits and land management practices than annual crops, such as longer growing seasons, higher productivity, and lower demands for nutrients and water. About 70% of biofuel aboveground biomass (leaf & livestem) is removed at harvest. Parameter values were developed by using observation data collected at the University of Illinois Energy Farm located in Central Midwestern United States (:ref:`Cheng et al., 2019<Chengetal2019>`).

In addition, CLM's default list of plant functional types (PFTs) includes an irrigated and unirrigated unmanaged C3 crop (:numref:`Table Crop plant functional types`) treated as a second C3 grass. The unmanaged C3 crop is only used when the crop model is not active and has grid cell coverage assigned from satellite data, and the unmanaged C3 irrigated crop type is currently not used since irrigation requires the crop model to be active. The default list of PFTs also includes twenty-one inactive crop PFTs that do not yet have associated parameters required for active management. Each of the inactive crop types is simulated using the parameters of the spatially closest associated crop type that is most similar to the functional type (e.g., C3 or C4), which is required to maintain similar phenological parameters based on temperature thresholds. Information detailing which parameters are used for each crop type is included in :numref:`Table Crop plant functional types`. It should be noted that PFT-level history output merges all crop types into the actively managed crop type, so analysis of crop-specific output will require use of the land surface dataset to remap the yields of each actively and inactively managed crop type. Otherwise, the actively managed crop type will include yields for that crop type and all inactively managed crop types that are using the same parameter set.

.. _Table Crop plant functional types:

.. table:: Crop plant functional types (PFTs) included in CLM5BGCCROP.

 ===  ===========================  ================  ===========================
 IVT  Plant function types (PFTs)  Management Class  Crop Parameters Used
 ===  ===========================  ================  ===========================
  15  c3 unmanaged rainfed crop    none              not applicable
  16  c3 unmanaged irrigated crop  none              not applicable
  17  rainfed temperate corn       active            rainfed temperate corn
  18  irrigated temperate corn     active            irrigated temperate corn
  19  rainfed spring wheat         active            rainfed spring wheat
  20  irrigated spring wheat       active            irrigated spring wheat
  21  rainfed winter wheat         inactive          rainfed spring wheat
  22  irrigated winter wheat       inactive          irrigated spring wheat
  23  rainfed temperate soybean    active            rainfed temperate soybean
  24  irrigated temperate soybean  active            irrigated temperate soybean
  25  rainfed barley               inactive          rainfed spring wheat
  26  irrigated barley             inactive          irrigated spring wheat
  27  rainfed winter barley        inactive          rainfed spring wheat
  28  irrigated winter barley      inactive          irrigated spring wheat
  29  rainfed rye                  inactive          rainfed spring wheat
  30  irrigated rye                inactive          irrigated spring wheat
  31  rainfed winter rye           inactive          rainfed spring wheat
  32  irrigated winter rye         inactive          irrigated spring wheat
  33  rainfed cassava              inactive          rainfed rice
  34  irrigated cassava            inactive          irrigated rice
  35  rainfed citrus               inactive          rainfed spring wheat
  36  irrigated citrus             inactive          irrigated spring wheat
  37  rainfed cocoa                inactive          rainfed rice
  38  irrigated cocoa              inactive          irrigated rice
  39  rainfed coffee               inactive          rainfed rice
  40  irrigated coffee             inactive          irrigated rice
  41  rainfed cotton               active            rainfed cotton
  42  irrigated cotton             active            irrigated cotton
  43  rainfed datepalm             inactive          rainfed cotton
  44  irrigated datepalm           inactive          irrigated cotton
  45  rainfed foddergrass          inactive          rainfed spring wheat
  46  irrigated foddergrass        inactive          irrigated spring wheat
  47  rainfed grapes               inactive          rainfed spring wheat
  48  irrigated grapes             inactive          irrigated spring wheat
  49  rainfed groundnuts           inactive          rainfed rice
  50  irrigated groundnuts         inactive          irrigated rice
  51  rainfed millet               inactive          rainfed tropical corn
  52  irrigated millet             inactive          irrigated tropical corn
  53  rainfed oilpalm              inactive          rainfed rice
  54  irrigated oilpalm            inactive          irrigated rice
  55  rainfed potatoes             inactive          rainfed spring wheat
  56  irrigated potatoes           inactive          irrigated spring wheat
  57  rainfed pulses               inactive          rainfed spring wheat
  58  irrigated pulses             inactive          irrigated spring wheat
  59  rainfed rapeseed             inactive          rainfed spring wheat
  60  irrigated rapeseed           inactive          irrigated spring wheat
  61  rainfed rice                 active            rainfed rice
  62  irrigated rice               active            irrigated rice
  63  rainfed sorghum              inactive          rainfed tropical corn
  64  irrigated sorghum            inactive          irrigated tropical corn
  65  rainfed sugarbeet            inactive          rainfed spring wheat
  66  irrigated sugarbeet          inactive          irrigated spring wheat
  67  rainfed sugarcane            active            rainfed sugarcane
  68  irrigated sugarcane          active            irrigated sugarcane
  69  rainfed sunflower            inactive          rainfed spring wheat
  70  irrigated sunflower          inactive          irrigated spring wheat
  71  rainfed miscanthus           active            rainfed miscanthus
  72  irrigated miscanthus         active            irrigated miscanthus
  73  rainfed switchgrass          active            rainfed switchgrass
  74  irrigated switchgrass        active            irrigated switchgrass
  75  rainfed tropical corn        active            rainfed tropical corn
  76  irrigated tropical corn      active            irrigated tropical corn
  77  rainfed tropical soybean     active            rainfed tropical soybean
  78  irrigated tropical soybean   active            irrigated tropical soybean
 ===  ===========================  ================  ===========================

.. _Phenology:

Phenology
^^^^^^^^^

CLM5-BGC includes evergreen, seasonally deciduous (responding to changes in day length), and stress deciduous (responding to changes in temperature and/or soil moisture) phenology algorithms (Chapter :numref:`rst_Vegetation Phenology and Turnover`). CLM5-BGC-crop uses the AgroIBIS crop phenology algorithm, consisting of three distinct phases.

Phase 1 starts at planting and ends with leaf emergence, phase 2 continues from leaf emergence to the beginning of grain fill, and phase 3 starts from the beginning of grain fill and ends with physiological maturity and harvest.

.. _Planting:

Planting
''''''''

All crops must meet the following requirements between the minimum planting date and the maximum planting date (for the northern hemisphere) in :numref:`Table Crop phenology parameters`:

.. math::
   :label: 25.1

   \begin{array}{c}
   {T_{10d} >T_{p} } \\
   {T_{10d}^{\min } >T_{p}^{\min } }  \\
   {GDD_{8} \ge GDD_{\min } }
   \end{array}

where :math:`{T}_{10d}` is the 10-day running mean of :math:`{T}_{2m}`, (the simulated 2-m air temperature during each model time step) and :math:`T_{10d}^{\min}` is the 10-day running mean of :math:`T_{2m}^{\min }` (the daily minimum of :math:`{T}_{2m}`). :math:`{T}_{p}` and :math:`T_{p}^{\min }` are crop-specific coldest planting temperatures (:numref:`Table Crop phenology parameters`), :math:`{GDD}_{8}` is the 20-year running mean growing degree-days (units are °C day) tracked from April through September (NH) above 8°C with maximum daily increments of 30 degree-days (see equation :eq:`25.3`), and :math:`{GDD}_{min }`\ is the minimum growing degree day requirement (:numref:`Table Crop phenology parameters`). :math:`{GDD}_{8}` does not change as quickly as :math:`{T}_{10d}` and :math:`T_{10d}^{\min }`, so it determines whether it is warm enough for the crop to be planted in a grid cell, while the 2-m air temperature variables determine the day when the crop may be planted if the :math:`{GDD}_{8}` threshold is met. If the requirements in equation :eq:`25.1` are not met by the maximum planting date, crops are still planted on the maximum planting date as long as :math:`{GDD}_{8} > 0`. In the southern hemisphere (SH) the NH requirements apply 6 months later.

At planting, each crop seed pool is assigned 3 gC m\ :sup:`-2` from its grain product pool. The seed carbon is transferred to the leaves upon leaf emergence. An equivalent amount of seed leaf N is assigned given the PFT's C to N ratio for leaves (:math:`{CN}_{leaf}` in :numref:`Table Crop allocation parameters`; this differs from AgroIBIS, which uses a seed leaf area index instead of seed C). The model updates the average growing degree-days necessary for the crop to reach vegetative and physiological maturity, :math:`{GDD}_{mat}`, according to the following AgroIBIS rules:

.. math::
   :label: 25.2

   \begin{array}{lll}
   GDD_{{\rm mat}}^{{\rm corn,sugarcane}} =0.85 GDD_{{\rm 8}} & {\rm \; \; \; and\; \; \; }& 950 <GDD_{{\rm mat}}^{{\rm corn,sugarcane}} <1850{}^\circ {\rm days} \\
   GDD_{{\rm mat}}^{{\rm spring\ wheat,cotton}} =GDD_{{\rm 0}} & {\rm \; \; \; and\; \; \; } & GDD_{{\rm mat}}^{{\rm spring\ wheat,cotton}} <1700{}^\circ {\rm days} \\
   GDD_{{\rm mat}}^{{\rm temp.soy}} =GDD_{{\rm 10}} & {\rm \; \; \; and\; \; \; } & GDD_{{\rm mat}}^{{\rm temp.soy}} <1900{}^\circ {\rm days} \\
   GDD_{{\rm mat}}^{{\rm rice}} =GDD_{{\rm 0}} & {\rm \; \; \; and\; \; \; } & GDD_{{\rm mat}}^{{\rm rice}} <2100{}^\circ {\rm days} \\
   GDD_{{\rm mat}}^{{\rm trop.soy}} =GDD_{{\rm 10}} & {\rm \; \; \; and\; \; \; } & GDD_{{\rm mat}}^{{\rm trop.soy}} <2100{}^\circ {\rm days}
   \end{array}

where :math:`{GDD}_{0}`, :math:`{GDD}_{8}`, and :math:`{GDD}_{10}` are the 20-year running mean growing degree-days tracked from April through September (NH) over 0°C, 8°C, and 10°C, respectively, with maximum daily increments of 26 degree-days (for :math:`{GDD}_{0}`) or 30 degree-days (for :math:`{GDD}_{8}` and :math:`{GDD}_{10}`). Equation :eq:`25.3` shows how we calculate :math:`{GDD}_{0}`, :math:`{GDD}_{8}`, and :math:`{GDD}_{10}` for each model timestep:

.. math::
   :label: 25.3

   \begin{array}{lll}
   GDD_{{\rm 0}} =GDD_{0} +T_{2{\rm m}} -T_{f} & \quad {\rm \; \; \; where\; \; \; } & 0 \le T_{2{\rm m}} -T_{f} \le 26{}^\circ {\rm days} \\
   GDD_{{\rm 8}} =GDD_{8} +T_{2{\rm m}} -T_{f} -8 & \quad {\rm \; \; \; where\; \; \; } & 0 \le T_{2{\rm m}} -T_{f} -8\le 30{}^\circ {\rm days} \\
   GDD_{{\rm 10}} =GDD_{10} +T_{2{\rm m}} -T_{f} -10 & \quad {\rm \; \; \; where\; \; \; } & 0 \le T_{2{\rm m}} -T_{f} -10\le 30{}^\circ {\rm days}
   \end{array}

where, if :math:`{T}_{2m}` - :math:`{T}_{f}` takes on values outside the above ranges within a day, then it equals the minimum or maximum value in the range for that day. :math:`{T}_{f}` is the freezing temperature of water and equals 273.15 K, :math:`{T}_{2m}` is the 2-m air temperature in units of K, and *GDD* is in units of degree-days.

.. _Leaf emergence:

Leaf emergence
''''''''''''''

According to AgroIBIS, leaves may emerge when the growing degree-days of soil temperature to 0.05 m depth (:math:`GDD_{T_{soi} }` ), which is tracked since planting, reaches 1 to 5% of :math:`{GDD}_{mat}` (see Phase 2 % :math:`{GDD}_{mat}` in :numref:`Table Crop phenology parameters`). The base temperature threshold values for :math:`GDD_{T_{soi} }` are listed in :numref:`Table Crop phenology parameters` (the same base temperature threshold values are also used for :math:`GDD_{T_{{\rm 2m}} }` in section :numref:`Grain Fill`), and leaf emergence (crop phenology phase 2) starts when this threshold is met. Leaf onset occurs in the first time step of phase 2, at which moment all seed C is transferred to leaf C. Subsequently, the leaf area index generally increases throughout phase 2 until it reaches a predetermined maximum value. Stem and root C also increase throughout phase 2 based on the carbon allocation algorithm in section :numref:`Leaf emergence to grain fill`.

.. _Grain fill:

Grain fill
''''''''''

The grain fill phase (phase 3) begins in one of two ways. The first potential trigger is based on temperature, similar to phase 2. A variable tracked since planting, similar to :math:`GDD_{T_{soi} }` but for 2-m air temperature, :math:`GDD_{T_{{\rm 2m}} }`, must reach a heat unit threshold, *h*, of of 40 to 65% of :math:`{GDD}_{mat}` (see Phase 3 % :math:`{GDD}_{mat}` in :numref:`Table Crop phenology parameters`). For crops with the C4 photosynthetic pathway (temperate and tropical corn, sugarcane), the :math:`{GDD}_{mat}` is based on an empirical function and ranges between 950 and 1850. The second potential trigger for phase 3 is based on leaf area index. When the maximum value of leaf area index is reached in phase 2 (:numref:`Table Crop allocation parameters`), phase 3 begins. In phase 3, the leaf area index begins to decline in response to a background litterfall rate calculated as the inverse of leaf longevity for the PFT as done in the BGC part of the model.

.. _Harvest:

Harvest
'''''''

Harvest is assumed to occur as soon as the crop reaches maturity. When :math:`GDD_{T_{{\rm 2m}} }` reaches 100% of :math:`{GDD}_{mat}` or the number of days past planting reaches a crop-specific maximum (:numref:`Table Crop phenology parameters`), then the crop is harvested. Harvest occurs in one time step using the BGC leaf offset algorithm.

.. _Table Crop phenology parameters:

.. list-table:: Crop phenology and morphology parameters for the active crop plant functional types (PFTs) in CLM5BGCCROP. Numbers in the first row correspond to the list of PFTs in :numref:`Table Crop plant functional types`.
   :header-rows: 1

   * - \
     - temperate corn
     - spring wheat
     - temperate soybean
     - cotton
     - rice
     - sugarcane
     - tropical corn
     - tropical soybean
     - miscanthus
     - switchgrass
   * - IVT
     - 17, 18
     - 19, 20
     - 23, 24
     - 41, 42
     - 61, 62
     - 67, 68
     - 75, 76
     - 77, 78
     - 71, 72
     - 73, 74
   * - :math:`Date_{planting}^{min}`
     - April 1
     - April 1
     - May 1
     - April 1
     - Janurary 1
     - Janurary 1
     - March 20
     - April 15
     - April 1
     - April 1
   * - :math:`Date_{planting}^{max}`
     - June 15
     - June 15
     - June 15
     - May 31
     - Feburary 28
     - March 31
     - April 15
     - June 31
     - June 15
     - June 15
   * - :math:`T_{p}`\(K)
     - 283.15
     - 280.15
     - 286.15
     - 294.15
     - 294.15
     - 294.15
     - 294.15
     - 294.15
     - 283.15
     - 283.15
   * - :math:`T_{p}^{ min }`\(K)
     - 279.15
     - 272.15
     - 279.15
     - 283.15
     - 283.15
     - 283.15
     - 283.15
     - 283.15
     - 279.15
     - 279.15
   * - :math:`{GDD}_{min}` (degree-days)
     - 50
     - 50
     - 50
     - 50
     - 50
     - 50
     - 50
     - 50
     - 50
     - 50
   * - base temperature for GDD (°C)
     - 8
     - 0
     - 10
     - 10
     - 10
     - 10
     - 10
     - 10
     - 8
     - 8
   * - :math:`{GDD}_{mat}` (degree-days)
     - 950-1850
     - ≤ 1700
     - ≤ 1900
     - ≤ 1700
     - ≤ 2100
     - 950-1850
     - 950-1850
     - ≤ 2100
     - 950-1850
     - 950-1850
   * - Phase 2 % :math:`{GDD}_{mat}`
     - 3%
     - 5%
     - 3%
     - 3%
     - 1%
     - 3%
     - 3%
     - 3%
     - 3%
     - 3%
   * - Phase 3 % :math:`{GDD}_{mat}`
     - 65%
     - 60%
     - 50%
     - 50%
     - 40%
     - 65%
     - 50%
     - 50%
     - 40%
     - 40%
   * - Max. growing season length (:math:`mxmat`)
     - 165
     - 150
     - 150
     - 160
     - 150
     - 300
     - 160
     - 150
     - 210
     - 210
   * - :math:`z_{top}^{\max }` (m)
     - 2.5
     - 1.2
     - 0.75
     - 1.5
     - 1.8
     - 4
     - 2.5
     - 1
     - 2.5
     - 2.5
   * - SLA (m :sup:`2` leaf g :sup:`-1` C)
     - 0.05
     - 0.035
     - 0.035
     - 0.035
     - 0.035
     - 0.05
     - 0.05
     - 0.035
     - 0.057
     - 0.049
   * - :math:`\chi _{L}` index
     - -0.5
     - -0.5
     - -0.5
     - -0.5
     - -0.5
     - -0.5
     - -0.5
     - -0.5
     - -0.5
     - -0.5
   * - grperc
     - 0.11
     - 0.11
     - 0.11
     - 0.11
     - 0.11
     - 0.11
     - 0.11
     - 0.11
     - 0.11
     - 0.11
   * - flnr
     - 0.293
     - 0.41
     - 0.41
     - 0.41
     - 0.41
     - 0.293
     - 0.293
     - 0.41
     - 0.293
     - 0.293
   * - fcur
     - 1
     - 1
     - 1
     - 1
     - 1
     - 1
     - 1
     - 1
     - 1
     - 1

Notes:

- :math:`Date_{planting}^{min}` and :math:`Date_{planting}^{max}` are the minimum and maximum planting dates (defining the "sowing window") in the Northern Hemisphere; the corresponding dates in the Southern Hemisphere are shifted by 6 months. (See Sect. :numref:`Planting`.) These parameters can also be set with more geographic variation via input map stream files ``stream_fldFileName_swindow_start`` and ``stream_fldFileName_swindow_end``.
- :math:`T_{p}` and :math:`T_{p}^{ min }` are crop-specific average and coldest planting temperatures, respectively. (See Sect. :numref:`Planting`.)
- :math:`GDD_{min}` is a threshold describing the coolest historical climate a patch can have had in order for a crop to be sown there; see Sect. :numref:`Planting` for details. 
- :math:`GDD_{mat}` is the heat unit index, in units of accumulated growing degree-days, a crop needs to reach maturity.
- :math:`mxmat` is the maximum growing season length (days past planting), at which harvest occurs even if heat unit index has not reached :math:`GDD_{mat}`.
- :math:`z_{top}^{\max }` is the maximum top-of-canopy height of a crop (see Sect. :numref:`Vegetation Structure`).
- SLA is specific leaf area (see Chapter :numref:`rst_Photosynthetic Capacity`).
- :math:`\chi _{L}` is the leaf orientation index, equals -1 for vertical, 0 for random, and 1 for horizontal leaf orientation. (See Sect. :numref:`Canopy Radiative Transfer`.)
- grperc is the growth respiration factor (see Sect. :numref:`Growth Respiration`). 
- flnr is the fraction of leaf N in the Rubisco enzyme (a.k.a. :math:`N_{cb}` in Sect. :numref:`Plant Nitrogen`).
- fcur is the fraction of allocation that goes to currently displayed growth (i.e., that is not sent to storage). See Sect. :numref:`Carbon Allocation to New Growth`.

.. _Allocation:

Allocation
^^^^^^^^^^

Allocation changes based on the crop phenology phases phenology (section :numref:`Phenology`). Simulated C assimilation begins every year upon leaf emergence in phase 2 and ends with harvest at the end of phase 3; therefore, so does the allocation of such C to the crop's leaf, live stem, fine root, and reproductive pools.

Typically, C:N ratios in plant tissue vary throughout the growing season and tend to be lower during early growth stages and higher in later growth stages. In order to account for this seasonal change, two sets of C:N ratios are established in CLM for the leaf, stem, and fine root of crops: one during the leaf emergence phase (phenology phase 2), and a second during grain fill phase (phenology phase 3). This modified C:N ratio approach accounts for the nitrogen retranslocation that occurs during the grain fill phase (phase 3) of crop growth. Leaf, stem, and root C:N ratios for phase 2 are calculated using the new CLM5 carbon and nitrogen allocation scheme (Chapter :numref:`rst_CN Allocation`), which provides a target C:N value (:numref:`Table Crop allocation parameters`) and allows C:N to vary through time. During grain fill (phase 3) of the crop growth cycle, a portion of the nitrogen in the plant tissues is moved to a storage pool to fulfill nitrogen demands of organ (reproductive pool) development, such that the resulting C:N ratio of the plant tissue is reflective of measurements at harvest. All C:N ratios were determined by calibration process, through comparisons of model output versus observations of plant carbon throughout the growing season.

The BGC part of the model keeps track of a term representing excess maintenance respiration, which supplies the carbon required for maintenance respiration during periods of low photosynthesis (Chapter :numref:`rst_Plant Respiration`). Carbon supply for excess maintenance respiration cannot continue to happen after harvest for annual crops, so at harvest the excess respiration pool is turned into a flux that extracts CO\ :sub:`2` directly from the atmosphere. This way any excess maintenance respiration remaining at harvest is eliminated as if such respiration had not taken place.

.. _Leaf emergence to grain fill:

Leaf emergence
''''''''''''''

During phase 2, the allocation coefficients (fraction of available C) to
each C pool are defined as:

.. math::
   :label: 25.4

   \begin{array}{l} {a_{repr} =0} \\ {a_{froot} =a_{froot}^{i} -(a_{froot}^{i} -a_{froot}^{f} )\frac{GDD_{T_{{\rm 2m}} } }{GDD_{{\rm mat}} } {\rm \; \; \; where\; \; \; }\frac{GDD_{T_{{\rm 2m}} } }{GDD_{{\rm mat}} } \le 1} \\ {a_{leaf} =(1-a_{froot} )\cdot \frac{a_{leaf}^{i} (e^{-b} -e^{-b\frac{GDD_{T_{{\rm 2m}} } }{h} } )}{e^{-b} -1} {\rm \; \; \; where\; \; \; }b=0.1} \\ {a_{livestem} =1-a_{repr} -a_{froot} -a_{leaf} } \end{array}

where :math:`a_{leaf}^{i}`, :math:`a_{froot}^{i}`, and :math:`a_{froot}^{f}` are initial and final values of these coefficients (:numref:`Table Crop allocation parameters`), and *h* is a heat unit threshold defined in section :numref:`Grain fill`. At a crop-specific maximum leaf area index, :math:`{L}_{max}` (:numref:`Table Crop allocation parameters`), carbon allocation is directed exclusively to the fine roots.

.. _Grain fill to harvest:

Grain fill
''''''''''

The calculation of :math:`a_{froot}` remains the same from phase 2 to phase 3. During grain fill (phase 3), other allocation coefficients change to:

.. math::
   :label: 25.5

   \begin{array}{ll}
   a_{leaf} =a_{leaf}^{i,3} & {\rm when} \quad a_{leaf}^{i,3} \le a_{leaf}^{f} \quad {\rm else} \\
   a_{leaf} =a_{leaf} \left(1-\frac{GDD_{T_{{\rm 2m}} } -h}{GDD_{{\rm mat}} d_{L} -h} \right)^{d_{alloc}^{leaf} } \ge a_{leaf}^{f} & {\rm where} \quad \frac{GDD_{T_{{\rm 2m}} } -h}{GDD_{{\rm mat}} d_{L} -h} \le 1 \\
    \\
   a_{livestem} =a_{livestem}^{i,3} & {\rm when} \quad a_{livestem}^{i,3} \le a_{livestem}^{f} \quad {\rm else} \\
   a_{livestem} =a_{livestem} \left(1-\frac{GDD_{T_{{\rm 2m}} } -h}{GDD_{{\rm mat}} d_{L} -h} \right)^{d_{alloc}^{stem} } \ge a_{livestem}^{f} & {\rm where} \quad \frac{GDD_{T_{{\rm 2m}} } -h}{GDD_{{\rm mat}} d_{L} -h} \le 1 \\
    \\
   a_{repr} =1-a_{froot} -a_{livestem} -a_{leaf}
   \end{array}

where :math:`a_{leaf}^{i,3}` and :math:`a_{livestem}^{i,3}` (initial values) equal the last :math:`a_{leaf}` and :math:`a_{livestem}` calculated in phase 2, :math:`d_{L}`, :math:`d_{alloc}^{leaf}` and :math:`d_{alloc}^{stem}` are leaf area index and leaf and stem allocation decline factors, and :math:`a_{leaf}^{f}` and :math:`a_{livestem}^{f}` are final values of these allocation coefficients (:numref:`Table Crop allocation parameters`).

.. _Nitrogen retranslocation for crops:

Nitrogen retranslocation for crops
''''''''''''''''''''''''''''''''''

Nitrogen retranslocation in crops occurs when nitrogen that was used for tissue growth of leaves, stems, and fine roots during the early growth season is remobilized and used for grain development (:ref:`Pollmer et al. 1979 <Pollmeretal1979>`, :ref:`Crawford et al. 1982 <Crawfordetal1982>`, :ref:`Simpson et al. 1983 <Simpsonetal1983>`, :ref:`Ta and Weiland 1992 <TaWeiland1992>`, :ref:`Barbottin et al. 2005 <Barbottinetal2005>`, :ref:`Gallais et al. 2006 <Gallaisetal2006>`, :ref:`Gallais et al. 2007 <Gallaisetal2007>`). Nitrogen allocation for crops follows that of natural vegetation, is supplied in CLM by the soil mineral nitrogen pool, and depends on C:N ratios for leaves, stems, roots, and organs. Nitrogen demand during organ development is fulfilled through retranslocation from leaves, stems, and roots. Nitrogen retranslocation is initiated at the beginning of the grain fill stage for all crops except soybean, for which retranslocation is after LAI decline. Nitrogen stored in the leaf and stem is moved into a storage retranslocation pool for all crops, and for wheat and rice, nitrogen in roots is also released into the retranslocation storage pool. The quantity of nitrogen mobilized depends on the C:N ratio of the plant tissue and is calculated as

.. math::
   :label: 25.6

   leaf\_ to\_ retransn=N_{leaf} -\frac{C_{leaf} }{CN_{leaf}^{f} }

.. math::
   :label: 25.7

   stemn\_ to\_ retransn=N_{stem} -\frac{C_{stem} }{CN_{stem}^{f} }

.. math::
   :label: 25.8

   frootn\_ to\_ retransn=N_{froot} -\frac{C_{froot} }{CN_{froot}^{f} }

where :math:`{C}_{leaf}`, :math:`{C}_{stem}`, and :math:`{C}_{froot}` is the carbon in the plant leaf, stem, and fine root, respectively, :math:`{N}_{leaf}`, :math:`{N}_{stem}`, and :math:`{N}_{froot}` is the nitrogen in the plant leaf, stem, and fine root, respectively, and :math:`CN^f_{leaf}`, :math:`CN^f_{stem}`, and :math:`CN^f_{froot}` is the post-grain fill C:N ratio of the leaf, stem, and fine root respectively (:numref:`Table Crop allocation parameters`). Since C:N measurements are often taken from mature crops, pre-grain development C:N ratios for leaves, stems, and roots in the model are optimized to allow maximum nitrogen accumulation for later use during organ development, and post-grain fill C:N ratios are assigned the same as crop residue. After nitrogen is moved into the retranslocated pool, the nitrogen in this pool is used to meet plant nitrogen demand by assigning the available nitrogen from the retranslocated pool equal to the plant nitrogen demand for each organ (:math:`{CN_{[organ]}^{f} }` in :numref:`Table Crop allocation parameters`). Once the retranslocation pool is depleted, soil mineral nitrogen pool is used to fulfill plant nitrogen demands.

.. _Harvest to food and seed:

Harvest
'''''''

Whereas live crop C and N in grain was formerly transferred to the litter pool upon harvest, CLM5 splits this between "food" and "seed" pools. In the former—more generally a "crop product" pool—C and N decay to the atmosphere over one year, similar to how the wood product pools work. The latter is used in the subsequent year to account for the C and N required for crop seeding.

Live leaf and stem biomass at harvest is transferred to biofuel, removed residue, and/or litter pools.

For the biofuel crops Miscanthus and switchgrass, 70% of live leaf and stem biomass at harvest is transferred to the crop product pool as described for "food" harvest above. This value can be changed for these crops—or set to something other than the default zero for any other crop—with the parameter :math:`biofuel\_harvfrac` (0-1).

50% of any remaining live leaf and stem biomass at harvest (after biofuel removal, if any) is removed to the crop product pool to represent off-field uses such as use for animal feed and bedding. This value can be changed with the parameter :math:`crop\_residue\_removal\_frac` (0–1). The default 50% is derived from :ref:`Smerald et al. 2023 <Smeraldetal2023>`, who found a global average of 50% of residues left on the field. This includes residues burned in the field, meaning that our implementation implictly assumes the CLM crop burning representation will handle those residues appropriately. 

The following equations illustrate how this works. Subscript :math:`p` refers to either the leaf or live stem biomass pool.

.. math::
   :label: 25.9

     CF_{p,biofuel} = \left({CS_{p} \mathord{\left/ {\vphantom {CS_{p}  \Delta t}} \right.} \Delta t}
     \right) * biofuel\_harvfrac

.. math::
   :label: harv_c_to_removed_residue

     CF_{p,removed\_residue} = \left({CS_{p} \mathord{\left/ {\vphantom {CS_{p}  \Delta t}} \right.} \Delta t}
     \right) * (1 - biofuel\_harvfrac) * crop\_residue\_removal\_frac

.. math::
   :label: 25.11

     CF_{p,litter} = \left({CS_{p} \mathord{\left/ {\vphantom {CS_{p}  \Delta t}} \right.} \Delta t}
     \right) * \left( 1-biofuel\_harvfrac  \right) * \left( 1-crop\_residue\_removal\_frac  \right) +CF_{p,alloc}

with corresponding nitrogen fluxes:

.. math::
   :label: 25.12

     NF_{p,biofuel} = \left({NS_{p} \mathord{\left/ {\vphantom {NS_{p}  \Delta t}} \right.} \Delta t}
     \right) * biofuel\_harvfrac

.. math::
   :label: harv_n_to_removed_residue

     NF_{p,removed\_residue} = \left({NS_{p} \mathord{\left/ {\vphantom {NS_{p}  \Delta t}} \right.} \Delta t}
     \right) * \left( 1 - biofuel\_harvfrac \right) * crop\_residue\_removal\_frac

.. math::
   :label: 25.14

     NF_{p,litter} = \left({NS_{p} \mathord{\left/ {\vphantom {NS_{p}  \Delta t}} \right.} \Delta t}
     \right) *  \left( 1-biofuel\_harvfrac  \right) *  \left( 1-crop\_residue\_removal\_frac  \right)

where CF is the carbon flux, CS is stored carbon, NF is the nitrogen flux, NS is stored nitrogen, and :math:`biofuel\_harvfrac` is the harvested fraction of leaf/livestem for biofuel feedstocks.

Annual food crop yields (g dry matter m\ :sup:`-2`) can be calculated by saving the GRAINC_TO_FOOD_ANN variable once per year, then postprocessing with Equation :eq:`25.15`. This calculation assumes that grain C is 45% of the total dry weight. Additionally, harvest is not typically 100% efficient, so analysis needs to assume that harvest efficiency is less---we use 85%.

.. math::
   :label: 25.15

     \text{Grain yield} = \frac{GRAINC\_TO\_FOOD\_ANN)*0.85}{0.45}

.. _Table Crop allocation parameters:

.. table:: Crop allocation parameters for the active crop plant functional types (PFTs) in CLM5BGCCROP. Numbers in the first row correspond to the list of PFTs in :numref:`Table Crop plant functional types`.

 ===========================================  ==============  ============  ==================  ======  ======  =========  =============  ================  ================  ================
 \                                            temperate corn  spring wheat  temperate soybean   cotton  rice    sugarcane  tropical corn  tropical soybean  miscanthus        switchgrass
 ===========================================  ==============  ============  ==================  ======  ======  =========  =============  ================  ================  ================
 IVT                                          17, 18          19, 20        23, 24              41, 42  61, 62  67, 68     75, 76         77, 78            71, 72            73, 74
 :math:`a_{leaf}^{i}`                         0.6             0.9           0.85                0.85    0.75    0.6        0.6            0.85              0.9               0.7
 :math:`{L}_{max}` (m :sup:`2`  m :sup:`-2`)  5               7             6                   6       7       5          5              6                 10                6.5
 :math:`a_{froot}^{i}`                        0.1             0.05          0.2                 0.2     0.1     0.1        0.1            0.2               0.11              0.14
 :math:`a_{froot}^{f}`                        0.05            0             0.2                 0.2     0       0.05       0.05           0.2               0.09              0.09
 :math:`a_{leaf}^{f}`                         0               0             0                   0       0       0          0              0                 0                 0
 :math:`a_{livestem}^{f}`                     0               0.05          0.3                 0.3     0.05    0          0              0.3               0                 0
 :math:`d_{L}`                                1.05            1.05          1.05                1.05    1.05    1.05       1.05           1.05              1.05              1.05
 :math:`d_{alloc}^{stem}`                     2               1             5                   5       1       2          2              5                 2                 2
 :math:`d_{alloc}^{leaf}`                     5               3             2                   2       3       5          5              2                 5                 5
 :math:`{CN}_{leaf}`                          25              20            20                  20      20      25         25             20                25                25
 :math:`{CN}_{stem}`                          50              50            50                  50      50      50         50             50                50                50
 :math:`{CN}_{froot}`                         42              42            42                  42      42      42         42             42                42                42
 :math:`CN^f_{leaf}`                          65              65            65                  65      65      65         65             65                65                65
 :math:`CN^f_{stem}`                          120             100           130                 130     100     120        120            130               120               120
 :math:`CN^f_{froot}`                         0               40            0                   0       40      0          0              0                 0                 0
 :math:`{CN}_{grain}`                         50              50            50                  50      50      50         50             50                50                50
 ===========================================  ==============  ============  ==================  ======  ======  =========  =============  ================  ================  ================

Notes: Crop growth phases and corresponding variables are described throughout the text. :math:`{CN}_{leaf}`, :math:`{CN}_{stem}`, and :math:`{CN}_{froot}` are the target C:N ratios used during the leaf emergence phase (phase 2).

.. _Other Features:

Other Features
^^^^^^^^^^^^^^

.. _Physical Crop Characteristics:

Physical Crop Characteristics
'''''''''''''''''''''''''''''
Leaf area index (*L*) is calculated as a function of specific leaf area (SLA, :numref:`Table Crop phenology parameters`) and leaf C. Stem area index (*S*) is equal to 0.1\ *L* for temperate and tropical corn, sugarcane, switchgrass, and miscanthus and 0.2\ *L* for other crops, as in AgroIBIS. All live C and N pools go to 0 after crop harvest, but the *S* is kept at 0.25 to simulate a post-harvest "stubble" on the ground.

Crop heights at the top and bottom of the canopy, :math:`{z}_{top}` and :math:`{z}_{bot}` (m), come from the AgroIBIS formulation:

.. math::
   :label: 25.16

   \begin{array}{l}
   {z_{top} =z_{top}^{\max } \left(\frac{L}{L_{\max } -1} \right)^{2} \ge 0.05{\rm \; where\; }\frac{L}{L_{\max } -1} \le 1} \\
   {z_{bot} =0.02{\rm m}}
   \end{array}

where :math:`z_{top}^{\max }` is the maximum top-of-canopy height of the crop (:numref:`Table Crop phenology parameters`) and :math:`L_{\max }` is the maximum leaf area index (:numref:`Table Crop allocation parameters`).

.. _Interactive fertilization:

Interactive Fertilization
'''''''''''''''''''''''''
CLM simulates fertilization by adding nitrogen directly to the soil mineral nitrogen pool to meet crop nitrogen demands using both industrial fertilizer and manure application. CLM's separate crop land unit ensures that natural vegetation will not access the fertilizer applied to crops. Fertilizer in CLM5BGCCROP is prescribed by crop functional types and varies spatially for each year based on the LUMIP land use and land cover change time series (LUH2 for historical and SSPs for future) (:ref:`Lawrence et al. 2016 <Lawrenceetal2016>`). One of two fields is used to prescribe industrial fertilizer based on the type of simulation. For non-transient simulations, annual fertilizer application in g N/m\ :sup:`2`/yr is specified on the land surface data set by the field CONST_FERTNITRO_CFT. In transient simulations, annual fertilizer application is specified on the land use time series file by the field FERTNITRO_CFT, which is also in g N/m\ :sup:`2`/yr. The values for both of these fields come from the LUMIP time series for each year. In addition to the industrial fertilizer, background manure fertilizer is specified on the parameter file by the field ``manunitro``. For perennial bioenergy crops, little fertilizer (56kg/ha/yr) is applied to switchgrass and no fertilizer is applied to Miscanthus. Note these rates are only based on local land management practices at the University of Illinois Energy Farm located in Central Midwestern United States :ref:`(Cheng et al., 2019)<Chengetal2019>` rather than the LUMIP timeseries. For the current CLM5BGCCROP, manure N is applied at a rate of 0.002 kg N/m\ :sup:`2`/yr. Because previous versions of CLM (e.g., CLM4) had rapid denitrification rates, fertilizer is applied slowly to minimize N loss (primarily through denitrification) and maximize plant uptake. The current implementation of CLM5 inherits this legacy, although denitrification rates are slower in the current version of the model (:ref:`Koven et al. 2013 <Kovenetal2013>`). As such, fertilizer application begins during the leaf emergence phase of crop development (phase 2) and continues for 20 days, which helps reduce large losses of nitrogen from leaching and denitrification during the early stage of crop development. The 20-day period is chosen as an optimization to limit fertilizer application to the emergence stage. A fertilizer counter in seconds, *f*, is set as soon as the leaf emergence phase for crops initiates:

.. math::
   :label: 25.17

    f = n \times 86400

where *n* is set to 20 fertilizer application days and 86400 is the number of seconds per day. When the crop enters phase 2 (leaf emergence) of its growth cycle, fertilizer application begins by initializing fertilizer amount to the total fertilizer at each column within the grid cell divided by the initialized *f*. Fertilizer is applied and *f* is decremented each time step until a zero balance on the counter is reached.

.. _Biological nitrogen fixation for soybeans:

Biological nitrogen fixation for soybeans
'''''''''''''''''''''''''''''''''''''''''
Biological N fixation for soybeans is calculated by the fixation and uptake of nitrogen module (Chapter :numref:`rst_FUN`) and is the same as N fixation in natural vegetation. Unlike natural vegetation, where a fraction of each PFT are N fixers, all soybeans are treated as N fixers.

.. _Latitude vary base tempereature for growing degree days:

Latitudinal variation in base growth tempereature
'''''''''''''''''''''''''''''''''''''''''''''''''
For most crops, :math:`GDD_{T_{{\rm 2m}} }` (growing degree days since planting) is the same in all locations. However, for both rainfed and irrigated spring wheat and sugarcane, the calculation of :math:`GDD_{T_{{\rm 2m}} }` allows for latitudinal variation:

.. math::
   :label: 25.18

   latitudinal\ variation\ in\ base\ T = \left\{
   \begin{array}{lr}
   baset +12 - 0.4 \times latitude &\qquad 0 \le latitude \le 30 \\
   baset +12 + 0.4 \times latitude &\qquad -30 \le latitude \le 0
   \end{array} \right\}

where :math:`baset` is the *base temperature for GDD* (7\ :sup:`th` row) in :numref:`Table Crop phenology parameters`. Such latitudinal variation in base temperature could slow :math:`GDD_{T_{{\rm 2m}} }` accumulation extend the growing season for regions within 30°S to 30°N for spring wheat and sugarcane.

.. _Separate reproductive pool:

Separate reproductive pool
''''''''''''''''''''''''''
One notable difference between natural vegetation and crops is the presence of reproductive carbon and nitrogen pools. Accounting for the reproductive pools helps determine whether crops are performing reasonably through yield calculations. The reproductive pool is maintained similarly to the leaf, stem, and fine root pools, but allocation of carbon and nitrogen does not begin until the grain fill stage of crop development. Equation :eq:`25.5` describes the carbon and nitrogen allocation coefficients to the reproductive pool. In CLM5BGCCROP, as allocation declines in stem, leaf, and root pools (see section :numref:`Grain fill to harvest`) during the grain fill stage of growth, increasing amounts of carbon and nitrogen are available for grain development.

.. _Tillage:

Tillage
'''''''
Tillage is represented as an enhancement of the decomposition rate coefficient; see section :numref:`decomp_mgmt_modifiers`.

.. _The irrigation model:

The irrigation model
--------------------

The CLM includes the option to irrigate cropland areas that are equipped for irrigation. The application of irrigation responds dynamically to the soil moisture conditions simulated by the CLM. This irrigation algorithm is based loosely on the implementation of :ref:`Ozdogan et al. (2010) <Ozdoganetal2010>`.

When irrigation is enabled, the crop areas of each grid cell are divided into irrigated and rainfed fractions according to a dataset of areas equipped for irrigation (:ref:`Portmann et al. 2010 <Portmannetal2010>`). Irrigated and rainfed crops are placed on separate soil columns, so that irrigation is only applied to the soil beneath irrigated crops.

In irrigated croplands, a check is made once per day to determine whether irrigation is required on that day. This check is made in the first time step after 6 AM local time. Irrigation is required if crop leaf area :math:`>` 0, and the available soil water is below a specified threshold.

The soil moisture deficit :math:`D_{irrig}` is

.. math::
   :label: 25.61

   D_{irrig} = \left\{
   \begin{array}{lr}
   w_{target} - w_{avail} &\qquad w_{thresh} > w_{avail} \\
   0 &\qquad w_{thresh} \le w_{avail}
   \end{array} \right\}

where :math:`w_{target}` is the irrigation target soil moisture (mm)

.. math::
   :label: 25.62

   w_{target} = \sum_{j=1}^{N_{irr}} \theta_{target} \Delta z_{j} \ .

The irrigation moisture threshold (mm) is

.. math::
   :label: 25.63

   w_{thresh} = f_{thresh} \left(w_{target} - w_{wilt}\right) + w_{wilt}

where :math:`w_{wilt}` is the wilting point soil moisture (mm)

.. math::
   :label: 25.64

   w_{wilt} = \sum_{j=1}^{N_{irr}} \theta_{wilt} \Delta z_{j} \ ,

and :math:`f_{thresh}` is a tuning parameter.  The available moisture in the soil (mm) is

.. math::
   :label: 25.65

   w_{avail} = \sum_{j=1}^{N_{irr}} \theta_{j} \Delta z_{j} \ ,

Note that :math:`w_{target}` is truly supposed to give the target soil moisture value that we're shooting for whenever irrigation happens; then the soil moisture deficit :math:`D_{irrig}` gives the difference between this target value and the current soil moisture. The irrigation moisture threshold :math:`w_{thresh}`, on the other hand, gives a threshold at which we decide to do any irrigation at all. The way this is written allows for the possibility that one may not want to irrigate every time there becomes even a tiny soil moisture deficit. Instead, one may want to wait until the deficit is larger before initiating irrigation; at that point, one doesn't want to just irrigate up to the "threshold" but instead up to the higher "target". The target should always be greater than or equal to the threshold.

:math:`N_{irr}` is the index of the soil layer corresponding to a specified depth :math:`z_{irrig}` (:numref:`Table Irrigation parameters`) and :math:`\Delta z_{j}` is the thickness of the soil layer in layer :math:`j` (section :numref:`Vertical Discretization`). :math:`\theta_{j}` is the volumetric soil moisture in layer :math:`j` (section :numref:`Soil Water`). :math:`\theta_{target}` and :math:`\theta_{wilt}` are the target and wilting point volumetric soil moisture values, respectively, and are determined by inverting :eq:`7.94` using soil matric potential parameters :math:`\Psi_{target}` and :math:`\Psi_{wilt}` (:numref:`Table Irrigation parameters`). After the soil moisture deficit :math:`D_{irrig}` is calculated, irrigation in an amount equal to :math:`\frac{D_{irrig}}{T_{irrig}}` (mm/s) is applied uniformly over the irrigation period :math:`T_{irrig}` (s). Irrigation water is applied directly to the ground surface, bypassing canopy interception (i.e., added to :math:`{q}_{grnd,liq}`: section :numref:`Canopy Water`).

To conserve mass, irrigation is removed from river water storage (Chapter :numref:`rst_River Transport Model (RTM)`). When river water storage is inadequate to meet irrigation demand, there are two options: 1) the additional water can be removed from the ocean model, or 2) the irrigation demand can be reduced such that river water storage is maintained above a specified threshold.

.. _Table Irrigation parameters:

.. table:: Irrigation parameters

 +--------------------------------------+-------------+
 | Parameter                            |             |
 +======================================+=============+
 | :math:`f_{thresh}`                   |  1.0        |
 +--------------------------------------+-------------+
 | :math:`z_{irrig}`       (m)          |  0.6        |
 +--------------------------------------+-------------+
 | :math:`\Psi_{target}`   (mm)         | -3400       |
 +--------------------------------------+-------------+
 | :math:`\Psi_{wilt}`     (mm)         | -150000     |
 +--------------------------------------+-------------+

.. add a reference to surface data in chapter2
 To accomplish this we downloaded data of percent irrigated and percent rainfed corn, soybean, and temperate cereals (wheat, barley, and rye) (:ref:`Portmann et al. 2010 <Portmannetal2010>`), available online from *ftp://ftp.rz.uni-frankfurt.de/pub/uni-frankfurt/physische\_geographie/hydrologie/public/data/MIRCA2000/harvested\_area\_grids.*
