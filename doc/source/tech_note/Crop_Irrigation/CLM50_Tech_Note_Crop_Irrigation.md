<!-- Define math macros here -->
$\newcommand{\gddthreshmat}{GDD_\textrm{*mat}}$
$\newcommand{\gddthreshmatbl}{\gddthreshmat^\textrm{bl}}$
$\newcommand{\gddaccsoil}{GDD_{T_\textrm{soi}}}$
$\newcommand{\ttwom}{T_\textrm{2m}}$
$\newcommand{\gddacctwom}{GDD_{\ttwom}}$
$\newcommand{\gddzero}{GDD_0}$
$\newcommand{\gddeight}{GDD_8}$
$\newcommand{\gddten}{GDD_{10}}$
$\newcommand{\gddx}{GDD_x}$
$\newcommand{\gddzerorun}{\overline{\gddzero}^\textrm{20yr}}$
$\newcommand{\gddeightrun}{\overline{\gddeight}^\textrm{20yr}}$
$\newcommand{\gddtenrun}{\overline{\gddten}^\textrm{20yr}}$
$\newcommand{\gddxrun}{\overline{\gddx}^\textrm{20yr}}$
$\newcommand{\gddxrunbl}{\overline{\gddx}^\textrm{20-yr,bl}}$
$\newcommand{\gddxdaymax}{\gddx^\textrm{daymax}}$
$\newcommand{\huithreshlfemerg}{h_\textrm{*lfemerg}}$
$\newcommand{\huithreshgrain}{h_\textrm{*grain}}$
$\newcommand{\parambaset}{T_\textrm{base}}$
$\newcommand{\paramztopmx}{z_\textrm{top}^\textrm{max}}$

(rst_crops and irrigation)=

# Crops and Irrigation

(the-crop-model)=

## The crop model: cash and bioenergy crops

### Introduction

Groups developing Earth System Models generally account for the human footprint on the landscape in simulations of historical and future climates. Traditionally we have represented this footprint with natural vegetation types and particularly grasses because they resemble many common crops. Most modeling efforts have not incorporated more explicit representations of land management such as crop type, planting, harvesting, tillage, fertilization, and irrigation, because global scale datasets of these factors have lagged behind vegetation mapping. As this begins to change, we increasingly find models that will simulate the biogeophysical and biogeochemical effects not only of natural but also human-managed land cover.

AgroIBIS is a state-of-the-art land surface model with options to simulate dynamic vegetation ({ref}`Kucharik et al. 2000 <Kuchariketal2000>`) and interactive crop management ({ref}`Kucharik and Brye 2003 <KucharikBrye2003>`). The interactive crop management parameterizations from AgroIBIS (March 2003 version) were coupled as a proof-of-concept to the Community Land Model version 3 \[CLM3.0, {ref}`Oleson et al. (2004) <Olesonetal2004>` \] (not published), then coupled to the CLM3.5 ({ref}`Levis et al. 2009 <Levisetal2009>`) and later released to the community with CLM4CN ({ref}`Levis et al. 2012 <Levisetal2012>`), and CLM4.5BGC. Additional updates after the release of CLM4.5 were available by request ({ref}`Levis et al. 2016 <Levisetal2016>`), and those are now incorporated into CLM5 and later.

With interactive crop management and, therefore, a more accurate representation of agricultural landscapes, we hope to improve CLM's simulated biogeophysics and biogeochemistry. These advances may improve fully coupled simulations with the Community Earth System Model (CESM), while helping human societies answer questions about changing food, energy, and water resources in response to climate, environmental, land use, and land management change (e.g., {ref}`Kucharik and Brye 2003 <KucharikBrye2003>`; {ref}`Lobell et al. 2006 <Lobelletal2006>`). As implemented here, the crop model uses the same physiology as the natural vegetation but with uses different crop-specific parameter values, phenology, and allocation, as well as fertilizer and irrigation management.

(crop-plant-functional-types)=

### Crop plant functional types

To allow crops to coexist with natural vegetation in a grid cell, the vegetated land unit is separated into a naturally vegetated land unit and a managed crop land unit. Unlike the plant functional types (PFTs) in the naturally vegetated land unit, the managed crop PFTs in the managed crop land unit do not share soil columns and thus permit for differences in the land management between crops. Each crop type has a rainfed and an irrigated PFT that are on independent soil columns. Crop area distributions are defined as explained in Sects. {numref}`Surface Data` and {numref}`rst_Transient Landcover Change`; see Sect. {numref}`Surface Heterogeneity and Data Structure` for more information on land units and soil columns.

CLM includes ten actively managed crop types (temperate soybean, tropical soybean, temperate corn, tropical corn, spring wheat, cotton, rice, sugarcane, miscanthus, and switchgrass) that are chosen based on the availability of corresponding algorithms in AgroIBIS and as developed by {ref}`Badger and Dirmeyer (2015)<BadgerandDirmeyer2015>` and described by {ref}`Levis et al. (2016)<Levisetal2016>`, or from available observations as described by {ref}`Cheng et al. (2019)<Chengetal2019>`. Sugarcane and tropical corn are both C4 plants and are therefore represented using the temperate corn functional form. Tropical soybean uses the temperate soybean functional form, while rice and cotton use the wheat functional form. In tropical regions, parameter values were developed for the Amazon Basin, and planting date window is shifted by six months relative to the Northern Hemisphere. Plantation areas of bioenergy crops are projected to expand throughout the 21st century as a major energy source to replace fossil fuels and mitigate climate change. Miscanthus and switchgrass are perennial bioenergy crops and have quite different physiological traits and land management practices than annual crops, such as longer growing seasons, higher productivity, and lower demands for nutrients and water. About 70% of biofuel aboveground biomass (leaf & livestem) is removed at harvest. Parameter values were developed by using observation data collected at the University of Illinois Energy Farm located in Central Midwestern United States ({ref}`Cheng et al., 2019<Chengetal2019>`).

In addition, CLM's default list of plant functional types (PFTs) includes an irrigated and unirrigated unmanaged C3 crop ({numref}`Table Crop plant functional types`) treated as a second C3 grass. The unmanaged C3 crop is only used when the crop model is not active and has grid cell coverage assigned from satellite data, and the unmanaged C3 irrigated crop type is currently not used since irrigation requires the crop model to be active. The default list of PFTs also includes twenty-one inactive crop PFTs that do not yet have associated parameters required for active management. Each of the inactive crop types is simulated using the parameters of the spatially closest associated crop type that is most similar to the functional type (e.g., C3 or C4), which is required to maintain similar phenological parameters based on temperature thresholds. Information detailing which parameters are used for each crop type is included in {numref}`Table Crop plant functional types`. It should be noted that PFT-level history output merges all crop types into the actively managed crop type, so analysis of crop-specific output will require use of the land surface dataset to remap the yields of each actively and inactively managed crop type. Otherwise, the actively managed crop type will include yields for that crop type and all inactively managed crop types that are using the same parameter set.

(table crop plant functional types)=

```{eval-rst}
.. table:: Crop plant functional types (PFTs) included in CLM with managed crops on (`BgcCrop` component sets).

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
```

(phenology)=

### Phenology

CLM uses the AgroIBIS crop phenology algorithm, consisting of three distinct phases.

Phase 1 starts at planting and ends with leaf emergence, phase 2 continues from leaf emergence to the beginning of grain fill, and phase 3 starts from the beginning of grain fill and ends with physiological maturity and harvest.

(planting)=

#### Planting

Each crop can be planted in each gridcell once per year, in the "sowing window." This is, by default, defined for each crop in each gridcell, with the start and end of the window drawn from the input files `stream_fldFileName_swindow_start` and `stream_fldFileName_swindow_end`, respectively. These sowing windows are centered on the sowing dates developed for phase 3 of the Inter-Sectoral Impacts Model Intercomparison Project's Agriculture sector runs ({ref}`Jägermeyr et al., 2021 <jägermeyretal2021>`, {ref}`Rabin et al., 2023 <rabinetal2023>`). The widths of the sowing windows are based on the original sowing windows from the `min_NH_planting_date`, `max_NH_planting_date`, `min_SH_planting_date`, and `max_SH_planting_date` parameters, which remain on the parameter file but are ignored by default (i.e., unless you set `cropcals_rx = .false.` in `user_nl_clm`). See {ref}`running-with-custom-crop-calendars`.

To be planted, a crop patch must meet the following requirements sometime within its sowing window:

$$
\begin{array}{c}
{T_{10d} >T_{p} } \\
{T_{10d}^{\min } >T_{p}^{\min } }  \\
{\gddeightrun \ge GDD_{\min } }
\end{array}
$$ (25.1)

where ${T}_{10d}$ is the 10-day running mean of $\ttwom$, (the simulated 2-m air temperature during each model time step) and $T_{10d}^{\min}$ is the 10-day running mean of $\ttwom^{\min }$ (the daily minimum of $\ttwom$). ${T}_{p}$ and $T_{p}^{\min }$ are crop-specific coldest planting temperatures ({numref}`Table Crop phenology parameters`), $\gddeightrun$ is the 20-year running mean growing degree-days (units are °C day) tracked from April through September (NH) above 8°C with maximum daily increments of 30 degree-days (see equation {eq}`25.3`), and ${GDD}_{min }$is the minimum growing degree day requirement ({numref}`Table Crop phenology parameters`). $\gddeightrun$ does not change as quickly as ${T}_{10d}$ and $T_{10d}^{\min }$, so it determines whether it is warm enough for the crop to be planted in a grid cell, while the 2-m air temperature variables determine the day when the crop may be planted if the $\gddeightrun$ threshold is met. If the requirements in equation {eq}`25.1` are not met by the maximum planting date, crops are still planted on the maximum planting date as long as $\gddeightrun > 0$.

At planting, each crop seed pool is assigned 3 gC m{sup}`-2` from its grain product pool. The seed carbon is transferred to the leaves upon leaf emergence. An equivalent amount of seed leaf N is assigned given the PFT's C to N ratio for leaves (${CN}_{leaf}$ in {numref}`Table Crop allocation parameters`; this differs from AgroIBIS, which uses a seed leaf area index instead of seed C).

#### Maturity requirement
At planting, CLM determines how many growing degree-days will be needed for the crop to reach maturity and thus be harvested. By default (i.e., `cropcals_rx_adapt = .true.`), this is set according to two input files with PFT-specific maps:
- `stream_fldfilename_cultivar_gdds` ($\gddthreshmatbl$): the average growing-degree days to reach maturity in the "baseline" period.
- `stream_fldFileName_gdd20_baseline` ($\gddxrunbl$): the means over the baseline period of $\gddzerorun$, $\gddeightrun$, and $\gddtenrun$.

Maturity requirement, $\gddthreshmat$, is then calculated as:

$$
\gddthreshmat = \max \left( 1,\ \gddthreshmatbl \times \frac{\gddxrun}{\gddxrunbl} \right),
$$ (gddmat-rx-adapt)

where $x$ is 0 (wheat, cotton, and rice), 8 (corn, sugarcane, _Miscanthus_, and switchgrass), or 10 (soybean). This allows the maturity requirement to "adapt" over time, being lower in cool periods and higher in warm periods. The baseline period is the 1980-2009 growing seasons (i.e., seasons where planting occurred in those calendar years, inclusive); baseline values were calculated based on a half-degree, land-only run with CRU-JRA climate forcings (Rabin et al., in prep.). The minimum value of 1 avoids numeric issues when $\gddthreshmat$ is in the denominator of a calculation.

- **Check baseline period**

If `cropcals_rx_adapt` is false but `cropcals_rx` is true, the calculation is just $\gddthreshmat = \gddthreshmatbl$.

If both `cropcals_rx_adapt` and `cropcals_rx` are false, or if $\gddthreshmatbl$ is negative, then CLM sets $\gddthreshmat$ according to various crop-specific rules based on $\gddx$, the PFT-specific parameter `hybgdd`, and hard-coded minimum and maximum values.

Equation {eq}`25.3` shows how we calculate $\gddzero$, $\gddeight$, and $\gddten$ for each model timestep:

$$
\gddx = \gddx + \frac{\max \left( \gddxdaymax,\ \min \left[ 0,\ \ttwom - 273.15 - x \right] \right)}{48} 
$$ (25.3)

where $\ttwom$ is the 2-m air temperature (K), 273.15 K is the freezing temperature of water, and $GDD$ is in units of °C-days. $\gddxdaymax$, the maximum daily growing degree-day accumulation, is 26°C for $x=0$ and 30°C for $x=8$ and $x=10$.

- **Is there a pre-existing symbol for number of timesteps in a day that we could use instead of 48?**

By default, the $\gddx$ values are set to zero at the beginning of the "$\gddx$ season" and then accumulated through its end: from April 1 through September 30 in the Northern Hemisphere and from October 1 through March 31 in the Southern Hemisphere. (Setting `stream_gdd20_seasons = .true.` would instead take those start and end dates from PFT-specific maps in the input files `stream_fldFileName_gdd20_season_start` and `stream_fldFileName_gdd20_season_end`, respectively; however, this is not scientifically supported.) At the end of each $\gddx$ season, the final value of $\gddx$ is incorporated into $\gddxrun$ like so:

$$
\gddxrun = \frac{\gddxrun \times \min(n-1, 19) + \gddx}{20},
$$ (update-gddxrun)

where $n$ is the number of years that $\gddxrun$ gas been calculated for. Note that this is not a true rolling 20-year mean, which would come with a memory cost in the simulation as a 20-member array would need to be saved for each PFT.

(leaf-emergence)=

#### Leaf emergence

The "leaf emergence" phase is the period of vegetative growth between when the leaves first emerge from the soil to when filling of the reproductive organ begins.

According to AgroIBIS, leaves may emerge when the growing degree-days of soil temperature to 0.05 m depth ($\gddaccsoil$ ), which is tracked since planting, reaches 1 to 5% of $\gddthreshmat$ (see $h_{lfemerg}$ in {numref}`Table Crop phenology parameters`). The base temperature threshold values for $\gddaccsoil$ are listed in {numref}`Table Crop phenology parameters` (the same base temperature threshold values are also used for $\gddacctwom$ in section {numref}`Grain Fill`), and leaf emergence (crop phenology phase 2) starts when this threshold is met. Leaf onset occurs in the first time step of phase 2, at which moment all seed C is transferred to leaf C. Subsequently, the leaf area index generally increases throughout phase 2 until it reaches a predetermined maximum value. Stem and root C also increase throughout phase 2 based on the carbon allocation algorithm in section {numref}`Leaf emergence to grain fill`.

(grain fill)=

#### Grain fill

The grain fill phase (phase 3) begins in one of two ways. The first potential trigger is based on temperature, similar to phase 2. A variable tracked since planting, similar to $\gddaccsoil$ but for 2-m air temperature, $\gddacctwom$, must reach a heat unit threshold, $\huithreshgrain$, of 40 to 65% of $\gddthreshmat$ (see {numref}`Table Crop phenology parameters`). The second potential trigger for phase 3 is based on leaf area index. When the maximum value of leaf area index is reached in phase 2 ({numref}`Table Crop allocation parameters`), phase 3 begins. In phase 3, the leaf area index begins to decline in response to a background litterfall rate calculated as the inverse of leaf longevity for the PFT as done in the BGC part of the model.

(harvest)=

#### Harvest

Harvest is assumed to occur as soon as the crop reaches maturity. When $\gddacctwom$ reaches 100% of $\gddthreshmat$ or the number of days past planting reaches a crop-specific maximum ({numref}`Table Crop phenology parameters`), then the crop is harvested. Harvest occurs in one time step using the BGC leaf offset algorithm.

(table crop phenology parameters)=

```{eval-rst}
.. list-table:: Crop phenology and morphology parameters for the active crop plant functional types (PFTs) in CLM with managed crops on (``BgcCrop`` component sets). Numbers in the first row correspond to the list of PFTs in :numref:`Table Crop plant functional types`. Where there are two values in a cell, they refer to the rainfed and irrigated functional types, respectively.
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
   * - :math:`\parambaset` (°C)
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
   * - :math:`\huithreshlfemerg` (% :math:`\gddthreshmat`)
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
   * - :math:`\huithreshgrain` (% :math:`\gddthreshmat`)
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
```

Notes:

- $T_{p}$ and $T_{p}^{ min }$ are crop-specific average and coldest planting temperatures, respectively. (See Sect. {numref}`Planting`.)
- $GDD_{min}$ is a threshold describing the coolest historical climate a patch can have had in order for a crop to be sown there; see Sect. {numref}`Planting` for details.
- $\parambaset$ is the minimum temperature for accumulating growing degree-days.
- $h_{lfemerg}$ and $h_{grainfill}$ are, respectively, the threshold fractions of $\gddthreshmat$ a crop must reach to enter the leaf-emergence phase (phase 2) and grain-filling phase (phase 3).
- $mxmat$ is the maximum growing season length (days past planting), at which harvest occurs even if heat unit index has not reached $\gddthreshmat$.
- $\paramztopmx$ is the maximum top-of-canopy height of a crop (see Sect. {numref}`Vegetation Structure`).
- SLA is specific leaf area (see Chapter {numref}`rst_Photosynthetic Capacity`).
- $\chi _{L}$ is the leaf orientation index, equals -1 for vertical, 0 for random, and 1 for horizontal leaf orientation. (See Sect. {numref}`Canopy Radiative Transfer`.)
- grperc is the growth respiration factor (see Sect. {numref}`Growth Respiration`).
- flnr is the fraction of leaf N in the Rubisco enzyme (a.k.a. $N_{cb}$ in Sect. {numref}`Plant Nitrogen`).
- fcur is the fraction of allocation that goes to currently displayed growth (i.e., that is not sent to storage). See Sect. {numref}`Carbon Allocation to New Growth`.

(allocation)=

### Allocation

Allocation changes based on the crop phenology phase (section {numref}`Phenology`). Simulated C assimilation begins every year upon leaf emergence in phase 2 and ends with harvest at the end of phase 3; therefore, so does the allocation of such C to the crop's leaf, live stem, fine root, and reproductive pools.

Typically, C:N ratios in plant tissue vary throughout the growing season and tend to be lower during early growth stages and higher in later growth stages. In order to account for this seasonal change, two sets of C:N ratios are established in CLM for the leaf, stem, and fine root of crops: one during the leaf emergence phase (phenology phase 2), and a second during grain fill phase (phenology phase 3). This modified C:N ratio approach accounts for the nitrogen retranslocation that occurs during the grain fill phase (phase 3) of crop growth. Leaf, stem, and root C:N ratios for phase 2 are calculated using the standard CLM carbon and nitrogen allocation scheme (Chapter {numref}`rst_CN Allocation`), which provides a target C:N value ({numref}`Table Crop allocation parameters`) and allows C:N to vary through time. During grain fill (phase 3) of the crop growth cycle, a portion of the nitrogen in the plant tissues is moved to a storage pool to fulfill nitrogen demands of organ (reproductive pool) development, such that the resulting C:N ratio of the plant tissue is reflective of measurements at harvest. All C:N ratios were determined by calibration process, through comparisons of model output versus observations of plant carbon throughout the growing season.

Carbon supply for excess maintenance respiration (Chapter {numref}`rst_Plant Respiration`) cannot continue to happen after harvest for annual crops, so at harvest the excess respiration pool is turned into a flux that extracts CO{sub}`2` directly from the atmosphere. This way any excess maintenance respiration remaining at harvest is eliminated as if such respiration had not taken place.

(leaf emergence to grain fill)=

#### Leaf emergence

During phase 2, the allocation coefficients (fraction of available C) to
each C pool are defined as:

$$
\begin{array}{l} {a_{repr} =0} \\
{a_{froot} =a_{froot}^{i} -(a_{froot}^{i} -a_{froot}^{f} ) \times {\rm min}\left(\frac{\gddacctwom }{\gddthreshmat }, 1\right)} \\
{a_{leaf} =(1-a_{froot} ) \times \frac{a_{leaf}^{i} (e^{-b} -e^{-b\frac{\gddacctwom }{\huithreshgrain} } )}{e^{-b} -1} {\rm \; \; \; where\; \; \; }b=0.1} \\
{a_{livestem} =1-a_{repr} -a_{froot} -a_{leaf} } \end{array}
$$ (eq-lfemerg-allocations)

where $a_{leaf}^{i}$, $a_{froot}^{i}$, and $a_{froot}^{f}$ are initial and final values of these coefficients, and $\huithreshgrain$ is the heat unit threshold to enter the grain-filling phase. At a crop-specific maximum leaf area index, ${L}_{max}$, carbon allocation is directed almost exclusively to the fine roots, with only 0.001% of carbon going to leaves. See {numref}`Table Crop allocation parameters` for parameter values.

(grain fill to harvest)=

#### Grain fill

The calculation of $a_{froot}$ remains the same from phase 2 (Eq. [](#eq-lfemerg-allocations)) to phase 3. During grain fill (phase 3), other allocation coefficients change to:

$$
\begin{array}{ll}
a_{leaf} =a_{leaf}^{i,3} & {\rm when} \quad a_{leaf}^{i,3} \le a_{leaf}^{f} \quad {\rm else} \\
a_{leaf} =a_{leaf} \left(1-\frac{\gddacctwom - \huithreshgrain}{\gddthreshmat d_{L} - \huithreshgrain} \right)^{d_{alloc}^{leaf} } \ge a_{leaf}^{f} & {\rm where} \quad \frac{\gddacctwom - \huithreshgrain}{\gddthreshmat d_{L} - \huithreshgrain} \le 1 \\
 \\
a_{livestem} =a_{livestem}^{i,3} & {\rm when} \quad a_{livestem}^{i,3} \le a_{livestem}^{f} \quad {\rm else} \\
a_{livestem} =a_{livestem} \left(1-\frac{\gddacctwom - \huithreshgrain}{\gddthreshmat d_{L} - \huithreshgrain} \right)^{d_{alloc}^{stem} } \ge a_{livestem}^{f} & {\rm where} \quad \frac{\gddacctwom - \huithreshgrain}{\gddthreshmat d_{L} - \huithreshgrain} \le 1 \\
 \\
a_{repr} =1-a_{froot} -a_{livestem} -a_{leaf}
\end{array}
$$ (25.5)

where $a_{leaf}^{i,3}$ and $a_{livestem}^{i,3}$ (initial values) equal the last $a_{leaf}$ and $a_{livestem}$ calculated in phase 2, $d_{L}$, $d_{alloc}^{leaf}$ and $d_{alloc}^{stem}$ are leaf area index and leaf and stem allocation decline factors, $a_{leaf}^{f}$ and $a_{livestem}^{f}$ are final values of these allocation coefficients, and $\huithreshgrain$ is the heat unit threshold to enter the grain-filling phase. See {numref}`Table Crop allocation parameters` for parameter values.

As in the leaf-emergence phase (Sect {numref}`leaf emergence to grain fill`), at a crop-specific maximum leaf area index, ${L}_{max}$, leaf allocation is reduced to 0.001%. The rest of the carbon that would have gone to leaves instead goes to the reproductive pool.

(nitrogen-retranslocation-for-crops)=

#### Nitrogen retranslocation for crops

Nitrogen retranslocation in crops occurs when nitrogen that was used for tissue growth of leaves, stems, and fine roots during the early growth season is remobilized and used for grain development ({ref}`Pollmer et al. 1979 <Pollmeretal1979>`, {ref}`Crawford et al. 1982 <Crawfordetal1982>`, {ref}`Simpson et al. 1983 <Simpsonetal1983>`, {ref}`Ta and Weiland 1992 <TaWeiland1992>`, {ref}`Barbottin et al. 2005 <Barbottinetal2005>`, {ref}`Gallais et al. 2006 <Gallaisetal2006>`, {ref}`Gallais et al. 2007 <Gallaisetal2007>`). Nitrogen allocation for crops follows that of natural vegetation, is supplied in CLM by the soil mineral nitrogen pool, and depends on C:N ratios for leaves, stems, roots, and organs. Nitrogen demand during organ development is fulfilled through retranslocation from leaves, stems, and roots. Nitrogen retranslocation is initiated at the beginning of the grain fill stage for all crops except soybean, for which retranslocation is after LAI decline. Nitrogen stored in the leaf and stem is moved into a storage retranslocation pool for all crops, and for wheat and rice, nitrogen in roots is also released into the retranslocation storage pool. The quantity of nitrogen mobilized depends on the C:N ratio of the plant tissue and is calculated as

$$
leaf\_ to\_ retransn=N_{leaf} -\frac{C_{leaf} }{CN_{leaf}^{f} }
$$ (25.6)

$$
stemn\_ to\_ retransn=N_{stem} -\frac{C_{stem} }{CN_{stem}^{f} }
$$ (25.7)

$$
frootn\_ to\_ retransn=N_{froot} -\frac{C_{froot} }{CN_{froot}^{f} }
$$ (25.8)

where ${C}_{leaf}$, ${C}_{stem}$, and ${C}_{froot}$ is the carbon in the plant leaf, stem, and fine root, respectively, ${N}_{leaf}$, ${N}_{stem}$, and ${N}_{froot}$ is the nitrogen in the plant leaf, stem, and fine root, respectively, and $CN^f_{leaf}$, $CN^f_{stem}$, and $CN^f_{froot}$ is the post-grain fill C:N ratio of the leaf, stem, and fine root respectively ({numref}`Table Crop allocation parameters`). Since C:N measurements are often taken from mature crops, pre-grain development C:N ratios for leaves, stems, and roots in the model are optimized to allow maximum nitrogen accumulation for later use during organ development, and post-grain fill C:N ratios are assigned the same as crop residue. After nitrogen is moved into the retranslocated pool, the nitrogen in this pool is used to meet plant nitrogen demand by assigning the available nitrogen from the retranslocated pool equal to the plant nitrogen demand for each organ (${CN_{[organ]}^{f} }$ in {numref}`Table Crop allocation parameters`). Once the retranslocation pool is depleted, soil mineral nitrogen pool is used to fulfill plant nitrogen demands.

(harvest to food and seed)=

#### Harvest

CLM splits live crop grain C and N between "food" and "seed" pools. In the former—more generally a "crop product" pool—C and N decay to the atmosphere over one year, similar to how the wood product pools work. The latter is used in the subsequent year to account for the C and N required for crop seeding.

Live leaf and stem biomass at harvest is transferred to biofuel, removed residue, and/or litter pools.

For the biofuel crops Miscanthus and switchgrass, 70% of live leaf and stem biomass at harvest is transferred to the crop product pool as described for "food" harvest above. This value can be changed for these crops—or set to something other than the default zero for any other crop—with the parameter $biofuel\_harvfrac$ (0-1).

50% of any remaining live leaf and stem biomass at harvest (after biofuel removal, if any) is removed to the crop product pool to represent off-field uses such as use for animal feed and bedding. This value can be changed with the parameter $crop\_residue\_removal\_frac$ (0–1). The default 50% is derived from {ref}`Smerald et al. 2023 <Smeraldetal2023>`, who found a global average of 50% of residues left on the field. This includes residues burned in the field, meaning that our implementation implictly assumes the CLM crop burning representation will handle those residues appropriately.

The following equations illustrate how this works. Subscript $p$ refers to either the leaf or live stem biomass pool.

$$
CF_{p,biofuel} = \left({CS_{p} \mathord{\left/ {\vphantom {CS_{p}  \Delta t}} \right.} \Delta t}
  \right) \times biofuel\_harvfrac
$$ (25.9)

$$
CF_{p,removed\_residue} = \left({CS_{p} \mathord{\left/ {\vphantom {CS_{p}  \Delta t}} \right.} \Delta t}
  \right) \times (1 - biofuel\_harvfrac) \times crop\_residue\_removal\_frac
$$ (harv_c_to_removed_residue)

$$
CF_{p,litter} = \left({CS_{p} \mathord{\left/ {\vphantom {CS_{p}  \Delta t}} \right.} \Delta t}
  \right) \times \left( 1-biofuel\_harvfrac  \right) \times \left( 1-crop\_residue\_removal\_frac  \right) +CF_{p,alloc}
$$ (25.11)

with corresponding nitrogen fluxes:

$$
NF_{p,biofuel} = \left({NS_{p} \mathord{\left/ {\vphantom {NS_{p}  \Delta t}} \right.} \Delta t}
  \right) \times biofuel\_harvfrac
$$ (25.12)

$$
NF_{p,removed\_residue} = \left({NS_{p} \mathord{\left/ {\vphantom {NS_{p}  \Delta t}} \right.} \Delta t}
  \right) \times \left( 1 - biofuel\_harvfrac \right) \times crop\_residue\_removal\_frac
$$ (harv_n_to_removed_residue)

$$
NF_{p,litter} = \left({NS_{p} \mathord{\left/ {\vphantom {NS_{p}  \Delta t}} \right.} \Delta t}
  \right) \times  \left( 1-biofuel\_harvfrac  \right) \times  \left( 1-crop\_residue\_removal\_frac  \right)
$$ (25.14)

where CF is the carbon flux, CS is stored carbon, NF is the nitrogen flux, NS is stored nitrogen, and $biofuel\_harvfrac$ is the harvested fraction of leaf/livestem for biofuel feedstocks.

Annual food crop yields (g dry matter m{sup}`-2`) can be calculated by saving the GRAINC_TO_FOOD_ANN variable once per year, then postprocessing with Equation {eq}`25.15`. This calculation assumes that grain C is 45% of the total dry weight. Additionally, harvest is not typically 100% efficient, so analysis needs to assume that harvest efficiency is less---we use 85%.

$$
\text{Grain yield} = \frac{GRAINC\_TO\_FOOD\_ANN) \times 0.85}{0.45}
$$ (25.15)

(table crop allocation parameters)=

```{eval-rst}
.. table:: Crop allocation parameters for the active crop plant functional types (PFTs) in CLM with managed crops on (`BgcCrop` component sets). Numbers in the first row correspond to the list of PFTs in :numref:`Table Crop plant functional types`.

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
```

Notes: Crop growth phases and corresponding variables are described throughout the text. ${CN}_{leaf}$, ${CN}_{stem}$, and ${CN}_{froot}$ are the target C:N ratios used during the leaf emergence phase (phase 2).

(other-features)=

### Other Features

(physical-crop-characteristics)=

#### Physical Crop Characteristics

Leaf area index ($L$) is calculated as a function of specific leaf area (SLA, {numref}`Table Crop phenology parameters`) and leaf C. Stem area index ($S$) is equal to 0.1$L$ for temperate and tropical corn, sugarcane, switchgrass, and miscanthus and 0.2$L$ for other crops, as in AgroIBIS. All live C and N pools go to 0 after crop harvest, but the $S$ is kept at 0.25 to simulate a post-harvest "stubble" on the ground.

Crop heights at the top and bottom of the canopy, ${z}_{top}$ and ${z}_{bot}$ (m), come from the AgroIBIS formulation:

$$
\begin{array}{l}
{z_{top} = \paramztopmx \left(\frac{L}{L_{\max } -1} \right)^{2} \ge 0.05{\rm \; where\; }\frac{L}{L_{\max } -1} \le 1} \\
{z_{bot} =0.02{\rm m}}
\end{array}
$$ (25.16)

where $\paramztopmx$ is the maximum top-of-canopy height of the crop ({numref}`Table Crop phenology parameters`) and $L_{\max }$ is the maximum leaf area index ({numref}`Table Crop allocation parameters`).

(interactive-fertilization)=

#### Interactive Fertilization

CLM simulates fertilization by adding nitrogen directly to the soil mineral nitrogen pool to meet crop nitrogen demands using both industrial fertilizer and manure application. CLM's separate crop land unit ensures that natural vegetation will not access the fertilizer applied to crops. Fertilizer in CLM is prescribed by crop functional types and varies spatially for each year based on the LUMIP land use and land cover change time series (LUH2 for historical and SSPs for future) ({ref}`Lawrence et al. 2016 <Lawrenceetal2016>`). One of two fields is used to prescribe industrial fertilizer based on the type of simulation. For non-transient simulations, annual fertilizer application in g N/m{sup}`2`/yr is specified on the land surface data set by the field CONST_FERTNITRO_CFT. In transient simulations, annual fertilizer application is specified on the land use time series file by the field FERTNITRO_CFT, which is also in g N/m{sup}`2`/yr. The values for both of these fields come from the LUMIP time series for each year. In addition to the industrial fertilizer, background manure fertilizer is specified on the parameter file by the field `manunitro`. For perennial bioenergy crops, little fertilizer (56kg/ha/yr) is applied to switchgrass and no fertilizer is applied to Miscanthus. Note these rates are only based on local land management practices at the University of Illinois Energy Farm located in Central Midwestern United States {ref}`(Cheng et al., 2019)<Chengetal2019>` rather than the LUMIP timeseries. Manure N is applied at a rate of 0.002 kg N/m{sup}`2`/yr. Because previous versions of CLM (e.g., CLM4) had rapid denitrification rates, fertilizer is applied slowly to minimize N loss (primarily through denitrification) and maximize plant uptake. The current implementation of CLM inherits this legacy, although denitrification rates are slower in the current version of the model ({ref}`Koven et al. 2013 <Kovenetal2013>`). As such, fertilizer application begins during the leaf emergence phase of crop development (phase 2) and continues for 20 days, which helps reduce large losses of nitrogen from leaching and denitrification during the early stage of crop development. The 20-day period is chosen as an optimization to limit fertilizer application to the emergence stage. A fertilizer counter in seconds, $f$, is set as soon as the leaf emergence phase for crops initiates:

$$
f = n \times 86400
$$ (25.17)

where $n$ is set to 20 fertilizer application days and 86400 is the number of seconds per day. When the crop enters phase 2 (leaf emergence) of its growth cycle, fertilizer application begins by initializing fertilizer amount to the total fertilizer at each column within the grid cell divided by the initialized $f$. Fertilizer is applied and $f$ is decremented each time step until a zero balance on the counter is reached.

(biological-nitrogen-fixation-for-soybeans)=

#### Biological nitrogen fixation for soybeans

Biological N fixation for soybeans is calculated by the fixation and uptake of nitrogen module (Chapter {numref}`rst_FUN`) and is the same as N fixation in natural vegetation. Unlike natural vegetation, where a fraction of each PFT are N fixers, all soybeans are treated as N fixers.

(latitude-vary-base-temperature-for-growing-degree-days)=

#### Latitudinal variation in base growth temperature

For most crops, $\gddacctwom$ (growing degree days since planting) is the same in all locations. However, for both rainfed and irrigated spring wheat and sugarcane, the calculation of $\gddacctwom$ allows for latitudinal variation:

$$
latitudinal\ variation\ in\ base\ T = \left\{
\begin{array}{lr}
\parambaset +12 - 0.4 \times latitude &\qquad 0 \le latitude \le 30 \\
\parambaset +12 + 0.4 \times latitude &\qquad -30 \le latitude \le 0
\end{array} \right\}
$$ (25.18)

where $\parambaset$ is the *base temperature for GDD* (7{sup}`th` row) in {numref}`Table Crop phenology parameters`. Such latitudinal variation in base temperature could slow $\gddacctwom$ accumulation extend the growing season for regions within 30°S to 30°N for spring wheat and sugarcane.

(separate-reproductive-pool)=

#### Separate reproductive pool

One notable difference between natural vegetation and crops is the presence of reproductive carbon and nitrogen pools. Accounting for the reproductive pools helps determine whether crops are performing reasonably through yield calculations. The reproductive pool is maintained similarly to the leaf, stem, and fine root pools, but allocation of carbon and nitrogen does not begin until the grain fill stage of crop development. Equation {eq}`25.5` describes the carbon and nitrogen allocation coefficients to the reproductive pool. In CLM, as allocation declines in stem, leaf, and root pools (see section {numref}`Grain fill to harvest`) during the grain fill stage of growth, increasing amounts of carbon and nitrogen are available for grain development.

(tillage)=

#### Tillage

Tillage is represented as an enhancement of the decomposition rate coefficient; see section {numref}`decomp_mgmt_modifiers`.

(the-irrigation-model)=

## The irrigation model

CLM includes the option to irrigate cropland areas that are equipped for irrigation. The application of irrigation responds dynamically to the simulated soil moisture conditions. This irrigation algorithm is based loosely on the implementation of {ref}`Ozdogan et al. (2010) <Ozdoganetal2010>`.

When irrigation is enabled, the crop areas of each grid cell are divided into irrigated and rainfed fractions according to a dataset of areas equipped for irrigation ({ref}`Portmann et al. 2010 <Portmannetal2010>`). Irrigated and rainfed crops are placed on separate soil columns, so that irrigation is only applied to the soil beneath irrigated crops.

In irrigated croplands, a check is made once per day to determine whether irrigation is required on that day. This check is made in the first time step after 6 AM local time. Irrigation is required if crop leaf area $>$ 0, and the available soil water is below a specified threshold.

The soil moisture deficit $D_{irrig}$ is

$$
D_{irrig} = \left\{
\begin{array}{lr}
w_{target} - w_{avail} &\qquad w_{thresh} > w_{avail} \\
0 &\qquad w_{thresh} \le w_{avail}
\end{array} \right\}
$$ (25.61)

where $w_{target}$ is the irrigation target soil moisture (mm)

$$
w_{target} = \sum_{j=1}^{N_{irr}} \theta_{target} \Delta z_{j} \ .
$$ (25.62)

The irrigation moisture threshold (mm) is

$$
w_{thresh} = f_{thresh} \left(w_{target} - w_{wilt}\right) + w_{wilt}
$$ (25.63)

where $w_{wilt}$ is the wilting point soil moisture (mm)

$$
w_{wilt} = \sum_{j=1}^{N_{irr}} \theta_{wilt} \Delta z_{j} \ ,
$$ (25.64)

and $f_{thresh}$ is a tuning parameter. The available moisture in the soil (mm) is

$$
w_{avail} = \sum_{j=1}^{N_{irr}} \theta_{j} \Delta z_{j} \ ,
$$ (25.65)

Note that $w_{target}$ is truly supposed to give the target soil moisture value that we're shooting for whenever irrigation happens; then the soil moisture deficit $D_{irrig}$ gives the difference between this target value and the current soil moisture. The irrigation moisture threshold $w_{thresh}$, on the other hand, gives a threshold at which we decide to do any irrigation at all. The way this is written allows for the possibility that one may not want to irrigate every time there becomes even a tiny soil moisture deficit. Instead, one may want to wait until the deficit is larger before initiating irrigation; at that point, one doesn't want to just irrigate up to the "threshold" but instead up to the higher "target". The target should always be greater than or equal to the threshold.

$N_{irr}$ is the index of the soil layer corresponding to a specified depth $z_{irrig}$ ({numref}`Table Irrigation parameters`) and $\Delta z_{j}$ is the thickness of the soil layer in layer $j$ (section {numref}`Vertical Discretization`). $\theta_{j}$ is the volumetric soil moisture in layer $j$ (section {numref}`Soil Water`). $\theta_{target}$ and $\theta_{wilt}$ are the target and wilting point volumetric soil moisture values, respectively, and are determined by inverting {eq}`7.94` using soil matric potential parameters $\Psi_{target}$ and $\Psi_{wilt}$ ({numref}`Table Irrigation parameters`). After the soil moisture deficit $D_{irrig}$ is calculated, irrigation in an amount equal to $\frac{D_{irrig}}{T_{irrig}}$ (mm/s) is applied uniformly over the irrigation period $T_{irrig}$ (s). Irrigation water is applied directly to the ground surface, bypassing canopy interception (i.e., added to ${q}_{grnd,liq}$: section {numref}`Canopy Water`).

To conserve mass, irrigation is removed from river water storage (Chapter {numref}`rst_MOSART`). When river water storage is inadequate to meet irrigation demand, there are two options: 1) the additional water can be removed from the ocean model, or 2) the irrigation demand can be reduced such that river water storage is maintained above a specified threshold.

(table irrigation parameters)=

```{eval-rst}
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
```

% add a reference to surface data in chapter2
% To accomplish this we downloaded data of percent irrigated and percent rainfed corn, soybean, and temperate cereals (wheat, barley, and rye) (:ref:`Portmann et al. 2010 <Portmannetal2010>`), available online from *ftp://ftp.rz.uni-frankfurt.de/pub/uni-frankfurt/physische\_geographie/hydrologie/public/data/MIRCA2000/harvested\_area\_grids.*
