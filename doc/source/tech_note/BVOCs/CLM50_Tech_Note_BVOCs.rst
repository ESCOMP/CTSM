.. _rst_Biogenic Volatile Organic Compounds (BVOCs):

Biogenic Volatile Organic Compounds (BVOCs)
===============================================

This chapter briefly describes the biogenic volatile organic compound (BVOC) emissions model implemented in CLM. The CLM3 version (Levis et al. 2003; Oleson et al. 2004) was based on Guenther et al. (1995). Heald et al. (2008) updated this scheme in CLM4 based on Guenther et al (2006). The current version was implemented in CLM4.5 and is based on MEGAN2.1 discussed in detail in Guenther et al. (2012). This update of MEGAN incorporates four main features: 1) expansion to 147 chemical compounds, 2) the treatment of the light-dependent fraction (LDF) for each compound, 3) inclusion of the inhibition of isoprene emission by atmospheric CO\ :sub:`2` and 4) emission factors mapped to the specific PFTs of the CLM.

MEGAN2.1 now describes the emissions of speciated monoterpenes, sesquiterpenes, oxygenated VOCs as well as isoprene. A flexible scheme has been implemented in the CLM to specify a subset of emissions. This allows for additional flexibility in grouping chemical compounds to form the lumped species frequently used in atmospheric chemistry. The mapping or grouping is therefore defined through a namelist parameter in drv\_flds\_in, e.g. megan\_specifier = 'ISOP = isoprene', 'BIGALK pentane + hexane + heptane + tricyclene'.

Terrestrial BVOC emissions from plants to the atmosphere are expressed as a flux, :math:`F_{i}` (:math:`\mu` \ g C m\ :sup:`-2` ground area h\ :sup:`-1`), for emission of chemical compound :math:`i`

.. math::
   :label: ZEqnNum964222

   F_{i} =\gamma _{i} \rho \sum _{j}\varepsilon _{i,j}  \left(wt\right)_{j}

where :math:`\gamma _{i}` is the emission activity factor accounting for responses to meteorological and phenological conditions, :math:`\rho` is the canopy loss and production factor also known as escape efficiency (set to 1), and :math:`\varepsilon _{i,\, j}` (:math:`\mu` \ g C m\ :sup:`-2` ground area h\ :sup:`-1`) is the emission factor at standard conditions of light, temperature, and leaf area for plant functional type *j* with fractional coverage :math:`\left(wt\right)_{j}` (Guenther et al. 2012). The emission activity factor :math:`\gamma _{i}` depends on plant functional type, temperature, LAI, leaf age, and soil moisture (Guenther et al. 2012) For isoprene only, the effect of CO\ :sub:`2` inhibition is now included as described by Heald et al. (2009). Previously, only isoprene was treated as a light-dependent emission. In MEGAN2.1, each chemical compound is assigned a LDF (ranging from 1.0 for isoprene to 0.2 for some monoterpenes, VOCs and acetone). The activity factor for the light response of emissions is therefore estimated as:

.. math::
   :label: 28.2)

   \gamma _{P,\, i} =\left(1-LDF_{i} \right)+\gamma _{P\_ LDF} LDF_{i}

where the LDF activity factor (:math:`\gamma _{P\_ LDF}` ) is specified as a function of PAR as in previous versions of MEGAN.

The values for each emission factor :math:`\epsilon _{i,\, j}` are now available for each of the plant functional types in the CLM and each chemical compound. This information is distributed through an external file, allowing for more frequent and easier updates.

The impact of drought on isoprene emissions is based on the theory proposed by Potosnak et al. (2014). Specifically, isoprene emissions are expected to increase under mild to moderate drought because drought raises leaf temperature, which stimulates isoprene emissions. Under severe drought, however, isoprene emissions are inhibited because substrate supply becomes constrained. Because the effect of leaf temperature is already represented by the leaf temperature module and its influence on isoprene emissions, only the inhibitory effect of severe drought is parameterized as (Wang et al., 2022):

.. math::
   :label: 28.2b)

   \gamma _{sub} = 1 / \left(1 + b_{1} exp(a1 (beta-0.2)) \right)

where a1 = -7.45 and b1=3.26 are empircal parameters.

Compared with Guenther et al. (2012), updates have been made to represent isoprene emissions from high-latitude plants, specifically boreal broadleaf deciduous shrubs (BBDS) and C3 Arctic grass (C3AG), in order to account for acclimation processes. These updates are based on leaf-enclosure and in situ measurements conducted at Toolik Field Station in Alaska, USA (Wang et al., 2024a, 2024b).
For BBDS, the isoprene emission factor is adjusted according to the mean temperature of the previous day as:

.. math::
   :label: 28.2c)

Eopt_bbds = 7.9 exp(0.22 (T_{24}-297.15) )

where T_{24} denotes the mean air temperature of the preceding day (Wang et al., 2024a).
For C3AG, the isoprene emission factor responds over a longer timescale of 10 days and is parameterized as a function of the mean air temperature over the preceding 10 days (T240) as (Wang et al., 2024b):

Eopt_C3AG = exp(0.22 (T_{240}-288.15) )

In addition, a dynamic temperature response curve for C3AG depends on recent temperature history as:

gamma_{T\_C3AG} = Eopt_C3AG exp((C_{C3AG}/R ( 1/303.15 - 1/Tleaf)))

where Tleaf denotes the leaf temperature and C_{C3AG} is the parameter controlling the isoprene temperature response of C3AG and changes varies with T_{240} as:

C_{C3AG} = 95 + 9.5 exp(0.53 (288.15-T_{240}))


References

Potosnak, M. J., LeStourgeon, L., Pallardy, S. G., Hosman, K. P., Gu, L., Karl, T., et al. (2014). Observed and modeled ecosystem isoprene fluxes from an oak-dominated temperate forest and the influence of drought stress. Atmospheric Environment, 84, 314–322. http://www.sciencedirect.com/science/article/pii/S1352231013009059
Wang, H., Lu, X., Seco, R., Stavrakou, T., Karl, T., Jiang, X., et al. (2022). Modeling Isoprene Emission Response to Drought and Heatwaves Within MEGAN Using Evapotranspiration Data and by Coupling With the Community Land Model. Journal of Advances in Modeling Earth Systems, 14(12), e2022MS003174. https://doi.org/10.1029/2022MS003174. https://doi.org/10.1029/2022MS003174
Wang, H., Welch, A., Nagalingam, S., Leong, C., Kittitananuvong, P., Barsanti, K. C., et al. (2024a). Arctic Heatwaves Could Significantly Influence the Isoprene Emissions From Shrubs. Geophysical Research Letters, 51(2), e2023GL107599. https://doi.org/10.1029/2023GL107599
Wang, H., Welch, A. M., Nagalingam, S., Leong, C., Czimczik, C. I., Tang, J., et al. (2024b). High temperature sensitivity of Arctic isoprene emissions explained by sedges. Nature Communications, 15(1), 6144. https://doi.org/10.1038/s41467-024-49960-0
