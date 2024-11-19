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
