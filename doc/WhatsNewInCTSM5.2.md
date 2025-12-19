# What's new in CTSM 5.2 (tag `ctsm5.2.005`)

## CTSM 5.2: New surface datasets and mksurfdata_esmf tool to create them

### Summary

- All new surface datasets, with updated input datasets.
- New mksurfdata_esmf tool to make global surface datasets.
- New tools to create inputs for regional simulations.
- ne0np4 grid: New 1979 surface dataset and 1979-2026 land use.
- Transient urban and lake by default turned on for transient cases.
- Ocean is run as baresoil rather than wetland (for clm6_0).
- The urban streams file was also updated.
- Update the README files.
- New FATES parameter file: Tree PFT allometry, allometric mode options, leaf maintenance scaling coefficients.
- New surface datasets
- The new surface datasets are incompatible with previous versions (for example the ctsm5.1 series)—ctsm5.2.0 and following versions can NOT use the previous ctsm5.1 datasets, and vice versa.

[!NOTE] See the section below about the new datasets used in their creation. Improvements in how landunits on coastal areas were also made.

## Fields added to the surface datasets in ctsm5.2:

- ORGC, BULK, CFRAG, PHAQ (soil data) (currently NOT used by CTSM)
- mapunits (map units from the soil dataset)
- LANDFRAC_MKSURFDATA (for reference NOT used by CTSM)
- PCT_OCEAN (previously PCT_WETLAND was used)
- Fields removed from the surface datasets in ctsm5.2:

   - AREA
   - PFTDATA_MASK

## New input data used for making surface datasets
- New soil dataset: ISRIC/WISE dataset (Batjes, 2016; https://doi.org/10.1016/j.geoderma.2016.01.034)
- New PFT, soil-color, LAI datasets: Created by Peter J. Lawrence (2022)
- New Glacier datasets: Glacier outlines from RGI version 6 (Arendt et al., 2017).
- vector data for GrIS and AIS retrieved from BedMachine version 4 and version 2 (Morlighem et al., 2017, 2020), respectively.
- 30-arcsec topography/land mask retrieved GMTED2010 (Danielson and Gesch, 2011).
- New urban datasets: Gao and O'Neill (2021) and Gao and Pesaresi (2022), Oleson and Feddema (2020)
- New lake datasets: HydroLake: Messager et. al. (2016)
- New mksurfdata_esmf Tool
- mksurfdata_esmf is a parallel version of mksurfdata that uses ESMF regridding directly so that offline mapping files don't have to be created as a separate step. This allows surface datasets to be created at much higher resolutions.

The build for the tool is based on the CESM/CIME build system and uses cmake. This allows the build to be kept up with changes in CESM. Currently it's only setup and working on Derecho, but this design will enable it to be built and run on any CESM-supported machine (or a machine that a user ports to).

Any input grid from ccs_config can be used, or the user can supply their own mesh file to define the output grid. The user no longer has to add to the list of valid resolutions (as in the now-deprecated mksurfdata_map).

Creation of supported single point datasets. These datasets are created through the use of subset_data.

Test datasets for dynUrban, dynLake, and dynPFT is done with a simple NCO script.

All datasets can be easily made by running make all in the tools/mksurfdata_esmf/ directory.

For detailed instructions, see tools/mksurfdata_esmf/README.md.