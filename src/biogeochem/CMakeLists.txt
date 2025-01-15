# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

list(APPEND clm_sources
  CNSharedParamsMod.F90
  CNPhenologyMod.F90
  CNSpeciesMod.F90
  CNDVType.F90
  DustEmisBase.F90
  DustEmisLeung2023.F90
  DustEmisZender2003.F90
  DustEmisFactory.F90
  CropReprPoolsMod.F90
  CropType.F90
  CNVegStateType.F90
  CNVegCarbonStateType.F90
  CNVegCarbonFluxType.F90
  CNVegMatrixMod.F90
  CNVegNitrogenStateType.F90
  CNVegNitrogenFluxType.F90
  CNCIsoAtmTimeSeriesReadMod.F90
  CNVegComputeSeedMod.F90
  SpeciesBaseType.F90
  SpeciesIsotopeType.F90
  SpeciesNonIsotopeType.F90
  )

sourcelist_to_parent(clm_sources)
