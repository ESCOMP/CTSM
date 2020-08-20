# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

list(APPEND clm_sources
  ctsm_CNSharedParamsMod.F90
  ctsm_CNSpeciesMod.F90
  ctsm_CNVegComputeSeedMod.F90
  ctsm_SpeciesBaseType.F90
  ctsm_SpeciesIsotopeType.F90
  ctsm_SpeciesNonIsotopeType.F90
  )

sourcelist_to_parent(clm_sources)
