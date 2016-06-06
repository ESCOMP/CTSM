# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

list(APPEND clm_sources
  CNSharedParamsMod.F90
  CNSpeciesMod.F90
  CNVegComputeSeedMod.F90
  SpeciesBaseType.F90
  SpeciesIsotopeType.F90
  SpeciesNonIsotopeType.F90
  )

sourcelist_to_parent(clm_sources)
