# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

list(APPEND clm_sources
  ctsm_InitInterpMindist.F90
  ctsm_InitInterpMultilevelBase.F90
  ctsm_InitInterpMultilevelCopy.F90
  ctsm_InitInterpMultilevelInterp.F90
  ctsm_InitInterpMultilevelSnow.F90
  ctsm_InitInterpMultilevelSplit.F90
  ctsm_InitInterpUtils.F90
  )

sourcelist_to_parent(clm_sources)
