# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

list(APPEND clm_sources
  ColumnType.F90
  pftconMod.F90
  FuncPedotransferMod.F90
  GridcellType.F90
  LandunitType.F90
  PatchType.F90
  TopoMod.F90
  abortutils.F90
  atm2lndMod.F90
  atm2lndType.F90
  accumulMod.F90
  clm_varcon.F90
  clm_varctl.F90
  clm_varpar.F90
  clm_varsur.F90 
  column_varcon.F90
  decompMod.F90
  filterColMod.F90
  glc2lndMod.F90
  glcBehaviorMod.F90
  initSubgridMod.F90
  landunit_varcon.F90
  lnd2glcMod.F90
  ncdio_utils.F90
  organicFileMod.F90
  paramUtilMod.F90
  subgridAveMod.F90
  subgridWeightsMod.F90
  surfrdUtilsMod.F90
  )

sourcelist_to_parent(clm_sources)
