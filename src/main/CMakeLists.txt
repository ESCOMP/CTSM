# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

list(APPEND clm_sources
  ctsm_ColumnType.F90
  ctsm_PftCon.F90
  ctsm_FuncPedotransfer.F90
  ctsm_GridcellType.F90
  ctsm_LandunitType.F90
  ctsm_PatchType.F90
  ctsm_Topo.F90
  ctsm_AbortUtils.F90
  ctsm_Atm2Lnd.F90
  ctsm_Atm2LndType.F90
  ctsm_Accumulators.F90
  ctsm_VarCon.F90
  ctsm_VarCtl.F90
  ctsm_VarPar.F90
  clm_varsur.F90 
  ctsm_ColumnVarCon.F90
  ctsm_Decomp.F90
  ctsm_FilterCol.F90
  ctsm_Glc2Lnd.F90
  ctsm_GlacierBehavior.F90
  ctsm_InitSubgrid.F90
  ctsm_InitVertical.F90
  ctsm_LandunitVarCon.F90
  ctsm_Lnd2Glc.F90
  ctsm_NcdIoUtils.F90
  ctsm_OrganicFile.F90
  ctsm_ParamUtil.F90
  ctsm_SubgridAve.F90
  ctsm_SubgridWeights.F90
  ctsm_SurfReadUtils.F90
  )

sourcelist_to_parent(clm_sources)
