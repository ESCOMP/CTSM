# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

list(APPEND clm_sources 
  ctsm_Aerosols.F90
  ctsm_DayLength.F90
  ctsm_BalanceCheck.F90
  ctsm_CanopyStateType.F90
  ctsm_EnergyFluxType.F90
  ctsm_GlacierSurfaceMassBalance.F90
  ctsm_HumanIndices.F90
  ctsm_InfiltrationExcessRunoff.F90
  ctsm_Irrigation.F90
  ctsm_LakeConstants.F90
  ctsm_QSat.F90
  ctsm_RootBiophys.F90
  ctsm_SaturatedExcessRunoff.F90
  ctsm_SnowCoverFractionBase.F90
  ctsm_SnowHydrology.F90
  ctsm_SnowSnicar.F90
  ctsm_SoilHydrology.F90
  ctsm_SoilHydrologyType.F90
  ctsm_SoilStateType.F90
  ctsm_SoilWaterRetentionCurve.F90
  ctsm_SolarAbsorbedType.F90
  ctsm_SurfaceAlbedoType.F90
  ctsm_TemperatureType.F90
  ctsm_TotalWaterAndHeat.F90
  ctsm_UrbanParamsType.F90
  ctsm_WaterAtm2LndBulkType.F90
  ctsm_WaterAtm2LndType.F90
  ctsm_WaterBalanceType.F90
  ctsm_WaterDiagnosticBulkType.F90
  ctsm_WaterDiagnosticType.F90
  ctsm_WaterFluxType.F90
  ctsm_WaterFluxBulkType.F90
  ctsm_WaterLnd2AtmBulkType.F90
  ctsm_WaterLnd2AtmType.F90
  ctsm_WaterInfoBaseType.F90
  ctsm_WaterInfoBulkType.F90
  ctsm_WaterInfoIsotopeType.F90
  ctsm_WaterInfoTracerType.F90
  ctsm_WaterStateBulkType.F90
  ctsm_WaterStateType.F90
  ctsm_WaterTracerContainerType.F90
  ctsm_WaterTracerUtils.F90
  ctsm_WaterType.F90
  )

sourcelist_to_parent(clm_sources)
