# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

list(APPEND clm_sources 
  AerosolMod.F90
  DaylengthMod.F90
  GlacierSurfaceMassBalanceMod.F90
  HumanIndexMod.F90
  IrrigationMod.F90
  LakeCon.F90
  QSatMod.F90
  RootBiophysMod.F90
  SnowHydrologyMod.F90
  SnowSnicarMod.F90
  SoilHydrologyType.F90
  SoilStateType.F90
  SoilWaterRetentionCurveMod.F90
  TemperatureType.F90
  TotalWaterAndHeatMod.F90
  UrbanParamsType.F90
  WaterDiagnosticBulkType.F90
  WaterDiagnosticType.F90
  WaterFluxType.F90
  WaterFluxBulkType.F90
  WaterInfoBaseType.F90
  WaterStateBulkType.F90
  WaterStateType.F90
  )

sourcelist_to_parent(clm_sources)
