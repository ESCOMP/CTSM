# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

list(APPEND clm_sources 
  AerosolMod.F90
  DaylengthMod.F90
  BalanceCheckMod.F90
  CanopyStateType.F90
  EnergyFluxType.F90
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
  SolarAbsorbedType.F90
  SurfaceAlbedoType.F90
  TemperatureType.F90
  TotalWaterAndHeatMod.F90
  UrbanParamsType.F90
  WaterfluxType.F90
  WaterStateType.F90
  )

sourcelist_to_parent(clm_sources)
