module readParamsMod

  !-----------------------------------------------------------------------
  !
  ! Read parameters
  ! module used to read parameters for individual modules and/or for some 
  ! well defined functionality (eg. ED).
  !
  ! ! USES:
  use clm_varctl , only : paramfile, iulog, use_fates, use_cn
  use spmdMod    , only : masterproc
  use fileutils  , only : getfil
  use ncdio_pio  , only : ncd_pio_closefile, ncd_pio_openfile
  use ncdio_pio  , only : file_desc_t , ncd_inqdid, ncd_inqdlen

  implicit none
  private
  !
  public :: readParameters

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readParameters (nutrient_competition_method, photosyns_inst)
    !
    ! ! USES:
    use CNSharedParamsMod                 , only : CNParamsReadShared
    use CNGapMortalityMod                 , only : readCNGapMortParams                    => readParams
    use CNMRespMod                        , only : readCNMRespParams                      => readParams
    use CNFUNMod                          , only : readCNFUNParams                        => readParams
    use CNPhenologyMod                    , only : readCNPhenolParams                     => readParams
    use SoilBiogeochemCompetitionMod      , only : readSoilBiogeochemCompetitionParams    => readParams
    use SoilBiogeochemNLeachingMod        , only : readSoilBiogeochemNLeachingParams      => readParams
    use SoilBiogeochemNitrifDenitrifMod   , only : readSoilBiogeochemNitrifDenitrifParams => readParams
    use SoilBiogeochemLittVertTranspMod   , only : readSoilBiogeochemLittVertTranspParams => readParams
    use SoilBiogeochemPotentialMod        , only : readSoilBiogeochemPotentialParams      => readParams
    use SoilBiogeochemDecompMod           , only : readSoilBiogeochemDecompParams         => readParams
    use SoilBiogeochemDecompCascadeBGCMod , only : readSoilBiogeochemDecompBgcParams      => readParams
    use SoilBiogeochemDecompCascadeCNMod  , only : readSoilBiogeochemDecompCnParams       => readParams
    use ch4Mod                            , only : readCH4Params                          => readParams
    use initVerticalMod                   , only : readParams_initVertical                => readParams
    use LunaMod                           , only : readParams_Luna                        => readParams
    use BareGroundFluxesMod               , only : readParams_BareGroundFluxes            => readParams
    use LakeFluxesMod                     , only : readParams_LakeFluxes                  => readParams
    use CanopyFluxesMod                   , only : readParams_CanopyFluxes                => readParams
    use CanopyTemperatureMod              , only : readParams_CanopyTemperature           => readParams
    use CanopyHydrologyMod                , only : readParams_CanopyHydrology             => readParams
    use SoilHydrologyMod                  , only : readParams_SoilHydrology               => readParams
    use SaturatedExcessRunoffMod          , only : readParams_SaturatedExcessRunoff       => readParams
    use SnowHydrologyMod                  , only : readParams_SnowHydrology               => readParams
    use SurfaceResistanceMod              , only : readParams_SurfaceResistance           => readParams
    use WaterDiagnosticBulkType           , only : readParams_WaterDiagnosticBulk         => readParams
    use NutrientCompetitionMethodMod      , only : nutrient_competition_method_type
    use clm_varctl,                         only : NLFilename_in
    use PhotosynthesisMod                 , only : photosyns_type
    !
    ! !ARGUMENTS:
    type(photosyns_type)                   , intent(in) :: photosyns_inst
    class(nutrient_competition_method_type), intent(in) :: nutrient_competition_method
    !
    ! !LOCAL VARIABLES:
    character(len=256) :: locfn ! local file name
    type(file_desc_t)  :: ncid  ! pio netCDF file id
    integer            :: dimid ! netCDF dimension id
    integer            :: npft  ! number of pfts on pft-physiology file
    character(len=32)  :: subname = 'readParameters'
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'paramMod.F90::'//trim(subname)//' :: reading CLM '//' parameters '
    end if

    call getfil (paramfile, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdid(ncid,'pft',dimid) 
    call ncd_inqdlen(ncid,dimid,npft) 

    !
    ! Above ground biogeochemistry...
    !
    if (use_cn) then
       call nutrient_competition_method%readParams(ncid)
       call readCNGapMortParams(ncid)
       call readCNMRespParams(ncid)
       call readCNFUNParams(ncid)
       call readCNPhenolParams(ncid)
    end if

    !
    ! Soil biogeochemistry...
    !
    if (use_cn .or. use_fates) then
       call readSoilBiogeochemCompetitionParams(ncid)
       call readSoilBiogeochemDecompBgcParams(ncid)
       call readSoilBiogeochemDecompCnParams(ncid)
       call readSoilBiogeochemDecompParams(ncid)
       call readSoilBiogeochemLittVertTranspParams(ncid)
       call readSoilBiogeochemNitrifDenitrifParams(ncid)
       call readSoilBiogeochemNLeachingParams(ncid)
       call readSoilBiogeochemPotentialParams(ncid)
       call CNParamsReadShared(ncid, NLFilename_in)  ! this is called CN params but really is for the soil biogeochem parameters

       call readCH4Params (ncid)
    end if

    !
    ! Main
    !
    call readParams_initVertical ( ncid )

    !
    ! Biogeophysics
    !
    call photosyns_inst%ReadParams( ncid )
    call readParams_Luna ( ncid )
    call readParams_BareGroundFluxes ( ncid )
    call readParams_LakeFluxes ( ncid )
    call readParams_CanopyFluxes ( ncid )
    call readParams_CanopyTemperature ( ncid )
    call readParams_CanopyHydrology ( ncid )
    call readParams_SoilHydrology ( ncid )
    call readParams_SaturatedExcessRunoff ( ncid )
    call readParams_SnowHydrology( ncid )
    call readParams_SurfaceResistance ( ncid )
    call readParams_WaterDiagnosticBulk ( ncid )

    !
    call ncd_pio_closefile(ncid)

  end subroutine readParameters

end module readParamsMod
