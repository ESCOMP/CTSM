module ctsm_ReadParams

  !-----------------------------------------------------------------------
  !
  ! Read parameters
  ! module used to read parameters for individual modules and/or for some 
  ! well defined functionality (eg. ED).
  !
  ! ! USES:
  use ctsm_VarCtl , only : paramfile, iulog, use_fates, use_cn
  use ctsm_Spmd    , only : masterproc
  use ctsm_FileUtils  , only : getfil
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
    use ctsm_CNSharedParamsMod                 , only : CNParamsReadShared
    use ctsm_CNGapMortalityMod                 , only : readCNGapMortParams                    => readParams
    use ctsm_CNMaintRespMod                        , only : readctsm_CNMaintRespParams                      => readParams
    use ctsm_CNFunMod                          , only : readctsm_CNFunParams                        => readParams
    use ctsm_CNPhenologyMod                    , only : readCNPhenolParams                     => readParams
    use ctsm_SoilBiogeochemCompetition      , only : readSoilBiogeochemCompetitionParams    => readParams
    use ctsm_SoilBiogeochemNLeaching        , only : readSoilBiogeochemNLeachingParams      => readParams
    use ctsm_SoilBiogeochemNitrifDenitrif   , only : readSoilBiogeochemNitrifDenitrifParams => readParams
    use ctsm_SoilBiogeochemLittVertTransp   , only : readSoilBiogeochemLittVertTranspParams => readParams
    use ctsm_SoilBiogeochemPotential        , only : readSoilBiogeochemPotentialParams      => readParams
    use ctsm_SoilBiogeochemDecomp           , only : readSoilBiogeochemDecompParams         => readParams
    use ctsm_SoilBiogeochemDecompCascadeBGC , only : readSoilBiogeochemDecompBgcParams      => readParams
    use ctsm_SoilBiogeochemDecompCascadeCN  , only : readSoilBiogeochemDecompCnParams       => readParams
    use ctsm_Methane                            , only : readCH4Params                          => readParams
    use ctsm_Luna                           , only : readParams_Luna                        => readParams
    use ctsm_BareGroundFluxes               , only : readParams_BareGroundFluxes            => readParams
    use ctsm_LakeFluxes                     , only : readParams_LakeFluxes                  => readParams
    use ctsm_CanopyFluxes                   , only : readParams_CanopyFluxes                => readParams
    use ctsm_UrbanFluxes                    , only : readParams_UrbanFluxes                 => readParams
    use ctsm_CanopyHydrology                , only : readParams_CanopyHydrology             => readParams
    use ctsm_SoilHydrology                  , only : readParams_SoilHydrology               => readParams
    use ctsm_SoilStateInitTimeConst         , only : readParams_SoilStateInitTimeConst      => readParams
    use ctsm_SoilWaterMovement              , only : readParams_SoilWaterMovement           => readParams
    use ctsm_SaturatedExcessRunoff          , only : readParams_SaturatedExcessRunoff       => readParams
    use ctsm_InfiltrationExcessRunoff       , only : readParams_InfiltrationExcessRunoff    => readParams
    use ctsm_SurfaceResistance              , only : readParams_SurfaceResistance           => readParams
    use ctsm_WaterDiagnosticBulkType           , only : readParams_WaterDiagnosticBulk         => readParams
    use ctsm_SnowHydrology                  , only : readParams_SnowHydrology               => readParams
    use ctsm_NutrientCompetitionMethodMod      , only : nutrient_competition_method_type
    use ctsm_VarCtl,                         only : NLFilename_in
    use ctsm_Photosynthesis                 , only : photosyns_type
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
       call readctsm_CNMaintRespParams(ncid)
       call readctsm_CNFunParams(ncid)
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
    ! Biogeophysics
    !
    call photosyns_inst%ReadParams( ncid )
    call readParams_Luna ( ncid )
    call readParams_BareGroundFluxes ( ncid )
    call readParams_LakeFluxes ( ncid )
    call readParams_CanopyFluxes ( ncid )
    call readParams_UrbanFluxes ( ncid )
    call readParams_CanopyHydrology ( ncid )
    call readParams_SoilHydrology ( ncid )
    call readParams_SoilStateInitTimeConst ( ncid )
    call readParams_SaturatedExcessRunoff ( ncid )
    call readParams_SoilWaterMovement ( ncid )
    call readParams_InfiltrationExcessRunoff ( ncid )
    call readParams_SurfaceResistance ( ncid )
    call readParams_WaterDiagnosticBulk ( ncid )
    call readParams_SnowHydrology (ncid)

    !
    call ncd_pio_closefile(ncid)

  end subroutine readParameters

end module ctsm_ReadParams
