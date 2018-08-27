module clm_instMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Instances and definitions of all data types
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use decompMod       , only : bounds_type
  use clm_varpar      , only : ndecomp_pools, nlevdecomp_full
  use clm_varctl      , only : use_cn, use_c13, use_c14, use_lch4, use_cndv, use_fates
  use clm_varctl      , only : use_century_decomp, use_crop
  use clm_varcon      , only : bdsno, c13ratio, c14ratio
  use landunit_varcon , only : istice_mec, istsoil
  use perf_mod        , only : t_startf, t_stopf
  use controlMod      , only : NLFilename

  !-----------------------------------------
  ! Constants
  !-----------------------------------------

  use UrbanParamsType                    , only : urbanparams_type   ! Constants 
  use UrbanParamsType                    , only : IsSimpleBuildTemp, IsProgBuildTemp
  use UrbanTimeVarType                   , only : urbantv_type
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use CNDVType                           , only : dgv_ecophyscon     ! Constants 

  !-----------------------------------------
  ! Definition of component types 
  !-----------------------------------------

  use AerosolMod                      , only : aerosol_type
  use CanopyStateType                 , only : canopystate_type
  use ch4Mod                          , only : ch4_type
  use CNVegetationFacade              , only : cn_vegetation_type
  use SoilBiogeochemStateType         , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonFluxType    , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemCarbonStateType   , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use CropType                        , only : crop_type
  use DryDepVelocity                  , only : drydepvel_type
  use DUSTMod                         , only : dust_type
  use EnergyFluxType                  , only : energyflux_type
  use FrictionVelocityMod             , only : frictionvel_type
  use GlacierSurfaceMassBalanceMod    , only : glacier_smb_type
  use InfiltrationExcessRunoffMod     , only : infiltration_excess_runoff_type
  use IrrigationMod                   , only : irrigation_type
  use LakeStateType                   , only : lakestate_type
  use OzoneBaseMod                    , only : ozone_base_type
  use OzoneFactoryMod                 , only : create_and_init_ozone_type
  use PhotosynthesisMod               , only : photosyns_type
  use SoilHydrologyType               , only : soilhydrology_type
  use SaturatedExcessRunoffMod        , only : saturated_excess_runoff_type
  use SoilStateType                   , only : soilstate_type
  use SolarAbsorbedType               , only : solarabs_type
  use SurfaceRadiationMod             , only : surfrad_type
  use SurfaceAlbedoType               , only : surfalb_type
  use TemperatureType                 , only : temperature_type
  use WaterType                       , only : water_type
  use UrbanParamsType                 , only : urbanparams_type
  use UrbanTimeVarType                , only : urbantv_type
  use HumanIndexMod                   , only : humanindex_type
  use VOCEmissionMod                  , only : vocemis_type
  use CNFireEmissionsMod              , only : fireemis_type
  use atm2lndType                     , only : atm2lnd_type
  use lnd2atmType                     , only : lnd2atm_type
  use lnd2glcMod                      , only : lnd2glc_type 
  use glc2lndMod                      , only : glc2lnd_type
  use glcBehaviorMod                  , only : glc_behavior_type
  use TopoMod                         , only : topo_type
  use GridcellType                    , only : grc
  use LandunitType                    , only : lun                
  use ColumnType                      , only : col                
  use PatchType                       , only : patch                
  use CLMFatesInterfaceMod            , only : hlm_fates_interface_type
  use SoilWaterRetentionCurveMod      , only : soil_water_retention_curve_type
  use NutrientCompetitionMethodMod    , only : nutrient_competition_method_type
  !
  use SoilStateInitTimeConstMod       , only : SoilStateInitTimeConst
  use SoilHydrologyInitTimeConstMod   , only : SoilHydrologyInitTimeConst
  use SurfaceAlbedoMod                , only : SurfaceAlbedoInitTimeConst 
  use LakeCon                         , only : LakeConInit 
  use SoilBiogeochemPrecisionControlMod, only: SoilBiogeochemPrecisionControlInit
  !
  implicit none
  public   ! By default everything is public 
  !
  !-----------------------------------------
  ! Instances of component types
  !-----------------------------------------

  ! Physics types 
  type(aerosol_type)                      :: aerosol_inst
  type(canopystate_type)                  :: canopystate_inst
  type(energyflux_type)                   :: energyflux_inst
  type(frictionvel_type)                  :: frictionvel_inst
  type(glacier_smb_type)                  :: glacier_smb_inst
  type(infiltration_excess_runoff_type)   :: infiltration_excess_runoff_inst
  type(irrigation_type)                   :: irrigation_inst
  type(lakestate_type)                    :: lakestate_inst
  class(ozone_base_type), allocatable     :: ozone_inst
  type(photosyns_type)                    :: photosyns_inst
  type(soilstate_type)                    :: soilstate_inst
  type(soilhydrology_type)                :: soilhydrology_inst
  type(saturated_excess_runoff_type)      :: saturated_excess_runoff_inst
  type(solarabs_type)                     :: solarabs_inst
  type(surfalb_type)                      :: surfalb_inst
  type(surfrad_type)                      :: surfrad_inst
  type(temperature_type)                  :: temperature_inst
  type(urbanparams_type)                  :: urbanparams_inst
  type(urbantv_type)                      :: urbantv_inst
  type(humanindex_type)                   :: humanindex_inst
  type(water_type)                        :: water_inst
  type(atm2lnd_type)                      :: atm2lnd_inst
  type(glc2lnd_type)                      :: glc2lnd_inst
  type(lnd2atm_type)                      :: lnd2atm_inst
  type(lnd2glc_type)                      :: lnd2glc_inst
  type(glc_behavior_type), target         :: glc_behavior
  type(topo_type)                         :: topo_inst
  class(soil_water_retention_curve_type) , allocatable :: soil_water_retention_curve

  ! CN vegetation types
  ! Eventually bgc_vegetation_inst will be an allocatable instance of an abstract
  ! interface
  type(cn_vegetation_type)                :: bgc_vegetation_inst

  class(nutrient_competition_method_type), allocatable :: nutrient_competition_method

  ! Soil biogeochem types 
  type(soilbiogeochem_state_type)         :: soilbiogeochem_state_inst
  type(soilbiogeochem_carbonstate_type)   :: soilbiogeochem_carbonstate_inst
  type(soilbiogeochem_carbonstate_type)   :: c13_soilbiogeochem_carbonstate_inst
  type(soilbiogeochem_carbonstate_type)   :: c14_soilbiogeochem_carbonstate_inst
  type(soilbiogeochem_carbonflux_type)    :: soilbiogeochem_carbonflux_inst
  type(soilbiogeochem_carbonflux_type)    :: c13_soilbiogeochem_carbonflux_inst
  type(soilbiogeochem_carbonflux_type)    :: c14_soilbiogeochem_carbonflux_inst
  type(soilbiogeochem_nitrogenstate_type) :: soilbiogeochem_nitrogenstate_inst
  type(soilbiogeochem_nitrogenflux_type)  :: soilbiogeochem_nitrogenflux_inst

  ! General biogeochem types
  type(ch4_type)                          :: ch4_inst
  type(crop_type)                         :: crop_inst
  type(dust_type)                         :: dust_inst
  type(vocemis_type)                      :: vocemis_inst
  type(fireemis_type)                     :: fireemis_inst
  type(drydepvel_type)                    :: drydepvel_inst

  ! FATES
  type(hlm_fates_interface_type)          :: clm_fates

  !
  public :: clm_instInit       ! Initialize
  public :: clm_instReadNML    ! Read in namelist
  public :: clm_instRest       ! Setup restart
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine clm_instReadNML( NLFilename )
    !
    ! !ARGUMENTS    
    implicit none
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
    ! Read in any namelists that must be read for any clm object instances that need it
    call canopystate_inst%ReadNML( NLFilename )
    call photosyns_inst%ReadNML(   NLFilename )
    if (use_cn .or. use_fates) then
       call crop_inst%ReadNML(     NLFilename )
    end if

  end subroutine clm_instReadNML

  !-----------------------------------------------------------------------
  subroutine clm_instInit(bounds)
    !
    ! !USES: 
    use clm_varpar                         , only : nlevsno, numpft
    use controlMod                         , only : nlfilename, fsurdat
    use domainMod                          , only : ldomain
    use SoilBiogeochemDecompCascadeBGCMod  , only : init_decompcascade_bgc
    use SoilBiogeochemDecompCascadeCNMod   , only : init_decompcascade_cn
    use SoilBiogeochemDecompCascadeContype , only : init_decomp_cascade_constants
    use SoilBiogeochemCompetitionMod       , only : SoilBiogeochemCompetitionInit
    
    use initVerticalMod                    , only : initVertical
    use accumulMod                         , only : print_accum_fields 
    use SoilWaterRetentionCurveFactoryMod  , only : create_soil_water_retention_curve
    use decompMod                          , only : get_proc_bounds
    !
    ! !ARGUMENTS    
    type(bounds_type), intent(in) :: bounds  ! processor bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: c,l,g
    integer               :: nclumps,nc
    integer               :: begp, endp
    integer               :: begc, endc
    integer               :: begl, endl
    type(bounds_type)     :: bounds_clump
    real(r8), allocatable :: h2osno_col(:)
    real(r8), allocatable :: snow_depth_col(:)

    integer :: dummy_to_make_pgi_happy
    !----------------------------------------------------------------------

    ! Note: h2osno_col and snow_depth_col are initialized as local variable 
    ! since they are needed to initialize vertical data structures  

    begp = bounds%begp; endp = bounds%endp 
    begc = bounds%begc; endc = bounds%endc 
    begl = bounds%begl; endl = bounds%endl

    allocate (h2osno_col(begc:endc))
    allocate (snow_depth_col(begc:endc))

    ! snow water
    do c = begc,endc
       l = col%landunit(c)
       g = col%gridcell(c)

       ! In areas that should be snow-covered, it can be problematic to start with 0 snow
       ! cover, because this can affect the long-term state through soil heating, albedo
       ! feedback, etc. On the other hand, we would introduce hysteresis by putting too
       ! much snow in places that are in a net melt regime, because the melt-albedo
       ! feedback may not activate on time (or at all). So, as a compromise, we start with
       ! a small amount of snow in places that are likely to be snow-covered for much or
       ! all of the year.
       if (lun%itype(l)==istice_mec) then
          h2osno_col(c) = 100._r8
       else if (lun%itype(l)==istsoil .and. abs(grc%latdeg(g)) >= 60._r8) then 
          h2osno_col(c) = 100._r8
       else
          h2osno_col(c) = 0._r8
       endif
       snow_depth_col(c)  = h2osno_col(c) / bdsno
    end do

    ! Initialize urban constants

    call urbanparams_inst%Init(bounds)
    call humanindex_inst%Init(bounds)

    ! Initialize urban time varying data
    call urbantv_inst%Init(bounds, NLFilename)

    ! Initialize vertical data components 

    call initVertical(bounds,               &
         glc_behavior, &
         snow_depth_col(begc:endc),              &
         urbanparams_inst%thick_wall(begl:endl), &
         urbanparams_inst%thick_roof(begl:endl))

    ! Initialize clm->drv and drv->clm data structures

    call atm2lnd_inst%Init( bounds, NLFilename )
    call lnd2atm_inst%Init( bounds, NLFilename )

    call glc2lnd_inst%Init( bounds, glc_behavior )
    call lnd2glc_inst%Init( bounds )

    ! Initialization of public data types

    call temperature_inst%Init(bounds,           &
         urbanparams_inst%em_roof(begl:endl),    &
         urbanparams_inst%em_wall(begl:endl),    &
         urbanparams_inst%em_improad(begl:endl), &
         urbanparams_inst%em_perroad(begl:endl), &
         IsSimpleBuildTemp(), IsProgBuildTemp() )

    call canopystate_inst%Init(bounds)

    call soilstate_inst%Init(bounds)
    call SoilStateInitTimeConst(bounds, soilstate_inst, nlfilename) ! sets hydraulic and thermal soil properties

    call water_inst%Init(bounds, &
         h2osno_col = h2osno_col(begc:endc), &
         snow_depth_col = snow_depth_col(begc:endc), &
         watsat_col = soilstate_inst%watsat_col(begc:endc, 1:), &
         t_soisno_col = temperature_inst%t_soisno_col(begc:endc, -nlevsno+1:))

    call glacier_smb_inst%Init(bounds)

    ! COMPILER_BUG(wjs, 2014-11-29, pgi 14.7) Without the following assignment, the
    ! assertion in energyflux_inst%Init fails with pgi 14.7 on yellowstone, presumably due
    ! to a compiler bug.
    dummy_to_make_pgi_happy = ubound(temperature_inst%t_grnd_col, 1)
    call energyflux_inst%Init(bounds, temperature_inst%t_grnd_col(begc:endc), &
         IsSimpleBuildTemp(), IsProgBuildTemp() )

    call aerosol_inst%Init(bounds, NLFilename)

    call frictionvel_inst%Init(bounds)

    call lakestate_inst%Init(bounds)
    call LakeConInit()

    allocate(ozone_inst, source = create_and_init_ozone_type(bounds))

    call photosyns_inst%Init(bounds)

    call soilhydrology_inst%Init(bounds, nlfilename)
    call SoilHydrologyInitTimeConst(bounds, soilhydrology_inst) ! sets time constant properties

    call saturated_excess_runoff_inst%Init(bounds)
    call infiltration_excess_runoff_inst%Init(bounds)

    call solarabs_inst%Init(bounds)

    call surfalb_inst%Init(bounds)
    call SurfaceAlbedoInitTimeConst(bounds)

    call surfrad_inst%Init(bounds)

    call dust_inst%Init(bounds)

    ! Once namelist options are added to control the soil water retention curve method,
    ! we'll need to either pass the namelist file as an argument to this routine, or pass
    ! the namelist value itself (if the namelist is read elsewhere).

    allocate(soil_water_retention_curve, &
         source=create_soil_water_retention_curve())

    call irrigation_inst%init(bounds, nlfilename, soilstate_inst, soil_water_retention_curve)

    call topo_inst%Init(bounds)

    ! Note - always initialize the memory for ch4_inst
    call ch4_inst%Init(bounds, soilstate_inst%cellorg_col(begc:endc, 1:), fsurdat, nlfilename)

    call vocemis_inst%Init(bounds)

    call fireemis_inst%Init(bounds)

    call drydepvel_inst%Init(bounds)

    if (use_cn .or. use_fates ) then

       ! Initialize soilbiogeochem_state_inst

       call soilbiogeochem_state_inst%Init(bounds)

       ! Initialize decompcascade constants
       ! Note that init_decompcascade_bgc and init_decompcascade_cn need 
       ! soilbiogeochem_state_inst to be initialized

       call init_decomp_cascade_constants()
       if (use_century_decomp) then
          call init_decompcascade_bgc(bounds, soilbiogeochem_state_inst, &
                                      soilstate_inst )
       else 
          call init_decompcascade_cn(bounds, soilbiogeochem_state_inst)
       end if

       ! Initalize soilbiogeochem carbon types

       call soilbiogeochem_carbonstate_inst%Init(bounds, carbon_type='c12', ratio=1._r8)
       if (use_c13) then
          call c13_soilbiogeochem_carbonstate_inst%Init(bounds, carbon_type='c13', ratio=c13ratio, &
               c12_soilbiogeochem_carbonstate_inst=soilbiogeochem_carbonstate_inst)
       end if
       if (use_c14) then
          call c14_soilbiogeochem_carbonstate_inst%Init(bounds, carbon_type='c14', ratio=c14ratio, &
               c12_soilbiogeochem_carbonstate_inst=soilbiogeochem_carbonstate_inst)
       end if

       call soilbiogeochem_carbonflux_inst%Init(bounds, carbon_type='c12') 
       if (use_c13) then
          call c13_soilbiogeochem_carbonflux_inst%Init(bounds, carbon_type='c13')
       end if
       if (use_c14) then
          call c14_soilbiogeochem_carbonflux_inst%Init(bounds, carbon_type='c14')
       end if

    end if

    if ( use_cn .or. use_fates) then 

       ! Initalize soilbiogeochem nitrogen types

       call soilbiogeochem_nitrogenstate_inst%Init(bounds, &
            soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools), &
            soilbiogeochem_carbonstate_inst%decomp_cpools_col(begc:endc,1:ndecomp_pools),  &
            soilbiogeochem_carbonstate_inst%decomp_cpools_1m_col(begc:endc, 1:ndecomp_pools))

       call soilbiogeochem_nitrogenflux_inst%Init(bounds) 

       ! Initialize precision control for soil biogeochemistry
       call SoilBiogeochemPrecisionControlInit( soilbiogeochem_carbonstate_inst, c13_soilbiogeochem_carbonstate_inst, &
                                                c14_soilbiogeochem_carbonstate_inst, soilbiogeochem_nitrogenstate_inst)

    end if ! end of if use_cn 

    ! Note - always call Init for bgc_vegetation_inst: some pieces need to be initialized always
    call bgc_vegetation_inst%Init(bounds, nlfilename)

    if (use_cn .or. use_fates) then
       call crop_inst%Init(bounds)
    end if

    
    ! Initialize the Functionaly Assembled Terrestrial Ecosystem Simulator (FATES)
    ! 
    if (use_fates) then
       call clm_fates%Init(bounds)
    end if

    deallocate (h2osno_col)
    deallocate (snow_depth_col)

    ! ------------------------------------------------------------------------
    ! Initialize accumulated fields
    ! ------------------------------------------------------------------------

    ! The time manager needs to be initialized before this called is made, since
    ! the step size is needed. 

    call t_startf('init_accflds')

    call atm2lnd_inst%InitAccBuffer(bounds)

    call temperature_inst%InitAccBuffer(bounds)
    
    call water_inst%InitAccBuffer(bounds)

    call energyflux_inst%InitAccBuffer(bounds)

    call canopystate_inst%InitAccBuffer(bounds)

    call bgc_vegetation_inst%InitAccBuffer(bounds)

    if (use_crop) then
       call crop_inst%InitAccBuffer(bounds)
    end if

    call print_accum_fields()

    call t_stopf('init_accflds')

  end subroutine clm_instInit

  !-----------------------------------------------------------------------
  subroutine clm_instRest(bounds, ncid, flag)
    !
    ! !USES:
    use ncdio_pio       , only : file_desc_t
    use UrbanParamsType , only : IsSimpleBuildTemp, IsProgBuildTemp
    use decompMod       , only : get_proc_bounds, get_proc_clumps, get_clump_bounds

    !
    ! !DESCRIPTION:
    ! Define/write/read CLM restart file.
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds          
    
    type(file_desc_t) , intent(inout) :: ncid ! netcdf id
    character(len=*)  , intent(in)    :: flag ! 'define', 'write', 'read' 

    ! Local variables
    integer                           :: nc, nclumps
    type(bounds_type)                 :: bounds_clump

    !-----------------------------------------------------------------------

    call atm2lnd_inst%restart (bounds, ncid, flag=flag)

    call canopystate_inst%restart (bounds, ncid, flag=flag)

    call energyflux_inst%restart (bounds, ncid, flag=flag, &
         is_simple_buildtemp=IsSimpleBuildTemp(), is_prog_buildtemp=IsProgBuildTemp())

    call frictionvel_inst% restart (bounds, ncid, flag=flag)

    call lakestate_inst%restart (bounds, ncid, flag=flag)

    call ozone_inst%restart (bounds, ncid, flag=flag)

    call photosyns_inst%restart (bounds, ncid, flag=flag)

    call soilhydrology_inst%restart (bounds, ncid, flag=flag)

    call solarabs_inst%restart (bounds, ncid, flag=flag)

    call temperature_inst%restart (bounds, ncid, flag=flag, &
         is_simple_buildtemp=IsSimpleBuildTemp(), is_prog_buildtemp=IsProgBuildTemp())

    call soilstate_inst%restart (bounds, ncid, flag=flag)

    call water_inst%restart(bounds, ncid, flag=flag, &
         watsat_col = soilstate_inst%watsat_col(bounds%begc:bounds%endc,:))

    call irrigation_inst%restart (bounds, ncid, flag=flag)

    call aerosol_inst%restart (bounds, ncid,  flag=flag, &
         h2osoi_ice_col=water_inst%waterstatebulk_inst%h2osoi_ice_col(bounds%begc:bounds%endc,:), &
         h2osoi_liq_col=water_inst%waterstatebulk_inst%h2osoi_liq_col(bounds%begc:bounds%endc,:))

    call surfalb_inst%restart (bounds, ncid, flag=flag, &
         tlai_patch=canopystate_inst%tlai_patch(bounds%begp:bounds%endp), &
         tsai_patch=canopystate_inst%tsai_patch(bounds%begp:bounds%endp))

    call topo_inst%restart (bounds, ncid, flag=flag)

    if (use_lch4) then
       call ch4_inst%restart(bounds, ncid, flag=flag)
    end if

    if ( use_cn ) then
       ! Need to do vegetation restart before soil bgc restart to get totvegc_col for purpose
       ! of resetting soil carbon at exit spinup when no vegetation is growing.
       call bgc_vegetation_inst%restart(bounds, ncid, flag=flag)

       call soilbiogeochem_nitrogenstate_inst%restart(bounds, ncid, flag=flag, &
            totvegc_col=bgc_vegetation_inst%get_totvegc_col(bounds))
       call soilbiogeochem_nitrogenflux_inst%restart(bounds, ncid, flag=flag)

       call crop_inst%restart(bounds, ncid, flag=flag)
    end if

    if (use_cn .or. use_fates) then

       call soilbiogeochem_state_inst%restart(bounds, ncid, flag=flag)
       call soilbiogeochem_carbonstate_inst%restart(bounds, ncid, flag=flag, carbon_type='c12', &
            totvegc_col=bgc_vegetation_inst%get_totvegc_col(bounds))

       if (use_c13) then
          call c13_soilbiogeochem_carbonstate_inst%restart(bounds, ncid, flag=flag, carbon_type='c13', &
               totvegc_col=bgc_vegetation_inst%get_totvegc_col(bounds), &
               c12_soilbiogeochem_carbonstate_inst=soilbiogeochem_carbonstate_inst)
       end if
       if (use_c14) then
          call c14_soilbiogeochem_carbonstate_inst%restart(bounds, ncid, flag=flag, carbon_type='c14', &
               totvegc_col=bgc_vegetation_inst%get_totvegc_col(bounds), &
               c12_soilbiogeochem_carbonstate_inst=soilbiogeochem_carbonstate_inst)
       end if
       call soilbiogeochem_carbonflux_inst%restart(bounds, ncid, flag=flag)
    endif

    if (use_fates) then

       call clm_fates%restart(bounds, ncid, flag=flag,  &
            waterdiagnosticbulk_inst=water_inst%waterdiagnosticbulk_inst, &
            canopystate_inst=canopystate_inst, &
            frictionvel_inst=frictionvel_inst)

    end if

 end subroutine clm_instRest

end module clm_instMod

