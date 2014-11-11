module CNAllocationMod
  
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines used in allocation model for coupled carbon
  ! nitrogen code.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg

  use abortutils          , only : endrun
  use decompMod           , only : bounds_type
  use CanopyStateType     , only : canopystate_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use CNStateType         , only : cnstate_type
  use PhotosynthesisType  , only : photosyns_type
  use CropType            , only : crop_type
  use EcophysConType      , only : ecophyscon
  use LandunitType        , only : lun                
  use ColumnType          , only : col                
  use PatchType           , only : pft                
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readCNAllocParams
  public :: CNAllocationInit         ! Initialization
  public :: CNAllocation             ! run method


  !
  !
  ! !PUBLIC DATA MEMBERS:
  character(len=*), parameter, public :: suplnAll='ALL'  ! Supplemental Nitrogen for all PFT's
  character(len=*), parameter, public :: suplnNon='NONE' ! No supplemental Nitrogen
  character(len=15), public :: suplnitro = suplnNon      ! Supplemental Nitrogen mode
  !
  !logical :: crop_supln  = .false.             !Prognostic crop receives supplemental Nitrogen
  
  !above ground
  integer, parameter :: clm_plant_bgc = 0       !
  integer       :: plant_bgc_method             !0: use clm's big leaf plant model
  
  !below ground
  integer, parameter :: clm45_soil_bgc = 0      !
  integer            :: soil_bgc_method         !0: use default soil bgc relased in clm45

  
  !-----------------------------------------------------------------------

contains

  subroutine readCNAllocParams ( ncid , nutrient_competition_method)
  !
  ! !USES:
  use ncdio_pio , only : file_desc_t
  use NutrientCompetitionMethodMod, only : nutrient_competition_method_type  
  
  implicit none
  ! !ARGUMENTS:  
  type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
  class(nutrient_competition_method_type), intent(in) :: nutrient_competition_method
  !
  ! !LOCAL VARIABLES:
  character(len=32)  :: subname = 'CNAllocParamsType'
  character(len=100) :: errCode = '-Error reading in parameters file:'
  logical            :: readv ! has variable been read in or not
  real(r8)           :: tempr ! temporary to read in parameter
  character(len=100) :: tString ! temp. var for reading
    
  call nutrient_competition_method%read_nutrient_competition_params(ncid)
  
  end subroutine readCNAllocParams
  !-----------------------------------------------------------------------
  subroutine CNAllocationInit ( bounds, nutrient_competition_method)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use clm_varcon      , only: secspday
    use clm_time_manager, only: get_step_size
    use clm_varpar      , only: crop_prog
    use clm_varctl      , only: iulog, cnallocate_carbon_only_set
    use shr_infnan_mod  , only: nan => shr_infnan_nan, assignment(=)
    use NutrientCompetitionMethodMod, only : nutrient_competition_method_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    class(nutrient_competition_method_type), intent(in) :: nutrient_competition_method
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'CNAllocationInit'
    logical :: carbon_only
    !-----------------------------------------------------------------------
  

    

    call nutrient_competition_method%init_nutrient_competition(bounds, crop_prog)

    ! Change namelist settings into private logical variables
    select case(suplnitro)
    case(suplnNon)
       Carbon_only = .false.
    case(suplnAll)
       Carbon_only = .true.
    case default
       write(iulog,*) 'Supplemental Nitrogen flag (suplnitro) can only be: ', &
            suplnNon, ' or ', suplnAll
       call endrun(msg='ERROR: supplemental Nitrogen flag is not correct'//&
            errMsg(__FILE__, __LINE__))
    end select

    call cnallocate_carbon_only_set(carbon_only)

  end subroutine CNAllocationInit

  !-----------------------------------------------------------------------
  subroutine CNAllocation (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       photosyns_vars, crop_vars, canopystate_vars, cnstate_vars,             &
       carbonstate_vars, carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars,  &
       nitrogenstate_vars, nitrogenflux_vars, nutrient_competition_method)
    !
    ! !USES:
    use shr_sys_mod      , only: shr_sys_flush
    use clm_varctl       , only: iulog, cnallocate_carbon_only
    use pftvarcon        , only: npcropmin, declfact, bfact, aleaff, arootf, astemf
    use pftvarcon        , only: arooti, fleafi, allconsl, allconss, grperc, grpnow, nsoybean
    use clm_varpar       , only: nlevsoi, nlevdecomp
    use clm_varcon       , only: nitrif_n2o_loss_frac, secspday
    use landunit_varcon  , only: istsoil, istcrop
    use clm_time_manager , only: get_step_size
    use NutrientCompetitionMethodMod, only : nutrient_competition_method_type
    implicit none
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc        ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)  ! filter for soil columns
    integer                  , intent(in)    :: num_soilp        ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)  ! filter for soil patches
    type(photosyns_type)     , intent(in)    :: photosyns_vars
    type(crop_type)          , intent(in)    :: crop_vars
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(carbonstate_type)   , intent(in)    :: carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    class(nutrient_competition_method_type), intent(in) :: nutrient_competition_method
    !Jinyun Tang: at this stage, the plant_nutrient_demand only calculates the plant ntirgeon demand.
    !I assume phosphorus dynamics will be included in the future. Also, I consider plant_nutrient_demand
    !as a generic interface to call actual nutrient calculation from different aboveground plantbgc. Right now
    !it is assumed the plant nutrient demand is summarized into columnwise demand, and the nutrient redistribution
    !after uptake is done by the plant bgc accordingly. When nutrient competition is required to be done at cohort level
    !both plant_nutrient_demand and do_nutrient_competition should be modified, but that modification should not significantly
    !change the current interface.
    
    call nutrient_competition_method%calc_plant_nutrient_demand(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
         photosyns_vars, canopystate_vars, crop_vars, carbonstate_vars, cnstate_vars, &
         carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars, c13_carbonflux_vars, c14_carbonflux_vars)
        
        
    call nutrient_competition_method%calc_nutrient_competition(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
        photosyns_vars,  cnstate_vars, carbonstate_vars,  carbonflux_vars, c13_carbonflux_vars, &
        c14_carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars)
      
      

 end subroutine CNAllocation
 
 
end module CNAllocationMod
