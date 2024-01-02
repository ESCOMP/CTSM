module dynHarvestMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle reading of the harvest data, as well as the state updates that happen as a
  ! result of harvest.
  !
  ! !USES:
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use decompMod               , only : bounds_type, bounds_level_proc
  use abortutils              , only : endrun
  use dynFileMod              , only : dyn_file_type
  use dynVarTimeUninterpMod   , only : dyn_var_time_uninterp_type
  use CNVegCarbonStateType    , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType     , only : cnveg_carbonflux_type
  use CNVegNitrogenStateType  , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType   , only : cnveg_nitrogenflux_type
  use SoilBiogeochemStateType , only : soilbiogeochem_state_type
  use pftconMod               , only : pftcon
  use clm_varcon              , only : grlnd
  use ColumnType              , only : col                
  use PatchType               , only : patch                
  use CNSharedParamsMod       , only : use_matrixcn
  use clm_varctl              , only : use_fates
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  !
  public :: dynHarvest_init    ! initialize data structures for harvest information, used by both FATES and non-FATES/CN
  public :: dynHarvest_interp  ! get harvest data for current time step, if needed, only used by non-FATES/CN
  public :: dynHarvest_interp_resolve_harvesttypes  ! get harvest data for current time step, if needed, harvest-type-resolved.  only used by FATES.
  public :: CNHarvest          ! harvest mortality routine for CN code, only used by non-FATES/CN
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: CNHarvestPftToColumn   ! gather patch-level harvest fluxes to the column level
  !
  ! !PRIVATE TYPES:

  ! Note that, since we have our own dynHarvest_file object (distinct from dynpft_file),
  ! we could theoretically have a separate file providing harvest data from that providing
  ! the pftdyn data
  type(dyn_file_type), target :: dynHarvest_file ! information for the file containing harvest data

  ! Define the underlying harvest variables
  integer, parameter, public :: num_harvest_inst = 5
  character(len=64), parameter, public :: harvest_varnames(num_harvest_inst) = &
       [character(len=64) :: 'HARVEST_VH1', 'HARVEST_VH2', 'HARVEST_SH1', 'HARVEST_SH2', 'HARVEST_SH3']
  
  type(dyn_var_time_uninterp_type) :: harvest_inst(num_harvest_inst)   ! value of each harvest variable

  real(r8) , allocatable   :: harvest(:) ! harvest rates
  logical, private         :: do_harvest ! whether we're in a period when we should do harvest
  character(len=*), parameter, private :: string_not_set = "not_set"  ! string to initialize with to indicate string wasn't set
  character(len=64), public, protected :: harvest_units = string_not_set ! units from harvest variables 
  character(len=64), parameter, public :: mass_units = "gC/m2/yr"
  character(len=64), parameter, public :: unitless_units = "unitless"
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dynHarvest_init(bounds, harvest_filename)
    !
    ! !DESCRIPTION:
    ! Initialize data structures for harvest information.
    ! This should be called once, during model initialization.
    ! 
    ! !USES:
    use dynVarTimeUninterpMod , only : dyn_var_time_uninterp_type
    use dynTimeInfoMod        , only : YEAR_POSITION_START_OF_TIMESTEP
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds           ! proc-level bounds
    character(len=*)  , intent(in) :: harvest_filename ! name of file containing harvest information
    !
    ! !LOCAL VARIABLES:
    integer :: varnum     ! counter for harvest variables
    integer :: num_points ! number of spatial points
    integer :: ier        ! error code
    character(len=64) :: units = string_not_set
    
    character(len=*), parameter :: subname = 'dynHarvest_init'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == bounds_level_proc, subname // ': argument must be PROC-level bounds')

    ! we only need to keep this summary variable in CN veg type
    if ( .not. use_fates ) then  
       allocate(harvest(bounds%begg:bounds%endg),stat=ier)
       if (ier /= 0) then
          call endrun(msg=' allocation error for harvest'//errMsg(sourcefile, __LINE__))
       end if
    end if

    ! Get the year from the START of the timestep for consistency with other dyn file
    ! stuff (though it shouldn't actually matter for harvest, given the way the
    ! once-per-year harvest is applied).
    dynHarvest_file = dyn_file_type(harvest_filename, YEAR_POSITION_START_OF_TIMESTEP)
    
    ! Get initial harvest data
    num_points = (bounds%endg - bounds%begg + 1)
    do varnum = 1, num_harvest_inst
       harvest_inst(varnum) = dyn_var_time_uninterp_type( &
            dyn_file=dynHarvest_file, varname=harvest_varnames(varnum), &
            dim1name=grlnd, conversion_factor=1.0_r8, &
            do_check_sums_equal_1=.false., data_shape=[num_points])
       call harvest_inst(varnum)%get_att("units",units)
       if ( trim(units) == string_not_set ) then
          units = unitless_units
       else if ( trim(units) == unitless_units ) then

       else if ( trim(units) /= mass_units ) then
          call endrun(msg=' bad units read in from file='//trim(units)//errMsg(sourcefile, __LINE__))
       end if
       if ( varnum > 1 .and. trim(units) /= trim(harvest_units) )then
          call endrun(msg=' harvest units are inconsitent on file ='// &
               trim(harvest_filename)//errMsg(sourcefile, __LINE__))
       end if
       harvest_units = units
       units = string_not_set
    end do

  end subroutine dynHarvest_init


  !-----------------------------------------------------------------------
  subroutine dynHarvest_interp(bounds)
    !
    ! !DESCRIPTION:
    ! Get harvest data for model time, when needed.
    !
    ! Note that harvest data are stored as rates (not weights) and so time interpolation
    ! is not necessary - the harvest rate is held constant through the year.  This is
    ! consistent with the treatment of changing PFT weights, where interpolation of the
    ! annual endpoint weights leads to a constant rate of change in PFT weight through the
    ! year, with abrupt changes in the rate at annual boundaries.
    !
    ! !USES:
    use dynTimeInfoMod , only : time_info_type
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: varnum       ! counter for harvest variables
    real(r8), allocatable :: this_data(:) ! data for a single harvest variable

    character(len=*), parameter :: subname = 'dynHarvest_interp'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == bounds_level_proc, subname // ': argument must be PROC-level bounds')

    call dynHarvest_file%time_info%set_current_year()

    ! Get total harvest for this time step
    harvest(bounds%begg:bounds%endg) = 0._r8

    if (dynHarvest_file%time_info%is_before_time_series()) then
       ! Turn off harvest before the start of the harvest time series
       do_harvest = .false.
    else
       ! Note that do_harvest stays true even past the end of the time series. This
       ! means that harvest rates will be maintained at the rate given in the last
       ! year of the file for all years past the end of this specified time series.
       do_harvest = .true.
       allocate(this_data(bounds%begg:bounds%endg))
       do varnum = 1, num_harvest_inst
          call harvest_inst(varnum)%get_current_data(this_data)
          harvest(bounds%begg:bounds%endg) = harvest(bounds%begg:bounds%endg) + &
               this_data(bounds%begg:bounds%endg)
       end do
       deallocate(this_data)
    end if

  end subroutine dynHarvest_interp


  !-----------------------------------------------------------------------
  subroutine dynHarvest_interp_resolve_harvesttypes(bounds, harvest_rates, after_start_of_harvest_ts)
    !
    ! !DESCRIPTION:
    ! Get harvest data for model time, when needed.
    !
    ! Note that harvest data are stored as rates (not weights) and so time interpolation
    ! is not necessary - the harvest rate is held constant through the year.  This is
    ! consistent with the treatment of changing PFT weights, where interpolation of the
    ! annual endpoint weights leads to a constant rate of change in PFT weight through the
    ! year, with abrupt changes in the rate at annual boundaries.
    !
    ! Note the difference between this and dynHarvest_interp is that here, we keep the different
    ! forcing sets distinct (e.g., for passing to FATES which has distinct primary and secondary lands)
    ! and thus store it in harvest_typeresolved
    !
    ! !USES:
    use dynTimeInfoMod , only : time_info_type
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8)         , intent(out):: harvest_rates(bounds%begg: , 1: ) ! output the harvest rates
    logical          , intent(out):: after_start_of_harvest_ts
    !
    ! !LOCAL VARIABLES:
    integer               :: varnum       ! counter for harvest variables
    real(r8), allocatable :: this_data(:) ! data for a single harvest variable

    character(len=*), parameter :: subname = 'dynHarvest_interp'
    !-----------------------------------------------------------------------

    call dynHarvest_file%time_info%set_current_year()

    if (dynHarvest_file%time_info%is_before_time_series()) then
       ! Turn off harvest before the start of the harvest time series
       after_start_of_harvest_ts = .false.
       harvest_rates(bounds%begg:bounds%endg,1:num_harvest_inst) = 0._r8
    else
       ! Note that do_harvest stays true even past the end of the time series. This
       ! means that harvest rates will be maintained at the rate given in the last
       ! year of the file for all years past the end of this specified time series.
       after_start_of_harvest_ts = .true.
       allocate(this_data(bounds%begg:bounds%endg))
       do varnum = 1, num_harvest_inst
          call harvest_inst(varnum)%get_current_data(this_data)
          harvest_rates(bounds%begg:bounds%endg,varnum) = this_data(bounds%begg:bounds%endg)
       end do
       deallocate(this_data)
    end if

  end subroutine dynHarvest_interp_resolve_harvesttypes

  !-----------------------------------------------------------------------
  subroutine CNHarvest (num_soilp, filter_soilp, &
       soilbiogeochem_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, &
       cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! Harvest mortality routine for coupled carbon-nitrogen code (CN)
    !
    ! !USES:
    use pftconMod       , only : noveg, nbrdlf_evr_shrub
    use clm_varcon      , only : secspday
    use clm_time_manager, only : get_step_size_real, is_beg_curr_year
    use CNVegMatrixMod  , only : matrix_update_gmc, matrix_update_gmn
    !
    ! !ARGUMENTS:
    integer                         , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type)  , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p                         ! patch index
    integer :: g                         ! gridcell index
    integer :: fp                        ! patch filter index
    real(r8):: thistreec                 ! carbon in this tree for calculating harvest fraction (gC/m2)
    real(r8):: cm                        ! rate for carbon harvest mortality (gC/m2/yr)
    real(r8):: am                        ! rate for fractional harvest mortality (1/yr)
    real(r8):: m                         ! rate for fractional harvest mortality (1/s)
    real(r8):: dtime                     ! model time step (s)
    !-----------------------------------------------------------------------

    associate(& 
         ivt                                 =>    patch%itype                                                      , & ! Input:  [integer (:)]  pft vegetation type                                
         
         leafc                               =>    cnveg_carbonstate_inst%leafc_patch                             , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C                                    
         frootc                              =>    cnveg_carbonstate_inst%frootc_patch                            , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C                               
         livestemc                           =>    cnveg_carbonstate_inst%livestemc_patch                         , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C                               
         deadstemc                           =>    cnveg_carbonstate_inst%deadstemc_patch                         , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C                               
         livecrootc                          =>    cnveg_carbonstate_inst%livecrootc_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C                        
         deadcrootc                          =>    cnveg_carbonstate_inst%deadcrootc_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C                        
         xsmrpool                            =>    cnveg_carbonstate_inst%xsmrpool_patch                          , & ! Input:  [real(r8) (:)]  (gC/m2) abstract C pool to meet excess MR demand  
         leafc_storage                       =>    cnveg_carbonstate_inst%leafc_storage_patch                     , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage                            
         frootc_storage                      =>    cnveg_carbonstate_inst%frootc_storage_patch                    , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage                       
         livestemc_storage                   =>    cnveg_carbonstate_inst%livestemc_storage_patch                 , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C storage                       
         deadstemc_storage                   =>    cnveg_carbonstate_inst%deadstemc_storage_patch                 , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C storage                       
         livecrootc_storage                  =>    cnveg_carbonstate_inst%livecrootc_storage_patch                , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage                
         deadcrootc_storage                  =>    cnveg_carbonstate_inst%deadcrootc_storage_patch                , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                
         gresp_storage                       =>    cnveg_carbonstate_inst%gresp_storage_patch                     , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration storage                
         leafc_xfer                          =>    cnveg_carbonstate_inst%leafc_xfer_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C transfer                           
         frootc_xfer                         =>    cnveg_carbonstate_inst%frootc_xfer_patch                       , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C transfer                      
         livestemc_xfer                      =>    cnveg_carbonstate_inst%livestemc_xfer_patch                    , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C transfer                      
         deadstemc_xfer                      =>    cnveg_carbonstate_inst%deadstemc_xfer_patch                    , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C transfer                      
         livecrootc_xfer                     =>    cnveg_carbonstate_inst%livecrootc_xfer_patch                   , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C transfer               
         deadcrootc_xfer                     =>    cnveg_carbonstate_inst%deadcrootc_xfer_patch                   , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
         gresp_xfer                          =>    cnveg_carbonstate_inst%gresp_xfer_patch                        , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration transfer               
         
         leafn                               =>    cnveg_nitrogenstate_inst%leafn_patch                           , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N                                    
         frootn                              =>    cnveg_nitrogenstate_inst%frootn_patch                          , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N                               
         livestemn                           =>    cnveg_nitrogenstate_inst%livestemn_patch                       , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N                               
         deadstemn                           =>    cnveg_nitrogenstate_inst%deadstemn_patch                       , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N                               
         livecrootn                          =>    cnveg_nitrogenstate_inst%livecrootn_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N                        
         deadcrootn                          =>    cnveg_nitrogenstate_inst%deadcrootn_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N                        
         retransn                            =>    cnveg_nitrogenstate_inst%retransn_patch                        , & ! Input:  [real(r8) (:)]  (gN/m2) plant pool of retranslocated N            
         leafn_storage                       =>    cnveg_nitrogenstate_inst%leafn_storage_patch                   , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N storage                            
         frootn_storage                      =>    cnveg_nitrogenstate_inst%frootn_storage_patch                  , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N storage                       
         livestemn_storage                   =>    cnveg_nitrogenstate_inst%livestemn_storage_patch               , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N storage                       
         deadstemn_storage                   =>    cnveg_nitrogenstate_inst%deadstemn_storage_patch               , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N storage                       
         livecrootn_storage                  =>    cnveg_nitrogenstate_inst%livecrootn_storage_patch              , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N storage                
         deadcrootn_storage                  =>    cnveg_nitrogenstate_inst%deadcrootn_storage_patch              , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N storage                
         leafn_xfer                          =>    cnveg_nitrogenstate_inst%leafn_xfer_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N transfer                           
         frootn_xfer                         =>    cnveg_nitrogenstate_inst%frootn_xfer_patch                     , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N transfer                      
         livestemn_xfer                      =>    cnveg_nitrogenstate_inst%livestemn_xfer_patch                  , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N transfer                      
         deadstemn_xfer                      =>    cnveg_nitrogenstate_inst%deadstemn_xfer_patch                  , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N transfer                      
         livecrootn_xfer                     =>    cnveg_nitrogenstate_inst%livecrootn_xfer_patch                 , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N transfer               
         deadcrootn_xfer                     =>    cnveg_nitrogenstate_inst%deadcrootn_xfer_patch                 , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N transfer               
         
         hrv_leafc_to_litter                 =>    cnveg_carbonflux_inst%hrv_leafc_to_litter_patch                , & ! Output: [real(r8) (:)]                                                    
         hrv_frootc_to_litter                =>    cnveg_carbonflux_inst%hrv_frootc_to_litter_patch               , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemc_to_litter             =>    cnveg_carbonflux_inst%hrv_livestemc_to_litter_patch            , & ! Output: [real(r8) (:)]                                                    
         wood_harvestc                       =>    cnveg_carbonflux_inst%wood_harvestc_patch                      , & ! Output: [real(r8) (:)]
         hrv_livecrootc_to_litter            =>    cnveg_carbonflux_inst%hrv_livecrootc_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootc_to_litter            =>    cnveg_carbonflux_inst%hrv_deadcrootc_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         hrv_xsmrpool_to_atm                 =>    cnveg_carbonflux_inst%hrv_xsmrpool_to_atm_patch                , & ! Output: [real(r8) (:)]                                                    
         hrv_leafc_storage_to_litter         =>    cnveg_carbonflux_inst%hrv_leafc_storage_to_litter_patch        , & ! Output: [real(r8) (:)]                                                    
         hrv_frootc_storage_to_litter        =>    cnveg_carbonflux_inst%hrv_frootc_storage_to_litter_patch       , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemc_storage_to_litter     =>    cnveg_carbonflux_inst%hrv_livestemc_storage_to_litter_patch    , & ! Output: [real(r8) (:)]                                                    
         hrv_deadstemc_storage_to_litter     =>    cnveg_carbonflux_inst%hrv_deadstemc_storage_to_litter_patch    , & ! Output: [real(r8) (:)]                                                    
         hrv_livecrootc_storage_to_litter    =>    cnveg_carbonflux_inst%hrv_livecrootc_storage_to_litter_patch   , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootc_storage_to_litter    =>    cnveg_carbonflux_inst%hrv_deadcrootc_storage_to_litter_patch   , & ! Output: [real(r8) (:)]                                                    
         hrv_gresp_storage_to_litter         =>    cnveg_carbonflux_inst%hrv_gresp_storage_to_litter_patch        , & ! Output: [real(r8) (:)]                                                    
         hrv_leafc_xfer_to_litter            =>    cnveg_carbonflux_inst%hrv_leafc_xfer_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         hrv_frootc_xfer_to_litter           =>    cnveg_carbonflux_inst%hrv_frootc_xfer_to_litter_patch          , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemc_xfer_to_litter        =>    cnveg_carbonflux_inst%hrv_livestemc_xfer_to_litter_patch       , & ! Output: [real(r8) (:)]                                                    
         hrv_deadstemc_xfer_to_litter        =>    cnveg_carbonflux_inst%hrv_deadstemc_xfer_to_litter_patch       , & ! Output: [real(r8) (:)]                                                    
         hrv_livecrootc_xfer_to_litter       =>    cnveg_carbonflux_inst%hrv_livecrootc_xfer_to_litter_patch      , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootc_xfer_to_litter       =>    cnveg_carbonflux_inst%hrv_deadcrootc_xfer_to_litter_patch      , & ! Output: [real(r8) (:)]                                                    
         hrv_gresp_xfer_to_litter            =>    cnveg_carbonflux_inst%hrv_gresp_xfer_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         
         hrv_leafn_to_litter                 =>    cnveg_nitrogenflux_inst%hrv_leafn_to_litter_patch              , & ! Output: [real(r8) (:)]                                                    
         hrv_frootn_to_litter                =>    cnveg_nitrogenflux_inst%hrv_frootn_to_litter_patch             , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemn_to_litter             =>    cnveg_nitrogenflux_inst%hrv_livestemn_to_litter_patch          , & ! Output: [real(r8) (:)]                                                    
         wood_harvestn                       =>    cnveg_nitrogenflux_inst%wood_harvestn_patch                    , & ! Output: [real(r8) (:)]
         hrv_livecrootn_to_litter            =>    cnveg_nitrogenflux_inst%hrv_livecrootn_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootn_to_litter            =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
         hrv_retransn_to_litter              =>    cnveg_nitrogenflux_inst%hrv_retransn_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         hrv_leafn_storage_to_litter         =>    cnveg_nitrogenflux_inst%hrv_leafn_storage_to_litter_patch      , & ! Output: [real(r8) (:)]                                                    
         hrv_frootn_storage_to_litter        =>    cnveg_nitrogenflux_inst%hrv_frootn_storage_to_litter_patch     , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemn_storage_to_litter     =>    cnveg_nitrogenflux_inst%hrv_livestemn_storage_to_litter_patch  , & ! Output: [real(r8) (:)]                                                    
         hrv_deadstemn_storage_to_litter     =>    cnveg_nitrogenflux_inst%hrv_deadstemn_storage_to_litter_patch  , & ! Output: [real(r8) (:)]                                                    
         hrv_livecrootn_storage_to_litter    =>    cnveg_nitrogenflux_inst%hrv_livecrootn_storage_to_litter_patch , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootn_storage_to_litter    =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_storage_to_litter_patch , & ! Output: [real(r8) (:)]                                                    
         hrv_leafn_xfer_to_litter            =>    cnveg_nitrogenflux_inst%hrv_leafn_xfer_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
         hrv_frootn_xfer_to_litter           =>    cnveg_nitrogenflux_inst%hrv_frootn_xfer_to_litter_patch        , & ! Output: [real(r8) (:)]                                                    
         hrv_livestemn_xfer_to_litter        =>    cnveg_nitrogenflux_inst%hrv_livestemn_xfer_to_litter_patch     , & ! Output: [real(r8) (:)]                                                    
         hrv_deadstemn_xfer_to_litter        =>    cnveg_nitrogenflux_inst%hrv_deadstemn_xfer_to_litter_patch     , & ! Output: [real(r8) (:)]                                                    
         hrv_livecrootn_xfer_to_litter       =>    cnveg_nitrogenflux_inst%hrv_livecrootn_xfer_to_litter_patch    , & ! Output: [real(r8) (:)]                                                    
         hrv_deadcrootn_xfer_to_litter       =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_xfer_to_litter_patch    , & ! Output: [real(r8) (:)]                                                    
         ileaf_to_iout_gmc                   =>    cnveg_carbonflux_inst%ileaf_to_iout_gm                         , &
         ileafst_to_iout_gmc                 =>    cnveg_carbonflux_inst%ileafst_to_iout_gm                       , &
         ileafxf_to_iout_gmc                 =>    cnveg_carbonflux_inst%ileafxf_to_iout_gm                       , &
         ifroot_to_iout_gmc                  =>    cnveg_carbonflux_inst%ifroot_to_iout_gm                        , &
         ifrootst_to_iout_gmc                =>    cnveg_carbonflux_inst%ifrootst_to_iout_gm                      , &
         ifrootxf_to_iout_gmc                =>    cnveg_carbonflux_inst%ifrootxf_to_iout_gm                      , &
         ilivestem_to_iout_gmc               =>    cnveg_carbonflux_inst%ilivestem_to_iout_gm                     , &
         ilivestemst_to_iout_gmc             =>    cnveg_carbonflux_inst%ilivestemst_to_iout_gm                   , &
         ilivestemxf_to_iout_gmc             =>    cnveg_carbonflux_inst%ilivestemxf_to_iout_gm                   , &
         ideadstem_to_iout_gmc               =>    cnveg_carbonflux_inst%ideadstem_to_iout_gm                     , &
         ideadstemst_to_iout_gmc             =>    cnveg_carbonflux_inst%ideadstemst_to_iout_gm                   , &
         ideadstemxf_to_iout_gmc             =>    cnveg_carbonflux_inst%ideadstemxf_to_iout_gm                   , &
         ilivecroot_to_iout_gmc              =>    cnveg_carbonflux_inst%ilivecroot_to_iout_gm                    , &
         ilivecrootst_to_iout_gmc            =>    cnveg_carbonflux_inst%ilivecrootst_to_iout_gm                  , &
         ilivecrootxf_to_iout_gmc            =>    cnveg_carbonflux_inst%ilivecrootxf_to_iout_gm                  , &
         ideadcroot_to_iout_gmc              =>    cnveg_carbonflux_inst%ideadcroot_to_iout_gm                    , &
         ideadcrootst_to_iout_gmc            =>    cnveg_carbonflux_inst%ideadcrootst_to_iout_gm                  , &
         ideadcrootxf_to_iout_gmc            =>    cnveg_carbonflux_inst%ideadcrootxf_to_iout_gm                  , &
         ileaf_to_iout_gmn                   =>    cnveg_nitrogenflux_inst%ileaf_to_iout_gm                       , &
         ileafst_to_iout_gmn                 =>    cnveg_nitrogenflux_inst%ileafst_to_iout_gm                     , &
         ileafxf_to_iout_gmn                 =>    cnveg_nitrogenflux_inst%ileafxf_to_iout_gm                     , &
         ifroot_to_iout_gmn                  =>    cnveg_nitrogenflux_inst%ifroot_to_iout_gm                      , &
         ifrootst_to_iout_gmn                =>    cnveg_nitrogenflux_inst%ifrootst_to_iout_gm                    , &
         ifrootxf_to_iout_gmn                =>    cnveg_nitrogenflux_inst%ifrootxf_to_iout_gm                    , &
         ilivestem_to_iout_gmn               =>    cnveg_nitrogenflux_inst%ilivestem_to_iout_gm                   , &
         ilivestemst_to_iout_gmn             =>    cnveg_nitrogenflux_inst%ilivestemst_to_iout_gm                 , &
         ilivestemxf_to_iout_gmn             =>    cnveg_nitrogenflux_inst%ilivestemxf_to_iout_gm                 , &
         ideadstem_to_iout_gmn               =>    cnveg_nitrogenflux_inst%ideadstem_to_iout_gm                   , &
         ideadstemst_to_iout_gmn             =>    cnveg_nitrogenflux_inst%ideadstemst_to_iout_gm                 , &
         ideadstemxf_to_iout_gmn             =>    cnveg_nitrogenflux_inst%ideadstemxf_to_iout_gm                 , &
         ilivecroot_to_iout_gmn              =>    cnveg_nitrogenflux_inst%ilivecroot_to_iout_gm                  , &
         ilivecrootst_to_iout_gmn            =>    cnveg_nitrogenflux_inst%ilivecrootst_to_iout_gm                , &
         ilivecrootxf_to_iout_gmn            =>    cnveg_nitrogenflux_inst%ilivecrootxf_to_iout_gm                , &
         ideadcroot_to_iout_gmn              =>    cnveg_nitrogenflux_inst%ideadcroot_to_iout_gm                  , &
         ideadcrootst_to_iout_gmn            =>    cnveg_nitrogenflux_inst%ideadcrootst_to_iout_gm                , &
         ideadcrootxf_to_iout_gmn            =>    cnveg_nitrogenflux_inst%ideadcrootxf_to_iout_gm                , &
         iretransn_to_iout_gmn               =>    cnveg_nitrogenflux_inst%iretransn_to_iout_gm                     &
         )

      dtime = get_step_size_real()

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         g = patch%gridcell(p)

         ! If this is a tree pft, then
         ! get the annual harvest "mortality" rate (am) from harvest array
         ! and convert to rate per second
         if (ivt(p) > noveg .and. ivt(p) < nbrdlf_evr_shrub) then

            if (do_harvest) then
               if (harvest_units == mass_units) then
                  thistreec = leafc(p) + frootc(p) + livestemc(p) + deadstemc(p) + livecrootc(p) + deadcrootc(p) + xsmrpool(p)
                  cm = harvest(g)
                  if (thistreec > 0.0_r8) then
                     am = min(0.98_r8,cm/thistreec)    ! Only harvest up to 98% so regrowth is possible PJL
                  else
                     am = 0._r8
                  end if
               else
                  am = harvest(g)
               end if

               ! Apply all harvest at the start of the year
               if (is_beg_curr_year()) then
                  m  = am/dtime
               else
                  m = 0._r8
               end if
            else
               m = 0._r8
            end if

            if(.not. use_matrixcn)then
               ! patch-level harvest carbon fluxes
               ! displayed pools
               hrv_leafc_to_litter(p)               = leafc(p)               * m
               hrv_frootc_to_litter(p)              = frootc(p)              * m
               hrv_livestemc_to_litter(p)           = livestemc(p)           * m
               wood_harvestc(p)                     = deadstemc(p)           * m
               hrv_livecrootc_to_litter(p)          = livecrootc(p)          * m
               hrv_deadcrootc_to_litter(p)          = deadcrootc(p)          * m
               hrv_xsmrpool_to_atm(p)               = xsmrpool(p)            * m

               ! storage pools
               hrv_leafc_storage_to_litter(p)       = leafc_storage(p)       * m
               hrv_frootc_storage_to_litter(p)      = frootc_storage(p)      * m
               hrv_livestemc_storage_to_litter(p)   = livestemc_storage(p)   * m
               hrv_deadstemc_storage_to_litter(p)   = deadstemc_storage(p)   * m
               hrv_livecrootc_storage_to_litter(p)  = livecrootc_storage(p)  * m
               hrv_deadcrootc_storage_to_litter(p)  = deadcrootc_storage(p)  * m
               hrv_gresp_storage_to_litter(p)       = gresp_storage(p)       * m

               ! transfer pools
               hrv_leafc_xfer_to_litter(p)          = leafc_xfer(p)          * m
               hrv_frootc_xfer_to_litter(p)         = frootc_xfer(p)         * m
               hrv_livestemc_xfer_to_litter(p)      = livestemc_xfer(p)      * m
               hrv_deadstemc_xfer_to_litter(p)      = deadstemc_xfer(p)      * m
               hrv_livecrootc_xfer_to_litter(p)     = livecrootc_xfer(p)     * m
               hrv_deadcrootc_xfer_to_litter(p)     = deadcrootc_xfer(p)     * m
               hrv_gresp_xfer_to_litter(p)          = gresp_xfer(p)          * m

               ! patch-level harvest mortality nitrogen fluxes
               ! displayed pools
               hrv_leafn_to_litter(p)               = leafn(p)               * m
               hrv_frootn_to_litter(p)              = frootn(p)              * m
               hrv_livestemn_to_litter(p)           = livestemn(p)           * m
               wood_harvestn(p)                     = deadstemn(p)           * m
               hrv_livecrootn_to_litter(p)          = livecrootn(p)          * m
               hrv_deadcrootn_to_litter(p)          = deadcrootn(p)          * m
               hrv_retransn_to_litter(p)            = retransn(p)            * m

               ! storage pools
               hrv_leafn_storage_to_litter(p)       = leafn_storage(p)       * m
               hrv_frootn_storage_to_litter(p)      = frootn_storage(p)      * m
               hrv_livestemn_storage_to_litter(p)   = livestemn_storage(p)   * m
               hrv_deadstemn_storage_to_litter(p)   = deadstemn_storage(p)   * m
               hrv_livecrootn_storage_to_litter(p)  = livecrootn_storage(p)  * m
               hrv_deadcrootn_storage_to_litter(p)  = deadcrootn_storage(p)  * m

               ! transfer pools
               hrv_leafn_xfer_to_litter(p)          = leafn_xfer(p)          * m
               hrv_frootn_xfer_to_litter(p)         = frootn_xfer(p)         * m
               hrv_livestemn_xfer_to_litter(p)      = livestemn_xfer(p)      * m
               hrv_deadstemn_xfer_to_litter(p)      = deadstemn_xfer(p)      * m
               hrv_livecrootn_xfer_to_litter(p)     = livecrootn_xfer(p)     * m
               hrv_deadcrootn_xfer_to_litter(p)     = deadcrootn_xfer(p)     * m

            ! NOTE: The non-matrix part of this update is in CNCStatUpdate2 CStateUpdate2h (EBK 11/25/2019)
            !   and for Nitrogen The non-matrix part of this update is in CNNStatUpdate2 NStateUpdate2h (EBK 11/25/2019)
            else
               ! patch-level harvest carbon fluxes
               ! displayed pools
               hrv_leafc_to_litter(p)      = matrix_update_gmc(p,ileaf_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.)      * &
                                             leafc(p)
               hrv_frootc_to_litter(p)     = matrix_update_gmc(p,ifroot_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.)     * &
                                             frootc(p)
               hrv_livestemc_to_litter(p)  = matrix_update_gmc(p,ilivestem_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.)  * &
                                             livestemc(p)
               wood_harvestc(p)            = matrix_update_gmc(p,ideadstem_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.)  * &
                                             deadstemc(p)
               hrv_livecrootc_to_litter(p) = matrix_update_gmc(p,ilivecroot_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             livecrootc(p)
               hrv_deadcrootc_to_litter(p) = matrix_update_gmc(p,ideadcroot_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             deadcrootc(p)
               hrv_xsmrpool_to_atm(p)      = xsmrpool(p)            * m

               ! storage pools
               hrv_leafc_storage_to_litter(p)      = matrix_update_gmc(p,ileafst_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.)      * &
                                                     leafc_storage(p)
               hrv_frootc_storage_to_litter(p)     = matrix_update_gmc(p,ifrootst_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.)     * &
                                                     frootc_storage(p)
               hrv_livestemc_storage_to_litter(p)  = matrix_update_gmc(p,ilivestemst_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.)  * &
                                                     livestemc_storage(p)
               hrv_deadstemc_storage_to_litter(p)  = matrix_update_gmc(p,ideadstemst_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.)  * &
                                                     deadstemc_storage(p)
               hrv_livecrootc_storage_to_litter(p) = matrix_update_gmc(p,ilivecrootst_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                                     livecrootc_storage(p)
               hrv_deadcrootc_storage_to_litter(p) = matrix_update_gmc(p,ideadcrootst_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                                     deadcrootc_storage(p)
               hrv_gresp_storage_to_litter(p)      = gresp_storage(p)       * m

               ! transfer pools
               hrv_leafc_xfer_to_litter(p)      = matrix_update_gmc(p,ileafxf_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.)      * &
                                                  leafc_xfer(p)
               hrv_frootc_xfer_to_litter(p)     = matrix_update_gmc(p,ifrootxf_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.)     * &
                                                  frootc_xfer(p)
               hrv_livestemc_xfer_to_litter(p)  = matrix_update_gmc(p,ilivestemxf_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.)  * &
                                                  livestemc_xfer(p)
               hrv_deadstemc_xfer_to_litter(p)  = matrix_update_gmc(p,ideadstemxf_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.)  * &
                                                  deadstemc_xfer(p)
               hrv_livecrootc_xfer_to_litter(p) = matrix_update_gmc(p,ilivecrootxf_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                                  livecrootc_xfer(p)
               hrv_deadcrootc_xfer_to_litter(p) = matrix_update_gmc(p,ideadcrootxf_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                                  deadcrootc_xfer(p)
               hrv_gresp_xfer_to_litter(p)      = gresp_xfer(p)          * m

               ! patch-level harvest mortality nitrogen fluxes
               ! displayed pools
               hrv_leafn_to_litter(p)      = matrix_update_gmn(p,ileaf_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.)      * &
                                             leafn(p)
               hrv_frootn_to_litter(p)     = matrix_update_gmn(p,ifroot_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.)     * &
                                             frootn(p)
               hrv_livestemn_to_litter(p)  = matrix_update_gmn(p,ilivestem_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.)  * &
                                             livestemn(p)
               wood_harvestn(p)            = matrix_update_gmn(p,ideadstem_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.)  * &
                                             deadstemn(p)
               hrv_livecrootn_to_litter(p) = matrix_update_gmn(p,ilivecroot_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             livecrootn(p)
               hrv_deadcrootn_to_litter(p) = matrix_update_gmn(p,ideadcroot_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             deadcrootn(p)
               hrv_retransn_to_litter(p)   = matrix_update_gmn(p,iretransn_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             retransn(p)

               ! storage pools
               hrv_leafn_storage_to_litter(p)      = matrix_update_gmn(p,ileafst_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.)     *&
                                                     leafn_storage(p)
               hrv_frootn_storage_to_litter(p)     = matrix_update_gmn(p,ifrootst_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.)    *&
                                                     frootn_storage(p)
               hrv_livestemn_storage_to_litter(p)  = matrix_update_gmn(p,ilivestemst_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) *&
                                                     livestemn_storage(p)
               hrv_deadstemn_storage_to_litter(p)  = matrix_update_gmn(p,ideadstemst_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) *&
                                                     deadstemn_storage(p)
               hrv_livecrootn_storage_to_litter(p) = matrix_update_gmn(p,ilivecrootst_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.)*&
                                                     livecrootn_storage(p)
               hrv_deadcrootn_storage_to_litter(p) = matrix_update_gmn(p,ideadcrootst_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.)*&
                                                     deadcrootn_storage(p)
               ! transfer pools
               hrv_leafn_xfer_to_litter(p)      = matrix_update_gmn(p,ileafxf_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.)      * &
                                                  leafn_xfer(p)
               hrv_frootn_xfer_to_litter(p)     = matrix_update_gmn(p,ifrootxf_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.)     * &
                                                  frootn_xfer(p)
               hrv_livestemn_xfer_to_litter(p)  = matrix_update_gmn(p,ilivestemxf_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.)  * &
                                                  livestemn_xfer(p)
               hrv_deadstemn_xfer_to_litter(p)  = matrix_update_gmn(p,ideadstemxf_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.)  * &
                                                  deadstemn_xfer(p)
               hrv_livecrootn_xfer_to_litter(p) = matrix_update_gmn(p,ilivecrootxf_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                                  livecrootn_xfer(p)
               hrv_deadcrootn_xfer_to_litter(p) = matrix_update_gmn(p,ideadcrootxf_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                                  deadcrootn_xfer(p)
            end if
         end if  ! end tree block
      end do ! end of pft loop

      ! gather all patch-level litterfall fluxes from harvest to the column
      ! for litter C and N inputs

      call CNHarvestPftToColumn(num_soilp, filter_soilp, &
           soilbiogeochem_state_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)

    end associate 

  end subroutine CNHarvest

 !-----------------------------------------------------------------------
 subroutine CNHarvestPftToColumn (num_soilp, filter_soilp, &
      soilbiogeochem_state_inst, CNVeg_carbonflux_inst, cnveg_nitrogenflux_inst)
   !
   ! !DESCRIPTION:
   ! called at the end of CNHarvest to gather all patch-level harvest litterfall fluxes
   ! to the column level and assign them to the three litter pools
   !
   ! !USES:
   use clm_varpar , only : nlevdecomp, maxsoil_patches, i_litr_min, i_litr_max, i_met_lit
   !
   ! !ARGUMENTS:
   integer                         , intent(in)    :: num_soilp       ! number of soil patches in filter
   integer                         , intent(in)    :: filter_soilp(:) ! patch filter for soil points
   type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
   type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
   type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
   !
   ! !LOCAL VARIABLES:
   integer :: fp,c,p,j,i  ! indices
   !-----------------------------------------------------------------------

   associate(                                                                                                   & 
        ivt                              =>    patch%itype                                                      , & ! Input:  [integer  (:)   ]  pft vegetation type                                
        wtcol                            =>    patch%wtcol                                                      , & ! Input:  [real(r8) (:)   ]  pft weight relative to column (0-1)               
        
        lf_f                             =>    pftcon%lf_f                                                    , & ! Input:  leaf litter fraction
        fr_f                             =>    pftcon%fr_f                                                    , & ! Input:  fine root litter fraction
        
        leaf_prof                        =>    soilbiogeochem_state_inst%leaf_prof_patch                      , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
        froot_prof                       =>    soilbiogeochem_state_inst%froot_prof_patch                     , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
        croot_prof                       =>    soilbiogeochem_state_inst%croot_prof_patch                     , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
        stem_prof                        =>    soilbiogeochem_state_inst%stem_prof_patch                      , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          
        
        hrv_leafc_to_litter              =>    cnveg_carbonflux_inst%hrv_leafc_to_litter_patch                , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootc_to_litter             =>    cnveg_carbonflux_inst%hrv_frootc_to_litter_patch               , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemc_to_litter          =>    cnveg_carbonflux_inst%hrv_livestemc_to_litter_patch            , & ! Input:  [real(r8) (:)   ]                                                    
        pwood_harvestc                   =>    cnveg_carbonflux_inst%wood_harvestc_patch                      , & ! Input:  [real(r8) (:)   ]
        hrv_livecrootc_to_litter         =>    cnveg_carbonflux_inst%hrv_livecrootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootc_to_litter         =>    cnveg_carbonflux_inst%hrv_deadcrootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafc_storage_to_litter      =>    cnveg_carbonflux_inst%hrv_leafc_storage_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootc_storage_to_litter     =>    cnveg_carbonflux_inst%hrv_frootc_storage_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemc_storage_to_litter  =>    cnveg_carbonflux_inst%hrv_livestemc_storage_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemc_storage_to_litter  =>    cnveg_carbonflux_inst%hrv_deadstemc_storage_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootc_storage_to_litter =>    cnveg_carbonflux_inst%hrv_livecrootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootc_storage_to_litter =>    cnveg_carbonflux_inst%hrv_deadcrootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_gresp_storage_to_litter      =>    cnveg_carbonflux_inst%hrv_gresp_storage_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafc_xfer_to_litter         =>    cnveg_carbonflux_inst%hrv_leafc_xfer_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootc_xfer_to_litter        =>    cnveg_carbonflux_inst%hrv_frootc_xfer_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemc_xfer_to_litter     =>    cnveg_carbonflux_inst%hrv_livestemc_xfer_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemc_xfer_to_litter     =>    cnveg_carbonflux_inst%hrv_deadstemc_xfer_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootc_xfer_to_litter    =>    cnveg_carbonflux_inst%hrv_livecrootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootc_xfer_to_litter    =>    cnveg_carbonflux_inst%hrv_deadcrootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_gresp_xfer_to_litter         =>    cnveg_carbonflux_inst%hrv_gresp_xfer_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        cwood_harvestc                   =>    cnveg_carbonflux_inst%wood_harvestc_col                        , & ! InOut:  [real(r8) (:)   ]
        harvest_c_to_litr_c              =>    cnveg_carbonflux_inst%harvest_c_to_litr_c_col                  , & ! InOut:  [real(r8) (:,:,:) ]  C fluxes associated with harvest to litter pools (gC/m3/s)
        harvest_c_to_cwdc                =>    cnveg_carbonflux_inst%harvest_c_to_cwdc_col                    , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with harvest to CWD pool (gC/m3/s)
        
        hrv_leafn_to_litter              =>    cnveg_nitrogenflux_inst%hrv_leafn_to_litter_patch              , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootn_to_litter             =>    cnveg_nitrogenflux_inst%hrv_frootn_to_litter_patch             , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemn_to_litter          =>    cnveg_nitrogenflux_inst%hrv_livestemn_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
        pwood_harvestn                   =>    cnveg_nitrogenflux_inst%wood_harvestn_patch                    , & ! Input:  [real(r8) (:)   ]
        hrv_livecrootn_to_litter         =>    cnveg_nitrogenflux_inst%hrv_livecrootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootn_to_litter         =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_retransn_to_litter           =>    cnveg_nitrogenflux_inst%hrv_retransn_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafn_storage_to_litter      =>    cnveg_nitrogenflux_inst%hrv_leafn_storage_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootn_storage_to_litter     =>    cnveg_nitrogenflux_inst%hrv_frootn_storage_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemn_storage_to_litter  =>    cnveg_nitrogenflux_inst%hrv_livestemn_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemn_storage_to_litter  =>    cnveg_nitrogenflux_inst%hrv_deadstemn_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootn_storage_to_litter =>    cnveg_nitrogenflux_inst%hrv_livecrootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootn_storage_to_litter =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_leafn_xfer_to_litter         =>    cnveg_nitrogenflux_inst%hrv_leafn_xfer_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_frootn_xfer_to_litter        =>    cnveg_nitrogenflux_inst%hrv_frootn_xfer_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livestemn_xfer_to_litter     =>    cnveg_nitrogenflux_inst%hrv_livestemn_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadstemn_xfer_to_litter     =>    cnveg_nitrogenflux_inst%hrv_deadstemn_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_livecrootn_xfer_to_litter    =>    cnveg_nitrogenflux_inst%hrv_livecrootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
        hrv_deadcrootn_xfer_to_litter    =>    cnveg_nitrogenflux_inst%hrv_deadcrootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
        cwood_harvestn                   =>    cnveg_nitrogenflux_inst%wood_harvestn_col                      , & ! InOut:  [real(r8) (:)   ]
        harvest_n_to_litr_n              =>    cnveg_nitrogenflux_inst%harvest_n_to_litr_n_col                , & ! InOut:  [real(r8) (:,:,:)]  N fluxes associated with harvest to litter pools (gN/m3/s)
        harvest_n_to_cwdn                =>    cnveg_nitrogenflux_inst%harvest_n_to_cwdn_col                    & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with harvest to CWD pool (gN/m3/s)
        )

     do j = 1, nlevdecomp
        do fp = 1,num_soilp
           p = filter_soilp(fp)
           c = patch%column(p)

           do i = i_litr_min, i_litr_max
              ! leaf harvest mortality carbon fluxes
              harvest_c_to_litr_c(c,j,i) = &
                 harvest_c_to_litr_c(c,j,i) + &
                 hrv_leafc_to_litter(p) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j)

              ! fine root harvest mortality carbon fluxes
              harvest_c_to_litr_c(c,j,i) = &
                 harvest_c_to_litr_c(c,j,i) + &
                 hrv_frootc_to_litter(p) * fr_f(ivt(p),i) * wtcol(p) * froot_prof(p,j)
           end do

           ! wood harvest mortality carbon fluxes
           harvest_c_to_cwdc(c,j)  = harvest_c_to_cwdc(c,j)  + &
                hrv_livestemc_to_litter(p)  * wtcol(p) * stem_prof(p,j) 
           harvest_c_to_cwdc(c,j) = harvest_c_to_cwdc(c,j) + &
                hrv_livecrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)
           harvest_c_to_cwdc(c,j) = harvest_c_to_cwdc(c,j) + &
                hrv_deadcrootc_to_litter(p) * wtcol(p) * croot_prof(p,j) 

           ! storage harvest mortality carbon fluxes
           ! Metabolic litter is treated differently than other types
           ! of litter, so it gets this additional line after the
           ! most recent loop over all litter types
           harvest_c_to_litr_c(c,j,i_met_lit) = &
              harvest_c_to_litr_c(c,j,i_met_lit) + &
              hrv_leafc_storage_to_litter(p) * wtcol(p) * leaf_prof(p,j) + &
              hrv_frootc_storage_to_litter(p) * wtcol(p) * froot_prof(p,j) + &
              hrv_livestemc_storage_to_litter(p) * wtcol(p) * stem_prof(p,j) + &
              hrv_deadstemc_storage_to_litter(p) * wtcol(p) * stem_prof(p,j) + &
              hrv_livecrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
              hrv_deadcrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
              hrv_gresp_storage_to_litter(p) * wtcol(p) * leaf_prof(p,j) + &

           ! transfer harvest mortality carbon fluxes
              hrv_leafc_xfer_to_litter(p) * wtcol(p) * leaf_prof(p,j) + &
              hrv_frootc_xfer_to_litter(p) * wtcol(p) * froot_prof(p,j) + &
              hrv_livestemc_xfer_to_litter(p) * wtcol(p) * stem_prof(p,j) + &
              hrv_deadstemc_xfer_to_litter(p) * wtcol(p) * stem_prof(p,j) + &
              hrv_livecrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
              hrv_deadcrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
              hrv_gresp_xfer_to_litter(p) * wtcol(p) * leaf_prof(p,j)

           do i = i_litr_min, i_litr_max
              harvest_n_to_litr_n(c,j,i) = &
                 harvest_n_to_litr_n(c,j,i) + &
                 ! leaf harvest mortality nitrogen fluxes
                 hrv_leafn_to_litter(p) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j) + &
                 ! fine root litter nitrogen fluxes
                 hrv_frootn_to_litter(p) * fr_f(ivt(p),i) * wtcol(p) * froot_prof(p,j)
           end do

           ! wood harvest mortality nitrogen fluxes
           harvest_n_to_cwdn(c,j)  = harvest_n_to_cwdn(c,j)  + &
                hrv_livestemn_to_litter(p)  * wtcol(p) * stem_prof(p,j)
           harvest_n_to_cwdn(c,j) = harvest_n_to_cwdn(c,j) + &
                hrv_livecrootn_to_litter(p) * wtcol(p) * croot_prof(p,j)
           harvest_n_to_cwdn(c,j) = harvest_n_to_cwdn(c,j) + &
                hrv_deadcrootn_to_litter(p) * wtcol(p) * croot_prof(p,j)

           ! Metabolic litter is treated differently than other types
           ! of litter, so it gets this additional line after the
           ! most recent loop over all litter types
           harvest_n_to_litr_n(c,j,i_met_lit) = &
              harvest_n_to_litr_n(c,j,i_met_lit) + &
              ! retranslocated N pool harvest mortality fluxes
              hrv_retransn_to_litter(p) * wtcol(p) * leaf_prof(p,j) + &
              ! storage harvest mortality nitrogen fluxes
              hrv_leafn_storage_to_litter(p) * wtcol(p) * leaf_prof(p,j) + &
              hrv_frootn_storage_to_litter(p) * wtcol(p) * froot_prof(p,j) + &
              hrv_livestemn_storage_to_litter(p) * wtcol(p) * stem_prof(p,j) + &
              hrv_deadstemn_storage_to_litter(p) * wtcol(p) * stem_prof(p,j) + &
              hrv_livecrootn_storage_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
              hrv_deadcrootn_storage_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
              ! transfer harvest mortality nitrogen fluxes
              hrv_leafn_xfer_to_litter(p) * wtcol(p) * leaf_prof(p,j) + &
              hrv_frootn_xfer_to_litter(p) * wtcol(p) * froot_prof(p,j) + &
              hrv_livestemn_xfer_to_litter(p) * wtcol(p) * stem_prof(p,j) + &
              hrv_deadstemn_xfer_to_litter(p) * wtcol(p) * stem_prof(p,j) + &
              hrv_livecrootn_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
              hrv_deadcrootn_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)

        end do
     end do
   
     do fp = 1,num_soilp
        p = filter_soilp(fp)
        c = patch%column(p)

        ! wood harvest mortality carbon fluxes to product pools
        cwood_harvestc(c)  = cwood_harvestc(c)  + &
             pwood_harvestc(p)  * wtcol(p)

        ! wood harvest mortality nitrogen fluxes to product pools
        cwood_harvestn(c)  = cwood_harvestn(c)  + &
             pwood_harvestn(p)  * wtcol(p)

     end do

   end associate 

 end subroutine CNHarvestPftToColumn

end module dynHarvestMod
