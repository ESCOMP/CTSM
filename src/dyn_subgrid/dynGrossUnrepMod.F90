module dynGrossUnrepMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle reading of the gross unrepresented landcover change data, as well as the 
  ! state updates that happen as a result.
  !
  ! !USES:
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use decompMod               , only : bounds_type, BOUNDS_LEVEL_PROC
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
  use clm_varpar              , only : natpft_size, i_litr_min, i_litr_max, i_met_lit
  use ColumnType              , only : col                
  use PatchType               , only : patch                
  use CNSharedParamsMod       , only : use_matrixcn
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  !
  public :: dynGrossUnrep_init    ! initialize data structures for grossunrep information
  public :: dynGrossUnrep_interp  ! get grossunrep data for current time step, if needed
  public :: CNGrossUnrep          ! grossunrep mortality routine for CN code
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: CNGrossUnrepPftToColumn   ! gather patch-level grossunrep fluxes to the column level
  !
  ! !PRIVATE TYPES:

  ! Note that, since we have our own dyngrossunrep_file object (distinct from dynpft_file),
  ! we could theoretically have a separate file providing grossunrep data from that providing
  ! the pftdyn data
  type(dyn_file_type), target :: dyngrossunrep_file ! information for the file containing grossunrep data

  ! Define the underlying grossunrep variables
  character(len=64), parameter :: grossunrep_varname = 'UNREPRESENTED_PFT_LULCC'
  
  type(dyn_var_time_uninterp_type) :: gru_inst   ! value of gross unrepresented variable

  real(r8) , allocatable   :: grossunrepfrac(:,:) ! gross unrepresented landcover change fraction
  logical                  :: do_grossunrep ! whether we're in a period when we should do gross unrepresented landcover change
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dynGrossUnrep_init(bounds, grossunrep_filename)
    !
    ! !DESCRIPTION:
    ! Initialize data structures for grossunrep information.
    ! This should be called once, during model initialization.
    ! 
    ! !USES:
    use dynVarTimeUninterpMod , only : dyn_var_time_uninterp_type
    use dynTimeInfoMod        , only : YEAR_POSITION_START_OF_TIMESTEP
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds              ! proc-level bounds
    character(len=*)  , intent(in) :: grossunrep_filename ! name of file containing grossunrep information
    !
    ! !LOCAL VARIABLES:
    integer :: num_points ! number of spatial points
    integer :: ier        ! error code
    
    character(len=*), parameter :: subname = 'dynGrossUnrep_init'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    allocate(grossunrepfrac(bounds%begg:bounds%endg,0:natpft_size-1),stat=ier)
    if (ier /= 0) then
       call endrun(msg=' allocation error for grossunrep'//errMsg(sourcefile, __LINE__))
    end if

    ! Get the year from the START of the timestep for consistency with other dyn file
    ! stuff (though it shouldn't actually matter for grossunrep, given the way the
    ! once-per-year grossunrep is applied).
    dyngrossunrep_file = dyn_file_type(grossunrep_filename, YEAR_POSITION_START_OF_TIMESTEP)
    
    ! Get initial grossunrep data
    num_points = (bounds%endg - bounds%begg + 1)
    gru_inst = dyn_var_time_uninterp_type( &
            dyn_file=dyngrossunrep_file, varname=grossunrep_varname, &
            dim1name=grlnd, conversion_factor=1.0_r8, &
            do_check_sums_equal_1=.false., data_shape=[num_points,natpft_size], &
	    allow_nodata=.true.)
       
  end subroutine dynGrossUnrep_init


  !-----------------------------------------------------------------------
  subroutine dynGrossUnrep_interp(bounds)
    !
    ! !DESCRIPTION:
    ! Get grossunrep data for model time, when needed.
    !
    ! Note that grossunrep data are stored as rates (not weights) and so time interpolation
    ! is not necessary - the grossunrep rate is held constant through the year.  This is
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

    character(len=*), parameter :: subname = 'dyngru_interp'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    call dyngrossunrep_file%time_info%set_current_year()

    ! Get total grossunrep for this time step
    grossunrepfrac(bounds%begg:bounds%endg,0:natpft_size-1) = 0._r8

    if (dyngrossunrep_file%time_info%is_before_time_series()) then
       ! Turn off grossunrep before the start of the grossunrep time series
       do_grossunrep = .false.
    else
       ! Note that do_grossunrep stays true even past the end of the time series. This
       ! means that grossunrep rates will be maintained at the rate given in the last
       ! year of the file for all years past the end of this specified time series.
       do_grossunrep = .true.
       call gru_inst%get_current_data(grossunrepfrac)
    end if

  end subroutine dynGrossUnrep_interp


  !-----------------------------------------------------------------------
  subroutine CNGrossUnrep (num_soilp, filter_soilp, &
       soilbiogeochem_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, &
       cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! GrossUnrep mortality routine for coupled carbon-nitrogen code (CN)
    !
    ! !USES:
    use pftconMod       , only : noveg, nbrdlf_evr_shrub, nc4_grass
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
    real(r8):: thistreec                 ! carbon in this tree for calculating grossunrep fraction (gC/m2)
    real(r8):: cm                        ! rate for carbon grossunrep mortality (gC/m2/yr)
    real(r8):: am                        ! rate for fractional grossunrep mortality (1/yr)
    real(r8):: m                         ! rate for fractional grossunrep mortality (1/s)
    real(r8):: dtime                     ! model time step (s)
    !-----------------------------------------------------------------------

    associate(& 
         ivt                                 =>    patch%itype                                                      , & ! Input:  [integer (:)]  pft vegetation type                                
         convfrac                            =>    pftcon%pconv                                                     , & ! Input:  [real (:)] pft conversion to atmosphere fraction
         
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
         
         gru_leafc_to_litter                 =>    cnveg_carbonflux_inst%gru_leafc_to_litter_patch                , & ! Output: [real(r8) (:)]                                                    
         gru_frootc_to_litter                =>    cnveg_carbonflux_inst%gru_frootc_to_litter_patch               , & ! Output: [real(r8) (:)]                                                    
         gru_livestemc_to_atm                =>    cnveg_carbonflux_inst%gru_livestemc_to_atm_patch               , & ! Output: [real(r8) (:)]                                                    
         gru_deadstemc_to_atm                =>    cnveg_carbonflux_inst%gru_deadstemc_to_atm_patch               , & ! Output: [real(r8) (:)]
         gru_wood_productc_gain              =>    cnveg_carbonflux_inst%gru_wood_productc_gain_patch             , & ! Output: [real(r8) (:)]
         gru_livecrootc_to_litter            =>    cnveg_carbonflux_inst%gru_livecrootc_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         gru_deadcrootc_to_litter            =>    cnveg_carbonflux_inst%gru_deadcrootc_to_litter_patch             , & ! Output: [real(r8) (:)]                                                    
         gru_xsmrpool_to_atm                 =>    cnveg_carbonflux_inst%gru_xsmrpool_to_atm_patch                , & ! Output: [real(r8) (:)]                                                    
         gru_leafc_storage_to_atm            =>    cnveg_carbonflux_inst%gru_leafc_storage_to_atm_patch           , & ! Output: [real(r8) (:)]                                                    
         gru_frootc_storage_to_atm           =>    cnveg_carbonflux_inst%gru_frootc_storage_to_atm_patch          , & ! Output: [real(r8) (:)]                                                    
         gru_livestemc_storage_to_atm        =>    cnveg_carbonflux_inst%gru_livestemc_storage_to_atm_patch       , & ! Output: [real(r8) (:)]                                                    
         gru_deadstemc_storage_to_atm        =>    cnveg_carbonflux_inst%gru_deadstemc_storage_to_atm_patch       , & ! Output: [real(r8) (:)]                                                    
         gru_livecrootc_storage_to_atm       =>    cnveg_carbonflux_inst%gru_livecrootc_storage_to_atm_patch      , & ! Output: [real(r8) (:)]                                                    
         gru_deadcrootc_storage_to_atm       =>    cnveg_carbonflux_inst%gru_deadcrootc_storage_to_atm_patch      , & ! Output: [real(r8) (:)]                                                    
         gru_gresp_storage_to_atm            =>    cnveg_carbonflux_inst%gru_gresp_storage_to_atm_patch           , & ! Output: [real(r8) (:)]                                                    
         gru_leafc_xfer_to_atm               =>    cnveg_carbonflux_inst%gru_leafc_xfer_to_atm_patch              , & ! Output: [real(r8) (:)]                                                    
         gru_frootc_xfer_to_atm              =>    cnveg_carbonflux_inst%gru_frootc_xfer_to_atm_patch             , & ! Output: [real(r8) (:)]                                                    
         gru_livestemc_xfer_to_atm           =>    cnveg_carbonflux_inst%gru_livestemc_xfer_to_atm_patch          , & ! Output: [real(r8) (:)]                                                    
         gru_deadstemc_xfer_to_atm           =>    cnveg_carbonflux_inst%gru_deadstemc_xfer_to_atm_patch          , & ! Output: [real(r8) (:)]                                                    
         gru_livecrootc_xfer_to_atm          =>    cnveg_carbonflux_inst%gru_livecrootc_xfer_to_atm_patch         , & ! Output: [real(r8) (:)]                                                    
         gru_deadcrootc_xfer_to_atm          =>    cnveg_carbonflux_inst%gru_deadcrootc_xfer_to_atm_patch         , & ! Output: [real(r8) (:)]                                                    
         gru_gresp_xfer_to_atm               =>    cnveg_carbonflux_inst%gru_gresp_xfer_to_atm_patch              , & ! Output: [real(r8) (:)]                                                    
         
         gru_leafn_to_litter                 =>    cnveg_nitrogenflux_inst%gru_leafn_to_litter_patch              , & ! Output: [real(r8) (:)]                                                    
         gru_frootn_to_litter                =>    cnveg_nitrogenflux_inst%gru_frootn_to_litter_patch             , & ! Output: [real(r8) (:)]                                                    
         gru_livestemn_to_atm                =>    cnveg_nitrogenflux_inst%gru_livestemn_to_atm_patch             , & ! Output: [real(r8) (:)]                                                    
         gru_deadstemn_to_atm                =>    cnveg_nitrogenflux_inst%gru_deadstemn_to_atm_patch             , & ! Output: [real(r8) (:)]
         gru_wood_productn_gain              =>    cnveg_nitrogenflux_inst%gru_wood_productn_gain_patch           , & ! Output: [real(r8) (:)]
         gru_livecrootn_to_litter            =>    cnveg_nitrogenflux_inst%gru_livecrootn_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
         gru_deadcrootn_to_litter            =>    cnveg_nitrogenflux_inst%gru_deadcrootn_to_litter_patch         , & ! Output: [real(r8) (:)]                                                    
         gru_retransn_to_litter              =>    cnveg_nitrogenflux_inst%gru_retransn_to_litter_patch           , & ! Output: [real(r8) (:)]                                                    
         gru_leafn_storage_to_atm            =>    cnveg_nitrogenflux_inst%gru_leafn_storage_to_atm_patch         , & ! Output: [real(r8) (:)]                                                    
         gru_frootn_storage_to_atm           =>    cnveg_nitrogenflux_inst%gru_frootn_storage_to_atm_patch        , & ! Output: [real(r8) (:)]                                                    
         gru_livestemn_storage_to_atm        =>    cnveg_nitrogenflux_inst%gru_livestemn_storage_to_atm_patch     , & ! Output: [real(r8) (:)]                                                    
         gru_deadstemn_storage_to_atm        =>    cnveg_nitrogenflux_inst%gru_deadstemn_storage_to_atm_patch     , & ! Output: [real(r8) (:)]                                                    
         gru_livecrootn_storage_to_atm       =>    cnveg_nitrogenflux_inst%gru_livecrootn_storage_to_atm_patch    , & ! Output: [real(r8) (:)]                                                    
         gru_deadcrootn_storage_to_atm       =>    cnveg_nitrogenflux_inst%gru_deadcrootn_storage_to_atm_patch    , & ! Output: [real(r8) (:)]                                                    
         gru_leafn_xfer_to_atm               =>    cnveg_nitrogenflux_inst%gru_leafn_xfer_to_atm_patch            , & ! Output: [real(r8) (:)]                                                    
         gru_frootn_xfer_to_atm              =>    cnveg_nitrogenflux_inst%gru_frootn_xfer_to_atm_patch           , & ! Output: [real(r8) (:)]                                                    
         gru_livestemn_xfer_to_atm           =>    cnveg_nitrogenflux_inst%gru_livestemn_xfer_to_atm_patch        , & ! Output: [real(r8) (:)]                                                    
         gru_deadstemn_xfer_to_atm           =>    cnveg_nitrogenflux_inst%gru_deadstemn_xfer_to_atm_patch        , & ! Output: [real(r8) (:)]                                                    
         gru_livecrootn_xfer_to_atm          =>    cnveg_nitrogenflux_inst%gru_livecrootn_xfer_to_atm_patch       , & ! Output: [real(r8) (:)]                                                    
         gru_deadcrootn_xfer_to_atm          =>    cnveg_nitrogenflux_inst%gru_deadcrootn_xfer_to_atm_patch       , & ! Output: [real(r8) (:)]
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

         ! If this is a natural pft, then
         ! get the annual grossunrep "mortality" rate (am) from grossunrep array
         ! and convert to rate per second
         if (ivt(p) > noveg .and. ivt(p) <= nc4_grass) then

            if (do_grossunrep) then
               am = grossunrepfrac(g,ivt(p))

               ! Apply all grossunrep at the start of the year
               if (is_beg_curr_year()) then
                  m  = am/dtime
               else
                  m = 0._r8
               end if
            else
               m = 0._r8
            end if

            if(.not. use_matrixcn)then
               ! patch-level gross unrepresented landcover change carbon fluxes
               ! displayed pools
               gru_leafc_to_litter(p)               = leafc(p)               * m
               gru_frootc_to_litter(p)              = frootc(p)              * m
               gru_livestemc_to_atm(p)              = livestemc(p)           * m
               gru_deadstemc_to_atm(p)              = deadstemc(p)           * m * convfrac(ivt(p))
               gru_wood_productc_gain(p)            = deadstemc(p)           * m * (1._r8 - convfrac(ivt(p)))
               gru_livecrootc_to_litter(p)          = livecrootc(p)          * m
               gru_deadcrootc_to_litter(p)          = deadcrootc(p)          * m
               gru_xsmrpool_to_atm(p)               = xsmrpool(p)            * m

               ! storage pools
               gru_leafc_storage_to_atm(p)          = leafc_storage(p)       * m
               gru_frootc_storage_to_atm(p)         = frootc_storage(p)      * m
               gru_livestemc_storage_to_atm(p)      = livestemc_storage(p)   * m
               gru_deadstemc_storage_to_atm(p)      = deadstemc_storage(p)   * m
               gru_livecrootc_storage_to_atm(p)     = livecrootc_storage(p)  * m
               gru_deadcrootc_storage_to_atm(p)     = deadcrootc_storage(p)  * m
               gru_gresp_storage_to_atm(p)          = gresp_storage(p)       * m

               ! transfer pools
               gru_leafc_xfer_to_atm(p)             = leafc_xfer(p)          * m
               gru_frootc_xfer_to_atm(p)            = frootc_xfer(p)         * m
               gru_livestemc_xfer_to_atm(p)         = livestemc_xfer(p)      * m
               gru_deadstemc_xfer_to_atm(p)         = deadstemc_xfer(p)      * m
               gru_livecrootc_xfer_to_atm(p)        = livecrootc_xfer(p)     * m
               gru_deadcrootc_xfer_to_atm(p)        = deadcrootc_xfer(p)     * m
               gru_gresp_xfer_to_atm(p)             = gresp_xfer(p)          * m

               ! patch-level gross unrepresented landcover change mortality nitrogen fluxes
               ! displayed pools
               gru_leafn_to_litter(p)               = leafn(p)               * m
               gru_frootn_to_litter(p)              = frootn(p)              * m
               gru_livestemn_to_atm(p)              = livestemn(p)           * m
               gru_deadstemn_to_atm(p)              = deadstemn(p)           * m * convfrac(ivt(p))
               gru_wood_productn_gain(p)            = deadstemn(p)           * m * (1._r8 - convfrac(ivt(p)))
               gru_livecrootn_to_litter(p)          = livecrootn(p)          * m
               gru_deadcrootn_to_litter(p)          = deadcrootn(p)          * m
               gru_retransn_to_litter(p)            = retransn(p)            * m

               ! storage pools
               gru_leafn_storage_to_atm(p)          = leafn_storage(p)       * m
               gru_frootn_storage_to_atm(p)         = frootn_storage(p)      * m
               gru_livestemn_storage_to_atm(p)      = livestemn_storage(p)   * m
               gru_deadstemn_storage_to_atm(p)      = deadstemn_storage(p)   * m
               gru_livecrootn_storage_to_atm(p)     = livecrootn_storage(p)  * m
               gru_deadcrootn_storage_to_atm(p)     = deadcrootn_storage(p)  * m

               ! transfer pools
               gru_leafn_xfer_to_atm(p)             = leafn_xfer(p)          * m
               gru_frootn_xfer_to_atm(p)            = frootn_xfer(p)         * m
               gru_livestemn_xfer_to_atm(p)         = livestemn_xfer(p)      * m
               gru_deadstemn_xfer_to_atm(p)         = deadstemn_xfer(p)      * m
               gru_livecrootn_xfer_to_atm(p)        = livecrootn_xfer(p)     * m
               gru_deadcrootn_xfer_to_atm(p)        = deadcrootn_xfer(p)     * m
            else  ! matrixcn solution
               ! patch-level gross unrepresented landcover change carbon fluxes
               ! displayed pools
               gru_leafc_to_litter(p) = matrix_update_gmc(p,ileaf_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             leafc(p)
               gru_frootc_to_litter(p) = matrix_update_gmc(p,ifroot_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             frootc(p)
               gru_livestemc_to_atm(p) = matrix_update_gmc(p,ilivestem_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             livestemc(p)
               gru_deadstemc_to_atm(p) = matrix_update_gmc(p,ideadstem_to_iout_gmc,m * convfrac(ivt(p)),dtime,cnveg_carbonflux_inst,.True.,.True.)  * &
                                             deadstemc(p)
               gru_wood_productc_gain(p) = matrix_update_gmc(p,ideadstem_to_iout_gmc,m * (1._r8 - convfrac(ivt(p))),dtime,cnveg_carbonflux_inst,.True.,.True.)  * &
                                             deadstemc(p)
               gru_livecrootc_to_litter(p) = matrix_update_gmc(p,ilivecroot_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             livecrootc(p)
               gru_deadcrootc_to_litter(p) = matrix_update_gmc(p,ideadcroot_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             deadcrootc(p)
               gru_xsmrpool_to_atm(p) = xsmrpool(p) * m

               ! storage pools
               gru_leafc_storage_to_atm(p) = matrix_update_gmc(p,ileafst_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             leafc_storage(p)
               gru_frootc_storage_to_atm(p) = matrix_update_gmc(p,ifrootst_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             frootc_storage(p)
               gru_livestemc_storage_to_atm(p) = matrix_update_gmc(p,ilivestemst_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             livestemc_storage(p)
               gru_deadstemc_storage_to_atm(p) = matrix_update_gmc(p,ideadstemst_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             deadstemc_storage(p)
               gru_livecrootc_storage_to_atm(p) = matrix_update_gmc(p,ilivecrootst_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             livecrootc_storage(p)
               gru_deadcrootc_storage_to_atm(p) = matrix_update_gmc(p,ideadcrootst_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             deadcrootc_storage(p)
               gru_gresp_storage_to_atm(p) = gresp_storage(p) * m

               ! transfer pools
               gru_leafc_xfer_to_atm(p) = matrix_update_gmc(p,ileafxf_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             leafc_xfer(p)
               gru_frootc_xfer_to_atm(p) = matrix_update_gmc(p,ifrootxf_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             frootc_xfer(p)
               gru_livestemc_xfer_to_atm(p) = matrix_update_gmc(p,ilivestemxf_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             livestemc_xfer(p)
               gru_deadstemc_xfer_to_atm(p) = matrix_update_gmc(p,ideadstemxf_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             deadstemc_xfer(p)
               gru_livecrootc_xfer_to_atm(p) = matrix_update_gmc(p,ilivecrootxf_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             livecrootc_xfer(p)
               gru_deadcrootc_xfer_to_atm(p) = matrix_update_gmc(p,ideadcrootxf_to_iout_gmc,m,dtime,cnveg_carbonflux_inst,.True.,.True.) * &
                                             deadcrootc_xfer(p)
               gru_gresp_xfer_to_atm(p) = gresp_xfer(p) * m

               ! patch-level gross unrepresented landcover change mortality nitrogen fluxes
               ! displayed pools
               gru_leafn_to_litter(p) = matrix_update_gmn(p,ileaf_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             leafn(p)
               gru_frootn_to_litter(p) = matrix_update_gmn(p,ifroot_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             frootn(p)
               gru_livestemn_to_atm(p) = matrix_update_gmn(p,ilivestem_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             livestemn(p)
               gru_deadstemn_to_atm(p) = matrix_update_gmn(p,ideadstem_to_iout_gmn,m * convfrac(ivt(p)),dtime,cnveg_nitrogenflux_inst,.True.,.True.)  * &
                                             deadstemn(p)
               gru_wood_productn_gain(p) = matrix_update_gmn(p,ideadstem_to_iout_gmn,m * (1._r8 - convfrac(ivt(p))),dtime,cnveg_nitrogenflux_inst,.True.,.True.)  * &
                                             deadstemn(p)
               gru_livecrootn_to_litter(p) = matrix_update_gmn(p,ilivecroot_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             livecrootn(p)
               gru_deadcrootn_to_litter(p) = matrix_update_gmn(p,ideadcroot_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             deadcrootn(p)
               gru_retransn_to_litter(p)   = matrix_update_gmn(p,iretransn_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             retransn(p)

               ! storage pools
               gru_leafn_storage_to_atm(p) = matrix_update_gmn(p,ileafst_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             leafn_storage(p)
               gru_frootn_storage_to_atm(p) = matrix_update_gmn(p,ifrootst_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             frootn_storage(p)
               gru_livestemn_storage_to_atm(p) = matrix_update_gmn(p,ilivestemst_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             livestemn_storage(p)
               gru_deadstemn_storage_to_atm(p) = matrix_update_gmn(p,ideadstemst_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             deadstemn_storage(p)
               gru_livecrootn_storage_to_atm(p) = matrix_update_gmn(p,ilivecrootst_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.)* &
                                             livecrootn_storage(p)
               gru_deadcrootn_storage_to_atm(p) = matrix_update_gmn(p,ideadcrootst_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.)* &
                                             deadcrootn_storage(p)
               ! transfer pools
               gru_leafn_xfer_to_atm(p) = matrix_update_gmn(p,ileafxf_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             leafn_xfer(p)
               gru_frootn_xfer_to_atm(p) = matrix_update_gmn(p,ifrootxf_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             frootn_xfer(p)
               gru_livestemn_xfer_to_atm(p) = matrix_update_gmn(p,ilivestemxf_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             livestemn_xfer(p)
               gru_deadstemn_xfer_to_atm(p) = matrix_update_gmn(p,ideadstemxf_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             deadstemn_xfer(p)
               gru_livecrootn_xfer_to_atm(p) = matrix_update_gmn(p,ilivecrootxf_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             livecrootn_xfer(p)
               gru_deadcrootn_xfer_to_atm(p) = matrix_update_gmn(p,ideadcrootxf_to_iout_gmn,m,dtime,cnveg_nitrogenflux_inst,.True.,.True.) * &
                                             deadcrootn_xfer(p)
            end if

         end if  ! end tree block

      end do ! end of pft loop

      ! gather all patch-level litterfall fluxes from grossunrep to the column
      ! for litter C and N inputs

      call CNGrossUnrepPftToColumn(num_soilp, filter_soilp, &
           soilbiogeochem_state_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)

    end associate 

  end subroutine CNGrossUnrep

 !-----------------------------------------------------------------------
 subroutine CNGrossUnrepPftToColumn (num_soilp, filter_soilp, &
      soilbiogeochem_state_inst, CNVeg_carbonflux_inst, cnveg_nitrogenflux_inst)
   !
   ! !DESCRIPTION:
   ! called at the end of CNGrossUnrep to gather all patch-level grossunrep litterfall fluxes
   ! to the column level and assign them to the three litter pools
   !
   ! !USES:
   use clm_varpar , only : maxsoil_patches, nlevdecomp
   !
   ! !ARGUMENTS:
   integer                         , intent(in)    :: num_soilp       ! number of soil patches in filter
   integer                         , intent(in)    :: filter_soilp(:) ! soil patch filter
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
        
        lf_f                             =>    pftcon%lf_f                                                    , & ! Input:  leaf litter fractions
        fr_f                             =>    pftcon%fr_f                                                    , & ! Input:  fine root litter fractions
        
        leaf_prof                        =>    soilbiogeochem_state_inst%leaf_prof_patch                      , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
        froot_prof                       =>    soilbiogeochem_state_inst%froot_prof_patch                     , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
        croot_prof                       =>    soilbiogeochem_state_inst%croot_prof_patch                     , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
        stem_prof                        =>    soilbiogeochem_state_inst%stem_prof_patch                      , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          
        
        gru_leafc_to_litter                 =>    cnveg_carbonflux_inst%gru_leafc_to_litter_patch                , & ! Input: [real(r8) (:)]                                                    
        gru_frootc_to_litter                =>    cnveg_carbonflux_inst%gru_frootc_to_litter_patch               , & ! Input: [real(r8) (:)]                                                    
        gru_wood_productc_gain              =>    cnveg_carbonflux_inst%gru_wood_productc_gain_patch             , & ! Input: [real(r8) (:)]
        gru_livecrootc_to_litter            =>    cnveg_carbonflux_inst%gru_livecrootc_to_litter_patch           , & ! Input: [real(r8) (:)]                                                    
        gru_deadcrootc_to_litter            =>    cnveg_carbonflux_inst%gru_deadcrootc_to_litter_patch           , & ! Input: [real(r8) (:)]                                                    
        gru_conv_cflux                      =>    cnveg_carbonflux_inst%gru_conv_cflux_patch                     , & ! Input: [real(r8) (:)   ]  C fluxes associated with conversion (immediate loss to atm)

        gru_c_to_litr_c                     =>    cnveg_carbonflux_inst%gru_c_to_litr_c_col                      , & ! Output:  [real(r8) (:,:,:) ]  C fluxes associated with gross unrepresented landcover change to litter pools (gC/m3/s)
        gru_c_to_cwdc_c                     =>    cnveg_carbonflux_inst%gru_c_to_cwdc_col                        , & ! Output:  [real(r8) (:,:) ]  C fluxes associated with grossunrep to CWD pool (gC/m3/s)

        gru_wood_productc_gain_c            =>    cnveg_carbonflux_inst%gru_wood_productc_gain_col               , & ! Input: [real(r8) (:)]
        
        gru_leafn_to_litter                 =>    cnveg_nitrogenflux_inst%gru_leafn_to_litter_patch              , & ! Input: [real(r8) (:)]                                                    
        gru_frootn_to_litter                =>    cnveg_nitrogenflux_inst%gru_frootn_to_litter_patch             , & ! Input: [real(r8) (:)]                                                    
        gru_wood_productn_gain              =>    cnveg_nitrogenflux_inst%gru_wood_productn_gain_patch           , & ! Input: [real(r8) (:)]
        gru_livecrootn_to_litter            =>    cnveg_nitrogenflux_inst%gru_livecrootn_to_litter_patch         , & ! Input: [real(r8) (:)]                                                    
        gru_deadcrootn_to_litter            =>    cnveg_nitrogenflux_inst%gru_deadcrootn_to_litter_patch         , & ! Input: [real(r8) (:)]                                                    
        gru_retransn_to_litter              =>    cnveg_nitrogenflux_inst%gru_retransn_to_litter_patch           , & ! Input: [real(r8) (:)]                                                    
        gru_n_to_litr_c                     =>    cnveg_nitrogenflux_inst%gru_n_to_litr_n_col                    , & ! Output:  [real(r8) (:,:,:) ]  N fluxes associated with grossunrep to litter pools (gN/m3/s)
        gru_n_to_cwdn_c                     =>    cnveg_nitrogenflux_inst%gru_n_to_cwdn_col                      , & ! Output:  [real(r8) (:,:) ]  N fluxes associated with grossunrep to CWD pool (gN/m3/s)

        gru_wood_productn_gain_c            =>    cnveg_nitrogenflux_inst%gru_wood_productn_gain_col               & ! Input: [real(r8) (:)]
        )

     do j = 1, nlevdecomp
        do fp = 1,num_soilp
           p = filter_soilp(fp)
           c = patch%column(p)

           do i = i_litr_min, i_litr_max
              gru_c_to_litr_c(c,j,i) = gru_c_to_litr_c(c,j,i) + &
                   ! leaf gross unrepresented landcover change mortality carbon fluxes
                   gru_leafc_to_litter(p) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j) + &
                   ! fine root gross unrepresented landcover change mortality carbon fluxes
                   gru_frootc_to_litter(p) * fr_f(ivt(p),i) * wtcol(p) * froot_prof(p,j)
              gru_n_to_litr_c(c,j,i) = gru_n_to_litr_c(c,j,i) + &
                   ! leaf gross unrepresented landcover change mortality nitrogen fluxes
                   gru_leafn_to_litter(p) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j) + &
                   ! fine root gross unrepresented landcover change mortality nitrogen fluxes
                   gru_frootn_to_litter(p) * fr_f(ivt(p),i) * wtcol(p) * froot_prof(p,j)
           end do

           ! coarse root gross unrepresented landcover change mortality carbon fluxes
           gru_c_to_cwdc_c(c,j) = gru_c_to_cwdc_c(c,j) + &
                gru_livecrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)
           gru_c_to_cwdc_c(c,j) = gru_c_to_cwdc_c(c,j) + &
                gru_deadcrootc_to_litter(p) * wtcol(p) * croot_prof(p,j) 

           ! coarse root gross unrepresented landcover change mortality nitrogen fluxes
           gru_n_to_cwdn_c(c,j) = gru_n_to_cwdn_c(c,j) + &
                gru_livecrootn_to_litter(p) * wtcol(p) * croot_prof(p,j)
           gru_n_to_cwdn_c(c,j) = gru_n_to_cwdn_c(c,j) + &
                gru_deadcrootn_to_litter(p) * wtcol(p) * croot_prof(p,j)

           ! retranslocated N pool gross unrepresented landcover change mortality fluxes
           ! process specific to i_met_lit, so we keep it outside
           ! the i_litr_min to i_litr_max loop above
           gru_n_to_litr_c(c,j,i_met_lit) =  &
                gru_n_to_litr_c(c,j,i_met_lit) + &
                gru_retransn_to_litter(p) * wtcol(p) * leaf_prof(p,j)

        end do
     end do
   
     do fp = 1,num_soilp
        p = filter_soilp(fp)
        c = patch%column(p)

        ! wood gross unrepresented landcover change mortality carbon fluxes to product pools
        gru_wood_productc_gain_c(c)  = gru_wood_productc_gain_c(c)  + &
             gru_wood_productc_gain(p)  * wtcol(p)

        ! wood gross unrepresented landcover change mortality nitrogen fluxes to product pools
        gru_wood_productn_gain_c(c)  = gru_wood_productn_gain_c(c)  + &
             gru_wood_productn_gain(p)  * wtcol(p)

     end do

   end associate 

 end subroutine CNGrossUnrepPftToColumn

end module dynGrossUnrepMod
