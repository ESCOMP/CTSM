module EDClmtypeInitMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Allocate EDClmtype components and initialize them to signaling NaN.
  !
  use shr_kind_mod, only: &
       r8 => shr_kind_r8
  use shr_infnan_mod, only: &
       nan => shr_infnan_nan, assignment(=)

  implicit none
  save
  private

  public :: EDInitClmType
  public :: EDinit_cohort_type

  private :: EDinit_pft_type
  private :: EDinit_pft_ecophys_constants
  private :: EDinit_pft_cflux_type

contains

  !------------------------------------------------------------------------
  subroutine EDInitClmType( bounds, EDpft, EDpftcon, EDpcf )
    !
    use decompMod , only : bounds_type
    use EDclmType , only : EDpft_type, EDpft_epc_type, EDpft_cflux_type

    implicit none

    type(bounds_type),       intent(in)    :: bounds  ! bounds
    type(EDpft_type),        intent(inout) :: EDpft
    type(EDpft_epc_type),    intent(inout) :: EDpftcon
    type (EDpft_cflux_type), intent(inout) :: EDpcf

    !------------------------------------------------------------------------

    call EDinit_pft_type( bounds%begp, bounds%endp, EDpft )
    call EDinit_pft_ecophys_constants( EDpftcon )
    call EDinit_pft_cflux_type( bounds%begp, bounds%endp, EDpcf )

  end subroutine EDInitClmType

  !------------------------------------------------------------------------
  subroutine EDinit_pft_type (beg, end, EDpft)
    !
    ! !DESCRIPTION:
    ! Initialize components of pft_type structure
    !
    use EDclmType , only : EDpft_type

    implicit none

    integer, intent(in) :: beg, end

    type(EDpft_type), intent(inout) :: EDpft
    !------------------------------------------------------------------------

    allocate(EDpft%ED_patch(beg:end)) 
    allocate(EDpft%ED_bareground(beg:end)) 
    allocate(EDpft%wtED(beg:end)) 

  end subroutine EDinit_pft_type

  !------------------------------------------------------------------------
  subroutine EDinit_cohort_type (beg, end, coh)
    !
    ! !DESCRIPTION:
    ! Initialize components of landunit_type structure
    !
    use EDtypesMod , only : cohort_type

    implicit none

    integer, intent(in) :: beg, end
    type(cohort_type), intent(inout):: coh
    !------------------------------------------------------------------------

    ! FIX(SPM,032414) pull this out and put in own ED source
    allocate(coh%gridcell(beg:end))

  end subroutine EDinit_cohort_type

  !------------------------------------------------------------------------
  subroutine EDinit_pft_ecophys_constants( EDpftcon )
    !
    use clm_varpar,  only : numpft
    use EDclmType,   only : EDpft_epc_type

    implicit none

    type(EDpft_epc_type), intent(inout) :: EDpftcon

    allocate(EDpftcon%max_dbh(0:numpft))
    allocate(EDpftcon%freezetol(0:numpft))
    allocate(EDpftcon%wood_density(0:numpft))
    allocate(EDpftcon%alpha_stem(0:numpft))
    allocate(EDpftcon%hgt_min(0:numpft))
    allocate(EDpftcon%cushion(0:numpft))
    allocate(EDpftcon%leaf_stor_priority(0:numpft))
    allocate(EDpftcon%leafwatermax(0:numpft))
    allocate(EDpftcon%rootresist(0:numpft))
    allocate(EDpftcon%soilbeta(0:numpft))
    allocate(EDpftcon%crown(0:numpft))
    allocate(EDpftcon%bark_scaler(0:numpft))
    allocate(EDpftcon%crown_kill(0:numpft))
    allocate(EDpftcon%initd(0:numpft))
    allocate(EDpftcon%sd_mort(0:numpft))
    allocate(EDpftcon%seed_rain(0:numpft))
    allocate(EDpftcon%BB_slope(0:numpft))
    allocate(EDpftcon%root_long(0:numpft))
    allocate(EDpftcon%seed_alloc(0:numpft))
    allocate(EDpftcon%clone_alloc(0:numpft))
    allocate(EDpftcon%sapwood_ratio(0:numpft))

    EDpftcon%max_dbh(:) = nan
    EDpftcon%freezetol(:) = nan
    EDpftcon%wood_density(:) = nan
    EDpftcon%alpha_stem(:) = nan
    EDpftcon%hgt_min(:) = nan
    EDpftcon%cushion(:) = nan
    EDpftcon%leaf_stor_priority(:) = nan
    EDpftcon%leafwatermax(:) = nan
    EDpftcon%rootresist(:) = nan
    EDpftcon%soilbeta(:) = nan
    EDpftcon%crown(:) = nan
    EDpftcon%bark_scaler(:) = nan
    EDpftcon%crown_kill(:) = nan
    EDpftcon%initd(:) = nan
    EDpftcon%sd_mort(:) = nan
    EDpftcon%seed_rain(:) = nan
    EDpftcon%BB_slope(:) = nan
    EDpftcon%root_long(:)=nan
    EDpftcon%seed_alloc(:)=nan
    EDpftcon%clone_alloc(:)=nan
    EDpftcon%sapwood_ratio(:)=nan

  end subroutine EDinit_pft_ecophys_constants

  subroutine EDinit_pft_cflux_type(beg, end, EDpcf)
    !
    ! !DESCRIPTION:
    ! Initialize pft carbon flux variables
    !
    use clm_varpar,  only : nlevgrnd
    use EDclmType,   only : EDpft_cflux_type

    implicit none

    integer, intent(in) :: beg, end
    type (EDpft_cflux_type), intent(inout) :: EDpcf
    !------------------------------------------------------------------------

    allocate(EDpcf%trimming(beg:end))
    allocate(EDpcf%canopy_spread(beg:end))
    allocate(EDpcf%GCcanopy(beg:end))
    allocate(EDpcf%area_plant(beg:end))
    allocate(EDpcf%area_trees(beg:end))    
    allocate(EDpcf%PFTbiomass(beg:end,1:nlevgrnd))
    allocate(EDpcf%PFTleafbiomass(beg:end,1:nlevgrnd))
    allocate(EDpcf%PFTstorebiomass(beg:end,1:nlevgrnd))
    allocate(EDpcf%PFTnindivs(beg:end,1:nlevgrnd)) 
    allocate(EDpcf%nesterov_fire_danger(beg:end))
    allocate(EDpcf%spitfire_ROS(beg:end))
    allocate(EDpcf%effect_wspeed(beg:end))
    allocate(EDpcf%TFC_ROS(beg:end))
    allocate(EDpcf%fire_intensity(beg:end))
    allocate(EDpcf%fire_area(beg:end))
    allocate(EDpcf%scorch_height(beg:end))
    allocate(EDpcf%fire_fuel_bulkd(beg:end))
    allocate(EDpcf%fire_fuel_eff_moist(beg:end))
    allocate(EDpcf%fire_fuel_sav(beg:end))
    allocate(EDpcf%fire_fuel_mef(beg:end))
    allocate(EDpcf%sum_fuel(beg:end))
    allocate(EDpcf%litter_in(beg:end))
    allocate(EDpcf%litter_out(beg:end))                   
    allocate(EDpcf%efpot(beg:end)) 
    allocate(EDpcf%rb(beg:end))    
    allocate(EDpcf%prec24  (beg:end))
    allocate(EDpcf%rh24  (beg:end))
    allocate(EDpcf%wind24  (beg:end))
    allocate(EDpcf%seed_bank  (beg:end))
    allocate(EDpcf%seed_decay  (beg:end))
    allocate(EDpcf%seeds_in  (beg:end))
    allocate(EDpcf%seed_germination  (beg:end))

    EDpcf%trimming(beg:end)                   = 0.0_r8 
    EDpcf%canopy_spread(beg:end)              = 0.0_r8 
    EDpcf%GCcanopy(beg:end)                   = 0.0_r8 
    EDpcf%area_plant(beg:end)                 = 0.0_r8 
    EDpcf%area_trees(beg:end)                 = 0.0_r8         
    EDpcf%PFTbiomass(beg:end,1:nlevgrnd)      = 0.0_r8
    EDpcf%PFTleafbiomass(beg:end,1:nlevgrnd)  = 0.0_r8
    EDpcf%PFTstorebiomass(beg:end,1:nlevgrnd) = 0.0_r8
    EDpcf%PFTnindivs(beg:end,1:nlevgrnd)      = 0.0_r8            
    EDpcf%nesterov_fire_danger(beg:end)       = 0.0_r8     
    EDpcf%spitfire_ROS(beg:end)               = 0.0_r8 
    EDpcf%effect_wspeed(beg:end)              = 0.0_r8 
    EDpcf%TFC_ROS(beg:end)                    = 0.0_r8 
    EDpcf%fire_intensity(beg:end)             = 0.0_r8 
    EDpcf%fire_area(beg:end)                  = 0.0_r8 
    EDpcf%scorch_height(beg:end)              = 0.0_r8 
    EDpcf%fire_fuel_bulkd(beg:end)            = 0.0_r8 
    EDpcf%fire_fuel_eff_moist(beg:end)        = 0.0_r8 
    EDpcf%fire_fuel_sav(beg:end)              = 0.0_r8 
    EDpcf%fire_fuel_mef(beg:end)              = 0.0_r8     
    EDpcf%sum_fuel(beg:end)                   = 0.0_r8     
    EDpcf%litter_in(beg:end)                  = 0.0_r8 
    EDpcf%litter_out(beg:end)                 = 0.0_r8       
    EDpcf%efpot(beg:end)                      = 0.0_r8       
    EDpcf%rb(beg:end)                         = 0.0_r8
    EDpcf%prec24(beg:end)                     = 0.0_r8
    EDpcf%rh24(beg:end)                       = 0.0_r8
    EDpcf%wind24(beg:end)                     = 0.0_r8
    EDpcf%seed_bank(beg:end)                  = 0.0_r8
    EDpcf%seeds_in(beg:end)                   = 0.0_r8
    EDpcf%seed_decay(beg:end)                 = 0.0_r8    
    EDpcf%seed_germination(beg:end)           = 0.0_r8    

  end subroutine EDinit_pft_cflux_type

end module EDClmtypeInitMod
