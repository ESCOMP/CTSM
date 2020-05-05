module CNProductsMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate loss fluxes from wood products pools, and update product pool state variables
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use shr_infnan_mod          , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use decompMod               , only : bounds_type
  use abortutils              , only : endrun
  use clm_time_manager        , only : get_step_size
  use SpeciesBaseType         , only : species_base_type
  use PatchType               , only : patch
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  type, public :: cn_products_type
     private
     ! ------------------------------------------------------------------------
     ! Public instance variables
     ! ------------------------------------------------------------------------

     real(r8), pointer, public :: product_loss_grc(:)   ! (g[C or N]/m2/s) total decomposition loss from ALL product pools

     ! ------------------------------------------------------------------------
     ! Private instance variables
     ! ------------------------------------------------------------------------

     class(species_base_type), allocatable :: species    ! C, N, C13, C14, etc.

     ! States
     real(r8), pointer :: cropprod1_grc(:)    ! (g[C or N]/m2) grain product pool, 1-year lifespan
     real(r8), pointer :: prod10_grc(:)       ! (g[C or N]/m2) wood product pool, 10-year lifespan
     real(r8), pointer :: prod100_grc(:)      ! (g[C or N]/m2) wood product pool, 100-year lifespan
     real(r8), pointer :: tot_woodprod_grc(:) ! (g[C or N]/m2) total wood product pool

     ! Fluxes: gains
     real(r8), pointer :: dwt_prod10_gain_grc(:)  ! (g[C or N]/m2/s) dynamic landcover addition to 10-year wood product pool
     real(r8), pointer :: dwt_prod100_gain_grc(:) ! (g[C or N]/m2/s) dynamic landcover addition to 100-year wood product pool
     real(r8), pointer :: dwt_woodprod_gain_grc(:) ! (g[C or N]/m2/s) dynamic landcover addition to wood product pools
     real(r8), pointer :: dwt_cropprod1_gain_grc(:) ! (g[C or N]/m2/s) dynamic landcover addition to 1-year crop product pool
     real(r8), pointer :: hrv_deadstem_to_prod10_patch(:)  ! (g[C or N]/m2/s) dead stem harvest to 10-year wood product pool
     real(r8), pointer :: hrv_deadstem_to_prod10_grc(:)  ! (g[C or N]/m2/s) dead stem harvest to 10-year wood product pool
     real(r8), pointer :: hrv_deadstem_to_prod100_patch(:) ! (g[C or N]/m2/s) dead stem harvest to 100-year wood product pool
     real(r8), pointer :: hrv_deadstem_to_prod100_grc(:) ! (g[C or N]/m2/s) dead stem harvest to 100-year wood product pool
     real(r8), pointer :: grain_to_cropprod1_patch(:) ! (g[C or N]/m2/s) grain to 1-year crop product pool
     real(r8), pointer :: grain_to_cropprod1_grc(:) ! (g[C or N]/m2/s) grain to 1-year crop product pool

     ! Fluxes: losses
     real(r8), pointer :: cropprod1_loss_grc(:)    ! (g[C or N]/m2/s) decomposition loss from 1-yr grain product pool
     real(r8), pointer :: prod10_loss_grc(:)       ! (g[C or N]/m2/s) decomposition loss from 10-yr wood product pool
     real(r8), pointer :: prod100_loss_grc(:)      ! (g[C or N]/m2/s) decomposition loss from 100-yr wood product pool
     real(r8), pointer :: tot_woodprod_loss_grc(:) ! (g[C or N]/m2/s) decompomposition loss from all wood product pools

   contains

     ! Infrastructure routines
     procedure, public  :: Init
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
     procedure, public  :: Restart

     ! Science routines
     procedure, public  :: UpdateProducts
     procedure, private :: PartitionWoodFluxes
     procedure, private :: PartitionGrainFluxes
     procedure, private :: ComputeSummaryVars

  end type cn_products_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds, species)
    ! !ARGUMENTS:
    class(cn_products_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds

    ! species tells whether this object is being used for C, N, C13, C14, etc. This is
    ! just used for naming history and restart fields
    class(species_base_type), intent(in) :: species

    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------

    allocate(this%species, source = species)

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    ! !ARGUMENTS:
    class(cn_products_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp,endp
    integer :: begg,endg

    character(len=*), parameter :: subname = 'InitAllocate'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp
    begg = bounds%begg
    endg = bounds%endg

    allocate(this%cropprod1_grc(begg:endg)) ; this%cropprod1_grc(:) = nan
    allocate(this%prod10_grc(begg:endg)) ; this%prod10_grc(:) = nan
    allocate(this%prod100_grc(begg:endg)) ; this%prod100_grc(:) = nan
    allocate(this%tot_woodprod_grc(begg:endg)) ; this%tot_woodprod_grc(:) = nan

    allocate(this%dwt_prod10_gain_grc(begg:endg)) ; this%dwt_prod10_gain_grc(:) = nan
    allocate(this%dwt_prod100_gain_grc(begg:endg)) ; this%dwt_prod100_gain_grc(:) = nan
    allocate(this%dwt_woodprod_gain_grc(begg:endg)) ; this%dwt_woodprod_gain_grc(:) = nan

    allocate(this%dwt_cropprod1_gain_grc(begg:endg)) ; this%dwt_cropprod1_gain_grc(:) = nan

    allocate(this%hrv_deadstem_to_prod10_patch(begp:endp)) ; this%hrv_deadstem_to_prod10_patch(:) = nan
    allocate(this%hrv_deadstem_to_prod10_grc(begg:endg)) ; this%hrv_deadstem_to_prod10_grc(:) = nan

    allocate(this%hrv_deadstem_to_prod100_patch(begp:endp)) ; this%hrv_deadstem_to_prod100_patch(:) = nan
    allocate(this%hrv_deadstem_to_prod100_grc(begg:endg)) ; this%hrv_deadstem_to_prod100_grc(:) = nan

    allocate(this%grain_to_cropprod1_patch(begp:endp)) ; this%grain_to_cropprod1_patch(:) = nan
    allocate(this%grain_to_cropprod1_grc(begg:endg)) ; this%grain_to_cropprod1_grc(:) = nan

    allocate(this%cropprod1_loss_grc(begg:endg)) ; this%cropprod1_loss_grc(:) = nan
    allocate(this%prod10_loss_grc(begg:endg)) ; this%prod10_loss_grc(:) = nan
    allocate(this%prod100_loss_grc(begg:endg)) ; this%prod100_loss_grc(:) = nan
    allocate(this%tot_woodprod_loss_grc(begg:endg)) ; this%tot_woodprod_loss_grc(:) = nan
    allocate(this%product_loss_grc(begg:endg)) ; this%product_loss_grc(:) = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    ! !USES:
    use histFileMod, only : hist_addfld1d
    use clm_varcon , only : spval
    !
    ! !ARGUMENTS:
    class(cn_products_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begg,endg
    character(len=:), allocatable :: active_if_non_isotope

    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------

    begg = bounds%begg
    endg = bounds%endg

    if (this%species%is_isotope()) then
       active_if_non_isotope = 'inactive'
    else
       active_if_non_isotope = 'active'
    end if

    this%cropprod1_grc(begg:endg) = spval
    call hist_addfld1d( &
         fname = this%species%hist_fname('CROPPROD1'), &
         units = 'g' // this%species%get_species() // '/m^2', &
         avgflag = 'A', &
         long_name = '1-yr grain product ' // this%species%get_species(), &
         ptr_gcell = this%cropprod1_grc, default=active_if_non_isotope)

    this%prod10_grc(begg:endg) = spval
    call hist_addfld1d( &
         fname = this%species%hist_fname('PROD10'), &
         units = 'g' // this%species%get_species() // '/m^2', &
         avgflag = 'A', &
         long_name = '10-yr wood product ' // this%species%get_species(), &
         ptr_gcell = this%prod10_grc, default='inactive')

    this%prod100_grc(begg:endg) = spval
    call hist_addfld1d( &
         fname = this%species%hist_fname('PROD100'), &
         units = 'g' // this%species%get_species() // '/m^2', &
         avgflag = 'A', &
         long_name = '100-yr wood product ' // this%species%get_species(), &
         ptr_gcell = this%prod100_grc, default='inactive')

    this%tot_woodprod_grc(begg:endg) = spval
    call hist_addfld1d( &
         fname = this%species%hist_fname('TOT_WOODPROD'), &
         units = 'g' // this%species%get_species() // '/m^2', &
         avgflag = 'A', &
         long_name = 'total wood product ' // this%species%get_species(), &
         ptr_gcell = this%tot_woodprod_grc, default=active_if_non_isotope)

    this%dwt_prod10_gain_grc(begg:endg) = spval
    call hist_addfld1d( &
         fname = this%species%hist_fname('DWT_PROD10', suffix='_GAIN'), &
         units = 'g' // this%species%get_species() // '/m^2/s', &
         avgflag = 'A', &
         long_name = 'landcover change-driven addition to 10-yr wood product pool', &
         ptr_gcell = this%dwt_prod10_gain_grc, default='inactive')

    this%dwt_prod100_gain_grc(begg:endg) = spval
    call hist_addfld1d( &
         fname = this%species%hist_fname('DWT_PROD100', suffix='_GAIN'), &
         units = 'g' // this%species%get_species() // '/m^2/s', &
         avgflag = 'A', &
         long_name = 'landcover change-driven addition to 100-yr wood product pool', &
         ptr_gcell = this%dwt_prod100_gain_grc, default='inactive')

    this%dwt_woodprod_gain_grc(begg:endg) = spval
    call hist_addfld1d( &
         fname = this%species%hist_fname('DWT_WOODPROD', suffix='_GAIN'), &
         units = 'g' // this%species%get_species() // '/m^2/s', &
         avgflag = 'A', &
         long_name = 'landcover change-driven addition to wood product pools', &
         ptr_gcell = this%dwt_woodprod_gain_grc, default=active_if_non_isotope)

    this%dwt_cropprod1_gain_grc(begg:endg) = spval
    call hist_addfld1d( &
         fname = this%species%hist_fname('DWT_CROPPROD1', suffix='_GAIN'), &
         units = 'g' // this%species%get_species() // '/m^2/s', &
         avgflag = 'A', &
         long_name = 'landcover change-driven addition to 1-year crop product pool', &
         ptr_gcell = this%dwt_cropprod1_gain_grc, default=active_if_non_isotope)

    this%cropprod1_loss_grc(begg:endg) = spval
    call hist_addfld1d( &
         fname = this%species%hist_fname('CROPPROD1', suffix='_LOSS'), &
         units = 'g' // this%species%get_species() // '/m^2/s', &
         avgflag = 'A', &
         long_name = 'loss from 1-yr grain product pool', &
         ptr_gcell = this%cropprod1_loss_grc, default=active_if_non_isotope)

    this%prod10_loss_grc(begg:endg) = spval
    call hist_addfld1d( &
         fname = this%species%hist_fname('PROD10', suffix='_LOSS'), &
         units = 'g' // this%species%get_species() // '/m^2/s', &
         avgflag = 'A', &
         long_name = 'loss from 10-yr wood product pool', &
         ptr_gcell = this%prod10_loss_grc, default='inactive')

    this%prod100_loss_grc(begg:endg) = spval
    call hist_addfld1d( &
         fname = this%species%hist_fname('PROD100', suffix='_LOSS'), &
         units = 'g' // this%species%get_species() // '/m^2/s', &
         avgflag = 'A', &
         long_name = 'loss from 100-yr wood product pool', &
         ptr_gcell = this%prod100_loss_grc, default='inactive')

    this%tot_woodprod_loss_grc(begg:endg) = spval
    call hist_addfld1d( &
         fname = this%species%hist_fname('TOT_WOODPROD', suffix='_LOSS'), &
         units = 'g' // this%species%get_species() // '/m^2/s', &
         avgflag = 'A', &
         long_name = 'total loss from wood product pools', &
         ptr_gcell = this%tot_woodprod_loss_grc, default=active_if_non_isotope)

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    ! !ARGUMENTS:
    class(cn_products_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g, p

    character(len=*), parameter :: subname = 'InitCold'
    !-----------------------------------------------------------------------

    do g = bounds%begg, bounds%endg
       this%cropprod1_grc(g) = 0._r8
       this%prod10_grc(g) = 0._r8
       this%prod100_grc(g) = 0._r8
       this%tot_woodprod_grc(g) = 0._r8
    end do

    ! Need to set these patch-level fluxes to 0 everywhere for the sake of special
    ! landunits (because they don't get set over special landunits in the run loop)
    do p = bounds%begp, bounds%endp
       this%hrv_deadstem_to_prod10_patch(p) = 0._r8
       this%hrv_deadstem_to_prod100_patch(p) = 0._r8
       this%grain_to_cropprod1_patch(p) = 0._r8
    end do

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag, &
       template_for_missing_fields, template_multiplier)
    ! !USES:
    use ncdio_pio  , only : file_desc_t, ncd_double
    use restUtilMod, only : restartvar, set_missing_from_template, set_grc_field_from_col_field
    !
    ! !ARGUMENTS:
    class(cn_products_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    type(file_desc_t), intent(inout) :: ncid
    character(len=*), intent(in) :: flag  ! 'read' or 'write'

    ! If template_for_missing_fields and template_multiplier are provided, then: When
    ! reading the restart file, for any field not present on the restart file, the field
    ! in this object is set equal to the corresponding field in
    ! template_for_missing_fields times template_multiplier.
    !
    ! The Restart routine must have been called on template_for_missing_fields before
    ! calling it on this object.
    ! 
    ! (Must provide both template_for_missing_fields and template_multiplier or neither)
    class(cn_products_type), optional, intent(in) :: template_for_missing_fields
    real(r8), optional, intent(in) :: template_multiplier

    !
    ! !LOCAL VARIABLES:
    logical :: template_provided
    logical :: readvar

    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------

    if (present(template_for_missing_fields) .and. present(template_multiplier)) then
       template_provided = .true.
    else if (present(template_for_missing_fields)) then
       call endrun(&
            msg='template_for_missing_fields provided; must also provide template_multiplier' // &
            errMsg(sourcefile, __LINE__))
    else if (present(template_multiplier)) then
       call endrun(&
            msg='template_multiplier provided; must also provide template_for_missing_fields' // &
            errMsg(sourcefile, __LINE__))
    else
       template_provided = .false.
    end if

    ! NOTE(wjs, 2016-03-29) Adding '_g' suffixes to the end of the restart field names to
    ! distinguish these gridcell-level restart fields from the obsolete column-level
    ! restart fields that are present on old restart files.

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%species%rest_fname('cropprod1', suffix='_g'), &
         xtype=ncd_double, dim1name='gridcell', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%cropprod1_grc)
    if (flag == 'read' .and. .not. readvar) then
       ! BACKWARDS_COMPATIBILITY(wjs, 2016-03-31) If the gridcell-level field isn't
       ! present, try to find a column-level field (which may be present on an older
       ! restart file).
       call set_grc_field_from_col_field( &
            bounds = bounds, &
            ncid = ncid, &
            varname = this%species%rest_fname('cropprod1'), &
            data_grc = this%cropprod1_grc, &
            readvar = readvar)

       ! If we still haven't found an appropriate field on the restart file, then set
       ! this field from the template, if provided
       if (.not. readvar .and. template_provided) then
          call set_missing_from_template(this%cropprod1_grc, &
               template_for_missing_fields%cropprod1_grc, &
               multiplier = template_multiplier)
       end if
    end if

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%species%rest_fname('prod10', suffix='_g'), &
         xtype=ncd_double, dim1name='gridcell', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod10_grc)
    if (flag == 'read' .and. .not. readvar) then
       ! BACKWARDS_COMPATIBILITY(wjs, 2016-03-31) If the gridcell-level field isn't
       ! present, try to find a column-level field (which may be present on an older
       ! restart file).
       call set_grc_field_from_col_field( &
            bounds = bounds, &
            ncid = ncid, &
            varname = this%species%rest_fname('prod10'), &
            data_grc = this%prod10_grc, &
            readvar = readvar)

       ! If we still haven't found an appropriate field on the restart file, then set
       ! this field from the template, if provided
       if (.not. readvar .and. template_provided) then
          call set_missing_from_template(this%prod10_grc, &
               template_for_missing_fields%prod10_grc, &
               multiplier = template_multiplier)
       end if
    end if

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%species%rest_fname('prod100', suffix='_g'), &
         xtype=ncd_double, dim1name='gridcell', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod100_grc)
    if (flag == 'read' .and. .not. readvar) then
       ! BACKWARDS_COMPATIBILITY(wjs, 2016-03-31) If the gridcell-level field isn't
       ! present, try to find a column-level field (which may be present on an older
       ! restart file).
       call set_grc_field_from_col_field( &
            bounds = bounds, &
            ncid = ncid, &
            varname = this%species%rest_fname('prod100'), &
            data_grc = this%prod100_grc, &
            readvar = readvar)

       ! If we still haven't found an appropriate field on the restart file, then set
       ! this field from the template, if provided
       if (.not. readvar .and. template_provided) then
          call set_missing_from_template(this%prod100_grc, &
               template_for_missing_fields%prod100_grc, &
               multiplier = template_multiplier)
       end if
    end if

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine UpdateProducts(this, bounds, &
       num_soilp, filter_soilp, &
       dwt_wood_product_gain_patch, &
       wood_harvest_patch, &
       dwt_crop_product_gain_patch, &
       grain_to_cropprod_patch)
    !
    ! !DESCRIPTION:
    ! Update all loss fluxes from wood and grain product pools, and update product pool
    ! state variables for both loss and gain terms
    !
    ! !ARGUMENTS:
    class(cn_products_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                 , intent(in)    :: filter_soilp(:) ! filter for soil patches

    ! dynamic landcover addition to wood product pools (g/m2/s) [patch]; although this is
    ! a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), intent(in) :: dwt_wood_product_gain_patch( bounds%begp: )

    ! wood harvest addition to wood product pools (g/m2/s) [patch]
    real(r8), intent(in) :: wood_harvest_patch( bounds%begp: )

    ! dynamic landcover addition to crop product pools (g/m2/s) [patch]; although this is
    ! a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), intent(in) :: dwt_crop_product_gain_patch( bounds%begp: )

    ! grain to crop product pool (g/m2/s) [patch]
    real(r8), intent(in) :: grain_to_cropprod_patch( bounds%begp: )
    !
    ! !LOCAL VARIABLES:
    integer  :: g        ! indices
    real(r8) :: dt       ! time step (seconds)
    real(r8) :: kprod1   ! decay constant for 1-year product pool
    real(r8) :: kprod10  ! decay constant for 10-year product pool
    real(r8) :: kprod100 ! decay constant for 100-year product pool
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(dwt_wood_product_gain_patch) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(wood_harvest_patch) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_crop_product_gain_patch) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(grain_to_cropprod_patch) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))

    call this%PartitionWoodFluxes(bounds, &
         num_soilp, filter_soilp, &
         dwt_wood_product_gain_patch(bounds%begp:bounds%endp), &
         wood_harvest_patch(bounds%begp:bounds%endp))

    call this%PartitionGrainFluxes(bounds, &
         num_soilp, filter_soilp, &
         dwt_crop_product_gain_patch(bounds%begp:bounds%endp), &
         grain_to_cropprod_patch(bounds%begp:bounds%endp))

    ! calculate losses from product pools
    ! the following (1/s) rate constants result in ~90% loss of initial state over 1, 10 and 100 years,
    ! respectively, using a discrete-time fractional decay algorithm.
    kprod1  = 7.2e-8
    kprod10 = 7.2e-9
    kprod100 = 7.2e-10

    do g = bounds%begg, bounds%endg
       ! calculate fluxes out of product pools (1/sec)
       this%cropprod1_loss_grc(g) = this%cropprod1_grc(g) * kprod1
       this%prod10_loss_grc(g)    = this%prod10_grc(g)    * kprod10
       this%prod100_loss_grc(g)   = this%prod100_grc(g)   * kprod100
    end do

    ! set time steps
    dt = real( get_step_size(), r8 )

    ! update product state variables
    do g = bounds%begg, bounds%endg

       ! fluxes into wood & grain product pools, from landcover change
       this%cropprod1_grc(g) = this%cropprod1_grc(g) + this%dwt_cropprod1_gain_grc(g)*dt
       this%prod10_grc(g)    = this%prod10_grc(g)    + this%dwt_prod10_gain_grc(g)*dt
       this%prod100_grc(g)   = this%prod100_grc(g)   + this%dwt_prod100_gain_grc(g)*dt

       ! fluxes into wood & grain product pools, from harvest
       this%cropprod1_grc(g) = this%cropprod1_grc(g) + this%grain_to_cropprod1_grc(g)*dt
       this%prod10_grc(g)    = this%prod10_grc(g)    + this%hrv_deadstem_to_prod10_grc(g)*dt
       this%prod100_grc(g)   = this%prod100_grc(g)   + this%hrv_deadstem_to_prod100_grc(g)*dt

       ! fluxes out of wood & grain product pools, from decomposition
       this%cropprod1_grc(g) = this%cropprod1_grc(g) - this%cropprod1_loss_grc(g)*dt
       this%prod10_grc(g)    = this%prod10_grc(g)    - this%prod10_loss_grc(g)*dt
       this%prod100_grc(g)   = this%prod100_grc(g)   - this%prod100_loss_grc(g)*dt

    end do

    call this%ComputeSummaryVars(bounds)

  end subroutine UpdateProducts

  !-----------------------------------------------------------------------
  subroutine PartitionWoodFluxes(this, bounds, &
       num_soilp, filter_soilp, &
       dwt_wood_product_gain_patch, &
       wood_harvest_patch)
    !
    ! !DESCRIPTION:
    ! Partition input wood fluxes into 10 and 100 year product pools
    !
    ! !USES:
    use pftconMod    , only : pftcon
    use subgridAveMod, only : p2g
    !
    ! !ARGUMENTS:
    class(cn_products_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                 , intent(in)    :: filter_soilp(:) ! filter for soil patches

    ! dynamic landcover addition to wood product pools (g/m2/s) [patch]; although this is
    ! a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), intent(in) :: dwt_wood_product_gain_patch( bounds%begp: )

    ! wood harvest addition to wood product pools (g/m2/s) [patch]
    real(r8), intent(in) :: wood_harvest_patch( bounds%begp: )

    !
    ! !LOCAL VARIABLES:
    integer :: fp
    integer :: p
    integer :: g
    real(r8) :: pprod10       ! PFT proportion of deadstem to 10-year product pool
    real(r8) :: pprod100      ! PFT proportion of deadstem to 100-year product pool
    real(r8) :: pprod_tot     ! PFT proportion of deadstem to any product pool
    real(r8) :: pprod10_frac  ! PFT fraction of deadstem to product pool that goes to 10-year product pool
    real(r8) :: pprod100_frac ! PFT fraction of deadstem to product pool that goes to 100-year product pool

    character(len=*), parameter :: subname = 'PartitionWoodFluxes'
    !-----------------------------------------------------------------------

    ! Partition patch-level harvest fluxes to 10 and 100-year product pools
    do fp = 1, num_soilp
       p = filter_soilp(fp)
       this%hrv_deadstem_to_prod10_patch(p)  = &
            wood_harvest_patch(p) * pftcon%pprodharv10(patch%itype(p))
       this%hrv_deadstem_to_prod100_patch(p) = &
            wood_harvest_patch(p) * (1.0_r8 - pftcon%pprodharv10(patch%itype(p)))
    end do

    ! Average harvest fluxes from patch to gridcell
    call p2g(bounds, &
         this%hrv_deadstem_to_prod10_patch(bounds%begp:bounds%endp), &
         this%hrv_deadstem_to_prod10_grc(bounds%begg:bounds%endg), &
         p2c_scale_type = 'unity', &
         c2l_scale_type = 'unity', &
         l2g_scale_type = 'unity')

    call p2g(bounds, &
         this%hrv_deadstem_to_prod100_patch(bounds%begp:bounds%endp), &
         this%hrv_deadstem_to_prod100_grc(bounds%begg:bounds%endg), &
         p2c_scale_type = 'unity', &
         c2l_scale_type = 'unity', &
         l2g_scale_type = 'unity')

    ! Zero the dwt gains
    do g = bounds%begg, bounds%endg
       this%dwt_prod10_gain_grc(g) = 0._r8
       this%dwt_prod100_gain_grc(g) = 0._r8
    end do

    ! Partition dynamic land cover fluxes to 10 and 100-year product pools.
    do p = bounds%begp, bounds%endp
       g = patch%gridcell(p)

       ! Note that pprod10 + pprod100 do NOT sum to 1: some fraction of the dwt changes
       ! was lost to other fluxes. dwt_wood_product_gain_patch gives the amount that goes
       ! to all product pools, so we need to determine the fraction of that flux that
       ! goes to each pool.
       pprod10 = pftcon%pprod10(patch%itype(p))
       pprod100 = pftcon%pprod100(patch%itype(p))
       pprod_tot = pprod10 + pprod100
       if (pprod_tot > 0) then
          pprod10_frac = pprod10 / pprod_tot
          pprod100_frac = pprod100 / pprod_tot
       else
          ! Avoid divide by 0
          pprod10_frac = 0._r8
          pprod100_frac = 0._r8
       end if

       ! Note that the patch-level fluxes are expressed per unit gridcell area. So, to go
       ! from patch-level fluxes to gridcell-level fluxes, we simply add up the various
       ! patch contributions, without having to multiply by any area weightings.
       this%dwt_prod10_gain_grc(g) = this%dwt_prod10_gain_grc(g) + &
            dwt_wood_product_gain_patch(p) * pprod10_frac
       this%dwt_prod100_gain_grc(g) = this%dwt_prod100_gain_grc(g) + &
            dwt_wood_product_gain_patch(p) * pprod100_frac
    end do

  end subroutine PartitionWoodFluxes

  !-----------------------------------------------------------------------
  subroutine PartitionGrainFluxes(this, bounds, &
       num_soilp, filter_soilp, &
       dwt_crop_product_gain_patch, &
       grain_to_cropprod_patch)
    !
    ! !DESCRIPTION:
    ! Partition input grain fluxes into crop product pools
    !
    ! For now this doesn't do much, since there is just a single (1-year) crop product
    ! pool. But this provides the capability to add different crop product pools in the
    ! future, without requiring any changes to code outside of this class. It also gives
    ! symmetry with the wood fluxes.
    !
    ! !USES:
    use subgridAveMod, only : p2g
    !
    ! !ARGUMENTS:
    class(cn_products_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                 , intent(in)    :: filter_soilp(:) ! filter for soil patches

    ! dynamic landcover addition to crop product pool (g/m2/s) [patch]; although this is
    ! a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), intent(in) :: dwt_crop_product_gain_patch( bounds%begp: )

    ! grain to crop product pool(s) (g/m2/s) [patch]
    real(r8)                , intent(in)    :: grain_to_cropprod_patch( bounds%begp: )
    !
    ! !LOCAL VARIABLES:
    integer :: fp
    integer :: p
    integer :: g

    character(len=*), parameter :: subname = 'PartitionGrainFluxes'
    !-----------------------------------------------------------------------

    ! Determine gains from crop harvest

    do fp = 1, num_soilp
       p = filter_soilp(fp)

       ! For now all crop product is put in the 1-year crop product pool
       this%grain_to_cropprod1_patch(p) = grain_to_cropprod_patch(p)
    end do

    call p2g(bounds, &
         this%grain_to_cropprod1_patch(bounds%begp:bounds%endp), &
         this%grain_to_cropprod1_grc(bounds%begg:bounds%endg), &
         p2c_scale_type = 'unity', &
         c2l_scale_type = 'unity', &
         l2g_scale_type = 'unity')

    ! Determine gains from dynamic landcover

    do g = bounds%begg, bounds%endg
       this%dwt_cropprod1_gain_grc(g) = 0._r8
    end do

    do p = bounds%begp, bounds%endp
       g = patch%gridcell(p)

       ! Note that the patch-level fluxes are expressed per unit gridcell area. So, to go
       ! from patch-level fluxes to gridcell-level fluxes, we simply add up the various
       ! patch contributions, without having to multiply by any area weightings.
       this%dwt_cropprod1_gain_grc(g) = this%dwt_cropprod1_gain_grc(g) + &
            dwt_crop_product_gain_patch(p)
    end do

  end subroutine PartitionGrainFluxes


  !-----------------------------------------------------------------------
  subroutine ComputeSummaryVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Compute summary variables in this object: sums across multiple product pools
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_products_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: g        ! indices

    character(len=*), parameter :: subname = 'ComputeSummaryVars'
    !-----------------------------------------------------------------------

    do g = bounds%begg, bounds%endg

       ! total wood products
       this%tot_woodprod_grc(g) = &
            this%prod10_grc(g) + &
            this%prod100_grc(g)

       ! total loss from wood products
       this%tot_woodprod_loss_grc(g) = &
            this%prod10_loss_grc(g) + &
            this%prod100_loss_grc(g)

       ! total loss from ALL products
       this%product_loss_grc(g) = &
            this%cropprod1_loss_grc(g) + &
            this%prod10_loss_grc(g) + &
            this%prod100_loss_grc(g)

       this%dwt_woodprod_gain_grc(g) = &
            this%dwt_prod100_gain_grc(g) + &
            this%dwt_prod10_gain_grc(g)
    end do

  end subroutine ComputeSummaryVars


end module CNProductsMod
