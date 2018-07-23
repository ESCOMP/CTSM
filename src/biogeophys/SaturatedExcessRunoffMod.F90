module SaturatedExcessRunoffMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Type and associated routines for calculating surface runoff due to saturated surface
  !
  ! This also includes calculations of fsat (fraction of each column that is saturated)
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use decompMod    , only : bounds_type
  use abortutils   , only : endrun
  use clm_varctl   , only : iulog, use_vichydro
  use clm_varcon   , only : spval
  use ColumnType   , only : column_type
  use SoilHydrologyType, only : soilhydrology_type
  use SoilStateType, only : soilstate_type
  use WaterFluxBulkType, only : waterfluxbulk_type

  implicit none
  save
  private

  ! !PUBLIC TYPES:

  type, public :: saturated_excess_runoff_type
     private
     ! Public data members
     ! Note: these should be treated as read-only by other modules
     real(r8), pointer, public :: fsat_col(:) ! fractional area with water table at surface

     ! Private data members
     integer :: fsat_method
     real(r8), pointer :: fcov_col(:) ! fractional impermeable area
   contains
     ! Public routines
     procedure, public :: Init

     procedure, public :: SaturatedExcessRunoff ! Calculate surface runoff due to saturated surface

     ! Private routines
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold

     procedure, private, nopass :: ComputeFsatTopmodel
     procedure, private, nopass :: ComputeFsatVic
  end type saturated_excess_runoff_type

  ! !PRIVATE DATA MEMBERS:

  integer, parameter :: FSAT_METHOD_TOPMODEL = 1
  integer, parameter :: FSAT_METHOD_VIC      = 2

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  ! ========================================================================
  ! Infrastructure routines
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize this saturated_excess_runoff_type object
    !
    ! !ARGUMENTS:
    class(saturated_excess_runoff_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Allocate memory for this saturated_excess_runoff_type object
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(saturated_excess_runoff_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc

    character(len=*), parameter :: subname = 'InitAllocate'
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc

    allocate(this%fsat_col(begc:endc))                 ; this%fsat_col(:)                 = nan
    allocate(this%fcov_col(begc:endc))                 ; this%fcov_col(:)                 = nan   

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize saturated_excess_runoff_type history variables
    !
    ! !USES:
    use histFileMod , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(saturated_excess_runoff_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc

    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc

    this%fcov_col(begc:endc) = spval
    call hist_addfld1d (fname='FCOV',  units='unitless',  &
         avgflag='A', long_name='fractional impermeable area', &
         ptr_col=this%fcov_col, l2g_scale_type='veg')

    this%fsat_col(begc:endc) = spval
    call hist_addfld1d (fname='FSAT',  units='unitless',  &
         avgflag='A', long_name='fractional area with water table at surface', &
         ptr_col=this%fsat_col, l2g_scale_type='veg')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Perform cold-start initialization for saturated_excess_runoff_type
    !
    ! !ARGUMENTS:
    class(saturated_excess_runoff_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitCold'
    !-----------------------------------------------------------------------

    ! TODO(wjs, 2017-07-12) We'll read fsat_method from namelist.
    if (use_vichydro) then
       this%fsat_method = FSAT_METHOD_VIC
    else
       this%fsat_method = FSAT_METHOD_TOPMODEL
    end if

  end subroutine InitCold

  ! ========================================================================
  ! Science routines
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine SaturatedExcessRunoff (this, bounds, num_hydrologyc, filter_hydrologyc, &
       col, soilhydrology_inst, soilstate_inst, waterfluxbulk_inst)
    !
    ! !DESCRIPTION:
    ! Calculate surface runoff due to saturated surface
    !
    ! Sets this%fsat_col and waterfluxbulk_inst%qflx_sat_excess_surf_col
    !
    ! !ARGUMENTS:
    class(saturated_excess_runoff_type), intent(inout) :: this
    type(bounds_type)        , intent(in)    :: bounds               
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(column_type)        , intent(in)    :: col
    type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    type(waterfluxbulk_type)     , intent(inout) :: waterfluxbulk_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c

    character(len=*), parameter :: subname = 'SaturatedExcessRunoff'
    !-----------------------------------------------------------------------

    associate(                                                        & 
         fcov                   =>    this%fcov_col                          , & ! Output: [real(r8) (:)   ]  fractional impermeable area
         fsat                   =>    this%fsat_col                          , & ! Output: [real(r8) (:)   ]  fractional area with water table at surface

         snl                    =>    col%snl                                , & ! Input:  [integer  (:)   ]  minus number of snow layers

         qflx_sat_excess_surf   =>    waterfluxbulk_inst%qflx_sat_excess_surf_col, & ! Output: [real(r8) (:)   ]  surface runoff due to saturated surface (mm H2O /s)
         qflx_floodc            =>    waterfluxbulk_inst%qflx_floodc_col         , & ! Input:  [real(r8) (:)   ]  column flux of flood water from RTM
         qflx_rain_plus_snomelt => waterfluxbulk_inst%qflx_rain_plus_snomelt_col , & ! Input: [real(r8) (:)   ] rain plus snow melt falling on the soil (mm/s)

         origflag               =>    soilhydrology_inst%origflag            , & ! Input:  logical
         fracice                =>    soilhydrology_inst%fracice_col           & ! Input:  [real(r8) (:,:) ]  fractional impermeability (-)
         )

    ! ------------------------------------------------------------------------
    ! Compute fsat
    ! ------------------------------------------------------------------------

    select case (this%fsat_method)
    case (FSAT_METHOD_TOPMODEL)
       call this%ComputeFsatTopmodel(bounds, num_hydrologyc, filter_hydrologyc, &
            soilhydrology_inst, soilstate_inst, &
            fsat = fsat(bounds%begc:bounds%endc))
    case (FSAT_METHOD_VIC)
       call this%ComputeFsatVic(bounds, num_hydrologyc, filter_hydrologyc, &
            soilhydrology_inst, &
            fsat = fsat(bounds%begc:bounds%endc))
    case default
       write(iulog,*) subname//' ERROR: Unrecognized fsat_method: ', this%fsat_method
       call endrun(subname//' ERROR: Unrecognized fsat_method')
    end select

    ! ------------------------------------------------------------------------
    ! Compute qflx_sat_excess_surf
    !
    ! assume qinmax (maximum infiltration rate) is large relative to
    ! qflx_rain_plus_snomelt in control
    ! ------------------------------------------------------------------------
    
    if (origflag == 1) then
       if (this%fsat_method == FSAT_METHOD_VIC) then
          ! NOTE(wjs, 2017-07-12) I'm not sure if it's the VIC fsat method per se that
          ! is incompatible with origflag, or some other aspect of VIC: The original
          ! check was for origflag == 1 and use_vichydro, which also appears in error
          ! checks elsewhere.
          call endrun(msg="VICHYDRO is not available for origflag=1"//errmsg(sourcefile, __LINE__))
       end if
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          fcov(c) = (1._r8 - fracice(c,1)) * fsat(c) + fracice(c,1)
          qflx_sat_excess_surf(c) = fcov(c) * qflx_rain_plus_snomelt(c)
       end do
    else
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          ! only send fast runoff directly to streams
          qflx_sat_excess_surf(c) = fsat(c) * qflx_rain_plus_snomelt(c)

          ! Set fcov just to have it on the history file
          fcov(c) = fsat(c)
       end do
    end if

    ! ------------------------------------------------------------------------
    ! For urban columns, send flood water flux to runoff
    ! ------------------------------------------------------------------------

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       if (col%urbpoi(c)) then
          ! send flood water flux to runoff for all urban columns
          qflx_sat_excess_surf(c) = qflx_sat_excess_surf(c) + qflx_floodc(c)
       end if
    end do

    end associate

  end subroutine SaturatedExcessRunoff

  !-----------------------------------------------------------------------
  subroutine ComputeFsatTopmodel(bounds, num_hydrologyc, filter_hydrologyc, &
       soilhydrology_inst, soilstate_inst, fsat)
    !
    ! !DESCRIPTION:
    ! Compute fsat using the TOPModel-based parameterization
    !
    ! This is the CLM default parameterization
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer, intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(in) :: soilhydrology_inst
    type(soilstate_type), intent(in) :: soilstate_inst
    real(r8), intent(inout) :: fsat( bounds%begc: ) ! fractional area with water table at surface
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c
    real(r8) :: fff ! decay factor (m-1)

    character(len=*), parameter :: subname = 'ComputeFsatTopmodel'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(fsat) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    associate( &
         frost_table      =>    soilhydrology_inst%frost_table_col  , & ! Input:  [real(r8) (:)   ]  frost table depth (m)
         zwt              =>    soilhydrology_inst%zwt_col          , & ! Input:  [real(r8) (:)   ]  water table depth (m)
         zwt_perched      =>    soilhydrology_inst%zwt_perched_col  , & ! Input:  [real(r8) (:)   ]  perched water table depth (m)

         wtfact           =>    soilstate_inst%wtfact_col             & ! Input:  [real(r8) (:)   ]  maximum saturated fraction for a gridcell
         )

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       fff = 0.5_r8
       if (frost_table(c) > zwt_perched(c) .and. frost_table(c) <= zwt(c)) then
          ! use perched water table to determine fsat (if present)
          fsat(c) = wtfact(c) * exp(-0.5_r8*fff*zwt_perched(c))
       else
          fsat(c) = wtfact(c) * exp(-0.5_r8*fff*zwt(c))
       end if
    end do

    end associate

  end subroutine ComputeFsatTopmodel

  !-----------------------------------------------------------------------
  subroutine ComputeFsatVic(bounds, num_hydrologyc, filter_hydrologyc, &
       soilhydrology_inst, fsat)
    !
    ! !DESCRIPTION:
    ! Compute fsat using the VIC-based parameterization
    !
    ! Citation: Wood et al. 1992, "A land-surface hydrology parameterization with subgrid
    ! variability for general circulation models", JGR 97(D3), 2717-2728.
    !
    ! This implementation gives a first-order approximation to saturated excess runoff.
    ! For now we're not including the more exact analytical solution, or even a better
    ! numerical approximation.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer, intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(in) :: soilhydrology_inst
    real(r8), intent(inout) :: fsat( bounds%begc: ) ! fractional area with water table at surface
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c
    real(r8) :: ex(bounds%begc:bounds%endc) ! exponent

    character(len=*), parameter :: subname = 'ComputeFsatVic'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(fsat) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    associate( &
         b_infil          =>    soilhydrology_inst%b_infil_col      , & ! Input:  [real(r8) (:)   ]  VIC b infiltration parameter
         top_max_moist    =>    soilhydrology_inst%top_max_moist_col, & ! Input:  [real(r8) (:)   ]  maximum soil moisture in top VIC layers
         top_moist_limited =>   soilhydrology_inst%top_moist_limited_col & ! Input:  [real(r8) (:) ]  soil moisture in top layers, limited to no greater than top_max_moist
         )

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       ex(c) = b_infil(c) / (1._r8 + b_infil(c))
       ! fsat is equivalent to A in VIC papers
       fsat(c) = 1._r8 - (1._r8 - top_moist_limited(c) / top_max_moist(c))**ex(c)
    end do

    end associate

  end subroutine ComputeFsatVic


end module SaturatedExcessRunoffMod
