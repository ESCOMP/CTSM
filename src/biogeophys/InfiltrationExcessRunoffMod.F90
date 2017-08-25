module InfiltrationExcessRunoffMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Type and associated routines for computing infiltration excess runoff and related
  ! variables
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use decompMod        , only : bounds_type
  use abortutils       , only : endrun
  use clm_varctl       , only : iulog, use_vichydro
  use clm_varcon       , only : spval, e_ice
  use clm_time_manager , only : get_step_size
  use SoilHydrologyType, only : soilhydrology_type
  use SoilStateType    , only : soilstate_type
  use SaturatedExcessRunoffMod, only : saturated_excess_runoff_type
  use WaterfluxType    , only : waterflux_type
  use WaterstateType   , only : waterstate_type

  implicit none
  save
  private

  ! !PUBLIC TYPES:

  type, public :: infiltration_excess_runoff_type
     private
     ! Public data members
     ! Note: these should be treated as read-only by other modules

     ! These are valid within the hydrology filter.
     !
     ! Both of these give averages over the entire column. However, qinmax is implicitly
     ! 0 over the fraction of the column given by fsat, and qflx_infl_excess is
     ! implicitly 0 over both fsat and frac_h2osfc.
     real(r8), pointer, public :: qinmax_col(:)  ! maximum infiltration capacity (mm H2O /s)
     real(r8), pointer, public :: qflx_infl_excess_col(:)  ! infiltration excess runoff (mm H2O /s)

     ! Private data members
     integer :: qinmax_method
   contains
     ! Public routines
     procedure, public :: Init

     procedure, public :: InfiltrationExcessRunoff ! Calculate surface runoff due to infiltration excess

     ! Private routines
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold

     procedure, private, nopass :: ComputeQinmaxHksat
     procedure, private, nopass :: ComputeQinmaxVic
  end type infiltration_excess_runoff_type

  ! !PRIVATE DATA MEMBERS:

  integer, parameter :: QINMAX_METHOD_HKSAT = 1
  integer, parameter :: QINMAX_METHOD_VIC  = 2

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
    ! Initialize this infiltration_excess_runoff_type object
    !
    ! !ARGUMENTS:
    class(infiltration_excess_runoff_type), intent(inout) :: this
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
    ! Allocate memory for this infiltration_excess_runoff_type object
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(infiltration_excess_runoff_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc

    character(len=*), parameter :: subname = 'InitAllocate'
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    allocate(this%qinmax_col          (begc:endc)); this%qinmax_col          (:) = nan
    allocate(this%qflx_infl_excess_col(begc:endc)); this%qflx_infl_excess_col(:) = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize infiltration_excess_runoff_type history variables
    !
    ! !USES:
    use histFileMod , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(infiltration_excess_runoff_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------

    ! Nothing to do for now

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Perform cold-start initialization for infiltration_excess_runoff_type
    !
    ! !ARGUMENTS:
    class(infiltration_excess_runoff_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitCold'
    !-----------------------------------------------------------------------

    ! TODO(wjs, 2017-08-14) We'll read qinmax_method from namelist.
    if (use_vichydro) then
       this%qinmax_method = QINMAX_METHOD_VIC
    else
       this%qinmax_method = QINMAX_METHOD_HKSAT
    end if

  end subroutine InitCold

  ! ========================================================================
  ! Science routines
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine InfiltrationExcessRunoff(this, bounds, num_hydrologyc, filter_hydrologyc, &
       soilhydrology_inst, soilstate_inst, saturated_excess_runoff_inst, waterflux_inst, waterstate_inst)
    !
    ! !DESCRIPTION:
    ! Calculate surface runoff due to infiltration excess
    !
    ! !ARGUMENTS:
    class(infiltration_excess_runoff_type) , intent(inout) :: this
    type(bounds_type)                      , intent(in)    :: bounds               
    integer                                , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                                , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type)               , intent(in)    :: soilhydrology_inst
    type(soilstate_type)                   , intent(in)    :: soilstate_inst
    type(saturated_excess_runoff_type)     , intent(in)    :: saturated_excess_runoff_inst
    type(waterflux_type)                   , intent(in)    :: waterflux_inst
    type(waterstate_type)                  , intent(in)    :: waterstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c
    real(r8) :: qinmax_on_unsaturated_area(bounds%begc:bounds%endc) ! maximum infiltration capacity on the unsaturated fraction of the column (mm H2O /s)

    character(len=*), parameter :: subname = 'InfiltrationExcessRunoff'
    !-----------------------------------------------------------------------

    associate( &
         qinmax           =>    this%qinmax_col                     , & ! Output: [real(r8) (:)   ]  maximum infiltration capacity (mm H2O /s)
         qflx_infl_excess =>    this%qflx_infl_excess_col           , & ! Output: [real(r8) (:)   ]  infiltration excess runoff (mm H2O /s)

         fsat             =>   saturated_excess_runoff_inst%fsat_col, & ! Input:  [real(r8) (:)   ]  fractional area with water table at surface       

         qflx_in_soil     =>    waterflux_inst%qflx_in_soil_col     , & ! Input:  [real(r8) (:)   ]  surface input to soil (mm/s)

         frac_h2osfc      =>    waterstate_inst%frac_h2osfc_col       & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
          )

    select case (this%qinmax_method)
    case (QINMAX_METHOD_HKSAT)
       call this%ComputeQinmaxHksat(bounds, num_hydrologyc, filter_hydrologyc, &
            soilhydrology_inst, soilstate_inst, &
            qinmax_on_unsaturated_area = qinmax_on_unsaturated_area(bounds%begc:bounds%endc))
    case (QINMAX_METHOD_VIC)
       call this%ComputeQinmaxVic(bounds, num_hydrologyc, filter_hydrologyc, &
            soilhydrology_inst, &
            fsat = fsat(bounds%begc:bounds%endc), &
            qflx_in_soil = qflx_in_soil(bounds%begc:bounds%endc), &
            qinmax_on_unsaturated_area = qinmax_on_unsaturated_area(bounds%begc:bounds%endc))
    case default
       write(iulog,*) subname//' ERROR: Unrecognized qinmax_method: ', this%qinmax_method
       call endrun(subname//' ERROR: Unrecognized qinmax_method')
    end select

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       qinmax(c) = (1._r8 - fsat(c)) * qinmax_on_unsaturated_area(c)
       qflx_infl_excess(c) = max(0._r8, &
            (qflx_in_soil(c) - (1.0_r8 - frac_h2osfc(c))*qinmax(c)))
    end do

    end associate

  end subroutine InfiltrationExcessRunoff

  !-----------------------------------------------------------------------
  subroutine ComputeQinmaxHksat(bounds, num_hydrologyc, filter_hydrologyc, &
       soilhydrology_inst, soilstate_inst, &
       qinmax_on_unsaturated_area)
    !
    ! !DESCRIPTION:
    ! Compute qinmax using a parameterization based on hksat
    !
    ! This is the CLM default parameterization
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer, intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(in) :: soilhydrology_inst
    type(soilstate_type), intent(in) :: soilstate_inst
    real(r8), intent(inout) :: qinmax_on_unsaturated_area( bounds%begc: ) ! maximum infiltration capacity on the unsaturated fraction of the column (mm H2O /s)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c

    character(len=*), parameter :: subname = 'ComputeQinmaxHksat'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(qinmax_on_unsaturated_area) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    associate( &
         icefrac          =>    soilhydrology_inst%icefrac_col      , & ! Input:  [real(r8) (:,:) ]  fraction of ice

         hksat            =>    soilstate_inst%hksat_col              & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
         )

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       qinmax_on_unsaturated_area(c) = minval(10._r8**(-e_ice*(icefrac(c,1:3)))*hksat(c,1:3))
    end do

    end associate

  end subroutine ComputeQinmaxHksat

  !-----------------------------------------------------------------------
  subroutine ComputeQinmaxVic(bounds, num_hydrologyc, filter_hydrologyc, &
       soilhydrology_inst, &
       fsat, qflx_in_soil, qinmax_on_unsaturated_area)
    !
    ! !DESCRIPTION:
    ! Compute qinmax using the VIC parameterization
    !
    ! Citation: Wood et al. 1992, "A land-surface hydrology parameterization with subgrid
    ! variability for general circulation models", JGR 97(D3), 2717-2728.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer, intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(in) :: soilhydrology_inst
    real(r8) , intent(in)    :: fsat( bounds%begc: )                       ! fractional area with water table at surface
    real(r8) , intent(in)    :: qflx_in_soil( bounds%begc: )               ! surface input to soil (mm/s)
    real(r8) , intent(inout) :: qinmax_on_unsaturated_area( bounds%begc: ) ! maximum infiltration capacity on the unsaturated fraction of the column (mm H2O /s)
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c
    real(r8) :: dtime       ! land model time step (sec)
    real(r8) :: top_icefrac ! ice fraction in top VIC layers
    real(r8) :: max_infil   ! max infiltration capacity in VIC (mm)
    real(r8) :: i_0         ! average soil moisture in top VIC layers (mm)
    real(r8) :: rsurf_vic   ! VIC surface runoff
    real(r8) :: basis       ! variable soil moisture holding capacity in top VIC layers for runoff calculation

    character(len=*), parameter :: subname = 'ComputeQinmaxVic'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(fsat) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(qflx_in_soil) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(qinmax_on_unsaturated_area) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    associate( &
          top_max_moist    =>    soilhydrology_inst%top_max_moist_col, & ! Input:  [real(r8) (:)   ]  maximum soil moisture in top VIC layers
          top_moist        =>    soilhydrology_inst%top_moist_col    , & ! Input:  [real(r8) (:)   ]  soil moisture in top VIC layers
          top_ice          =>    soilhydrology_inst%top_ice_col      , & ! Input:  [real(r8) (:)   ]  ice len in top VIC layers
          b_infil          =>    soilhydrology_inst%b_infil_col        & ! Input:  [real(r8) (:)   ]  VIC b infiltration parameter
          )

    dtime = get_step_size()

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       top_icefrac = min(1._r8,top_ice(c)/top_max_moist(c))
       max_infil = (1._r8+b_infil(c)) * top_max_moist(c)
       i_0 = max_infil * (1._r8 - (1._r8 - fsat(c))**(1._r8/b_infil(c)))
       if(qflx_in_soil(c) <= 0._r8) then
          rsurf_vic = 0._r8
       else if(max_infil <= 0._r8) then
          rsurf_vic = qflx_in_soil(c)
       else if((i_0 + qflx_in_soil(c)*dtime) > max_infil) then             !(Eq.(3a) Wood et al. 1992)
          rsurf_vic = (qflx_in_soil(c)*dtime - top_max_moist(c) + top_moist(c))/dtime
       else                                                                !(Eq.(3b) Wood et al. 1992)
          basis = 1._r8 - (i_0 + qflx_in_soil(c)*dtime)/max_infil
          rsurf_vic = (qflx_in_soil(c)*dtime - top_max_moist(c) + top_moist(c)    &
               + top_max_moist(c) * basis**(1._r8 + b_infil(c)))/dtime
       end if
       rsurf_vic = min(qflx_in_soil(c), rsurf_vic)
       qinmax_on_unsaturated_area(c) = 10._r8**(-e_ice*top_icefrac)*(qflx_in_soil(c) - rsurf_vic)
    end do

    end associate

  end subroutine ComputeQinmaxVic

end module InfiltrationExcessRunoffMod
