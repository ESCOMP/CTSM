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
  use WaterFluxBulkType    , only : waterfluxbulk_type
  use WaterDiagnosticBulkType, only : waterdiagnosticbulk_type

  implicit none
  save
  private

  ! !PUBLIC TYPES:

  type, public :: infiltration_excess_runoff_type
     private
     ! Public data members
     ! Note: these should be treated as read-only by other modules
     real(r8), pointer, public :: qinmax_col(:)  ! maximum infiltration rate (mm H2O /s)

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
  end type infiltration_excess_runoff_type

  ! !PRIVATE DATA MEMBERS:

  ! For methods that don't generate any infiltration excess runoff, we get this end result
  ! by specifying a huge qinmax value - i.e., an effectively infinite max infiltration
  ! rate.
  !
  ! 1e200 mm H2O /s seems large enough to be effectively infinite, while not being so
  ! large as to cause floating point overflows elsewhere.
  real(r8), parameter :: qinmax_unlimited = 1.e200_r8  ! mm H2O /s

  ! The 'none' option avoids generating any infiltration excess runoff by setting qinmax
  ! to a huge value
  integer, parameter :: QINMAX_METHOD_NONE = 0
  integer, parameter :: QINMAX_METHOD_HKSAT = 1

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
       this%qinmax_method = QINMAX_METHOD_NONE
    else
       this%qinmax_method = QINMAX_METHOD_HKSAT
    end if

  end subroutine InitCold

  ! ========================================================================
  ! Science routines
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine InfiltrationExcessRunoff(this, bounds, num_hydrologyc, filter_hydrologyc, &
       soilhydrology_inst, soilstate_inst, saturated_excess_runoff_inst, waterfluxbulk_inst, &
       waterdiagnosticbulk_inst)
    !
    ! !DESCRIPTION:
    ! Calculate surface runoff due to infiltration excess
    !
    ! Sets waterfluxbulk_inst%qflx_infl_excess_col and this%qinmax_col. These are valid within
    ! the hydrology filter. Both of these give averages over the entire column. However,
    ! qinmax is implicitly 0 over the fraction of the column given by fsat, and
    ! qflx_infl_excess is implicitly 0 over both fsat and frac_h2osfc.
    !
    ! !ARGUMENTS:
    class(infiltration_excess_runoff_type) , intent(in)    :: this
    type(bounds_type)                      , intent(in)    :: bounds               
    integer                                , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                                , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type)               , intent(in)    :: soilhydrology_inst
    type(soilstate_type)                   , intent(in)    :: soilstate_inst
    type(saturated_excess_runoff_type)     , intent(in)    :: saturated_excess_runoff_inst
    type(waterfluxbulk_type)                   , intent(in)    :: waterfluxbulk_inst
    type(waterdiagnosticbulk_type)         , intent(in)    :: waterdiagnosticbulk_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c
    real(r8) :: qinmax_on_unsaturated_area(bounds%begc:bounds%endc) ! maximum infiltration rate on the unsaturated fraction of the column (mm H2O /s)

    character(len=*), parameter :: subname = 'InfiltrationExcessRunoff'
    !-----------------------------------------------------------------------

    associate( &
         qinmax           =>    this%qinmax_col                     , & ! Output: [real(r8) (:)   ]  maximum infiltration rate (mm H2O /s)

         fsat             =>   saturated_excess_runoff_inst%fsat_col, & ! Input:  [real(r8) (:)   ]  fractional area with water table at surface       

         qflx_infl_excess =>    waterfluxbulk_inst%qflx_infl_excess_col , & ! Output: [real(r8) (:)   ]  infiltration excess runoff (mm H2O /s)
         qflx_in_soil     =>    waterfluxbulk_inst%qflx_in_soil_col     , & ! Input:  [real(r8) (:)   ]  surface input to soil (mm/s)

         frac_h2osfc      =>    waterdiagnosticbulk_inst%frac_h2osfc_col & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
          )

    select case (this%qinmax_method)
    case (QINMAX_METHOD_NONE)
       ! NOTE(wjs, 2017-09-01) I'm putting this here for clarity, though it could be
       ! moved to initialization for efficiency
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          qinmax_on_unsaturated_area(c) = qinmax_unlimited
       end do
    case (QINMAX_METHOD_HKSAT)
       call this%ComputeQinmaxHksat(bounds, num_hydrologyc, filter_hydrologyc, &
            soilhydrology_inst, soilstate_inst, &
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
    real(r8), intent(inout) :: qinmax_on_unsaturated_area( bounds%begc: ) ! maximum infiltration rate on the unsaturated fraction of the column (mm H2O /s)
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

end module InfiltrationExcessRunoffMod
