module WaterType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Container for derived types relating to water, both for bulk water and for isotopes
  ! and other tracers.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use decompMod               , only : bounds_type
  use clm_varpar              , only : nlevsno
  use ncdio_pio               , only : file_desc_t
  use WaterFluxBulkType       , only : waterfluxbulk_type
  use WaterStateBulkType      , only : waterstatebulk_type
  use WaterDiagnosticBulkType , only : waterdiagnosticbulk_type
  use WaterBalanceType        , only : waterbalance_type

  implicit none
  private

  !
  ! !PUBLIC TYPES:
  type, public :: water_type
     type(waterfluxbulk_type)       :: waterfluxbulk_inst
     type(waterstatebulk_type)      :: waterstatebulk_inst
     type(waterdiagnosticbulk_type) :: waterdiagnosticbulk_inst
     type(waterbalance_type)        :: waterbalancebulk_inst

   contains
     procedure, public :: Init
     procedure, public :: InitAccBuffer
     procedure, public :: InitAccVars
     procedure, public :: UpdateAccVars
     procedure, public :: Restart
  end type water_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds, &
       h2osno_col, snow_depth_col, watsat_col, t_soisno_col)
    !
    ! !DESCRIPTION:
    ! Initialize all water variables
    !
    ! !ARGUMENTS:
    class(water_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    real(r8)          , intent(in) :: h2osno_col(bounds%begc:)
    real(r8)          , intent(in) :: snow_depth_col(bounds%begc:)
    real(r8)          , intent(in) :: watsat_col(bounds%begc:, 1:)          ! volumetric soil water at saturation (porosity)
    real(r8)          , intent(in) :: t_soisno_col(bounds%begc:, -nlevsno+1:) ! col soil temperature (Kelvin)
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc

    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------

    begc = bounds%begc
    endc = bounds%endc

    SHR_ASSERT_ALL((ubound(h2osno_col) == [endc]), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(snow_depth_col) == [endc]), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(watsat_col, 1) == endc), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(t_soisno_col, 1) == endc), errMsg(sourcefile, __LINE__))

    call this%waterstatebulk_inst%InitBulk(bounds, &
         h2osno_input_col = h2osno_col(begc:endc),       &
         watsat_col = watsat_col(begc:endc, 1:),   &
         t_soisno_col = t_soisno_col(begc:endc, -nlevsno+1:) )

    call this%waterdiagnosticbulk_inst%InitBulk(bounds, &
         snow_depth_input_col = snow_depth_col(begc:endc),    &
         waterstatebulk_inst = this%waterstatebulk_inst )

    call this%waterbalancebulk_inst%Init(bounds)

    call this%waterfluxbulk_inst%InitBulk(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAccBuffer(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for all water variables
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    class(water_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds

    character(len=*), parameter :: subname = 'InitAccBuffer'
    !-----------------------------------------------------------------------

    call this%waterfluxbulk_inst%InitAccBuffer(bounds)

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize variables that are associated with accumulated fields
    !
    ! !ARGUMENTS:
    class(water_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitAccVars'
    !-----------------------------------------------------------------------

    call this%waterfluxbulk_inst%initAccVars(bounds)

  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine UpdateAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Update accumulated variables
    !
    ! Should be called every time step
    !
    ! !ARGUMENTS:
    class(water_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'UpdateAccVars'
    !-----------------------------------------------------------------------

    call this%waterfluxbulk_inst%UpdateAccVars(bounds)

  end subroutine UpdateAccVars


  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag, &
       watsat_col)
    !
    ! !DESCRIPTION:
    ! Read/write information to/from restart file for all water variables
    !
    ! !ARGUMENTS:
    class(water_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read', 'write' or 'define'
    real(r8)         , intent(in)    :: watsat_col (bounds%begc:, 1:)  ! volumetric soil water at saturation (porosity)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(watsat_col, 1) == bounds%endc), errMsg(sourcefile, __LINE__))

    call this%waterfluxbulk_inst%restartBulk (bounds, ncid, flag=flag)

    call this%waterstatebulk_inst%restartBulk (bounds, ncid, flag=flag, &
         watsat_col=watsat_col(bounds%begc:bounds%endc,:))

    call this%waterdiagnosticbulk_inst%restartBulk (bounds, ncid, flag=flag, &
         waterstatebulk_inst=this%waterstatebulk_inst)

  end subroutine Restart

end module WaterType
