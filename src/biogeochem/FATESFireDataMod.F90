module FATESFireDataMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! module for FATES to obtain fire inputs from data
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_log_mod, only: errmsg => shr_log_errMsg
  use abortutils, only: endrun
  use clm_varctl, only: iulog
  use decompMod, only: bounds_type
  use CNFireMethodMod, only: cnfire_method_type
  use CNFireBaseMod, only: cnfire_base_type
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: fates_fire_data_type
  !
  type, extends(cnfire_base_type) :: fates_fire_data_type
      ! !PRIVATE MEMBER DATA:
      real(r8), public, pointer :: lnfm24(:)  ! Daily avg lightning by grid cell (#/km2/hr)
    contains
      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public :: InitAccBuffer  ! Initialize accumulation processes
      procedure, public :: InitAccVars  ! Initialize accumulation variables
      procedure, public :: UpdateAccVars  ! Update/extract accumulations vars

  end type fates_fire_data_type

  character(len=*), parameter, private :: sourcefile = __FILE__
  !
  ! !PRIVATE MEMBER DATA:
  !-----------------------------------------------------------------------

  interface fates_fire_data_type
     ! initialize a new cnfire_base object
     module procedure constructor
  end interface fates_fire_data_type
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  type(fates_fire_data_type) function constructor()
    !
    ! !DESCRIPTION:
    ! Creates an object of type fates_fire_data_type
    ! !ARGUMENTS:
    constructor%need_lightning_and_popdens = .true.
  end function constructor

  !-----------------------------------------------------------------------
  subroutine InitAccBuffer (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for all required module accumulated fields
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    !
    ! !USES
    use clm_varcon, only : spval
    use accumulMod, only : init_accum_field
    !
    ! !ARGUMENTS:
    class(fates_fire_data_type) :: this
    type(bounds_type), intent(in) :: bounds

    ! !LOCAL VARIABLES:
    integer  :: begg, endg
    integer  :: ier
    !---------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    allocate(this%lnfm24(begg:endg), stat=ier)
    if (ier/=0) then
       call endrun(msg="allocation error for lnfm24"//&
            errmsg(sourcefile, __LINE__))
    endif
    this%lnfm24(:) = spval
    call init_accum_field (name='lnfm24', units='strikes/km2/hr', &
         desc='24hr average of lightning strikes',  accum_type='runmean', &
         accum_period=-1, subgrid_type='gridcell', numlev=1, init_value=0._r8)

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart
    ! file is read in and the accumulation buffer is obtained)
    !
    ! !USES
    use accumulMod       , only : extract_accum_field
    use clm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(fates_fire_data_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begg, endg
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslg(:)  ! temporary
    !---------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    ! Allocate needed dynamic memory for single level patch field
    allocate(rbufslg(begg:endg), stat=ier)
    if (ier/=0) then
       call endrun(msg="allocation error for rbufslg"//&
            errmsg(sourcefile, __LINE__))
    endif

    ! Determine time step
    nstep = get_nstep()

    call extract_accum_field ('lnfm24', rbufslg, nstep)
    this%lnfm24(begg:endg) = rbufslg(begg:endg)

    deallocate(rbufslg)

  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine UpdateAccVars (this, bounds)
    !
    ! USES
    use clm_time_manager, only: get_nstep
    use accumulMod, only: update_accum_field, extract_accum_field
    use abortutils, only: endrun
    !
    ! !ARGUMENTS:
    class(fates_fire_data_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: dtime                 ! timestep size [seconds]
    integer :: nstep                 ! timestep number
    integer :: ier                   ! error status
    integer :: begg, endg
    real(r8), pointer :: rbufslg(:)  ! temporary single level - gridcell level
    !---------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    nstep = get_nstep()

    ! Allocate needed dynamic memory for single level gridcell field

    allocate(rbufslg(begg:endg), stat=ier)
    if (ier/=0) then
       write(iulog,*)'UpdateAccVars allocation error for rbuf1dg'
       call endrun(msg=errmsg(sourcefile, __LINE__))
    endif

    ! Accumulate and extract lnfm24
    rbufslg(begg:endg) = this%forc_lnfm(begg:endg)
    call update_accum_field  ('lnfm24', rbufslg, nstep)
    call extract_accum_field ('lnfm24', this%lnfm24, nstep)

    deallocate(rbufslg)

  end subroutine UpdateAccVars

end module FATESFireDataMod
