module FATESFireNoDataMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! module for FATES when not obtaining fire inputs from data
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_log_mod, only: errmsg => shr_log_errMsg
  use abortutils, only: endrun
  use clm_varctl, only: iulog
  use decompMod, only: bounds_type
  use FatesFireBase, only: fates_fire_base_type
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: fates_fire_no_data_type
  !
  type, extends(fates_fire_base_type) :: fates_fire_no_data_type
      private

      ! !PRIVATE MEMBER DATA:

    contains
      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public :: FATESNoFireInit! Initialization
      procedure, public :: FireInit => FATESNoFireInit
      procedure, public :: need_lightning_and_popdens
      procedure, public :: GetLight24     ! Return the 24-hour averaged lightning data
      procedure, public :: GetGDP         ! Return the global gdp data
      procedure, public :: InitAccBuffer  ! Initialize accumulation processes
      procedure, public :: InitAccVars    ! Initialize accumulation variables
      procedure, public :: UpdateAccVars  ! Update/extract accumulations vars

  end type fates_fire_no_data_type

  character(len=*), parameter, private :: sourcefile = __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine FATESNoFireInit( this, bounds )
    !
    ! !DESCRIPTION:
    ! Initialize No Fire data module for FATES
    use shr_fire_emis_mod, only : shr_fire_emis_mechcomps_n
    use shr_log_mod      , only : errMsg => shr_log_errMsg
    use clm_varctl       , only : fates_spitfire_mode
    ! !ARGUMENTS:
    class(fates_fire_no_data_type) :: this
    type(bounds_type), intent(in) :: bounds

    if ( (shr_fire_emis_mechcomps_n > 0) .and. (fates_spitfire_mode == 0) ) then
      write(iulog,*) "Fire emissions can NOT be active for fates_spitfire_mode=0 (no_fire)",  &
                  errMsg(sourcefile, __LINE__)
      call endrun(msg="Having fire emissions on requires fates_spitfire_mode to be something besides no_fire (0)" )
      return
    end if
    call this%CNFireInit( bounds )

  end subroutine FATESNoFireInit

  !------------------------------------------------------------------------
  function need_lightning_and_popdens(this)
    ! !ARGUMENTS:
    class(fates_fire_no_data_type), intent(in) :: this
    logical :: need_lightning_and_popdens  ! function result
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'need_lightning_and_popdens'
    !-----------------------------------------------------------------------

    need_lightning_and_popdens = .false.
  end function need_lightning_and_popdens

  !-----------------------------------------------------------------------
  function GetLight24( this ) result(lnfm24)
    !
    ! !DESCRIPTION: Get the 24-hour averaged lightning data
    ! !USES
    !
    ! !ARGUMENTS:
    class(fates_fire_no_data_type) :: this
    real(r8), pointer :: lnfm24(:)
    !---------------------------------------------------------------------
    call endrun( "GetLight24 should NOT be called for the FATES No-Data case" )
    !---------------------------------------------------------------------
  end function
  
  !-----------------------------------------------------------------------
  function GetGDP( this ) result(gdp)
    !
    ! !DESCRIPTION: Get the global gross domestic product data
    ! !USES
    !
    ! !ARGUMENTS:
    class(fates_fire_no_data_type) :: this
    real(r8), pointer :: gdp(:)
    !---------------------------------------------------------------------
    call endrun( "GetGDP should NOT be called for the FATES No-Data case" )
    !---------------------------------------------------------------------
  end function

  !-----------------------------------------------------------------------
  subroutine InitAccBuffer (this, bounds)
    !
    ! !DESCRIPTION:
    ! EMPTY subroutine for the no_data case.
    ! Initialize accumulation buffer for all required module accumulated fields
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    !
    ! !USES
    !
    ! !ARGUMENTS:
    class(fates_fire_no_data_type) :: this
    type(bounds_type), intent(in) :: bounds

    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------


  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! EMPTY subroutine for the no_data case.
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart
    ! file is read in and the accumulation buffer is obtained)
    !
    ! !USES
    !
    ! !ARGUMENTS:
    class(fates_fire_no_data_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------


  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine UpdateAccVars (this, bounds)
    !
    ! !DESCRIPTION:
    ! EMPTY subroutine for the no_data case.
    !
    ! !USES
    !
    ! !ARGUMENTS:
    class(fates_fire_no_data_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------


  end subroutine UpdateAccVars

end module FATESFireNoDataMod
