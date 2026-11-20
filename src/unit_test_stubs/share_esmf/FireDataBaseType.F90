module FireDataBaseType

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! module for handling of fire data
  ! UNIT-TEST STUB for fire data Streams
  !       This just allows the fire code to be tested without
  !       reading in the streams data, by faking it and setting it to a
  !       constant value.
  !
  ! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use clm_varctl       , only : iulog
  use spmdMod          , only : masterproc, mpicom, iam
  use abortutils       , only : endrun
  use decompMod        , only : bounds_type
  use FireMethodType   , only : fire_method_type
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: fire_base_type
  !
  type, abstract, extends(fire_method_type) :: fire_base_type
    private
      ! !PRIVATE MEMBER DATA:
      real(r8), public, pointer :: forc_hdm(:)  ! Human population density
      real(r8), public, pointer :: forc_lnfm(:) ! Lightning frequency
      real(r8), public, pointer :: gdp_lf_col(:)   ! col global real gdp data (k US$/capita)
      real(r8), public, pointer :: peatf_lf_col(:) ! col global peatland fraction data (0-1)
      integer , public, pointer :: abm_lf_col(:)   ! col global peak month of crop fire emissions
      
    contains
      !
      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public :: BaseFireInit             ! Initialization of Fire
      procedure, public :: FireInit => BaseFireInit ! Initialization of Fire
      procedure, public :: FireInterp               ! Interpolate fire data
      procedure, public :: BaseFireReadNML          ! Read in the namelist
      procedure, public :: FireReadNML => BaseFireReadNML  ! Read in the namelist
      procedure(need_lightning_and_popdens_interface), public, deferred :: &
           need_lightning_and_popdens               ! Returns true if need lightning & popdens

  end type fire_base_type

  abstract interface
     !-----------------------------------------------------------------------
     function need_lightning_and_popdens_interface(this) result(need_lightning_and_popdens)
       !
       ! !DESCRIPTION:
       ! Returns true if need lightning and popdens, false otherwise
       !
       ! USES
       import :: fire_base_type
       !
       ! !ARGUMENTS:
       class(fire_base_type), intent(in) :: this
       logical :: need_lightning_and_popdens  ! function result
       !-----------------------------------------------------------------------
     end function need_lightning_and_popdens_interface
  end interface

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine BaseFireReadNML( this, bounds, NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for Fire
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(fire_base_type) :: this
    type(bounds_type), intent(in) :: bounds
    character(len=*), intent(in) :: NLFilename ! Namelist filename
  end subroutine BaseFireReadNML

  !================================================================
  subroutine BaseFireInit( this, bounds )
    !
    ! !DESCRIPTION:
    ! Initialize CN Fire module
    ! !USES:
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(fire_base_type) :: this
    type(bounds_type), intent(in) :: bounds
    !-----------------------------------------------------------------------

    if ( this%need_lightning_and_popdens() ) then

    end if

  end subroutine BaseFireInit

  !================================================================
  subroutine FireInterp(this,bounds)
    !
    ! !DESCRIPTION:
    ! Interpolate CN Fire datasets
    !
    ! !ARGUMENTS:
    class(fire_base_type) :: this
    type(bounds_type), intent(in) :: bounds
    !-----------------------------------------------------------------------

    if ( this%need_lightning_and_popdens() ) then

    end if

  end subroutine FireInterp

end module FireDataBaseType
