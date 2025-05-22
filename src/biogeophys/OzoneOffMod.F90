module OzoneOffMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Provides an implementatio of ozone_base_type for the ozone-off case. Note that very
  ! little needs to be done in this case, so this module mainly provides empty
  ! implementations to satisfy the interface.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod, only : r8 => shr_kind_r8
  use decompMod   , only : bounds_type
  use OzoneBaseMod, only : ozone_base_type

  implicit none
  save
  private

  ! !PUBLIC TYPES:
  type, extends(ozone_base_type), public :: ozone_off_type
     private
   contains
     procedure, public :: Init
     procedure, public :: Restart
     procedure, public :: CalcOzoneUptake
     procedure, public :: CalcOzoneStress
  end type ozone_off_type

  interface ozone_off_type
     module procedure constructor
  end interface ozone_off_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  function constructor() result(ozone_off)
    !
    ! !DESCRIPTION:
    ! Return an instance of ozone_off_type
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ozone_off_type) :: ozone_off  ! function result
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    ! DO NOTHING (simply return a variable of the appropriate type)

    ! Eventually this should call the Init routine (or replace the Init routine
    ! entirely). But I think it would be confusing to do that until we switch everything
    ! to use a constructor rather than the init routine.

  end function constructor


  subroutine Init(this,  bounds, o3_veg_stress_method)
    !
    ! !DESCRIPTION:
    ! Initialize ozone data structure
    !
    ! ! USES:
    use abortutils,only : endrun
    !
    ! !ARGUMENTS:
    class(ozone_off_type) , intent(inout) :: this
    type(bounds_type)     , intent(in)           :: bounds
    character(len=*), intent(in)                      :: o3_veg_stress_method
    !-----------------------------------------------------------------------

    if (o3_veg_stress_method /= 'unset' ) call endrun(' unconsistent choice of o3_veg_stress_method in init OzoneOffMod.')

    call this%InitAllocateBase(bounds)
    call this%InitColdBase(bounds)

  end subroutine Init

  subroutine Restart(this, bounds, ncid, flag)
    use ncdio_pio , only : file_desc_t

    class(ozone_off_type) :: this
    type(bounds_type), intent(in) :: bounds
    type(file_desc_t) , intent(inout) :: ncid ! netcdf id
    character(len=*)  , intent(in)    :: flag ! 'read', 'write' or 'define'

    ! DO NOTHING

  end subroutine Restart

  subroutine CalcOzoneUptake(this, bounds, num_exposedvegp, filter_exposedvegp, &
       forc_pbot, forc_th, rssun, rssha, rb, ram, tlai, forc_o3)

    class(ozone_off_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    integer  , intent(in) :: num_exposedvegp           ! number of points in filter_exposedvegp
    integer  , intent(in) :: filter_exposedvegp(:)     ! patch filter for non-snow-covered veg
    real(r8) , intent(in) :: forc_pbot( bounds%begc: ) ! atmospheric pressure (Pa)
    real(r8) , intent(in) :: forc_th( bounds%begc: )   ! atmospheric potential temperature (K)
    real(r8) , intent(in) :: rssun( bounds%begp: )     ! leaf stomatal resistance, sunlit leaves (s/m)
    real(r8) , intent(in) :: rssha( bounds%begp: )     ! leaf stomatal resistance, shaded leaves (s/m)
    real(r8) , intent(in) :: rb( bounds%begp: )        ! boundary layer resistance (s/m)
    real(r8) , intent(in) :: ram( bounds%begp: )       ! aerodynamical resistance (s/m)
    real(r8) , intent(in) :: tlai( bounds%begp: )      ! one-sided leaf area index, no burying by snow
    real(r8) , intent(in) :: forc_o3( bounds%begg: )   ! atmospheric ozone (mol/mol)

    ! Enforce expected array sizes (mainly so that a debug-mode threaded test with
    ! ozone-off can pick up problems with the call to this routine)
    SHR_ASSERT_ALL_FL((ubound(forc_pbot) == (/bounds%endc/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(forc_th) == (/bounds%endc/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(rssun) == (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(rssha) == (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(rb) == (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(ram) == (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(tlai) == (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(forc_o3) == (/bounds%endg/)), sourcefile, __LINE__)

    ! Do nothing: In the ozone off case, we don't need to track ozone uptake

  end subroutine CalcOzoneUptake

  subroutine CalcOzoneStress(this, bounds, &
       num_exposedvegp, filter_exposedvegp, &
       num_noexposedvegp, filter_noexposedvegp)
    class(ozone_off_type), intent(inout) :: this
    type(bounds_type)    , intent(in) :: bounds
    integer              , intent(in) :: num_exposedvegp
    integer              , intent(in) :: filter_exposedvegp(:)
    integer              , intent(in) :: num_noexposedvegp
    integer              , intent(in) :: filter_noexposedvegp(:)

    ! Outputs (stress terms) are already fixed at 1 from cold start initialization

  end subroutine CalcOzoneStress

end module OzoneOffMod
