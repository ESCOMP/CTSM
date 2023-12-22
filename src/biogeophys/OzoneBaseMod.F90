module OzoneBaseMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Define the interface for ozone_type, which calculates ozone-induced stress. The type
  ! defined here is abstract; it will get instantiated as a concrete type that extends
  ! this base type (e.g., an ozone-off or ozone-on version).
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use decompMod   , only : bounds_type

  implicit none
  save
  private

  ! !PUBLIC TYPES:
  type, abstract, public :: ozone_base_type
     private

     ! Public data members
     ! These should be treated as read-only by other modules (except that they can be
     ! modified by extensions of the ozone_base_type)
     real(r8), pointer, public :: o3coefvsha_patch(:)          ! ozone coefficient for photosynthesis, shaded leaves (0 - 1)
     real(r8), pointer, public :: o3coefvsun_patch(:)         ! ozone coefficient for photosynthesis, sunlit leaves (0 - 1)
     real(r8), pointer, public :: o3coefgsha_patch(:)         ! ozone coefficient for conductance, shaded leaves (0 - 1)
     real(r8), pointer, public :: o3coefgsun_patch(:)         ! ozone coefficient for conductance, sunlit leaves (0 - 1)
     real(r8), pointer, public :: o3coefjmaxsha_patch(:)  ! ozone coefficient for max electron transport rate, shaded leaves (0 - 1)
     real(r8), pointer, public :: o3coefjmaxsun_patch(:)  ! ozone coefficient for max electron transport rate, sunlit leaves (0 - 1)
   contains
     ! The following routines need to be implemented by all type extensions
     procedure(Init_interface)            , public, deferred :: Init
     procedure(Restart_interface)         , public, deferred :: Restart
     procedure(CalcOzoneUptake_interface) , public, deferred :: CalcOzoneUptake
     procedure(CalcOzoneStress_interface) , public, deferred :: CalcOzoneStress

     ! The following routines should only be called by extensions of the ozone_base_type
     procedure, public :: InitAllocateBase
     procedure, public :: InitColdBase

  end type ozone_base_type

  abstract interface

     subroutine Init_interface(this, bounds, o3_veg_stress_method)
       use decompMod, only : bounds_type
       import :: ozone_base_type

       class(ozone_base_type), intent(inout) :: this
       type(bounds_type), intent(in) :: bounds
       character(len=*),        intent(in) :: o3_veg_stress_method

     end subroutine Init_interface

     subroutine Restart_interface(this, bounds, ncid, flag)
       use decompMod , only : bounds_type
       use ncdio_pio , only : file_desc_t
       import :: ozone_base_type

       class(ozone_base_type)            :: this
       type(bounds_type) , intent(in)    :: bounds
       type(file_desc_t) , intent(inout) :: ncid ! netcdf id
       character(len=*)  , intent(in)    :: flag ! 'read', 'write' or 'define'
     end subroutine Restart_interface

     subroutine CalcOzoneUptake_interface(this, bounds, num_exposedvegp, filter_exposedvegp, &
          forc_pbot, forc_th, rssun, rssha, rb, ram, tlai, forc_o3)
       use decompMod    , only : bounds_type
       use shr_kind_mod , only : r8 => shr_kind_r8
       import :: ozone_base_type

       class(ozone_base_type) , intent(inout) :: this
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
     end subroutine CalcOzoneUptake_interface

     subroutine CalcOzoneStress_interface(this, bounds, &
          num_exposedvegp, filter_exposedvegp, &
          num_noexposedvegp, filter_noexposedvegp)
       use decompMod, only : bounds_type
       import :: ozone_base_type

       class(ozone_base_type) , intent(inout) :: this
       type(bounds_type)      , intent(in)    :: bounds
       integer                , intent(in)    :: num_exposedvegp         ! number of points in filter_exposedvegp
       integer                , intent(in)    :: filter_exposedvegp(:)   ! patch filter for non-snow-covered veg
       integer                , intent(in)    :: num_noexposedvegp       ! number of points in filter_noexposedvegp
       integer                , intent(in)    :: filter_noexposedvegp(:) ! patch filter for veg where frac_veg_nosno is 0
     end subroutine CalcOzoneStress_interface
  end interface

contains

  !-----------------------------------------------------------------------
  subroutine InitAllocateBase(this, bounds)
    !
    ! !DESCRIPTION:
    ! Allocate variables in the base class
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(ozone_base_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp

    character(len=*), parameter :: subname = 'InitAllocateBase'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    allocate(this%o3coefvsha_patch(begp:endp))  ; this%o3coefvsha_patch(:)                = nan
    allocate(this%o3coefvsun_patch(begp:endp))  ; this%o3coefvsun_patch(:)                = nan
    allocate(this%o3coefgsha_patch(begp:endp))  ; this%o3coefgsha_patch(:)                = nan
    allocate(this%o3coefgsun_patch(begp:endp))  ; this%o3coefgsun_patch(:)               = nan
    allocate(this%o3coefjmaxsha_patch(begp:endp))  ; this%o3coefjmaxsha_patch(:) = nan
    allocate(this%o3coefjmaxsun_patch(begp:endp))  ; this%o3coefjmaxsun_patch(:) = nan


  end subroutine InitAllocateBase


  !-----------------------------------------------------------------------
  subroutine InitColdBase(this, bounds)
    !
    ! !DESCRIPTION:
    ! Do cold start initialization for variables in the base class. Note that this
    ! initialization will be the same for all ozone implementations, including the
    ! ozone-off implementation.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(ozone_base_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp

    character(len=*), parameter :: subname = 'InitColdBase'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    this%o3coefvsha_patch(begp:endp)        = 1._r8
    this%o3coefvsun_patch(begp:endp)        = 1._r8
    this%o3coefgsha_patch(begp:endp)        = 1._r8
    this%o3coefgsun_patch(begp:endp)        = 1._r8
    this%o3coefjmaxsha_patch(begp:endp) = 1._r8
    this%o3coefjmaxsun_patch(begp:endp) = 1._r8

  end subroutine InitColdBase

end module OzoneBaseMod
