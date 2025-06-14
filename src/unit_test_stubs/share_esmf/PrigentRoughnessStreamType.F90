module PrigentRoughnessStreamType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! UNIT-TEST STUB for Prigent Roughtness Stream
  !       This just allows the Leung dust emission code to be tested without
  !       reading in the streams data, by faking it and setting it to a
  !       constant value.
  ! NOTE: THIS SHOULD BE REMOVED WHEN THE MAIN VERSION ALLOWS IT TO BE OFF.
  !       https://github.com/ESCOMP/CTSM/issues/2381
  ! Removing this version will remove testing code duplication.
  ! !USES
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use clm_varctl       , only : iulog
  use abortutils       , only : endrun
  use decompMod        , only : bounds_type

  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: prigent_roughness_stream_type
    real(r8), pointer, public :: prigent_rghn  (:)         ! Prigent et al. (1997) roughness length (m)
  contains

      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public :: Init            ! Initialize and read data in
      procedure, public :: UseStreams      ! If Prigent rougness streams will be used
      procedure, public :: IsStreamInit    ! If the streams have been initialized and read in, so data can be used
      procedure, public :: Clean           ! Clean and deallocate

  end type prigent_roughness_stream_type

  logical, private :: InitDone = .false.  ! If initialization of streams has been done

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine Init(this, bounds, NLFilename)
   !
   ! Initialize the prigent roughness stream object
   !
   ! Uses:
   !
   ! arguments
   implicit none
   class(prigent_roughness_stream_type) :: this
   type(bounds_type), intent(in)   :: bounds
   character(len=*),  intent(in)   :: NLFilename   ! Namelist filename
   !
   ! local variables
   !-----------------------------------------------------------------------
   allocate(this%prigent_rghn(bounds%begg:bounds%endg))
   this%prigent_rghn(:) = 1.0_r8     ! This should likely be a smaller value like 0.1 based on the Prigent dataset average values
   InitDone = .true.

  end subroutine Init

  !==============================================================================
  logical function UseStreams(this)
    !
    ! !DESCRIPTION:
    ! Return true if the Prigent Roughness stream is being used
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    class(prigent_roughness_stream_type) :: this
    !
    !
    if ( .not. this%IsStreamInit() )then
       UseStreams = .false.
       call endrun( "Not initialized" )
    else
       UseStreams = .true.
    end if
  end function UseStreams

  !==============================================================================
  logical function IsStreamInit(this)
    !
    ! !DESCRIPTION:
    ! Return true if the streams have been initialized
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    class(prigent_roughness_stream_type) :: this
    !
    character(len=*), parameter :: subname = 'PrigentRoughnessStream::IsStreamInit'
    !
    if ( InitDone )then
      IsStreamInit = .true.
    else
      IsStreamInit = .false.
    end if
  end function IsStreamInit

  !==============================================================================
  subroutine Clean(this)
    !
    ! Deallocate and clean the object
    !
    ! Uses:
    !
    ! arguments
    implicit none
    class(prigent_roughness_stream_type) :: this
    !
    ! local variables
    !-----------------------------------------------------------------------
    deallocate(this%prigent_rghn)
    this%prigent_rghn => NULL()
    InitDone = .false.
 
   end subroutine Clean

end module PrigentRoughnessStreamType
