module ZenderSoilErodStreamType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! UNIT-TEST STUB for ZenderSoilErodStreamType
  !
  ! !USES
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use abortutils       , only : endrun
  use decompMod        , only : bounds_type

  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: soil_erod_stream_type
  contains

      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public :: Init            ! Initialize and read data in
      procedure, public :: CalcDustSource ! Calculate dust source spatial filter (basically truncating stream data value smaller than 0.1 following CAM's practice) based on input streams
      procedure, public :: UseStreams      ! If streams will be used

  end type soil_erod_stream_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine Init(this, bounds, NLFilename)
   !
   ! Initialize the Zender soil eroditability stream object
   !
   ! Uses:
   !
   ! arguments
   implicit none
   class(soil_erod_stream_type) :: this
   type(bounds_type), intent(in)   :: bounds
   character(len=*),  intent(in)   :: NLFilename   ! Namelist filename
   !
   ! local variables
   !-----------------------------------------------------------------------

  end subroutine Init

  !==============================================================================
  logical function UseStreams(this)
    !
    ! !DESCRIPTION:
    ! Return true if the Zender method is being used and the soil erodability
    ! file is being used with it
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    class(soil_erod_stream_type) :: this
    !
    ! !LOCAL VARIABLES:
    UseStreams = .false.
  end function UseStreams

  !==============================================================================
  subroutine CalcDustSource(this, bounds, soil_erod)
    !
    ! !DESCRIPTION:
    !  Calculate the soil eroditability for the Zender dust method.
    !
    ! !USES:
    !USES
    !
    ! !ARGUMENTS:
    implicit none
    class(soil_erod_stream_type)             :: this
    type(bounds_type)              , intent(in)    :: bounds
    real(r8)                       , intent(inout) :: soil_erod(bounds%begc:)      ! [fraction] rock drag partition factor (roughness effect)
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

  end subroutine CalcDustSource

end module ZenderSoilErodStreamType
