module ExcessIceStreamType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Stub module for Excess ice streams
  !
  ! !USES
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use abortutils       , only : endrun
  use decompMod        , only : bounds_type

  ! !PUBLIC TYPES:
  implicit none
  private

  public  :: UseExcessIceStreams      ! If streams will be used

  type, public :: excessicestream_type
  contains

      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public  :: Init            ! Initialize and read data in
      procedure, public  :: CalcExcessIce   ! Calculate excess ice ammount
  end type excessicestream_type
    ! ! PRIVATE DATA:

  logical :: namelist_read = .false.

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine Init(this, bounds, NLFilename)
    !
    ! arguments
    implicit none
    class(excessicestream_type) :: this
    type(bounds_type), intent(in)   :: bounds
    character(len=*),  intent(in)   :: NLFilename   ! Namelist filename

    namelist_read = .true.

  end subroutine Init

  subroutine CalcExcessIce(this,bounds,exice_bulk_init)
  
  ! only transfers grid values to columns
   implicit none
   class(excessicestream_type)        :: this
   type(bounds_type),  intent(in)     :: bounds
   real(r8)         ,  intent(inout)  :: exice_bulk_init(bounds%begc:bounds%endc) 
   !
   ! !LOCAL VARIABLES:
   call endrun(msg=' ERROR CalcExcessIce stub is being called and should NOT be'//errMsg(sourcefile, __LINE__))

  end subroutine CalcExcessIce

  logical function UseExcessIceStreams()
  !
  ! !DESCRIPTION:
  ! Return true if
  !
  ! !USES:
  !
  ! !ARGUMENTS:
  implicit none
  !
  ! !LOCAL VARIABLES:
  if ( .not. namelist_read ) then
      call endrun(msg=' ERROR UseExcessIceStreams being called, but namelist has not been read yet'//errMsg(sourcefile, __LINE__))
  end if
  UseExcessIceStreams = .false.
end function UseExcessIceStreams

end module ExcessIceStreamType
