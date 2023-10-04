module ExcessIceStreamType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Stub for ExcessIceStreams for the MCT driver. So that MCT can be used
  ! without excess ice streams.
  !
  ! !USES
  use shr_kind_mod     , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use spmdMod          , only : mpicom, masterproc
  use clm_varctl       , only : iulog
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

      ! !PRIVATE MEMBER FUNCTIONS:
      procedure, private :: ReadNML     ! Read in namelist

  end type excessicestream_type
  ! ! PRIVATE DATA:

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine Init(this, bounds, NLFilename)
    !
    !
    ! arguments
    implicit none
    class(excessicestream_type) :: this
    type(bounds_type), intent(in)   :: bounds
    character(len=*),  intent(in)   :: NLFilename   ! Namelist filename

    !
    ! local variables

    call this%ReadNML( bounds, NLFileName )
  end subroutine Init

  subroutine CalcExcessIce(this,bounds,exice_bulk_init)
  
  ! only transfers grid values to columns
   implicit none
   class(excessicestream_type)        :: this
   type(bounds_type),  intent(in)     :: bounds
   real(r8)         ,  intent(inout)  :: exice_bulk_init(bounds%begc:bounds%endc) 
   !
   ! !LOCAL VARIABLES:
   
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
  UseExcessIceStreams = .false.
end function UseExcessIceStreams

subroutine ReadNML(this, bounds, NLFilename)
  !
  ! Read the namelist data stream information.
  !
  ! Uses:
  use shr_nl_mod       , only : shr_nl_find_group_name
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use shr_mpi_mod      , only : shr_mpi_bcast
  !
  ! arguments
  implicit none
  class(excessicestream_type)   :: this
  type(bounds_type), intent(in) :: bounds
  character(len=*),  intent(in) :: NLFilename   ! Namelist filename
  !
  ! local variables
  integer            :: nu_nml    ! unit for namelist file
  integer            :: nml_error ! namelist i/o error flag
  logical            :: use_excess_ice_streams = .false.         ! logical to turn on use of excess ice streams
  character(len=CL)  :: stream_fldFileName_exice = ' '
  character(len=CL)  :: stream_mapalgo_exice = 'none'
  character(len=*), parameter :: namelist_name = 'exice_streams'    ! MUST agree with name in namelist and read
  character(len=*), parameter :: subName = "('exice_streams::ReadNML')"
  !-----------------------------------------------------------------------

  namelist /exice_streams/ &               ! MUST agree with namelist_name above
       stream_mapalgo_exice,  stream_fldFileName_exice, use_excess_ice_streams
  !-----------------------------------------------------------------------
  ! Default values for namelist

  ! Read excess ice namelist
  if (masterproc) then
     open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
     call shr_nl_find_group_name(nu_nml, namelist_name, status=nml_error)
     if (nml_error == 0) then
        read(nu_nml, nml=exice_streams,iostat=nml_error)   ! MUST agree with namelist_name above
        if (nml_error /= 0) then
           call endrun(msg=' ERROR reading '//namelist_name//' namelist'//errMsg(sourcefile, __LINE__))
        end if
     else
        call endrun(msg=' ERROR finding '//namelist_name//' namelist'//errMsg(sourcefile, __LINE__))
     end if
     close(nu_nml)
  endif

  call shr_mpi_bcast(use_excess_ice_streams   , mpicom)

  if (masterproc) then
     if ( use_excess_ice_streams ) then
        call endrun(msg=' ERROR excess ice streams can NOT be on for the MCT driver'//errMsg(sourcefile, __LINE__))
     end if
     if ( trim(stream_fldFileName_exice) /= ''  ) then
        call endrun(msg=' ERROR stream_fldFileName_exice can NOT be set for the MCT driver'//errMsg(sourcefile, __LINE__))
     end if
     if ( trim(stream_mapalgo_exice) /= 'none'  ) then
        call endrun(msg=' ERROR stream_mapalgo_exice can only be none for the MCT driver'//errMsg(sourcefile, __LINE__))
     end if
  endif

end subroutine ReadNML

end module ExcessIceStreamType
