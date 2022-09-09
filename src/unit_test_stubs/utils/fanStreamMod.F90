module FanStreamMod

  !----------------------------------------------------------------------- 
  ! Contains methods for reading in FAN nitrogen deposition (in the form of
  ! manure) data file
  ! Also includes functions for fan stream file handling and 
  ! interpolation.
  !
  ! uses:
  use shr_kind_mod     , only : r8 => shr_kind_r8, CL => shr_kind_cl
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use abortutils       , only : endrun
  use decompMod        , only : bounds_type
  !
  implicit none
  private

  ! Public interfaces
  public :: fanstream_init            ! Initialize FAN streams data
  public :: fanstream_interp          ! interpolates between two years of FAN file data
  public :: set_bcast_fanstream_pars  ! Set teh namelist parameters for the FAN streams

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !==============================================================================

contains

  !==============================================================================

  subroutine set_bcast_fanstream_pars(str_yr_first, str_yr_last, mdl_yr_align, mapalgo, str_filename, crop_man_is_percrop)
    !-----------------------------------------------------------------------
    !    
    ! Set the FAN stream namelist parameters
    !
    !-----------------------------------------------------------------------
    ! uses:
    ! Arguments:
    integer, intent(in) :: str_yr_first, str_yr_last, mdl_yr_align
    ! whether manure_sgrz and manure_ngrz are per crop or land area:
    logical, intent(in) :: crop_man_is_percrop 
    character(len=*), intent(in) :: str_filename, mapalgo
    ! Local variables:
    character(*), parameter :: subName = "('set_bcast_fanstream_pars')"
    !-----------------------------------------------------------------------
    call endrun(msg=subName//'ERROR this is a stub and should not be called'//errMsg(sourcefile, __LINE__))

  end subroutine set_bcast_fanstream_pars

  !************************************************************************************
  
  subroutine fanstream_init(bounds, NLFilename)
   !-----------------------------------------------------------------------
   !    
   ! Initialize data stream information for FAN
   !
   !-----------------------------------------------------------------------
   ! uses:
   !
   ! Arguments:
   implicit none
   type(bounds_type), intent(in) :: bounds  
   character(len=*),  intent(in) :: NLFilename   ! Namelist filename
   !
   ! Local variables:
   character(*), parameter :: subName = "('fanstream_init')"

   call endrun(msg=subName//'ERROR this is a stub and should not be called'//errMsg(sourcefile, __LINE__))

 end subroutine fanstream_init
  
 !================================================================

 subroutine fanstream_interp(bounds, atm2lnd_inst)
   !-----------------------------------------------------------------------
   !
   ! Interpoalte the FAN data to the current simulation time
   !
   !-----------------------------------------------------------------------
   use atm2lndType     , only : atm2lnd_type
   !
   ! Arguments
   type(bounds_type) , intent(in)    :: bounds  
   type(atm2lnd_type), intent(inout) :: atm2lnd_inst
   !
   ! Local variables
   character(*), parameter :: subName = "('fanstream_interp')"
   !-----------------------------------------------------------------------

   call endrun(msg=subName//'ERROR this is a stub and should not be called'//errMsg(sourcefile, __LINE__))

 end subroutine fanstream_interp
    
end module FanStreamMod
