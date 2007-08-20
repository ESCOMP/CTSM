module system_messages

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: system_messages
!
! !DESCRIPTION:
! Contains general purpose routines for checking system messages.
!
! !USES:
  use abortutils, only : endrun
  use clm_varctl, only : iulog
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: allocation_err ! allocation error message
  public :: netcdf_err     ! netCD return error message
  public :: iobin_err      ! binary i/o  error message
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: allocation_err
!
! !INTERFACE:
  subroutine allocation_err( ier, routine_name, array_name, nsize )
!
! !DESCRIPTION:
! Issue error message after non-zero return from an allocate statement.
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ier                   ! status from allocate statement
    character(len=*), intent(in) :: routine_name ! routine name that invoked allocate
    character(len=*), intent(in), optional :: array_name   ! array name
    integer, intent(in), optional :: nsize       ! size of attempted allocation
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------
    if ( ier /= 0 ) then
       write(iulog,*)'ERROR trying to allocate memory in routine: ' // trim(routine_name)
       write(iulog,*)'  Variable name: ' // trim(array_name)
       write(iulog,*)'  Size of allocation ', nsize
       call endrun
    end if
    return
  end subroutine allocation_err


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: netcdf_err
!
! !INTERFACE:
  subroutine netcdf_err (ier, routine_name)
!
! !DESCRIPTION:
! If error detected in netCDF call, issue error message and abort.
!
! !ARGUMENTS:
    implicit none
#include <netcdf.inc>
    integer, intent(in) :: ier                   ! return code from netCDF call
    character(len=*), intent(in) :: routine_name ! calling routine name
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------
    if ( ier /= NF_NOERR ) then
       write(iulog,*)' ERROR from netcdf library call: '
       write(iulog,*) nf_strerror( ier )
       write(iulog,*)' called from routine ', trim(routine_name)
       call endrun
    endif
    return
  end subroutine netcdf_err


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: iobin_err
!
! !INTERFACE:
  subroutine iobin_err (ier, routine_name)
!
! !DESCRIPTION:
! If error detected during i/o binary read/write
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ier              ! return code from netCDF call
    character(len=*), intent(in) :: routine_name ! calling routine name
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------
    if (ier /= 0) then
       write(iulog,*)' i/o binary error ',ier
       write(iulog,*)' called from routine ', trim(routine_name)
       call endrun
    end if

  end subroutine iobin_err

end module system_messages
