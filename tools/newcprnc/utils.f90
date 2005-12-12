module utils
!
! Utility routines
!
   implicit none

   private
   save

   public :: endrun
   public :: lenchr

CONTAINS

   subroutine endrun (msg)
!-------------------------------------------------------------------------------------------
! Purpose: Print an optional string and abort
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in), optional :: msg    ! string to be printed

      if (present (msg)) then
         write(6,*)'ENDRUN:', msg
      else
         write(6,*)'ENDRUN: called without a message string'
      end if
      
      call abort ()
   end subroutine endrun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function lenchr (string)
!-------------------------------------------------------------------------------------------
! Purpose: Determine the position of the last non-null, non-blank character in the string.
!          This function is needed (e.g. instead of len_trim or trim) because various
!          netcdf routines null-fill output strings.
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in) :: string       !  Input character string
!
! Local workspace
!
      integer :: i
!
!-----------------------------------------------------------------------
!
      lenchr = 0
      do i=len(string),1,-1
         if (string(i:i) /= ' ' .and. string(i:i) /= char(0)) then
            lenchr = i
            exit
         end if
      end do

      return
   end function lenchr
end module utils
