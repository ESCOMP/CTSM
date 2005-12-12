#include <misc.h>
#include <preproc.h>

module fileutils

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fileutils
!
! !DESCRIPTION:
! Module containing file I/O utilities
!
! !USES:
  use abortutils, only : endrun
!
! !PUBLIC TYPES:
  implicit none
  save
  logical, public :: lsmiou(99)  !I/O file unit numbers (1 to 99)
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: get_filename  !Returns filename given full pathname
  public :: set_filename  !Set remote full path filename
  public :: opnfil        !Open local unformatted or formatted file
  public :: getfil        !Obtain local copy of file
  public :: putfil        !Dispose file to Mass Store
  public :: relavu        !Close and release Fortran unit no longer in use
  public :: getavu        !Get next available Fortran unit number
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
  private:: shell_cmd
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_filename
!
! !INTERFACE:
  character(len=256) function get_filename (fulpath)
!
! !DESCRIPTION:
! Returns filename given full pathname
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: fulpath !full pathname
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer i               !loop index
    integer klen            !length of fulpath character string
!------------------------------------------------------------------------

    klen = len_trim(fulpath)
    do i = klen, 1, -1
       if (fulpath(i:i) == '/') go to 10
    end do
    i = 0
10  get_filename = fulpath(i+1:klen)

    return
  end function get_filename

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_filename
!
! !INTERFACE:
  character(len=256) function set_filename (rem_dir, loc_fn)
!
! !DESCRIPTION:
!
! !ARGUMENTS:
!
    implicit none
    character(len=*), intent(in)  :: rem_dir !remote directory
    character(len=*), intent(in)  :: loc_fn  !local full path filename
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i   !integer
!------------------------------------------------------------------------

    set_filename = ' '
    do i = len_trim(loc_fn), 1, -1
       if (loc_fn(i:i)=='/') go to 10
    end do
    i = 0
10  set_filename = trim(rem_dir) // loc_fn(i+1:len_trim(loc_fn))

  end function set_filename

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getfil
!
! !INTERFACE:
   subroutine getfil (fulpath, locfn, iflag)
!
! !DESCRIPTION:
! Obtain local copy of file
! First check current working directory
! Next check full pathname[fulpath] on disk
! Finally check full pathname[fulpath] on mass store
!
! !ARGUMENTS:
     implicit none
     character(len=*), intent(in)  :: fulpath !MSS or permanent disk full pathname
     character(len=*), intent(out) :: locfn   !output local file name
     integer, optional, intent(in) :: iflag   !0=>abort if file not found 1=>do not abort
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
     integer i               !loop index
     integer klen            !length of fulpath character string
     integer ierr            !error status
     logical lexist          !true if local file exists
     character(len=256) text !mswrite command
!------------------------------------------------------------------------

     ! get local file name from full name: start at end. look for first "/"

     klen = len_trim(fulpath)
     do i = klen, 1, -1
        if (fulpath(i:i).eq.'/') go to 100
     end do
     i = 0
100  locfn = fulpath(i+1:klen)
     if (len_trim(locfn) == 0) then
        write(6,*)'(GETFIL): local filename has zero length'
        call endrun
     else
        write(6,*)'(GETFIL): attempting to find local file ',  &
             trim(locfn)
     endif

     ! first check if file is in current working directory.

     inquire (file=locfn,exist=lexist)
     if (lexist) then
        write (6,*) '(GETFIL): using ',trim(locfn), &
             ' in current working directory'
        RETURN
     endif

     ! second check for full pathname on disk

     inquire(file=fulpath,exist=lexist)
     if (lexist) then
        locfn = trim(fulpath)
        write(6,*)'(GETFIL): using ',trim(fulpath)
        return
     endif

     ! finally check on mass store

     text='msread '//trim(locfn)//' '//trim(fulpath)
     call shell_cmd(text, ierr)
     if (ierr==0) then
        write(6,*)'(GETFIL): File ',trim(locfn),' read from MSS'
     else  ! all tries to get file have been unsuccessful
        write(6,*)'(GETFIL): failed cmd=',trim(text)
        if (present(iflag) .and. iflag==0) then
           call endrun
        else
           RETURN
        endif
     end if

   end subroutine getfil

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: putfil
!
! !INTERFACE:
   subroutine putfil(locfn, mssfpn, pass, irt, lremov)
!
! !DESCRIPTION:
! Dispose to Mass Store only if nonzero retention period.
! Put mswrite command in background for asynchronous behavior.
! The string put into 'cmd' below needs to be changed to
! the appropriate archival command for the users system
! if a shell command 'mswrite' does not exist.
!
! !ARGUMENTS:
     implicit none
     character(len=*), intent(in) :: locfn   ! Local filename
     character(len=*), intent(in) :: mssfpn  ! Mass Store full pathname
     character(len=*), intent(in) :: pass    ! write password
     integer, intent(in) :: irt              ! Mass Store retention time
     logical, intent(in) :: lremov           ! true=>remove local file
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
     character(len=256) cmd     ! Command string
     character(len=256) cmdtem  ! Temporary for command string
     character(len=  4) crt     ! Retention time as characters
     character(len= 16) wpass   ! Write password
     integer ier                ! error number
!------------------------------------------------------------------------

     if (irt/=0) then
        wpass = ' '
        if (pass(1:1) /= ' ') wpass = ' -w ' // trim(pass)
        write (crt,'(i4)') irt
        write (cmd,'(100a)') 'mswrite ',' -t ',crt,trim(wpass),' ',&
             trim(locfn),' ',trim(mssfpn)
        if (lremov) then
           cmdtem = '('//trim(cmd)//' && /bin/rm '//trim(locfn)//' )&'
        else
           cmdtem = '('//trim(cmd)//' )&'
        end if
        write(6,*)'(PUTFIL): Issuing shell cmd:',trim(cmdtem)
        call shell_cmd(cmdtem, ier)
        if (ier /= 0) then
           write(6,*)'(PUTFIL): Error from shell cmd'
           call endrun
        end if
     endif

   end subroutine putfil

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: opnfil
!
! !INTERFACE:
   subroutine opnfil (locfn, iun, form)
!
! !DESCRIPTION:
! Open file locfn in unformatted or formatted form on unit iun
!
! !ARGUMENTS:
!
     implicit none
     character(len=*), intent(in):: locfn  !file name
     integer, intent(in):: iun             !fortran unit number
     character(len=1), intent(in):: form   !file format: u = unformatted,
                                           !f = formatted
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
     integer ioe             !error return from fortran open
     character(len=11) ft    !format type: formatted. unformatted
!------------------------------------------------------------------------

     if (len_trim(locfn) == 0) then
        write(6,*)'(OPNFIL): local filename has zero length'
        call endrun
     endif
     if (form=='u' .or. form=='U') then
        ft = 'unformatted'
     else
        ft = 'formatted  '
     end if
     open (unit=iun,file=locfn,status='unknown',form=ft,iostat=ioe)
     if (ioe /= 0) then
        write(6,*)'(OPNFIL): failed to open file ',trim(locfn),        &
             &     ' on unit ',iun,' ierr=',ioe
        call endrun
     else
        write(6,*)'(OPNFIL): Successfully opened file ',trim(locfn),   &
             &     ' on unit= ',iun
     end if

   end subroutine opnfil

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getavu
!
! !INTERFACE:
  integer function getavu()
!
! !DESCRIPTION:
! Get next available Fortran unit number.
! If COUP_CSM of OFFLINE is defined,  get next available Fortran unit
! number itst. Set lsmiou(itst). If COUP_CAM is defined, use CAM function
! navu to get available unit number, in which case lsmiou is not needed.
!
! !USES:
#if (defined COUP_CAM)
    use units     !CAM units module
#endif
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Gordon Bonan
! Modified for clm2 by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer itst  !Fortran unit number
!------------------------------------------------------------------------

#if (defined COUP_CAM)
    getavu = getunit()
    RETURN
#else
    do itst = 1, 99
       if (.not.lsmiou(itst)) then
          getavu = itst
          lsmiou(itst) = .true.
          RETURN
       end if
    end do
    write (6,*) 'GETAVU error: ran out of Fortran unit numbers'
    call endrun
#endif
  end function getavu

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: relavu
!
! !INTERFACE:
  subroutine relavu (iunit)
!
! !DESCRIPTION:
! Close and release Fortran unit no longer in use!
! If COUP_CSM or OFFLINE is defined, close and release Fortran unit
! number iunit and set lsmiou(iunit) to false.
! If COUP_CAM is defined, use CAM function relunit to close/release
! unit number.
!
! !USES:
#if (defined COUP_CAM)
    use units     !CAM units module
#endif
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: iunit    !Fortran unit number
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!EOP
!------------------------------------------------------------------------

#if (defined COUP_CAM)
    close(iunit)
    call freeunit(iunit)
#else
    if (.not.lsmiou(iunit)) then
       write (6,*) 'RELAVU eror: unit ',iunit,' is not flagged as in use'
       call endrun
    end if
    if (iunit<1 .or. iunit>99) then
       write (6,*) 'RELAVU error: attempt to return out of range unit'
       call endrun
    end if
    close(iunit)
    lsmiou(iunit) = .false.
#endif

  end subroutine relavu

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: shell_cmd
!
! !INTERFACE:
  subroutine shell_cmd(text, ier)
!
! !DESCRIPTION:
! Invoke shell command
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: text
    integer         , intent(out):: ier
!
! !REVISION HISTORY:
! Created by CAM core group
! Modified for clm2 by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
#if ( defined UNICOSMP )
   integer, external :: ishell ! System routine, execute shell command
#elif ( defined IRIX64 )
   integer, external :: system_cmd  ! Wrapper to "C" system command
#elif (!defined AIX)
   integer, external :: system ! System routine, execute shell command
#endif
!------------------------------------------------------------------------

#if ( defined UNICOSMP )
   ier = ishell(trim(text))
#elif ( defined AIX )
   call system(trim(text), ier)
#elif ( defined IRIX64 )
   ier = system_cmd(trim(text))
#else
   ier = system(trim(text))
#endif

 end subroutine shell_cmd

end module fileutils
