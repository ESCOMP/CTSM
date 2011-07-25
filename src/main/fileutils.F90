
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
  use clm_varctl, only : iulog
  use spmdMod   , only : masterproc
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: get_filename  !Returns filename given full pathname
  public :: set_filename  !Set remote full path filename
  public :: opnfil        !Open local unformatted or formatted file
  public :: getfil        !Obtain local copy of file
  public :: putfil        !Dispose file to archival system
  public :: relavu        !Close and release Fortran unit no longer in use
  public :: getavu        !Get next available Fortran unit number
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !PRIVATE MEMBER FUNCTIONS: None
!EOP
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
!
! !LOCAL VARIABLES:
!EOP
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
!
! !LOCAL VARIABLES:
!EOP
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
! Finally check full pathname[fulpath] on archival system
! 
! !USES:
     use shr_file_mod, only: shr_file_get
     use clm_varctl  , only: fget_archdev
! !ARGUMENTS:
     implicit none
     character(len=*), intent(in)  :: fulpath !Archival or permanent disk full pathname
     character(len=*), intent(out) :: locfn   !output local file name
     integer, optional, intent(in) :: iflag   !0=>abort if file not found 1=>do not abort
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
     integer i               !loop index
     integer klen            !length of fulpath character string
     integer ierr            !error status
     logical lexist          !true if local file exists
     character(len=len(fulpath)+5)  :: fulpath2 !Archival full pathname
!------------------------------------------------------------------------

     ! get local file name from full name: start at end. look for first "/"

     klen = len_trim(fulpath)
     do i = klen, 1, -1
        if (fulpath(i:i).eq.'/') go to 100
     end do
     i = 0
100  locfn = fulpath(i+1:klen)
     if (len_trim(locfn) == 0) then
	if (masterproc) write(iulog,*)'(GETFIL): local filename has zero length'
        call endrun
     else
        if (masterproc) write(iulog,*)'(GETFIL): attempting to find local file ',  &
             trim(locfn)
     endif

     ! first check if file is in current working directory.

     inquire (file=locfn,exist=lexist)
     if (lexist) then
        if (masterproc) write(iulog,*) '(GETFIL): using ',trim(locfn), &
             ' in current working directory'
        RETURN
     endif

     ! second check for full pathname on disk if no prepended "type:"

     if ( index(fulpath,":") == 0 )then
        inquire(file=fulpath,exist=lexist)
        if (lexist) then
           locfn = trim(fulpath)
           if (masterproc) write(iulog,*)'(GETFIL): using ',trim(fulpath)
           return
        endif
        fulpath2 = trim(fget_archdev)//trim(fulpath)
     else
        fulpath2 = trim(fulpath)
     end if

     ! finally check on full archive path location

     call shr_file_get( ierr, locfn, fulpath2 )
     if (ierr==0) then
        if (masterproc) write(iulog,*)'(GETFIL): File ',trim(locfn),' read in from: ', fulpath2
     else  ! all tries to get file have been unsuccessful
        if (masterproc) write(iulog,*)'(GETFIL): failed getting file from full path: ', fulpath2
        if (present(iflag) .and. iflag==0) then
           call endrun ('GETFIL: FAILED to get '//trim(fulpath2))
        else
           RETURN
        endif
     end if

     ! And now make sure file was successfully transfered

     inquire (file=locfn,exist=lexist)
     if ( .not. lexist) then
        if (masterproc) write(iulog,*)'(GETFIL): failed transferring file to local path: ', locfn
        if (present(iflag) .and. iflag==0) then
           call endrun ('GETFIL: file not transfered to local path' )
        end if
     endif

   end subroutine getfil

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: putfil
!
! !INTERFACE:
   subroutine putfil(locfn, fulpath, pass, irt, lremov)
!
! !DESCRIPTION:
! Dispose to archival system only if nonzero retention period.
! Put mswrite command in background for asynchronous behavior.
! The string put into 'cmd' below needs to be changed to
! the appropriate archival command for the users system
! if a shell command 'mswrite' does not exist.
!
! !USES:
     use shr_file_mod, only: shr_file_put
! !ARGUMENTS:
     implicit none
     character(len=*), intent(in) :: locfn   ! Local filename
     character(len=*), intent(in) :: fulpath ! archive full pathname
     character(len=*), intent(in) :: pass    ! write password
     integer, intent(in) :: irt              ! Archival system retention time
     logical, intent(in) :: lremov           ! true=>remove local file
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
     integer ier                ! error number
!------------------------------------------------------------------------

     if (irt/=0) then
        if (pass(1:1) /= ' ')then
           call shr_file_put( ier, locfn, fulpath, passwd=pass, rtpd=irt, &
                              async=.false., remove=lremov )
        else
           call shr_file_put( ier, locfn, fulpath, rtpd=irt, async=.false., &
                              remove=lremov )
        end if

        if (ier /= 0) then
           call endrun ('PUTFIL: Error from shell shr_file_put')
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
!
! !LOCAL VARIABLES:
!EOP
     integer ioe             !error return from fortran open
     character(len=11) ft    !format type: formatted. unformatted
!------------------------------------------------------------------------

     if (len_trim(locfn) == 0) then
        write(iulog,*)'(OPNFIL): local filename has zero length'
        call endrun
     endif
     if (form=='u' .or. form=='U') then
        ft = 'unformatted'
     else
        ft = 'formatted  '
     end if
     open (unit=iun,file=locfn,status='unknown',form=ft,iostat=ioe)
     if (ioe /= 0) then
        write(iulog,*)'(OPNFIL): failed to open file ',trim(locfn),        &
             &     ' on unit ',iun,' ierr=',ioe
        call endrun
     else if ( masterproc )then
        write(iulog,*)'(OPNFIL): Successfully opened file ',trim(locfn),   &
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
!
! !USES:
   use shr_file_mod, only : shr_file_getUnit
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Gordon Bonan
! Modified for clm2 by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
!------------------------------------------------------------------------

    getavu = shr_file_getunit()

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
!
! !USES:
   use shr_file_mod, only : shr_file_freeUnit
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

    close(iunit)
    call shr_file_freeUnit(iunit)

  end subroutine relavu

end module fileutils
