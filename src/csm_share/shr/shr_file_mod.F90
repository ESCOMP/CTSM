!===============================================================================
! CVS $Id$
! CVS $Source$
! CVS $Name$
!===============================================================================
! put/get local files into/from archival location
!===============================================================================

MODULE shr_file_mod

   use shr_kind_mod, only: SHR_KIND_IN
   use shr_sys_mod, only: shr_sys_system   ! System calls
   
   IMPLICIT none

CONTAINS

!===============================================================================

SUBROUTINE shr_file_put(rcode,loc_fn,rem_fn,passwd,rtpd,async,remove)

   !----- arguments -----
   integer(SHR_KIND_IN),intent(out) :: rcode   ! return code
   character(len=*)    ,intent(in)  :: loc_fn  ! local filename
   character(len=*)    ,intent(in)  :: rem_fn  ! remote filename
   character(len=*)    ,optional    :: passwd  ! password
   integer(SHR_KIND_IN),optional    :: rtpd    ! MSS retention period
   logical             ,optional    :: remove  ! true <=> rm after put
   logical             ,optional    :: async   ! true <=> asynchronous put

   !----- local -----  
   character(len=256)   :: passwd2  ! password
   integer(SHR_KIND_IN) :: rtpd2    ! MSS retention period
   logical              :: remove2  ! true <=> rm after put
   logical              :: async2   ! true <=> asynchronous put
   character(len=256)   :: rfn      ! rem_fn without the destination prefix
   character(len=256)   :: cmd      ! command sent to system call

   !----- formats -----
   character(len=*),parameter :: F00 = "('(shr_file_put) ',4a)"
   character(len=*),parameter :: F01 = "('(shr_file_put) ',a,i3,2a)"
   character(len=*),parameter :: F02 = "(a,i4)"

!-------------------------------------------------------------------------------
! PURPOSE: 
!    a generic, extensible put-local-file-into-archive routine
! USAGE:
!    call shr_file_put(rcode,"foo","/home/user/foo")
!    call shr_file_put(rcode,"foo","cp:/home/user/foo",remove=.true.)
!    call shr_file_put(rcode,"foo","mss:/USER/foo",rtpd=3650)
!-------------------------------------------------------------------------------

   remove2 =.false. ; if ( PRESENT(remove )) remove2 = remove
   async2  =.true.  ; if ( PRESENT(async  )) async2  = async
   passwd2 = " "    ; if ( PRESENT(passwd )) passwd2 = passwd
   rtpd2   = 365    ; if ( PRESENT(rtpd   )) rtpd2   = rtpd
   rcode = 0

   if ( index(rem_fn,":") == 0 .or. rem_fn(1:3) == "cp:" ) then
      !------------------------------------------------------
      ! put via unix cp
      !------------------------------------------------------
      rfn = rem_fn
      if ( rem_fn(1:3) == "cp:") rfn = rem_fn(4:len_trim(rem_fn))
      cmd = 'cp -f '//trim(loc_fn)//' '//trim(rfn)
      if (remove2) cmd = trim(cmd)//' && rm '//trim(loc_fn)
      if (async2 ) cmd = trim(cmd)//' & '
      call shr_sys_system(trim(cmd),rcode)
   else if ( rem_fn(1:4) == "mss:" ) then
      !------------------------------------------------------
      ! put onto NCAR's MSS
      !------------------------------------------------------
      if (rtpd2 > 9999) rtpd2 = 9999
      write(cmd,F02) '/usr/local/bin/msrcp -period ',rtpd2
      if (async2 )                cmd = trim(cmd)//' -async '
      if (len_trim(passwd2) > 0 ) cmd = trim(cmd)//' -wpwd '//trim(passwd)
      cmd = trim(cmd)//' '//trim(loc_fn)//' '//trim(rem_fn)
      if (remove2) cmd = trim(cmd)//' && rm '//trim(loc_fn)
      call shr_sys_system(trim(cmd),rcode)
   else if ( rem_fn(1:5) == "hpss:") then
      !------------------------------------------------------
      ! put onto LANL's hpss
      !------------------------------------------------------
      rcode = -1
      cmd = 'rem_fn='//trim(rem_fn)//' loc_fn='//trim(loc_fn)
      write(6,F00) 'ERROR: hpss option not yet implemented'
   else if ( rem_fn(1:5) == "rcp:" ) then
      !------------------------------------------------------
      ! put via rcp
      !------------------------------------------------------
      rcode = -1
      cmd = 'rem_fn='//trim(rem_fn)//' loc_fn='//trim(loc_fn)
      write(6,F00) 'ERROR: rcp option not yet implemented'
   else if ( rem_fn(1:5) == "null:" ) then
      ! do nothing
      cmd = "null prefix => no file archival, do nothing"
      rcode = 0
   else 
      !------------------------------------------------------
      ! unrecognized remote file location
      !------------------------------------------------------
      rcode = -1
      cmd = 'rem_fn='//trim(rem_fn)//' loc_fn='//trim(loc_fn)
      write(6,F00) 'ERROR: unrecognized archive device = ',trim(rem_fn)
   end if

   write(6,F01) 'rcode =',rcode,' cmd = ', trim(cmd)

END SUBROUTINE shr_file_put

!===============================================================================

SUBROUTINE shr_file_get(rcode,loc_fn,rem_fn,passwd,async,clobber)

   !----- arguments -----
   character(len=*)    ,intent(in)  :: loc_fn  ! local filename
   character(len=*)    ,intent(in)  :: rem_fn  ! remote filename
   integer(SHR_KIND_IN),intent(out) :: rcode   ! return code
   character(len=*)    ,optional    :: passwd  ! password
   logical             ,optional    :: async   ! true <=> asynchronous get
   logical             ,optional    :: clobber ! true <=> clobber existing file

   !----- local -----
   character(len=256) :: passwd2  ! password
   logical            :: async2   ! true <=> asynchronous get
   logical            :: clobber2 ! true <=> clobber existing file
   character(len=256) :: rfn      ! rem_fn without the destination prefix
   character(len=256) :: cmd      ! command sent to system call
   logical            :: exists   ! true <=> local file aready exists

   !----- formats -----
   character(len=*),parameter :: F00 = "('(shr_file_get) ',4a)"
   character(len=*),parameter :: F01 = "('(shr_file_get) ',a,i3,2a)"

!-------------------------------------------------------------------------------
! PURPOSE: 
!    a generic, extensible get-local-file-from-archive routine
! USAGE:
!    call shr_file_get(rcode,"foo","/home/user/foo")
!    call shr_file_get(rcode,"foo","cp:/home/user/foo",remove=.true.)
!    call shr_file_get(rcode,"foo","mss:/USER/foo",clobber=.true.)
!-------------------------------------------------------------------------------

   passwd2  = " "     ; if (PRESENT(passwd )) passwd2  = passwd
   async2   = .false. ; if (PRESENT(async  )) async2   = async
   clobber2 = .false. ; if (PRESENT(clobber)) clobber2 = clobber
   rcode = 0

   inquire(file=trim(loc_fn),exist=exists)

   if ( index(rem_fn,":") == 0 .or. rem_fn(1:3) == "cp:" ) then
      !------------------------------------------------------
      ! get via unix cp
      !------------------------------------------------------
      rfn = rem_fn
      if ( rem_fn(1:3) == "cp:") rfn = rem_fn(4:len_trim(rem_fn))
      cmd = 'cp '
      if (clobber2) cmd = 'cp -f '
      cmd = trim(cmd)//' '//trim(rfn)//' '//trim(loc_fn)
      if (async2) cmd = trim(cmd)//' & '
      call shr_sys_system(trim(cmd),rcode)
   else if ( rem_fn(1:4) == "mss:" ) then
      !------------------------------------------------------
      ! get from NCAR's MSS
      !------------------------------------------------------
      cmd = '/usr/local/bin/msrcp '
      if (.not. clobber2) cmd = trim(cmd)//' -noreplace '
      if (async2)   cmd = trim(cmd)//' -async '
      cmd = trim(cmd)//' '//trim(rem_fn)//' '//trim(loc_fn)
      call shr_sys_system(trim(cmd),rcode)
   else if ( rem_fn(1:5) == "hpss:") then
      !------------------------------------------------------
      ! get from LANL's hpss
      !------------------------------------------------------
      rcode = -1
      cmd = 'rem_fn='//trim(rem_fn)//' loc_fn='//trim(loc_fn)
      write(6,F00) 'ERROR: hpss option not yet implemented'
   else if ( rem_fn(1:5) == "rcp:" ) then
      !------------------------------------------------------
      ! get via rcp
      !------------------------------------------------------
      rcode = -1
      cmd = 'rem_fn='//trim(rem_fn)//' loc_fn='//trim(loc_fn)
      write(6,F00) 'ERROR: rcp option not yet implemented'
   else if ( rem_fn(1:5) == "null:" ) then
      ! do nothing
      cmd = "null prefix => no file retrieval, do nothing"
      rcode = 0
   else 
      !------------------------------------------------------
      ! unrecognized remote file location
      !------------------------------------------------------
      rcode = -1
      cmd = 'rem_fn='//trim(rem_fn)//' loc_fn='//trim(loc_fn)
      write(6,F00) 'ERROR: unrecognized archive device = ',trim(rem_fn)
   end if

   write(6,F01) 'rcode =',rcode,' cmd = ', trim(cmd)

END SUBROUTINE shr_file_get

!===============================================================================

END MODULE shr_file_mod
