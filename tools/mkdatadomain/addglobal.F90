!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: addglobal -- add global attributes to file.
!
! !INTERFACE:
   subroutine addglobal (ncid)
!
! !DESCRIPTION:
! Add several global attributes to the output file, including: version, revision_id,
! history with date, hostname, and user creating file.
!
! !USES:
  use shr_kind_mod  , only : SHR_KIND_CL
  implicit none
  include 'netcdf.inc'
!
! !ARGUMENTS:
  integer, intent(in) :: ncid
!
! !LOCAL VARIABLES:
  integer :: numchars
  integer :: values(8)
  integer :: hnum
  integer :: hlen
  character(len= 8) :: date
  character(len=10) :: time
  character(len= 5) :: zone
  character(len=18) :: datetime
  character(len=SHR_KIND_CL) :: version = &
  "$HeadURL$"
  character(len=SHR_KIND_CL) :: revision_id = &
  "$Id$"
  character(len=16) :: logname
  character(len=16) :: hostname
  character(len=SHR_KIND_CL) :: str
  character(len=1500) :: hist
  character(len=32) :: subname = "addglobal"
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

  call date_and_time (date, time, zone, values)

  datetime(1:8) =        date(5:6) // '/' // date(7:8) // '/' // date(3:4)
  datetime(9:)  = ' ' // time(1:2) // ':' // time(3:4) // ':' // time(5:6) // ' '

  call getenv ('LOGNAME', logname)
  call getenv ('HOST', hostname)

  str = 'NCAR-CESM:CF-1.0'
  call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
       'Conventions', len_trim(str), trim(str)), subname)

  str = 'CESM domain data:'
  call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
       'title', len_trim(str), trim(str)), subname)

  str = 'Standard CESM1.0 domain specification file created from CLM inputdata files:'
  call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
       'user_comment', len_trim(str), trim(str)), subname)

  str = 'from CLM fraction and griddata files'
  call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
       'source', len_trim(str), trim(str)), subname)

  hlen = 0
  hist = ' '
  if (nf_inq_attid (ncid, nf_global, 'history', hnum) == nf_noerr) then
     call check_ret(  nf_inq_attlen (ncid, nf_global, 'history', hlen), subname )
     call check_ret(  nf_get_att_text (ncid, nf_global, 'history', hist), subname )
  end if

  hist = trim (hist) // char(10) // datetime // trim (logname) // ':' // trim (hostname)

  ! Add "3" to account for first newline and colons between each of 2 trimmed strings

  hlen = hlen + len(datetime) + len_trim(logname) + len_trim(hostname) + 3

  if (hlen > len(hist)) then
     write(6,*)'Warning: history attribute too long: truncating'
     hlen = len(hist)
  end if

  numchars = len_trim (hist)
  call check_ret(  nf_put_att_text (ncid, nf_global, 'history', numchars, hist), subname )

  write(6,*) "Add SVN_version and Id to global file attributes"
  numchars = len_trim (version)
  call check_ret(  nf_put_att_text (ncid, nf_global, 'mkdatadomain_version', numchars, version), subname )
  call check_ret(  nf_put_att_text (ncid, nf_global, 'SVN_url',              numchars, version), subname )
  numchars = len_trim (revision_id)
  call check_ret(  nf_put_att_text (ncid, nf_global, 'mkdatadomain_version_Id', numchars, revision_id), subname )
  call check_ret(  nf_put_att_text (ncid, nf_global, 'source_code',             numchars, revision_id), subname )

  write(6,*) "Done adding global attributes"

  return
end subroutine addglobal
  
 
