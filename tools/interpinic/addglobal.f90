subroutine addglobal (ncid, cmdline)

  implicit none
  include 'netcdf.inc'

  ! ------------------------ arguments------------------------------------
  integer, intent(in) :: ncid
  character(len=*), intent(in) ::  cmdline
  !-----------------------------------------------------------------------

  ! ------------------------ local variables -----------------------------
  integer :: ret
  integer :: numchars
  integer :: values(8)
  integer :: hnum
  integer :: hlen
  character(len= 8) :: date
  character(len=10) :: time
  character(len= 5) :: zone
  character(len=18) :: datetime
  character(len=16) :: logname
  character(len=16) :: hostname
  character(len=1500) :: hist
  !-----------------------------------------------------------------------

  call date_and_time (date, time, zone, values)

  datetime(1:8) =        date(5:6) // '/' // date(7:8) // '/' // date(3:4)
  datetime(9:)  = ' ' // time(1:2) // ':' // time(3:4) // ':' // time(5:6) // ' '

  call getenv ('LOGNAME', logname)
  call getenv ('HOST', hostname)

  hlen = 0
  hist = ' '
  if (nf_inq_attid (ncid, nf_global, 'history', hnum) == nf_noerr) then
     ret = nf_inq_attlen (ncid, nf_global, 'history', hlen)
     call wrap_get_att_text (ncid, nf_global, 'history', hist)
  end if

  hist = trim (hist) // char(10) // datetime // trim (logname) // ':' // &
       trim (hostname) // ':' // trim (cmdline)

  ! Add "3" to account for first newline and colons between each of 2 trimmed strings

  hlen = hlen + len(datetime) + len_trim(logname) + len_trim(hostname) + len_trim(cmdline) + 3

  if (hlen > len(hist)) then
     write(6,*)'Warning: history attribute too long: truncating'
     hlen = len(hist)
  end if

  numchars = len_trim (hist)
  ret = nf_put_att_text (ncid, nf_global, 'history', numchars, hist)

  return
end subroutine addglobal
  
 
