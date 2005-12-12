subroutine wrap_create (path, cmode, ncid)
  implicit none
  include 'netcdf.inc'
      
  character*(*) path
  integer cmode, ncid
  
  integer ret

  ret = nf_create (path, cmode, ncid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_inq_dimid (nfid, dimname, dimid)
  implicit none
  include 'netcdf.inc'
  
  integer nfid, dimid
  character*(*) dimname
  
  integer ret
  
  ret = nf_inq_dimid (nfid, dimname, dimid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_inq_nvars (nfid, nvars)
  implicit none
  include 'netcdf.inc'
  
  integer nfid, nvars
  
  integer ret
  
  ret = nf_inq_nvars (nfid, nvars)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end

subroutine wrap_inq_dimlen (nfid, dimid, dimlen)
  implicit none
  include 'netcdf.inc'
  
  integer nfid, dimid, dimlen
      
  integer ret
      
  ret = nf_inq_dimlen (nfid, dimid, dimlen)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_inq_unlimdim (nfid, unlimdimid)
  implicit none
  include 'netcdf.inc'

  integer nfid, unlimdimid

  integer ret

  ret = nf_inq_unlimdim (nfid, unlimdimid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_inq_var (nfid, varid, varname, xtype, ndims, dimids, natts)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, xtype, ndims, dimids(*), natts
  character*(*) varname

  integer ret

  ret = nf_inq_var (nfid, varid, varname, xtype, ndims, dimids, natts)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_inq_varid (nfid, varname, varid)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  character*(*) varname

  integer ret

  ret = nf_inq_varid (nfid, varname, varid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_inq_varname (nfid, varid, varname)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  character*(*) varname

  integer ret
  
  ret = nf_inq_varname (nfid, varid, varname)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_get_att_text (nfid, varid, attname, atttext)
  implicit none
      
  include 'netcdf.inc'

  integer nfid, varid
  character*(*) attname, atttext

  integer ret

  ret = nf_get_att_text (nfid, varid, attname, atttext)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_get_var_realx (nfid, varid, arr)
  use precision
  implicit none

  include 'netcdf.inc'

  integer nfid, varid
  real(r8) arr(*)

  integer ret
  ret = nf_get_var_double (nfid, varid, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_get_vara_realx (nfid, varid, start, count, arr)
  use precision
  implicit none

  include 'netcdf.inc'

  integer nfid, varid
  integer start(*), count(*)
  real(r8) arr(*)

  integer ret
  ret = nf_get_vara_double (nfid, varid, start, count, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_get_vara_int (nfid, varid, start, count, arr)
  use precision
  implicit none

  include 'netcdf.inc'

  integer nfid, varid
  integer start(*), count(*)
  integer arr(*)

  integer ret
  ret = nf_get_vara_int (nfid, varid, start, count, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_get_var_int (nfid, varid, arr)
  use precision
  implicit none

  include 'netcdf.inc'

  integer nfid, varid
  integer arr(*)

  integer ret
  ret = nf_get_var_int (nfid, varid, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_get_vara_text (nfid, varid, start, count, text)
  use precision
  implicit none

  include 'netcdf.inc'

  integer nfid, varid
  integer start(*), count(*)
  character*(*) text(*)

  integer ret
  ret = nf_get_vara_text (nfid, varid, start, count, text)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_open (path, omode, ncid)
  implicit none
  include 'netcdf.inc'
      
  character*(*) path
  integer omode
  integer ncid
  integer ret

  ret = nf_open (path, omode, ncid)
  if (ret.ne.NF_NOERR) then
    write(6,*)'WRAP_OPEN: nf_open failed for file ',path
    call handle_error (ret)
  end if
end subroutine

subroutine wrap_close (ncid)
  implicit none
  include 'netcdf.inc'

  integer ncid
  integer ret

  ret = nf_close (ncid)
  if (ret.ne.NF_NOERR) then
    write(6,*)'WRAP_CLOSE: nf_close failed for id ',ncid
    call handle_error (ret)
  end if
end subroutine

subroutine handle_error(ret)
  implicit none
  include 'netcdf.inc'
  
  integer ret
      
  write(6,*)nf_strerror(ret)
  call abort
end subroutine

