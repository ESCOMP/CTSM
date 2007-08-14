subroutine wrap_inq_varid (nfid, varname, varid)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  character*(*) varname

  integer ret

  ret = nf_inq_varid (nfid, varname, varid)
  if (ret /= NF_NOERR) then
    write(6,*)'nf_inq_varid: ',trim(varname),' not found'
    call handle_error (ret)
  end if
end subroutine wrap_inq_varid

subroutine wrap_inq_var (nfid, varid, varname, xtype, ndims, dimids, natts)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, xtype, ndims, dimids(nf_max_dims), natts
  character*(*) varname

  integer ret

  ret = nf_inq_var (nfid, varid, varname, xtype, ndims, dimids, natts)
  if (ret /= NF_NOERR) call handle_error (ret)
end subroutine wrap_inq_var

subroutine wrap_get_var_double (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  real*8 arr(*)

  integer ret

  ret = nf_get_var_double (nfid, varid, arr)
  if (ret /= NF_NOERR) call handle_error (ret)
end subroutine wrap_get_var_double

subroutine wrap_get_var_int (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer arr(*)

  integer ret

  ret = nf_get_var_int (nfid, varid, arr)
  if (ret /= NF_NOERR) call handle_error (ret)
end subroutine wrap_get_var_int

subroutine wrap_get_vara_double (nfid, varid, start, count, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, start(nf_max_dims), count(nf_max_dims)
  real*8 arr(*)

  integer ret

  ret = nf_get_vara_double (nfid, varid, start, count, arr)
  if (ret /= NF_NOERR) call handle_error (ret)
end subroutine wrap_get_vara_double

subroutine wrap_get_vara_int (nfid, varid, start, count, arr)
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, start(*), count(*)
  integer arr(*)

  integer ret

  ret = nf_get_vara_int (nfid, varid, start, count, arr)
  if (ret /= NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_put_vara_text (nfid, varid, start, count, text)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer start(*), count(*)
  character*(*) text

  integer ret

  ret =      nf_put_vara_text (nfid, varid, start, count, text)
  if (ret /= NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_get_vara_text (nfid, varid, start, count, text)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer start(*), count(*)
  character*(*) text

  integer ret

  ret =      nf_get_vara_text (nfid, varid, start, count, text)
  if (ret /= NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_put_var_double (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  real*8 arr(*)

  integer ret
  ret = nf_put_var_double (nfid, varid, arr)
  if (ret /= NF_NOERR) call handle_error (ret)
end subroutine wrap_put_var_double

subroutine wrap_put_vara_double (nfid, varid, start, count, arr)
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer start(*), count(*)
  real(r8) arr(*)
  
  integer ret
  ret =      nf_put_vara_double (nfid, varid, start, count, arr)
  if (ret /= NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_put_var_int (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer arr(*)

  integer ret
  ret = nf_put_var_int (nfid, varid, arr)
  if (ret /= NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_put_vara_int (nfid, varid, start, count, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, start(*), count(*)
  integer arr(*)

  integer ret
  ret =      nf_put_vara_int (nfid, varid, start, count, arr)
  if (ret /= NF_NOERR) call handle_error (ret)
end subroutine

subroutine handle_error (ret)
  implicit none
  include 'netcdf.inc'

  integer ret

  write(6,*) "NetCDF error code = ", ret
  write(6,*) nf_strerror (ret)
  call abort
end subroutine handle_error
