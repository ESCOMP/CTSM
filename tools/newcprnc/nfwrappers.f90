module nfwrappers
!
! Wrapper routines for netcdf library calls.  Abort on failure.
! Arguments are exactly as the library routines called, except where
! an optional argument is allowed (inq_varid and get_att_text). When
! present, if the library call fails then the output variable is set
! to the value of the optional argument.
!
   use utils, only: endrun

   implicit none

PUBLIC

   include 'netcdf.inc'

   integer, private :: ret

CONTAINS

subroutine wrap_create (path, cmode, ncid)
   character(len=*), intent(in) :: path
   integer, intent(in) :: cmode, ncid
  
   ret = nf_create (path, cmode, ncid)
   if (ret /= NF_NOERR) call endrun('wrap_create:'//path//':'//nf_strerror(ret))
   return
end subroutine wrap_create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_open (path, omode, ncid)
   character(len=*), intent(in) :: path
   integer, intent(in) :: omode
   integer, intent(out):: ncid

   ret = nf_open (path, omode, ncid)
   if (ret /= NF_NOERR) then
      call endrun ('wrap_open:'//path//':'//nf_strerror(ret))
   end if
   return
end subroutine wrap_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_close (ncid)
   integer, intent(in) :: ncid

   ret = nf_close (ncid)
   if (ret /= NF_NOERR) then
      write(6,*)'WRAP_CLOSE: nf_close failed for id ', ncid
      call endrun ('wrap_close:'//nf_strerror(ret))
   end if
   return
end subroutine wrap_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_inq_dimid (nfid, dimname, dimid)
   integer, intent(in) :: nfid
   integer, intent(out) :: dimid
   character(len=*), intent(in) :: dimname
  
   ret = nf_inq_dimid (nfid, dimname, dimid)
   if (ret /= NF_NOERR) then
      write(6,*) "WRAP_INQ_DIMID: error "//trim(nf_strerror(ret))
!     call endrun ('wrap_inq_dimid:'//nf_strerror (ret))
   end if
   return
end subroutine wrap_inq_dimid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_inq_ndims (nfid, ndims)
   integer, intent(in) :: nfid
   integer, intent(out) :: ndims
  
   ret = nf_inq_ndims (nfid, ndims)
   if (ret /= NF_NOERR) then
      call endrun ('wrap_inq_ndims:'//nf_strerror (ret))
   end if
   return
end subroutine wrap_inq_ndims

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_inq_dimname (nfid, dimid, dimname)
   integer, intent(in) :: nfid
   integer, intent(in) :: dimid
   character(len=*), intent(out) :: dimname
  
   ret = nf_inq_dimname (nfid, dimid, dimname)
   if (ret /= NF_NOERR) then
      call endrun ('wrap_inq_dimname:'//nf_strerror (ret))
   end if
   return
end subroutine wrap_inq_dimname

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_inq_nvars (nfid, nvars)
   integer, intent(in) :: nfid
   integer, intent(out) :: nvars
  
   ret = nf_inq_nvars (nfid, nvars)
   if (ret /= NF_NOERR) call endrun ('wrap_inq_nvars:'//nf_strerror (ret))
   return
end subroutine wrap_inq_nvars

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_inq_dimlen (nfid, dimid, dimlen)
   integer, intent(in) :: nfid, dimid
   integer, intent(out) :: dimlen
      
   ret = nf_inq_dimlen (nfid, dimid, dimlen)
   if (ret /= NF_NOERR) then
      write(6,*) "WRAP_INQ_DIMLEN: error "//trim(nf_strerror(ret))
!      call endrun ('wrap_inq_dimlen:'//nf_strerror (ret))
   end if
   return
end subroutine wrap_inq_dimlen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_inq_unlimdim (nfid, unlimdimid)
   integer, intent(in) :: nfid
   integer, intent(out) :: unlimdimid
   
   ret = nf_inq_unlimdim (nfid, unlimdimid)
   if (ret /= NF_NOERR) call endrun ('wrap_inq_unlimdim:'//nf_strerror (ret))
   return
end subroutine wrap_inq_unlimdim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_inq_var (nfid, varid, varname, xtype, ndims, dimids, natts)
   integer, intent(in) :: nfid, varid
   integer, intent(out) :: xtype, ndims, dimids(*), natts
   character(len=*), intent(out) :: varname
   
   ret = nf_inq_var (nfid, varid, varname, xtype, ndims, dimids, natts)
   if (ret /= NF_NOERR) call endrun ('wrap_inq_var:'//nf_strerror (ret))
   return
end subroutine wrap_inq_var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_inq_vartype (nfid, varid, xtype)
   integer, intent(in) :: nfid, varid
   integer, intent(out) :: xtype
   
   ret = nf_inq_vartype (nfid, varid, xtype)
   if (ret /= NF_NOERR) call endrun ('wrap_inq_vartype:'//nf_strerror (ret))
   return
end subroutine wrap_inq_vartype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_inq_vardimid (nfid, varid, dimids)
   integer, intent(in) :: nfid, varid
   integer, intent(out) :: dimids(NF_MAX_DIMS)
   
   ret = nf_inq_vardimid (nfid, varid, dimids)
   if (ret /= NF_NOERR) call endrun ('wrap_inq_vardimids:'//nf_strerror (ret))
   return
end subroutine wrap_inq_vardimid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_inq_varid (nfid, varname, varid, defval)
   integer, intent(in) :: nfid
   integer, intent(out) :: varid
   character(len=*), intent(in) :: varname
   integer, intent(in), optional :: defval
   
   ret = nf_inq_varid (nfid, varname, varid)
   if (ret /= NF_NOERR) then
      if (present (defval)) then
         write(6,*)'wrap_inq_varid: setting varid for ', varname, ' to ', defval
         varid = defval
      else
         call endrun ('wrap_inq_varid: varname='//varname//':'//nf_strerror (ret))
      end if
   end if
   return
end subroutine wrap_inq_varid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_inq_varname (nfid, varid, varname)
   integer, intent(in) :: nfid, varid
   character(len=*), intent(out) :: varname

   ret = nf_inq_varname (nfid, varid, varname)
   if (ret /= NF_NOERR) call endrun ('wrap_inq_varname:'//nf_strerror (ret))
   return
end subroutine wrap_inq_varname

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_get_att_text (nfid, varid, attname, atttext, defval)
   integer, intent(in) :: nfid, varid
   character(len=*), intent(in) :: attname
   character(len=*), intent(out) :: atttext
   character(len=*), intent(in), optional :: defval
   
   ret = nf_get_att_text (nfid, varid, attname, atttext)
   if (ret /= NF_NOERR) then
      if (present (defval)) then
         write(6,*)'wrap_get_att_text: setting attname ', attname, ' to default=', defval
         atttext = defval
      else
         write(6,*)'wrap_get_att_text:attname=', attname
         call endrun ('wrap_get_att_text:'//nf_strerror (ret))
      end if
   end if
   return
end subroutine wrap_get_att_text

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_get_var_double (nfid, varid, arr)
   integer, intent(in) :: nfid, varid
   double precision, intent(out) :: arr(:)

   ret = nf_get_var_double (nfid, varid, arr)
   if (ret /= NF_NOERR) then
      call endrun ('wrap_get_var_double:'//nf_strerror (ret))
   end if
   return
end subroutine wrap_get_var_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_get_vara_double (nfid, varid, start, kount, arr)
   integer, intent(in) :: nfid, varid
   integer, intent(in) :: start(*), kount(*)
   double precision, intent(out) :: arr(*)

   ret = nf_get_vara_double (nfid, varid, start, kount, arr)
   if (ret /= NF_NOERR) call endrun ('wrap_get_vara_double:'//nf_strerror (ret))
   return
end subroutine wrap_get_vara_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_get_vara_int (nfid, varid, start, kount, arr)
   integer, intent(in) :: nfid, varid
   integer, intent(in) :: start(*), kount(*)
   integer, intent(out) :: arr(*)
   
   ret = nf_get_vara_int (nfid, varid, start, kount, arr)
   if (ret /= NF_NOERR) call endrun ('wrap_get_vara_int:'//nf_strerror (ret))
   return
end subroutine wrap_get_vara_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_get_var_int (nfid, varid, arr)
   integer, intent(in) :: nfid, varid
   integer, intent(out) :: arr(*)

   ret = nf_get_var_int (nfid, varid, arr)
   if (ret /= NF_NOERR) then
      call endrun ('wrap_get_var_int:'//nf_strerror (ret))
   end if
   return
end subroutine wrap_get_var_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrap_get_vara_text (nfid, varid, start, kount, text)
   integer, intent(in) :: nfid, varid
   integer, intent(in) :: start(*), kount(*)
   character(len=*), intent(out) :: text

   integer ret
   ret = nf_get_vara_text (nfid, varid, start, kount, text)
   if (ret /= NF_NOERR) call endrun ('wrap_get_vara_text:'//nf_strerror (ret))
   return
end subroutine wrap_get_vara_text
end module nfwrappers
