#!/bin/csh

set filename = ncdio_global_subs.inc

touch $filename
mv -f $filename $filename.old
touch $filename

foreach DIM (var 1d 2d 3d)
foreach TYPE (int real char)

  set DLEN1 = "size(data)"
  set DLEN2 = "size(data,dim=n)"
  if ($DIM == var) then
     set CDIMS = "       "
     set NDIMS = 0
     set DLEN1 = "1"
     set DLEN2 = "1"
  endif
  if ($DIM == 1d) then
     set CDIMS = "(:)    "
     set NDIMS = 1
  endif
  if ($DIM == 2d) then
     set CDIMS = "(:,:)  "
     set NDIMS = 2
  endif
  if ($DIM == 3d) then
     set CDIMS = "(:,:,:)"
     set NDIMS = 3
  endif

  if ($TYPE == real) then
     set ATYPE = "real(r8)        "
     set PTYPE = "real(r8)        "
     set NCTYPE = "double"
     set MPTYPE = "MPI_REAL8"
  endif
  if ($TYPE == int) then
     set ATYPE = "integer         "
     set PTYPE = "integer         "
     set NCTYPE = "int"
     set MPTYPE = "MPI_INTEGER"
  endif
  if ($TYPE == char) then
     set ATYPE = "character(len=*)"
     set PTYPE = "character"
     set NCTYPE = "text"
     set MPTYPE = "MPI_CHARACTER"
     set CDIMS = " "
     set DLEN1 = "len(data)"
     set DLEN2 = "len(data)"
  endif

  set gen_sub = "true"
  if ($TYPE == char && $DIM != 1d) then
     set gen_sub = "false"
  endif

# if then doesn't work, use goto
#  if ($gen_sub == true) then
  if ($gen_sub != "true") goto skipsub

echo $DIM $TYPE

cat >> $filename <<EOF

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_${TYPE}_${DIM}
!
! !INTERFACE:
  subroutine ncd_ioglobal_${TYPE}_${DIM}(varname, data, flag, ncid, readvar, nt, bcast, usepio)
!
! !DESCRIPTION:
! netcdf I/O of global ${DIM} ${TYPE} array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    ${ATYPE}, intent(inout) :: data ${CDIMS}     ! raw data
    logical         , optional, intent(out):: readvar   ! was var read?
    integer         , optional, intent(in) :: nt        ! time sample index
    logical         , optional, intent(in) :: bcast     ! bcast on read?
    logical         , optional, intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    ${PTYPE},pointer :: piodata(:)  ! copy of data in 1d
    integer :: n,n1,n2,n3           ! local counter
    integer :: varid                ! netCDF variable id
    integer :: vdnum                ! vardesc number
#if (defined BUILDPIO)
    type(var_desc_t),pointer :: pio_vardesc   ! local vardesc pointer
#endif
    integer :: ier                  ! error code
    integer :: start(4), count(4)   ! output bounds
    integer :: nd,did(4),ld(4)      ! var/dim error checking
    character(len=32) :: vname      ! variable error checking
    character(len=32) :: dname      ! dimension error checking
    logical :: varpresent           ! if true, variable is on tape
    logical :: lbcast               ! local copy of bcast flag
    logical :: lusepio              ! local usepio variable
    integer,parameter :: ndims = ${NDIMS}  ! data dims
    character(len=*),parameter :: subname='ncd_ioglobal_${TYPE}_${DIM}'
!-----------------------------------------------------------------------

    lusepio = pio_def
    if (present(usepio)) then
       lusepio = usepio
    endif

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(varname),lusepio
    endif

    if (lusepio) then
       n = ${DLEN1}
       allocate(piodata(n))
EOF

if ($TYPE == char) then
cat >> $filename <<EOF
       do n1 = 1,n
          piodata(n1) = data(n1:n1)
       enddo
EOF
endif

if ($TYPE != char && $NDIMS == 0) then
cat >> $filename <<EOF
       piodata(1) = data
EOF
endif

if ($TYPE != char && $NDIMS == 1) then
cat >> $filename <<EOF
       piodata = data
EOF
endif

if ($TYPE != char && $NDIMS == 2) then
cat >> $filename <<EOF
       n = 0
       do n2 = 1, size(data,dim=2)
       do n1 = 1, size(data,dim=1)
          n = n + 1
          piodata(n) = data(n1,n2)
       enddo
       enddo
EOF
endif

if ($TYPE != char && $NDIMS == 3) then
cat >> $filename <<EOF
       n = 0
       do n3 = 1, size(data,dim=3)
       do n2 = 1, size(data,dim=2)
       do n1 = 1, size(data,dim=1)
          n = n + 1
          piodata(n) = data(n1,n2,n3)
       enddo
       enddo
       enddo
EOF
endif

if ($TYPE != char && $NDIMS >= 4) then
cat >> $filename <<EOF
          if (masterproc) write(iulog,*) trim(subname),' ERROR: allocate of piodata'
          call endrun()
EOF
endif

cat >> $filename <<EOF
    endif

       start = 1
       count = 1
       lbcast = lbcast_def
       if (present(bcast)) then
          lbcast = bcast
       endif

       if (flag == 'write') then
          if (lusepio) then
#if (defined BUILDPIO)
             call ncd_inqvdesc(varname, vdnum, subname, usepio=lusepio)
             pio_vardesc => pio_vardesc_list(vdnum)%pio_varDesc
             if (vdnum < 1 .or. vdnum > pio_num_vardesc) then
                write(iulog,*) trim(subname),' ERROR in vdnum from inqvdesc ',trim(varname),vdnum
                call endrun()
             endif
             call ncd_inqvid(ncid, varname, varid, subname, usepio=lusepio, pio_varDesc=pio_varDesc)
             call ncd_inqvname(ncid, varid, vname, subname, usepio=lusepio, pio_varDesc=pio_varDesc)
             call ncd_inqvdims(ncid, varid, nd, subname, usepio=lusepio, pio_varDesc=pio_varDesc)
             call ncd_inqvdids(ncid, varid, did, subname, usepio=lusepio, pio_varDesc=pio_varDesc)
#endif
          else
             call ncd_inqvid(ncid, varname, varid, subname, usepio=lusepio)
             call ncd_inqvname(ncid, varid, vname, subname, usepio=lusepio)
             call ncd_inqvdims(ncid, varid, nd, subname, usepio=lusepio)
             call ncd_inqvdids(ncid, varid, did, subname, usepio=lusepio)
          endif
          if (masterproc .and. (trim(varname) /= trim(vname))) then
             write(iulog,*) trim(subname),' ERROR: varnames do not match ',trim(varname),' ',trim(vname)
             call endrun()
          endif
          if (masterproc .and. (nd - ndims > 1 .or. nd - ndims < 0)) then
             write(iulog,*) trim(subname),' ERROR: array ndims ne cdf var ndims ',trim(varname),ndims,nd
             call endrun()
          endif
          do n = 1, ndims
             call ncd_inqdlen(ncid, did(n), ld(n), subname, usepio=lusepio)
             call ncd_inqdname(ncid, did(n), dname, subname, usepio=lusepio)
             count(n) = ${DLEN2}
             if (masterproc .and. count(n) /= ld(n)) then
                write(iulog,*) trim(subname),' ERROR: array size ne cdf var size ',trim(varname),n,trim(dname),count(n),ld(n)
                call endrun()
             endif
          enddo
          if (present(nt)) then
             start(ndims+1) = nt
          endif
          if (lusepio) then
#if (defined BUILDPIO)
             call check_ret_pio(PIO_put_var(pio_file, varid, start, count, piodata), subname)
#endif
          else
             if (masterproc) call check_ret(nf_put_vara_${NCTYPE}(ncid, varid, start, count, data), subname)
          endif

       else if (flag == 'read') then
          call ncd_inqvid(ncid, varname, varid, subname, readvar=varpresent, usepio=lusepio)
          if (varpresent) then
             if (single_column) then
                call scam_field_offsets(ncid,'undefined',start,count)
                if (lusepio) then
#if (defined BUILDPIO)
!                   call check_ret_pio(PIO_get_var(ncid, varid, start, count, piodata), subname)
                    if (masterproc) write(iulog,*) trim(subname),' pio not implemented'
                    call endrun()
#endif
                else
                   if (masterproc) call check_ret(nf_get_vara_${NCTYPE}(ncid, varid, start, count, data), subname)
                endif
             else
                if (lusepio) then
#if (defined BUILDPIO)
!                   call check_ret_pio(PIO_get_var(ncid, varid, piodata), subname)
                    if (masterproc) write(iulog,*) trim(subname),' pio not implemented'
                    call endrun()
#endif
                else
                   if (masterproc) call check_ret(nf_get_var_${NCTYPE}(ncid, varid, data), subname)
                endif
             endif
          endif

          if (lbcast) then
             call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
             if (ier /= 0) then
                write(iulog,*) trim(subname), &
                               ' ERROR from mpi_bcast for varpresent'
                call endrun()
             endif
             if (varpresent) then
                call mpi_bcast(data, ${DLEN1}, ${MPTYPE}, 0, mpicom, ier)
                if (ier /= 0) then
                   write(iulog,*) trim(subname), &
                                  ' ERROR from mpi_bcast for data'
                   call endrun()
                endif
             endif
          endif
          if (present(readvar)) readvar = varpresent

       endif   ! flag

    if (lusepio) then
       deallocate(piodata)
    endif

  end subroutine ncd_ioglobal_${TYPE}_${DIM}

EOF

skipsub:

end
end

exit
