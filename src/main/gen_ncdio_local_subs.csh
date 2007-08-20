#!/bin/csh

set filename = ncdio_local_subs.inc

touch ncdio.F90

touch $filename
mv -f $filename $filename.old
touch $filename

foreach TYPE (int real)

  if ($TYPE == real) then
     set ATYPE = "real(r8)"
     set NCTYPE = "double"
     set MPTYPE = "MPI_REAL8"
     set SPVAL = "spval"
  endif
  if ($TYPE == int) then
     set ATYPE = "integer "
     set NCTYPE = "int"
     set MPTYPE = "MPI_INTEGER"
     set SPVAL = "ispval"
  endif

echo $TYPE

cat >> $filename <<EOF
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_${TYPE}_1d
!
! !INTERFACE:
  subroutine ncd_iolocal_${TYPE}_1d(varname, data, dim1name, &
       flag, ncid, nlonxy, nlatxy, nt, readvar, missing, usepio)
!
! !DESCRIPTION:
! I/O for 1d int field
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
    integer         , intent(in)  :: ncid               ! input unit
    character(len=*), intent(in)  :: varname            ! variable name
    ${ATYPE}        , pointer     :: data(:)            ! local decomposition data
    character(len=*), intent(in)  :: dim1name           ! dimension name
    integer         , optional, intent(in) :: nlonxy    ! 2d longitude size
    integer         , optional, intent(in) :: nlatxy    ! 2d latitude size
    integer         , optional, intent(in) :: nt        ! time sample index
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    ${ATYPE}        , optional, intent(in) :: missing   ! missing value
    logical         , optional, intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: dims                     ! dimensions
    integer :: gsize                    ! size of global array
    integer :: ier                      ! error status
    integer :: start(3)                 ! starting indices for netcdf field
    integer :: count(3)                 ! count values for netcdf field
    ${ATYPE}:: lmissing                 ! local missing value
    logical :: lusepio           ! local usepio variable
    character(len=8) :: clmlevel        ! clmlevel
    character(len=*),parameter :: subname='ncd_iolocal_${TYPE}_1d' ! subroutine name
!-----------------------------------------------------------------------

    lusepio = pio_def
    if (present(usepio)) then
       lusepio = usepio
    endif

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(flag),' ',trim(varname),lusepio
    endif

    if ((present(nlonxy) .and. .not.present(nlatxy)) .or. &
        (present(nlatxy) .and. .not.present(nlonxy))) then
       write(iulog,*) trim(subname),' error nlonxy/nlatxy must be both or neither present '
       call endrun()
    endif

    if (present(missing)) then
       lmissing = missing
    else
       lmissing = ${SPVAL}
    endif

    clmlevel = dim1name
    if (present(nlonxy) .and. present(nlatxy)) then
       if (dim1name == nameg .or. dim1name == grlnd) then
          clmlevel = grlnd
       elseif (dim1name == allrof .or. dim1name == gratm) then
          ! continue, acceptable and default behavior for now
       else
          if (masterproc) write(iulog,*) trim(subname),' warning incorrect use of dim1name and nlonxy/nlatxy ',trim(dim1name),nlonxy,nlatxy
       endif
    endif

    gsize = get_clmlevel_gsize(clmlevel)
    start = 1
    count = 1
    call get_clmlevel_dsize(clmlevel,dims,count(1),count(2))
    if (dims == 1) then
       if (present(nt)) then
          start(2) = nt
       endif
    elseif (dims == 2) then
       if (present(nt)) then
          start(3) = nt
       endif
    else
       write(iulog,*) trim(subname),' error dims incorrect ',clmlevel,dims
       call endrun()
    endif

    call ncd_iolocal_gs_${TYPE}1d(ncid, varname, flag, data, clmlevel, start, count, ier, lmissing, usepio=lusepio)

    if (present(readvar)) then
       readvar = .false.
       if (ier == 0) readvar = .true.
    endif

  end subroutine ncd_iolocal_${TYPE}_1d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_${TYPE}_2d
!
! !INTERFACE:
  subroutine ncd_iolocal_${TYPE}_2d(varname, data, dim1name, dim2name, &
             lowerb2, upperb2, flag, ncid, nlonxy, nlatxy, nt, readvar, missing, usepio)
!
! !DESCRIPTION:
! Netcdf i/o of 2d initial integer field out to netCDF file
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
    integer         , intent(in)  :: ncid               ! input unit
    character(len=*), intent(in)  :: varname            ! variable name
    ${ATYPE}        , pointer     :: data(:,:)          ! local decomposition input data
    character(len=*), intent(in)  :: dim1name           ! dimension 1 name
    character(len=*), intent(in)  :: dim2name           ! dimension 2 name
    integer         , optional, intent(in) :: nlonxy    ! 2d longitude size
    integer         , optional, intent(in) :: nlatxy    ! 2d latitude size
    integer         , optional, intent(in) :: nt        ! time sample index
    integer         , optional, intent(in) :: lowerb2,upperb2 ! lower and upper bounds of second dimension
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    ${ATYPE}        , optional, intent(in) :: missing   ! missing value
    logical         , optional, intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: k                        ! index
    integer :: dims                     ! dimensions
    integer :: gsize                    ! size of global array
    integer :: ier                      ! error status
    integer :: start(4)                 ! starting indices for netcdf field
    integer :: count(4)                 ! count values for netcdf field
    integer :: lb1,ub1                  ! lower/upper bound of dim 1
    integer :: lb2,ub2                  ! lower/upper bound of dim 2
    ${ATYPE}:: lmissing                 ! local missing value
    logical :: lusepio           ! local usepio variable
    ${ATYPE},pointer :: data1d(:)       ! 1 level data
    character(len=8) :: clmlevel        ! clmlevel
    character(len=*),parameter :: subname='ncd_iolocal_${TYPE}_2d' ! subroutine name
!-----------------------------------------------------------------------

    lusepio = pio_def
    if (present(usepio)) then
       lusepio = usepio
    endif

    if (lusepio .and. lowmem2d) then
       write(iulog,*) trim(subname),' ERROR usepio and lowmem2d are both true'
       call endrun()
    endif

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(flag),' ',trim(varname),lusepio
    endif

    if ((present(nlonxy) .and. .not.present(nlatxy)) .or. &
        (present(nlatxy) .and. .not.present(nlonxy))) then
       write(iulog,*) trim(subname),' error nlonxy/nlatxy must be both or neither present '
       call endrun()
    endif

    lb1 = lbound(data, dim=1)
    ub1 = ubound(data, dim=1)

    if (present(lowerb2)) then
       lb2 = lowerb2
    else
       lb2 = lbound(data, dim=2)
    end if
    if (present(upperb2)) then
       ub2 = upperb2
    else
       ub2 = ubound(data, dim=2)
    end if

    if (present(missing)) then
       lmissing = missing
    else
       lmissing = ${SPVAL}
    endif

    clmlevel = dim1name
    if (present(nlonxy) .and. present(nlatxy)) then
       if (dim1name == nameg .or. dim1name == grlnd) then
          clmlevel = grlnd
       else
          write(iulog,*) trim(subname),' error in dim1name and nlonxy/nlatxy ',trim(dim1name),nlonxy,nlatxy
          call endrun()
       endif
    endif

    gsize = get_clmlevel_gsize(clmlevel)
    start = 1
    count = 1  
    call get_clmlevel_dsize(clmlevel,dims,count(1),count(2))

    if (lowmem2d) then
       allocate(data1d(lb1:ub1))
       do k = lb2,ub2
          if (dims == 1) then
#if (defined SWITCH_DIMS)
             start(1) = k-lb2+1
             count(1) = 1
             count(2) = gsize
#else
             count(1) = gsize
             start(2) = k-lb2+1
             count(2) = 1
#endif
             if (present(nt)) then
                start(3) = nt
             endif
          elseif (dims == 2) then
             ! count set by dsize ok
             start(3) = k-lb2+1
             if (present(nt)) then
                start(4) = nt
             endif
          else
             write(iulog,*) trim(subname),' error dims incorrect ',clmlevel,dims
             call endrun()
          endif
          if (flag == 'write') data1d(:) = data(:,k)
          call ncd_iolocal_gs_${TYPE}1d(ncid, varname, flag, data1d, clmlevel, start, count, ier, lmissing, usepio=lusepio)
          if (flag == 'read') data(:,k) = data1d(:)
       enddo
       deallocate(data1d)
    else
       if (dims == 1) then
#if (defined SWITCH_DIMS)
          start(1) = 1
          count(1) = ub2-lb2+1
          count(2) = gsize
#else
          count(1) = gsize
          start(2) = 1
          count(2) = ub2-lb2+1
#endif
          if (present(nt)) then
             start(3) = nt
          endif
       elseif (dims == 2) then
          ! count set by dsize ok
          start(3) = 1
          count(3) = ub2-lb2+1
          if (present(nt)) then
             start(4) = nt
          endif
       else
          write(iulog,*) trim(subname),' error dims incorrect ',clmlevel,dims
          call endrun()
       endif
       call ncd_iolocal_gs_${TYPE}2d(ncid, varname, flag, data, clmlevel, start, count, ier, lmissing, usepio=lusepio)
    endif

    if (present(readvar)) then
       readvar = .false.
       if (ier == 0) readvar = .true.
    endif

  end subroutine ncd_iolocal_${TYPE}_2d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_gs_${TYPE}1d
!
! !INTERFACE:
  subroutine ncd_iolocal_gs_${TYPE}1d(ncid, varname, flag, data, clmlevel, start, count, status, missing, usepio)
!
! !DESCRIPTION:
! Netcdf i/o of 2d initial real field out to netCDF file
!
! !USES:
  use spmdGathScatMod, only : scatter_data_from_master, gather_data_to_master
!
! !ARGUMENTS:
    implicit none
    integer          ,intent(in)  :: ncid       ! input unit
    character(len=*) ,intent(in)  :: varname    ! variable name
    character(len=*) ,intent(in)  :: flag       ! 'read' or 'write'
    ${ATYPE},pointer              :: data(:)    ! local decomposition input data (out)
    character(len=*) ,intent(in)  :: clmlevel   ! type of grid
    integer, optional,intent(in)  :: start(:)   ! netcdf start index
    integer, optional,intent(in)  :: count(:)   ! netcdf count index
    integer, optional,intent(out) :: status     ! return code
    ${ATYPE},optional,intent(in)  :: missing    ! missing value
    logical, optional,intent(in)  :: usepio     ! use pio lib
    !--- rcodes:
    !      0  : success
    !    -99  : general error
    !     -5  : var not found on read
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer           :: n
    ${ATYPE}, pointer :: arrayg(:)
    integer           :: gsize      ! array global size from gsmap
    integer           :: lstart(4),lcount(4)  ! local start/count arrays
    integer           :: varid      ! varid
    integer           :: ndims      ! ndims for var
    integer           :: dids(4)    ! dim ids
    character(len=32) :: dname(4)   ! dim names
    integer           :: dlen(4)    ! dim lens
#if (defined BUILDPIO)
    type(pio_vardesc_plus_type),pointer  :: pio_vardesc_plus
    type(pio_iodesc_plus_type) ,pointer  :: pio_iodesc_plus
    integer           :: vdnum      ! vardesc num in list
    integer           :: basetype   ! pio initdecomp info
    integer           :: lenblocks  ! pio initdecomp info
    integer           :: dims(4)    ! pio initdecomp info
    integer           :: iodnum     ! iodesc num in list
    integer,pointer   :: compDOF(:)
    integer,pointer   :: ioDOF(:)
    integer(pio_offset),pointer :: pstart(:),pcount(:)
    
#endif
    logical           :: varpresent ! if true, variable is on tape
    integer           :: rcode      ! local return code
    integer           :: ier        ! error code
    integer :: data_offset              ! offset to single grid point for column model
    integer :: ndata                    ! count of pft's or columns to read
    logical :: lusepio           ! local usepio variable
    character(len=*),parameter :: subname='ncd_iolocal_gs_${TYPE}1d' ! subroutine name
!-----------------------------------------------------------------------

    lusepio = pio_def
    if (present(usepio)) then
       lusepio = usepio
    endif

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(flag),' ',trim(varname),' ',trim(clmlevel),lusepio
    endif

! tcx remove
! if (lusepio) then
!#if (defined BUILDPIO)
!    if (masterproc) write(iulog,*) trim(subname),' pio not implemented'
!!    call endrun('pio error')
!#endif
! else

   rcode = 0
   lstart = 1
   lcount = 1
   if (present(start).and.present(count)) then
      lstart(1:size(start)) = start(1:size(start))
      lcount(1:size(count)) = count(1:size(count))
   endif
   gsize = get_clmlevel_gsize(clmlevel)
   if (masterproc) then
      allocate(arrayg(gsize))
   endif

   if (flag == 'read') then
      if (masterproc) then
         call check_var(ncid, varname, varid, varpresent)
         if (varpresent) then
            if (single_column) then
               call scam_field_offsets(ncid,clmlevel,lstart,lcount)
               call check_ret(nf_get_vara_${NCTYPE}(ncid, varid, lstart, lcount, arrayg), subname)
            else
               if (present(start).and.present(count)) then
                  call check_ret(nf_get_vara_${NCTYPE}(ncid, varid, start, count, arrayg), subname)
               else
                  call check_ret(nf_get_var_${NCTYPE}(ncid, varid, arrayg), subname)
               endif
            endif
         else
            rcode = -5
         endif
      endif
      call scatter_data_from_master(data,arrayg,clmlevel)
   elseif (flag == 'write') then
      if (lusepio) then
#if (defined BUILDPIO) 
         call ncd_inqvdesc(varname,vdnum,subname,usepio=lusepio)
         if (vdnum < 1 .or. vdnum > pio_num_vardesc) then
            write(iulog,*) trim(subname),' ERROR in vdnum from inqvdesc ',trim(varname),vdnum
            call endrun()
         endif
         pio_vardesc_plus => pio_vardesc_list(vdnum)

         iodnum = pio_vardesc_plus%iodnum
         if (iodnum < 1 .or. iodnum > pio_num_iodesc) then
            write(iulog,*) trim(subname),' ERROR in iodnum from vardesc ',trim(varname),iodnum
            call endrun()
         endif
         pio_iodesc_plus => pio_iodesc_list(iodnum)

         !------------------------
         !--- setup iodesc if it's not set yet ---
         !------------------------
         if (.not. pio_iodesc_plus%set) then
            baseTYPE = ${MPTYPE}
            dims(:) = 1
            ndims = pio_iodesc_plus%ndims
            do n = 1,ndims
               call ncd_inqdlen(ncid,pio_iodesc_plus%dimids(n),dims(n),trim(subname),usepio=lusepio)
            enddo
            lenBLOCKS = 1

            call ncd_setDOF(clmlevel,dims,compDOF,ioDOF,pstart,pcount)

            !--- pio call ---
            call pio_initDecomp(pio_File,baseTYPE,dims,lenBLOCKS,compDOF,ioDOF,pstart,pcount,pio_iodesc_plus%pio_ioDesc)

            deallocate(compDOF)
            deallocate(IODOF)
            deallocate(pstart)
            deallocate(pcount)

            pio_iodesc_plus%set = .true.
         endif
         !------------------------
         !--- end setup iodesc ---
         !------------------------

         if (masterproc) then
            write(iulog,*) trim(subname),' write_darray1 ',vdnum,iodnum,pio_iodesc_plus%ndims,pio_iodesc_plus%dimids
            write(iulog,*) trim(subname),' write_darray2 ',vdnum,iodnum,dims,size(data)
            call shr_sys_flush(iulog)
         endif
         call mpi_barrier(mpicom,ier)

         call pio_setVarDesc(pio_iodesc_plus%pio_iodesc,pio_vardesc_plus%pio_vardesc)
         call PIO_write_darray(pio_File,pio_vardesc_plus%pio_varDesc,data,ier)
#endif 
      else
         if (present(missing)) then
            call gather_data_to_master(data,arrayg,clmlevel,missing)
         else
            call gather_data_to_master(data,arrayg,clmlevel)
         endif
         if (masterproc) then
            call check_ret(nf_inq_varid(ncid, varname, varid), subname)
            if (present(start).and.present(count)) then
               call check_ret(nf_put_vara_${NCTYPE}(ncid, varid, start, count, arrayg), subname)
            else
               call check_ret(nf_put_var_${NCTYPE}(ncid, varid, arrayg), subname)
            endif
         endif
      endif
   else
      if (masterproc) then
         write(iulog,*) subname,' error: unsupported flag ',trim(flag)
         call endrun()
      endif
   endif

   if (masterproc) then
      deallocate(arrayg)
   endif

   if (present(status)) then
      call mpi_bcast(rcode, 1, MPI_INTEGER, 0, mpicom, ier)
      status = rcode
   endif

!tcx remove
! endif

  end subroutine ncd_iolocal_gs_${TYPE}1d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_gs_${TYPE}2d
!
! !INTERFACE:
  subroutine ncd_iolocal_gs_${TYPE}2d(ncid, varname, flag, data, clmlevel, start, count, status, missing, usepio)
!
! !DESCRIPTION:
! Netcdf i/o of 2d initial real field out to netCDF file
!
! !USES:
  use spmdGathScatMod, only : scatter_data_from_master, gather_data_to_master
!
! !ARGUMENTS:
    implicit none
    integer          ,intent(in)  :: ncid       ! input unit
    character(len=*) ,intent(in)  :: varname    ! variable name
    character(len=*) ,intent(in)  :: flag       ! 'read' or 'write'
    ${ATYPE},pointer              :: data(:,:)  ! local decomposition input data (out)
    character(len=*) ,intent(in)  :: clmlevel   ! type of grid
    integer, optional,intent(in)  :: start(:)   ! netcdf start index
    integer, optional,intent(in)  :: count(:)   ! netcdf count index
    integer, optional,intent(out) :: status     ! return code
    ${ATYPE},optional,intent(in)  :: missing    ! missing value
    logical, optional,intent(in)  :: usepio     ! use pio lib
    !--- rcodes:
    !      0  : success
    !    -99  : general error
    !     -5  : var not found on read
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer           :: n
    ${ATYPE}, pointer :: arrayg(:,:)
    integer           :: gsize      ! array global size from gsmap
    integer           :: ksize      ! level ndims
    integer           :: lstart(4),lcount(4)  ! local start/count arrays
    logical           :: varpresent ! if true, variable is on tape
    integer           :: varid      ! varid
    integer           :: ndims      ! ndims for var
    integer           :: dids(4)    ! dim ids
    character(len=32) :: dname(4)   ! dim names
    integer           :: dlen(4)    ! dim lens
    integer           :: rcode      ! local return code
    integer           :: ier        ! error code
    integer :: data_offset              ! offset to single grid point for column model
    integer :: ndata                    ! count of pft's or columns to read
    logical :: lusepio           ! local usepio variable
    character(len=*),parameter :: subname='ncd_iolocal_gs_${TYPE}2d' ! subroutine name
!-----------------------------------------------------------------------

    lusepio = pio_def
    if (present(usepio)) then
       lusepio = usepio
    endif

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(flag),' ',trim(varname),' ',trim(clmlevel),lusepio
    endif

 if (lusepio) then
#if (defined BUILDPIO)
    if (masterproc) write(iulog,*) trim(subname),' pio not implemented'
!    call endrun('pio error')
#endif
 else

   rcode = 0
   lstart = 1
   lcount = 1
   if (present(start).and.present(count)) then
      lstart(1:size(start)) = start(1:size(start))
      lcount(1:size(count)) = count(1:size(count))
   endif
   gsize = get_clmlevel_gsize(clmlevel)
   ksize = size(data,dim=2)
   if (masterproc) then
      allocate(arrayg(gsize,ksize))
   endif

   if (flag == 'read') then
      if (masterproc) then
         call check_var(ncid, varname, varid, varpresent)
         if (varpresent) then
            if (single_column) then
               call scam_field_offsets(ncid,clmlevel,lstart,lcount)
               call check_ret(nf_get_vara_${NCTYPE}(ncid, varid, lstart, lcount, arrayg), subname)
            else
               if (present(start).and.present(count)) then
                  call check_ret(nf_get_vara_${NCTYPE}(ncid, varid, start, count, arrayg), subname)
               else
                  call check_ret(nf_get_var_${NCTYPE}(ncid, varid, arrayg), subname)
               endif
            endif
         else
            rcode = -5
         endif
      endif
      call scatter_data_from_master(data,arrayg,clmlevel)
   elseif (flag == 'write') then
      if (present(missing)) then
         call gather_data_to_master(data,arrayg,clmlevel,missing)
      else
         call gather_data_to_master(data,arrayg,clmlevel)
      endif
      if (masterproc) then
         call check_ret(nf_inq_varid(ncid, varname, varid), subname)
         if (present(start).and.present(count)) then
            call check_ret(nf_put_vara_${NCTYPE}(ncid, varid, start, count, arrayg), subname)
         else
            call check_ret(nf_put_var_${NCTYPE}(ncid, varid, arrayg), subname)
         endif
      endif
   else
      if (masterproc) then
         write(iulog,*) subname,' error: unsupported flag ',trim(flag)
         call endrun()
      endif
   endif

   if (masterproc) then
      deallocate(arrayg)
   endif

   if (present(status)) then
      call mpi_bcast(rcode, 1, MPI_INTEGER, 0, mpicom, ier)
      status = rcode
   endif

 endif

  end subroutine ncd_iolocal_gs_${TYPE}2d
EOF

end

exit


