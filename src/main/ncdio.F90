#include <misc.h>
#include <preproc.h>

module ncdio

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdioMod
!
! !DESCRIPTION:
! Generic interfaces to write fields to netcdf files
!
! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use spmdMod        , only : masterproc, mpicom, MPI_REAL8, MPI_INTEGER, &
                              MPI_LOGICAL
  use clm_varcon     , only : spval,ispval
  use shr_sys_mod    , only : shr_sys_flush
  use abortutils     , only : endrun
  use clm_varctl     , only : scmlon,scmlat, single_column
  use clm_mct_mod
  use spmdGathScatMod
!
! !PUBLIC TYPES:
  implicit none
  include 'netcdf.inc'
  save
  public :: check_ret   ! checks return status of netcdf calls
  public :: check_var   ! determine if variable is on netcdf file
  public :: check_dim   ! validity check on dimension
!
! !REVISION HISTORY:
!
!EOP
!
! !PRIVATE METHODS:
!
  interface ncd_iolocal
     module procedure ncd_iolocal_int_1d
     module procedure ncd_iolocal_real_1d
     module procedure ncd_iolocal_int_2d
     module procedure ncd_iolocal_real_2d
     module procedure ncd_iolocal_gs_real
     module procedure ncd_iolocal_gs_int
  end interface
  interface ncd_ioglobal
     module procedure ncd_ioglobal_int_var
     module procedure ncd_ioglobal_real_var
     module procedure ncd_ioglobal_int_1d
     module procedure ncd_ioglobal_real_1d
     module procedure ncd_ioglobal_int_2d
     module procedure ncd_ioglobal_real_2d
     module procedure ncd_ioglobal_int_3d
     module procedure ncd_ioglobal_real_3d
  end interface
  private :: get_size_dim1      ! obtain size of first dimension
  private :: scam_field_offsets ! get offset to proper lat/lon gridcell for SCAM
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_dim
!
! !INTERFACE:
  subroutine check_dim(ncid, dimname, value)
!
! !DESCRIPTION:
! Validity check on dimension
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: dimname
    integer, intent(in) :: value
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: dimid, dimlen    ! temporaries
!-----------------------------------------------------------------------

    call check_ret(nf_inq_dimid (ncid, trim(dimname), dimid), 'check_dim')
    call check_ret(nf_inq_dimlen (ncid, dimid, dimlen), 'check_dim')
    if (dimlen /= value) then
       write (6,*) 'CHECK_DIM error: mismatch of input dimension ',dimlen, &
            ' with expected value ',value,' for variable ',trim(dimname)
       call endrun()
    end if

  end subroutine check_dim

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_var
!
! !INTERFACE:
  subroutine check_var(ncid, varname, varid, readvar)
!
! !DESCRIPTION:
! Check if variable is on netcdf file
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)          :: ncid
    character(len=*), intent(in) :: varname
    integer, intent(out)         :: varid
    logical, intent(out)         :: readvar
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ret     ! return value
!-----------------------------------------------------------------------

    readvar = .true.
    if (masterproc) then
       ret = nf_inq_varid (ncid, varname, varid)
       if (ret/=NF_NOERR) then
          write(6,*)'CHECK_VAR: variable ',trim(varname),' is not on initial dataset'
#ifndef UNICOSMP
          call shr_sys_flush(6)
#endif
          readvar = .false.
       end if
    end if
  end subroutine check_var

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_ret
!
! !INTERFACE:
  subroutine check_ret(ret, calling)
!
! !DESCRIPTION:
! Check return status from netcdf call
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ret
    character(len=*) :: calling
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------

    if (ret /= NF_NOERR) then
       write(6,*)'netcdf error from ',trim(calling)
       write(6,*)'netcdf strerror = ',trim(NF_STRERROR(ret))
       call endrun()
    end if

  end subroutine check_ret

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_defvar
!
! !INTERFACE:
  subroutine ncd_defvar(ncid, varname, xtype, &
       dim1name, dim2name, dim3name, dim4name, dim5name, &
       long_name, units, cell_method, missing_value, fill_value, &
       imissing_value, ifill_value)
!
! !DESCRIPTION:
!  Define a netcdf variable
!
! !ARGUMENTS:
    implicit none
    integer         , intent(in)  :: ncid                    ! input unit
    character(len=*), intent(in)  :: varname                 ! variable name
    integer         , intent(in)  :: xtype                   ! external type
    character(len=*), intent(in), optional :: dim1name       ! dimension name
    character(len=*), intent(in), optional :: dim2name       ! dimension name
    character(len=*), intent(in), optional :: dim3name       ! dimension name
    character(len=*), intent(in), optional :: dim4name       ! dimension name
    character(len=*), intent(in), optional :: dim5name       ! dimension name
    character(len=*), intent(in), optional :: long_name      ! attribute
    character(len=*), intent(in), optional :: units          ! attribute
    character(len=*), intent(in), optional :: cell_method    ! attribute
    real(r8)        , intent(in), optional :: missing_value  ! attribute for real
    real(r8)        , intent(in), optional :: fill_value     ! attribute for real
    integer         , intent(in), optional :: imissing_value ! attribute for int
    integer         , intent(in), optional :: ifill_value    ! attribute for int
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: n              ! indices
    integer :: ndims          ! dimension counter
    integer :: dimid(5)       ! dimension ids
    integer :: varid          ! variable id
    integer :: itmp           ! temporary
    logical :: switchdim      ! true=> permute dim1 and dim2 for output
    character(len=256) :: str ! temporary
    character(len=32) :: subname='NCD_DEFVAR_REAL' ! subroutine name
!-----------------------------------------------------------------------

    if (.not. masterproc) return

    ! Determine dimension ids for variable

    dimid(:) = 0

    if (present(dim1name)) then
       call check_ret(nf_inq_dimid(ncid, dim1name, dimid(1)), subname)
    end if
    if (present(dim2name)) then
       call check_ret(nf_inq_dimid(ncid, dim2name, dimid(2)), subname)
    end if
    if (present(dim3name)) then
       call check_ret(nf_inq_dimid(ncid, dim3name, dimid(3)), subname)
    end if
    if (present(dim4name)) then
       call check_ret(nf_inq_dimid(ncid, dim4name, dimid(4)), subname)
    end if
    if (present(dim5name)) then
       call check_ret(nf_inq_dimid(ncid, dim5name, dimid(5)), subname)
    end if

    ! Permute dim1 and dim2 if necessary
    ! (If first dimension corresponds to a clmtype 1d type and
    ! if second dimension is a level dimension)

    if (present(dim1name) .and. present(dim2name)) then
       if (dim1name=='gridcell' .or. dim1name=='landunit' .or. &
           dim1name=='column'   .or. dim1name=='pft') then
          switchdim = .false.
          if (dim2name(1:3)=='lev' .or. dim2name(1:3)=='num') then
             switchdim = .true.
          end if
#if (defined CASA)   
          if (dim2name=='npools' .or. dim2name=='nlive') then
             switchdim = .true.
          end if
#endif
          if (switchdim) then
             itmp = dimid(2)
             dimid(2) = dimid(1)
             dimid(1) = itmp
          endif
       end if
    end if

    ! Define variable

    if (present(dim1name)) then
       ndims = 0
       do n = 1, size(dimid)
          if (dimid(n) /= 0) ndims = ndims + 1
       end do
       call check_ret(nf_def_var(ncid, trim(varname), xtype, ndims, dimid(1:ndims), varid), subname)
    else
       call check_ret(nf_def_var(ncid, varname, xtype, 0, 0, varid), subname)
    end if
    if (present(long_name)) then
       call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
    end if
    if (present(units)) then
       call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
    end if
    if (present(cell_method)) then
       str = 'time: ' // trim(cell_method)
       call check_ret(nf_put_att_text(ncid, varid, 'cell_method', len_trim(str), trim(str)), subname)
    end if
    if (present(fill_value)) then
       call check_ret(nf_put_att_double(ncid, varid, '_FillValue', xtype, 1, fill_value), subname)
    end if
    if (present(missing_value)) then
       call check_ret(nf_put_att_double(ncid, varid, 'missing_value', xtype, 1, missing_value), subname)
    end if
    if (present(ifill_value)) then
       call check_ret(nf_put_att_int(ncid, varid, '_FillValue', xtype, 1, ifill_value), subname)
    end if
    if (present(imissing_value)) then
       call check_ret(nf_put_att_int(ncid, varid, 'missing_value', xtype, 1, imissing_value), subname)
    end if

  end subroutine ncd_defvar

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_int_1d
!
! !INTERFACE:

  subroutine ncd_iolocal_int_1d(varname, data, dim1name, &
       flag, ncid, nlonxy, nlatxy, nt, readvar, imissing)
!
! !DESCRIPTION:
! I/O for 1d int field
!
! !USES:
  use decompMod, only : map_dc2sn, map_sn2dc, ldecomp
  use spmdGathScatMod, only : scatter_data_from_master, gather_data_to_master
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
    integer         , intent(in)  :: ncid               ! input unit
    character(len=*), intent(in)  :: varname            ! variable name
    integer         , pointer     :: data(:)            ! local decomposition data
    character(len=*), intent(in)  :: dim1name           ! dimension name
    integer         , optional, intent(in) :: nlonxy    ! 2d longitude size
    integer         , optional, intent(in) :: nlatxy    ! 2d latitude size
    integer         , optional, intent(in) :: nt        ! time sample index
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: imissing  ! value to set missing data to

! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,k,n,ixy,jxy          ! indices
    integer :: ndims                    ! dimension counter
    integer :: dimid(3)                 ! dimension ids
    integer :: varid                    ! variable id
    integer :: nsize                    ! size of global array
    integer :: ier                      ! error status
    integer :: start(3)                 ! starting indices for netcdf field
    integer :: count(3)                 ! count values for netcdf field
    integer, pointer :: iglobdc(:)      ! global decomposition initial data
    integer, pointer :: iglobsn(:)      ! global s->n initial data
    integer, pointer :: fldxy(:,:)      ! grid-average single-level field
    integer :: data_offset              ! offset to single grid point for column model
    integer :: ndata                    ! count of pft's or columns to read
    character(len=256):: str            ! temporary
    character(len=32) :: subname='NCD_IOLOCAL_INT_1D' ! subroutine name
    logical :: varpresent               ! if true, variable is on tape
!-----------------------------------------------------------------------

    ! Get size of dim1name

    if (masterproc) then
       nsize = get_size_dim1(dim1name)
       allocate (iglobdc(nsize), iglobsn(nsize), stat=ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' allocation error'; call endrun()
       end if
    end if

    ! Write field either as 1d field or as xy field

    if (flag == 'write') then

       call gather_data_to_master (data, iglobdc, clmlevel=dim1name)
       if (masterproc) then

          call check_ret(nf_inq_varid(ncid, varname, varid), subname)

          if (present(nlonxy) .and. present(nlatxy)) then

             ! Write xy field

             start(1) = 1;  count(1) = nlonxy
             start(2) = 1;  count(2) = nlatxy
             if (present(nt)) then
                start(3) = nt;  count(3) = 1
             end if
             allocate(fldxy(nlonxy,nlatxy), stat=ier)
             if (ier /= 0) then
                write(6,*)subname,' allocation error for fldxy'; call endrun()
             end if
             if (dim1name /= 'gridcell') then
                write(6,*)subname,' error: 1d clm output type must be ',&
                     'at gridcell level if 2d xy output is requested'; call endrun()
             end if
             if ( .not. present(imissing) )then
                fldxy(:,:) = ispval
             else
                fldxy(:,:) = imissing
             end if
!dir$ concurrent
!cdir nodep
             do k = 1,nsize
                ixy = ldecomp%gdc2i(k)
                jxy = ldecomp%gdc2j(k)
                fldxy(ixy,jxy) = iglobdc(k)
             end do
             call check_ret(nf_put_vara_int(ncid, varid, start, count, fldxy), subname)
             deallocate(fldxy)

          else

             ! Write 1d field

             start(1) = 1; count(1) = nsize
             if (present(nt)) then
                start(2) = nt; count(2) = 1
             end if
             call map_dc2sn(iglobdc, iglobsn, dim1name)
             call check_ret(nf_put_vara_int(ncid, varid, start, count, iglobsn), subname)

          end if

       end if   ! end of if-masterproc block

    else if (flag == 'read') then

      if (masterproc) then
        call check_var(ncid, varname, varid, varpresent)
        if (varpresent) then
        if (single_column) then
           call scam_field_offsets(ncid,dim1name,data_offset,ndata)
           start(1) = data_offset; count(1) = ndata
           call check_ret(nf_get_vara_int(ncid, varid, start, count, iglobsn), subname)
        else
           call check_ret(nf_get_var_int(ncid, varid, iglobsn), subname)
           call map_sn2dc(iglobsn, iglobdc, dim1name)
        end if
        end if
     end if

       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)subname,' error from mpi_bcast'; call endrun()
       end if
       if (varpresent) call scatter_data_from_master(data, iglobdc, clmlevel=dim1name)
       if (present(readvar)) readvar = varpresent

    end if

    if (masterproc) deallocate (iglobdc, iglobsn)

  end subroutine ncd_iolocal_int_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_real_1d
!
! !INTERFACE:
  subroutine ncd_iolocal_real_1d(varname, data, dim1name, &
       flag, ncid, nlonxy, nlatxy, nt, readvar)
!
! !DESCRIPTION:
! I/O for 1d int field
!
! !USES:
  use decompMod, only : map_dc2sn, map_sn2dc, ldecomp
#ifdef RTM
  use RunoffMod      , only : runoff
#endif
  use spmdGathScatMod, only : scatter_data_from_master, gather_data_to_master
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
    integer         , intent(in)  :: ncid               ! input unit
    character(len=*), intent(in)  :: varname            ! variable name
    real(r8)        , pointer     :: data(:)            ! local decomposition data
    character(len=*), intent(in)  :: dim1name           ! dimension name
    integer         , optional, intent(in) :: nlonxy    ! 2d longitude size
    integer         , optional, intent(in) :: nlatxy    ! 2d latitude size
    integer         , optional, intent(in) :: nt        ! time sample index
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,k,n,ixy,jxy          ! indices
    integer :: ndims                    ! dimension counter
    integer :: dimid(3)                 ! dimension ids
    integer :: varid                    ! variable id
    integer :: nsize                    ! size of global array
    integer :: ier                      ! error status
    integer :: start(3)                 ! starting indices for netcdf field
    integer :: count(3)                 ! count values for netcdf field
    integer :: data_offset              ! offset to single grid point for column model
    integer :: ndata                    ! count of pft's or columns to read
    real(r8), pointer :: rglobdc(:)     ! global decomposition initial data
    real(r8), pointer :: rglobsn(:)     ! global s->n initial data
    real(r8), pointer :: fldxy(:,:)     ! grid-average single-level field
    character(len=256):: str            ! temporary
    character(len=32) :: subname='NCD_IOLOCAL_REAL_1D' ! subroutine name
    logical :: varpresent               ! if true, variable is on tape
!-----------------------------------------------------------------------

    ! Get size of dim1name

    if (masterproc) then
       nsize = get_size_dim1(dim1name)
       allocate (rglobdc(nsize), rglobsn(nsize), stat=ier)
       if (ier /= 0) then
          write(6,*)subname,' allocation error'; call endrun()
       end if
    end if

    ! Write field either as 1d field as or xy field

    if (flag == 'write') then

       call gather_data_to_master (data, rglobdc, clmlevel=dim1name)
       if (masterproc) then

          ! Define variable if it has not already been defined on the tape

          call check_ret(nf_inq_varid(ncid, varname, varid), subname)

          if (present(nlonxy) .and. present(nlatxy)) then

             ! Write xy output

             start(1) = 1;  count(1) = nlonxy
             start(2) = 1;  count(2) = nlatxy
             if (present(nt)) then
                start(3) = nt;  count(3) = 1
             end if
             allocate(fldxy(nlonxy,nlatxy), stat=ier)
             if (ier /= 0) then
                write(6,*)subname,' error: allocation error for fldxy'; call endrun()
             end if
             fldxy(:,:) = spval
             select case (dim1name)
#ifdef RTM
             case('allrof')
!dir$ concurrent
!cdir nodep
                do k = 1,nsize
                      ixy = runoff%gdc2i(k)
                      jxy = runoff%gdc2j(k)
                      fldxy(ixy,jxy) = rglobdc(k)
                end do
             case('lndrof')
!dir$ concurrent
!cdir nodep
                do k = 1,nsize
                   if (runoff%mask(k) == 1) then
                      ixy = runoff%gdc2i(k)
                      jxy = runoff%gdc2j(k)
                      fldxy(ixy,jxy) = rglobdc(k)
                   endif
                end do
             case('ocnrof')
!dir$ concurrent
!cdir nodep
                do k = 1,nsize
                   if (runoff%mask(k) == 2) then
                      ixy = runoff%gdc2i(k)
                      jxy = runoff%gdc2j(k)
                      fldxy(ixy,jxy) = rglobdc(k)
                   endif
                end do
#endif
             case default
                if (dim1name /= 'gridcell') then
                   write(6,*)subname,' error: 1d clm output type must be ',&
                        'at gridcell level if 2d xy output is requested'; call endrun()
                end if
!dir$ concurrent
!cdir nodep
                do k = 1,nsize
                   ixy = ldecomp%gdc2i(k)
                   jxy = ldecomp%gdc2j(k)
                   fldxy(ixy,jxy) = rglobdc(k)
                end do
             end select
             call check_ret(nf_put_vara_double(ncid, varid, start, count, fldxy), subname)
             deallocate(fldxy)

          else

             ! Write one-dimensional output

             start(1) = 1; count(1) = nsize
             if (present(nt)) then
                start(2) = nt; count(2) = 1
             end if
             call map_dc2sn(rglobdc, rglobsn, dim1name)
             call check_ret(nf_put_vara_double(ncid, varid, start, count, rglobsn), subname)

          end if

       end if   ! end of if-masterproc block

    else if (flag == 'read') then

       if (masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) then
          if (single_column) then
             call scam_field_offsets(ncid,dim1name,data_offset,ndata)
             start(1) = data_offset; count(1) = ndata
             call check_ret(nf_get_vara_double(ncid, varid, start, count, rglobsn), subname)
          else
          if (present(nlonxy) .and. present(nlatxy)) then

             ! Write xy output

             start(1) = 1;  count(1) = nlonxy
             start(2) = 1;  count(2) = nlatxy
             if (present(nt)) then
                start(3) = nt;  count(3) = 1
             end if
             allocate(fldxy(nlonxy,nlatxy), stat=ier)
             if (ier /= 0) then
                write(6,*)subname,' error: allocation error for fldxy'; call endrun()
             end if
             call check_ret(nf_get_vara_double(ncid, varid, start, count, fldxy), subname)
             rglobdc(:) = 0._r8
             select case (dim1name)
#ifdef RTM
             case('allrof')
!dir$ concurrent
!cdir nodep
                do k = 1,nsize
                      ixy = runoff%gdc2i(k)
                      jxy = runoff%gdc2j(k)
                      rglobdc(k) = fldxy(ixy,jxy)
                end do
             case('lndrof')
!dir$ concurrent
!cdir nodep
                do k = 1,nsize
                   if (runoff%mask(k) == 1) then
                      ixy = runoff%gdc2i(k)
                      jxy = runoff%gdc2j(k)
                      rglobdc(k) = fldxy(ixy,jxy)
                   endif
                end do
             case('ocnrof')
!dir$ concurrent
!cdir nodep
                do k = 1,nsize
                   if (runoff%mask(k) == 2) then
                      ixy = runoff%gdc2i(k)
                      jxy = runoff%gdc2j(k)
                      rglobdc(k) = fldxy(ixy,jxy)
                   endif
                end do
#endif
             case default
                if (dim1name /= 'gridcell') then
                   write(6,*)subname,' error: 1d clm output type must be ',&
                        'at gridcell level if 2d xy output is requested'; call endrun()
                end if
!dir$ concurrent
!cdir nodep
                do k = 1,nsize
                   ixy = ldecomp%gdc2i(k)
                   jxy = ldecomp%gdc2j(k)
                   rglobdc(k) = fldxy(ixy,jxy)
                end do
             end select
             deallocate(fldxy)

          else
             call check_ret(nf_get_var_double(ncid, varid, rglobsn), subname)
             call map_sn2dc(rglobsn, rglobdc, dim1name)
          end if
          end if !  (single_column)
          endif
       end if
       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)subname,' error from mpi_bcast'; call endrun()
       end if
       if (varpresent) call scatter_data_from_master(data, rglobdc, clmlevel=dim1name)
       if (present(readvar)) readvar = varpresent
    end if

    if (masterproc) deallocate (rglobdc, rglobsn)

  end subroutine ncd_iolocal_real_1d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_int_2d
!
! !INTERFACE:
  subroutine ncd_iolocal_int_2d(varname, data, dim1name, dim2name, &
             lowerb2, upperb2, flag, ncid, nlonxy, nlatxy, nt, readvar)
!
! !DESCRIPTION:
! Netcdf i/o of 2d initial integer field out to netCDF file
!
! !USES:
  use decompMod, only : map_dc2sn, map_sn2dc, ldecomp
  use spmdGathScatMod, only : scatter_data_from_master, gather_data_to_master
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
    integer         , intent(in)  :: ncid               ! input unit
    character(len=*), intent(in)  :: varname            ! variable name
    integer         , pointer     :: data(:,:)          ! local decomposition input data
    character(len=*), intent(in)  :: dim1name           ! dimension 1 name
    character(len=*), intent(in)  :: dim2name           ! dimension 2 name
    integer         , optional, intent(in) :: nlonxy    ! 2d longitude size
    integer         , optional, intent(in) :: nlatxy    ! 2d latitude size
    integer         , optional, intent(in) :: nt        ! time sample index
    integer         , optional, intent(in) :: lowerb2,upperb2 ! lower and upper bounds of second dimension
    logical         , optional, intent(out):: readvar  ! true => variable is on initial dataset (read only)
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: j,k,ixy,jxy              ! indices
    integer :: ier                      ! error status
    integer :: ndims                    ! dimension counter
    integer :: start(4), count(4)       ! bounds for io
    integer :: dimid(4)                 ! dimension ids
    integer :: varid                    ! variable id
    integer :: lb1,ub1,lb2,ub2          ! bounds of data
    integer :: nsize                    ! size of global array second dimension
    integer, pointer :: datap(:,:)      ! permutted 2d data
    integer, pointer :: iglobdc(:,:)    ! global decomposition initial data
    integer, pointer :: iglobsn(:,:)    ! global s->n initial data
    integer, pointer :: fldxy(:,:,:)    ! temporary
    integer :: data_offset              ! offset to single grid point for column model
    integer :: ndata                    ! count of pft's or columns to read
    character(len=256):: str            ! temporary
    character(len=32) :: subname='NCD_IOLOCAL_INT_2D' ! subroutine name
    logical :: varpresent               ! if true, variable is on tape
!-----------------------------------------------------------------------

    ! Dynamic memory allocation (note that data is permutted)

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
    allocate (datap(lb2:ub2,lb1:ub1), stat=ier)
    if (ier /= 0) then
       write(6,*)subname,' allocation error'; call endrun()
    end if
    if (masterproc) then
       nsize = get_size_dim1(dim1name)
       allocate (iglobsn(lb2:ub2,nsize), iglobdc(lb2:ub2,nsize), stat=ier)
       if (ier /= 0) then
          write(6,*)subname,' allocation error'; call endrun()
       end if
    end if

    if (flag == 'write') then

       ! Permute 2d data for output and write out permuted data

       do j = lb2,ub2
!dir$ concurrent
!cdir nodep
          do k = lb1,ub1
             datap(j,k) = data(k,j)
          end do
       end do
       call gather_data_to_master (datap, iglobdc, clmlevel=dim1name)

       if (masterproc) then

          call check_ret(nf_inq_varid(ncid, varname, varid), subname)

          if (present(nlonxy) .and. present(nlatxy)) then

             start(1) = 1;  count(1) = nlonxy
             start(2) = 1;  count(2) = nlatxy
             start(3) = 1;  count(3) = ub2-lb2+1
             if (present(nt)) then
                start(4) = nt;  count(4) = 1
             end if
             allocate(fldxy(nlonxy,nlatxy,lb2:ub2), stat=ier)
             if (ier /= 0) then
                write (6,*)subname,' allocation error for fldxy'; call endrun()
             end if
             if (dim1name /= 'gridcell') then
                write(6,*)subname,' error: 1d clm output type must be ',&
                     'at gridcell level if 2d xy output is requested'; call endrun()
             end if
             fldxy(:,:,:) = ispval
             do j = lb2,ub2
!dir$ concurrent
!cdir nodep
                do k = 1,nsize
                   ixy = ldecomp%gdc2i(k)
                   jxy = ldecomp%gdc2j(k)
                   fldxy(ixy,jxy,j) = iglobdc(j,k)
                end do
             end do
             call check_ret(nf_put_vara_int(ncid, varid, start, count, fldxy), subname)
             deallocate(fldxy)

          else

             start(1) = 1;  count(1) = ub2-lb2+1
             start(2) = 1;  count(2) = nsize
             if (present(nt)) then
                start(3) = nt;  count(3) = 1
             end if
             call map_dc2sn(iglobdc, iglobsn, dim1name, lb2, ub2)
             call check_ret(nf_put_vara_int(ncid, varid, start, count, iglobsn), subname)

          end if

       end if   ! end of if-masterproc block

    else if (flag == 'read') then

       ! Determine if will read variable, read permutted variable data
       ! from netcdf file and unpermute the data

       if (masterproc) then
         call check_var(ncid, varname, varid, varpresent)
         if (varpresent) then
         if (single_column)then
            call scam_field_offsets(ncid,dim1name,data_offset,ndata)
            start(1) = 1;  count(1) = ub2-lb2+1
            start(2) = data_offset;  count(2) = ndata
            call check_ret(nf_get_vara_int(ncid, varid, start, count, iglobsn), subname)
         else
            call check_ret(nf_get_var_int(ncid, varid, iglobsn), subname)
            call map_sn2dc(iglobsn, iglobdc, dim1name, lb2, ub2)
         end if
         end if
       end if

       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast'; call endrun()
       end if
       if (varpresent) call scatter_data_from_master(datap, iglobdc, clmlevel=dim1name)
       if (varpresent) then
          do j = lb2,ub2
!dir$ concurrent
!cdir nodep
             do k = lb1,ub1
                data(k,j) = datap(j,k)
             end do
          end do
       end if
       if (present(readvar)) readvar = varpresent

    end if

    deallocate(datap)
    if (masterproc) deallocate (iglobdc, iglobsn)

  end subroutine ncd_iolocal_int_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_real_2d
!
! !INTERFACE:
  subroutine ncd_iolocal_real_2d(varname, data, dim1name, dim2name, &
             lowerb2, upperb2, flag, ncid, nlonxy, nlatxy, nt, readvar)
!
! !DESCRIPTION:
! Netcdf i/o of 2d initial integer field out to netCDF file
!
! !USES:
  use decompMod, only : map_dc2sn, map_sn2dc, ldecomp
  use spmdGathScatMod, only : scatter_data_from_master, gather_data_to_master
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
    integer         , intent(in)  :: ncid               ! input unit
    character(len=*), intent(in)  :: varname            ! variable name
    real(r8)        , pointer     :: data(:,:)          ! local decomposition input data
    character(len=*), intent(in)  :: dim1name           ! dimension 1 name
    character(len=*), intent(in)  :: dim2name           ! dimension 2 name
    integer         , optional, intent(in) :: nlonxy    ! 2d longitude size
    integer         , optional, intent(in) :: nlatxy    ! 2d latitude size
    integer         , optional, intent(in) :: nt        ! time sample index
    integer         , optional, intent(in) :: lowerb2,upperb2 ! lower and upper bounds of second dimension
    logical         , optional, intent(out):: readvar  ! true => variable is on initial dataset (read only)
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: j,k,ixy,jxy              ! indices
    integer :: ier                      ! error status
    integer :: ndims                    ! dimension counter
    integer :: start(4), count(4)       ! bounds for io
    integer :: dimid(4)                 ! dimension ids
    integer :: varid                    ! variable id
    integer :: lb1,ub1,lb2,ub2          ! bounds of data
    integer :: nsize                    ! size of global array second dimension         
    integer :: data_offset              ! offset to single grid point for column model
    integer :: ndata                    ! count of pft's or columns to read
    real(r8), pointer :: datap(:,:)     ! permutted 2d data
    real(r8), pointer :: rglobdc(:,:)   ! global decomposition initial data
    real(r8), pointer :: rglobsn(:,:)   ! global s->n initial data
    real(r8), pointer :: fldxy(:,:,:)   ! temporary
    character(len=256):: str            ! temporary
    logical :: varpresent               ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOLOCAL_REAL_2D' ! subroutine name
!-----------------------------------------------------------------------

    ! Dynamic memory allocation (note that data is permutted)

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
    allocate (datap(lb2:ub2,lb1:ub1), stat=ier)
    if (ier /= 0) then
       write(6,*)subname,' allocation error'; call endrun()
    end if
    if (masterproc) then
       nsize = get_size_dim1(dim1name)
       allocate (rglobsn(lb2:ub2,nsize), rglobdc(lb2:ub2,nsize), stat=ier)
       if (ier /= 0) then
          write(6,*)subname,' allocation error'; call endrun()
       end if
    end if

    if (flag == 'write') then

       ! Permute 2d data for output and write out permuted data

       do j = lb2,ub2
!dir$ concurrent
!cdir nodep
          do k = lb1,ub1
             datap(j,k) = data(k,j)
          end do
       end do
       call gather_data_to_master (datap, rglobdc, clmlevel=trim(dim1name))

       if (masterproc) then

          call check_ret(nf_inq_varid(ncid, varname, varid), subname)

          if (present(nlonxy) .and. present(nlatxy)) then

             start(1) = 1;  count(1) = nlonxy
             start(2) = 1;  count(2) = nlatxy
             start(3) = 1;  count(3) = ub2-lb2+1
             if (present(nt)) then
                start(4) = nt;  count(4) = 1
             end if
             allocate(fldxy(nlonxy,nlatxy,lb2:ub2), stat=ier)
             if (ier /= 0) then
                write (6,*)subname,' error : allocation error for fldxy'; call endrun()
             end if
             if (dim1name /= 'gridcell') then
                write(6,*)subname,' error: 1d clm output type must be ',&
                     'at gridcell level if 2d xy output is requested'; call endrun()
             end if
             fldxy(:,:,:) = spval
             do j = lb2,ub2
!dir$ concurrent
!cdir nodep
                do k = 1,nsize
                   ixy = ldecomp%gdc2i(k)
                   jxy = ldecomp%gdc2j(k)
                   fldxy(ixy,jxy,j) = rglobdc(j,k)
                end do
             end do
             call check_ret(nf_put_vara_double(ncid, varid, start, count, fldxy), subname)
             deallocate(fldxy)

          else

             start(1) = 1;  count(1) = ub2-lb2+1
             start(2) = 1;  count(2) = nsize
             if (present(nt)) then
                start(3) = nt;  count(3) = 1
             end if
             call map_dc2sn(rglobdc, rglobsn, dim1name, lb2, ub2)
             call check_ret(nf_put_vara_double(ncid, varid, start, count, rglobsn), subname)

          end if
       end if   ! end of if-masterproc block

    else if (flag == 'read') then

       ! Determine if will read variable, read permutted variable data
       ! from netcdf file and unpermute the data

       if (masterproc) then
         call check_var(ncid, varname, varid, varpresent)
         if (varpresent) then
         if (single_column)then
            call scam_field_offsets(ncid,dim1name,data_offset,ndata)
            start(1) = 1;  count(1) = ub2-lb2+1
            start(2) = data_offset;  count(2) = ndata
            call check_ret(nf_get_vara_double(ncid, varid, start, count, rglobsn), subname)
         else
            call check_ret(nf_get_var_double(ncid, varid, rglobsn), subname)
            call map_sn2dc(rglobsn, rglobdc, dim1name, lb2, ub2)
         end if
         end if
       end if

       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast'; call endrun()
       end if
       if (varpresent) call scatter_data_from_master(datap, rglobdc, clmlevel=dim1name)
       if (varpresent) then
          do j = lb2,ub2
!dir$ concurrent
!cdir nodep
             do k = lb1,ub1
                data(k,j) = datap(j,k)
             end do
          end do
       end if
       if (present(readvar)) readvar = varpresent

    end if

    deallocate(datap)
    if (masterproc) deallocate (rglobdc, rglobsn)

  end subroutine ncd_iolocal_real_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_gs_real
!
! !INTERFACE:
  subroutine ncd_iolocal_gs_real(ncid, varname, flag, data, beg, end, gsmap, perm, start, count)
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
    real(r8),pointer              :: data(:)    ! local decomposition input data (out)
    integer          ,intent(in)  :: beg        ! local start index
    integer          ,intent(in)  :: end        ! local end index
    type(mct_gsMap)  ,intent(in)  :: gsmap      ! gsmap associate with data decomp
    integer,pointer               :: perm(:)    ! permute array assoicated with gsmap
    integer, optional,intent(in)  :: start(:)   ! netcdf start index
    integer, optional,intent(in)  :: count(:)   ! netcdf count index

!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
  integer varid
  real(r8), pointer :: arrayg(:)
  integer           :: gsize      ! array global size from gsmap
  character(len=32) :: subname='NCD_IOGLOBAL_GS_REAL' ! subroutine name
!-----------------------------------------------------------------------

   gsize = mct_gsmap_gsize(gsmap)
   if (flag == 'read') then
      if (masterproc) then
         allocate(arrayg(gsize))
         call check_ret(nf_inq_varid(ncid, varname, varid), subname)
         if (present(start).and.present(count)) then
            call check_ret(nf_get_vara_double(ncid, varid, start, count, arrayg), subname)
         else
            call check_ret(nf_get_var_double(ncid, varid, arrayg), subname)
         endif
      endif
      call scatter_data_from_master(data,arrayg,gsMap,perm,beg,end)
      if (masterproc) then
         deallocate(arrayg)
      endif
   else
      if (masterproc) then
         write(6,*) subname,' error: unsupported flag ',trim(flag)
         call endrun()
      endif
   endif

  end subroutine ncd_iolocal_gs_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_gs_int
!
! !INTERFACE:
  subroutine ncd_iolocal_gs_int(ncid, varname, flag, data, beg, end, gsmap, perm, start, count)
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
    integer,pointer               :: data(:)    ! local decomposition input data
    integer          ,intent(in)  :: beg        ! local start index
    integer          ,intent(in)  :: end        ! local end index
    type(mct_gsMap)  ,intent(in)  :: gsmap      ! gsmap associate with data decomp
    integer, pointer              :: perm(:)    ! permute array assoicated with gsmap
    integer, optional,intent(in)  :: start(:)   ! netcdf start index
    integer, optional,intent(in)  :: count(:)   ! netcdf count index

!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
  integer varid
  integer, pointer  :: arrayg(:)
  integer           :: gsize      ! array global size from gsmap
  character(len=32) :: subname='NCD_IOGLOBAL_GS_INT' ! subroutine name
!-----------------------------------------------------------------------

   gsize = mct_gsmap_gsize(gsmap)
   if (flag == 'read') then
      if (masterproc) then
         allocate(arrayg(gsize))
         call check_ret(nf_inq_varid(ncid, varname, varid), subname)
         if (present(start).and.present(count)) then
            call check_ret(nf_get_vara_int(ncid, varid, start, count, arrayg), subname)
         else
            call check_ret(nf_get_var_int(ncid, varid, arrayg), subname)
         endif
      endif
      call scatter_data_from_master(data,arrayg,gsMap,perm,beg,end)
      if (masterproc) then
         deallocate(arrayg)
      endif
   else
      if (masterproc) then
         write(6,*) subname,' error: unsupported flag ',trim(flag)
         call endrun()
      endif
   endif

  end subroutine ncd_iolocal_gs_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_var
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_var(varname, data, flag, ncid, readvar, nt)
!
! !DESCRIPTION:
! I/O of integer variable
!

! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    integer         , intent(inout) :: data             ! local decomposition data
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                            ! error status
    integer :: dimid(1)                       ! dimension id
    integer :: start(1), count(1)             ! output bounds
    integer :: varid                          ! variable id
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_INT_VAR' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = nt; count(1) = 1
             call check_ret(nf_put_vara_int(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_int(ncid, varid, data), subname)
          end if
       end if

    else if (flag == 'read') then

       if (masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) call check_ret(nf_get_var_int(ncid, varid, data), subname)
       end if
       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast for varpresent'; call endrun()
       end if
       if (varpresent) then
          call mpi_bcast(data, 1, MPI_INTEGER, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*)trim(subname),' error from mpi_bcast for data'; call endrun()
          end if
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_int_var

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_var
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_var(varname, data, flag, ncid, readvar, nt)
!
! !DESCRIPTION:
! I/O of real variable
!

! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    real(r8)        , intent(inout) :: data             ! local decomposition data
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                            ! error status
    integer :: dimid(1)                       ! dimension id
    integer :: start(1), count(1)             ! output bounds
    integer :: varid                          ! variable id
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_REAL_VAR' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = nt; count(1) = 1
             call check_ret(nf_put_vara_double(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_double(ncid, varid, data), subname)
          end if
       end if

    else if (flag == 'read') then

       if (masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) call check_ret(nf_get_var_double(ncid, varid, data), subname)
       end if
       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast for varpresent'; call endrun()
       end if
       if (varpresent) then
          call mpi_bcast(data, 1, MPI_REAL8, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*)trim(subname),' error from mpi_bcast for data'; call endrun()
          end if
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_real_var

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_1d
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_1d(varname, data, flag, ncid, readvar, nt)
!
! !DESCRIPTION:
! Master I/O for 1d integer data
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    integer         , intent(inout) :: data(:)          ! local decomposition data
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: dimid(2), ndims                ! dimension ids
    integer :: start(2), count(2)             ! output bounds
    integer :: ier                            ! error code
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_INT_1D' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data)
             start(2) = nt; count(2) = 1
             call check_ret(nf_put_vara_int(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_int(ncid, varid, data), subname)
          end if
       end if

    else if (flag == 'read') then

       if (masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) call check_ret(nf_get_var_int(ncid, varid, data), subname)
       end if
       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast for varpresent'; call endrun()
       endif
       if (varpresent) then
          call mpi_bcast(data, size(data), MPI_INTEGER, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*)trim(subname),' error from mpi_bcast for data'; call endrun()
          end if
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_int_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_1d
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_1d(varname, data, flag, ncid, readvar, nt)
!
! !DESCRIPTION:
! Master I/O for 1d real data
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    real(r8)        , intent(inout) :: data(:)          ! local decomposition input data
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: ier                            ! error code
    integer :: dimid(2), ndims                ! dimension ids
    integer :: start(2), count(2)             ! output bounds
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_REAL_1D' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data)
             start(2) = nt; count(2) = 1
             call check_ret(nf_put_vara_double(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_double(ncid, varid, data), subname)
          end if
       end if

    else if (flag == 'read') then

       if (masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) call check_ret(nf_get_var_double(ncid, varid, data), subname)
       end if
       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast for varpresent'; call endrun()
       endif
       if (varpresent) then
          call mpi_bcast(data, size(data), MPI_REAL8, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*)trim(subname),' error from mpi_bcast for data'; call endrun()
          end if
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_real_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_2d
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_2d(varname, data, flag, ncid, readvar, nt)
!
! !DESCRIPTION:
! netcdf I/O of global 2d integer array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    integer         , intent(inout) :: data(:,:)        ! local decomposition input data
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: dimid(3), ndims                ! dimension ids
    integer :: start(3), count(3)             ! output bounds
    integer :: ier                            ! error code
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_2D_INT_IO' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = nt; count(3) = 1
             call check_ret(nf_put_vara_int(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_int(ncid, varid, data), subname)
          end if
       end if

    else if (flag == 'read') then

       if (masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) call check_ret(nf_get_var_int(ncid, varid, data), subname)
       end if
       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast for varpresent'; call endrun()
       endif
       if (varpresent) then
          call mpi_bcast(data, size(data), MPI_INTEGER, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*)trim(subname),' error from mpi_bcast for data'; call endrun()
          end if
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_int_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_2d
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_2d(varname, data, long_name, units, flag, &
                                  ncid, readvar, nt)
!
! !DESCRIPTION:
! netcdf I/O of global 2d real array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    real(r8)        , intent(inout) :: data(:,:)        ! local decomposition input data
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: ier                            ! error code
    integer :: dimid(3), ndims                ! dimension ids
    integer :: start(3), count(3)             ! output bounds
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_REAL_2D' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = nt; count(3) = 1
             call check_ret(nf_put_vara_double(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_double(ncid, varid, data), subname)
          end if
       end if

    else if (flag == 'read') then

       if (masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) call check_ret(nf_get_var_double(ncid, varid, data), subname)
       end if
       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast for varpresent'; call endrun()
       endif
       if (varpresent) then
          call mpi_bcast(data, size(data), MPI_REAL8, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*)trim(subname),' error from mpi_bcast for data'; call endrun()
          end if
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_real_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_3d
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_3d(varname, data, long_name, units, flag, &
                                 ncid, readvar, nt)
!
! !DESCRIPTION:
! netcdf I/O of global 3d integer array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    integer         , intent(inout) :: data(:,:,:)      ! local decomposition input data
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                    ! netCDF variable id
    integer :: dimid(4), ndims          ! dimension ids
    integer :: start(4), count(4)       ! output bounds
    integer :: ier                      ! error code
    logical :: varpresent               ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_3D_INT_IO' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = 1;  count(3) = size(data, dim=3)
             start(4) = nt; count(4) = 1
             call check_ret(nf_put_vara_int(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_int(ncid, varid, data), subname)
          end if
       end if

    else if (flag == 'read') then

       if (masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) call check_ret(nf_get_var_int(ncid, varid, data), subname)
       end if
       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast for varpresent'; call endrun()
       endif
       if (varpresent) then
          call mpi_bcast(data, size(data), MPI_INTEGER, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*)trim(subname),' error from mpi_bcast for data'; call endrun()
          end if
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_int_3d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_3d
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_3d(varname, data, long_name, units, flag, &
                                  ncid, readvar, nt)
!
! !DESCRIPTION:
! netcdf I/O of global 3d real array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    real(r8)        , intent(inout) :: data(:,:,:)      ! local decomposition input data
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                    ! netCDF variable id
    integer :: ier                      ! error code
    integer :: dimid(4), ndims          ! dimension ids
    integer :: start(4), count(4)       ! output bounds
    logical :: varpresent               ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_REAL_3D' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = 1;  count(3) = size(data, dim=3)
             start(4) = nt; count(4) = 1
             call check_ret(nf_put_vara_double(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_double(ncid, varid, data), subname)
          end if
       end if

    else if (flag == 'read') then

       if (masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) call check_ret(nf_get_var_double(ncid, varid, data), subname)
       end if
       call mpi_bcast(varpresent, 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= 0) then
          write(6,*)trim(subname),' error from mpi_bcast for varpresent'; call endrun()
       endif
       if (varpresent) then
          call mpi_bcast(data, size(data), MPI_REAL8, 0, mpicom, ier)
          if (ier /= 0) then
             write(6,*)trim(subname),' error from mpi_bcast for data'; call endrun()
          end if
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_real_3d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_size_dim1
!
! !INTERFACE:
  integer function get_size_dim1 (dim1name)
!
! !DESCRIPTION:
! Determine 1d size from dim1name
!
! !USES:
  use decompMod, only : get_proc_global
#ifdef RTM
  use RunoffMod, only : get_proc_rof_global
#endif
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: dim1name    !type of clm 1d array
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: nump        ! total number of pfts across all processors
    integer :: numc        ! total number of columns across all processors
    integer :: numl        ! total number of landunits across all processors
    integer :: numg        ! total number of gridcells across all processors
#ifdef RTM
    integer :: num_rtm     ! total number of all rtm cells on all procs
    integer :: num_lndrof  ! total number of land runoff across all procs
    integer :: num_ocnrof  ! total number of ocean runoff across all procs
#endif
!-----------------------------------------------------------------------
    ! Determine necessary indices

    call get_proc_global(numg, numl, numc, nump)
#ifdef RTM
    call get_proc_rof_global(num_rtm, num_lndrof, num_ocnrof)
#endif

    select case (dim1name)
    case('gridcell')
       get_size_dim1 = numg
    case('landunit')
       get_size_dim1 = numl
    case('column')
       get_size_dim1 = numc
    case('pft')
       get_size_dim1 = nump
#ifdef RTM
    case('allrof')
       get_size_dim1 = num_rtm
    case('lndrof')
       get_size_dim1 = num_rtm
    case('ocnrof')
       get_size_dim1 = num_rtm
#endif
    case default
       write(6,*) 'GET1DSIZE does not match dim1 type: ', trim(dim1name)
       call endrun()
    end select

  end function get_size_dim1
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine scam_field_offsets
!
! !INTERFACE:
  subroutine scam_field_offsets(ncid,dim1name,data_offset, ndata)
!
! !DESCRIPTION: 
! Read/Write initial data from/to netCDF instantaneous initial data file 
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use decompMod   , only : get_proc_bounds
    use clm_varpar  , only : maxpatch
    use nanMod      , only : nan
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: dim1name ! dimension 1 name
    integer, intent(in)  :: ncid             ! netCDF dataset id
    integer, intent(out) :: data_offset      ! offset into land array 
                                             ! 1st column 
    integer, intent(out) :: ndata            ! number of column (or 
                                             ! pft points to read)
!
! !CALLED FROM: subroutine inicfields
!
! !REVISION HISTORY:
! Created by John Truesdale

!EOP
!
! !LOCAL VARIABLES:
    real(r8) , pointer :: cols1dlon(:)       ! holds cols1d_ixy var
    real(r8) , pointer :: cols1dlat(:)       ! holds cols1d_jxy var
    real(r8) , pointer :: pfts1dlon(:)       ! holds pfts1d_ixy var
    real(r8) , pointer :: pfts1dlat(:)       ! holds pfts1d_jxy var
    integer cols(maxpatch)                   ! grid cell columns for scam
    integer pfts(maxpatch)                   ! grid cell pfts for scam
    integer :: cc,i                          ! index variable
    integer :: totpfts                       ! total number of pfts
    integer :: totcols                       ! total number of columns
    integer, save :: col_offset              ! offset into land array of 
                                             ! starting column
    integer, save :: pi_offset               ! offset into land array of 
                                             ! starting pft needed for scam
    integer :: dimid                         ! netCDF dimension id
    integer :: varid                         ! netCDF variable id
    integer, save :: begp, endp              ! per-proc beg/end pft indices
    integer, save :: begc, endc              ! per-proc beg/end col indices 
    integer, save :: begl, endl              ! per-proc beg/end land indices
    integer, save :: begg, endg              ! per-proc beg/end gridcell ind
    logical, save :: firsttime = .true.      ! determine offsets once.
    integer latidx,lonidx,ret
    real(r8) closelat,closelon

!------------------------------------------------------------------------
    if ( firsttime) then
       call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,latidx,lonidx)

       call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

       ret = nf_inq_dimid (ncid, 'column', dimid)
       if(ret==NF_NOERR) return
       if (ret/=NF_EBADDIM) print *, 'NETCDF ERROR: ', NF_STRERROR(ret)

       ret = nf_inq_dimlen (ncid, dimid, totcols)
       if (ret/=NF_NOERR) print *, 'NETCDF ERROR: ', NF_STRERROR(ret)

       ret = nf_inq_dimid (ncid, 'pft', dimid)
       if(ret==NF_NOERR) return
       if (ret/=NF_EBADDIM) print *, 'NETCDF ERROR: ', NF_STRERROR(ret)

       ret = nf_inq_dimlen (ncid, dimid, totpfts)
       if (ret/=NF_NOERR) print *, 'NETCDF ERROR: ', NF_STRERROR(ret)

       allocate (pfts1dlon(totpfts))
       allocate (pfts1dlat(totpfts))
       
       ret =  nf_inq_varid (ncid, 'pfts1d_ixy', varid)
       if (ret/=NF_NOERR) then
         write(6,*)'inq_varid: id for pfts1d_ixy not found'
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
       end if

       ret = nf_get_var_double (ncid, varid, pfts1dlon)
       if (ret/=NF_NOERR) then
         write(6,*)'GET_VAR_REALX: error reading pfts1dlon, varid =', varid
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
       end if

       ret =  nf_inq_varid (ncid, 'pfts1d_jxy', varid)
       if (ret/=NF_NOERR) then
         write(6,*)'inq_varid: id for pfts1d_jxy not found'
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
       end if

       ret = nf_get_var_double (ncid, varid, pfts1dlat)
       if (ret/=NF_NOERR) then
         write(6,*)'GET_VAR_REALX: error reading pfts1dlat, varid =', varid
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
       end if


       allocate (cols1dlon(totcols))
       allocate (cols1dlat(totcols))

       ret =  nf_inq_varid (ncid, 'cols1d_ixy', varid)
       if (ret/=NF_NOERR) then
         write(6,*)'inq_varid: id for cols1d_ixy not found'
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
       end if

       ret = nf_get_var_double (ncid, varid, cols1dlon)
       if (ret/=NF_NOERR) then
         write(6,*)'GET_VAR_REALX: error reading cols1dlon, varid =', varid
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
       end if

       ret =  nf_inq_varid (ncid, 'cols1d_jxy', varid)
       if (ret/=NF_NOERR) then
         write(6,*)'inq_varid: id for cols1d_jxy not found'
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
       end if

       ret = nf_get_var_double (ncid, varid, cols1dlat)
       if (ret/=NF_NOERR) then
         write(6,*)'GET_VAR_REALX: error reading cols1dlat, varid =', varid
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
       end if

       cols(:)=nan
       pfts(:)=nan
       col_offset=nan
       pi_offset=nan
       i=1
       do cc = 1, totcols
          if (cols1dlon(cc).eq.lonidx.and.cols1dlat(cc).eq.latidx) then
             cols(i)=cc
             i=i+1
          end if
       end do
       if (endc-begc+1.ne.i-1) then
          write(6,*)'error in number of columns read for this gridcell',endc,begc,i
          call endrun
       end if
       if (i.eq.1) then
          write(6,*)'couldnt find any columns for this latitude ',latidx,' and longitude ',lonidx
          call endrun
       else
          col_offset=cols(1)
       end if

       i=1
       do cc = 1, totpfts
          if (pfts1dlon(cc).eq.lonidx.and.pfts1dlat(cc).eq.latidx) then
             pfts(i)=cc
             i=i+1
          end if
       end do
       if (endp-begp+1.ne.i-1) then
          write(6,*)'error in number of pfts read for this gridcell',endp,begp,i
          call endrun
       end if
       if (i.eq.1) then
          write(6,*)'couldnt find any pfts for this latitude ',latidx,' and longitude ',lonidx
          call endrun
       else
          pi_offset=pfts(1)
       end if

       deallocate (pfts1dlon)
       deallocate (pfts1dlat)
       deallocate (cols1dlon)
       deallocate (cols1dlat)
       firsttime = .false.
    endif

    
    if (dim1name == 'pft') then
       data_offset = pi_offset
       ndata = endp-begp+1
    else if (dim1name == 'column') then
       data_offset = col_offset
       ndata = endc-begc+1
    else
       write(6,*)'error calculation array offsets for SCAM'; call endrun()
    endif
  end subroutine scam_field_offsets

end module ncdio
