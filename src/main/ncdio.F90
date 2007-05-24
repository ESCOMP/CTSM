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
  use clmtype        , only : gratm, grlnd, nameg, namel, namec, namep, allrof
  use clm_varcon     , only : spval,ispval
  use shr_sys_mod    , only : shr_sys_flush
  use abortutils     , only : endrun
  use clm_varctl     , only : single_column
  use clm_mct_mod
  use spmdGathScatMod
  use decompMod      , only : get_clmlevel_gsize
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
     module procedure ncd_iolocal_gs_real1d
     module procedure ncd_iolocal_gs_int1d
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
  private :: scam_field_offsets ! get offset to proper lat/lon gridcell for SCAM
  logical,parameter,private :: lbcast_def = .false.  ! lbcast default
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
          write(6,*)'CHECK_VAR: variable ',trim(varname),' is not on dataset'
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
       write(6,*)'netcdf error from ',trim(calling),':',trim(NF_STRERROR(ret))
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
    character(len=*),parameter :: subname='ncd_defvar_real' ! subroutine name
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
    integer         , optional, intent(in) :: imissing  ! missing value
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: gsize                    ! size of global array
    integer :: ier                      ! error status
    integer :: start(3)                 ! starting indices for netcdf field
    integer :: count(3)                 ! count values for netcdf field
    integer :: lmissing                 ! local missing value
    character(len=8) :: clmlevel        ! clmlevel
    character(len=*),parameter :: subname='ncd_iolocal_int_1d' ! subroutine name
!-----------------------------------------------------------------------

    if ((present(nlonxy) .and. .not.present(nlatxy)) .or. &
        (present(nlatxy) .and. .not.present(nlonxy))) then
       write(6,*) trim(subname),' error nlonxy/nlatxy must be both or neither present '
       call endrun()
    endif

    if (present(imissing)) then
       lmissing = imissing
    else
       lmissing = ispval
    endif

    clmlevel = dim1name
    if (present(nlonxy) .and. present(nlatxy)) then
       if (dim1name == nameg .or. dim1name == grlnd) then
          clmlevel = grlnd
       elseif (dim1name == allrof .or. dim1name == gratm) then
          ! continue, acceptable and default behavior for now
       else
          if (masterproc) write(6,*) trim(subname),' warning use dim1name and nlonxy/nlatxy ',trim(dim1name),nlonxy,nlatxy
       endif
    endif

    gsize = get_clmlevel_gsize(clmlevel)

    start = 1
    count = 1
    if (present(nlonxy) .and. present(nlatxy)) then
       count(1) = nlonxy
       count(2) = nlatxy
       if (present(nt)) then
          start(3) = nt
       endif
    else
       count(1) = gsize
       if (present(nt)) then
          start(2) = nt
       endif
    endif

    call ncd_iolocal_gs_int1d(ncid, varname, flag, data, clmlevel, start, count, ier, lmissing)

    if (present(readvar)) then
       readvar = .false.
       if (ier == 0) readvar = .true.
    endif

  end subroutine ncd_iolocal_int_1d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_real_1d
!
! !INTERFACE:
  subroutine ncd_iolocal_real_1d(varname, data, dim1name, &
       flag, ncid, nlonxy, nlatxy, nt, readvar, missing)
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
    real(r8)        , pointer     :: data(:)            ! local decomposition data
    character(len=*), intent(in)  :: dim1name           ! dimension name
    integer         , optional, intent(in) :: nlonxy    ! 2d longitude size
    integer         , optional, intent(in) :: nlatxy    ! 2d latitude size
    integer         , optional, intent(in) :: nt        ! time sample index
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    real(r8)        , optional, intent(in) :: missing   ! missing value
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: gsize                    ! size of global array
    integer :: ier                      ! error status
    integer :: start(3)                 ! starting indices for netcdf field
    integer :: count(3)                 ! count values for netcdf field
    real(r8):: lmissing                 ! local missing value
    character(len=8) :: clmlevel        ! clmlevel
    character(len=*),parameter :: subname='ncd_iolocal_real_1d' ! subroutine name
!-----------------------------------------------------------------------

    if ((present(nlonxy) .and. .not.present(nlatxy)) .or. &
        (present(nlatxy) .and. .not.present(nlonxy))) then
       write(6,*) trim(subname),' error nlonxy/nlatxy must be both or neither present '
       call endrun()
    endif

    if (present(missing)) then
       lmissing = missing
    else
       lmissing = spval
    endif

    clmlevel = dim1name
    if (present(nlonxy) .and. present(nlatxy)) then
       if (dim1name == nameg .or. dim1name == grlnd) then
          clmlevel = grlnd
       elseif (dim1name == allrof .or. dim1name == gratm) then
          ! continue, acceptable and default behavior for now
       else
          write(6,*) trim(subname),' warning use dim1name and nlonxy/nlatxy ',trim(dim1name),nlonxy,nlatxy
       endif
    endif

    gsize = get_clmlevel_gsize(clmlevel)

    start = 1
    count = 1
    if (present(nlonxy) .and. present(nlatxy)) then
       count(1) = nlonxy
       count(2) = nlatxy
       if (present(nt)) then
          start(3) = nt
       endif
    else
       count(1) = gsize
       if (present(nt)) then
          start(2) = nt
       endif
    endif

    call ncd_iolocal_gs_real1d(ncid, varname, flag, data, clmlevel, start, count, ier, lmissing)

    if (present(readvar)) then
       readvar = .false.
       if (ier == 0) readvar = .true.
    endif

  end subroutine ncd_iolocal_real_1d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_int_2d
!
! !INTERFACE:
  subroutine ncd_iolocal_int_2d(varname, data, dim1name, dim2name, &
             lowerb2, upperb2, flag, ncid, nlonxy, nlatxy, nt, readvar, imissing)
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
    integer         , pointer     :: data(:,:)          ! local decomposition input data
    character(len=*), intent(in)  :: dim1name           ! dimension 1 name
    character(len=*), intent(in)  :: dim2name           ! dimension 2 name
    integer         , optional, intent(in) :: nlonxy    ! 2d longitude size
    integer         , optional, intent(in) :: nlatxy    ! 2d latitude size
    integer         , optional, intent(in) :: nt        ! time sample index
    integer         , optional, intent(in) :: lowerb2,upperb2 ! lower and upper bounds of second dimension
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: imissing  ! missing value
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: k                        ! index
    integer :: gsize                    ! size of global array
    integer :: ier                      ! error status
    integer :: start(4)                 ! starting indices for netcdf field
    integer :: count(4)                 ! count values for netcdf field
    integer :: lb1,ub1                  ! lower/upper bound of dim 1
    integer :: lb2,ub2                  ! lower/upper bound of dim 2
    integer :: lmissing                 ! local missing value
    integer,pointer  :: data1d(:)       ! 1 level data
    character(len=8) :: clmlevel        ! clmlevel
    character(len=*),parameter :: subname='ncd_iolocal_int_2d' ! subroutine name
!-----------------------------------------------------------------------

    if ((present(nlonxy) .and. .not.present(nlatxy)) .or. &
        (present(nlatxy) .and. .not.present(nlonxy))) then
       write(6,*) trim(subname),' error nlonxy/nlatxy must be both or neither present '
       call endrun()
    endif

    lb1 = lbound(data, dim=1)
    ub1 = ubound(data, dim=1)
    allocate(data1d(lb1:ub1))

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

    if (present(imissing)) then
       lmissing = imissing
    else
       lmissing = ispval
    endif

    clmlevel = dim1name
    if (present(nlonxy) .and. present(nlatxy)) then
       if (dim1name == nameg .or. dim1name == grlnd) then
          clmlevel = grlnd
       else
          write(6,*) trim(subname),' error in dim1name and nlonxy/nlatxy ',trim(dim1name),nlonxy,nlatxy
          call endrun()
       endif
    endif

    gsize = get_clmlevel_gsize(clmlevel)

    do k = lb2,ub2
       start = 1
       count = 1
       if (present(nlonxy) .and. present(nlatxy)) then
          count(1) = nlonxy
          count(2) = nlatxy
          start(3) = k-lb2+1
          if (present(nt)) then
             start(4) = nt
          endif
       else
          start(1) = k-lb2+1
          count(2) = gsize
          if (present(nt)) then
             start(3) = nt
          endif
       endif
       if (flag == 'write') data1d(:) = data(:,k)
       call ncd_iolocal_gs_int1d(ncid, varname, flag, data1d, clmlevel, start, count, ier, lmissing)
       if (flag == 'read') data(:,k) = data1d(:)
    enddo

    deallocate(data1d)

    if (present(readvar)) then
       readvar = .false.
       if (ier == 0) readvar = .true.
    endif

  end subroutine ncd_iolocal_int_2d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_real_2d
!
! !INTERFACE:
  subroutine ncd_iolocal_real_2d(varname, data, dim1name, dim2name, &
             lowerb2, upperb2, flag, ncid, nlonxy, nlatxy, nt, readvar, missing)
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
    real(r8)        , pointer     :: data(:,:)          ! local decomposition input data
    character(len=*), intent(in)  :: dim1name           ! dimension 1 name
    character(len=*), intent(in)  :: dim2name           ! dimension 2 name
    integer         , optional, intent(in) :: nlonxy    ! 2d longitude size
    integer         , optional, intent(in) :: nlatxy    ! 2d latitude size
    integer         , optional, intent(in) :: nt        ! time sample index
    integer         , optional, intent(in) :: lowerb2,upperb2 ! lower and upper bounds of second dimension
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    real(r8)        , optional, intent(in) :: missing   ! missing value
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: k                        ! index
    integer :: gsize                    ! size of global array
    integer :: ier                      ! error status
    integer :: start(4)                 ! starting indices for netcdf field
    integer :: count(4)                 ! count values for netcdf field
    integer :: lb1,ub1                  ! lower/upper bound of dim 1
    integer :: lb2,ub2                  ! lower/upper bound of dim 2
    real(r8):: lmissing                 ! local missing value
    real(r8),pointer :: data1d(:)       ! 1 level data
    character(len=8) :: clmlevel        ! clmlevel
    character(len=*),parameter :: subname='ncd_iolocal_real_2d' ! subroutine name
!-----------------------------------------------------------------------

    if ((present(nlonxy) .and. .not.present(nlatxy)) .or. &
        (present(nlatxy) .and. .not.present(nlonxy))) then
       write(6,*) trim(subname),' error nlonxy/nlatxy must be both or neither present '
       call endrun()
    endif

    lb1 = lbound(data, dim=1)
    ub1 = ubound(data, dim=1)
    allocate(data1d(lb1:ub1))

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
       lmissing = spval
    endif

    clmlevel = dim1name
    if (present(nlonxy) .and. present(nlatxy)) then
       if (dim1name == nameg .or. dim1name == grlnd) then
          clmlevel = grlnd
       else
          write(6,*) trim(subname),' error in dim1name and nlonxy/nlatxy ',trim(dim1name),nlonxy,nlatxy
          call endrun()
       endif
    endif

    gsize = get_clmlevel_gsize(clmlevel)

    do k = lb2,ub2
       start = 1
       count = 1
       if (present(nlonxy) .and. present(nlatxy)) then
          count(1) = nlonxy
          count(2) = nlatxy
          start(3) = k-lb2+1
          if (present(nt)) then
             start(4) = nt
          endif
       else
          start(1) = k-lb2+1
          count(2) = gsize
          if (present(nt)) then
             start(3) = nt
          endif
       endif
       if (flag == 'write') data1d(:) = data(:,k)
       call ncd_iolocal_gs_real1d(ncid, varname, flag, data1d, clmlevel, start, count, ier, lmissing)
       if (flag == 'read') data(:,k) = data1d(:)
    enddo

    deallocate(data1d)

    if (present(readvar)) then
       readvar = .false.
       if (ier == 0) readvar = .true.
    endif

  end subroutine ncd_iolocal_real_2d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_gs_real1d
!
! !INTERFACE:
  subroutine ncd_iolocal_gs_real1d(ncid, varname, flag, data, clmlevel, start, count, status, missing)
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
    character(len=*) ,intent(in)  :: clmlevel   ! type of grid
    integer, optional,intent(in)  :: start(:)   ! netcdf start index
    integer, optional,intent(in)  :: count(:)   ! netcdf count index
    integer, optional,intent(out) :: status    ! return code
    real(r8),optional,intent(in)  :: missing    ! missing value
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
  integer varid
  real(r8), pointer :: arrayg(:)
  integer           :: gsize      ! array global size from gsmap
  integer           :: lstart(4),lcount(4)  ! local start/count arrays
  logical           :: varpresent ! if true, variable is on tape
  integer           :: rcode      ! local return code
  integer           :: ier        ! error code
  integer :: data_offset              ! offset to single grid point for column model
  integer :: ndata                    ! count of pft's or columns to read
  character(len=*),parameter :: subname='ncd_iolocal_gs_real1d' ! subroutine name
!-----------------------------------------------------------------------

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
!tcx
!               call scam_field_offsets(ncid,clmlevel,data_offset,ndata)
!               lstart(1) = data_offset; lcount(1) = ndata
!               call check_ret(nf_get_vara_double(ncid, varid, lstart, lcount, arrayg), subname)
               call scam_field_offsets(ncid,clmlevel,lstart,lcount)
               call check_ret(nf_get_vara_double(ncid, varid, lstart, lcount, arrayg), subname)
            else
               if (present(start).and.present(count)) then
                  call check_ret(nf_get_vara_double(ncid, varid, start, count, arrayg), subname)
               else
                  call check_ret(nf_get_var_double(ncid, varid, arrayg), subname)
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
            call check_ret(nf_put_vara_double(ncid, varid, start, count, arrayg), subname)
         else
            call check_ret(nf_put_var_double(ncid, varid, arrayg), subname)
         endif
      endif
   else
      if (masterproc) then
         write(6,*) subname,' error: unsupported flag ',trim(flag)
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

  end subroutine ncd_iolocal_gs_real1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_gs_int1d
!
! !INTERFACE:
  subroutine ncd_iolocal_gs_int1d(ncid, varname, flag, data, clmlevel, start, count, status, imissing)
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
    character(len=*) ,intent(in)  :: clmlevel   ! type of grid
    integer,pointer               :: data(:)    ! local decomposition input data
    integer, optional,intent(in)  :: start(:)   ! netcdf start index
    integer, optional,intent(in)  :: count(:)   ! netcdf count index
    integer, optional,intent(out) :: status    ! return code
    integer, optional,intent(in)  :: imissing   ! missing value
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
  integer varid
  integer, pointer  :: arrayg(:)
  integer           :: gsize      ! array global size from gsmap
  integer           :: lstart(4),lcount(4)  ! local start/count arrays
  logical           :: varpresent ! if true, variable is on tape
  integer           :: rcode      ! local return code
  integer           :: ier        ! error code
  integer           :: data_offset! offset to single grid point for column model
  integer           :: ndata      ! count of pft's or columns to read
  character(len=*),parameter :: subname='ncd_iolocal_gs_int1d' ! subroutine name
!-----------------------------------------------------------------------

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
!tcx
!               call scam_field_offsets(ncid,clmlevel,data_offset,ndata)
!               lstart(1) = data_offset; lcount(1) = ndata
!               call check_ret(nf_get_vara_int(ncid, varid, lstart, lcount, arrayg), subname)
               call scam_field_offsets(ncid,clmlevel,lstart,lcount)
               call check_ret(nf_get_vara_int(ncid, varid, lstart, lcount, arrayg), subname)
            else
               if (present(start).and.present(count)) then
                  call check_ret(nf_get_vara_int(ncid, varid, start, count, arrayg), subname)
               else
                  call check_ret(nf_get_var_int(ncid, varid, arrayg), subname)
               endif
            endif
         else
            rcode = -5
         endif
      endif
      call scatter_data_from_master(data,arrayg,clmlevel)
   elseif (flag == 'write') then
      if (present(imissing)) then
         call gather_data_to_master(data,arrayg,clmlevel,imissing)
      else
         call gather_data_to_master(data,arrayg,clmlevel)
      endif
      if (masterproc) then
         call check_ret(nf_inq_varid(ncid, varname, varid), subname)
         if (present(start).and.present(count)) then
            call check_ret(nf_put_vara_int(ncid, varid, start, count, arrayg), subname)
         else
            call check_ret(nf_put_var_int(ncid, varid, arrayg), subname)
         endif
      endif
   else
      if (masterproc) then
         write(6,*) subname,' error: unsupported flag ',trim(flag)
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

  end subroutine ncd_iolocal_gs_int1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_var
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_var(varname, data, flag, ncid, readvar, nt, bcast)
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
    logical         , optional, intent(in) :: bcast     ! bcast on read?
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                            ! error status
    integer :: start(4), count(4)             ! output bounds
    integer :: varid                          ! variable id
    logical :: varpresent                     ! if true, variable is on tape
    logical :: lbcast                   ! local copy of bcast flag
    character(len=*),parameter :: subname='ncd_ioglobal_int_var' ! subroutine name
!-----------------------------------------------------------------------

    start = 1
    count = 1
    lbcast = lbcast_def
    if (present(bcast)) then
       lbcast = bcast
    endif

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
          if (varpresent) then
             if (single_column) then
                call scam_field_offsets(ncid,'undefined',start,count)
                call check_ret(nf_get_vara_int(ncid, varid, start, count, data), subname)
             else
                call check_ret(nf_get_var_int(ncid, varid, data), subname)
             endif
          endif
       end if
       if (lbcast) then
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
       endif
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_int_var

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_var
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_var(varname, data, flag, ncid, readvar, nt, bcast)
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
    logical         , optional, intent(in) :: bcast     ! bcast on read?
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                            ! error status
    integer :: start(4), count(4)             ! output bounds
    integer :: varid                          ! variable id
    logical :: varpresent                     ! if true, variable is on tape
    logical :: lbcast                   ! local copy of bcast flag
    character(len=*),parameter :: subname='ncd_ioglobal_real_var' ! subroutine name
!-----------------------------------------------------------------------

    start = 1
    count = 1
    lbcast = lbcast_def
    if (present(bcast)) then
       lbcast = bcast
    endif

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
          if (varpresent) then
             if (single_column) then
                call scam_field_offsets(ncid,'undefined',start,count)
                call check_ret(nf_get_vara_double(ncid, varid, start, count, data), subname)
             else
                call check_ret(nf_get_var_double(ncid, varid, data), subname)
             endif
          endif
       end if
       if (lbcast) then
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
       endif
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_real_var

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_1d
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_1d(varname, data, flag, ncid, readvar, nt, bcast)
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
    logical         , optional, intent(in) :: bcast     ! bcast on read?
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: start(4), count(4)             ! output bounds
    integer :: ier                            ! error code
    logical :: varpresent                     ! if true, variable is on tape
    logical :: lbcast                   ! local copy of bcast flag
    character(len=*),parameter :: subname='ncd_ioglobal_int_1d' ! subroutine name
!-----------------------------------------------------------------------

    start = 1
    count = 1
    lbcast = lbcast_def
    if (present(bcast)) then
       lbcast = bcast
    endif

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
          if (varpresent) then
             if (single_column) then
                call scam_field_offsets(ncid,'undefined',start,count)
                call check_ret(nf_get_vara_int(ncid, varid, start, count, data), subname)
             else
                call check_ret(nf_get_var_int(ncid, varid, data), subname)
             endif
          endif
       end if
       if (lbcast) then
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
       endif
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_int_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_1d
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_1d(varname, data, flag, ncid, readvar, nt, bcast)
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
    logical         , optional, intent(in) :: bcast     ! bcast on read?
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: ier                            ! error code
    integer :: start(4), count(4)             ! output bounds
    logical :: varpresent                     ! if true, variable is on tape
    logical :: lbcast                   ! local copy of bcast flag
    character(len=*),parameter :: subname='ncd_ioglobal_real_1d' ! subroutine name
!-----------------------------------------------------------------------

    start = 1
    count = 1
    lbcast = lbcast_def
    if (present(bcast)) then
       lbcast = bcast
    endif

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
          if (varpresent) then
             if (single_column) then
                call scam_field_offsets(ncid,'undefined',start,count)
                call check_ret(nf_get_vara_double(ncid, varid, start, count, data), subname)
             else
                call check_ret(nf_get_var_double(ncid, varid, data), subname)
             endif
          endif
       end if
       if (lbcast) then
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
       endif
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_real_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_2d
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_2d(varname, data, flag, ncid, readvar, nt, bcast)
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
    logical         , optional, intent(in) :: bcast     ! bcast on read?
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: start(4), count(4)             ! output bounds
    integer :: ier                            ! error code
    logical :: varpresent                     ! if true, variable is on tape
    logical :: lbcast                   ! local copy of bcast flag
    character(len=*),parameter :: subname='ncd_ioglobal_2d_int_io' ! subroutine name
!-----------------------------------------------------------------------

    start = 1
    count = 1
    lbcast = lbcast_def
    if (present(bcast)) then
       lbcast = bcast
    endif

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
          if (single_column) then
             call scam_field_offsets(ncid,'undefined',start,count)
             call check_ret(nf_get_vara_int(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_get_var_int(ncid, varid, data), subname)
          endif
       end if
       if (lbcast) then
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
       endif
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_int_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_2d
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_2d(varname, data, flag, ncid, readvar, nt, bcast)
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
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
    logical         , optional, intent(in) :: bcast     ! bcast on read?
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: ier                            ! error code
    integer :: start(4), count(4)             ! output bounds
    logical :: varpresent                     ! if true, variable is on tape
    logical :: lbcast                   ! local copy of bcast flag
    character(len=*),parameter :: subname='ncd_ioglobal_real_2d' ! subroutine name
!-----------------------------------------------------------------------

    start = 1
    count = 1
    lbcast = lbcast_def
    if (present(bcast)) then
       lbcast = bcast
    endif

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
          if (single_column) then
             call scam_field_offsets(ncid,'undefined',start,count)
             call check_ret(nf_get_vara_double(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_get_var_double(ncid, varid, data), subname)
          endif
       end if
       if (lbcast) then
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
       endif
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_real_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_3d
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_3d(varname, data, flag, ncid, readvar, nt, bcast)
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
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
    logical         , optional, intent(in) :: bcast     ! bcast on read?
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                    ! netCDF variable id
    integer :: start(4), count(4)       ! output bounds
    integer :: ier                      ! error code
    logical :: varpresent               ! if true, variable is on tape
    logical :: lbcast                   ! local copy of bcast flag
    character(len=*),parameter :: subname='ncd_ioglobal_3d_int_io' ! subroutine name
!-----------------------------------------------------------------------

    start = 1
    count = 1
    lbcast = lbcast_def
    if (present(bcast)) then
       lbcast = bcast
    endif

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
          if (varpresent) then
             if (single_column) then
                call scam_field_offsets(ncid,'undefined',start,count)
                call check_ret(nf_get_vara_int(ncid, varid, start, count, data), subname)
             else
                call check_ret(nf_get_var_int(ncid, varid, data), subname)
             endif
          endif
       end if
       if (lbcast) then
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
       endif
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_int_3d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_3d
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_3d(varname, data, flag, ncid, readvar, nt, bcast)
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
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
    logical         , optional, intent(in) :: bcast     ! bcast on read?
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                    ! netCDF variable id
    integer :: ier                      ! error code
    integer :: start(4), count(4)       ! output bounds
    logical :: varpresent               ! if true, variable is on tape
    logical :: lbcast                   ! local copy of bcast flag
    character(len=*),parameter :: subname='ncd_ioglobal_real_3d' ! subroutine name
!-----------------------------------------------------------------------

    start = 1
    count = 1
    lbcast = lbcast_def
    if (present(bcast)) then
       lbcast = bcast
    endif

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
          if (varpresent) then
             if (single_column) then
                call scam_field_offsets(ncid,'undefined',start,count)
                call check_ret(nf_get_vara_double(ncid, varid, start, count, data), subname)
             else
                call check_ret(nf_get_var_double(ncid, varid, data), subname)
             endif
          endif
       end if
       if (lbcast) then
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
       endif
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_real_3d
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine scam_field_offsets
!
! !INTERFACE:
  subroutine scam_field_offsets(ncid,dim1name,start,count)
!
! !DESCRIPTION: 
! Read/Write initial data from/to netCDF instantaneous initial data file 
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use decompMod   , only : get_proc_bounds
    use clm_varpar  , only : maxpatch
    use nanMod      , only : nan
    use clm_varctl  , only : scmlon,scmlat,single_column
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: dim1name ! dimension 1 name
    integer, intent(in)    :: ncid             ! netCDF dataset id
    integer, intent(inout) :: start(:)
    integer, intent(inout) :: count(:)
!
! !CALLED FROM: subroutine inicfields
!
! !REVISION HISTORY:
! Created by John Truesdale

!EOP
!
! !LOCAL VARIABLES:
    integer :: data_offset      ! offset into land array 1st column 
    integer :: ndata            ! number of column (or pft points to read)
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
    character(len=32) :: subname = 'scam_field_offsets'

!------------------------------------------------------------------------
    call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,latidx,lonidx)
    start(1) = lonidx
    count(1) = 1
    start(2) = latidx
    count(2) = 1
    write(6,*) trim(subname),' scam_setlatlonidx ',lonidx,latidx

    if ( firsttime) then
       write(6,*) trim(subname),' firsttime=',firsttime
       call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,latidx,lonidx)

       call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
       write(6,*) trim(subname),' beg,end=',begg, endg, begl, endl, begc, endc, begp, endp

       ret = nf_inq_dimid (ncid, 'column', dimid)
       write(6,*) trim(subname),' column ret=',ret,nf_noerr
       if(ret==NF_NOERR) return
       if (ret/=NF_EBADDIM) print *, 'NETCDF ERROR: ', NF_STRERROR(ret)

       ret = nf_inq_dimlen (ncid, dimid, totcols)
       if (ret/=NF_NOERR) print *, 'NETCDF ERROR: ', NF_STRERROR(ret)

       ret = nf_inq_dimid (ncid, 'pft', dimid)
       write(6,*) trim(subname),' pft ret=',ret,nf_noerr
       if(ret==NF_NOERR) return
       if (ret/=NF_EBADDIM) print *, 'NETCDF ERROR: ', NF_STRERROR(ret)

       ret = nf_inq_dimlen (ncid, dimid, totpfts)
       if (ret/=NF_NOERR) print *, 'NETCDF ERROR: ', NF_STRERROR(ret)

       write(6,*) trim(subname),' totals ',totcols,totpfts

       allocate (pfts1dlon(totpfts))
       allocate (pfts1dlat(totpfts))
       
       ret =  nf_inq_varid (ncid, 'pfts1d_ixy', varid)
       if (ret/=NF_NOERR) then
         write(6,*)'inq_varid: id for pfts1d_ixy not found'
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
         return
       end if

       ret = nf_get_var_double (ncid, varid, pfts1dlon)
       if (ret/=NF_NOERR) then
         write(6,*)'GET_VAR_REALX: error reading pfts1dlon, varid =', varid
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
         return
       end if

       ret =  nf_inq_varid (ncid, 'pfts1d_jxy', varid)
       if (ret/=NF_NOERR) then
         write(6,*)'inq_varid: id for pfts1d_jxy not found'
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
         return
       end if

       ret = nf_get_var_double (ncid, varid, pfts1dlat)
       if (ret/=NF_NOERR) then
         write(6,*)'GET_VAR_REALX: error reading pfts1dlat, varid =', varid
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
         return
       end if


       allocate (cols1dlon(totcols))
       allocate (cols1dlat(totcols))

       ret =  nf_inq_varid (ncid, 'cols1d_ixy', varid)
       if (ret/=NF_NOERR) then
         write(6,*)'inq_varid: id for cols1d_ixy not found'
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
         return
       end if

       ret = nf_get_var_double (ncid, varid, cols1dlon)
       if (ret/=NF_NOERR) then
         write(6,*)'GET_VAR_REALX: error reading cols1dlon, varid =', varid
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
         return
       end if

       ret =  nf_inq_varid (ncid, 'cols1d_jxy', varid)
       if (ret/=NF_NOERR) then
         write(6,*)'inq_varid: id for cols1d_jxy not found'
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
         return
       end if

       ret = nf_get_var_double (ncid, varid, cols1dlat)
       if (ret/=NF_NOERR) then
         write(6,*)'GET_VAR_REALX: error reading cols1dlat, varid =', varid
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
         return
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
!          call endrun
       end if
       if (i.eq.1) then
          write(6,*)'couldnt find any columns for this latitude ',latidx,' and longitude ',lonidx
!          call endrun
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
!          call endrun
       end if
       if (i.eq.1) then
          write(6,*)'couldnt find any pfts for this latitude ',latidx,' and longitude ',lonidx
!          call endrun
       else
          pi_offset=pfts(1)
       end if

       deallocate (pfts1dlon)
       deallocate (pfts1dlat)
       deallocate (cols1dlon)
       deallocate (cols1dlat)
       firsttime = .false.
    endif

    write(6,*) trim(subname),' offsets ',pi_offset,col_offset
    
    if (dim1name == 'pft') then
       data_offset = pi_offset
       ndata = endp-begp+1
    else if (dim1name == 'column') then
       data_offset = col_offset
       ndata = endc-begc+1
    else
       write(6,*)'error calculation array offsets for SCAM'
!       call endrun()
    endif
  end subroutine scam_field_offsets
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine scam_field_offsets_old
!
! !INTERFACE:
  subroutine scam_field_offsets_old(ncid,dim1name,data_offset, ndata)
!
! !DESCRIPTION: 
! Read/Write initial data from/to netCDF instantaneous initial data file 
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use decompMod   , only : get_proc_bounds
    use clm_varpar  , only : maxpatch
    use nanMod      , only : nan
    use clm_varctl  , only : scmlon,scmlat,single_column
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
  end subroutine scam_field_offsets_old

end module ncdio
