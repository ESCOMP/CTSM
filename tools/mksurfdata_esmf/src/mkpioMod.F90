module mkpioMod

  use ESMF
  use pio
  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_kind_mod , only : i2 => shr_kind_i2, i4 => shr_kind_i4
  use shr_kind_mod , only : cl => shr_kind_cl, cs => shr_kind_cs
  use shr_sys_mod  , only : shr_sys_abort
  use mkutilsMod   , only : chkerr
  use mkvarctl     , only : root_task, ndiag, mpicom, outnc_1d

  implicit none
  private

#include <mpif.h>

  public :: mkpio_get_rawdata
  public :: mkpio_get_rawdata_level
  public :: mkpio_iodesc_rawdata
  public :: mkpio_iodesc_output
  public :: mkpio_wopen
  public :: mkpio_close
  public :: mkpio_defvar
  public :: mkpio_def_spatial_var
  public :: mkpio_get_dimlengths
  public :: mkpio_put_time_slice

  interface mkpio_get_rawdata_level
     module procedure mkpio_get_rawdata1d_level_real4
     module procedure mkpio_get_rawdata1d_level_real8
     module procedure mkpio_get_rawdata2d_level_real8
  end interface mkpio_get_rawdata_level

  interface mkpio_get_rawdata
     module procedure mkpio_get_rawdata1d_int
     module procedure mkpio_get_rawdata1d_real4
     module procedure mkpio_get_rawdata1d_real8
     module procedure mkpio_get_rawdata2d_real4
     module procedure mkpio_get_rawdata2d_real8
  end interface mkpio_get_rawdata

  interface mkpio_def_spatial_var
     module procedure mkpio_def_spatial_var_0lev
     module procedure mkpio_def_spatial_var_1lev
     module procedure mkpio_def_spatial_var_2lev
  end interface mkpio_def_spatial_var

  interface mkpio_put_time_slice
     module procedure mkpio_put_time_slice_1d
     module procedure mkpio_put_time_slice_2d
  end interface mkpio_put_time_slice

  integer               , public :: pio_iotype
  integer               , public :: pio_ioformat
  type(iosystem_desc_t) , public :: pio_iosystem

  logical :: debug = .false.

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mkpio_get_rawdata1d_int(pioid, varname, mesh_i, data_i, rc)

    ! input/output variables
    type(file_desc_t), intent(inout) :: pioid
    character(len=*) , intent(in)    :: varname   ! field name in rawdata file
    type(ESMF_Mesh)  , intent(in)    :: mesh_i
    integer          , intent(inout) :: data_i(:) ! input raw data
    integer          , intent(out)   :: rc

    ! local variables
    type(var_desc_t)       :: pio_varid
    integer                :: pio_vartype
    type(io_desc_t)        :: pio_iodesc
    type(io_desc_t)        :: pio_iodesc_mask
    integer                :: lsize
    integer                :: rcode
    integer                :: n
    integer(i2), allocatable :: data_short(:)
    character(len=*), parameter :: subname = 'mkpio_get_rawdata1d_int'
    !-------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get data_i - Read in varname from filename
    lsize = size(data_i)

    ! Create io descriptor for input raw data
    ! This will query the raw data file for the dimensions of the variable varname and
    ! create iodesc for either single or multi level input data
    call mkpio_iodesc_rawdata(mesh_i, trim(varname), pioid, pio_varid, pio_vartype, pio_iodesc, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_iodesc for "//trim(varname)//" in "//trim(subname))

    ! Read the input raw data
    if (pio_vartype == PIO_INT) then
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_i, rcode)
    else if (pio_vartype == PIO_SHORT) then
       allocate(data_short(lsize))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_short, rcode)
       data_i = int(data_short, i4)
       deallocate(data_short)
    else
       call shr_sys_abort(subName//" ERROR: only int and short is supported for "//trim(varname))
    end if
    call ESMF_VMLogMemInfo("After pio_read_darray for "//trim(varname)//" in "//trim(subname))

    ! Free the memory of the io descriptor
    call pio_freedecomp(pioid, pio_iodesc)
    call ESMF_VMLogMemInfo("After call to pio_freedecomp for "//trim(varname))

  end subroutine mkpio_get_rawdata1d_int

  !===============================================================
  subroutine mkpio_get_rawdata1d_real4(pioid, varname, mesh_i, data_i, rc)

    ! input/output variables
    type(file_desc_t), intent(inout) :: pioid
    character(len=*) , intent(in)    :: varname   ! field name in rawdata file
    type(ESMF_Mesh)  , intent(in)    :: mesh_i
    real(r4)         , intent(inout) :: data_i(:) ! input raw data
    integer          , intent(out)   :: rc

    ! local variables
    type(var_desc_t)          :: pio_varid
    integer                   :: pio_vartype
    type(io_desc_t)           :: pio_iodesc
    type(io_desc_t)           :: pio_iodesc_mask
    integer(i2) , allocatable :: data_short(:)
    integer(i4) , allocatable :: data_int(:)
    real(r8)    , allocatable :: data_double(:)
    integer                   :: lsize
    integer                   :: rcode
    integer                   :: n
    character(len=*), parameter :: subname = 'mkpio_get_rawdata1d_real4'
    !-------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get data_i - Read in varname from filename
    lsize = size(data_i)

    ! Create io descriptor for input raw data
    ! This will query the raw data file for the dimensions of the variable varname and
    ! create iodesc for either single or multi level input data
    call mkpio_iodesc_rawdata(mesh_i, trim(varname), pioid, pio_varid, pio_vartype, pio_iodesc, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_iodesc for varname "//trim(varname)//" in "//trim(subname))

    ! Read the input raw data
    if (pio_vartype == PIO_SHORT) then
       allocate(data_short(lsize))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_short, rcode)
       data_i(:) = real(data_short(:), kind=r8)
       deallocate(data_short)
    else if (pio_vartype == PIO_INT) then
       allocate(data_int(lsize))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_int, rcode)
       data_i(:) = real(data_int(:), kind=r4)
       deallocate(data_int)
    else if (pio_vartype == PIO_REAL) then
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_i, rcode)
    else if (pio_vartype == PIO_DOUBLE) then
       allocate(data_double(lsize))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_double, rcode)
       data_i(:) = real(data_double(:), kind=r4)
       deallocate(data_double)
    else
       call shr_sys_abort(subName//" ERROR: only real and double types are supported for "//trim(varname))
    end if
    call ESMF_VMLogMemInfo("After call to pio_read_darray for varname "//trim(varname))

    call pio_freedecomp(pioid, pio_iodesc)
    call ESMF_VMLogMemInfo("After call to pio_freedecomp for "//trim(varname))

  end subroutine mkpio_get_rawdata1d_real4

  !===============================================================
  subroutine mkpio_get_rawdata1d_real8(pioid, varname, mesh_i, data_i, nt, rc)

    ! input/output variables
    type(file_desc_t), intent(inout) :: pioid
    character(len=*) , intent(in)    :: varname   ! field name in rawdata file
    type(ESMF_Mesh)  , intent(in)    :: mesh_i
    real(r8)         , intent(inout) :: data_i(:) ! input raw data
    integer, optional, intent(in)    :: nt
    integer          , intent(out)   :: rc


    ! local variables
    type(var_desc_t)          :: pio_varid
    integer                   :: pio_vartype
    type(io_desc_t)           :: pio_iodesc
    type(io_desc_t)           :: pio_iodesc_mask
    integer(i2) , allocatable :: data_short(:)
    integer(i4) , allocatable :: data_int(:)
    real(r4)    , allocatable :: data_real(:)
    real(r8)    , allocatable :: data_double(:)
    integer                   :: lsize
    integer                   :: rcode
    integer                   :: n
    character(len=*), parameter :: subname = 'mkpio_get_rawdata1d_real8'
    !-------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get data_i - Read in varname from filename
    lsize = size(data_i)

    ! Create io descriptor for input raw data
    ! This will query the raw data file for the dimensions of the variable varname and
    ! create iodesc for either single or multi level input data
    call mkpio_iodesc_rawdata(mesh_i, trim(varname), pioid, pio_varid, pio_vartype, pio_iodesc, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (present(nt))  then
       call pio_setframe(pioid, pio_varid, int(nt,kind=Pio_Offset_Kind))
    end if

    ! Read the input raw data
    call ESMF_VMLogMemInfo("After mkpio_iodesc for varname for "//trim(varname)//" in "//trim(subname))
    if (pio_vartype == PIO_SHORT) then
       allocate(data_short(lsize))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_short, rcode)
       data_i(:) = real(data_short(:), kind=r8)
       deallocate(data_short)
    else if (pio_vartype == PIO_INT) then
       allocate(data_int(lsize))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_int, rcode)
       data_i(:) = real(data_int(:), kind=r8)
       deallocate(data_int)
    else if (pio_vartype == PIO_REAL) then
       allocate(data_real(lsize))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_real, rcode)
       data_i(:) = real(data_real(:), kind=r8)
       deallocate(data_real)
    else if (pio_vartype == PIO_DOUBLE) then
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_i, rcode)
    else
       call shr_sys_abort(subName//" ERROR: supported variable type not found for "//trim(varname))
    end if
    call ESMF_VMLogMemInfo("After call to pio_read_darrayy for varname "//trim(varname))

    call pio_freedecomp(pioid, pio_iodesc)
    call ESMF_VMLogMemInfo("After call to pio_freedecomp for "//trim(varname))

  end subroutine mkpio_get_rawdata1d_real8

  !===============================================================
  subroutine mkpio_get_rawdata2d_real4(pioid, varname, mesh_i, data_i, rc)

    ! input/output variables
    type(file_desc_t) , intent(inout) :: pioid
    character(len=*)  , intent(in)    :: varname   ! field name in rawdata file
    type(ESMF_Mesh)   , intent(in)    :: mesh_i
    real(r4)          , intent(inout) :: data_i(:,:) ! input raw data
    integer           , intent(out)   :: rc

    ! local variables
    type(var_desc_t)       :: pio_varid
    integer                :: pio_vartype
    type(io_desc_t)        :: pio_iodesc
    real(r4), allocatable  :: data_real1d(:)
    real(r4), allocatable  :: data_real2d(:,:)
    real(r8), allocatable  :: data_double1d(:)
    real(r8), allocatable  :: data_double2d(:,:)
    real(r4), allocatable  :: landmask(:)
    integer                :: lsize, nlev
    integer                :: n,l
    integer                :: rcode
    character(len=*), parameter :: subname = ' mkpio_get_rawdata_2d'
    !-------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get data_i - Read in varname from filename
    ! Note that data_i is coming in as (nlev,lsize) in terms of dimensions
    nlev  = size(data_i, dim=1)
    lsize = size(data_i, dim=2)

    ! Create io descriptor for input raw data
    ! This will query the raw data file for the dimensions of the variable varname and
    ! create iodesc for either single or multi level input data
    call ESMF_VMLogMemInfo("Before mkpio_iodesc in "//trim(subname))
    call mkpio_iodesc_rawdata(mesh_i, trim(varname), pioid, pio_varid, pio_vartype, pio_iodesc, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_iodesc in "//trim(subname))

    ! Read the input raw data (all levels read at once)
    ! - levels are the innermost dimension for esmf fields
    ! - levels are the outermost dimension in pio reads
    ! Input data is read into (lsize,nlev) array and then transferred to data_i(nlev,lsize)
    if (pio_vartype == PIO_REAL) then
       allocate(data_real2d(lsize,nlev))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_real2d, rcode)
       do l = 1,nlev
          do n = 1,lsize
             data_i(l,n) = data_real2d(n,l)
          end do
       end do
       deallocate(data_real2d)
    else if (pio_vartype == PIO_DOUBLE) then
       allocate(data_double2d(lsize,nlev))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_double2d, rcode)
       do l = 1,nlev
          do n = 1,lsize
             data_i(l,n) = real(data_double2d(n,l), kind=r4)
          end do
       end do
       deallocate(data_double2d)
    else
       call shr_sys_abort(subName//" ERROR: only real and double types are supported for "//trim(varname))
    end if
    call pio_freedecomp(pioid, pio_iodesc)

  end subroutine mkpio_get_rawdata2d_real4

  !===============================================================
  subroutine mkpio_get_rawdata2d_real8(pioid, varname, mesh_i, data_i, setframe, rc)

    ! input/output variables
    type(file_desc_t) , intent(inout) :: pioid
    character(len=*)  , intent(in)    :: varname   ! field name in rawdata file
    type(ESMF_Mesh)   , intent(in)    :: mesh_i
    real(r8)          , intent(inout) :: data_i(:,:) ! input raw data
    integer, optional , intent(in)    :: setframe
    integer           , intent(out)   :: rc

    ! local variables
    type(var_desc_t)       :: pio_varid
    integer                :: pio_vartype
    type(io_desc_t)        :: pio_iodesc
    type(io_desc_t)        :: pio_iodesc_mask
    real(r4), allocatable  :: data_real1d(:)
    real(r4), allocatable  :: data_real2d(:,:)
    real(r8), allocatable  :: data_double1d(:)
    real(r8), allocatable  :: data_double2d(:,:)
    real(r4), allocatable  :: landmask(:)
    integer                :: lsize, nlev
    integer                :: n,l
    integer                :: rcode
    character(len=*), parameter :: subname = ' mkpio_get_rawdata_2d'
    !-------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get data_i - Read in varname from filename
    ! Note that data_i is coming in as (nlev,lsize) in terms of dimensions
    nlev  = size(data_i, dim=1)
    lsize = size(data_i, dim=2)

    ! Create io descriptor for input raw data
    ! This will query the raw data file for the dimensions of the variable varname and
    ! create iodesc for either single or multi level input data
    call ESMF_VMLogMemInfo("Before mkpio_iodesc in "//trim(subname))
    call mkpio_iodesc_rawdata(mesh_i, trim(varname), pioid, pio_varid, pio_vartype, pio_iodesc, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_iodesc in "//trim(subname))

    ! Read the input raw data (all levels read at once)
    ! - levels are the innermost dimension for esmf fields
    ! - levels are the outermost dimension in pio reads
    ! Input data is read into (lsize,nlev) array and then transferred to data_i(nlev,lsize)
    if (pio_vartype == PIO_REAL) then
       allocate(data_real2d(lsize,nlev))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_real2d, rcode)
       do l = 1,nlev
          do n = 1,lsize
             data_i(l,n) = real(data_real2d(n,l), kind=r8)
          end do
       end do
       deallocate(data_real2d)
    else if (pio_vartype == PIO_DOUBLE) then
       allocate(data_double2d(lsize,nlev))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_double2d, rcode)
       do l = 1,nlev
          do n = 1,lsize
             data_i(l,n) = data_double2d(n,l)
          end do
       end do
       deallocate(data_double2d)
    else
       call shr_sys_abort(subName//"ERROR: only real and double types are supported")
    end if
    call pio_freedecomp(pioid, pio_iodesc)

  end subroutine mkpio_get_rawdata2d_real8

  !===============================================================
  subroutine mkpio_iodesc_rawdata( mesh, varname, pioid, pio_varid,  pio_vartype, pio_iodesc, rc)

    ! Determine pio io descriptor for variable on rawdata file

    ! input/output variables
    type(ESMF_Mesh)   , intent(in)    :: mesh
    character(len=*)  , intent(in)    :: varname
    type(file_desc_t) , intent(inout) :: pioid
    type(var_desc_t)  , intent(out)   :: pio_varid
    integer           , intent(out)   :: pio_vartype
    type(io_desc_t)   , intent(inout) :: pio_iodesc
    integer           , intent(out)   :: rc

    ! local variables
    type(ESMF_DistGrid)     :: distGrid
    integer                 :: n, ndims
    integer, allocatable    :: compdof(:)
    integer, allocatable    :: compdof3d(:)
    integer, allocatable    :: dimids(:)
    integer, allocatable    :: dimlens(:)
    character(len=cs)       :: dimname
    integer                 :: lsize
    integer                 :: nlev
    integer                 :: cnt, m
    integer                 :: offset
    integer                 :: rCode ! pio return code (only used when pio error handling is PIO_BCAST_ERROR)
    integer                 :: unlimdim
    logical                 :: unlimited_dim
    character(*), parameter :: subname = '(mkpio_iodesc_rawdata) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMLogMemInfo("Beginning setting compdof for "//trim(varname))
    call ESMF_MeshGet(mesh, elementdistGrid=distGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_DistGridGet(distGrid, localDe=0, elementCount=lsize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(compdof(lsize))
    call ESMF_DistGridGet(distGrid, localDe=0, seqIndexList=compdof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("Ending setting compdof for "//trim(varname))

    ! get pio variable id, type and number of dimensions
    rcode = pio_inq_varid(pioid, trim(varname), pio_varid)
    rcode = pio_inq_vartype(pioid, pio_varid, pio_vartype)
    rcode = pio_inq_varndims(pioid, pio_varid, ndims)

    ! get variable dimension sizes
    allocate(dimids(ndims))
    allocate(dimlens(ndims))
    rcode = pio_inq_vardimid(pioid, pio_varid, dimids(1:ndims))
    do n = 1, ndims
       rcode = pio_inq_dimlen(pioid, dimids(n), dimlens(n))
    end do
    rcode = pio_inq_unlimdim(pioid, unlimdim)
    unlimited_dim = (dimids(ndims) == unlimdim)

    ! Create compdof3d if needed
    ! Assume that input data is always lon,lat as first two dimensions
    nlev = 0
    if (ndims == 3 .and. .not. unlimited_dim) then
       nlev = dimlens(3)
    else if (ndims == 3 .and. unlimited_dim) then
       ! do nothing - keep nlev at 0
    else if (ndims == 4 .and. .not. unlimited_dim) then
       nlev = dimlens(3)*dimlens(4)
    else if (ndims == 4 .and. unlimited_dim) then
       nlev = dimlens(3)
    end if

    if (nlev > 0) then
       offset = dimlens(1)*dimlens(2)
       allocate(compdof3d(nlev*lsize))
       cnt = 0
       do n = 1,nlev
          do m = 1,size(compdof)
             cnt = cnt + 1
             compdof3d(cnt) = (n-1)*offset + compdof(m)
          enddo
       enddo
    end if

    ! determine io descriptor for this variable
    if (ndims == 1) then
       call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1)/), compdof, pio_iodesc)
       if (root_task .and. debug) then
          write(ndiag,'(a,i20)') ' set iodesc for rawdata: '//trim(varname)//' with dim(1) = ',&
               dimlens(1)
       end if
    else if (ndims == 2) then
       call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1),dimlens(2)/), compdof, pio_iodesc)
       if (root_task .and. debug) then
          write(ndiag,'(a,i8,i8)') ' set iodesc for rawdata: '//trim(varname)//' with dim(1),dim(2) = ',&
               dimlens(1),dimlens(2)
       end if
    else if (ndims == 3) then
       if (unlimited_dim) then
          call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1),dimlens(2)/), compdof, pio_iodesc)
          if (root_task .and. debug) then
             write(ndiag,'(a,i8,i8)') ' set iodesc for rawdata: '//trim(varname)//' with dim(1),dim(2) = ',&
                  dimlens(1),dimlens(2)
          end if
       else
          call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1),dimlens(2),dimlens(3)/), compdof3d, pio_iodesc)
          if (root_task .and. debug) then
             write(ndiag,'(a,i8,i8,i8)') ' set iodesc for rawdata: '//trim(varname)//' with dim(1),dim(2),dim(3) = ',&
                  dimlens(1),dimlens(2),dimlens(3)
          end if
       end if
    else if (ndims == 4) then
       if (unlimited_dim) then
          call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1),dimlens(2),dimlens(3)/), compdof3d, pio_iodesc)
          if (root_task .and. debug) then
             write(ndiag,'(a,i8,i8,i8)') ' set iodesc for rawdata: '//trim(varname)//' with dim(1),dim(2),dim(3) = ',&
                  dimlens(1),dimlens(2),dimlens(3)
          end if
       else
          write(6,*)' ndims = ',ndims
          write(6,*)' unlimited_dim = ',unlimited_dim
          call shr_sys_abort('for lon/lat support up to 3 input spatial dims plus a time dim')
       end if
    else
       call shr_sys_abort('rawdata input for variable '//trim(varname)//'  must have ndims either 1,2,3 or 4')
    end if

    ! deallocate memory
    deallocate(compdof)
    if (allocated(compdof3d)) deallocate(compdof3d)
    call ESMF_VMLogMemInfo("Finished setting iodesc for "//trim(varname)//" in "//trim(subname))

  end subroutine mkpio_iodesc_rawdata

  !===============================================================
  subroutine mkpio_iodesc_output(pioid, mesh, varname, pio_iodesc, rc)

    ! Create pio_iodesc for varname

    ! input/output variables
    type(file_desc_t) , intent(inout) :: pioid
    character(len=*)  , intent(in)    :: varname
    type(ESMF_Mesh)   , intent(in)    :: mesh
    type(io_desc_t)   , intent(inout) :: pio_iodesc
    integer           , intent(out)   :: rc

    ! local variables
    type(ESMF_DistGrid)     :: distGrid
    integer                 :: n, ndims
    integer, allocatable    :: compdof(:)
    integer, allocatable    :: compdof3d(:)
    integer, allocatable    :: dimids(:)
    integer, allocatable    :: dimlens(:)
    character(len=cs)       :: dimname
    integer                 :: lsize
    integer                 :: nlev
    integer                 :: cnt, m
    integer                 :: offset
    type(var_desc_t)        :: pio_varid
    integer                 :: pio_vartype
    integer                 :: rCode ! pio return code (only used when pio error handling is PIO_BCAST_ERROR)
    integer                 :: unlimdim
    logical                 :: unlimited_dim
    character(*), parameter :: subname = '(shr_strdata_set_stream_iodesc) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! get pio variable id, type and dimension information
    rcode = pio_inq_varid(pioid, trim(varname), pio_varid)
    rcode = pio_inq_vartype(pioid, pio_varid, pio_vartype)
    rcode = pio_inq_varndims(pioid, pio_varid, ndims)
    allocate(dimids(ndims))
    allocate(dimlens(ndims))
    rcode = pio_inq_vardimid(pioid, pio_varid, dimids(1:ndims))
    do n = 1, ndims
       rcode = pio_inq_dimlen(pioid, dimids(n), dimlens(n))
    end do
    rcode = pio_inq_unlimdim(pioid, unlimdim)
    unlimited_dim = (dimids(ndims) == unlimdim)

    ! Get compdof from mesh
    call ESMF_MeshGet(mesh, elementdistGrid=distGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_DistGridGet(distGrid, localDe=0, elementCount=lsize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(compdof(lsize))
    call ESMF_DistGridGet(distGrid, localDe=0, seqIndexList=compdof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Create compdof3d if needed
    nlev = 0
    if (outnc_1d) then
       offset = dimlens(1)
       if (ndims == 2 .and. .not. unlimited_dim) then
          nlev = dimlens(2)
       else if (ndims == 3 .and. unlimited_dim) then
          nlev = dimlens(2)
       else if (ndims == 3 .and. .not. unlimited_dim) then
          nlev = dimlens(2)*dimlens(3)
       end if
    else
       offset = dimlens(1)*dimlens(2)
       if (ndims == 3 .and. .not. unlimited_dim) then
          nlev = dimlens(3)
       else if (ndims == 3 .and. unlimited_dim) then
          ! do nothing - keep nlev at 0
       else if (ndims == 4 .and. .not. unlimited_dim) then
          nlev = dimlens(3)*dimlens(4)
       else if (ndims == 4 .and. unlimited_dim) then
          nlev = dimlens(3)
       end if
    end if

    if (nlev > 0) then
       allocate(compdof3d(nlev*lsize))
       cnt = 0
       do n = 1,nlev
          do m = 1,size(compdof)
             cnt = cnt + 1
             compdof3d(cnt) = (n-1)*offset + compdof(m)
          enddo
       enddo
    end if

    ! determine io descriptor for this variable
    if (outnc_1d) then
       ! Assume that can have (gridcell), (gridcell,lev), (gridcell,lev,time)
       ! Where lev would correspond to an undistributed dimension in esmf
       if (ndims == 1)  then
          call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1)/), compdof, pio_iodesc)
          if (root_task .and. debug) then
             write(ndiag,'(a,i20)') ' set iodesc for output data: '//trim(varname)//' with dim(1) = ',&
                  dimlens(1)
          end if
       else if (ndims == 2) then
          if (unlimited_dim) then
             call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1)/), compdof, pio_iodesc)
             if (root_task .and. debug) then
                write(ndiag,'(a,i8)') ' set iodesc for output data with time dim: '//trim(varname)//&
                     ' with dim(1) = ',dimlens(1)
             end if
          else
             call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1),dimlens(2)/), compdof3d, pio_iodesc)
             if (root_task .and. debug) then
                write(ndiag,'(a,i8,i8)') ' set iodesc for output data: '//trim(varname)//&
                     ' with dim(1),dim(2) = ',dimlens(1),dimlens(2)
             end if
          end if
       else if (ndims == 3) then
          if (unlimited_dim) then
             call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1),dimlens(2)/), compdof3d, pio_iodesc)
             if (root_task .and. debug) then
                write(ndiag,'(a,i8,i8)') ' set iodesc for output data with time dim: '//trim(varname)//&
                     ' with dim(1),dim(2) = ',dimlens(1),dimlens(2)
             end if
          else
             call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1),dimlens(2),dimlens(3)/), compdof3d, pio_iodesc)
             if (root_task .and. debug) then
                write(ndiag,'(a,i8,i8,i8)') ' set iodesc for output data: '//trim(varname)//&
                     ' with dim(1),dim(2),dim(3) = ',dimlens(1),dimlens(2),dimlens(3)
             end if
          end if
       end if
    else
       ! Assume that can have (lon,lat), (lon,lat,lev1), (lon,lat,lev1,lev2), (lon,lat,time) or (lon,lat,lev1,time)
       if (ndims == 2) then
          call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1),dimlens(2)/), compdof, pio_iodesc)
          if (root_task .and. debug) then
             write(ndiag,'(a,i8,i8)') ' set iodesc for output data: '//trim(varname)//&
                  ' with dim(1),dim(2)= ',dimlens(1),dimlens(2)
          end if
       else if (ndims == 3) then
          if (unlimited_dim) then
             call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1),dimlens(2)/), compdof, pio_iodesc)
             if (root_task .and. debug) then
                write(ndiag,'(a,i8,i8)') ' set iodesc for output data with time dim : '//trim(varname)//&
                     ' with dim(1),dim(2)= ', dimlens(1),dimlens(2)
             end if
          else
             call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1),dimlens(2),dimlens(3)/), compdof3d, pio_iodesc)
             if (root_task .and. debug) then
                write(ndiag,'(a,i8,i8,i8)') ' set iodesc for output data: '//trim(varname)//&
                     ' with dim(1),dim(2),dim3(3)= ',dimlens(1),dimlens(2),dimlens(3)
             end if
          end if
       else if (ndims == 4) then
          if (unlimited_dim) then
             call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1),dimlens(2),dimlens(3)/), compdof3d, pio_iodesc)
             if (root_task .and. debug) then
                write(ndiag,'(a,i8,i8,i8)') ' set iodesc for output data with time dim : '//trim(varname)//&
                     ' with dim(1),dim(2),dimlens(3)= ', dimlens(1),dimlens(2),dimlens(3)
             end if
          else
             call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1),dimlens(2),dimlens(3),dimlens(4)/), compdof3d, pio_iodesc)
             if (root_task .and. debug) then
                write(ndiag,'(a,i8,i8,i8,i8)') ' set iodesc for output data: '//trim(varname)//&
                     ' with dim(1),dim(2),dimlens(3),dimlens(4)= ', dimlens(1),dimlens(2),dimlens(3),dimlens(4)
             end if
          end if
       end if
    end if

    ! deallocate memory
    deallocate(compdof)
    if (allocated(compdof3d)) deallocate(compdof3d)

  end subroutine mkpio_iodesc_output

  !===============================================================================
  logical function mkpio_file_exists(filename)

    !---------------
    ! inquire if i/o file exists
    !---------------

    ! input/output variables
    character(len=*) , intent(in) :: filename

    ! local variables
    integer :: tmp(1)
    integer :: ier
    !-------------------------------------------------------------------------------

    tmp(1) = 0
    mkpio_file_exists = .false.
    if (root_task) then
       inquire(file=trim(filename), exist=mkpio_file_exists)
       if (mkpio_file_exists) tmp(1) = 1
    end if
    call mpi_bcast(tmp(1), 1, MPI_INTEGER, 0, mpicom, ier)
    if (tmp(1) == 1) mkpio_file_exists = .true.

  end function mkpio_file_exists

  !===============================================================================
  subroutine mkpio_wopen(filename, clobber, pioid)

    !---------------
    ! open netcdf file
    !---------------

    use pio , only : PIO_IOTYPE_PNETCDF, PIO_IOTYPE_NETCDF, PIO_BCAST_ERROR, PIO_INTERNAL_ERROR
    use pio , only : pio_openfile, pio_createfile, PIO_GLOBAL, pio_enddef
    use pio , only : pio_put_att, pio_get_att
    use pio , only : pio_seterrorhandling, pio_file_is_open, pio_clobber, pio_write, pio_noclobber

    ! input/output arguments
    character(len=*)  , intent(in)    :: filename
    logical           , intent(in)    :: clobber
    type(file_desc_t) , intent(inout) :: pioid

    ! local variables
    integer           :: rcode
    integer           :: nmode
    character(*),parameter :: subName = '(mkpio_wopen) '
    !-------------------------------------------------------------------------------

    ! filename not open
    if (root_task) then
       write(ndiag,'(a)') "opening output file "//trim(filename)
    end if

    if (mkpio_file_exists(filename)) then
       if (clobber) then
          nmode = pio_clobber
          ! only applies to classic NETCDF files.
          if(pio_iotype == PIO_IOTYPE_NETCDF .or. pio_iotype == PIO_IOTYPE_PNETCDF) then
             nmode = ior(nmode,pio_ioformat)
          endif
          rcode = pio_createfile(pio_iosystem, pioid, pio_iotype, trim(filename), nmode)
          if (root_task) write(ndiag,'(a)') trim(subname)//' creating file '//trim(filename)
       else
          rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(filename), pio_write)
          if (root_task) write(ndiag,'(a)') trim(subname)//' opening file '//trim(filename)
       endif
    else
       ! only applies to classic NETCDF files.
       nmode = pio_noclobber
       if (pio_iotype == PIO_IOTYPE_NETCDF .or. pio_iotype == PIO_IOTYPE_PNETCDF) then
          nmode = ior(nmode,pio_ioformat)
       endif
       rcode = pio_createfile(pio_iosystem, pioid, pio_iotype, trim(filename), nmode)
       if (root_task) write(ndiag,'(a)') trim(subname) //' creating file '// trim(filename)
    endif
    call ESMF_LogWrite("successfully opened output file "//trim(filename), ESMF_LOGMSG_INFO)

  end subroutine mkpio_wopen

  !===============================================================================
  subroutine mkpio_close(pioid, filename, rc)

    !---------------
    ! close netcdf file
    !---------------

    use pio, only: pio_file_is_open, pio_closefile

    ! input/output variables
    type(file_desc_t) , intent(in)  :: pioid
    character(*)      , intent(in)  :: filename
    integer           , intent(out) :: rc

    ! local variables
    character(len=CL) :: wfilename
    character(*),parameter :: subName = '(mkpio_close) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (.not. pio_file_is_open(pioid)) then
       ! filename not open, just return
    elseif (trim(wfilename) == trim(filename)) then
       ! filename matches, close it
       call pio_closefile(pioid)
    else
       ! different filename is open, abort
       if (root_task) then
          write(ndiag,'(a)') trim(subname)//' different  wfilename and filename currently open, aborting '
          write(ndiag,'(a)') 'filename  = ',trim(filename)
          write(ndiag,'(a)') 'wfilename = ',trim(wfilename)
       end if
       call ESMF_LogWrite(subname//'different file currently open, aborting '//trim(filename), ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
    endif

  end subroutine mkpio_close

  !===============================================================================
  subroutine mkpio_defvar(pioid, varname, xtype, &
       dim1name, dim2name, dim3name, dim4name, dim5name, &
       long_name, units, missing_value, fill_value, imissing_value, ifill_value)

    !  Define a pio variable

    ! input/output variables
    type(file_desc_t) , intent(in)           :: pioid
    character(len=*)  , intent(in)           :: varname        ! variable name
    integer           , intent(in)           :: xtype          ! external type
    character(len=*)  , intent(in), optional :: dim1name       ! dimension name
    character(len=*)  , intent(in), optional :: dim2name       ! dimension name
    character(len=*)  , intent(in), optional :: dim3name       ! dimension name
    character(len=*)  , intent(in), optional :: dim4name       ! dimension name
    character(len=*)  , intent(in), optional :: dim5name       ! dimension name
    character(len=*)  , intent(in), optional :: long_name      ! attribute
    character(len=*)  , intent(in), optional :: units          ! attribute
    real(r8)          , intent(in), optional :: missing_value  ! attribute for real
    real(r8)          , intent(in), optional :: fill_value     ! attribute for real
    integer           , intent(in), optional :: imissing_value ! attribute for int
    integer           , intent(in), optional :: ifill_value    ! attribute for int

    ! !LOCAL VARIABLES:
    integer            :: n         ! indices
    integer            :: ndims     ! dimension counter
    integer            :: dimid(5)  ! dimension ids
    type(var_desc_t)   :: pio_varid ! variable id
    integer            :: itmp      ! temporary
    character(len=CS)  :: str       ! temporary
    integer            :: rcode
    character(*), parameter :: subname='mkpio_defvar_real' ! subroutine name
    !-----------------------------------------------------------------------

    ! Determine dimension ids for variable

    dimid(:) = 0

    if (present(dim1name)) then
       rcode = pio_inq_dimid(pioid, dim1name, dimid(1))
    end if
    if (present(dim2name)) then
       rcode = pio_inq_dimid(pioid, dim2name, dimid(2))
    end if
    if (present(dim3name)) then
       rcode = pio_inq_dimid(pioid, dim3name, dimid(3))
    end if
    if (present(dim4name)) then
       rcode = pio_inq_dimid(pioid, dim4name, dimid(4))
    end if
    if (present(dim5name)) then
       rcode = pio_inq_dimid(pioid, dim5name, dimid(5))
    end if

    ! Define variable
    if (present(dim1name)) then
       ndims = 0
       do n = 1, size(dimid)
          if (dimid(n) /= 0) ndims = ndims + 1
       end do
       rcode = pio_def_var(pioid, trim(varname), xtype, dimid(1:ndims), pio_varid)
    else
       rcode = pio_def_var(pioid, trim(varname), xtype, pio_varid)
    end if

    ! Add attributes to variable
    if (present(long_name)) then
       rcode = pio_put_att(pioid, pio_varid, 'long_name', trim(long_name))
    end if
    if (present(units)) then
       rcode = pio_put_att(pioid, pio_varid, 'units', trim(units))
    end if
    if (present(fill_value)) then
       rcode = pio_put_att(pioid, pio_varid, '_FillValue', fill_value)
    end if
    if (present(missing_value)) then
       rcode = pio_put_att(pioid, pio_varid, 'missing_value', missing_value)
    end if
    if (present(ifill_value)) then
       rcode = pio_put_att(pioid, pio_varid, '_FillValue', ifill_value)
    end if
    if (present(imissing_value)) then
       rcode = pio_put_att(pioid, pio_varid, 'missing_value', imissing_value)
    end if

  end subroutine mkpio_defvar

  ! ========================================================================
  ! mkpio_def_spatial_var routines: define a spatial pio variable
  ! (convenience wrapper to mkpio_defvar)
  ! ========================================================================

  subroutine mkpio_def_spatial_var_0lev(pioid, varname, xtype, long_name, units)

    ! Define a spatial netCDF variable (convenience wrapper to mkpio_defvar)
    ! The variable in question has ONLY spatial dimensions (no level or time dimensions)

    ! !ARGUMENTS:
    type(file_desc_t) , intent(in) :: pioid
    character(len=*)  , intent(in) :: varname   ! variable name
    integer           , intent(in) :: xtype     ! external type
    character(len=*)  , intent(in) :: long_name ! attribute
    character(len=*)  , intent(in) :: units     ! attribute

    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'mkpio_def_spatial_var_0lev'
    !-----------------------------------------------------------------------

    if (outnc_1d) then
       call mkpio_defvar(pioid=pioid, varname=varname, xtype=xtype, &
            dim1name='gridcell', long_name=long_name, units=units)
    else
       call mkpio_defvar(pioid=pioid, varname=varname, xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', long_name=long_name, units=units)
    end if

  end subroutine mkpio_def_spatial_var_0lev

  !-----------------------------------------------------------------------
  subroutine mkpio_def_spatial_var_1lev(pioid, varname, xtype, lev1name, long_name, units)
    !
    ! Define a spatial netCDF variable (convenience wrapper to mkpio_defvar)
    ! The variable in question has one level (or time) dimension in addition to its
    ! spatial dimensions

    ! input/output variables
    type(file_desc_t) , intent(in)           :: pioid
    character(len=*) , intent(in) :: varname   ! variable name
    integer          , intent(in) :: xtype     ! external type
    character(len=*) , intent(in) :: lev1name  ! name of level (or time) dimension
    character(len=*) , intent(in) :: long_name ! attribute
    character(len=*) , intent(in) :: units     ! attribute

    ! local variables
    character(len=*), parameter :: subname = 'mkpio_def_spatial_var_1lev'
    !-----------------------------------------------------------------------

    if (outnc_1d) then
       call mkpio_defvar(pioid, varname=varname, xtype=xtype, &
            dim1name='gridcell', dim2name=lev1name, &
            long_name=long_name, units=units)
    else
       call mkpio_defvar(pioid, varname=varname, xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat',dim3name=lev1name, &
            long_name=long_name, units=units)
    end if

  end subroutine mkpio_def_spatial_var_1lev

  !-----------------------------------------------------------------------
  subroutine mkpio_def_spatial_var_2lev(pioid, varname, xtype, lev1name, lev2name, long_name, units)
    !
    ! Define a spatial netCDF variable (convenience wrapper to mkpio_defvar)
    !
    ! The variable in question has two level (or time) dimensions in addition to its
    ! spatial dimensions
    !
    ! input/output variables
    type(file_desc_t) , intent(in) :: pioid
    character(len=*)  , intent(in) :: varname   ! variable name
    integer           , intent(in) :: xtype     ! external type
    character(len=*)  , intent(in) :: lev1name  ! name of first level (or time) dimension
    character(len=*)  , intent(in) :: lev2name  ! name of second level (or time) dimension
    character(len=*)  , intent(in) :: long_name ! attribute
    character(len=*)  , intent(in) :: units     ! attribute

    ! local variables:
    character(len=*), parameter :: subname = 'mkpio_def_spatial_var_2lev'
    !-----------------------------------------------------------------------

    if (outnc_1d) then
       call mkpio_defvar(pioid=pioid, varname=varname, xtype=xtype, &
            dim1name='gridcell', dim2name=lev1name, dim3name=lev2name, &
            long_name=long_name, units=units)
    else
       call mkpio_defvar(pioid=pioid, varname=varname, xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', dim3name=lev1name, dim4name=lev2name, &
            long_name=long_name, units=units)
    end if

  end subroutine mkpio_def_spatial_var_2lev

  ! ========================================================================
  subroutine mkpio_put_time_slice_1d(pioid, pio_varid, pio_iodesc, time_index, data)

    ! Write a single time slice of a 1-d variable

    ! input/output variables
    type(file_desc_t) , intent(inout) :: pioid
    type(var_desc_t)  , intent(inout) :: pio_varid
    type(io_desc_t)   , intent(inout) :: pio_iodesc
    integer           , intent(in)    :: time_index ! time index in file
    real(r8)          , intent(in)    :: data(:)    ! data to write (a single time slice)
    !
    ! local variables:
    integer :: rcode
    character(len=*), parameter :: subname = 'mkpio_put_time_slice_1d'
    !-----------------------------------------------------------------------

    call pio_setframe(pioid, pio_varid, int(time_index, kind=Pio_Offset_Kind))
    call pio_write_darray(pioid, pio_varid, pio_iodesc, data, rcode)

  end subroutine mkpio_put_time_slice_1d

  ! ========================================================================
  subroutine mkpio_put_time_slice_2d(pioid, pio_varid, pio_iodesc, time_index, data)

    ! Write a single time slice of a 2-d variable

    ! input/output variables
    type(file_desc_t) , intent(inout) :: pioid
    type(var_desc_t)  , intent(inout) :: pio_varid
    type(io_desc_t)   , intent(inout) :: pio_iodesc
    integer           , intent(in)    :: time_index ! time index in file
    real(r8)          , intent(in)    :: data(:,:)  ! data to write (a single time slice)
    !
    ! local variables:
    integer :: rcode
    character(len=*), parameter :: subname = 'mkpio_put_time_slice_2d'
    !-----------------------------------------------------------------------

    call pio_setframe(pioid, pio_varid, int(time_index, kind=Pio_Offset_Kind))
    call pio_write_darray(pioid, pio_varid, pio_iodesc, data, rcode)

  end subroutine mkpio_put_time_slice_2d

  !===============================================================
  subroutine mkpio_get_dimlengths(pioid, varname, ndims, dimlengths)

    ! Returns the number of dimensions and an array containing the dimension lengths of a
    ! variable in an open netcdf file.
    ! Entries 1:ndims in the returned dim_lengths array contain the dimension lengths; the
    ! remaining entries in that vector are meaningless. The dim_lengths array must be large
    ! enough to hold all ndims values; if not, the code aborts (this can be ensured by passing
    ! in an array of length nf_max_var_dims).
    !
    ! input/otuput variables
    type(file_desc_t) , intent(in)  :: pioid
    character(len=*)  , intent(in)  :: varname        ! name of variable of interest
    integer           , intent(out) :: ndims          ! number of dimensions of variable
    integer           , intent(out) :: dimlengths(:)  ! lengths of dimensions of variable
    !
    ! local variables
    type(var_desc_t)     :: pio_varid
    integer, allocatable :: dimids(:)
    integer              :: i
    integer              :: rcode
    character(len=*), parameter :: subname = 'mkpio_get_dimlengths'
    !------------------------------------------------------------------------------

    rcode = pio_inq_varid(pioid, trim(varname), pio_varid)
    rcode = pio_inq_varndims(pioid, pio_varid, ndims)

    if (ndims > size(dimlengths)) then
       call shr_sys_abort(trim(subname)//' ERROR: dimlengths too small')
    end if

    allocate(dimids(ndims))
    rcode = pio_inq_vardimid(pioid, pio_varid, dimids(1:ndims))
    dimlengths(:) = 0  ! pre-fill with 0 so we won't have garbage in elements past ndims
    do i = 1, ndims
       rcode = pio_inq_dimlen(pioid, dimids(i), dimlengths(i))
    end do
    deallocate(dimids)

  end subroutine mkpio_get_dimlengths

  ! ========================================================================
  subroutine mkpio_get_rawdata1d_level_real4(pioid, pio_iodesc, unlimited_index, varname, data_i)

    ! input/output variables
    type(file_desc_t), intent(inout) :: pioid
    type(io_desc_t)  , intent(inout) :: pio_iodesc
    integer          , intent(in)    :: unlimited_index
    character(len=*) , intent(in)    :: varname   ! field name in rawdata file
    real(r4)         , intent(inout) :: data_i(:) ! input raw data

    ! local variables
    type(var_desc_t)          :: pio_varid
    integer                   :: pio_vartype
    integer(i2) , allocatable :: data_short(:)
    integer(i4) , allocatable :: data_int(:)
    real(r8)    , allocatable :: data_double(:)
    integer                   :: ns_i
    integer                   :: rcode
    character(len=*), parameter :: subname = 'mkpio_get_rawdata_level_real4'
    !-------------------------------------------------

    ! Get variable id and type
    rcode = pio_inq_varid(pioid, trim(varname), pio_varid)
    rcode = pio_inq_vartype(pioid, pio_varid, pio_vartype)

    ! Set unlimited frame index
    call pio_setframe(pioid, pio_varid, int(unlimited_index, kind=Pio_Offset_Kind))

    ! Read the input raw data
    ns_i = size(data_i)
    if (pio_vartype == PIO_SHORT) then
       allocate(data_short(ns_i))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_short, rcode)
       data_i(:) = real(data_short(:), kind=r8)
       deallocate(data_short)
    else if (pio_vartype == PIO_INT) then
       allocate(data_int(ns_i))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_int, rcode)
       data_i(:) = real(data_int(:), kind=r4)
       deallocate(data_int)
    else if (pio_vartype == PIO_REAL) then
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_i, rcode)
    else if (pio_vartype == PIO_DOUBLE) then
       ns_i = size(data_i)
       allocate(data_double(ns_i))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_double, rcode)
       data_i(:) = real(data_double(:), kind=r4)
       deallocate(data_double)
    else
       call shr_sys_abort(subName//" ERROR: vartype not supported for "//trim(varname))
    end if
    call ESMF_VMLogMemInfo("After call to pio_read_darray for varname "//trim(varname))

  end subroutine mkpio_get_rawdata1d_level_real4

  ! ========================================================================
  subroutine mkpio_get_rawdata1d_level_real8(pioid, pio_iodesc, unlimited_index, varname, data_i)

    ! input/output variables
    type(file_desc_t), intent(inout) :: pioid
    type(io_desc_t)  , intent(inout) :: pio_iodesc
    integer          , intent(in)    :: unlimited_index
    character(len=*) , intent(in)    :: varname   ! field name in rawdata file
    real(r8)         , intent(inout) :: data_i(:) ! input raw data

    ! local variables
    type(var_desc_t)          :: pio_varid
    integer                   :: pio_vartype
    integer(i2) , allocatable :: data_short(:)
    integer(i4) , allocatable :: data_int(:)
    real(r4)    , allocatable :: data_real(:)
    integer                   :: ns_i
    integer                   :: rcode
    character(len=*), parameter :: subname = 'mkpio_get_rawdata1d_real4'
    !-------------------------------------------------

    ! Get variable id and type
    rcode = pio_inq_varid(pioid, trim(varname), pio_varid)
    rcode = pio_inq_vartype(pioid, pio_varid, pio_vartype)

    ! Set unlimited frame index
    call pio_setframe(pioid, pio_varid, int(unlimited_index, kind=Pio_Offset_Kind))

    ! Read the input raw data
    ns_i = size(data_i)
    if (pio_vartype == PIO_SHORT) then
       allocate(data_short(ns_i))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_short, rcode)
       data_i(:) = real(data_short(:), kind=r8)
       deallocate(data_short)
    else if (pio_vartype == PIO_INT) then
       allocate(data_int(ns_i))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_int, rcode)
       data_i(:) = real(data_int(:), kind=r4)
       deallocate(data_int)
    else if (pio_vartype == PIO_REAL) then
       allocate(data_real(ns_i))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_real, rcode)
       data_i(:) = real(data_real(:), kind=r8)
       deallocate(data_real)
    else if (pio_vartype == PIO_DOUBLE) then
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_i, rcode)
    else
       call shr_sys_abort(subName//" ERROR: vartype not supported for "//trim(varname))
    end if
    call ESMF_VMLogMemInfo("After call to pio_read_darray for varname "//trim(varname))

  end subroutine mkpio_get_rawdata1d_level_real8

  ! ========================================================================
  subroutine mkpio_get_rawdata2d_level_real8(pioid, pio_iodesc, unlimited_index, varname, data_i)

    ! input/output variables
    type(file_desc_t), intent(inout) :: pioid
    type(io_desc_t)  , intent(inout) :: pio_iodesc
    integer          , intent(in)    :: unlimited_index
    character(len=*) , intent(in)    :: varname     ! field name in rawdata file
    real(r8)         , intent(inout) :: data_i(:,:) ! input raw data

    ! local variables
    type(var_desc_t)          :: pio_varid
    integer                   :: pio_vartype
    integer(i2) , allocatable :: data_short2d(:,:)
    integer(i4) , allocatable :: data_int2d(:,:)
    real(r4)    , allocatable :: data_real2d(:,:)
    real(r8)    , allocatable :: data_double2d(:,:)
    integer                   :: ns_i, nlev
    integer                   :: n, l
    integer                   :: rcode
    character(len=*), parameter :: subname = 'mkpio_get_rawdata1d_real4'
    !-------------------------------------------------

    ! Get variable id and type
    rcode = pio_inq_varid(pioid, trim(varname), pio_varid)
    rcode = pio_inq_vartype(pioid, pio_varid, pio_vartype)

    ! Set unlimited frame index
    call pio_setframe(pioid, pio_varid, int(unlimited_index, kind=Pio_Offset_Kind))

    ! Read the input raw data
    ! - levels are the innermost dimension for esmf fields
    ! - levels are the outermost dimension in pio reads
    ! Input data is read into (ns_i,nlev) array and then transferred to data_i(nlev,ns_i)
    ! which sill be used for esmf regridding
    nlev = size(data_i, dim=1)
    ns_i = size(data_i, dim=2)
    if (pio_vartype == PIO_SHORT) then
       allocate(data_short2d(ns_i,nlev))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_short2d, rcode)
       do l = 1,nlev
          do n = 1,ns_i
             data_i(l,n) = real(data_short2d(n,l), kind=r8)
          end do
       end do
       deallocate(data_short2d)
    else if (pio_vartype == PIO_INT) then
       allocate(data_int2d(ns_i,nlev))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_int2d, rcode)
       do l = 1,nlev
          do n = 1,ns_i
             data_i(l,n) = real(data_int2d(n,l), kind=r8)
          end do
       end do
       deallocate(data_int2d)
    else if (pio_vartype == PIO_REAL) then
       allocate(data_real2d(ns_i,nlev))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_real2d, rcode)
       do l = 1,nlev
          do n = 1,ns_i
             data_i(l,n) = real(data_real2d(n,l), kind=r8)
          end do
       end do
       deallocate(data_real2d)
    else if (pio_vartype == PIO_DOUBLE) then
       allocate(data_double2d(ns_i,nlev))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_double2d, rcode)
       do l = 1,nlev
          do n = 1,ns_i
             data_i(l,n) = data_double2d(n,l)
          end do
       end do
       deallocate(data_double2d)
    else
       call shr_sys_abort(subName//" ERROR: vartype not supported for "//trim(varname))
    end if
    call ESMF_VMLogMemInfo("After call to pio_read_darray for varname "//trim(varname))

  end subroutine mkpio_get_rawdata2d_level_real8

end module mkpioMod
