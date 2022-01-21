module mkpioMod

  use ESMF ! TODO: put in only statements
  use pio  ! TODO: put in only statements
  use shr_kind_mod , only : r8 => shr_kind_r8, cl => shr_kind_cl, cs => shr_kind_cs, r4 => shr_kind_r4
  use shr_sys_mod  , only : shr_sys_abort
  use mkutilsMod   , only : chkerr
  use mkvarctl     , only : root_task, ndiag, mpicom

  implicit none
  private

#include <mpif.h>

  public :: mkpio_get_rawdata
  public :: mkpio_iodesc_rawdata
  public :: mkpio_iodesc_output
  public :: mkpio_wopen
  public :: mkpio_close
  public :: mkpio_defvar
  public :: mkpio_def_spatial_var
  public :: mkpio_get_dim_lengths
  public :: mkpio_put_time_slice

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

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mkpio_get_rawdata(filename, varname, mesh_i, data_i, rc)

    ! input/output variables
    character(len=*)       , intent(in)    :: filename  ! file name of rawdata file
    character(len=*)       , intent(in)    :: varname   ! field name in rawdata file
    type(ESMF_Mesh)        , intent(in)    :: mesh_i
    real(r8)               , intent(inout) :: data_i(:) ! input raw data
    integer                , intent(out)   :: rc

    ! local variables
    type(file_desc_t)      :: pioid
    type(var_desc_t)       :: pio_varid
    type(io_desc_t)        :: pio_iodesc 
    integer                :: pio_vartype
    real(r4), allocatable  :: data_real(:)
    real(r8), allocatable  :: data_double(:)
    real(r8), pointer      :: dataptr(:)
    integer                :: lsize 
    integer                :: rcode
    character(len=*), parameter :: subname = 'mklakwat'
    !-------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get data_i - Read in varname from filename 
    lsize = size(data_i)
    call ESMF_VMLogMemInfo("Before pio_openfile in regrid_data for "//trim(filename))
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(filename), pio_nowrite)
    call ESMF_VMLogMemInfo("After pio_openfile in regrid_data")
    call ESMF_VMLogMemInfo("After field get")
    call mkpio_iodesc_rawdata(mesh_i, trim(varname), pioid, pio_varid, pio_vartype, pio_iodesc, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_iodesc in regrid_data")
    if (pio_vartype == PIO_REAL) then
       allocate(data_real(lsize))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_real, rcode)
       data_i(:) = real(data_real(:), kind=r8)
       deallocate(data_real)
    else if (pio_vartype == PIO_DOUBLE) then
       allocate(data_double(lsize))
       call pio_read_darray(pioid, pio_varid, pio_iodesc, data_double, rcode)
       data_i(:) = data_double(:)
       deallocate(data_double)
    else
       call shr_sys_abort(subName//"ERROR: only real and double types are supported")
    end if
    call pio_freedecomp(pioid, pio_iodesc)
    call pio_closefile(pioid)
    call ESMF_VMLogMemInfo("After pio_read_darry in regrid_data")

  end subroutine mkpio_get_rawdata

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
    integer, pointer        :: compdof(:)
    integer, allocatable    :: dimids(:)
    integer, allocatable    :: dimlens(:)
    integer                 :: lsize
    integer                 :: rCode ! pio return code (only used when pio error handling is PIO_BCAST_ERROR)
    character(*), parameter :: subname = '(shr_strdata_set_stream_iodesc) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !if (root_task) then
       !write(ndiag,'(a,i8, i8)') 'setting iodesc for : '//trim(varname)// &
       !' with dimlens(1), dimlens2 = ',dimlens(1),dimlens(2)
    !end if

    call ESMF_VMLogMemInfo("Beginning setting compdof")
    call ESMF_MeshGet(mesh, elementdistGrid=distGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_DistGridGet(distGrid, localDe=0, elementCount=lsize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(compdof(lsize))
    call ESMF_DistGridGet(distGrid, localDe=0, seqIndexList=compdof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("Ending setting compdof") 

    ! get pio variable id, type and number of dimensions
    call ESMF_VMLogMemInfo("Beginning getting variable id")
    call PIO_seterrorhandling(pioid, PIO_BCAST_ERROR)
    rcode = pio_inq_varid(pioid, trim(varname), pio_varid)
    rcode = pio_inq_vartype(pioid, pio_varid, pio_vartype)
    rcode = pio_inq_varndims(pioid, pio_varid, ndims)
    call ESMF_VMLogMemInfo("Ending getting variable id")

    ! get variable dimension sizes
    call ESMF_VMLogMemInfo("Beginning getting dims")
    allocate(dimids(ndims))
    allocate(dimlens(ndims))
    rcode = pio_inq_vardimid(pioid, pio_varid, dimids(1:ndims))
    do n = 1, ndims
       rcode = pio_inq_dimlen(pioid, dimids(n), dimlens(n))
    end do
    call ESMF_VMLogMemInfo("End getting dims")

    ! determine io descriptor for this variable
    if (ndims == 1) then
       call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1)/), compdof, pio_iodesc)
    else if (ndims == 2) then
       call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1),dimlens(2)/), compdof, pio_iodesc)
    else
       call shr_sys_abort('rawdata input must have ndims either 1 or 2')
    end if

    ! deallocate memory
    deallocate(compdof)

  end subroutine mkpio_iodesc_rawdata

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
  subroutine mkpio_wopen(pioid, filename, clobber)

    !---------------
    ! open netcdf file
    !---------------

    use pio , only : PIO_IOTYPE_PNETCDF, PIO_IOTYPE_NETCDF, PIO_BCAST_ERROR, PIO_INTERNAL_ERROR
    use pio , only : pio_openfile, pio_createfile, PIO_GLOBAL, pio_enddef
    use pio , only : pio_put_att, pio_get_att
    use pio , only : pio_seterrorhandling, pio_file_is_open, pio_clobber, pio_write, pio_noclobber

    ! input/output arguments
    type(file_desc_t) , intent(inout) :: pioid
    character(len=*)  , intent(in)    :: filename
    logical           , intent(in)    :: clobber

    ! local variables
    integer           :: rcode
    integer           :: nmode
    character(*),parameter :: subName = '(mkpio_wopen) '
    !-------------------------------------------------------------------------------

    ! filename not open
    call ESMF_LogWrite("opening output file "//trim(filename), ESMF_LOGMSG_INFO)

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

    use mkvarctl, only : outnc_1d

    ! !ARGUMENTS:
    type(file_desc_t) , intent(in)           :: pioid
    character(len=*) , intent(in) :: varname   ! variable name
    integer          , intent(in) :: xtype     ! external type
    character(len=*) , intent(in) :: long_name ! attribute
    character(len=*) , intent(in) :: units     ! attribute

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
    ! !DESCRIPTION:
    ! Define a spatial netCDF variable (convenience wrapper to mkpio_defvar)
    !
    ! The variable in question has one level (or time) dimension in addition to its
    ! spatial dimensions

    use mkvarctl, only : outnc_1d

    ! !ARGUMENTS:
    type(file_desc_t) , intent(in)           :: pioid
    character(len=*) , intent(in) :: varname   ! variable name
    integer          , intent(in) :: xtype     ! external type
    character(len=*) , intent(in) :: lev1name  ! name of level (or time) dimension
    character(len=*) , intent(in) :: long_name ! attribute
    character(len=*) , intent(in) :: units     ! attribute

    ! !LOCAL VARIABLES:
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
    ! !DESCRIPTION:
    ! Define a spatial netCDF variable (convenience wrapper to mkpio_defvar)
    !
    ! The variable in question has two level (or time) dimensions in addition to its
    ! spatial dimensions
    !
    ! !USES:
    use mkvarctl, only : outnc_1d
    !
    ! !ARGUMENTS:
    type(file_desc_t) , intent(in)           :: pioid
    character(len=*) , intent(in) :: varname   ! variable name
    integer          , intent(in) :: xtype     ! external type
    character(len=*) , intent(in) :: lev1name  ! name of first level (or time) dimension
    character(len=*) , intent(in) :: lev2name  ! name of second level (or time) dimension
    character(len=*) , intent(in) :: long_name ! attribute
    character(len=*) , intent(in) :: units     ! attribute

    ! !LOCAL VARIABLES:
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

  !-----------------------------------------------------------------------
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
  subroutine mkpio_get_dim_lengths(pioid, varname, ndims, dim_lengths)

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
    integer           , intent(out) :: dim_lengths(:) ! lengths of dimensions of variable
    !
    ! local variables
    type(var_desc_t) :: pio_varid
    integer          :: dimids(size(dim_lengths))
    integer          :: i
    integer          :: rcode
    character(len=*), parameter :: subname = 'get_dim_lengths'
    !------------------------------------------------------------------------------

    rcode = pio_inq_varid(pioid, trim(varname), pio_varid)
    rcode = pio_inq_varndims(pioid, pio_varid, ndims)

    if (ndims > size(dim_lengths)) then
       call shr_sys_abort(trim(subname)//' ERROR: dim_lengths too small')
    end if

    dim_lengths(:) = 0  ! pre-fill with 0 so we won't have garbage in elements past ndims
    do i = 1, ndims
       rcode = pio_inq_dimlen(pioid, dimids(i), dim_lengths(i))
    end do

  end subroutine mkpio_get_dim_lengths

  !===============================================================
  subroutine mkpio_iodesc_output(pioid, mesh, varname, pio_iodesc, rc)

    ! Write variable varname to output file

    ! input/output variables
    type(file_desc_t) , intent(inout) :: pioid
    character(len=*)  , intent(in)    :: varname
    type(ESMF_Mesh)   , intent(in)    :: mesh
    type(io_desc_t)   , intent(inout) :: pio_iodesc 
    
    integer           , intent(out)   :: rc

    ! local variables
    type(ESMF_DistGrid)     :: distGrid
    integer                 :: n, ndims
    integer, pointer        :: compdof(:)
    integer, allocatable    :: dimids(:)
    integer, allocatable    :: dimlens(:)
    integer                 :: lsize
    type(var_desc_t)        :: pio_varid
    integer                 :: pio_vartype
    character(len=cs)       :: dimname
    integer                 :: rCode ! pio return code (only used when pio error handling is PIO_BCAST_ERROR)
    character(*), parameter :: subname = '(shr_strdata_set_stream_iodesc) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       !write(ndiag,'(a,i8, i8)') 'setting iodesc for : '//trim(varname)//' with dimlens(1), dimlens2 = ',dimlens(1),dimlens(2)
       write(ndiag,'(a)') 'setting iodesc for : '//trim(varname)
    end if

    call ESMF_MeshGet(mesh, elementdistGrid=distGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_DistGridGet(distGrid, localDe=0, elementCount=lsize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(compdof(lsize))
    call ESMF_DistGridGet(distGrid, localDe=0, seqIndexList=compdof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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

    ! determine io descriptor for this variable
    if (ndims == 1)  then
       call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1)/), compdof, pio_iodesc)
    else if (ndims == 2) then
       if (trim(dimname) == 'time' .or. trim(dimname) == 'nt') then
          call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1)/), compdof, pio_iodesc)
       else
          call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1),dimlens(2)/), compdof, pio_iodesc)
       end if
    else if (ndims == 3) then
       if (trim(dimname) == 'time' .or. trim(dimname) == 'nt') then
          call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1),dimlens(2)/), compdof, pio_iodesc)
       else
          call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1),dimlens(2),dimlens(3)/), compdof, pio_iodesc)
       end if
    else if (ndims == 4) then
    end if

    ! deallocate memory
    deallocate(compdof)

  end subroutine mkpio_iodesc_output

end module mkpioMod
