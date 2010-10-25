module ncdio_pio

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdio_pioMod
!
! !DESCRIPTION:
! Generic interfaces to write fields to netcdf files
!
! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8, i8=>shr_kind_i8, shr_kind_cl
  use spmdMod        , only : masterproc, mpicom, iam, npes,  &
	                      MPI_REAL8, MPI_INTEGER, MPI_LOGICAL
  use clmtype        , only : gratm, grlnd, nameg, namel, namec, namep, allrof
  use clm_varcon     , only : spval,ispval
  use clm_varctl     , only : single_column, iulog
  use shr_sys_mod    , only : shr_sys_flush
  use shr_file_mod   , only : shr_file_getunit, shr_file_freeunit
  use shr_string_mod , only : shr_string_toUpper
  use abortutils     , only : endrun
  use decompMod      , only : get_clmlevel_gsize,get_clmlevel_gsmap
  use perf_mod       , only : t_startf, t_stopf
  use fileutils      , only : getavu, relavu
  use clm_mct_mod
  use pio
!
! !PUBLIC TYPES:
  implicit none
  include 'netcdf.inc'
  save

!
! !PUBLIC MEMBER FUNCTIONS:
!
  public :: check_ret     ! checks return status of netcdf calls
  public :: check_var     ! determine if variable is on netcdf file
  public :: check_dim     ! validity check on dimension
  public :: ncd_pio_openfile
  public :: ncd_pio_createfile
  public :: ncd_pio_init  ! called from clm_comp
  public :: ncd_enddef    ! end define mode
  public :: ncd_putatt    ! put attribute
  public :: ncd_defdim    ! define dimension
  public :: ncd_inqdid    ! inquire dimension id
  public :: ncd_inqdname  ! inquire dimension name
  public :: ncd_inqdlen   ! inquire dimension length
  public :: ncd_defvar    ! define variables
  public :: ncd_inqvid    ! inquire variable id
  public :: ncd_inqvname  ! inquire variable name
  public :: ncd_inqvdims  ! inquire variable ndims
  public :: ncd_inqvdids  ! inquire variable dimids
  public :: ncd_io        ! write local data

  integer,parameter,public :: ncd_int       = nf_int
  integer,parameter,public :: ncd_float     = nf_float
  integer,parameter,public :: ncd_double    = nf_double
  integer,parameter,public :: ncd_char      = nf_char
  integer,parameter,public :: ncd_global    = nf_global
  integer,parameter,public :: ncd_write     = nf_write
  integer,parameter,public :: ncd_nowrite   = nf_nowrite
  integer,parameter,public :: ncd_clobber   = nf_clobber
  integer,parameter,public :: ncd_noclobber = nf_noclobber
  integer,parameter,public :: ncd_share     = nf_share
  integer,parameter,public :: ncd_fill      = nf_fill
  integer,parameter,public :: ncd_nofill    = nf_nofill
  integer,parameter,public :: ncd_unlimited = nf_unlimited
!
! !REVISION HISTORY:
!
!
! !PRIVATE MEMBER FUNCTIONS:
!

  interface ncd_putatt
     module procedure ncd_putatt_int
     module procedure ncd_putatt_real
     module procedure ncd_putatt_char
  end interface

  interface ncd_defvar
     module procedure ncd_defvar_bynf
     module procedure ncd_defvar_bygrid
  end interface

  interface ncd_io 
     module procedure ncd_io_int_var0_nf
     module procedure ncd_io_real_var0_nf
     module procedure ncd_io_int_var1_nf
     module procedure ncd_io_real_var1_nf
     module procedure ncd_io_int_var2_nf
     module procedure ncd_io_real_var2_nf
     module procedure ncd_io_int_var1
     module procedure ncd_io_real_var1
     module procedure ncd_io_int_var2
     module procedure ncd_io_real_var2
     module procedure ncd_io_int_var3
     module procedure ncd_io_real_var3
  end interface

  private :: ncd_inqiodesc      ! inquire variable descriptor
  private :: ncd_getiodesc      ! obtain iodesc
  private :: scam_field_offsets ! get offset to proper lat/lon gridcell for SCAM

  integer,parameter,private :: debug = 0             ! local debug level

  integer , parameter  , public  :: max_string_len = 256     ! length of strings
  real(r8), parameter  , public  :: fillvalue = 1.e36_r8     ! fill value for netcdf fields

  integer, private :: io_stride, num_iotasks, io_type

  type(iosystem_desc_t), public  :: pio_subsystem

  type iodesc_plus_type
     character(len=32) :: name
     type(IO_desc_t)   :: iodesc
     integer           :: type
     integer           :: ndims
     integer           :: dims(4)
  end type iodesc_plus_type
  integer,parameter      ,private :: max_iodesc = 100
  integer                ,private :: num_iodesc = 0
  type(iodesc_plus_type) ,private, target :: iodesc_list(max_iodesc)

!EOP
!-----------------------------------------------------------------------

contains

  subroutine ncd_pio_init(nlfilename)
    character(len=*), intent(in) :: nlfilename

    call read_namelist_pio(nlfilename)

    call PIO_init(iam, mpicom, num_iotasks, 1, io_stride, PIO_Rearr_box, PIO_subsystem)

  end subroutine ncd_pio_init

  !=======================================================================

  subroutine read_namelist_pio(nlfilename)
    character(len=*), intent(in) :: nlfilename

    character(len=80):: iotype_name

    namelist /clmpio_nml/ io_stride, num_iotasks, iotype_name
    integer  :: nu_nml    ! unit for namelist file
    integer  :: nml_error ! namelist i/o error flag
    integer  :: ierr      ! error status

    io_stride   = -1   ! set based on num_iotasks value when initialized < 0
    num_iotasks = -1   ! set based on io_stride value when initialized < 0
    iotype_name = 'netcdf'

    if (masterproc) then
       nu_nml = getavu()
       open( nu_nml, file=trim(nlfilename), status='old', iostat=nml_error )
       if (nml_error /= 0) then
          nml_error = -1
       else
          nml_error =  1
       endif
       do while (nml_error > 0)
          read(nu_nml, nml=clmpio_nml,iostat=nml_error)
          if (nml_error > 0) read(nu_nml,*)  ! for Nagware compiler
       end do
       call relavu( nu_nml )
       call set_pio_parameters(iotype_name)
       
       write(iulog,*) 'CLM PIO parameter settings...'
       write(iulog,*) '  io_stride   = ',io_stride
       write(iulog,*) '  iotype_name = ',iotype_name
       write(iulog,*) '  num_iotasks = ',num_iotasks
    end if

    call mpi_bcast(io_type    , 1, mpi_integer, 0, mpicom, ierr)
    call mpi_bcast(io_stride  , 1, mpi_integer, 0, mpicom, ierr)
    call mpi_bcast(num_iotasks, 1, mpi_integer, 0, mpicom, ierr)

  end subroutine read_namelist_pio

  !=======================================================================

  subroutine ncd_pio_openfile(file, fname, mode)
    type(file_desc_t), intent(inout) :: file
    character(len=*) , intent(in)    :: fname
    integer          , intent(in)    :: mode

    integer :: ierr

    ierr = pio_openfile(pio_subsystem, file, io_type, fname, mode)

    if(ierr/= PIO_NOERR) then
       call endrun('Failed to open file')
    else if(pio_subsystem%io_rank==0) then
       write(iulog,*) 'Opened existing file ', trim(fname), file%fh
    end if

  end subroutine ncd_pio_openfile

  !=======================================================================

  subroutine ncd_pio_createfile(file, fname)
    type(file_desc_t), intent(inout) :: file
    character(len=*),  intent(in)    :: fname

    integer :: ierr

    ierr = pio_createfile(pio_subsystem, file, io_type, fname, ior(PIO_CLOBBER,PIO_64BIT_OFFSET))

    if(ierr/= PIO_NOERR) then
       call endrun('Failed to open restart file to write')
    else if(pio_subsystem%io_rank==0) then
       write(iulog,*) 'Opened file ', trim(fname),  ' to write', file%fh
    end if

  end subroutine ncd_pio_createfile

  !=======================================================================

  subroutine set_pio_parameters(io_type_name)
    character(len=*), intent(in) :: io_type_name

    character(len=16) :: tmpname

    tmpname = shr_string_toupper(io_type_name)

    if(tmpname .eq. 'NETCDF') then
       io_type = iotype_netcdf
    else if(tmpname .eq. 'PNETCDF') then
       io_type = iotype_pnetcdf
    else if(tmpname .eq. 'NETCDF4P') then
       io_type = pio_iotype_netcdf4p
    else if(tmpname .eq. 'NETCDF4C') then
       io_type = pio_iotype_netcdf4c
    else
       write(iulog,*) 'Bad io_type argument - using iotype_netcdf'
       io_type=iotype_netcdf
    end if

    if(io_stride>0.and.num_iotasks<0) then
       num_iotasks = npes/io_stride
    else if(num_iotasks>0 .and. io_stride<0) then
       io_stride = npes/num_iotasks
    else if(num_iotasks<0 .and. io_stride<0) then
       io_stride = max(1,npes/4)
       num_iotasks = npes/io_stride
    end if

    if(io_stride*num_iotasks> npes .or. io_stride<=0 .or. num_iotasks<=0) then
       write(iulog,*) 'io_stride or num_iotasks out of bounds - resetting to defaults ',io_stride, num_iotasks 
       io_stride = max(1,npes/4)
       num_iotasks = npes/io_stride
    end if
    write(iulog,*) 'Using io_type ',tmpname,' with stride ',io_stride,' and ',num_iotasks,' io tasks'

  end subroutine set_pio_parameters

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_var
!
! !INTERFACE:
  subroutine check_var(ncid, varname, vardesc, readvar)
!
! !DESCRIPTION:
! Check if variable is on netcdf file
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout)  :: ncid
    character(len=*) , intent(in)     :: varname
    type(Var_desc_t) , intent(out)    :: vardesc
    logical          , intent(out)    :: readvar
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: ret     ! return value
    character(len=*),parameter :: subname='check_var' ! subroutine name
!-----------------------------------------------------------------------

    readvar = .true.
    call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
    ret = PIO_inq_varid (ncid, varname, vardesc)
    if (ret /= PIO_noerr) then
       readvar = .false.
       if (masterproc) write(iulog,*) trim(subname),': variable ',trim(varname),' is not on dataset'
    end if
    call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
    
  end subroutine check_var

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
    type(file_desc_t),intent(in) :: ncid
    character(len=*), intent(in) :: dimname
    integer, intent(in)          :: value
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: dimid, dimlen    ! temporaries
    integer :: status           ! error code	
    character(len=*),parameter :: subname='check_dim' ! subroutine name
!-----------------------------------------------------------------------

    status = pio_inq_dimid (ncid, trim(dimname), dimid)
    status = pio_inq_dimlen (ncid, dimid, dimlen)
    if (dimlen /= value) then
       write(iulog,*) trim(subname),' ERROR: mismatch of input dimension ',dimlen, &
            ' with expected value ',value,' for variable ',trim(dimname)
       call endrun()
    end if

  end subroutine check_dim

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_ret
!
! !INTERFACE:
  subroutine check_ret(ret, cstring)
!
! !DESCRIPTION:
! Check return status from netcdf call
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ret
    character(len=*) :: cstring
!
! !REVISION HISTORY:
!
!EOP
    character(len=*),parameter :: subname='check_ret' ! subroutine name
!-----------------------------------------------------------------------

    if (ret /= NF_NOERR) then
       write(iulog,*)'netcdf error from ',trim(subname),':',trim(cstring),':',trim(NF_STRERROR(ret))
       call endrun()
    end if

  end subroutine check_ret

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_enddef
!
! !INTERFACE:
  subroutine ncd_enddef(ncid)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    status = PIO_enddef(ncid)

  end subroutine ncd_enddef

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqdid
!
! !INTERFACE:
  subroutine ncd_inqdid(ncid,name,dimid)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(in) :: ncid      ! netcdf file id
    character(len=*), intent(in) :: name      ! dimension name
    integer         , intent(out):: dimid     ! dimension id
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    status = PIO_inq_dimid(ncid,name,dimid)

  end subroutine ncd_inqdid

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqdlen
!
! !INTERFACE:
  subroutine ncd_inqdlen(ncid,dimid,len)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(in) :: ncid      ! netcdf file id
    integer         , intent(in) :: dimid     ! dimension id
    integer         , intent(out):: len       ! dimension len
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    len = -1
    status = PIO_inq_dimlen(ncid,dimid,len)

  end subroutine ncd_inqdlen

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqdname
!
! !INTERFACE:
  subroutine ncd_inqdname(ncid,dimid,dname)
!
! !DESCRIPTION:
! inquire dim name
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(in) :: ncid      ! netcdf file id
    integer          , intent(in) :: dimid     ! dimension id
    character(len=*) , intent(out):: dname     ! dimension name
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    status = PIO_inq_dimname(ncid,dimid,dname)

  end subroutine ncd_inqdname

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqiodesc
!
! !INTERFACE:
  subroutine ncd_inqiodesc(ndims,dims,type,iodnum)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    integer         ,intent(in) :: ndims     ! number of dims
    integer         ,intent(in) :: dims(:)   ! dims
    integer         ,intent(in) :: type      ! data type
    integer         ,intent(out):: iodnum    ! iodesc num
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    logical :: found             ! search flag
    integer :: n,m
    character(len=*),parameter :: subname='ncd_inqiodesc' ! subroutine name

!-----------------------------------------------------------------------

    iodnum = -1

    n = 1
    found = .false.
    do while (n <= num_iodesc .and. .not.found)
       if (ndims == iodesc_list(n)%ndims .and. type == iodesc_list(n)%type) then
          found = .true.
          do m = 1,ndims
             if (dims(m) /= iodesc_list(n)%dims(m)) then
                found = .false.
             endif
          enddo
       endif
       if (found) then
          iodnum = n
       endif
       n = n + 1
    enddo

  end subroutine ncd_inqiodesc

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqvid
!
! !INTERFACE:
  subroutine ncd_inqvid(ncid,name,varid,vardesc,readvar)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid      ! netcdf file id
    character(len=*) , intent(in)    :: name      ! variable name
    integer          , intent(out)   :: varid     ! variable id
    type(Var_desc_t) , intent(out)   :: vardesc   ! variable descriptor
    logical, optional, intent(out)   :: readvar   ! does variable exist
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: ret               ! return code
    character(len=*),parameter :: subname='ncd_inqvid' ! subroutine name
!-----------------------------------------------------------------------

    if (present(readvar)) then
       readvar = .false.
       call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
       ret = PIO_inq_varid(ncid,name,vardesc)
       if (ret /= PIO_noerr) then
          if (masterproc) write(iulog,*) trim(subname),': variable ',trim(name),' is not on dataset'
          readvar = .false.
       else
          readvar = .true.
       end if
       call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
    else
       ret = PIO_inq_varid(ncid,name,vardesc)
    endif
    varid = vardesc%varid
 
  end subroutine ncd_inqvid

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqvdims
!
! !INTERFACE:
  subroutine ncd_inqvdims(ncid,varid,ndims,vardesc)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(in)   :: ncid      ! netcdf file id
    integer          , intent(in)   :: varid     ! variable id
    integer          , intent(out)  :: ndims     ! variable ndims
    type(Var_desc_t) , intent(inout):: vardesc
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    ndims = -1
    status = PIO_inq_varndims(ncid,vardesc,ndims)

  end subroutine ncd_inqvdims

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqvname
!
! !INTERFACE:
  subroutine ncd_inqvname(ncid,varid,vname,vardesc)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(in)   :: ncid      ! netcdf file id
    integer          , intent(in)   :: varid     ! variable id
    character(len=*) , intent(out)  :: vname     ! variable vname
    type(Var_desc_t) , intent(inout):: vardesc
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    vname = ''
    status = PIO_inq_varname(ncid,vardesc,vname)

  end subroutine ncd_inqvname

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqvdids
!
! !INTERFACE:
  subroutine ncd_inqvdids(ncid,dids,vardesc)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(in)  :: ncid      ! netcdf file id
    integer         ,intent(out)  :: dids(:)   ! variable dids
    type(Var_desc_t),intent(inout):: vardesc
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    dids = -1
    status = PIO_inq_vardimid(ncid,vardesc,dids)

  end subroutine ncd_inqvdids

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_putatt_int
!
! !INTERFACE:
  subroutine ncd_putatt_int(ncid,varid,attrib,value,xtype)
!
! !DESCRIPTION:
! Check return status from netcdf call
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
    integer          ,intent(in)    :: varid     ! netcdf var id
    character(len=*) ,intent(in)    :: attrib    ! netcdf attrib
    integer          ,intent(in)    :: value     ! netcdf attrib value
    integer,optional ,intent(in)    :: xtype     ! netcdf data type
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    status = PIO_put_att(ncid,varid,trim(attrib),value)

  end subroutine ncd_putatt_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_putatt_char
!
! !INTERFACE:
  subroutine ncd_putatt_char(ncid,varid,attrib,value,xtype)
!
! !DESCRIPTION:
! Check return status from netcdf call
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
    integer          ,intent(in)    :: varid     ! netcdf var id
    character(len=*) ,intent(in)    :: attrib    ! netcdf attrib
    character(len=*) ,intent(in)    :: value     ! netcdf attrib value
    integer,optional ,intent(in)    :: xtype     ! netcdf data type
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    status = PIO_put_att(ncid,varid,trim(attrib),value)

  end subroutine ncd_putatt_char

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_putatt_real
!
! !INTERFACE:
  subroutine ncd_putatt_real(ncid,varid,attrib,value,xtype)
!
! !DESCRIPTION:
! Check return status from netcdf call
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
    integer          ,intent(in)    :: varid     ! netcdf var id
    character(len=*) ,intent(in)    :: attrib    ! netcdf attrib
    real(r8)         ,intent(in)    :: value     ! netcdf attrib value
    integer          ,intent(in)    :: xtype     ! netcdf data type
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
    real*4  :: value4
!-----------------------------------------------------------------------

    value4 = value

    if (xtype == nf_double) then
       status = PIO_put_att(ncid,varid,trim(attrib),value)
    else
       status = PIO_put_att(ncid,varid,trim(attrib),value4)
    endif

  end subroutine ncd_putatt_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_defdim
!
! !INTERFACE:
  subroutine ncd_defdim(ncid,attrib,value,dimid)
!
! !DESCRIPTION:
! Check return status from netcdf call
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(in) :: ncid      ! netcdf file id
    character(len=*) , intent(in) :: attrib    ! netcdf attrib
    integer          , intent(in) :: value     ! netcdf attrib value
    integer          , intent(out):: dimid     ! netcdf dimension id
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    status = pio_def_dim(ncid,attrib,value,dimid)

  end subroutine ncd_defdim

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_defvar_bynf
!
! !INTERFACE:
  subroutine ncd_defvar_bynf(ncid, varname, xtype, ndims, dimid, varid, &
                             long_name, units, cell_method, missing_value, fill_value, &
                             imissing_value, ifill_value)
!
! !DESCRIPTION:
!  Define a netcdf variable
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid                  ! netcdf file id
    character(len=*) , intent(in)  :: varname                 ! variable name
    integer          , intent(in)  :: xtype                   ! external type
    integer          , intent(in)  :: ndims                   ! number of dims
    integer          , intent(inout) :: varid                 ! returned var id
    integer          , intent(in), optional :: dimid(:)       ! dimids
    character(len=*) , intent(in), optional :: long_name      ! attribute
    character(len=*) , intent(in), optional :: units          ! attribute
    character(len=*) , intent(in), optional :: cell_method    ! attribute
    real(r8)         , intent(in), optional :: missing_value  ! attribute for real
    real(r8)         , intent(in), optional :: fill_value     ! attribute for real
    integer          , intent(in), optional :: imissing_value ! attribute for int
    integer          , intent(in), optional :: ifill_value    ! attribute for int
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: n                   ! indices
    integer :: ldimid(4)           ! local dimid
    integer :: dimid0(1)           ! local dimid
    integer :: status              ! error status 
    type(var_desc_t)   :: vardesc  ! local vardesc
    character(len=128) :: dimname  ! temporary
    character(len=256) :: str      ! temporary
    character(len=*),parameter :: subname='ncd_defvar_bynf' ! subroutine name
!-----------------------------------------------------------------------

    varid = -1

    dimid0 = 0
    ldimid = 0
    if (present(dimid)) then
       ldimid(1:ndims) = dimid(1:ndims)
    else   ! ndims must be zero if dimid not present
       if (ndims /= 0) then
          write(iulog,*) trim(subname),' ERROR: dimid not supplied and ndims ne 0 ',trim(varname),ndims
          call endrun()
       endif
    endif
       
    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(varname),xtype,ndims,ldimid(1:ndims)
    endif
    
    if (ndims >  0) then 
       status = pio_inq_dimname(ncid,ldimid(ndims),dimname)
    end if

    ! Define variable
    if (present(dimid)) then
       status = PIO_def_var(ncid,trim(varname),xtype,dimid(1:ndims),vardesc)
    else
       status = PIO_def_var(ncid,trim(varname),xtype,dimid0        ,vardesc)
    endif
    varid = vardesc%varid

    if (present(long_name)) then
       call ncd_putatt(ncid, varid, 'long_name', trim(long_name))
    end if
    if (present(units)) then
       call ncd_putatt(ncid, varid, 'units', trim(units))
    end if
    if (present(cell_method)) then
       str = 'time: ' // trim(cell_method)
       call ncd_putatt(ncid, varid, 'cell_methods', trim(str))
    end if
    if (present(fill_value)) then
       call ncd_putatt(ncid, varid, '_FillValue', fill_value, xtype)
    end if
    if (present(missing_value)) then
       call ncd_putatt(ncid, varid, 'missing_value', missing_value, xtype)
    end if
    if (present(ifill_value)) then
       call ncd_putatt(ncid, varid, '_FillValue', ifill_value, xtype)
    end if
    if (present(imissing_value)) then
       call ncd_putatt(ncid, varid, 'missing_value', imissing_value, xtype)
    end if

  end subroutine ncd_defvar_bynf

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_defvar_bygrid
!
! !INTERFACE:
  subroutine ncd_defvar_bygrid(ncid, varname, xtype, &
                               dim1name, dim2name, dim3name, dim4name, dim5name, &
                               long_name, units, cell_method, missing_value, fill_value, &
                               imissing_value, ifill_value, switchdim)
!
! !DESCRIPTION:
!  Define a netcdf variable
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid                 ! netcdf file id
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
    logical         , intent(in), optional :: switchdim      ! true=> permute dim1 and dim2 for output
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: n              ! indices
    integer :: ndims          ! dimension counter
    integer :: dimid(5)       ! dimension ids
    integer :: varid          ! variable id
    integer :: itmp           ! temporary
    character(len=256) :: str ! temporary
    character(len=*),parameter :: subname='ncd_defvar_bygrid' ! subroutine name
!-----------------------------------------------------------------------

    dimid(:) = 0

    ! Determine dimension ids for variable

    if (present(dim1name)) call ncd_inqdid(ncid, dim1name, dimid(1))
    if (present(dim2name)) call ncd_inqdid(ncid, dim2name, dimid(2))
    if (present(dim3name)) call ncd_inqdid(ncid, dim3name, dimid(3))
    if (present(dim4name)) call ncd_inqdid(ncid, dim4name, dimid(4))
    if (present(dim5name)) call ncd_inqdid(ncid, dim5name, dimid(5))

    ! Permute dim1 and dim2 if necessary

    if (present(switchdim)) then
       itmp = dimid(2)
       dimid(2) = dimid(1)
       dimid(1) = itmp
    end if

    ! Define variable

    ndims = 0
    if (present(dim1name)) then
       do n = 1, size(dimid)
          if (dimid(n) /= 0) ndims = ndims + 1
       end do
    end if

    call ncd_defvar_bynf(ncid,varname,xtype,ndims,dimid,varid)

    if (present(long_name)) then
       call ncd_putatt(ncid, varid, 'long_name', trim(long_name))
    end if
    if (present(units)) then
       call ncd_putatt(ncid, varid, 'units', trim(units))
    end if
    if (present(cell_method)) then
       str = 'time: ' // trim(cell_method)
       call ncd_putatt(ncid, varid, 'cell_methods', trim(str))
    end if
    if (present(fill_value)) then
       call ncd_putatt(ncid, varid, '_FillValue', fill_value, xtype)
    end if
    if (present(missing_value)) then
       call ncd_putatt(ncid, varid, 'missing_value', missing_value, xtype)
    end if
    if (present(ifill_value)) then
       call ncd_putatt(ncid, varid, '_FillValue', ifill_value, xtype)
    end if
    if (present(imissing_value)) then
       call ncd_putatt(ncid, varid, 'missing_value', imissing_value, xtype)
    end if

  end subroutine ncd_defvar_bygrid

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_int_var0_nf
!
! !INTERFACE:
  subroutine ncd_io_int_var0_nf(varname, data, flag, ncid, readvar, nt)
!
! !DESCRIPTION:
! netcdf I/O of global var int array
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid      ! netcdf file id
    character(len=*) , intent(in)    :: flag      ! 'read' or 'write'
    character(len=*) , intent(in)    :: varname   ! variable name
    integer          , intent(inout) :: data      ! raw data
    logical, optional, intent(out)   :: readvar   ! was var read?
    integer, optional, intent(in)    :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                ! netCDF variable id
    integer :: start(1), count(1)   ! output bounds
    integer :: status               ! error code
    logical :: varpresent           ! if true, variable is on tape
    integer :: temp(1)              ! temporary
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_int_var0_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             call scam_field_offsets(ncid,'undefined',start,count)
             status = pio_get_var(ncid, varid, start, count, temp)
	     data = temp(1)
          else
             status = pio_get_var(ncid, varid, data)
          endif
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       if (present(nt))	then
          start(1) = nt
          count(1) = 1
       else
          start(1) = 1
          count(1) = 1
       end if
       call ncd_inqvid  (ncid, varname, varid, vardesc)
       temp(1) = data
       status = pio_put_var(ncid, varid, start, count, temp)

    endif   ! flag

  end subroutine ncd_io_int_var0_nf

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_real_var0_nf
!
! !INTERFACE:
  subroutine ncd_io_real_var0_nf(varname, data, flag, ncid, readvar, nt)
!
! !DESCRIPTION:
! netcdf I/O of global var int array
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid      ! netcdf file id
    character(len=*) , intent(in)    :: flag      ! 'read' or 'write'
    character(len=*) , intent(in)    :: varname   ! variable name
    real(r8)         , intent(inout) :: data      ! raw data
    logical, optional, intent(out)   :: readvar   ! was var read?
    integer, optional, intent(in)    :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                ! netCDF variable id
    integer :: start(1), count(1)   ! output bounds
    integer :: status               ! error code
    logical :: varpresent           ! if true, variable is on tape
    real(r8):: temp(1)              ! temporary          	
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_real_var0_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             call scam_field_offsets(ncid, 'undefined', start, count)
             status = pio_get_var(ncid, vardesc, start, count, temp)
	     data = temp(1)
          else
             status = pio_get_var(ncid, vardesc, data)
          endif
       endif
       if (present(readvar)) readvar = varpresent

    else if (flag == 'write') then

       if (present(nt))	then
          start(1) = nt
          count(1) = 1
       else
          start(1) = 1
          count(1) = 1
       end if
       call ncd_inqvid  (ncid, varname, varid, vardesc)
       temp(1) = data
       status = pio_put_var(ncid, varid, start, count, temp)

    endif   ! flag

  end subroutine ncd_io_real_var0_nf

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_int_var1_nf
!
! !INTERFACE:
  subroutine ncd_io_int_var1_nf(varname, data, flag, ncid, readvar, nt)
!
! !DESCRIPTION:
! netcdf I/O of global var int array
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid      ! netcdf file id
    character(len=*) , intent(in)    :: flag      ! 'read' or 'write'
    character(len=*) , intent(in)    :: varname   ! variable name
    integer          , intent(inout) :: data(:)   ! raw data
    logical, optional, intent(out)   :: readvar   ! was var read?
    integer, optional, intent(in)    :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                ! netCDF variable id
    integer :: start(2), count(2)   ! output bounds
    integer :: status               ! error code
    logical :: varpresent           ! if true, variable is on tape
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_int_var1_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             call scam_field_offsets(ncid,'undefined',start,count)
             status = pio_get_var(ncid, varid, start, count, data)
          else
             status = pio_get_var(ncid, varid, data)
          endif
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       if (present(nt))	then
          start(1) = 1
          count(1) = size(data)
	  start(2) = nt
	  count(2) = 1
       else
          start(1) = 1
          count(1) = size(data)
	  start(2) = 1
	  count(2) = 1
       end if
       call ncd_inqvid  (ncid, varname, varid, vardesc)
       status = pio_put_var(ncid, varid, start, count, data)

    endif   ! flag

  end subroutine ncd_io_int_var1_nf

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_real_var1_nf
!
! !INTERFACE:
  subroutine ncd_io_real_var1_nf(varname, data, flag, ncid, readvar, nt)
!
! !DESCRIPTION:
! netcdf I/O of global var int array
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid             ! netcdf file id
    character(len=*) , intent(in)    :: flag             ! 'read' or 'write'
    character(len=*) , intent(in)    :: varname          ! variable name
    real(r8)         , intent(inout) :: data(:)          ! raw data
    logical          , optional, intent(out):: readvar   ! was var read?
    integer          , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                ! netCDF variable id
    integer :: start(2), count(2)   ! output bounds
    integer :: status               ! error code
    logical :: varpresent           ! if true, variable is on tape
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_real_var1_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             call scam_field_offsets(ncid,'undefined',start,count)
             status = pio_get_var(ncid, varid, start, count, data)
          else
             status = pio_get_var(ncid, varid, data)
          endif
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       if (present(nt))	then
          start(1) = 1
	  start(2) = nt
          count(1) = size(data)
	  count(2) = 1
       else
          start(1) = 1
	  start(2) = 1
          count(1) = size(data)
	  count(2) = 1
       end if
       call ncd_inqvid  (ncid, varname, varid, vardesc)
       status = pio_put_var(ncid, varid, start, count, data)

    endif   ! flag

  end subroutine ncd_io_real_var1_nf

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_int_var2_nf
!
! !INTERFACE:
  subroutine ncd_io_int_var2_nf(varname, data, flag, ncid, readvar, nt)
!
! !DESCRIPTION:
! netcdf I/O of global var int array
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid             ! netcdf file id
    character(len=*) , intent(in)    :: flag             ! 'read' or 'write'
    character(len=*) , intent(in)    :: varname          ! variable name
    integer          , intent(inout) :: data(:,:)        ! raw data
    logical          , optional, intent(out):: readvar   ! was var read?
    integer          , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                ! netCDF variable id
    integer :: start(3), count(3)   ! output bounds
    integer :: status               ! error code
    logical :: varpresent           ! if true, variable is on tape
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_int_var2_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             call scam_field_offsets(ncid,'undefined',start,count)
             status = pio_get_var(ncid, varid, start, count, data)
          else
             status = pio_get_var(ncid, varid, data)
          endif
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       if (present(nt))	then
          start(1) = 1
	  start(2) = 1
	  start(3) = nt
          count(1) = size(data, dim=1)
	  count(2) = size(data, dim=2)
	  count(3) = 1
       else
          start(1) = 1
	  start(2) = 1
	  start(3) = 1
          count(1) = size(data, dim=1)
	  count(2) = size(data, dim=2)
	  count(3) = 1
       end if
       call ncd_inqvid(ncid, varname, varid, vardesc)
       status = pio_put_var(ncid, varid, start, count, data)

    endif   

  end subroutine ncd_io_int_var2_nf

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_real_var2_nf
!
! !INTERFACE:
  subroutine ncd_io_real_var2_nf(varname, data, flag, ncid, readvar, nt)
!
! !DESCRIPTION:
! netcdf I/O of global var int array
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid             ! netcdf file id
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    character(len=*), intent(in)    :: varname          ! variable name
    real(r8)        , intent(inout) :: data(:,:)        ! raw data
    logical         , optional, intent(out):: readvar   ! was var read?
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                ! netCDF variable id
    integer :: start(3), count(3)   ! output bounds
    integer :: status               ! error code
    logical :: varpresent           ! if true, variable is on tape
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_real_var2_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             call scam_field_offsets(ncid,'undefined',start,count)
             status = pio_get_var(ncid, varid, start, count, data)
          else
             status = pio_get_var(ncid, varid, data)
          endif
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       if (present(nt))	then
          start(1) = 1
	  start(2) = 1
	  start(3) = nt
          count(1) = size(data, dim=1)
	  count(2) = size(data, dim=2)
	  count(3) = 1
       else
          start(1) = 1
	  start(2) = 1
	  start(3) = 1
          count(1) = size(data, dim=1)
	  count(2) = size(data, dim=2)
	  count(3) = 1
       end if
       call ncd_inqvid  (ncid, varname, varid, vardesc)
       status = pio_put_var(ncid, varid, start, count, data)

    endif   

  end subroutine ncd_io_real_var2_nf

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_int_var1
!
! !INTERFACE:
  subroutine ncd_io_int_var1(varname, data, dim1name, &
                             flag, ncid, nt, readvar)
!
! !DESCRIPTION:
! I/O for 1d int field
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid             ! netcdf file id
    character(len=*) , intent(in)    :: flag             ! 'read' or 'write'
    character(len=*) , intent(in)    :: varname          ! variable name
    integer          , pointer       :: data(:)          ! local decomposition data
    character(len=*) , intent(in)    :: dim1name         ! dimension name
    integer          , optional, intent(in) :: nt        ! time sample index
    logical          , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
!EOP
    character(len=8)  :: clmlevel   ! clmlevel
    character(len=32) :: dimname    ! temporary
    integer           :: n          ! index	
    integer           :: iodnum     ! iodesc num in list
    integer           :: varid      ! varid
    integer           :: ndims      ! ndims for var
    integer           :: ndims_iod  ! ndims iodesc for var
    integer           :: dims(4)    ! dim sizes	 
    integer           :: dids(4)    ! dim ids
    integer           :: start(3)   ! netcdf start index
    integer           :: count(3)   ! netcdf count index
    integer           :: status     ! error code  
    logical           :: varpresent ! if true, variable is on tape
    integer           :: xtype      ! netcdf data type
    integer                , pointer  :: compDOF(:)
    type(iodesc_plus_type) , pointer  :: iodesc_plus
    type(var_desc_t)                  :: vardesc
    character(len=*),parameter :: subname='ncd_io_int_var1' ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(flag),' ',trim(varname),' ',trim(clmlevel)
    end if

    clmlevel = dim1name

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             start(:) = 1
             count(:) = 1
             call scam_field_offsets(ncid,clmlevel,start,count)
             if (trim(clmlevel) == gratm .or. trim(clmlevel) == grlnd) then
                if (present(nt)) start(3) = nt
             else
                if (present(nt)) start(2) = nt
             end if
             status = pio_get_var(ncid, varid, start, count, data)
          else
             status = pio_inq_varndims(ncid, vardesc, ndims)
	     status = pio_inq_vardimid(ncid, vardesc, dids)
             status = pio_inq_vartype (ncid, vardesc, xtype)
	     status = pio_inq_dimname(ncid,dids(ndims),dimname)
             if ('time' == trim(dimname)) then
                ndims_iod = ndims - 1
             else
                ndims_iod = ndims
             end if
             do n = 1,ndims_iod
                status = pio_inq_dimlen(ncid,dids(n),dims(n))
             enddo
             call ncd_getiodesc(clmlevel, ndims_iod, dims(1:ndims_iod), xtype, iodnum)
             iodesc_plus => iodesc_list(iodnum)
             if (present(nt)) then
                call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
             end if
             call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, data, status)
          end if
       end if
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       call ncd_inqvid(ncid, varname ,varid, vardesc)
       status = pio_inq_varndims(ncid, vardesc, ndims)
       status = pio_inq_vardimid(ncid, vardesc, dids)
       status = pio_inq_vartype (ncid, vardesc, xtype)
       status = pio_inq_dimname(ncid,dids(ndims),dimname)
       if ('time' == trim(dimname)) then
          ndims_iod = ndims - 1
       else
          ndims_iod = ndims
       end if
       do n = 1,ndims_iod
          status = pio_inq_dimlen(ncid,dids(n),dims(n))
       enddo
       call ncd_getiodesc(clmlevel, ndims_iod, dims(1:ndims_iod), xtype, iodnum)
       iodesc_plus => iodesc_list(iodnum)
       if (present(nt)) then
          call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
       end if
       call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, data, status, fillval=0)

    else

       if (masterproc) then
          write(iulog,*) subname,' error: unsupported flag ',trim(flag)
          call endrun()
       endif

    endif

  end subroutine ncd_io_int_var1

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_real_var1
!
! !INTERFACE:
  subroutine ncd_io_real_var1(varname, data, dim1name, &
                              flag, ncid, nt, readvar)
!
! !DESCRIPTION:
! I/O for 1d int field
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid             ! netcdf file id
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
    character(len=*), intent(in)  :: varname            ! variable name
    real(r8)        , pointer     :: data(:)            ! local decomposition data
    character(len=*), intent(in)  :: dim1name           ! dimension name
    integer         , optional, intent(in) :: nt        ! time sample index
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=8)  :: clmlevel   ! clmlevel
    character(len=32) :: dimname    ! temporary
    integer           :: iodnum     ! iodesc num in list
    integer           :: varid      ! varid
    integer           :: ndims      ! ndims for var
    integer           :: ndims_iod  ! ndims iodesc for var
    integer           :: n          ! index	
    integer           :: dims(4)    ! dim sizes	 
    integer           :: dids(4)    ! dim ids
    integer           :: start(3)   ! netcdf start index
    integer           :: count(3)   ! netcdf count index
    integer           :: status     ! error code  
    logical           :: varpresent ! if true, variable is on tape
    integer           :: xtype      ! netcdf data type
    integer                , pointer  :: compDOF(:)
    type(iodesc_plus_type) , pointer  :: iodesc_plus
    type(var_desc_t)                  :: vardesc
    character(len=*),parameter :: subname='ncd_io_real_var1' ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(flag),' ',trim(varname),' ',trim(clmlevel)
    endif

    clmlevel = dim1name

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             start(:) = 1
             count(:) = 1
             call scam_field_offsets(ncid,clmlevel,start,count)
             if (trim(clmlevel) == gratm .or. trim(clmlevel) == grlnd) then
                if (present(nt)) start(3) = nt
             else
                if (present(nt)) start(2) = nt
             end if
             status = pio_get_var(ncid, varid, start, count, data)
          else
             status = pio_inq_varndims(ncid, vardesc, ndims)
             status = pio_inq_vardimid(ncid,vardesc, dids)
             status = pio_inq_vartype(ncid, vardesc, xtype)
             status = pio_inq_dimname(ncid,dids(ndims),dimname)
             if ('time' == trim(dimname)) then
                ndims_iod = ndims - 1
             else
                ndims_iod = ndims
             end if
             do n = 1,ndims_iod
                status = pio_inq_dimlen(ncid,dids(n),dims(n))
             enddo
             call ncd_getiodesc(clmlevel, ndims_iod, dims(1:ndims_iod), xtype, iodnum)
             iodesc_plus => iodesc_list(iodnum)
             if (present(nt)) then
                call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
             end if
             call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, data, status)
          end if
       end if
       if (present(readvar)) readvar = varpresent
       
    elseif (flag == 'write') then
       
       call ncd_inqvid(ncid, varname ,varid, vardesc)
       status = pio_inq_varndims(ncid, vardesc, ndims)
       status = pio_inq_vardimid(ncid, vardesc, dids)
       status = pio_inq_vartype (ncid, vardesc, xtype)
       status = pio_inq_dimname(ncid,dids(ndims),dimname)
       if ('time' == trim(dimname)) then
          ndims_iod = ndims - 1
       else
          ndims_iod = ndims
       end if
       do n = 1,ndims_iod
          status = pio_inq_dimlen(ncid,dids(n),dims(n))
       enddo
       call ncd_getiodesc(clmlevel, ndims_iod, dims(1:ndims_iod), xtype, iodnum)
       iodesc_plus => iodesc_list(iodnum)
       if (present(nt)) then
          call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
       end if
       call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, data, status, fillval=spval)
       
    else
       
       if (masterproc) then
          write(iulog,*) subname,' error: unsupported flag ',trim(flag)
          call endrun()
       endif
       
    endif

  end subroutine ncd_io_real_var1

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_int_var2
!
! !INTERFACE:
  subroutine ncd_io_int_var2(varname, data, dim1name, &
	                   lowerb2, upperb2, flag, ncid, nt, readvar, switchdim)
!
! !DESCRIPTION:
! Netcdf i/o of 2d initial integer field out to netCDF file
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid          ! netcdf file id
    character(len=*) , intent(in)  :: flag            ! 'read' or 'write'
    character(len=*) , intent(in)  :: varname         ! variable name
    integer          , pointer     :: data(:,:)       ! local decomposition input data
    character(len=*) , intent(in)  :: dim1name        ! dimension 1 name
    integer, optional, intent(in)  :: nt              ! time sample index
    integer, optional, intent(in)  :: lowerb2,upperb2 ! lower and upper bounds of second dimension
    logical, optional, intent(out) :: readvar         ! true => variable is on initial dataset (read only)
    logical, optional, intent(in)  :: switchdim	      ! true=> permute dim1 and dim2 for output
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer, pointer  :: temp(:,:)
    integer           :: ndim1,ndim2 	
    character(len=8)  :: clmlevel   ! clmlevel
    character(len=32) :: dimname    ! temporary	
    integer           :: status     ! error status
    integer           :: ndims      ! ndims total for var
    integer           :: ndims_iod  ! ndims iodesc for var
    integer           :: varid      ! varid
    integer           :: n          ! index	
    integer           :: dims(4)    ! dim sizes	 
    integer           :: dids(4)    ! dim ids
    integer           :: iodnum     ! iodesc num in list
    integer           :: start(4)   ! netcdf start index
    integer           :: count(4)   ! netcdf count index
    logical           :: varpresent ! if true, variable is on tape
    integer           :: xtype      ! netcdf data type
    integer           :: i,j	
    integer           :: lb1,lb2
    integer           :: ub1,ub2
    type(iodesc_plus_type) , pointer  :: iodesc_plus
    type(var_desc_t)                  :: vardesc
    character(len=*),parameter :: subname='ncd_io_int_var2' ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(flag),' ',trim(varname),' ',trim(clmlevel)
    end if

    clmlevel = dim1name	

    if (present(switchdim)) then
       lb1 = lbound(data, dim=1)
       ub1 = ubound(data, dim=1)
       lb2 = lbound(data, dim=2)
       ub2 = ubound(data, dim=2)
       if (present(lowerb2)) lb2 = lowerb2
       if (present(upperb2)) ub2 = upperb2
       allocate(temp(lb2:ub2,lb1:ub1))
    end if

    if (flag == 'read') then
       
       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             start(:) = 1
             count(:) = 1
             call scam_field_offsets(ncid, clmlevel, start, count)
             status = pio_get_var(ncid, vardesc, start, count, data)
             if (trim(clmlevel) == gratm .or. trim(clmlevel) == grlnd) then
                count(3) = size(data,dim=2)
                if (present(nt)) start(4) = nt
             else
                count(2) = size(data,dim=2)
                if (present(nt)) start(3) = nt
             end if
          else
             status = pio_inq_varndims(ncid, vardesc, ndims)
             status = pio_inq_vardimid(ncid, vardesc, dids)
             status = pio_inq_vartype (ncid, vardesc, xtype)
             status = pio_inq_dimname(ncid, dids(ndims), dimname)
             if (ndims == 0) then
                write(iulog,*) trim(subname),' ERROR: ndims must be greater than 0'
                call endrun()
             end if
             if ('time' == trim(dimname)) then
                ndims_iod = ndims - 1
             else
                ndims_iod = ndims
             end if
             do n = 1,ndims_iod
                status = pio_inq_dimlen(ncid,dids(n),dims(n))
             enddo
             if (present(switchdim)) then
                call ncd_getiodesc(clmlevel, ndims_iod, dims(1:ndims_iod), xtype, iodnum, switchdim=.true.)
             else
                call ncd_getiodesc(clmlevel, ndims_iod, dims(1:ndims_iod), xtype, iodnum)
             end if
             iodesc_plus => iodesc_list(iodnum)
             if (present(nt)) then
                call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
             end if
             if (present(switchdim)) then
                call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, temp, status)
                do j = lb2,ub2
                do i = lb1,ub1
                   data(i,j) = temp(j,i) 
                end do
                end do
             else
                call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, data, status)
             end if
          end if
       end if
       if (present(readvar)) readvar = varpresent

    else if (flag == 'write') then

       call ncd_inqvid(ncid, varname ,varid, vardesc)
       status = pio_inq_varndims(ncid, vardesc, ndims)
       status = pio_inq_vardimid(ncid, vardesc , dids)
       status = pio_inq_vartype (ncid, vardesc, xtype)
       if (ndims == 0) then
          write(iulog,*) trim(subname),' ERROR: ndims must be greater than 0'
          call endrun()
       end if
       status = pio_inq_dimname(ncid,dids(ndims),dimname)
       if ('time' == trim(dimname)) then
          ndims_iod = ndims - 1
       else
          ndims_iod = ndims
       end if
       do n = 1,ndims_iod
          status = pio_inq_dimlen(ncid,dids(n),dims(n))
       enddo
       if (present(switchdim)) then
          call ncd_getiodesc(clmlevel, ndims_iod, dims(1:ndims_iod), xtype, iodnum, switchdim=.true.)
       else
          call ncd_getiodesc(clmlevel, ndims_iod, dims(1:ndims_iod), xtype, iodnum)
       end if
       iodesc_plus => iodesc_list(iodnum)
       if (present(nt)) then
          call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
       end if
       if (present(switchdim)) then
	  do j = lb2,ub2
	  do i = lb1,ub1
             temp(j,i) = data(i,j)
          end do
          end do
          call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, temp, status, fillval=0)
       else
          call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, data, status, fillval=0)
       end if

    else
       
       if (masterproc) then
          write(iulog,*) subname,' error: unsupported flag ',trim(flag)
          call endrun()
       endif
       
    endif

  end subroutine ncd_io_int_var2

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_real_var2
!
! !INTERFACE:
  subroutine ncd_io_real_var2(varname, data, dim1name, &
   	                      lowerb2, upperb2, flag, ncid, nt, readvar, switchdim)
!
! !DESCRIPTION:
! Netcdf i/o of 2d initial integer field out to netCDF file
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid          ! netcdf file id
    character(len=*) , intent(in)  :: flag            ! 'read' or 'write'
    character(len=*) , intent(in)  :: varname         ! variable name
    real(r8)         , pointer     :: data(:,:)       ! local decomposition input data
    character(len=*) , intent(in)  :: dim1name        ! dimension 1 name
    integer, optional, intent(in)  :: nt              ! time sample index
    integer, optional, intent(in)  :: lowerb2,upperb2 ! lower and upper bounds of second dimension
    logical, optional, intent(out) :: readvar         ! true => variable is on initial dataset (read only)
    logical, optional, intent(in)  :: switchdim	      ! true=> permute dim1 and dim2 for output
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    real(r8), pointer :: temp(:,:)
    integer           :: ndim1,ndim2 	
    character(len=8)  :: clmlevel   ! clmlevel
    character(len=32) :: dimname    ! temporary	
    integer           :: status     ! error status
    integer           :: ndims      ! ndims total for var
    integer           :: ndims_iod  ! ndims iodesc for var
    integer           :: varid      ! varid
    integer           :: n          ! index	
    integer           :: dims(4)    ! dim sizes	 
    integer           :: dids(4)    ! dim ids
    integer           :: iodnum     ! iodesc num in list
    integer           :: start(4)   ! netcdf start index
    integer           :: count(4)   ! netcdf count index
    logical           :: varpresent ! if true, variable is on tape
    integer           :: xtype      ! netcdf data type
    integer           :: i,j	
    integer           :: lb1,lb2
    integer           :: ub1,ub2
    type(iodesc_plus_type) , pointer  :: iodesc_plus
    type(var_desc_t)                  :: vardesc
    character(len=*),parameter :: subname='ncd_io_real_var2' ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(flag),' ',trim(varname),' ',trim(clmlevel)
    end if

    clmlevel = dim1name	

    if (present(switchdim)) then
       lb1 = lbound(data, dim=1)
       ub1 = ubound(data, dim=1)
       lb2 = lbound(data, dim=2)
       ub2 = ubound(data, dim=2)
       if (present(lowerb2)) lb2 = lowerb2
       if (present(upperb2)) ub2 = upperb2
       allocate(temp(lb2:ub2,lb1:ub1))
    end if

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             start(:) = 1
             count(:) = 1
             call scam_field_offsets(ncid, clmlevel, start, count)
             if (trim(clmlevel) == gratm .or. trim(clmlevel) == grlnd) then
                count(3) = size(data,dim=2)
                if (present(nt)) start(4) = nt
             else
                count(2) = size(data,dim=2)
                if (present(nt)) start(3) = nt
             end if
             status = pio_get_var(ncid, vardesc, start, count, data)
          else
             status = pio_inq_varndims(ncid, vardesc, ndims)
             status = pio_inq_vartype(ncid, vardesc, xtype)
	     status = pio_inq_vardimid(ncid,vardesc, dids)
             status = pio_inq_dimname(ncid, dids(ndims), dimname)
             if (ndims == 0) then
                write(iulog,*) trim(subname),' ERROR: ndims must be greater than 0'
                call endrun()
             end if
             if ('time' == trim(dimname)) then
                ndims_iod = ndims - 1
             else
                ndims_iod = ndims
             end if
             do n = 1,ndims_iod
                status = pio_inq_dimlen(ncid,dids(n),dims(n))
             enddo
             if (present(switchdim)) then
                call ncd_getiodesc(clmlevel, ndims_iod, dims(1:ndims_iod), xtype, iodnum, switchdim=.true.)
             else
                call ncd_getiodesc(clmlevel, ndims_iod, dims(1:ndims_iod), xtype, iodnum)
             end if
             iodesc_plus => iodesc_list(iodnum)
             if (present(nt)) then
                call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
             end if
             if (present(switchdim)) then
                call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, temp, status)
                do j = lb2,ub2
                do i = lb1,ub1
                   data(i,j) = temp(j,i) 
                end do
                end do
             else
                call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, data, status)
             end if
          end if
       end if
       if (present(readvar)) readvar = varpresent

    else if (flag == 'write') then

       call ncd_inqvid(ncid, varname ,varid, vardesc)
       status = pio_inq_varndims(ncid, vardesc, ndims)
       status = pio_inq_vardimid(ncid, vardesc , dids)
       status = pio_inq_vartype (ncid, vardesc, xtype)
       if (ndims == 0) then
          write(iulog,*) trim(subname),' ERROR: ndims must be greater than 0'
          call endrun()
       end if
       status = pio_inq_dimname(ncid,dids(ndims),dimname)
       if ('time' == trim(dimname)) then
          ndims_iod = ndims - 1
       else
          ndims_iod = ndims
       end if
       do n = 1,ndims_iod
          status = pio_inq_dimlen(ncid,dids(n),dims(n))
       enddo
       if (present(switchdim)) then
          call ncd_getiodesc(clmlevel, ndims_iod, dims(1:ndims_iod), xtype, iodnum, switchdim=.true.)
       else
          call ncd_getiodesc(clmlevel, ndims_iod, dims(1:ndims_iod), xtype, iodnum)
       end if
       iodesc_plus => iodesc_list(iodnum)
       if (present(nt)) then
          call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
       end if
       if (present(switchdim)) then
	  do j = lb2,ub2
	  do i = lb1,ub1
             temp(j,i) = data(i,j)
          end do
          end do
          call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, temp, status, fillval=spval)
       else
          call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, data, status, fillval=spval)
       end if

    else

       if (masterproc) then
          write(iulog,*) subname,' error: unsupported flag ',trim(flag)
          call endrun()
       endif

    endif

    if (present(switchdim)) then
       deallocate(temp)
    end if

  end subroutine ncd_io_real_var2

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_int_var3
!
! !INTERFACE:
  subroutine ncd_io_int_var3(varname, data, dim1name, &
   	                      flag, ncid, nt, readvar)
!
! !DESCRIPTION:
! Netcdf i/o of 2d initial integer field out to netCDF file
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid          ! netcdf file id
    character(len=*) , intent(in)  :: flag            ! 'read' or 'write'
    character(len=*) , intent(in)  :: varname         ! variable name
    integer          , pointer     :: data(:,:,:)     ! local decomposition input data
    character(len=*) , intent(in)  :: dim1name        ! dimension 1 name
    integer, optional, intent(in)  :: nt              ! time sample index
    logical, optional, intent(out) :: readvar         ! true => variable is on initial dataset (read only)
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer           :: ndim1,ndim2 	
    character(len=8)  :: clmlevel   ! clmlevel
    character(len=32) :: dimname    ! temporary	
    integer           :: status     ! error status
    integer           :: ndims      ! ndims total for var
    integer           :: ndims_iod  ! ndims iodesc for var
    integer           :: varid      ! varid
    integer           :: n          ! index	
    integer           :: dims(4)    ! dim sizes	 
    integer           :: dids(4)    ! dim ids
    integer           :: iodnum     ! iodesc num in list
    integer           :: start(5)   ! netcdf start index
    integer           :: count(5)   ! netcdf count index
    logical           :: varpresent ! if true, variable is on tape
    integer           :: xtype      ! netcdf data type
    type(iodesc_plus_type) , pointer  :: iodesc_plus
    type(var_desc_t)                  :: vardesc
    character(len=*),parameter :: subname='ncd_io_int_var3' ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(flag),' ',trim(varname),' ',trim(clmlevel)
    end if

    clmlevel = dim1name	

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             start(:) = 1
             count(:) = 1
             call scam_field_offsets(ncid, clmlevel, start, count)
             if (trim(clmlevel) == gratm .or. trim(clmlevel) == grlnd) then
                count(3) = size(data,dim=2)
                count(4) = size(data,dim=3)
                if (present(nt)) start(5) = nt
             else
                count(2) = size(data,dim=2)
                count(3) = size(data,dim=3)
                if (present(nt)) start(4) = nt
             end if
             status = pio_get_var(ncid, vardesc, start, count, data)
          else
             status = pio_inq_varndims(ncid, vardesc, ndims)
	     status = pio_inq_vardimid(ncid,vardesc, dids)
             status = pio_inq_vartype(ncid, vardesc, xtype)
             status = pio_inq_dimname(ncid, dids(ndims), dimname)
             if (ndims == 0) then
                write(iulog,*) trim(subname),' ERROR: ndims must be greater than 0'
                call endrun()
             end if
             if ('time' == trim(dimname)) then
                ndims_iod = ndims - 1
             else
                ndims_iod = ndims
             end if
             do n = 1,ndims_iod
                status = pio_inq_dimlen(ncid,dids(n),dims(n))
             enddo
             call ncd_getiodesc(clmlevel, ndims_iod, dims(1:ndims_iod), xtype, iodnum)
             iodesc_plus => iodesc_list(iodnum)
             if (present(nt)) then
                call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
             end if
             call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, data, status)
          end if
       end if
       if (present(readvar)) readvar = varpresent

    else if (flag == 'write') then

       call ncd_inqvid(ncid, varname ,varid, vardesc)
       status = pio_inq_varndims(ncid, vardesc, ndims)
       status = pio_inq_vardimid(ncid, vardesc , dids)
       status = pio_inq_vartype (ncid, vardesc, xtype)
       if (ndims == 0) then
          write(iulog,*) trim(subname),' ERROR: ndims must be greater than 0'
          call endrun()
       end if
       status = pio_inq_dimname(ncid,dids(ndims),dimname)
       if ('time' == trim(dimname)) then
          ndims_iod = ndims - 1
       else
          ndims_iod = ndims
       end if
       do n = 1,ndims_iod
          status = pio_inq_dimlen(ncid,dids(n),dims(n))
       enddo
       call ncd_getiodesc(clmlevel, ndims_iod, dims(1:ndims_iod), xtype, iodnum)
       iodesc_plus => iodesc_list(iodnum)
       if (present(nt)) then
          call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
       end if
       call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, data, status)

    else

       if (masterproc) then
          write(iulog,*) subname,' error: unsupported flag ',trim(flag)
          call endrun()
       endif

    endif

  end subroutine ncd_io_int_var3

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_real_var3
!
! !INTERFACE:
  subroutine ncd_io_real_var3(varname, data, dim1name, &
   	                      flag, ncid, nt, readvar)
!
! !DESCRIPTION:
! Netcdf i/o of 2d initial integer field out to netCDF file
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid          ! netcdf file id
    character(len=*) , intent(in)  :: flag            ! 'read' or 'write'
    character(len=*) , intent(in)  :: varname         ! variable name
    real(r8)         , pointer     :: data(:,:,:)     ! local decomposition input data
    character(len=*) , intent(in)  :: dim1name        ! dimension 1 name
    integer, optional, intent(in)  :: nt              ! time sample index
    logical, optional, intent(out) :: readvar         ! true => variable is on initial dataset (read only)
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    real(r8), allocatable :: temp(:,:)
    integer               :: ndim1,ndim2 	
    character(len=8)  :: clmlevel   ! clmlevel
    character(len=32) :: dimname    ! temporary	
    integer           :: status     ! error status
    integer           :: ndims      ! ndims total for var
    integer           :: ndims_iod  ! ndims iodesc for var
    integer           :: varid      ! varid
    integer           :: n          ! index	
    integer           :: dims(4)    ! dim sizes	 
    integer           :: dids(4)    ! dim ids
    integer           :: iodnum     ! iodesc num in list
    integer           :: start(5)   ! netcdf start index
    integer           :: count(5)   ! netcdf count index
    logical           :: varpresent ! if true, variable is on tape
    integer           :: xtype      ! netcdf data type
    type(iodesc_plus_type) , pointer  :: iodesc_plus
    type(var_desc_t)                  :: vardesc
    character(len=*),parameter :: subname='ncd_io_real_var3' ! subroutine name
!-----------------------------------------------------------------------


    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(flag),' ',trim(varname),' ',trim(clmlevel)
    end if

    clmlevel = dim1name	

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             start(:) = 1
             count(:) = 1
             call scam_field_offsets(ncid, clmlevel, start, count)
             if (trim(clmlevel) == gratm .or. trim(clmlevel) == grlnd) then
                count(3) = size(data,dim=2)
                count(4) = size(data,dim=3)
                if (present(nt)) start(5) = nt
             else
                count(2) = size(data,dim=2)
                count(3) = size(data,dim=3)
                if (present(nt)) start(4) = nt
             end if
             status = pio_get_var(ncid, vardesc, start, count, data)
          else
             status = pio_inq_varndims(ncid, vardesc, ndims)
	     status = pio_inq_vardimid(ncid, vardesc, dids)
             status = pio_inq_vartype(ncid, vardesc, xtype)
             status = pio_inq_dimname(ncid, dids(ndims), dimname)
             if (ndims == 0) then
                write(iulog,*) trim(subname),' ERROR: ndims must be greater than 0'
                call endrun()
             end if
             if ('time' == trim(dimname)) then
                ndims_iod = ndims - 1
             else
                ndims_iod = ndims
             end if
             do n = 1,ndims_iod
                status = pio_inq_dimlen(ncid,dids(n),dims(n))
             enddo
             call ncd_getiodesc(clmlevel, ndims_iod, dims(1:ndims_iod), xtype, iodnum)
             iodesc_plus => iodesc_list(iodnum)
             if (present(nt)) then
                call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
             end if
             call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, data, status)
          end if
       end if
       if (present(readvar)) readvar = varpresent

    else if (flag == 'write') then

       call ncd_inqvid(ncid, varname ,varid, vardesc)
       status = pio_inq_varndims(ncid, vardesc, ndims)
       status = pio_inq_vardimid(ncid, vardesc , dids)
       status = pio_inq_vartype (ncid, vardesc, xtype)
       if (ndims == 0) then
          write(iulog,*) trim(subname),' ERROR: ndims must be greater than 0'
          call endrun()
       end if
       status = pio_inq_dimname(ncid,dids(ndims),dimname)
       if ('time' == trim(dimname)) then
          ndims_iod = ndims - 1
       else
          ndims_iod = ndims
       end if
       do n = 1,ndims_iod
          status = pio_inq_dimlen(ncid,dids(n),dims(n))
       enddo
       call ncd_getiodesc(clmlevel, ndims_iod, dims(1:ndims_iod), xtype, iodnum)
       iodesc_plus => iodesc_list(iodnum)
       if (present(nt)) then
          call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
       end if
       call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, data, status, fillval=spval)

    else

       if (masterproc) then
          write(iulog,*) subname,' error: unsupported flag ',trim(flag)
          call endrun()
       endif

    endif

  end subroutine ncd_io_real_var3

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
    use clm_varpar  , only : maxpatch
    use clm_varctl  , only : scmlon,scmlat,single_column
    use nanMod      , only : nan
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t) , intent(inout) :: ncid     ! netcdf file id
    character(len=*)  , intent(in)    :: dim1name ! dimension 1 name
    integer           , intent(inout) :: start(:) ! start index
    integer           , intent(inout) :: count(:) ! count to retrieve
!
! !CALLED FROM: subroutine inicfields
!
! !REVISION HISTORY:
! Created by John Truesdale

!
! !LOCAL VARIABLES:
!EOP
    integer :: data_offset                   ! offset into land array 1st column 
    integer :: ndata                         ! number of column (or pft points to read)
    real(r8) , pointer :: cols1dlon(:)       ! holds cols1d_ixy var
    real(r8) , pointer :: cols1dlat(:)       ! holds cols1d_jxy var
    real(r8) , pointer :: pfts1dlon(:)       ! holds pfts1d_ixy var
    real(r8) , pointer :: pfts1dlat(:)       ! holds pfts1d_jxy var
    integer :: cols(maxpatch)                ! grid cell columns for scam
    integer :: pfts(maxpatch)                ! grid cell pfts for scam
    integer :: cc,i                          ! index variable
    integer :: totpfts                       ! total number of pfts
    integer :: totcols                       ! total number of columns
    integer :: dimid                         ! netCDF dimension id
    integer :: varid                         ! netCDF variable id
    integer :: status                        ! return code
    integer :: latidx,lonidx                 ! latitude/longitude indices
    real(r8):: closelat,closelon             ! closest latitude and longitude indices
    character(len=32) :: subname = 'scam_field_offsets'
!------------------------------------------------------------------------

    ! find closest land grid cell for this point

    call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,latidx,lonidx)
     
    if (dim1name == 'column') then

       status = pio_inq_dimid(ncid, 'column', dimid)
       status = pio_inq_dimlen(ncid, dimid, totcols)

       allocate (cols1dlon(totcols))
       allocate (cols1dlat(totcols))
       
       status = pio_inq_varid(ncid, 'cols1d_lon', varid)
       status = pio_get_var(ncid, varid, cols1dlon)
       status = pio_inq_varid(ncid, 'cols1d_lat', varid)
       status = pio_get_var(ncid, varid, cols1dlat)
       
       cols(:)     = nan
       data_offset = nan
       i = 1
       ndata = 0
       do cc = 1, totcols
          if (cols1dlon(cc) == closelon.and.cols1dlat(cc) == closelat) then
             cols(i)=cc
             ndata  =i
             i=i+1
          end if
       end do
       if (ndata == 0) then
          write(iulog,*)'couldnt find any columns for this latitude ',latidx,' and longitude ',lonidx
          call endrun
       else
          data_offset=cols(1)
       end if
       
       deallocate (cols1dlon)
       deallocate (cols1dlat)
       
       start(1) = data_offset
       count(1) = ndata
       
    else if (dim1name == 'pft') then
       
       status = pio_inq_dimid(ncid, 'pft', dimid)
       status = pio_inq_dimlen(ncid, dimid, totpfts)

       allocate (pfts1dlon(totpfts))
       allocate (pfts1dlat(totpfts))

       status = pio_inq_varid(ncid, 'pfts1d_lon', varid)
       status = pio_get_var(ncid, varid, pfts1dlon)

       status = pio_inq_varid(ncid, 'pfts1d_lat', varid)
       status = pio_get_var(ncid, varid, pfts1dlat)
       
       pfts(:)     = nan
       data_offset = nan
       i     = 1
       ndata = 0
       do cc = 1, totpfts
          if (pfts1dlon(cc) == closelon.and.pfts1dlat(cc) == closelat) then
             pfts(i)=cc
             ndata  =i
             i=i+1
          end if
       end do
       if (ndata == 0) then
          write(iulog,*)'couldnt find any pfts for this latitude ',closelat,' and longitude ',closelon
          call endrun
       else
          data_offset=pfts(1)
       end if
       
       deallocate (pfts1dlon)
       deallocate (pfts1dlat)
       
       start(1) = data_offset
       count(1) = ndata
       
    else
       
       start(1) = lonidx
       count(1) = 1
       start(2) = latidx
       count(2) = 1
       write(iulog,*) trim(subname),' scam_setlatlonidx ',lonidx,latidx
       
    endif

  end subroutine scam_field_offsets

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine ncd_getiodesc
!
! !INTERFACE:
  subroutine ncd_getiodesc(clmlevel, ndims, dims, xtype, iodnum, switchdim)
!
! !DESCRIPTION: 
! Returns an index to an io descriptor
!
! !USES:
!
! !ARGUMENTS:
    character(len=8), intent(in)  :: clmlevel   ! clmlevel
    integer         , intent(in)  :: ndims      ! ndims for var	
    integer         , intent(in)  :: dims(:)    ! dim sizes
    integer         , intent(in)  :: xtype      ! file external type
    integer         , intent(out) :: iodnum     ! iodesc num in list
    logical,optional, intent(in)  :: switchdim  ! switch level dimension and first dim 

! !REVISION HISTORY:
! Created by Mariana Vertenstein

!
! !LOCAL VARIABLES:
!EOP
    integer :: k,m,n,cnt                     ! indices
    integer :: basetype                      ! pio basetype
    integer :: gsmap_lsize                   ! local size of gsmap
    integer :: gsmap_gsize                   ! global size of gsmap
    integer :: fullsize                      ! size of entire array on cdf
    integer :: gsize                         ! global size of clmlevel
    integer :: vsize                         ! other dimensions
    integer :: vsize1, vsize2                ! other dimensions
    integer :: status                        ! error status
    type(mct_gsMap),pointer       :: gsmap   ! global seg map
    integer, pointer,dimension(:) :: gsmOP   ! gsmap ordered points
    integer, pointer  :: compDOF(:)
    character(len=32) :: subname = 'ncd_getiodesc'
!------------------------------------------------------------------------

    call ncd_inqiodesc(ndims, dims, xtype, iodnum)

    if (iodnum < 1) then

       if (ndims > 0) then 
          num_iodesc = num_iodesc + 1
          if (num_iodesc > max_iodesc) then
             write(iulog,*) trim(subname),' ERROR num_iodesc gt max_iodesc ',max_iodesc
             call endrun()
          endif
          iodnum = num_iodesc
          if (masterproc .and. debug > 1) then
             write(iulog,*) trim(subname),' creating iodesc at iodnum,ndims,dims(1:ndims),xtype',&
                  iodnum,ndims,dims(1:ndims),xtype
          endif
       end if

       if (xtype == nf_double ) then
          basetype = PIO_DOUBLE
       else if (xtype == nf_float) then
          basetype  = PIO_DOUBLE
       else if (xtype == nf_int) then
          basetype = PIO_INT
       end if

       call get_clmlevel_gsmap(clmlevel,gsmap)
       gsize = get_clmlevel_gsize(clmlevel)
       gsmap_lsize = mct_gsmap_lsize(gsmap,mpicom)
       gsmap_gsize = mct_gsmap_gsize(gsmap)

       call mct_gsmap_OP(gsmap,iam,gsmOP)

       fullsize = 1
       do n = 1,ndims
          fullsize = fullsize*dims(n)
       enddo

       vsize = fullsize / gsize
       if (mod(fullsize,gsize) /= 0) then
          write(iulog,*) subname,' ERROR in vsize ',fullsize,gsize,vsize
          call endrun()
       endif

       allocate(compDOF(gsmap_lsize*vsize))

       if (present(switchdim)) then
          allocate(compDOF(gsmap_lsize*vsize))
          if (switchdim) then
             cnt = 0
             do m = 1,gsmap_lsize
                do n = 1,vsize
                   cnt = cnt + 1
                   compDOF(cnt) = (gsmOP(m)-1)*vsize + n
                enddo
             enddo
          else
             write(iulog,*) subname,' ERROR switch dims present must have switchdim true' 
             call endrun()
          end if
       else         ! currently allow for up to two vertical dimensions
          if (vsize /= 1 .and. vsize /= dims(ndims)) then
             vsize1 = vsize/dims(ndims)
             vsize2 = dims(ndims)
             if (vsize1*vsize2 /= vsize) then
                write(iulog,*)'vsize1= ',vsize1,' vsize2= ',vsize2,' vsize= ',vsize
                call endrun('error in vsize1 and vsize2 computation')
             end if
             cnt = 0
             do k = 1,vsize2
                do n = 1,vsize1
                   do m = 1,gsmap_lsize
                      cnt = cnt + 1
                      compDOF(cnt) = (k-1)*vsize1*gsmap_gsize + (n-1)*gsmap_gsize +  gsmOP(m) 
                   enddo
                enddo
             end do
          else
             cnt = 0
             do n = 1,vsize
                do m = 1,gsmap_lsize
                   cnt = cnt + 1
                   compDOF(cnt) = (n-1)*gsmap_gsize + gsmOP(m)
                enddo
             enddo
          end if
       end if

       if (debug > 1) then
          do m = 0,npes-1
             if (iam == m) then
                write(iulog,*) trim(subname),' sizes1  = ',iam,gsize,gsmap_gsize,gsmap_lsize
                write(iulog,*) trim(subname),' sizes2  = ',iam,fullsize,npes,vsize
                write(iulog,*) trim(subname),' compDOF = ',iam,size(compDOF),minval(compDOF),maxval(compDOF)
                call shr_sys_flush(iulog)
             endif
             call mpi_barrier(mpicom,status)
          enddo
       endif

       deallocate(gsmOP)

       call pio_initdecomp(pio_subsystem, baseTYPE, dims(1:ndims), compDOF, iodesc_list(iodnum)%iodesc)

       deallocate(compDOF)

       iodesc_list(iodnum)%type  = xtype
       iodesc_list(iodnum)%ndims = ndims
       iodesc_list(iodnum)%dims  = 0
       iodesc_list(iodnum)%dims(1:ndims) = dims(1:ndims)

    else

       if (iodnum > num_iodesc) then
          write(iulog,*) trim(subname),' ERROR: iodnum out of range ',iodnum,num_iodesc
          call endrun()
       endif

    endif

  end subroutine ncd_getiodesc

end module ncdio_pio
