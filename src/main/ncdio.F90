#include <misc.h>
#include <preproc.h>
#define BUILDPIO
#undef  BUILDPIO
#undef  SWITCH_DIMS
#define SWITCH_DIMS

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
  use clm_varctl     , only : iulog
  use clm_mct_mod
  use spmdGathScatMod
  use decompMod      , only : get_clmlevel_gsize,get_clmlevel_dsize
#if (defined BUILDPIO)
  use piolib_mod     ! _EXTERNAL
  use pio_types      ! _EXTERNAL
  use pio_kinds      , only : pio_offset
#endif
!
! !PUBLIC TYPES:
  implicit none
  include 'netcdf.inc'
  save

  public :: check_ret   ! checks return status of netcdf calls
#if (defined BUILDPIO)
  public :: check_ret_pio   ! checks return status of pio calls
#endif
  public :: check_var   ! determine if variable is on netcdf file
  public :: check_dim   ! validity check on dimension
  public :: ncd_open    ! open file
  public :: ncd_close   ! close file
  public :: ncd_redef   ! enter define mode
  public :: ncd_enddef  ! end define mode
  public :: ncd_setfill ! set file value
  public :: ncd_putatt  ! put attribute
  public :: ncd_defdim  ! define dimension
  public :: ncd_inqdid  ! inquire dimension id
  public :: ncd_inqdname ! inquire dimension name
  public :: ncd_inqdlen ! inquire dimension length
  public :: ncd_defvar  ! define variables
  public :: ncd_inqvid  ! inquire variable id
  public :: ncd_inqvname ! inquire variable name
  public :: ncd_inqvdims ! inquire variable ndims
  public :: ncd_inqvdids ! inquire variable dimids
  public :: ncd_iolocal ! write local data
  public :: ncd_ioglobal! write global data

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
  integer,parameter,public :: ncd_64bit_offset = nf_64bit_offset
  integer,parameter,public :: ncd_fill      = nf_fill
  integer,parameter,public :: ncd_nofill    = nf_nofill
  integer,parameter,public :: ncd_unlimited = nf_unlimited

!
! !REVISION HISTORY:
!
!EOP
!
! !PRIVATE METHODS:
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
     module procedure ncd_ioglobal_char_1d
     module procedure ncd_ioglobal_int_2d
     module procedure ncd_ioglobal_real_2d
     module procedure ncd_ioglobal_int_3d
     module procedure ncd_ioglobal_real_3d
  end interface
  private :: ncd_inqvdesc  ! inquire variable descriptor
  private :: ncd_inqiodesc ! inquire variable descriptor
#if (defined BUILDPIO)
  private :: ncd_setDOF    ! set DOF arrays for pio
#endif
  private :: scam_field_offsets ! get offset to proper lat/lon gridcell for SCAM

  logical,parameter,private :: lbcast_def = .false.  ! lbcast default
  integer,parameter,private :: debug = 0             ! local debug level
  logical,private,parameter :: lowmem2d = .true.
  logical,parameter,private :: pio_def = .false.

#if (defined BUILDPIO)
  logical,parameter,private :: pio_UseRearranger = .true.
  logical,parameter,private :: pio_SerialCDF = .false.
  integer,private           :: pio_num_iotasks
  integer,private           :: pio_num_aggregator
  integer,private           :: pio_DebugLevel = 2

  type(File_desc_t)         :: pio_File

  type pio_iodesc_plus_type
     character(len=32) :: name
     logical           :: set
     integer           :: ndims
     integer           :: dimids(4)
     type(IO_desc_t),pointer   :: pio_ioDesc
  end type pio_iodesc_plus_type
  integer,parameter     ,private :: pio_max_iodesc = 50
  integer               ,private :: pio_num_iodesc = 0
  type(pio_iodesc_plus_type) ,private, target :: pio_iodesc_list(pio_max_iodesc)

  type pio_vardesc_plus_type
     character(len=64) :: name
     integer           :: iodnum   ! iodesc associated with vardesc
     type(Var_desc_t)  :: pio_varDesc
  end type pio_vardesc_plus_type
  integer,parameter     ,private :: pio_max_vardesc = 500
  integer               ,private :: pio_num_vardesc = 0
  type(pio_vardesc_plus_type),private, target :: pio_varDesc_list(pio_max_vardesc)

#endif
  
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_dim
!
! !INTERFACE:
  subroutine check_dim(ncid, dimname, value, usepio)
!
! !DESCRIPTION:
! Validity check on dimension
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: dimname
    integer, intent(in) :: value
    logical,optional,intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: dimid, dimlen    ! temporaries
    logical :: lusepio           ! local usepio variable
    character(len=*),parameter :: subname='check_dim' ! subroutine name
!-----------------------------------------------------------------------

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
    if (masterproc) write(iulog,*) trim(subname),' WARNING: pio not implemented'
 else
    if (.not. masterproc) return

    call check_ret(nf_inq_dimid (ncid, trim(dimname), dimid), 'check_dim')
    call check_ret(nf_inq_dimlen (ncid, dimid, dimlen), 'check_dim')
    if (dimlen /= value) then
       write(iulog,*) trim(subname),' ERROR: mismatch of input dimension ',dimlen, &
            ' with expected value ',value,' for variable ',trim(dimname)
       call endrun()
    end if
 endif

  end subroutine check_dim

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_var
!
! !INTERFACE:
  subroutine check_var(ncid, varname, varid, readvar, usepio)
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
    logical,optional,intent(in)  :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ret     ! return value
    logical :: lusepio           ! local usepio variable
#if (defined BUILDPIO)
    type(Var_desc_t) :: pio_varDesc
#endif
    character(len=*),parameter :: subname='check_var' ! subroutine name
!-----------------------------------------------------------------------

 varid = -1
 readvar = .true.

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
#if (defined BUILDPIO)
    ret = PIO_inq_varid (pio_File, varname, pio_varDesc)
    if (ret/=PIO_noerr) then
       if (masterproc) write(iulog,*) trim(subname),': variable ',trim(varname),' is not on dataset'
       readvar = .false.
    end if
#endif
 else
    if (.not. masterproc) return
    ret = nf_inq_varid (ncid, varname, varid)
    if (ret/=NF_NOERR) then
       write(iulog,*) trim(subname),': variable ',trim(varname),' is not on dataset'
       readvar = .false.
    end if
 endif

  end subroutine check_var

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

#if (defined BUILDPIO)
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_ret_pio
!
! !INTERFACE:
  subroutine check_ret_pio(ret, cstring)
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
    character(len=*),parameter :: subname='check_ret_pio' ! subroutine name
!-----------------------------------------------------------------------

    if (ret /= PIO_noerr) then
       write(iulog,*)'pio error from ',trim(subname),':',trim(cstring),': ERROR = ',ret
       call endrun()
    end if

  end subroutine check_ret_pio
#endif
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_create
!
! !INTERFACE:
  subroutine ncd_create(filename,mode,ncid,cstring,usepio)
!
! !DESCRIPTION:
! create netcdf file
!
! !ARGUMENTS:
    implicit none
    character(len=*),intent(in) :: filename  ! file to open
    integer         ,intent(in) :: mode      ! nf_create mode
    integer         ,intent(out):: ncid      ! netcdf file id
    character(len=*),intent(in) :: cstring   ! comment string
    logical,optional,intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    integer :: stride        ! pe stride for num_iotasks
    integer :: iotype        ! type of output file
    character(len=*),parameter :: subname='ncd_create' ! subroutine name

!-----------------------------------------------------------------------

 ncid = -1

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
#if (defined BUILDPIO)
    call PIO_setDebugLevel(pio_DebugLevel)
    if (pio_UseRearranger) then
       if (pio_SerialCDF) then
          iotype = iotype_netcdf_rearrange      ! serial netcdf no rearrange
       else
          iotype = iotype_pnetcdf_rearrange     ! parallel netcdf no rearrange
       endif
    else
       if (pio_SerialCDF) then
          iotype = iotype_netcdf                ! serial netcdf w/ rearrange
       else
          iotype = iotype_pnetcdf               ! parallel netcdf w/ rearrange
       endif
    endif
    pio_num_iotasks = npes
    pio_num_aggregator = npes
    stride = max(npes/pio_num_iotasks,1)        ! regular pe stride
!tcx    if (debug > 1) then
       write(iulog,*) trim(subname),' iam = ',iam,pio_num_iotasks,pio_num_aggregator,stride,iotype
!    endif
    call PIO_initFile(iam, mpicom, &
             pio_num_iotasks,pio_num_aggregator, stride,iotype,pio_File)
    call check_ret_pio(PIO_CreateFile(pio_File,trim(filename)), &
       trim(subname)//':'//trim(cstring))
#endif
 else

    if (.not. masterproc) return

    call check_ret(nf_create(trim(filename),mode,ncid), &
       trim(subname)//':'//trim(cstring))

 endif

  end subroutine ncd_create

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_open
!
! !INTERFACE:
  subroutine ncd_open(filename,mode,ncid,cstring,usepio)
!
! !DESCRIPTION:
! open netcdf file
!
! !ARGUMENTS:
    implicit none
    character(len=*),intent(in) :: filename  ! file to open
    integer         ,intent(in) :: mode      ! nf_create mode
    integer         ,intent(out):: ncid      ! netcdf file id
    character(len=*),intent(in) :: cstring   ! comment string
    logical,optional,intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    character(len=*),parameter :: subname='ncd_open' ! subroutine name

!-----------------------------------------------------------------------

 ncid = -1

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then

#if (defined BUILDPIO)
    call PIO_setDebugLevel(pio_DebugLevel)
    call check_ret_pio(PIO_OpenFile(pio_File,trim(filename)), &
       trim(subname)//':'//trim(cstring))
#endif

 else

    if (.not. masterproc) return

    call check_ret(nf_open(trim(filename),mode,ncid), &
       trim(subname)//':'//trim(cstring))

 endif

  end subroutine ncd_open

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_close
!
! !INTERFACE:
  subroutine ncd_close(ncid,cstring,usepio)
!
! !DESCRIPTION:
! close netcdf file
!
! !ARGUMENTS:
    implicit none
    integer         ,intent(in) :: ncid      ! netcdf file id
    character(len=*),intent(in) :: cstring   ! comment string
    logical,optional,intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    character(len=*),parameter :: subname='ncd_close' ! subroutine name

!-----------------------------------------------------------------------

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
#if (defined BUILDPIO)
    call PIO_CloseFile(pio_File)
#endif
 else

    if (.not. masterproc) return

    call check_ret(nf_close(ncid), &
       trim(subname)//':'//trim(cstring))

 endif

  end subroutine ncd_close

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_redef
!
! !INTERFACE:
  subroutine ncd_redef(ncid,cstring,usepio)
!
! !DESCRIPTION:
! redef netcdf file
!
! !ARGUMENTS:
    implicit none
    integer         ,intent(in) :: ncid      ! netcdf file id
    character(len=*),intent(in) :: cstring   ! comment string
    logical,optional,intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    character(len=*),parameter :: subname='ncd_redef' ! subroutine name

!-----------------------------------------------------------------------

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
#if (defined BUILDPIO)
    if (masterproc) write(iulog,*) trim(subname),':',trim(cstring),' WARNING: pio not implemented'
#endif
 else

    if (.not. masterproc) return

    call check_ret(nf_redef(ncid), &
       trim(subname)//':'//trim(cstring))

 endif

  end subroutine ncd_redef

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_enddef
!
! !INTERFACE:
  subroutine ncd_enddef(ncid,cstring,usepio)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    integer         ,intent(in) :: ncid      ! netcdf file id
    character(len=*),intent(in) :: cstring   ! comment string
    logical,optional,intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    character(len=*),parameter :: subname='ncd_enddef' ! subroutine name

!-----------------------------------------------------------------------

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
#if (defined BUILDPIO)
    call check_ret_pio(PIO_enddef(pio_File),&
       trim(subname)//':'//trim(cstring))
#endif
 else

    if (.not. masterproc) return

    call check_ret(nf_enddef(ncid), &
       trim(subname)//':'//trim(cstring))

 endif

  end subroutine ncd_enddef

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_setfill
!
! !INTERFACE:
  subroutine ncd_setfill(ncid,mode,old_mode,cstring,usepio)
!
! !DESCRIPTION:
! setfill netcdf file
!
! !ARGUMENTS:
    implicit none
    integer         ,intent(in) :: ncid      ! netcdf file id
    integer         ,intent(in) :: mode      ! fill mode
    integer         ,intent(out):: old_mode  ! old mode
    character(len=*),intent(in) :: cstring   ! comment string
    logical,optional,intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    character(len=*),parameter :: subname='ncd_setfill' ! subroutine name

!-----------------------------------------------------------------------

 old_mode = -1

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
#if (defined BUILDPIO)
    if (masterproc) write(iulog,*) trim(subname),':',trim(cstring),' WARNING: pio not implemented'
#endif
 else
    if (.not. masterproc) return
    call check_ret(nf_set_fill(ncid,mode,old_mode), &
       trim(subname)//':'//trim(cstring))
 endif

  end subroutine ncd_setfill

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqdid
!
! !INTERFACE:
  subroutine ncd_inqdid(ncid,name,dimid,cstring,usepio)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    integer         ,intent(in) :: ncid      ! netcdf file id
    character(len=*),intent(in) :: name      ! dimension name
    integer         ,intent(out):: dimid     ! dimension id
    character(len=*),intent(in) :: cstring   ! comment string
    logical,optional,intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    character(len=*),parameter :: subname='ncd_inqdid' ! subroutine name

!-----------------------------------------------------------------------

 dimid = -1

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
#if (defined BUILDPIO)
    call check_ret_pio(PIO_inq_dimid(pio_File,name,dimid),&
       trim(subname)//':'//trim(cstring))
#endif
 else
    if (.not. masterproc) return
    call check_ret(nf_inq_dimid(ncid,name,dimid), &
       trim(subname)//':'//trim(cstring))
 endif

  end subroutine ncd_inqdid

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqdlen
!
! !INTERFACE:
  subroutine ncd_inqdlen(ncid,dimid,len,cstring,usepio)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    integer         ,intent(in) :: ncid      ! netcdf file id
    integer         ,intent(in) :: dimid     ! dimension id
    integer         ,intent(out):: len       ! dimension len
    character(len=*),intent(in) :: cstring   ! comment string
    logical,optional,intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    character(len=*),parameter :: subname='ncd_inqdlen' ! subroutine name

!-----------------------------------------------------------------------

 len = -1

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
#if (defined BUILDPIO)
    call check_ret_pio(PIO_inq_dimlen(pio_File,dimid,len),&
       trim(subname)//':'//trim(cstring))
#endif
 else
    if (.not. masterproc) return
    call check_ret(nf_inq_dimlen(ncid,dimid,len), &
       trim(subname)//':'//trim(cstring))
 endif

  end subroutine ncd_inqdlen

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqdname
!
! !INTERFACE:
  subroutine ncd_inqdname(ncid,dimid,dname,cstring,usepio)
!
! !DESCRIPTION:
! inquire dim name
!
! !ARGUMENTS:
    implicit none
    integer         ,intent(in) :: ncid      ! netcdf file id
    integer         ,intent(in) :: dimid     ! dimension id
    character(len=*),intent(out):: dname     ! dimension name
    character(len=*),intent(in) :: cstring   ! comment string
    logical,optional,intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    character(len=*),parameter :: subname='ncd_inqdname' ! subroutine name

!-----------------------------------------------------------------------

 dname = ''

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
#if (defined BUILDPIO)
    call check_ret_pio(PIO_inq_dimname(pio_File,dimid,dname),&
       trim(subname)//':'//trim(cstring))
#endif
 else
    if (.not. masterproc) return
    call check_ret(nf_inq_dimname(ncid,dimid,dname), &
       trim(subname)//':'//trim(cstring))
 endif

  end subroutine ncd_inqdname

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqvdesc
!
! !INTERFACE:
  subroutine ncd_inqvdesc(name,vdnum,cstring,usepio)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    character(len=*),intent(in) :: name      ! variable name
    integer         ,intent(out):: vdnum     ! vardesc num
    character(len=*),intent(in) :: cstring   ! comment string
    logical,optional,intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    logical :: found             ! search flag
    integer :: n
    character(len=*),parameter :: subname='ncd_inqvdesc' ! subroutine name

!-----------------------------------------------------------------------

    vdnum = -1

    lusepio = pio_def
    if (present(usepio)) then
       lusepio = usepio
    endif

    if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

    if (lusepio) then
#if (defined BUILDPIO)
       n = 1
       found = .false.
       do while (n <= pio_num_vardesc .and. .not.found)
          if (trim(name) == trim(pio_vardesc_list(n)%name)) then
             found = .true.
             vdnum = n
          endif
          n = n + 1
       enddo
#endif
    else
       ! do nothing
    endif

  end subroutine ncd_inqvdesc

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqiodesc
!
! !INTERFACE:
  subroutine ncd_inqiodesc(ndims,dimids,iodnum,cstring,usepio)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    integer         ,intent(in) :: ndims     ! number of dims
    integer         ,intent(in) :: dimids(:) ! dimids
    integer         ,intent(out):: iodnum    ! iodesc num
    character(len=*),intent(in) :: cstring   ! comment string
    logical,optional,intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    logical :: found             ! search flag
    integer :: n,m
    character(len=*),parameter :: subname='ncd_inqiodesc' ! subroutine name

!-----------------------------------------------------------------------

    iodnum = -1

    lusepio = pio_def
    if (present(usepio)) then
       lusepio = usepio
    endif

    if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

    if (lusepio) then
#if (defined BUILDPIO)
       n = 1
       found = .false.
       do while (n <= pio_num_iodesc .and. .not.found)
          if (ndims == pio_iodesc_list(n)%ndims) then
             found = .true.
             do m = 1,ndims
                if (dimids(m) /= pio_iodesc_list(n)%dimids(m)) then
                   found = .false.
                endif
             enddo
          endif
          if (found) then
             iodnum = n
          endif
          n = n + 1
       enddo
#endif
    else
       ! do nothing
    endif

  end subroutine ncd_inqiodesc

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqvid
!
! !INTERFACE:
  subroutine ncd_inqvid(ncid,name,varid,cstring,readvar,usepio,pio_vardesc)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    integer         ,intent(in) :: ncid      ! netcdf file id
    character(len=*),intent(in) :: name      ! variable name
    integer         ,intent(out):: varid     ! variable id
    character(len=*),intent(in) :: cstring   ! comment string
    logical,optional,intent(out):: readvar   ! does variable exist
    logical,optional,intent(in) :: usepio    ! use pio lib
#if (defined BUILDPIO)
    type(Var_desc_t),optional,intent(inout):: pio_varDesc
#else
    integer,optional,intent(out) :: pio_varDesc       ! dummy
#endif
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    integer :: ret               ! return code
    character(len=*),parameter :: subname='ncd_inqvid' ! subroutine name

!-----------------------------------------------------------------------

 varid = -1
 if (present(readvar)) then
    readvar = .false.
 endif

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
#if (defined BUILDPIO)
    if (present(readvar)) then
       ret = PIO_inq_varid(pio_File,name,pio_varDesc)
       if (ret/=PIO_noerr) then
          if (masterproc) write(iulog,*) trim(subname),': variable ',trim(name),' is not on dataset'
          readvar = .false.
       else
          readvar = .true.
       end if
    else
       call check_ret_pio(PIO_inq_varid(pio_File,name,pio_varDesc), &
          trim(subname)//':'//trim(cstring))
    endif
    varid = pio_varDesc%varid
#endif
 else
    if (.not. masterproc) return
    if (present(readvar)) then
       ret = nf_inq_varid (ncid, name, varid)
       if (ret/=NF_NOERR) then
          if (masterproc) write(iulog,*) trim(subname),': variable ',trim(name),' is not on dataset'
          readvar = .false.
       else
          readvar = .true.
       end if
    else
       call check_ret(nf_inq_varid(ncid,name,varid), &
          trim(subname)//':'//trim(cstring))
    endif
 endif

  end subroutine ncd_inqvid

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqvdims
!
! !INTERFACE:
  subroutine ncd_inqvdims(ncid,varid,ndims,cstring,usepio,pio_varDesc)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    integer         ,intent(in) :: ncid      ! netcdf file id
    integer         ,intent(in) :: varid     ! variable id
    integer         ,intent(out):: ndims     ! variable ndims
    character(len=*),intent(in) :: cstring   ! comment string
    logical,optional,intent(in) :: usepio    ! use pio lib
#if (defined BUILDPIO)
    type(Var_desc_t),optional,intent(inout):: pio_varDesc
#else
    integer,optional,intent(in) :: pio_varDesc       ! dummy
#endif
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    character(len=*),parameter :: subname='ncd_inqvdims' ! subroutine name

!-----------------------------------------------------------------------

 ndims = -1

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
#if (defined BUILDPIO)
    if (.not.present(pio_varDesc)) then
       write(iulog,*) trim(subname),' ERROR pio_varDesc must be an argument'
       call endrun()
    endif
    call check_ret_pio(PIO_inq_varndims(pio_File,pio_varDesc,ndims),&
       trim(subname)//':'//trim(cstring))
#endif
 else
    if (.not. masterproc) return
    call check_ret(nf_inq_varndims(ncid,varid,ndims), &
       trim(subname)//':'//trim(cstring))
 endif

  end subroutine ncd_inqvdims

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqvname
!
! !INTERFACE:
  subroutine ncd_inqvname(ncid,varid,vname,cstring,usepio,pio_varDesc)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    integer         ,intent(in) :: ncid      ! netcdf file id
    integer         ,intent(in) :: varid     ! variable id
    character(len=*),intent(out):: vname     ! variable vname
    character(len=*),intent(in) :: cstring   ! comment string
    logical,optional,intent(in) :: usepio    ! use pio lib
#if (defined BUILDPIO)
    type(Var_desc_t),optional,intent(inout):: pio_varDesc
#else
    integer,optional,intent(in) :: pio_varDesc       ! dummy
#endif
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    character(len=*),parameter :: subname='ncd_inqvname' ! subroutine name

!-----------------------------------------------------------------------

 vname = ''

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
#if (defined BUILDPIO)
    if (.not.present(pio_varDesc)) then
       write(iulog,*) trim(subname),' ERROR pio_varDesc must be an argument'
       call endrun()
    endif
    call check_ret_pio(PIO_inq_varname(pio_File,pio_varDesc,vname),&
       trim(subname)//':'//trim(cstring))
#endif
 else
    if (.not. masterproc) return
    call check_ret(nf_inq_varname(ncid,varid,vname), &
       trim(subname)//':'//trim(cstring))
 endif

  end subroutine ncd_inqvname

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqvdids
!
! !INTERFACE:
  subroutine ncd_inqvdids(ncid,varid,dids,cstring,usepio,pio_varDesc)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    integer         ,intent(in) :: ncid      ! netcdf file id
    integer         ,intent(in) :: varid     ! variable id
    integer         ,intent(out):: dids(:)   ! variable dids
    character(len=*),intent(in) :: cstring   ! comment string
    logical,optional,intent(in) :: usepio    ! use pio lib
#if (defined BUILDPIO)
    type(Var_desc_t),optional,intent(inout):: pio_varDesc
#else
    integer,optional,intent(in) :: pio_varDesc       ! dummy
#endif
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    character(len=*),parameter :: subname='ncd_inqvdids' ! subroutine name

!-----------------------------------------------------------------------

 dids = -1

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
#if (defined BUILDPIO)
    if (.not.present(pio_varDesc)) then
       write(iulog,*) trim(subname),' ERROR pio_varDesc must be an argument'
       call endrun()
    endif
    call check_ret_pio(PIO_inq_vardimid(pio_File,pio_varDesc,dids),&
       trim(subname)//':'//trim(cstring))
#endif
 else
    if (.not. masterproc) return
    call check_ret(nf_inq_vardimid(ncid,varid,dids), &
       trim(subname)//':'//trim(cstring))
 endif

  end subroutine ncd_inqvdids

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_putatt_int
!
! !INTERFACE:
  subroutine ncd_putatt_int(ncid,varid,attrib,value,cstring,xtype,usepio)
!
! !DESCRIPTION:
! Check return status from netcdf call
!
! !ARGUMENTS:
    implicit none
    integer         ,intent(in) :: ncid      ! netcdf file id
    integer         ,intent(in) :: varid     ! netcdf var id
    character(len=*),intent(in) :: attrib    ! netcdf attrib
    integer         ,intent(in) :: value     ! netcdf attrib value
    character(len=*),intent(in) :: cstring   ! comment string
    integer,optional,intent(in) :: xtype     ! netcdf data type
    logical,optional,intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    integer :: lxtype
    character(len=*),parameter :: subname='ncd_putatt_int' ! subroutine name

!-----------------------------------------------------------------------

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
#if (defined BUILDPIO)
    call check_ret_pio(PIO_put_att(pio_File,varid,trim(attrib),value), &
       trim(subname)//':'//trim(cstring))
#endif
 else

    if (.not. masterproc) return

    if (present(xtype)) then
       lxtype = xtype
    else
       lxtype = nf_int
    endif

    call check_ret(nf_put_att_int(ncid, varid, trim(attrib), lxtype, 1, value), &
       trim(subname)//':'//trim(cstring))

 endif

  end subroutine ncd_putatt_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_putatt_char
!
! !INTERFACE:
  subroutine ncd_putatt_char(ncid,varid,attrib,value,cstring,xtype,usepio)
!
! !DESCRIPTION:
! Check return status from netcdf call
!
! !ARGUMENTS:
    implicit none
    integer         ,intent(in) :: ncid      ! netcdf file id
    integer         ,intent(in) :: varid     ! netcdf var id
    character(len=*),intent(in) :: attrib    ! netcdf attrib
    character(len=*),intent(in) :: value     ! netcdf attrib value
    character(len=*),intent(in) :: cstring   ! comment string
    integer,optional,intent(in) :: xtype     ! netcdf data type
    logical,optional,intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    character(len=*),parameter :: subname='ncd_putatt_char' ! subroutine name

!-----------------------------------------------------------------------

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
#if (defined BUILDPIO)
    call check_ret_pio(PIO_put_att(pio_File,varid,trim(attrib),value), &
       trim(subname)//':'//trim(cstring))
#endif
 else

    if (.not. masterproc) return

    call check_ret(nf_put_att_text(ncid, varid, trim(attrib), len_trim(value), value), &
       trim(subname)//':'//trim(cstring))

 endif

  end subroutine ncd_putatt_char

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_putatt_real
!
! !INTERFACE:
  subroutine ncd_putatt_real(ncid,varid,attrib,value,cstring,xtype,usepio)
!
! !DESCRIPTION:
! Check return status from netcdf call
!
! !ARGUMENTS:
    implicit none
    integer         ,intent(in) :: ncid      ! netcdf file id
    integer         ,intent(in) :: varid     ! netcdf var id
    character(len=*),intent(in) :: attrib    ! netcdf attrib
    real(r8)        ,intent(in) :: value     ! netcdf attrib value
    character(len=*),intent(in) :: cstring   ! comment string
    integer,optional,intent(in) :: xtype     ! netcdf data type
    logical,optional,intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    integer :: lxtype
    real*4  :: value4
    character(len=*),parameter :: subname='ncd_putatt_real' ! subroutine name

!-----------------------------------------------------------------------

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 value4 = value

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
#if (defined BUILDPIO)
    if (lxtype == nf_double) then
       call check_ret_pio(PIO_put_att(pio_File,varid,trim(attrib),value), &
       trim(subname)//':'//trim(cstring))
    else
       call check_ret_pio(PIO_put_att(pio_File,varid,trim(attrib),value4), &
       trim(subname)//':'//trim(cstring))
    endif
#endif
 else

    if (.not. masterproc) return

    if (present(xtype)) then
       lxtype = xtype
    else
       lxtype = nf_double
    endif

    call check_ret(nf_put_att_double(ncid, varid, trim(attrib), lxtype, 1, value), &
       trim(subname)//':'//trim(cstring))

 endif

  end subroutine ncd_putatt_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_defdim
!
! !INTERFACE:
  subroutine ncd_defdim(ncid,attrib,value,dimid,cstring,usepio)
!
! !DESCRIPTION:
! Check return status from netcdf call
!
! !ARGUMENTS:
    implicit none
    integer         ,intent(in) :: ncid      ! netcdf file id
    character(len=*),intent(in) :: attrib    ! netcdf attrib
    integer         ,intent(in) :: value     ! netcdf attrib value
    integer         ,intent(out):: dimid     ! netcdf dimension id
    character(len=*),intent(in) :: cstring   ! comment string
    logical,optional,intent(in) :: usepio    ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: lusepio           ! local usepio variable
    character(len=*),parameter :: subname='ncd_defdim' ! subroutine name

!-----------------------------------------------------------------------

 dimid = -1

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

 if (lusepio) then
#if (defined BUILDPIO)
    call check_ret_pio(PIO_def_dim(pio_File,attrib,value,dimid), &
       trim(subname)//':'//trim(cstring))
#endif
 else

    if (.not. masterproc) return

    call check_ret(nf_def_dim(ncid, trim(attrib), value, dimid), &
       trim(subname)//':'//trim(cstring))

 endif

  end subroutine ncd_defdim

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_defvar_bynf
!
! !INTERFACE:
  subroutine ncd_defvar_bynf(ncid, varname, xtype, ndims, dimid, varid, &
       cstring, long_name, units, cell_method, missing_value, fill_value, &
       imissing_value, ifill_value, usepio)
!
! !DESCRIPTION:
!  Define a netcdf variable
!
! !ARGUMENTS:
    implicit none
    integer         , intent(in)  :: ncid                    ! input unit
    character(len=*), intent(in)  :: varname                 ! variable name
    integer         , intent(in)  :: xtype                   ! external type
    integer         , intent(in)  :: ndims                   ! number of dims
    integer         , intent(in), optional :: dimid(:)       ! dimids
    integer         , intent(out) :: varid                   ! returned var id
    character(len=*), intent(in), optional :: cstring        ! caller string
    character(len=*), intent(in), optional :: long_name      ! attribute
    character(len=*), intent(in), optional :: units          ! attribute
    character(len=*), intent(in), optional :: cell_method    ! attribute
    real(r8)        , intent(in), optional :: missing_value  ! attribute for real
    real(r8)        , intent(in), optional :: fill_value     ! attribute for real
    integer         , intent(in), optional :: imissing_value ! attribute for int
    integer         , intent(in), optional :: ifill_value    ! attribute for int
    logical         , intent(in), optional :: usepio         ! use pio lib
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: n              ! indices
    character(len=256) :: str ! temporary
    character(len=256) :: lsubname ! temporary
    integer :: ldimid(4)      ! local dimid
    integer :: dimid0(1)      ! local dimid
    integer :: vdnum          ! vardesc num
    integer :: iodnum         ! iodesc num
    logical :: lusepio        ! local usepio variable
    character(len=*),parameter :: subname='ncd_defvar_bynf' ! subroutine name
!-----------------------------------------------------------------------

 varid = -1

 lusepio = pio_def
 if (present(usepio)) then
    lusepio = usepio
 endif

 if (present(cstring)) then
    lsubname = trim(subname)//':'//trim(cstring)
 else
    lsubname = trim(subname)
 endif

 dimid0 = 0
 ldimid = 0
 if (present(dimid)) then
    ldimid(1:ndims) = dimid(1:ndims)
 else   ! ndims must be zero if dimid not present
    if (ndims /= 0) then
       write(iulog,*) trim(lsubname),' ERROR: dimid not supplied and ndims ne 0 ',trim(varname),ndims
       call endrun()
    endif
 endif

 if (masterproc .and. debug > 1) then
    write(iulog,*) trim(subname),' ',trim(varname),lusepio
    write(iulog,*) trim(subname),' ',trim(varname),xtype,ndims,ldimid(1:ndims)
 endif

 if (lusepio) then
#if (defined BUILDPIO)
    call ncd_inqvdesc(trim(varname),vdnum,lsubname,usepio=lusepio)
    if (vdnum < 1) then
       pio_num_vardesc = pio_num_vardesc + 1
       if (pio_num_vardesc > pio_max_vardesc) then
          write(iulog,*) trim(subname),' ERROR num_vardesc gt max_vardesc ',trim(varname),pio_max_vardesc
          call endrun()
       endif
       vdnum = pio_num_vardesc
       pio_varDesc_list(vdnum)%name = trim(varname)

       call ncd_inqiodesc(ndims,ldimid,iodnum,lsubname,usepio=lusepio)

       if (iodnum < 1) then
          ! generate iodesc
          pio_num_iodesc = pio_num_iodesc + 1
          if (pio_num_iodesc > pio_max_iodesc) then
             write(iulog,*) trim(subname),' ERROR num_iodesc gt max_iodesc ',trim(varname),pio_max_iodesc
             call endrun()
          endif
          iodnum = pio_num_iodesc
          pio_iodesc_list(iodnum)%ndims = ndims
          pio_iodesc_list(iodnum)%dimids = 0
          pio_iodesc_list(iodnum)%dimids(1:ndims) = ldimid(1:ndims)
          pio_iodesc_list(iodnum)%set = .false.
!tcx          if (masterproc .and. debug > 1) then
          if (masterproc) then
             write(iulog,*) trim(subname),' creating iodesc at ',trim(varname),vdnum,iodnum,ndims,ldimid(1:ndims)
         endif
       else
          if (iodnum > pio_num_iodesc) then
             write(iulog,*) trim(lsubname),' ERROR: iodnum out of range ',iodnum,pio_num_iodesc
             call endrun()
          endif
       endif
       if (masterproc .and. debug > 1) then
          write(iulog,*) trim(subname),' creating vardesc for ',trim(varname),vdnum,iodnum
       endif
       pio_varDesc_list(vdnum)%iodnum = iodnum
       call pio_setVarDesc(pio_iodesc_list(iodnum)%pio_iodesc,pio_vardesc_list(vdnum)%pio_vardesc)

    else  ! vardesc already exists
       if (vdnum > pio_num_vardesc) then
          write(iulog,*) trim(lsubname),' ERROR: vdnum out of range ',vdnum,pio_num_vardesc
          call endrun()
       endif
       if (masterproc .and. debug > 1) then
          write(iulog,*) trim(subname),' vardesc exists for ',trim(varname),vdnum
       endif
    endif

    if (present(dimid)) then
       call check_ret_pio(PIO_def_var(pio_File,trim(varname),xtype, &
          dimid(1:ndims),  pio_varDesc_list(vdnum)%pio_varDesc), &
          trim(lsubname))
    else
       call check_ret_pio(PIO_def_var(pio_File,trim(varname),xtype, &
          dimid0        ,  pio_varDesc_list(vdnum)%pio_varDesc), &
          trim(lsubname))
    endif
    varid = pio_varDesc_list(vdnum)%pio_varDesc%varid
#endif
 else

    if (.not. masterproc) return

    ! Define variable
    call check_ret(nf_def_var(ncid, trim(varname), xtype, ndims, ldimid, varid), lsubname)

 endif

    if (present(long_name)) then
       call ncd_putatt(ncid, varid, 'long_name', trim(long_name), lsubname, usepio=lusepio)
    end if
    if (present(units)) then
       call ncd_putatt(ncid, varid, 'units', trim(units), lsubname, usepio=lusepio)
    end if
    if (present(cell_method)) then
       str = 'time: ' // trim(cell_method)
       call ncd_putatt(ncid, varid, 'cell_method', trim(str), lsubname, usepio=lusepio)
    end if
    if (present(fill_value)) then
       call ncd_putatt(ncid, varid, '_FillValue', fill_value, lsubname, xtype, usepio=lusepio)
    end if
    if (present(missing_value)) then
       call ncd_putatt(ncid, varid, 'missing_value', missing_value, lsubname, xtype, usepio=lusepio)
    end if
    if (present(ifill_value)) then
       call ncd_putatt(ncid, varid, '_FillValue', ifill_value, lsubname, xtype, usepio=lusepio)
    end if
    if (present(imissing_value)) then
       call ncd_putatt(ncid, varid, 'missing_value', imissing_value, lsubname, xtype, usepio=lusepio)
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
       imissing_value, ifill_value, usepio)
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
    logical         , intent(in), optional :: usepio         ! use pio lib
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
    logical :: lusepio           ! local usepio variable
    character(len=256) :: str ! temporary
    character(len=*),parameter :: subname='ncd_defvar_bygrid' ! subroutine name
!-----------------------------------------------------------------------

    lusepio = pio_def
    if (present(usepio)) then
       lusepio = usepio
    endif

    dimid(:) = 0

    if (masterproc .and. debug > 1) write(iulog,*) trim(subname),lusepio

    if (lusepio) then
       ! continue
    else
       if (.not. masterproc) return
    endif

    ! Determine dimension ids for variable

    if (present(dim1name)) then
       call ncd_inqdid(ncid, dim1name, dimid(1), subname, lusepio)
    end if
    if (present(dim2name)) then
       call ncd_inqdid(ncid, dim2name, dimid(2), subname, lusepio)
    end if
    if (present(dim3name)) then
       call ncd_inqdid(ncid, dim3name, dimid(3), subname, lusepio)
    end if
    if (present(dim4name)) then
       call ncd_inqdid(ncid, dim4name, dimid(4), subname, lusepio)
    end if
    if (present(dim5name)) then
       call ncd_inqdid(ncid, dim5name, dimid(5), subname, lusepio)
    end if


#if (defined SWITCH_DIMS)
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
#endif

 ! Define variable

    ndims = 0
    if (present(dim1name)) then
       do n = 1, size(dimid)
          if (dimid(n) /= 0) ndims = ndims + 1
       end do
    endif

    call ncd_defvar_bynf(ncid,varname,xtype,ndims,dimid,varid,subname,usepio=lusepio)

    if (present(long_name)) then
       call ncd_putatt(ncid, varid, 'long_name', trim(long_name), subname, usepio=lusepio)
    end if
    if (present(units)) then
       call ncd_putatt(ncid, varid, 'units', trim(units), subname, usepio=lusepio)
    end if
    if (present(cell_method)) then
       str = 'time: ' // trim(cell_method)
       call ncd_putatt(ncid, varid, 'cell_method', trim(str), subname, usepio=lusepio)
    end if
    if (present(fill_value)) then
       call ncd_putatt(ncid, varid, '_FillValue', fill_value, subname, xtype, usepio=lusepio)
    end if
    if (present(missing_value)) then
       call ncd_putatt(ncid, varid, 'missing_value', missing_value, subname, xtype, usepio=lusepio)
    end if
    if (present(ifill_value)) then
       call ncd_putatt(ncid, varid, '_FillValue', ifill_value, subname, xtype, usepio=lusepio)
    end if
    if (present(imissing_value)) then
       call ncd_putatt(ncid, varid, 'missing_value', imissing_value, subname, xtype, usepio=lusepio)
    end if

  end subroutine ncd_defvar_bygrid

!-----------------------------------------------------------------------
!***** begin include ncdio_local_subs.inc *****
!      this next section has been autogenerated from 
!      --- gen_ncdio_local_subs.csh ---
!      and manually included.  this section ends at "end include ..."
!-----------------------------------------------------------------------
!  #include <ncdio_local_subs.inc>
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_int_1d
!
! !INTERFACE:
  subroutine ncd_iolocal_int_1d(varname, data, dim1name, &
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
    integer         , pointer     :: data(:)            ! local decomposition data
    character(len=*), intent(in)  :: dim1name           ! dimension name
    integer         , optional, intent(in) :: nlonxy    ! 2d longitude size
    integer         , optional, intent(in) :: nlatxy    ! 2d latitude size
    integer         , optional, intent(in) :: nt        ! time sample index
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: missing   ! missing value
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
    integer :: lmissing                 ! local missing value
    logical :: lusepio           ! local usepio variable
    character(len=8) :: clmlevel        ! clmlevel
    character(len=*),parameter :: subname='ncd_iolocal_int_1d' ! subroutine name
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
       lmissing = ispval
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

    call ncd_iolocal_gs_int1d(ncid, varname, flag, data, clmlevel, start, count, ier, lmissing, usepio=lusepio)

    if (present(readvar)) then
       readvar = .false.
       if (ier == 0) readvar = .true.
    endif

  end subroutine ncd_iolocal_int_1d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_int_2d
!
! !INTERFACE:
  subroutine ncd_iolocal_int_2d(varname, data, dim1name, dim2name, &
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
    integer         , pointer     :: data(:,:)          ! local decomposition input data
    character(len=*), intent(in)  :: dim1name           ! dimension 1 name
    character(len=*), intent(in)  :: dim2name           ! dimension 2 name
    integer         , optional, intent(in) :: nlonxy    ! 2d longitude size
    integer         , optional, intent(in) :: nlatxy    ! 2d latitude size
    integer         , optional, intent(in) :: nt        ! time sample index
    integer         , optional, intent(in) :: lowerb2,upperb2 ! lower and upper bounds of second dimension
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: missing   ! missing value
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
    integer :: lmissing                 ! local missing value
    logical :: lusepio           ! local usepio variable
    integer ,pointer :: data1d(:)       ! 1 level data
    character(len=8) :: clmlevel        ! clmlevel
    character(len=*),parameter :: subname='ncd_iolocal_int_2d' ! subroutine name
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
       lmissing = ispval
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
          call ncd_iolocal_gs_int1d(ncid, varname, flag, data1d, clmlevel, start, count, ier, lmissing, usepio=lusepio)
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
       call ncd_iolocal_gs_int2d(ncid, varname, flag, data, clmlevel, start, count, ier, lmissing, usepio=lusepio)
    endif

    if (present(readvar)) then
       readvar = .false.
       if (ier == 0) readvar = .true.
    endif

  end subroutine ncd_iolocal_int_2d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_gs_int1d
!
! !INTERFACE:
  subroutine ncd_iolocal_gs_int1d(ncid, varname, flag, data, clmlevel, start, count, status, missing, usepio)
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
    integer ,pointer              :: data(:)    ! local decomposition input data (out)
    character(len=*) ,intent(in)  :: clmlevel   ! type of grid
    integer, optional,intent(in)  :: start(:)   ! netcdf start index
    integer, optional,intent(in)  :: count(:)   ! netcdf count index
    integer, optional,intent(out) :: status     ! return code
    integer ,optional,intent(in)  :: missing    ! missing value
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
    integer , pointer :: arrayg(:)
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
    character(len=*),parameter :: subname='ncd_iolocal_gs_int1d' ! subroutine name
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
            baseTYPE = MPI_INTEGER
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
               call check_ret(nf_put_vara_int(ncid, varid, start, count, arrayg), subname)
            else
               call check_ret(nf_put_var_int(ncid, varid, arrayg), subname)
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

  end subroutine ncd_iolocal_gs_int1d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_gs_int2d
!
! !INTERFACE:
  subroutine ncd_iolocal_gs_int2d(ncid, varname, flag, data, clmlevel, start, count, status, missing, usepio)
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
    integer ,pointer              :: data(:,:)  ! local decomposition input data (out)
    character(len=*) ,intent(in)  :: clmlevel   ! type of grid
    integer, optional,intent(in)  :: start(:)   ! netcdf start index
    integer, optional,intent(in)  :: count(:)   ! netcdf count index
    integer, optional,intent(out) :: status     ! return code
    integer ,optional,intent(in)  :: missing    ! missing value
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
    integer , pointer :: arrayg(:,:)
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
    character(len=*),parameter :: subname='ncd_iolocal_gs_int2d' ! subroutine name
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
      if (present(missing)) then
         call gather_data_to_master(data,arrayg,clmlevel,missing)
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

  end subroutine ncd_iolocal_gs_int2d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_real_1d
!
! !INTERFACE:
  subroutine ncd_iolocal_real_1d(varname, data, dim1name, &
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
    real(r8)        , pointer     :: data(:)            ! local decomposition data
    character(len=*), intent(in)  :: dim1name           ! dimension name
    integer         , optional, intent(in) :: nlonxy    ! 2d longitude size
    integer         , optional, intent(in) :: nlatxy    ! 2d latitude size
    integer         , optional, intent(in) :: nt        ! time sample index
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    real(r8)        , optional, intent(in) :: missing   ! missing value
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
    real(r8):: lmissing                 ! local missing value
    logical :: lusepio           ! local usepio variable
    character(len=8) :: clmlevel        ! clmlevel
    character(len=*),parameter :: subname='ncd_iolocal_real_1d' ! subroutine name
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
       lmissing = spval
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

    call ncd_iolocal_gs_real1d(ncid, varname, flag, data, clmlevel, start, count, ier, lmissing, usepio=lusepio)

    if (present(readvar)) then
       readvar = .false.
       if (ier == 0) readvar = .true.
    endif

  end subroutine ncd_iolocal_real_1d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_real_2d
!
! !INTERFACE:
  subroutine ncd_iolocal_real_2d(varname, data, dim1name, dim2name, &
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
    real(r8)        , pointer     :: data(:,:)          ! local decomposition input data
    character(len=*), intent(in)  :: dim1name           ! dimension 1 name
    character(len=*), intent(in)  :: dim2name           ! dimension 2 name
    integer         , optional, intent(in) :: nlonxy    ! 2d longitude size
    integer         , optional, intent(in) :: nlatxy    ! 2d latitude size
    integer         , optional, intent(in) :: nt        ! time sample index
    integer         , optional, intent(in) :: lowerb2,upperb2 ! lower and upper bounds of second dimension
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    real(r8)        , optional, intent(in) :: missing   ! missing value
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
    real(r8):: lmissing                 ! local missing value
    logical :: lusepio           ! local usepio variable
    real(r8),pointer :: data1d(:)       ! 1 level data
    character(len=8) :: clmlevel        ! clmlevel
    character(len=*),parameter :: subname='ncd_iolocal_real_2d' ! subroutine name
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
       lmissing = spval
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
          call ncd_iolocal_gs_real1d(ncid, varname, flag, data1d, clmlevel, start, count, ier, lmissing, usepio=lusepio)
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
       call ncd_iolocal_gs_real2d(ncid, varname, flag, data, clmlevel, start, count, ier, lmissing, usepio=lusepio)
    endif

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
  subroutine ncd_iolocal_gs_real1d(ncid, varname, flag, data, clmlevel, start, count, status, missing, usepio)
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
    integer, optional,intent(out) :: status     ! return code
    real(r8),optional,intent(in)  :: missing    ! missing value
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
    real(r8), pointer :: arrayg(:)
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
    character(len=*),parameter :: subname='ncd_iolocal_gs_real1d' ! subroutine name
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
            baseTYPE = MPI_REAL8
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
               call check_ret(nf_put_vara_double(ncid, varid, start, count, arrayg), subname)
            else
               call check_ret(nf_put_var_double(ncid, varid, arrayg), subname)
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

  end subroutine ncd_iolocal_gs_real1d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_gs_real2d
!
! !INTERFACE:
  subroutine ncd_iolocal_gs_real2d(ncid, varname, flag, data, clmlevel, start, count, status, missing, usepio)
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
    real(r8),pointer              :: data(:,:)  ! local decomposition input data (out)
    character(len=*) ,intent(in)  :: clmlevel   ! type of grid
    integer, optional,intent(in)  :: start(:)   ! netcdf start index
    integer, optional,intent(in)  :: count(:)   ! netcdf count index
    integer, optional,intent(out) :: status     ! return code
    real(r8),optional,intent(in)  :: missing    ! missing value
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
    real(r8), pointer :: arrayg(:,:)
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
    character(len=*),parameter :: subname='ncd_iolocal_gs_real2d' ! subroutine name
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

  end subroutine ncd_iolocal_gs_real2d
!------------------------------------------------------------------------
!***** end include ncdio_local_subs.inc *****
!------------------------------------------------------------------------

!-----------------------------------------------------------------------
!***** begin include ncdio_global_subs.inc *****
!      this next section has been autogenerated from 
!      --- gen_ncdio_global_subs.csh ---
!      and manually included.  this section ends at "end include ..."
!-----------------------------------------------------------------------
!  #include <ncdio_global_subs.inc>

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_var
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_var(varname, data, flag, ncid, readvar, nt, bcast, usepio)
!
! !DESCRIPTION:
! netcdf I/O of global var int array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    integer         , intent(inout) :: data             ! raw data
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
    integer         ,pointer :: piodata(:)  ! copy of data in 1d
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
    integer,parameter :: ndims = 0  ! data dims
    character(len=*),parameter :: subname='ncd_ioglobal_int_var'
!-----------------------------------------------------------------------

    lusepio = pio_def
    if (present(usepio)) then
       lusepio = usepio
    endif

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(varname),lusepio
    endif

    if (lusepio) then
       n = 1
       allocate(piodata(n))
       piodata(1) = data
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
             count(n) = 1
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
             if (masterproc) call check_ret(nf_put_vara_int(ncid, varid, start, count, data), subname)
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
                   if (masterproc) call check_ret(nf_get_vara_int(ncid, varid, start, count, data), subname)
                endif
             else
                if (lusepio) then
#if (defined BUILDPIO)
!                   call check_ret_pio(PIO_get_var(ncid, varid, piodata), subname)
                    if (masterproc) write(iulog,*) trim(subname),' pio not implemented'
                    call endrun()
#endif
                else
                   if (masterproc) call check_ret(nf_get_var_int(ncid, varid, data), subname)
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
                call mpi_bcast(data, 1, MPI_INTEGER, 0, mpicom, ier)
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

  end subroutine ncd_ioglobal_int_var


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_var
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_var(varname, data, flag, ncid, readvar, nt, bcast, usepio)
!
! !DESCRIPTION:
! netcdf I/O of global var real array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    real(r8)        , intent(inout) :: data             ! raw data
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
    real(r8)        ,pointer :: piodata(:)  ! copy of data in 1d
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
    integer,parameter :: ndims = 0  ! data dims
    character(len=*),parameter :: subname='ncd_ioglobal_real_var'
!-----------------------------------------------------------------------

    lusepio = pio_def
    if (present(usepio)) then
       lusepio = usepio
    endif

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(varname),lusepio
    endif

    if (lusepio) then
       n = 1
       allocate(piodata(n))
       piodata(1) = data
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
             count(n) = 1
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
             if (masterproc) call check_ret(nf_put_vara_double(ncid, varid, start, count, data), subname)
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
                   if (masterproc) call check_ret(nf_get_vara_double(ncid, varid, start, count, data), subname)
                endif
             else
                if (lusepio) then
#if (defined BUILDPIO)
!                   call check_ret_pio(PIO_get_var(ncid, varid, piodata), subname)
                    if (masterproc) write(iulog,*) trim(subname),' pio not implemented'
                    call endrun()
#endif
                else
                   if (masterproc) call check_ret(nf_get_var_double(ncid, varid, data), subname)
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
                call mpi_bcast(data, 1, MPI_REAL8, 0, mpicom, ier)
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

  end subroutine ncd_ioglobal_real_var


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_1d
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_1d(varname, data, flag, ncid, readvar, nt, bcast, usepio)
!
! !DESCRIPTION:
! netcdf I/O of global 1d int array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    integer         , intent(inout) :: data (:)         ! raw data
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
    integer         ,pointer :: piodata(:)  ! copy of data in 1d
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
    integer,parameter :: ndims = 1  ! data dims
    character(len=*),parameter :: subname='ncd_ioglobal_int_1d'
!-----------------------------------------------------------------------

    lusepio = pio_def
    if (present(usepio)) then
       lusepio = usepio
    endif

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(varname),lusepio
    endif

    if (lusepio) then
       n = size(data)
       allocate(piodata(n))
       piodata = data
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
             count(n) = size(data,dim=n)
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
             if (masterproc) call check_ret(nf_put_vara_int(ncid, varid, start, count, data), subname)
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
                   if (masterproc) call check_ret(nf_get_vara_int(ncid, varid, start, count, data), subname)
                endif
             else
                if (lusepio) then
#if (defined BUILDPIO)
!                   call check_ret_pio(PIO_get_var(ncid, varid, piodata), subname)
                    if (masterproc) write(iulog,*) trim(subname),' pio not implemented'
                    call endrun()
#endif
                else
                   if (masterproc) call check_ret(nf_get_var_int(ncid, varid, data), subname)
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
                call mpi_bcast(data, size(data), MPI_INTEGER, 0, mpicom, ier)
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

  end subroutine ncd_ioglobal_int_1d


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_1d
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_1d(varname, data, flag, ncid, readvar, nt, bcast, usepio)
!
! !DESCRIPTION:
! netcdf I/O of global 1d real array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    real(r8)        , intent(inout) :: data (:)         ! raw data
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
    real(r8)        ,pointer :: piodata(:)  ! copy of data in 1d
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
    integer,parameter :: ndims = 1  ! data dims
    character(len=*),parameter :: subname='ncd_ioglobal_real_1d'
!-----------------------------------------------------------------------

    lusepio = pio_def
    if (present(usepio)) then
       lusepio = usepio
    endif

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(varname),lusepio
    endif

    if (lusepio) then
       n = size(data)
       allocate(piodata(n))
       piodata = data
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
             count(n) = size(data,dim=n)
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
             if (masterproc) call check_ret(nf_put_vara_double(ncid, varid, start, count, data), subname)
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
                   if (masterproc) call check_ret(nf_get_vara_double(ncid, varid, start, count, data), subname)
                endif
             else
                if (lusepio) then
#if (defined BUILDPIO)
!                   call check_ret_pio(PIO_get_var(ncid, varid, piodata), subname)
                    if (masterproc) write(iulog,*) trim(subname),' pio not implemented'
                    call endrun()
#endif
                else
                   if (masterproc) call check_ret(nf_get_var_double(ncid, varid, data), subname)
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
                call mpi_bcast(data, size(data), MPI_REAL8, 0, mpicom, ier)
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

  end subroutine ncd_ioglobal_real_1d


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_char_1d
!
! !INTERFACE:
  subroutine ncd_ioglobal_char_1d(varname, data, flag, ncid, readvar, nt, bcast, usepio)
!
! !DESCRIPTION:
! netcdf I/O of global 1d char array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    character(len=*), intent(inout) :: data       ! raw data
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
    character,pointer :: piodata(:)  ! copy of data in 1d
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
    integer,parameter :: ndims = 1  ! data dims
    character(len=*),parameter :: subname='ncd_ioglobal_char_1d'
!-----------------------------------------------------------------------

    lusepio = pio_def
    if (present(usepio)) then
       lusepio = usepio
    endif

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(varname),lusepio
    endif

    if (lusepio) then
       n = len(data)
       allocate(piodata(n))
       do n1 = 1,n
          piodata(n1) = data(n1:n1)
       enddo
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
             count(n) = len(data)
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
             if (masterproc) call check_ret(nf_put_vara_text(ncid, varid, start, count, data), subname)
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
                   if (masterproc) call check_ret(nf_get_vara_text(ncid, varid, start, count, data), subname)
                endif
             else
                if (lusepio) then
#if (defined BUILDPIO)
!                   call check_ret_pio(PIO_get_var(ncid, varid, piodata), subname)
                    if (masterproc) write(iulog,*) trim(subname),' pio not implemented'
                    call endrun()
#endif
                else
                   if (masterproc) call check_ret(nf_get_var_text(ncid, varid, data), subname)
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
                call mpi_bcast(data, len(data), MPI_CHARACTER, 0, mpicom, ier)
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

  end subroutine ncd_ioglobal_char_1d


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_2d
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_2d(varname, data, flag, ncid, readvar, nt, bcast, usepio)
!
! !DESCRIPTION:
! netcdf I/O of global 2d int array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    integer         , intent(inout) :: data (:,:)       ! raw data
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
    integer         ,pointer :: piodata(:)  ! copy of data in 1d
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
    integer,parameter :: ndims = 2  ! data dims
    character(len=*),parameter :: subname='ncd_ioglobal_int_2d'
!-----------------------------------------------------------------------

    lusepio = pio_def
    if (present(usepio)) then
       lusepio = usepio
    endif

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(varname),lusepio
    endif

    if (lusepio) then
       n = size(data)
       allocate(piodata(n))
       n = 0
       do n2 = 1, size(data,dim=2)
       do n1 = 1, size(data,dim=1)
          n = n + 1
          piodata(n) = data(n1,n2)
       enddo
       enddo
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
             count(n) = size(data,dim=n)
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
             if (masterproc) call check_ret(nf_put_vara_int(ncid, varid, start, count, data), subname)
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
                   if (masterproc) call check_ret(nf_get_vara_int(ncid, varid, start, count, data), subname)
                endif
             else
                if (lusepio) then
#if (defined BUILDPIO)
!                   call check_ret_pio(PIO_get_var(ncid, varid, piodata), subname)
                    if (masterproc) write(iulog,*) trim(subname),' pio not implemented'
                    call endrun()
#endif
                else
                   if (masterproc) call check_ret(nf_get_var_int(ncid, varid, data), subname)
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
                call mpi_bcast(data, size(data), MPI_INTEGER, 0, mpicom, ier)
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

  end subroutine ncd_ioglobal_int_2d


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_2d
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_2d(varname, data, flag, ncid, readvar, nt, bcast, usepio)
!
! !DESCRIPTION:
! netcdf I/O of global 2d real array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    real(r8)        , intent(inout) :: data (:,:)       ! raw data
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
    real(r8)        ,pointer :: piodata(:)  ! copy of data in 1d
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
    integer,parameter :: ndims = 2  ! data dims
    character(len=*),parameter :: subname='ncd_ioglobal_real_2d'
!-----------------------------------------------------------------------

    lusepio = pio_def
    if (present(usepio)) then
       lusepio = usepio
    endif

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(varname),lusepio
    endif

    if (lusepio) then
       n = size(data)
       allocate(piodata(n))
       n = 0
       do n2 = 1, size(data,dim=2)
       do n1 = 1, size(data,dim=1)
          n = n + 1
          piodata(n) = data(n1,n2)
       enddo
       enddo
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
             count(n) = size(data,dim=n)
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
             if (masterproc) call check_ret(nf_put_vara_double(ncid, varid, start, count, data), subname)
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
                   if (masterproc) call check_ret(nf_get_vara_double(ncid, varid, start, count, data), subname)
                endif
             else
                if (lusepio) then
#if (defined BUILDPIO)
!                   call check_ret_pio(PIO_get_var(ncid, varid, piodata), subname)
                    if (masterproc) write(iulog,*) trim(subname),' pio not implemented'
                    call endrun()
#endif
                else
                   if (masterproc) call check_ret(nf_get_var_double(ncid, varid, data), subname)
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
                call mpi_bcast(data, size(data), MPI_REAL8, 0, mpicom, ier)
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

  end subroutine ncd_ioglobal_real_2d


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_3d
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_3d(varname, data, flag, ncid, readvar, nt, bcast, usepio)
!
! !DESCRIPTION:
! netcdf I/O of global 3d int array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    integer         , intent(inout) :: data (:,:,:)     ! raw data
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
    integer         ,pointer :: piodata(:)  ! copy of data in 1d
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
    integer,parameter :: ndims = 3  ! data dims
    character(len=*),parameter :: subname='ncd_ioglobal_int_3d'
!-----------------------------------------------------------------------

    lusepio = pio_def
    if (present(usepio)) then
       lusepio = usepio
    endif

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(varname),lusepio
    endif

    if (lusepio) then
       n = size(data)
       allocate(piodata(n))
       n = 0
       do n3 = 1, size(data,dim=3)
       do n2 = 1, size(data,dim=2)
       do n1 = 1, size(data,dim=1)
          n = n + 1
          piodata(n) = data(n1,n2,n3)
       enddo
       enddo
       enddo
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
             count(n) = size(data,dim=n)
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
             if (masterproc) call check_ret(nf_put_vara_int(ncid, varid, start, count, data), subname)
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
                   if (masterproc) call check_ret(nf_get_vara_int(ncid, varid, start, count, data), subname)
                endif
             else
                if (lusepio) then
#if (defined BUILDPIO)
!                   call check_ret_pio(PIO_get_var(ncid, varid, piodata), subname)
                    if (masterproc) write(iulog,*) trim(subname),' pio not implemented'
                    call endrun()
#endif
                else
                   if (masterproc) call check_ret(nf_get_var_int(ncid, varid, data), subname)
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
                call mpi_bcast(data, size(data), MPI_INTEGER, 0, mpicom, ier)
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

  end subroutine ncd_ioglobal_int_3d


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_3d
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_3d(varname, data, flag, ncid, readvar, nt, bcast, usepio)
!
! !DESCRIPTION:
! netcdf I/O of global 3d real array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    real(r8)        , intent(inout) :: data (:,:,:)     ! raw data
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
    real(r8)        ,pointer :: piodata(:)  ! copy of data in 1d
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
    integer,parameter :: ndims = 3  ! data dims
    character(len=*),parameter :: subname='ncd_ioglobal_real_3d'
!-----------------------------------------------------------------------

    lusepio = pio_def
    if (present(usepio)) then
       lusepio = usepio
    endif

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(varname),lusepio
    endif

    if (lusepio) then
       n = size(data)
       allocate(piodata(n))
       n = 0
       do n3 = 1, size(data,dim=3)
       do n2 = 1, size(data,dim=2)
       do n1 = 1, size(data,dim=1)
          n = n + 1
          piodata(n) = data(n1,n2,n3)
       enddo
       enddo
       enddo
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
             count(n) = size(data,dim=n)
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
             if (masterproc) call check_ret(nf_put_vara_double(ncid, varid, start, count, data), subname)
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
                   if (masterproc) call check_ret(nf_get_vara_double(ncid, varid, start, count, data), subname)
                endif
             else
                if (lusepio) then
#if (defined BUILDPIO)
!                   call check_ret_pio(PIO_get_var(ncid, varid, piodata), subname)
                    if (masterproc) write(iulog,*) trim(subname),' pio not implemented'
                    call endrun()
#endif
                else
                   if (masterproc) call check_ret(nf_get_var_double(ncid, varid, data), subname)
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
                call mpi_bcast(data, size(data), MPI_REAL8, 0, mpicom, ier)
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

  end subroutine ncd_ioglobal_real_3d

!------------------------------------------------------------------------
!***** end include ncdio_global_subs.inc *****
!------------------------------------------------------------------------

#if (defined BUILDPIO)
!BOP
!
! !IROUTINE: subroutine ncd_setDOF
!
! !INTERFACE:
!------------------------------------------------------------------------
  subroutine ncd_setDOF(clmlevel,dims,compDOF,ioDOF,start,count)
!
! !DESCRIPTION: 
! Setup DOF and start/count arrays for pio
!
! !USES:
    use decompMod      , only : get_clmlevel_gsize,get_clmlevel_dsize
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: clmlevel   ! dimension 1 name
    integer         , intent(in)    :: dims(:)
    integer,pointer , intent(inout) :: compDOF(:)
    integer,pointer , intent(inout) :: ioDOF(:)
    integer(pio_offset),pointer , intent(inout) :: start(:)
    integer(pio_offset),pointer , intent(inout) :: count(:)
!
! !REVISION HISTORY:
! Created by T Craig, Aug 2007

!EOP
!
! !LOCAL VARIABLES:
    integer :: m,n,cnt,n1,n2,n3,n4,n1s,n2s,n3s,n4s
    integer :: ndims                         ! size of dims
    integer :: gsize                         ! global size of clmlevel
    integer :: gsmap_lsize                   ! local size of gsmap
    integer :: gsmap_gsize                   ! global size of gsmap
    integer :: fullsize                      ! size of entire array on cdf
    integer :: vsize                         ! other dimensions
    integer :: iosize                        ! local iodof size
    integer :: ier                           ! error status
    logical :: found                         ! found flag
    type(mct_gsMap),pointer       :: gsmap   ! global seg map
    integer, pointer,dimension(:) :: perm    ! mct permuter
    integer, pointer,dimension(:) :: gsmOP   ! gsmap ordered points
    character(len=32) :: subname = 'ncd_setDOF'
!-----------------------------------------------

    call get_clmlevel_gsmap(clmlevel,gsmap,perm)

    ndims = size(dims)
    gsize = get_clmlevel_gsize(clmlevel)
    gsmap_lsize = mct_gsmap_lsize(gsmap,mpicom)
    gsmap_gsize = mct_gsmap_gsize(gsmap)

    fullsize = 1
    do n = 1,ndims
       fullsize = fullsize*dims(n)
    enddo

    vsize = fullsize / gsize
    if (mod(fullsize,gsize) /= 0) then
       write(iulog,*) subname,' ERROR in vsize ',fullsize,gsize,vsize
       call endrun()
    endif

    allocate(start(4),count(4))
    start = 1
    count = 1
    count(1:ndims) = dims(1:ndims)

    found = .false.
    do n = ndims,1,-1
       if (.not.found .and. dims(n) >= npes) then
          found = .true.
          count(n) = dims(n) / npes
          start(n) = iam*count(n) + 1
          m = mod(dims(n),npes)
          start(n) = start(n) + min(iam,m)
          if (iam < m) count(n) = count(n) + 1
       endif
    enddo

    if (.not.found) then
       write(iulog,*) trim(subname),' ERROR: start,count not computed for npes = ', &
          npes,' and dims = ',dims
       call endrun()
    endif

    iosize = 1
    do n = 1,ndims
       iosize = iosize*count(n)
    enddo

    allocate(compDOF(gsmap_lsize*vsize))
    allocate(ioDOF(iosize))

    call mct_gsmap_OP(gsmap,iam,gsmOP)

    cnt = 0
    do n = 1,vsize
       do m = 1,gsmap_lsize
          cnt = cnt + 1
          compDOF(cnt) = (n-1)*gsmap_gsize + gsmOP(m)
       enddo
    enddo

    cnt = 0
    do n4 = start(4),start(4)+count(4)-1
       n4s = (n4-1)*dims(3)*dims(2)*dims(1)
       do n3 = start(3),start(3)+count(3)-1
          n3s = n4s + (n3-1)*dims(2)*dims(1)
          do n2 = start(2),start(2)+count(2)-1
             n2s = n3s + (n2-1)*dims(1)
             do n1 = start(1),start(1)+count(1)-1
                n1s = n2s + (n1)
                cnt = cnt + 1
                ioDOF(cnt) = n1s
             enddo
          enddo
       enddo
    enddo

    do m = 0,npes-1
    if (iam == m) then
!tcx    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' sizes1 = ',iam,gsize,gsmap_gsize,gsmap_lsize
       write(iulog,*) trim(subname),' sizes2 = ',iam,fullsize,npes,vsize,iosize
       write(iulog,*) trim(subname),' dims = ',iam,(start(n),count(n),dims(n),n=1,4)
       write(iulog,*) trim(subname),' compDOF = ',iam,size(compDOF),minval(compDOF),maxval(compDOF)
       write(iulog,*) trim(subname),' ioDOF = ',iam,size(ioDOF),minval(ioDOF),maxval(ioDOF)
       call shr_sys_flush(iulog)
!    endif
    endif
    call mpi_barrier(mpicom,ier)
    enddo

    deallocate(gsmOP)

 end subroutine ncd_setDOF
#endif

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
    write(iulog,*) trim(subname),' scam_setlatlonidx ',lonidx,latidx

    if ( firsttime) then
       write(iulog,*) trim(subname),' firsttime=',firsttime
       call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,latidx,lonidx)

       call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
       write(iulog,*) trim(subname),' beg,end=',begg, endg, begl, endl, begc, endc, begp, endp

       ret = nf_inq_dimid (ncid, 'column', dimid)
       write(iulog,*) trim(subname),' column ret=',ret,nf_noerr
       if(ret==NF_NOERR) return
       if (ret/=NF_EBADDIM) print *, 'NETCDF ERROR: ', NF_STRERROR(ret)

       ret = nf_inq_dimlen (ncid, dimid, totcols)
       if (ret/=NF_NOERR) print *, 'NETCDF ERROR: ', NF_STRERROR(ret)

       ret = nf_inq_dimid (ncid, 'pft', dimid)
       write(iulog,*) trim(subname),' pft ret=',ret,nf_noerr
       if(ret==NF_NOERR) return
       if (ret/=NF_EBADDIM) print *, 'NETCDF ERROR: ', NF_STRERROR(ret)

       ret = nf_inq_dimlen (ncid, dimid, totpfts)
       if (ret/=NF_NOERR) print *, 'NETCDF ERROR: ', NF_STRERROR(ret)

       write(iulog,*) trim(subname),' totals ',totcols,totpfts

       allocate (pfts1dlon(totpfts))
       allocate (pfts1dlat(totpfts))
       
       ret =  nf_inq_varid (ncid, 'pfts1d_ixy', varid)
       if (ret/=NF_NOERR) then
         write(iulog,*)'inq_varid: id for pfts1d_ixy not found'
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
         return
       end if

       ret = nf_get_var_double (ncid, varid, pfts1dlon)
       if (ret/=NF_NOERR) then
         write(iulog,*)'GET_VAR_REALX: error reading pfts1dlon, varid =', varid
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
         return
       end if

       ret =  nf_inq_varid (ncid, 'pfts1d_jxy', varid)
       if (ret/=NF_NOERR) then
         write(iulog,*)'inq_varid: id for pfts1d_jxy not found'
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
         return
       end if

       ret = nf_get_var_double (ncid, varid, pfts1dlat)
       if (ret/=NF_NOERR) then
         write(iulog,*)'GET_VAR_REALX: error reading pfts1dlat, varid =', varid
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
         return
       end if


       allocate (cols1dlon(totcols))
       allocate (cols1dlat(totcols))

       ret =  nf_inq_varid (ncid, 'cols1d_ixy', varid)
       if (ret/=NF_NOERR) then
         write(iulog,*)'inq_varid: id for cols1d_ixy not found'
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
         return
       end if

       ret = nf_get_var_double (ncid, varid, cols1dlon)
       if (ret/=NF_NOERR) then
         write(iulog,*)'GET_VAR_REALX: error reading cols1dlon, varid =', varid
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
         return
       end if

       ret =  nf_inq_varid (ncid, 'cols1d_jxy', varid)
       if (ret/=NF_NOERR) then
         write(iulog,*)'inq_varid: id for cols1d_jxy not found'
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
         return
       end if

       ret = nf_get_var_double (ncid, varid, cols1dlat)
       if (ret/=NF_NOERR) then
         write(iulog,*)'GET_VAR_REALX: error reading cols1dlat, varid =', varid
         print *, 'NETCDF ERROR: ', NF_STRERROR(ret)
         return
       end if

       cols(:)=-9999
       pfts(:)=-9999
       col_offset=-9999
       pi_offset=-9999
       i=1
       do cc = 1, totcols
          if (cols1dlon(cc).eq.lonidx.and.cols1dlat(cc).eq.latidx) then
             cols(i)=cc
             i=i+1
          end if
       end do
       if (endc-begc+1.ne.i-1) then
          write(iulog,*)'error in number of columns read for this gridcell',endc,begc,i
!          call endrun
       end if
       if (i.eq.1) then
          write(iulog,*)'couldnt find any columns for this latitude ',latidx,' and longitude ',lonidx
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
          write(iulog,*)'error in number of pfts read for this gridcell',endp,begp,i
!          call endrun
       end if
       if (i.eq.1) then
          write(iulog,*)'couldnt find any pfts for this latitude ',latidx,' and longitude ',lonidx
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

    write(iulog,*) trim(subname),' offsets ',pi_offset,col_offset
    
    if (dim1name == 'pft') then
       data_offset = pi_offset
       ndata = endp-begp+1
    else if (dim1name == 'column') then
       data_offset = col_offset
       ndata = endc-begc+1
    else
       write(iulog,*)'error calculation array offsets for SCAM'
!       call endrun()
    endif
  end subroutine scam_field_offsets
!------------------------------------------------------------------------

end module ncdio
