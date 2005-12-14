#include <misc.h>
#include <preproc.h>

module histFileMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: histFileMod
!
! !DESCRIPTION:
! Module containing methods to for CLM history file handling.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_getenv, shr_sys_flush
  use abortutils  , only : endrun
  use clm_varcon  , only : spval
  use shr_sys_mod , only : shr_sys_flush 
  implicit none
  save
  private
  include "netcdf.inc"
!
! !PUBLIC TYPES:
!
! Constants
!
  integer , public, parameter :: max_tapes = 6          ! max number of history tapes
  integer , public, parameter :: max_flds = 1000        ! max number of history fields
  integer , public, parameter :: max_namlen = 32        ! maximum number of characters for field name
!
! Counters
!
  integer , public :: ntapes = 0         ! index of max history file requested
!
! Namelist
!
  integer :: ni                          ! implicit index below
  logical, public :: &
       hist_empty_htapes  = .false.      ! namelist: flag indicates no default history fields
  integer, public :: &
       hist_ndens(max_tapes) = 2         ! namelist: output density of netcdf history files
  integer, public :: &
       hist_mfilt(max_tapes) = 30        ! namelist: number of time samples per tape
  logical, public :: &
       hist_dov2xy(max_tapes) = (/.true.,(.true.,ni=2,max_tapes)/) ! namelist: true=> do grid averaging
  integer, public :: &
       hist_nhtfrq(max_tapes) = (/0, (-24, ni=2,max_tapes)/)        ! namelist: history write freq(0=monthly)
  character(len=1), public :: &
       hist_avgflag_pertape(max_tapes) = (/(' ',ni=1,max_tapes)/)   ! namelist: per tape averaging flag
  character(len=max_namlen), public :: &
       hist_type1d_pertape(max_tapes)  = (/(' ',ni=1,max_tapes)/)   ! namelist: per tape type1d

  character(len=max_namlen+2), public :: &
       fincl(max_flds,max_tapes) = ' '   ! namelist-equivalence list of fields to add

  character(len=max_namlen+2), public :: &
       hist_fincl1(max_flds)        ! namelist: list of fields to add
  character(len=max_namlen+2), public :: &
       hist_fincl2(max_flds)        ! namelist: list of fields to add
  character(len=max_namlen+2), public :: &
       hist_fincl3(max_flds)        ! namelist: list of fields to add
  character(len=max_namlen+2), public :: &
       hist_fincl4(max_flds)        ! namelist: list of fields to add
  character(len=max_namlen+2), public :: &
       hist_fincl5(max_flds)        ! namelist: list of fields to add
  character(len=max_namlen+2), public :: &
       hist_fincl6(max_flds)        ! namelist: list of fields to add

  character(len=max_namlen), public :: &
       fexcl(max_flds,max_tapes) = ' '   ! namelist-equivalence list of fields to remove

  character(len=max_namlen), public :: &
       hist_fexcl1(max_flds)  ! namelist: list of fields to remove
  character(len=max_namlen), public :: &
       hist_fexcl2(max_flds)  ! namelist: list of fields to remove
  character(len=max_namlen), public :: &
       hist_fexcl3(max_flds)  ! namelist: list of fields to remove
  character(len=max_namlen), public :: &
       hist_fexcl4(max_flds)  ! namelist: list of fields to remove
  character(len=max_namlen), public :: &
       hist_fexcl5(max_flds)  ! namelist: list of fields to remove
  character(len=max_namlen), public :: &
       hist_fexcl6(max_flds)  ! namelist: list of fields to remove
!
! Equivalence used to satisfy namelist on a wide variety of platforms
! NOTE: It is *ASSUMED* that max_tapes is 6
!
  equivalence (hist_fincl1,fincl(1,1))
  equivalence (hist_fincl2,fincl(1,2))
  equivalence (hist_fincl3,fincl(1,3))
  equivalence (hist_fincl4,fincl(1,4))
  equivalence (hist_fincl5,fincl(1,5))
  equivalence (hist_fincl6,fincl(1,6))

  equivalence (hist_fexcl1,fexcl(1,1))
  equivalence (hist_fexcl2,fexcl(1,2))
  equivalence (hist_fexcl3,fexcl(1,3))
  equivalence (hist_fexcl4,fexcl(1,4))
  equivalence (hist_fexcl5,fexcl(1,5))
  equivalence (hist_fexcl6,fexcl(1,6))
!
! Restart
!
  logical, public :: if_writrest    ! true=> write restart file now
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: add_fld1d               ! Add a 1d single-level field to the master field list
  public :: add_fld2d               ! Add a 2d multi-level field to the master field list
  public :: add_subscript           ! Add a 2d subscript dimension
  public :: masterlist_make_active  ! Add a field to the given history file default "on" list
  public :: masterlist_printflds    ! Print summary of master field list
  public :: htapes_build            ! Initialize history file handler for initial or continue run
  public :: update_hbuf             ! Updates history buffer for all fields and tapes
  public :: htapes_wrapup           ! Write and/or dispose history tape(s)
  public :: restart_history         ! Read/write history file restart data
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! PRIVATE MEMBER FUNCTIONS:
  private :: masterlist_addfld         ! Add a field to the master field list
  private :: masterlist_change_timeavg ! Override default history tape contents for specific tape
  private :: htapes_fieldlist          ! Define the contents of each history file based on namelist
  private :: htape_addfld              ! Add a field to the active list for a history tape
  private :: htape_create              ! Define contents of history file t
  private :: htape_timeconst           ! Write time constant values to primary history tape
  private :: hfields_normalize         ! Normalize history file fields by number of accumulations
  private :: hfields_zero              ! Zero out accumulation and hsitory buffers for a tape
  private :: hfields_write             ! Write a variable to a history tape
  private :: hfields_1dinfo            ! Define/output 1d subgrid info if appropriate
  private :: update_hbuf_field_1d      ! Updates history buffer for specific field and tape
  private :: update_hbuf_field_2d      ! Updates history buffer for specific field and tape 
  private :: list_index                ! Find index of field in exclude list
  private :: getname                   ! Retrieve name portion of input "inname"
  private :: set_hist_filename         ! Determine history dataset filenames

! PRIVATE TYPES:
! Constants
!
  integer, parameter :: max_chars = 128        ! max chars for char variables
!
! Subscript dimensions
!
  integer, parameter :: max_subs = 100         ! max number of subscripts
  integer            :: num_subs = 0
  character(len=32)  :: subs_name(max_subs)
  integer            :: subs_dim(max_subs)
!
! Derived types
!
  type field_info
     character(len=max_namlen) :: name         ! field name
     character(len=max_chars)  :: long_name    ! long name
     character(len=max_chars)  :: units        ! units
     character(len=8) :: type1d                ! clm pointer first dimension type type ["gridcell","landunit","column","pft","roflnd","rofocn"]
     character(len=8) :: type1d_out            ! hbuf first dimension type ["gridcell","landunit","column","pft","roflnd","rofocn"]
     character(len=8) :: type2d                ! hbuf second dimension type ["levsoi","levlak","numrad","subname(n)"]
     character(len=8) :: typexy                ! xy grid type ["clm", "rof"]
     integer :: beg1d                          ! on-node 1d clm pointer start index
     integer :: end1d                          ! on-node 1d clm pointer end index
     integer :: num1d                          ! size of clm pointer first dimension (overall all nodes)
     integer :: beg1d_out                      ! on-node 1d hbuf pointer start index
     integer :: end1d_out                      ! on-node 1d hbuf pointer end index
     integer :: num1d_out                      ! size of hbuf first dimension (overall all nodes)
     integer :: num2d                          ! size of hbuf second dimension (e.g. number of vertical levels)
     integer :: nlonxy                         ! xy grid longitude size
     integer :: nlatxy                         ! xy grid latitude size
     integer :: hpindex                        ! history pointer index selected out of possible list
     character(len=8) :: p2c_scale_type        ! scale factor when averaging pft to column
     character(len=8) :: c2l_scale_type        ! scale factor when averaging column to landunit
     character(len=8) :: l2g_scale_type        ! scale factor when averaging landunit to gridcell
  end type field_info

  type master_entry
     type (field_info)  :: field               ! field information
     logical            :: actflag(max_tapes)  ! active/inactive flag
     character(len=1)   :: avgflag(max_tapes)  ! time averaging flag ("X","A","M" or "I",)
  end type master_entry

  type history_entry
     type (field_info) :: field                ! field information
     character(len=1)  :: avgflag              ! time averaging flag
     real(r8), pointer :: hbuf(:,:)            ! history buffer (dimensions: dim1d x num2d)
     integer , pointer :: nacs(:,:)            ! accumulation counter (dimensions: dim1d x num2d)
  end type history_entry

  type history_tape
     integer  :: nflds                         ! number of active fields on tape
     integer  :: ntimes                        ! current number of time samples on tape
     integer  :: mfilt                         ! maximum number of time samples per tape
     integer  :: nhtfrq                        ! number of time samples per tape
     integer  :: ncprec                        ! netcdf output precision
     logical  :: dov2xy                        ! true => do xy average for all fields
     logical  :: is_endhist                    ! true => current time step is end of history interval
     real(r8) :: begtime                       ! time at beginning of history averaging interval
     type (history_entry) :: hlist(max_flds)   ! array of active history tape entries
  end type history_tape

  type clmpoint_is
     integer, pointer :: ptr(:)
  end type clmpoint_is
  type clmpoint_ia
     integer, pointer :: ptr(:,:)
  end type clmpoint_ia
  type clmpoint_rs
     real(r8), pointer :: ptr(:)
  end type clmpoint_rs
  type clmpoint_ra
     real(r8), pointer :: ptr(:,:)
  end type clmpoint_ra
!
! Pointers into clmtype arrays
!
  integer, parameter :: max_mapflds = 1000
  type (clmpoint_rs) :: clmptr_rs(max_mapflds)
  type (clmpoint_ra) :: clmptr_ra(max_mapflds)
!
! Master list: an array of master_entry entities
!
  type (master_entry) :: masterlist(max_flds)  ! master field list
!
! History tape: an array of history_tape entities (only active fields)
!
  type (history_tape) :: tape(max_tapes)       ! array history tapes
!
! Namelist input
!
! Counters
!
  integer :: nfmaster = 0                        ! number of fields in master field list
!
! Other variables
!
  character(len=16) :: host                      ! host name
  character(len= 8) :: logname                   ! user name
  character(len=max_chars) :: locfnh(max_tapes)  ! local history file names
  logical :: htapes_defined = .false.            ! flag indicates history contents have been defined
!
! NetCDF  Id's
!
  integer :: nfid(max_tapes)                 ! file ids
  integer :: time_dimid                      ! time dimension id
  integer :: hist_interval_dimid             ! time bounds dimension id
  integer :: strlen_dimid                    ! string dimension id
  integer :: mcdate_id(max_tapes)            ! current date, yyyymmdd format (mcdate) ids
  integer :: mcsec_id(max_tapes)             ! current seconds in day (mcsec) ids
  integer :: mdcur_id(max_tapes)             ! current day (from base day) (mdcur) ids
  integer :: mscur_id(max_tapes)             ! current seconds of current day (mdcur) ids
  integer :: nstep_id(max_tapes)             ! current nstep ids
  integer :: time_var_id(max_tapes)          ! time coordinate variable ids
  integer :: time_bounds_id(max_tapes)       ! time bounds ids
  integer :: date_written_id(max_tapes)      ! current system date ids
  integer :: time_written_id(max_tapes)      ! current system time ids
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: masterlist_printflds
!
! !INTERFACE:
  subroutine masterlist_printflds()
!
! !DESCRIPTION:
! Print summary of master field list.
!
! !USES:
    use spmdMod, only : masterproc
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 03/2003
!
!EOP
!
! LOCAL VARIABLES:
    integer nf
!-----------------------------------------------------------------------

    if (masterproc) then
       write(6,*)'MASTERLIST_BUILD: number of master fields = ',nfmaster
       write(6,*)' ******* MASTER FIELD LIST *******'
       do nf = 1,nfmaster
          write(6,9000)nf, masterlist(nf)%field%name, masterlist(nf)%field%units
9000      format (i5,1x,a32,1x,a16)
       end do
    end if

  end subroutine masterlist_printflds

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: masterlist_addfld
!
! !INTERFACE:
  subroutine masterlist_addfld (fname, type1d, type2d, num2d,  &
        typexy, nlonxy, nlatxy, &
        units, avgflag, long_name, hpindex, &
        p2c_scale_type, c2l_scale_type, l2g_scale_type)
!
! !DESCRIPTION:
! Add a field to the master field list. Put input arguments of
! field name, units, number of levels, averaging flag, and long name
! into a type entry in the global master field list (masterlist).
!
! !USES:
    use clmtype  , only : nameg, namel, namec, namep, lndrof, ocnrof
    use decompMod, only : get_proc_bounds, get_proc_global
#if (defined RTM)
    use RunoffMod, only : get_proc_rof_bounds, get_proc_rof_global
#endif
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fname        ! field name
    character(len=*), intent(in) :: type1d       ! 1d output type
    character(len=*), intent(in) :: type2d       ! 2d output type
    character(len=*), intent(in) :: typexy       ! xy grid output type
    integer         , intent(in) :: nlonxy       ! xy grid longitude dimension
    integer         , intent(in) :: nlatxy       ! xy grid latitude dimension
    integer         , intent(in) :: num2d        ! size of second dimension (e.g. number of vertical levels)
    character(len=*), intent(in) :: units        ! units of field
    character(len=1), intent(in) :: avgflag      ! time averaging flag
    character(len=*), intent(in) :: long_name    ! long name of field
    integer         , intent(in) :: hpindex      ! clmtype index for history buffer output
    character(len=*), intent(in) :: p2c_scale_type ! scale type for subgrid averaging of pfts to column
    character(len=*), intent(in) :: c2l_scale_type ! scale type for subgrid averaging of columns to landunits
    character(len=*), intent(in) :: l2g_scale_type ! scale type for subgrid averaging of landunits to gridcells
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: n            ! loop index
    integer :: f            ! masterlist index
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer :: numg         ! total number of gridcells across all processors
    integer :: numl         ! total number of landunits across all processors
    integer :: numc         ! total number of columns across all processors
    integer :: nump         ! total number of pfts across all processors
    integer :: beg_lndrof, end_lndrof, num_lndrof  ! land  runoff bounds
    integer :: beg_ocnrof, end_ocnrof, num_ocnrof  ! ocean runoff bounds
!------------------------------------------------------------------------

    ! Determine bounds

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)
#if (defined RTM)
    call get_proc_rof_bounds(beg_lndrof, end_lndrof, beg_ocnrof, end_ocnrof)
    call get_proc_rof_global(num_lndrof, num_ocnrof)
#endif
    ! Ensure that new field is not all blanks

    if (fname == ' ') then
       write(6,*)'error masterlist_addfld: blank field name not allowed'
       call endrun()
    end if

    ! Ensure that new field doesn't already exist

    do n = 1,nfmaster
       if (masterlist(n)%field%name == fname) then
          write(6,*)'error masterlist_addfld:', fname, ' already on list'
          call endrun()
       end if
    end do

    ! Increase number of fields on master field list

    nfmaster = nfmaster + 1
    f = nfmaster

    ! Check number of fields in master list against maximum number for master list

    if (nfmaster > max_flds) then
       write(6,*)'error masterlist_addfld: too many fields for primary history file ', &
            '-- max_flds,nfmaster=', max_flds, nfmaster
       call endrun()
    end if

    ! Add field to master list

    masterlist(f)%field%name           = fname
    masterlist(f)%field%long_name      = long_name
    masterlist(f)%field%units          = units
    masterlist(f)%field%type1d         = type1d
    masterlist(f)%field%type1d_out     = type1d
    masterlist(f)%field%type2d         = type2d
    masterlist(f)%field%num2d          = num2d
    masterlist(f)%field%typexy         = typexy
    masterlist(f)%field%nlonxy         = nlonxy
    masterlist(f)%field%nlatxy         = nlatxy
    masterlist(f)%field%hpindex        = hpindex
    masterlist(f)%field%p2c_scale_type = p2c_scale_type
    masterlist(f)%field%c2l_scale_type = c2l_scale_type
    masterlist(f)%field%l2g_scale_type = l2g_scale_type

    select case (type1d)
    case (nameg)
       masterlist(f)%field%beg1d = begg
       masterlist(f)%field%end1d = endg
       masterlist(f)%field%num1d = numg
    case (namel)
       masterlist(f)%field%beg1d = begl
       masterlist(f)%field%end1d = endl
       masterlist(f)%field%num1d = numl
    case (namec)
       masterlist(f)%field%beg1d = begc
       masterlist(f)%field%end1d = endc
       masterlist(f)%field%num1d = numc
    case (namep)
       masterlist(f)%field%beg1d = begp
       masterlist(f)%field%end1d = endp
       masterlist(f)%field%num1d = nump
#if (defined RTM)
    case (lndrof)
       masterlist(f)%field%beg1d = beg_lndrof
       masterlist(f)%field%end1d = end_lndrof
       masterlist(f)%field%num1d = num_lndrof
    case (ocnrof)
       masterlist(f)%field%beg1d = beg_ocnrof
       masterlist(f)%field%end1d = end_ocnrof
       masterlist(f)%field%num1d = num_ocnrof
#endif
    case default
       write(6,*)'MASTERLIST_ADDFLD: unknown 1d output type= ',type1d
       call endrun()
    end select

    ! The following two fields are used only in master field list,
    ! NOT in the runtime active field list
    ! ALL FIELDS IN THE MASTER LIST ARE INITIALIZED WITH THE ACTIVE
    ! FLAG SET TO FALSE

    masterlist(f)%avgflag(:) = avgflag
    masterlist(f)%actflag(:) = .false.

  end subroutine masterlist_addfld

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: htapes_build
!
! !INTERFACE:
  subroutine htapes_build ()
!
! !DESCRIPTION:
! Initialize history file for initial or continuation run.  For example,
! on an initial run, this routine initializes ``ntapes'' history files.
! On a restart run, this routine only initializes history files declared
! beyond what existed on the previous run.  Files which already existed on
! the previous run have already been initialized (i.e. named and opened)
! in routine restart\_history.  Loop over tapes and fields per tape setting
! appropriate variables and calling appropriate routines
!
! !USES:
    use spmdMod, only : masterproc
    use time_manager, only: get_curr_date, get_curr_time
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: i                   ! index
    integer :: ier                 ! error code
    integer :: t, f                ! tape, field indices
    integer :: day, sec            ! day and seconds from base date
!-----------------------------------------------------------------------

    if (masterproc) then
       write(6,*) 'Initializing clm2 history files'
       write(6,'(72a1)') ("-",i=1,60)
    endif

    ! Get users logname and machine hostname

    if ( masterproc )then
       logname = ' '
       call shr_sys_getenv('LOGNAME', logname, ier)
       if (ier /= 0 .or. logname=='        ') then
          write(6,*)'HTAPES_BUILD error: Cannot find LOGNAME environment variable'
          call endrun()
       end if
       host = ' '
       call shr_sys_getenv('HOST', host, ier)
    end if

    ! Override averaging flag for all fields on a particular tape
    ! if namelist input so specifies

    do t=1,max_tapes
       if (hist_avgflag_pertape(t) /= ' ') then
          call masterlist_change_timeavg (t)
       end if
    end do

    ! Define field list information for all history files.
    ! Update ntapes to reflect number of active history files
    ! (note, branch runs can have additional auxiliary history files
    ! declared).

    call htapes_fieldlist()

    ! Determine elapased time since reference date

    call get_curr_time(day, sec)

    ! Set number of time samples in each history file and
    ! time of a beginning of current averaging interval.
    ! (note - the following entries will be overwritten by history restart)
    ! Determine if xy averaging is done for all fields on the tape, etc
    ! Note - with netcdf, only 1 (nf_double) and 2 (nf_float) are allowed

    do t=1,ntapes
       tape(t)%ntimes = 0
       tape(t)%begtime = day + sec/86400._r8
       tape(t)%dov2xy = hist_dov2xy(t)
       tape(t)%nhtfrq = hist_nhtfrq(t)
       if (hist_nhtfrq(t) == 0) hist_mfilt(t) = 1
       tape(t)%mfilt = hist_mfilt(t)
       if (hist_ndens(t) == 1) then
          tape(t)%ncprec = nf_double
       else
          tape(t)%ncprec = nf_float
       endif
    end do

    if (masterproc) then
       write(6,*) 'Successfully initialized clm2 history files'
       write(6,'(72a1)') ("-",i=1,60)
    endif

  end subroutine htapes_build

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: masterlist_make_active
!
! !INTERFACE:
  subroutine masterlist_make_active (name, tape_index, avgflag)
!
! !DESCRIPTION:
! Add a field to the default ``on'' list for a given history file.
! Also change the default time averaging flag if requested.
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: name          ! field name
    integer, intent(in) :: tape_index             ! history tape index
    character(len=1), intent(in), optional :: avgflag  ! time averaging flag
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: f            ! field index
    logical :: found        ! flag indicates field found in masterlist
!-----------------------------------------------------------------------

    ! Check validity of input arguments

    if (tape_index > max_tapes) then
       write(6,*)'MASTERLIST_MAKE_ACTIVE error: tape index=', tape_index, ' is too big'
       call endrun()
    end if

    if (present(avgflag)) then
       if ( avgflag /= ' ' .and. &
            avgflag /= 'A' .and. avgflag /= 'I' .and. &
            avgflag /= 'X' .and. avgflag /= 'M') then
          write(6,*)'MASTERLIST_MAKE_ACTIVE error: unknown averaging flag=', avgflag
          call endrun()
       endif
    end if

    ! Look through master list for input field name.
    ! When found, set active flag for that tape to true.
    ! Also reset averaging flag if told to use other than default.

    found = .false.
    do f = 1,nfmaster
       if (trim(name) == trim(masterlist(f)%field%name)) then
          masterlist(f)%actflag(tape_index) = .true.
          if (present(avgflag)) then
             if (avgflag/= ' ') masterlist(f)%avgflag(tape_index) = avgflag
          end if
          found = .true.
          exit
       end if
    end do
    if (.not. found) then
       write(6,*)'MASTERLIST_MAKE_ACTIVE error: field=', name, ' not found'
       call endrun()
    end if

  end subroutine masterlist_make_active

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: masterlist_change_timeavg
!
! !INTERFACE:
  subroutine masterlist_change_timeavg (t)
!
! !DESCRIPTION:
! Override default history tape contents for a specific tape.
! Copy the flag into the master field list.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t         ! history tape index
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: f                     ! field index
    character(len=1) :: avgflag      ! lcl equiv of hist_avgflag_pertape(t)
!-----------------------------------------------------------------------

    avgflag = hist_avgflag_pertape(t)

    do f = 1,nfmaster
       select case (avgflag)
       case ('A')
          masterlist(f)%avgflag(t) = avgflag
       case ('I')
          masterlist(f)%avgflag(t) = avgflag
       case ('X')
          masterlist(f)%avgflag(t) = avgflag
       case ('M')
          masterlist(f)%avgflag(t) = avgflag
       case default
          write(6,*)'MASTERLIST_CHANGE_TIMEAVG error: unknown avgflag=',avgflag
          call endrun ()
       end select
    end do

  end subroutine masterlist_change_timeavg

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: htapes_fieldlist
!
! !INTERFACE:
  subroutine htapes_fieldlist()
!
! !DESCRIPTION:
! Define the contents of each history file based on namelist
! input for initial or branch run, and restart data if a restart run.
! Use arrays fincl and fexcl to modify default history tape contents.
! Then sort the result alphanumerically.
!
! !USES:
    use spmdMod, only : masterproc
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: t, f                ! tape, field indices
    integer :: ff                  ! index into include, exclude and fprec list
    character(len=max_namlen) :: name       ! field name portion of fincl (i.e. no avgflag separator)
    character(len=max_namlen) :: mastername ! name from masterlist field
    character(len=1)  :: avgflag    ! averaging flag
    character(len=1)  :: prec_acc   ! history buffer precision flag
    character(len=1)  :: prec_wrt   ! history buffer write precision flag
    type (history_entry) :: tmp     ! temporary used for swapping
!-----------------------------------------------------------------------

    ! First ensure contents of fincl and fexcl are valid names

    do t = 1,max_tapes
       f = 1
       do while (f < max_flds .and. fincl(f,t) /= ' ')
          name = getname (fincl(f,t))
          do ff = 1,nfmaster
             mastername = masterlist(ff)%field%name
             if (name == mastername) exit
          end do
          if (name /= mastername) then
             write(6,*)'HTAPES_FIELDLIST error: ', trim(name), ' in fincl(', f, ') ',&
                  'for history tape ',t,' not found'
             call endrun()
          end if
          f = f + 1
       end do

       f = 1
       do while (f < max_flds .and. fexcl(f,t) /= ' ')
          do ff = 1,nfmaster
             mastername = masterlist(ff)%field%name
             if (fexcl(f,t) == mastername) exit
          end do
          if (fexcl(f,t) /= mastername) then
             write(6,*)'HTAPES_FIELDLIST error: ', fexcl(f,t), ' in fexcl(', f, ') ', &
                  'for history tape ',t,' not found'
             call endrun()
          end if
          f = f + 1
       end do
    end do

    tape(:)%nflds = 0
    do t = 1,max_tapes

       ! Loop through the masterlist set of field names and determine if any of those
       ! are in the FINCL or FEXCL arrays
       ! The call to list_index determines the index in the FINCL or FEXCL arrays
       ! that the masterlist field corresponds to
       ! Add the field to the tape if specified via namelist (FINCL[1-max_tapes]),
       ! or if it is on by default and was not excluded via namelist (FEXCL[1-max_tapes]).

       do f = 1,nfmaster
          mastername = masterlist(f)%field%name
          call list_index (fincl(1,t), mastername, ff)

          if (ff > 0) then

             ! if field is in include list, ff > 0 and htape_addfld
             ! will not be called for field

             avgflag = getflag (fincl(ff,t))
             call htape_addfld (t, f, avgflag)

          else if (.not. hist_empty_htapes) then

             ! find index of field in exclude list

             call list_index (fexcl(1,t), mastername, ff)

             ! if field is in exclude list, ff > 0 and htape_addfld
             ! will not be called for field
             ! if field is not in exclude list, ff =0 and htape_addfld
             ! will be called for field (note that htape_addfld will be
             ! called below only if field is not in exclude list OR in
             ! include list

             if (ff == 0 .and. masterlist(f)%actflag(t)) then
                call htape_addfld (t, f, ' ')
             end if

          end if
       end do

       ! Specification of tape contents now complete.
       ! Sort each list of active entries

       do f = tape(t)%nflds-1,1,-1
          do ff = 1,f
             if (tape(t)%hlist(ff)%field%name > tape(t)%hlist(ff+1)%field%name) then

                tmp = tape(t)%hlist(ff)
                tape(t)%hlist(ff  ) = tape(t)%hlist(ff+1)
                tape(t)%hlist(ff+1) = tmp

             else if (tape(t)%hlist(ff  )%field%name == tape(t)%hlist(ff+1)%field%name) then

                write(6,*)'HTAPES_FIELDLIST error: Duplicate field ', tape(t)%hlist(ff  )%field%name
                write(6,*)'t,ff,name=',t,ff,tape(t)%hlist(ff  )%field%name
                call endrun()

             end if
          end do
       end do

       if (masterproc) then
          if (tape(t)%nflds > 0) then
             write(6,*)'HTAPES_FIELDLIST: Included fields tape ',t,'=',tape(t)%nflds
          end if
          do f = 1,tape(t)%nflds
             write(6,*) f,' ',tape(t)%hlist(f)%field%name, &
                  tape(t)%hlist(f)%field%num2d,' ',tape(t)%hlist(f)%avgflag
          end do
       end if
    end do

    ! Determine total number of active history tapes

    ntapes = 0
    do t = max_tapes,1,-1
       if (tape(t)%nflds > 0) then
          ntapes = t
          exit
       end if
    end do

    ! Ensure there are no "holes" in tape specification, i.e. empty tapes.
    ! Enabling holes should not be difficult if necessary.

    do t = 1,ntapes
       if (tape(t)%nflds  ==  0) then
          write(6,*)'HTAPES_FIELDLIST error: Tape ',t,' is empty'
          call endrun()
       end if
    end do

    ! Check that the number of history files declared does not exceed
    ! the maximum allowed.

    if (ntapes > max_tapes) then
       write(6,*)'HTAPES_FIELDLIST error: Too many history files declared, max=',max_tapes
       write(6,*)'To increase, change parameter max_tapes.'
       call endrun()
    end if

    ! Change 1d output per tape output flag if requested - only for history
    ! tapes greater than 1 where 2d xy averaging is not enabled

    do t = 2,ntapes
       if (hist_type1d_pertape(t) /= ' ' .and. (.not. hist_dov2xy(t))) then
          select case (trim(hist_type1d_pertape(t)))
          case ('PFTS','COLS', 'LAND', 'GRID')
             write(6,*)'history tape ',t,' will have 1d output type of ',hist_type1d_pertape(t)
          case default
             write(6,*)'HTAPES_FIELDLIST error : unknown namelist type1d per tape=',hist_type1d_pertape(t)
             call endrun()
          end select
       end if
    end do

    if (masterproc) then
       write(6,*) 'There will be a total of ',ntapes,' history tapes'
       do t=1,ntapes
          write(6,*)
          if (hist_nhtfrq(t) == 0) then
             write(6,*)'History tape ',t,' write frequency is MONTHLY'
          else
             write(6,*)'History tape ',t,' write frequency = ',hist_nhtfrq(t)
          endif
          if (hist_dov2xy(t)) then
             write(6,*)'All fields on history tape ',t,' are grid averaged'
          else
             write(6,*)'All fields on history tape ',t,' are not grid averaged'
          end if
          write(6,*)'Number of time samples on history tape ',t,' is ',hist_mfilt(t)
          write(6,*)'Output precision on history tape ',t,'=',hist_ndens(t)
          write(6,*)
       end do
    end if

    ! Set flag indicating h-tape contents are now defined (needed by masterlist_addfld)

    htapes_defined = .true.

  end subroutine htapes_fieldlist

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: htape_addfld
!
! !INTERFACE:
  subroutine htape_addfld (t, f, avgflag)
!
! !DESCRIPTION:
! Add a field to the active list for a history tape. Copy the data from
! the master field list to the active list for the tape.
!
! !USES:
    use decompMod , only : get_proc_bounds, get_proc_global
    use clmtype   , only : nameg, namel, namec, namep, lndrof, ocnrof
    use decompMod , only : get_proc_bounds, get_proc_global
#if (defined RTM)
    use RunoffMod , only : get_proc_rof_bounds, get_proc_rof_global
#endif
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t                 ! history tape index
    integer, intent(in) :: f                 ! field index from master field list
    character(len=1), intent(in) :: avgflag  ! time averaging flag
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: n                    ! field index on defined tape
    character(len=8) :: type1d      ! clm pointer 1d type
    character(len=8) :: type1d_out  ! history buffer 1d type
    integer :: begp, endp           ! per-proc beginning and ending pft indices
    integer :: begc, endc           ! per-proc beginning and ending column indices
    integer :: begl, endl           ! per-proc beginning and ending landunit indices
    integer :: begg, endg           ! per-proc gridcell ending gridcell indices
    integer :: numg                 ! total number of gridcells across all processors
    integer :: numl                 ! total number of landunits across all processors
    integer :: numc                 ! total number of columns across all processors
    integer :: nump                 ! total number of pfts across all processors
    integer :: num2d                ! size of second dimension (e.g. .number of vertical levels)
    integer :: beg1d_out,end1d_out  ! history output per-proc 1d beginning and ending indices
    integer :: num1d_out            ! history output 1d size
#if (defined RTM)
    integer :: beg_roflnd           ! per-proc beginning land runoff index
    integer :: end_roflnd           ! per-proc ending land runoff index
    integer :: beg_rofocn           ! per-proc beginning ocean runoff index
    integer :: end_rofocn           ! per-proc ending ocean runoff index
    integer :: num_roflnd           ! total number of land runoff points across all procs
    integer :: num_rofocn           ! total number of ocean runoff points across all procs
#endif
!-----------------------------------------------------------------------

    ! Ensure that it is not to late to add a field to the history tape

    if (htapes_defined) then
       write(6,*)'HTAPE_ADDFLD error : attempt to add field ', &
            masterlist(f)%field%name, ' after history files are set'
       call endrun()
    end if

    tape(t)%nflds = tape(t)%nflds + 1
    n = tape(t)%nflds

    ! Copy field information

    tape(t)%hlist(n)%field = masterlist(f)%field

    ! Determine bounds

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)
#if (defined RTM)
    call get_proc_rof_bounds(beg_roflnd, end_roflnd, beg_rofocn, end_rofocn)
    call get_proc_rof_global(num_roflnd, num_rofocn)
#endif

    ! Modify type1d_out if necessary

    if (hist_dov2xy(t)) then

       ! If xy output averaging is requested, set output 1d type to gridcell

       type1d = tape(t)%hlist(n)%field%type1d
       if (type1d == namel .or. type1d == namec .or. type1d == namep) then
          tape(t)%hlist(n)%field%type1d_out = nameg
       end if

    else if (t > 1 .and. hist_type1d_pertape(t) /= ' ') then

       ! Set output 1d type  based on namelist setting of  hist_type1d_pertape
       ! Only applies to tapes other than primary and when xy output is not required

       select case (trim(hist_type1d_pertape(t)))
       case('GRID')
          tape(t)%hlist(n)%field%type1d_out = nameg
       case('LAND')
          tape(t)%hlist(n)%field%type1d_out = namel
       case('COLS')
          tape(t)%hlist(n)%field%type1d_out = namec
       case ('PFTS')
          tape(t)%hlist(n)%field%type1d_out = namep
       case default
          write(6,*)'HTAPE_ADDFLD error: unknown input hist_type1d_pertape= ', hist_type1d_pertape(t)
          call endrun()
       end select

    endif

    ! Determine output 1d dimensions

    type1d_out = tape(t)%hlist(n)%field%type1d_out
    if      (type1d_out == nameg) then
       beg1d_out = begg
       end1d_out = endg
       num1d_out = numg
    else if (type1d_out == namel) then
       beg1d_out = begl
       end1d_out = endl
       num1d_out = numl
    else if (type1d_out == namec) then
       beg1d_out = begc
       end1d_out = endc
       num1d_out = numc
    else if (type1d_out == namep) then
       beg1d_out = begp
       end1d_out = endp
       num1d_out = nump
#ifdef RTM
    else if (type1d_out == lndrof) then
       beg1d_out = beg_roflnd
       end1d_out = end_roflnd
       num1d_out = num_roflnd
    else if (type1d_out == ocnrof) then
       beg1d_out = beg_rofocn
       end1d_out = end_rofocn
       num1d_out = num_rofocn
#endif
    else
       write(6,*)'HTAPE_ADDFLD error: incorrect value of type1d_out= ',type1d_out
       call endrun()
    end if

    tape(t)%hlist(n)%field%beg1d_out = beg1d_out
    tape(t)%hlist(n)%field%end1d_out = end1d_out
    tape(t)%hlist(n)%field%num1d_out = num1d_out
    
    ! Alloccate and initialize history buffer and related info

    num2d = tape(t)%hlist(n)%field%num2d
    allocate (tape(t)%hlist(n)%hbuf(beg1d_out:end1d_out,num2d))
    allocate (tape(t)%hlist(n)%nacs(beg1d_out:end1d_out,num2d))
    tape(t)%hlist(n)%hbuf(:,:) = 0._r8
    tape(t)%hlist(n)%nacs(:,:) = 0

    ! Set time averaging flag based on masterlist setting or
    ! override the default averaging flag with namelist setting

    select case (avgflag)
    case (' ')
       tape(t)%hlist(n)%avgflag = masterlist(f)%avgflag(t)
    case ('A','I','X','M')
       tape(t)%hlist(n)%avgflag = avgflag
    case default
       write(6,*)'HTAPE_ADDFLD error: unknown avgflag=', avgflag
       call endrun()
    end select

  end subroutine htape_addfld

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: update_hbuf
!
! !INTERFACE:
  subroutine update_hbuf()
!
! !DESCRIPTION:
! Accumulate (or take min, max, etc. as appropriate) input field
! into its history buffer for appropriate tapes.
!
! !USES:
    use decompMod , only : get_proc_bounds, get_proc_global
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: t                   ! tape index
    integer :: f                   ! field index
    integer :: begp, endp          ! per-proc beginning and ending pft indices
    integer :: begc, endc          ! per-proc beginning and ending column indices
    integer :: begl, endl          ! per-proc beginning and ending landunit indices
    integer :: begg, endg          ! per-proc gridcell ending gridcell indices
    integer :: num2d               ! size of second dimension (e.g. number of vertical levels)
!-----------------------------------------------------------------------

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    do t = 1,ntapes
!$OMP PARALLEL DO PRIVATE (f, num2d)
#if !defined (USE_OMP)
!CSD$ PARALLEL DO PRIVATE (f, num2d)
#endif
       do f = 1,tape(t)%nflds
          num2d = tape(t)%hlist(f)%field%num2d
          if ( num2d == 1) then
             call update_hbuf_field_1d (t, f, begp, endp, begc, endc, begl, endl, begg, endg)
          else
             call update_hbuf_field_2d (t, f, begp, endp, begc, endc, begl, endl, begg, endg, num2d)
          end if
       end do
#if !defined (USE_OMP)
!CSD$ END PARALLEL DO
#endif
!$OMP END PARALLEL DO
    end do

  end subroutine update_hbuf

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: update_hbuf_field_1d
!
! !INTERFACE:
  subroutine update_hbuf_field_1d (t, f, begp, endp, begc, endc, begl, endl, begg, endg)
!
! !DESCRIPTION:
! Accumulate (or take min, max, etc. as appropriate) input field
! into its history buffer for appropriate tapes.
!
! !USES:
    use clmtype
    use subgridAveMod, only : p2g, c2g, l2g
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t            ! tape index
    integer, intent(in) :: f            ! field index
    integer, intent(in) :: begp, endp   ! per-proc beginning and ending pft indices
    integer, intent(in) :: begc, endc   ! per-proc beginning and ending column indices
    integer, intent(in) :: begl, endl   ! per-proc beginning and ending landunit indices
    integer, intent(in) :: begg, endg   ! per-proc gridcell ending gridcell indices
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: hpindex                 ! history pointer index
    integer  :: k                       ! gridcell, landunit, column or pft index
    integer  :: beg1d,end1d             ! beginning and ending indices
    logical  :: checkwt                 ! true => check weight of pft relative to gridcell
    logical  :: valid                   ! true => history operation is valid
    logical  :: map2gcell               ! true => map clm pointer field to gridcell
    character(len=8)  :: type1d         ! 1d clm pointerr type   ["gridcell","landunit","column","pft","roflnd","rofocn"]
    character(len=8)  :: type1d_out     ! 1d history buffer type ["gridcell","landunit","column","pft","roflnd","rofocn"]
    character(len=1)  :: avgflag        ! time averaging flag
    character(len=8)  :: p2c_scale_type ! scale type for subgrid averaging of pfts to column
    character(len=8)  :: c2l_scale_type ! scale type for subgrid averaging of columns to landunits
    character(len=8)  :: l2g_scale_type ! scale type for subgrid averaging of landunits to gridcells
    real(r8), pointer :: hbuf(:,:)      ! history buffer
    integer , pointer :: nacs(:,:)      ! accumulation counter
    real(r8), pointer :: pwtgcell(:)    ! weight of pft relative to corresponding gridcell
    real(r8), pointer :: field(:)       ! clm 1d pointer field
    real(r8) :: field_gcell(begg:endg)  ! gricell level field (used if mapping to gridcell is done)
    integer j
!-----------------------------------------------------------------------

    avgflag        =  tape(t)%hlist(f)%avgflag
    nacs           => tape(t)%hlist(f)%nacs
    hbuf           => tape(t)%hlist(f)%hbuf
    beg1d          =  tape(t)%hlist(f)%field%beg1d
    end1d          =  tape(t)%hlist(f)%field%end1d
    type1d         =  tape(t)%hlist(f)%field%type1d
    type1d_out     =  tape(t)%hlist(f)%field%type1d_out
    p2c_scale_type =  tape(t)%hlist(f)%field%p2c_scale_type
    c2l_scale_type =  tape(t)%hlist(f)%field%c2l_scale_type
    l2g_scale_type =  tape(t)%hlist(f)%field%l2g_scale_type
    hpindex        =  tape(t)%hlist(f)%field%hpindex
    field          => clmptr_rs(hpindex)%ptr

    ! set variables to check weights when allocate all pfts

    map2gcell = .false.
    if (type1d_out == nameg) then
       if (type1d == namep) then
          call p2g(begp, endp, begc, endc, begl, endl, begg, endg, field, field_gcell, &
               p2c_scale_type, c2l_scale_type, l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namec) then
          call c2g(begc, endc, begl, endl, begg, endg, field, field_gcell, &
               c2l_scale_type, l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namel) then
          call l2g(begl, endl, begg, endg, field, field_gcell, &
               l2g_scale_type)
          map2gcell = .true.
       end if
    end if

    if (map2gcell) then  ! Map to gridcell

       ! note that in this case beg1d = begg and end1d=endg
       select case (avgflag)
       case ('I') ! Instantaneous
!dir$ concurrent
!cdir nodep
          do k = begg,endg
             if (field_gcell(k) /= spval) then
                hbuf(k,1) = field_gcell(k)
             else
                hbuf(k,1) = spval
             end if
             nacs(k,1) = 1
          end do
       case ('A') ! Time average
!dir$ concurrent
!cdir nodep
          do k = begg,endg
             if (field_gcell(k) /= spval) then
                if (nacs(k,1) == 0) hbuf(k,1) = 0._r8
                hbuf(k,1) = hbuf(k,1) + field_gcell(k)
                nacs(k,1) = nacs(k,1) + 1
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
          end do
       case ('X') ! Maximum over time
!dir$ concurrent
!cdir nodep
          do k = begg,endg
             if (field_gcell(k) /= spval) then
                if (nacs(k,1) == 0) hbuf(k,1) = -1.e50_r8
                hbuf(k,1) = max( hbuf(k,1), field_gcell(k) )
             else
                hbuf(k,1) = spval
             endif
             nacs(k,1) = 1
          end do
       case ('M') ! Minimum over time
!dir$ concurrent
!cdir nodep
          do k = begg,endg
             if (field_gcell(k) /= spval) then
                if (nacs(k,1) == 0) hbuf(k,1) = +1.e50_r8
                hbuf(k,1) = min( hbuf(k,1), field_gcell(k) )
             else
                hbuf(k,1) = spval
             endif
             nacs(k,1) = 1
          end do
       case default
          write(6,*)'UPDATE_HBUF_FIELD_1d error: invalid time averaging flag ', avgflag
          call endrun()
       end select

    else  ! Do not map to gridcell

       pwtgcell => clm3%g%l%c%p%wtgcell
       checkwt = .false.
       if (type1d == 'pft') checkwt = .true.

       select case (avgflag)
       case ('I') ! Instantaneous
!dir$ concurrent
!cdir nodep
          do k = beg1d,end1d
             valid = .true.
             if (checkwt) then
                if (pwtgcell(k) == 0._r8) valid = .false.
             end if
             if (valid) then
                if (field(k) /= spval) then
                   hbuf(k,1) = field(k)
                else
                   hbuf(k,1) = spval
                end if
             else
                hbuf(k,1) = spval
             end if
             nacs(k,1) = 1
          end do
       case ('A') ! Time average
!dir$ concurrent
!cdir nodep
          do k = beg1d,end1d
             valid = .true.
             if (checkwt) then
                if (pwtgcell(k) == 0._r8) valid = .false.
             end if
             if (valid) then
                if (field(k) /= spval) then
                   if (nacs(k,1) == 0) hbuf(k,1) = 0._r8
                   hbuf(k,1) = hbuf(k,1) + field(k)
                   nacs(k,1) = nacs(k,1) + 1
                else
                   if (nacs(k,1) == 0) hbuf(k,1) = spval
                end if
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
          end do
       case ('X') ! Maximum over time
!dir$ concurrent
!cdir nodep
          do k = beg1d,end1d
             valid = .true.
             if (checkwt) then
                if (pwtgcell(k) == 0._r8) valid = .false.
             end if
             if (valid) then
                if (field(k) /= spval) then
                   if (nacs(k,1) == 0) hbuf(k,1) = -1.e50_r8
                   hbuf(k,1) = max( hbuf(k,1), field(k) )
                else
                   if (nacs(k,1) == 0) hbuf(k,1) = spval
                end if
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
             nacs(k,1) = 1
          end do
       case ('M') ! Minimum over time
!dir$ concurrent
!cdir nodep
          do k = beg1d,end1d
             valid = .true.
             if (checkwt) then
                if (pwtgcell(k) == 0._r8) valid = .false.
             end if
             if (valid) then
                if (field(k) /= spval) then
                   if (nacs(k,1) == 0) hbuf(k,1) = +1.e50_r8
                   hbuf(k,1) = min( hbuf(k,1), field(k) )
                else
                   if (nacs(k,1) == 0) hbuf(k,1) = spval
                end if
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
             nacs(k,1) = 1
          end do
       case default
          write(6,*)'UPDATE_HBUF_FIELD_1d error: invalid time averaging flag ', avgflag
          call endrun()
       end select
    end if

  end subroutine update_hbuf_field_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: update_hbuf_field_2d
!
! !INTERFACE:
  subroutine update_hbuf_field_2d (t, f, begp, endp, begc, endc, begl, endl, begg, endg, num2d)
!
! !DESCRIPTION:
! Accumulate (or take min, max, etc. as appropriate) input field
! into its history buffer for appropriate tapes.
!
! !USES:
    use clmtype
    use subgridAveMod, only : p2g, c2g, l2g
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t            ! tape index
    integer, intent(in) :: f            ! field index
    integer, intent(in) :: begp, endp   ! per-proc beginning and ending pft indices
    integer, intent(in) :: begc, endc   ! per-proc beginning and ending column indices
    integer, intent(in) :: begl, endl   ! per-proc beginning and ending landunit indices
    integer, intent(in) :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer, intent(in) :: num2d        ! size of second dimension
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: hpindex                 ! history pointer index
    integer  :: k                       ! gridcell, landunit, column or pft index
    integer  :: j                       ! level index
    integer  :: beg1d,end1d             ! beginning and ending indices
    logical  :: checkwt                 ! true => check weight of pft relative to gridcell
    logical  :: valid                   ! true => history operation is valid
    logical  :: map2gcell               ! true => map clm pointer field to gridcell
    character(len=8)  :: type1d         ! 1d clm pointerr type   ["gridcell","landunit","column","pft","roflnd","rofocn"]
    character(len=8)  :: type1d_out     ! 1d history buffer type ["gridcell","landunit","column","pft","roflnd","rofocn"]
    character(len=1)  :: avgflag        ! time averaging flag
    character(len=8)  :: p2c_scale_type ! scale type for subgrid averaging of pfts to column
    character(len=8)  :: c2l_scale_type ! scale type for subgrid averaging of columns to landunits
    character(len=8)  :: l2g_scale_type ! scale type for subgrid averaging of landunits to gridcells
    real(r8), pointer :: hbuf(:,:)      ! history buffer
    integer , pointer :: nacs(:,:)      ! accumulation counter
    real(r8), pointer :: pwtgcell(:)    ! weight of pft relative to corresponding gridcell
    real(r8), pointer :: field(:,:)     ! clm 2d pointer field
    real(r8) :: field_gcell(begg:endg,num2d) ! gricell level field (used if mapping to gridcell is done)
!-----------------------------------------------------------------------

    avgflag        =  tape(t)%hlist(f)%avgflag
    nacs           => tape(t)%hlist(f)%nacs
    hbuf           => tape(t)%hlist(f)%hbuf
    beg1d          =  tape(t)%hlist(f)%field%beg1d
    end1d          =  tape(t)%hlist(f)%field%end1d
    type1d         =  tape(t)%hlist(f)%field%type1d
    type1d_out     =  tape(t)%hlist(f)%field%type1d_out
    p2c_scale_type =  tape(t)%hlist(f)%field%p2c_scale_type
    c2l_scale_type =  tape(t)%hlist(f)%field%c2l_scale_type
    l2g_scale_type =  tape(t)%hlist(f)%field%l2g_scale_type
    hpindex        =  tape(t)%hlist(f)%field%hpindex
    field          => clmptr_ra(hpindex)%ptr(:,1:num2d)

    ! set variables to check weights when allocate all pfts

    map2gcell = .false.
    if (type1d_out == nameg) then
       if (type1d == namep) then
          call p2g(begp, endp, begc, endc, begl, endl, begg, endg, num2d, field, field_gcell, &
               p2c_scale_type, c2l_scale_type, l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namec) then
          call c2g(begc, endc, begl, endl, begg, endg, num2d, field, field_gcell, &
               c2l_scale_type, l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namel) then
          call l2g(begl, endl, begg, endg, num2d, field, field_gcell, &
               l2g_scale_type)
          map2gcell = .true.
       end if
    end if

    if (map2gcell) then  ! Map to gridcell

       ! note that in this case beg1d = begg and end1d=endg
       select case (avgflag)
       case ('I') ! Instantaneous
          do j = 1,num2d
!dir$ concurrent
!cdir nodep
             do k = begg,endg
                if (field_gcell(k,j) /= spval) then
                   hbuf(k,j) = field_gcell(k,j)
                else
                   hbuf(k,j) = spval
                end if
                nacs(k,j) = 1
             end do
          end do
       case ('A') ! Time average
          do j = 1,num2d
!dir$ concurrent
!cdir nodep
             do k = begg,endg
                if (field_gcell(k,j) /= spval) then
                   if (nacs(k,j) == 0) hbuf(k,j) = 0._r8
                   hbuf(k,j) = hbuf(k,j) + field_gcell(k,j)
                   nacs(k,j) = nacs(k,j) + 1
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                endif
             end do
          end do
       case ('X') ! Maximum over time
          do j = 1,num2d
!dir$ concurrent
!cdir nodep
             do k = begg,endg
                if (field_gcell(k,j) /= spval) then
                   if (nacs(k,j) == 0) hbuf(k,j) = -1.e50_r8
                   hbuf(k,j) = max( hbuf(k,j), field_gcell(k,j) )
                else
                   hbuf(k,j) = spval
                endif
                nacs(k,j) = 1
             end do
          end do
       case ('M') ! Minimum over time
          do j = 1,num2d
!dir$ concurrent
!cdir nodep
             do k = begg,endg
                if (field_gcell(k,j) /= spval) then
                   if (nacs(k,j) == 0) hbuf(k,j) = +1.e50_r8
                   hbuf(k,j) = min( hbuf(k,j), field_gcell(k,j) )
                else
                   hbuf(k,j) = spval
                endif
                nacs(k,j) = 1
             end do
          end do
       case default
          write(6,*)'UPDATE_HBUF_FIELD_2d error: invalid time averaging flag ', avgflag
          call endrun()
       end select

    else  ! Do not map to gridcell

       ! Note that since field points to an array section the
       ! bounds are field(1:end1d-beg1d+1, num2d) - therefore
       ! need to do the shifting below

       pwtgcell => clm3%g%l%c%p%wtgcell
       checkwt = .false.
       if (type1d == 'pft') checkwt = .true.

       select case (avgflag)
       case ('I') ! Instantaneous
          do j = 1,num2d
!dir$ concurrent
!cdir nodep
             do k = beg1d,end1d
                valid = .true.
                if (checkwt) then
                   if (pwtgcell(k) == 0._r8) valid = .false.
                end if
                if (valid) then
                   if (field(k-beg1d+1,j) /= spval) then
                      hbuf(k,j) = field(k-beg1d+1,j)
                   else
                      hbuf(k,j) = spval
                   end if
                else
                   hbuf(k,j) = spval
                end if
                nacs(k,j) = 1
             end do
          end do
       case ('A') ! Time average
          do j = 1,num2d
!dir$ concurrent
!cdir nodep
             do k = beg1d,end1d
                valid = .true.
                if (checkwt) then
                   if (pwtgcell(k) == 0._r8) valid = .false.
                end if
                if (valid) then
                   if (field(k-beg1d+1,j) /= spval) then
                      if (nacs(k,j) == 0) hbuf(k,j) = 0._r8
                      hbuf(k,j) = hbuf(k,j) + field(k-beg1d+1,j)
                      nacs(k,j) = nacs(k,j) + 1
                   else
                      if (nacs(k,j) == 0) hbuf(k,j) = spval
                   end if
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
             end do
          end do
       case ('X') ! Maximum over time
          do j = 1,num2d
!dir$ concurrent
!cdir nodep
             do k = beg1d,end1d
                valid = .true.
                if (checkwt) then
                   if (pwtgcell(k) == 0._r8) valid = .false.
                end if
                if (valid) then
                   if (field(k-beg1d+1,j) /= spval) then
                      if (nacs(k,j) == 0) hbuf(k,j) = -1.e50_r8
                      hbuf(k,j) = max( hbuf(k,j), field(k-beg1d+1,j) )
                   else
                      if (nacs(k,j) == 0) hbuf(k,j) = spval
                   end if
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
                nacs(k,j) = 1
             end do
          end do
       case ('M') ! Minimum over time
          do j = 1,num2d
!dir$ concurrent
!cdir nodep
             do k = beg1d,end1d
                valid = .true.
                if (checkwt) then
                   if (pwtgcell(k) == 0._r8) valid = .false.
                end if
                if (valid) then
                   if (field(k-beg1d+1,j) /= spval) then
                      if (nacs(k,j) == 0) hbuf(k,j) = +1.e50_r8
                      hbuf(k,j) = min( hbuf(k,j), field(k-beg1d+1,j))
                   else
                      if (nacs(k,j) == 0) hbuf(k,j) = spval
                   end if
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
                nacs(k,j) = 1
             end do
          end do
       case default
          write(6,*)'UPDATE_HBUF_FIELD_2d error: invalid time averaging flag ', avgflag
          call endrun()
       end select
    end if

  end subroutine update_hbuf_field_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: hfields_normalize
!
! !INTERFACE:
  subroutine hfields_normalize (t)
!
! !DESCRIPTION:
! Normalize fields on a history file by the number of accumulations.
! Loop over fields on the tape.  Need averaging flag and number of
! accumulations to perform normalization.
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t       ! tape index
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: f                   ! field index
    integer :: k                   ! 1d index
    integer :: j                   ! 2d index
    logical :: aflag               ! averaging flag
    integer :: beg1d_out,end1d_out ! hbuf 1d beginning and ending indices
    integer :: num2d               ! hbuf size of second dimension (e.g. number of vertical levels)
    character(len=1)  :: avgflag   ! averaging flag
    real(r8), pointer :: hbuf(:,:) ! history buffer
    integer , pointer :: nacs(:,:) ! accumulation counter
!-----------------------------------------------------------------------
!dir$ inlinenever hfields_normalize

    ! Normalize by number of accumulations for time averaged case

    do f = 1,tape(t)%nflds
       avgflag   =  tape(t)%hlist(f)%avgflag
       beg1d_out =  tape(t)%hlist(f)%field%beg1d_out
       end1d_out =  tape(t)%hlist(f)%field%end1d_out
       num2d     =  tape(t)%hlist(f)%field%num2d
       nacs      => tape(t)%hlist(f)%nacs
       hbuf      => tape(t)%hlist(f)%hbuf

       if (avgflag == 'A') then
          aflag = .true.
       else
          aflag = .false.
       end if

       do j = 1, num2d
          do k = beg1d_out, end1d_out
             if (aflag .and. nacs(k,j) /= 0) then
                hbuf(k,j) = hbuf(k,j) / float(nacs(k,j))
             end if
          end do
       end do
    end do

  end subroutine hfields_normalize

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: hfields_zero
!
! !INTERFACE:
  subroutine hfields_zero (t)
!
! !DESCRIPTION:
! Zero out accumulation and history buffers for a given history tape.
! Loop through fields on the tape.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t     ! tape index
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: f                 ! field index
!-----------------------------------------------------------------------

    do f = 1,tape(t)%nflds
       tape(t)%hlist(f)%hbuf(:,:) = 0._r8
       tape(t)%hlist(f)%nacs(:,:) = 0
    end do

  end subroutine hfields_zero

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: htape_create
!
! !INTERFACE:
  subroutine htape_create (t)
!
! !DESCRIPTION:
! Define contents of history file t. Issue the required netcdf
! wrapper calls to define the history file contents.
!
! !USES:
    use clmtype
    use ncdio
    use decompMod   , only : get_proc_global
    use clm_varpar  , only : lsmlon, lsmlat, nlevsoi
    use clm_varctl  , only : caseid, ctitle, frivinp_rtm, fsurdat, finidat, fpftcon
    use clm_varsur  , only : fullgrid, longxy, latixy, area, landfrac, landmask, numlon, lsmedge
    use clm_varcon  , only : zsoi, zlak
    use fileutils   , only : get_filename
    use time_manager, only : get_ref_date
#if (defined RTM)
    use RtmMod      , only : latixy_r, longxy_r
    use RunoffMod   , only : get_proc_rof_global
#endif	
#if (defined CASA)
    use CASAMod,    only : nlive, npools
#endif
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t   ! tape index
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: f                   ! field index
    integer :: p,c,l,n             ! indices
    integer :: ier                 ! error code
    integer :: yr,mon,day,nbsec    ! year,month,day,seconds components of a date
    integer :: hours,minutes,secs  ! hours,minutes,seconds of hh:mm:ss
    integer :: num2d               ! size of second dimension (e.g. number of vertical levels)
    integer :: dimid               ! dimension id temporary
    integer :: dim1id(1)           ! netCDF dimension id
    integer :: dim2id(2)           ! netCDF dimension id
    integer :: ndims               ! dimension counter
    integer :: omode               ! returned mode from netCDF call
    integer :: ncprec              ! output netCDF write precision
    integer :: ret                 ! netCDF error status
    integer :: nump                ! total number of pfts across all processors
    integer :: numc                ! total number of columns across all processors
    integer :: numl                ! total number of landunits across all processors
    integer :: numg                ! total number of gridcells across all processors
    real(r8):: lonvar(lsmlon)      ! only used for full grid
    real(r8):: latvar(lsmlat)      ! only used for full grid
#ifdef RTM
    real(r8):: lonvar_rof(rtmlon)  ! only used for full grid
    real(r8):: latvar_rof(rtmlat)  ! only used for full grid
    integer :: num_lndrof          ! total number of land runoff across all procs
    integer :: num_ocnrof          ! total number of ocean runoff across all procs
#endif
    character(len=  8) :: curdate  ! current date
    character(len=  8) :: curtime  ! current time
    character(len= 10) :: basedate ! base date (yyyymmdd)
    character(len=  8) :: basesec  ! base seconds
    character(len=256) :: name     ! name of attribute
    character(len=256) :: units    ! units of attribute
    character(len=256) :: str      ! global attribute string
    character(len=  1) :: avgflag  ! time averaging flag
    character(len=32) :: subname = 'htape_create'
!-----------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_global(numg, numl, numc, nump)
#ifdef RTM
    call get_proc_rof_global(num_lndrof, num_ocnrof)
#endif

    ! define output write precsion for tape

    ncprec = tape(t)%ncprec

    ! Create new netCDF file. File will be in define mode

    write(6,*)'Opening netcdf htape ', trim(locfnh(t))
    call check_ret(nf_create (trim(locfnh(t)), nf_clobber, nfid(t)), subname)
    call check_ret(nf_set_fill (nfid(t), nf_nofill, omode), subname)

    ! Create global attributes. Attributes are used to store information
    ! about the data set. Global attributes are information about the
    ! data set as a whole, as opposed to a single variable

    str = 'CF-1.0'
    call check_ret(&
         nf_put_att_text (nfid(t), NF_GLOBAL, 'conventions', len_trim(str), trim(str)), subname)

    call getdatetime(curdate, curtime)
    str = 'created on ' // curdate // ' ' // curtime
    call check_ret(&
         nf_put_att_text(nfid(t), NF_GLOBAL, 'history', len_trim(str), trim(str)), subname)

    call shr_sys_getenv('LOGNAME', str, ier)
    call check_ret(&
         nf_put_att_text (nfid(t), NF_GLOBAL, 'logname',len_trim(str), trim(str)), subname)

    call shr_sys_getenv('HOST', str, ier)
    call check_ret(&
         nf_put_att_text (nfid(t), NF_GLOBAL, 'host', len_trim(str), trim(str)), subname)

    str = 'Community Land Model: CLM2'
    call check_ret(&
         nf_put_att_text (nfid(t), NF_GLOBAL, 'source', len_trim(str), trim(str)), subname)

    str = '$Name$'
    call check_ret(&
         nf_put_att_text (nfid(t), NF_GLOBAL, 'version', len_trim(str), trim(str)), subname)
     
    str = '$Id$'
    call check_ret(&
         nf_put_att_text (nfid(t), NF_GLOBAL, 'revision_id', len_trim(str), trim(str)), subname)

    str = ctitle
    call check_ret(&
         nf_put_att_text (nfid(t), NF_GLOBAL, 'case_title', len_trim(str), trim(str)), subname)

    str = caseid
    call check_ret(&
         nf_put_att_text (nfid(t), NF_GLOBAL, 'case_id', len_trim(str), trim(str)), subname)

    str = get_filename(fsurdat)
    call check_ret(&
         nf_put_att_text(nfid(t), NF_GLOBAL, 'Surface_dataset', len_trim(str), trim(str)), subname)

    if (finidat == ' ') then
       str = 'arbitrary initialization'
    else
       str = get_filename(finidat)
    endif
    call check_ret(&
         nf_put_att_text(nfid(t), NF_GLOBAL, 'Initial_conditions_dataset', len_trim(str), trim(str)), subname)

    str = get_filename(fpftcon)
    call check_ret(&
         nf_put_att_text(nfid(t), NF_GLOBAL, 'PFT_physiological_constants_dataset', len_trim(str), trim(str)), subname)

    if (frivinp_rtm /= ' ') then
       str = get_filename(frivinp_rtm)
       call check_ret(&
            nf_put_att_text(nfid(t), NF_GLOBAL, 'RTM_input_datset', len_trim(str), trim(str)), subname)
    endif

    ! Define dimensions.
    ! Time is an unlimited dimension. Character string is treated as an array of characters.

    call check_ret(nf_def_dim (nfid(t), 'gridcell', numg   , dimid), subname)
    call check_ret(nf_def_dim (nfid(t), 'landunit', numl   , dimid), subname)
    call check_ret(nf_def_dim (nfid(t), 'column'  , numc   , dimid), subname)
    call check_ret(nf_def_dim (nfid(t), 'pft'     , nump   , dimid), subname)
#ifdef RTM
    call check_ret(nf_def_dim (nfid(t), 'ocnrof', num_ocnrof, dimid), subname)
    call check_ret(nf_def_dim (nfid(t), 'lndrof', num_lndrof, dimid), subname)
#endif
    call check_ret(nf_def_dim (nfid(t), 'levsoi', nlevsoi, dimid), subname)
    call check_ret(nf_def_dim (nfid(t), 'levlak', nlevlak, dimid), subname)
    call check_ret(nf_def_dim (nfid(t), 'numrad', numrad , dimid), subname)
#if (defined CASA)
    call check_ret(nf_def_dim (nfid(t), 'nlive', nlive , dimid), subname)
    call check_ret(nf_def_dim (nfid(t), 'npools', npools , dimid), subname)
#endif
    do n = 1,num_subs
       call check_ret(nf_def_dim (nfid(t), subs_name(n), subs_dim(n), dimid), subname)
    end do
    call check_ret(nf_def_dim (nfid(t), 'lon'   , lsmlon, dimid), subname)
    call check_ret(nf_def_dim (nfid(t), 'lat'   , lsmlat, dimid), subname)
#ifdef RTM
    call check_ret(nf_def_dim (nfid(t), 'lonrof', rtmlon, dimid), subname)
    call check_ret(nf_def_dim (nfid(t), 'latrof', rtmlat, dimid), subname)
#endif
    call check_ret(nf_def_dim (nfid(t), 'time', nf_unlimited, time_dimid), subname)
    call check_ret(nf_def_dim (nfid(t), 'hist_interval', 2, hist_interval_dimid), subname)
    call check_ret(nf_def_dim (nfid(t), 'string_length', 8, strlen_dimid), subname)

    write(6,*)'HTAPE_CREATE: Successfully defined netcdf history file ',t

    ! Define coordinate variables (including time)

    if (fullgrid) then
       call ncd_defvar(varname='lon', xtype=ncprec, dim1name='lon', &
            long_name='coordinate longitude', units='degrees_east', ncid=nfid(t))
       call ncd_defvar(varname='lat', xtype=ncprec, dim1name='lat', &
            long_name='coordinate latitude', units='degrees_north', ncid=nfid(t))
    endif
#ifdef RTM
    call ncd_defvar(varname='lonrof', xtype=ncprec, dim1name='lonrof', &
         long_name='runoff coordinate longitude', units='degrees_east', ncid=nfid(t))
    call ncd_defvar(varname='latrof', xtype=ncprec, dim1name='latrof', &
         long_name='runoff coordinate latitude', units='degrees_north', ncid=nfid(t))
#endif

    call get_ref_date(yr, mon, day, nbsec)
    hours   = nbsec / 3600
    minutes = (nbsec - hours*3600) / 60
    secs    = (nbsec - hours*3600 - minutes*60)
    write(basedate,80) yr,mon,day
80  format(i4.4,'-',i2.2,'-',i2.2)
    write(basesec ,90) hours, minutes, secs
90  format(i2.2,':',i2.2,':',i2.2)
    dim1id(1) = time_dimid
    call check_ret(nf_def_var(nfid(t), 'time', ncprec, 1, dim1id, time_var_id(t)), subname)
    str = 'time'
    call check_ret(nf_put_att_text (nfid(t), time_var_id(t), 'logname', len_trim(str), trim(str)), subname)
    str = 'days since ' // basedate // " " // basesec
    call check_ret(nf_put_att_text (nfid(t), time_var_id(t), 'units', len_trim(str), trim(str)), subname)
    str = 'noleap'
    call check_ret(nf_put_att_text (nfid(t), time_var_id(t), 'calendar', len_trim(str), trim(str)), subname)

    call ncd_defvar(varname='levsoi', xtype=ncprec, dim1name='levsoi', &
         long_name='coordinate soil levels', units='m', ncid=nfid(t))

    call ncd_defvar(varname='levlak', xtype=ncprec, dim1name='levlak', &
         long_name='coordinate lake levels', units='m', ncid=nfid(t))

    ! Define time-invariant grid variables

    call ncd_defvar(varname='edgen', xtype=ncprec, &
         long_name='northern edge of surface grid', units='degrees_north', ncid=nfid(t))
    
    call ncd_defvar(varname='edgee', xtype=ncprec, &
         long_name='eastern edge of surface grid', units='degrees_east', ncid=nfid(t))
    
    call ncd_defvar(varname='edges', xtype=ncprec, &
         long_name='southern edge of surface grid', units='degrees_north', ncid=nfid(t))
    
    call ncd_defvar(varname='edgew', xtype=ncprec, &
         long_name='western edge of surface grid', units='degrees_east', ncid=nfid(t))

    if (fullgrid) then
       call ncd_defvar(varname='longxy', xtype=ncprec, dim1name='lon', dim2name='lat', &
            long_name='longitude', units='degrees_east',  ncid=nfid(t))
    else
       call ncd_defvar(varname='rlongxy', xtype=ncprec, dim1name='lon', dim2name='lat',&
            long_name='rlongitude', units='degrees_east', ncid=nfid(t))
    endif

    call ncd_defvar(varname='latixy', xtype=ncprec, dim1name='lon', dim2name='lat',&
         long_name='latitude', units='degrees_north', ncid=nfid(t))

    call ncd_defvar(varname='area', xtype=ncprec, dim1name='lon', dim2name='lat',&
         long_name='grid cell areas', units='km^2', ncid=nfid(t))

    call ncd_defvar(varname='numlon', xtype=nf_int, dim1name='lat', &
         long_name='number of longitudes at each latitude', ncid=nfid(t))

    call ncd_defvar(varname='landfrac', xtype=ncprec, dim1name='lon', dim2name='lat', &
         long_name='land fraction', ncid=nfid(t))

    call ncd_defvar(varname='landmask', xtype=nf_int, dim1name='lon', dim2name='lat', &
         long_name='land/ocean mask (0.=ocean and 1.=land)', ncid=nfid(t))

    ! Define time information

    dim1id(1) = time_dimid
    call check_ret(nf_def_var(nfid(t) , 'mcdate', nf_int, 1, dim1id  , mcdate_id(t)), subname)
    str = 'current date (YYYYMMDD)'
    call check_ret(nf_put_att_text (nfid(t), mcdate_id(t), 'logname', len_trim(str), trim(str)), subname)

    call check_ret(nf_def_var(nfid(t) , 'mcsec' , nf_int, 1, dim1id , mcsec_id(t)), subname)
    str = 'current seconds of current date'
    call check_ret(nf_put_att_text (nfid(t), mcsec_id(t), 'logname', len_trim(str), trim(str)), subname)
    str = 's'
    call check_ret(nf_put_att_text (nfid(t), mcsec_id(t), 'units', len_trim(str), trim(str)), subname)

    call check_ret(nf_def_var(nfid(t) , 'mdcur' , nf_int, 1, dim1id , mdcur_id(t)), subname)
    str = 'current day (from base day)'
    call check_ret(nf_put_att_text (nfid(t), mdcur_id(t), 'logname', len_trim(str), trim(str)), subname)

    call check_ret(nf_def_var(nfid(t) , 'mscur' , nf_int, 1, dim1id , mscur_id(t)), subname)
    str = 'current seconds of current day'
    call check_ret(nf_put_att_text (nfid(t), mscur_id(t), 'logname', len_trim(str), trim(str)), subname)

    call check_ret(nf_def_var(nfid(t) , 'nstep' , nf_int, 1, dim1id , nstep_id(t)), subname)
    str = 'time step'
    call check_ret(nf_put_att_text (nfid(t), nstep_id(t), 'logname', len_trim(str), trim(str)), subname)

    dim2id(1) = hist_interval_dimid;  dim2id(2) = time_dimid
    call check_ret(nf_def_var(nfid(t), 'time_bounds', nf_double, 2, dim2id, time_bounds_id(t)), subname)
    str = 'history time interval endpoints'
    call check_ret(nf_put_att_text (nfid(t), time_bounds_id(t), 'logname', len_trim(str), trim(str)), subname)

    dim2id(1) = strlen_dimid;  dim2id(2) = time_dimid
    call check_ret(nf_def_var(nfid(t), 'date_written', nf_char, 2, dim2id, date_written_id(t)), subname)
    call check_ret(nf_def_var(nfid(t), 'time_written', nf_char, 2, dim2id, time_written_id(t)), subname)

    ! Define time constant variables

    if (t == 1) call htape_timeconst(t, mode='define')

    ! Define time variant variables

    call hfields_write(t, mode='define')

    ! Exit define model

    call check_ret(nf_enddef(nfid(t)), subname)

    ! Write out variables

    if (fullgrid) then
       lonvar(:) = longxy(1:lsmlon,1)
       latvar(:) = latixy(1,1:lsmlat)
       call ncd_ioglobal(varname='lon', data=lonvar, ncid=nfid(t), flag='write')
       call ncd_ioglobal(varname='lat', data=latvar, ncid=nfid(t), flag='write')
    endif
#ifdef RTM
    lonvar_rof = longxy_r(1:rtmlon,1)
    latvar_rof = latixy_r(1,1:rtmlat)
    call ncd_ioglobal(varname='lonrof', data=lonvar_rof, ncid=nfid(t), flag='write')
    call ncd_ioglobal(varname='latrof', data=latvar_rof, ncid=nfid(t), flag='write')
#endif
    call ncd_ioglobal(varname='levsoi', data=zsoi, ncid=nfid(t), flag='write')
    call ncd_ioglobal(varname='levlak', data=zlak, ncid=nfid(t), flag='write')

    call ncd_ioglobal(varname='edgen', data=lsmedge(1), ncid=nfid(t), flag='write')
    call ncd_ioglobal(varname='edgee', data=lsmedge(2), ncid=nfid(t), flag='write')
    call ncd_ioglobal(varname='edges', data=lsmedge(3), ncid=nfid(t), flag='write')
    call ncd_ioglobal(varname='edgew', data=lsmedge(4), ncid=nfid(t), flag='write')

    if (fullgrid) then
       call ncd_ioglobal(varname='longxy' , data=longxy, ncid=nfid(t), flag='write')
    else
       call ncd_ioglobal(varname='rlongxy', data=longxy, ncid=nfid(t), flag='write')
    endif
    call ncd_ioglobal(varname='latixy'  , data=latixy  , ncid=nfid(t), flag='write')
    call ncd_ioglobal(varname='area'    , data=area    , ncid=nfid(t), flag='write')
    call ncd_ioglobal(varname='landfrac', data=landfrac, ncid=nfid(t), flag='write')
    call ncd_ioglobal(varname='landmask', data=landmask, ncid=nfid(t), flag='write')
    call ncd_ioglobal(varname='numlon'  , data=numlon  , ncid=nfid(t), flag='write')

  end subroutine htape_create

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: htape_timeconst
!
! !INTERFACE:
  subroutine htape_timeconst(t, mode)
!
! !DESCRIPTION:
! Write time constant values to primary history tape.
! Issue the required netcdf wrapper calls to define the history file
! contents.
!
! !USES:
    use clmtype
    use subgridAveMod, only : c2g
    use decompMod    , only : get_proc_bounds, get_proc_global
    use clm_varpar   , only : lsmlon, lsmlat, nlevsoi
    use ncdio
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t              ! tape index
    character(len=*), intent(in) :: mode  ! 'define' or 'write'
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: c,l,lev,ifld,vid           ! indices
    integer :: ier                        ! error status
    integer :: begp, endp                 ! per-proc beginning and ending pft indices
    integer :: begc, endc                 ! per-proc beginning and ending column indices
    integer :: begl, endl                 ! per-proc beginning and ending landunit indices
    integer :: begg, endg                 ! per-proc gridcell ending gridcell indices
    integer :: numg                       ! total number of gridcells across all processors
    integer :: numl                       ! total number of landunits across all processors
    integer :: numc                       ! total number of columns across all processors
    integer :: nump                       ! total number of pfts across all processors
    character(len=max_chars) :: long_name ! variable long name
    character(len=max_namlen):: varname   ! variable name
    character(len=max_namlen):: units     ! variable units
    real(r8), pointer :: histi(:,:)       ! temporary
    real(r8), pointer :: histo(:,:)       ! temporary
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
!-----------------------------------------------------------------------

    if (mode == 'define') then

       if (masterproc) then
          do ifld = 1,5
             if (ifld == 1) then
                varname='ZSOI'
                long_name='soil depth'; units = 'm'
             else if (ifld == 2) then
                varname='DZSOI'
                long_name='soil thickness'; units = 'm'
             else if (ifld == 3) then
                varname='WATSAT'
                long_name='saturated soil water content (porosity)';  units = 'mm3/mm3'
             else if (ifld == 4) then
                varname='SUCSAT'
                long_name='saturated soil matric potential'; units = 'mm'
             else if (ifld == 5) then
                varname = 'BSW'
                long_name='slope of soil water retention curve'; units = 'unitless'
             end if
             if (tape(t)%dov2xy) then
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name='lon', dim2name='lat', dim3name='levsoi', &
                     long_name=long_name, units=units, missing_value=spval, fill_value=spval)
             else
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name='column', dim2name='levsoi', &
                     long_name=long_name, units=units, missing_value=spval, fill_value=spval)
             end if
          end do
       end if

    else if (mode == 'write') then

       ! Set pointers into derived type and get necessary bounds

       lptr => clm3%g%l
       cptr => clm3%g%l%c

       call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
       call get_proc_global(numg, numl, numc, nump)

       allocate(histi(begc:endc,nlevsoi), stat=ier)
       if (ier /= 0) then
          write(6,*)'HTAPE_TIMECONST error: allocation error for histi'; call endrun()
       end if

       ! Write time constant fields

       if (tape(t)%dov2xy) then
          allocate(histo(begg:endg,nlevsoi), stat=ier)
          if (ier /= 0) then
             write(6,*) "HTAPE_TIMECONST error: allocation error for histo"; call endrun()
          end if
       end if

       do ifld = 1,5
          if (ifld == 1) then
             varname='ZSOI'
          else if (ifld == 2) then
             varname='DZSOI'
          else if (ifld == 3) then
             varname='WATSAT'
          else if (ifld == 4) then
             varname='SUCSAT'
          else if (ifld == 5) then
             varname = 'BSW'
          end if
          histi(:,:) = spval
          do lev = 1,nlevsoi
!dir$ concurrent
!cdir nodep
             do c = begc, endc
                l = cptr%landunit(c)
                if (.not. lptr%lakpoi(l)) then
                   if (ifld ==1) histi(c,lev) = cptr%cps%z(c,lev)
                   if (ifld ==2) histi(c,lev) = cptr%cps%dz(c,lev)
                   if (ifld ==3) histi(c,lev) = cptr%cps%watsat(c,lev)
                   if (ifld ==4) histi(c,lev) = cptr%cps%sucsat(c,lev)
                   if (ifld ==5) histi(c,lev) = cptr%cps%bsw(c,lev)
                end if
             end do
          end do
          if (tape(t)%dov2xy) then
             histo(:,:) = spval
             call c2g(begc, endc, begl, endl, begg, endg, nlevsoi, histi, histo, &
                  c2l_scale_type='unity', l2g_scale_type='unity')

             call ncd_iolocal(varname=varname, dim1name='gridcell', dim2name='levsoi', data=histo, &
                  ncid=nfid(t), nlonxy=lsmlon, nlatxy=lsmlat, flag='write')
          else
             call ncd_iolocal(varname=varname, dim1name='column'  , dim2name='levsoi', data=histi, &
                  ncid=nfid(t), flag='write')
          end if
       end do

       if (tape(t)%dov2xy) deallocate(histo)
       deallocate(histi)

    end if

  end subroutine htape_timeconst

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: hfields_write
!
! !INTERFACE:
  subroutine hfields_write(t, mode)
!
! !DESCRIPTION:
! Write history tape. If SPMD, first gather the data to the master
! processor. Issue the netcdf call to write the variable.
!
! !USES:
    use ncdio
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t                ! tape index
    character(len=*), intent(in) :: mode    ! 'define' or 'write'
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: f                         ! field index
    integer :: k                         ! 1d index
    integer :: c,l,p                     ! indices
    integer :: nlonxy, nlatxy            ! xy grid dimensions
    integer :: beg1d_out                 ! on-node 1d hbuf pointer start index
    integer :: end1d_out                 ! on-node 1d hbuf pointer end index
    integer :: num1d_out                 ! size of hbuf first dimension (overall all nodes)
    integer :: num2d                     ! hbuf second dimension size
    integer :: nt                        ! time index
    integer :: ier                       ! error status
    character(len=1)         :: avgflag  ! time averaging flag
    character(len=max_chars) :: long_name! long name
    character(len=max_chars) :: units    ! units
    character(len=max_namlen):: varname  ! variable name
    character(len=32) :: avgstr          ! time averaging type
    character(len=8)  :: type1d_out      ! history output 1d type
    character(len=8)  :: type2d          ! history output 2d type
    character(len=32) :: dim1name        ! temporary
    character(len=32) :: dim2name        ! temporary
    real(r8), pointer :: histo(:,:)      ! temporary
    real(r8), pointer :: hist1do(:)      ! temporary
    character(len=32) :: subname = 'hfields_write'
!-----------------------------------------------------------------------

    ! Write/define 1d topological info

    if (.not. tape(t)%dov2xy) then
       if (mode == 'define') then
          call hfields_1dinfo(t, mode='define')
       else if (mode == 'write') then
          call hfields_1dinfo(t, mode='write')
       end if
    end if

    ! Define time-dependent variables create variables and attributes for field list

    do f = 1,tape(t)%nflds

       ! Set history field variables

       varname    = tape(t)%hlist(f)%field%name
       long_name  = tape(t)%hlist(f)%field%long_name
       units      = tape(t)%hlist(f)%field%units
       avgflag    = tape(t)%hlist(f)%avgflag
       type1d_out = tape(t)%hlist(f)%field%type1d_out
       beg1d_out  = tape(t)%hlist(f)%field%beg1d_out
       end1d_out  = tape(t)%hlist(f)%field%end1d_out
       num1d_out  = tape(t)%hlist(f)%field%num1d_out
       type2d     = tape(t)%hlist(f)%field%type2d
       num2d      = tape(t)%hlist(f)%field%num2d
       nlonxy     = tape(t)%hlist(f)%field%nlonxy
       nlatxy     = tape(t)%hlist(f)%field%nlatxy
       nt         = tape(t)%ntimes

       if (mode == 'define') then

          select case (avgflag)
          case ('A')
             avgstr = 'mean'
          case ('I')
             avgstr = 'instantaneous'
          case ('X')
             avgstr = 'maximum'
          case ('M')
             avgstr = 'minimum'
          case default
             write(6,*)subname,' error: unknown time averaging flag (avgflag)=',avgflag; call endrun()
          end select

          if (type1d_out == 'lndrof' .or. type1d_out == 'ocnrof') then

             dim1name = 'lonrof'; dim2name = 'latrof'
             call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                  dim1name=dim1name, dim2name=dim2name, dim3name='time', &
                  long_name=long_name, units=units, cell_method=avgstr, missing_value=spval, fill_value=spval)


          else if (tape(t)%dov2xy) then

             dim1name = 'lon'; dim2name = 'lat'
             if (num2d == 1) then
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name=dim1name, dim2name=dim2name, dim3name='time', &
                     long_name=long_name, units=units, cell_method=avgstr, missing_value=spval, fill_value=spval)
             else
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name=dim1name, dim2name=dim2name, dim3name=type2d, dim4name='time', &
                     long_name=long_name, units=units, cell_method=avgstr, missing_value=spval, fill_value=spval)
             end if

          else

             dim1name = type1d_out
             if (num2d == 1) then
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name=dim1name, dim2name='time', &
                     long_name=long_name, units=units, cell_method=avgstr, missing_value=spval, fill_value=spval)
             else
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name=dim1name, dim2name=type2d, dim3name='time', &
                     long_name=long_name, units=units, cell_method=avgstr, missing_value=spval, fill_value=spval)
             end if

          end if

       else if (mode == 'write') then

          ! Determine output buffer

          histo => tape(t)%hlist(f)%hbuf

          ! Allocate dynamic memory

          if (num2d == 1) then
             allocate(hist1do(beg1d_out:end1d_out), stat=ier)
             if (ier /= 0) then
                write(6,*)'hfields_write allocation error'; call endrun()
             end if
             hist1do(beg1d_out:end1d_out) = histo(beg1d_out:end1d_out,1)
          end if

          ! Write history output.  Always output land and ocean runoff on xy grid.

          if (type1d_out == 'lndrof' .or. type1d_out == 'ocnrof') then

             call ncd_iolocal(flag='write', varname=varname, dim1name=type1d_out, &
                  data=hist1do, ncid=nfid(t), nlonxy=nlonxy, nlatxy=nlatxy, nt=nt)

          else if (tape(t)%dov2xy) then

             if (num2d == 1) then
                call ncd_iolocal(flag='write', varname=varname, dim1name=type1d_out, &
                     data=hist1do, ncid=nfid(t), nlonxy=nlonxy, nlatxy=nlatxy, nt=nt)
             else
                call ncd_iolocal(flag='write', varname=varname, dim1name=type1d_out, dim2name=type2d, &
                     data=histo, ncid=nfid(t), nlonxy=nlonxy, nlatxy=nlatxy, nt=nt)
             end if

          else

             if (num2d == 1) then
                call ncd_iolocal(flag='write', varname=varname, dim1name=type1d_out, &
                     data=hist1do, ncid=nfid(t), nt=nt)
             else
                call ncd_iolocal(flag='write', varname=varname, dim1name=type1d_out, dim2name=type2d, &
                     data=histo, ncid=nfid(t), nt=nt)
             end if

          end if

          ! Deallocate dynamic memory

          if (num2d == 1) then
             deallocate(hist1do)
          end if

       end if

    end do

  end subroutine hfields_write

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: hfields_1dinfo
!
! !INTERFACE:
  subroutine hfields_1dinfo(t, mode)
!
! !DESCRIPTION:
! Write/define 1d info for history tape.
!
! !USES:
    use ncdio
    use clmtype         , only : nameg, namel, namec, namep
    use decompMod       , only : get_proc_bounds, get_proc_global, map_dc2sn
    use initGridCellsMod, only : get_sn_land1d, get_sn_cols1d, get_sn_pfts1d
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t                ! tape index
    character(len=*), intent(in) :: mode    ! 'define' or 'write'
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: f                         ! field index
    integer :: k                         ! 1d index
    integer :: c,l,p                     ! indices
    integer :: begp, endp                ! per-proc beginning and ending pft indices
    integer :: begc, endc                ! per-proc beginning and ending column indices
    integer :: begl, endl                ! per-proc beginning and ending landunit indices
    integer :: begg, endg                ! per-proc gridcell ending gridcell indices
    integer :: numg                      ! total number of gridcells across all processors
    integer :: numl                      ! total number of landunits across all processors
    integer :: numc                      ! total number of columns across all processors
    integer :: nump                      ! total number of pfts across all processors
    integer :: num2d                     ! 2d size (e.g. number of vertical levels)
    integer :: ncid                      ! netcdf file id
    integer :: ier                       ! errir status
    integer , pointer :: iltemp(:)       ! temporary
    integer , pointer :: ictemp(:)       ! temporary
    integer , pointer :: iptemp(:)       ! temporary
    integer , pointer :: ictype(:)       ! temporary
    integer , pointer :: iptype(:)       ! temporary
    type(gridcell_type), pointer :: gptr ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr ! pointer to pft derived subtype
!-----------------------------------------------------------------------

    ncid = nfid(t)

    if (mode == 'define') then

       if (masterproc) then

          ! Define gridcell info

          call ncd_defvar(varname='grid1d_lon', xtype=nf_double, dim1name='gridcell', &
               long_name='gridcell longitude', units='degrees_east', ncid=ncid)

          call ncd_defvar(varname='grid1d_lat', xtype=nf_double,  dim1name='gridcell', &
               long_name='gridcell latitude', units='degrees_north', ncid=ncid)

          call ncd_defvar(varname='grid1d_ixy', xtype=nf_int, dim1name='gridcell', &
               long_name='2d longitude index of corresponding gridcell', ncid=ncid)

          call ncd_defvar(varname='grid1d_jxy', xtype=nf_int, dim1name='gridcell', &
               long_name='2d latitude index of corresponding gridcell', ncid=ncid)

          ! Define landunit info

          call ncd_defvar(varname='land1d_lon', xtype=nf_double, dim1name='landunit', &
               long_name='landunit longitude', units='degrees_east', ncid=ncid)

          call ncd_defvar(varname='land1d_lat', xtype=nf_double, dim1name='landunit', &
               long_name='landunit latitude', units='degrees_north', ncid=ncid)

          call ncd_defvar(varname='land1d_ixy', xtype=nf_int, dim1name='landunit', &
               long_name='2d longitude index of corresponding landunit', ncid=ncid)

          call ncd_defvar(varname='land1d_jxy', xtype=nf_int, dim1name='landunit', &
               long_name='2d latitude index of corresponding landunit', ncid=ncid)

          call ncd_defvar(varname='land1d_gi', xtype=nf_int, dim1name='landunit', &
               long_name='1d grid index of corresponding landunit', ncid=ncid)

          call ncd_defvar(varname='land1d_wtgcell', xtype=nf_double, dim1name='landunit', &
               long_name='landunit weight relative to corresponding gridcell', ncid=ncid)

          call ncd_defvar(varname='land1d_ityplunit', xtype=nf_int, dim1name='landunit', &
               long_name='landunit type (vegetated,urban,lake,wetland or glacier)', ncid=ncid)

          ! Define column info

          call ncd_defvar(varname='cols1d_lon', xtype=nf_double, dim1name='column', &
               long_name='column longitude', units='degrees_east', ncid=ncid)

          call ncd_defvar(varname='cols1d_lat', xtype=nf_double, dim1name='column', &
               long_name='column latitude', units='degrees_north', ncid=ncid)

          call ncd_defvar(varname='cols1d_ixy', xtype=nf_int, dim1name='column', &
               long_name='2d longitude index of corresponding column', ncid=ncid)

          call ncd_defvar(varname='cols1d_jxy', xtype=nf_int, dim1name='column', &
               long_name='2d latitude index of corresponding column', ncid=ncid)

          call ncd_defvar(varname='cols1d_gi', xtype=nf_int, dim1name='column', &
               long_name='1d grid index of corresponding column', ncid=ncid)

          call ncd_defvar(varname='cols1d_li', xtype=nf_int, dim1name='column', &
               long_name='1d landunit index of corresponding column', ncid=ncid)

          call ncd_defvar(varname='cols1d_wtgcell', xtype=nf_double, dim1name='column', &
               long_name='column weight relative to corresponding gridcell', ncid=ncid)

          call ncd_defvar(varname='cols1d_wtlunit', xtype=nf_double, dim1name='column', &
               long_name='column weight relative to corresponding landunit', ncid=ncid)

          call ncd_defvar(varname='cols1d_itype_lunit', xtype=nf_int, dim1name='column', &
               long_name='column landunit type (vegetated,urban,lake,wetland or glacier)', ncid=ncid)

          ! Define pft info

          call ncd_defvar(varname='pfts1d_lon', xtype=nf_double, dim1name='pft', &
               long_name='pft longitude', units='degrees_east', ncid=ncid)

          call ncd_defvar(varname='pfts1d_lat', xtype=nf_double, dim1name='pft', &
               long_name='pft latitude', units='degrees_north', ncid=ncid)

          call ncd_defvar(varname='pfts1d_ixy', xtype=nf_int, dim1name='pft', &
               long_name='2d longitude index of corresponding pft', ncid=ncid)

          call ncd_defvar(varname='pfts1d_jxy', xtype=nf_int, dim1name='pft', &
               long_name='2d latitude index of corresponding pft', ncid=ncid)

          call ncd_defvar(varname='pfts1d_gi', xtype=nf_int, dim1name='pft', &
               long_name='1d grid index of corresponding pft', ncid=ncid)

          call ncd_defvar(varname='pfts1d_li', xtype=nf_int, dim1name='pft', &
               long_name='1d landunit index of corresponding pft', ncid=ncid)

          call ncd_defvar(varname='pfts1d_ci', xtype=nf_int, dim1name='pft', &
               long_name='1d column index of corresponding pft', ncid=ncid)

          call ncd_defvar(varname='pfts1d_wtgcell', xtype=nf_double, dim1name='pft', &
               long_name='pft weight relative to corresponding gridcell', ncid=ncid)

          call ncd_defvar(varname='pfts1d_wtlunit', xtype=nf_double, dim1name='pft', &
               long_name='pft weight relative to corresponding landunit', ncid=ncid)

          call ncd_defvar(varname='pfts1d_wtcol', xtype=nf_double, dim1name='pft', &
               long_name='pft weight relative to corresponding column', ncid=ncid)

          call ncd_defvar(varname='pfts1d_itype_veg', xtype=nf_int, dim1name='pft', &
               long_name='pft vegetation type', ncid=ncid)

          call ncd_defvar(varname='pfts1d_itype_lunit', xtype=nf_int, dim1name='pft', &
               long_name='pft landunit type (vegetated,urban,lake,wetland or glacier)', ncid=ncid)

       end if

    else if (mode == 'write') then

       ! Set pointers into derived type

       gptr => clm3%g
       lptr => clm3%g%l
       cptr => clm3%g%l%c
       pptr => clm3%g%l%c%p

       ! Determine bounds

       call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
       call get_proc_global(numg, numl, numc, nump)

       if (masterproc) then
          allocate(iltemp(numl), ictemp(numc), iptemp(nump), ictype(numc), iptype(nump), stat=ier)
          if (ier /= 0) call endrun('hfields_1dinfo allocation error')
       end if

       ! Write gridcell info

       call ncd_iolocal(varname='grid1d_lon', data=gptr%londeg, dim1name='gridcell', ncid=ncid, flag='write')
       call ncd_iolocal(varname='grid1d_lat', data=gptr%latdeg, dim1name='gridcell', ncid=ncid, flag='write')
       call ncd_iolocal(varname='grid1d_ixy', data=gptr%ixy   , dim1name='gridcell', ncid=ncid, flag='write')
       call ncd_iolocal(varname='grid1d_jxy', data=gptr%jxy   , dim1name='gridcell', ncid=ncid, flag='write')

       ! Write landunit info

       call ncd_iolocal(varname='land1d_lon', data=lptr%londeg, dim1name='landunit', ncid=ncid, flag='write')
       call ncd_iolocal(varname='land1d_lat', data=lptr%latdeg, dim1name='landunit', ncid=ncid, flag='write')
       call ncd_iolocal(varname='land1d_ixy', data=lptr%ixy   , dim1name='landunit', ncid=ncid, flag='write')
       call ncd_iolocal(varname='land1d_jxy', data=lptr%jxy   , dim1name='landunit', ncid=ncid, flag='write')
       if (masterproc) then
          call get_sn_land1d(iltemp, type1d='gridcell',numl=numl)
          call ncd_ioglobal(varname='land1d_gi', data=iltemp, ncid=ncid, flag='write')
       end if
       call ncd_iolocal(varname='land1d_wtgcell'  , data=lptr%wtgcell, dim1name='landunit', ncid=ncid, flag='write')
       call ncd_iolocal(varname='land1d_ityplunit', data=lptr%itype  , dim1name='landunit', ncid=ncid, flag='write')

       ! Write column info

       call ncd_iolocal(varname='cols1d_lon', data=cptr%londeg, dim1name='column', ncid=ncid, flag='write')
       call ncd_iolocal(varname='cols1d_lat', data=cptr%latdeg, dim1name='column', ncid=ncid, flag='write')
       call ncd_iolocal(varname='cols1d_ixy', data=cptr%ixy   , dim1name='column', ncid=ncid, flag='write')
       call ncd_iolocal(varname='cols1d_jxy', data=cptr%jxy   , dim1name='column', ncid=ncid, flag='write')
       if (masterproc) then
          call get_sn_cols1d(ictemp, type1d='gridcell',numc=numc)
          call ncd_ioglobal(varname='cols1d_gi', data=ictemp, ncid=ncid, flag='write')
          call get_sn_cols1d(ictemp, type1d='landunit',numc=numc)
          call ncd_ioglobal(varname='cols1d_li', data=ictemp, ncid=ncid, flag='write')
       end if
       call ncd_iolocal(varname='cols1d_wtgcell', data=cptr%wtgcell, dim1name='column', ncid=ncid, flag='write')
       call ncd_iolocal(varname='cols1d_wtlunit', data=cptr%wtlunit, dim1name='column', ncid=ncid, flag='write')
       if (masterproc) then
!dir$ concurrent
!cdir nodep
          do c = 1,numc
             l = cptr%landunit(c)
             ictype(c) = lptr%itype(l)
          end do
          call map_dc2sn(ictype, ictemp, type1d=namec)
          call ncd_ioglobal(varname='cols1d_itype_lunit', data=ictemp, ncid=ncid, flag='write')
       end if

       ! Write pft info

       call ncd_iolocal(varname='pfts1d_lon', data=pptr%londeg, dim1name='pft', ncid=ncid, flag='write')
       call ncd_iolocal(varname='pfts1d_lat', data=pptr%latdeg, dim1name='pft', ncid=ncid, flag='write')
       call ncd_iolocal(varname='pfts1d_ixy', data=pptr%ixy   , dim1name='pft', ncid=ncid, flag='write')
       call ncd_iolocal(varname='pfts1d_jxy', data=pptr%jxy   , dim1name='pft', ncid=ncid, flag='write')
       if (masterproc) then
          call get_sn_pfts1d(iptemp, type1d='gridcell',nump=nump)
          call ncd_ioglobal(varname='pfts1d_gi', data=iptemp, ncid=ncid, flag='write')
          call get_sn_pfts1d(iptemp, type1d='landunit',nump=nump)
          call ncd_ioglobal(varname='pfts1d_li', data=iptemp, ncid=ncid, flag='write')
          call get_sn_pfts1d(iptemp, type1d='column',nump=nump)
          call ncd_ioglobal(varname='pfts1d_ci', data=iptemp, ncid=ncid, flag='write')
       end if
       call ncd_iolocal(varname='pfts1d_wtgcell'  , data=pptr%wtgcell, dim1name='pft', ncid=ncid, flag='write')
       call ncd_iolocal(varname='pfts1d_wtlunit'  , data=pptr%wtlunit, dim1name='pft', ncid=ncid, flag='write')
       call ncd_iolocal(varname='pfts1d_wtcol'    , data=pptr%wtcol  , dim1name='pft', ncid=ncid, flag='write')
       call ncd_iolocal(varname='pfts1d_itype_veg', data=pptr%itype  , dim1name='pft', ncid=ncid, flag='write')
       if (masterproc) then
!dir$ concurrent
!cdir nodep
          do p = 1,nump
             l = pptr%landunit(p)
             iptype(p) = lptr%itype(l)
          end do
          call map_dc2sn(iptype, iptemp, type1d=namep)
          call ncd_ioglobal(varname='pfts1d_itype_lunit', data=iptemp, ncid=ncid, flag='write')
       end if

       if (masterproc) then
          deallocate(iltemp, ictemp, iptemp, ictype, iptype)
       end if

    end if

  end subroutine hfields_1dinfo

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: htapes_wrapup
!
! !INTERFACE:
  subroutine htapes_wrapup ()
!
! !DESCRIPTION:
! Write and/or dispose history tape(s)
! Determine if next time step is beginning of history interval and if so:
!   increment the current time sample counter, open a new history file
!   and if needed (i.e., when ntim = 1), write history data to current
!   history file, reset field accumulation counters to zero.
! If primary history file is full or at the last time step of the simulation,
!   write restart dataset and close and dispose all history fiels.
! If history file is full or at the last time step of the simulation:
!   close history file and dispose to mass store (only if file is open)
!   and reset time sample counter to zero if file is full.
! Daily-averaged data for the first day in September are written on
!   date = 00/09/02 with mscur = 0.
! Daily-averaged data for the first day in month mm are written on
!   date = yyyy/mm/02 with mscur = 0.
! Daily-averaged data for the 30th day (last day in September) are written
!   on date = 0000/10/01 mscur = 0.
! Daily-averaged data for the last day in month mm are written on
!   date = yyyy/mm+1/01 with mscur = 0.
!
! !USES:
    use clm_varctl   , only : archive_dir, mss_wpass, mss_irt
    use fileutils    , only : set_filename, putfil
    use spmdMod      , only : masterproc
    use time_manager , only : get_nstep, get_curr_date, get_curr_time, get_prev_date
    use shr_const_mod, only : SHR_CONST_CDAY
    use ncdio
    use shr_sys_mod  , only : shr_sys_flush
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: t                      ! tape index
    integer :: f                      ! field index
    integer :: ier                    ! error code
    integer :: nstep                  ! current step
    integer :: day                    ! current day (1 -> 31)
    integer :: mon                    ! current month (1 -> 12)
    integer :: yr                     ! current year (0 -> ...)
    integer :: mcsec                  ! seconds of current date
    integer :: mdcur                  ! current day
    integer :: mscur                  ! seconds of current day
    integer :: mcdate                 ! current date
    integer :: daym1                  ! nstep-1 day (1 -> 31)
    integer :: monm1                  ! nstep-1 month (1 -> 12)
    integer :: yrm1                   ! nstep-1 year (0 -> ...)
    integer :: mcsecm1                ! nstep-1 time of day [seconds]
    integer :: start1                 ! start 1d value
    integer :: count1                 ! count 1d value
    integer :: start2(2)              ! start 2d values
    integer :: count2(2)              ! count 2d values
    character(len=256) :: rem_dir     ! remote (archive) directory
    character(len=256) :: rem_fn      ! remote (archive) filename
    character(len=  8) :: cdate       ! system date
    character(len=  8) :: ctime       ! system time
    character(len=256) :: str         ! global attribute string
    real(r8):: time                   ! current time
    real(r8):: timedata(2)            ! time interval boundaries
    logical :: if_stop                ! true => last time step of run
    logical :: if_disphist(max_tapes) ! true => save and dispose history file
    logical :: if_remvhist(max_tapes) ! true => remove local history file after dispose
    character(len=32) :: subname = 'htapes_wrapup'
!-----------------------------------------------------------------------

    ! get current step

    nstep = get_nstep()

    ! Set calendar for current time step

    call get_curr_date (yr, mon, day, mcsec)
    mcdate = yr*10000 + mon*100 + day

    ! Set calendar for current for previous time step

    call get_prev_date (yrm1, monm1, daym1, mcsecm1)

    ! Set elapased time since reference date

    call get_curr_time (mdcur, mscur)

    ! Loop over active history tapes, create new history files if necessary
    ! and write data to history files if end of history interval.
    do t = 1, ntapes

       ! Skip nstep=0 if monthly average

       if (nstep==0 .and. tape(t)%nhtfrq==0) cycle

       ! Determine if end of history interval
       tape(t)%is_endhist = .false.
       if (tape(t)%nhtfrq==0) then   !monthly average
          if (mon /= monm1) tape(t)%is_endhist = .true.
       else
          if (mod(nstep,tape(t)%nhtfrq) == 0) tape(t)%is_endhist = .true.
       end if

       ! If end of history interval

       if (tape(t)%is_endhist) then

          ! Normalize history buffer if time averaged

          call hfields_normalize(t)

          ! Increment current time sample counter.

          tape(t)%ntimes = tape(t)%ntimes + 1

          ! Create history file if appropriate and build time comment

          ! If first time sample, generate unique history file name, open file
          ! and write time constant fields to file

          if (tape(t)%ntimes == 1) then
             if (masterproc) then
                locfnh(t) = set_hist_filename (hist_freq=tape(t)%nhtfrq, hist_file=t)
                write(6,*)'HTAPES_WRAPUP: Creating history file ', trim(locfnh(t)), &
                     ' at nstep = ',get_nstep()
                call htape_create (t)
             endif

             ! Write time constant history variables to primary tape
             if (t==1) call htape_timeconst(t, mode='write')

          endif

          ! Write time information to history file

          if (masterproc) then

             start1 = tape(t)%ntimes
             count1 = 1

             call get_curr_time (mdcur, mscur)
             call get_curr_date (yr, mon, day, mcsec)
             mcdate = yr*10000 + mon*100 + day

             call check_ret(nf_put_vara_int (nfid(t), mcdate_id(t), start1, count1, mcdate), subname)
             call check_ret(nf_put_vara_int (nfid(t), mcsec_id(t), start1, count1, mcsec), subname)
             call check_ret(nf_put_vara_int (nfid(t), mdcur_id(t), start1, count1, mdcur), subname)
             call check_ret(nf_put_vara_int (nfid(t), mscur_id(t), start1, count1, mscur), subname)
             call check_ret(nf_put_vara_int (nfid(t), nstep_id(t), start1, count1, nstep), subname)

             time = mdcur + mscur/86400._r8
             call check_ret(nf_put_vara_double (nfid(t), time_var_id(t), start1, count1, time), subname)

             start2(1) = 1
             start2(2) = tape(t)%ntimes
             count2(1) = 2
             count2(2) = 1
             timedata(1) = tape(t)%begtime
             timedata(2) = time
             call check_ret(nf_put_vara_double(nfid(t), time_bounds_id(t), start2, count2, timedata), subname)

             start2(1) = 1
             start2(2) = tape(t)%ntimes
             count2(1) = 8
             count2(2) = 1
             call getdatetime (cdate, ctime)
             call check_ret(nf_put_vara_text(nfid(t), date_written_id(t), start2, count2, cdate), subname)
             call check_ret(nf_put_vara_text(nfid(t), time_written_id(t), start2, count2, ctime), subname)

             write(6,*)
             write(6,*)'HTAPES_WRAPUP: Writing current time sample to local history file ', &
                  trim(locfnh(t)),' at nstep = ',get_nstep(), &
                  ' for history time interval beginning at ', tape(t)%begtime, &
                  ' and ending at ',time
             write(6,*)
             call shr_sys_flush(6)

             ! Update beginning time of next interval

             tape(t)%begtime = time

          endif

          ! Write history time samples

          call hfields_write(t, mode='write')

          ! Zero necessary history buffers

          call hfields_zero(t)

       end if

    end do  ! end loop over history tapes

    ! Determine if file needs to be closed and disposed

    call do_close_dispose (ntapes, tape(:)%ntimes, tape(:)%mfilt, if_stop, if_disphist, if_remvhist)

    ! Close open history file and dispose to mass store
    ! Auxilary files may have been closed and saved off without being full,
    ! must reopen the files

    do t = 1, ntapes
       if (if_disphist(t)) then
          if (masterproc) then
             if (tape(t)%ntimes /= 0) then
                write(6,*)
                write(6,*) 'HTAPES_WRAPUP: Closing local history file ',&
                     trim(locfnh(t)),' at nstep = ', get_nstep()
                write(6,*)
                call check_ret(nf_close(nfid(t)), subname)
                rem_dir = trim(archive_dir) // '/hist/'
                rem_fn = set_filename(rem_dir, locfnh(t))
                call putfil (locfnh(t), rem_fn, mss_wpass, mss_irt, if_remvhist(t))
                if (.not.if_stop .and. (tape(t)%ntimes/=tape(t)%mfilt)) then
                   call check_ret(nf_open (trim(locfnh(t)), NF_WRITE, nfid(t)), subname)
                end if
             else
                write(6,*)'HTAPES_WRAPUP: history tape ',t,': no open file to close'
             end if
          endif
       endif
    end do

    ! Determine if time to write restarts before resetting number of
    ! time samples below

#if (defined OFFLINE) || (defined COUP_CAM)
    if_writrest = .false.
    if (tape(1)%ntimes == tape(1)%mfilt) if_writrest = .true.
#endif

    ! Reset number of time samples to zero if file is full

    do t = 1, ntapes
       if (if_disphist(t) .and. tape(t)%ntimes==tape(t)%mfilt) then
          tape(t)%ntimes = 0
       end if
    end do

  end subroutine htapes_wrapup

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_history
!
! !INTERFACE:
  subroutine restart_history (nio, flag)
!
! !DESCRIPTION:
! Read/write history file restart data.
! If the current history file(s) are not full, file(s) are opened
! so that subsequent time samples are added until the file is full.
! A new history file is used on a branch run.
!
! !USES:
    use iobinary
    use ncdio
    use clmtype   , only : nameg, namel, namec, namep, ocnrof, lndrof
    use decompMod , only : get_proc_bounds, get_proc_global
#if (defined RTM)
    use RunoffMod , only : get_proc_rof_bounds, get_proc_rof_global
#endif
    use clm_varctl, only : archive_dir, nsrest, mss_irt
    use fileutils , only : set_filename, getfil
#if (defined SPMD)
    use spmdMod   , only : masterproc, mpicom, MPI_REAL8, MPI_INTEGER, MPI_CHARACTER
#else
    use spmdMod   , only : masterproc
#endif
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nio               !restart unit
    character(len=*), intent(in) :: flag     !'read' or 'write'
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    character(len=256) :: fnameh(max_tapes)      ! full name of history file
    character(len=256) :: rem_dir                ! remote (archive) directory
    character(len=256) :: filename               ! generic filename
    character(len=  8) :: type1d                 ! clm pointer 1d type
    character(len=  8) :: type1d_out             ! hbuf 1d type
    integer :: k                                 ! 1d index
    integer :: lev                               ! level index
    integer :: t                                 ! tape index
    integer :: f                                 ! field index
    integer :: ier                               ! error status
    integer :: num1d,beg1d,end1d                 ! 1d size, beginning and ending indices
    integer :: num1d_out,beg1d_out,end1d_out     ! 1d size, beginning and ending indices
    integer :: num2d                             ! 2d size (e.g. number of vertical levels)
    real(r8), pointer :: hbuf(:,:)               ! history buffer
    integer , pointer :: nacs(:,:)               ! accumulation counter
    integer , pointer :: ibuf1d(:)               ! temporary
    real(r8), pointer :: rbuf1d(:)               ! temporary
    integer , pointer :: ibuf2d(:,:)             ! temporary
    real(r8), pointer :: rbuf2d(:,:)             ! temporary
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer :: numg         ! total number of gridcells across all processors
    integer :: numl         ! total number of landunits across all processors
    integer :: numc         ! total number of columns across all processors
    integer :: nump         ! total number of pfts across all processors
#if (defined RTM)
    integer :: beg_roflnd   ! per-proc beginning land runoff index
    integer :: end_roflnd   ! per-proc ending land runoff index
    integer :: beg_rofocn   ! per-proc beginning ocean runoff index
    integer :: end_rofocn   ! per-proc ending ocean runoff index
    integer :: num_roflnd   ! total number of land runoff points across all procs
    integer :: num_rofocn   ! total number of ocean runoff points across all procs
#endif
    character(len=32) :: subname = 'restart_history'
    integer :: varid(max_flds,max_tapes)   ! no longer used (here for backwards compat)
    integer :: is_endhist_int, dov2xy_int
!------------------------------------------------------------------------

    ! If branch run, initialize file times and return

    if (flag == 'read') then
       if (nsrest == 3) then
          do t = 1,ntapes
             tape(t)%ntimes = 0
          end do
          RETURN
       end if
    endif

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)
#if (defined RTM)
    call get_proc_rof_bounds(beg_roflnd, end_roflnd, beg_rofocn, end_rofocn)
    call get_proc_rof_global(num_roflnd, num_rofocn)
#endif

    ! Read/write history file data only for restart run (not for branch run)

    if (masterproc) then

        if (flag == 'write') then
           varid(:,:) = -999._r8 ! this is done for backwards compatibility of restarts
           write (nio, iostat=ier) ntapes, varid, fincl, fexcl
        else if (flag == 'read') then
           read (nio, iostat=ier) ntapes, varid, fincl, fexcl
        endif
        if (ier /= 0) then
           write (6,*)' RESTART_HISTORY error: read/write error 1',ier,' on i/o unit = ',nio
           call endrun()
        end if

        do t = 1,ntapes
           if (flag == 'write') then
              write (nio, iostat=ier)  &
                   tape(t)%nflds,      &
                   tape(t)%ntimes,     &
                   tape(t)%nhtfrq,     &
                   tape(t)%mfilt,      &
                   tape(t)%ncprec,     &
                   tape(t)%begtime,    &
                   tape(t)%is_endhist, &
                   tape(t)%dov2xy
           else if (flag == 'read') then
              read (nio, iostat=ier)   &
                   tape(t)%nflds,      &
                   tape(t)%ntimes,     &
                   tape(t)%nhtfrq,     &
                   tape(t)%mfilt,      &
                   tape(t)%ncprec,     &
                   tape(t)%begtime,    &
                   tape(t)%is_endhist, &
                   tape(t)%dov2xy
           endif
           if (ier /= 0) then
              write (6,*)' RESTART_HISTORY error: read/write error 2',ier,' on i/o unit = ',nio
              call endrun()
           end if
           if (flag == 'write') then
              if ( tape(t)%is_endhist ) then
                 is_endhist_int = 1
              else
                 is_endhist_int = 0
              end if
              if ( tape(t)%dov2xy ) then
                 dov2xy_int = 1
              else
                 dov2xy_int = 0
              end if
              write (nio, iostat=ier) is_endhist_int, dov2xy_int
           else if (flag == 'read') then
              read (nio, iostat=ier) is_endhist_int, dov2xy_int
              if ( is_endhist_int /= 0 ) then
                 tape(t)%is_endhist = .true.
              else
                 tape(t)%is_endhist = .false.
              end if
              if ( dov2xy_int /= 0 ) then
                 tape(t)%dov2xy = .true.
              else
                 tape(t)%dov2xy = .false.
              end if
           endif
           if (ier /= 0) then
              write (6,*)' RESTART_HISTORY error: read/write error 3',ier,' on i/o unit = ',nio
              call endrun()
           end if
           do f=1,tape(t)%nflds
              if (flag == 'write') then
                 write (nio,iostat=ier) tape(t)%hlist(f)%field%name,           &
                                        tape(t)%hlist(f)%field%long_name,      &
                                        tape(t)%hlist(f)%field%units,          &
                                        tape(t)%hlist(f)%field%num2d,          &
                                        tape(t)%hlist(f)%field%hpindex,        &
                                        tape(t)%hlist(f)%field%type1d,         &
                                        tape(t)%hlist(f)%field%type1d_out,     &
                                        tape(t)%hlist(f)%field%type2d,         &
                                        tape(t)%hlist(f)%field%typexy,         &
                                        tape(t)%hlist(f)%field%nlonxy,         &
                                        tape(t)%hlist(f)%field%nlatxy,         &
                                        tape(t)%hlist(f)%field%p2c_scale_type, &
                                        tape(t)%hlist(f)%field%c2l_scale_type, &
                                        tape(t)%hlist(f)%field%l2g_scale_type, &
                                        tape(t)%hlist(f)%avgflag
              else if (flag == 'read') then
                 read  (nio,iostat=ier) tape(t)%hlist(f)%field%name,           &
                                        tape(t)%hlist(f)%field%long_name,      &
                                        tape(t)%hlist(f)%field%units,          &
                                        tape(t)%hlist(f)%field%num2d,          &
                                        tape(t)%hlist(f)%field%hpindex,        &
                                        tape(t)%hlist(f)%field%type1d,         &
                                        tape(t)%hlist(f)%field%type1d_out,     &
                                        tape(t)%hlist(f)%field%type2d,         &
                                        tape(t)%hlist(f)%field%typexy,         &
                                        tape(t)%hlist(f)%field%nlonxy,         &
                                        tape(t)%hlist(f)%field%nlatxy,         &
                                        tape(t)%hlist(f)%field%p2c_scale_type, &
                                        tape(t)%hlist(f)%field%c2l_scale_type, &
                                        tape(t)%hlist(f)%field%l2g_scale_type, &
                                        tape(t)%hlist(f)%avgflag
              end if
              if (ier /= 0) then
                 write(6,*)'RESTART_HISTORY error: read/write error 4: End or error condition writing ', &
                      'history restart field ',f,' from tape ',t
                 call endrun()
              end if
           end do
        end do

     endif  ! end of if-masterproc block

#if ( defined SPMD)
     if (flag == 'read') then

        ! Broadcast history information from masterprocessor

        call mpi_bcast (ntapes, 1, mpi_integer, 0, mpicom, ier)
        do t = 1,ntapes
           call mpi_bcast (tape(t)%nflds     , 1, MPI_INTEGER, 0, mpicom, ier)
           call mpi_bcast (tape(t)%ntimes    , 1, MPI_INTEGER, 0, mpicom, ier)
           call mpi_bcast (tape(t)%nhtfrq    , 1, MPI_INTEGER, 0 ,mpicom, ier)
           call mpi_bcast (tape(t)%mfilt     , 1, MPI_INTEGER, 0, mpicom, ier)
           call mpi_bcast (tape(t)%ncprec    , 1, MPI_INTEGER, 0, mpicom, ier)
           call mpi_bcast (tape(t)%is_endhist, 1, MPI_LOGICAL, 0, mpicom, ier)
           call mpi_bcast (tape(t)%dov2xy    , 1, MPI_LOGICAL, 0, mpicom, ier)
           call mpi_bcast (tape(t)%begtime   , 1, MPI_REAL8  , 0, mpicom, ier)
           do f = 1,tape(t)%nflds
              call mpi_bcast (tape(t)%hlist(f)%field%name          , max_namlen, MPI_CHARACTER, 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%field%units         , max_chars , MPI_CHARACTER, 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%field%num2d         , 1         , MPI_INTEGER  , 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%field%hpindex       , 1         , MPI_INTEGER  , 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%field%type1d        , 8         , MPI_CHARACTER, 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%field%type1d_out    , 8         , MPI_CHARACTER, 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%field%type2d        , 8         , MPI_CHARACTER, 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%field%typexy        , 8         , MPI_CHARACTER, 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%field%nlonxy        , 1         , MPI_INTEGER  , 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%field%p2c_scale_type, 8         , MPI_CHARACTER, 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%field%c2l_scale_type, 8         , MPI_CHARACTER, 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%field%l2g_scale_type, 8         , MPI_CHARACTER, 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%field%nlatxy        , 1         , MPI_INTEGER  , 0, mpicom, ier)
              call mpi_bcast (tape(t)%hlist(f)%avgflag             , 1         , MPI_CHARACTER, 0, mpicom, ier)
           end do
        end do
     endif
#endif

     ! Allocate memory for history buffers - read only

     if (flag == 'read') then
        do t = 1,ntapes
           do f = 1,tape(t)%nflds

              type1d_out = tape(t)%hlist(f)%field%type1d_out
              select case (type1d_out)
              case (nameg)
                 num1d_out = numg
                 beg1d_out = begg
                 end1d_out = endg
              case (namel)
                 num1d_out = numl
                 beg1d_out = begl
                 end1d_out = endl
              case (namec)
                 num1d_out = numc
                 beg1d_out = begc
                 end1d_out = endc
              case (namep)
                 num1d_out = nump
                 beg1d_out = begp
                 end1d_out = endp
#if (defined RTM)
              case (lndrof)
                 num1d_out = num_roflnd
                 beg1d_out = beg_roflnd
                 end1d_out = end_roflnd
              case (ocnrof)
                 num1d_out = num_rofocn
                 beg1d_out = beg_rofocn
                 end1d_out = end_rofocn
#endif
              case default
                 write(6,*)'RESTART_HISTORY error: read unknown 1d output type=',type1d_out
                 call endrun ()
              end select

              tape(t)%hlist(f)%field%num1d_out = num1d_out
              tape(t)%hlist(f)%field%beg1d_out = beg1d_out
              tape(t)%hlist(f)%field%end1d_out = end1d_out

              num2d  = tape(t)%hlist(f)%field%num2d
              allocate (tape(t)%hlist(f)%hbuf(beg1d_out:end1d_out,num2d), &
                        tape(t)%hlist(f)%nacs(beg1d_out:end1d_out,num2d), stat=ier)
              if (ier /= 0) then
                 write(6,*)' RESTART_HISTORY error: allocation error for hbuf,nacs at t,f=',t,f
                 call endrun()
              endif
              tape(t)%hlist(f)%hbuf(:,:) = 0._r8
              tape(t)%hlist(f)%nacs(:,:) = 0

              type1d = tape(t)%hlist(f)%field%type1d
              select case (type1d)
              case (nameg)
                 num1d = numg
                 beg1d = begg
                 end1d = endg
              case (namel)
                 num1d = numl
                 beg1d = begl
                 end1d = endl
              case (namec)
                 num1d = numc
                 beg1d = begc
                 end1d = endc
              case (namep)
                 num1d = nump
                 beg1d = begp
                 end1d = endp
#if (defined RTM)
              case (lndrof)
                 num1d = num_roflnd
                 beg1d = beg_roflnd
                 end1d = end_roflnd
              case (ocnrof)
                 num1d = num_rofocn
                 beg1d = beg_rofocn
                 end1d = end_rofocn
#endif
              case default
                 write(6,*)'RESTART_HISTORY error: read unknown 1d type=',type1d
                 call endrun ()
              end select

              tape(t)%hlist(f)%field%num1d = num1d
              tape(t)%hlist(f)%field%beg1d = beg1d
              tape(t)%hlist(f)%field%end1d = end1d

           end do
        end do
     endif

     ! Loop over tapes - only read/write accumulators and counters if needed
     ! (if not end of history interval)

     do t = 1,ntapes
        if (.not. tape(t)%is_endhist) then

           if (masterproc) then
              if (flag == 'write') then
                 write(6,*)'RESTART_HISTORY: Writing history restart information for tape ',t
              else
                 write(6,*)'RESTART_HISTORY: Reading history restart information for tape ',t
              endif
           endif

           do f=1,tape(t)%nflds

              type1d_out =  tape(t)%hlist(f)%field%type1d_out
              num1d_out  =  tape(t)%hlist(f)%field%num1d_out
              beg1d_out  =  tape(t)%hlist(f)%field%beg1d_out
              end1d_out  =  tape(t)%hlist(f)%field%end1d_out
              num2d      =  tape(t)%hlist(f)%field%num2d
              nacs       => tape(t)%hlist(f)%nacs
              hbuf       => tape(t)%hlist(f)%hbuf

              if (num2d == 1) then
                 allocate(ibuf1d(num1d_out))
                 allocate(rbuf1d(num1d_out))
                 if (flag == 'read') then
                    call readin (nio, ibuf1d, clmlevel=type1d_out)
                    call readin (nio, rbuf1d, clmlevel=type1d_out)
!dir$ concurrent
!cdir nodep
                    do k = beg1d_out,end1d_out
                       nacs(k,1) = ibuf1d(k)
                       hbuf(k,1) = rbuf1d(k)
                    end do
                 else if (flag == 'write') then
!dir$ concurrent
!cdir nodep
                    do k = beg1d_out,end1d_out
                       ibuf1d(k) = nacs(k,1)
                       rbuf1d(k) = hbuf(k,1)
                    end do
                    call wrtout (nio, ibuf1d, clmlevel=type1d_out)
                    call wrtout (nio, rbuf1d, clmlevel=type1d_out)
                 endif
                 deallocate(ibuf1d)
                 deallocate(rbuf1d)
              else
                 allocate(ibuf2d(num2d,num1d_out))
                 allocate(rbuf2d(num2d,num1d_out))
                 if (flag == 'read') then
                    call readin (nio, ibuf2d, clmlevel=type1d_out)
                    call readin (nio, rbuf2d, clmlevel=type1d_out)
                    do lev = 1,num2d
!dir$ concurrent
!cdir nodep
                       do k = beg1d_out,end1d_out
                          nacs(k,lev) = ibuf2d(lev,k)
                          hbuf(k,lev) = rbuf2d(lev,k)
                       end do
                    end do
                 else if (flag == 'write') then
                    do lev = 1,num2d
!dir$ concurrent
!cdir nodep
                       do k = beg1d_out,end1d_out
                          ibuf2d(lev,k) = nacs(k,lev)
                          rbuf2d(lev,k) = hbuf(k,lev)
                       end do
                    end do
                    call wrtout (nio, ibuf2d, clmlevel=type1d_out)
                    call wrtout (nio, rbuf2d, clmlevel=type1d_out)
                 endif
                 deallocate(ibuf2d)
                 deallocate(rbuf2d)
              endif

           end do   ! end of fields do-loop
        end if   ! end of if-end-of-history-interval block
     end do   ! end of tape do-loop

     ! Read names of history files. If history file is not full,
     ! open history file and obtain time dependent netCDF variable id's

     if (flag == 'read') then
        if (masterproc) then
           do t = 1,ntapes
              read (nio) fnameh(t)
              if (tape(t)%ntimes /= 0) then
                 call getfil (fnameh(t), locfnh(t), 0)
                 call check_ret(nf_open (locfnh(t), nf_write, nfid(t)), subname)
                 call check_ret(nf_inq_varid(nfid(t), 'mcdate', mcdate_id(t)), subname)
                 call check_ret(nf_inq_varid(nfid(t), 'mcsec' , mcsec_id(t)), subname)
                 call check_ret(nf_inq_varid(nfid(t), 'mdcur' , mdcur_id(t)), subname)
                 call check_ret(nf_inq_varid(nfid(t), 'mscur' , mscur_id(t)), subname)
                 call check_ret(nf_inq_varid(nfid(t), 'nstep' , nstep_id(t)), subname)
                 call check_ret(nf_inq_varid(nfid(t), 'time'  , time_var_id(t)), subname)
                 call check_ret(nf_inq_varid(nfid(t), 'time_bounds' , time_bounds_id(t)), subname)
                 call check_ret(nf_inq_varid(nfid(t), 'date_written', date_written_id(t)), subname)
                 call check_ret(nf_inq_varid(nfid(t), 'time_written', time_written_id(t)), subname)
              end if
           end do
        endif
     endif

     ! Write name of current history file(s)

     if (flag == 'write') then
        if (masterproc) then
           do t = 1,ntapes
              if (mss_irt == 0) then
                 filename = locfnh(t)
              else
                 rem_dir = trim(archive_dir) // '/hist/'
                 filename = set_filename(rem_dir, locfnh(t))
              endif
              write(nio) filename
           end do
        endif
     endif

   end subroutine restart_history

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: getname
!
! !INTERFACE:
   character(len=max_namlen) function getname (inname)
!
! !DESCRIPTION:
! Retrieve name portion of inname. If an averaging flag separater character
! is present (:) in inname, lop it off.
!
! !ARGUMENTS:
     implicit none
     character(len=*), intent(in) :: inname
!
! !REVISION HISTORY:
! Created by Jim Rosinski
!
!EOP
!
! LOCAL VARIABLES:
     integer :: length
     integer :: i
!-----------------------------------------------------------------------

     length = len (inname)

     if (length < max_namlen .or. length > max_namlen+2) then
        write(6,*) 'getname: bad length=',length
        call endrun()
     end if

     getname = ' '
     do i = 1,max_namlen
        if (inname(i:i) == ':') exit
        getname(i:i) = inname(i:i)
     end do

   end function getname

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: getflag
!
! !INTERFACE:
   character(len=1) function getflag (inname)
!
! !DESCRIPTION:
! Retrieve flag portion of inname. If an averaging flag separater character
! is present (:) in inname, return the character after it as the flag
!
! !ARGUMENTS:
     implicit none
     character(len=*) inname   ! character string
!
! !REVISION HISTORY:
! Created by Jim Rosinski
!
!EOP
!
! LOCAL VARIABLES:
     integer :: length         ! length of inname
     integer :: i              ! loop index
!-----------------------------------------------------------------------

     length = len (inname)

     if (length < max_namlen .or. length > max_namlen+2) then
        write(6,*) 'getname: bad length=',length
        call endrun()
     end if

     getflag = ' '
     do i = 1,length
        if (inname(i:i) == ':') then
           getflag = inname(i+1:i+1)
           exit
        end if
     end do

   end function getflag

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: list_index
!
! !INTERFACE:
   subroutine list_index (list, name, index)
!
! !DESCRIPTION:
!
! !USES:
!
! !ARGUMENTS:
     implicit none
     character(len=*), intent(in) :: list(max_flds)  ! input list of names, possibly ":" delimited
     character(len=max_namlen), intent(in) :: name   ! name to be searched for
     integer, intent(out) :: index                   ! index of "name" in "list"
!
! !REVISION HISTORY:
! Created by Jim Rosinski
!
!EOP
!
! LOCAL VARIABLES:
     character(len=max_namlen) :: listname           ! input name with ":" stripped off.
     integer f                                       ! field index
!-----------------------------------------------------------------------

     ! Only list items

     index = 0
     do f=1,max_flds
        listname = getname (list(f))
        if (listname == ' ') exit
        if (listname == name) then
           index = f
           exit
        end if
     end do

   end subroutine list_index

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_hist_filename
!
! !INTERFACE:
  character(len=256) function set_hist_filename (hist_freq, hist_file)
!
! !DESCRIPTION:
! Determine history dataset filenames.
!
! !USES:
    use clm_varctl, only : caseid
    use time_manager, only : get_curr_date, get_prev_date
!
! !ARGUMENTS:
   implicit none
   integer, intent(in)  :: hist_freq   !history file frequency
   integer, intent(in)  :: hist_file   !history file index
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
   character(len=256) :: cdate       !date char string
   character(len=  1) :: hist_index  !p,1 or 2 (currently)
   integer :: day                    !day (1 -> 31)
   integer :: mon                    !month (1 -> 12)
   integer :: yr                     !year (0 -> ...)
   integer :: sec                    !seconds into current day
!-----------------------------------------------------------------------

   if (hist_freq == 0 ) then   !monthly
      call get_prev_date (yr, mon, day, sec)
      write(cdate,'(i4.4,"-",i2.2)') yr,mon
   else                        !other
      call get_curr_date (yr, mon, day, sec)
      write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
   endif
   write(hist_index,'(i1.1)') hist_file - 1
   set_hist_filename = "./"//trim(caseid)//".clm2.h"//hist_index//"."//&
        trim(cdate)//".nc"

  end function set_hist_filename

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: add_fld1d
!
! !INTERFACE:
  subroutine add_fld1d (fname, units, avgflag, long_name, &
                        ptr_gcell, ptr_lunit, ptr_col, ptr_pft, ptr_roflnd, &
                        ptr_rofocn, p2c_scale_type, c2l_scale_type, &
                        l2g_scale_type, set_lake, default, typexy)
!
! !DESCRIPTION:
! Initialize a single level history field. The pointer, ptrhist,
! is a pointer to the clmtype array that the history buffer will use.
! The value of type1d passed to masterlist\_add\_fld determines which of the
! 1d type of the output and the beginning and ending indices the history
! buffer field). Default history contents for given field on all tapes
! are set by calling [masterlist\_make\_active] for the appropriate tape.
! After the masterlist is built, routine [htapes\_build] is called for an
! initial or branch run to initialize the actual history tapes.
!
! !USES:
    use clmtype
    use decompMod , only : get_proc_bounds
    use clm_varpar, only : lsmlon, lsmlat, rtmlon, rtmlat
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)           :: fname          ! field name
    character(len=*), intent(in)           :: units          ! units of field
    character(len=1), intent(in)           :: avgflag        ! time averaging flag
    character(len=*), intent(in)           :: long_name      ! long name of field
    real(r8)        , optional, pointer    :: ptr_gcell(:)   ! pointer to gridcell array
    real(r8)        , optional, pointer    :: ptr_lunit(:)   ! pointer to landunit array
    real(r8)        , optional, pointer    :: ptr_col(:)     ! pointer to column array
    real(r8)        , optional, pointer    :: ptr_pft(:)     ! pointer to pft array
    real(r8)        , optional, pointer    :: ptr_roflnd(:)  ! pointer to land channel runoff
    real(r8)        , optional, pointer    :: ptr_rofocn(:)  ! pointer to ocean runoff
    real(r8)        , optional, intent(in) :: set_lake       ! value to set lakes to
    character(len=*), optional, intent(in) :: p2c_scale_type ! scale type for subgrid averaging of pfts to column
    character(len=*), optional, intent(in) :: c2l_scale_type ! scale type for subgrid averaging of columns to landunits
    character(len=*), optional, intent(in) :: l2g_scale_type ! scale type for subgrid averaging of landunits to gridcells
    character(len=*), optional, intent(in) :: default        ! if set to 'inactive, field will not appear on primary tape
    character(len=*), optional, intent(in) :: typexy         ! by default set to 'clm', for other values need to specify
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: p,c,l,g                 ! indices
    integer :: begp, endp              ! per-proc beginning and ending pft indices
    integer :: begc, endc              ! per-proc beginning and ending column indices
    integer :: begl, endl              ! per-proc beginning and ending landunit indices
    integer :: begg, endg              ! per-proc gridcell ending gridcell indices
    integer :: hpindex                 ! history buffer pointer index
    integer :: xynlon                  ! xy grid latitude size
    integer :: xynlat                  ! xy grid latitude size
    character(len=8) :: xytype         ! xy grid output type
    character(len=8) :: type1d         ! 1d type
    character(len=8) :: scale_type_p2c ! scale type for subgrid averaging of pfts to column
    character(len=8) :: scale_type_c2l ! scale type for subgrid averaging of columns to landunits
    character(len=8) :: scale_type_l2g ! scale type for subgrid averaging of landunits to gridcells
!------------------------------------------------------------------------

    ! Determine processor bounds

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! History buffer pointer

    hpindex = pointer_index()

    if (present(ptr_gcell)) then

       type1d = nameg
       clmptr_rs(hpindex)%ptr => ptr_gcell

    else if (present(ptr_lunit)) then

       type1d = namel
       clmptr_rs(hpindex)%ptr => ptr_lunit
       if (present(set_lake)) then
          do l = begl,endl
             if (clm3%g%l%lakpoi(l)) ptr_lunit(l) = set_lake
          end do
       end if

    else if (present(ptr_col)) then

       type1d = namec
       clmptr_rs(hpindex)%ptr => ptr_col
       if (present(set_lake)) then
          do c = begc,endc
             l = clm3%g%l%c%landunit(c)
             if (clm3%g%l%lakpoi(l)) ptr_col(c) = set_lake
          end do
       end if

    else if (present(ptr_pft)) then

       type1d = namep
       clmptr_rs(hpindex)%ptr => ptr_pft
       if (present(set_lake)) then
          do p = begp,endp
             l = clm3%g%l%c%p%landunit(p)
             if (clm3%g%l%lakpoi(l)) ptr_pft(p) = set_lake
          end do
       end if

#if (defined RTM)
    else if (present(ptr_roflnd)) then

       type1d = lndrof
       clmptr_rs(hpindex)%ptr => ptr_roflnd

    else if (present(ptr_rofocn)) then

       type1d = ocnrof
       clmptr_rs(hpindex)%ptr => ptr_rofocn
#endif
    else

       write(6,*)'ADDFLD_1D error: must specify a valid pointer index'
       write(6,*)'choices are [ptr_pft, ptr_lunit, ptr_col, ptr_pft] ', &
             'and if RTM is defined also [ptr_roflnd.ptr_rofocn]'
       call endrun()

    end if

    ! Set output grid(xy) type and number of longitudes and latitudes of output grid

    xytype = 'clm'
    xynlon = lsmlon
    xynlat = lsmlat

    if (present(typexy)) then
       if (typexy == 'rof') then
          xytype = typexy
          xynlon = rtmlon
          xynlat = rtmlat
       else
          write(6,*)'ADD_FLD1D error: currently only an optional xy grid type of [rof] is accepted'
          call endrun()
       end if
    end if

    ! Set scaling factor

    scale_type_p2c = 'unity'
    scale_type_c2l = 'unity'
    scale_type_l2g = 'unity'

    if (present(p2c_scale_type)) scale_type_p2c = p2c_scale_type
    if (present(p2c_scale_type)) scale_type_c2l = c2l_scale_type
    if (present(p2c_scale_type)) scale_type_l2g = l2g_scale_type

    ! Add field to masterlist

    call masterlist_addfld (fname=trim(fname), type1d=type1d, type2d='unset', num2d=1, &
         typexy=xytype, nlonxy=xynlon, nlatxy=xynlat, &
         units=units, avgflag=avgflag, long_name=long_name, hpindex=hpindex, &
         p2c_scale_type=scale_type_p2c, c2l_scale_type=scale_type_c2l, l2g_scale_type=scale_type_l2g)

    if (present(default)) then
       if (trim(default) == 'inactive') return
    else
       call masterlist_make_active (name=trim(fname), tape_index=1)
    end if

  end subroutine add_fld1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: add_fld2d
!
! !INTERFACE:
  subroutine add_fld2d (fname, type2d, units, avgflag, long_name, &
                        ptr_gcell, ptr_lunit, ptr_col, ptr_pft, &
                        p2c_scale_type, c2l_scale_type, l2g_scale_type, &
                        set_lake, default, typexy)
!
! !DESCRIPTION:
! Initialize a single level history field. The pointer, ptrhist,
! is a pointer to the clmtype array that the history buffer will use.
! The value of type1d passed to masterlist\_add\_fld determines which of the
! 1d type of the output and the beginning and ending indices the history
! buffer field). Default history contents for given field on all tapes
! are set by calling [masterlist\_make\_active] for the appropriatae tape.
! After the masterlist is built, routine [htapes\_build] is called for an
! initial or branch run to initialize the actual history tapes.
!
! !USES:
    use clmtype
    use decompMod , only : get_proc_bounds
    use clm_varpar, only : lsmlon, lsmlat, nlevsoi, nlevlak, numrad 
#if (defined CASA)
    use CASAMod,    only : nlive, npools
#endif
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fname                    ! field name
    character(len=*), intent(in) :: type2d                   ! 2d output type
    character(len=*), intent(in) :: units                    ! units of field
    character(len=1), intent(in) :: avgflag                  ! time averaging flag
    character(len=*), intent(in) :: long_name                ! long name of field
    real(r8)        , optional, pointer    :: ptr_gcell(:,:) ! pointer to gridcell array
    real(r8)        , optional, pointer    :: ptr_lunit(:,:) ! pointer to landunit array
    real(r8)        , optional, pointer    :: ptr_col(:,:)   ! pointer to column array
    real(r8)        , optional, pointer    :: ptr_pft(:,:)   ! pointer to pft array
    real(r8)        , optional, intent(in) :: set_lake       ! value to set lakes to
    character(len=*), optional, intent(in) :: p2c_scale_type ! scale type for subgrid averaging of pfts to column
    character(len=*), optional, intent(in) :: c2l_scale_type ! scale type for subgrid averaging of columns to landunits
    character(len=*), optional, intent(in) :: l2g_scale_type ! scale type for subgrid averaging of landunits to gridcells
    character(len=*), optional, intent(in) :: default        ! if set to 'inactive, field will not appear on primary tape
    character(len=*), optional, intent(in) :: typexy         ! by default set to 'clm', for other values need to specify
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: p,c,l,g                 ! indices
    integer :: num2d                   ! size of second dimension (e.g. number of vertical levels)
    integer :: begp, endp              ! per-proc beginning and ending pft indices
    integer :: begc, endc              ! per-proc beginning and ending column indices
    integer :: begl, endl              ! per-proc beginning and ending landunit indices
    integer :: begg, endg              ! per-proc gridcell ending gridcell indices
    integer :: hpindex                 ! history buffer index
    integer :: xynlon                  ! xy grid latitude size
    integer :: xynlat                  ! xy grid latitude size
    character(len=8) :: xytype         ! xy grid output type
    character(len=8) :: type1d         ! 1d output type
    character(len=8) :: scale_type_p2c ! scale type for subgrid averaging of pfts to column
    character(len=8) :: scale_type_c2l ! scale type for subgrid averaging of columns to landunits
    character(len=8) :: scale_type_l2g ! scale type for subgrid averaging of landunits to gridcells
!------------------------------------------------------------------------

    ! Determine second dimension size

    select case (type2d)
    case ('levsoi')
       num2d = nlevsoi
    case ('levlak')
       num2d = nlevlak
    case ('numrad')
       num2d = numrad
#if (defined CASA)
    case ('nlive')
       num2d = nlive
    case ('npools')
       num2d = npools
#endif
    case default
       write(6,*)'ADD_FLD2D error: unsupported 2d type ',type2d
       write(6,*)'currently supported types for multi level fields are [levsoi,levlak,numrad', &
#if (defined CASA)
          ',nlive,npools', &
#endif
          ']'
       call endrun()
    end select

    ! Determine processor bounds

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! History buffer pointer

    hpindex = pointer_index()

    if (present(ptr_gcell)) then

       type1d = nameg
       clmptr_ra(hpindex)%ptr => ptr_gcell

    else if (present(ptr_lunit)) then

       type1d = namel
       clmptr_ra(hpindex)%ptr => ptr_lunit
       if (present(set_lake)) then
          do l = begl,endl
             if (clm3%g%l%lakpoi(l)) ptr_lunit(l,:) = set_lake
          end do
       end if

    else if (present(ptr_col)) then

       type1d = namec
       clmptr_ra(hpindex)%ptr => ptr_col
       if (present(set_lake)) then
          do c = begc,endc
             l = clm3%g%l%c%landunit(c)
             if (clm3%g%l%lakpoi(l)) ptr_col(c,:) = set_lake
          end do
       end if

    else if (present(ptr_pft)) then

       type1d = namep
       clmptr_ra(hpindex)%ptr => ptr_pft
       if (present(set_lake)) then
          do p = begp,endp
             l = clm3%g%l%c%p%landunit(p)
             if (clm3%g%l%lakpoi(l)) ptr_pft(p,:) = set_lake
          end do
       end if

    else

       write(6,*)'addfld_2d error: must specify a valid pointer index'
       write(6,*)'choices are ptr_pft, ptr_lunit, ptr_col, ptr_pft '
       call endrun()

    end if

    ! Set output grid(xy) type and number of longitudes and latitudes of output grid

    xytype = 'clm'
    xynlon = lsmlon
    xynlat = lsmlat

    if (present(typexy)) then
       write(6,*)'addfld_2d error: currently no optional grid type is accepted for 2d fields'
       call endrun()
    end if

    ! Set scaling factor

    scale_type_p2c = 'unity'
    scale_type_c2l = 'unity'
    scale_type_l2g = 'unity'

    if (present(p2c_scale_type)) scale_type_p2c = p2c_scale_type
    if (present(p2c_scale_type)) scale_type_c2l = c2l_scale_type
    if (present(p2c_scale_type)) scale_type_l2g = l2g_scale_type

    ! Add field to masterlist

    call masterlist_addfld (fname=trim(fname), type1d=type1d, type2d=type2d, num2d=num2d, &
         typexy=xytype, nlonxy=xynlon, nlatxy=xynlat, &
         units=units, avgflag=avgflag, long_name=long_name, hpindex=hpindex, &
         p2c_scale_type=scale_type_p2c, c2l_scale_type=scale_type_c2l, l2g_scale_type=scale_type_l2g)

    if (present(default)) then
       if (trim(default) == 'inactive') return
    else
       call masterlist_make_active (name=trim(fname), tape_index=1)
    end if

  end subroutine add_fld2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: pointer_index
!
! !INTERFACE:
  integer function pointer_index ()
!
! !DESCRIPTION:
! Set the current pointer index and increment the value of the index.
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
    integer, save :: lastindex = 1
!-----------------------------------------------------------------------

    pointer_index = lastindex
    lastindex = lastindex + 1
    if (lastindex > max_mapflds) then
       write(6,*)'error pointer_index: ',&
            ' lastindex = ',lastindex,' greater than max_mapflds= ',max_mapflds
       call endrun()
    endif

  end function pointer_index

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: add_subscript
!
! !INTERFACE:
  subroutine add_subscript(subname, subdim)
!
! !DESCRIPTION:
! Add a history variable to the output history tape.
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: subname ! name of subscript
    integer         , intent(in) :: subdim  ! dimension of subscript
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
!-----------------------------------------------------------------------

    num_subs = num_subs + 1
    if (num_subs > max_subs) then
       write(6,*)'error add_subscript: ',&
            ' num_subs = ',num_subs,' greater than max_subs= ',max_subs
       call endrun()
    endif
    subs_name(num_subs) = subname
    subs_dim(num_subs) =  subdim

  end subroutine add_subscript

end module histFileMod
