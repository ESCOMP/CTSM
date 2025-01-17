module histFileMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing methods to for CLM history file handling.
  ! See 'history_tape' type for more details.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_sys_mod    , only : shr_sys_flush
  use spmdMod        , only : masterproc
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog, use_fates, compname, use_cn, use_crop
  use clm_varcon     , only : spval, ispval
  use clm_varcon     , only : grlnd, nameg, namel, namec, namep
  use decompMod      , only : get_proc_bounds, get_proc_global, bounds_type, get_global_index, get_global_index_array
  use decompMod      , only : subgrid_level_gridcell, subgrid_level_landunit, subgrid_level_column
  use GridcellType   , only : grc
  use LandunitType   , only : lun
  use ColumnType     , only : col
  use PatchType      , only : patch
  use EDParamsMod    , only : nclmax
  use EDParamsMod    , only : nlevleaf
  use FatesInterfaceTypesMod , only : nlevsclass, nlevage, nlevcoage
  use FatesInterfaceTypesMod , only : nlevheight
  use FatesInterfaceTypesMod , only : nlevdamage
  use FatesConstantsMod      , only : n_landuse_cats
  use FatesFuelClassesMod    , only : num_fuel_classes
  use FatesLitterMod         , only : ncwd
  use PRTGenericMod          , only : num_elements_fates  => num_elements
  use FatesInterfaceTypesMod , only : numpft_fates => numpft
  use ncdio_pio

  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  !
  ! Constants
  !
  integer , public, parameter :: max_tapes = 10         ! max number of history tapes
  integer , public, parameter :: max_flds = 2500        ! max number of history fields
  integer , public, parameter :: max_namlen = 64        ! maximum number of characters for field name
  integer , public, parameter :: scale_type_strlen = 32 ! maximum number of characters for scale types
  integer , private, parameter :: avgflag_strlen = 10   ! maximum number of characters for avgflag
  integer , private, parameter :: hist_dim_name_length = 16 ! lenngth of character strings in dimension names

  ! Possible ways to treat multi-layer snow fields at times when no snow is present in a
  ! given layer. Note that the public parameters are the only ones that can be used by
  ! calls to hist_addfld2d; the private parameters are just used internally by the
  ! histFile implementation.
  integer , private, parameter :: no_snow_MIN = 1                 ! minimum valid value for this flag
  integer , public , parameter :: no_snow_normal = 1              ! normal treatment, which should be used for most fields (use spval when snow layer not present)
  integer , public , parameter :: no_snow_zero = 2                ! average in a 0 value for times when the snow layer isn't present
  integer , private, parameter :: no_snow_MAX = 2                 ! maximum valid value for this flag
  integer , private, parameter :: no_snow_unset = no_snow_MIN - 1 ! flag specifying that field is NOT a multi-layer snow field
  !
  ! Counters
  !
  ! ntapes gives the index of the max history file requested. There can be "holes" in the
  ! numbering - e.g., we can have h0, h1 and h3 tapes, but no h2 tape (because there are
  ! no fields on the h2 tape). In this case, ntapes will be 4 (for h0, h1, h2 and h3,
  ! since h3 is the last requested file), not 3 (the number of files actually produced).
  integer , private :: ntapes = 0        ! index of max history file requested
  !
  ! Namelist
  !
  integer :: ni                          ! implicit index below
  integer, public :: &
       hist_ndens(max_tapes) = 2         ! namelist: output density of netcdf history files
  integer, public :: &
       hist_mfilt(max_tapes) = (/ 1, (30, ni=2, max_tapes)/)        ! namelist: number of time samples per tape
  logical, public :: &
       hist_dov2xy(max_tapes) = (/.true.,(.true.,ni=2,max_tapes)/) ! namelist: true=> do grid averaging
  integer, public :: &
       hist_nhtfrq(max_tapes) = (/0, (-24, ni=2,max_tapes)/)        ! namelist: history write freq(0=monthly)
  character(len=avgflag_strlen), public :: &
       hist_avgflag_pertape(max_tapes) = (/(' ',ni=1,max_tapes)/)   ! namelist: per tape averaging flag
  character(len=max_namlen), public :: &
       hist_type1d_pertape(max_tapes)  = (/(' ',ni=1,max_tapes)/)   ! namelist: per tape type1d

  logical, public :: &
       hist_empty_htapes  = .false.      ! namelist: disable default-active history fields (which
                                         ! only exist on history tape 1). Use hist_fincl1 to enable
                                         ! select fields on top of this.

  character(len=max_namlen+2), public :: &
       hist_fincl1(max_flds) = ' '       ! namelist: list of fields to include in history tape 1
                                         !            aka 'h0' history file.
  character(len=max_namlen+2), public :: &
       hist_fincl2(max_flds) = ' '       ! namelist: list of fields to include in history tape 2
  character(len=max_namlen+2), public :: &
       hist_fincl3(max_flds) = ' '       ! namelist: list of fields to include in history tape 3
  character(len=max_namlen+2), public :: &
       hist_fincl4(max_flds) = ' '       ! namelist: list of fields to include in history tape 4
  character(len=max_namlen+2), public :: &
       hist_fincl5(max_flds) = ' '       ! namelist: list of fields to include in history tape 5
  character(len=max_namlen+2), public :: &
       hist_fincl6(max_flds) = ' '       ! namelist: list of fields to include in history tape 6
  character(len=max_namlen+2), public :: &
       hist_fincl7(max_flds) = ' '       ! namelist: list of fields to include in history tape 7
  character(len=max_namlen+2), public :: &
       hist_fincl8(max_flds) = ' '       ! namelist: list of fields to include in history tape 8
  character(len=max_namlen+2), public :: &
       hist_fincl9(max_flds) = ' '       ! namelist: list of fields to include in history tape 9
  character(len=max_namlen+2), public :: &
       hist_fincl10(max_flds) = ' '      ! namelist: list of fields to include in history tape 10

  character(len=max_namlen+2), public :: &
       fincl(max_flds,max_tapes)         ! copy of hist_fincl* fields in 2-D format. Note Fortran
                                         ! used to have a bug in 2-D namelists, thus this workaround.

  character(len=max_namlen+2), public :: &
       hist_fexcl1(max_flds) = ' ' ! namelist: list of fields to exclude from history tape 1
                                   !           aka 'h0' history file.
  character(len=max_namlen+2), public :: &
       hist_fexcl2(max_flds) = ' ' ! namelist: list of fields to exclude from history tape 2
  character(len=max_namlen+2), public :: &
       hist_fexcl3(max_flds) = ' ' ! namelist: list of fields to exclude from history tape 3
  character(len=max_namlen+2), public :: &
       hist_fexcl4(max_flds) = ' ' ! namelist: list of fields to exclude from history tape 4
  character(len=max_namlen+2), public :: &
       hist_fexcl5(max_flds) = ' ' ! namelist: list of fields to exclude from history tape 5
  character(len=max_namlen+2), public :: &
       hist_fexcl6(max_flds) = ' ' ! namelist: list of fields to exclude from history tape 6
  character(len=max_namlen+2), public :: &
       hist_fexcl7(max_flds) = ' ' ! namelist: list of fields to exclude from history tape 7
  character(len=max_namlen+2), public :: &
       hist_fexcl8(max_flds) = ' ' ! namelist: list of fields to exclude from history tape 8
  character(len=max_namlen+2), public :: &
       hist_fexcl9(max_flds) = ' ' ! namelist: list of fields to exclude from history tape 9
  character(len=max_namlen+2), public :: &
       hist_fexcl10(max_flds) = ' ' ! namelist: list of fields to exclude from history tape 10

  character(len=max_namlen+2), public :: &
       fexcl(max_flds,max_tapes)         ! copy of hist_fexcl* fields in 2-D format. Note Fortran
                                         ! used to have a bug in 2-D namelists, thus this workaround.

  logical, private :: if_disphist(max_tapes)   ! restart, true => save history file
  !
  ! !PUBLIC MEMBER FUNCTIONS:  (in rough call order)
  public :: hist_addfld1d        ! Add a 1d single-level field to the list of all history fields
  public :: hist_addfld2d        ! Add a 2d multi-level field to the list of all history fields
  public :: hist_addfld_decomp   ! Add a 1d/2d field based on patch or column data


  public :: hist_printflds       ! Print summary of list of all history fields
  public :: htapes_fieldlist     ! Finalize history file field lists, intersecting allhistfldlist with
                                 ! namelist params.

  public :: hist_htapes_build    ! Initialize history file handler (for initial or continued run)
  public :: hist_update_hbuf     ! Accumulate into history buffer (all fields and tapes)
  public :: hist_htapes_wrapup   ! Write history tape(s)

  public :: hist_restart_ncd     ! Read/write history file restart data
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: is_mapping_upto_subgrid   ! Is this field being mapped up to a higher subgrid level?
  private :: allhistfldlist_make_active    ! Declare a single field active for a single tape
  private :: allhistfldlist_addfld         ! Add a field to the list of all history fields
  private :: allhistfldlist_change_timeavg ! Override default history tape contents for specific tape
  private :: htape_addfld              ! Transfer field metadata from allhistfldlist to a history tape.
  private :: htape_create              ! Define netcdf metadata of history file t
  private :: htape_add_ltype_metadata  ! Add global metadata defining landunit types
  private :: htape_add_ctype_metadata  ! Add global metadata defining column types
  private :: htape_add_natpft_metadata ! Add global metadata defining natpft types
  private :: htape_add_cft_metadata    ! Add global metadata defining cft types
  private :: htape_timeconst           ! Write time constant values to history tape
  private :: htape_timeconst3D         ! Write time constant 3D values to primary history tape
  private :: hfields_normalize         ! Normalize history file fields by number of accumulations
  private :: hfields_zero              ! Zero out accumulation and hsitory buffers for a tape
  private :: hfields_write             ! Write a variable to a history tape
  private :: hfields_1dinfo            ! Define/output 1d subgrid info if appropriate
  private :: hist_update_hbuf_field_1d ! Updates history buffer for specific field and tape
  private :: hist_update_hbuf_field_2d ! Updates history buffer for specific field and tape
  private :: calc_weight_local_time    ! Calculate weight for time interpolation for local time flag
  private :: hist_set_snow_field_2d    ! Set values in history field dimensioned by levsno
  private :: list_index                ! Find index of field in exclude list
  private :: set_hist_filename         ! Determine history dataset filenames
  public  :: getname                   ! Retrieve name portion of input "inname" (PUBLIC FOR FATES)
  private :: getflag                   ! Retrieve flag
  private :: next_history_pointer_index ! Latest index into raw history data (clmptr_r*) arrays
  private :: max_nFields               ! The max number of fields on any tape
  private :: avgflag_valid             ! Whether a given avgflag is a valid option
  private :: add_landunit_mask_metadata ! Add landunit_mask metadata for the given history field
  !
  ! !PRIVATE TYPES:
  ! Constants
  !
  integer, parameter :: max_length_filename = 199 ! max length of a filename. on most linux systems this
                                                  ! is 255. But this can't be increased until all hard
                                                  ! coded values throughout the i/o stack are updated.
  integer, parameter :: max_chars = 199        ! max chars for char variables

  ! type2d value for a field without a level dimension. This value is important for the
  ! following reasons (as of 2023-08-21):
  ! - type2d is used to determine the sort order of history fields both within the history
  !   file (e.g., what you see from 'ncdump -h') and in the documentation that lists all
  !   history fields. For these purposes, it is important that variables with
  !   type2d_unset appear before variables with a real type2d, so type2d_unset should
  !   appear early in alphabetical sort order. (If type2d_unset were changed to something
  !   that appeared later in alphabetical sort order, then sort_hist_list should be
  !   changed to have some special handling of fields with type2d_unset, forcing them to
  !   appear first.)
  ! - This will soon be added to the history field documentation, so should be a sensible
  !   value for the type2d column in that output.
  character(len=*), parameter :: type2d_unset = '-'
  !
  type field_info
     character(len=max_namlen) :: name         ! field name
     character(len=max_chars)  :: long_name    ! long name
     character(len=max_chars)  :: units        ! units
     character(len=hist_dim_name_length) :: type1d                ! pointer to first dimension type from data type (nameg, etc)
     character(len=hist_dim_name_length) :: type1d_out            ! hbuf first dimension type from data type (nameg, etc)
     character(len=hist_dim_name_length) :: type2d                ! hbuf second dimension type ["levgrnd","levlak","numrad","ltype","natpft","cft","glc_nec","elevclas","subname(n)","mxsowings","mxharvests"]
     integer :: beg1d                          ! on-node 1d clm pointer start index
     integer :: end1d                          ! on-node 1d clm pointer end index
     integer :: num1d                          ! size of clm pointer first dimension (all nodes)
     integer :: beg1d_out                      ! on-node 1d hbuf pointer start index
     integer :: end1d_out                      ! on-node 1d hbuf pointer end index
     integer :: num1d_out                      ! size of hbuf first dimension (all nodes)
     integer :: numdims                        ! the actual number of dimensions, this allows
                                               ! for 2D arrays, where the second dimension is allowed
                                               ! to be 1
     integer :: num2d                          ! size of hbuf second dimension (e.g. number of vertical levels)
     integer :: hpindex                        ! index into raw history data (clmptr_r*) arrays
     character(len=scale_type_strlen) :: p2c_scale_type       ! scale factor when averaging patch to column
     character(len=scale_type_strlen) :: c2l_scale_type       ! scale factor when averaging column to landunit
     character(len=scale_type_strlen) :: l2g_scale_type       ! scale factor when averaging landunit to gridcell
     integer :: no_snow_behavior               ! for multi-layer snow fields, flag saying how to treat times when a given snow layer is absent
  end type field_info

  ! Metadata about a single history field.
  type, abstract :: entry_base
     type (field_info) :: field                ! field information
  contains
     procedure(copy_entry_interface), deferred :: copy
  end type entry_base

  abstract interface
     subroutine copy_entry_interface(this, other)
        ! set this = other
        import :: entry_base
        class(entry_base), intent(out) :: this
        class(entry_base), intent(in) :: other
     end subroutine copy_entry_interface
  end interface

  ! Additional per-field metadata. See also history_entry. 
  ! For the primary history tape, some fields are enabled here (inside hist_addfld* 
  ! call)  but then can be overridden by namelist params (like hist_fincl1). The
  ! fields for other history tapes are theoretically settable here but in
  ! practice are all disabled.  Fields for those tapes have to be specified
  ! explicitly and manually via hist_fincl2 et al.
  type, extends(entry_base) :: allhistfldlist_entry
     logical :: actflag(max_tapes)  ! which history tapes to write to. 
     character(len=avgflag_strlen) :: avgflag(max_tapes)  ! type of time averaging
  contains
     procedure :: copy => copy_allhistfldlist_entry
  end type allhistfldlist_entry

  ! Actual per-field history data, accumulated from clmptr_r* vars. See also allhistfldlist_entry.
  type, extends(entry_base) :: history_entry
     character(len=avgflag_strlen) :: avgflag  ! time averaging flag ("X","A","M","I","SUM")
     real(r8), pointer :: hbuf(:,:)            ! history buffer (dimensions: dim1d x num2d)
     integer , pointer :: nacs(:,:)            ! accumulation counter (dimensions: dim1d x num2d)
  contains
     procedure :: copy => copy_history_entry
  end type history_entry

  ! Each 'history tape' accumulates output values for a set of fields marked 'active' for this run,
  ! at a given time frequency and precision.  The first ('primary') tape defaults to a non-empty set
  ! of active fields (see hist_addfld* methods), overridable by namelist flags,  while the other 
  ! tapes are entirely manually configured via namelist flags. The set of active fields across all
  ! tapes is assembled in the 'allhistfldlist' variable. Note that the first history tape is index 1 in
  ! the code but contains 'h0' in its output filenames (see set_hist_filename method).
  type history_tape
     integer  :: nflds                         ! number of active fields on tape
     integer  :: ntimes                        ! current number of time samples on tape
     integer  :: mfilt                         ! maximum number of time samples per tape
     integer  :: nhtfrq                        ! number of time samples per tape
     integer  :: ncprec                        ! netcdf output precision
     logical  :: dov2xy                        ! true => do xy average for all fields
     logical  :: is_endhist                    ! true => current time step is end of history interval
     real(r8) :: begtime                       ! time at beginning of history averaging interval
     type (history_entry) :: hlist(max_flds)   ! array of active history tape entries.
                                               ! The ordering matches the allhistfldlist's.
  end type history_tape

  type clmpoint_rs                             ! Pointer to real scalar data (1D)
     real(r8), pointer :: ptr(:)
  end type clmpoint_rs
  type clmpoint_ra                             ! Pointer to real array data (2D)
     real(r8), pointer :: ptr(:,:)
  end type clmpoint_ra

  ! Raw history field data (not accumulated). One entry per history field, indexed by 'hpindex'
  ! aka the history pointer index. For accumulated values see 'tape'.  
  integer, parameter :: max_mapflds = 2500     ! Maximum number of fields to track
  type (clmpoint_rs) :: clmptr_rs(max_mapflds) ! Real scalar data (1D)
  type (clmpoint_ra) :: clmptr_ra(max_mapflds) ! Real array data (2D)
  !
  ! History field metadata including which history tapes (if any) it should be output to, and
  ! type of accumulation to perform. This list contains all possible fields, and their field ordering 
  ! is arbitrary, as it depends on the order of hist_addfld* calls in the code.
  ! For the field data itself, see 'tape'.
  !
  type (allhistfldlist_entry) :: allhistfldlist(max_flds)  ! list of all history fields
  !
  ! Whether each history tape is in use in this run. If history_tape_in_use(i) is false,
  ! then data in tape(i) is undefined and should not be referenced.
  !
  logical :: history_tape_in_use(max_tapes)  ! whether each history tape is in use in this run
  !
  ! The actual (accumulated) history data for all active fields in each in-use tape. See
  ! 'history_tape_in_use' for in-use tapes, and 'allhistfldlist' for active fields. See also
  ! clmptr_r* variables for raw history data.
  ! 
  type (history_tape) :: tape(max_tapes)       ! array of history tapes
  !
  ! Namelist input
  !
  ! Counters
  !
  integer :: nallhistflds = 0                        ! number of fields in list of all history fields
  !
  ! Other variables
  !
  character(len=max_length_filename) :: locfnh(max_tapes)  ! local history file names
  character(len=max_length_filename) :: locfnhr(max_tapes) ! local history restart file names
  logical :: htapes_defined = .false.        ! flag indicates history output fields have been defined
  !
  ! NetCDF  Id's
  !
  type(file_desc_t), target :: nfid(max_tapes)       ! file ids
  type(file_desc_t), target :: ncid_hist(max_tapes)  ! file ids for history restart files
  integer :: time_dimid                      ! time dimension id
  integer :: nbnd_dimid                      ! time bounds dimension id
  integer :: strlen_dimid                    ! string dimension id
  !
  ! Time Constant variable names and filename
  !
  character(len=max_chars) :: TimeConst3DVars_Filename = ' '
  !
  ! time_period_freq variable
  !
  character(len=max_chars) :: time_period_freq         = ' '

  character(len=max_chars) :: TimeConst3DVars          = ' '

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine hist_printflds()
    !
    ! !DESCRIPTION:
    ! Print summary of list of all history fields.
    !
    ! !USES:
    use clm_varctl, only: hist_fields_list_file
    use fileutils, only: getavu, relavu
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer, parameter :: ncol = 5  ! number of table columns
    integer nf, i, j  ! do-loop counters
    integer hist_fields_file  ! file unit number
    integer width_col(ncol)  ! widths of table columns
    integer width_col_sum  ! widths of columns summed, including spaces
    character(len=3) str_width_col(ncol)  ! string version of width_col
    character(len=3) str_w_col_sum  ! string version of width_col_sum
    character(len=7) file_identifier  ! fates identifier used in file_name
    character(len=26) file_name  ! hist_fields_file.rst with or without fates
    character(len=99) fmt_txt  ! format statement
    character(len=*),parameter :: subname = 'CLM_hist_printflds'
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) trim(subname),' : number of history fields = ',nallhistflds
       write(iulog,*)' ******* LIST OF ALL HISTORY FIELDS *******'
       do nf = 1,nallhistflds
          write(iulog,9000)nf, allhistfldlist(nf)%field%name, allhistfldlist(nf)%field%units
9000      format (i5,1x,a32,1x,a16)
       end do
       call shr_sys_flush(iulog)
    end if

    ! Print list of all history fields in separate text file when namelist
    ! variable requests it. Text file is formatted in the .rst
    ! (reStructuredText) format for easy introduction of the file to
    ! the CTSM's web-based documentation.

    ! First sort the list to be in alphabetical order
    call sort_hist_list(1, nallhistflds, allhistfldlist)

    if (masterproc .and. hist_fields_list_file) then
       ! Hardwired table column widths to fit the table on a computer
       ! screen. Some strings will be truncated as a result of the
       ! current choices (35, 16, 94, 65, 7). In sphinx (ie the web-based
       ! documentation), text that has not been truncated will wrap
       ! around in the available space.
       width_col(1) = 35  ! variable name column
       width_col(2) = hist_dim_name_length  ! level dimension column
       width_col(3) = 94  ! long description column
       width_col(4) = 65  ! units column
       width_col(5) = 7  ! active (T or F) column
       width_col_sum = sum(width_col) + ncol - 1  ! sum of widths & blank spaces

       ! Convert integer widths to strings for use in format statements
       ! These write statements are not outputting to files
       do i = 1, ncol
          write(str_width_col(i),'(i0)') width_col(i)
       end do
       write(str_w_col_sum,'(i0)') width_col_sum

       ! Open hist_fields_file
       hist_fields_file = getavu()  ! get next available file unit number
       if (use_fates) then
          file_identifier = 'fates'
       else
          file_identifier = 'nofates'
       end if
       file_name = 'history_fields_' // trim(file_identifier) // '.rst'
       open(unit = hist_fields_file, file = file_name,  &
            status = 'replace', action = 'write', form = 'formatted')

       ! File title
       fmt_txt = '(a)'
       write(hist_fields_file,fmt_txt) '============================='
       write(hist_fields_file,fmt_txt) 'CTSM History Fields (' // trim(file_identifier) // ')'
       write(hist_fields_file,fmt_txt) '============================='
       write(hist_fields_file,*)

       ! A warning message and flags from the current CTSM case
       write(hist_fields_file,fmt_txt) 'CAUTION: Not all variables are relevant / present for all CTSM cases.'
       write(hist_fields_file,fmt_txt) 'Key flags used in this CTSM case:'
       fmt_txt = '(a,l)'
       write(hist_fields_file,fmt_txt) 'use_cn = ', use_cn
       write(hist_fields_file,fmt_txt) 'use_crop = ', use_crop
       write(hist_fields_file,fmt_txt) 'use_fates = ', use_fates
       write(hist_fields_file,*)

       ! Table header
       ! Concatenate strings needed in format statement
       do i = 1, ncol
          fmt_txt = '('//str_width_col(i)//'a,x)'
          write(hist_fields_file,fmt_txt,advance='no') ('=', j=1,width_col(i))
       end do
       write(hist_fields_file,*)  ! next write statement will now appear in new line

       ! Table title
       fmt_txt = '(a)'
       write(hist_fields_file,fmt_txt) 'CTSM History Fields'

       ! Sub-header
       ! Concatenate strings needed in format statement
       fmt_txt = '('//str_w_col_sum//'a)'
       write(hist_fields_file,fmt_txt) ('-', i=1, width_col_sum)
       ! Concatenate strings needed in format statement
       fmt_txt = '(a'//str_width_col(1)//',x,a'//str_width_col(2)//',x,a'//str_width_col(3)//',x,a'//str_width_col(4)//',x,a'//str_width_col(5)//')'
       write(hist_fields_file,fmt_txt) 'Variable Name',  &
                           'Level Dim.', 'Long Description', 'Units', 'Active?'

       ! End header, same as header
       ! Concatenate strings needed in format statement
       do i = 1, ncol
          fmt_txt = '('//str_width_col(i)//'a,x)'
          write(hist_fields_file,fmt_txt,advance='no') ('=', j=1,width_col(i))
       end do
       write(hist_fields_file,*)  ! next write statement will now appear in new line

       ! Main table
       ! Concatenate strings needed in format statement
       fmt_txt = '(a'//str_width_col(1)//',x,a'//str_width_col(2)//',x,a'//str_width_col(3)//',x,a'//str_width_col(4)//',l'//str_width_col(5)//')'
       do nf = 1,nallhistflds
          write(hist_fields_file,fmt_txt) &
             allhistfldlist(nf)%field%name,  &
             allhistfldlist(nf)%field%type2d,  &
             allhistfldlist(nf)%field%long_name,  &
             allhistfldlist(nf)%field%units,  &
             allhistfldlist(nf)%actflag(1)
       end do

       ! Table footer, same as header
       ! Concatenate strings needed in format statement
       do i = 1, ncol
          fmt_txt = '('//str_width_col(i)//'a,x)'
          write(hist_fields_file,fmt_txt,advance='no') ('=', j=1,width_col(i))
       end do

       call shr_sys_flush(hist_fields_file)
       close(unit = hist_fields_file)
       call relavu(hist_fields_file)  ! close and release file unit number
    end if

  end subroutine hist_printflds

  !-----------------------------------------------------------------------
  subroutine allhistfldlist_addfld (fname, numdims, type1d, type1d_out, &
        type2d, num2d, units, avgflag, long_name, hpindex, &
        p2c_scale_type, c2l_scale_type, l2g_scale_type, &
        no_snow_behavior)
    !
    ! !DESCRIPTION:
    ! Add a field to the list of all history fields. Put input arguments of
    ! field name, units, number of levels, averaging flag, and long name
    ! into a type entry in the global list of all history fields (allhistfldlist).
    !
    ! The optional argument no_snow_behavior should be given when this is a multi-layer
    ! snow field, and should be absent otherwise. It should take on one of the no_snow_*
    ! parameters defined above
    !
    ! !ARGUMENTS:
    character(len=*), intent(in)  :: fname            ! field name
    integer         , intent(in)  :: numdims          ! number of dimensions
    character(len=*), intent(in)  :: type1d           ! 1d data type
    character(len=*), intent(in)  :: type1d_out       ! 1d output type
    character(len=*), intent(in)  :: type2d           ! 2d output type
    integer         , intent(in)  :: num2d            ! size of second dimension (e.g. number of vertical levels)
    character(len=*), intent(in)  :: units            ! units of field
    character(len=*), intent(in)  :: avgflag          ! time averaging flag
    character(len=*), intent(in)  :: long_name        ! long name of field
    integer         , intent(in)  :: hpindex          ! index into raw history data (clmptr_r*) arrays
    character(len=*), intent(in)  :: p2c_scale_type   ! scale type for subgrid averaging of pfts to column
    character(len=*), intent(in)  :: c2l_scale_type   ! scale type for subgrid averaging of columns to landunits
    character(len=*), intent(in)  :: l2g_scale_type   ! scale type for subgrid averaging of landunits to gridcells
    integer, intent(in), optional :: no_snow_behavior ! if a multi-layer snow field, behavior to use for absent snow layers
    !
    ! !LOCAL VARIABLES:
    integer :: n            ! loop index
    integer :: f            ! allhistfldlist index
    integer :: numa         ! total number of atm cells across all processors
    integer :: numg         ! total number of gridcells across all processors
    integer :: numl         ! total number of landunits across all processors
    integer :: numc         ! total number of columns across all processors
    integer :: nump         ! total number of pfts across all processors
    type(bounds_type) :: bounds
    character(len=*),parameter :: subname = 'allhistfldlist_addfld'
    !------------------------------------------------------------------------

    if (.not. avgflag_valid(avgflag, blank_valid=.true.)) then
       write(iulog,*) trim(subname),' ERROR: unknown averaging flag=', avgflag
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! Determine bounds

    call get_proc_bounds(bounds)
    call get_proc_global(ng=numg, nl=numl, nc=numc, np=nump)

    ! Ensure that new field is not all blanks

    if (fname == ' ') then
       write(iulog,*) trim(subname),' ERROR: blank field name not allowed'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! Ensure that new field name isn't too long

    if (len_trim(fname) > max_namlen ) then
       write(iulog,*) trim(subname),' ERROR: field name too long: ', trim(fname)
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    ! Ensure that new field doesn't already exist

    do n = 1,nallhistflds
       if (allhistfldlist(n)%field%name == fname) then
          write(iulog,*) trim(subname),' ERROR:', fname, ' already on list'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end do

    ! Increase number of fields on list of all history fields

    nallhistflds = nallhistflds + 1
    f = nallhistflds

    ! Check number of fields in list against maximum number

    if (nallhistflds > max_flds) then
       write(iulog,*) trim(subname),' ERROR: too many fields for primary history file ', &
            '-- max_flds,nallhistflds=', max_flds, nallhistflds
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! Add field to list of all history fields

    allhistfldlist(f)%field%name           = fname
    allhistfldlist(f)%field%long_name      = long_name
    allhistfldlist(f)%field%units          = units
    allhistfldlist(f)%field%type1d         = type1d
    allhistfldlist(f)%field%type1d_out     = type1d_out
    allhistfldlist(f)%field%type2d         = type2d
    allhistfldlist(f)%field%numdims        = numdims
    allhistfldlist(f)%field%num2d          = num2d
    allhistfldlist(f)%field%hpindex        = hpindex
    allhistfldlist(f)%field%p2c_scale_type = p2c_scale_type
    allhistfldlist(f)%field%c2l_scale_type = c2l_scale_type
    allhistfldlist(f)%field%l2g_scale_type = l2g_scale_type

    select case (type1d)
    case (grlnd)
       allhistfldlist(f)%field%beg1d = bounds%begg
       allhistfldlist(f)%field%end1d = bounds%endg
       allhistfldlist(f)%field%num1d = numg
    case (nameg)
       allhistfldlist(f)%field%beg1d = bounds%begg
       allhistfldlist(f)%field%end1d = bounds%endg
       allhistfldlist(f)%field%num1d = numg
    case (namel)
       allhistfldlist(f)%field%beg1d = bounds%begl
       allhistfldlist(f)%field%end1d = bounds%endl
       allhistfldlist(f)%field%num1d = numl
    case (namec)
       allhistfldlist(f)%field%beg1d = bounds%begc
       allhistfldlist(f)%field%end1d = bounds%endc
       allhistfldlist(f)%field%num1d = numc
    case (namep)
       allhistfldlist(f)%field%beg1d = bounds%begp
       allhistfldlist(f)%field%end1d = bounds%endp
       allhistfldlist(f)%field%num1d = nump
    case default
       write(iulog,*) trim(subname),' ERROR: unknown 1d output type= ',type1d
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end select

    if (present(no_snow_behavior)) then
       allhistfldlist(f)%field%no_snow_behavior = no_snow_behavior
    else
       allhistfldlist(f)%field%no_snow_behavior = no_snow_unset
    end if

    ! The following two fields are used only in list of all history fields,
    ! NOT in the runtime active field list
    ! ALL FIELDS IN THE FORMER ARE INITIALIZED WITH THE ACTIVE
    ! FLAG SET TO FALSE

    allhistfldlist(f)%avgflag(:) = avgflag
    allhistfldlist(f)%actflag(:) = .false.

  end subroutine allhistfldlist_addfld

  !-----------------------------------------------------------------------
  subroutine hist_htapes_build ()
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
    use clm_time_manager, only: get_prev_time
    use clm_varcon      , only: secspday
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: i                   ! index
    integer :: ier                 ! error code
    integer :: t, f                ! tape, field indices
    integer :: day, sec            ! day and seconds from base date
    character(len=*),parameter :: subname = 'hist_htapes_build'
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*)  trim(subname),' Initializing ', trim(compname), ' history files'
       write(iulog,'(72a1)') ("-",i=1,60)
       call shr_sys_flush(iulog)
    endif

    ! Define field list information for all history files.
    ! Update ntapes to reflect number of active history files
    ! Note - branch runs can have additional auxiliary history files
    ! declared).

    call htapes_fieldlist()

    ! Determine if gridcell (xy) averaging is done for all fields on tape

    do t=1,ntapes
       tape(t)%dov2xy = hist_dov2xy(t)
       if (masterproc) then
          write(iulog,*)trim(subname),' hist tape = ',t,&
               ' written with dov2xy= ',tape(t)%dov2xy
       end if
    end do

    ! Set number of time samples in each history file and
    ! Note - the following entries will be overwritten by history restart
    ! Note - with netcdf, only 1 (ncd_double) and 2 (ncd_float) are allowed

    do t=1,ntapes
       tape(t)%ntimes = 0
       tape(t)%dov2xy = hist_dov2xy(t)
       tape(t)%nhtfrq = hist_nhtfrq(t)
       tape(t)%mfilt = hist_mfilt(t)
       if (hist_ndens(t) == 1) then
          tape(t)%ncprec = ncd_double
       else
          tape(t)%ncprec = ncd_float
       endif
    end do

    ! Set time of beginning of current averaging interval
    ! First etermine elapased time since reference date

    call get_prev_time(day, sec)
    do t=1,ntapes
       tape(t)%begtime = day + sec/secspday
    end do

    if (masterproc) then
       write(iulog,*)  trim(subname),' Successfully initialized ', trim(compname), ' history files'
       write(iulog,'(72a1)') ("-",i=1,60)
       call shr_sys_flush(iulog)
    endif

  end subroutine hist_htapes_build

  !-----------------------------------------------------------------------
  subroutine allhistfldlist_make_active (name, tape_index, avgflag)
    !
    ! !DESCRIPTION:
    ! Add a field to the default ``on'' list for a given history file.
    ! Also change the default time averaging flag if requested.
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: name          ! field name
    integer, intent(in) :: tape_index             ! history tape index
    character(len=*), intent(in), optional :: avgflag  ! time averaging flag
    !
    ! !LOCAL VARIABLES:
    integer :: f            ! field index
    logical :: found        ! flag indicates field found in allhistfldlist
    character(len=*),parameter :: subname = 'allhistfldlist_make_active'
    !-----------------------------------------------------------------------

    ! Check validity of input arguments

    if (tape_index > max_tapes) then
       write(iulog,*) trim(subname),' ERROR: tape index=', tape_index, ' is too big'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    if (present(avgflag)) then
       if (.not. avgflag_valid(avgflag, blank_valid=.true.)) then
          write(iulog,*) trim(subname),' ERROR: unknown averaging flag=', avgflag
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
    end if

    ! Look through list of all history fields for input field name.
    ! When found, set active flag for that tape to true.
    ! Also reset averaging flag if told to use other than default.

    found = .false.
    do f = 1,nallhistflds
       if (trim(name) == trim(allhistfldlist(f)%field%name)) then
          allhistfldlist(f)%actflag(tape_index) = .true.
          if (present(avgflag)) then
             if (avgflag/= ' ') allhistfldlist(f)%avgflag(tape_index) = avgflag
          end if
          found = .true.
          exit
       end if
    end do
    if (.not. found) then
       write(iulog,*) trim(subname),' ERROR: field=', name, ' not found'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

  end subroutine allhistfldlist_make_active

  !-----------------------------------------------------------------------
  subroutine allhistfldlist_change_timeavg (t)
    !
    ! !DESCRIPTION:
    ! Override default history tape contents for a specific tape.
    ! Copy the flag into the list of all history fields.
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t         ! history tape index
    !
    ! !LOCAL VARIABLES:
    integer :: f                     ! field index
    character(len=avgflag_strlen) :: avgflag      ! local equiv of hist_avgflag_pertape(t)
    character(len=*),parameter :: subname = 'allhistfldlist_change_timeavg'
    !-----------------------------------------------------------------------

    avgflag = hist_avgflag_pertape(t)
    if (.not. avgflag_valid(avgflag, blank_valid = .false.)) then
       write(iulog,*) trim(subname),' ERROR: unknown avgflag=',avgflag
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    do f = 1,nallhistflds
       allhistfldlist(f)%avgflag(t) = avgflag
    end do

  end subroutine allhistfldlist_change_timeavg

  !-----------------------------------------------------------------------
  subroutine htapes_fieldlist()
    !
    ! !DESCRIPTION:
    ! Define the contents of each history file based on namelist
    ! input for initial or branch run, and restart data if a restart run.
    ! Fill and use arrays fincl and fexcl to modify default history tape contents.
    ! Then sort the result alphanumerically.
    !
    ! Sets history_tape_in_use and htapes_defined. Fills fields in 'tape' array.
    ! Optionally updates allhistfldlist avgflag.
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: t, f                         ! tape, field indices
    integer :: ff                           ! index into include, exclude and fprec list
    character(len=max_namlen) :: name       ! field name portion of fincl (i.e. no avgflag separator)
    character(len=max_namlen) :: allhistfldname ! name from allhistfldlist field
    character(len=avgflag_strlen) :: avgflag ! averaging flag
    character(len=1)  :: prec_acc           ! history buffer precision flag
    character(len=1)  :: prec_wrt           ! history buffer write precision flag
    character(len=*),parameter :: subname = 'htapes_fieldlist'
    !-----------------------------------------------------------------------

    ! Override averaging flag for all fields on a particular tape
    ! if namelist input so specifies

    do t=1,max_tapes
       if (hist_avgflag_pertape(t) /= ' ') then
          call allhistfldlist_change_timeavg (t)
       end if
    end do

    fincl(:,1)  = hist_fincl1(:)
    fincl(:,2)  = hist_fincl2(:)
    fincl(:,3)  = hist_fincl3(:)
    fincl(:,4)  = hist_fincl4(:)
    fincl(:,5)  = hist_fincl5(:)
    fincl(:,6)  = hist_fincl6(:)
    fincl(:,7)  = hist_fincl7(:)
    fincl(:,8)  = hist_fincl8(:)
    fincl(:,9)  = hist_fincl9(:)
    fincl(:,10) = hist_fincl10(:)

    fexcl(:,1)  = hist_fexcl1(:)
    fexcl(:,2)  = hist_fexcl2(:)
    fexcl(:,3)  = hist_fexcl3(:)
    fexcl(:,4)  = hist_fexcl4(:)
    fexcl(:,5)  = hist_fexcl5(:)
    fexcl(:,6)  = hist_fexcl6(:)
    fexcl(:,7)  = hist_fexcl7(:)
    fexcl(:,8)  = hist_fexcl8(:)
    fexcl(:,9)  = hist_fexcl9(:)
    fexcl(:,10) = hist_fexcl10(:)


    ! First ensure contents of fincl and fexcl are valid names

    do t = 1,max_tapes
       f = 1
       do while (f < max_flds .and. fincl(f,t) /= ' ')
          name = getname (fincl(f,t))
          do ff = 1,nallhistflds
             allhistfldname = allhistfldlist(ff)%field%name
             if (name == allhistfldname) exit
          end do
          if (name /= allhistfldname) then
             write(iulog,*) trim(subname),' ERROR: ', trim(name), ' in fincl(', f, ') ',&
                  'for history tape ',t,' not found'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
          f = f + 1
       end do

       f = 1
       do while (f < max_flds .and. fexcl(f,t) /= ' ')
          do ff = 1,nallhistflds
             allhistfldname = allhistfldlist(ff)%field%name
             if (fexcl(f,t) == allhistfldname) exit
          end do
          if (fexcl(f,t) /= allhistfldname) then
             write(iulog,*) trim(subname),' ERROR: ', fexcl(f,t), ' in fexcl(', f, ') ', &
                  'for history tape ',t,' not found'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
          f = f + 1
       end do
    end do

    history_tape_in_use(:) = .false.
    tape(:)%nflds = 0
    do t = 1,max_tapes

       ! Loop through the allhistfldlist set of field names and determine if any of those
       ! are in the FINCL or FEXCL arrays
       ! The call to list_index determines the index in the FINCL or FEXCL arrays
       ! that the allhistfldlist field corresponds to
       ! Add the field to the tape if specified via namelist (FINCL[1-max_tapes]),
       ! or if it is on by default and was not excluded via namelist (FEXCL[1-max_tapes]).

       do f = 1,nallhistflds
          allhistfldname = allhistfldlist(f)%field%name
          call list_index (fincl(1,t), allhistfldname, ff)

          if (ff > 0) then

             ! if field is in include list, ff > 0 and htape_addfld
             ! will be called for field

             avgflag = getflag (fincl(ff,t))
             call htape_addfld (t, f, avgflag)

          else if (.not. hist_empty_htapes) then

             ! find index of field in exclude list

             call list_index (fexcl(1,t), allhistfldname, ff)

             ! if field is in exclude list, ff > 0 and htape_addfld
             ! will not be called for field
             ! if field is not in exclude list, ff =0 and htape_addfld
             ! will be called for field (note that htape_addfld will be
             ! called below only if field is not in exclude list OR in
             ! include list

             if (ff == 0 .and. allhistfldlist(f)%actflag(t)) then
                call htape_addfld (t, f, ' ')
             end if

          end if
       end do

       ! Specification of tape contents now complete.
       ! Sort each list of active entries
       call sort_hist_list(t, tape(t)%nflds, tape(t)%hlist)

       if (masterproc) then
          if (tape(t)%nflds > 0) then
             write(iulog,*) trim(subname),' : Included fields tape ',t,'=',tape(t)%nflds
          end if
          do f = 1,tape(t)%nflds
             write(iulog,*) f,' ',tape(t)%hlist(f)%field%name, &
                  tape(t)%hlist(f)%field%num2d,' ',tape(t)%hlist(f)%avgflag
          end do
          call shr_sys_flush(iulog)
       end if
    end do

    ! Determine index of max active history tape, and whether each tape is in use

    ntapes = 0
    do t = max_tapes,1,-1
       if (tape(t)%nflds > 0) then
          ntapes = t
          exit
       end if
    end do

    do t = 1, ntapes
       if (tape(t)%nflds > 0) then
          history_tape_in_use(t) = .true.
       end if
    end do

    ! Change 1d output per tape output flag if requested - only for history
    ! tapes where 2d xy averaging is not enabled

    do t = 1,ntapes
       if (hist_type1d_pertape(t) /= ' ' .and. (.not. hist_dov2xy(t))) then
          select case (trim(hist_type1d_pertape(t)))
          case ('PFTS','COLS', 'LAND', 'GRID')
             if ( masterproc ) &
             write(iulog,*)'history tape ',t,' will have 1d output type of ',hist_type1d_pertape(t)
          case default
             write(iulog,*) trim(subname),' ERROR: unknown namelist type1d per tape=',hist_type1d_pertape(t)
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end select
       end if
    end do

    if (masterproc) then
       write(iulog,*) 'There will be a total of ',ntapes,' history tapes'
       do t=1,ntapes
          write(iulog,*)
          if (hist_nhtfrq(t) == 0) then
             write(iulog,*)'History tape ',t,' write frequency is MONTHLY'
          else
             write(iulog,*)'History tape ',t,' write frequency = ',hist_nhtfrq(t)
          endif
          if (hist_dov2xy(t)) then
             write(iulog,*)'All fields on history tape ',t,' are grid averaged'
          else
             write(iulog,*)'All fields on history tape ',t,' are not grid averaged'
          end if
          write(iulog,*)'Number of time samples on history tape ',t,' is ',hist_mfilt(t)
          write(iulog,*)'Output precision on history tape ',t,'=',hist_ndens(t)
          if (.not. history_tape_in_use(t)) then
             write(iulog,*) 'History tape ',t,' does not have any fields,'
             write(iulog,*) 'so it will not be written!'
          end if
          write(iulog,*)
       end do
       call shr_sys_flush(iulog)
    end if

    ! Set flag indicating h-tape contents are now defined (needed by allhistfldlist_addfld)

    htapes_defined = .true.


  end subroutine htapes_fieldlist

  !-----------------------------------------------------------------------
  subroutine copy_allhistfldlist_entry(this, other)
    ! set this = other
    class(allhistfldlist_entry), intent(out) :: this
    class(entry_base), intent(in) :: other

    select type(this)
    type is (allhistfldlist_entry)
       select type(other)
       type is (allhistfldlist_entry)
          this = other
       class default
          call endrun('Unexpected type of "other" in copy_allhistfldlist_entry')
       end select
    class default
       call endrun('Unexpected type of "this" in copy_allhistfldlist_entry')
    end select
  end subroutine copy_allhistfldlist_entry

  !-----------------------------------------------------------------------
  subroutine copy_history_entry(this, other)
    ! set this = other
    class(history_entry), intent(out) :: this
    class(entry_base), intent(in) :: other

    select type(this)
    type is (history_entry)
       select type(other)
       type is (history_entry)
          this = other
       class default
          call endrun('Unexpected type of "other" in copy_history_entry')
       end select
    class default
       call endrun('Unexpected type of "this" in copy_history_entry')
    end select
  end subroutine copy_history_entry

  !-----------------------------------------------------------------------
  subroutine sort_hist_list(t, n_fields, hist_list)

    ! !DESCRIPTION:
    ! Sort list of history variable names hist_list in alphabetical
    ! order.

    ! !ARGUMENTS:
    integer, intent(in) :: t  ! tape index
    integer, intent(in) :: n_fields  ! number of fields
    class(entry_base), intent(inout) :: hist_list(:)

    ! !LOCAL VARIABLES:
    integer :: f, ff  ! field indices
    class(entry_base), allocatable :: tmp

    character(len=*), parameter :: subname = 'sort_hist_list'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL(size(hist_list) >= n_fields, sourcefile, __LINE__)

    if (n_fields < 2) then
       return
    end if

    allocate(tmp, source = hist_list(1))

    do f = n_fields-1, 1, -1
       do ff = 1, f
          ! First sort by the name of the level dimension; then, within the list of
          ! fields with the same level dimension, sort by field name. Sorting first by
          ! the level dimension gives a significant performance improvement especially
          ! notable on lustre file systems such as on derecho.
          if (hist_list(ff)%field%type2d > hist_list(ff+1)%field%type2d .or. &
               (hist_list(ff)%field%type2d == hist_list(ff+1)%field%type2d .and. &
               hist_list(ff)%field%name > hist_list(ff+1)%field%name)) then

             call tmp%copy(hist_list(ff))
             call hist_list(ff  )%copy(hist_list(ff+1))
             call hist_list(ff+1)%copy(tmp)

          end if
       end do
    end do

  end subroutine sort_hist_list

  !-----------------------------------------------------------------------
  logical function is_mapping_upto_subgrid( type1d, type1d_out ) result ( mapping)
    !
    ! !DESCRIPTION:
    !
    ! Return true if this field will be mapped into a higher subgrid level
    ! If false it will be output on it's native grid
    !
    ! !ARGUMENTS:
    implicit none
    character(len=8), intent(in) :: type1d      ! clm pointer 1d type
    character(len=8), intent(in) :: type1d_out  ! history buffer 1d type
    !
    mapping = .false.
    if (type1d_out == nameg .or. type1d_out == grlnd) then
       if (type1d == namep) then
          mapping = .true.
       else if (type1d == namec) then
          mapping = .true.
       else if (type1d == namel) then
          mapping = .true.
       end if
    else if (type1d_out == namel ) then
       if (type1d == namep) then
          mapping = .true.
       else if (type1d == namec) then
          mapping = .true.
       end if
    else if (type1d_out == namec ) then
       if (type1d == namep) then
          mapping = .true.
       end if
    end if
  end function is_mapping_upto_subgrid

  !-----------------------------------------------------------------------
  subroutine htape_addfld (t, f, avgflag)
    !
    ! !DESCRIPTION:
    ! Add a field to a history tape, copying metadata from the list of all history fields
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t                 ! history tape index
    integer, intent(in) :: f                 ! field index from list of all history fields
    character(len=*), intent(in) :: avgflag  ! time averaging flag
    !
    ! !LOCAL VARIABLES:
    integer :: n                    ! field index on defined tape
    character(len=hist_dim_name_length) :: type1d      ! clm pointer 1d type
    character(len=hist_dim_name_length) :: type1d_out  ! history buffer 1d type
    integer :: numa                 ! total number of atm cells across all processors
    integer :: numg                 ! total number of gridcells across all processors
    integer :: numl                 ! total number of landunits across all processors
    integer :: numc                 ! total number of columns across all processors
    integer :: nump                 ! total number of pfts across all processors
    integer :: num2d                ! size of second dimension (e.g. .number of vertical levels)
    integer :: beg1d_out,end1d_out  ! history output per-proc 1d beginning and ending indices
    integer :: beg1d,end1d          ! beginning and ending indices for this field (assume already set)
    integer :: num1d_out            ! history output 1d size
    type(bounds_type) :: bounds
    character(len=avgflag_strlen) :: avgflag_temp  ! local copy of hist_avgflag_pertape(t)
    character(len=*),parameter :: subname = 'htape_addfld'
    !-----------------------------------------------------------------------

    ! Ensure that it is not too late to add a field to the history tape

    if (htapes_defined) then
       write(iulog,*) trim(subname),' ERROR: attempt to add field ', &
            allhistfldlist(f)%field%name, ' after history files are set'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    tape(t)%nflds = tape(t)%nflds + 1
    n = tape(t)%nflds

    ! Copy field information

    tape(t)%hlist(n)%field = allhistfldlist(f)%field

    ! Determine bounds

    call get_proc_bounds(bounds)
    call get_proc_global(ng=numg, nl=numl, nc=numc, np=nump)

    ! Modify type1d_out if necessary

    if (hist_dov2xy(t)) then

       ! If xy output averaging is requested, set output 1d type to grlnd
       ! ***NOTE- the following logic is what permits non lat/lon grids to
       ! be written to clm history file

       type1d = tape(t)%hlist(n)%field%type1d

       if (type1d == nameg .or. &
           type1d == namel .or. &
           type1d == namec .or. &
           type1d == namep) then
          tape(t)%hlist(n)%field%type1d_out = grlnd
       end if
       if (type1d == grlnd) then
          tape(t)%hlist(n)%field%type1d_out = grlnd
       end if

    else if (hist_type1d_pertape(t) /= ' ') then

       ! Set output 1d type  based on namelist setting of  hist_type1d_pertape
       ! Only applies to tapes when xy output is not required

       type1d = tape(t)%hlist(n)%field%type1d

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
          write(iulog,*) trim(subname),' ERROR: unknown input hist_type1d_pertape= ', hist_type1d_pertape(t)
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select

    endif

    ! Determine output 1d dimensions

    type1d_out = tape(t)%hlist(n)%field%type1d_out
    if (type1d_out == grlnd) then
       beg1d_out = bounds%begg
       end1d_out = bounds%endg
       num1d_out = numg
    else if (type1d_out == nameg) then
       beg1d_out = bounds%begg
       end1d_out = bounds%endg
       num1d_out = numg
    else if (type1d_out == namel) then
       beg1d_out = bounds%begl
       end1d_out = bounds%endl
       num1d_out = numl
    else if (type1d_out == namec) then
       beg1d_out = bounds%begc
       end1d_out = bounds%endc
       num1d_out = numc
    else if (type1d_out == namep) then
       beg1d_out = bounds%begp
       end1d_out = bounds%endp
       num1d_out = nump
    else
       write(iulog,*) trim(subname),' ERROR: incorrect value of type1d_out= ',type1d_out
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! Output bounds for the field
    tape(t)%hlist(n)%field%beg1d_out = beg1d_out
    tape(t)%hlist(n)%field%end1d_out = end1d_out
    tape(t)%hlist(n)%field%num1d_out = num1d_out

    ! Fields native bounds
    beg1d = allhistfldlist(f)%field%beg1d
    end1d = allhistfldlist(f)%field%end1d

    ! Alloccate and initialize history buffer and related info

    num2d = tape(t)%hlist(n)%field%num2d
    if ( is_mapping_upto_subgrid( type1d, type1d_out ) ) then
       allocate (tape(t)%hlist(n)%hbuf(beg1d_out:end1d_out,num2d))
       allocate (tape(t)%hlist(n)%nacs(beg1d_out:end1d_out,num2d))
    else
       allocate (tape(t)%hlist(n)%hbuf(beg1d:end1d,num2d))
       allocate (tape(t)%hlist(n)%nacs(beg1d:end1d,num2d))
    end if
    tape(t)%hlist(n)%hbuf(:,:) = 0._r8
    tape(t)%hlist(n)%nacs(:,:) = 0

    ! Set time averaging flag based on allhistfldlist setting or
    ! override the default averaging flag with namelist setting

    if (.not. avgflag_valid(avgflag, blank_valid=.true.)) then
       write(iulog,*) trim(subname),' ERROR: unknown avgflag=', avgflag
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    if (avgflag == ' ') then
       tape(t)%hlist(n)%avgflag = allhistfldlist(f)%avgflag(t)
    else
       tape(t)%hlist(n)%avgflag = avgflag
    end if

    ! Override this tape's avgflag if nhtfrq == 1
    if (tape(t)%nhtfrq == 1) then  ! output is instantaneous
       hist_avgflag_pertape(t) = 'I'
    end if
    ! Override this field's avgflag if the namelist or the previous line
    ! has set this tape to
    ! - instantaneous (I) or
    ! - local time (L)
    avgflag_temp = hist_avgflag_pertape(t)
    if (avgflag_temp == 'I' .or. avgflag_temp(1:1) == 'L') then
       tape(t)%hlist(n)%avgflag = avgflag_temp
    end if

  end subroutine htape_addfld

  !-----------------------------------------------------------------------
  subroutine hist_update_hbuf(bounds)
    !
    ! !DESCRIPTION:
    ! Accumulate (or take min, max, etc. as appropriate) input field
    ! into its history buffer for appropriate tapes.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: t                   ! tape index
    integer :: f                   ! field index
    integer :: num2d               ! size of second dimension (e.g. number of vertical levels)
    integer :: numdims             ! number of dimensions
    character(len=*),parameter :: subname = 'hist_update_hbuf'
    character(len=hist_dim_name_length) :: type2d     ! hbuf second dimension type ["levgrnd","levlak","numrad","ltype","natpft","cft","glc_nec","elevclas","subname(n)","mxsowings","mxharvests"]
    !-----------------------------------------------------------------------

    do t = 1,ntapes
!$OMP PARALLEL DO PRIVATE (f, num2d, numdims)
       do f = 1,tape(t)%nflds

          numdims = tape(t)%hlist(f)%field%numdims

          if ( numdims == 1) then
             call hist_update_hbuf_field_1d (t, f, bounds)
          else
             num2d = tape(t)%hlist(f)%field%num2d
             call hist_update_hbuf_field_2d (t, f, bounds, num2d)
          end if
       end do
!$OMP END PARALLEL DO
    end do

  end subroutine hist_update_hbuf

  !-----------------------------------------------------------------------
  subroutine hist_update_hbuf_field_1d (t, f, bounds)
    !
    ! !DESCRIPTION:
    ! Accumulate (or take min, max, etc. as appropriate) input field
    ! into its history buffer for appropriate tapes.
    !
    ! This canNOT be called from within a threaded region (see comment below regarding the
    ! call to p2g, and the lack of explicit bounds on its arguments; see also bug 1786)
    !
    ! !USES:
    use subgridAveMod   , only : p2g, c2g, l2g, p2l, c2l, p2c
    use decompMod       , only : bounds_level_proc
    use clm_varcon      , only : degpsec, isecspday
    use clm_time_manager, only : get_curr_date
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t            ! tape index
    integer, intent(in) :: f            ! field index
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: hpindex                 ! index into raw history data (clmptr_r*) arrays
    integer  :: k                       ! gridcell, landunit, column or patch index
    integer  :: beg1d,end1d             ! beginning and ending indices
    integer  :: beg1d_out,end1d_out     ! beginning and ending indices on output grid
    logical  :: check_active            ! true => check 'active' flag of each point (this refers to a point being active, NOT a history field being active)
    logical  :: valid                   ! true => history operation is valid
    logical  :: map2gcell               ! true => map clm pointer field to gridcell
    character(len=hist_dim_name_length)  :: type1d         ! 1d clm pointerr type   ["gridcell","landunit","column","pft"]
    character(len=hist_dim_name_length)  :: type1d_out     ! 1d history buffer type ["gridcell","landunit","column","pft"]
    character(len=avgflag_strlen) :: avgflag ! time averaging flag
    character(len=scale_type_strlen)  :: p2c_scale_type ! scale type for subgrid averaging of pfts to column
    character(len=scale_type_strlen)  :: c2l_scale_type ! scale type for subgrid averaging of columns to landunits
    character(len=scale_type_strlen) :: l2g_scale_type ! scale type for subgrid averaging of landunits to gridcells
    real(r8), pointer :: hbuf(:,:)      ! history buffer
    integer , pointer :: nacs(:,:)      ! accumulation counter
    real(r8), pointer :: field(:)       ! clm 1d pointer field
    logical , pointer :: active(:)      ! flag saying whether each point is active (used for type1d = landunit/column/pft) (this refers to a point being active, NOT a history field being active)
    real(r8), allocatable :: field_gcell(:)  ! gricell level field (used if mapping to gridcell is done)
    integer j
    character(len=*),parameter :: subname = 'hist_update_hbuf_field_1d'
    integer k_offset                    ! offset for mapping sliced subarray pointers when outputting variables in PFT/col vector form
    integer :: year                      ! year (0, ...) for nstep
    integer :: month                     ! month (1, ..., 12) for nstep
    integer :: day                       ! day of month (1, ..., 31) for nstep
    integer :: secs                      ! seconds into current date for nstep
    integer :: local_secpl               ! seconds into current date in local time
    integer :: tod                       ! Desired local solar time of output in seconds
    integer :: weight                    ! Weight for linear interpolation in time for local time avgflag
    integer, allocatable :: grid_index(:)             ! Grid cell index for longitude
    integer, allocatable :: tods(:)
    character(len=1) :: avgflag_trim     ! first character of avgflag

    !-----------------------------------------------------------------------

    SHR_ASSERT_FL(bounds%level == bounds_level_proc, sourcefile, __LINE__)

    avgflag        =  tape(t)%hlist(f)%avgflag
    nacs           => tape(t)%hlist(f)%nacs
    hbuf           => tape(t)%hlist(f)%hbuf
    beg1d          =  tape(t)%hlist(f)%field%beg1d
    end1d          =  tape(t)%hlist(f)%field%end1d
    beg1d_out      =  tape(t)%hlist(f)%field%beg1d_out
    end1d_out      =  tape(t)%hlist(f)%field%end1d_out
    type1d         =  tape(t)%hlist(f)%field%type1d
    type1d_out     =  tape(t)%hlist(f)%field%type1d_out
    p2c_scale_type =  tape(t)%hlist(f)%field%p2c_scale_type
    c2l_scale_type =  tape(t)%hlist(f)%field%c2l_scale_type
    l2g_scale_type =  tape(t)%hlist(f)%field%l2g_scale_type
    hpindex        =  tape(t)%hlist(f)%field%hpindex
    field          => clmptr_rs(hpindex)%ptr

    call get_curr_date (year, month, day, secs)

    ! set variables to check weights when allocate all pfts

    map2gcell = .false.
    if (type1d_out == nameg .or. type1d_out == grlnd) then
       SHR_ASSERT_FL(beg1d_out == bounds%begg, sourcefile, __LINE__)
       SHR_ASSERT_FL(end1d_out == bounds%endg, sourcefile, __LINE__)
       if (type1d == namep) then
          ! In this and the following calls, we do NOT explicitly subset field using
          ! bounds (e.g., we do NOT do field(bounds%begp:bounds%endp). This is because,
          ! for some fields, the lower bound has been reset to 1 due to taking a pointer
          ! to an array slice. Thus, this code will NOT work properly if done within a
          ! threaded region! (See also bug 1786)
          allocate( field_gcell(beg1d_out:end1d_out) )
          call p2g(bounds, &
               field, &
               field_gcell(bounds%begg:bounds%endg), &
               p2c_scale_type, c2l_scale_type, l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namec) then
          allocate( field_gcell(beg1d_out:end1d_out) )
          call c2g(bounds, &
               field, &
               field_gcell(bounds%begg:bounds%endg), &
               c2l_scale_type, l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namel) then
          allocate( field_gcell(beg1d_out:end1d_out) )
          call l2g(bounds, &
               field, &
               field_gcell(bounds%begg:bounds%endg), &
               l2g_scale_type)
          map2gcell = .true.
       end if
    end if
    if (type1d_out == namel ) then
       SHR_ASSERT_FL(beg1d_out == bounds%begl, sourcefile, __LINE__)
       SHR_ASSERT_FL(end1d_out == bounds%endl, sourcefile, __LINE__)
       if (type1d == namep) then
          ! In this and the following calls, we do NOT explicitly subset field using
          ! bounds (e.g., we do NOT do field(bounds%begp:bounds%endp). This is because,
          ! for some fields, the lower bound has been reset to 1 due to taking a pointer
          ! to an array slice. Thus, this code will NOT work properly if done within a
          ! threaded region! (See also bug 1786)
          allocate( field_gcell(beg1d_out:end1d_out) )
          call p2l(bounds, &
               field, &
               field_gcell(beg1d_out:end1d_out), &
               p2c_scale_type, c2l_scale_type)
          map2gcell = .true.
       else if (type1d == namec) then
          allocate( field_gcell(beg1d_out:end1d_out) )
          call c2l(bounds, &
               field, &
               field_gcell(beg1d_out:end1d_out), &
               c2l_scale_type)
          map2gcell = .true.
       end if
    end if
    if (type1d_out == namec ) then
       SHR_ASSERT_FL(beg1d_out == bounds%begc, sourcefile, __LINE__)
       SHR_ASSERT_FL(end1d_out == bounds%endc, sourcefile, __LINE__)
       if (type1d == namep) then
          ! In this and the following calls, we do NOT explicitly subset field using
          ! bounds (e.g., we do NOT do field(bounds%begp:bounds%endp). This is because,
          ! for some fields, the lower bound has been reset to 1 due to taking a pointer
          ! to an array slice. Thus, this code will NOT work properly if done within a
          ! threaded region! (See also bug 1786)
          allocate( field_gcell(beg1d_out:end1d_out) )
          call p2c(bounds, &
               field, &
               field_gcell(beg1d_out:end1d_out), &
               p2c_scale_type)
          map2gcell = .true.
       end if
    end if
    if ( map2gcell .and. .not. is_mapping_upto_subgrid(type1d, type1d_out) )then
       call endrun(msg=trim(subname)//' ERROR: mapping upto subgrid level is inconsistent'//errMsg(sourcefile, __LINE__))
    end if
    if ( .not. map2gcell .and. is_mapping_upto_subgrid(type1d, type1d_out) )then
       call endrun(msg=trim(subname)//' ERROR: mapping upto subgrid level is inconsistent'//errMsg(sourcefile, __LINE__))
    end if

    if (map2gcell) then  ! Map to gridcell

       ! note that in this case beg1d = begg and end1d=endg
       avgflag_trim = avgflag(1:1)
       select case (avgflag_trim)
       case ('I') ! Instantaneous
          do k = beg1d_out, end1d_out
             if (field_gcell(k) /= spval) then
                hbuf(k,1) = field_gcell(k)
             else
                hbuf(k,1) = spval
             end if
             nacs(k,1) = 1
          end do
       case ('A', 'S') ! Time average / sum
          do k = beg1d_out, end1d_out
             if (field_gcell(k) /= spval) then
                if (nacs(k,1) == 0) hbuf(k,1) = 0._r8
                hbuf(k,1) = hbuf(k,1) + field_gcell(k)
                nacs(k,1) = nacs(k,1) + 1
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
          end do
       case ('X') ! Maximum over time
          do k = beg1d_out, end1d_out
             if (field_gcell(k) /= spval) then
                if (nacs(k,1) == 0) hbuf(k,1) = -1.e50_r8
                hbuf(k,1) = max( hbuf(k,1), field_gcell(k) )
             else
                hbuf(k,1) = spval
             endif
             nacs(k,1) = 1
          end do
       case ('M') ! Minimum over time
          do k = beg1d_out, end1d_out
             if (field_gcell(k) /= spval) then
                if (nacs(k,1) == 0) hbuf(k,1) = +1.e50_r8
                hbuf(k,1) = min( hbuf(k,1), field_gcell(k) )
             else
                hbuf(k,1) = spval
             endif
             nacs(k,1) = 1
          end do
       case ('L') ! Local solar time
          read(avgflag(2:6), *) tod
          do k = beg1d_out, end1d_out
             if (field_gcell(k) /= spval) then
                local_secpl = secs + grc%londeg(k)/degpsec
                local_secpl = mod(local_secpl,isecspday)
                weight = calc_weight_local_time(local_secpl, tod)
                if (weight > 0) then
                   if (nacs(k,1) == 0) hbuf(k,1) = 0._r8
                   hbuf(k,1) = hbuf(k,1) + field_gcell(k)*real(weight)
                   nacs(k,1) = nacs(k,1) + weight
                end if
              else
                 if (nacs(k,1) == 0) hbuf(k,1) = spval
              end if
          end do
       case default
          write(iulog,*) trim(subname),' ERROR: invalid time averaging flag ', avgflag
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select
       deallocate( field_gcell )

    else  ! Do not map to gridcell

       allocate( grid_index(beg1d:end1d) )

       ! For data defined on the pft, col or landunit, we need to check if a point is active
       ! to determine whether that point should be assigned spval
       if (type1d == namep) then
          check_active = .true.
          active => patch%active
          grid_index = patch%gridcell
       else if (type1d == namec) then
          check_active = .true.
          active => col%active
          grid_index = col%gridcell
       else if (type1d == namel) then
          check_active = .true.
          active =>lun%active
          grid_index = lun%gridcell
       else
          check_active = .false.
       end if

       avgflag_trim = avgflag(1:1)

       select case (avgflag_trim)
       case ('I') ! Instantaneous
          do k = beg1d,end1d
             valid = .true.
             if (check_active) then
                if (.not. active(k)) valid = .false.
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
       case ('A', 'S') ! Time average / sum
          ! create mappings for array slice pointers (which go from 1 to size(field) rather than beg1d to end1d)
          if ( end1d .eq. ubound(field,1) ) then
             k_offset = 0
          else
             k_offset = 1 - beg1d
          endif
          do k = beg1d,end1d
             valid = .true.
             if (check_active) then
                if (.not. active(k)) valid = .false.
             end if
             if (valid) then
                if (field(k+k_offset) /= spval) then   ! add k_offset
                   if (nacs(k,1) == 0) hbuf(k,1) = 0._r8
                   hbuf(k,1) = hbuf(k,1) + field(k+k_offset)   ! add k_offset
                   nacs(k,1) = nacs(k,1) + 1
                else
                   if (nacs(k,1) == 0) hbuf(k,1) = spval
                end if
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
          end do
       case ('X') ! Maximum over time
          do k = beg1d,end1d
             valid = .true.
             if (check_active) then
                if (.not. active(k)) valid = .false.
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
          do k = beg1d,end1d
             valid = .true.
             if (check_active) then
                if (.not. active(k)) valid = .false.
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
       case ('L') ! Local solar time
          read(avgflag(2:6), *) tod
          if ( end1d .eq. ubound(field,1) ) then
             k_offset = 0
          else
             k_offset = 1 - beg1d
          endif
             do k = beg1d, end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) then
                      valid = .false.
                   else
                      local_secpl = secs + grc%londeg(grid_index(k))/degpsec
                   end if
                else
                   local_secpl = secs + grc%londeg(k)/degpsec
                end if
                local_secpl = mod(local_secpl,isecspday)
                if (valid) then
                   weight = calc_weight_local_time(local_secpl, tod)
                   if (weight > 0 .and. field(k+k_offset) /= spval) then
                      if (nacs(k,1) == 0) hbuf(k,1) = 0._r8
                      hbuf(k,1) = hbuf(k,1) + field(k+k_offset)*real(weight)
                      nacs(k,1) = nacs(k,1) + weight
                   end if
                else
                   if (nacs(k,1) == 0) hbuf(k,1) = spval
                end if
             end do
       case default
          write(iulog,*) trim(subname),' ERROR: invalid time averaging flag ', avgflag
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select
    end if

  end subroutine hist_update_hbuf_field_1d

  !-----------------------------------------------------------------------
  subroutine hist_update_hbuf_field_2d (t, f, bounds, num2d)
    !
    ! !DESCRIPTION:
    ! Accumulate (or take min, max, etc. as appropriate) input field
    ! into its history buffer for appropriate tapes.
    !
    ! This canNOT be called from within a threaded region (see comment below regarding the
    ! call to p2g, and the lack of explicit bounds on its arguments; see also bug 1786)
    !
    ! !USES:
    use subgridAveMod   , only : p2g, c2g, l2g, p2l, c2l, p2c
    use decompMod       , only : bounds_level_proc
    use clm_varctl      , only : iulog
    use clm_varcon      , only : degpsec, isecspday
    use clm_time_manager, only : get_curr_date
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t            ! tape index
    integer, intent(in) :: f            ! field index
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num2d        ! size of second dimension
    !
    ! !LOCAL VARIABLES:
    integer  :: hpindex                 ! index into raw history data (clmptr_r*) arrays
    integer  :: k                       ! gridcell, landunit, column or patch index
    integer  :: j                       ! level index
    integer  :: beg1d,end1d             ! beginning and ending indices
    integer  :: beg1d_out,end1d_out     ! beginning and ending indices for output level
    logical  :: check_active            ! true => check 'active' flag of each point (this refers to a point being active, NOT a history field being active)
    logical  :: valid                   ! true => history operation is valid
    logical  :: map2gcell               ! true => map clm pointer field to gridcell
    character(len=hist_dim_name_length)  :: type1d         ! 1d clm pointerr type   ["gridcell","landunit","column","pft"]
    character(len=hist_dim_name_length)  :: type1d_out     ! 1d history buffer type ["gridcell","landunit","column","pft"]
    character(len=avgflag_strlen) :: avgflag ! time averaging flag
    character(len=scale_type_strlen) :: p2c_scale_type ! scale type for subgrid averaging of pfts to column
    character(len=scale_type_strlen) :: c2l_scale_type ! scale type for subgrid averaging of columns to landunits
    character(len=scale_type_strlen) :: l2g_scale_type ! scale type for subgrid averaging of landunits to gridcells
    integer  :: no_snow_behavior        ! for multi-layer snow fields, behavior to use when a given layer is absent
    real(r8), pointer :: hbuf(:,:)      ! history buffer
    integer , pointer :: nacs(:,:)      ! accumulation counter
    real(r8), pointer :: field(:,:)     ! clm 2d pointer field
    logical           :: field_allocated! whether 'field' was allocated here
    logical , pointer :: active(:)      ! flag saying whether each point is active (used for type1d = landunit/column/pft)
                                        !(this refers to a point being active, NOT a history field being active)
    real(r8), allocatable :: field_gcell(:,:) ! gridcell level field (used if mapping to gridcell is done)
    character(len=*),parameter :: subname = 'hist_update_hbuf_field_2d'
    integer :: year                      ! year (0, ...) for nstep
    integer :: month                     ! month (1, ..., 12) for nstep
    integer :: day                       ! day of month (1, ..., 31) for nstep
    integer :: secs                      ! seconds into current date for nstep
    integer :: local_secpl               ! seconds into current date in local time
    integer :: tod                       ! Desired local solar time of output in seconds
    integer :: weight                    ! Weight for linear interpolation in time for local time avgflag
    integer, allocatable :: grid_index(:)             ! Grid cell index for longitude
    integer, allocatable :: tods(:)
    character(len=1) :: avgflag_trim     ! first character of avgflag

    !-----------------------------------------------------------------------

    SHR_ASSERT_FL(bounds%level == bounds_level_proc, sourcefile, __LINE__)

    avgflag             =  tape(t)%hlist(f)%avgflag
    nacs                => tape(t)%hlist(f)%nacs
    hbuf                => tape(t)%hlist(f)%hbuf
    beg1d               =  tape(t)%hlist(f)%field%beg1d
    end1d               =  tape(t)%hlist(f)%field%end1d
    beg1d_out           =  tape(t)%hlist(f)%field%beg1d_out
    end1d_out           =  tape(t)%hlist(f)%field%end1d_out
    type1d              =  tape(t)%hlist(f)%field%type1d
    type1d_out          =  tape(t)%hlist(f)%field%type1d_out
    p2c_scale_type      =  tape(t)%hlist(f)%field%p2c_scale_type
    c2l_scale_type      =  tape(t)%hlist(f)%field%c2l_scale_type
    l2g_scale_type      =  tape(t)%hlist(f)%field%l2g_scale_type
    no_snow_behavior    =  tape(t)%hlist(f)%field%no_snow_behavior
    hpindex             =  tape(t)%hlist(f)%field%hpindex

    call get_curr_date (year, month, day, secs)

    if (no_snow_behavior /= no_snow_unset) then
       ! For multi-layer snow fields, build a special output variable that handles
       ! missing snow layers appropriately

       ! Note, regarding bug 1786: The following allocation is not what we would want if
       ! this routine were operating in a threaded region (or, more generally, within a
       ! loop over nclumps) - in that case we would want to use the bounds information for
       ! this clump. But currently that's not possible because the bounds of some fields
       ! have been reset to 1 - see also bug 1786. Similarly, if we wanted to allow
       ! operation within a loop over clumps, we would need to pass 'bounds' to
       ! hist_set_snow_field_2d rather than relying on beg1d & end1d (which give the proc,
       ! bounds not the clump bounds)

       allocate(field(lbound(clmptr_ra(hpindex)%ptr, 1) : ubound(clmptr_ra(hpindex)%ptr, 1), 1:num2d))
       field_allocated = .true.

       call hist_set_snow_field_2d(field, clmptr_ra(hpindex)%ptr, no_snow_behavior, type1d, &
            beg1d, end1d)
    else

       field => clmptr_ra(hpindex)%ptr(:,1:num2d)
       field_allocated = .false.
    end if

    ! set variables to check weights when allocate all pfts

    map2gcell = .false.
    if (type1d_out == nameg .or. type1d_out == grlnd) then
       SHR_ASSERT_FL(beg1d_out == bounds%begg, sourcefile, __LINE__)
       SHR_ASSERT_FL(end1d_out == bounds%endg, sourcefile, __LINE__)
       if (type1d == namep) then
          ! In this and the following calls, we do NOT explicitly subset field using
          ! (e.g., we do NOT do field(bounds%begp:bounds%endp). This is because,
          ! for some fields, the lower bound has been reset to 1 due to taking a pointer
          ! to an array slice. Thus, this code will NOT work properly if done within a
          ! threaded region! (See also bug 1786)
          allocate(field_gcell(bounds%begg:bounds%endg,num2d) )
          call p2g(bounds, num2d, &
               field, &
               field_gcell(bounds%begg:bounds%endg, :), &
               p2c_scale_type, c2l_scale_type, l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namec) then
          allocate(field_gcell(bounds%begg:bounds%endg,num2d) )
          call c2g(bounds, num2d, &
               field, &
               field_gcell(bounds%begg:bounds%endg, :), &
               c2l_scale_type, l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namel) then
          allocate(field_gcell(bounds%begg:bounds%endg,num2d) )
          call l2g(bounds, num2d, &
               field, &
               field_gcell(bounds%begg:bounds%endg, :), &
               l2g_scale_type)
          map2gcell = .true.
       end if
    else if ( type1d_out == namel )then
       SHR_ASSERT_FL(beg1d_out == bounds%begl, sourcefile, __LINE__)
       SHR_ASSERT_FL(end1d_out == bounds%endl, sourcefile, __LINE__)
       if (type1d == namep) then
          ! In this and the following calls, we do NOT explicitly subset field using
          ! (e.g., we do NOT do field(bounds%begp:bounds%endp). This is because,
          ! for some fields, the lower bound has been reset to 1 due to taking a pointer
          ! to an array slice. Thus, this code will NOT work properly if done within a
          ! threaded region! (See also bug 1786)
          allocate(field_gcell(beg1d_out:end1d_out,num2d))
          call p2l(bounds, num2d, &
               field, &
               field_gcell(beg1d_out:end1d_out, :), &
               p2c_scale_type, c2l_scale_type)
          map2gcell = .true.
       else if (type1d == namec) then
          allocate(field_gcell(beg1d_out:end1d_out,num2d))
          call c2l(bounds, num2d, &
               field, &
               field_gcell(beg1d_out:end1d_out, :), &
               c2l_scale_type)
          map2gcell = .true.
       end if
    else if ( type1d_out == namec )then
       SHR_ASSERT_FL(beg1d_out == bounds%begc, sourcefile, __LINE__)
       SHR_ASSERT_FL(end1d_out == bounds%endc, sourcefile, __LINE__)
       if (type1d == namep) then
          ! In this and the following calls, we do NOT explicitly subset field using
          ! (e.g., we do NOT do field(bounds%begp:bounds%endp). This is because,
          ! for some fields, the lower bound has been reset to 1 due to taking a pointer
          ! to an array slice. Thus, this code will NOT work properly if done within a
          ! threaded region! (See also bug 1786)
          allocate(field_gcell(beg1d_out:end1d_out,num2d))
          call p2c(bounds, num2d, &
               field, &
               field_gcell(beg1d_out:end1d_out, :), &
               p2c_scale_type)
          map2gcell = .true.
       end if
    end if
    if ( map2gcell .and. .not. is_mapping_upto_subgrid(type1d, type1d_out) )then
       call endrun(msg=trim(subname)//' ERROR: mapping upto subgrid level is inconsistent'//errMsg(sourcefile, __LINE__))
    end if
    if ( .not. map2gcell .and. is_mapping_upto_subgrid(type1d, type1d_out) )then
       call endrun(msg=trim(subname)//' ERROR: mapping upto subgrid level is inconsistent'//errMsg(sourcefile, __LINE__))
    end if

    if (map2gcell) then  ! Map to gridcell

       avgflag_trim = avgflag(1:1)
       ! note that in this case beg1d = begg and end1d=endg
       select case (avgflag_trim)
       case ('I') ! Instantaneous
          do j = 1,num2d
             do k = beg1d_out, end1d_out
                if (field_gcell(k,j) /= spval) then
                   hbuf(k,j) = field_gcell(k,j)
                else
                   hbuf(k,j) = spval
                end if
                nacs(k,j) = 1
             end do
          end do
       case ('A', 'S') ! Time average / sum
          do j = 1,num2d
             do k = beg1d_out, end1d_out
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
             do k = beg1d_out, end1d_out
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
             do k = beg1d_out, end1d_out
                if (field_gcell(k,j) /= spval) then
                   if (nacs(k,j) == 0) hbuf(k,j) = +1.e50_r8
                   hbuf(k,j) = min( hbuf(k,j), field_gcell(k,j) )
                else
                   hbuf(k,j) = spval
                endif
                nacs(k,j) = 1
             end do
          end do
       case ('L') ! Local solar time
          read(avgflag(2:6), *) tod
          do j = 1,num2d
             do k = beg1d_out, end1d_out
                if (field_gcell(k,j) /= spval) then
                   local_secpl = secs + grc%londeg(k)/degpsec
                   local_secpl = mod(local_secpl,isecspday)
                   weight = calc_weight_local_time(local_secpl, tod)
                   if (weight > 0) then
                      if (nacs(k,j) == 0) hbuf(k,j) = 0._r8
                      hbuf(k,j) = hbuf(k,j) + field_gcell(k,j)*real(weight)
                      nacs(k,j) = nacs(k,j) + weight
                   end if
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
             end do
          end do
       case default
          write(iulog,*) trim(subname),' ERROR: invalid time averaging flag ', avgflag
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select
       deallocate( field_gcell )

    else  ! Do not map to gridcell

       ! For data defined on the pft, col or landunit, we need to check if a point is active
       ! to determine whether that point should be assigned spval
       if (type1d == namep) then
          check_active = .true.
          active => patch%active
          allocate(grid_index(bounds%begg:bounds%endg) )
          grid_index = patch%gridcell
       else if (type1d == namec) then
          check_active = .true.
          active => col%active
          allocate(grid_index(bounds%begg:bounds%endg) )
          grid_index = col%gridcell
       else if (type1d == namel) then
          check_active = .true.
          active =>lun%active
          allocate(grid_index(bounds%begg:bounds%endg) )
          grid_index = lun%gridcell
       else
          check_active = .false.
       end if

       ! Note that since field points to an array section the
       ! bounds are field(1:end1d-beg1d+1, num2d) - therefore
       ! need to do the shifting below
       avgflag_trim = avgflag(1:1)
       select case (avgflag_trim)
       case ('I') ! Instantaneous
          do j = 1,num2d
             do k = beg1d,end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) valid = .false.
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
       case ('A', 'S') ! Time average / sum
          do j = 1,num2d
             do k = beg1d,end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) valid = .false.
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
             do k = beg1d,end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) valid = .false.
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
             do k = beg1d,end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) valid = .false.
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
       case ('L') ! Local solar time
          read(avgflag(2:6), *) tod
          do j = 1,num2d
             do k = beg1d, end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) then
                      valid = .false.
                   else
                      local_secpl = secs + grc%londeg(grid_index(k))/degpsec
                   end if
                else
                   local_secpl = secs + grc%londeg(k)/degpsec
                end if
                local_secpl = mod(local_secpl,isecspday)
                if (valid) then
                   weight = calc_weight_local_time(local_secpl, tod)
                   if (weight > 0 .and. field(k-beg1d+1,j) /= spval) then
                      if (nacs(k,j) == 0) hbuf(k,j) = 0._r8
                      hbuf(k,j) = hbuf(k,j) + field(k-beg1d+1,j)*real(weight)
                      nacs(k,j) = nacs(k,j) + weight
                   end if
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
             end do
          end do
       case default
          write(iulog,*) trim(subname),' ERROR: invalid time averaging flag ', avgflag
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select
    end if

    if (field_allocated) then
       deallocate(field)
    end if

  end subroutine hist_update_hbuf_field_2d

  !-----------------------------------------------------------------------
  function calc_weight_local_time(local_secpl, tod) result(weight)
    !
    ! !DESCRIPTION:
    ! Calculates weight for linear intepolation in time for local time
    ! average flag
    !
    ! !USES:
    use clm_varcon      , only : isecspday
    use clm_time_manager, only : get_step_size
    !
    ! !ARGUMENTS:
    integer                :: weight       ! function result
    integer, intent(inout) :: local_secpl  ! seconds into current date in local time
    integer, intent(in)    :: tod          ! Desired local solar time of output in seconds

    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'calc_weight_local_time'
    integer :: dtime                       ! timestep size [seconds]
    !-----------------------------------------------------------------------

    weight = 0
    dtime = get_step_size()

    if (tod < dtime .and. local_secpl > isecspday-dtime ) then
       local_secpl = local_secpl - isecspday
    end if
    if (local_secpl >= tod - dtime .and. local_secpl < tod ) then
       weight = dtime-tod+local_secpl
    else if (local_secpl >= tod .and. local_secpl < tod + dtime ) then
       weight = dtime+tod-local_secpl
    end if

  end function calc_weight_local_time

  !-----------------------------------------------------------------------
  subroutine hist_set_snow_field_2d (field_out, field_in, no_snow_behavior, type1d, beg1d, end1d)
    !
    ! !DESCRIPTION:
    ! Set values in history field dimensioned by levsno.
    !
    ! This routine handles what to do when a given snow layer doesn't exist for a given
    ! point, based on the no_snow_behavior argument. Options are:
    !
    ! - no_snow_normal: This is the normal behavior, which applies to most snow fields:
    !   Use spval (missing value flag). This means that temporal averages will just
    !   consider times when a particular snow layer actually existed
    !
    ! - no_snow_zero: Average in a 0 value for times when the snow layer isn't present
    !
    ! Input and output fields can be defined at the patch or column level
    !
    ! !ARGUMENTS:
    integer         , intent(in)  :: beg1d                    ! beginning spatial index
    integer         , intent(in)  :: end1d                    ! ending spatial index
    real(r8)        , intent(out) :: field_out( beg1d: , 1: ) ! output field [point, lev]
    real(r8)        , intent(in)  :: field_in ( beg1d: , 1: ) ! input field [point, lev]
    integer         , intent(in)  :: no_snow_behavior         ! behavior to use when a snow layer is absent
    character(len=*), intent(in)  :: type1d                   ! 1d clm pointer type ("column" or "pft")
    !
    ! !LOCAL VARIABLES:
    integer :: num_levels             ! total number of possible snow layers
    integer :: point
    integer :: level
    integer :: num_snow_layers        ! number of snow layers that exist at a point
    integer :: num_nonexistent_layers
    integer :: c                      ! column index
    real(r8):: no_snow_val            ! value to use when a snow layer is missing
    character(len=*), parameter :: subname = 'hist_set_snow_field_2d'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(field_out, 1) == end1d), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(field_in , 1) == end1d), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(field_out, 2) == ubound(field_in, 2)), sourcefile, __LINE__)

    associate(&
    snl            => col%snl  &   ! Input: [integer (:)] number of snow layers (negative)
    )

    num_levels = ubound(field_in, 2)

    ! Determine no_snow_val
    select case (no_snow_behavior)
    case (no_snow_normal)
       no_snow_val = spval
    case (no_snow_zero)
       no_snow_val = 0._r8
    case default
       write(iulog,*) trim(subname), ' ERROR: unrecognized no_snow_behavior: ', &
            no_snow_behavior
       call endrun()
    end select

    do point = beg1d, end1d

       ! Get number of snow layers at this point

       if (type1d == namec) then
          c = point
       else if (type1d == namep) then
          c = patch%column(point)
       else
          write(iulog,*) trim(subname), ' ERROR: Only implemented for patch and col-level fields'
          write(iulog,*) 'type1d = ', trim(type1d)
          call endrun()
       end if

       num_snow_layers = abs(snl(c))
       num_nonexistent_layers = num_levels - num_snow_layers

       ! Fill output field appropriately for each layer
       ! When only a subset of snow layers exist, it is the LAST num_snow_layers that exist
       ! Levels are rearranged such that the top snow layer (surface layer) becomes level 1, etc.

       do level = num_levels, (num_levels-num_nonexistent_layers+1), -1
          field_out(point, level) = no_snow_val
       end do
       do level = (num_levels-num_nonexistent_layers), 1, -1
          field_out(point, level) = field_in(point, level+num_nonexistent_layers)
       end do

    end do

    end associate

  end subroutine hist_set_snow_field_2d


  !-----------------------------------------------------------------------
  subroutine hfields_normalize (t)
    !
    ! !DESCRIPTION:
    ! Normalize fields on a history file by the number of accumulations.
    ! Loop over fields on the tape.  Need averaging flag and number of
    ! accumulations to perform normalization.
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t       ! tape index
    !
    ! !LOCAL VARIABLES:
    integer :: f                   ! field index
    integer :: k                   ! 1d index
    integer :: j                   ! 2d index
    logical :: aflag               ! averaging flag
    integer :: beg1d,end1d         ! hbuf 1d beginning and ending indices
    integer :: num2d               ! hbuf size of second dimension (e.g. number of vertical levels)
    character(len=avgflag_strlen)  :: avgflag   ! averaging flag
    real(r8), pointer :: hbuf(:,:) ! history buffer
    integer , pointer :: nacs(:,:) ! accumulation counter
    character(len=*),parameter :: subname = 'hfields_normalize'
    !-----------------------------------------------------------------------

    ! Normalize by number of accumulations for time averaged case

    do f = 1,tape(t)%nflds
       avgflag   =  tape(t)%hlist(f)%avgflag
       if ( is_mapping_upto_subgrid(tape(t)%hlist(f)%field%type1d, tape(t)%hlist(f)%field%type1d_out) )then
          beg1d =  tape(t)%hlist(f)%field%beg1d_out
          end1d =  tape(t)%hlist(f)%field%end1d_out
       else
          beg1d =  tape(t)%hlist(f)%field%beg1d
          end1d =  tape(t)%hlist(f)%field%end1d
       end if
       num2d     =  tape(t)%hlist(f)%field%num2d
       nacs      => tape(t)%hlist(f)%nacs
       hbuf      => tape(t)%hlist(f)%hbuf

       if (avgflag == 'A' .or. avgflag(1:1) == 'L') then
          aflag = .true.
       else
          aflag = .false.
       end if

       do j = 1, num2d
          do k = beg1d, end1d
             if (aflag .and. nacs(k,j) /= 0) then
                hbuf(k,j) = hbuf(k,j) / float(nacs(k,j))
             elseif (avgflag(1:1) == 'L' .and. nacs(k,j) == 0) then
                hbuf(k,j) = spval
             end if
          end do
       end do
    end do

  end subroutine hfields_normalize

  !-----------------------------------------------------------------------
  subroutine hfields_zero (t)
    !
    ! !DESCRIPTION:
    ! Zero out accumulation and history buffers for a given history tape.
    ! Loop through fields on the tape.
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t     ! tape index
    !
    ! !LOCAL VARIABLES:
    integer :: f                 ! field index
    character(len=*),parameter :: subname = 'hfields_zero'
    !-----------------------------------------------------------------------

    do f = 1,tape(t)%nflds
       tape(t)%hlist(f)%hbuf(:,:) = 0._r8
       tape(t)%hlist(f)%nacs(:,:) = 0
    end do

  end subroutine hfields_zero

  !-----------------------------------------------------------------------
  subroutine htape_create (t, histrest)
    !
    ! !DESCRIPTION:
    ! Define netcdf metadata of history file t.
    !
    ! !USES:
    use clm_varpar      , only : nlevgrnd, nlevsno, nlevlak, nlevurb, nlevmaxurbgrnd, numrad, nlevcan, nvegwcs,nlevsoi
    use clm_varpar      , only : natpft_size, cft_size, maxpatch_glc, nlevdecomp_full, mxsowings, mxharvests
    use landunit_varcon , only : max_lunit
    use clm_varctl      , only : caseid, ctitle, fsurdat, finidat, paramfile
    use clm_varctl      , only : hillslope_file
    use clm_varctl      , only : version, hostname, username, conventions, source
    use clm_varctl      , only : use_hillslope,nhillslope,max_columns_hillslope
    use domainMod       , only : ldomain
    use fileutils       , only : get_filename
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t                   ! tape index
    logical, intent(in), optional :: histrest  ! if creating the history restart file
    !
    ! !LOCAL VARIABLES:
    integer :: f                   ! field index
    integer :: p,c,l,n             ! indices
    integer :: ier                 ! error code
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
    integer :: numa                ! total number of atm cells across all processors
    logical :: lhistrest           ! local history restart flag
    type(file_desc_t), pointer :: lnfid     ! local file id
    character(len=  8) :: curdate  ! current date
    character(len=  8) :: curtime  ! current time
    character(len=256) :: name     ! name of attribute
    character(len=256) :: units    ! units of attribute
    character(len=256) :: str      ! global attribute string
    character(len=*),parameter :: subname = 'htape_create'
    !-----------------------------------------------------------------------

    if ( present(histrest) )then
       lhistrest = histrest
    else
       lhistrest = .false.
    end if

    ! Determine necessary indices

    call get_proc_global(ng=numg, nl=numl, nc=numc, np=nump)

    ! define output write precsion for tape

    ncprec = tape(t)%ncprec
    if (lhistrest) then
       lnfid => ncid_hist(t)
    else
       lnfid => nfid(t)
    endif

    ! Create new netCDF file. It will be in define mode

    if ( .not. lhistrest )then
       if (masterproc) then
          write(iulog,*) trim(subname),' : Opening netcdf htape ', &
                                      trim(locfnh(t))
          call shr_sys_flush(iulog)
       end if
       call ncd_pio_createfile(lnfid, trim(locfnh(t)))
       call ncd_putatt(lnfid, ncd_global, 'title', 'CLM History file information' )
       call ncd_putatt(lnfid, ncd_global, 'comment', &
          "NOTE: None of the variables are weighted by land fraction!" )
    else
       if (masterproc) then
          write(iulog,*) trim(subname),' : Opening netcdf rhtape ', &
                                      trim(locfnhr(t))
          call shr_sys_flush(iulog)
       end if
       call ncd_pio_createfile(lnfid, trim(locfnhr(t)))
       call ncd_putatt(lnfid, ncd_global, 'title', &
          'CLM Restart History information, required to continue a simulation' )
       call ncd_putatt(lnfid, ncd_global, 'comment', &
                       "This entire file NOT needed for startup or branch simulations")
    end if

    ! Create global attributes. Attributes are used to store information
    ! about the data set. Global attributes are information about the
    ! data set as a whole, as opposed to a single variable

    call ncd_putatt(lnfid, ncd_global, 'Conventions', trim(conventions))
    call getdatetime(curdate, curtime)
    str = 'created on ' // curdate // ' ' // curtime
    call ncd_putatt(lnfid, ncd_global, 'history' , trim(str))
    call ncd_putatt(lnfid, ncd_global, 'source'  , trim(source))
    call ncd_putatt(lnfid, ncd_global, 'hostname', trim(hostname))
    call ncd_putatt(lnfid, ncd_global, 'username', trim(username))
    call ncd_putatt(lnfid, ncd_global, 'version' , trim(version))

    str = &
    '$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $'
    call ncd_putatt(lnfid, ncd_global, 'revision_id', trim(str))
    call ncd_putatt(lnfid, ncd_global, 'case_title', trim(ctitle))
    call ncd_putatt(lnfid, ncd_global, 'case_id', trim(caseid))
    str = get_filename(fsurdat)
    call ncd_putatt(lnfid, ncd_global, 'Surface_dataset', trim(str))
    str = get_filename(hillslope_file)
    call ncd_putatt(lnfid, ncd_global, 'Hillslope_dataset', trim(str))
    if (finidat == ' ') then
       str = 'arbitrary initialization'
    else
       str = get_filename(finidat)
    endif
    call ncd_putatt(lnfid, ncd_global, 'Initial_conditions_dataset', trim(str))
    str = get_filename(paramfile)
    call ncd_putatt(lnfid, ncd_global, 'PFT_physiological_constants_dataset', trim(str))

    ! Define dimensions.
    ! Time is an unlimited dimension. Character string is treated as an array of characters.

    ! Global uncompressed dimensions (including non-land points)
    if (ldomain%isgrid2d) then
       call ncd_defdim(lnfid, 'lon'   , ldomain%ni, dimid)
       call ncd_defdim(lnfid, 'lat'   , ldomain%nj, dimid)
    else
       call ncd_defdim(lnfid, trim(grlnd), ldomain%ns, dimid)
    end if

    ! Global compressed dimensions (not including non-land points)
    call ncd_defdim(lnfid, trim(nameg), numg, dimid)
    call ncd_defdim(lnfid, trim(namel), numl, dimid)
    call ncd_defdim(lnfid, trim(namec), numc, dimid)
    call ncd_defdim(lnfid, trim(namep), nump, dimid)

    ! "level" dimensions
    call ncd_defdim(lnfid, 'levgrnd', nlevgrnd, dimid)
    call ncd_defdim(lnfid, 'levsoi', nlevsoi, dimid)
    if (nlevurb > 0) then
       call ncd_defdim(lnfid, 'levurb' , nlevurb, dimid)
    end if
    call ncd_defdim(lnfid, 'levmaxurbgrnd' , nlevmaxurbgrnd, dimid)
    call ncd_defdim(lnfid, 'levlak' , nlevlak, dimid)
    call ncd_defdim(lnfid, 'numrad' , numrad , dimid)
    call ncd_defdim(lnfid, 'levsno' , nlevsno , dimid)
    call ncd_defdim(lnfid, 'ltype', max_lunit, dimid)
    call ncd_defdim(lnfid, 'nlevcan',nlevcan, dimid)
    call ncd_defdim(lnfid, 'nvegwcs',nvegwcs, dimid)
    if (use_hillslope) then
       call ncd_defdim(lnfid, 'nhillslope',nhillslope, dimid)
       call ncd_defdim(lnfid, 'max_columns_hillslope',max_columns_hillslope, dimid)
    endif
    call ncd_defdim(lnfid, 'mxsowings' , mxsowings , dimid)
    call ncd_defdim(lnfid, 'mxharvests' , mxharvests , dimid)
    call htape_add_ltype_metadata(lnfid)
    call htape_add_ctype_metadata(lnfid)
    call ncd_defdim(lnfid, 'natpft', natpft_size, dimid)
    if (cft_size > 0) then
       call ncd_defdim(lnfid, 'cft', cft_size, dimid)
       call htape_add_cft_metadata(lnfid)
    end if
    call ncd_defdim(lnfid, 'glc_nec' , maxpatch_glc , dimid)
    ! elevclas (in contrast to glc_nec) includes elevation class 0 (bare land)
    ! (although on the history file it will go 1:(nec+1) rather than 0:nec)
    call ncd_defdim(lnfid, 'elevclas' , maxpatch_glc + 1, dimid)

    call ncd_defdim(lnfid, 'string_length', hist_dim_name_length, strlen_dimid)
    call ncd_defdim(lnfid, 'scale_type_string_length', scale_type_strlen, dimid)
    call ncd_defdim( lnfid, 'levdcmp', nlevdecomp_full, dimid)

    if(use_fates)then
       call ncd_defdim(lnfid, 'fates_levscag', nlevsclass * nlevage, dimid)
       call ncd_defdim(lnfid, 'fates_levscagpf', nlevsclass * nlevage * numpft_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levagepft', nlevage * numpft_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levscls', nlevsclass, dimid)
       call ncd_defdim(lnfid, 'fates_levcacls', nlevcoage, dimid)
       call ncd_defdim(lnfid, 'fates_levpft', numpft_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levage', nlevage, dimid)
       call ncd_defdim(lnfid, 'fates_levheight', nlevheight, dimid)
       call ncd_defdim(lnfid, 'fates_levfuel', num_fuel_classes, dimid)
       call ncd_defdim(lnfid, 'fates_levcwdsc', ncwd, dimid)
       call ncd_defdim(lnfid, 'fates_levscpf', nlevsclass*numpft_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levcapf', nlevcoage*numpft_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levcan', nclmax, dimid)
       call ncd_defdim(lnfid, 'fates_levleaf', nlevleaf, dimid)
       call ncd_defdim(lnfid, 'fates_levcnlf', nlevleaf * nclmax, dimid)
       call ncd_defdim(lnfid, 'fates_levcnlfpf', nlevleaf * nclmax * numpft_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levcdsc', nlevdamage * nlevsclass, dimid)
       call ncd_defdim(lnfid, 'fates_levcdpf', nlevdamage * nlevsclass * numpft_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levcdam', nlevdamage, dimid)
       call ncd_defdim(lnfid, 'fates_levelem', num_elements_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levelpft', num_elements_fates * numpft_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levelcwd', num_elements_fates * ncwd, dimid)
       call ncd_defdim(lnfid, 'fates_levelage', num_elements_fates * nlevage, dimid)
       call ncd_defdim(lnfid, 'fates_levagefuel', nlevage * num_fuel_classes, dimid)
       call ncd_defdim(lnfid, 'fates_levclscpf', nclmax*nlevsclass*numpft_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levlanduse', n_landuse_cats, dimid)
       call ncd_defdim(lnfid, 'fates_levlulu', n_landuse_cats * n_landuse_cats, dimid)
    end if

    if ( .not. lhistrest )then
       call ncd_defdim(lnfid, 'nbnd', 2, nbnd_dimid)
       call ncd_defdim(lnfid, 'time', ncd_unlimited, time_dimid)
       if (masterproc)then
          write(iulog,*) trim(subname), &
                          ' : Successfully defined netcdf history file ',t
          call shr_sys_flush(iulog)
       end if
    else
       if (masterproc)then
          write(iulog,*) trim(subname), &
                          ' : Successfully defined netcdf restart history file ',t
          call shr_sys_flush(iulog)
       end if
    end if

  end subroutine htape_create

  !-----------------------------------------------------------------------
  subroutine htape_add_ltype_metadata(lnfid)
    !
    ! !DESCRIPTION:
    ! Add global metadata defining landunit types
    !
    ! !USES:
    use landunit_varcon, only : max_lunit, landunit_names, landunit_name_length
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: lnfid ! local file id
    !
    ! !LOCAL VARIABLES:
    integer :: ltype  ! landunit type
    character(len=*), parameter :: att_prefix = 'ltype_'  ! prefix for attributes
    character(len=len(att_prefix)+landunit_name_length) :: attname ! attribute name

    character(len=*), parameter :: subname = 'htape_add_ltype_metadata'
    !-----------------------------------------------------------------------

    do ltype = 1, max_lunit
       attname = att_prefix // landunit_names(ltype)
       call ncd_putatt(lnfid, ncd_global, attname, ltype)
    end do

  end subroutine htape_add_ltype_metadata

  !-----------------------------------------------------------------------
  subroutine htape_add_ctype_metadata(lnfid)
    !
    ! !DESCRIPTION:
    ! Add global metadata defining column types
    !
    ! !USES:
    use column_varcon, only : write_coltype_metadata
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: lnfid ! local file id
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: att_prefix = 'ctype_'  ! prefix for attributes

    character(len=*), parameter :: subname = 'htape_add_ctype_metadata'
    !-----------------------------------------------------------------------

    call write_coltype_metadata(att_prefix, lnfid)

  end subroutine htape_add_ctype_metadata

  !-----------------------------------------------------------------------
  subroutine htape_add_natpft_metadata(lnfid)
    !
    ! !DESCRIPTION:
    ! Add global metadata defining natpft types
    !
    ! !USES:
    use clm_varpar, only : natpft_lb, natpft_ub
    use pftconMod , only : pftname_len, pftname
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: lnfid ! local file id
    !
    ! !LOCAL VARIABLES:
    integer :: ptype  ! patch type
    integer :: ptype_1_indexing ! patch type, translated to 1 indexing
    character(len=*), parameter :: att_prefix = 'natpft_' ! prefix for attributes
    character(len=len(att_prefix)+pftname_len) :: attname ! attribute name

    character(len=*), parameter :: subname = 'htape_add_natpft_metadata'
    !-----------------------------------------------------------------------

    do ptype = natpft_lb, natpft_ub
       ptype_1_indexing = ptype + (1 - natpft_lb)
       attname = att_prefix // pftname(ptype)
       call ncd_putatt(lnfid, ncd_global, attname, ptype_1_indexing)
    end do

  end subroutine htape_add_natpft_metadata

  !-----------------------------------------------------------------------
  subroutine htape_add_cft_metadata(lnfid)
    !
    ! !DESCRIPTION:
    ! Add global metadata defining natpft types
    !
    ! !USES:
    use clm_varpar, only : cft_lb, cft_ub
    use pftconMod , only : pftname_len, pftname
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: lnfid ! local file id
    !
    ! !LOCAL VARIABLES:
    integer :: ptype  ! patch type
    integer :: ptype_1_indexing ! patch type, translated to 1 indexing
    character(len=*), parameter :: att_prefix = 'cft_'    ! prefix for attributes
    character(len=len(att_prefix)+pftname_len) :: attname ! attribute name

    character(len=*), parameter :: subname = 'htape_add_cft_metadata'
    !-----------------------------------------------------------------------

    do ptype = cft_lb, cft_ub
       ptype_1_indexing = ptype + (1 - cft_lb)
       attname = att_prefix // pftname(ptype)
       call ncd_putatt(lnfid, ncd_global, attname, ptype_1_indexing)
    end do

  end subroutine htape_add_cft_metadata

  !-----------------------------------------------------------------------
  subroutine htape_timeconst3D(t, &
       bounds, watsat_col, sucsat_col, bsw_col, hksat_col, &
       cellsand_col, cellclay_col, mode)
    !
    ! !DESCRIPTION:
    ! Write time constant 3D variables to history tapes.
    ! Only write out when this subroutine is called (normally only for
    ! primary history files at very first time-step, nstep=1).
    ! Issue the required netcdf wrapper calls to define the history file
    ! contents.
    !
    ! !USES:
    use subgridAveMod  , only : c2g
    use clm_varpar     , only : nlevgrnd ,nlevlak, nlevmaxurbgrnd, nlevsoi
    use shr_string_mod , only : shr_string_listAppend
    use domainMod      , only : ldomain
    !
    ! !ARGUMENTS:
    integer           , intent(in) :: t    ! tape index
    type(bounds_type) , intent(in) :: bounds
    real(r8)          , intent(in) :: watsat_col( bounds%begc:,1: )
    real(r8)          , intent(in) :: sucsat_col( bounds%begc:,1: )
    real(r8)          , intent(in) :: bsw_col( bounds%begc:,1: )
    real(r8)          , intent(in) :: hksat_col( bounds%begc:,1: )
    real(r8)          , intent(in) :: cellsand_col( bounds%begc:,1: )
    real(r8)          , intent(in) :: cellclay_col( bounds%begc:,1: )
    character(len=*)  , intent(in) :: mode ! 'define' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: c,l,lev,ifld               ! indices
    integer :: ier                        ! error status
    character(len=max_chars) :: long_name ! variable long name
    character(len=max_namlen):: varname   ! variable name
    character(len=max_namlen):: units     ! variable units
    integer :: varid                      ! variable id
    !
    real(r8), pointer :: histi(:,:)       ! temporary
    real(r8), pointer :: histo(:,:)       ! temporary
    integer, parameter :: nflds = 6       ! Number of 3D time-constant fields
    character(len=*),parameter :: subname = 'htape_timeconst3D'
    character(len=*),parameter :: varnames(nflds) = (/ &
                                                        'ZSOI  ', &
                                                        'DZSOI ', &
                                                        'WATSAT', &
                                                        'SUCSAT', &
                                                        'BSW   ', &
                                                        'HKSAT '  &
                                                    /)
    ! Scale type for subgrid averaging of landunits to grid cells
    ! WJS (10-25-11): Note about l2g_scale_type in the following: ZSOI & DZSOI are
    ! currently constant in space, except for urban points, so their scale type
    ! doesn't matter at the moment as long as it excludes urban points. I am using
    ! 'nonurb' so that the values are output everywhere where the fields are
    ! constant (i.e., everywhere except urban points). For the other fields, I am
    ! using 'veg' to be consistent with the l2g_scale_type that is now used for many
    ! of the 3-d time-variant fields; in theory, though, one might want versions of
    ! these variables output for different landunits.
    character(len=scale_type_strlen) :: l2g_scale_type(nflds) = [ &
         'nonurb', &  ! ZSOI
         'nonurb', &  ! DZSOI
         'veg   ', &  ! WATSAT
         'veg   ', &  ! SUCSAT
         'veg   ', &  ! BSW
         'veg   '  &  ! HKSAT
         ]
    real(r8), pointer :: histil(:,:)      ! temporary
    real(r8), pointer :: histol(:,:)
    integer, parameter :: nfldsl = 2
    character(len=*),parameter :: varnamesl(nfldsl) = (/ &
                                                          'ZLAKE ', &
                                                          'DZLAKE' &
                                                      /)
    real(r8), pointer :: histit(:,:)      ! temporary
    real(r8), pointer :: histot(:,:)
    integer, parameter :: nfldst = 2
    character(len=*),parameter :: varnamest(nfldst) = (/ &
                                                          'PCT_SAND ', &
                                                          'PCT_CLAY '  &
                                                      /)
    ! Scale type for subgrid averaging of landunits to grid cells, for lake fields
    character(len=scale_type_strlen) :: l2g_scale_typel(nfldsl) = [ &
         'lake', &  ! ZLAKE
         'lake'  &  ! DZLAKE
         ]

    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(watsat_col)   == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sucsat_col)   == (/bounds%endc, nlevgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(bsw_col)      == (/bounds%endc, nlevgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(hksat_col)    == (/bounds%endc, nlevgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(cellsand_col) == (/bounds%endc, nlevsoi/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(cellclay_col) == (/bounds%endc, nlevsoi/)), sourcefile, __LINE__)

    !-------------------------------------------------------------------------------
    !***      Non-time varying 3D fields                    ***
    !***      Only write out when this subroutine is called ***
    !***       Normally only called once for primary tapes  ***
    !-------------------------------------------------------------------------------

    if (mode == 'define') then

       do ifld = 1,nflds
          ! Field indices MUST match varnames array order above!
          if (ifld == 1) then
             long_name='soil depth'; units = 'm'
          else if (ifld == 2) then
             long_name='soil thickness'; units = 'm'
          else if (ifld == 3) then
             long_name='saturated soil water content (porosity)';  units = 'mm3/mm3'
          else if (ifld == 4) then
             long_name='saturated soil matric potential'; units = 'mm'
          else if (ifld == 5) then
             long_name='slope of soil water retention curve'; units = 'unitless'
          else if (ifld == 6) then
             long_name='saturated hydraulic conductivity'; units = 'mm s-1'
          else
             call endrun(msg=' ERROR: bad 3D time-constant field index'//errMsg(sourcefile, __LINE__))
          end if
          if (tape(t)%dov2xy) then
             if (ldomain%isgrid2d) then
                call ncd_defvar(ncid=nfid(t), varname=trim(varnames(ifld)), xtype=tape(t)%ncprec,&
                     dim1name='lon', dim2name='lat', dim3name='levgrnd', &
                     long_name=long_name, units=units, missing_value=spval, fill_value=spval, &
                     varid=varid)
             else
                call ncd_defvar(ncid=nfid(t), varname=trim(varnames(ifld)), xtype=tape(t)%ncprec, &
                     dim1name=grlnd, dim2name='levgrnd', &
                     long_name=long_name, units=units, missing_value=spval, fill_value=spval, &
                     varid=varid)
             end if

             call add_landunit_mask_metadata(nfid(t), varid, l2g_scale_type(ifld))
          else
             call ncd_defvar(ncid=nfid(t), varname=trim(varnames(ifld)), xtype=tape(t)%ncprec, &
                  dim1name=namec, dim2name='levgrnd', &
                  long_name=long_name, units=units, missing_value=spval, fill_value=spval)
          end if
          call shr_string_listAppend(TimeConst3DVars,varnames(ifld))
       end do

    else if (mode == 'write') then

       allocate(histi(bounds%begc:bounds%endc,nlevgrnd), stat=ier)
       if (ier /= 0) then
          write(iulog,*) trim(subname),' ERROR: allocation error for histi'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

       ! Write time constant fields

       if (tape(t)%dov2xy) then
          allocate(histo(bounds%begg:bounds%endg,nlevgrnd), stat=ier)
          if (ier /= 0) then
             write(iulog,*)  trim(subname),' ERROR: allocation error for histo'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       end if

       do ifld = 1,nflds

          histi(:,:) = spval
          do lev = 1,nlevgrnd
             do c = bounds%begc,bounds%endc
                l = col%landunit(c)
                   ! Field indices MUST match varnames array order above!
                   if (ifld ==1) histi(c,lev) = col%z(c,lev)
                   if (ifld ==2) histi(c,lev) = col%dz(c,lev)
                   if (ifld ==3) histi(c,lev) = watsat_col(c,lev)
                   if (ifld ==4) histi(c,lev) = sucsat_col(c,lev)
                   if (ifld ==5) histi(c,lev) = bsw_col(c,lev)
                   if (ifld ==6) histi(c,lev) = hksat_col(c,lev)
             end do
          end do
          if (tape(t)%dov2xy) then
             histo(:,:) = spval

             call c2g(bounds, nlevgrnd, &
                  histi(bounds%begc:bounds%endc, :), &
                  histo(bounds%begg:bounds%endg, :), &
                  c2l_scale_type='unity', l2g_scale_type=l2g_scale_type(ifld))

             if (ldomain%isgrid2d) then
                call ncd_io(varname=trim(varnames(ifld)), dim1name=grlnd, &
                     data=histo, ncid=nfid(t), flag='write')
             else
                call ncd_io(varname=trim(varnames(ifld)), dim1name=grlnd, &
                     data=histo, ncid=nfid(t), flag='write')
             end if
          else
             call ncd_io(varname=trim(varnames(ifld)), dim1name=namec, &
                  data=histi, ncid=nfid(t), flag='write')
          end if
       end do

       if (tape(t)%dov2xy) deallocate(histo)
       deallocate(histi)

    end if  ! (define/write mode

    if (mode == 'define') then
       do ifld = 1,nfldsl
          ! Field indices MUST match varnamesl array order above!
          if (ifld == 1) then
             long_name='lake layer node depth'; units = 'm'
          else if (ifld == 2) then
             long_name='lake layer thickness'; units = 'm'
          else
             call endrun(msg=' ERROR: bad 3D time-constant field index'//errMsg(sourcefile, __LINE__))
          end if
          if (tape(t)%dov2xy) then
             if (ldomain%isgrid2d) then
                call ncd_defvar(ncid=nfid(t), varname=trim(varnamesl(ifld)), xtype=tape(t)%ncprec,&
                     dim1name='lon', dim2name='lat', dim3name='levlak', &
                     long_name=long_name, units=units, missing_value=spval, fill_value=spval, &
                     varid=varid)
             else
                call ncd_defvar(ncid=nfid(t), varname=trim(varnamesl(ifld)), xtype=tape(t)%ncprec, &
                     dim1name=grlnd, dim2name='levlak', &
                     long_name=long_name, units=units, missing_value=spval, fill_value=spval, &
                     varid=varid)
             end if

             call add_landunit_mask_metadata(nfid(t), varid, l2g_scale_typel(ifld))
          else
             call ncd_defvar(ncid=nfid(t), varname=trim(varnamesl(ifld)), xtype=tape(t)%ncprec, &
                  dim1name=namec, dim2name='levlak', &
                  long_name=long_name, units=units, missing_value=spval, fill_value=spval)
          end if
          call shr_string_listAppend(TimeConst3DVars,varnamesl(ifld))
       end do

    else if (mode == 'write') then

       allocate(histil(bounds%begc:bounds%endc,nlevlak), stat=ier)
       if (ier /= 0) then
          write(iulog,*) trim(subname),' ERROR: allocation error for histil'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

       ! Write time constant fields

       if (tape(t)%dov2xy) then
          allocate(histol(bounds%begg:bounds%endg,nlevlak), stat=ier)
          if (ier /= 0) then
             write(iulog,*)  trim(subname),' ERROR: allocation error for histol'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       end if

       do ifld = 1,nfldsl
          histil(:,:) = spval
          do lev = 1,nlevlak
             do c = bounds%begc,bounds%endc
                l = col%landunit(c)
                if (lun%lakpoi(l)) then
                   ! Field indices MUST match varnamesl array order above!
                   if (ifld ==1) histil(c,lev) = col%z_lake(c,lev)
                   if (ifld ==2) histil(c,lev) = col%dz_lake(c,lev)
                end if
             end do
          end do
          if (tape(t)%dov2xy) then
             histol(:,:) = spval
             call c2g(bounds, nlevlak, &
                  histil(bounds%begc:bounds%endc, :), &
                  histol(bounds%begg:bounds%endg, :), &
                  c2l_scale_type='unity', l2g_scale_type=l2g_scale_typel(ifld))
             if (ldomain%isgrid2d) then
                call ncd_io(varname=trim(varnamesl(ifld)), dim1name=grlnd, &
                     data=histol, ncid=nfid(t), flag='write')
             else
                call ncd_io(varname=trim(varnamesl(ifld)), dim1name=grlnd, &
                     data=histol, ncid=nfid(t), flag='write')
             end if
          else
             call ncd_io(varname=trim(varnamesl(ifld)), dim1name=namec,  &
                  data=histil, ncid=nfid(t), flag='write')
          end if
       end do

       if (tape(t)%dov2xy) deallocate(histol)
       deallocate(histil)

    end if  ! (define/write mode

    if (mode == 'define') then
       do ifld = 1,nfldst
          ! Field indices MUST match varnamest array order above!
          if (ifld == 1) then
             long_name='percent sand'; units = 'percent'
          else if (ifld == 2) then
             long_name='percent clay'; units = 'percent'
          else
             call endrun(msg=' ERROR: bad 3D time-constant field index'//errMsg(sourcefile, __LINE__))
          end if
          if (tape(t)%dov2xy) then
             if (ldomain%isgrid2d) then
                call ncd_defvar(ncid=nfid(t), varname=trim(varnamest(ifld)), xtype=tape(t)%ncprec,&
                     dim1name='lon', dim2name='lat', dim3name='levsoi', &
                     long_name=long_name, units=units, missing_value=spval, fill_value=spval)
             else
                call ncd_defvar(ncid=nfid(t), varname=trim(varnamest(ifld)), xtype=tape(t)%ncprec, &
                     dim1name=grlnd, dim2name='levsoi', &
                     long_name=long_name, units=units, missing_value=spval, fill_value=spval)
             end if
          else
             call ncd_defvar(ncid=nfid(t), varname=trim(varnamest(ifld)), xtype=tape(t)%ncprec, &
                  dim1name=namec, dim2name='levsoi', &
                  long_name=long_name, units=units, missing_value=spval, fill_value=spval)
          end if
          call shr_string_listAppend(TimeConst3DVars,varnamest(ifld))
       end do

    else if (mode == 'write') then

       allocate(histit(bounds%begc:bounds%endc,nlevsoi), stat=ier)
       if (ier /= 0) then
          write(iulog,*) trim(subname),' ERROR: allocation error for histit'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

       ! Write time constant fields

       if (tape(t)%dov2xy) then
          allocate(histot(bounds%begg:bounds%endg,nlevsoi), stat=ier)
          if (ier /= 0) then
             write(iulog,*)  trim(subname),' ERROR: allocation error for histot'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       end if

       do ifld = 1,nfldst
          histit(:,:) = spval
          do lev = 1,nlevsoi
             do c = bounds%begc,bounds%endc
                ! Field indices MUST match varnamesl array order above!
                if (ifld ==1) histit(c,lev) = cellsand_col(c,lev)
                if (ifld ==2) histit(c,lev) = cellclay_col(c,lev)
             end do
          end do
          if (tape(t)%dov2xy) then
             histot(:,:) = spval
             call c2g(bounds, nlevsoi, &
                  histit(bounds%begc:bounds%endc, :), &
                  histot(bounds%begg:bounds%endg, :), &
                  c2l_scale_type='unity', l2g_scale_type='veg')
             if (ldomain%isgrid2d) then
                call ncd_io(varname=trim(varnamest(ifld)), dim1name=grlnd, &
                     data=histot, ncid=nfid(t), flag='write')
             else
                call ncd_io(varname=trim(varnamest(ifld)), dim1name=grlnd, &
                     data=histot, ncid=nfid(t), flag='write')
             end if
          else
             call ncd_io(varname=trim(varnamest(ifld)), dim1name=namec,  &
                  data=histit, ncid=nfid(t), flag='write')
          end if
       end do

       if (tape(t)%dov2xy) deallocate(histot)
       deallocate(histit)

    end if  ! (define/write mode

  end subroutine htape_timeconst3D

  !-----------------------------------------------------------------------
  subroutine htape_timeconst(t, mode)
    !
    ! !DESCRIPTION:
    ! Write time constant values to primary history tape.
    use clm_time_manager, only : get_step_size
    use SoilBiogeochemDecompCascadeConType, only : decomp_method, no_soil_decomp
    ! Issue the required netcdf wrapper calls to define the history file
    ! contents.
    !
    ! !USES:
    use clm_varpar      , only : nlevsoi
    use clm_varctl      , only : use_hillslope
    use clm_varcon      , only : zsoi, zlak, secspday, isecspday, isecsphr, isecspmin, ispval
    use domainMod       , only : ldomain, lon1d, lat1d
    use clm_time_manager, only : get_nstep, get_curr_date, get_curr_time
    use clm_time_manager, only : get_ref_date, get_calendar, NO_LEAP_C, GREGORIAN_C
    use FatesInterfaceTypesMod, only : fates_hdim_levsclass
    use FatesInterfaceTypesMod, only : fates_hdim_pfmap_levscpf
    use FatesInterfaceTypesMod, only : fates_hdim_scmap_levscpf
    use FatesInterfaceTypesMod, only : fates_hdim_levcoage
    use FatesInterfaceTypesMod, only : fates_hdim_pfmap_levcapf
    use FatesInterfaceTypesMod, only : fates_hdim_camap_levcapf
    use FatesInterfaceTypesMod, only : fates_hdim_levage
    use FatesInterfaceTypesMod, only : fates_hdim_levheight
    use FatesInterfaceTypesMod, only : fates_hdim_levpft
    use FatesInterfaceTypesMod, only : fates_hdim_scmap_levscag
    use FatesInterfaceTypesMod, only : fates_hdim_agmap_levscag
    use FatesInterfaceTypesMod, only : fates_hdim_scmap_levscagpft
    use FatesInterfaceTypesMod, only : fates_hdim_agmap_levscagpft
    use FatesInterfaceTypesMod, only : fates_hdim_pftmap_levscagpft
    use FatesInterfaceTypesMod, only : fates_hdim_agmap_levagepft
    use FatesInterfaceTypesMod, only : fates_hdim_pftmap_levagepft
    use FatesInterfaceTypesMod, only : fates_hdim_levfuel
    use FatesInterfaceTypesMod, only : fates_hdim_levdamage
    use FatesInterfaceTypesMod, only : fates_hdim_levcwdsc
    use FatesInterfaceTypesMod, only : fates_hdim_levcan
    use FatesInterfaceTypesMod, only : fates_hdim_levleaf
    use FatesInterfaceTypesMod, only : fates_hdim_canmap_levcnlf
    use FatesInterfaceTypesMod, only : fates_hdim_lfmap_levcnlf
    use FatesInterfaceTypesMod, only : fates_hdim_canmap_levcnlfpf
    use FatesInterfaceTypesMod, only : fates_hdim_lfmap_levcnlfpf
    use FatesInterfaceTypesMod, only : fates_hdim_pftmap_levcnlfpf
    use FatesInterfaceTypesMod, only : fates_hdim_levelem
    use FatesInterfaceTypesMod, only : fates_hdim_elmap_levelpft
    use FatesInterfaceTypesMod, only : fates_hdim_pftmap_levelpft
    use FatesInterfaceTypesMod, only : fates_hdim_elmap_levelcwd
    use FatesInterfaceTypesMod, only : fates_hdim_cwdmap_levelcwd
    use FatesInterfaceTypesMod, only : fates_hdim_elmap_levelage
    use FatesInterfaceTypesMod, only : fates_hdim_agemap_levelage
    use FatesInterfaceTypesMod, only : fates_hdim_agmap_levagefuel
    use FatesInterfaceTypesMod, only : fates_hdim_fscmap_levagefuel
    use FatesInterfaceTypesMod, only : fates_hdim_scmap_levcdsc
    use FatesInterfaceTypesMod, only : fates_hdim_cdmap_levcdsc
    use FatesInterfaceTypesMod, only : fates_hdim_scmap_levcdpf
    use FatesInterfaceTypesMod, only : fates_hdim_cdmap_levcdpf
    use FatesInterfaceTypesMod, only : fates_hdim_pftmap_levcdpf
    use FatesInterfaceTypesMod, only : fates_hdim_levlanduse


    !
    ! !ARGUMENTS:
    integer, intent(in) :: t              ! tape index
    integer :: dtime                      ! timestep size
    character(len=*), intent(in) :: mode  ! 'define' or 'write'
    !
    integer :: sec_hist_nhtfrq            ! hist_nhtfrq converted to seconds
    ! !LOCAL VARIABLES:
    integer :: vid,n,i,j,m,c              ! indices
    integer :: nstep                      ! current step
    integer :: mcsec                      ! seconds of current date
    integer :: mdcur                      ! current day
    integer :: mscur                      ! seconds of current day
    integer :: mcdate                     ! current date
    integer :: yr,mon,day,nbsec           ! year,month,day,seconds components of a date
    integer :: hours,minutes,secs         ! hours,minutes,seconds of hh:mm:ss
    character(len= 12) :: step_or_bounds  ! string used in long_name of several time variables
    character(len= 10) :: basedate        ! base date (yyyymmdd)
    character(len=  8) :: basesec         ! base seconds
    character(len=  8) :: cdate           ! system date
    character(len=  8) :: ctime           ! system time
    real(r8):: time                       ! current time
    real(r8):: timedata(2)                ! time interval boundaries
    integer :: dim1id(1)                  ! netCDF dimension id
    integer :: dim2id(2)                  ! netCDF dimension id
    integer :: varid                      ! netCDF variable id
    character(len=max_chars) :: long_name ! variable long name
    character(len=max_namlen):: varname   ! variable name
    character(len=max_namlen):: units     ! variable units
    character(len=max_namlen):: cal       ! calendar from the time-manager
    character(len=max_namlen):: caldesc   ! calendar description to put on file
    character(len=256):: str              ! global attribute string
    real(r8), pointer :: histo(:,:)       ! temporary
    integer :: status
    real(r8) :: zsoi_1d(1)
    type(bounds_type) :: bounds
    integer :: ier                        ! error status
    integer, pointer :: icarr(:)          ! temporary
    character(len=*),parameter :: subname = 'htape_timeconst'
    !-----------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    !***     Time constant grid variables only on first time-sample of file ***
    !-------------------------------------------------------------------------------

    call get_proc_bounds(bounds)


    if (tape(t)%ntimes == 1) then
       if (mode == 'define') then
          call ncd_defvar(varname='levgrnd', xtype=tape(t)%ncprec, &
               dim1name='levgrnd', &
               long_name='coordinate ground levels', units='m', ncid=nfid(t))
          call ncd_defvar(varname='levsoi', xtype=tape(t)%ncprec, &
               dim1name='levsoi', &
               long_name='coordinate soil levels (equivalent to top nlevsoi levels of levgrnd)', units='m', ncid=nfid(t))
          call ncd_defvar(varname='levlak', xtype=tape(t)%ncprec, &
               dim1name='levlak', &
               long_name='coordinate lake levels', units='m', ncid=nfid(t))
          call ncd_defvar(varname='levdcmp', xtype=tape(t)%ncprec, dim1name='levdcmp', &
               long_name='coordinate levels for soil decomposition variables', units='m', ncid=nfid(t))

          if (use_hillslope .and. .not.tape(t)%dov2xy)then
             call ncd_defvar(varname='hillslope_distance', xtype=ncd_double, &
                  dim1name=namec, long_name='hillslope column distance', &
                  units='m', ncid=nfid(t))             
             call ncd_defvar(varname='hillslope_width', xtype=ncd_double, &
                  dim1name=namec, long_name='hillslope column width', &
                  units='m', ncid=nfid(t))             
             call ncd_defvar(varname='hillslope_area', xtype=ncd_double, &
                  dim1name=namec, long_name='hillslope column area', &
                  units='m2', ncid=nfid(t))
             call ncd_defvar(varname='hillslope_elev', xtype=ncd_double, &
                  dim1name=namec, long_name='hillslope column elevation', &
                  units='m', ncid=nfid(t))             
             call ncd_defvar(varname='hillslope_slope', xtype=ncd_double, &
                  dim1name=namec, long_name='hillslope column slope', &
                  units='m/m', ncid=nfid(t))
             call ncd_defvar(varname='hillslope_aspect', xtype=ncd_double, &
                  dim1name=namec, long_name='hillslope column aspect', &
                  units='radians', ncid=nfid(t))
             call ncd_defvar(varname='hillslope_index', xtype=ncd_int, &
                  dim1name=namec, long_name='hillslope index', &
                  ncid=nfid(t))             
             call ncd_defvar(varname='hillslope_cold', xtype=ncd_int, &
                  dim1name=namec, long_name='hillslope downhill column index', &
                  ncid=nfid(t))             
             call ncd_defvar(varname='hillslope_colu', xtype=ncd_int, &
                  dim1name=namec, long_name='hillslope uphill column index', &
                  ncid=nfid(t))             
          end if

          if(use_fates)then

             call ncd_defvar(varname='fates_levscls', xtype=tape(t)%ncprec, dim1name='fates_levscls', &
                  long_name='FATES diameter size class lower bound', units='cm', ncid=nfid(t))
             call ncd_defvar(varname='fates_scmap_levscag', xtype=ncd_int, dim1name='fates_levscag', &
                   long_name='FATES size-class map into size x patch age', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_agmap_levscag', xtype=ncd_int, dim1name='fates_levscag', &
                   long_name='FATES age-class map into size x patch age', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_pftmap_levscpf',xtype=ncd_int, dim1name='fates_levscpf', &
                  long_name='FATES pft index of the combined pft-size class dimension', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_scmap_levscpf',xtype=ncd_int, dim1name='fates_levscpf', &
                  long_name='FATES size index of the combined pft-size class dimension', units='-', ncid=nfid(t))
             ! Units are dash here with units of yr added to the long name so
             ! that postprocessors (like ferret) won't get confused with what
             ! the time coordinate is. EBK Nov/3/2021 (see #1540)
             call ncd_defvar(varname='fates_levcacls', xtype=tape(t)%ncprec, dim1name='fates_levcacls', &
                  long_name='FATES cohort age class lower bound (yr)', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_pftmap_levcapf',xtype=ncd_int, dim1name='fates_levcapf', &
                  long_name='FATES pft index of the combined pft-cohort age class dimension', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_camap_levcapf',xtype=ncd_int, dim1name='fates_levcapf', &
                  long_name='FATES cohort age index of the combined pft-cohort age dimension', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_levage',xtype=tape(t)%ncprec, dim1name='fates_levage', &
                  long_name='FATES patch age (yr)', ncid=nfid(t))
             call ncd_defvar(varname='fates_levheight',xtype=tape(t)%ncprec, dim1name='fates_levheight', &
                  long_name='FATES height (m)', ncid=nfid(t))
             call ncd_defvar(varname='fates_levpft',xtype=ncd_int, dim1name='fates_levpft', &
                  long_name='FATES pft number', ncid=nfid(t))
             call ncd_defvar(varname='fates_levfuel',xtype=ncd_int, dim1name='fates_levfuel', &
                  long_name='FATES fuel index', ncid=nfid(t))
             call ncd_defvar(varname='fates_levcwdsc',xtype=ncd_int, dim1name='fates_levcwdsc', &
                  long_name='FATES cwd size class', ncid=nfid(t))
             call ncd_defvar(varname='fates_levcan',xtype=ncd_int, dim1name='fates_levcan', &
                  long_name='FATES canopy level', ncid=nfid(t))
             call ncd_defvar(varname='fates_levleaf',xtype=ncd_int, dim1name='fates_levleaf', &
                  long_name='FATES leaf+stem level', units='VAI', ncid=nfid(t))
             call ncd_defvar(varname='fates_canmap_levcnlf',xtype=ncd_int, dim1name='fates_levcnlf', &
                  long_name='FATES canopy level of combined canopy-leaf dimension', ncid=nfid(t))
             call ncd_defvar(varname='fates_lfmap_levcnlf',xtype=ncd_int, dim1name='fates_levcnlf', &
                  long_name='FATES leaf level of combined canopy-leaf dimension', ncid=nfid(t))
             call ncd_defvar(varname='fates_canmap_levcnlfpf',xtype=ncd_int, dim1name='fates_levcnlfpf', &
                  long_name='FATES canopy level of combined canopy x leaf x pft dimension', ncid=nfid(t))
             call ncd_defvar(varname='fates_lfmap_levcnlfpf',xtype=ncd_int, dim1name='fates_levcnlfpf', &
                  long_name='FATES leaf level of combined canopy x leaf x pft dimension', ncid=nfid(t))
             call ncd_defvar(varname='fates_pftmap_levcnlfpf',xtype=ncd_int, dim1name='fates_levcnlfpf', &
                  long_name='FATES PFT level of combined canopy x leaf x pft dimension', ncid=nfid(t))
             call ncd_defvar(varname='fates_scmap_levscagpft', xtype=ncd_int, dim1name='fates_levscagpf', &
                  long_name='FATES size-class map into size x patch age x pft', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_agmap_levscagpft', xtype=ncd_int, dim1name='fates_levscagpf', &
                   long_name='FATES age-class map into size x patch age x pft', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_pftmap_levscagpft', xtype=ncd_int, dim1name='fates_levscagpf', &
                   long_name='FATES pft map into size x patch age x pft', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_pftmap_levagepft', xtype=ncd_int, dim1name='fates_levagepft', &
                   long_name='FATES pft map into patch age x pft', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_agmap_levagepft', xtype=ncd_int, dim1name='fates_levagepft', &
                   long_name='FATES age-class map into patch age x pft', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_agmap_levagefuel', xtype=ncd_int, dim1name='fates_levagefuel', &
                   long_name='FATES age-class map into patch age x fuel size', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_fscmap_levagefuel', xtype=ncd_int, dim1name='fates_levagefuel', &
                   long_name='FATES fuel size-class map into patch age x fuel size', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_cdmap_levcdsc',xtype=ncd_int, dim1name='fates_levcdsc', &
                  long_name='FATES damage index of the combined damage-size dimension', ncid=nfid(t))
             call ncd_defvar(varname='fates_scmap_levcdsc',xtype=ncd_int, dim1name='fates_levcdsc', &
                  long_name='FATES size index of the combined damage-size dimension', ncid=nfid(t))
             call ncd_defvar(varname='fates_cdmap_levcdpf',xtype=ncd_int, dim1name='fates_levcdpf', &
                  long_name='FATES damage index of the combined damage-size-PFT dimension', ncid=nfid(t))
             call ncd_defvar(varname='fates_scmap_levcdpf',xtype=ncd_int, dim1name='fates_levcdpf', &
                  long_name='FATES size index of the combined damage-size-PFT dimension', ncid=nfid(t))
             call ncd_defvar(varname='fates_pftmap_levcdpf',xtype=ncd_int, dim1name='fates_levcdpf', &
                  long_name='FATES pft index of the combined damage-size-PFT dimension', ncid=nfid(t))
             call ncd_defvar(varname='fates_levcdam', xtype=tape(t)%ncprec, dim1name='fates_levcdam', &
                  long_name='FATES damage class lower bound', units='unitless', ncid=nfid(t))
             call ncd_defvar(varname='fates_levlanduse',xtype=ncd_int, dim1name='fates_levlanduse', &
                   long_name='FATES land use label', ncid=nfid(t))

          end if


       elseif (mode == 'write') then
          if ( masterproc ) write(iulog, *) ' zsoi:',zsoi
          call ncd_io(varname='levgrnd', data=zsoi, ncid=nfid(t), flag='write')
          call ncd_io(varname='levsoi', data=zsoi(1:nlevsoi), ncid=nfid(t), flag='write')
          call ncd_io(varname='levlak' , data=zlak, ncid=nfid(t), flag='write')
          if ( decomp_method /= no_soil_decomp )then
             call ncd_io(varname='levdcmp', data=zsoi, ncid=nfid(t), flag='write')
          else
             zsoi_1d(1) = 1._r8
             call ncd_io(varname='levdcmp', data=zsoi_1d, ncid=nfid(t), flag='write')
          end if

          if (use_hillslope .and. .not.tape(t)%dov2xy) then
             call ncd_io(varname='hillslope_distance' , data=col%hill_distance, dim1name=namec, ncid=nfid(t), flag='write')
             call ncd_io(varname='hillslope_width' , data=col%hill_width, dim1name=namec, ncid=nfid(t), flag='write')
             call ncd_io(varname='hillslope_area' , data=col%hill_area, dim1name=namec, ncid=nfid(t), flag='write')
             call ncd_io(varname='hillslope_elev' , data=col%hill_elev, dim1name=namec, ncid=nfid(t), flag='write')
             call ncd_io(varname='hillslope_slope' , data=col%hill_slope, dim1name=namec, ncid=nfid(t), flag='write')
             call ncd_io(varname='hillslope_aspect' , data=col%hill_aspect, dim1name=namec, ncid=nfid(t), flag='write')
             call ncd_io(varname='hillslope_index' , data=col%hillslope_ndx, dim1name=namec, ncid=nfid(t), flag='write')

             ! write global indices rather than local indices
             allocate(icarr(bounds%begc:bounds%endc),stat=ier)
             if (ier /= 0) then
                call endrun(msg=' allocation error of icarr'//errMsg(sourcefile, __LINE__))
             end if

             do c = bounds%begc,bounds%endc
                if (col%cold(c) /= ispval) then 
                   icarr(c)= get_global_index(subgrid_index=col%cold(c), subgrid_level=subgrid_level_column)
                else
                   icarr(c)= col%cold(c)
                endif
             enddo
             
             call ncd_io(varname='hillslope_cold' , data=icarr, dim1name=namec, ncid=nfid(t), flag='write')

             do c = bounds%begc,bounds%endc
                if (col%colu(c) /= ispval) then 
                   icarr(c)= get_global_index(subgrid_index=col%colu(c), subgrid_level=subgrid_level_column)
                else
                   icarr(c)= col%colu(c)
                endif
             enddo

             call ncd_io(varname='hillslope_colu' , data=icarr, dim1name=namec, ncid=nfid(t), flag='write')
             deallocate(icarr)
          endif

          if(use_fates)then
             call ncd_io(varname='fates_scmap_levscag',data=fates_hdim_scmap_levscag, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_agmap_levscag',data=fates_hdim_agmap_levscag, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levscls',data=fates_hdim_levsclass, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levcacls',data=fates_hdim_levcoage, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_pftmap_levscpf',data=fates_hdim_pfmap_levscpf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_scmap_levscpf',data=fates_hdim_scmap_levscpf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_pftmap_levcapf',data=fates_hdim_pfmap_levcapf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_camap_levcapf',data=fates_hdim_camap_levcapf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levage',data=fates_hdim_levage, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levheight',data=fates_hdim_levheight, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levpft',data=fates_hdim_levpft, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levfuel',data=fates_hdim_levfuel, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levcdam',data=fates_hdim_levdamage, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levcwdsc',data=fates_hdim_levcwdsc, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levcan',data=fates_hdim_levcan, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levleaf',data=fates_hdim_levleaf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_canmap_levcnlf',data=fates_hdim_canmap_levcnlf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_lfmap_levcnlf',data=fates_hdim_lfmap_levcnlf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_canmap_levcnlfpf',data=fates_hdim_canmap_levcnlfpf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_lfmap_levcnlfpf',data=fates_hdim_lfmap_levcnlfpf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_pftmap_levcnlfpf',data=fates_hdim_pftmap_levcnlfpf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_scmap_levscagpft',data=fates_hdim_scmap_levscagpft, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_agmap_levscagpft',data=fates_hdim_agmap_levscagpft, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_pftmap_levscagpft',data=fates_hdim_pftmap_levscagpft, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_pftmap_levagepft',data=fates_hdim_pftmap_levagepft, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_agmap_levagepft',data=fates_hdim_agmap_levagepft, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_agmap_levagefuel',data=fates_hdim_agmap_levagefuel, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_fscmap_levagefuel',data=fates_hdim_fscmap_levagefuel, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_scmap_levcdsc',data=fates_hdim_scmap_levcdsc, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_cdmap_levcdsc',data=fates_hdim_cdmap_levcdsc, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_scmap_levcdpf',data=fates_hdim_scmap_levcdpf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_cdmap_levcdpf',data=fates_hdim_cdmap_levcdpf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_pftmap_levcdpf',data=fates_hdim_pftmap_levcdpf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levlanduse',data=fates_hdim_levlanduse, ncid=nfid(t), flag='write')
          end if

       endif
    endif

    !-------------------------------------------------------------------------------
    !***     Time definition variables ***
    !-------------------------------------------------------------------------------

    ! For define mode -- only do this for first time-sample
    if (mode == 'define' .and. tape(t)%ntimes == 1) then
       call get_ref_date(yr, mon, day, nbsec)
       nstep = get_nstep()
       hours   = nbsec / 3600
       minutes = (nbsec - hours*3600) / 60
       secs    = (nbsec - hours*3600 - minutes*60)
       write(basedate,80) yr,mon,day
80     format(i4.4,'-',i2.2,'-',i2.2)
       write(basesec ,90) hours, minutes, secs
90     format(i2.2,':',i2.2,':',i2.2)

       dim1id(1) = time_dimid
       str = 'days since ' // basedate // " " // basesec
       if (hist_avgflag_pertape(t) /= 'I') then  ! NOT instantaneous fields tape
          step_or_bounds = 'time_bounds'
          long_name = 'time at exact middle of ' // step_or_bounds
          call ncd_defvar(nfid(t), 'time', tape(t)%ncprec, 1, dim1id, varid, &
               long_name=long_name, units=str)
          call ncd_putatt(nfid(t), varid, 'bounds', 'time_bounds')
       else  ! instantaneous fields tape
          step_or_bounds = 'time step'
          long_name = 'time at end of ' // step_or_bounds
          call ncd_defvar(nfid(t), 'time', tape(t)%ncprec, 1, dim1id, varid, &
               long_name=long_name, units=str)
       end if
       cal = get_calendar()
       if (      trim(cal) == NO_LEAP_C   )then
          caldesc = "noleap"
       else if ( trim(cal) == GREGORIAN_C )then
          caldesc = "gregorian"
       end if
       call ncd_putatt(nfid(t), varid, 'calendar', caldesc)

       dim1id(1) = time_dimid
       long_name = 'current date (YYYYMMDD) at end of ' // step_or_bounds
       call ncd_defvar(nfid(t) , 'mcdate', ncd_int, 1, dim1id , varid, &
          long_name = long_name)
       call ncd_putatt(nfid(t), varid, 'calendar', caldesc)
       !
       ! add global attribute time_period_freq
       !
       if (hist_nhtfrq(t) < 0) then !hour need to convert to seconds
          sec_hist_nhtfrq = abs(hist_nhtfrq(t))*3600
       else
          sec_hist_nhtfrq = hist_nhtfrq(t)
       end if

       dtime = get_step_size()
       if (sec_hist_nhtfrq == 0) then !month
          time_period_freq = 'month_1'
       else if (mod(sec_hist_nhtfrq*dtime,isecspday) == 0) then ! day
          write(time_period_freq,999) 'day_',sec_hist_nhtfrq*dtime/isecspday
       else if (mod(sec_hist_nhtfrq*dtime,isecsphr) == 0) then ! hour
          write(time_period_freq,999) 'hour_',(sec_hist_nhtfrq*dtime)/isecsphr
       else if (mod(sec_hist_nhtfrq*dtime,isecspmin) == 0) then ! minute
          write(time_period_freq,999) 'minute_',(sec_hist_nhtfrq*dtime)/isecspmin
       else                     ! second
          write(time_period_freq,999) 'second_',sec_hist_nhtfrq*dtime
       end if
999    format(a,i0)

       call ncd_putatt(nfid(t), ncd_global, 'time_period_freq',          &
                          trim(time_period_freq))

       long_name = 'current seconds of current date at end of ' // step_or_bounds
       call ncd_defvar(nfid(t) , 'mcsec' , ncd_int, 1, dim1id , varid, &
          long_name = long_name, units='s')
       call ncd_putatt(nfid(t), varid, 'calendar', caldesc)
       long_name = 'current day (from base day) at end of ' // step_or_bounds
       call ncd_defvar(nfid(t) , 'mdcur' , ncd_int, 1, dim1id , varid, &
          long_name = long_name)
       call ncd_putatt(nfid(t), varid, 'calendar', caldesc)
       long_name = 'current seconds of current day at end of ' // step_or_bounds
       call ncd_defvar(nfid(t) , 'mscur' , ncd_int, 1, dim1id , varid, &
          long_name = long_name)
       call ncd_putatt(nfid(t), varid, 'calendar', caldesc)
       call ncd_defvar(nfid(t) , 'nstep' , ncd_int, 1, dim1id , varid, &
          long_name = 'time step')

       dim2id(1) = nbnd_dimid;  dim2id(2) = time_dimid
       if (hist_avgflag_pertape(t) /= 'I') then  ! NOT instantaneous fields tape
          call ncd_defvar(nfid(t), 'time_bounds', ncd_double, 2, dim2id, varid, &
             long_name = 'time interval endpoints', &
             units = str)
          call ncd_putatt(nfid(t), varid, 'calendar', caldesc)
       end if

       dim2id(1) = strlen_dimid;  dim2id(2) = time_dimid
       call ncd_defvar(nfid(t), 'date_written', ncd_char, 2, dim2id, varid)
       call ncd_defvar(nfid(t), 'time_written', ncd_char, 2, dim2id, varid)

       if ( len_trim(TimeConst3DVars_Filename) > 0 )then
          call ncd_putatt(nfid(t), ncd_global, 'Time_constant_3Dvars_filename', &
                          trim(TimeConst3DVars_Filename))
       end if
       if ( len_trim(TimeConst3DVars)          > 0 )then
          call ncd_putatt(nfid(t), ncd_global, 'Time_constant_3Dvars',          &
                          trim(TimeConst3DVars))
       end if

    elseif (mode == 'write') then

       call get_curr_time (mdcur, mscur)
       call get_curr_date (yr, mon, day, mcsec)
       mcdate = yr*10000 + mon*100 + day
       nstep = get_nstep()

       call ncd_io('mcdate', mcdate, 'write', nfid(t), nt=tape(t)%ntimes)
       call ncd_io('mcsec' , mcsec , 'write', nfid(t), nt=tape(t)%ntimes)
       call ncd_io('mdcur' , mdcur , 'write', nfid(t), nt=tape(t)%ntimes)
       call ncd_io('mscur' , mscur , 'write', nfid(t), nt=tape(t)%ntimes)
       call ncd_io('nstep' , nstep , 'write', nfid(t), nt=tape(t)%ntimes)

       timedata(1) = tape(t)%begtime  ! beginning time
       timedata(2) = mdcur + mscur/secspday  ! end time
       if (hist_avgflag_pertape(t) /= 'I') then  ! NOT instantaneous fields tape
          time = (timedata(1) + timedata(2)) * 0.5_r8
          call ncd_io('time_bounds', timedata, 'write', nfid(t), nt=tape(t)%ntimes)
       else
          time = timedata(2)
       end if
       call ncd_io('time'  , time  , 'write', nfid(t), nt=tape(t)%ntimes)

       call getdatetime (cdate, ctime)
       call ncd_io('date_written', cdate, 'write', nfid(t), nt=tape(t)%ntimes)

       call ncd_io('time_written', ctime, 'write', nfid(t), nt=tape(t)%ntimes)

    endif

    !-------------------------------------------------------------------------------
    !***     Grid definition variables ***
    !-------------------------------------------------------------------------------
    ! For define mode -- only do this for first time-sample
    if (mode == 'define' .and. tape(t)%ntimes == 1) then

       if (ldomain%isgrid2d) then
          call ncd_defvar(varname='lon', xtype=tape(t)%ncprec, dim1name='lon', &
              long_name='coordinate longitude', units='degrees_east', &
              ncid=nfid(t), missing_value=spval, fill_value=spval)
       else
          call ncd_defvar(varname='lon', xtype=tape(t)%ncprec, &
              dim1name=grlnd, &
              long_name='coordinate longitude', units='degrees_east', ncid=nfid(t), &
              missing_value=spval, fill_value=spval)
       end if
       if (ldomain%isgrid2d) then
          call ncd_defvar(varname='lat', xtype=tape(t)%ncprec, dim1name='lat', &
              long_name='coordinate latitude', units='degrees_north', &
              ncid=nfid(t), missing_value=spval, fill_value=spval)
       else
          call ncd_defvar(varname='lat', xtype=tape(t)%ncprec, &
              dim1name=grlnd, &
              long_name='coordinate latitude', units='degrees_north', ncid=nfid(t), &
              missing_value=spval, fill_value=spval)
       end if
       if (ldomain%isgrid2d) then
          call ncd_defvar(varname='area', xtype=tape(t)%ncprec, &
              dim1name='lon', dim2name='lat',&
              long_name='grid cell areas', units='km^2', ncid=nfid(t), &
              missing_value=spval, fill_value=spval)
       else
          call ncd_defvar(varname='area', xtype=tape(t)%ncprec, &
              dim1name=grlnd, &
              long_name='grid cell areas', units='km^2', ncid=nfid(t), &
              missing_value=spval, fill_value=spval)
       end if
       if (ldomain%isgrid2d) then
          call ncd_defvar(varname='landfrac', xtype=tape(t)%ncprec, &
              dim1name='lon', dim2name='lat', &
              long_name='land fraction', ncid=nfid(t), &
              missing_value=spval, fill_value=spval)
       else
          call ncd_defvar(varname='landfrac', xtype=tape(t)%ncprec, &
              dim1name=grlnd, &
              long_name='land fraction', ncid=nfid(t), &
              missing_value=spval, fill_value=spval)
       end if
       if (ldomain%isgrid2d) then
          call ncd_defvar(varname='landmask', xtype=ncd_int, &
              dim1name='lon', dim2name='lat', &
              long_name='land/ocean mask (0.=ocean and 1.=land)', ncid=nfid(t), &
              imissing_value=ispval, ifill_value=ispval)
       else
          call ncd_defvar(varname='landmask', xtype=ncd_int, &
              dim1name=grlnd, &
              long_name='land/ocean mask (0.=ocean and 1.=land)', ncid=nfid(t), &
              imissing_value=ispval, ifill_value=ispval)
       end if
       if (ldomain%isgrid2d) then
          call ncd_defvar(varname='nbedrock' , xtype=ncd_int, &
              dim1name='lon', dim2name='lat', &
              long_name='index of shallowest bedrock layer', ncid=nfid(t), &
              imissing_value=ispval, ifill_value=ispval)
       else
          call ncd_defvar(varname='nbedrock' , xtype=ncd_int, &
              dim1name=grlnd, &
              long_name='index of shallowest bedrock layer', ncid=nfid(t), &
              imissing_value=ispval, ifill_value=ispval)
       end if

    else if (mode == 'write') then

       ! Most of this is constant and only needs to be done on tape(t)%ntimes=1
       ! But, some may change for dynamic PATCH mode for example

       if (ldomain%isgrid2d) then
          call ncd_io(varname='lon', data=lon1d, ncid=nfid(t), flag='write')
          call ncd_io(varname='lat', data=lat1d, ncid=nfid(t), flag='write')
       else
          call ncd_io(varname='lon', data=ldomain%lonc, dim1name=grlnd, ncid=nfid(t), flag='write')
          call ncd_io(varname='lat', data=ldomain%latc, dim1name=grlnd, ncid=nfid(t), flag='write')
       end if
       call ncd_io(varname='area'    , data=ldomain%area, dim1name=grlnd, ncid=nfid(t), flag='write')
       call ncd_io(varname='landfrac', data=ldomain%frac, dim1name=grlnd, ncid=nfid(t), flag='write')
       call ncd_io(varname='landmask', data=ldomain%mask, dim1name=grlnd, ncid=nfid(t), flag='write')
       call ncd_io(varname='nbedrock' , data=grc%nbedrock, dim1name=grlnd, ncid=nfid(t), flag='write')

    end if  ! (define/write mode

  end subroutine htape_timeconst

  !-----------------------------------------------------------------------
  subroutine hfields_write(t, mode)
    !
    ! !DESCRIPTION:
    ! Write history tape.  Issue the call to write the variable.
    !
    ! !USES:
    use domainMod , only : ldomain
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t                ! tape index
    character(len=*), intent(in) :: mode    ! 'define' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: f                         ! field index
    integer :: k                         ! 1d index
    integer :: c,l,p                     ! indices
    integer :: beg1d                     ! on-node 1d field pointer start index
    integer :: end1d                     ! on-node 1d field pointer end index
    integer :: beg1d_out                 ! on-node 1d hbuf pointer start index
    integer :: end1d_out                 ! on-node 1d hbuf pointer end index
    integer :: num1d_out                 ! size of hbuf first dimension (overall all nodes)
    integer :: num2d                     ! hbuf second dimension size
    integer :: nt                        ! time index
    integer :: ier                       ! error status
    integer :: numdims                   ! number of dimensions
    integer :: varid                     ! variable id
    character(len=avgflag_strlen) :: avgflag  ! time averaging flag
    character(len=max_chars) :: long_name! long name
    character(len=max_chars) :: units    ! units
    character(len=max_namlen):: varname  ! variable name
    character(len=32) :: avgstr          ! time averaging type
    character(len=hist_dim_name_length)  :: type1d          ! field 1d type
    character(len=hist_dim_name_length)  :: type1d_out      ! history output 1d type
    character(len=hist_dim_name_length)  :: type2d          ! history output 2d type
    character(len=scale_type_strlen)     :: l2g_scale_type
    character(len=32) :: dim1name        ! temporary
    character(len=32) :: dim2name        ! temporary
    real(r8), pointer :: histo(:,:)      ! temporary
    real(r8), pointer :: hist1do(:)      ! temporary
    character(len=*),parameter :: subname = 'hfields_write'
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

       varname        = tape(t)%hlist(f)%field%name
       long_name      = tape(t)%hlist(f)%field%long_name
       units          = tape(t)%hlist(f)%field%units
       avgflag        = tape(t)%hlist(f)%avgflag
       type1d         = tape(t)%hlist(f)%field%type1d
       type1d_out     = tape(t)%hlist(f)%field%type1d_out
       beg1d          = tape(t)%hlist(f)%field%beg1d
       end1d          = tape(t)%hlist(f)%field%end1d
       beg1d_out      = tape(t)%hlist(f)%field%beg1d_out
       end1d_out      = tape(t)%hlist(f)%field%end1d_out
       num1d_out      = tape(t)%hlist(f)%field%num1d_out
       type2d         = tape(t)%hlist(f)%field%type2d
       numdims        = tape(t)%hlist(f)%field%numdims
       num2d          = tape(t)%hlist(f)%field%num2d
       l2g_scale_type = tape(t)%hlist(f)%field%l2g_scale_type
       nt             = tape(t)%ntimes

       if (mode == 'define') then

          select case (avgflag(1:1))
          case ('A')
             avgstr = 'mean'
          case ('I')
             avgstr = 'instantaneous'
          case ('X')
             avgstr = 'maximum'
          case ('M')
             avgstr = 'minimum'
          case ('S')
             avgstr = 'sum'
          case ('L')
             avgstr = 'local solar time'
          case default
             write(iulog,*) trim(subname),' ERROR: unknown time averaging flag (avgflag)=',avgflag
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end select

          if (type1d_out == grlnd) then
             if (ldomain%isgrid2d) then
                dim1name = 'lon'      ; dim2name = 'lat'
             else
                dim1name = trim(grlnd); dim2name = 'undefined'
             end if
          else
             dim1name = type1d_out ; dim2name = 'undefined'
          endif

          if (dim2name == 'undefined') then
             if (numdims == 1) then
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name=dim1name, dim2name='time', &
                     long_name=long_name, units=units, cell_method=avgstr, &
                     missing_value=spval, fill_value=spval, &
                     varid=varid)
             else
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name=dim1name, dim2name=type2d, dim3name='time', &
                     long_name=long_name, units=units, cell_method=avgstr, &
                     missing_value=spval, fill_value=spval, &
                     varid=varid)
             end if
          else
             if (numdims == 1) then
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name=dim1name, dim2name=dim2name, dim3name='time', &
                     long_name=long_name, units=units, cell_method=avgstr, &
                     missing_value=spval, fill_value=spval, &
                     varid=varid)
             else
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name=dim1name, dim2name=dim2name, dim3name=type2d, dim4name='time', &
                     long_name=long_name, units=units, cell_method=avgstr, &
                     missing_value=spval, fill_value=spval, &
                     varid=varid)
             end if
          endif

          if (type1d_out == nameg .or. type1d_out == grlnd) then
             call add_landunit_mask_metadata(nfid(t), varid, l2g_scale_type)
          end if

       else if (mode == 'write') then

          ! Determine output buffer

          histo => tape(t)%hlist(f)%hbuf

          ! Allocate dynamic memory

          if (numdims == 1) then
             allocate(hist1do(beg1d_out:end1d_out), stat=ier)
             if (ier /= 0) then
                write(iulog,*) trim(subname),' ERROR: allocation'
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if
             hist1do(beg1d_out:end1d_out) = histo(beg1d_out:end1d_out,1)
          end if

          ! Write history output.  Always output land and ocean runoff on xy grid.

          if (numdims == 1) then
             call ncd_io(flag='write', varname=varname, &
                  dim1name=type1d_out, data=hist1do, ncid=nfid(t), nt=nt)
          else
             call ncd_io(flag='write', varname=varname, &
                  dim1name=type1d_out, data=histo, ncid=nfid(t), nt=nt)
          end if


          ! Deallocate dynamic memory

          if (numdims == 1) then
             deallocate(hist1do)
          end if

       end if

    end do

  end subroutine hfields_write

  !-----------------------------------------------------------------------
  subroutine hfields_1dinfo(t, mode)
    !
    ! !DESCRIPTION:
    ! Write/define 1d info for history tape.
    !
    ! !USES:
    use decompMod   , only : gindex_global
    use domainMod   , only : ldomain, ldomain
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t                ! tape index
    character(len=*), intent(in) :: mode    ! 'define' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: f                         ! field index
    integer :: k                         ! 1d index
    integer :: g,c,l,p                   ! indices
    integer :: ier                       ! errir status
    integer :: gindex                    ! global gridcell index
    real(r8), pointer :: rgarr(:)        ! temporary
    real(r8), pointer :: rcarr(:)        ! temporary
    real(r8), pointer :: rlarr(:)        ! temporary
    real(r8), pointer :: rparr(:)        ! temporary
    integer , pointer :: igarr(:)        ! temporary
    integer , pointer :: icarr(:)        ! temporary
    integer , pointer :: ilarr(:)        ! temporary
    integer , pointer :: iparr(:)        ! temporary
    type(file_desc_t), pointer :: ncid   ! netcdf file
    type(bounds_type) :: bounds
    character(len=*),parameter :: subname = 'hfields_1dinfo'
!-----------------------------------------------------------------------

    call get_proc_bounds(bounds)

    ncid => nfid(t)

    if (mode == 'define') then

          ! Define gridcell info

          call ncd_defvar(varname='grid1d_lon', xtype=ncd_double, dim1name=nameg, &
               long_name='gridcell longitude', units='degrees_east', fill_value=spval, ncid=ncid)

          call ncd_defvar(varname='grid1d_lat', xtype=ncd_double,  dim1name=nameg, &
               long_name='gridcell latitude', units='degrees_north', fill_value=spval, ncid=ncid)

          call ncd_defvar(varname='grid1d_ixy', xtype=ncd_int, dim1name=nameg, &
               long_name='2d longitude index of corresponding gridcell', ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='grid1d_jxy', xtype=ncd_int, dim1name=nameg, &
               long_name='2d latitude index of corresponding gridcell', ifill_value=ispval, ncid=ncid)

          ! Define landunit info

          call ncd_defvar(varname='land1d_lon', xtype=ncd_double, dim1name=namel, &
               long_name='landunit longitude', units='degrees_east', fill_value=spval, ncid=ncid)

          call ncd_defvar(varname='land1d_lat', xtype=ncd_double, dim1name=namel, &
               long_name='landunit latitude', units='degrees_north', fill_value=spval, ncid=ncid)

          call ncd_defvar(varname='land1d_ixy', xtype=ncd_int, dim1name=namel, &
               long_name='2d longitude index of corresponding landunit', ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='land1d_jxy', xtype=ncd_int, dim1name=namel, &
               long_name='2d latitude index of corresponding landunit', ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='land1d_gi', xtype=ncd_int, dim1name=namel, &
               long_name='1d grid index of corresponding landunit', ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='land1d_wtgcell', xtype=ncd_double, dim1name=namel, &
               long_name='landunit weight relative to corresponding gridcell', fill_value=spval, ncid=ncid)

          call ncd_defvar(varname='land1d_ityplunit', xtype=ncd_int, dim1name=namel, &
               long_name='landunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)', &
                  ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='land1d_active', xtype=ncd_log, dim1name=namel, &
               long_name='true => do computations on this landunit', ifill_value=0, ncid=ncid)

          ! Define column info

          call ncd_defvar(varname='cols1d_lon', xtype=ncd_double, dim1name=namec, &
               long_name='column longitude', units='degrees_east', fill_value=spval, ncid=ncid)

          call ncd_defvar(varname='cols1d_lat', xtype=ncd_double, dim1name=namec, &
               long_name='column latitude', units='degrees_north', fill_value=spval, ncid=ncid)

          call ncd_defvar(varname='cols1d_ixy', xtype=ncd_int, dim1name=namec, &
               long_name='2d longitude index of corresponding column', ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='cols1d_jxy', xtype=ncd_int, dim1name=namec, &
               long_name='2d latitude index of corresponding column', ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='cols1d_gi', xtype=ncd_int, dim1name=namec, &
               long_name='1d grid index of corresponding column', ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='cols1d_li', xtype=ncd_int, dim1name=namec, &
               long_name='1d landunit index of corresponding column', ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='cols1d_wtgcell', xtype=ncd_double, dim1name=namec, &
               long_name='column weight relative to corresponding gridcell', fill_value=spval, ncid=ncid)

          call ncd_defvar(varname='cols1d_wtlunit', xtype=ncd_double, dim1name=namec, &
               long_name='column weight relative to corresponding landunit', fill_value=spval, ncid=ncid)

          call ncd_defvar(varname='cols1d_itype_col', xtype=ncd_int, dim1name=namec, &
               long_name='column type (see global attributes)', ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='cols1d_itype_lunit', xtype=ncd_int, dim1name=namec, &
               long_name='column landunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)', &
                  ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='cols1d_active', xtype=ncd_log, dim1name=namec, &
               long_name='true => do computations on this column', ifill_value=0, ncid=ncid)

          call ncd_defvar(varname='cols1d_nbedrock', xtype=ncd_int, dim1name=namec, &
               long_name='column bedrock depth index', ifill_value=ispval, ncid=ncid)

          ! Define patch info

          call ncd_defvar(varname='pfts1d_lon', xtype=ncd_double, dim1name=namep, &
               long_name='pft longitude', units='degrees_east', fill_value=spval, ncid=ncid)

          call ncd_defvar(varname='pfts1d_lat', xtype=ncd_double, dim1name=namep, &
               long_name='pft latitude', units='degrees_north', fill_value=spval, ncid=ncid)

          call ncd_defvar(varname='pfts1d_ixy', xtype=ncd_int, dim1name=namep, &
               long_name='2d longitude index of corresponding pft', ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='pfts1d_jxy', xtype=ncd_int, dim1name=namep, &
               long_name='2d latitude index of corresponding pft', ifill_value=ispval,  ncid=ncid)

          call ncd_defvar(varname='pfts1d_gi', xtype=ncd_int, dim1name=namep, &
               long_name='1d grid index of corresponding pft', ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='pfts1d_li', xtype=ncd_int, dim1name=namep, &
               long_name='1d landunit index of corresponding pft', ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='pfts1d_ci', xtype=ncd_int, dim1name=namep, &
               long_name='1d column index of corresponding pft', ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='pfts1d_wtgcell', xtype=ncd_double, dim1name=namep, &
               long_name='pft weight relative to corresponding gridcell', fill_value=spval, ncid=ncid)

          call ncd_defvar(varname='pfts1d_wtlunit', xtype=ncd_double, dim1name=namep, &
               long_name='pft weight relative to corresponding landunit', fill_value=spval, ncid=ncid)

          call ncd_defvar(varname='pfts1d_wtcol', xtype=ncd_double, dim1name=namep, &
               long_name='pft weight relative to corresponding column', fill_value=spval, ncid=ncid)

          call ncd_defvar(varname='pfts1d_itype_veg', xtype=ncd_int, dim1name=namep, &
               long_name='pft vegetation type', ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='pfts1d_itype_col', xtype=ncd_int, dim1name=namep, &
               long_name='pft column type (see global attributes)', ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='pfts1d_itype_lunit', xtype=ncd_int, dim1name=namep, &
               long_name='pft landunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)',  &
                  ifill_value=ispval, ncid=ncid)

          call ncd_defvar(varname='pfts1d_active', xtype=ncd_log, dim1name=namep, &
               ifill_value=0, long_name='true => do computations on this pft', ncid=ncid)

    else if (mode == 'write') then

       ! Determine bounds

       allocate(&
            rgarr(bounds%begg:bounds%endg),&
            rlarr(bounds%begl:bounds%endl),&
            rcarr(bounds%begc:bounds%endc),&
            rparr(bounds%begp:bounds%endp),&
            stat=ier)
       if (ier /= 0) then
          call endrun(msg=' hfields_1dinfo allocation error of rarrs'//errMsg(sourcefile, __LINE__))
       end if

       allocate(&
            igarr(bounds%begg:bounds%endg),&
            ilarr(bounds%begl:bounds%endl),&
            icarr(bounds%begc:bounds%endc),&
            iparr(bounds%begp:bounds%endp),stat=ier)
       if (ier /= 0) then
          call endrun(msg=' hfields_1dinfo allocation error of iarrs'//errMsg(sourcefile, __LINE__))
       end if

       ! Write gridcell info

       call ncd_io(varname='grid1d_lon', data=grc%londeg, dim1name=nameg, ncid=ncid, flag='write')
       call ncd_io(varname='grid1d_lat', data=grc%latdeg, dim1name=nameg, ncid=ncid, flag='write')
       do g = bounds%begg,bounds%endg
         gindex = gindex_global(g-bounds%begg+1)
         igarr(g)= mod(gindex-1,ldomain%ni) + 1
       enddo
       call ncd_io(varname='grid1d_ixy', data=igarr      , dim1name=nameg, ncid=ncid, flag='write')
       do g = bounds%begg,bounds%endg
         gindex = gindex_global(g-bounds%begg+1)
         igarr(g)= (gindex-1)/ldomain%ni + 1
       enddo
       call ncd_io(varname='grid1d_jxy', data=igarr      , dim1name=nameg, ncid=ncid, flag='write')

       ! Write landunit info

       do l = bounds%begl,bounds%endl
         rlarr(l) = grc%londeg(lun%gridcell(l))
       enddo
       call ncd_io(varname='land1d_lon', data=rlarr, dim1name=namel, ncid=ncid, flag='write')
       do l = bounds%begl,bounds%endl
         rlarr(l) = grc%latdeg(lun%gridcell(l))
       enddo
       call ncd_io(varname='land1d_lat', data=rlarr, dim1name=namel, ncid=ncid, flag='write')
       do l= bounds%begl,bounds%endl
         gindex = gindex_global(lun%gridcell(l)-bounds%begg+1)
         ilarr(l) = mod(gindex-1,ldomain%ni) + 1
       enddo
       call ncd_io(varname='land1d_ixy', data=ilarr, dim1name=namel, ncid=ncid, flag='write')
       do l = bounds%begl,bounds%endl
         gindex = gindex_global(lun%gridcell(l)-bounds%begg+1)
         ilarr(l) = (gindex-1)/ldomain%ni + 1
       enddo
       call ncd_io(varname='land1d_jxy'      , data=ilarr        , dim1name=namel, ncid=ncid, flag='write')
       ilarr = get_global_index_array(lun%gridcell(bounds%begl:bounds%endl), bounds%begl, bounds%endl, &
            subgrid_level=subgrid_level_gridcell)
       call ncd_io(varname='land1d_gi'       , data=ilarr, dim1name=namel, ncid=ncid, flag='write')
       call ncd_io(varname='land1d_wtgcell'  , data=lun%wtgcell , dim1name=namel, ncid=ncid, flag='write')
       call ncd_io(varname='land1d_ityplunit', data=lun%itype   , dim1name=namel, ncid=ncid, flag='write')
       call ncd_io(varname='land1d_active'   , data=lun%active  , dim1name=namel, ncid=ncid, flag='write')

       ! Write column info

       do c = bounds%begc,bounds%endc
         rcarr(c) = grc%londeg(col%gridcell(c))
       enddo
       call ncd_io(varname='cols1d_lon', data=rcarr, dim1name=namec, ncid=ncid, flag='write')
       do c = bounds%begc,bounds%endc
         rcarr(c) = grc%latdeg(col%gridcell(c))
       enddo
       call ncd_io(varname='cols1d_lat', data=rcarr, dim1name=namec, ncid=ncid, flag='write')
       do c = bounds%begc,bounds%endc
         gindex = gindex_global(col%gridcell(c)-bounds%begg+1)
         icarr(c) = mod(gindex-1,ldomain%ni) + 1
       enddo
       call ncd_io(varname='cols1d_ixy', data=icarr, dim1name=namec, ncid=ncid, flag='write')
       do c = bounds%begc,bounds%endc
         gindex = gindex_global(col%gridcell(c)-bounds%begg+1)
         icarr(c) = (gindex-1)/ldomain%ni + 1
       enddo
       call ncd_io(varname='cols1d_jxy'    , data=icarr         ,dim1name=namec, ncid=ncid, flag='write')
       icarr = get_global_index_array(col%gridcell(bounds%begc:bounds%endc), bounds%begc, bounds%endc, &
            subgrid_level=subgrid_level_gridcell)
       call ncd_io(varname='cols1d_gi'     , data=icarr, dim1name=namec, ncid=ncid, flag='write')
       icarr = get_global_index_array(col%landunit(bounds%begc:bounds%endc), bounds%begc, bounds%endc, &
            subgrid_level=subgrid_level_landunit)
       call ncd_io(varname='cols1d_li', data=icarr            , dim1name=namec, ncid=ncid, flag='write')

       call ncd_io(varname='cols1d_wtgcell', data=col%wtgcell , dim1name=namec, ncid=ncid, flag='write')
       call ncd_io(varname='cols1d_wtlunit', data=col%wtlunit , dim1name=namec, ncid=ncid, flag='write')
       call ncd_io(varname='cols1d_itype_col', data=col%itype , dim1name=namec, ncid=ncid, flag='write')

       do c = bounds%begc,bounds%endc
         icarr(c) = lun%itype(col%landunit(c))
       enddo
       call ncd_io(varname='cols1d_itype_lunit', data=icarr    , dim1name=namec, ncid=ncid, flag='write')

       call ncd_io(varname='cols1d_active' , data=col%active  , dim1name=namec, ncid=ncid, flag='write')
       call ncd_io(varname='cols1d_nbedrock', data=col%nbedrock , dim1name=namec, ncid=ncid, flag='write')

       ! Write patch info

       do p = bounds%begp,bounds%endp
         rparr(p) = grc%londeg(patch%gridcell(p))
       enddo
       call ncd_io(varname='pfts1d_lon', data=rparr, dim1name=namep, ncid=ncid, flag='write')
       do p = bounds%begp,bounds%endp
         rparr(p) = grc%latdeg(patch%gridcell(p))
       enddo
       call ncd_io(varname='pfts1d_lat', data=rparr, dim1name=namep, ncid=ncid, flag='write')
       do p = bounds%begp,bounds%endp
         gindex = gindex_global(patch%gridcell(p)-bounds%begg+1)
         iparr(p) = mod(gindex-1,ldomain%ni) + 1
       enddo
       call ncd_io(varname='pfts1d_ixy', data=iparr, dim1name=namep, ncid=ncid, flag='write')
       do p = bounds%begp,bounds%endp
         gindex = gindex_global(patch%gridcell(p)-bounds%begg+1)
         iparr(p) = (gindex-1)/ldomain%ni + 1
       enddo
       call ncd_io(varname='pfts1d_jxy'      , data=iparr        , dim1name=namep, ncid=ncid, flag='write')

       iparr = get_global_index_array(patch%gridcell(bounds%begp:bounds%endp), bounds%begp, bounds%endp, &
            subgrid_level=subgrid_level_gridcell)
       call ncd_io(varname='pfts1d_gi'       , data=iparr, dim1name=namep, ncid=ncid, flag='write')
       iparr = get_global_index_array(patch%landunit(bounds%begp:bounds%endp), bounds%begp, bounds%endp, &
            subgrid_level=subgrid_level_landunit)
       call ncd_io(varname='pfts1d_li'       , data=iparr, dim1name=namep, ncid=ncid, flag='write')
       iparr = get_global_index_array(patch%column(bounds%begp:bounds%endp), bounds%begp, bounds%endp, &
            subgrid_level=subgrid_level_column)
       call ncd_io(varname='pfts1d_ci'  , data=iparr              , dim1name=namep, ncid=ncid, flag='write')

       call ncd_io(varname='pfts1d_wtgcell'  , data=patch%wtgcell , dim1name=namep, ncid=ncid, flag='write')
       call ncd_io(varname='pfts1d_wtlunit'  , data=patch%wtlunit , dim1name=namep, ncid=ncid, flag='write')
       call ncd_io(varname='pfts1d_wtcol'    , data=patch%wtcol   , dim1name=namep, ncid=ncid, flag='write')
       call ncd_io(varname='pfts1d_itype_veg', data=patch%itype   , dim1name=namep, ncid=ncid, flag='write')

       do p = bounds%begp,bounds%endp
          iparr(p) = col%itype(patch%column(p))
       end do
       call ncd_io(varname='pfts1d_itype_col', data=iparr         , dim1name=namep, ncid=ncid, flag='write')

       do p = bounds%begp,bounds%endp
          iparr(p) = lun%itype(patch%landunit(p))
       enddo
       call ncd_io(varname='pfts1d_itype_lunit', data=iparr      , dim1name=namep, ncid=ncid, flag='write')

       call ncd_io(varname='pfts1d_active'   , data=patch%active  , dim1name=namep, ncid=ncid, flag='write')

       deallocate(rgarr,rlarr,rcarr,rparr)
       deallocate(igarr,ilarr,icarr,iparr)

    end if

  end subroutine hfields_1dinfo

  !-----------------------------------------------------------------------
  subroutine hist_htapes_wrapup( rstwr, nlend, bounds, &
       watsat_col, sucsat_col, bsw_col, hksat_col, cellsand_col, cellclay_col)
    !
    ! !DESCRIPTION:
    ! Write history tape(s)
    ! Determine if next time step is beginning of history interval and if so:
    !   increment the current time sample counter, open a new history file
    !   and if needed (i.e., when ntim = 1), write history data to current
    !   history file, reset field accumulation counters to zero.
    ! If primary history file is full or at the last time step of the simulation,
    !   write restart dataset and close all history fiels.
    ! If history file is full or at the last time step of the simulation:
    !   close history file
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
    use clm_time_manager, only : get_nstep, get_curr_date, get_curr_time, get_prev_date
    use clm_varcon      , only : secspday
    use perf_mod        , only : t_startf, t_stopf
    use clm_varpar      , only : nlevgrnd, nlevmaxurbgrnd, nlevsoi
    !
    ! !ARGUMENTS:
    logical, intent(in) :: rstwr    ! true => write restart file this step
    logical, intent(in) :: nlend    ! true => end of run on this step
    type(bounds_type) , intent(in) :: bounds
    real(r8)          , intent(in) :: watsat_col( bounds%begc:,1: )
    real(r8)          , intent(in) :: sucsat_col( bounds%begc:,1: )
    real(r8)          , intent(in) :: bsw_col( bounds%begc:,1: )
    real(r8)          , intent(in) :: hksat_col( bounds%begc:,1: )
    real(r8)          , intent(in) :: cellsand_col( bounds%begc:,1: )
    real(r8)          , intent(in) :: cellclay_col( bounds%begc:,1: )
    !
    ! !LOCAL VARIABLES:
    integer :: t                          ! tape index
    integer :: f                          ! field index
    integer :: ier                        ! error code
    integer :: nstep                      ! current step
    integer :: day                        ! current day (1 -> 31)
    integer :: mon                        ! current month (1 -> 12)
    integer :: yr                         ! current year (0 -> ...)
    integer :: mdcur                      ! current day
    integer :: mscur                      ! seconds of current day
    integer :: mcsec                      ! current time of day [seconds]
    integer :: daym1                      ! nstep-1 day (1 -> 31)
    integer :: monm1                      ! nstep-1 month (1 -> 12)
    integer :: yrm1                       ! nstep-1 year (0 -> ...)
    integer :: mcsecm1                    ! nstep-1 time of day [seconds]
    real(r8):: time                       ! current time
    character(len=256) :: str             ! global attribute string
    logical :: if_stop                    ! true => last time step of run
    logical, save :: do_3Dtconst = .true. ! true => write out 3D time-constant data
    character(len=*),parameter :: subname = 'hist_htapes_wrapup'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(watsat_col)    == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sucsat_col)    == (/bounds%endc, nlevgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(bsw_col)       == (/bounds%endc, nlevgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(hksat_col)     == (/bounds%endc, nlevgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(cellsand_col)  == (/bounds%endc, nlevsoi/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(cellclay_col)  == (/bounds%endc, nlevsoi/)), sourcefile, __LINE__)

    ! get current step

    nstep = get_nstep()

    ! Set calendar for current time step

    call get_curr_date (yr, mon, day, mcsec)
    call get_curr_time (mdcur, mscur)
    time = mdcur + mscur/secspday

    ! Set calendar for current for previous time step

    call get_prev_date (yrm1, monm1, daym1, mcsecm1)

    ! Loop over active history tapes, create new history files if necessary
    ! and write data to history files if end of history interval.
    do t = 1, ntapes

       if (.not. history_tape_in_use(t)) then
          cycle
       end if

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

          ! If first time sample, generate unique history file name, open file,
          ! define dims, vars, etc.


          if (tape(t)%ntimes == 1) then
             call t_startf('hist_htapes_wrapup_define')
             locfnh(t) = set_hist_filename (hist_freq=tape(t)%nhtfrq, &
                                            hist_mfilt=tape(t)%mfilt, hist_file=t)
             if (masterproc) then
                write(iulog,*) trim(subname),' : Creating history file ', trim(locfnh(t)), &
                     ' at nstep = ',get_nstep()
                write(iulog,*)'calling htape_create for file t = ',t
             endif
             call htape_create (t)

             ! Define time-constant field variables
             call htape_timeconst(t, mode='define')

             ! Define 3D time-constant field variables on first history tapes
             if ( do_3Dtconst .and. t == 1) then
                call htape_timeconst3D(t, &
                     bounds, watsat_col, sucsat_col, bsw_col, hksat_col, &
                     cellsand_col, cellclay_col, mode='define')
                TimeConst3DVars_Filename = trim(locfnh(t))
             end if

             ! Define model field variables
             call hfields_write(t, mode='define')

             ! Exit define model
             call ncd_enddef(nfid(t))
             call t_stopf('hist_htapes_wrapup_define')
          endif

          call t_startf('hist_htapes_wrapup_tconst')
          ! Write time constant history variables
          call htape_timeconst(t, mode='write')

          ! Write 3D time constant history variables to first history tapes
          if ( do_3Dtconst .and. t == 1 .and. tape(t)%ntimes == 1 )then
             call htape_timeconst3D(t, &
                  bounds, watsat_col, sucsat_col, bsw_col, hksat_col, &
                  cellsand_col, cellclay_col, mode='write')
             do_3Dtconst = .false.
          end if

          if (masterproc) then
             write(iulog,*)
             write(iulog,*) trim(subname),' : Writing current time sample to local history file ', &
                  trim(locfnh(t)),' at nstep = ',get_nstep(), &
                  ' for history time interval beginning at ', tape(t)%begtime, &
                  ' and ending at ',time
             write(iulog,*)
             call shr_sys_flush(iulog)
          endif

          ! Update beginning time of next interval
          tape(t)%begtime = time
          call t_stopf('hist_htapes_wrapup_tconst')

          ! Write history time samples
          call t_startf('hist_htapes_wrapup_write')
          call hfields_write(t, mode='write')
          call t_stopf('hist_htapes_wrapup_write')

          ! Zero necessary history buffers
          call hfields_zero(t)

       end if

    end do  ! end loop over history tapes

    ! Determine if file needs to be closed

    call hist_do_disp (ntapes, tape(:)%ntimes, tape(:)%mfilt, if_stop, if_disphist, rstwr, nlend)

    ! Close open history file
    ! Auxilary files may have been closed and saved off without being full,
    ! must reopen the files

    do t = 1, ntapes
       if (.not. history_tape_in_use(t)) then
          cycle
       end if

       if (if_disphist(t)) then
          if (tape(t)%ntimes /= 0) then
             if (masterproc) then
                write(iulog,*)
                write(iulog,*)  trim(subname),' : Closing local history file ',&
                     trim(locfnh(t)),' at nstep = ', get_nstep()
                write(iulog,*)
             endif

            call ncd_pio_closefile(nfid(t))

             if (.not.if_stop .and. (tape(t)%ntimes/=tape(t)%mfilt)) then
                call ncd_pio_openfile (nfid(t), trim(locfnh(t)), ncd_write)
             end if
          else
             if (masterproc) then
                write(iulog,*) trim(subname),' : history tape ',t,': no open file to close'
             end if
          endif
       endif
    end do

    ! Reset number of time samples to zero if file is full

    do t = 1, ntapes
       if (.not. history_tape_in_use(t)) then
          cycle
       end if

       if (if_disphist(t) .and. tape(t)%ntimes==tape(t)%mfilt) then
          tape(t)%ntimes = 0
       end if
    end do

  end subroutine hist_htapes_wrapup

  !-----------------------------------------------------------------------
  subroutine hist_restart_ncd (bounds, ncid, flag, rdate)
    !
    ! !DESCRIPTION:
    ! Read/write history file restart data.
    ! If the current history file(s) are not full, file(s) are opened
    ! so that subsequent time samples are added until the file is full.
    ! A new history file is used on a branch run.
    !
    ! !USES:
    use clm_varctl      , only : nsrest, caseid, inst_suffix, nsrStartup, nsrBranch
    use fileutils       , only : getfil
    use domainMod       , only : ldomain
    use clm_varpar      , only : nlevgrnd, nlevlak, numrad, nlevdecomp_full, mxsowings, mxharvests
    use clm_time_manager, only : is_restart
    use restUtilMod     , only : iflag_skip
    use pio
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid     ! netcdf file
    character(len=*) , intent(in)    :: flag     !'read' or 'write'
    character(len=*) , intent(in), optional :: rdate    ! restart file time stamp for name
    !
    ! !LOCAL VARIABLES:
    integer :: max_nflds                     ! Max number of fields
    integer :: num1d,beg1d,end1d             ! 1d size, beginning and ending indices
    integer :: num1d_out,beg1d_out,end1d_out ! 1d size, beginning and ending indices
    integer :: num2d                         ! 2d size (e.g. number of vertical levels)
    integer :: numa                 ! total number of atm cells across all processors
    integer :: numg                 ! total number of gridcells across all processors
    integer :: numl                 ! total number of landunits across all processors
    integer :: numc                 ! total number of columns across all processors
    integer :: nump                 ! total number of pfts across all processors
    character(len=max_namlen) :: name            ! variable name
    character(len=max_namlen) :: name_acc        ! accumulator variable name
    character(len=max_namlen) :: long_name       ! long name of variable
    character(len=max_chars)  :: long_name_acc   ! long name for accumulator
    character(len=max_chars)  :: units           ! units of variable
    character(len=max_chars)  :: units_acc       ! accumulator units
    character(len=max_chars)  :: fname           ! full name of history file
    character(len=max_chars)  :: locrest(max_tapes)  ! local history restart file names
    character(len=max_length_filename) :: my_locfnh  ! temporary version of locfnh
    character(len=max_length_filename) :: my_locfnhr ! temporary version of locfnhr

    character(len=max_namlen),allocatable :: tname(:)
    character(len=max_chars), allocatable :: tunits(:),tlongname(:)
    character(len=hist_dim_name_length), allocatable :: tmpstr(:,:)
    character(len=scale_type_strlen), allocatable :: p2c_scale_type(:)
    character(len=scale_type_strlen), allocatable :: c2l_scale_type(:)
    character(len=scale_type_strlen), allocatable :: l2g_scale_type(:)
    character(len=avgflag_strlen), allocatable :: tavgflag(:)
    integer :: start(2)

    character(len=1)   :: hnum                   ! history file index
    character(len=hist_dim_name_length)   :: type1d                 ! clm pointer 1d type
    character(len=hist_dim_name_length)   :: type1d_out             ! history buffer 1d type
    character(len=hist_dim_name_length)   :: type2d                 ! history buffer 2d type
    character(len=32)  :: dim1name               ! temporary
    character(len=32)  :: dim2name               ! temporary
    type(var_desc_t)   :: name_desc              ! variable descriptor for name
    type(var_desc_t)   :: longname_desc          ! variable descriptor for long_name
    type(var_desc_t)   :: units_desc             ! variable descriptor for units
    type(var_desc_t)   :: type1d_desc            ! variable descriptor for type1d
    type(var_desc_t)   :: type1d_out_desc        ! variable descriptor for type1d_out
    type(var_desc_t)   :: type2d_desc            ! variable descriptor for type2d
    type(var_desc_t)   :: avgflag_desc           ! variable descriptor for avgflag
    type(var_desc_t)   :: p2c_scale_type_desc    ! variable descriptor for p2c_scale_type
    type(var_desc_t)   :: c2l_scale_type_desc    ! variable descriptor for c2l_scale_type
    type(var_desc_t)   :: l2g_scale_type_desc    ! variable descriptor for l2g_scale_type
    integer :: status                            ! error status
    integer :: dimid                             ! dimension ID
    integer :: k                                 ! 1d index
    integer :: ntapes_onfile                     ! number of history tapes on the restart file
    logical, allocatable :: history_tape_in_use_onfile(:) ! whether a given history tape is in use, according to the restart file
    integer :: nflds_onfile                      ! number of history fields on the restart file
    logical :: readvar                           ! whether a variable was read successfully
    integer :: t                                 ! tape index
    integer :: f                                 ! field index
    integer :: varid                             ! variable id
    integer, allocatable :: itemp(:)             ! temporary
    real(r8), pointer :: hbuf(:,:)               ! history buffer
    real(r8), pointer :: hbuf1d(:)               ! 1d history buffer
    integer , pointer :: nacs(:,:)               ! accumulation counter
    integer , pointer :: nacs1d(:)               ! 1d accumulation counter
    integer           :: ier                     ! error code
    type(Var_desc_t)  :: vardesc                 ! netCDF variable description
    character(len=*),parameter :: subname = 'hist_restart_ncd'
!------------------------------------------------------------------------

    call get_proc_global(ng=numg, nl=numl, nc=numc, np=nump)

    ! If branch run, initialize file times and return

    if (flag == 'read') then
       if (nsrest == nsrBranch) then
          do t = 1,ntapes
             tape(t)%ntimes = 0
          end do
          return
       end if
       ! If startup run just return
       if (nsrest == nsrStartup) then
          RETURN
       end if
    endif

    ! Read history file data only for restart run (not for branch run)

    !
    ! First when writing out and in define mode, create files and define all variables
    !
    !================================================
    if (flag == 'define') then
    !================================================

       if (.not. present(rdate)) then
          call endrun(msg=' variable rdate must be present for writing restart files'//&
               errMsg(sourcefile, __LINE__))
       end if

       !
       ! On master restart file add ntapes/max_chars dimension
       ! and then add the history and history restart filenames
       !
       call ncd_defdim( ncid, 'ntapes'       , ntapes      , dimid)
       call ncd_defdim( ncid, 'max_chars'    , max_chars   , dimid)

       call ncd_defvar(ncid=ncid, varname='history_tape_in_use', xtype=ncd_log, &
            long_name="Whether this history tape is in use", &
            dim1name="ntapes")
       ier = PIO_inq_varid(ncid, 'history_tape_in_use', vardesc)
       ier = PIO_put_att(ncid, vardesc%varid, 'interpinic_flag', iflag_skip)

       call ncd_defvar(ncid=ncid, varname='locfnh', xtype=ncd_char, &
            long_name="History filename",     &
            comment="This variable NOT needed for startup or branch simulations", &
            dim1name='max_chars', dim2name="ntapes" )
       ier = PIO_inq_varid(ncid, 'locfnh', vardesc)
       ier = PIO_put_att(ncid, vardesc%varid, 'interpinic_flag', iflag_skip)

       call ncd_defvar(ncid=ncid, varname='locfnhr', xtype=ncd_char, &
            long_name="Restart history filename",     &
            comment="This variable NOT needed for startup or branch simulations", &
            dim1name='max_chars', dim2name="ntapes" )
       ier = PIO_inq_varid(ncid, 'locfnhr', vardesc)
       ier = PIO_put_att(ncid, vardesc%varid, 'interpinic_flag', iflag_skip)

       ! max_nflds is the maximum number of fields on any tape
       ! max_flds is the maximum number possible number of fields

       max_nflds = max_nFields()

       ! Loop over tapes - write out namelist information to each restart-history tape
       ! only read/write accumulators and counters if needed

       do t = 1,ntapes
          if (.not. history_tape_in_use(t)) then
             cycle
          end if

          ! Create the restart history filename and open it
          write(hnum,'(i1.1)') t-1
          locfnhr(t) = "./" // trim(caseid) //"."// trim(compname) // trim(inst_suffix) &
                        // ".rh" // hnum //"."// trim(rdate) //".nc"

          call htape_create( t, histrest=.true. )

          ! Add read/write accumultators and counters if needed
          if (.not. tape(t)%is_endhist) then
             do f = 1,tape(t)%nflds
                name           =  tape(t)%hlist(f)%field%name
                long_name      =  tape(t)%hlist(f)%field%long_name
                units          =  tape(t)%hlist(f)%field%units
                name_acc       =  trim(name) // "_acc"
                units_acc      =  "unitless positive integer"
                long_name_acc  =  trim(long_name) // " accumulator number of samples"
                type1d_out     =  tape(t)%hlist(f)%field%type1d_out
                type2d         =  tape(t)%hlist(f)%field%type2d
                num2d          =  tape(t)%hlist(f)%field%num2d
                nacs           => tape(t)%hlist(f)%nacs
                hbuf           => tape(t)%hlist(f)%hbuf

                if (type1d_out == grlnd) then
                   if (ldomain%isgrid2d) then
                      dim1name = 'lon'      ; dim2name = 'lat'
                   else
                      dim1name = trim(grlnd); dim2name = 'undefined'
                   end if
                else
                   dim1name = type1d_out ; dim2name = 'undefined'
                endif

                if (dim2name == 'undefined') then
                   if (num2d == 1) then
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name), xtype=ncd_double, &
                           dim1name=dim1name, &
                           long_name=trim(long_name), units=trim(units))
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name_acc), xtype=ncd_int,  &
                           dim1name=dim1name, &
                           long_name=trim(long_name_acc), units=trim(units_acc))
                   else
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name), xtype=ncd_double, &
                           dim1name=dim1name, dim2name=type2d, &
                           long_name=trim(long_name), units=trim(units))
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name_acc), xtype=ncd_int,  &
                           dim1name=dim1name, dim2name=type2d, &
                           long_name=trim(long_name_acc), units=trim(units_acc))
                   end if
                else
                   if (num2d == 1) then
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name), xtype=ncd_double, &
                           dim1name=dim1name, dim2name=dim2name, &
                           long_name=trim(long_name), units=trim(units))
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name_acc), xtype=ncd_int,  &
                           dim1name=dim1name, dim2name=dim2name, &
                           long_name=trim(long_name_acc), units=trim(units_acc))
                   else
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name), xtype=ncd_double, &
                           dim1name=dim1name, dim2name=dim2name, dim3name=type2d, &
                           long_name=trim(long_name), units=trim(units))
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name_acc), xtype=ncd_int,  &
                           dim1name=dim1name, dim2name=dim2name, dim3name=type2d, &
                           long_name=trim(long_name_acc), units=trim(units_acc))
                   end if
                endif
             end do
          endif

          !
          ! Add namelist information to each restart history tape
          !
          call ncd_defdim( ncid_hist(t), 'fname_lenp2'  , max_namlen+2, dimid)
          call ncd_defdim( ncid_hist(t), 'fname_len'    , max_namlen  , dimid)
          call ncd_defdim( ncid_hist(t), 'avgflag_len'  , avgflag_strlen, dimid)
          call ncd_defdim( ncid_hist(t), 'scalar'       , 1           , dimid)
          call ncd_defdim( ncid_hist(t), 'max_chars'    , max_chars   , dimid)
          call ncd_defdim( ncid_hist(t), 'max_nflds'    , max_nflds   ,  dimid)
          call ncd_defdim( ncid_hist(t), 'max_flds'     , max_flds    , dimid)

          call ncd_defvar(ncid=ncid_hist(t), varname='nhtfrq', xtype=ncd_int, &
               long_name="Frequency of history writes",               &
               comment="Namelist item", &
               units="absolute value of negative is in hours, 0=monthly, positive is time-steps",     &
               dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='mfilt', xtype=ncd_int, &
               long_name="Number of history time samples on a file", units="unitless",     &
               comment="Namelist item", &
               dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='ncprec', xtype=ncd_int, &
               long_name="Flag for data precision", flag_values=(/1,2/), &
               comment="Namelist item", &
               nvalid_range=(/1,2/), &
               flag_meanings=(/"single-precision", "double-precision"/), &
               dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='dov2xy', xtype=ncd_log, &
               long_name="Output on 2D grid format (TRUE) or vector format (FALSE)", &
               comment="Namelist item", &
               dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='fincl', xtype=ncd_char, &
               comment="Namelist item", &
               long_name="Fieldnames to include", &
               dim1name='fname_lenp2', dim2name='max_flds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='fexcl', xtype=ncd_char, &
               comment="Namelist item", &
               long_name="Fieldnames to exclude",  &
               dim1name='fname_lenp2', dim2name='max_flds' )

          call ncd_defvar(ncid=ncid_hist(t), varname='nflds', xtype=ncd_int, &
               long_name="Number of fields on file", units="unitless",        &
               dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='ntimes', xtype=ncd_int, &
               long_name="Number of time steps on file", units="time-step",     &
               dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='is_endhist', xtype=ncd_log, &
               long_name="End of history file", dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='begtime', xtype=ncd_double, &
               long_name="Beginning time", units="time units",     &
               dim1name='scalar')

          call ncd_defvar(ncid=ncid_hist(t), varname='num2d', xtype=ncd_int, &
               long_name="Size of second dimension", units="unitless",     &
               dim1name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='hpindex', xtype=ncd_int, &
               long_name="History pointer index", units="unitless",     &
               dim1name='max_nflds' )

          call ncd_defvar(ncid=ncid_hist(t), varname='avgflag', xtype=ncd_char, &
               long_name="Averaging flag", &
               units="A=Average, X=Maximum, M=Minimum, I=Instantaneous, SUM=Sum", &
               dim1name='avgflag_len', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='name', xtype=ncd_char, &
               long_name="Fieldnames",  &
               dim1name='fname_len', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='long_name', xtype=ncd_char, &
               long_name="Long descriptive names for fields", &
               dim1name='max_chars', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='units', xtype=ncd_char, &
               long_name="Units for each history field output", &
               dim1name='max_chars', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='type1d', xtype=ncd_char, &
               long_name="1st dimension type", &
               dim1name='string_length', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='type1d_out', xtype=ncd_char, &
               long_name="1st output dimension type", &
               dim1name='string_length', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='type2d', xtype=ncd_char, &
               long_name="2nd dimension type", &
               dim1name='string_length', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='p2c_scale_type', xtype=ncd_char, &
               long_name="PFT to column scale type", &
               dim1name='scale_type_string_length', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='c2l_scale_type', xtype=ncd_char, &
               long_name="column to landunit scale type", &
               dim1name='scale_type_string_length', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='l2g_scale_type', xtype=ncd_char, &
               long_name="landunit to gridpoint scale type", &
               dim1name='scale_type_string_length', dim2name='max_nflds' )

          call ncd_enddef(ncid_hist(t))

       end do   ! end of ntapes loop

       RETURN

    !
    ! First write out namelist information to each restart history file
    !
    !================================================
    else if (flag == 'write') then
    !================================================

       ! Add history filenames to master restart file
       do t = 1,ntapes
          call ncd_io('history_tape_in_use', history_tape_in_use(t), 'write', ncid, nt=t)
          if (history_tape_in_use(t)) then
             my_locfnh  = locfnh(t)
             my_locfnhr = locfnhr(t)
          else
             my_locfnh  = 'non_existent_file'
             my_locfnhr = 'non_existent_file'
          end if
          call ncd_io('locfnh',  my_locfnh,  'write', ncid, nt=t)
          call ncd_io('locfnhr', my_locfnhr, 'write', ncid, nt=t)
       end do

       fincl(:,1)  = hist_fincl1(:)
       fincl(:,2)  = hist_fincl2(:)
       fincl(:,3)  = hist_fincl3(:)
       fincl(:,4)  = hist_fincl4(:)
       fincl(:,5)  = hist_fincl5(:)
       fincl(:,6)  = hist_fincl6(:)
       fincl(:,7)  = hist_fincl7(:)
       fincl(:,8)  = hist_fincl8(:)
       fincl(:,9)  = hist_fincl9(:)
       fincl(:,10) = hist_fincl10(:)

       fexcl(:,1)  = hist_fexcl1(:)
       fexcl(:,2)  = hist_fexcl2(:)
       fexcl(:,3)  = hist_fexcl3(:)
       fexcl(:,4)  = hist_fexcl4(:)
       fexcl(:,5)  = hist_fexcl5(:)
       fexcl(:,6)  = hist_fexcl6(:)
       fexcl(:,7)  = hist_fexcl7(:)
       fexcl(:,8)  = hist_fexcl8(:)
       fexcl(:,9)  = hist_fexcl9(:)
       fexcl(:,10) = hist_fexcl10(:)

       max_nflds = max_nFields()

       start(1)=1

       !
       ! Add history namelist data to each history restart tape
       !
       allocate(itemp(max_nflds))

       do t = 1,ntapes
          if (.not. history_tape_in_use(t)) then
             cycle
          end if

          call ncd_io(varname='fincl', data=fincl(:,t), ncid=ncid_hist(t), flag='write')

          call ncd_io(varname='fexcl', data=fexcl(:,t), ncid=ncid_hist(t), flag='write')

          call ncd_io(varname='is_endhist', data=tape(t)%is_endhist, ncid=ncid_hist(t), flag='write')

          call ncd_io(varname='dov2xy', data=tape(t)%dov2xy, ncid=ncid_hist(t), flag='write')

          itemp(:) = 0
          do f=1,tape(t)%nflds
             itemp(f) = tape(t)%hlist(f)%field%num2d
          end do
          call ncd_io(varname='num2d', data=itemp(:), ncid=ncid_hist(t), flag='write')

          itemp(:) = 0
          do f=1,tape(t)%nflds
             itemp(f) = tape(t)%hlist(f)%field%hpindex
          end do
          call ncd_io(varname='hpindex', data=itemp(:), ncid=ncid_hist(t), flag='write')

          call ncd_io('nflds',        tape(t)%nflds,   'write', ncid_hist(t) )
          call ncd_io('ntimes',       tape(t)%ntimes,  'write', ncid_hist(t) )
          call ncd_io('nhtfrq',  tape(t)%nhtfrq,  'write', ncid_hist(t) )
          call ncd_io('mfilt',   tape(t)%mfilt,   'write', ncid_hist(t) )
          call ncd_io('ncprec',  tape(t)%ncprec,  'write', ncid_hist(t) )
          call ncd_io('begtime',      tape(t)%begtime, 'write', ncid_hist(t) )
          allocate(tmpstr(tape(t)%nflds,3 ),tname(tape(t)%nflds), &
               tavgflag(tape(t)%nflds),tunits(tape(t)%nflds),tlongname(tape(t)%nflds), &
               p2c_scale_type(tape(t)%nflds), c2l_scale_type(tape(t)%nflds), &
               l2g_scale_type(tape(t)%nflds))
          do f=1,tape(t)%nflds
             tname(f)  = tape(t)%hlist(f)%field%name
             tunits(f) = tape(t)%hlist(f)%field%units
             tlongname(f) = tape(t)%hlist(f)%field%long_name
             tmpstr(f,1) = tape(t)%hlist(f)%field%type1d
             tmpstr(f,2) = tape(t)%hlist(f)%field%type1d_out
             tmpstr(f,3) = tape(t)%hlist(f)%field%type2d
             tavgflag(f) = tape(t)%hlist(f)%avgflag
             p2c_scale_type(f) = tape(t)%hlist(f)%field%p2c_scale_type
             c2l_scale_type(f) = tape(t)%hlist(f)%field%c2l_scale_type
             l2g_scale_type(f) = tape(t)%hlist(f)%field%l2g_scale_type
          end do
          call ncd_io( 'name', tname, 'write',ncid_hist(t))
          call ncd_io('long_name', tlongname, 'write', ncid_hist(t))
          call ncd_io('units', tunits, 'write',ncid_hist(t))
          call ncd_io('type1d', tmpstr(:,1), 'write', ncid_hist(t))
          call ncd_io('type1d_out', tmpstr(:,2), 'write', ncid_hist(t))
          call ncd_io('type2d', tmpstr(:,3), 'write', ncid_hist(t))
          call ncd_io('avgflag',tavgflag , 'write', ncid_hist(t))
          call ncd_io('p2c_scale_type', p2c_scale_type, 'write', ncid_hist(t))
          call ncd_io('c2l_scale_type', c2l_scale_type, 'write', ncid_hist(t))
          call ncd_io('l2g_scale_type', l2g_scale_type, 'write', ncid_hist(t))
          deallocate(tname,tlongname,tunits,tmpstr,tavgflag)
          deallocate(p2c_scale_type, c2l_scale_type, l2g_scale_type)
       enddo
       deallocate(itemp)

    !
    ! Read in namelist information
    !
    !================================================
    else if (flag == 'read') then
    !================================================

       call ncd_inqdlen(ncid,dimid,ntapes_onfile, name='ntapes')
       if (is_restart()) then
          if (ntapes_onfile /= ntapes) then
             write(iulog,*) 'ntapes = ', ntapes, ' ntapes_onfile = ', ntapes_onfile
             call endrun(msg=' ERROR: number of ntapes differs from restart file. '// &
                  'You can NOT change history options on restart.', &
                  additional_msg=errMsg(sourcefile, __LINE__))
          end if

          if (ntapes > 0) then
             allocate(history_tape_in_use_onfile(ntapes))
             call ncd_io('history_tape_in_use', history_tape_in_use_onfile, 'read', ncid, &
                  readvar=readvar)
             if (.not. readvar) then
                ! BACKWARDS_COMPATIBILITY(wjs, 2018-10-06) Old restart files do not have
                ! 'history_tape_in_use'. However, before now, this has implicitly been
                ! true for all tapes <= ntapes.
                history_tape_in_use_onfile(:) = .true.
             end if
             do t = 1, ntapes
                if (history_tape_in_use_onfile(t) .neqv. history_tape_in_use(t)) then
                   write(iulog,*) subname//' ERROR: history_tape_in_use on restart file'
                   write(iulog,*) 'disagrees with current run: For tape ', t
                   write(iulog,*) 'On restart file: ', history_tape_in_use_onfile(t)
                   write(iulog,*) 'In current run : ', history_tape_in_use(t)
                   write(iulog,*) 'This suggests that this tape was empty in one case,'
                   write(iulog,*) 'but non-empty in the other. (history_tape_in_use .false.'
                   write(iulog,*) 'means that history tape is empty.)'
                   call endrun(msg=' ERROR: history_tape_in_use differs from restart file. '// &
                        'You can NOT change history options on restart.', &
                        additional_msg=errMsg(sourcefile, __LINE__))
                end if
             end do

             call ncd_io('locfnh',  locfnh(1:ntapes),  'read', ncid )
             call ncd_io('locfnhr', locrest(1:ntapes), 'read', ncid )
             do t = 1,ntapes
                call strip_null(locrest(t))
                call strip_null(locfnh(t))
             end do
          end if
       end if

       ! Determine necessary indices - the following is needed if model decomposition is different on restart

       start(1)=1

       if ( is_restart() )then
          do t = 1,ntapes
             if (.not. history_tape_in_use(t)) then
                cycle
             end if

             call getfil( locrest(t), locfnhr(t), 0 )
             call ncd_pio_openfile (ncid_hist(t), trim(locfnhr(t)), ncd_nowrite)

             if ( t == 1 )then

                call ncd_inqdlen(ncid_hist(1),dimid,max_nflds,name='max_nflds')

                allocate(itemp(max_nflds))
             end if

             call ncd_inqvid(ncid_hist(t), 'name',           varid, name_desc)
             call ncd_inqvid(ncid_hist(t), 'long_name',      varid, longname_desc)
             call ncd_inqvid(ncid_hist(t), 'units',          varid, units_desc)
             call ncd_inqvid(ncid_hist(t), 'type1d',         varid, type1d_desc)
             call ncd_inqvid(ncid_hist(t), 'type1d_out',     varid, type1d_out_desc)
             call ncd_inqvid(ncid_hist(t), 'type2d',         varid, type2d_desc)
             call ncd_inqvid(ncid_hist(t), 'avgflag',        varid, avgflag_desc)
             call ncd_inqvid(ncid_hist(t), 'p2c_scale_type', varid, p2c_scale_type_desc)
             call ncd_inqvid(ncid_hist(t), 'c2l_scale_type', varid, c2l_scale_type_desc)
             call ncd_inqvid(ncid_hist(t), 'l2g_scale_type', varid, l2g_scale_type_desc)

             call ncd_io(varname='fincl', data=fincl(:,t), ncid=ncid_hist(t), flag='read')

             call ncd_io(varname='fexcl', data=fexcl(:,t), ncid=ncid_hist(t), flag='read')

             call ncd_io('nflds',   nflds_onfile, 'read', ncid_hist(t) )
             if ( nflds_onfile /= tape(t)%nflds )then
                write(iulog,*) 'nflds = ', tape(t)%nflds, ' nflds_onfile = ', nflds_onfile
                call endrun(msg=' ERROR: number of fields different than on restart file!,'// &
                     ' you can NOT change history options on restart!' //&
                     errMsg(sourcefile, __LINE__))
             end if
             call ncd_io('ntimes',  tape(t)%ntimes, 'read', ncid_hist(t) )
             call ncd_io('nhtfrq',  tape(t)%nhtfrq, 'read', ncid_hist(t) )
             call ncd_io('mfilt',   tape(t)%mfilt, 'read', ncid_hist(t) )
             call ncd_io('ncprec',  tape(t)%ncprec, 'read', ncid_hist(t) )
             call ncd_io('begtime', tape(t)%begtime, 'read', ncid_hist(t) )

             call ncd_io(varname='is_endhist', data=tape(t)%is_endhist, ncid=ncid_hist(t), flag='read')
             call ncd_io(varname='dov2xy', data=tape(t)%dov2xy, ncid=ncid_hist(t), flag='read')
             call ncd_io(varname='num2d', data=itemp(:), ncid=ncid_hist(t), flag='read')
             do f=1,tape(t)%nflds
                tape(t)%hlist(f)%field%num2d = itemp(f)
             end do

             call ncd_io(varname='hpindex', data=itemp(:), ncid=ncid_hist(t), flag='read')
             do f=1,tape(t)%nflds
                tape(t)%hlist(f)%field%hpindex = itemp(f)
             end do

             do f=1,tape(t)%nflds
                start(2) = f
                call ncd_io( name_desc,           tape(t)%hlist(f)%field%name,       &
                             'read', ncid_hist(t), start )
                call ncd_io( longname_desc,       tape(t)%hlist(f)%field%long_name,  &
                             'read', ncid_hist(t), start )
                call ncd_io( units_desc,          tape(t)%hlist(f)%field%units,      &
                             'read', ncid_hist(t), start )
                call ncd_io( type1d_desc,         tape(t)%hlist(f)%field%type1d,     &
                             'read', ncid_hist(t), start )
                call ncd_io( type1d_out_desc,     tape(t)%hlist(f)%field%type1d_out, &
                             'read', ncid_hist(t), start )
                call ncd_io( type2d_desc,         tape(t)%hlist(f)%field%type2d,     &
                             'read', ncid_hist(t), start )
                call ncd_io( avgflag_desc,        tape(t)%hlist(f)%avgflag,          &
                             'read', ncid_hist(t), start )
                call ncd_io( p2c_scale_type_desc, tape(t)%hlist(f)%field%p2c_scale_type,   &
                             'read', ncid_hist(t), start )
                call ncd_io( c2l_scale_type_desc, tape(t)%hlist(f)%field%c2l_scale_type,   &
                             'read', ncid_hist(t), start )
                call ncd_io( l2g_scale_type_desc, tape(t)%hlist(f)%field%l2g_scale_type,   &
                             'read', ncid_hist(t), start )
                call strip_null(tape(t)%hlist(f)%field%name)
                call strip_null(tape(t)%hlist(f)%field%long_name)
                call strip_null(tape(t)%hlist(f)%field%units)
                call strip_null(tape(t)%hlist(f)%field%type1d)
                call strip_null(tape(t)%hlist(f)%field%type1d_out)
                call strip_null(tape(t)%hlist(f)%field%type2d)
                call strip_null(tape(t)%hlist(f)%field%p2c_scale_type)
                call strip_null(tape(t)%hlist(f)%field%c2l_scale_type)
                call strip_null(tape(t)%hlist(f)%field%l2g_scale_type)
                call strip_null(tape(t)%hlist(f)%avgflag)

                type1d_out = trim(tape(t)%hlist(f)%field%type1d_out)
                select case (trim(type1d_out))
                case (grlnd)
                   num1d_out = numg
                   beg1d_out = bounds%begg
                   end1d_out = bounds%endg
                case (nameg)
                   num1d_out = numg
                   beg1d_out = bounds%begg
                   end1d_out = bounds%endg
                case (namel)
                   num1d_out = numl
                   beg1d_out = bounds%begl
                   end1d_out = bounds%endl
                case (namec)
                   num1d_out = numc
                   beg1d_out = bounds%begc
                   end1d_out = bounds%endc
                case (namep)
                   num1d_out = nump
                   beg1d_out = bounds%begp
                   end1d_out = bounds%endp
                case default
                   write(iulog,*) trim(subname),' ERROR: read unknown 1d output type=',trim(type1d_out)
                   call endrun(msg=errMsg(sourcefile, __LINE__))
                end select

                tape(t)%hlist(f)%field%num1d_out = num1d_out
                tape(t)%hlist(f)%field%beg1d_out = beg1d_out
                tape(t)%hlist(f)%field%end1d_out = end1d_out

                num2d  = tape(t)%hlist(f)%field%num2d
                allocate (tape(t)%hlist(f)%hbuf(beg1d_out:end1d_out,num2d), &
                          tape(t)%hlist(f)%nacs(beg1d_out:end1d_out,num2d), &
                          stat=status)
                if (status /= 0) then
                   write(iulog,*) trim(subname),' ERROR: allocation error for hbuf,nacs at t,f=',t,f
                   call endrun(msg=errMsg(sourcefile, __LINE__))
                endif
                tape(t)%hlist(f)%hbuf(:,:) = 0._r8
                tape(t)%hlist(f)%nacs(:,:) = 0

                type1d = tape(t)%hlist(f)%field%type1d
                select case (type1d)
                case (grlnd)
                   num1d = numg
                   beg1d = bounds%begg
                   end1d = bounds%endg
                case (nameg)
                   num1d = numg
                   beg1d = bounds%begg
                   end1d = bounds%endg
                case (namel)
                   num1d = numl
                   beg1d = bounds%begl
                   end1d = bounds%endl
                case (namec)
                   num1d = numc
                   beg1d = bounds%begc
                   end1d = bounds%endc
                case (namep)
                   num1d = nump
                   beg1d = bounds%begp
                   end1d = bounds%endp
                case default
                   write(iulog,*) trim(subname),' ERROR: read unknown 1d type=',type1d
                   call endrun(msg=errMsg(sourcefile, __LINE__))
                end select

                tape(t)%hlist(f)%field%num1d = num1d
                tape(t)%hlist(f)%field%beg1d = beg1d
                tape(t)%hlist(f)%field%end1d = end1d

             end do   ! end of flds loop

             ! If history file is not full, open it

             if (tape(t)%ntimes /= 0) then
                call ncd_pio_openfile (nfid(t), trim(locfnh(t)), ncd_write)
             end if

          end do  ! end of tapes loop

          hist_fincl1(:)  = fincl(:,1)
          hist_fincl2(:)  = fincl(:,2)
          hist_fincl3(:)  = fincl(:,3)
          hist_fincl4(:)  = fincl(:,4)
          hist_fincl5(:)  = fincl(:,5)
          hist_fincl6(:)  = fincl(:,6)
          hist_fincl7(:)  = fincl(:,7)
          hist_fincl8(:)  = fincl(:,8)
          hist_fincl9(:)  = fincl(:,9)
          hist_fincl10(:) = fincl(:,10)

          hist_fexcl1(:)  = fexcl(:,1)
          hist_fexcl2(:)  = fexcl(:,2)
          hist_fexcl3(:)  = fexcl(:,3)
          hist_fexcl4(:)  = fexcl(:,4)
          hist_fexcl5(:)  = fexcl(:,5)
          hist_fexcl6(:)  = fexcl(:,6)
          hist_fexcl7(:)  = fexcl(:,7)
          hist_fexcl8(:)  = fexcl(:,8)
          hist_fexcl9(:)  = fexcl(:,9)
          hist_fexcl10(:) = fexcl(:,10)

       end if

       if ( allocated(itemp) ) deallocate(itemp)

    end if

    !======================================================================
    ! Read/write history file restart data.
    ! If the current history file(s) are not full, file(s) are opened
    ! so that subsequent time samples are added until the file is full.
    ! A new history file is used on a branch run.
    !======================================================================

    if (flag == 'write') then

       do t = 1,ntapes
          if (.not. history_tape_in_use(t)) then
             cycle
          end if

          if (.not. tape(t)%is_endhist) then

             do f = 1,tape(t)%nflds
                name       =  tape(t)%hlist(f)%field%name
                name_acc   =  trim(name) // "_acc"
                type1d_out =  tape(t)%hlist(f)%field%type1d_out
                type2d     =  tape(t)%hlist(f)%field%type2d
                num2d      =  tape(t)%hlist(f)%field%num2d
                beg1d_out  =  tape(t)%hlist(f)%field%beg1d_out
                end1d_out  =  tape(t)%hlist(f)%field%end1d_out
                nacs       => tape(t)%hlist(f)%nacs
                hbuf       => tape(t)%hlist(f)%hbuf

                if (num2d == 1) then
                   allocate(hbuf1d(beg1d_out:end1d_out), &
                            nacs1d(beg1d_out:end1d_out), stat=status)
                   if (status /= 0) then
                      write(iulog,*) trim(subname),' ERROR: allocation'
                      call endrun(msg=errMsg(sourcefile, __LINE__))
                   end if

                   hbuf1d(beg1d_out:end1d_out) = hbuf(beg1d_out:end1d_out,1)
                   nacs1d(beg1d_out:end1d_out) = nacs(beg1d_out:end1d_out,1)

                   call ncd_io(ncid=ncid_hist(t), flag='write', varname=trim(name), &
                        dim1name=type1d_out, data=hbuf1d)
                   call ncd_io(ncid=ncid_hist(t), flag='write', varname=trim(name_acc), &
                        dim1name=type1d_out, data=nacs1d)

                   deallocate(hbuf1d)
                   deallocate(nacs1d)
                else
                   call ncd_io(ncid=ncid_hist(t), flag='write', varname=trim(name), &
                        dim1name=type1d_out, data=hbuf)
                   call ncd_io(ncid=ncid_hist(t), flag='write', varname=trim(name_acc), &
                        dim1name=type1d_out, data=nacs)
                end if

             end do

          end if  ! end of is_endhist block

          call ncd_pio_closefile(ncid_hist(t))

       end do   ! end of ntapes loop

    else if (flag == 'read') then

       ! Read history restart information if history files are not full

       do t = 1,ntapes
          if (.not. history_tape_in_use(t)) then
             cycle
          end if

          if (.not. tape(t)%is_endhist) then

             do f = 1,tape(t)%nflds
                name       =  tape(t)%hlist(f)%field%name
                name_acc   =  trim(name) // "_acc"
                type1d_out =  tape(t)%hlist(f)%field%type1d_out
                type2d     =  tape(t)%hlist(f)%field%type2d
                num2d      =  tape(t)%hlist(f)%field%num2d
                beg1d_out  =  tape(t)%hlist(f)%field%beg1d_out
                end1d_out  =  tape(t)%hlist(f)%field%end1d_out
                nacs       => tape(t)%hlist(f)%nacs
                hbuf       => tape(t)%hlist(f)%hbuf

                if (num2d == 1) then
                   allocate(hbuf1d(beg1d_out:end1d_out), &
                        nacs1d(beg1d_out:end1d_out), stat=status)
                   if (status /= 0) then
                      write(iulog,*) trim(subname),' ERROR: allocation'
                      call endrun(msg=errMsg(sourcefile, __LINE__))
                   end if

                   call ncd_io(ncid=ncid_hist(t), flag='read', varname=trim(name), &
                        dim1name=type1d_out, data=hbuf1d)
                   call ncd_io(ncid=ncid_hist(t), flag='read', varname=trim(name_acc), &
                        dim1name=type1d_out, data=nacs1d)

                   hbuf(beg1d_out:end1d_out,1) = hbuf1d(beg1d_out:end1d_out)
                   nacs(beg1d_out:end1d_out,1) = nacs1d(beg1d_out:end1d_out)

                   deallocate(hbuf1d)
                   deallocate(nacs1d)
                else
                   call ncd_io(ncid=ncid_hist(t), flag='read', varname=trim(name), &
                        dim1name=type1d_out, data=hbuf)
                   call ncd_io(ncid=ncid_hist(t), flag='read', varname=trim(name_acc), &
                        dim1name=type1d_out, data=nacs)
                end if
             end do

          end if

          call ncd_pio_closefile(ncid_hist(t))

       end do

    end if

  end subroutine hist_restart_ncd

  !-----------------------------------------------------------------------
  integer function max_nFields()
    !
    ! !DESCRIPTION:
    ! Get the maximum number of fields on all tapes.
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: t  ! index
    character(len=*),parameter :: subname = 'max_nFields'
    !-----------------------------------------------------------------------

    max_nFields = 0
    do t = 1,ntapes
       max_nFields = max(max_nFields, tape(t)%nflds)
    end do
    return
  end function max_nFields

  !-----------------------------------------------------------------------
  character(len=max_namlen) function getname (inname)
    !
    ! !DESCRIPTION:
    ! Retrieve name portion of inname. If an averaging flag separater character
    ! is present (:) in inname, lop it off.
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: inname
    !
    ! !LOCAL VARIABLES:
    integer :: length
    integer :: i
    character(len=*),parameter :: subname = 'getname'
    !-----------------------------------------------------------------------

     length = len (inname)

     if (length < max_namlen .or. length > max_namlen+2) then
        write(iulog,*) trim(subname),' ERROR: bad length=',length
        call endrun(msg=errMsg(sourcefile, __LINE__))
     end if

     getname = ' '
     do i = 1,max_namlen
        if (inname(i:i) == ':') exit
        getname(i:i) = inname(i:i)
     end do

   end function getname

   !-----------------------------------------------------------------------
   character(len=avgflag_strlen) function getflag (inname)
     !
     ! !DESCRIPTION:
     ! Retrieve flag portion of inname. If an averaging flag separater character
     ! is present (:) in inname, return the character after it as the flag
     !
     ! !ARGUMENTS:
     character(len=*) inname   ! character string
     !
     ! !LOCAL VARIABLES:
     integer :: length         ! length of inname
     integer :: i              ! loop index
     character(len=*),parameter :: subname = 'getflag'
     !-----------------------------------------------------------------------

     length = len (inname)

     if (length < max_namlen .or. length > max_namlen+2) then
        write(iulog,*) trim(subname),' ERROR: bad length=',length
        call endrun(msg=errMsg(sourcefile, __LINE__))
     end if

     getflag = ' '
     do i = 1,length
        if (inname(i:i) == ':') then
           getflag = trim(inname(i+1:length))
           exit
        end if
     end do

   end function getflag

   !-----------------------------------------------------------------------
   subroutine list_index (list, name, index)
     !
     ! !ARGUMENTS:
     character(len=*), intent(in) :: list(max_flds)  ! input list of names, possibly ":" delimited
     character(len=max_namlen), intent(in) :: name   ! name to be searched for
     integer, intent(out) :: index                   ! index of "name" in "list"
     !
     ! !LOCAL VARIABLES:
     !EOP
     character(len=max_namlen) :: listname           ! input name with ":" stripped off.
     integer f                                       ! field index
     character(len=*),parameter :: subname = 'list_index'
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
   character(len=max_length_filename) function set_hist_filename (hist_freq, hist_mfilt, hist_file)
     !
     ! !DESCRIPTION:
     ! Determine history dataset filenames.
     ! Note that the first history tape is index 1 in the code but contains 'h0' in its output
     ! filenames.
     !
     ! !USES:
     use clm_varctl, only : caseid, inst_suffix
     use clm_time_manager, only : get_curr_date, get_prev_date
     !
     ! !ARGUMENTS:
     integer, intent(in)  :: hist_freq   !history file frequency
     integer, intent(in)  :: hist_mfilt  !history file number of time-samples
     integer, intent(in)  :: hist_file   !history file index
     !
     ! !LOCAL VARIABLES:
     !EOP
     character(len=max_chars) :: cdate !date char string
     character(len=  1) :: hist_index  !p,1 or 2 (currently)
     integer :: day                    !day (1 -> 31)
     integer :: mon                    !month (1 -> 12)
     integer :: yr                     !year (0 -> ...)
     integer :: sec                    !seconds into current day
     integer :: filename_length
     character(len=*),parameter :: subname = 'set_hist_filename'
     !-----------------------------------------------------------------------

   if (hist_freq == 0 .and. hist_mfilt == 1) then   !monthly
      call get_prev_date (yr, mon, day, sec)
      write(cdate,'(i4.4,"-",i2.2)') yr,mon
   else                        !other
      call get_curr_date (yr, mon, day, sec)
      write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
   endif
   write(hist_index,'(i1.1)') hist_file - 1
   set_hist_filename = "./"//trim(caseid)//"."//trim(compname)//trim(inst_suffix)//&
                       ".h"//hist_index//"."//trim(cdate)//".nc"

   ! check to see if the concatenated filename exceeded the
   ! length. Simplest way to do this is ensure that the file
   ! extension is '.nc'.
   filename_length = len_trim(set_hist_filename)
   if (set_hist_filename(filename_length-2:filename_length) /= '.nc') then
      write(iulog, '(a,a,a,a,a)') 'ERROR: ', subname, &
           ' : expected file extension ".nc", received extension "', &
           set_hist_filename(filename_length-2:filename_length), '"'
      write(iulog, '(a,a,a,a,a)') 'ERROR: ', subname, &
           ' : filename : "', set_hist_filename, '"'
      write(iulog, '(a,a,a,i3,a,i3)') 'ERROR: ', subname, &
           ' Did the constructed filename exceed the maximum length? : filename length = ', &
           filename_length, ', max length = ', max_length_filename
      call endrun(msg=errMsg(sourcefile, __LINE__))
   end if
  end function set_hist_filename

  !-----------------------------------------------------------------------
  subroutine hist_addfld1d (fname, units, avgflag, long_name, type1d_out, &
                        ptr_gcell, ptr_lunit, ptr_col, ptr_patch, ptr_lnd, &
                        ptr_atm, p2c_scale_type, c2l_scale_type, &
                        l2g_scale_type, set_lake, set_nolake, set_urb, set_nourb, &
                        set_noglc, set_spec, default)
    !
    ! !DESCRIPTION:
    ! Initialize a single level history field. The pointer inputs, ptr\_*,
    ! point to the appropriate-type array storing the raw history data points.
    ! The value of type1d passed to allhistfldlist\_add\_fld determines which of the
    ! 1d type of the output and the beginning and ending indices the history
    ! buffer field). All fields default to being written to the first history tape
    ! unless 'default' is set to 'inactive'.
    !
    ! !ARGUMENTS:
    character(len=*), intent(in)           :: fname          ! field name
    character(len=*), intent(in)           :: units          ! units of field
    character(len=*), intent(in)           :: avgflag        ! time averaging flag
    character(len=*), intent(in)           :: long_name      ! long name of field
    character(len=*), optional, intent(in) :: type1d_out     ! output type (from data type)
    real(r8)        , optional, pointer    :: ptr_gcell(:)   ! pointer to gridcell array
    real(r8)        , optional, pointer    :: ptr_lunit(:)   ! pointer to landunit array
    real(r8)        , optional, pointer    :: ptr_col(:)     ! pointer to column array
    real(r8)        , optional, pointer    :: ptr_patch(:)   ! pointer to patch array
    real(r8)        , optional, pointer    :: ptr_lnd(:)     ! pointer to lnd array
    real(r8)        , optional, pointer    :: ptr_atm(:)     ! pointer to atm array
    real(r8)        , optional, intent(in) :: set_lake       ! value to set lakes to
    real(r8)        , optional, intent(in) :: set_nolake     ! value to set non-lakes to
    real(r8)        , optional, intent(in) :: set_urb        ! value to set urban to
    real(r8)        , optional, intent(in) :: set_nourb      ! value to set non-urban to
    real(r8)        , optional, intent(in) :: set_noglc      ! value to set non-glacier to
    real(r8)        , optional, intent(in) :: set_spec       ! value to set special to
    character(len=*), optional, intent(in) :: p2c_scale_type ! scale type for subgrid averaging of pfts to column
    character(len=*), optional, intent(in) :: c2l_scale_type ! scale type for subgrid averaging of columns to landunits
    character(len=*), optional, intent(in) :: l2g_scale_type ! scale type for subgrid averaging of landunits to gridcells
    character(len=*), optional, intent(in) :: default        ! if set to 'inactive, field will not appear on primary tape
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l,g                 ! indices
    integer :: hpindex                 ! history buffer pointer index
    character(len=hist_dim_name_length) :: l_type1d       ! 1d data type
    character(len=hist_dim_name_length) :: l_type1d_out   ! 1d output type
    character(len=scale_type_strlen) :: scale_type_p2c ! scale type for subgrid averaging of pfts to column
    character(len=scale_type_strlen) :: scale_type_c2l ! scale type for subgrid averaging of columns to landunits
    character(len=scale_type_strlen) :: scale_type_l2g ! scale type for subgrid averaging of landunits to gridcells
    type(bounds_type):: bounds         ! boudns
    character(len=16):: l_default      ! local version of 'default'
    character(len=*),parameter :: subname = 'hist_addfld1d'
!------------------------------------------------------------------------

    ! Determine processor bounds

    call get_proc_bounds(bounds)

    ! History buffer pointer

    hpindex = next_history_pointer_index()

    if (present(ptr_lnd)) then
       l_type1d = grlnd
       l_type1d_out = grlnd
       clmptr_rs(hpindex)%ptr => ptr_lnd

    else if (present(ptr_gcell)) then
       l_type1d = nameg
       l_type1d_out = nameg
       clmptr_rs(hpindex)%ptr => ptr_gcell

    else if (present(ptr_lunit)) then
       l_type1d = namel
       l_type1d_out = namel
       clmptr_rs(hpindex)%ptr => ptr_lunit
       if (present(set_lake)) then
          do l = bounds%begl,bounds%endl
             if (lun%lakpoi(l)) ptr_lunit(l) = set_lake
          end do
       end if
       if (present(set_nolake)) then
          do l = bounds%begl,bounds%endl
             if (.not.(lun%lakpoi(l))) ptr_lunit(l) = set_nolake
          end do
       end if
       if (present(set_urb)) then
          do l = bounds%begl,bounds%endl
             if (lun%urbpoi(l)) ptr_lunit(l) = set_urb
          end do
       end if
       if (present(set_nourb)) then
          do l = bounds%begl,bounds%endl
             if (.not.(lun%urbpoi(l))) ptr_lunit(l) = set_nourb
          end do
       end if
       if (present(set_spec)) then
          do l = bounds%begl,bounds%endl
             if (lun%ifspecial(l)) ptr_lunit(l) = set_spec
          end do
       end if

    else if (present(ptr_col)) then
       l_type1d = namec
       l_type1d_out = namec
       clmptr_rs(hpindex)%ptr => ptr_col
       if (present(set_lake)) then
          do c = bounds%begc,bounds%endc
             l =col%landunit(c)
             if (lun%lakpoi(l)) ptr_col(c) = set_lake
          end do
       end if
       if (present(set_nolake)) then
          do c = bounds%begc,bounds%endc
             l =col%landunit(c)
             if (.not.(lun%lakpoi(l))) ptr_col(c) = set_nolake
          end do
       end if
       if (present(set_urb)) then
          do c = bounds%begc,bounds%endc
             l =col%landunit(c)
             if (lun%urbpoi(l)) ptr_col(c) = set_urb
          end do
       end if
       if (present(set_nourb)) then
          do c = bounds%begc,bounds%endc
             l =col%landunit(c)
             if (.not.(lun%urbpoi(l))) ptr_col(c) = set_nourb
          end do
       end if
       if (present(set_spec)) then
          do c = bounds%begc,bounds%endc
             l =col%landunit(c)
             if (lun%ifspecial(l)) ptr_col(c) = set_spec
          end do
       end if
       if (present(set_noglc)) then
          do c = bounds%begc,bounds%endc
             l =col%landunit(c)
             if (.not.(lun%glcpoi(l))) ptr_col(c) = set_noglc
          end do
       endif

    else if (present(ptr_patch)) then
       l_type1d = namep
       l_type1d_out = namep
       clmptr_rs(hpindex)%ptr => ptr_patch
       if (present(set_lake)) then
          do p = bounds%begp,bounds%endp
             l =patch%landunit(p)
             if (lun%lakpoi(l)) ptr_patch(p) = set_lake
          end do
       end if
       if (present(set_nolake)) then
          do p = bounds%begp,bounds%endp
             l =patch%landunit(p)
             if (.not.(lun%lakpoi(l))) ptr_patch(p) = set_nolake
          end do
       end if
       if (present(set_urb)) then
          do p = bounds%begp,bounds%endp
             l =patch%landunit(p)
             if (lun%urbpoi(l)) ptr_patch(p) = set_urb
          end do
       end if
       if (present(set_nourb)) then
          do p = bounds%begp,bounds%endp
             l =patch%landunit(p)
             if (.not.(lun%urbpoi(l))) ptr_patch(p) = set_nourb
          end do
       end if
       if (present(set_spec)) then
          do p = bounds%begp,bounds%endp
             l =patch%landunit(p)
             if (lun%ifspecial(l)) ptr_patch(p) = set_spec
          end do
       end if
       if (present(set_noglc)) then
          do p = bounds%begp,bounds%endp
             l =patch%landunit(p)
             if (.not.(lun%glcpoi(l))) ptr_patch(p) = set_noglc
          end do
       end if
    else
       write(iulog,*) trim(subname),' ERROR: must specify a valid pointer index,', &
          ' choices are [ptr_atm, ptr_lnd, ptr_gcell, ptr_lunit, ptr_col, ptr_patch] '
       call endrun(msg=errMsg(sourcefile, __LINE__))

    end if

    ! Set scaling factor

    scale_type_p2c = 'unity'
    scale_type_c2l = 'unity'
    scale_type_l2g = 'unity'

    if (present(p2c_scale_type)) scale_type_p2c = p2c_scale_type
    if (present(c2l_scale_type)) scale_type_c2l = c2l_scale_type
    if (present(l2g_scale_type)) scale_type_l2g = l2g_scale_type
    if (present(type1d_out)) l_type1d_out = type1d_out

    ! Add field to allhistfldlist

    call allhistfldlist_addfld (fname=trim(fname), numdims=1, type1d=l_type1d, &
          type1d_out=l_type1d_out, type2d=type2d_unset, num2d=1, &
          units=units, avgflag=avgflag, long_name=long_name, hpindex=hpindex, &
          p2c_scale_type=scale_type_p2c, c2l_scale_type=scale_type_c2l, &
          l2g_scale_type=scale_type_l2g)

    l_default = 'active'
    if (present(default)) then
       l_default = default
    end if
    if (trim(l_default) == 'inactive') then
       return
    else
       call allhistfldlist_make_active (name=trim(fname), tape_index=1)
    end if

  end subroutine hist_addfld1d

  !-----------------------------------------------------------------------
  subroutine hist_addfld2d (fname, type2d, units, avgflag, long_name, type1d_out, &
                        ptr_gcell, ptr_lunit, ptr_col, ptr_patch, ptr_lnd, ptr_atm, &
                        p2c_scale_type, c2l_scale_type, l2g_scale_type, &
                        set_lake, set_nolake, set_urb, set_nourb, set_spec, &
                        no_snow_behavior, default)
    !
    ! !DESCRIPTION:
    ! Initialize a single level history field. The pointer inputs, ptr\_*,
    ! point to the appropriate-type array storing the raw history data points.
    ! The value of type1d passed to allhistfldlist\_add\_fld determines which of the
    ! 1d type of the output and the beginning and ending indices the history
    ! buffer field). All fields default to being written to the first history tape
    ! unless 'default' is set to 'inactive'.
    !
    ! !USES:
    use clm_varpar      , only : nlevgrnd, nlevsno, nlevlak, numrad, nlevdecomp_full, nlevcan, nvegwcs,nlevsoi
    use clm_varpar      , only : natpft_size, cft_size, maxpatch_glc, mxsowings, mxharvests
    use landunit_varcon , only : max_lunit
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: fname                      ! field name
    character(len=*), intent(in) :: type2d                     ! 2d output type
    character(len=*), intent(in) :: units                      ! units of field
    character(len=*), intent(in) :: avgflag                    ! time averaging flag
    character(len=*), intent(in) :: long_name                  ! long name of field
    character(len=*), optional, intent(in) :: type1d_out       ! output type (from data type)
    real(r8)        , optional, pointer    :: ptr_atm(:,:)     ! pointer to atm array
    real(r8)        , optional, pointer    :: ptr_lnd(:,:)     ! pointer to lnd array
    real(r8)        , optional, pointer    :: ptr_gcell(:,:)   ! pointer to gridcell array
    real(r8)        , optional, pointer    :: ptr_lunit(:,:)   ! pointer to landunit array
    real(r8)        , optional, pointer    :: ptr_col(:,:)     ! pointer to column array
    real(r8)        , optional, pointer    :: ptr_patch(:,:)     ! pointer to patch array
    real(r8)        , optional, intent(in) :: set_lake         ! value to set lakes to
    real(r8)        , optional, intent(in) :: set_nolake       ! value to set non-lakes to
    real(r8)        , optional, intent(in) :: set_urb          ! value to set urban to
    real(r8)        , optional, intent(in) :: set_nourb        ! value to set non-urban to
    real(r8)        , optional, intent(in) :: set_spec         ! value to set special to
    integer         , optional, intent(in) :: no_snow_behavior ! if a multi-layer snow field, behavior to use for absent snow layers (should be one of the public no_snow_* parameters defined above)
    character(len=*), optional, intent(in) :: p2c_scale_type   ! scale type for subgrid averaging of pfts to column
    character(len=*), optional, intent(in) :: c2l_scale_type   ! scale type for subgrid averaging of columns to landunits
    character(len=*), optional, intent(in) :: l2g_scale_type   ! scale type for subgrid averaging of landunits to gridcells
    character(len=*), optional, intent(in) :: default          ! if set to 'inactive, field will not appear on primary tape
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l,g                 ! indices
    integer :: num2d                   ! size of second dimension (e.g. number of vertical levels)
    integer :: hpindex                 ! history buffer index
    character(len=hist_dim_name_length) :: l_type1d         ! 1d data type
    character(len=hist_dim_name_length) :: l_type1d_out     ! 1d output type
    character(len=scale_type_strlen) :: scale_type_p2c ! scale type for subgrid averaging of pfts to column
    character(len=scale_type_strlen) :: scale_type_c2l ! scale type for subgrid averaging of columns to landunits
    character(len=scale_type_strlen) :: scale_type_l2g ! scale type for subgrid averaging of landunits to gridcells
    type(bounds_type):: bounds
    character(len=16):: l_default      ! local version of 'default'
    character(len=*),parameter :: subname = 'hist_addfld2d'
!------------------------------------------------------------------------

    call get_proc_bounds(bounds)

    ! Error-check no_snow_behavior optional argument: It should be present if and only if
    ! type2d is 'levsno', and its value should be one of the public no_snow_* parameters
    ! defined above.
    if (present(no_snow_behavior)) then
       if (type2d /= 'levsno') then
          write(iulog,*) trim(subname), &
               ' ERROR: Only specify no_snow_behavior for fields with dimension levsno'
          call endrun()
       end if

       if (no_snow_behavior < no_snow_MIN .or. no_snow_behavior > no_snow_MAX) then
          write(iulog,*) trim(subname), &
               ' ERROR: Invalid value for no_snow_behavior: ', no_snow_behavior
          call endrun()
       end if

    else  ! no_snow_behavior is absent
       if (type2d == 'levsno') then
          write(iulog,*) trim(subname), &
               ' ERROR: must specify no_snow_behavior for fields with dimension levsno'
          call endrun()
       end if
    end if

    ! Determine second dimension size

    select case (type2d)
    case ('levgrnd')
       num2d = nlevgrnd
    case ('levsoi')
       num2d = nlevsoi
    case ('levlak')
       num2d = nlevlak
    case ('numrad')
       num2d = numrad
    case ('levdcmp')
       num2d = nlevdecomp_full
    case ('mxsowings')
       num2d = mxsowings
    case ('mxharvests')
       num2d = mxharvests
    case ('fates_levscls')
       num2d = nlevsclass
    case('fates_levcacls')
       num2d = nlevcoage
    case ('fates_levpft')
       num2d = numpft_fates
    case ('fates_levage')
       num2d = nlevage
    case ('fates_levheight')
       num2d = nlevheight
    case ('fates_levfuel')
       num2d = num_fuel_classes
    case ('fates_levcwdsc')
       num2d = ncwd
    case ('fates_levscpf')
       num2d = nlevsclass*numpft_fates
    case ('fates_levcapf')
       num2d = nlevcoage*numpft_fates
    case ('fates_levscag')
       num2d = nlevsclass*nlevage
    case ('fates_levscagpf')
       num2d = nlevsclass*nlevage*numpft_fates
    case ('fates_levagepft')
       num2d = nlevage*numpft_fates
    case ('fates_levcan')
       num2d = nclmax
    case ('fates_levleaf')
       num2d = nlevleaf
    case ('fates_levcnlf')
       num2d = nlevleaf * nclmax
    case ('fates_levcnlfpf')
       num2d = nlevleaf * nclmax * numpft_fates
    case ('fates_levcdsc')
       num2d = nlevdamage * nlevsclass
    case ('fates_levcdpf')
       num2d = nlevdamage * nlevsclass * numpft_fates
    case ('fates_levcdam')
       num2d = nlevdamage
    case ('ltype')
       num2d = max_lunit
    case ('natpft')
       num2d = natpft_size
    case ('fates_levelem')
       num2d = num_elements_fates
    case ('fates_levelpft')
       num2d = num_elements_fates*numpft_fates
    case ('fates_levelcwd')
       num2d = num_elements_fates*ncwd
    case ('fates_levelage')
       num2d = num_elements_fates*nlevage
    case ('fates_levagefuel')
       num2d = nlevage*num_fuel_classes
    case('fates_levclscpf')
       num2d = nclmax * nclmax * numpft_fates
    case ('fates_levlanduse')
       num2d = n_landuse_cats
    case ('fates_levlulu')
       num2d = n_landuse_cats * n_landuse_cats
    case('cft')
       if (cft_size > 0) then
          num2d = cft_size
       else
          write(iulog,*) trim(subname),' ERROR: 2d type =', trim(type2d), &
               ' only valid for cft_size > 0'
          call endrun()
       end if
    case ('glc_nec')
       num2d = maxpatch_glc
    case ('elevclas')
       ! add one because indexing starts at 0 (elevclas, unlike glc_nec, includes the
       ! bare ground "elevation class")
       num2d = maxpatch_glc + 1
    case ('levsno')
       num2d = nlevsno
    case ('nlevcan')
        num2d = nlevcan
    case ('nvegwcs')
        num2d = nvegwcs
    case default
       write(iulog,*) trim(subname),' ERROR: unsupported 2d type ',type2d, &
          ' currently supported types for multi level fields are: ', &
          '[levgrnd,levsoi,levlak,numrad,levdcmp,levtrc,ltype,natpft,cft,glc_nec,elevclas,levsno,nvegwcs,mxsowings,mxharvests]'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end select

    ! History buffer pointer
    hpindex = next_history_pointer_index()


    if (present(ptr_lnd)) then
       l_type1d = grlnd
       l_type1d_out = grlnd
       clmptr_ra(hpindex)%ptr => ptr_lnd

    else if (present(ptr_gcell)) then
       l_type1d = nameg
       l_type1d_out = nameg
       clmptr_ra(hpindex)%ptr => ptr_gcell

    else if (present(ptr_lunit)) then
       l_type1d = namel
       l_type1d_out = namel
       clmptr_ra(hpindex)%ptr => ptr_lunit

       if (present(set_lake)) then
          do l = bounds%begl,bounds%endl
             if (lun%lakpoi(l)) ptr_lunit(l,:) = set_lake
          end do
       end if
       if (present(set_nolake)) then
          do l = bounds%begl,bounds%endl
             if (.not.(lun%lakpoi(l))) ptr_lunit(l,:) = set_nolake
          end do
       end if
       if (present(set_urb)) then
          do l = bounds%begl,bounds%endl
             if (lun%urbpoi(l)) ptr_lunit(l,:) = set_urb
          end do
       end if
       if (present(set_nourb)) then
          do l = bounds%begl,bounds%endl
             if (.not.(lun%urbpoi(l))) ptr_lunit(l,:) = set_nourb
          end do
       end if
       if (present(set_spec)) then
          do l = bounds%begl,bounds%endl
             if (lun%ifspecial(l)) ptr_lunit(l,:) = set_spec
          end do
       end if

    else if (present(ptr_col)) then
       l_type1d = namec
       l_type1d_out = namec
       clmptr_ra(hpindex)%ptr => ptr_col
       if (present(set_lake)) then
          do c = bounds%begc,bounds%endc
             l =col%landunit(c)
             if (lun%lakpoi(l)) ptr_col(c,:) = set_lake
          end do
       end if
       if (present(set_nolake)) then
          do c = bounds%begc,bounds%endc
             l =col%landunit(c)
             if (.not.(lun%lakpoi(l))) ptr_col(c,:) = set_nolake
          end do
       end if
       if (present(set_urb)) then
          do c = bounds%begc,bounds%endc
             l =col%landunit(c)
             if (lun%urbpoi(l)) ptr_col(c,:) = set_urb
          end do
       end if
       if (present(set_nourb)) then
          do c = bounds%begc,bounds%endc
             l =col%landunit(c)
             if (.not.(lun%urbpoi(l))) ptr_col(c,:) = set_nourb
          end do
       end if
       if (present(set_spec)) then
          do c = bounds%begc,bounds%endc
             l =col%landunit(c)
             if (lun%ifspecial(l)) ptr_col(c,:) = set_spec
          end do
       end if

    else if (present(ptr_patch)) then
       l_type1d = namep
       l_type1d_out = namep
       clmptr_ra(hpindex)%ptr => ptr_patch
       
       if (present(set_lake)) then
          do p = bounds%begp,bounds%endp
             l =patch%landunit(p)
             if (lun%lakpoi(l)) ptr_patch(p,:) = set_lake
          end do
       end if
       if (present(set_nolake)) then
          do p = bounds%begp,bounds%endp
             l =patch%landunit(p)
             if (.not.(lun%lakpoi(l))) ptr_patch(p,:) = set_nolake
          end do
       end if
       if (present(set_urb)) then
          do p = bounds%begp,bounds%endp
             l =patch%landunit(p)
             if (lun%urbpoi(l)) ptr_patch(p,:) = set_urb
          end do
       end if
       if (present(set_nourb)) then
          do p = bounds%begp,bounds%endp
             l =patch%landunit(p)
             if (.not.(lun%urbpoi(l))) ptr_patch(p,:) = set_nourb
          end do
       end if
       if (present(set_spec)) then
          do p = bounds%begp,bounds%endp
             l =patch%landunit(p)
             if (lun%ifspecial(l)) ptr_patch(p,:) = set_spec
          end do
       end if

    else
       write(iulog,*) trim(subname),' ERROR: must specify a valid pointer index,', &
          ' choices are ptr_atm, ptr_lnd, ptr_gcell, ptr_lunit, ptr_col, ptr_patch'
       call endrun(msg=errMsg(sourcefile, __LINE__))

    end if

    ! Set scaling factor

    scale_type_p2c = 'unity'
    scale_type_c2l = 'unity'
    scale_type_l2g = 'unity'

    if (present(p2c_scale_type)) scale_type_p2c = p2c_scale_type
    if (present(c2l_scale_type)) scale_type_c2l = c2l_scale_type
    if (present(l2g_scale_type)) scale_type_l2g = l2g_scale_type
    if (present(type1d_out)) l_type1d_out = type1d_out

    ! Add field to allhistfldlist

    call allhistfldlist_addfld (fname=trim(fname), numdims=2, type1d=l_type1d, &
          type1d_out=l_type1d_out, type2d=type2d, num2d=num2d, &
          units=units, avgflag=avgflag, long_name=long_name, hpindex=hpindex, &
          p2c_scale_type=scale_type_p2c, c2l_scale_type=scale_type_c2l, &
          l2g_scale_type=scale_type_l2g, no_snow_behavior=no_snow_behavior)

    l_default = 'active'
    if (present(default)) then
       l_default = default
    end if
    if (trim(l_default) == 'inactive') then
       return
    else
       call allhistfldlist_make_active (name=trim(fname), tape_index=1)
    end if

  end subroutine hist_addfld2d

  !-----------------------------------------------------------------------
  subroutine hist_addfld_decomp (fname, type2d, units, avgflag, long_name, ptr_col, &
       ptr_patch, l2g_scale_type, default)

    !
    ! !DESCRIPTION:
    ! Adds 1-D or 2-D history field based on column data (if ptr_col is present),
    ! patch data (otherwise).
    ! 
    ! !USES:
    use clm_varpar  , only : nlevdecomp_full
    use clm_varctl  , only : iulog
    use abortutils  , only : endrun
    use shr_log_mod , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: fname                    ! field name
    character(len=*), intent(in) :: type2d                   ! 2d output type, if 2d output is chosen
    character(len=*), intent(in) :: units                    ! units of field
    character(len=*), intent(in) :: avgflag                  ! time averaging flag
    character(len=*), intent(in) :: long_name                ! long name of field
    real(r8)        , optional, pointer    :: ptr_col(:,:)   ! pointer to column array
    real(r8)        , optional, pointer    :: ptr_patch(:,:)   ! pointer to patch array
    character(len=*), optional, intent(in) :: l2g_scale_type ! scale type for subgrid averaging of landunits to gridcells
    character(len=*), optional, intent(in) :: default        ! if set to 'inactive, field will not appear on primary tape
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer  :: ptr_1d(:)
    !-----------------------------------------------------------------------

    if (present(ptr_col)) then

       ! column-level data
       if (present(default)) then
          if ( nlevdecomp_full > 1 ) then
             call hist_addfld2d (fname=trim(fname), units=units, type2d=type2d, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_col=ptr_col, l2g_scale_type=l2g_scale_type, default=default)
          else
             ptr_1d => ptr_col(:,1)
             call hist_addfld1d (fname=trim(fname), units=units, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_col=ptr_1d, l2g_scale_type=l2g_scale_type, default=default)
          endif
       else
          if ( nlevdecomp_full > 1 ) then
             call hist_addfld2d (fname=trim(fname), units=units, type2d=type2d, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_col=ptr_col, l2g_scale_type=l2g_scale_type)
          else
             ptr_1d => ptr_col(:,1)
             call hist_addfld1d (fname=trim(fname), units=units, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_col=ptr_1d, l2g_scale_type=l2g_scale_type)
          endif
       endif

    else if (present(ptr_patch)) then

       ! patch-level data
       if (present(default)) then
          if ( nlevdecomp_full > 1 ) then
             call hist_addfld2d (fname=trim(fname), units=units, type2d=type2d, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_patch=ptr_patch, l2g_scale_type=l2g_scale_type, default=default)
          else
             ptr_1d => ptr_patch(:,1)
             call hist_addfld1d (fname=trim(fname), units=units, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_patch=ptr_1d, l2g_scale_type=l2g_scale_type, default=default)
          endif
       else
          if ( nlevdecomp_full > 1 ) then
             call hist_addfld2d (fname=trim(fname), units=units, type2d=type2d, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_patch=ptr_patch, l2g_scale_type=l2g_scale_type)
          else
             ptr_1d => ptr_patch(:,1)
             call hist_addfld1d (fname=trim(fname), units=units, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_patch=ptr_1d, l2g_scale_type=l2g_scale_type)
          endif
       endif

    else
       write(iulog, *) ' error: hist_addfld_decomp needs either patch or column level pointer'
       write(iulog, *) fname
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

  end subroutine hist_addfld_decomp

  !-----------------------------------------------------------------------
  integer function next_history_pointer_index ()
    !
    ! !DESCRIPTION:
    ! Return the next free index in clmptr_r* arrays (for a new history field to write to)
    ! Aka 'hpindex', e.g. field_info.hpindex.
    !
    ! !ARGUMENTS:
    !
    integer, save :: lastindex = 1
    character(len=*),parameter :: subname = 'next_history_pointer_index'
    !-----------------------------------------------------------------------

    next_history_pointer_index = lastindex
    lastindex = lastindex + 1
    if (lastindex > max_mapflds) then
       write(iulog,*) trim(subname),' ERROR: ',&
            ' lastindex = ',lastindex,' greater than max_mapflds= ',max_mapflds
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

  end function next_history_pointer_index

  !-----------------------------------------------------------------------

  subroutine strip_null(str)
    character(len=*), intent(inout) :: str
    integer :: i
    do i=1,len(str)
       if(ichar(str(i:i))==0) str(i:i)=' '
    end do
  end subroutine strip_null

  !------------------------------------------------------------------------
  subroutine hist_do_disp (ntapes, hist_ntimes, hist_mfilt, if_stop, if_disphist, rstwr, nlend)
    !
    ! !DESCRIPTION:
    ! Determine logic for closing and/or disposing history file
    ! Sets values for if_disphist, if_stop (arguments)
    ! Remove history files unless this is end of run or
    ! history file is not full.
    !
    ! !ARGUMENTS:
    integer, intent(in)  :: ntapes              !actual number of history tapes
    integer, intent(in)  :: hist_ntimes(ntapes) !current numbers of time samples on history tape
    integer, intent(in)  :: hist_mfilt(ntapes)  !maximum number of time samples per tape
    logical, intent(out) :: if_stop             !true => last time step of run
    logical, intent(out) :: if_disphist(ntapes) !true => save and dispose history file
    logical, intent(in)  :: rstwr
    logical, intent(in)  :: nlend
    !
    ! !LOCAL VARIABLES:
    integer :: t                   ! history tape index
    logical :: rest_now            ! temporary
    logical :: stop_now            ! temporary
    !------------------------------------------------------------------------

    rest_now = .false.
    stop_now = .false.

    if (nlend) stop_now = .true.
    if (rstwr) rest_now = .true.

    if_stop = stop_now

    if (stop_now) then
       ! End of run -  dispose all history files

       if_disphist(1:ntapes) = .true.

    else if (rest_now) then
       ! Restart - dispose all history files

       do t = 1,ntapes
          if_disphist(t) = .true.
       end do
    else
       ! Dispose

       if_disphist(1:ntapes) = .false.
       do t = 1,ntapes
          if (hist_ntimes(t) ==  hist_mfilt(t)) then
             if_disphist(t) = .true.
          endif
       end do
    endif

  end subroutine hist_do_disp

  !-----------------------------------------------------------------------
  function avgflag_valid(avgflag, blank_valid) result(valid)
    !
    ! !DESCRIPTION:
    ! Returns true if the given avgflag is a valid option, false if not
    !
    ! !USES:
    use clm_varcon      , only : isecspday
    use clm_time_manager, only : get_step_size
    !
    ! !ARGUMENTS:
    logical :: valid  ! function result
    character(len=*), intent(in) :: avgflag
    logical, intent(in) :: blank_valid  ! whether ' ' is a valid avgflag in this context
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'avgflag_valid'
    integer :: tod                      ! Desired local solar time of output in seconds
    integer :: dtime                    ! timestep size [seconds]
    !-----------------------------------------------------------------------

    ! This initial check is mainly here to catch the possibility that someone has added a
    ! new "valid" avgflag option that exceeds avgflag_strlen
    if (len_trim(avgflag) > avgflag_strlen) then
       valid = .false.

    else if (avgflag == ' ' .and. blank_valid) then
       valid = .true.
    else if (avgflag == 'A' .or. avgflag == 'I' .or. &
         avgflag == 'X' .or. avgflag == 'M' .or. &
         avgflag == 'SUM') then
       valid = .true.
    else if (avgflag(1:1) == 'L') then
       dtime = get_step_size()
       if ( len_trim(avgflag) < 6 )then
          valid = .false.
       else
          read(avgflag(2:6), *) tod
          if (tod >= 0 .and. tod <= isecspday) then
             valid = .true.
             if(tod < dtime .or. isecspday - tod <= dtime) then
                write(iulog,*) 'Warning: Local time history output ', avgflag, ' is closer than ', &
                   'dtime to midnight! This is problematic particularly for daily output.'
             end if
          else
             valid = .false.
          end if
       end if
    else
       valid = .false.
    end if

  end function avgflag_valid

  !-----------------------------------------------------------------------
  subroutine add_landunit_mask_metadata(ncid, varid, l2g_scale_type)
    !
    ! !DESCRIPTION:
    ! Add landunit_mask metadata for the given history field
    !
    ! !ARGUMENTS:
    class(file_desc_t), intent(inout) :: ncid  ! netcdf file id
    integer           , intent(in)    :: varid ! netcdf var id
    character(len=*)  , intent(in)    :: l2g_scale_type ! l2g_scale_type for this variable
    !
    ! !LOCAL VARIABLES:
    character(len=:), allocatable :: landunit_mask_string

    character(len=*), parameter :: subname = 'add_landunit_mask_metadata'
    !-----------------------------------------------------------------------

    if (l2g_scale_type == 'unity') then
       ! BUG(wjs, 2021-04-19, ESCOMP/CTSM#1347) Once we consistently set l2g_scale_type
       ! for all variables, and have stopped using other mechanisms (particularly the
       ! setting of variables to spval everywhere) then we can stop setting this to
       ! 'unknown': we can instead set this to something like 'all', with reasonable
       ! confidence that the field truly applies over all landunits.
       landunit_mask_string = 'unknown'
    else
       landunit_mask_string = l2g_scale_type
    end if

    call ncd_putatt(ncid, varid, 'landunit_mask', landunit_mask_string)

  end subroutine add_landunit_mask_metadata

end module histFileMod
