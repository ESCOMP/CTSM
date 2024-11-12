module CropType

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing variables needed for the crop model
  !
  ! TODO(wjs, 2014-08-05) Move more crop-specific variables into here - many are
  ! currently in CNVegStateType
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use spmdMod             , only : masterproc
  use abortutils          , only : endrun
  use decompMod           , only : bounds_type
  use clm_varcon          , only : spval
  use clm_varctl          , only : iulog, use_crop
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC DATA TYPES:
  public :: latbaset
  !

  ! Possible values of cphase
  real(r8), parameter, public :: cphase_not_planted = 0._r8
  real(r8), parameter, public :: cphase_planted     = 1._r8
  real(r8), parameter, public :: cphase_leafemerge  = 2._r8
  real(r8), parameter, public :: cphase_grainfill   = 3._r8
  real(r8), parameter, public :: cphase_harvest     = 4._r8

  ! Crop state variables structure
  type, public :: crop_type

     integer , pointer :: nyrs_crop_active_patch  (:)   ! number of years this crop patch has been active (0 for non-crop patches)
     logical , pointer :: croplive_patch          (:)   ! patch Flag, true if planted, not harvested
     integer , pointer :: harvdate_patch          (:)   ! most recent patch harvest date; 999 if currently (or never) planted
     real(r8), pointer :: fertnitro_patch         (:)   ! patch fertilizer nitrogen
     real(r8), pointer :: gddtsoi_patch           (:)   ! patch growing degree-days from planting (top two soil layers) (ddays)
     real(r8), pointer :: vf_patch                (:)   ! patch vernalization factor for cereal
     real(r8), pointer :: cphase_patch            (:)   ! phenology phase (see cphase_* constants above for possible values)
     integer , pointer :: sowing_reason_patch     (:)   ! reason for most recent sowing of this patch
     real(r8), pointer :: latbaset_patch          (:)   ! Latitude vary baset for hui (degree C)
     character(len=20) :: baset_mapping
     real(r8) :: baset_latvary_intercept
     real(r8) :: baset_latvary_slope
     logical , pointer :: sown_in_this_window           (:)   ! patch flag. True if the crop has already been sown during the current sowing window. False otherwise or if not in a sowing window.
     integer , pointer :: rx_swindow_starts_thisyr_patch(:,:) ! all prescribed sowing window start dates for this patch this year (day of year) [patch, mxsowings]
     integer , pointer :: rx_swindow_ends_thisyr_patch  (:,:) ! all prescribed sowing window end   dates for this patch this year (day of year) [patch, mxsowings]
     real(r8), pointer :: rx_cultivar_gdds_thisyr_patch (:,:) ! all cultivar GDD targets for this patch this year (ddays) [patch, mxsowings]
     real(r8), pointer :: gdd20_baseline_patch          (:)   ! GDD20 baseline for this patch (ddays) [patch]
     real(r8), pointer :: gdd20_season_start_patch(:) ! gdd20 season start date for this patch (day of year) [patch]. Real to enable history field.
     real(r8), pointer :: gdd20_season_end_patch  (:) ! gdd20 season end   date for this patch (day of year) [patch]. Real to enable history field.
     real(r8), pointer :: sdates_thisyr_patch     (:,:) ! all actual sowing dates for this patch this year (day of year) [patch, mxsowings]
     real(r8), pointer :: swindow_starts_thisyr_patch(:,:) ! all sowing window start dates for this patch this year (day of year) [patch, mxsowings]
     real(r8), pointer :: swindow_ends_thisyr_patch  (:,:) ! all sowing window end   dates for this patch this year (day of year) [patch, mxsowings]
     real(r8), pointer :: sdates_perharv_patch    (:,:) ! all actual sowing dates for crops *harvested* this year (day of year) [patch, mxharvests]
     real(r8), pointer :: syears_perharv_patch    (:,:) ! all actual sowing years for crops *harvested* this year (day of year) [patch, mxharvests]
     real(r8), pointer :: hdates_thisyr_patch     (:,:) ! all actual harvest dates for this patch this year (day of year) [patch, mxharvests]
     real(r8), pointer :: gddaccum_thisyr_patch   (:,:) ! accumulated GDD at harvest for this patch this year (ddays) [patch, mxharvests]
     real(r8), pointer :: hui_thisyr_patch        (:,:) ! accumulated heat unit index at harvest for this patch this year (ddays) [patch, mxharvests]
     real(r8), pointer :: sowing_reason_thisyr_patch  (:,:) ! reason for each sowing for this patch this year [patch, mxsowings]
     real(r8), pointer :: sowing_reason_perharv_patch (:,:) ! reason for each sowing of crops *harvested* this year [patch, mxharvests]
     real(r8), pointer :: harvest_reason_thisyr_patch (:,:) ! reason for each harvest for this patch this year [patch, mxharvests]
     integer , pointer :: sowing_count            (:)   ! number of sowing events this year for this patch
     integer , pointer :: harvest_count           (:)   ! number of sowing events this year for this patch
     ! gddaccum tracks the actual growing degree-days accumulated over the growing season.
     ! hui also accumulates growing degree-days, but can be boosted if full leafout is
     ! achieved before the GDD threshold for grain fill has been reached; see CropPhenology().
     real(r8), pointer :: hui_patch               (:)   ! crop patch heat unit index (ddays)
     real(r8), pointer :: gddaccum_patch          (:)   ! patch growing degree-days from planting (air) (ddays)

   contains
     ! Public routines
     procedure, public  :: Init               ! Initialize the crop type
     procedure, public  :: InitAccBuffer
     procedure, public  :: InitAccVars
     procedure, public  :: Restart
     procedure, public  :: ReadNML            ! Read in the crop_inparm namelist

     ! NOTE(wjs, 2014-09-29) need to rename this from UpdateAccVars to CropUpdateAccVars
     ! to prevent cryptic error messages with pgi (v. 13.9 on yellowstone)
     ! This is probably related to this bug
     ! <http://www.pgroup.com/userforum/viewtopic.php?t=4285>, which was fixed in pgi 14.7.
     procedure, public  :: CropUpdateAccVars

     procedure, public  :: CropIncrementYear

     ! Private routines
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
     procedure, private, nopass :: checkDates

  end type crop_type

  character(len=*), parameter, private :: baset_map_constant = 'constant'
  character(len=*), parameter, private :: baset_map_latvary  = 'varytropicsbylat'
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !ARGUMENTS:
    class(crop_type) , intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------

    call this%InitAllocate(bounds)

    if (use_crop) then
       call this%InitHistory(bounds)
       call this%InitCold(bounds)
    end if

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine ReadNML(this, NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for CropType
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    class(crop_type) , intent(inout) :: this
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'Crop::ReadNML'
    character(len=*), parameter :: nmlname = 'crop_inparm'
    !-----------------------------------------------------------------------
    character(len=20) :: baset_mapping
    real(r8) :: baset_latvary_intercept
    real(r8) :: baset_latvary_slope
    namelist /crop_inparm/ baset_mapping, baset_latvary_intercept, baset_latvary_slope

    ! Initialize options to default values, in case they are not specified in
    ! the namelist

    baset_mapping           = 'constant'
    baset_latvary_intercept = 12._r8
    baset_latvary_slope     = 0.4_r8
    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=crop_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (baset_mapping           , mpicom)
    call shr_mpi_bcast (baset_latvary_intercept , mpicom)
    call shr_mpi_bcast (baset_latvary_slope     , mpicom)

    this%baset_mapping           = baset_mapping
    this%baset_latvary_intercept = baset_latvary_intercept
    this%baset_latvary_slope     = baset_latvary_slope
    if (      trim(this%baset_mapping) == baset_map_constant ) then
       if ( masterproc ) write(iulog,*) 'baset mapping for ALL crops are constant'
    else if ( trim(this%baset_mapping) == baset_map_latvary ) then
       if ( masterproc ) write(iulog,*) 'baset mapping for crops vary with latitude'
    else
       call endrun(msg="Bad value for baset_mapping in "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
    end if

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=crop_inparm)
       write(iulog,*) ' '
    end if

    !-----------------------------------------------------------------------

  end subroutine ReadNML

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    ! !USES:
    !
    use clm_varpar, only : mxsowings, mxharvests
    !
    ! !ARGUMENTS:
    class(crop_type) , intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp

    character(len=*), parameter :: subname = 'InitAllocate'
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    allocate(this%nyrs_crop_active_patch(begp:endp)) ; this%nyrs_crop_active_patch(:) = 0
    allocate(this%croplive_patch (begp:endp)) ; this%croplive_patch (:) = .false.
    allocate(this%harvdate_patch (begp:endp)) ; this%harvdate_patch (:) = huge(1)
    allocate(this%fertnitro_patch (begp:endp)) ; this%fertnitro_patch (:) = spval
    allocate(this%hui_patch (begp:endp))      ; this%hui_patch      (:) = spval
    allocate(this%gddaccum_patch (begp:endp)) ; this%gddaccum_patch (:) = spval
    allocate(this%gddtsoi_patch  (begp:endp)) ; this%gddtsoi_patch  (:) = spval
    allocate(this%vf_patch       (begp:endp)) ; this%vf_patch       (:) = 0.0_r8
    allocate(this%cphase_patch   (begp:endp)) ; this%cphase_patch   (:) = cphase_not_planted
    allocate(this%sowing_reason_patch (begp:endp)) ; this%sowing_reason_patch (:) = -1
    allocate(this%latbaset_patch (begp:endp)) ; this%latbaset_patch (:) = spval
    allocate(this%sown_in_this_window(begp:endp)) ; this%sown_in_this_window(:) = .false.
    allocate(this%rx_swindow_starts_thisyr_patch(begp:endp,1:mxsowings)); this%rx_swindow_starts_thisyr_patch(:,:) = -1
    allocate(this%rx_swindow_ends_thisyr_patch(begp:endp,1:mxsowings))  ; this%rx_swindow_ends_thisyr_patch  (:,:) = -1
    allocate(this%rx_cultivar_gdds_thisyr_patch(begp:endp,1:mxsowings)) ; this%rx_cultivar_gdds_thisyr_patch(:,:) = spval
    allocate(this%gdd20_baseline_patch(begp:endp)) ; this%gdd20_baseline_patch(:) = spval
    allocate(this%gdd20_season_start_patch(begp:endp)); this%gdd20_season_start_patch(:) = spval
    allocate(this%gdd20_season_end_patch(begp:endp))  ; this%gdd20_season_end_patch  (:) = spval
    allocate(this%sdates_thisyr_patch(begp:endp,1:mxsowings)) ; this%sdates_thisyr_patch(:,:) = spval
    allocate(this%swindow_starts_thisyr_patch(begp:endp,1:mxsowings)) ; this%swindow_starts_thisyr_patch(:,:) = spval
    allocate(this%swindow_ends_thisyr_patch  (begp:endp,1:mxsowings)) ; this%swindow_ends_thisyr_patch  (:,:) = spval
    allocate(this%sdates_perharv_patch(begp:endp,1:mxharvests)) ; this%sdates_perharv_patch(:,:) = spval
    allocate(this%syears_perharv_patch(begp:endp,1:mxharvests)) ; this%syears_perharv_patch(:,:) = spval
    allocate(this%hdates_thisyr_patch(begp:endp,1:mxharvests)) ; this%hdates_thisyr_patch(:,:) = spval
    allocate(this%gddaccum_thisyr_patch(begp:endp,1:mxharvests)) ; this%gddaccum_thisyr_patch(:,:) = spval
    allocate(this%hui_thisyr_patch(begp:endp,1:mxharvests)) ; this%hui_thisyr_patch(:,:) = spval
    allocate(this%sowing_reason_thisyr_patch(begp:endp,1:mxsowings)) ; this%sowing_reason_thisyr_patch(:,:) = spval
    allocate(this%sowing_reason_perharv_patch(begp:endp,1:mxharvests)) ; this%sowing_reason_perharv_patch(:,:) = spval
    allocate(this%harvest_reason_thisyr_patch(begp:endp,1:mxharvests)) ; this%harvest_reason_thisyr_patch(:,:) = spval
    allocate(this%sowing_count(begp:endp)) ; this%sowing_count(:) = 0
    allocate(this%harvest_count(begp:endp)) ; this%harvest_count(:) = 0

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod    , only : hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(crop_type),  intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp

    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    this%fertnitro_patch(begp:endp) = spval
    call hist_addfld1d (fname='FERTNITRO', units='gN/m2/yr', &
         avgflag='A', long_name='Nitrogen fertilizer for each crop', &
         ptr_patch=this%fertnitro_patch, default='inactive')

    this%hui_patch(begp:endp) = spval
    call hist_addfld1d (fname='HUI', units='ddays', &
         avgflag='A', long_name='Crop patch heat unit index', &
         ptr_patch=this%hui_patch, default='inactive')

    this%gddaccum_patch(begp:endp) = spval
    call hist_addfld1d (fname='GDDACCUM', units='ddays', &
         avgflag='A', long_name='Accumulated growing degree days past planting date for crop', &
         ptr_patch=this%gddaccum_patch, default='inactive')

    this%gddtsoi_patch(begp:endp) = spval
    call hist_addfld1d (fname='GDDTSOI', units='ddays', &
         avgflag='A', long_name='Growing degree-days from planting (top two soil layers)', &
         ptr_patch=this%gddtsoi_patch, default='inactive')

    this%cphase_patch(begp:endp) = spval
    call hist_addfld1d (fname='CPHASE', units='0-not planted, 1-planted, 2-leaf emerge, 3-grain fill, 4-harvest', &
         avgflag='A', long_name='crop phenology phase', &
         ptr_patch=this%cphase_patch, default='active')

    if ( (trim(this%baset_mapping) == baset_map_latvary) )then
       this%latbaset_patch(begp:endp) = spval
       call hist_addfld1d (fname='LATBASET', units='degree C', &
            avgflag='A', long_name='latitude vary base temperature for hui', &
            ptr_patch=this%latbaset_patch, default='inactive')
    end if

    this%sdates_thisyr_patch(begp:endp,:) = spval
    call hist_addfld2d (fname='SDATES', units='day of year', type2d='mxsowings', &
         avgflag='I', long_name='actual crop sowing dates; should only be output annually', &
         ptr_patch=this%sdates_thisyr_patch, default='inactive')

    this%swindow_starts_thisyr_patch(begp:endp,:) = spval
    call hist_addfld2d (fname='SWINDOW_STARTS', units='day of year', type2d='mxsowings', &
         avgflag='I', long_name='crop sowing window start dates; should only be output annually', &
         ptr_patch=this%swindow_starts_thisyr_patch, default='inactive')

    this%swindow_ends_thisyr_patch(begp:endp,:) = spval
    call hist_addfld2d (fname='SWINDOW_ENDS', units='day of year', type2d='mxsowings', &
         avgflag='I', long_name='crop sowing window end dates; should only be output annually', &
         ptr_patch=this%swindow_ends_thisyr_patch, default='inactive')

    this%sdates_perharv_patch(begp:endp,:) = spval
    call hist_addfld2d (fname='SDATES_PERHARV', units='day of year', type2d='mxharvests', &
         avgflag='I', long_name='actual sowing dates for crops harvested this year; should only be output annually', &
         ptr_patch=this%sdates_perharv_patch, default='inactive')

    this%syears_perharv_patch(begp:endp,:) = spval
    call hist_addfld2d (fname='SYEARS_PERHARV', units='year', type2d='mxharvests', &
         avgflag='I', long_name='actual sowing years for crops harvested this year; should only be output annually', &
         ptr_patch=this%syears_perharv_patch, default='inactive')

    this%hdates_thisyr_patch(begp:endp,:) = spval
    call hist_addfld2d (fname='HDATES', units='day of year', type2d='mxharvests', &
         avgflag='I', long_name='actual crop harvest dates; should only be output annually', &
         ptr_patch=this%hdates_thisyr_patch, default='inactive')

    this%gddaccum_thisyr_patch(begp:endp,:) = spval
    call hist_addfld2d (fname='GDDACCUM_PERHARV', units='ddays', type2d='mxharvests', &
         avgflag='I', long_name='At-harvest accumulated growing degree days past planting date for crop; should only be output annually', &
         ptr_patch=this%gddaccum_thisyr_patch, default='inactive')

    this%hui_thisyr_patch(begp:endp,:) = spval
    call hist_addfld2d (fname='HUI_PERHARV', units='ddays', type2d='mxharvests', &
         avgflag='I', long_name='At-harvest accumulated heat unit index for crop; should only be output annually', &
         ptr_patch=this%hui_thisyr_patch, default='inactive')

    this%sowing_reason_thisyr_patch(begp:endp,:) = spval
    call hist_addfld2d (fname='SOWING_REASON', units='unitless', type2d='mxsowings', &
         avgflag='I', long_name='Reason for each crop sowing; should only be output annually', &
         ptr_patch=this%sowing_reason_thisyr_patch, default='inactive')

    this%sowing_reason_perharv_patch(begp:endp,:) = spval
    call hist_addfld2d (fname='SOWING_REASON_PERHARV', units='unitless', type2d='mxharvests', &
         avgflag='I', long_name='Reason for sowing of each crop harvested this year; should only be output annually', &
         ptr_patch=this%sowing_reason_perharv_patch, default='inactive')

    this%harvest_reason_thisyr_patch(begp:endp,:) = spval
    call hist_addfld2d (fname='HARVEST_REASON_PERHARV', units='1 = mature; 2 = max season length; 3 = incorrect Dec. 31 '// &
    'sowing; 4 = sowing today; 5 = sowing tomorrow; 6 = tomorrow == idop; 7 = killed by cold temperature during vernalization', &
         type2d='mxharvests', &
         avgflag='I', long_name='Reason for each crop harvest; should only be output annually', &
         ptr_patch=this%harvest_reason_thisyr_patch, default='inactive')

    this%gdd20_baseline_patch(begp:endp) = spval
    call hist_addfld1d (fname='GDD20_BASELINE', units='ddays', &
         avgflag='I', long_name='Baseline mean growing-degree days accumulated during accumulation period (from input)', &
         ptr_patch=this%gdd20_baseline_patch, default='inactive')
    this%gdd20_season_start_patch(begp:endp) = spval
    call hist_addfld1d (fname='GDD20_SEASON_START', units='day of year', &
         avgflag='I', long_name='Start of the GDD20 accumulation season (from input)', &
         ptr_patch=this%gdd20_season_start_patch, default='inactive')
    this%gdd20_season_end_patch(begp:endp) = spval
    call hist_addfld1d (fname='GDD20_SEASON_END', units='day of year', &
         avgflag='I', long_name='End of the GDD20 accumulation season (from input)', &
         ptr_patch=this%gdd20_season_end_patch, default='inactive')

  end subroutine InitHistory

  subroutine InitCold(this, bounds)
    ! !USES:
    use LandunitType, only : lun
    use landunit_varcon, only : istcrop
    use PatchType, only : patch
    use clm_instur, only : fert_cft
    use pftconMod        , only : pftcon
    use GridcellType     , only : grc
    use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
    ! !ARGUMENTS:
    class(crop_type),  intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: l, g, p, ivt ! indices
    logical :: latvary_baset

    character(len=*), parameter :: subname = 'InitCold'
    !-----------------------------------------------------------------------

    latvary_baset = trim(this%baset_mapping) == baset_map_latvary
    if (.not. latvary_baset) then
        this%latbaset_patch(bounds%begp:bounds%endp) = nan
    end if

    do p= bounds%begp,bounds%endp
       l = patch%landunit(p)

       this%nyrs_crop_active_patch(p) = 0

       if (lun%itype(l) == istcrop) then
          g = patch%gridcell(p)
          ivt = patch%itype(p)
          this%fertnitro_patch(p) = fert_cft(g,ivt)

          if (latvary_baset) then
              this%latbaset_patch(p) = latbaset(pftcon%baset(ivt), grc%latdeg(g), this%baset_latvary_intercept, this%baset_latvary_slope)
          end if
       end if
    end do

  end subroutine InitCold

  !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
  subroutine InitAccBuffer (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for all required module accumulated fields
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    ! Each interval and accumulation type is unique to each field processed.
    ! Routine [initAccBuffer] defines the fields to be processed
    ! and the type of accumulation.
    ! Routine [updateAccVars] does the actual accumulation for a given field.
    ! Fields are accumulated by calls to subroutine [update_accum_field].
    ! To accumulate a field, it must first be defined in subroutine [initAccVars]
    ! and then accumulated by calls to [updateAccVars].
    !
    ! Should only be called if use_crop is true
    !
    ! !USES
    use accumulMod       , only : init_accum_field
    !
    ! !ARGUMENTS:
    class(crop_type) , intent(in) :: this
    type(bounds_type), intent(in) :: bounds

    !
    ! !LOCAL VARIABLES:
    integer, parameter :: not_used = huge(1)

    !---------------------------------------------------------------------

    ! BACKWARDS_COMPATIBILITY (ssr/wjs, 2022-03-16): old_name specification
    ! allows correct restarting from files that used GDDPLANT.
    call init_accum_field (name='HUI', units='K', &
         desc='heat unit index accumulated since planting', accum_type='runaccum', accum_period=not_used,  &
         subgrid_type='pft', numlev=1, init_value=0._r8, old_name='GDDPLANT')

    call init_accum_field (name='GDDACCUM', units='K', &
         desc='growing degree-days from planting', accum_type='runaccum', accum_period=not_used,  &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field (name='GDDTSOI', units='K', &
         desc='growing degree-days from planting (top two soil layers)', accum_type='runaccum', accum_period=not_used,  &
         subgrid_type='pft', numlev=1, init_value=0._r8)

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES:
    use accumulMod       , only : extract_accum_field
    use clm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(crop_type),  intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begp, endp
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslp(:)  ! temporary

    character(len=*), parameter :: subname = 'InitAccVars'
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    ! Allocate needed dynamic memory for single level patch field
    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)' in '
       call endrun(msg=" allocation error for rbufslp"//&
            errMsg(sourcefile, __LINE__))
    endif

    nstep = get_nstep()

    call extract_accum_field ('HUI', rbufslp, nstep)
    this%hui_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('GDDACCUM', rbufslp, nstep)
    this%gddaccum_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('GDDTSOI', rbufslp, nstep)
    this%gddtsoi_patch(begp:endp)  = rbufslp(begp:endp)

    deallocate(rbufslp)

  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, cnveg_state_inst, flag)
    !
    ! !USES:
    use restUtilMod
    use ncdio_pio
    use PatchType, only : patch
    use pftconMod, only : npcropmin, npcropmax
    use clm_varpar, only : mxsowings, mxharvests
    ! BACKWARDS_COMPATIBILITY(wjs/ssr, 2023-01-09)
    use CNVegstateType, only : cnveg_state_type
    use clm_time_manager , only : get_curr_calday, get_curr_date
    !
    ! !ARGUMENTS:
    class(crop_type), intent(inout)  :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid
    type(cnveg_state_type) , intent(inout) :: cnveg_state_inst ! BACKWARDS_COMPATIBILITY(wjs/ssr, 2023-01-09)
    character(len=*) , intent(in)    :: flag
    !
    ! !LOCAL VARIABLES:
    integer, pointer :: temp1d(:) ! temporary
    integer :: restyear
    integer :: p
    logical :: readvar   ! determine if variable is on initial file
    integer :: seasons_found, seasons_loopvar      ! getting number of sowings/harvests in patch
    ! BACKWARDS_COMPATIBILITY(wjs/ssr, 2023-01-09)
    integer jday      ! julian day of the year
    integer kyr       ! current year
    integer kmo       ! month of year  (1, ..., 12)
    integer kda       ! day of month   (1, ..., 31)
    integer mcsec     ! seconds of day (0, ..., seconds/day)
    ! BACKWARDS_COMPATIBILITY(ssr, 2023-01-13)
    logical read_hdates_thisyr_patch

    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------

    if (use_crop) then
       call restartvar(ncid=ncid, flag=flag, varname='nyrs_crop_active', xtype=ncd_int, &
            dim1name='pft', &
            long_name='Number of years this crop patch has been active (0 for non-crop patches)', &
            units='years', &
            interpinic_flag='interp', readvar=readvar, data=this%nyrs_crop_active_patch)
       if (flag == 'read' .and. .not. readvar) then
          ! BACKWARDS_COMPATIBILITY(wjs, 2017-02-17) Old restart files did not have this
          ! patch-level variable. Instead, they had a single scalar tracking the number
          ! of years the crop model ran. Copy this scalar onto all *active* crop patches.

          ! Some arguments in the following restartvar call are irrelevant, because we
          ! only call this for 'read'. I'm simply maintaining the old restartvar call.
          call restartvar(ncid=ncid, flag=flag,  varname='restyear', xtype=ncd_int,  &
               long_name='Number of years prognostic crop ran', units="years", &
               interpinic_flag='copy', readvar=readvar, data=restyear)
          if (readvar) then
             do p = bounds%begp, bounds%endp
                if (patch%itype(p) >= npcropmin .and. patch%itype(p) <= npcropmax .and. &
                     patch%active(p)) then
                   this%nyrs_crop_active_patch(p) = restyear
                end if
             end do
          end if
       end if

       allocate(temp1d(bounds%begp:bounds%endp))
       if (flag == 'write') then
          do p= bounds%begp,bounds%endp
             if (this%croplive_patch(p)) then
                temp1d(p) = 1
             else
                temp1d(p) = 0
             end if
          end do
       end if
       call restartvar(ncid=ncid, flag=flag,  varname='croplive', xtype=ncd_log,  &
            dim1name='pft', &
            long_name='Flag that crop is alive, but not harvested', &
            interpinic_flag='interp', readvar=readvar, data=temp1d)
       if (flag == 'read') then
          do p= bounds%begp,bounds%endp
             if (temp1d(p) == 1) then
                this%croplive_patch(p) = .true.
             else
                this%croplive_patch(p) = .false.
             end if
          end do
       end if
       deallocate(temp1d)

       allocate(temp1d(bounds%begp:bounds%endp))
       if (flag == 'write') then
          do p= bounds%begp,bounds%endp
             if (this%sown_in_this_window(p)) then
                temp1d(p) = 1
             else
                temp1d(p) = 0
             end if
          end do
       end if
       call restartvar(ncid=ncid, flag=flag,  varname='sown_in_this_window', xtype=ncd_log,  &
            dim1name='pft', &
            long_name='Flag that patch was sown already during the current sowing window', &
            interpinic_flag='interp', readvar=readvar, data=temp1d)
       if (flag == 'read') then
          do p= bounds%begp,bounds%endp
             if (temp1d(p) == 1) then
                this%sown_in_this_window(p) = .true.
             else
                this%sown_in_this_window(p) = .false.
             end if
          end do
       end if
       deallocate(temp1d)

       call restartvar(ncid=ncid, flag=flag,  varname='harvdate', xtype=ncd_int,  &
            dim1name='pft', long_name='harvest date', units='jday', nvalid_range=(/1,366/), &
            interpinic_flag='interp', readvar=readvar, data=this%harvdate_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='vf', xtype=ncd_double,  &
            dim1name='pft', long_name='vernalization factor', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%vf_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='cphase',xtype=ncd_double, &
            dim1name='pft', long_name='crop phenology phase', &
            units='0-not planted, 1-planted, 2-leaf emerge, 3-grain fill, 4-harvest', &
            interpinic_flag='interp', readvar=readvar, data=this%cphase_patch)
       if (flag=='read' )then
          call this%checkDates( )  ! Check that restart date is same calendar date (even if year is different)
                                   ! This is so that it properly goes through
                                   ! the crop phases
       end if

       call restartvar(ncid=ncid, flag=flag,  varname='sowing_reason_patch',xtype=ncd_int, &
            dim1name='pft', long_name='sowing reason for this patch', &
            units='none', &
            interpinic_flag='interp', readvar=readvar, data=this%sowing_reason_patch)

       ! Read or write variable(s) with mxsowings dimension
       ! BACKWARDS_COMPATIBILITY(wjs/ssr, 2022-02-02) See note in CallRestartvarDimOK()
       if (CallRestartvarDimOK(ncid, flag, 'mxsowings')) then
           call restartvar(ncid=ncid, flag=flag, varname='sdates_thisyr_patch', xtype=ncd_double,  &
                dim1name='pft', dim2name='mxsowings', switchdim=.true., &
                long_name='crop sowing dates for this patch this year', units='day of year', &
                scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=this%sdates_thisyr_patch)
           call restartvar(ncid=ncid, flag=flag, varname='swindow_starts_thisyr_patch', xtype=ncd_double,  &
                dim1name='pft', dim2name='mxsowings', switchdim=.true., &
                long_name='sowing window start dates for this patch this year', units='day of year', &
                scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=this%swindow_starts_thisyr_patch)
           call restartvar(ncid=ncid, flag=flag, varname='swindow_ends_thisyr_patch', xtype=ncd_double,  &
                dim1name='pft', dim2name='mxsowings', switchdim=.true., &
                long_name='sowing window end dates for this patch this year', units='day of year', &
                scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=this%swindow_ends_thisyr_patch)
           ! Fill variable(s) derived from read-in variable(s)
           if (flag == 'read' .and. readvar) then
             do p = bounds%begp,bounds%endp
                seasons_found = 0
                do seasons_loopvar = 1,mxsowings
                   if (this%sdates_thisyr_patch(p,seasons_loopvar) >= 1 .and. this%sdates_thisyr_patch(p,seasons_loopvar) <= 366) then
                      seasons_found = seasons_loopvar
                   else
                      exit
                   end if
                end do ! loop through possible sowings
                this%sowing_count(p) = seasons_found
             end do ! loop through patches
           end if
       end if

       ! Read or write variable(s) with mxharvests dimension
       ! BACKWARDS_COMPATIBILITY(wjs/ssr, 2022-02-02) See note in CallRestartvarDimOK()
       if (CallRestartvarDimOK(ncid, flag, 'mxharvests')) then
           call restartvar(ncid=ncid, flag=flag, varname='sdates_perharv_patch', xtype=ncd_double,  &
                dim1name='pft', dim2name='mxharvests', switchdim=.true., &
                long_name='sowing dates for crops harvested in this patch this year', units='day of year', &
                scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=this%sdates_perharv_patch)
           call restartvar(ncid=ncid, flag=flag, varname='syears_perharv_patch', xtype=ncd_double,  &
                dim1name='pft', dim2name='mxharvests', switchdim=.true., &
                long_name='sowing years for crops harvested in this patch this year', units='year', &
                scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=this%syears_perharv_patch)
           call restartvar(ncid=ncid, flag=flag, varname='hdates_thisyr_patch', xtype=ncd_double,  &
                dim1name='pft', dim2name='mxharvests', switchdim=.true., &
                long_name='crop harvest dates for this patch this year', units='day of year', &
                scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=read_hdates_thisyr_patch, data=this%hdates_thisyr_patch)
           call restartvar(ncid=ncid, flag=flag, varname='gddaccum_thisyr_patch', xtype=ncd_double,  &
                dim1name='pft', dim2name='mxharvests', switchdim=.true., &
                long_name='accumulated GDD at harvest for this patch this year', units='ddays', &
                scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=this%gddaccum_thisyr_patch)
           call restartvar(ncid=ncid, flag=flag, varname='hui_thisyr_patch', xtype=ncd_double,  &
                dim1name='pft', dim2name='mxharvests', switchdim=.true., &
                long_name='accumulated heat unit index at harvest for this patch this year', units='ddays', &
                scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=this%hui_thisyr_patch)
           call restartvar(ncid=ncid, flag=flag, varname='sowing_reason_thisyr_patch', xtype=ncd_double,  &
                dim1name='pft', dim2name='mxsowings', switchdim=.true., &
                long_name='reason for each sowing for this patch this year', units='unitless', &
                scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=this%sowing_reason_thisyr_patch)
           call restartvar(ncid=ncid, flag=flag, varname='sowing_reason_perharv_patch', xtype=ncd_double,  &
                dim1name='pft', dim2name='mxharvests', switchdim=.true., &
                long_name='reason for sowing of each crop harvested this year', units='unitless', &
                scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=this%sowing_reason_perharv_patch)
           call restartvar(ncid=ncid, flag=flag, varname='harvest_reason_thisyr_patch', xtype=ncd_double,  &
                dim1name='pft', dim2name='mxharvests', switchdim=.true., &
                long_name='reason for each harvest for this patch this year', units='unitless', &
                scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=this%harvest_reason_thisyr_patch)

           ! Fill variable(s) derived from read-in variable(s)
           if (flag == 'read') then
             jday = get_curr_calday()
             call get_curr_date(kyr, kmo, kda, mcsec)
             do p = bounds%begp,bounds%endp

                ! Harvest count
                if (read_hdates_thisyr_patch) then
                   seasons_found = 0
                   do seasons_loopvar = 1,mxharvests
                      if (this%hdates_thisyr_patch(p,seasons_loopvar) >= 1 .and. this%hdates_thisyr_patch(p,seasons_loopvar) <= 366) then
                         seasons_found = seasons_loopvar
                      else
                         exit
                      end if
                   end do ! loop through possible harvests
                   this%harvest_count(p) = seasons_found
                end if

                ! Year of planting
                ! Calculating this here instead of saving in restart file to allow for
                ! sensible iyop values in startup/hybrid runs.
                ! * Assumes no growing season is longer than 364 days (or 365 days if
                !   spanning a leap day).
                if (cnveg_state_inst%idop_patch(p) <= jday) then
                    cnveg_state_inst%iyop_patch(p) = kyr
                else
                    cnveg_state_inst%iyop_patch(p) = kyr - 1
                end if
             end do ! loop through patches
           end if
       end if

       ! BACKWARDS_COMPATIBILITY(wjs/ssr, 2023-01-09)
       if (flag == 'read') then
           do p = bounds%begp,bounds%endp
               ! Will be needed until we can rely on all restart files including sowing_reason_patch.
               if (this%croplive_patch(p) .and. this%sowing_reason_patch(p) < 0) then
                  this%sowing_reason_patch(p) = 0
               end if
           end do ! loop through patches
       end if
    end if

  end subroutine Restart


  !-----------------------------------------------------------------------
  subroutine CropUpdateAccVars(this, bounds, t_ref2m_patch, t_soisno_col)
    !
    ! !DESCRIPTION:
    ! Update accumulated variables. Should be called every time step.
    ! Should only be called if use_crop is true.
    !
    ! !USES:
    use accumulMod       , only : update_accum_field, extract_accum_field, markreset_accum_field
    use shr_const_mod    , only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use clm_time_manager , only : get_step_size, get_nstep
    use clm_varpar       , only : nlevsno, nlevmaxurbgrnd
    use pftconMod        , only : nswheat, nirrig_swheat, pftcon
    use pftconMod        , only : nwwheat, nirrig_wwheat
    use pftconMod        , only : nsugarcane, nirrig_sugarcane
    use ColumnType       , only : col
    use PatchType        , only : patch
    !
    ! !ARGUMENTS:
    implicit none
    class(crop_type)       , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    real(r8)               , intent(in)    :: t_ref2m_patch( bounds%begp:)
    real(r8)               , intent(inout) :: t_soisno_col(bounds%begc:, -nlevsno+1:)
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,g ! indices
    integer :: ivt   ! vegetation type
    integer :: dtime ! timestep size [seconds]
    integer :: nstep ! timestep number
    integer :: ier   ! error status
    integer :: begp, endp
    integer :: begc, endc
    real(r8), pointer :: rbufslp(:)      ! temporary single level - patch level
    character(len=*), parameter :: subname = 'CropUpdateAccVars'
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(t_ref2m_patch)  == (/endp/))          , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(t_soisno_col)   == (/endc,nlevmaxurbgrnd/)) , sourcefile, __LINE__)

    dtime = get_step_size()
    nstep = get_nstep()

    ! Allocate needed dynamic memory for single level patch field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)'update_accum_hist allocation error for rbuf1dp'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

    ! Update HUI. This is not standard for accumulation fields,
    ! but HUI needs it because it can be changed outside this
    ! accumulation routine (see CropPhenology). This requires
    ! the accumulation buffer to be reset.
    call extract_accum_field ('HUI', rbufslp, nstep)
    do p = begp,endp
      rbufslp(p) = max(0.0_r8,this%hui_patch(p)-rbufslp(p))
    end do
    call update_accum_field  ('HUI', rbufslp, nstep)

    ! Accumulate and extract HUI and GDDACCUM
    do p = begp,endp
       if (this%croplive_patch(p)) then ! relative to planting date
          ivt = patch%itype(p)
          if ( (trim(this%baset_mapping) == baset_map_latvary) .and. &
             ((ivt == nswheat) .or. (ivt == nirrig_swheat) .or. &
              (ivt == nsugarcane) .or. (ivt == nirrig_sugarcane)) ) then
             rbufslp(p) = max(0._r8, min(pftcon%mxtmp(ivt), &
             t_ref2m_patch(p)-(SHR_CONST_TKFRZ + this%latbaset_patch(p)))) &
             * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = max(0._r8, min(pftcon%mxtmp(ivt), &
             t_ref2m_patch(p)-(SHR_CONST_TKFRZ + pftcon%baset(ivt)))) &
             * dtime/SHR_CONST_CDAY
          end if
          if (ivt == nwwheat .or. ivt == nirrig_wwheat) then
             rbufslp(p) = rbufslp(p) * this%vf_patch(p)
          end if
       else
          call markreset_accum_field('HUI', p)
          call markreset_accum_field('GDDACCUM', p)
       end if
    end do
    call update_accum_field  ('HUI', rbufslp, nstep)
    call extract_accum_field ('HUI', this%hui_patch, nstep)
    call update_accum_field  ('GDDACCUM', rbufslp, nstep)
    call extract_accum_field ('GDDACCUM', this%gddaccum_patch, nstep)

    ! Accumulate and extract GDDTSOI
    ! In agroibis this variable is calculated
    ! to 0.05 m, so here we use the top two soil layers

    do p = begp,endp
       if (this%croplive_patch(p)) then ! relative to planting date
          ivt = patch%itype(p)
          c   = patch%column(p)
          rbufslp(p) = max(0._r8, min(pftcon%mxtmp(ivt), &
               ((t_soisno_col(c,1)*col%dz(c,1) + &
               t_soisno_col(c,2)*col%dz(c,2))/(col%dz(c,1)+col%dz(c,2))) - &
               (SHR_CONST_TKFRZ + pftcon%baset(ivt)))) * dtime/SHR_CONST_CDAY
          if (ivt == nwwheat .or. ivt == nwwheat) then
             rbufslp(p) = rbufslp(p) * this%vf_patch(p)
          end if
       else
          call markreset_accum_field('GDDTSOI', p)
       end if
    end do
    call update_accum_field  ('GDDTSOI', rbufslp, nstep)
    call extract_accum_field ('GDDTSOI', this%gddtsoi_patch, nstep)

    deallocate(rbufslp)

  end subroutine CropUpdateAccVars

  !-----------------------------------------------------------------------
  subroutine CropIncrementYear (this, num_pcropp, filter_pcropp)
    !
    ! !DESCRIPTION:
    ! Increment the crop year, if appropriate
    !
    ! This routine should be called every time step
    !
    ! !USES:
    use clm_time_manager , only : get_curr_date, is_first_step
    !
    ! !ARGUMENTS:
    class(crop_type) :: this
    integer , intent(in) :: num_pcropp       ! number of prog. crop patches in filter
    integer , intent(in) :: filter_pcropp(:) ! filter for prognostic crop patches
    !
    ! !LOCAL VARIABLES:
    integer kyr   ! current year
    integer kmo   ! month of year  (1, ..., 12)
    integer kda   ! day of month   (1, ..., 31)
    integer mcsec ! seconds of day (0, ..., seconds/day)
    integer :: fp, p
    !-----------------------------------------------------------------------

    call get_curr_date (   kyr, kmo, kda, mcsec)
    ! Update nyrs when it's the end of the year (unless it's the very start of the
    ! run). This assumes that, if this patch is active at the end of the year, then it was
    ! active for the whole year.
    if ((kmo == 1 .and. kda == 1 .and. mcsec == 0)) then
       do fp = 1, num_pcropp
          p = filter_pcropp(fp)

          this%nyrs_crop_active_patch(p) = this%nyrs_crop_active_patch(p) + 1
       end do
    end if

  end subroutine CropIncrementYear

  !-----------------------------------------------------------------------
  subroutine checkDates( )
    !
    ! !DESCRIPTION:
    ! Make sure the dates are compatible. The date given to startup the model
    ! and the date on the restart file must be the same although years can be
    ! different. The dates need to be checked when the restart file is being
    ! read in for a startup or branch case (they are NOT allowed to be different
    ! for a restart case).
    !
    ! For the prognostic crop model the date of planting is tracked and growing
    ! degree days is tracked (with a 20 year mean) -- so shifting the start dates
    ! messes up these bits of saved information.
    !
    ! !ARGUMENTS:
    use clm_time_manager, only : get_driver_start_ymd, get_start_date
    use clm_varctl      , only : iulog
    use clm_varctl      , only : nsrest, nsrBranch, nsrStartup
    !
    ! !LOCAL VARIABLES:
    integer :: stymd       ! Start date YYYYMMDD from driver
    integer :: styr        ! Start year from driver
    integer :: stmon_day   ! Start date MMDD from driver
    integer :: rsmon_day   ! Restart date MMDD from restart file
    integer :: rsyr        ! Restart year from restart file
    integer :: rsmon       ! Restart month from restart file
    integer :: rsday       ! Restart day from restart file
    integer :: tod         ! Restart time of day from restart file
    character(len=*), parameter :: formDate = '(A,i4.4,"/",i2.2,"/",i2.2)' ! log output format
    character(len=32) :: subname = 'CropRest::checkDates'
    !-----------------------------------------------------------------------
    !
    ! If branch or startup make sure the startdate is compatible with the date
    ! on the restart file.
    !
    if ( nsrest == nsrBranch .or. nsrest == nsrStartup )then
       stymd       = get_driver_start_ymd()
       styr        = stymd / 10000
       stmon_day   = stymd - styr*10000
       call get_start_date( rsyr, rsmon, rsday, tod )
       rsmon_day = rsmon*100 + rsday
       if ( masterproc ) &
            write(iulog,formDate) 'Date on the restart file is: ', rsyr, rsmon, rsday
       if ( stmon_day /= rsmon_day )then
          write(iulog,formDate) 'Start date is: ', styr, stmon_day/100, &
               (stmon_day - stmon_day/100)
          call endrun(msg=' ERROR: For prognostic crop to work correctly, the start date (month and day)'// &
               ' and the date on the restart file needs to match (years can be different)'//&
               errMsg(sourcefile, __LINE__))
       end if
    end if

  end subroutine checkDates

  real(r8) function latbaset(baset, latdeg, baset_latvary_intercept, baset_latvary_slope)
    ! !ARGUMENTS:
    real(r8), intent(in) :: baset
    real(r8), intent(in) :: latdeg
    real(r8), intent(in) :: baset_latvary_intercept
    real(r8), intent(in) :: baset_latvary_slope

    ! Was originally
    !     maxlat = baset_latvary_intercept / baset_latvary_slope
    !     if (abs(latdeg) > maxlat) then
    !         latbaset = baset
    !     else
    !         latbaset = baset + baset_latvary_intercept - baset_latvary_slope*abs(latdeg)
    !     end if
    ! But the one-liner below should improve efficiency, at least marginally.

    latbaset = baset + baset_latvary_intercept - min(baset_latvary_intercept, baset_latvary_slope * abs(latdeg))

  end function latbaset

end module CropType
