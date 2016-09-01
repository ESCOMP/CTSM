module subgridMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! sub-grid data and mapping types and modules
  !
  ! TODO(wjs, 2015-12-08) Much of the logic here duplicates (in some sense) logic in
  ! initGridCellsMod. The duplication should probably be extracted into routines shared
  ! between these modules (or the two modules should be combined into one).
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use spmdMod        , only : masterproc
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog
  use clm_instur     , only : wt_lunit, urban_valid, wt_cft
  use glcBehaviorMod , only : glc_behavior_type
  use EDtypesMod, only : cohorts_per_col

  implicit none
  private   
  save

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: subgrid_get_gcellinfo   ! Obtain gridcell properties, summed across all landunits

  ! Routines to get info for each landunit:
  public :: subgrid_get_info_natveg
  public :: subgrid_get_info_urban_tbd
  public :: subgrid_get_info_urban_hd
  public :: subgrid_get_info_urban_md
  public :: subgrid_get_info_lake
  public :: subgrid_get_info_wetland
  public :: subgrid_get_info_glacier
  public :: subgrid_get_info_glacier_mec
  public :: subgrid_get_info_crop
  public :: crop_patch_exists ! returns true if the given crop patch should be created in memory

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: subgrid_get_info_urban

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine subgrid_get_gcellinfo (gi, glc_behavior, &
       nlunits, ncols, npatches, ncohorts)
    !
    ! !DESCRIPTION:
    ! Obtain gridcell properties, aggregated across all landunits
    !
    ! !ARGUMENTS
    integer , intent(in)  :: gi       ! grid cell index
    type(glc_behavior_type), intent(in) :: glc_behavior
    integer , intent(out) :: nlunits  ! number of landunits
    integer , intent(out) :: ncols    ! number of columns 
    integer , intent(out) :: npatches ! number of patchs 
    integer , intent(out) :: ncohorts ! number of cohorts 
    !
    ! !LOCAL VARIABLES:
    ! Counts from a single landunit:
    integer :: ncohorts_temp
    integer :: npatches_temp
    integer :: ncols_temp
    integer :: nlunits_temp

    ! atm_topo is arbitrary for the sake of getting these counts. We don't have a true
    ! atm_topo value at the point of this call, so use 0.
    real(r8), parameter :: atm_topo = 0._r8
    !------------------------------------------------------------------------------

    npatches = 0
    ncols    = 0
    nlunits  = 0
    ncohorts = 0

    call subgrid_get_info_natveg(gi, ncohorts, npatches_temp, ncols_temp, nlunits_temp)
    call accumulate_counters()

    call subgrid_get_info_urban_tbd(gi, npatches_temp, ncols_temp, nlunits_temp)
    call accumulate_counters()

    call subgrid_get_info_urban_hd(gi, npatches_temp, ncols_temp, nlunits_temp)
    call accumulate_counters()

    call subgrid_get_info_urban_md(gi, npatches_temp, ncols_temp, nlunits_temp)
    call accumulate_counters()

    call subgrid_get_info_lake(gi, npatches_temp, ncols_temp, nlunits_temp)
    call accumulate_counters()

    call subgrid_get_info_wetland(gi, npatches_temp, ncols_temp, nlunits_temp)
    call accumulate_counters()

    call subgrid_get_info_glacier(gi, npatches_temp, ncols_temp, nlunits_temp)
    call accumulate_counters()

    call subgrid_get_info_glacier_mec(gi, atm_topo, glc_behavior, &
         npatches_temp, ncols_temp, nlunits_temp)
    call accumulate_counters()

    call subgrid_get_info_crop(gi, npatches_temp, ncols_temp, nlunits_temp)
    call accumulate_counters()
   

  contains
    subroutine accumulate_counters
      ! Accumulate running sums of patches, columns and landunits.
      !
      ! This uses local variables in the parent subroutine as both inputs and outputs

      npatches = npatches + npatches_temp
      ncols = ncols + ncols_temp
      nlunits = nlunits + nlunits_temp

    end subroutine accumulate_counters

  end subroutine subgrid_get_gcellinfo

  !-----------------------------------------------------------------------
  subroutine subgrid_get_info_natveg(gi, ncohorts, npatches, ncols, nlunits)
    !
    ! !DESCRIPTION:
    ! Obtain properties for natural vegetated landunit in this grid cell
    !
    ! !USES
    use clm_varpar, only : natpft_size
    !
    ! !ARGUMENTS:
    integer, intent(in)  :: gi        ! grid cell index
    integer, intent(out) :: ncohorts  ! number of nat veg cohorts in this grid cell
    integer, intent(out) :: npatches  ! number of nat veg patches in this grid cell
    integer, intent(out) :: ncols     ! number of nat veg columns in this grid cell
    integer, intent(out) :: nlunits   ! number of nat veg landunits in this grid cell
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'subgrid_get_info_natveg'
    !-----------------------------------------------------------------------

    ! To support dynamic landunits, we have a naturally vegetated landunit in every grid
    ! cell, because it might need to come into existence even if its weight is 0 at the
    ! start of the run. And to support transient patches or dynamic vegetation, we always
    ! allocate space for ALL patches on this landunit.

    npatches = natpft_size

    ! Assume that the vegetated landunit has one column
    nlunits = 1
    ncols = 1

    ! -------------------------------------------------------------------------
    ! Number of cohorts is set here
    ! ED cohorts (via FATES) populate all natural vegetation columns.
    ! Current implementations mostly assume that only one column contains
    ! natural vegetation, which is synonomous with the soil column. 
    ! For restart output however, we will allocate the cohort vector space
    ! based on all columns.
    ! -------------------------------------------------------------------------

    ncohorts = ncols*cohorts_per_col

  end subroutine subgrid_get_info_natveg

  !-----------------------------------------------------------------------
  subroutine subgrid_get_info_urban_tbd(gi, npatches, ncols, nlunits)
    !
    ! !DESCRIPTION:
    ! Obtain properties for urban tbd landunit in this grid cell
    !
    ! !ARGUMENTS:
    integer, intent(in)  :: gi        ! grid cell index
    integer, intent(out) :: npatches  ! number of urban tbd patches in this grid cell
    integer, intent(out) :: ncols     ! number of urban tbd columns in this grid cell
    integer, intent(out) :: nlunits   ! number of urban tbd landunits in this grid cell
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'subgrid_get_info_urban_tbd'
    !-----------------------------------------------------------------------

    call subgrid_get_info_urban(gi, npatches, ncols, nlunits)

  end subroutine subgrid_get_info_urban_tbd

  !-----------------------------------------------------------------------
  subroutine subgrid_get_info_urban_hd(gi, npatches, ncols, nlunits)
    !
    ! !DESCRIPTION:
    ! Obtain properties for urban hd landunit in this grid cell
    !
    ! !ARGUMENTS:
    integer, intent(in)  :: gi        ! grid cell index
    integer, intent(out) :: npatches  ! number of urban hd patches in this grid cell
    integer, intent(out) :: ncols     ! number of urban hd columns in this grid cell
    integer, intent(out) :: nlunits   ! number of urban hd landunits in this grid cell
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'subgrid_get_info_urban_hd'
    !-----------------------------------------------------------------------

    call subgrid_get_info_urban(gi, npatches, ncols, nlunits)

  end subroutine subgrid_get_info_urban_hd

  !-----------------------------------------------------------------------
  subroutine subgrid_get_info_urban_md(gi, npatches, ncols, nlunits)
    !
    ! !DESCRIPTION:
    ! Obtain properties for urban md landunit in this grid cell
    !
    ! !ARGUMENTS:
    integer, intent(in)  :: gi        ! grid cell index
    integer, intent(out) :: npatches  ! number of urban md patches in this grid cell
    integer, intent(out) :: ncols     ! number of urban md columns in this grid cell
    integer, intent(out) :: nlunits   ! number of urban md landunits in this grid cell
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'subgrid_get_info_urban_md'
    !-----------------------------------------------------------------------

    call subgrid_get_info_urban(gi, npatches, ncols, nlunits)

  end subroutine subgrid_get_info_urban_md

  !-----------------------------------------------------------------------
  subroutine subgrid_get_info_urban(gi, npatches, ncols, nlunits)
    !
    ! !DESCRIPTION:
    ! Obtain properties for one of the urban landunits in this grid cell
    !
    ! This is shared for all urban landunits, because currently they are all treated the same.
    !
    ! !USES
    use clm_varpar, only : maxpatch_urb
    !
    ! !ARGUMENTS:
    integer, intent(in)  :: gi        ! grid cell index
    integer, intent(out) :: npatches  ! number of urban patches in this grid cell, for one urban landunit
    integer, intent(out) :: ncols     ! number of urban columns in this grid cell, for one urban landunit
    integer, intent(out) :: nlunits   ! number of urban landunits in this grid cell, for one urban landunit
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'subgrid_get_info_urban'
    !-----------------------------------------------------------------------

    ! To support dynamic landunits, we have all urban landunits in every grid cell that
    ! has valid urban parameters, because they might need to come into existence even if
    ! their weight is 0 at the start of the run. And for simplicity, we always allocate
    ! space for ALL columns on the urban landunits.

    if (urban_valid(gi)) then
       npatches = maxpatch_urb
       ncols = npatches
       nlunits = 1
    else
       npatches = 0
       ncols = 0
       nlunits = 0
    end if

  end subroutine subgrid_get_info_urban

  !-----------------------------------------------------------------------
  subroutine subgrid_get_info_lake(gi, npatches, ncols, nlunits)
    !
    ! !DESCRIPTION:
    ! Obtain properties for lake landunit in this grid cell
    !
    ! !USES:
    use landunit_varcon, only : istdlak
    !
    ! !ARGUMENTS:
    integer, intent(in)  :: gi        ! grid cell index
    integer, intent(out) :: npatches  ! number of lake patches in this grid cell
    integer, intent(out) :: ncols     ! number of lake columns in this grid cell
    integer, intent(out) :: nlunits   ! number of lake landunits in this grid cell
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'subgrid_get_info_lake'
    !-----------------------------------------------------------------------

    ! We currently do NOT allow the lake landunit to expand via dynamic landunits, so we
    ! only need to allocate space for it where its weight is currently non-zero.

    if (wt_lunit(gi, istdlak) > 0.0_r8) then
       npatches = 1
       ncols = 1
       nlunits = 1
    else
       npatches = 0
       ncols = 0
       nlunits = 0
    end if

  end subroutine subgrid_get_info_lake

  !-----------------------------------------------------------------------
  subroutine subgrid_get_info_wetland(gi, npatches, ncols, nlunits)
    !
    ! !DESCRIPTION:
    ! Obtain properties for wetland landunit in this grid cell
    !
    ! !USES:
    use landunit_varcon, only : istwet
    !
    ! !ARGUMENTS:
    integer, intent(in)  :: gi        ! grid cell index
    integer, intent(out) :: npatches  ! number of wetland patches in this grid cell
    integer, intent(out) :: ncols     ! number of wetland columns in this grid cell
    integer, intent(out) :: nlunits   ! number of wetland landunits in this grid cell
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'subgrid_get_info_wetland'
    !-----------------------------------------------------------------------

    ! We currently do NOT allow the wetland landunit to expand via dynamic landunits, so we
    ! only need to allocate space for it where its weight is currently non-zero.

    if (wt_lunit(gi, istwet) > 0.0_r8) then
       npatches = 1
       ncols = 1
       nlunits = 1
    else
       npatches = 0
       ncols = 0
       nlunits = 0
    end if

  end subroutine subgrid_get_info_wetland
  
  !-----------------------------------------------------------------------
  subroutine subgrid_get_info_glacier(gi, npatches, ncols, nlunits)
    !
    ! !DESCRIPTION:
    ! Obtain properties for glacier landunit in this grid cell
    !
    ! !USES:
    use landunit_varcon, only : istice
    !
    ! !ARGUMENTS:
    integer, intent(in)  :: gi        ! grid cell index
    integer, intent(out) :: npatches  ! number of glacier patches in this grid cell
    integer, intent(out) :: ncols     ! number of glacier columns in this grid cell
    integer, intent(out) :: nlunits   ! number of glacier landunits in this grid cell
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'subgrid_get_info_glacier'
    !-----------------------------------------------------------------------

    ! We currently do NOT allow the glacier landunit to expand via dynamic landunits, so
    ! we only need to allocate space for it where its weight is currently non-zero. (If we
    ! have dynamic glacier area, we will be using glacier_mec landunits rather than
    ! glacier landunits.)

    if (wt_lunit(gi, istice) > 0.0_r8) then
       npatches = 1
       ncols = 1
       nlunits = 1
    else
       npatches = 0
       ncols = 0
       nlunits = 0
    end if

  end subroutine subgrid_get_info_glacier

    !-----------------------------------------------------------------------
  subroutine subgrid_get_info_glacier_mec(gi, atm_topo, glc_behavior, npatches, ncols, nlunits)
    !
    ! !DESCRIPTION:
    ! Obtain properties for glacier_mec landunit in this grid cell
    !
    ! !ARGUMENTS:
    integer, intent(in)  :: gi        ! grid cell index
    real(r8), intent(in) :: atm_topo  ! atmosphere's topographic height for this grid cell (m)
    type(glc_behavior_type), intent(in) :: glc_behavior
    integer, intent(out) :: npatches  ! number of glacier_mec patches in this grid cell
    integer, intent(out) :: ncols     ! number of glacier_mec columns in this grid cell
    integer, intent(out) :: nlunits   ! number of glacier_mec landunits in this grid cell
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'subgrid_get_info_glacier_mec'
    !-----------------------------------------------------------------------

    call glc_behavior%get_num_glc_mec_subgrid(gi, atm_topo, npatches, ncols, nlunits)

  end subroutine subgrid_get_info_glacier_mec

  !-----------------------------------------------------------------------
  subroutine subgrid_get_info_crop(gi, npatches, ncols, nlunits)
    !
    ! !DESCRIPTION:
    ! Obtain properties for crop landunit in this grid cell
    !
    ! !USES:
    use clm_varpar, only : cft_lb, cft_ub
    !
    ! !ARGUMENTS:
    integer, intent(in)  :: gi        ! grid cell index
    integer, intent(out) :: npatches  ! number of nat veg patches in this grid cell
    integer, intent(out) :: ncols     ! number of nat veg columns in this grid cell
    integer, intent(out) :: nlunits   ! number of nat veg landunits in this grid cell
    !
    ! !LOCAL VARIABLES:
    integer :: cft  ! crop functional type index

    character(len=*), parameter :: subname = 'subgrid_get_info_crop'
    !-----------------------------------------------------------------------

    npatches = 0

    do cft = cft_lb, cft_ub
       if (crop_patch_exists(gi, cft)) then
          npatches = npatches + 1
       end if
    end do

    if (npatches > 0) then
       ncols = npatches
       nlunits = 1
    else
       ncols = 0
       nlunits = 0
    end if

  end subroutine subgrid_get_info_crop

  !-----------------------------------------------------------------------
  function crop_patch_exists(gi, cft) result(exists)
    !
    ! !DESCRIPTION:
    ! Returns true if a patch should be created in memory for the given crop functional
    ! type in this grid cell.
    !
    ! This just applies to the crop landunit: it always returns .false. if
    ! create_crop_landunit is .false.
    !
    ! !USES:
    use clm_varpar           , only : cft_lb, cft_ub
    use clm_varctl           , only : create_crop_landunit
    use pftconmod            , only : pftcon
    use landunit_varcon      , only : istcrop
    use dynSubgridControlMod , only : get_do_transient_crops
    !
    ! !ARGUMENTS:
    logical :: exists  ! function result
    integer, intent(in) :: gi  ! grid cell index
    integer, intent(in) :: cft ! crop functional type
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'crop_patch_exists'
    !-----------------------------------------------------------------------

    if (create_crop_landunit) then
       SHR_ASSERT(cft >= cft_lb, errMsg(sourcefile, __LINE__))
       SHR_ASSERT(cft <= cft_ub, errMsg(sourcefile, __LINE__))

       if (get_do_transient_crops()) then
          ! To support dynamic landunits, we have all possible crop columns in every grid
          ! cell, because they might need to come into existence even if their weight is 0 at
          ! the start of the run.
          if (pftcon%is_pft_known_to_model(cft)) then
             exists = .true.
          else
             exists = .false.
          end if

       else
          ! For a run without transient crops, only allocate memory for crops that are
          ! actually present in this run. (This will require running init_interp when
          ! changing between a transient crop run and a non-transient run.)
          if (wt_lunit(gi, istcrop) > 0.0_r8 .and. wt_cft(gi, cft) > 0.0_r8) then
             exists = .true.
          else
             exists = .false.
          end if
       end if

    else  ! create_crop_landunit false
       exists = .false.
    end if

  end function crop_patch_exists



end module subgridMod
