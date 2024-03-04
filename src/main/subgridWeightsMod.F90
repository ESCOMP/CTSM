module subgridWeightsMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handles modifications, error-checks and diagnostics related to changing subgrid weights
  !
  ! ----- Requirements for subgrid weights that are enforced here -----
  !
  ! (These requirements are checked in check_weights/weights_okay)
  !
  ! Note: in the following, 'active' refers to a pft, column, landunit or grid cell over
  ! which computations are performed, and 'inactive' refers to a pft, column or landunit
  ! where computations are NOT performed (grid cells are always active).
  !
  ! (1) For all columns, landunits and grid cells, the sum of all subgrid weights of its
  !     children (or grandchildren, etc.) is equal to 1. For example:
  !     - For all columns, the sum of all patch weights on the column equals 1
  !     - For all landunits, the sum of all col weights on the landunit equals 1
  !     - For all grid cells, the sum of all patch weights on the grid cell equals 1
  !     - etc.
  !
  ! (2) For all ACTIVE columns, landunits and grid cells, the sum of all subgrid weights of
  !     its ACTIVE children (or grandchildren, etc.) is equal to 1. For example:
  !     - For all active columns, the sum of all patch weights on the column equals 1 when
  !       just considering active pfts
  !     - For all active landunits, the sum of all col weights on the landunit equals 1 when
  !       just considering active cols
  !     - For ALL grid cells, the sum of all patch weights on the grid cell equals 1 when
  !       just considering active pfts -- note that all grid cells are considered active!
  !     - etc.
  !
  ! (3) For all INACTIVE columns, landunits and grid cells, the sum of all subgrid weights of
  !     its ACTIVE children, grandchildren, etc. are equal to either 0 or 1. For example:
  !     - For all inactive columns, the sum of all patch weights on the column equals either 0
  !       or 1 when just considering active pfts
  !     - For all inactive landunits, the sum of all col weights on the landunit equals
  !       either 0 or 1 when just considering active cols
  !     - etc.
  !
  ! Another way of stating (2) and (3) is that the sum of weights of all ACTIVE pfts, cols
  ! or landunits on their parent/grandparent/etc. is always equal to either 0 or 1 -- and
  ! must be equal to 1 if this parent/grandparent, etc. is itself active.
  !
  ! Note that, together, conditions (1) and (2) imply that any pft, col or landunit whose
  ! weight on the grid cell is non-zero must be active. In addition, these conditions imply
  ! that any patch whose weight on the column is non-zero must be active if the column is
  ! active (and similarly for any patch on an active landunit, and any col on an active
  ! landunit).
  !
  !
  ! ----- Implications of these requirements for computing subgrid averages -----
  !
  ! The preferred way to average from, say, patch to col is:
  !    colval(c) = 0
  !    do p = pfti(c), pftf(c)
  !       if (active(p)) colval(c) = colval(c) + pftval(p) * wtcol(p)
  ! (where wtcol(p) is the weight of the patch on the column)
  ! If column c is active, then the above conditions guarantee that the pwtcol values
  ! included in the above sum will sum to 1. If column c is inactive, then the above
  ! conditions guarantee that the pwtcol values included in the above sum will sum to
  ! either 1 or 0; if they sum to 0, then colval(c) will remain 0.
  !
  ! Another acceptable method is the following; this method accommodates some unknown
  ! fraction of pftval's being set to spval, and leaves colval set at spval if there are no
  ! valid patch values:
  !    colval(c) = spval
  !    sumwt(c) = 0
  !    do p = pfti(c), pftf(c)
  !       if (active(p) .and. wtcol(p) /= 0) then
  !          if (pftval(p) /= spval) then
  !             if (sumwt(c) == 0) colval(c) = 0
  !             colval(c) = colval(c) + pftval(p) * wtcol(p)
  !             sumwt(c) = sumwt(c) + wtcol(p)
  !          end if
  !       end if
  !    end do
  !    if (sumwt(c) /= 0) then
  !       colval(c) = colval(c) / sumwt(c)
  !    end if
  ! Note that here we check the condition (active(p) .and. wtcol(p) /= 0). We need to
  ! include a check for wtcol(p) /= 0 because we don't want to set colval(c) = 0 for zero-
  ! weight pfts in this line:
  !             if (sumwt(c) == 0) colval(c) = 0
  ! And we include a check for active(p) because we don't want to assume that pftval(p) has
  ! been set to spval for inactive pfts -- we want to allow for the possibility that
  ! pftval(p) will be NaN for inactive pfts.
  !
  !
  ! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use abortutils   , only : endrun
  use clm_varctl   , only : iulog, all_active, run_zero_weight_urban, use_fates, use_fates_sp
  use decompMod    , only : bounds_type, subgrid_level_landunit, subgrid_level_column, subgrid_level_patch
  use GridcellType , only : grc
  use LandunitType , only : lun
  use ColumnType   , only : col
  use PatchType    , only : patch
  use glcBehaviorMod , only : glc_behavior_type
  !
  ! PUBLIC TYPES:
  implicit none
  save

  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: init_subgrid_weights_mod      ! initialize stuff in this module
  public :: compute_higher_order_weights  ! given p2c, c2l and l2g weights, compute other weights
  public :: set_active                    ! set 'active' flags at pft, column & landunit level
  public :: check_weights                 ! check subgrid weights
  public :: get_landunit_weight           ! get the weight of a given landunit on a single grid cell
  public :: set_landunit_weight           ! set the weight of a given landunit on a single grid cell
  public :: is_gcell_all_ltypeX           ! determine whether a grid cell is 100% covered by the given landunit type
  public :: set_subgrid_diagnostic_fields ! set all subgrid weights diagnostic fields
  !
  ! !REVISION HISTORY:
  ! Created by Bill Sacks
  !
  ! !PRIVATE TYPES:
  type subgrid_weights_diagnostics_type
     ! This type contains diagnostics on subgrid weights, for output to the history file
     real(r8), pointer :: pct_landunit(:,:)  ! % of each landunit on the grid cell [begg:endg, 1:max_lunit]
     real(r8), pointer :: pct_nat_pft(:,:)   ! % of each pft, as % of landunit [begg:endg, natpft_lb:natpft_ub]
     real(r8), pointer :: pct_cft(:,:)       ! % of each crop functional type, as % of landunit [begg:endg, cft_lb:cft_ub]
     real(r8), pointer :: pct_glc_mec(:,:)   ! % of each glacier elevation class, as % of landunit [begg:endg, 1:maxpatch_glc]
  end type subgrid_weights_diagnostics_type

  type(subgrid_weights_diagnostics_type) :: subgrid_weights_diagnostics

  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: is_active_l                  ! determine whether the given landunit is active
  private :: is_active_c                  ! determine whether the given column is active
  private :: is_active_p                  ! determine whether the given patch is active
  private :: weights_okay                 ! determine if sum of weights satisfies requirements laid out above
  private :: set_pct_landunit_diagnostics ! set pct_landunit diagnostic field
  private :: set_pct_glc_mec_diagnostics  ! set pct_glc_mec diagnostic field
  private :: set_pct_pft_diagnostics      ! set pct_nat_pft & pct_cft diagnostic fields

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine init_subgrid_weights_mod(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize stuff in this module
    !
    ! !USES:
    use landunit_varcon, only : max_lunit
    use clm_varpar     , only : maxpatch_glc, natpft_size, cft_size
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use decompMod      , only : bounds_level_proc
    use histFileMod    , only : hist_addfld2d
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'init_subgrid_weights_mod'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL(bounds%level == bounds_level_proc, sourcefile, __LINE__)

    ! ------------------------------------------------------------------------
    ! Allocate variables in subgrid_weights_diagnostics
    ! ------------------------------------------------------------------------

    ! Note that, because these variables are output to the history file, it appears that
    ! their lower bounds need to start at 1 (e.g., 1:natpft_size rather than
    ! natpft_lb:natpft_ub)
    allocate(subgrid_weights_diagnostics%pct_landunit(bounds%begg:bounds%endg, 1:max_lunit))
    subgrid_weights_diagnostics%pct_landunit(:,:) = nan
    allocate(subgrid_weights_diagnostics%pct_nat_pft(bounds%begg:bounds%endg, 1:natpft_size))
    subgrid_weights_diagnostics%pct_nat_pft(:,:) = nan
    allocate(subgrid_weights_diagnostics%pct_cft(bounds%begg:bounds%endg, 1:cft_size))
    subgrid_weights_diagnostics%pct_cft(:,:) = nan
    allocate(subgrid_weights_diagnostics%pct_glc_mec(bounds%begg:bounds%endg, 1:maxpatch_glc))
    subgrid_weights_diagnostics%pct_glc_mec(:,:) = nan

    ! ------------------------------------------------------------------------
    ! Add history fields
    ! ------------------------------------------------------------------------

    call hist_addfld2d (fname='PCT_LANDUNIT', units='%', type2d='ltype', &
         avgflag='A', long_name='% of each landunit on grid cell', &
         ptr_lnd=subgrid_weights_diagnostics%pct_landunit)

    if(.not.use_fates.or.use_fates_sp) then
       call hist_addfld2d (fname='PCT_NAT_PFT', units='%', type2d='natpft', &
             avgflag='A', long_name='% of each PFT on the natural vegetation (i.e., soil) landunit', &
             ptr_lnd=subgrid_weights_diagnostics%pct_nat_pft)
    end if

    if (cft_size > 0) then
       call hist_addfld2d (fname='PCT_CFT', units='%', type2d='cft', &
            avgflag='A', long_name='% of each crop on the crop landunit', &
            ptr_lnd=subgrid_weights_diagnostics%pct_cft)
    end if

    call hist_addfld2d (fname='PCT_GLC_MEC', units='%', type2d='glc_nec', &
         avgflag='A', long_name='% of each GLC elevation class on the glacier landunit', &
         ptr_lnd=subgrid_weights_diagnostics%pct_glc_mec)

  end subroutine init_subgrid_weights_mod


  !-----------------------------------------------------------------------
  subroutine compute_higher_order_weights(bounds)
    !
    ! !DESCRIPTION:
    ! Assuming patch%wtcol, col%wtlunit and lun%wtgcell have already been computed, compute
    ! the "higher-order" weights: patch%wtlunit, patch%wtgcell and col%wtgcell, for all p and c
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! clump bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p, c, l      ! indices for pft, col & landunit
    !------------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       col%wtgcell(c) = col%wtlunit(c) * lun%wtgcell(l)
    end do

    do p = bounds%begp, bounds%endp
       c = patch%column(p)
       patch%wtlunit(p) = patch%wtcol(p) * col%wtlunit(c)
       patch%wtgcell(p) = patch%wtcol(p) * col%wtgcell(c)
    end do
  end subroutine compute_higher_order_weights

  !-----------------------------------------------------------------------
  subroutine set_active(bounds, glc_behavior)
    !
    ! !DESCRIPTION:
    ! Set 'active' flags at the pft, column and landunit level
    ! (note that grid cells are always active)
    !
    ! This should be called whenever any weights change (e.g., patch weights on the column,
    ! landunit weights on the grid cell, etc.).
    !
    ! Ensures that we don't have any active patch on an inactive column, or an active column on an
    ! inactive landunit (since these conditions could lead to garbage data)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    type(glc_behavior_type), intent(in) :: glc_behavior
    !
    ! !LOCAL VARIABLES:
    integer :: l,c,p       ! loop counters

    character(len=*), parameter :: subname = 'set_active'
    !------------------------------------------------------------------------

    do l = bounds%begl,bounds%endl
       lun%active(l) = is_active_l(l, glc_behavior)
    end do

    do c = bounds%begc,bounds%endc
       l = col%landunit(c)
       col%active(c) = is_active_c(c, glc_behavior)
       if (col%active(c) .and. .not. lun%active(l)) then
          write(iulog,*) trim(subname),' ERROR: active column found on inactive landunit', &
                         'at c = ', c, ', l = ', l
          call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, msg=errMsg(sourcefile, __LINE__))
       end if
    end do

    do p = bounds%begp,bounds%endp
       c = patch%column(p)
       patch%active(p) = is_active_p(p)
       if (patch%active(p) .and. .not. col%active(c)) then
          write(iulog,*) trim(subname),' ERROR: active patch found on inactive column', &
                         'at p = ', p, ', c = ', c
          call endrun(subgrid_index=p, subgrid_level=subgrid_level_patch, msg=errMsg(sourcefile, __LINE__))
       end if
    end do

  end subroutine set_active

  !-----------------------------------------------------------------------
  logical function is_active_l(l, glc_behavior)
    !
    ! !DESCRIPTION:
    ! Determine whether the given landunit is active
    !
    ! !USES:
    use landunit_varcon, only : istsoil, istice, isturb_MIN, isturb_MAX, istdlak
    use clm_instur     , only : pct_urban_max
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: l   ! landunit index
    type(glc_behavior_type), intent(in) :: glc_behavior
    !
    ! !LOCAL VARIABLES:
    integer :: g  ! grid cell index
    integer :: dens_index ! urban density index
    !------------------------------------------------------------------------

    if (all_active) then
       is_active_l = .true.

    else
       g =lun%gridcell(l)

       is_active_l = .false.

       ! ------------------------------------------------------------------------
       ! General conditions under which is_active_l NEEDS to be true in order to satisfy
       ! the requirements laid out at the top of this module:
       ! ------------------------------------------------------------------------
       if (lun%wtgcell(l) > 0) is_active_l = .true.

       ! ------------------------------------------------------------------------
       ! Conditions under which is_active_p is set to true because we want extra virtual landunits:
       ! ------------------------------------------------------------------------

       if (lun%itype(l) == istice .and. &
            glc_behavior%has_virtual_columns_grc(g)) then
          is_active_l = .true.
       end if

       ! Set urban land units to active, as long as memory has been allocated for such land units, either
       ! through the run_zero_weight_urban setting which runs all urban landunits in each grid cell or
       ! through pct_urban_max which is the maximum percent urban for each density type in a transient run.
       ! (See subgridMod.F90 for this logic).
       ! By doing this, urban land units are also run virtually in grid cells which will grow
       ! urban during the transient run.

       if ( (lun%itype(l) >= isturb_MIN .and. lun%itype(l) <= isturb_MAX) ) then
          is_active_l = .true.
       end if

       ! In general, include a virtual natural vegetation landunit. This aids
       ! initialization of a new landunit; and for runs that are coupled to CISM, this
       ! provides bare land SMB forcing even if there is no vegetated area.
       !
       ! Also (echoing the similar comment in glcBehaviorMod): We need all glacier and
       ! vegetated points to be active in the icemask region for the sake of init_interp -
       ! since we only interpolate onto active points, and we don't know which points will
       ! have non-zero area until after initialization (as long as we can't send
       ! information from glc to clm in initialization). (If we had an inactive vegetated
       ! point in the icemask region, according to the weights on the surface dataset, and
       ! ran init_interp, this point would keep its cold start initialization
       ! values. Then, in the first time step of the run loop, it's possible that this
       ! point would become active because, according to glc, there is actually > 0% bare
       ! ground in that grid cell. We don't do any state / flux adjustments in the first
       ! time step after init_interp due to glacier area changes, so this vegetated column
       ! would remain at its cold start initialization values, which would be a Bad
       ! Thing. Ensuring that all vegetated points within the icemask are active gets
       ! around this problem - as well as having other benefits, as noted above.)
       if (lun%itype(l) == istsoil) then
          is_active_l = .true.
       end if

       ! Set all lake land units to active
       ! By doing this, lakes are also run virtually in grid cells which will grow
       ! lakes during the transient run.

       if (lun%itype(l) == istdlak) then
            is_active_l = .true.
        end if

    end if

  end function is_active_l

  !-----------------------------------------------------------------------
  logical function is_active_c(c, glc_behavior)
    !
    ! !DESCRIPTION:
    ! Determine whether the given column is active
    !
    ! !USES:
    use landunit_varcon, only : istice, isturb_MIN, isturb_MAX
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: c   ! column index
    type(glc_behavior_type), intent(in) :: glc_behavior
    !
    ! !LOCAL VARIABLES:
    integer :: l  ! landunit index
    integer :: g  ! grid cell index
    !------------------------------------------------------------------------

    if (all_active) then
       is_active_c = .true.

    else
       l =col%landunit(c)
       g =col%gridcell(c)

       is_active_c = .false.

       ! ------------------------------------------------------------------------
       ! General conditions under which is_active_c NEEDS to be true in order to satisfy
       ! the requirements laid out at the top of this module:
       ! ------------------------------------------------------------------------
       if (lun%active(l) .and. col%wtlunit(c) > 0._r8) is_active_c = .true.

       ! ------------------------------------------------------------------------
       ! Conditions under which is_active_c is set to true because we want extra virtual columns:
       ! ------------------------------------------------------------------------

       if (lun%itype(l) == istice .and. &
            glc_behavior%has_virtual_columns_grc(g)) then
          is_active_c = .true.
       end if

       ! We don't really need to run over 0-weight urban columns. But because of some
       ! messiness in the urban code (many loops are over the landunit filter, then drill
       ! down to columns - so we would need to add 'col%active(c)' conditionals in many
       ! places) it keeps the code cleaner to run over 0-weight urban columns. This generally
       ! shouldn't add much computation time, since in most places, all urban columns are
       ! non-zero weight if the landunit is non-zero weight.
       if (lun%active(l) .and. (lun%itype(l) >= isturb_MIN .and. lun%itype(l) <= isturb_MAX)) then
          is_active_c = .true.
       end if
    end if

  end function is_active_c

  !-----------------------------------------------------------------------
  logical function is_active_p(p)
    !
    ! !DESCRIPTION:
    ! Determine whether the given patch is active
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p   ! patch index
    !
    ! !LOCAL VARIABLES:
    integer :: c  ! column index
    !------------------------------------------------------------------------

    if (all_active) then
       is_active_p = .true.

    else
       c =patch%column(p)

       is_active_p = .false.

       ! ------------------------------------------------------------------------
       ! General conditions under which is_active_p NEEDS to be true in order to satisfy
       ! the requirements laid out at the top of this module:
       ! ------------------------------------------------------------------------
       if (col%active(c) .and. patch%wtcol(p) > 0._r8) is_active_p = .true.

    end if

  end function is_active_p

  !-----------------------------------------------------------------------
  function get_landunit_weight(g, ltype) result(weight)
    !
    ! !DESCRIPTION:
    ! Get the subgrid weight of a given landunit type on a single grid cell
    !
    ! !USES:
    use clm_varcon, only : ispval
    !
    ! !ARGUMENTS:
    real(r8) :: weight  ! function result
    integer , intent(in) :: g     ! grid cell index
    integer , intent(in) :: ltype ! landunit type of interest
    !
    ! !LOCAL VARIABLES:
    integer :: l ! landunit index

    character(len=*), parameter :: subname = 'get_landunit_weight'
    !-----------------------------------------------------------------------

    l = grc%landunit_indices(ltype, g)
    if (l == ispval) then
       weight = 0._r8
    else
       weight = lun%wtgcell(l)
    end if

  end function get_landunit_weight

  !-----------------------------------------------------------------------
  subroutine set_landunit_weight(g, ltype, weight)
    !
    ! !DESCRIPTION:
    ! Set the subgrid weight of a given landunit type on a single grid cell
    !
    ! !USES:
    use clm_varcon, only : ispval
    !
    ! !ARGUMENTS:
    integer , intent(in) :: g      ! grid cell index
    integer , intent(in) :: ltype  ! landunit type of interest
    real(r8), intent(in) :: weight ! new weight of this landunit
    !
    ! !LOCAL VARIABLES:
    integer :: l ! landunit index

    character(len=*), parameter :: subname = 'set_landunit_weight'
    !-----------------------------------------------------------------------

    l = grc%landunit_indices(ltype, g)
    if (l /= ispval) then
       lun%wtgcell(l) = weight
    else if (weight > 0._r8) then
       write(iulog,*) subname//' ERROR: Attempt to assign non-zero weight to a non-existent landunit'
       write(iulog,*) 'g, l, ltype, weight = ', g, l, ltype, weight
       call endrun(subgrid_index=l, subgrid_level=subgrid_level_landunit, msg=errMsg(sourcefile, __LINE__))
    end if

  end subroutine set_landunit_weight


  !-----------------------------------------------------------------------
  function is_gcell_all_ltypeX(g, ltype) result(all_ltypeX)
    !
    ! !DESCRIPTION:
    ! Determine if the given grid cell is 100% covered by the landunit type given by ltype
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    logical :: all_ltypeX        ! function result
    integer, intent(in) :: g     ! grid cell index
    integer, intent(in) :: ltype ! landunit type of interest
    !
    ! !LOCAL VARIABLES:
    real(r8) :: wt_lunit ! subgrid weight of the given landunit

    real(r8), parameter :: tolerance = 1.e-13_r8  ! tolerance for checking whether landunit's weight is 1
    character(len=*), parameter :: subname = 'is_gcell_all_ltypeX'
    !------------------------------------------------------------------------------

    wt_lunit = get_landunit_weight(g, ltype)
    if (wt_lunit >= (1._r8 - tolerance)) then
       all_ltypeX = .true.
    else
       all_ltypeX = .false.
    end if

  end function is_gcell_all_ltypeX

  !------------------------------------------------------------------------------
  subroutine check_weights (bounds, active_only)
    !
    ! !DESCRIPTION:
    ! Check subgrid weights.
    !
    ! This routine operates in two different modes, depending on the value of active_only. If
    ! active_only is true, then we check the sum of weights of the ACTIVE children,
    ! grandchildren, etc. of a given point. If active_only is false, then we check the sum of
    ! weights of ALL children, grandchildren, etc. of a given point.
    !
    ! Normally this routine will be called twice: once with active_only=false, and once with
    ! active_only=true.
    !
    ! !USES
    !
    ! !ARGUMENTS
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    logical, intent(in) :: active_only ! true => check sum of weights just of ACTIVE children, grandchildren, etc.
    !
    ! !LOCAL VARIABLES:
    integer :: g,l,c,p     ! loop counters
    real(r8), allocatable :: sumwtcol(:), sumwtlunit(:), sumwtgcell(:)
    logical :: error_found                ! true if we find an error
    character(len=*), parameter :: subname = 'check_weights'
    !------------------------------------------------------------------------------

    allocate(sumwtcol(bounds%begc:bounds%endc))
    allocate(sumwtlunit(bounds%begl:bounds%endl))
    allocate(sumwtgcell(bounds%begg:bounds%endg))

    error_found = .false.

    ! Check patch-level weights
    sumwtcol(bounds%begc : bounds%endc) = 0._r8
    sumwtlunit(bounds%begl : bounds%endl) = 0._r8
    sumwtgcell(bounds%begg : bounds%endg) = 0._r8

    do p = bounds%begp,bounds%endp
       c = patch%column(p)
       l = patch%landunit(p)
       g = patch%gridcell(p)

       if ((active_only .and. patch%active(p)) .or. .not. active_only) then
          sumwtcol(c) = sumwtcol(c) + patch%wtcol(p)
          sumwtlunit(l) = sumwtlunit(l) + patch%wtlunit(p)
          sumwtgcell(g) = sumwtgcell(g) + patch%wtgcell(p)
       end if
    end do

    do c = bounds%begc,bounds%endc
       if (.not. weights_okay(sumwtcol(c), active_only, col%active(c))) then
          write(iulog,*) trim(subname),' ERROR: at c = ',c,'total PFT weight is ',sumwtcol(c), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do

    do l = bounds%begl,bounds%endl
       if (.not. weights_okay(sumwtlunit(l), active_only, lun%active(l))) then
          write(iulog,*) trim(subname),' ERROR: at l = ',l,'total PFT weight is ',sumwtlunit(l), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do

    do g = bounds%begg,bounds%endg
       if (.not. weights_okay(sumwtgcell(g), active_only, i_am_active=.true.)) then
          write(iulog,*) trim(subname),' ERROR: at g = ',g,'total PFT weight is ',sumwtgcell(g), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do

    ! Check col-level weights
    sumwtlunit(bounds%begl : bounds%endl) = 0._r8
    sumwtgcell(bounds%begg : bounds%endg) = 0._r8

    do c = bounds%begc,bounds%endc
       l = col%landunit(c)
       g = col%gridcell(c)

       if ((active_only .and. col%active(c)) .or. .not. active_only) then
          sumwtlunit(l) = sumwtlunit(l) + col%wtlunit(c)
          sumwtgcell(g) = sumwtgcell(g) + col%wtgcell(c)
       end if
    end do

    do l = bounds%begl,bounds%endl
       if (.not. weights_okay(sumwtlunit(l), active_only, lun%active(l))) then
          write(iulog,*) trim(subname),' ERROR: at l = ',l,'total col weight is ',sumwtlunit(l), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do

    do g = bounds%begg,bounds%endg
       if (.not. weights_okay(sumwtgcell(g), active_only, i_am_active=.true.)) then
          write(iulog,*) trim(subname),' ERROR: at g = ',g,'total col weight is ',sumwtgcell(g), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do

    ! Check landunit-level weights
    sumwtgcell(bounds%begg : bounds%endg) = 0._r8

    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       if ((active_only .and. lun%active(l)) .or. .not. active_only) then
          sumwtgcell(g) = sumwtgcell(g) + lun%wtgcell(l)
       end if
    end do

    do g = bounds%begg,bounds%endg
       if (.not. weights_okay(sumwtgcell(g), active_only, i_am_active=.true.)) then
          write(iulog,*) trim(subname),' ERROR: at g = ',g,'total lunit weight is ',sumwtgcell(g), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do

    deallocate(sumwtcol, sumwtlunit, sumwtgcell)

    if (error_found) then
       write(iulog,*) ' '
       write(iulog,*) 'If you are seeing this message at the beginning of a run with'
       write(iulog,*) 'use_init_interp = .true. and init_interp_method = "use_finidat_areas",'
       write(iulog,*) 'and you are seeing weights less than 1, then a likely cause is:'
       write(iulog,*) 'For the above-mentioned grid cell(s):'
       write(iulog,*) 'The matching input grid cell had some non-zero-weight subgrid type'
       write(iulog,*) 'that is not present in memory in the new run.'
       write(iulog,*) ' '
       write(iulog,*) 'If you are using a ctsm5.2 or later fsurdat file containing'
       write(iulog,*) 'PCT_OCEAN > 0, then you may solve the error by setting'
       write(iulog,*) 'convert_ocean_to_land = .true.'
       write(iulog,*) ' '
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! Success

  end subroutine check_weights

  !-----------------------------------------------------------------------
  logical function weights_okay(sumwts, active_weights_only, i_am_active)
    !
    ! !DESCRIPTION:
    ! Determine if sumwts (the sum of weights of children, grandchildren or
    ! great-grandchildren of a column, landunit or grid cell) satisfies the requirements laid
    ! out above.
    !
    ! The way this is determined depends on the values of two other variables:
    ! - active_weights_only: does sumwts just include weights of active children,
    !   grandchildren or great-grandchilden? (alternative is that it includes weights of ALL
    !   children, grandchildren or great-grandchildren)
    ! - i_am_active: true if the column, landunit or grid cell of interest is active
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: sumwts              ! sum of weights of children, grandchildren or great-grandchildren
    logical , intent(in) :: active_weights_only ! true if sumwts just includes active children, etc.
    logical , intent(in) :: i_am_active         ! true if the current point is active
    !
    ! !LOCAL VARIABLES:
    logical :: weights_equal_1
    real(r8), parameter :: tolerance = 1.e-12_r8  ! tolerance for checking whether weights sum to 1
    !------------------------------------------------------------------------

    weights_equal_1 = (abs(sumwts - 1._r8) <= tolerance)

    if (active_weights_only) then
       if (i_am_active) then        ! condition (2) above
          weights_okay = weights_equal_1
       else                         ! condition (3) above
          weights_okay = (sumwts == 0._r8 .or. weights_equal_1)
       end if
    else                            ! condition (1) above
       ! (note that i_am_active is irrelevant in this case)
       weights_okay = weights_equal_1
    end if

  end function weights_okay

  !-----------------------------------------------------------------------
  subroutine set_subgrid_diagnostic_fields(bounds)
    !
    ! !DESCRIPTION:
    ! Set history fields giving diagnostics about subgrid weights
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'set_subgrid_diagnostic_fields'
    !-----------------------------------------------------------------------

    call set_pct_landunit_diagnostics(bounds)

    ! Note: (MV, 10-17-14): The following has an use_fates if-block around it since
    ! the pct_pft_diagnostics referens to patch%itype(p) which is not used by ED
    ! Note: (SPM, 10-20-15): If this isn't set then debug mode with intel and
    ! yellowstone will fail when trying to write pct_nat_pft since it contains
    ! all NaN's.
    call set_pct_pft_diagnostics(bounds)

    call set_pct_glc_mec_diagnostics(bounds)

  end subroutine set_subgrid_diagnostic_fields

  !-----------------------------------------------------------------------
  subroutine set_pct_landunit_diagnostics(bounds)
    !
    ! !DESCRIPTION:
    ! Set pct_landunit diagnostic field: % of each landunit on the grid cell
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g, l  ! grid cell & landunit indices
    integer :: ltype ! landunit type

    character(len=*), parameter :: subname = 'set_pct_landunit_diagnostics'
    !-----------------------------------------------------------------------

    subgrid_weights_diagnostics%pct_landunit(bounds%begg:bounds%endg, :) = 0._r8

    do l = bounds%begl, bounds%endl
       g = lun%gridcell(l)
       ltype = lun%itype(l)
       subgrid_weights_diagnostics%pct_landunit(g, ltype) = lun%wtgcell(l) * 100._r8
    end do

  end subroutine set_pct_landunit_diagnostics

  !-----------------------------------------------------------------------
  subroutine set_pct_glc_mec_diagnostics(bounds)
    !
    ! !DESCRIPTION:
    ! Set pct_glc_mec diagnostic field: % of each glc_mec column on the glc landunit
    !
    ! Note that pct_glc_mec will be 0 for all elevation classes in a grid cell that does
    ! not have a glc landunit. However, it will still sum to 100% for a grid cell
    ! that has a 0-weight (i.e., virtual) glc landunit.
    !
    ! !USES:
    use landunit_varcon, only : istice
    use column_varcon, only : col_itype_to_ice_class
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: c,l,g          ! indices
    integer :: ice_class      ! ice class (1..maxpatch_glc)

    character(len=*), parameter :: subname = 'set_pct_glc_mec_diagnostics'
    !-----------------------------------------------------------------------

    subgrid_weights_diagnostics%pct_glc_mec(bounds%begg:bounds%endg, :) = 0._r8

    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
       l = col%landunit(c)
       if (lun%itype(l) == istice) then
          ice_class = col_itype_to_ice_class(col%itype(c))
          subgrid_weights_diagnostics%pct_glc_mec(g, ice_class) = col%wtlunit(c) * 100._r8
       end if
    end do

  end subroutine set_pct_glc_mec_diagnostics

  !-----------------------------------------------------------------------
  subroutine set_pct_pft_diagnostics(bounds)
    !
    ! !DESCRIPTION:
    ! Set pct_nat_pft & pct_cft diagnostic fields: % of PFTs on their landunit
    !
    ! !USES:
    use landunit_varcon, only : istsoil, istcrop
    use clm_varpar, only : natpft_lb, cft_lb
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p,l,g           ! indices
    integer :: ptype           ! patch itype
    integer :: ptype_1indexing ! patch itype, translated into 1-indexing for the given landunit type

    character(len=*), parameter :: subname = 'set_pct_pft_diagnostics'
    !-----------------------------------------------------------------------

    subgrid_weights_diagnostics%pct_nat_pft(bounds%begg:bounds%endg, :) = 0._r8

    ! Note that pct_cft will be 0-size if cft_size is 0 (which can happen if we don't
    ! have a crop landunit). But it doesn't hurt to have this line setting all elements
    ! to 0, and doing this always allows us to avoid extra logic which could be a
    ! maintenance problem.
    subgrid_weights_diagnostics%pct_cft(bounds%begg:bounds%endg, :) = 0._r8

    do p = bounds%begp,bounds%endp
       g = patch%gridcell(p)
       l = patch%landunit(p)
       ptype = patch%itype(p)
       if (lun%itype(l) == istsoil .and. (.not.use_fates.or.use_fates_sp) ) then
          ptype_1indexing = ptype + (1 - natpft_lb)
          subgrid_weights_diagnostics%pct_nat_pft(g, ptype_1indexing) = patch%wtlunit(p) * 100._r8
       else if (lun%itype(l) == istcrop) then
          ptype_1indexing = ptype + (1 - cft_lb)
          subgrid_weights_diagnostics%pct_cft(g, ptype_1indexing) = patch%wtlunit(p) * 100._r8
       end if
    end do

  end subroutine set_pct_pft_diagnostics

end module subgridWeightsMod
