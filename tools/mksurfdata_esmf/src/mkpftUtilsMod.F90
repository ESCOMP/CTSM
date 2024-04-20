module mkpftUtilsMod

  !-----------------------------------------------------------------------
  ! Lower-level utilities used in making PFT data.
  !
  ! These are separated out from mkpftMod mainly as an aid to testing.
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_abort

  implicit none
  private

  public :: convert_from_p2g      ! Convert a p2g array into pct_pft_type objects

  private :: get_default_natpft   ! Get the default natural pft breakdown, for a 0-area natural veg. landunit
  private :: get_default_cft      ! Get the default cft breakdown, for a 0-area crop landunit

  interface convert_from_p2g
     module procedure convert_from_p2g_default
     module procedure convert_from_p2g_missing_crops
  end interface convert_from_p2g

!===============================================================
contains
!===============================================================

  subroutine convert_from_p2g_default(pct_p2g, pctnatpft, pctcft)
    !
    ! !DESCRIPTION:
    ! Given the % of each pft on the grid cell, create pct_pft_type objects that give % of
    ! each pft on the landunit and % of each landunit on the grid cell.
    !
    ! !USES:
    use mkpctPftTypeMod  , only : pct_pft_type
    use mkpftConstantsMod, only : natpft_lb, natpft_ub, num_cft, cft_lb, cft_ub
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: pct_p2g(natpft_lb:) ! % of each pft on the grid cell (includes crops as well as natural veg types)
    type(pct_pft_type), intent(out) :: pctnatpft ! natural PFT cover
    type(pct_pft_type), intent(out) :: pctcft    ! crop (CFT) COVER
    !
    ! !LOCAL VARIABLES:
    real(r8), allocatable :: default_natpft(:) ! default p2l for natural PFTs, for grid cells where the current size of the natural veg landunit is 0
    real(r8), allocatable :: default_cft(:)    ! default p2l for CFTs, for grid cells where the current size of the crop landunit is 0

    character(len=*), parameter :: subname = 'convert_from_p2g_default'
    !-----------------------------------------------------------------------

    if (ubound(pct_p2g, 1) /= cft_ub) then
       write(6,*) subname, ' ERROR: upper bound of pct_p2g should be cft_ub'
       write(6,*) 'ubound(pct_p2g), cft_ub = ', ubound(pct_p2g), cft_ub
       call shr_sys_abort()
    end if

    allocate(default_natpft(natpft_lb:natpft_ub))
    default_natpft = get_default_natpft()
    pctnatpft = pct_pft_type(pct_p2g(natpft_lb:natpft_ub), natpft_lb, default_natpft)
    deallocate(default_natpft)

    if (num_cft > 0) then
       allocate(default_cft(cft_lb:cft_ub))
       default_cft = get_default_cft()
       pctcft = pct_pft_type(pct_p2g(cft_lb:cft_ub), cft_lb, default_cft)
       deallocate(default_cft)
    else
       ! create an empty placeholder, with 0 area on the grid cell
       pctcft = pct_pft_type()
    end if

  end subroutine convert_from_p2g_default

  !-----------------------------------------------------------------------
  subroutine convert_from_p2g_missing_crops(pct_p2g, pctcft_saved, pctnatpft, pctcft)
    !
    ! !DESCRIPTION:
    ! Given the % of each pft on the grid cell, create pct_pft_type objects that give %
    ! of each pft on the landunit and % of each landunit on the grid cell.
    !
    ! This version of the routine assumes that pct_p2g only includes the standard PFTs -
    ! not prognostic crops. It takes the relative crop cover from pctcft_saved, and uses
    ! the % cover of the generic c3 crop in pct_p2g to specify the total crop landunit
    ! area.
    !
    ! Typically, pct_p2g will have an upper bound of numstdpft; however, this is not
    ! assumed. Any upper bound is fine as long as the upper bound is greater than
    ! natpft_ub and includes c3cropindex.
    !
    ! Assumptions:
    ! - We are running with prognostic crops (i.e., NOT an empty crop landunit - although
    !   it's fine for the crop landunit area to be 0%)
    ! - In pct_p2g, the only non-zero areas should be:
    !   - Areas of PFTs on the natural veg landunit
    !   - The area of the generic c3 crop
    !
    ! !USES:
    use mkpctPftTypeMod   , only : pct_pft_type
    use mkpftConstantsMod , only : c3cropindex, natpft_lb, natpft_ub, num_cft
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: pct_p2g(natpft_lb:) ! % of each pft on the grid cell (includes crops as well as natural veg types)
    type(pct_pft_type), intent(in)  :: pctcft_saved ! saved crop cover information, used to specify the relative cover of each crop
    type(pct_pft_type), intent(out) :: pctnatpft ! natural PFT cover
    type(pct_pft_type), intent(out) :: pctcft    ! crop (CFT) COVER
    !
    ! !LOCAL VARIABLES:
    real(r8), allocatable :: default_natpft(:) ! default p2l for natural PFTs, for grid cells where the current size of the natural veg landunit is 0
    integer  :: pft_index
    real(r8) :: crop_area  ! area of the crop landunit on the grid cell

    character(len=*), parameter :: subname = 'convert_from_p2g_missing_crops'
    !-----------------------------------------------------------------------
    
    ! Error checking on inputs

    if (num_cft == 0) then
       write(6,*) subname, ' ERROR: this routine should only be called when running with prognostic crops'
       write(6,*) '(i.e., with num_cft > 0)'
       call shr_sys_abort()
    end if

    do pft_index = natpft_ub + 1, ubound(pct_p2g, 1)
       if (pft_index /= c3cropindex .and. pct_p2g(pft_index) > 0._r8) then
          write(6,*) subname, ' ERROR: in pct_p2g, the only non-zero areas should be:'
          write(6,*) '  - areas of PFTs on the natural veg landunit'
          write(6,*) '  - the area of the generic c3 crop'
          write(6,*) '(we do not currently handle the case where the transient input dataset'
          write(6,*) 'has non-zero areas for both pft 15 and pft 16)'
          write(6,*) 'pft_index, area = ', pft_index, pct_p2g(pft_index)
          call shr_sys_abort()
       end if
    end do

    ! Done error checking on inputs

    allocate(default_natpft(natpft_lb:natpft_ub))
    default_natpft = get_default_natpft()
    pctnatpft = pct_pft_type(pct_p2g(natpft_lb:natpft_ub), natpft_lb, default_natpft)
    deallocate(default_natpft)

    pctcft = pctcft_saved
    crop_area = pct_p2g(c3cropindex)
    call pctcft%set_pct_l2g(crop_area)

  end subroutine convert_from_p2g_missing_crops

  !-----------------------------------------------------------------------
  function get_default_natpft() result(default_natpft)
    !
    ! !DESCRIPTION:
    ! Get the default natural pft breakdown, for a 0-area natural veg. landunit.
    !
    ! Currently we use the same default everywhere. In the future, we could change this
    ! to compute default_natpft based on some function of location (e.g., different
    ! values for high latitudes than low latitudes, etc.).
    !
    ! !USES:
    use mkpftConstantsMod, only : baregroundindex, natpft_lb, natpft_ub
    !
    ! !ARGUMENTS:
    real(r8), allocatable :: default_natpft(:)  ! function result
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'get_default_natpft'
    !-----------------------------------------------------------------------
    
    allocate(default_natpft(natpft_lb:natpft_ub))
    default_natpft(:) = 0._r8
    default_natpft(baregroundindex) = 100._r8

  end function get_default_natpft

  !-----------------------------------------------------------------------
  function get_default_cft() result(default_cft)
    !
    ! !DESCRIPTION:
    ! Get the default cft breakdown, for a 0-area crop landunit.
    !
    ! !USES:
    use mkpftConstantsMod, only : c3cropindex, cft_lb, cft_ub
    !
    ! !ARGUMENTS:
    real(r8), allocatable :: default_cft(:)  ! function result
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'get_default_cft'
    !-----------------------------------------------------------------------
    
    allocate(default_cft(cft_lb:cft_ub))
    default_cft(:) = 0._r8
    default_cft(c3cropindex) = 100._r8

  end function get_default_cft


end module mkpftUtilsMod


