module surfrdUtilsMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains utility methods that can be used when reading surface datasets or similar
  ! datasets (such as the landuse_timeseries dataset)
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod , only : r8 => shr_kind_r8
  use clm_varctl   , only : iulog
  use abortutils   , only : endrun
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use spmdMod      , only : masterproc
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: check_sums_equal_1  ! Confirm that sum(arr(n,:)) == 1 for all n
  public :: renormalize         ! Renormalize an array
  public :: convert_cft_to_pft  ! Conversion of crop CFT to natural veg PFT:w
  public :: collapse_crop_types ! Collapse unused crop types into types used in this run
  public :: collapse_all_pfts  ! Collapse to dominant plant types
  public :: collapse_crop_var  ! Collapse crop variables according to cft weights determined in previous "collapse" subroutines

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !-----------------------------------------------------------------------

contains
  
  !-----------------------------------------------------------------------
  subroutine check_sums_equal_1(arr, lb, name, caller, ier, sumto)
    !
    ! !DESCRIPTION:
    ! Confirm that sum(arr(n,:)) == 1 for all n. If this isn't true for any n, abort with a message.
    !
    ! !ARGUMENTS:
    integer         , intent(in) :: lb           ! lower bound of the first dimension of arr
    real(r8)        , intent(in) :: arr(lb:,:)   ! array to check
    character(len=*), intent(in) :: name         ! name of array
    character(len=*), intent(in) :: caller       ! identifier of caller, for more meaningful error messages
    integer, optional, intent(out):: ier         ! Return an error code rather than abort
    real(r8), optional, intent(out):: sumto      ! The value the array should sum to (1.0 if not provided)
    !
    ! !LOCAL VARIABLES:
    logical :: found
    integer :: nl
    integer :: nindx
    real(r8), parameter :: eps = 1.e-13_r8
    real(r8) :: TotalSum
    !-----------------------------------------------------------------------

    TotalSum = 1._r8
    if ( present(sumto) ) TotalSum = sumto
    if( present(ier) ) ier = 0
    found = .false.

    do nl = lbound(arr, 1), ubound(arr, 1)
       if (abs(sum(arr(nl,:)) - TotalSum) > eps) then
          found = .true.
          nindx = nl
          exit
       end if
    end do

    if (found) then
       write(iulog,*) trim(caller), ' ERROR: sum of ', trim(name), ' not ', TotalSum, ' at nl=', nindx
       write(iulog,*) 'sum is: ', sum(arr(nindx,:))
       if( present(ier) ) then
          ier = -10
       else
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

  end subroutine check_sums_equal_1

  !-----------------------------------------------------------------------
  subroutine renormalize(arr, lb, normal)
    !
    ! !DESCRIPTION:
    ! Re normalize an array so that it sums to the input value
    !
    ! !ARGUMENTS:
    integer         , intent(in) :: lb            ! lower bound of the first dimension of arr
    real(r8)        , intent(inout) :: arr(lb:,:) ! array to check
    real(r8)        , intent(in) :: normal        ! normal to sum to
    !
    ! !LOCAL VARIABLES:
    integer :: nl        ! Array index
    real(r8) :: arr_sum  ! sum of array
    real(r8) :: ratio    ! ratio to multiply by
    !-----------------------------------------------------------------------

    do nl = lbound(arr, 1), ubound(arr, 1)
       arr_sum = sum(arr(nl,:))
       if ( arr_sum /= 0.0_r8 )then
          ratio     = normal / arr_sum
          arr(nl,:) = arr(nl,:) * ratio
       end if
    end do

  end subroutine renormalize

!-----------------------------------------------------------------------
  subroutine convert_cft_to_pft( begg, endg, cftsize, wt_cft )
    !
    ! !DESCRIPTION:
    !        Convert generic crop types that were read in as seperate CFT's on
    !        a crop landunit, and put them on the vegetated landunit.
    ! !USES:
    use clm_instur      , only : wt_lunit, wt_nat_patch
    use clm_varpar      , only : cft_size, natpft_size
    use pftconMod       , only : nc3crop
    use landunit_varcon , only : istsoil, istcrop
    ! !ARGUMENTS:
    implicit none
    integer          , intent(in)    :: begg, endg
    integer          , intent(in)    :: cftsize          ! CFT size
    real(r8)         , intent(inout) :: wt_cft(begg:,:)  ! CFT weights
    !
    ! !LOCAL VARIABLES:
    integer :: g    ! index
!-----------------------------------------------------------------------
    SHR_ASSERT_ALL((ubound(wt_cft) == (/endg, cftsize/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(wt_nat_patch) == (/endg, nc3crop+cftsize-1/)), errMsg(sourcefile, __LINE__))

    do g = begg, endg
       if ( wt_lunit(g,istcrop) > 0.0_r8 )then
          ! Move CFT over to PFT and do weighted average of the crop and soil parts
          wt_nat_patch(g,:) = wt_nat_patch(g,:) * wt_lunit(g,istsoil)
          wt_cft(g,:)       = wt_cft(g,:) * wt_lunit(g,istcrop)
          wt_nat_patch(g,nc3crop:) = wt_cft(g,:)      ! Add crop CFT's to end of natural veg PFT's
          wt_lunit(g,istsoil) = (wt_lunit(g,istsoil) + wt_lunit(g,istcrop)) ! Add crop landunit to soil landunit
          wt_nat_patch(g,:)   =  wt_nat_patch(g,:) / wt_lunit(g,istsoil)
          wt_lunit(g,istcrop) = 0.0_r8                ! Zero out crop CFT's
       else
          wt_nat_patch(g,nc3crop:) = 0.0_r8    ! Make sure generic crops are zeroed out
       end if
    end do

  end subroutine convert_cft_to_pft

  !-----------------------------------------------------------------------
  subroutine collapse_all_pfts(wt_lunit, wt_nat_patch, natpft_size, wt_cft, cft_size, begg, endg, n_dom_soil_patches)
    !
    ! DESCRIPTION
    ! Collapse to the top N dominant soil patches (n_dom_soil_patches) among all
    ! present pfts, cfts, & bare ground.
    ! - Bare ground could be up to 1 patch before collapsing.
    ! - Pfts could be up to 14 before collapsing.
    ! - Cfts could be up to 2 or 78 with use_crop = .F. or .T., respectively.
    !
    ! !USES:
    use clm_varpar, only: natpft_lb, natpft_ub, cft_lb, cft_ub, maxveg
    use landunit_varcon, only: istsoil, istcrop, max_lunit
    use pftconMod, only: noveg, nc3crop, pftcon
    use array_utils, only: find_k_max_indices
    !
    ! !ARGUMENTS:
    ! Use begg and endg rather than 'bounds', because bounds may not be
    ! available yet when this is called
    integer, intent(in) :: begg  ! Beginning grid cell index
    integer, intent(in) :: endg  ! Ending grid cell index
    integer, intent(in) :: natpft_size  ! CFT size
    integer, intent(in) :: cft_size  ! CFT size
    integer, intent(in) :: n_dom_soil_patches  ! # dominant soil patches
    ! These arrays modified in-place
    ! Weights of landunits
    real(r8), intent(inout) :: wt_lunit(begg:,:)
    ! Weights of managed (cft) and unmanaged (nat_patch) soil patches per
    ! grid cell
    ! Dimensioned [g, natpft_lb:natpft_ub]
    real(r8), intent(inout) :: wt_nat_patch(begg:, natpft_lb:)
    ! Dimensioned [g, cft_lb:cft_lb+cft_size-1]
    real(r8), intent(inout) :: wt_cft(begg:, cft_lb:)
    !
    ! !LOCAL VARIABLES:
    integer :: g  ! gridcell index
    integer :: m  ! soil patch index
    integer :: n  ! index of the order of the dominant soil patches
    integer, allocatable :: max_indices(:)  ! array of dominant soil patch index values
    real(r8) :: wt_all_patch(natpft_lb:cft_ub)
    real(r8) :: wt_dom_nat_sum
    real(r8) :: wt_dom_cft_sum
    character(len=*), parameter :: subname = 'collapse_all_pfts'

    !-----------------------------------------------------------------------

    if (cft_size > 0) then  ! The opposite applies only if use_fates

       SHR_ASSERT_ALL((ubound(wt_lunit) == (/endg, max_lunit/)), errMsg(sourcefile, __LINE__))
       SHR_ASSERT_ALL((ubound(wt_nat_patch) == (/endg, natpft_lb+natpft_size-1/)), errMsg(sourcefile, __LINE__))
       SHR_ASSERT_ALL((ubound(wt_cft) == (/endg, cft_lb+cft_size-1/)), errMsg(sourcefile, __LINE__))
       SHR_ASSERT_ALL((ubound(wt_all_patch) == (/natpft_lb+natpft_size+cft_size-1/)), errMsg(sourcefile, __LINE__))

       ! Find the top N dominant soil patches to collapse pft data to
       ! n_dom_soil_patches < 0 is not allowed (error check in controlMod.F90)
       ! n_dom_soil_patches == 0 or 78 currently mean "do not collapse pfts"
       ! so skip over this subroutine's work
       if (n_dom_soil_patches > 0 .and. n_dom_soil_patches < maxveg) then
          allocate(max_indices(n_dom_soil_patches))
          do g = begg, endg
             ! Normalize nat and cft weights by their landunit weights in
             ! the gridcell and concatenate into a single array before calling
             ! find_k_max_indices
             wt_all_patch(:natpft_ub) = wt_nat_patch(g,:) * &
                                        wt_lunit(g,istsoil)
             wt_all_patch(cft_lb:) = wt_cft(g,:) * &
                                     wt_lunit(g,istcrop)
             max_indices = 0._r8  ! initialize
             call find_k_max_indices(wt_all_patch(:), &
                                     natpft_lb, &
                                     n_dom_soil_patches, &
                                     max_indices)

             ! Adjust wt_nat_patch, wt_cft, and their respective wt_lunit by
             ! normalizing the dominant weights to 1 (currently they sum to <= 1)
             ! and setting the remaining weights to 0.
             ! TODO Possible to use existing function renormalize?
             wt_dom_nat_sum = 0._r8  ! initialize the dominant pft sum
             wt_dom_cft_sum = 0._r8  ! initialize the dominant cft sum
             do n = 1, n_dom_soil_patches
                m = max_indices(n)
                if (m <= natpft_ub) then  ! calculate the dominant pft sum
                   wt_dom_nat_sum = wt_nat_patch(g,m) + wt_dom_nat_sum
                else if (m >= cft_lb) then  ! calculate the dominant cft sum
                   wt_dom_cft_sum = wt_cft(g,m) + wt_dom_cft_sum
                else  ! error message
                   call endrun(subname//': m must range [natpft_lb, cft_ub]')
                end if
             end do
             ! Normalize dominant pft weights to 1; however,
             ! if non-existent, then merge the nat and crop landunit weights
             ! and set the pft weights and pft landunit weight to 0.
             ! Don't bother merging if they're not both > 0.
             if (wt_dom_nat_sum <= 0._r8) then
                if (wt_lunit(g,istcrop) > 0._r8 .and. wt_lunit(g,istsoil) > 0._r8) then
                   wt_lunit(g,istcrop) = wt_lunit(g,istcrop) + wt_lunit(g,istsoil)
                end if
                wt_lunit(g,istsoil) = 0._r8
                wt_nat_patch(g,:) = 0._r8
                wt_nat_patch(g,noveg) = 1._r8  ! to pass error check below
             else
                do n = 1, n_dom_soil_patches
                   m = max_indices(n)
                   if (m <= natpft_ub) then
                      wt_nat_patch(g,m) = wt_nat_patch(g,m) / wt_dom_nat_sum
                   end if
                end do
             end if
             ! Normalize dominant cft weights to 1; however,
             ! if non-existent, then merge the nat and crop landunit weights
             ! and set the cft weights and cft landunit weight to 0.
             ! Don't bother merging if they're not both > 0.
             if (wt_dom_cft_sum <= 0._r8) then
                if (wt_lunit(g,istcrop) > 0._r8 .and. wt_lunit(g,istsoil) > 0._r8) then
                   wt_lunit(g,istsoil) = wt_lunit(g,istsoil) + wt_lunit(g,istcrop)
                end if
                wt_lunit(g,istcrop) = 0._r8
                wt_cft(g,:) = 0._r8
                wt_cft(g,nc3crop) = 1._r8  ! to pass error check below
             else
                do n = 1, n_dom_soil_patches
                   m = max_indices(n)
                   if (m >= cft_lb) then
                      wt_cft(g,m) = wt_cft(g,m) / wt_dom_cft_sum
                   end if
                end do
             end if
             ! Set non-dominant pft and cft weights to 0 for cases where dominant
             ! pfts or cfts were found
             do m = 0, maxveg
                if (.not. any(max_indices == m)) then
                   if (m <= natpft_ub .and. wt_dom_nat_sum > 0._r8) then
                      wt_nat_patch(g,m) = 0._r8
                   else if (m >= cft_lb .and. wt_dom_cft_sum > 0._r8) then
                      wt_cft(g,m) = 0._r8
                   end if
                   ! Set to .F. regardless of whether dominant pfts and cfts were
                   ! found
                   pftcon%is_pft_known_to_model(m) = .false.
                end if
             end do

          end do

          ! Error checks
          call check_sums_equal_1(wt_lunit, begg, 'wt_lunit', subname)
          call check_sums_equal_1(wt_cft, begg, 'wt_cft', subname)
          call check_sums_equal_1(wt_nat_patch, begg, 'wt_nat_patch', subname)

          deallocate(max_indices)
       end if
    end if

  end subroutine collapse_all_pfts

  !-----------------------------------------------------------------------
  subroutine collapse_crop_var(crop_var, cft_size, begg, endg)
    !
    ! DESCRIPTION
    ! After collapsing pft & cft weights in collapse_all_pfts, now ensure that
    ! crop-related variables are consistent with the new crop weights (wt_cft).
    !
    ! Crop-related variables (locally named crop_var in this subroutine):
    !
    ! !USES:
    use clm_varpar, only: cft_lb, cft_ub
    use pftconMod, only: pftcon
    !
    ! !ARGUMENTS:
    ! Use begg and endg rather than 'bounds', because bounds may not be
    ! available yet when this is called
    integer, intent(in) :: begg  ! Beginning grid cell index
    integer, intent(in) :: endg  ! Ending grid cell index
    integer, intent(in) :: cft_size  ! CFT size

    ! Crop variable dimensioned [g, cft_lb:cft_lb+cft_size-1] modified in-place
    real(r8), intent(inout) :: crop_var(begg:, cft_lb:)

    ! !LOCAL VARIABLES:
    integer :: g  ! gridcell index
    integer :: m  ! cft index

    character(len=*), parameter :: subname = 'collapse_crop_var'
    !-----------------------------------------------------------------------

    if (cft_size > 0) then  ! The opposite applies only if use_fates

       SHR_ASSERT_ALL((ubound(crop_var) == (/endg, cft_lb+cft_size-1/)), errMsg(sourcefile, __LINE__))

       do g = begg, endg
          do m = cft_lb, cft_ub
             if (.not. pftcon%is_pft_known_to_model(m)) then
                crop_var(g,m) = 0._r8
             end if
          end do
       end do

    end if

  end subroutine collapse_crop_var

  !-----------------------------------------------------------------------
  subroutine collapse_crop_types(wt_cft, fert_cft, cftsize, begg, endg, verbose, sumto)
    !
    ! !DESCRIPTION:
    ! Collapse unused crop types into types used in this run.
    !
    ! !USES:
    use clm_varctl , only : irrigate, use_crop
    use clm_varpar , only : cft_lb, cft_ub, maxveg
    use pftconMod  , only : nc3crop, nc3irrig, pftcon
    !
    ! !ARGUMENTS:

    ! Note that we use begg and endg rather than 'bounds', because bounds may not be
    ! available yet when this is called
    integer, intent(in) :: begg     ! Beginning grid cell index
    integer, intent(in) :: endg     ! Ending grid cell index
    integer, intent(in) :: cftsize  ! CFT size

    ! Weight and fertilizer of each CFT in each grid cell; dimensioned [g, cft_lb:cft_lb+cftsize-1]
    ! This array is modified in-place
    real(r8), intent(inout) :: wt_cft(begg:, cft_lb:)
    real(r8), intent(inout) :: fert_cft(begg:, cft_lb:)
    real(r8), intent(in), optional :: sumto                   ! What weights should sum to for one grid-cell

    logical, intent(in) :: verbose  ! If true, print some extra information
    !
    ! !LOCAL VARIABLES:
    integer :: g
    integer :: m
    real(r8) :: wt_cft_to
    real(r8) :: wt_cft_from
    real(r8) :: wt_cft_merge
    real(r8) :: TotalSum       ! What the total is expected to sum to

    character(len=*), parameter :: subname = 'collapse_crop_types'
    !-----------------------------------------------------------------------

    if (cftsize > 0) then  ! The opposite applies only if use_fates

       SHR_ASSERT_ALL((ubound(wt_cft)   == (/endg, cft_lb+cftsize-1/)), errMsg(sourcefile, __LINE__))
       SHR_ASSERT_ALL((ubound(fert_cft) == (/endg, cft_lb+cftsize-1/)), errMsg(sourcefile, __LINE__))

       TotalSum = 1.0_r8
       if ( present(sumto) ) TotalSum = sumto  ! e.g. sumto may = 100._r8

       ! -----------------------------------------------------------------------
       ! If not using irrigation, merge irrigated CFTs into rainfed CFTs
       ! -----------------------------------------------------------------------

       if (.not. irrigate) then
          if (verbose .and. masterproc) then
             write(iulog,*) trim(subname)//' irrigate=.F., so merging irrigated pfts with rainfed'
          end if

          do g = begg, endg
             ! Left Hand Side: merged rainfed+irrigated crop pfts from nc3crop
             !                 to maxveg-1, stride 2
             ! Right Hand Side: rainfed crop pfts from nc3crop to maxveg-1,
             !                  stride 2
             ! plus             irrigated crop pfts from nc3irrig to maxveg,
             !                  stride 2
             ! where stride 2 means "every other"
             wt_cft(g, nc3crop:maxveg-1:2) = &
                  wt_cft(g, nc3crop:maxveg-1:2) + wt_cft(g, nc3irrig:maxveg:2)
             wt_cft(g, nc3irrig:maxveg:2)  = 0._r8
          end do

          call check_sums_equal_1(wt_cft, begg, 'wt_cft', subname//': irrigation', sumto=TotalSum)
       end if

       ! -----------------------------------------------------------------------
       ! Merge CFTs into the list of crops that CLM knows how to model
       ! -----------------------------------------------------------------------

       if (verbose .and. masterproc .and. use_crop) then
          write(iulog, *) trim(subname) // ' merging wheat, barley, and rye into temperate cereals'
          write(iulog, *) trim(subname) // ' clm knows how to model corn, temperate cereals, and soybean'
          write(iulog, *) trim(subname) // ' all other crops are lumped with the generic crop pft'
       else if (verbose .and. masterproc .and. .not. use_crop) then
          write(iulog, *) trim(subname) // ' merging crops into C3 generic crops'
       end if

       do g = begg, endg
          do m = 1, maxveg
             if (m /= pftcon%mergetoclmpft(m)) then
                wt_cft_to = wt_cft(g, pftcon%mergetoclmpft(m))
                wt_cft_from = wt_cft(g, m)
                wt_cft_merge = wt_cft_to + wt_cft_from
                wt_cft(g, pftcon%mergetoclmpft(m)) = wt_cft_merge
                wt_cft(g, m) = 0._r8
                if (wt_cft_merge > 0._r8) then
                   fert_cft(g,pftcon%mergetoclmpft(m)) = (wt_cft_to * fert_cft(g,pftcon%mergetoclmpft(m)) + &
                                                         wt_cft_from * fert_cft(g,m)) / wt_cft_merge
                end if
                pftcon%is_pft_known_to_model(m) = .false.
             end if
          end do

       end do

       call check_sums_equal_1(wt_cft, begg, 'wt_cft', subname//': mergetoclmpft', sumto=TotalSum)
       if ( .not. use_crop )then
          if ( any(wt_cft(begg:endg,cft_ub+1:) /= 0.0_r8) )then
             call endrun(msg = subname//' without prognostic crops (use_crop=F) and weight of CFT of prognostic crop'//&
                  ' is not zero as expected' // errMsg(sourcefile, __LINE__))
          end if
          if ( any(fert_cft(begg:endg,cft_ub+1:) /= 0.0_r8) )then
             call endrun(msg = subname//' without prognostic crops (use_crop=F) and fertilizer of prognostic crop'// &
                  ' is not zero as expected' // errMsg(sourcefile, __LINE__))
          end if
       end if

    end if

  end subroutine collapse_crop_types


end module surfrdUtilsMod
