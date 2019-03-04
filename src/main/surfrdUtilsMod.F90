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
  public :: collapse_to_dominant ! Collapse to dominant pfts or landunits
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
    use clm_varpar      , only : cft_size
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
  subroutine collapse_to_dominant(weight, lower_bound, upper_bound, begg, endg, n_dominant)
    !
    ! DESCRIPTION
    ! Collapse to the top N dominant pfts or landunits (n_dominant)
    !
    ! !USES:
    use array_utils, only: find_k_max_indices
    !
    ! !ARGUMENTS:
    ! Use begg and endg rather than 'bounds', because bounds may not be
    ! available yet when this is called
    integer, intent(in) :: begg  ! Beginning grid cell index
    integer, intent(in) :: endg  ! Ending grid cell index
    integer, intent(in) :: lower_bound  ! lower bound of pft or landunit indices
    integer, intent(in) :: upper_bound  ! upper bound of pft or landunit indices
    integer, intent(in) :: n_dominant  ! # dominant pfts or landunits
    ! This array modified in-place
    ! Weights of pfts or landunits per grid cell
    ! Dimensioned [g, lower_bound:upper_bound]
    real(r8), intent(inout) :: weight(begg:endg, lower_bound:upper_bound)
    !
    ! !LOCAL VARIABLES:
    integer :: g  ! gridcell index
    integer :: m  ! pft or landunit index
    integer :: n  ! index of the order of the dominant pfts or landunits
    integer, allocatable :: max_indices(:)  ! array of dominant pft or landunit index values
    real(r8) :: wt_dom_sum

    character(len=*), parameter :: subname = 'collapse_to_dominant'

    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(weight) == (/endg, upper_bound/)), errMsg(sourcefile, __LINE__))

    ! Find the top N dominant pfts or landunits to collapse the data to
    ! n_dominant < 0 is not allowed (error check in controlMod.F90)
    ! Default value n_dominant = 0 or a user-selected n_dominant = upper_bound
    ! means "do not collapse pfts" and skip over this subroutine's work
    if (n_dominant > 0 .and. n_dominant < upper_bound) then
       allocate(max_indices(n_dominant))
       do g = begg, endg
          max_indices = 0  ! initialize
          call find_k_max_indices(weight(g,:), lower_bound, n_dominant, &
                                  max_indices)

          ! Adjust weight by normalizing the dominant weights to 1
          ! (currently they sum to <= 1) and setting the remaining weights to 0.
          wt_dom_sum = 0._r8  ! initialize the dominant pft or landunit sum
          do n = 1, n_dominant
             m = max_indices(n)
             wt_dom_sum = weight(g,m) + wt_dom_sum
          end do
          ! Normalize dominant pft or landunit weights to 1; if non-existent,
          ! set the weights to 0.
          if (wt_dom_sum <= 0._r8) then
             call endrun(msg = subname//' wt_dom_sum should never be <= 0'//&
                  ' but it is here' // errMsg(sourcefile, __LINE__))
          else
             do n = 1, n_dominant
                m = max_indices(n)
                weight(g,m) = weight(g,m) / wt_dom_sum
             end do
          end if
          ! Set non-dominant weights to 0
          do m = lower_bound, upper_bound
             if (.not. any(max_indices == m)) then
                weight(g,m) = 0._r8
             end if
          end do

       end do

       ! Error checks
       call check_sums_equal_1(weight, begg, 'weight', subname)

       deallocate(max_indices)
    end if

  end subroutine collapse_to_dominant

  !-----------------------------------------------------------------------
  subroutine collapse_crop_var(crop_var, cft_size, begg, endg)
    !
    ! DESCRIPTION
    ! After collapse_crop_types, ensure that
    ! crop-related variables are consistent with the new crop weights (wt_cft).
    !
    ! List of crop-related variables (locally named crop_var):
    ! -
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
