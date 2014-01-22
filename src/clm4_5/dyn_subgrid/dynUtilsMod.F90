module dynUtilsMod

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Miscellaneous lower-level utilities needed by the dyn subgrid code
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dyn_time_weights   ! gives the weight of time1 in the current annual interval

contains

  !-----------------------------------------------------------------------
  subroutine dyn_time_weights(wt1, offset)
    !
    ! !DESCRIPTION:
    ! Assumes we have data at annual intervals, provided on Jan. 1.
    ! 
    ! Returns the weight of time1 in the current annual interval, as determined by the
    ! model's current time (with optional offset). The weight of time2 could be derived
    ! as (1.0 - wt1).
    !
    ! !USES:
    use clm_time_manager, only : get_curr_calday, get_days_per_year
    !
    ! !ARGUMENTS:
    real(r8), intent(out)          :: wt1       ! weight of time1 (the left-hand time point)
    integer , intent(in), optional :: offset    ! offset (in sec.) from the current model time
    !
    ! !LOCAL VARIABLES:
    integer  :: l_offset           ! local version of offset
    real(r8) :: cday               ! current calendar day (1.0 = 0Z on Jan 1)
    real(r8) :: days_per_year      ! days per year
    !-----------------------------------------------------------------------

    if (present(offset)) then
       l_offset = offset
    else
       l_offset = 0._r8
    end if

    cday          = get_curr_calday(offset=l_offset)
    days_per_year = get_days_per_year()

    wt1 = ((days_per_year + 1._r8) - cday)/days_per_year

  end subroutine dyn_time_weights

end module dynUtilsMod
