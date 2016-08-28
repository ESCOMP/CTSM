module EDVecCohortType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! cohortype.  mimics CLM vector subgrid types.  For now this holds ED data that is
  ! necessary in the rest of CLM
  !
  ! !USES:

  ! !PUBLIC TYPES:
  implicit none
  public
  !
  type, public :: ed_vec_cohort_type
     integer :: cohorts_per_column
     integer , pointer :: column(:) !index into column level quantities
   contains
     procedure, public :: Init
  end type ed_vec_cohort_type

  type(ed_vec_cohort_type), public :: ed_vec_cohort
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, beg, end)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(ed_vec_cohort_type) :: this
    integer, intent(in) :: beg, end
    !------------------------------------------------------------------------

    ! FIX(SPM,032414) pull this out and put in own ED source

    allocate(this%column(beg:end))

  end subroutine Init

end module EDVecCohortType
