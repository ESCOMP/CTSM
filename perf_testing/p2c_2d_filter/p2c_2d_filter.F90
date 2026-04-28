module p2c_2d_filter_mod

  !-----------------------------------------------------------------------
  ! Standalone extraction of p2c_2d_filter from src/main/subgridAveMod.F90.
  !
  ! Stage 2: col/patch lookups replaced with array arguments. Argument
  ! attributes (pointer, etc.) are preserved from the original signature
  ! so the only change in behavior is where the arrays come from. The
  ! module now depends only on shr_kind_mod for r8.
  !-----------------------------------------------------------------------

  use shr_kind_mod , only : r8 => shr_kind_r8

  implicit none
  private

  public :: p2c_2d_filter

contains

  !-----------------------------------------------------------------------
  subroutine p2c_2d_filter (lev, numfc, filterc, &
                            patchi, patchf,      &  ! were col%patchi, col%patchf
                            active, wtcol,       &  ! were patch%active, patch%wtcol
                            patcharr, colarr)
    !
    ! !DESCRIPTION:
    ! perform patch to column averaging for multi level patch arrays
    !
    ! !ARGUMENTS:
    integer , intent(in)  :: lev
    integer , intent(in)  :: numfc
    integer , intent(in)  :: filterc(numfc)
    integer , pointer     :: patchi(:)        ! beginning patch index for each column
    integer , pointer     :: patchf(:)        ! ending patch index for each column
    logical , pointer     :: active(:)        ! true=>do computations on this patch
    real(r8), pointer     :: wtcol(:)         ! patch weight relative to column
    real(r8), pointer     :: patcharr(:,:)
    real(r8), pointer     :: colarr(:,:)
    !
    ! !LOCAL VARIABLES:
    integer :: fc,c,p,j    ! indices
    !-----------------------------------------------------------------------

    do j = 1,lev
       do fc = 1,numfc
          c = filterc(fc)
          colarr(c,j) = 0._r8
          do p = patchi(c), patchf(c)
             if (active(p)) colarr(c,j) = colarr(c,j) + patcharr(p,j) * wtcol(p)
          end do
       end do
    end do

  end subroutine p2c_2d_filter

end module p2c_2d_filter_mod
