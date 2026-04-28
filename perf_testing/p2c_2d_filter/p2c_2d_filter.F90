module p2c_2d_filter_mod

  !-----------------------------------------------------------------------
  ! Standalone extraction of p2c_2d_filter from src/main/subgridAveMod.F90.
  !
  ! Stage 1: verbatim copy. Still depends on ColumnType (col) and PatchType
  ! (patch); does not yet compile outside the CTSM build. Subsequent stages
  ! convert the col/patch lookups into array arguments and then drop the
  ! shr_kind_mod dependency.
  !-----------------------------------------------------------------------

  use shr_kind_mod , only : r8 => shr_kind_r8
  use ColumnType   , only : col
  use PatchType    , only : patch

  implicit none
  private

  public :: p2c_2d_filter

contains

  !-----------------------------------------------------------------------
  subroutine p2c_2d_filter (lev, numfc, filterc, patcharr, colarr)
    !
    ! !DESCRIPTION:
    ! perform patch to column averaging for multi level patch arrays
    !
    ! !ARGUMENTS:
    integer , intent(in)  :: lev
    integer , intent(in)  :: numfc
    integer , intent(in)  :: filterc(numfc)
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
          do p = col%patchi(c), col%patchf(c)
             if (patch%active(p)) colarr(c,j) = colarr(c,j) + patcharr(p,j) * patch%wtcol(p)
          end do
       end do
    end do

  end subroutine p2c_2d_filter

end module p2c_2d_filter_mod
