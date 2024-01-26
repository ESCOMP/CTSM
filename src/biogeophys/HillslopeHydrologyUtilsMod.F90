module HillslopeHydrologyUtilsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Utilities used in HillslopeHydrologyMod
  !
  ! !USES:
#include "shr_assert.h"
  use decompMod      , only : bounds_type
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use spmdMod        , only : masterproc, iam
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog

  ! !PUBLIC TYPES:
  implicit none

  private
  save

  real(r8), parameter :: toosmall_distance_default  = 1e-6

  ! !PUBLIC MEMBER FUNCTIONS:
  public HillslopeSoilThicknessProfile_linear

contains

  !------------------------------------------------------------------------
  subroutine HillslopeSoilThicknessProfile_linear(bounds, soil_depth_lowland, soil_depth_upland, toosmall_distance_in)
    !
    ! !DESCRIPTION:
    ! Modify soil thickness across hillslope by changing
    ! nbedrock according to the "Linear" method
    !
    ! !USES:
    use LandunitType    , only : lun
    use ColumnType      , only : col
    use clm_varpar      , only : nlevsoi
    use clm_varcon      , only : zisoi
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in) :: soil_depth_lowland, soil_depth_upland
    real(r8), intent(in), optional :: toosmall_distance_in
    !
    ! !LOCAL VARIABLES
    real(r8) :: min_hill_dist, max_hill_dist
    real(r8) :: toosmall_distance
    real(r8) :: soil_depth_col
    real(r8) :: m, b
    integer :: c, j, l

    if (present(toosmall_distance_in)) then
       toosmall_distance = toosmall_distance_in
    else
      toosmall_distance = toosmall_distance_default
    end if

    do l = bounds%begl,bounds%endl
       min_hill_dist = minval(col%hill_distance(lun%coli(l):lun%colf(l)))
       max_hill_dist = maxval(col%hill_distance(lun%coli(l):lun%colf(l)))

       if (abs(max_hill_dist - min_hill_dist) > toosmall_distance) then
          m = (soil_depth_lowland - soil_depth_upland)/ &
               (max_hill_dist - min_hill_dist)
       else
          m = 0._r8
       end if
       b = soil_depth_upland

       do c =  lun%coli(l), lun%colf(l)
          if (col%is_hillslope_column(c) .and. col%active(c)) then
             soil_depth_col = m*(max_hill_dist - col%hill_distance(c)) + b
             do j = 1,nlevsoi
               if ((zisoi(j-1) <  soil_depth_col) .and. (zisoi(j) >= soil_depth_col)) then
                  col%nbedrock(c) = j
                  exit
               end if
             enddo
          end if
       enddo
    enddo
   end subroutine HillslopeSoilThicknessProfile_linear
end module HillslopeHydrologyUtilsMod