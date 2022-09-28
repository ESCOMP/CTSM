module circle

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  real(r8), parameter, public :: pi = 3.14159265358979323846_r8

  public :: circle_area

contains

  function circle_area(r)
    real(r8), intent(in) :: r
    real(r8) :: circle_area

    circle_area = pi*r*r

  end function circle_area

end module circle
