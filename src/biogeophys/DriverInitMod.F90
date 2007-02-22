#include <misc.h>
#include <preproc.h>

module DriverInitMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: DriverInitMod
!
! !DESCRIPTION:
! Initialization of driver variables needed from previous timestep
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: DriverInit
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: DriverInit
!
! !INTERFACE:
  subroutine DriverInit(lbc, ubc, lbp, ubp, &
             num_nolakec, filter_nolakec, num_lakec, filter_lakec)
!
! !DESCRIPTION:
! Initialization of driver variables needed from previous timestep
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype
    use clm_varpar   , only : nlevsoi, nlevsno
    use subgridAveMod, only : p2c

!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                    ! column-index bounds
    integer, intent(in) :: lbp, ubp                    ! pft-index bounds
    integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
    integer, intent(in) :: num_lakec                   ! number of column non-lake points in column filter
    integer, intent(in) :: filter_lakec(ubc-lbc+1)     ! column filter for non-lake points
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in variables
!
    real(r8), pointer :: pwtgcell(:)           ! weight of pft wrt corresponding gridcell
    integer , pointer :: snl(:)                ! number of snow layers
    real(r8), pointer :: h2osno(:)             ! snow water (mm H2O)
    integer , pointer :: frac_veg_nosno_alb(:) ! fraction of vegetation not covered by snow (0 OR 1) [-]
    integer , pointer :: frac_veg_nosno(:)     ! fraction of vegetation not covered by snow (0 OR 1 now) [-] (pft-level)
    real(r8), pointer :: h2osoi_ice(:,:)       ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)       ! liquid water (kg/m2)
!
! local pointers to original implicit out variables
!
    logical , pointer :: do_capsnow(:)         ! true => do snow capping
    real(r8), pointer :: h2osno_old(:)         ! snow water (mm H2O) at previous time step
    real(r8), pointer :: frac_iceold(:,:)      ! fraction of ice relative to the tot water
!
! !OTHER LOCAL VARIABLES:
!
    integer :: c, p, f, j, fc            ! indices
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (column-level)

    snl                => clm3%g%l%c%cps%snl
    h2osno             => clm3%g%l%c%cws%h2osno
    h2osno_old         => clm3%g%l%c%cws%h2osno_old
    do_capsnow         => clm3%g%l%c%cps%do_capsnow
    frac_iceold        => clm3%g%l%c%cps%frac_iceold
    h2osoi_ice         => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq         => clm3%g%l%c%cws%h2osoi_liq
    frac_veg_nosno_alb => clm3%g%l%c%p%pps%frac_veg_nosno_alb
    frac_veg_nosno     => clm3%g%l%c%p%pps%frac_veg_nosno

    ! Assign local pointers to derived type members (pft-level)

    pwtgcell           => clm3%g%l%c%p%wtgcell

!dir$ concurrent
!cdir nodep
    do c = lbc, ubc

      ! Save snow mass at previous time step
      h2osno_old(c) = h2osno(c)

      ! Decide whether to cap snow
      if (h2osno(c) > 1000._r8) then
         do_capsnow(c) = .true.
      else
         do_capsnow(c) = .false.
      end if

    end do

    ! Initialize fraction of vegetation not covered by snow (pft-level)

!dir$ concurrent
!cdir nodep
    do p = lbp,ubp
       if (pwtgcell(p)>0._r8) then
          frac_veg_nosno(p) = frac_veg_nosno_alb(p)
       else
          frac_veg_nosno(p) = 0._r8
       end if
    end do

    ! Initialize set of previous time-step variables
    ! Ice fraction of snow at previous time step

    do j = -nlevsno+1,0
!dir$ concurrent
!cdir nodep
      do f = 1, num_nolakec
         c = filter_nolakec(f)
         if (j >= snl(c) + 1) then
            frac_iceold(c,j) = h2osoi_ice(c,j)/(h2osoi_liq(c,j)+h2osoi_ice(c,j))
         end if
      end do
    end do

  end subroutine DriverInit

end module DriverInitMod
