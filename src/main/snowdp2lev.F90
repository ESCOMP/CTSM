#include <misc.h>
#include <preproc.h>

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: snowdp2lev
!
! !INTERFACE:
subroutine snowdp2lev(lbc, ubc)
!
! !DESCRIPTION:
! Create snow layers and interfaces given snow depth.
! Note that cps%zi(0) is set in routine iniTimeConst.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clmtype
  use clm_varpar  , only : nlevsno
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: lbc, ubc                    ! column bounds
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
  integer , pointer :: clandunit(:)  ! landunit index associated with each column
  real(r8), pointer :: snowdp(:)     ! snow height (m)
  logical , pointer :: lakpoi(:)     ! true => landunit is a lake point
!
! local pointers to implicit out arguments
!
  integer , pointer :: snl(:)        ! number of snow layers
  real(r8), pointer :: z(:,:)        ! layer depth  (m) over snow only
  real(r8), pointer :: dz(:,:)       ! layer thickness depth (m) over snow only
  real(r8), pointer :: zi(:,:)       ! interface depth (m) over snow only
!
!
! !LOCAL VARIABLES:
!EOP
  integer :: c,l,j      !indices
!-----------------------------------------------------------------------

  ! Assign local pointers to derived subtypes components (landunit-level)

  lakpoi => clm3%g%l%lakpoi

  ! Assign local pointers to derived type members (column-level)

  clandunit => clm3%g%l%c%landunit
  snowdp    => clm3%g%l%c%cps%snowdp
  snl       => clm3%g%l%c%cps%snl
  zi        => clm3%g%l%c%cps%zi
  dz        => clm3%g%l%c%cps%dz
  z         => clm3%g%l%c%cps%z

  ! Initialize snow levels and interfaces (lake and non-lake points)

!dir$ concurrent
!cdir nodep
  do c = lbc, ubc
     dz(c,-nlevsno+1: 0) = 1.e36_r8
     z (c,-nlevsno+1: 0) = 1.e36_r8
     zi(c,-nlevsno  :-1) = 1.e36_r8
  end do

  ! Determine snow levels and interfaces for non-lake points

!dir$ concurrent
!cdir nodep
  do c = lbc,ubc
     l = clandunit(c)
     if (.not. lakpoi(l)) then
        if (snowdp(c) < 0.01_r8) then
           snl(c) = 0
           dz(c,-nlevsno+1:0) = 0._r8
           z (c,-nlevsno+1:0) = 0._r8
           zi(c,-nlevsno+0:0) = 0._r8
        else
           if ((snowdp(c) >= 0.01_r8) .and. (snowdp(c) <= 0.03_r8)) then
              snl(c) = -1
              dz(c,0)  = snowdp(c)
           else if ((snowdp(c) > 0.03_r8) .and. (snowdp(c) <= 0.04_r8)) then
              snl(c) = -2
              dz(c,-1) = snowdp(c)/2._r8
              dz(c, 0) = dz(c,-1)
           else if ((snowdp(c) > 0.04_r8) .and. (snowdp(c) <= 0.07_r8)) then
              snl(c) = -2
              dz(c,-1) = 0.02_r8
              dz(c, 0) = snowdp(c) - dz(c,-1)
           else if ((snowdp(c) > 0.07_r8) .and. (snowdp(c) <= 0.12_r8)) then
              snl(c) = -3
              dz(c,-2) = 0.02_r8
              dz(c,-1) = (snowdp(c) - 0.02_r8)/2._r8
              dz(c, 0) = dz(c,-1)
           else if ((snowdp(c) > 0.12_r8) .and. (snowdp(c) <= 0.18_r8)) then
              snl(c) = -3
              dz(c,-2) = 0.02_r8
              dz(c,-1) = 0.05_r8
              dz(c, 0) = snowdp(c) - dz(c,-2) - dz(c,-1)
           else if ((snowdp(c) > 0.18_r8) .and. (snowdp(c) <= 0.29_r8)) then
              snl(c) = -4
              dz(c,-3) = 0.02_r8
              dz(c,-2) = 0.05_r8
              dz(c,-1) = (snowdp(c) - dz(c,-3) - dz(c,-2))/2._r8
              dz(c, 0) = dz(c,-1)
           else if ((snowdp(c) > 0.29_r8) .and. (snowdp(c) <= 0.41_r8)) then
              snl(c) = -4
              dz(c,-3) = 0.02_r8
              dz(c,-2) = 0.05_r8
              dz(c,-1) = 0.11_r8
              dz(c, 0) = snowdp(c) - dz(c,-3) - dz(c,-2) - dz(c,-1)
           else if ((snowdp(c) > 0.41_r8) .and. (snowdp(c) <= 0.64_r8)) then
              snl(c) = -5
              dz(c,-4) = 0.02_r8
              dz(c,-3) = 0.05_r8
              dz(c,-2) = 0.11_r8
              dz(c,-1) = (snowdp(c) - dz(c,-4) - dz(c,-3) - dz(c,-2))/2._r8
              dz(c, 0) = dz(c,-1)
           else if (snowdp(c) > 0.64_r8) then
              snl(c) = -5
              dz(c,-4) = 0.02_r8
              dz(c,-3) = 0.05_r8
              dz(c,-2) = 0.11_r8
              dz(c,-1) = 0.23_r8
              dz(c, 0)=snowdp(c)-dz(c,-4)-dz(c,-3)-dz(c,-2)-dz(c,-1)
           endif
        end if
     end if
  end do

  ! The following loop is currently not vectorized

  do c = lbc,ubc
     l = clandunit(c)
     if (.not. lakpoi(l)) then
        do j = 0, snl(c)+1, -1
           z(c,j)    = zi(c,j) - 0.5_r8*dz(c,j)
           zi(c,j-1) = zi(c,j) - dz(c,j)
        end do
     end if
  end do

  ! Determine snow levels and interfaces for lake points

!dir$ concurrent
!cdir nodep
  do c = lbc,ubc
     l = clandunit(c)
     if (lakpoi(l)) then
        snl(c) = 0
        dz(c,-nlevsno+1:0) = 0._r8
        z (c,-nlevsno+1:0) = 0._r8
        zi(c,-nlevsno+0:0) = 0._r8
     end if
  end do

end subroutine snowdp2lev
