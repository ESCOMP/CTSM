#include <misc.h>
#include <preproc.h>

module initGridIndexMod

!------------------------------------------------------------------------------
!BOP
!
! !MODULE: initGridIndexMod
!
! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
  public set_index_xy  ! Obtain s->n index foreach each decomposition index
  public set_index_sn  ! Obtain gridcell i,j indices for model decomposition 
                       ! and s->n ordering     
!
! !DESCRIPTION:
! Obtain i,j indices for each land gridcell for model decomposition
! and s->n ordering     
!
! !REVISION HISTORY:
! 7/2004 Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------------

contains
  
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_index_xy
!
! !INTERFACE:
   subroutine set_index_xy(numg, ixy, jxy, ixy_sn, jxy_sn)
!
! !DESCRIPTION:
! Obtain i,j indices for each land gridcell for model decomposition and s->n ordering     
!
! !USES
     use decompMod, only : gcelldc, gcellsn
!
! !ARGUMENTS
     implicit none
     integer, intent(in)  :: numg
     integer, intent(out) :: ixy(numg)
     integer, intent(out) :: jxy(numg)
     integer, intent(out) :: ixy_sn(numg)
     integer, intent(out) :: jxy_sn(numg)
! !REVISION HISTORY:
! 2003.09.12  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES
     integer :: g     ! indices
!------------------------------------------------------------------------------

!dir$ concurrent
!cdir nodep
     do g = 1,numg
        ixy(g)    = gcelldc(g)%ixy
        jxy(g)    = gcelldc(g)%jxy
        ixy_sn(g) = gcellsn(g)%ixy
        jxy_sn(g) = gcellsn(g)%jxy
     end do

   end subroutine set_index_xy

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_index_sn
!
! !INTERFACE:
   subroutine set_index_sn(numg, snindex, type1d)
!
! !DESCRIPTION:
! Obtain s->n index foreach each corresponding decomposition index
!
! !USES
     use clmtype  , only : nameg, namel, namec, namep
     use decompMod, only : gcelldc, gcellsn
!
! !ARGUMENTS
     implicit none
     integer         , intent(in) :: numg
     character(len=*), intent(in) :: type1d
     integer         , pointer    :: snindex(:)
! !REVISION HISTORY:
! 2003.09.12  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES
     integer :: g,l,c,p,gsn    ! incides
!------------------------------------------------------------------------------

     if (type1d == nameg) then

        do g = 1,numg
           gsn = gcelldc(g)%gsn
           snindex(g) = gsn
        end do

     else if (type1d == namel) then

        do g = 1,numg
           gsn = gcelldc(g)%gsn
           do l = 1, gcelldc(g)%lf - gcelldc(g)%li + 1
              snindex(l) = gcellsn(gsn)%li + l - 1
           end do
        end do

     else if (type1d == namec) then

        do g = 1,numg
           gsn = gcelldc(g)%gsn
           do c = 1,gcelldc(g)%cf - gcelldc(g)%ci + 1
              snindex(c) = gcellsn(gsn)%ci + c - 1
           end do
        end do

     else if (type1d == namep) then

        do g = 1,numg
           gsn = gcelldc(g)%gsn
           do p = 1,gcelldc(g)%pf - gcelldc(g)%pi + 1
              snindex(p) = gcellsn(gsn)%pi + p - 1
           end do
        end do

     end if

   end subroutine set_index_sn

end module initGridIndexMod
