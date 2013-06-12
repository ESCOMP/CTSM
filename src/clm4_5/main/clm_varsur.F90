
module clm_varsur

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varsur
!
! !DESCRIPTION:
! Module containing 2-d surface boundary data information
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! land model grid - moved to domainMod
!
! surface boundary data, these are all "gdc" local 
! Note that some of these need to be pointers (as opposed to just allocatable arrays) to
! match the ncd_io interface; for consistency, we make them all pointers
!
  real(r8), pointer :: wt_lunit(:,:)     ! weight of each landunit on the grid cell

  real(r8), pointer :: wt_nat_pft(:,:)   ! for natural veg landunit, weight of each
                                         ! pft on the landunit (adds to 1.0 on the
                                         ! landunit for all all grid cells, even
                                         ! those without any natural pft)
                                         ! (second dimension goes natpft_lb:natpft_ub)

  real(r8), pointer :: wt_cft(:,:)       ! for crop landunit, weight of each
                                         ! cft on the landunit (adds to 1.0 on the
                                         ! landunit for all all grid cells, even
                                         ! those without any crop)
                                         ! (second dimension goes cft_lb:cft_ub)

  real(r8), pointer :: wt_glc_mec(:,:)   ! for glc_mec landunits, weight of glacier
                                         ! in each elevation class (adds to 1.0 on the
                                         ! landunit for all grid cells, even those
                                         ! without any glacier)

  real(r8), pointer :: pctspec(:)        ! percent of spec lunits wrt gcell

  real(r8), pointer :: topo_glc_mec(:,:) ! subgrid glacier_mec sfc elevation

!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 2005-11-01 Moved grid to domainMod, T Craig
!
!EOP
!-----------------------------------------------------------------------

end module clm_varsur
