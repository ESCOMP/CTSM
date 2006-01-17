#include <misc.h>
#include <preproc.h>

module domainMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: domainMod
!
! !DESCRIPTION:
! Module containing 2-d global surface boundary data information
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use nanMod
  use abortutils  , only : endrun
!
! !PUBLIC TYPES:
  implicit none
  private
!
  public :: domain_type

  type domain_type
     integer          :: ni,nj         ! size of global arrays (lsmlon,lsmlat)
     real(r8)         :: edges(4)      ! edges (N,E,S,W)
     integer ,pointer :: mask(:,:)     ! land mask: 1 = land. 0 = ocean
     real(r8),pointer :: frac(:,:)     ! fractional land
     real(r8),pointer :: latc(:,:)     ! latitude of grid cell (deg)
     real(r8),pointer :: lonc(:,:)     ! longitude of grid cell (deg)
     real(r8),pointer :: area(:,:)     ! grid cell area (km**2)
     real(r8),pointer :: lats(:,:)     ! grid cell latitude, S edge (deg)
     real(r8),pointer :: latn(:,:)     ! grid cell latitude, N edge (deg)
     real(r8),pointer :: lonw(:,:)     ! grid cell longitude, W edge (deg)
     real(r8),pointer :: lone(:,:)     ! grid cell longitude, E edge (deg)
     character*16     :: domain_set    ! flag to check if domain is set
  end type domain_type

  type(domain_type),public :: ldomain
  type(domain_type),public :: adomain

!
! !PUBLIC MEMBER FUNCTIONS:
  public domain_init          ! allocates/nans domain types
  public domain_setptrs       ! sets external pointer arrays into domain
!
!
! !REVISION HISTORY:
! Originally clm_varsur by Mariana Vertenstein
! Migrated from clm_varsur to domainMod by T Craig
!
  character*16,parameter :: domain_set   = 'domain_set      '
  character*16,parameter :: domain_unset = 'NOdomain_unsetNO'
!
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_init
!
! !INTERFACE:
  subroutine domain_init(domain,ni,nj)
!
! !DESCRIPTION:
! This subroutine allocates and nans the domain type
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type) :: domain        ! domain datatype
    integer           :: ni,nj      ! grid size, 2d
!
! !REVISION HISTORY:
!   Created by T Craig
!
!EOP
!
! LOCAL VARIABLES:
    integer ier
!
!------------------------------------------------------------------------------
    if (domain%domain_set == domain_set) then
       call domain_clean(domain)
    endif

    allocate(domain%mask(ni,nj),domain%frac(ni,nj),domain%latc(ni,nj), &
             domain%lonc(ni,nj),domain%area(ni,nj),stat=ier)
    if (ier /= 0) then
       write(6,*) 'domain_init ERROR: allocate mask, frac, lat, lon, area '
       call endrun()
    endif
    allocate(domain%lats(ni,nj),domain%latn(ni,nj),domain%lonw(ni,nj),domain%lone(ni,nj), &
       stat=ier)
    if (ier /= 0) then
       write(6,*) 'domain_init ERROR: allocate lats, latn, lonw, lone'
       call endrun()
    endif

    domain%ni       = ni
    domain%nj       = nj
    domain%edges    = nan
    domain%mask     = bigint
    domain%frac     = nan
    domain%latc     = nan
    domain%lonc     = nan
    domain%area     = nan
    domain%lats     = nan
    domain%latn     = nan
    domain%lonw     = nan
    domain%lone     = nan

    domain%domain_set = domain_set

end subroutine domain_init
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_clean
!
! !INTERFACE:
  subroutine domain_clean(domain)
!
! !DESCRIPTION:
! This subroutine deallocates the domain type
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type) :: domain        ! domain datatype
!
! !REVISION HISTORY:
!   Created by T Craig
!
!EOP
!
! LOCAL VARIABLES:
    integer ier
!
!------------------------------------------------------------------------------
    if (domain%domain_set == domain_set) then
       write(6,*) 'domain_clean: cleaning ',domain%ni,domain%nj
       deallocate(domain%mask,domain%frac,domain%latc, &
              domain%lonc,domain%area,stat=ier)
       if (ier /= 0) then
          write(6,*) 'domain_clean ERROR: deallocate mask, frac, lat, lon, area '
          call endrun()
       endif
       deallocate(domain%lats,domain%latn,domain%lonw,domain%lone, &
          stat=ier)
       if (ier /= 0) then
          write(6,*) 'domain_clean ERROR: deallocate lats, latn, lonw, lone'
          call endrun()
       endif
    else
       write(6,*) 'domain_clean WARN: clean domain unecessary '
    endif

    domain%ni         = bigint
    domain%nj         = bigint
    domain%domain_set = domain_unset

end subroutine domain_clean
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_setptrs
!
! !INTERFACE:
  subroutine domain_setptrs(domain,ni,nj,mask,frac,latc,lonc,area, &
     lats,latn,lonw,lone)
!
! !DESCRIPTION:
! This subroutine sets external pointer arrays to arrays in domain
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(in)  :: domain        ! domain datatype
    integer ,optional :: ni,nj      ! grid size, 2d
    integer ,optional,pointer  :: mask(:,:)
    real(r8),optional,pointer  :: frac(:,:)
    real(r8),optional,pointer  :: latc(:,:)
    real(r8),optional,pointer  :: lonc(:,:)
    real(r8),optional,pointer  :: area(:,:)
    real(r8),optional,pointer  :: lats(:,:)
    real(r8),optional,pointer  :: latn(:,:)
    real(r8),optional,pointer  :: lonw(:,:)
    real(r8),optional,pointer  :: lone(:,:)
!
! !REVISION HISTORY:
!   Created by T Craig
!
!EOP
!
! LOCAL VARIABLES:
!
!------------------------------------------------------------------------------
    if (present(ni)) then
      ni = domain%ni
    endif
    if (present(nj)) then
      nj = domain%nj
    endif
    if (present(mask)) then
      mask => domain%mask
    endif
    if (present(frac)) then
      frac => domain%frac
    endif
    if (present(latc)) then
      latc => domain%latc
    endif
    if (present(lonc)) then
      lonc => domain%lonc
    endif
    if (present(area)) then
      area => domain%area
    endif
    if (present(lats)) then
      lats => domain%lats
    endif
    if (present(latn)) then
      latn => domain%latn
    endif
    if (present(lonw)) then
      lonw => domain%lonw
    endif
    if (present(lone)) then
      lone => domain%lone
    endif

end subroutine domain_setptrs
!------------------------------------------------------------------------------

end module domainMod
