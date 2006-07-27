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
     integer          :: ns         ! size of arrays
     integer          :: ni,nj      ! axis lenght if 2d (nj=1 if unstructured)
     real(r8)         :: edges(4)   ! edges (N,E,S,W)
     integer ,pointer :: mask(:)    ! land mask: 1 = land, 0 = ocean
     real(r8),pointer :: frac(:)    ! fractional land
     real(r8),pointer :: topo(:)    ! topography
     real(r8),pointer :: latc(:)    ! latitude of grid cell (deg)
     real(r8),pointer :: lonc(:)    ! longitude of grid cell (deg)
     real(r8),pointer :: area(:)    ! grid cell area (km**2)
     real(r8),pointer :: lats(:)    ! grid cell latitude, S edge (deg)
     real(r8),pointer :: latn(:)    ! grid cell latitude, N edge (deg)
     real(r8),pointer :: lonw(:)    ! grid cell longitude, W edge (deg)
     real(r8),pointer :: lone(:)    ! grid cell longitude, E edge (deg)
     character*16     :: domain_set ! flag to check if domain is set
     !--- following are valid only for land domain ---
     integer ,pointer :: pftm(:)    ! pft  mask: 1=real, 0=fake, -1=notset
     real(r8),pointer :: nara(:)    ! normalized area in upscaling (km**2),
     real(r8),pointer :: ntop(:)    ! normalized topo for downscaling (m)
     integer ,pointer :: gatm(:)    ! overlapping atm gridcell, 1d glo
  end type domain_type

  type(domain_type),public,pointer :: ldomain
  type(domain_type),public,pointer :: ndomains(:)
  type(domain_type),public,target  :: adomain

!
! !PUBLIC MEMBER FUNCTIONS:
  public domain_init          ! allocates/nans domain types
  public domain_setptrs       ! sets external pointer arrays into domain
  public domain_check         ! write out domain info
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

    allocate(domain%mask(ni*nj),domain%frac(ni*nj),domain%latc(ni*nj), &
             domain%pftm(ni*nj),domain%area(ni*nj),domain%lonc(ni*nj), &
             domain%nara(ni*nj),domain%topo(ni*nj),domain%ntop(ni*nj), &
             domain%gatm(ni*nj), &
             stat=ier)
    if (ier /= 0) then
       write(6,*) 'domain_init ERROR: allocate mask, frac, lat, lon, area '
       call endrun()
    endif
    allocate(domain%lats(ni*nj),domain%latn(ni*nj),domain%lonw(ni*nj),domain%lone(ni*nj), &
       stat=ier)
    if (ier /= 0) then
       write(6,*) 'domain_init ERROR: allocate lats, latn, lonw, lone'
       call endrun()
    endif

    domain%ns       = ni*nj
    domain%ni       = ni
    domain%nj       = nj
    domain%edges    = nan
    domain%mask     = bigint
    domain%frac     = -1.0e36
    domain%topo     = 0._r8
    domain%latc     = nan
    domain%lonc     = nan
    domain%area     = nan
    domain%lats     = nan
    domain%latn     = nan
    domain%lonw     = nan
    domain%lone     = nan

    domain%domain_set = domain_set

    domain%pftm     = bigint
    domain%nara     = 0._r8
    domain%ntop     = -1.0e36
    domain%gatm     = bigint

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
                  domain%lonc,domain%area,domain%pftm, &
                  domain%nara,domain%topo,domain%ntop, &
                  domain%gatm, &
                  stat=ier)
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

    domain%ns         = bigint
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
  subroutine domain_setptrs(domain,ns,ni,nj,mask,pftm, &
     frac,topo,latc,lonc,area,nara,ntop,gatm,lats,latn,lonw,lone)
!
! !DESCRIPTION:
! This subroutine sets external pointer arrays to arrays in domain
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(in)  :: domain    ! domain datatype
    integer ,optional :: ns,ni,nj              ! grid size, 2d
    integer ,optional,pointer  :: mask(:)  
    integer ,optional,pointer  :: pftm(:)  
    real(r8),optional,pointer  :: frac(:)  
    real(r8),optional,pointer  :: topo(:)  
    real(r8),optional,pointer  :: latc(:)  
    real(r8),optional,pointer  :: lonc(:)  
    real(r8),optional,pointer  :: area(:)  
    real(r8),optional,pointer  :: nara(:)  
    real(r8),optional,pointer  :: ntop(:)  
    integer ,optional,pointer  :: gatm(:)  
    real(r8),optional,pointer  :: lats(:)  
    real(r8),optional,pointer  :: latn(:)  
    real(r8),optional,pointer  :: lonw(:)  
    real(r8),optional,pointer  :: lone(:)  
!
! !REVISION HISTORY:
!   Created by T Craig
!
!EOP
!
! LOCAL VARIABLES:
!
!------------------------------------------------------------------------------
    if (present(ns)) then
      ns = domain%ns
    endif
    if (present(ni)) then
      ni = domain%ni
    endif
    if (present(nj)) then
      nj = domain%nj
    endif
    if (present(mask)) then
      mask => domain%mask
    endif
    if (present(pftm)) then
      pftm => domain%pftm
    endif
    if (present(frac)) then
      frac => domain%frac
    endif
    if (present(topo)) then
      topo => domain%topo
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
    if (present(nara)) then
      nara => domain%nara
    endif
    if (present(ntop)) then
      ntop => domain%ntop
    endif
    if (present(gatm)) then
      gatm => domain%gatm
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
!BOP
!
! !IROUTINE: domain_check
!
! !INTERFACE:
  subroutine domain_check(domain)
!
! !DESCRIPTION:
! This subroutine write domain info
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(in)  :: domain        ! domain datatype
!
! !REVISION HISTORY:
!   Created by T Craig
!
!EOP
!
! LOCAL VARIABLES:
!
!------------------------------------------------------------------------------

    write(6,*) '  domain_check domain_set= ',trim(domain_set)
    write(6,*) '  domain_check ns        = ',domain%ns
    write(6,*) '  domain_check ni,nj     = ',domain%ni,domain%nj
    write(6,*) '  domain_check edgeNESW  = ',domain%edges
    write(6,*) '  domain_check lonc = ',minval(domain%lonc),maxval(domain%lonc)
    write(6,*) '  domain_check latc = ',minval(domain%latc),maxval(domain%latc)
    write(6,*) '  domain_check mask = ',minval(domain%mask),maxval(domain%mask)
    write(6,*) '  domain_check frac = ',minval(domain%frac),maxval(domain%frac)
    write(6,*) '  domain_check topo = ',minval(domain%topo),maxval(domain%topo)
    write(6,*) '  domain_check area = ',minval(domain%area),maxval(domain%area)
    write(6,*) '  domain_check latn = ',minval(domain%latn),maxval(domain%latn)
    write(6,*) '  domain_check lone = ',minval(domain%lone),maxval(domain%lone)
    write(6,*) '  domain_check lats = ',minval(domain%lats),maxval(domain%lats)
    write(6,*) '  domain_check lonw = ',minval(domain%lonw),maxval(domain%lonw)
    write(6,*) ' '

end subroutine domain_check
!------------------------------------------------------------------------------

end module domainMod
