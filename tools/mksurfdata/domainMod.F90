
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
  use nanMod      , only : nan, bigint
!
! !PUBLIC TYPES:
  implicit none
  private
!
  public :: domain_type

  character*16,parameter, private :: domain_set   = 'domain_set      '
  character*16,parameter, private :: domain_unset = 'NOdomain_unsetNO'
  real(r8), parameter, private :: area_init = -9999._r8
  real(r8), parameter, private :: frac_init = -9999._r8
  real(r8), parameter, private :: topo_init = 0.0_r8
  integer,  parameter, private :: mask_init = bigint

  type domain_type
     integer          :: ni,nj                        ! size of global arrays (lsmlon,lsmlat)
     logical          :: fullgrid                     ! fullgrid
     real(r8)         :: edgen                        ! lsmedge north
     real(r8)         :: edgee                        ! lsmedge east
     real(r8)         :: edges                        ! lsmedge south
     real(r8)         :: edgew                        ! lsmedge west
     integer ,pointer :: numlon(:)                    ! numlon
     integer ,pointer :: mask(:,:)                    ! land mask: 1 = land. 0 = ocean
     real(r8),pointer :: frac(:,:)                    ! fractional land
     real(r8),pointer :: topo(:,:)                    ! topography,elevation (m)
     real(r8),pointer :: latixy(:,:)                  ! latitude of grid cell (deg)
     real(r8),pointer :: longxy(:,:)                  ! longitude of grid cell (deg)
     real(r8),pointer :: area(:,:)                    ! grid cell area (km**2)
     real(r8),pointer :: lats(:,:)                    ! grid cell latitude, S edge (deg)
     real(r8),pointer :: latn(:,:)                    ! grid cell latitude, N edge (deg)
     real(r8),pointer :: lonw(:,:)                    ! grid cell longitude, W edge (deg)
     real(r8),pointer :: lone(:,:)                    ! grid cell longitude, E edge (deg)
     character*16     :: domain_set = domain_unset    ! flag to check if domain is set
  end type domain_type

!
! !PUBLIC MEMBER FUNCTIONS:
  public domain_init          ! allocates/nans domain types
  public domain_isSet         ! if domain is set
  public domain_clean         ! deallocate domain
  public domain_setptrs       ! sets external pointer arrays into domain
  public domain_check         ! write out domain stats
!
!
! !REVISION HISTORY:
! Originally clm_varsur by Mariana Vertenstein
! Migrated from clm_varsur to domainMod by T Craig
!
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
!
! !LOCAL VARIABLES:
    integer ier
!
!EOP
!------------------------------------------------------------------------------
    if ( domain_IsSet(domain) )then
        return
!       call domain_clean(domain)
    endif

    allocate(domain%mask(ni,nj),domain%frac(ni,nj),domain%latixy(ni,nj), &
             domain%longxy(ni,nj),domain%area(ni,nj),domain%topo(ni,nj), &
             stat=ier)
    if (ier /= 0) then
       write(6,*) 'domain_init ERROR: allocate mask, frac, lat, lon, area '
       stop
    endif
    allocate(domain%lats(ni,nj),domain%latn(ni,nj),domain%lonw(ni,nj),domain%lone(ni,nj), &
       stat=ier)
    if (ier /= 0) then
       write(6,*) 'domain_init ERROR: allocate lats, latn, lonw, lone'
       stop
    endif
    allocate(domain%numlon(nj),stat=ier)
    if (ier /= 0) then
       write(6,*) 'domain_init ERROR: allocate numlon'
       stop
    endif

    domain%ni       = ni
    domain%nj       = nj
    domain%fullgrid = .true.
    domain%numlon   = bigint
    domain%edgen    = nan
    domain%edgee    = nan
    domain%edges    = nan
    domain%edgew    = nan
    domain%mask     = bigint
    domain%frac     = frac_init
    domain%topo     = topo_init
    domain%area     = area_init
    domain%latixy   = nan
    domain%longxy   = nan
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
!
! !LOCAL VARIABLES:
    integer ier
!
!EOP
!------------------------------------------------------------------------------
    if ( domain_IsSet(domain) )then
       write(6,*) 'domain_clean: cleaning ',domain%ni,domain%nj
       deallocate(domain%mask,domain%frac,domain%latixy, &
              domain%longxy,domain%area,domain%topo,stat=ier)
       if (ier /= 0) then
          write(6,*) 'domain_clean ERROR: deallocate mask, frac, lat, lon, area '
          stop
       endif
       deallocate(domain%lats,domain%latn,domain%lonw,domain%lone, &
          stat=ier)
       if (ier /= 0) then
          write(6,*) 'domain_clean ERROR: deallocate lats, latn, lonw, lone'
          stop
       endif
       deallocate(domain%numlon,stat=ier)
       if (ier /= 0) then
          write(6,*) 'domain_clean ERROR: deallocate numlon'
          stop
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
  subroutine domain_setptrs(domain,ni,nj,mask,frac,topo,latixy,longxy,area, &
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
    real(r8),optional,pointer  :: topo(:,:)
    real(r8),optional,pointer  :: latixy(:,:)
    real(r8),optional,pointer  :: longxy(:,:)
    real(r8),optional,pointer  :: area(:,:)
    real(r8),optional,pointer  :: lats(:,:)
    real(r8),optional,pointer  :: latn(:,:)
    real(r8),optional,pointer  :: lonw(:,:)
    real(r8),optional,pointer  :: lone(:,:)
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
!
!EOP
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
    if (present(topo)) then
      topo => domain%topo
    endif
    if (present(latixy)) then
      latixy => domain%latixy
    endif
    if (present(longxy)) then
      longxy => domain%longxy
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
!BOP
!
! !IROUTINE: domain_isSet
!
! !INTERFACE:
  logical function domain_isSet(domain, areaset, fracset, maskset, toposet)
!
! !DESCRIPTION:
! This returns .true. if the domain is set. Also optionally returns
! logical variables on the status of if specific variables are set in the domain.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type), intent(IN)  :: domain        ! domain datatype
    logical, optional, intent(OUT) :: areaset       ! area is set
    logical, optional, intent(OUT) :: fracset       ! land fraction is set
    logical, optional, intent(OUT) :: maskset       ! land mask is set
    logical, optional, intent(OUT) :: toposet       ! topography is set
!
! !REVISION HISTORY:
!   Created by E. Kluzek
!
!
! !LOCAL VARIABLES:
    integer :: i, j   ! indices
!
!EOP
!------------------------------------------------------------------------------
    if ( index(domain%domain_set,domain_set) == 1 ) then
        domain_isSet = .true.
        if ( present(areaset) ) areaset = .not. all( domain%area == area_init )
        if ( present(fracset) ) fracset = .not. all( domain%frac == frac_init )
        if ( present(maskset) ) maskset = .not. all( domain%mask == mask_init )
        if ( present(toposet) ) toposet = .not. all( domain%topo == topo_init )
    else
        domain_isSet = .false.
        if ( present(areaset) ) areaset = domain_isSet
        if ( present(fracset) ) fracset = domain_isSet
        if ( present(maskset) ) maskset = domain_isSet
        if ( present(toposet) ) toposet = domain_isSet
    endif

end function domain_isSet
!------------------------------------------------------------------------------


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
!
! !LOCAL VARIABLES:
!
!EOP
!------------------------------------------------------------------------------

    write(6,*) ' '
    write(6,*) 'domain_check domain_set= ',trim(domain_set)
    write(6,*) 'domain_check ni,nj     = ',domain%ni,domain%nj
    write(6,*) 'domain_check edgeNESW  = ',domain%edgen,domain%edgee,domain%edges,domain%edgew
    write(6,*) 'domain_check longxy = ',minval(domain%longxy),maxval(domain%longxy)
    write(6,*) 'domain_check latixy = ',minval(domain%latixy),maxval(domain%latixy)
    write(6,*) 'domain_check mask   = ',minval(domain%mask),maxval(domain%mask)
    write(6,*) 'domain_check frac   = ',minval(domain%frac),maxval(domain%frac)
    write(6,*) 'domain_check topo   = ',minval(domain%topo),maxval(domain%topo)
    write(6,*) 'domain_check area   = ',minval(domain%area),maxval(domain%area)
    write(6,*) 'domain_check latn   = ',minval(domain%latn),maxval(domain%latn)
    write(6,*) 'domain_check lone   = ',minval(domain%lone),maxval(domain%lone)
    write(6,*) 'domain_check lats   = ',minval(domain%lats),maxval(domain%lats)
    write(6,*) 'domain_check lonw   = ',minval(domain%lonw),maxval(domain%lonw)

end subroutine domain_check
!------------------------------------------------------------------------------

end module domainMod
