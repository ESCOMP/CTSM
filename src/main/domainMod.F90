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
  use spmdMod     , only : masterproc
  use abortutils  , only : endrun
  use clm_varctl  , only : iulog
!
! !PUBLIC TYPES:
  implicit none
  private
!
  public :: domain_type
  public :: latlon_type

  !--- this typically contains local domain info with arrays dim begg:endg ---
  type domain_type
     integer          :: ns         ! global size of domain
     integer          :: ni,nj      ! global axis if 2d (nj=1 if unstructured)
     integer          :: nbeg,nend  ! local beg/end indices
     character(len=8) :: clmlevel   ! grid type
     logical          :: decomped   ! decomposed locally or global copy
     logical          :: regional   ! regional or global grid
     logical          :: areaset    ! has area been set
     integer ,pointer :: mask(:)    ! land mask: 1 = land, 0 = ocean
     real(r8),pointer :: frac(:)    ! fractional land
     real(r8),pointer :: topo(:)    ! topography
     real(r8),pointer :: latc(:)    ! latitude of grid cell (deg)
     real(r8),pointer :: lonc(:)    ! longitude of grid cell (deg)
     real(r8),pointer :: area(:)    ! grid cell area (km**2)
     real(r8),pointer :: asca(:)    ! area scaling from CESM driver
     character*16     :: set        ! flag to check if domain is set
     integer ,pointer :: glcmask(:) ! glc mask: 1=sfc mass balance required by GLC component
                                    !           0=SMB not required (default)
     !--- following are valid only for land domain ---
     integer ,pointer :: pftm(:)    ! pft mask: 1=real, 0=fake, -1=notset
     real(r8),pointer :: nara(:)    ! normalized area in upscaling (km**2),
     real(r8),pointer :: ntop(:)    ! normalized topo for downscaling (m)
  end type domain_type

  !--- this contains global info about a grid, lats and lons are 1d
  !--- global arrays of size ni or nj which assume regular lat/lon grids only
  type latlon_type
     integer          :: ns         ! global size of domain
     integer          :: ni,nj      ! global axis if 2d (nj=1 if unstructured)
     character*16     :: set        ! flag to check if domain is set
     logical          :: regional   ! regional or global grid
     real(r8)         :: edges(4)   ! global edges (N,E,S,W)
     real(r8),pointer :: latc(:)    ! latitude of 1d grid cell (deg)
     real(r8),pointer :: lonc(:)    ! longitude of 1d grid cell (deg)
     real(r8),pointer :: lats(:)    ! latitude of 1d south grid cell edge (deg)
     real(r8),pointer :: latn(:)    ! latitude of 1d north grid cell edge (deg)
     real(r8),pointer :: lonw(:)    ! longitude of 1d west grid cell edge (deg)
     real(r8),pointer :: lone(:)    ! longitude of 1d east grid cell edge (deg)
  end type latlon_type

  type(domain_type),public    :: adomain
  type(domain_type),public    :: ldomain

  type(latlon_type),public    :: alatlon
  type(latlon_type),public    :: llatlon

  integer ,pointer,public     :: gatm(:)   ! gatm pulled out of domain
  integer ,pointer,public     :: amask(:)  ! global atm mask
  integer, pointer,public     :: pftm(:)   ! pft mask for lnd grid
!
! !PUBLIC MEMBER FUNCTIONS:
  public domain_init          ! allocates/nans domain types
  public domain_clean         ! deallocates domain types
  public domain_setptrs       ! sets external pointer arrays into domain
  public domain_check         ! write out domain info
  public latlon_init          ! allocates/nans domain types
  public latlon_check         ! write out domain info
  public latlon_clean         ! deallocate domain info
  public latlon_setsame          ! copy one domain to another
!
!
! !REVISION HISTORY:
! Originally clm_varsur by Mariana Vertenstein
! Migrated from clm_varsur to domainMod by T Craig
!
  character*16,parameter :: set   = 'domain_set      '
  character*16,parameter :: unset = 'NOdomain_unsetNO'
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
  subroutine domain_init(domain,ni,nj,nbeg,nend,clmlevel)
!
! !DESCRIPTION:
! This subroutine allocates and nans the domain type
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type) :: domain        ! domain datatype
    integer           :: ni,nj         ! grid size, 2d
    integer,optional  :: nbeg,nend     ! beg/end indices
    character(len=*),optional:: clmlevel      ! grid type
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer ier
    integer nb,ne
!
!------------------------------------------------------------------------------

    nb = 1
    ne = ni*nj
    if (present(nbeg)) then
       if (present(nend)) then
          nb = nbeg
          ne = nend
       endif
    endif

    if (domain%set == set) then
       call domain_clean(domain)
    endif
    allocate(domain%mask(nb:ne),domain%frac(nb:ne),domain%latc(nb:ne), &
             domain%pftm(nb:ne),domain%area(nb:ne),domain%lonc(nb:ne), &
             domain%nara(nb:ne),domain%topo(nb:ne),domain%ntop(nb:ne), &
             domain%asca(nb:ne),domain%glcmask(nb:ne),stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'domain_init ERROR: allocate mask, frac, lat, lon, area '
       call endrun()
    endif

    if (present(clmlevel)) then
       domain%clmlevel = clmlevel
    endif

    domain%ns       = ni*nj
    domain%ni       = ni
    domain%nj       = nj
    domain%nbeg     = nb
    domain%nend     = ne
    domain%mask     = -9999
    domain%frac     = -1.0e36
    domain%topo     = 0._r8
    domain%latc     = nan
    domain%lonc     = nan
    domain%area     = nan

    domain%set      = set
    domain%regional = .false.
    domain%areaset  = .false.
    if (domain%nbeg == 1 .and. domain%nend == domain%ns) then
       domain%decomped = .false.
    else
       domain%decomped = .true.
    endif

    domain%pftm     = -9999
    domain%glcmask  = 0  
    domain%nara     = 0._r8
    domain%ntop     = -1.0e36
    domain%asca     = 1._r8

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
!EOP
    integer ier
!
!------------------------------------------------------------------------------
    if (domain%set == set) then
       if (masterproc) then
          write(iulog,*) 'domain_clean: cleaning ',domain%ni,domain%nj
       endif
       deallocate(domain%mask,domain%frac,domain%latc, &
                  domain%lonc,domain%area,domain%pftm, &
                  domain%nara,domain%topo,domain%ntop, &
                  domain%asca,domain%glcmask,stat=ier)
       if (ier /= 0) then
          write(iulog,*) 'domain_clean ERROR: deallocate mask, frac, lat, lon, area '
          call endrun()
       endif
    else
       if (masterproc) then
          write(iulog,*) 'domain_clean WARN: clean domain unecessary '
       endif
    endif

    domain%clmlevel   = unset
    domain%ns         = bigint
    domain%ni         = bigint
    domain%nj         = bigint
    domain%nbeg       = bigint
    domain%nend       = bigint
    domain%set        = unset
    domain%decomped   = .true.
    domain%regional   = .false.
    domain%areaset    = .false.

end subroutine domain_clean
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_setptrs
!
! !INTERFACE:
  subroutine domain_setptrs(domain,ns,ni,nj,nbeg,nend,decomped,regional, &
                            mask,pftm,glcmask,clmlevel, &
                            frac,topo,latc,lonc,area,nara,ntop,asca)
!
! !DESCRIPTION:
! This subroutine sets external pointer arrays to arrays in domain
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(in)  :: domain    ! domain datatype
    integer ,optional :: ns,ni,nj,nbeg,nend    ! grid size, 2d, beg/end
    character(len=*),optional :: clmlevel      ! grid type
    logical, optional :: decomped              ! decomped or global
    logical, optional :: regional              ! regional or global
    integer ,optional,pointer  :: mask(:)  
    integer ,optional,pointer  :: pftm(:)  
    integer ,optional,pointer  :: glcmask(:)  
    real(r8),optional,pointer  :: frac(:)  
    real(r8),optional,pointer  :: topo(:)  
    real(r8),optional,pointer  :: latc(:)  
    real(r8),optional,pointer  :: lonc(:)  
    real(r8),optional,pointer  :: area(:)  
    real(r8),optional,pointer  :: nara(:)  
    real(r8),optional,pointer  :: ntop(:)  
    real(r8),optional,pointer  :: asca(:)  
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
!
!EOP
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
    if (present(nbeg)) then
      nbeg = domain%nbeg
    endif
    if (present(nend)) then
      nend = domain%nend
    endif
    if (present(clmlevel)) then
      clmlevel = domain%clmlevel
    endif
    if (present(decomped)) then
      decomped = domain%decomped
    endif
    if (present(regional)) then
      regional = domain%regional
    endif
    if (present(mask)) then
      mask => domain%mask
    endif
    if (present(pftm)) then
      pftm => domain%pftm
    endif
    if (present(glcmask)) then
      glcmask => domain%glcmask
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
    if (present(asca)) then
      asca => domain%asca
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
!
! !LOCAL VARIABLES:
!
!EOP
!------------------------------------------------------------------------------

  if (masterproc) then
    write(iulog,*) '  domain_check set       = ',trim(domain%set)
    write(iulog,*) '  domain_check decomped  = ',domain%decomped
    write(iulog,*) '  domain_check regional  = ',domain%regional
    write(iulog,*) '  domain_check areaset   = ',domain%areaset
    write(iulog,*) '  domain_check ns        = ',domain%ns
    write(iulog,*) '  domain_check ni,nj     = ',domain%ni,domain%nj
    write(iulog,*) '  domain_check clmlevel  = ',trim(domain%clmlevel)
    write(iulog,*) '  domain_check nbeg,nend = ',domain%nbeg,domain%nend
    write(iulog,*) '  domain_check lonc = ',minval(domain%lonc),maxval(domain%lonc)
    write(iulog,*) '  domain_check latc = ',minval(domain%latc),maxval(domain%latc)
    write(iulog,*) '  domain_check mask = ',minval(domain%mask),maxval(domain%mask)
    write(iulog,*) '  domain_check frac = ',minval(domain%frac),maxval(domain%frac)
    write(iulog,*) '  domain_check topo = ',minval(domain%topo),maxval(domain%topo)
    write(iulog,*) '  domain_check area = ',minval(domain%area),maxval(domain%area)
    write(iulog,*) '  domain_check pftm = ',minval(domain%pftm),maxval(domain%pftm)
    write(iulog,*) '  domain_check nara = ',minval(domain%nara),maxval(domain%nara)
    write(iulog,*) '  domain_check ntop = ',minval(domain%ntop),maxval(domain%ntop)
    write(iulog,*) '  domain_check asca = ',minval(domain%asca),maxval(domain%asca)
    write(iulog,*) '  domain_check glcmask = ',minval(domain%glcmask),maxval(domain%glcmask)
    write(iulog,*) ' '
  endif

end subroutine domain_check
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: latlon_init
!
! !INTERFACE:
  subroutine latlon_init(latlon,ni,nj)
!
! !DESCRIPTION:
! This subroutine allocates and nans the latlon type
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(latlon_type) :: latlon        ! latlon datatype
    integer           :: ni,nj         ! grid size, 2d
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

    if (latlon%set == set) then
       call latlon_clean(latlon)
    endif

    allocate(latlon%latc(nj),latlon%lonc(ni), &
             latlon%lats(nj),latlon%latn(nj), &
             latlon%lonw(ni),latlon%lone(ni), &
             stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'latlon_init ERROR: allocate '
       call endrun()
    endif

    latlon%ns       = ni*nj
    latlon%ni       = ni
    latlon%nj       = nj
    latlon%latc     = nan
    latlon%lonc     = nan
    latlon%lats     = nan
    latlon%latn     = nan
    latlon%lonw     = nan
    latlon%lone     = nan
    latlon%edges    = nan

    latlon%set      = set
    latlon%regional = .false.

end subroutine latlon_init
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: latlon_clean
!
! !INTERFACE:
  subroutine latlon_clean(latlon)
!
! !DESCRIPTION:
! This subroutine allocates and nans the latlon type
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(latlon_type) :: latlon        ! latlon datatype
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

    if (latlon%set == unset) then
       return
    endif

    deallocate(latlon%latc,latlon%lonc, &
               latlon%lats,latlon%latn, &
               latlon%lonw,latlon%lone, &
               stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'latlon_clean ERROR: deallocate '
       call endrun()
    endif

    latlon%ns       = bigint
    latlon%ni       = bigint
    latlon%nj       = bigint
    latlon%edges    = nan

    latlon%set      = unset
    latlon%regional = .false.

end subroutine latlon_clean
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: latlon_setsame
!
! !INTERFACE:
  subroutine latlon_setsame(latlon1,latlon2)
!
! !DESCRIPTION:
! This subroutine copies parts of latlon2 = latlon1 specifically for
! setting a finemesh lats/lons to coarsemesh grid
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(latlon_type),intent(in)    :: latlon1        ! latlon datatype
    type(latlon_type),intent(inout) :: latlon2        ! latlon datatype
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

    if (latlon1%ni /= latlon2%ni .or. latlon1%nj /= latlon2%nj) then
       write(iulog,*) 'latlon_setsame: error on size',latlon1%ni,latlon1%nj,latlon2%ni,latlon2%nj
       call endrun()
    endif

    if (masterproc) then
       write(iulog,*) 'latlon_setsame: copying ',latlon1%ni,latlon1%nj
    endif
    latlon2%edges    = latlon1%edges
    latlon2%latc     = latlon1%latc
    latlon2%lonc     = latlon1%lonc
    latlon2%lats     = latlon1%lats
    latlon2%latn     = latlon1%latn
    latlon2%lonw     = latlon1%lonw
    latlon2%lone     = latlon1%lone

    latlon2%set      = latlon1%set
    latlon2%regional = latlon1%regional

end subroutine latlon_setsame
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: latlon_check
!
! !INTERFACE:
  subroutine latlon_check(latlon)
!
! !DESCRIPTION:
! This subroutine write latlon info
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(latlon_type),intent(in)  :: latlon        ! latlon datatype
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
!
!EOP
!------------------------------------------------------------------------------

  if (masterproc .and. latlon%set == set) then
    write(iulog,*) '  latlon_check ns        = ',latlon%ns
    write(iulog,*) '  latlon_check ni,nj     = ',latlon%ni,latlon%nj
    write(iulog,*) '  latlon_check set       = ',latlon%set
    write(iulog,*) '  latlon_check regional  = ',latlon%regional
    write(iulog,*) '  latlon_check edgeNESW  = ',latlon%edges
    write(iulog,*) '  latlon_check lonc = ',minval(latlon%lonc),maxval(latlon%lonc)
    write(iulog,*) '  latlon_check latc = ',minval(latlon%latc),maxval(latlon%latc)
    write(iulog,*) '  latlon_check lonw = ',minval(latlon%lonw),maxval(latlon%lonw)
    write(iulog,*) '  latlon_check lone = ',minval(latlon%lone),maxval(latlon%lone)
    write(iulog,*) '  latlon_check lats = ',minval(latlon%lats),maxval(latlon%lats)
    write(iulog,*) '  latlon_check latn = ',minval(latlon%latn),maxval(latlon%latn)
    write(iulog,*) ' '
  endif

end subroutine latlon_check
!------------------------------------------------------------------------------

end module domainMod
