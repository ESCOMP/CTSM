module downscaleMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: downscaleMod
!
! !DESCRIPTION:
! Area averaging routines
! These routines are used for area-average mapping of a field from one
! grid to another for the purposes of downscaling and using the fine mesh grid
!
! !USES:
  use clm_varpar   , only : numrad
  use clm_varctl   , only : iulog
  use domainMod    , only : domain_type, domain_setptrs, latlon_type
  use shr_const_mod, only : SHR_CONST_PI
  use shr_kind_mod , only : r8 => shr_kind_r8
  use spmdMod      , only : iam,masterproc
  use abortutils   , only : endrun
  use nanMod      
!
! !PUBLIC TYPES:
  implicit none
  private

  type map_type
     private
     ! lower level in hierarchy
     character(len=32) :: name
     character(len=16) :: type        ! global, dst, src, etc
     integer           :: ni_i,nj_i   ! size of src grid ni,nj
     integer           :: ni_o,nj_o   ! size of dst grid ni,nj
     integer           :: nwts        ! size of row, col, S (local)
     integer , pointer :: src(:)      ! src index (COL)
     integer , pointer :: dst(:)      ! dst index (ROW)
     real(r8), pointer :: S(:)        ! wt of overlap input cell
     integer           :: dstmo       ! max num of overlaps of dst cell
  end type map_type
  public map_type

  type(map_type),public        :: map1dl_a2l   ! a2l mapping 1d loc glo to gdc
  type(map_type),public        :: map1dl_l2a   ! l2a mapping 1d loc gdc to glo

  character(len=16),parameter,public :: map_typelocal  = 'local'
  character(len=16),parameter,public :: map_typeglobal = 'global'

!
! !PUBLIC MEMBER FUNCTIONS:
  public :: map_init
  public :: map_setptrs
  public :: map_setmapsFM
  public :: map_setgatmFM
  public :: map_maparrayl
  public :: map_maparrayg

  interface map_maparrayl
     module procedure map_maparrayl_arr
     module procedure map_maparrayl_rev
     module procedure map_maparrayl_av
  end interface
!
! !REVISION HISTORY:
! Created by Sam Levis
! Updated to clm2.1 data structures by Mariana Vertenstein
! 2005.11.01 Updated and cleaned by T Craig
!
!
! !PRIVATE MEMBER FUNCTIONS:
!EOP
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_init
!
! !INTERFACE:
  subroutine map_init(map,domain_i,domain_o,nwts,name,type)
!
! !DESCRIPTION:
! This subroutine initializes the map datatype
!
! !USES:
!
! !ARGUMENTS:
  implicit none
  type(map_type)    , intent(inout)       :: map
  type(domain_type) , intent(in),target   :: domain_i
  type(domain_type) , intent(in),target   :: domain_o
  integer           , intent(in)          :: nwts     ! number of wts
  character(len=*)  , intent(in),optional :: name
  character(len=*)  , intent(in),optional :: type
!
! !REVISION HISTORY:
! 2006.06.28  T Craig  Creation.
!
!
! !LOCAL VARIABLES:
!EOP
  integer ier    ! error flag
!------------------------------------------------------------------------------

  map%ni_i = domain_i%ni
  map%nj_i = domain_i%nj
  map%ni_o = domain_o%ni
  map%nj_o = domain_o%nj
  if (present(name)) then
    map%name = trim(name)
  else
    map%name = 'unset'
  endif
  if (present(type)) then
    map%type = trim(type)
  else
    map%type = 'unset'
  endif

  map%nwts = nwts
  allocate(map%src(nwts),map%dst(nwts),map%S(nwts),stat=ier)
  if (ier /= 0) then
     write(iulog,*) 'map_init ERROR: allocate map'
     call endrun()
  endif

  map%src = -1
  map%dst = -1
  map%S   = nan

end subroutine map_init

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_setptrs
!
! !INTERFACE:
  subroutine map_setptrs(map,name,type,ni_i,nj_i,ni_o,nj_o, &
     nwts,src,dst,wts,dstmo)
!
! !DESCRIPTION:
! This subroutine sets external pointer arrays to arrays in map
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(map_type)   ,intent(in)        :: map
    character(len=*) ,optional          :: name     
    character(len=*) ,optional          :: type
    integer          ,optional          :: ni_i,nj_i
    integer          ,optional          :: ni_o,nj_o
    integer          ,optional          :: nwts
    integer          ,optional,pointer  :: src(:)
    integer          ,optional,pointer  :: dst(:)
    real(r8)         ,optional,pointer  :: wts(:)
    integer          ,optional          :: dstmo
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
!
!EOP
!------------------------------------------------------------------------------
    if (present(name)) then
      name = map%name
    endif
    if (present(type)) then
      type = map%type
    endif
    if (present(ni_i)) then
      ni_i = map%ni_i
    endif
    if (present(nj_i)) then
      nj_i = map%nj_i
    endif
    if (present(ni_o)) then
      ni_o = map%ni_o
    endif
    if (present(nj_o)) then
      nj_o = map%nj_o
    endif
    if (present(nwts)) then
      nwts = map%nwts
    endif
    if (present(src)) then
      src => map%src
    endif
    if (present(dst)) then
      dst => map%dst
    endif
    if (present(wts)) then
      wts => map%S
    endif
    if (present(dstmo)) then
      dstmo = map%dstmo
    endif

end subroutine map_setptrs
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_maparrayl_arr
!
! !INTERFACE:
  subroutine map_maparrayl_arr(begg_i, endg_i, begg_o, endg_o, nflds, &
                               fld_i, fld_o, map)
!
! !DESCRIPTION:
! This subroutine maps arrays, local 1d with 1d map type
!
! !USES:
!
! !ARGUMENTS:
  implicit none
  integer                         :: begg_i,endg_i         !beg,end of input grid
  integer                         :: begg_o,endg_o         !beg,end of output grid
  integer                         :: nflds                 !number of fields being mapped
  real(r8), intent(in)            :: fld_i(begg_i:endg_i,nflds)
  real(r8), intent(out)           :: fld_o(begg_o:endg_o,nflds)
  type(map_type), intent(in)      :: map
!
! !REVISION HISTORY:
! 2005.11.15  T Craig  Creation.
! 2006.3.30   P Worley Restructuring for improved vector performance
!

! !LOCAL VARIABLES:
!EOP
  integer :: g_o ,g_i              ! gridcell indices
  integer :: n,ifld                ! loop counters
  real(r8):: wtx                   ! wt for map

!
!------------------------------------------------------------------------------

    if (trim(map%type) /= trim(map_typelocal)) then
       write(iulog,*) 'map_maparrayl_arr WARNING: map type not correct, ', &
                   map%name,map%type
    endif

    ! check for errors in overlap

    if (minval(map%src) < begg_i .or. maxval(map%src) > endg_i) then
       write(iulog,*) 'map_maparrayl_arr ERROR: src out of bounds:', &
                   minval(map%src),maxval(map%src),begg_i,endg_i
       call endrun()
    endif

    if (minval(map%dst) < begg_o .or. maxval(map%dst) > endg_o) then
       write(iulog,*) 'map_maparrayl_arr ERROR: dst out of bounds:', &
                   minval(map%dst),maxval(map%dst),begg_o,endg_o
       call endrun()
    endif

    ! initialize field on output grid to zero everywhere

    fld_o(:,:) = 0._r8

    ! map flds

    do ifld = 1,nflds
    do n = 1,map%nwts
       g_i = map%src(n)
       g_o = map%dst(n)
       wtx = map%S(n)
       fld_o(g_o,ifld) = fld_o(g_o,ifld) + wtx * fld_i(g_i,ifld)
    enddo
    enddo

end subroutine map_maparrayl_arr

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_maparrayl_rev
!
! !INTERFACE:
  subroutine map_maparrayl_rev(begg_i, endg_i, begg_o, endg_o, nflds, &
                              fld_i, fld_o, map, reverse_order)
!
! !DESCRIPTION:
! This subroutine maps arrays, local 1d with 1d map type
!
! !USES:
  use mct_mod
!
! !ARGUMENTS:
  implicit none
  integer                         :: begg_i,endg_i         !beg,end of input grid
  integer                         :: begg_o,endg_o         !beg,end of output grid
  integer                         :: nflds                 !number of fields being mapped
  real(r8), intent(in)            :: fld_i(nflds,endg_i-begg_i+1)
  real(r8), intent(out)           :: fld_o(nflds,endg_o-begg_o+1)
  type(map_type), intent(in)      :: map
  logical, intent(in)             :: reverse_order
!
! !REVISION HISTORY:
! 2005.11.15  T Craig  Creation.
! 2006.3.30   P Worley Restructuring for improved vector performance
! 2010.5.1    M Vertenstein Added new interfaces
!

! !LOCAL VARIABLES:
!EOP
  integer :: g_o ,g_i              ! gridcell indices
  integer :: n,ifld                ! loop counters
  real(r8):: wtx                   ! wt for map
!
!------------------------------------------------------------------------------

    if (trim(map%type) /= trim(map_typelocal)) then
       write(iulog,*) 'map_maparrayl_rev WARNING: map type not correct, ', &
                   map%name,map%type
    endif

    ! check for errors in overlap

    if (minval(map%src) < begg_i .or. maxval(map%src) > endg_i) then
       write(iulog,*) 'map_maparrayl_rev ERROR: src out of bounds:', &
                   minval(map%src),maxval(map%src),begg_i,endg_i
       call endrun()
    endif

    if (minval(map%dst) < begg_o .or. maxval(map%dst) > endg_o) then
       write(iulog,*) 'map_maparrayl_rev ERROR: dst out of bounds:', &
                   minval(map%dst),maxval(map%dst),begg_o,endg_o
       call endrun()
    endif

    ! initialize field on output grid to zero everywhere

    fld_o(:,:) = 0._r8

    ! map flds

    do n = 1,map%nwts
       g_i = map%src(n) - begg_i + 1
       g_o = map%dst(n) - begg_o + 1
       wtx = map%S(n)
       do ifld = 1,nflds
          fld_o(ifld,g_o) = fld_o(ifld,g_o) + wtx * fld_i(ifld,g_i)
       enddo
    enddo

end subroutine map_maparrayl_rev

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_maparrayl_av
!
! !INTERFACE:
  subroutine map_maparrayl_av(begg_i, endg_i, begg_o, endg_o, nflds, &
                              fld_i, fld_o, map)
!
! !DESCRIPTION:
! This subroutine maps arrays, local 1d with 1d map type
!
! !USES:
  use mct_mod
!
! !ARGUMENTS:
  implicit none
  integer                         :: begg_i,endg_i         !beg,end of input grid
  integer                         :: begg_o,endg_o         !beg,end of output grid
  integer                         :: nflds                 !number of fields being mapped
  type(mct_aVect), intent(inout)  :: fld_i                 !land model import and export states
  type(mct_aVect), intent(inout)  :: fld_o                 !land model import and export states
  type(map_type), intent(in)      :: map
!
! !REVISION HISTORY:
! 2005.11.15  T Craig  Creation.
! 2006.3.30   P Worley Restructuring for improved vector performance
!

! !LOCAL VARIABLES:
!EOP
  integer :: g_o ,g_i              ! gridcell indices
  integer :: n,ifld                ! loop counters
  real(r8):: wtx                   ! wt for map
!
!------------------------------------------------------------------------------

    if (trim(map%type) /= trim(map_typelocal)) then
       write(iulog,*) 'map_maparrayl_av WARNING: map type not correct, ', &
                   map%name,map%type
    endif

    ! check for errors in overlap

    if (minval(map%src) < begg_i .or. maxval(map%src) > endg_i) then
       write(iulog,*) 'map_maparrayl_av ERROR: src out of bounds:', &
                   minval(map%src),maxval(map%src),begg_i,endg_i
       call endrun()
    endif

    if (minval(map%dst) < begg_o .or. maxval(map%dst) > endg_o) then
       write(iulog,*) 'map_maparrayl_av ERROR: dst out of bounds:', &
                   minval(map%dst),maxval(map%dst),begg_o,endg_o
       call endrun()
    endif

    ! initialize field on output grid to zero everywhere

    call mct_aVect_zero(fld_o)

    ! map flds

    do n = 1,map%nwts
       g_i = map%src(n) - begg_i + 1
       g_o = map%dst(n) - begg_o + 1
       wtx = map%S(n)
       do ifld = 1,nflds
          fld_o%rAttr(ifld,g_o) = fld_o%rAttr(ifld,g_o) + wtx * fld_i%rAttr(ifld,g_i)
       enddo
    enddo

end subroutine map_maparrayl_av

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_maparrayg
!
! !INTERFACE:
  subroutine map_maparrayg(fld_i, fld_o, map)
!
! !DESCRIPTION:
! This subroutine maps arrays, global 1d with 1d map type
!
! !USES:
!
! !ARGUMENTS:
  implicit none
  real(r8), intent(in)            :: fld_i(:,:)
  real(r8), intent(out)           :: fld_o(:,:)
  type(map_type), intent(in)      :: map
!
! !REVISION HISTORY:
! 2005.11.15  T Craig  Creation.
! 2006.3.30   P Worley Restructuring for improved vector performance
!

! !LOCAL VARIABLES:
!EOP
  integer :: g_o ,g_i              ! gridcell indices
  integer :: n,ifld,nflds          ! loop counters
  real(r8):: wtx                   ! wt for map

!
!------------------------------------------------------------------------------

    if (trim(map%type) /= trim(map_typeglobal)) then
       write(iulog,*) 'map_maparrayg WARNING: map type not correct, ', &
                   map%name,map%type
    endif

    nflds = size(fld_i,dim=2)
    ! initialize field on output grid to zero everywhere

    fld_o(:,:) = 0._r8

    ! map flds

    do ifld = 1,nflds
    do n = 1,map%nwts
       g_i = map%src(n)
       g_o = map%dst(n)
       wtx = map%S(n)
       fld_o(g_o,ifld) = fld_o(g_o,ifld) + wtx * fld_i(g_i,ifld)
    enddo
    enddo

end subroutine map_maparrayg

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_setgatm
!
! !INTERFACE:
  subroutine map_setgatmFM(gatm, alatlon, llatlon, amask, pftm)
!
! !DESCRIPTION:
! Set gatm index in ldomain.  unique for finemesh.
!
! !USES:
!
  use domainMod, only : latlon_type
!
! !ARGUMENTS:
  implicit none
  integer, pointer               :: gatm(:)
  type(latlon_type), intent(in)  :: alatlon
  type(latlon_type), intent(in)  :: llatlon
  integer          , intent(in)  :: amask(:)
  integer          , intent(in)  :: pftm(:)
!
! !REVISION HISTORY:
! 2006.06.28  T Craig  Creation.
! 2006.08.23  P Worley Performance optimizations
!
!
! !LOCAL VARIABLES:
!EOP
    integer          :: nlon_a       !input  grid: max number of longitude pts
    integer          :: nlat_a       !input  grid: number of latitude  points
    integer          :: ns_a         !input  grid: total number of cells
!    integer ,pointer :: mask_a(:)    !input grid: mask
    real(r8),pointer :: lon_a(:)     !input grid: longitude (degrees)
    real(r8),pointer :: lat_a(:)     !input grid: latitude  (degrees)
    real(r8),pointer :: lone_a(:)    !input grid: longitude, E edge (degrees)
    real(r8),pointer :: lonw_a(:)    !input grid: longitude, W edge (degrees)
    real(r8),pointer :: latn_a(:)    !input grid: latitude , N edge (degrees)
    real(r8),pointer :: lats_a(:)    !input grid: latitude , S edge (degrees)
    integer          :: nlon_l       !output grid: max number of longitude pts
    integer          :: nlat_l       !output grid: number of latitude  points
    integer          :: ns_l         !output grid: total number of cells
!    integer ,pointer :: gatm_l(:)    !output grid: atm grid cell overlapping
    real(r8),pointer :: lon_l(:)     !output grid: longitude (degrees)
    real(r8),pointer :: lat_l(:)     !output grid: latitude  (degrees)
    real(r8),pointer :: lone_l(:)    !output grid: longitude, E edge (degrees)
    real(r8),pointer :: lonw_l(:)    !output grid: longitude, W edge (degrees)
    real(r8),pointer :: latn_l(:)    !output grid: latitude , N edge (degrees)
    real(r8),pointer :: lats_l(:)    !output grid: latitude , S edge (degrees)
!    integer ,pointer :: pftm_l(:)    !output grid: cell frac
    integer          :: n            !loop counters
    integer          :: ia,ja,na     !indices for atm grid
    integer          :: il,jl,nl     !indices for lnd grid
    integer          :: if,jf,nf     !found indices
    integer          :: noffset=1    !noffset
    real(r8)         :: doffset =360_r8 !offset value
    real(r8)         :: offset       !offset*n_offset
    logical          :: found        !local logical 
    logical          :: overlapgrid  !are atm and lnd grids 1:1
    logical          :: latlongrid   !are atm and lnd grids regular lat/lon
    real(r8),parameter :: relerr = 1.0e-6    ! error limit
    integer          :: cnt,cnt0,cmax  !counters
!tcx fix, lat_l_local should be removed when limit no longer needed
    real(r8)         :: lat_l_local  !local copy of lat_l(il,jl), adjusted
    real(r8),parameter:: eps = 1.0e-8  ! eps for check

    integer, pointer :: ilfound(:)     ! list over overlap i indices
    integer, pointer :: jlfound(:)     ! list over overlap j indices
    integer, pointer :: nocnt(:)      ! lnd cell count per atm cell
    integer, pointer :: nooff(:)      ! atm cell offset in nomap
    integer, pointer :: nomap(:)      ! map from atm cell to lnd cells
    integer :: noidx
!
!------------------------------------------------------------------------

    !--- set pointers into domains, initialize gridmap a2l gridmap ---

    ns_a   =  alatlon%ns
    nlon_a =  alatlon%ni
    nlat_a =  alatlon%nj
    lat_a  => alatlon%latc
    lon_a  => alatlon%lonc
    latn_a => alatlon%latn
    lats_a => alatlon%lats
    lone_a => alatlon%lone
    lonw_a => alatlon%lonw

    ns_l   =  llatlon%ns
    nlon_l =  llatlon%ni
    nlat_l =  llatlon%nj
    lat_l  => llatlon%latc
    lon_l  => llatlon%lonc
    latn_l => llatlon%latn
    lats_l => llatlon%lats
    lone_l => llatlon%lone
    lonw_l => llatlon%lonw

    allocate(gatm(ns_l))

    !--- search for the overlap, input to output, 1:1 "disaggregation" ---
    !--- this is the coarse to fine map where there should be exactly
    !--- one coarse overlap point for each fine point 
    !--- three possible search algorithms, overlapgrid, latlongrid, neither
    !--- figure out which search scheme to use

    !--- overlapgrid means both grids are identical to relerr
    overlapgrid = .false.
    if (nlon_a == nlon_l .and. nlat_a == nlat_l) then
       overlapgrid = .true.
       do jl = 1, nlat_l
          if (abs( lat_l(jl)- lat_a(jl)) > relerr .or. &
              abs(latn_l(jl)-latn_a(jl)) > relerr .or. &
              abs(lats_l(jl)-lats_a(jl)) > relerr) then
             overlapgrid = .false.
          endif
       enddo
       do il = 1, nlon_l
          if (abs( lon_l(il)- lon_a(il)) > relerr .or. &
              abs(lone_l(il)-lone_a(il)) > relerr .or. &
              abs(lonw_l(il)-lonw_a(il)) > relerr) then
             overlapgrid = .false.
          endif
       enddo
    endif

    !--- latlongrid means atm grid is regular latlon grid
    !--- always true for now, may not be in the future
    latlongrid = .true.

    if (masterproc) write(iulog,*) 'setgatm overlapgrid,latlongrid = ',overlapgrid,latlongrid

!pw major restructuring follows
    if (overlapgrid) then
       gatm = -1
       do nl = 1,ns_l
          if (pftm(nl) >= 0) then       ! only real or fake points
             if (amask(nl) /= 0) then
                gatm(nl)=nl
             endif
          endif
       enddo

    elseif (latlongrid) then
!pw Still need restructuring for vectorization.

       allocate(ilfound(nlon_l),jlfound(nlat_l))
       ilfound = 0
       jlfound = 0

       do jl = 1, nlat_l
          found  = .false.
          lat_l_local = min(max(lat_l(jl),-90.0_r8),90.0_r8)  !limit [-90,90]
          do ja = 1,nlat_a
             if ((ja == 1 .and. lat_l_local <= latn_a(ja) .and.   &
                                lat_l_local >= lats_a(ja)) .or.   &
                 (ja >  1 .and. lat_l_local <= latn_a(ja) .and.   &
                                lat_l_local >  lats_a(ja))) then
                if (found) then
                   write(iulog,*) 'map_setgatm WARNING: found > 1 pt j ', &
                      jl,ja,jlfound(jl),lat_l(jl),lat_l_local
                   call endrun()
                endif
                jlfound(jl) = ja
                found = .true.
             endif
          enddo  ! ja
       enddo !jl

       do il = 1, nlon_l
          found  = .false.
          do ia = 1,nlon_a
          do n = -noffset,noffset
             offset = n*doffset
             if ((ia == 1 .and. lon_l(il)+offset <= lone_a(ia) .and.   &
                                lon_l(il)+offset >  lonw_a(ia)) .or.   &
                 (ia >  1 .and. lon_l(il)+offset <= lone_a(ia) .and.   &
                                lon_l(il)+offset >  lonw_a(ia))) then
                if (found) then
                   write(iulog,*) 'map_setgatm WARNING: found > 1 pt i ', &
                      il,ia,ilfound(il),lon_l(il)+offset
                   call endrun()
                endif
                ilfound(il) = ia
                found = .true.
             endif
          enddo  ! n, offset
          enddo  ! ia
       enddo

       gatm = -1
       do jl = 1,nlat_l
       do il = 1,nlon_l
          nl = (jl-1)*nlon_l + il
          if (pftm(nl) >= 0) then       ! only real or fake points
             if = ilfound(il)
             jf = jlfound(jl)
             nf = (jf-1)*nlon_a + if
             if (if == 0 .or. jf == 0) then
                write(iulog,*) 'map_setgatm ERROR: pt not found, ', &
                   il,lon_l(il), '_l', &
                   minval(lon_l),maxval(lon_l),           &
                   minval(lat_l),maxval(lat_l),'_awe',    &
                   minval(lonw_a),maxval(lonw_a),         &
                   minval(lone_a),maxval(lone_a),'_asn',  &
                   minval(lats_a),maxval(lats_a),         &
                   minval(latn_a),maxval(latn_a)
                call endrun()
             endif
             if (amask(nf) /= 0) then
                gatm(nl)=nf
             endif
          endif
       enddo
       enddo

       deallocate(ilfound,jlfound)

    else
!pw Still need restructuring for vectorization

       write(iulog,*) 'map_setgatm ERROR: irregular lat lon grid not supported'
       call endrun()

    endif
!pw major restructuring ends

    !--- remove fake land if there is at least one real land overlap point ---
    allocate(nocnt(ns_a),nooff(ns_a),nomap(ns_l))

    nocnt = 0
    do nl = 1,ns_l
       na = gatm(nl)
       if (na > 0) then
          nocnt(na) = nocnt(na) + 1
       endif
    enddo

    nooff(1) = 1
    do na = 2,ns_a
       nooff(na) = nooff(na-1) + nocnt(na-1)
    enddo

    nocnt = 0
    nomap = -1
    do nl = 1,ns_l
       na = gatm(nl)
       if (na > 0) then
         nomap(nooff(na)+nocnt(na)) = nl
         nocnt(na) = nocnt(na) + 1
       endif
    enddo

    do na = 1,ns_a
       found = .false.        ! check if any points are real land
       do noidx = 0,nocnt(na)-1
          nl = nomap(nooff(na)+noidx)
          if (pftm(nl) > 0 ) then
             found = .true.
          endif
       enddo
       if (found) then        ! if so, keep only real land points
          do noidx = 0,nocnt(na)-1
             nl = nomap(nooff(na)+noidx)
             if (pftm(nl) <= 0 ) then
                gatm(nl) = -1
             endif
       enddo
       endif
    enddo

    !--- check that valid fine grid points have coarse mapping gridpoints
    if (masterproc) then
       nocnt = 0
       do nl = 1,ns_l
          na = gatm(nl)
          if (na > 0) then
             nocnt(na) = nocnt(na) + 1
          endif
       enddo

       found = .true.
       do na = 1,ns_a
          if ((amask(na) /= 0) .and. (nocnt(na) == 0)) then
             found = .false.
          endif
       enddo
       if (.not. found) then
          do na = 1,ns_a
             if ((amask(na) /= 0) .and. (nocnt(na) == 0)) then
                write(iulog,*) 'map_setgatm ERROR: invalid f->c index ', &
                   na,amask(na)
                call endrun()
             endif
          enddo
       endif

    endif

    deallocate(nocnt,nooff,nomap)

end subroutine map_setgatmFM

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_setmapsFM
!
! !INTERFACE:
  subroutine map_setmapsFM(domain_a, domain_l, gatm_l, name)
!
! !DESCRIPTION:
! Set course to fine mesh maps and reverse.  domain_a should be coarse
! (atm) mesh, domain_l is fine (land) mesh.
! Simple overlap algorithm
!   - Find every fine gridcell within coarse gridcell
!   - Keep "real" cells unless there are none
!   - Weights based on areas of cells used, sum(weights)==1
!   - Use a2l mapping to set l2a mapping
!
! !USES:
  use decompMod, only : ldecomp, adecomp
  use decompMod, only : get_proc_bounds, get_proc_bounds_atm
  use clm_varsur, only : wtxy
!
! !ARGUMENTS:
  implicit none
  type(domain_type),intent(in)    :: domain_a
  type(domain_type),intent(in)    :: domain_l
  integer          ,intent(in)    :: gatm_l(:)
  character(len=*) ,optional,intent(in) :: name
!
! !REVISION HISTORY:
! 2005.12.01  T Craig  Creation.
!
!
! !LOCAL VARIABLES:
!EOP
    integer           :: ns_a
    integer           :: ns_l
    integer           :: na,nag
    integer           :: nl,nlg
    real(r8),pointer :: area_a(:)    !input grid: cell area
    real(r8),pointer :: topo_a(:)    !input grid: cell topo/elevation
    real(r8),pointer :: frac_a(:)    !input grid: cell frac
    real(r8),pointer :: area_l(:)    !output grid: cell area
    real(r8),pointer :: nara_l(:)    !output grid: cell equiv upscale area
    real(r8),pointer :: topo_l(:)    !output grid: cell topo/elevation
    real(r8),pointer :: ntop_l(:)    !output grid: cell equiv downscale topo
    real(r8),pointer :: frac_l(:)    !output grid: cell frac

    integer          :: n_a2l       !a2l nwts
    integer          :: mo_a2l      !a2l dstmo
    integer ,pointer :: src_a2l(:)  !a2l src indices
    integer ,pointer :: dst_a2l(:)  !a2l dst indices
    real(r8),pointer :: wts_a2l(:)  !a2l wts indices
    integer          :: n_l2a       !l2a nwts
    integer          :: mo_l2a      !l2a dstmo
    integer ,pointer :: src_l2a(:)  !l2a src indices
    integer ,pointer :: dst_l2a(:)  !l2a dst indices
    real(r8),pointer :: wts_l2a(:)  !l2a wts indices

    integer          :: na2l,nl2a
    integer          :: begg,endg
    integer          :: abegg,aendg
    integer ,pointer :: ncnta(:)    !number of overlap points in l2a
    integer ,pointer :: cnta(:)     !number of overlap points in l2a
    real(r8),pointer :: asum(:)     !number of overlap points in l2a
    real(r8),pointer :: tsum(:)     !number of overlap points in l2a
    real(r8),parameter:: eps = 1.0e-8  ! eps for check
!
!------------------------------------------------------------------------

    !--- get local beg/end

    call get_proc_bounds    ( begg,  endg)
    call get_proc_bounds_atm(abegg, aendg)

    !--- set pointers into domains, initialize gridmap a2l gridmap ---

    call domain_setptrs(domain_a,ns=ns_a,area=area_a,frac=frac_a,topo=topo_a)
    call domain_setptrs(domain_l,ns=ns_l,area=area_l,frac=frac_l,topo=topo_l, &
                        nara=nara_l,ntop=ntop_l)

    !--- allocate temporaries

    allocate(ncnta(abegg:aendg),asum(abegg:aendg),tsum(abegg:aendg),cnta(abegg:aendg))

    !--- count local number of wts needed, na2l, nl2a
    !--- also set global arrays ncnta and asum for later computations

    na2l = 0
    nl2a = 0
    ncnta = 0
    asum = 0.0_r8
    do nl = begg,endg
       nlg = ldecomp%gdc2glo(nl)
       nag = gatm_l(nlg)
       if (nag <= 0) then
          write(iulog,*) 'map_setmapsFM ERROR nag0 <= 0',nl,nlg,nag
          call endrun()
       endif
       na  = adecomp%glo2gdc(nag)
       if (nlg < 1 .or. nlg > ns_l .or. &
           nag < 1 .or. nag > ns_a .or. &
           na  < abegg .or. na > aendg) then
          write(iulog,*) 'map_setmapsFM ERROR in index0 ',nl,nlg,nag,na,ns_l,ns_a,abegg,aendg
          call endrun()
       endif
       ncnta(na) = ncnta(na) + 1
       asum(na) = asum(na) + area_l(nl)
       na2l = na2l + 1
       nl2a = nl2a + 1
    enddo

    !--- initialize and allocate maps
    call map_init(map1dl_a2l,domain_a,domain_l,na2l,name='setmapsFM_a2l', &
                  type=map_typelocal)
    call map_init(map1dl_l2a,domain_l,domain_a,nl2a,name='setmapsFM_l2a', &
                  type=map_typelocal)

    !--- set pointers to the maps

    call map_setptrs(map1dl_a2l,src=src_a2l,dst=dst_a2l,wts=wts_a2l, &
                     nwts=n_a2l,dstmo=mo_a2l)
    call map_setptrs(map1dl_l2a,src=src_l2a,dst=dst_l2a,wts=wts_l2a, &
                     nwts=n_l2a,dstmo=mo_l2a)

    !--- set dstmo in maps

    map1dl_a2l%dstmo = 1
    map1dl_l2a%dstmo = maxval(ncnta)

    !--- set src,dst,wts in maps

    na2l = 0
    nl2a = 0
    cnta = 0
    do nl = begg,endg
       nlg = ldecomp%gdc2glo(nl)
       nag = gatm_l(nlg)
       if (nag <= 0) then
          write(iulog,*) 'map_setmapsFM ERROR nag1 <= 0',nl,nlg,nag
          call endrun()
       endif
       na  = adecomp%glo2gdc(nag)
       if (nlg < 1 .or. nlg > ns_l .or. &
           nag < 1 .or. nag > ns_a .or. &
           na  < abegg .or. na > aendg) then
          write(iulog,*) 'map_setmapsFM ERROR in index1 ',nl,nlg,nag,na,ns_l,ns_a,abegg,aendg
          call endrun()
       endif

       na2l = na2l + 1
       if (na2l > n_a2l) then
          write(iulog,*) 'map_setmapsFM ERROR na2l > n_a2l ',na2l,n_a2l
          call endrun()
       endif
       src_a2l(na2l) = na
       dst_a2l(na2l) = nl
       wts_a2l(na2l) = 1.0_r8

       nl2a = nl2a + 1
       cnta(na) = cnta(na) + 1
       if (nl2a > n_l2a .or. cnta(na) > ncnta(na)) then
          write(iulog,*) 'map_setmapsFM ERROR nl2a > n_l2a ',nl2a,n_l2a,cnta(na),ncnta(na)
          call endrun()
       endif
       src_l2a(nl2a) = nl
       dst_l2a(nl2a) = na
       if (ncnta(na) == 1) then
          wts_l2a(nl2a) = 1.0_r8
       else
          if (asum(na) <= 0.0_r8) then
             write(iulog,*) 'map_setmapsFM ERROR asum <= 0 ',nlg,nag,asum(na)
             call endrun()
          endif
          wts_l2a(nl2a) = (area_l(nl)/asum(na))
       endif
    enddo

    !--- update ldomain arrays based on upscale/downscale stuff

    tsum = 0.0_r8
    do nl = begg,endg
       nlg = ldecomp%gdc2glo(nl)
       nag = gatm_l(nlg)
       if (nag <= 0) then
          write(iulog,*) 'map_setmapsFM ERROR nag2 <= 0',nl,nlg,nag
          call endrun()
       endif
       na  = adecomp%glo2gdc(nag)
       if (nlg < 1 .or. nlg > ns_l .or. &
           nag < 1 .or. nag > ns_a .or. &
           na  < abegg .or. na > aendg) then
          write(iulog,*) 'map_setmapsFM ERROR in index2 ',nl,nlg,nag,na,ns_l,ns_a,abegg,aendg
          call endrun()
       endif
       !??? MV: Note that frac_l is not set yet - so it is -1.0e36 - so first if will never be
       ! exercised - seems like this is a bug???
       if (frac_l(nl) > 0.0_r8) then
          ntop_l(nl) = topo_l(nl)                      ! set topo ovr lnd
       else
          ntop_l(nl) = max(0.0_r8,topo_l(nl))          ! set topo ovr ocn 
       endif
       tsum(na) = tsum(na) + ntop_l(nl)*area_l(nl)   ! area wt topo avg
    enddo

    do nl = begg,endg
       nlg = ldecomp%gdc2glo(nl)
       nag = gatm_l(nlg)
       if (nag <= 0) then
          write(iulog,*) 'map_setmapsFM ERROR nag3 <= 0',nl,nlg,nag
          call endrun()
       endif
       na  = adecomp%glo2gdc(nag)
       if (nlg < 1 .or. nlg > ns_l .or. &
           nag < 1 .or. nag > ns_a .or. &
           na  < abegg .or. na > aendg) then
          write(iulog,*) 'map_setmapsFM ERROR in index3 ',nl,nlg,nag,na,ns_l,ns_a,abegg,aendg
          call endrun()
       endif
       if (ncnta(na) == 1) then
          nara_l(nl) = area_a(na)
          ntop_l(nl) = topo_a(na)
       else
          if (asum(na) <= 0.0_r8) then
             write(iulog,*) 'map_setmapsFM ERROR2 asum <= 0 ',nl,na,asum(na)
             call endrun()
          endif
          if (tsum(na) <  0.0_r8) then
             write(iulog,*) 'map_setmapsFM ERROR2 tsum < 0 ',nl,na,tsum(na)
             call endrun()
          endif
          nara_l(nl)   = (area_l(nl)/asum(na))*area_a(na)
! DOWNSCALING
!-v-v-v-v-v- land topo elevation adjustment for downscaling -v-v-v-v-
!----------- want avg land topo to be equal to atm topo and want the
!----------- variability in finemesh land topo to be preserved
          if (tsum(na) == 0.) then
             ntop_l(nl)   = ntop_l(nl)+topo_a(na)
          elseif (topo_a(na) > 0.) then
             ntop_l(nl)   = (ntop_l(nl)/(tsum(na)/asum(na)))*topo_a(na)
          else
             ntop_l(nl)   = (ntop_l(nl)-(tsum(na)/asum(na)))+topo_a(na)
          endif
!-^-^-^-^-^- land topo elevation adjustment for downscaling -^-^-^-^-
       endif
    enddo

    !--- check that areas and topos match up

    asum = 0.0_r8
    tsum = 0.0_r8
    do nl = begg,endg
       nlg = ldecomp%gdc2glo(nl)
       nag = gatm_l(nlg)
       if (nag <= 0) then
          write(iulog,*) 'map_setmapsFM ERROR nag4 <= 0',nl,nlg,nag
          call endrun()
       endif
       na  = adecomp%glo2gdc(nag)
       if (nlg < 1 .or. nlg > ns_l .or. &
           nag < 1 .or. nag > ns_a .or. &
           na  < abegg .or. na > aendg) then
          write(iulog,*) 'map_setmapsFM ERROR in index4 ',nl,nlg,nag,na,ns_l,ns_a,abegg,aendg
          call endrun()
       endif
       asum(na) = asum(na) + nara_l(nl)
       tsum(na) = tsum(na) + nara_l(nl)*ntop_l(nl)
    enddo

    do na = abegg, aendg
       if (asum(na) > 0.0_r8) then
          if (abs(asum(na)-area_a(na)) > eps .or. &
              abs(tsum(na)/area_a(na) - topo_a(na)) > eps) then
                 write(iulog,*) ' map_setmapsFM ERROR: ERROR in nara,ntop, ',asum(na),area_a(na),tsum(na),topo_a(na),eps
                 call endrun()
          endif
       endif
    enddo

    !--- clean up
    deallocate(ncnta,cnta,asum,tsum)

    call map_checkmap(map1dl_a2l)
    call map_checkmap(map1dl_l2a)

end subroutine map_setmapsFM

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_checkmap
!
! !INTERFACE:
  subroutine map_checkmap(map)
!
! !DESCRIPTION:
! Checks the map for consistency
!
! !USES:
  use spmdMod
  use decompMod, only : get_proc_bounds, get_proc_bounds_atm
!
! !ARGUMENTS:
  implicit none
  type(map_type),intent(in) :: map
!
! !REVISION HISTORY:
! 2005.12.01  T Craig  Creation.
!
!
! !LOCAL VARIABLES:
!EOP
    integer          :: nlon_i       !input  grid: max number of longitude pts
    integer          :: nlat_i       !input  grid: number of latitude  points
    integer          :: nlon_o       !output grid: max number of longitude pts
    integer          :: nlat_o       !output grid: number of latitude  points
    integer          :: nwts         !local num of weights
    integer          :: dstmo        !max num of overlapping cells
    integer ,pointer :: src(:)       !src index
    integer ,pointer :: dst(:)       !dst index
    real(r8),pointer :: wts(:)       !weight
    integer          :: n            !loop counters
    real(r8),pointer :: rsum(:)      !local array for deriving values
    real(r8)         :: rmin,rmax    !local min/max values
    real(r8)         :: smin,smax    !local min/max values
    integer          :: imin,imax    !local min/max values
    integer          :: nmin,nmax,nsum  !local min/max/sum values
    integer          :: ier          ! error flag
!
!------------------------------------------------------------------------

    !--- get general info ---
    call map_setptrs(map,ni_i=nlon_i,nj_i=nlat_i, &
                         ni_o=nlon_o,nj_o=nlat_o)
    call map_setptrs(map,nwts=nwts,src=src,dst=dst,wts=wts,dstmo=dstmo)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'map_checkmap name          = ',trim(map%name)
       write(iulog,*) 'map_checkmap type          = ',trim(map%type)
       write(iulog,*) 'map_checkmap src grid      = ',nlon_i,nlat_i
       write(iulog,*) 'map_checkmap dst grid      = ',nlon_o,nlat_o
       write(iulog,*) 'map_checkmap dstmo         = ',dstmo
    endif

    call mpi_reduce(nwts,nmin,1,MPI_INTEGER,MPI_MIN,0,mpicom,ier)
    if (ier /= 0) then
       write(iulog,*)'mpi_reduce nwts error: ',ier
    endif
    call mpi_reduce(nwts,nmax,1,MPI_INTEGER,MPI_MAX,0,mpicom,ier)
    if (ier /= 0) then
       write(iulog,*)'mpi_reduce nwts error: ',ier
    endif
    call mpi_reduce(nwts,nsum,1,MPI_INTEGER,MPI_SUM,0,mpicom,ier)
    if (ier /= 0) then
       write(iulog,*)'mpi_reduce nwts error: ',ier
    endif
    if (masterproc) then
       write(iulog,*) 'map_checkmap nwts gmin/max = ',nmin,nmax
       write(iulog,*) 'map_checkmap nwts gsum     = ',nsum
    endif

    imin = minval(src)
    imax = maxval(src)
    call mpi_reduce(imin,nmin,1,MPI_INTEGER,MPI_MIN,0,mpicom,ier)
    if (ier /= 0) then
       write(iulog,*)'mpi_reduce src min error: ',ier
    endif
    call mpi_reduce(imax,nmax,1,MPI_INTEGER,MPI_MAX,0,mpicom,ier)
    if (ier /= 0) then
       write(iulog,*)'mpi_reduce src max error: ',ier
    endif
    if (masterproc) then
       write(iulog,*) 'map_checkmap src gmin/max  = ',nmin,nmax
    endif

    imin = minval(dst)
    imax = maxval(dst)
    call mpi_reduce(imin,nmin,1,MPI_INTEGER,MPI_MIN,0,mpicom,ier)
    if (ier /= 0) then
       write(iulog,*)'mpi_reduce dst min error: ',ier
    endif
    call mpi_reduce(imax,nmax,1,MPI_INTEGER,MPI_MAX,0,mpicom,ier)
    if (ier /= 0) then
       write(iulog,*)'mpi_reduce dst max error: ',ier
    endif
    if (masterproc) then
       write(iulog,*) 'map_checkmap dst gmin/max  = ',nmin,nmax
    endif

    rmin = minval(wts)
    rmax = maxval(wts)
    call mpi_reduce(rmin,smin,1,MPI_REAL8,MPI_MIN,0,mpicom,ier)
    if (ier /= 0) then
       write(iulog,*)'mpi_reduce wts min error: ',ier
    endif
    call mpi_reduce(rmax,smax,1,MPI_REAL8,MPI_MAX,0,mpicom,ier)
    if (ier /= 0) then
       write(iulog,*)'mpi_reduce wts max error: ',ier
    endif
    if (masterproc) then
       write(iulog,*) 'map_checkmap wts gmin/max  = ',smin,smax
    endif

    imin = minval(dst)
    imax = maxval(dst)
    allocate(rsum(imin:imax))
    rsum = 0.0_r8
    do n = 1,nwts
       if (dst(n) < imin .or. dst(n) > imax) then
          write(iulog,*) 'map_checkmap dst index error ',n,dst(n),imin,imax
          call endrun()
       endif
       rsum(dst(n)) = rsum(dst(n)) + wts(n)
    enddo
    rmin = minval(rsum)
    rmax = maxval(rsum)
    deallocate(rsum)
    call mpi_reduce(rmin,smin,1,MPI_REAL8,MPI_MIN,0,mpicom,ier)
    if (ier /= 0) then
       write(iulog,*)'mpi_reduce swts min error: ',ier
    endif
    call mpi_reduce(rmax,smax,1,MPI_REAL8,MPI_MAX,0,mpicom,ier)
    if (ier /= 0) then
       write(iulog,*)'mpi_reduce swts max error: ',ier
    endif
    if (masterproc) then
       write(iulog,*) 'map_checkmap swts gmin/max  = ',smin,smax
    endif

end subroutine map_checkmap

end module downscaleMod



















