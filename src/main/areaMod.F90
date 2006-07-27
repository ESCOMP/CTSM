#include <misc.h>
#include <preproc.h>

module areaMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: areaMod
!
! !DESCRIPTION:
! Area averaging routines
! Thes routines are used for area-average mapping of a field from one
! grid to another.
!
! !USES:
  use clm_varcon   , only : re
  use clm_varpar   , only : numrad
  use domainMod    , only : domain_type, domain_setptrs
  use shr_const_mod, only : SHR_CONST_PI
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_sys_mod  , only : shr_sys_flush
  use spmdMod      , only : masterproc
  use nanMod      
  use abortutils   , only : endrun
!
! !PUBLIC TYPES:
  implicit none
  private

  type gridmap_type
     private
     ! lower level in hierarchy
     character(len=32)         :: name
     character(len=16)         :: type        ! global, dst, src, etc
     type(domain_type),pointer :: domain_i    ! domain_i
     type(domain_type),pointer :: domain_o    ! domain_o
     integer                   :: mx_ovr       ! max num of overlapping cells
     integer          ,pointer :: n_ovr(:,:)   ! number of overlapping cells
     integer          ,pointer :: i_ovr(:,:,:) ! i index of overlap input cell
     integer          ,pointer :: j_ovr(:,:,:) ! j index of overlap input cell
     real(r8)         ,pointer :: w_ovr(:,:,:) ! wt of overlap input cell
  end type gridmap_type
  public gridmap_type

  type map_type
     private
     ! lower level in hierarchy
     character(len=32)         :: name
     character(len=16)         :: type        ! global, dst, src, etc
     type(domain_type),pointer :: domain_i    ! domain_i
     type(domain_type),pointer :: domain_o    ! domain_o
     integer                   :: nwts          ! size of row, col, S (local)
     integer          ,pointer :: src(:)      ! src index (COL)
     integer          ,pointer :: dst(:)      ! dst index (ROW)
     real(r8)         ,pointer :: S(:)        ! wt of overlap input cell
     integer                   :: dstmo       ! max num of overlaps of dst cell
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
  public :: map_maparray
  public :: map_setgatm
  public :: gridmap_checkmap
  public :: areaini        ! area averaging initialization
  public :: areaave        ! area averaging of field from input to output grids
  interface celledge
     module procedure celledge_regional
     module procedure celledge_global  
  end interface
  interface cellarea
     module procedure cellarea_regional
     module procedure cellarea_global
  end interface
  public :: celledge
  public :: cellarea
!
! !REVISION HISTORY:
! Created by Sam Levis
! Updated to clm2.1 data structures by Mariana Vertenstein
! 2005.11.01 Updated and cleaned by T Craig
!
!EOP
!
! PRIVATE MEMBER FUNCTIONS:
  private :: areamap   ! weights and indices for area of overlap between grids
  private :: areaovr   ! area of overlap between grid cells
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gridmap_init
!
! !INTERFACE:
  subroutine gridmap_init(gridmap,domain_i,domain_o,mwts,name,type)
!
! !DESCRIPTION:
! This subroutine initializes the gridmap datatype
!
! !USES:
!
! !ARGUMENTS:
  implicit none
  type(gridmap_type), intent(inout)       :: gridmap
  type(domain_type) , intent(in),target   :: domain_i
  type(domain_type) , intent(in),target   :: domain_o
  integer           , intent(in)          :: mwts     ! max number of wts
  character(len=*)  , intent(in),optional :: name
  character(len=*)  , intent(in),optional :: type
!
! !REVISION HISTORY:
! 2005.11.15  T Craig  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
  integer ni,nj  ! size of domain_o
  integer ier    ! error flag
!------------------------------------------------------------------------------

  gridmap%domain_i => domain_i
  gridmap%domain_o => domain_o
  if (present(name)) then
    gridmap%name = trim(name)
  else
    gridmap%name = 'unset'
  endif
  if (present(type)) then
    gridmap%type = trim(type)
  else
    gridmap%type = 'unset'
  endif

  ni = domain_o%ni
  nj = domain_o%nj
  gridmap%mx_ovr = mwts
  allocate(gridmap%n_ovr(ni,nj)     , gridmap%i_ovr(ni,nj,mwts), &
           gridmap%j_ovr(ni,nj,mwts), gridmap%w_ovr(ni,nj,mwts),stat=ier)
  if (ier /= 0) then
     write(6,*) 'gridmap_init ERROR: allocate gridmap'
     call endrun()
  endif

  gridmap%n_ovr = bigint
  gridmap%i_ovr = -1
  gridmap%j_ovr = -1
  gridmap%w_ovr = nan

end subroutine gridmap_init
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
!EOP
!
! !LOCAL VARIABLES:
  integer ier    ! error flag
!------------------------------------------------------------------------------

  map%domain_i => domain_i
  map%domain_o => domain_o
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
     write(6,*) 'map_init ERROR: allocate map'
     call endrun()
  endif

  map%src = -1
  map%dst = -1
  map%S   = nan

end subroutine map_init
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gridmap_setptrs
!
! !INTERFACE:
  subroutine gridmap_setptrs(gridmap,name,type,domain_i,domain_o, &
     mx_ovr,n_ovr,i_ovr,j_ovr,w_ovr)
!
! !DESCRIPTION:
! This subroutine sets external pointer arrays to arrays in gridmap
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(gridmap_type),intent(in)       :: gridmap
    character(len=*) ,optional          :: name     
    character(len=*) ,optional          :: type
    type(domain_type),optional,pointer  :: domain_i
    type(domain_type),optional,pointer  :: domain_o
    integer          ,optional          :: mx_ovr
    integer          ,optional,pointer  :: n_ovr(:,:)
    integer          ,optional,pointer  :: i_ovr(:,:,:)
    integer          ,optional,pointer  :: j_ovr(:,:,:)
    real(r8)         ,optional,pointer  :: w_ovr(:,:,:)
!
! !REVISION HISTORY:
!   Created by T Craig
!
!EOP
!
! LOCAL VARIABLES:
!
!------------------------------------------------------------------------------
    if (present(name)) then
      name = gridmap%name
    endif
    if (present(type)) then
      type = gridmap%type
    endif
    if (present(domain_i)) then
      domain_i => gridmap%domain_i
    endif
    if (present(domain_o)) then
      domain_o => gridmap%domain_o
    endif
    if (present(mx_ovr)) then
      mx_ovr = gridmap%mx_ovr
    endif
    if (present(n_ovr)) then
      n_ovr => gridmap%n_ovr
    endif
    if (present(i_ovr)) then
      i_ovr => gridmap%i_ovr
    endif
    if (present(j_ovr)) then
      j_ovr => gridmap%j_ovr
    endif
    if (present(w_ovr)) then
      w_ovr => gridmap%w_ovr
    endif

end subroutine gridmap_setptrs
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_setptrs
!
! !INTERFACE:
  subroutine map_setptrs(map,name,type,domain_i,domain_o, &
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
    type(domain_type),optional,pointer  :: domain_i
    type(domain_type),optional,pointer  :: domain_o
    integer          ,optional          :: nwts
    integer          ,optional,pointer  :: src(:)
    integer          ,optional,pointer  :: dst(:)
    real(r8)         ,optional,pointer  :: wts(:)
    integer          ,optional          :: dstmo
!
! !REVISION HISTORY:
!   Created by T Craig
!
!EOP
!
! LOCAL VARIABLES:
!
!------------------------------------------------------------------------------
    if (present(name)) then
      name = map%name
    endif
    if (present(type)) then
      type = map%type
    endif
    if (present(domain_i)) then
      domain_i => map%domain_i
    endif
    if (present(domain_o)) then
      domain_o => map%domain_o
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
! !IROUTINE: map_maparray
!
! !INTERFACE:
  subroutine map_maparray(begg_i, endg_i, begg_o, endg_o, nflds, &
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
!EOP

! !LOCAL VARIABLES:
  integer :: g_o ,g_i              ! gridcell indices
  integer :: n,ifld                ! loop counters
  real(r8):: wtx                   ! wt for map

!
!------------------------------------------------------------------------------

    if (trim(map%type) /= trim(map_typelocal)) then
       write(6,*) 'map_maparray WARNING: map type not correct, ', &
                   map%name,map%type
    endif

    ! check for errors in overlap

    if (minval(map%src) < begg_i .or. maxval(map%src) > endg_i) then
       write(6,*) 'map_maparray ERROR: src out of bounds:', &
                   minval(map%src),maxval(map%src),begg_i,endg_i
       call endrun()
    endif

    if (minval(map%dst) < begg_o .or. maxval(map%dst) > endg_o) then
       write(6,*) 'map_maparray ERROR: dst out of bounds:', &
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

end subroutine map_maparray

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_setgatm
!
! !INTERFACE:
  subroutine map_setgatm(domain_a, domain_l)
!
! !DESCRIPTION:
! Set gatm index in ldomain.  unique for finemesh.
!
! !USES:
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(in)    :: domain_a
  type(domain_type), intent(in)    :: domain_l
!
! !REVISION HISTORY:
! 2006.06.28  T Craig  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
    integer          :: nlon_a       !input  grid: max number of longitude pts
    integer          :: nlat_a       !input  grid: number of latitude  points
    integer          :: ns_a         !input  grid: total number of cells
    integer ,pointer :: mask_a(:)    !input grid: mask
    real(r8),pointer :: lon_a(:)     !input grid: longitude (degrees)
    real(r8),pointer :: lat_a(:)     !input grid: latitude  (degrees)
    real(r8),pointer :: lone_a(:)    !input grid: longitude, E edge (degrees)
    real(r8),pointer :: lonw_a(:)    !input grid: longitude, W edge (degrees)
    real(r8),pointer :: latn_a(:)    !input grid: latitude , N edge (degrees)
    real(r8),pointer :: lats_a(:)    !input grid: latitude , S edge (degrees)
    integer          :: nlon_l       !output grid: max number of longitude pts
    integer          :: nlat_l       !output grid: number of latitude  points
    integer          :: ns_l         !output grid: total number of cells
    integer ,pointer :: gatm_l(:)    !output grid: atm grid cell overlapping
    real(r8),pointer :: lon_l(:)     !output grid: longitude (degrees)
    real(r8),pointer :: lat_l(:)     !output grid: latitude  (degrees)
    real(r8),pointer :: lone_l(:)    !output grid: longitude, E edge (degrees)
    real(r8),pointer :: lonw_l(:)    !output grid: longitude, W edge (degrees)
    real(r8),pointer :: latn_l(:)    !output grid: latitude , N edge (degrees)
    real(r8),pointer :: lats_l(:)    !output grid: latitude , S edge (degrees)
    integer ,pointer :: pftm_l(:)    !output grid: cell frac
    integer          :: n            !loop counters
    integer          :: ii,ji,ni     !indices for input grid
    integer          :: io,jo,no     !indices for output grid
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
    real(r8)         :: lat_l_local  !local copy of lat_l(io,jo), adjusted
    real(r8),parameter:: eps = 1.0e-8  ! eps for check
!
!------------------------------------------------------------------------

    !--- set pointers into domains, initialize gridmap a2l gridmap ---

    call domain_setptrs(domain_a,ns=ns_a,ni=nlon_a,nj=nlat_a, &
       latc=lat_a,lonc=lon_a,mask=mask_a,  &
       latn=latn_a,lats=lats_a,lone=lone_a,lonw=lonw_a)
    call domain_setptrs(domain_l,ns=ns_l,ni=nlon_l,nj=nlat_l, &
       pftm=pftm_l,latc=lat_l,lonc=lon_l, &
       gatm=gatm_l,latn=latn_l,lats=lats_l,lone=lone_l,lonw=lonw_l)

    !--- search for the overlap, input to output, 1:1 "disaggregation" ---
    !--- this is the coarse to fine map where there should be exactly
    !--- one coarse overlap point for each fine point 
    !--- three possible search algorithms, overlapgrid, latlongrid, neither
    !--- figure out which search scheme to use

    !--- overlapgrid means both grids are identical
    overlapgrid = .false.
    if (nlon_a == nlon_l .and. nlat_a == nlat_l) then
       overlapgrid = .true.
       do no = 1, nlat_l*nlon_l
          if (abs( lat_l(no)- lat_a(no)) > relerr .or. &
              abs( lon_l(no)- lon_a(no)) > relerr .or. &
              abs(lone_l(no)-lone_a(no)) > relerr .or. &
              abs(lonw_l(no)-lonw_a(no)) > relerr .or. &
              abs(latn_l(no)-latn_a(no)) > relerr .or. &
              abs(lats_l(no)-lats_a(no)) > relerr) then
             overlapgrid = .false.
          endif
       enddo
    endif

    !--- latlongrid means atm grid is regular latlon grid
    !--- assume true and then set false if not
    latlongrid = .true.
    do ni = 1,nlat_a*nlon_a
       ii = mod((ni-1),nlon_a)+1
       ji = (int((ni-1)/nlon_a))*nlon_a + 1
       if (abs( lat_a(ni) -  lat_a(ji)) > relerr .or. &
           abs(latn_a(ni) - latn_a(ji)) > relerr .or. &
           abs(lats_a(ni) - lats_a(ji)) > relerr) then
           latlongrid = .false.
       endif
       if (abs( lon_a(ni) -  lon_a(ii)) > relerr .or. &
           abs(lone_a(ni) - lone_a(ii)) > relerr .or. &
           abs(lonw_a(ni) - lonw_a(ii)) > relerr) then
           latlongrid = .false.
       endif
    enddo

    if (masterproc) write(6,*) 'setgatm overlapgrid,latlongrid = ',overlapgrid,latlongrid

    do jo = 1, nlat_l
    do io = 1, nlon_l
    no = (jo-1)*nlon_l + io
    gatm_l(no) = -1
    if (pftm_l(no) >= 0) then       ! only real or fake points
       found  = .false.
       lat_l_local = min(max(lat_l(no),-90.0_r8),90.0_r8)  !limit [-90,90]

       if (overlapgrid) then
          found = .true.
          if = io
          jf = jo
       elseif (latlongrid) then
          do ji = 1,nlat_a
             ni = (ji-1)*nlon_a + 1
             if ((ji == 1 .and. lat_l_local <= latn_a(ni) .and.   &
                                lat_l_local >= lats_a(ni)) .or.   &
                 (ji >  1 .and. lat_l_local <= latn_a(ni) .and.   &
                                lat_l_local >  lats_a(ni))) then
                if (found) then
                   write(6,*) 'map_setgatm WARNING: found > 1 pt j ', &
                      io,jo,no,ji,ni,jf,lon_l(no),lat_l(no),lat_l_local
                   call endrun()
                endif
                jf = ji
                found = .true.
             endif
          enddo  ! ji
          if (found) then     ! move on to i
             found = .false.
          else                ! stop
             write(6,*) 'map_setgatm ERROR: pt not found, ', &
                io,jo,lon_l(no),lat_l(no),lat_l_local, '_o', &
                minval(lon_l),maxval(lon_l),           &
                minval(lat_l),maxval(lat_l),'_iwe',    &
                minval(lonw_a),maxval(lonw_a),         &
                minval(lone_a),maxval(lone_a),'_isn',  &
                minval(lats_a),maxval(lats_a),         &
                minval(latn_a),maxval(latn_a)
             call endrun()
          endif
          do ii = 1,nlon_a
          ni = (jf-1)*nlon_a + ii
          do n = -noffset,noffset
             offset = n*doffset
             if (lon_l(no)+offset <= lone_a(ni) .and.   &
                 lon_l(no)+offset >= lonw_a(ni)) then
                if (found) then
                   write(6,*) 'map_setgatm WARNING: found > 1 pt i ', &
                      io,jo,no,lon_l(no),lat_l(no),lat_l_local
                   call endrun()
                endif
                if = ii
                found = .true.
             endif
          enddo  ! n, offset
          enddo  ! ii
       else
          do n = -noffset,noffset
          offset = n*doffset
          do ji = 1,nlat_a
          do ii = 1,nlon_a
             ni = (ji-1)*nlon_a + ii
             if (lon_l(no)+offset <= lone_a(ni) .and.   &
                 lon_l(no)+offset >= lonw_a(ni) .and.   &
                 lat_l_local      <= latn_a(ni) .and.   &
                 lat_l_local      >= lats_a(ni)) then
                if (found) then
                   write(6,*) 'map_setgatm WARNING: found > 1 pt', &
                      io,jo,no,lon_l(no),lat_l(no),lat_l_local
                   call endrun()
                endif
                found = .true.
                if = ii
                jf = ji
             endif
          enddo  ! ii
          enddo  ! ji
          enddo  ! n, offset
       endif

       nf = (jf-1)*nlon_a + if
       if (found) then
          if (mask_a(nf) /= 0) then
             gatm_l(no)=nf
          endif
       else
          write(6,*) 'map_setgatm ERROR: pt not found, ', &
             io,jo,no,lon_l(no),lat_l(no),lat_l_local, '_o', &
             minval(lon_l),maxval(lon_l),           &
             minval(lat_l),maxval(lat_l),'_iwe',    &
             minval(lonw_a),maxval(lonw_a),         &
             minval(lone_a),maxval(lone_a),'_isn',  &
             minval(lats_a),maxval(lats_a),         &
             minval(latn_a),maxval(latn_a)
          call endrun()
       endif

    endif
    enddo
    enddo

    !--- remove fake land if there is at least one real land overlap point ---

    do ni = 1,ns_a
       found = .false.        ! check if any points are real land
       do no = 1,ns_l
          if (gatm_l(no) == ni .and. pftm_l(no) > 0 ) then
             found = .true.
          endif
       enddo
       if (found) then        ! if so, keep only real land points
       do no = 1,ns_l
          if (gatm_l(no) == ni) then
             gatm_l(no) = -1
             if (pftm_l(no) > 0 ) then
                gatm_l(no) = ni
             endif
          endif
       enddo
       endif
    enddo

    !--- check that valid fine grid points have coarse mapping gridpoints
    if (masterproc) then
    do ni = 1,ns_a
       if (mask_a(ni) /= 0) then
          found = .false.
          do no = 1,ns_l
             if (gatm_l(no) == ni) found = .true.
          enddo
          if (.not.found) then
             write(6,*) 'map_setgatm ERROR: invalid f->c index ', &
                ni,mask_a(ni)
             call endrun()
          endif
       endif
    enddo
    endif

end subroutine map_setgatm
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_setmapsFM
!
! !INTERFACE:
  subroutine map_setmapsFM(domain_a, domain_l, map1dl_a2l, map1dl_l2a, name)
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
!
! !ARGUMENTS:
  implicit none
  type(domain_type),intent(in)    :: domain_a
  type(domain_type),intent(in)    :: domain_l
  type(map_type)   ,intent(inout) :: map1dl_a2l
  type(map_type)   ,intent(inout) :: map1dl_l2a
  character(len=*) ,optional,intent(in) :: name
!
! !REVISION HISTORY:
! 2005.12.01  T Craig  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
    integer           :: ns_a
    integer           :: ns_l
    integer           :: na,nag
    integer           :: nl,nlg
!    integer          :: nlon_a       !input  grid: max number of longitude pts
!    integer          :: nlat_a       !input  grid: number of latitude  points
    real(r8),pointer :: area_a(:)    !input grid: cell area
    real(r8),pointer :: topo_a(:)    !input grid: cell topo/elevation
!    real(r8),pointer :: frac_a(:)    !input grid: cell frac
!    integer ,pointer :: mask_a(:)    !input grid: mask
!    real(r8),pointer :: lon_a(:)     !input grid: longitude (degrees)
!    real(r8),pointer :: lat_a(:)     !input grid: latitude  (degrees)
!    real(r8),pointer :: lone_a(:)    !input grid: longitude, E edge (degrees)
!    real(r8),pointer :: lonw_a(:)    !input grid: longitude, W edge (degrees)
!    real(r8),pointer :: latn_a(:)    !input grid: latitude , N edge (degrees)
!    real(r8),pointer :: lats_a(:)    !input grid: latitude , S edge (degrees)
!    integer          :: nlon_l       !output grid: max number of longitude pts
!    integer          :: nlat_l       !output grid: number of latitude  points
    real(r8),pointer :: area_l(:)    !output grid: cell area
    real(r8),pointer :: nara_l(:)    !output grid: cell equiv upscale area
    real(r8),pointer :: topo_l(:)    !output grid: cell topo/elevation
    real(r8),pointer :: ntop_l(:)    !output grid: cell equiv downscale topo
    integer ,pointer :: gatm_l(:)    !output grid: atm grid cell overlapping
!    integer ,pointer :: gatm_l2(:)   !output grid: atm grid cell overlapping
    real(r8),pointer :: frac_l(:)    !output grid: cell frac
!    real(r8),pointer :: lon_l(:)     !output grid: longitude (degrees)
!    real(r8),pointer :: lat_l(:)     !output grid: latitude  (degrees)
!    real(r8),pointer :: lone_l(:)    !output grid: longitude, E edge (degrees)
!    real(r8),pointer :: lonw_l(:)    !output grid: longitude, W edge (degrees)
!    real(r8),pointer :: latn_l(:)    !output grid: latitude , N edge (degrees)
!    real(r8),pointer :: lats_l(:)    !output grid: latitude , S edge (degrees)
!    integer ,pointer :: pftm_l(:)    !output grid: cell frac
!    integer          :: mx_a2l       !max overlapping cells
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
!    integer          :: mx_l2a       !max overlapping cells
!    integer ,pointer :: n_l2a(:,:)   !lon index, overlapping input cell
!    integer ,pointer :: i_l2a(:,:,:) !lon index, overlapping input cell
!    integer ,pointer :: j_l2a(:,:,:) !lat index, overlapping input cell
!    real(r8),pointer :: w_l2a(:,:,:) !overlap weights for input cells
!    character(len=32):: lname        !gridmap name, local variable
!    integer          :: n            !loop counters
!    integer          :: ii,ji,ni     !indices for input grid
!    integer          :: io,jo,no     !indices for output grid
!    integer          :: if,jf,nf     !found indices
!    integer          :: noffset=1    !noffset
!    real(r8)         :: doffset =360_r8 !offset value
!    real(r8)         :: offset       !offset*n_offset
!    logical          :: found        !local logical 
!    logical          :: overlapgrid  !are atm and lnd grids 1:1
!    logical          :: latlongrid   !are atm and lnd grids regular lat/lon
    integer          :: na2l,nl2a
    integer          :: begg,endg
    integer          :: abegg,aendg
    integer ,pointer :: ncnta(:)    !number of overlap points in l2a
    integer ,pointer :: cnta(:)     !number of overlap points in l2a
    real(r8),pointer :: asum(:)     !number of overlap points in l2a
    real(r8),pointer :: tsum(:)     !number of overlap points in l2a
!    real(r8)         :: sum1,sum2    !temporary sums
!    real(r8)         :: max1,max2    !temporary maxs
!    integer          :: nold         !temporary for n
!    real(r8),parameter :: relerr = 1.0e-6    ! error limit
    real(r8),parameter:: eps = 1.0e-8  ! eps for check
!
!------------------------------------------------------------------------

    !--- get local beg/end

    call get_proc_bounds    ( begg,  endg)
    call get_proc_bounds_atm(abegg, aendg)

    !--- set pointers into domains, initialize gridmap a2l gridmap ---

    call domain_setptrs(domain_a,ns=ns_a,area=area_a, &
       topo=topo_a)
    call domain_setptrs(domain_l,ns=ns_l, &
       area=area_l, frac=frac_l, &
       nara=nara_l,topo=topo_l,ntop=ntop_l, &
       gatm=gatm_l)

    !--- allocate temporaries

    allocate(ncnta(ns_a),asum(ns_a),tsum(ns_a),cnta(ns_a))

    !--- count local number of wts needed, na2l, nl2a
    !--- also set global arrays ncnta and asum for later computations

    na2l = 0
    nl2a = 0
    ncnta = 0
    asum = 0.0_r8
    do nlg = 1, ns_l
       nag = gatm_l(nlg)
       if (nag > ns_a) then
          write(6,*) 'map_setmapsFM ERROR nag > ns_a ',nag,ns_a
          call endrun()
       endif
       if (nag > 0) then
          ncnta(nag) = ncnta(nag) + 1
          asum(nag) = asum(nag) + area_l(nlg)
          if (ldecomp%glo2gdc(nlg) >=  begg .and. &
              ldecomp%glo2gdc(nlg) <=  endg) na2l = na2l + 1
          if (adecomp%glo2gdc(nag) >= abegg .and. &
              adecomp%glo2gdc(nag) <= aendg) nl2a = nl2a + 1
       endif
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
       na  = adecomp%glo2gdc(nag)

       if (nlg < 1 .or. nlg > ns_l .or. &
           nag < 1 .or. nag > ns_a .or. &
           na  < abegg .or. na > aendg) then
          write(6,*) 'map_setmapsFM ERROR in index ',nl,nlg,nag,na,ns_l,ns_a,abegg,aendg
          call endrun()
       endif

       na2l = na2l + 1
       if (na2l > n_a2l) then
          write(6,*) 'map_setmapsFM ERROR na2l > n_a2l ',na2l,n_a2l
          call endrun()
       endif
       src_a2l(na2l) = na
       dst_a2l(na2l) = nl
       wts_a2l(na2l) = 1.0_r8

       nl2a = nl2a + 1
       cnta(nag) = cnta(nag) + 1
       if (nl2a > n_l2a .or. cnta(nag) > ncnta(nag)) then
          write(6,*) 'map_setmapsFM ERROR nl2a > n_l2a ',nl2a,n_l2a,cnta(nag),ncnta(nag)
          call endrun()
       endif
       src_l2a(nl2a) = nl
       dst_l2a(nl2a) = na
       if (ncnta(nag) == 1) then
          wts_l2a(nl2a) = 1.0_r8
       else
          if (asum(nag) <= 0.0_r8) then
             write(6,*) 'map_setmapsFM ERROR asum <= 0 ',nlg,nag,asum(nag)
             call endrun()
          endif
          wts_l2a(nl2a) = (area_l(nlg)/asum(nag))
       endif
    enddo

    !--- update ldomain arrays based on upscale/downscale stuff

    tsum = 0.0_r8
    do nlg = 1, ns_l
       nag = gatm_l(nlg)
       if (frac_l(nlg) > 0.0_r8) then
          ntop_l(nlg) = topo_l(nlg)                      ! set topo ovr lnd
       else
          ntop_l(nlg) = max(0.0_r8,topo_l(nlg))          ! set topo ovr ocn 
       endif
       if (nag > 0) then
          tsum(nag) = tsum(nag) + ntop_l(nlg)*area_l(nlg)   ! area wt topo avg
       endif
    enddo

    do nlg = 1, ns_l
       nag = gatm_l(nlg)
       if (nag > 0) then
       if (ncnta(nag) == 1) then
          nara_l(nlg) = area_a(nag)
          ntop_l(nlg) = topo_a(nag)
       else
          if (asum(nag) <= 0.0_r8) then
             write(6,*) 'map_setmapsFM ERROR2 asum <= 0 ',nlg,nag,asum(nag)
             call endrun()
          endif
          if (tsum(nag) <  0.0_r8) then
             write(6,*) 'map_setmapsFM ERROR2 tsum < 0 ',nlg,nag,tsum(nag)
             call endrun()
          endif
          nara_l(nlg)   = (area_l(nlg)/asum(nag))*area_a(nag)
!-v-v-v-v-v- land topo elevation adjustment for downscaling -v-v-v-v-
          if (topo_a(nag) > 0.) then
             ntop_l(nlg)   = (ntop_l(nlg)/(tsum(nag)/asum(nag)))*topo_a(nag)
          else
             ntop_l(nlg)   = (ntop_l(nlg)-(tsum(nag)/asum(nag)))+topo_a(nag)
          endif
!-^-^-^-^-^- land topo elevation adjustment for downscaling -^-^-^-^-
       endif
       endif
    enddo

    !--- check that areas and topos match up

    asum = 0.0_r8
    tsum = 0.0_r8
    do nlg = 1, ns_l
       nag = gatm_l(nlg)
       if (nag > 0) then
       asum(nag) = asum(nag) + nara_l(nlg)
       tsum(nag) = tsum(nag) + nara_l(nlg)*ntop_l(nlg)
       endif
    enddo

    do nag = 1, ns_a
       if (asum(nag) > 0.0_r8) then
          if (abs(asum(nag)-area_a(nag)) > eps .or. &
              abs(tsum(nag)/area_a(nag) - topo_a(nag)) > eps) then
                 write(6,*) ' map_setmapsFM ERROR: ERROR in nara,ntop, ',asum(nag),area_a(nag),tsum(nag),topo_a(nag),eps
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
! !IROUTINE: gridmap_checkmap
!
! !INTERFACE:
  subroutine gridmap_checkmap(gridmap)
!
! !DESCRIPTION:
! Checks the gridmap for consistency
!
! !USES:
!
! !ARGUMENTS:
  implicit none
  type(gridmap_type),intent(in) :: gridmap
!
! !REVISION HISTORY:
! 2005.12.01  T Craig  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
    integer          :: nlon_i       !input  grid: max number of longitude pts
    integer          :: nlat_i       !input  grid: number of latitude  points
    integer          :: nlon_o       !output grid: max number of longitude pts
    integer          :: nlat_o       !output grid: number of latitude  points
    integer          :: mx_ovr       !max overlapping cells
    integer ,pointer :: n_ovr(:,:)   !lon index, overlapping input cell
    integer ,pointer :: i_ovr(:,:,:) !lon index, overlapping input cell
    integer ,pointer :: j_ovr(:,:,:) !lat index, overlapping input cell
    real(r8),pointer :: w_ovr(:,:,:) !overlap weights for input cells
    integer          :: i,j,n        !loop counters
    real(r8)         :: sum          !running sum
    real(r8)         :: rmin,rmax    !local min/max values
    integer          :: imin,imax    !local min/max values
!
!------------------------------------------------------------------------

    !--- set pointers into domains ---
    call domain_setptrs(gridmap%domain_i,ni=nlon_i,nj=nlat_i)
    call domain_setptrs(gridmap%domain_o,ni=nlon_o,nj=nlat_o)
    call gridmap_setptrs(gridmap,mx_ovr=mx_ovr,n_ovr=n_ovr,i_ovr=i_ovr, &
       j_ovr=j_ovr,w_ovr=w_ovr)

    if (masterproc) then
       write(6,*) ' '
       write(6,*) 'gridmap_checkmap name          = ',trim(gridmap%name)
       write(6,*) 'gridmap_checkmap type          = ',trim(gridmap%type)
       write(6,*) 'gridmap_checkmap src grid      = ',nlon_i,nlat_i
       write(6,*) 'gridmap_checkmap dst grid      = ',nlon_o,nlat_o
       write(6,*) 'gridmap_checkmap mx_ovr        = ',mx_ovr
       write(6,*) 'gridmap_checkmap n_ovr min/max = ',minval(n_ovr),maxval(n_ovr)
    endif

    imin = bigint
    imax = -1
    do j = 1,nlat_o
    do i = 1,nlon_o
    if (n_ovr(i,j) > 0) then
       imin = min(imin,n_ovr(i,j))
       imax = max(imax,n_ovr(i,j))
    endif
    enddo
    enddo
    if (masterproc) then
       write(6,*) 'gridmap_checkmap n_ovr nonzero = ',imin,imax
    endif

    imin = bigint
    imax = -1
    do j = 1,nlat_o
    do i = 1,nlon_o
    do n = 1,n_ovr(i,j)
       imin = min(imin,i_ovr(i,j,n))
       imax = max(imax,i_ovr(i,j,n))
    enddo
    enddo
    enddo
    if (masterproc) then
       write(6,*) 'gridmap_checkmap i_ovr min/max = ',imin,imax
    endif

    imin = bigint
    imax = -1
    do j = 1,nlat_o
    do i = 1,nlon_o
    do n = 1,n_ovr(i,j)
       imin = min(imin,j_ovr(i,j,n))
       imax = max(imax,j_ovr(i,j,n))
    enddo
    enddo
    enddo
    if (masterproc) then
       write(6,*) 'gridmap_checkmap j_ovr min/max = ',imin,imax
    endif

    rmin =  1.0e30
    rmax = -1.0e30
    do j = 1,nlat_o
    do i = 1,nlon_o
    do n = 1,n_ovr(i,j)
       rmin = min(rmin,w_ovr(i,j,n))
       rmax = max(rmax,w_ovr(i,j,n))
    enddo
    enddo
    enddo
    if (masterproc) then
       write(6,*) 'gridmap_checkmap w_ovr min/max = ',rmin,rmax
    endif

    rmin =  1.0e30
    rmax = -1.0e30
    do j = 1,nlat_o
    do i = 1,nlon_o
       sum = 0.0_r8
       do n = 1,n_ovr(i,j)
          sum = sum + w_ovr(i,j,n)
       enddo
       rmin = min(rmin,sum)
       rmax = max(rmax,sum)
    enddo
    enddo
    if (masterproc) then
       write(6,*) 'gridmap_checkmap wsum min/max = ',rmin,rmax
    endif

end subroutine gridmap_checkmap

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
!EOP
!
! !LOCAL VARIABLES:
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
    integer ,pointer :: isum(:)      !local array for deriving values
    real(r8),pointer :: rloc(:)      !local array for deriving values
    integer ,pointer :: iloc(:)      !local array for deriving values
    real(r8)         :: rmin,rmax    !local min/max values
    integer          :: imin,imax    !local min/max values
    integer          :: begg,endg    !local beg/end
    integer ,pointer :: igv1(:)      ! integer gather vector
    integer ,pointer :: igv2(:)      ! integer gather vector
    real(r8),pointer :: rgv1(:)      ! real gather vector
    real(r8),pointer :: rgv2(:)      ! real gather vector
    integer          :: ier          ! error flag
!
!------------------------------------------------------------------------

    !--- get general info ---
    call get_proc_bounds    ( begg,  endg)
    call domain_setptrs(map%domain_i,ni=nlon_i,nj=nlat_i)
    call domain_setptrs(map%domain_o,ni=nlon_o,nj=nlat_o)
    call map_setptrs(map,nwts=nwts,src=src,dst=dst,wts=wts,dstmo=dstmo)

    !--- allocate temporaries ---
    allocate(igv1(0:npes-1),igv2(0:npes-1))
    allocate(rgv1(0:npes-1),rgv2(0:npes-1))
    allocate(rsum(begg:endg),isum(begg:endg))
    allocate(rloc(nwts),iloc(nwts))

    if (masterproc) then
       write(6,*) ' '
       write(6,*) 'map_checkmap name          = ',trim(map%name)
       write(6,*) 'map_checkmap type          = ',trim(map%type)
       write(6,*) 'map_checkmap src grid      = ',nlon_i,nlat_i
       write(6,*) 'map_checkmap dst grid      = ',nlon_o,nlat_o
!       write(6,*) 'map_checkmap begg/endg     = ',begg,endg
!       write(6,*) 'map_checkmap nwts          = ',iam,nwts
       write(6,*) 'map_checkmap dstmo         = ',dstmo
!       write(6,*) 'map_checkmap src min/max   = ',minval(src),maxval(src)
!       write(6,*) 'map_checkmap dst min/max   = ',minval(dst),maxval(dst)
!       write(6,*) 'map_checkmap wts min/max   = ',minval(wts),maxval(wts)
    endif

#if (defined SPMD)
    call mpi_gather(nwts,1,MPI_INTEGER,igv1,1,MPI_INTEGER,0,mpicom,ier)
    if (ier /= 0) then
       write(6,*)'mpi_gather nwts error: ',ier
    endif
#else
    igv1(0) = nwts
#endif
    if (masterproc) then
!       write(6,*) 'map_checkmap nwts list     = ',igv1
       write(6,*) 'map_checkmap nwts gmin/max = ',minval(igv1),maxval(igv1)
       do n = 1,npes-1
          igv1(0) = igv1(0) + igv1(n)
       enddo
       write(6,*) 'map_checkmap nwts gsum     = ',igv1(0)
    endif

    imin = minval(src)
    imax = maxval(src)
#if (defined SPMD)
    call mpi_gather(imin,1,MPI_INTEGER,igv1,1,MPI_INTEGER,0,mpicom,ier)
    if (ier /= 0) then
       write(6,*)'mpi_gather src min error: ',ier
    endif
    call mpi_gather(imax,1,MPI_INTEGER,igv2,1,MPI_INTEGER,0,mpicom,ier)
    if (ier /= 0) then
       write(6,*)'mpi_gather src max error: ',ier
    endif
#else
    igv1(0) = imin
    igv2(0) = imax
#endif
    if (masterproc) then
       write(6,*) 'map_checkmap src gmin/max  = ',minval(igv1),maxval(igv2)
    endif

    imin = minval(dst)
    imax = maxval(dst)
#if (defined SPMD)
    call mpi_gather(imin,1,MPI_INTEGER,igv1,1,MPI_INTEGER,0,mpicom,ier)
    if (ier /= 0) then
       write(6,*)'mpi_gather dst min error: ',ier
    endif
    call mpi_gather(imax,1,MPI_INTEGER,igv2,1,MPI_INTEGER,0,mpicom,ier)
    if (ier /= 0) then
       write(6,*)'mpi_gather dst max error: ',ier
    endif
#else
    igv1(0) = imin
    igv2(0) = imax
#endif
    if (masterproc) then
       write(6,*) 'map_checkmap dst gmin/max  = ',minval(igv1),maxval(igv2)
    endif

    rmin = minval(wts)
    rmax = maxval(wts)
#if (defined SPMD)
    call mpi_gather(rmin,1,MPI_REAL8,rgv1,1,MPI_REAL8,0,mpicom,ier)
    if (ier /= 0) then
       write(6,*)'mpi_gather wts min error: ',ier
    endif
    call mpi_gather(rmax,1,MPI_REAL8,rgv2,1,MPI_REAL8,0,mpicom,ier)
    if (ier /= 0) then
       write(6,*)'mpi_gather wts max error: ',ier
    endif
#else
    rgv1(0) = rmin
    rgv2(0) = rmax
#endif
    if (masterproc) then
       write(6,*) 'map_checkmap wts gmin/max  = ',minval(rgv1),maxval(rgv2)
    endif

    rsum = 0.0_r8
    do n = 1,nwts
       rsum(dst(n)) = rsum(dst(n)) + wts(n)
    enddo
    rmin = minval(rsum)
    rmax = maxval(rsum)
#if (defined SPMD)
    call mpi_gather(rmin,1,MPI_REAL8,rgv1,1,MPI_REAL8,0,mpicom,ier)
    if (ier /= 0) then
       write(6,*)'mpi_gather swts min error: ',ier
    endif
    call mpi_gather(rmax,1,MPI_REAL8,rgv2,1,MPI_REAL8,0,mpicom,ier)
    if (ier /= 0) then
       write(6,*)'mpi_gather swts max error: ',ier
    endif
#else
    rgv1(0) = rmin
    rgv2(0) = rmax
#endif
    if (masterproc) then
       write(6,*) 'map_checkmap swts gmin/max = ',minval(rgv1),maxval(rgv2)
    endif

    deallocate(rloc,iloc,igv1,igv2,rgv1,rgv2)

end subroutine map_checkmap

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: areaini
!
! !INTERFACE:
  subroutine areaini (domain_i, domain_o, gridmap, &
                      fracin, fracout, name )
!
! !DESCRIPTION:
! area averaging initialization
! This subroutine is used for area-average mapping of a field from one
! grid to another.
!
!    areaini  - initializes indices and weights for area-averaging from
!               input grid to output grid
!    areamap  - called by areaini: finds indices and weights
!    areaovr  - called by areamap: finds if cells overlap and area of overlap
!    areaave  - does area-averaging from input grid to output grid
!
! To map from one grid to another, must first call areaini to build
! the indices and weights (iovr_ovr, jovr_ovr, wovr_ovr). Then must
! call areaave to get new field on output grid.
!
! Not all grid cells on the input grid will be used in the area-averaging
! of a field to the output grid. Only input grid cells with [fland_i] = 1
! contribute to output grid cell average. If [fland_i] = 0, input grid cell
! does not contribute to output grid cell. This distinction is not usually
! required for atm -> land mapping, because all cells on the atm grid have
! data. But when going from land -> atm, only land grid cells have data.
! Non-land grid cells on surface grid do not have data. So if output grid cell
! overlaps with land and non-land cells (input grid), can only use land
! grid cells when computing area-average.
!
! o Input and output grids can be ANY resolution BUT:
!
!   a. Grid orientation -- Grids can be oriented south to north
!      (i.e. cell(lat+1) is north of cell(lat)) or from north to
!      south (i.e. cell(lat+1) is south of cell(lat)). Both grids must be
!      oriented from west to east, i.e., cell(lon+1) must be east of cell(lon)
!
!   b. Grid domain -- Grids do not have to be global. Both grids are defined
!      by their north, east, south, and west edges (edge_i and edge_o in
!      this order, i.e., edge_i(1) is north and edge_i(4) is west).
!
!      For partial grids, northern and southern edges are any latitude
!      between 90 (North Pole) and -90 (South Pole). Western and eastern
!      edges are any longitude between -180 and 180, with longitudes
!      west of Greenwich negative.
!
!      For global grids, northern and southern edges are 90 (North Pole)
!      and -90 (South Pole). The grids do not have to start at the
!      same longitude, i.e., one grid can start at Dateline and go east;
!      the other grid can start at Greenwich and go east. Longitudes for
!      the western edge of the cells must increase continuously and span
!      360 degrees. Examples
!
!                              West edge    East edge
!                            ---------------------------------------------------
!      Dateline            :        -180 to 180        (negative W of Greenwich)
!      Greenwich (centered):    0 - dx/2 to 360 - dx/2
!
!   c. Both grids can have variable number of longitude points for each
!      latitude strip. However, the western edge of the first point in each
!      latitude must be the same for all latitudes. Likewise, for the
!      eastern edge of the last point. That is, each latitude strip must span
!      the same longitudes, but the number of points to do this can be different
!
!   d. One grid can be a sub-set (i.e., smaller domain) than the other grid.
!      In this way, an atmospheric dataset for the entire globe can be
!      used in a simulation for a region 30N to 50N and 130W to 70W -- the
!      code will extract the appropriate data. The two grids do not have to
!      be the same resolution. Area-averaging will work for full => partial
!      grid but obviously will not work for partial => full grid.
!
! o Field values fld_i on an  input grid with dimensions nlon_i and nlat_i =>
!   field values fld_o on an output grid with dimensions nlon_o and nlat_o as
!
!   fld_o(io,jo) =
!   fld_i(i_ovr(io,jo,    1),j_ovr(io,jo,    1)) * w_ovr(io,jo,   1)
!                             ... + ... +
!   fld_i(i_ovr(io,jo,mx_ovr),j_ovr(io,jo,mx_ovr)) * w_ovr(io,jo,mx_ovr)
!
! o Error checks:
!
!   Overlap weights of input cells sum to 1 for each output cell.
!   Global sum of dummy field is conserved for input => output area-average.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type) ,intent(in)        :: domain_i   ! input domain
    type(domain_type) ,intent(in)        :: domain_o   ! output domain
    type(gridmap_type),intent(inout)     :: gridmap    ! gridmap
    real(r8), intent(in),optional,target :: fracin(:)
    real(r8), intent(in),optional,target :: fracout(:)
    character(len=*),intent(in),optional :: name
!
! !REVISION HISTORY:
! Created by Gordon Bonan
! 2005.11.20 Updated by T Craig
!
!EOP
!
! LOCAL VARIABLES:
    integer          :: nlon_i       !input  grid: max number of longitude pts
    integer          :: nlat_i       !input  grid: number of latitude  points
    real(r8),pointer :: area_i(:)    !input grid: cell area
    real(r8),pointer :: fland_i(:)   !input grid: cell frac
    real(r8),pointer :: lone_i(:)    !input grid: longitude, E edge (degrees)
    real(r8),pointer :: lonw_i(:)    !input grid: longitude, W edge (degrees)
    real(r8),pointer :: latn_i(:)    !input grid: latitude, N edge (degrees)
    real(r8),pointer :: lats_i(:)    !input grid: latitude, S edge (degrees)
    integer          :: nlon_o       !output grid: max number of longitude pts
    integer          :: nlat_o       !output grid: number of latitude  points
    real(r8),pointer :: area_o(:)    !output grid: cell area
    real(r8),pointer :: fland_o(:)   !output grid: cell frac
    real(r8),pointer :: lone_o(:)    !output grid: longitude, E edge  (degrees)
    real(r8),pointer :: lonw_o(:)    !output grid: longitude, W edge  (degrees)
    real(r8),pointer :: latn_o(:)    !output grid: latitude, N edge (degrees)
    real(r8),pointer :: lats_o(:)    !output grid: latitude, S edge (degrees)

    real(r8),allocatable :: fld_o(:,:) !output grid: dummy field
    real(r8),allocatable :: fld_i(:,:) !input grid: dummy field
    real(r8) :: sum_fldo               !global sum of dummy output field
    real(r8) :: sum_fldi               !global sum of dummy input field
    real(r8) :: relerr = 0.00001_r8    !relative error for error checks
    integer  :: mwts                   !max number of wts per cell in map
    integer  :: ii,ji,ni               !input  grid loop index
    integer  :: io,jo,no               !output grid loop index
    real(r8) :: dx_i                   !input grid  longitudinal range
    real(r8) :: dy_i                   !input grid  latitudinal  range
    real(r8) :: dx_o                   !output grid longitudinal range
    real(r8) :: dy_o                   !output grid latitudinal  range
    character(len=32) :: lname         !gridmap name, local variable
    integer  :: ier                    !error status
!------------------------------------------------------------------------

    !--- set pointers into domain ---
    call domain_setptrs(domain_i,ni=nlon_i,nj=nlat_i,area=area_i, &
       latn=latn_i,lats=lats_i,lone=lone_i,lonw=lonw_i,frac=fland_i)
    call domain_setptrs(domain_o,ni=nlon_o,nj=nlat_o,area=area_o, &
       latn=latn_o,lats=lats_o,lone=lone_o,lonw=lonw_o,frac=fland_o)

    lname = 'areaini'
    if (present(name)) then
       lname = trim(name)
    endif

    !--- get mx_ovr, allocate gridmap, set local pointers to gridmap ---
    call areaovr(domain_i,domain_o,noffset=1,mx_ovr=mwts)
    call gridmap_init(gridmap,domain_i,domain_o,mwts,name=lname,type=map_typeglobal)

    if (present(fracin)) then
       fland_i => fracin
    else
       write(6,*) 'areaini ERROR: fracin required'
       call endrun()
    endif
    if (present(fracout)) then
       fland_o => fracout
    else
       write(6,*) 'areaini ERROR: fracout required'
       call endrun()
    endif

    ! Dynamically allocate memory

    allocate (fld_o(nlon_o,nlat_o), fld_i(nlon_i,nlat_i), stat=ier)
    if (ier /= 0) then
       write (6,*) 'areaini(): allocation error'
       call endrun
    end if

    ! Get indices and weights for mapping from input grid to output grid

    call areamap (domain_i , domain_o , gridmap, &
                  fland_i  , fland_o )

    ! Error check: global sum fld_o = global sum fld_i.
    ! This true only if both grids span the same domain.

    dx_i = lone_i(nlon_i) - lonw_i(1)
    dx_o = lone_o(nlon_o) - lonw_o(1)

    if (latn_i((nlat_i-1)*nlon_i + 1) > latn_i(1)) then      !South to North grid
       dy_i = latn_i((nlat_i-1)*nlon_i + 1) - lats_i(1)
    else                                      !North to South grid
       dy_i = latn_i(1) - lats_i((nlat_i-1)*nlon_i + 1)
    end if
    if (latn_o((nlat_o-1)*nlon_o + 1) > latn_o(1)) then      !South to North grid
       dy_o = latn_o((nlat_o-1)*nlon_o + 1) - lats_o(1)
    else                                      !North to South grid
       dy_o = latn_o(1) - lats_o((nlat_o-1)*nlon_o + 1)
    end if

    if (abs(dx_i-dx_o)>relerr .or. abs(dy_i-dy_o)>relerr) then
       write (6,*) 'AREAINI warning: conservation check not valid for'
       write (6,*) '   input  grid of ',nlon_i,' x ',nlat_i
       write (6,*) '   output grid of ',nlon_o,' x ',nlat_o
       return
    end if

    ! make dummy input field and sum globally

    sum_fldi = 0._r8
    do ji = 1, nlat_i
       do ii = 1, nlon_i
          ni = (ji-1)*nlon_i + ii
          fld_i(ii,ji) = ((ji-1)*nlon_i + ii) * fland_i(ni)
          sum_fldi = sum_fldi + area_i(ni)*fld_i(ii,ji)
       end do
    end do

    ! area-average output field from input field

    call areaave (fld_i , fld_o , gridmap)

    ! global sum of output field -- must multiply by fraction of output
    ! grid that is land as determined by input grid

    sum_fldo = 0._r8
    do jo = 1, nlat_o
       do io = 1, nlon_o
          no = (jo-1)*nlon_o + io
          sum_fldo = sum_fldo + area_o(no)*fld_o(io,jo) * fland_o(no)
       end do
    end do

    ! check for conservation

    if ( abs(sum_fldo/sum_fldi-1._r8) > relerr ) then
       write (6,*) 'AREAINI error: input field not conserved'
       write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
       write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
       call endrun
    end if

    deallocate (fld_o, fld_i)

  end subroutine areaini

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: areaave
!
! !INTERFACE:
  subroutine areaave (fld_i , fld_o , gridmap, scale_i)
!
! !DESCRIPTION:
! Mapping of field from input to output grids, 2d global fields
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8)          ,intent(in) :: fld_i(:,:)   !input grid : field
    real(r8)          ,intent(out):: fld_o(:,:)   !field for output grid
    type(gridmap_type),intent(in) :: gridmap      ! gridmap
    real(r8),optional ,intent(in) :: scale_i(:)   !input scale field
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: nlat_i    !input grid : number of latitude points
    integer  :: nlon_i    !input grid : max number longitude points
    integer  :: nlat_o    !output grid: number of latitude points
    integer  :: nlon_o    !output grid: max number of longitude points
    integer          :: mx_ovr       !max overlapping cells
    integer ,pointer :: n_ovr(:,:)   !lon index, overlapping input cell
    integer ,pointer :: i_ovr(:,:,:) !lon index, overlapping input cell
    integer ,pointer :: j_ovr(:,:,:) !lat index, overlapping input cell
    real(r8),pointer :: w_ovr(:,:,:) !overlap weights for input cells
    integer io,jo,no      !index for output grid
    integer ii,ji,ni      !index for input grid
    integer n             !overlapping cell index
!------------------------------------------------------------------------
!dir$ inlinenever areaave

    call domain_setptrs(gridmap%domain_i,ni=nlon_i,nj=nlat_i)
    call domain_setptrs(gridmap%domain_o,ni=nlon_o,nj=nlat_o)
    call gridmap_setptrs(gridmap,mx_ovr=mx_ovr,n_ovr=n_ovr,i_ovr=i_ovr,j_ovr=j_ovr,w_ovr=w_ovr)

    if (trim(gridmap%type) /= trim(map_typeglobal)) then
       write(6,*) 'areaave WARNING: gridmap type not global, ',gridmap%name,gridmap%type
    endif

    ! initialize field on output grid to zero everywhere

!$OMP PARALLEL DO PRIVATE (jo,io)
    do jo = 1, nlat_o
       do io = 1, nlon_o
          fld_o(io,jo) = 0._r8
       end do
    end do
!$OMP END PARALLEL DO

    ! loop through overlapping cells on input grid to make area-average

    do n = 1, mx_ovr
!$OMP PARALLEL DO PRIVATE (jo,io,ii,ji)
       do jo = 1, nlat_o
          do io =1, nlon_o
                ii = i_ovr(io,jo,n)
                ji = j_ovr(io,jo,n)
                ni = (ji-1)*nlon_i + ii
                if (present(scale_i)) then
                   fld_o(io,jo) = fld_o(io,jo) + &
                                  w_ovr(io,jo,n)*fld_i(ii,ji)*scale_i(ni)
                else
                   fld_o(io,jo) = fld_o(io,jo) + &
                                  w_ovr(io,jo,n)*fld_i(ii,ji)
                endif
          end do
       end do
!$OMP END PARALLEL DO
    end do

    return
  end subroutine areaave

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: areamap
!
! !INTERFACE:
  subroutine areamap (domain_i, domain_o, gridmap, &
                      fland_i  , fland_o )
!
! !DESCRIPTION:
! Weights and indices for area of overlap between grids
! Get indices and weights for area-averaging between input and output grids.
! For each output grid cell find:
!    o number of input grid cells that overlap with output grid cell (n_ovr)
!    o longitude index (1 <= i_ovr <= nlon_i) of the overlapping input grid cell
!    o latitude index  (1 <= j_ovr <= nlat_i) of the overlapping input grid cell
!    o fractional overlap of input grid cell (w_ovr)
! so that for
! field values fld_i on an  input grid with dimensions nlon_i and nlat_i
! field values fld_o on an output grid with dimensions nlon_o and nlat_o are
! fld_o(io,jo) =
! fld_i(i_ovr(io,jo,     1),j_ovr(io,jo,     1)) * w_ovr(io,jo,     1) +
!                             ... + ... +
! fld_i(i_ovr(io,jo,mx_ovr),j_ovr(io,jo,mx_ovr)) * w_ovr(io,jo,mx_ovr)
!
! Note: mx_ovr is some number greater than n_ovr. Weights of zero are
! used for the excess points
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type), intent(in) :: domain_i
    type(domain_type), intent(in) :: domain_o
    type(gridmap_type),intent(inout) :: gridmap

    real(r8),intent(in) :: fland_i(:)     ! input grid : mask (0, 1)
    real(r8),intent(in) :: fland_o(:)     ! output grid: fraction that is land
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!EOP
!
! LOCAL VARIABLES:
    integer          :: mx_ovr      ! max num input cells that overlap 
    integer ,pointer :: n_ovr(:,:)  ! number of overlapping input cells
    integer ,pointer :: i_ovr(:,:,:)! lon index, overlapping input cell
    integer ,pointer :: j_ovr(:,:,:)! lat index, overlapping input cell
    real(r8),pointer :: w_ovr(:,:,:)! overlap weights for input cells
    integer :: io,jo,no             !output grid loop index
    integer :: ii,ji,ni             !input  grid loop index
    integer :: n                    !weights index loop
    real(r8) :: f_ovr               !sum of overlap weights
    real(r8) :: relerr = 0.00001_r8 !max error: sum overlap weights ne 1
    real(r8) :: dx_i                !input grid  longitudinal range
    real(r8) :: dy_i                !input grid  latitudinal  range
    real(r8) :: dx_o                !output grid longitudinal range
    real(r8) :: dy_o                !output grid latitudinal  range
    integer  :: nlon_i              !input size, i
    integer  :: nlat_i              !input size, j
    real(r8),pointer :: lone_i(:)   !input grid: longitude, E edge (degrees)
    real(r8),pointer :: lonw_i(:)   !input grid: longitude, W edge (degrees)
    real(r8),pointer :: latn_i(:)   !input grid: latitude, N edge (degrees)
    real(r8),pointer :: lats_i(:)   !input grid: latitude, S edge (degrees)
    integer  :: nlon_o              !output size, i
    integer  :: nlat_o              !output size, j
    real(r8),pointer :: lone_o(:)   !output grid: longitude, E edge (degrees)
    real(r8),pointer :: lonw_o(:)   !output grid: longitude, W edge (degrees)
    real(r8),pointer :: latn_o(:)   !output grid: latitude, N edge (degrees)
    real(r8),pointer :: lats_o(:)   !output grid: latitude, S edge (degrees)
    real(r8),pointer :: area_o(:)   ! output grid: cell area
!------------------------------------------------------------------------

    call domain_setptrs(domain_i,ni=nlon_i,nj=nlat_i, &
       latn=latn_i,lats=lats_i,lone=lone_i,lonw=lonw_i)
    call domain_setptrs(domain_o,ni=nlon_o,nj=nlat_o, &
       latn=latn_o,lats=lats_o,lone=lone_o,lonw=lonw_o,area=area_o)
    call gridmap_setptrs(gridmap,mx_ovr=mx_ovr,n_ovr=n_ovr, &
       i_ovr=i_ovr,j_ovr=j_ovr,w_ovr=w_ovr)

    ! --------------------------------------------------------------------
    ! Initialize overlap weights on output grid to zero for maximum
    ! number of overlapping points. Set lat and lon indices of overlapping
    ! input cells to dummy values. Set number of overlapping cells to zero
    ! --------------------------------------------------------------------

    do n = 1, mx_ovr
       do jo = 1, nlat_o
          do io = 1, nlon_o
             i_ovr(io,jo,n) = 1
             j_ovr(io,jo,n) = 1
             w_ovr(io,jo,n) = 0._r8
          end do
       end do
    end do

    do jo = 1, nlat_o
       do io = 1, nlon_o
          n_ovr(io,jo) = 0
       end do
    end do

    call areaovr (domain_i, domain_o, &
                  noffset=1, &
                  n_ovr=n_ovr, i_ovr=i_ovr, j_ovr=j_ovr, w_ovr=w_ovr  )

    ! --------------------------------------------------------------------
    ! Normalize areas of overlap to get fractional contribution of each
    ! overlapping grid cell (input grid) to grid cell average on output grid.
    ! Normally, do this by dividing area of overlap by area of output grid cell.
    ! But, only have data for land cells on input grid. So if output grid cell
    ! overlaps with land and non-land cells (input grid), do not have valid
    ! non-land data for area-average. Instead, weight by area of land using
    ! [fland_i], which has a value of one for land and zero for ocean. If
    ! [fland_i] = 1, input grid cell contributes to output grid cell average.
    ! If [fland_i] = 0, input grid cell does not contribute to output grid cell
    ! average.
    ! --------------------------------------------------------------------

    do jo = 1, nlat_o
       do io = 1, nlon_o
          no = (jo-1)*nlon_o + io

          ! find total land area of overlapping input cells

          f_ovr = 0._r8
          do n = 1, n_ovr(io,jo)
             ii = i_ovr(io,jo,n)
             ji = j_ovr(io,jo,n)
             ni = (ji-1)*nlon_i + ii
             f_ovr = f_ovr + w_ovr(io,jo,n)*fland_i(ni)
          end do

          ! make sure area of overlap is less than or equal to output grid cell area

          if ((f_ovr-area_o(no))/area_o(no) > relerr) then
             write (6,*) 'AREAMAP error: area not conserved for lon,lat = ',io,jo,no
             write (6,'(a30,e20.10)') 'sum of overlap area = ',f_ovr
             write (6,'(a30,e20.10)') 'area of output grid = ',area_o(no)
             call endrun
          end if

          ! make weights

          do n = 1, n_ovr(io,jo)
             ii = i_ovr(io,jo,n)
             ji = j_ovr(io,jo,n)
             ni = (ji-1)*nlon_i + ii
             if (f_ovr > 0._r8) then
                w_ovr(io,jo,n) = w_ovr(io,jo,n)*fland_i(ni) / f_ovr
             else
                w_ovr(io,jo,n) = 0._r8
             end if
          end do

       end do
    end do

    ! --------------------------------------------------------------------
    ! Error check: overlap weights for input grid cells must sum to 1. This
    ! is always true if both grids span the same domain. However, if one
    ! grid is a subset of the other grid, this is only true when mapping
    ! from the full grid to the subset. When input grid covers a smaller
    ! domain than the output grid, this test is not valid.
    ! --------------------------------------------------------------------

    dx_i = lone_i(nlon_i) - lonw_i(1)
    dx_o = lone_o(nlon_o) - lonw_o(1)

    if (latn_i((nlat_i-1)*nlon_i + 1) > latn_i(1)) then      !South to North grid
       dy_i = latn_i((nlat_i-1)*nlon_i + 1) - lats_i(1)
    else                                      !North to South grid
       dy_i = latn_i(1) - lats_i((nlat_i-1)*nlon_i + 1)
    end if
    if (latn_o((nlat_o-1)*nlon_o + 1) > latn_o(1)) then      !South to North grid
       dy_o = latn_o((nlat_o-1)*nlon_o + 1) - lats_o(1)
    else                                      !North to South grid
       dy_o = latn_o(1) - lats_o((nlat_o-1)*nlon_o + 1)
    end if

    if (abs(dx_i-dx_o)>relerr .or. abs(dy_i-dy_o)>relerr) then
       if (dx_i<dx_o .or. dy_i<dy_o) then
          write (6,*) 'AREAMAP warning: area-average not valid for '
          write (6,*) '   input  grid of ',nlon_i,' x ',nlat_i
          write (6,*) '   output grid of ',nlon_o,' x ',nlat_o
          return
       end if
    end if

    do jo = 1, nlat_o
       do io = 1, nlon_o
          no = (jo-1)*nlon_o + io
          f_ovr = 0._r8

          do n = 1, mx_ovr
             f_ovr = f_ovr + w_ovr(io,jo,n)
          end do

          ! error check only valid if output grid cell has land. non-land cells
          ! will have weights equal to zero

          if (fland_o(no) > 0._r8) then
             if (abs(f_ovr-1._r8) > relerr) then
                write (6,*) 'AREAMAP error: area not conserved for lon,lat = ',io,jo,no
                write (6,'(a30,e20.10)') 'sum of overlap weights = ',f_ovr
                call endrun
             end if
          end if

       end do
    end do

    return
  end subroutine areamap

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: areaovr
!
! !INTERFACE:
  subroutine areaovr (domain_i, domain_o, &
                      noffset, &
                      mx_ovr , n_ovr  , i_ovr    , j_ovr , w_ovr )
!
! !DESCRIPTION:
! Find area of overlap between grid cells
! For each output grid cell: find overlapping input grid cell and area of
! input grid cell that overlaps with output grid cell. Cells overlap if:
!
! southern edge of input grid < northern edge of output grid AND
! northern edge of input grid > southern edge of output grid
!
! western edge of input grid < eastern edge of output grid AND
! eastern edge of input grid > western edge of output grid
!
!           lon_o(io,jo)      lon_o(io+1,jo)
!
!              |                   |
!              --------------------- lat_o(jo+1)
!              |                   |
!              |                   |
!    xxxxxxxxxxxxxxx lat_i(ji+1)   |
!    x         |   x               |
!    x  input  |   x   output      |
!    x  cell   |   x    cell       |
!    x  ii,ji  |   x   io,jo       |
!    x         |   x               |
!    x         ----x---------------- lat_o(jo  )
!    x             x
!    xxxxxxxxxxxxxxx lat_i(ji  )
!    x             x
! lon_i(ii,ji) lon_i(ii+1,ji)
!
!
! The above diagram assumes both grids are oriented South to North. Other
! combinations of North to South and South to North grids are possible:
!
!     Input Grid    Output Grid
!     -------------------------
! (1)   S to N        S to N
! (2)   N to S        N to S
! (3)   S to N        N to S
! (4)   N to S        S to N
!
! The code has been modified to allow for North to South grids. Verification
! that these changes work are:
!    o (1) and (4) give same results for output grid
!    o (2) and (3) give same results for output grid
!    o (2) and (4) give same results for output grid when output grid inverted
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(in) :: domain_i
    type(domain_type),intent(in) :: domain_o

    integer , intent(in)   ,optional :: noffset      ! number of offsets to 
                                                     ! try to find overlaps
    integer , intent(inout),optional :: mx_ovr       !max num of overlap points
    integer , intent(inout),optional :: n_ovr(:,:)   !number of overlap pts
    integer , intent(inout),optional :: i_ovr(:,:,:) !lon index of overlap pts
    integer , intent(inout),optional :: j_ovr(:,:,:) !lat index of overlap pts
    real(r8), intent(inout),optional :: w_ovr(:,:,:) !weight of overlap pts
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!EOP
!
! LOCAL VARIABLES:
    integer           :: nlon_i       !input grid : max lon points
    integer           :: nlat_i       !input grid : lat points
    real(r8), pointer :: lone_i(:)    !input grid : cell E edge lon (deg)
    real(r8), pointer :: lonw_i(:)    !input grid : cell W edge lon (deg)
    real(r8), pointer :: latn_i(:)    !input grid : cell N edge lat (deg)
    real(r8), pointer :: lats_i(:)    !input grid : cell S edge lat (deg)
    integer           :: nlon_o       !output grid: max lon points
    integer           :: nlat_o       !output grid: lat points
    real(r8), pointer :: lone_o(:)    !output grid: cell E edge lon (deg)
    real(r8), pointer :: lonw_o(:)    !output grid: cell W edge lon (deg)
    real(r8), pointer :: latn_o(:)    !output grid: cell N edge lat (deg)
    real(r8), pointer :: lats_o(:)    !output grid: cell S edge lat (deg)

    integer, parameter :: mx_ovr_ceiling = 100000 !emergency limit
    integer io,jo,no       !output grid loop index
    integer indexo         !output grid lat. index according to orientn
    integer ii,ji,ni       !input  grid loop index
    integer indexi         !input grid lat. index according to orientn
    real(r8) lonw          !west longitudes of overlap
    real(r8) lone          !east longitudes of overlap
    real(r8) dx            !difference in longitudes
    real(r8) lats          !south latitudes of overlap
    real(r8) latn          !north latitudes of overlap
    real(r8) dy            !difference in latitudes
    integer  size3         !size of 3rd dim in map arrays
    real(r8) deg2rad       !pi/180
    real(r8) a_ovr         !area of overlap
    integer  n             !overlapping cell index
    integer  noffsetl      !local, number of offsets to test, 0=none, default=1
    real(r8) offset        !value of offset used to shift x-grid 360 degrees
    integer ,allocatable :: n_ovrl(:,:)   ! local copy of novr
    real(r8),allocatable :: lonw_il(:)    ! local copy of lonw_i with offset
    real(r8),allocatable :: lone_il(:)    ! local copy of lone_i with offset
    integer ier            ! error flag
!------------------------------------------------------------------------
 
    call domain_setptrs(domain_i,ni=nlon_i,nj=nlat_i, &
                        lats=lats_i,latn=latn_i,lonw=lonw_i,lone=lone_i)
    call domain_setptrs(domain_o,ni=nlon_o,nj=nlat_o, &
                        lats=lats_o,latn=latn_o,lonw=lonw_o,lone=lone_o)

    allocate(n_ovrl(nlon_o,nlat_o),lonw_il(nlon_i*nlat_i), &
                                  lone_il(nlon_i*nlat_i),stat=ier)
    if (ier /= 0) then
       write (6,*) 'areaovr(): allocation error'
       call endrun
    end if

    deg2rad = SHR_CONST_PI / 180._r8
    noffsetl = 1
    if (present(noffset)) then
       noffsetl = noffset
    endif
    size3 = mx_ovr_ceiling
    if (present(w_ovr)) then
      size3 = size(w_ovr,3)
    endif
    n_ovrl(:,:) = 0

    do n = 0,noffsetl   ! loop through offsets

       if (lonw_i(1) < lonw_o(1)) then
          offset = (n*360)
       else
          offset = -(n*360)
       end if
       do ji = 1, nlat_i
       do ii = 1, nlon_i
          ni = (ji-1)*nlon_i + ii
          lonw_il(ni) = lonw_i(ni) + offset
          lone_il(ni) = lone_i(ni) + offset
       end do
       end do

       ! for all output grid cells-

       do jo = 1, nlat_o

       if (latn_o((nlat_o-1)*nlon_o + 1) > latn_o(1)) then
          indexo  = jo          !south to north at the center of cell
       else
          indexo  = nlat_o+1-jo !north to south at the center of cell
       end if

       do io = 1, nlon_o

          ! loop through all input grid cells to find overlap with output grid

          do ji = 1, nlat_i

          if (latn_i((nlat_i-1)*nlon_i + 1) > latn_i(1)) then
             indexi  = ji          !south to north at the center of cell
          else
             indexi  = nlat_i+1-ji !north to south at the center of cell
          end if

          ! lats overlap

          if ( lats_i((indexi-1)*nlon_i + 1)<latn_o((indexo-1)*nlon_o + 1) .and. &
               latn_i((indexi-1)*nlon_i + 1)>lats_o((indexo-1)*nlon_o + 1) ) then

          do ii = 1, nlon_i

             ! lons overlap

             ni = (indexi-1)*nlon_i + ii
             no = (indexo-1)*nlon_o + io

             if (lonw_il(ni)<lone_o(no) .and. &
                 lone_il(ni)>lonw_o(no)) then

                ! increment number of overlapping cells.
                ! make sure 0 < n_ovrl < size3, not bigger than dimension

                n_ovrl(io,indexo) = n_ovrl(io,indexo) + 1
                if (n_ovrl(io,indexo) > mx_ovr_ceiling) then
                   write (6,*) 'AREAOVR error: n_ovr= ', &
                      n_ovrl(io,indexo),' exceeded mx_ovr_ceiling = ', &
                      mx_ovr_ceiling,' for output lon,lat = ',io,indexo
                   call endrun
                end if
                if (n_ovrl(io,indexo) > size3) then
                   write (6,*) 'AREAOVR error: n_ovr= ', &
                      n_ovrl(io,indexo),' exceeded size of arrays = ', &
                      size3,' for output lon,lat = ',io,indexo
                   call endrun
                end if

                if (present(i_ovr).and.present(j_ovr).and.present(w_ovr)) then
                   ! determine area of overlap

                   lone = min(lone_o(no),lone_il(ni))*deg2rad 
                   lonw = max(lonw_o(no),lonw_il(ni))*deg2rad 
                   dx = max(0.0_r8,(lone-lonw))
                   latn = min(latn_o(no),latn_i(ni))*deg2rad 
                   lats = max(lats_o(no),lats_i(ni))*deg2rad 
                   dy = max(0.0_r8,(sin(latn)-sin(lats)))
                   a_ovr = dx*dy*re*re

                   ! save lat, lon, area

                   i_ovr(io,indexo,n_ovrl(io,indexo)) = ii
                   j_ovr(io,indexo,n_ovrl(io,indexo)) = indexi
                   w_ovr(io,indexo,n_ovrl(io,indexo)) = a_ovr
                endif

             end if
          end do
          end if
          end do

       end do
       end do

    enddo   ! offset loop

    if (present(n_ovr)) then
       n_ovr = n_ovrl
    endif
    if (present(mx_ovr)) then
       mx_ovr = maxval(n_ovrl)
    endif

    deallocate(n_ovrl,lonw_il,lone_il)

    return
  end subroutine areaovr

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cellarea_regional
!
! !INTERFACE:
  subroutine cellarea_regional (domain, edgen, edgee, edges, edgew)
!
! !DESCRIPTION:
! Comute area of grid cells (square kilometers) - regional grid
! Verify total area from grid cells is same as area of grid
! as defined by its edges
! (can become global as special case)
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type), intent(inout) :: domain
    real(r8), intent(in) :: edges            
    real(r8), intent(in) :: edgen            
    real(r8), intent(in) :: edgew            
    real(r8), intent(in) :: edgee            
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: nlat           
    integer  :: nlon           
    real(r8), pointer :: lats(:)
    real(r8), pointer :: latn(:)
    real(r8), pointer :: lonw(:)
    real(r8), pointer :: lone(:)
    real(r8), pointer :: area(:)  
    integer i,j,n               !indices
    real(r8) deg2rad            !pi/180
    real(r8) global             !summed area
    real(r8) dx                 !cell width: E-W
    real(r8) dy                 !cell width: N-S
    real(r8) garea              !true area for error check
!------------------------------------------------------------------------

    call domain_setptrs(domain,ni=nlon,nj=nlat,area=area, &
                        lats=lats,latn=latn,lonw=lonw,lone=lone)

    !--- compute area from lats/lons ---
    call cellarea_global(domain)

    !--- sum local areas ---
    global = 0._r8
    do n = 1, nlat*nlon
       global = global + area(n)
    end do

    !--- compute global area ---
    deg2rad = (SHR_CONST_PI) / 180._r8
    dx = (edgee - edgew) * deg2rad
    dy = sin(edgen*deg2rad) - sin(edges*deg2rad)
    garea =  dx*dy*re*re

    if (abs(global-garea)/garea > 0.00001_r8) then
       write (6,*) 'CELLAREA error: correct area is ',garea, &
            ' but summed area of grid cells is ',global
       call endrun
    end if

    return
  end subroutine cellarea_regional

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cellarea_global
!
! !INTERFACE:
  subroutine cellarea_global (domain)
!
! !DESCRIPTION:
! Area of grid cells (square kilometers)- global grid
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type), intent(inout) :: domain
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: nlat           
    integer  :: nlon           
    real(r8), pointer :: lats(:)
    real(r8), pointer :: latn(:)
    real(r8), pointer :: lonw(:)
    real(r8), pointer :: lone(:)
    real(r8), pointer :: area(:)  
    integer i,j,n               !indices
    real(r8) deg2rad            !pi/180
    real(r8) dx                 !cell width: E-W
    real(r8) dy                 !cell width: N-S
!------------------------------------------------------------------------

    ! Note: supports general lat/lon grids

    call domain_setptrs(domain,ni=nlon,nj=nlat,area=area, &
                        lats=lats,latn=latn,lonw=lonw,lone=lone)

    deg2rad = (SHR_CONST_PI) / 180._r8
    do n = 1, nlat*nlon
       dx = (lone(n) - lonw(n)) * deg2rad
       dy = sin(latn(n)*deg2rad) - sin(lats(n)*deg2rad) 
       area(n) = dx*dy*re*re
    end do

    return
  end subroutine cellarea_global

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: celledge_regional
!
! !INTERFACE:
  subroutine celledge_regional (domain, edgen, edgee, edges, edgew)
!
! !DESCRIPTION:
! Southern and western edges of grid cells - regional grid
! (can become global as special case)
! Latitudes -- southern/northern edges for each latitude strip.
! For grids oriented South to North, the southern
! and northern edges of latitude strip [j] are:
!        southern = lats(j  )
!        northern = lats(j+1)
! For grids oriented North to South: the southern
! and northern edges of latitude strip [j] are:
!        northern = lats(j  )
!        southern = lats(j+1)
! In both cases, [lats] must be dimensioned lats(lat+1)
! Longitudes -- western edges. Longitudes for the western edge of the
! cells must increase continuously and span 360 degrees. Assume that
! grid starts at Dateline with western edge on Dateline Western edges
! correspond to [lonc] (longitude at center of cell) and range from
! -180 to 180 with negative longitudes west of Greenwich.
! Partial grids that do not span 360 degrees are allowed so long as they
! have the convention of Grid 1 with
!      western edge of grid: >= -180 and < 180
!      eastern edge of grid: > western edge  and <= 180
! [lonw] must be dimensioned lonw(lon+1,lat) because each latitude
! strip can have variable longitudinal resolution
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout) :: domain
    real(r8), intent(in) :: edgen             !northern edge of grid (degrees)
    real(r8), intent(in) :: edgee             !eastern edge of grid (degrees)
    real(r8), intent(in) :: edges             !southern edge of grid (degrees)
    real(r8), intent(in) :: edgew             !western edge of grid (degrees)

!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 2005.11.20 Updated to domain datatype by T Craig
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: nlon              
    integer  :: nlat              
    real(r8),pointer :: lonc(:) 
    real(r8),pointer :: latc(:) 
    real(r8),pointer :: lats(:)   
    real(r8),pointer :: latn(:)   
    real(r8),pointer :: lonw(:)   
    real(r8),pointer :: lone(:)   
    integer i,j,n,nim1,njm1 !indices
    real(r8) dx             !cell width
!------------------------------------------------------------------------

    nlon = domain%ni
    nlat = domain%nj
    lonc => domain%lonc
    latc => domain%latc
    lats => domain%lats
    latn => domain%latn
    lonw => domain%lonw
    lone => domain%lone

    ! Latitudes
    ! Assumes lats are constant on an i line

    if (nlat == 1) then                      ! single latitude
       lats(:)    = edges
       latn(:)    = edgen
    elseif (latc(nlon+1) > latc(1)) then  ! South to North grid
       lats(:) = edges
       latn(:) = edgen
       do j = 2, nlat
       do i = 1, nlon
          n    = (j-1)*nlon + i
          njm1 = (j-2)*nlon + i
          lats(n) = (latc(njm1) + latc(n)) / 2._r8
          latn(njm1) = lats(n)
       end do
       end do
    else                                     ! North to South grid
       lats(:) = edges
       latn(:) = edgen
       do j = 2, nlat
       do i = 1, nlon
          n    = (j-1)*nlon + i
          njm1 = (j-2)*nlon + i
          latn(n)    = (latc(njm1) + latc(n)) / 2._r8
          lats(njm1) = latn(n)
       end do
       end do
    end if

    ! Longitudes
    ! Western edge of first grid cell -- since grid starts with western
    ! edge on Dateline, lonw(1,j)=-180. This is the same as [edgew].

    lonw(:) = edgew
    lone(:) = edgee
    dx = (edgee - edgew) / nlon
    do j = 1, nlat
    do i = 2, nlon
       n    = (j-1)*nlon + i
       nim1 = (j-1)*nlon + i-1
       lonw(n)    = lonw(n) + (i-1)*dx
       lone(nim1) = lonw(n)
    end do
    end do

    return

  end subroutine celledge_regional

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: celledge_global
!
! !INTERFACE:
  subroutine celledge_global (domain)
!
! !DESCRIPTION:
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout) :: domain

!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 2005.11.20 Updated to domain datatype by T Craig
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: nlon              
    integer  :: nlat              
    real(r8),pointer :: lonc(:) 
    real(r8),pointer :: latc(:) 
    real(r8),pointer :: lats(:)   
    real(r8),pointer :: latn(:)   
    real(r8),pointer :: lonw(:)   
    real(r8),pointer :: lone(:)   
    integer i,j,n,nim1,njm1 !indices
    real(r8) dx             !cell width
!------------------------------------------------------------------------

    nlon = domain%ni
    nlat = domain%nj
    lonc => domain%lonc
    latc => domain%latc
    lats => domain%lats
    latn => domain%latn
    lonw => domain%lonw
    lone => domain%lone

    ! Latitudes
    lats(:) = -90._r8
    latn(:) =  90._r8
    do j = 2, nlat   
    do i = 1, nlon
       n    = (j-1)*nlon + i
       njm1 = (j-2)*nlon + i
       lats(n) = (latc(njm1) + latc(n)) / 2._r8
       latn(njm1) = lats(n)
    end do
    end do

    ! Longitudes

    if (lonc(1) >= 0._r8) then
       dx = 360._r8/(nlon)
       lonw(:) = -dx/2._r8
       lone(:) = -dx/2._r8 + (nlon)*dx
       do j = 1, nlat
       do i = 2, nlon
          n    = (j-1)*nlon + i
          nim1 = (j-1)*nlon + i-1
          lonw(n)    = -dx/2._r8 + (i-1)*dx
          lone(nim1) = lonw(n)
       end do
       end do
    else
       write(6,*)'global non-regional grids currently only supported ', &
            'for grids starting at greenwich and centered on Greenwich'
       call endrun
    endif

    return
  end subroutine celledge_global

!-----------------------------------------------------------------------

end module areaMod



















