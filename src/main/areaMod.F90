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

  character(len=16),parameter,public :: gridmap_typelocal  = 'local'
  character(len=16),parameter,public :: gridmap_typeglobal = 'global'
  character(len=16),parameter,public :: gridmap_typedst    = 'dst'
  character(len=16),parameter,public :: gridmap_typesrc    = 'src'

!
! !PUBLIC MEMBER FUNCTIONS:
  public :: gridmap_init
  public :: gridmap_setptrs
  public :: gridmap_maparray
  public :: gridmap_setmapsFM
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
  gridmap%i_ovr = bigint
  gridmap%j_ovr = bigint
  gridmap%w_ovr = nan

end subroutine gridmap_init
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
! !IROUTINE: gridmap_maparray
!
! !INTERFACE:
  subroutine gridmap_maparray(fld_i,fld_o,gridmap,type)
!
! !DESCRIPTION:
! This subroutine maps arrays, local 1d
! Need a decomp, use type to set decomps for now, could be improved (tcx fix)
!
! !USES:
  use decompMod, only : ldecomp, adecomp, get_proc_bounds, get_proc_bounds_atm, decomp_type
!
! !ARGUMENTS:
  implicit none
  real(r8),pointer                :: fld_i(:)
  real(r8),pointer                :: fld_o(:)
  type(gridmap_type), intent(in)  :: gridmap
  character(len=*)  , intent(in)  :: type
!
! !REVISION HISTORY:
! 2005.11.15  T Craig  Creation.
!
!EOP

! !LOCAL VARIABLES:
  integer :: begg_o,endg_o         !beg,end of output grid
  integer :: begg_i,endg_i         !beg,end of input grid
  integer :: g_o ,g_i              !gridcell indices
  integer          :: mx_ovr       !max overlapping cells
  integer ,pointer :: n_ovr(:,:)   !lon index, overlapping input cell
  integer ,pointer :: i_ovr(:,:,:) !lon index, overlapping input cell
  integer ,pointer :: j_ovr(:,:,:) !lat index, overlapping input cell
  real(r8),pointer :: w_ovr(:,:,:) !overlap weights for input cells
  integer          :: n            !loop counters
  integer          :: ii,ji        !indices for input grid
  integer          :: io,jo        !indices for output grid
  logical          :: a2ltype      !a2l or l2a type
  type(decomp_type),pointer :: decomp_o
  type(decomp_type),pointer :: decomp_i

!
!------------------------------------------------------------------------------

    if (trim(type) == 'a2l') then
       a2ltype = .true.
    elseif (trim(type) == 'l2a') then
       a2ltype = .false.
    else
       write(6,*) 'gridmap_maparry ERROR type must be a2l or l2a:',trim(type)
       call endrun()
    endif

    if (a2ltype) then
       decomp_i => adecomp
       call get_proc_bounds_atm(begg_i, endg_i)
       decomp_o => ldecomp
       call get_proc_bounds    (begg_o, endg_o)
    else
       decomp_i => ldecomp
       call get_proc_bounds    (begg_i, endg_i)
       decomp_o => adecomp
       call get_proc_bounds_atm(begg_o, endg_o)
    endif

    call gridmap_setptrs(gridmap,mx_ovr=mx_ovr,n_ovr=n_ovr,i_ovr=i_ovr, &
       j_ovr=j_ovr,w_ovr=w_ovr)

    if (trim(gridmap%type) /= trim(gridmap_typelocal)) then
       write(6,*) 'gridmap_maparray WARNING: gridmap type not correct, ', &
                   gridmap%name,gridmap%type
    endif

    ! initialize field on output grid to zero everywhere

    do g_o = begg_o,endg_o
       fld_o(g_o) = 0._r8
    end do

    ! loop through overlapping cells on input grid to make area-average

    do g_o = begg_o,endg_o
       io = decomp_o%gdc2i(g_o)
       jo = decomp_o%gdc2j(g_o)
       do n = 1, n_ovr(io,jo)
          ii = i_ovr(io,jo,n)
          ji = j_ovr(io,jo,n)
          g_i = decomp_i%ij2gdc(ii,ji)
          if (g_i < begg_i .or. g_i > endg_i) then
             write(6,*) 'gridmap_maparray ERROR: g_i out of bounds ',g_i, &
                         begg_i,endg_i,ii,ji,io,jo,g_o,n
             call endrun()
          endif
          fld_o(g_o) = fld_o(g_o) + w_ovr(io,jo,n)*fld_i(g_i)
       end do
    end do

end subroutine gridmap_maparray

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gridmap_setmapsFM
!
! !INTERFACE:
  subroutine gridmap_setmapsFM(domain_i, domain_o, gridmap_i2o, gridmap_o2i, name)
!
! !DESCRIPTION:
! Set course to fine mesh maps and reverse.  domain_i should be coarse
! (atm) mesh, domain_o is fine (land) mesh.
! Simple overlap algorithm
!   - Find every fine gridcell within coarse gridcell
!   - Keep "real" cells unless there are none
!   - Weights based on areas of cells used, sum(weights)==1
!   - Use i2o mapping to set o2i mapping
!
! !USES:
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(in)    :: domain_i
  type(domain_type), intent(in)    :: domain_o
  type(gridmap_type),intent(inout) :: gridmap_i2o
  type(gridmap_type),intent(inout) :: gridmap_o2i
  character(len=*) ,optional,intent(in) :: name
!
! !REVISION HISTORY:
! 2005.12.01  T Craig  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
    integer          :: nlon_i       !input  grid: max number of longitude pts
    integer          :: nlat_i       !input  grid: number of latitude  points
    real(r8),pointer :: area_i(:,:)  !input grid: cell area
    real(r8),pointer :: fland_i(:,:) !input grid: cell frac
    integer ,pointer :: mask_i(:,:)  !input grid: mask
    real(r8),pointer :: lon_i(:,:)   !input grid: longitude (degrees)
    real(r8),pointer :: lat_i(:,:)   !input grid: latitude  (degrees)
    real(r8),pointer :: lone_i(:,:)  !input grid: longitude, E edge (degrees)
    real(r8),pointer :: lonw_i(:,:)  !input grid: longitude, W edge (degrees)
    real(r8),pointer :: latn_i(:,:)  !input grid: latitude , N edge (degrees)
    real(r8),pointer :: lats_i(:,:)  !input grid: latitude , S edge (degrees)
    integer          :: nlon_o       !output grid: max number of longitude pts
    integer          :: nlat_o       !output grid: number of latitude  points
    real(r8),pointer :: area_o(:,:)  !output grid: cell area
    real(r8),pointer :: nara_o(:,:)  !output grid: cell equiv upscale area
    real(r8),pointer :: fland_o(:,:) !output grid: cell frac
    real(r8),pointer :: lon_o(:,:)   !output grid: longitude (degrees)
    real(r8),pointer :: lat_o(:,:)   !output grid: latitude  (degrees)
    real(r8),pointer :: lone_o(:,:)  !output grid: longitude, E edge (degrees)
    real(r8),pointer :: lonw_o(:,:)  !output grid: longitude, W edge (degrees)
    real(r8),pointer :: latn_o(:,:)  !output grid: latitude , N edge (degrees)
    real(r8),pointer :: lats_o(:,:)  !output grid: latitude , S edge (degrees)
    integer ,pointer :: pftm_o(:,:) !output grid: cell frac
    integer          :: mx_i2o       !max overlapping cells
    integer ,pointer :: n_i2o(:,:)   !lon index, overlapping input cell
    integer ,pointer :: i_i2o(:,:,:) !lon index, overlapping input cell
    integer ,pointer :: j_i2o(:,:,:) !lat index, overlapping input cell
    real(r8),pointer :: w_i2o(:,:,:) !overlap weights for input cells
    integer          :: mx_o2i       !max overlapping cells
    integer ,pointer :: n_o2i(:,:)   !lon index, overlapping input cell
    integer ,pointer :: i_o2i(:,:,:) !lon index, overlapping input cell
    integer ,pointer :: j_o2i(:,:,:) !lat index, overlapping input cell
    real(r8),pointer :: w_o2i(:,:,:) !overlap weights for input cells
    character(len=32):: lname        !gridmap name, local variable
    integer          :: n            !loop counters
    integer          :: ii,ji        !indices for input grid
    integer          :: io,jo        !indices for output grid
    integer          :: if,jf        !found indices
    integer          :: noffset=1    !noffset
    real(r8)         :: doffset =360_r8 !offset value
    real(r8)         :: offset       !offset*n_offset
    logical          :: found        !local logical 
    logical          :: overlapgrid  !are atm and lnd grids 1:1
    logical          :: latlongrid   !are atm and lnd grids regular lat/lon
    integer ,pointer :: ncnta(:,:)   !number of overlap points in o2i
    real(r8)         :: sum          !sum of weights
    integer          :: nold         !temporary for n
    real(r8),parameter :: relerr = 1.0e-6    ! error limit
!tcx fix, lat_o_local should be removed when limit no longer needed
    real(r8)         :: lat_o_local  !local copy of lat_o(io,jo), adjusted
!
!------------------------------------------------------------------------

    !--- set pointers into domains, initialize gridmap i2o gridmap ---

    call domain_setptrs(domain_i,ni=nlon_i,nj=nlat_i,area=area_i, &
       latc=lat_i,lonc=lon_i,frac=fland_i,mask=mask_i, &
       latn=latn_i,lats=lats_i,lone=lone_i,lonw=lonw_i)
    call domain_setptrs(domain_o,ni=nlon_o,nj=nlat_o,pftm=pftm_o,area=area_o, &
       latc=lat_o,lonc=lon_o,frac=fland_o,nara=nara_o, &
       latn=latn_o,lats=lats_o,lone=lone_o,lonw=lonw_o)

    mx_i2o = 1
    lname = 'setmapsFM_a2l'
    if (present(name)) then
      lname = trim(name)//'_a2l'
    endif
    call gridmap_init(gridmap_i2o,domain_i,domain_o,mx_i2o,name=lname, &
       type=gridmap_typelocal)
    call gridmap_setptrs(gridmap_i2o,mx_ovr=mx_i2o,n_ovr=n_i2o,i_ovr=i_i2o, &
       j_ovr=j_i2o,w_ovr=w_i2o)

    !--- search for the overlap, input to output, 1:1 "disaggregation" ---
    !--- this is the coarse to fine map where there should be exactly
    !--- one coarse overlap point for each fine point 
    !--- three possible search algorithms, overlapgrid, latlongrid, neither
    !--- figure out which search scheme to use

    !--- overlapgrid means both grids are identical
    overlapgrid = .false.
    if (nlon_i == nlon_o .and. nlat_i == nlat_o) then
       overlapgrid = .true.
       do jo = 1, nlat_o
       do io = 1, nlon_o
          if (abs( lat_o(io,jo)- lat_i(io,jo)) > relerr .or. &
              abs( lon_o(io,jo)- lon_i(io,jo)) > relerr .or. &
              abs(lone_o(io,jo)-lone_i(io,jo)) > relerr .or. &
              abs(lonw_o(io,jo)-lonw_i(io,jo)) > relerr .or. &
              abs(latn_o(io,jo)-latn_i(io,jo)) > relerr .or. &
              abs(lats_o(io,jo)-lats_i(io,jo)) > relerr) then
             overlapgrid = .false.
          endif
       enddo
       enddo
    endif

    !--- latlongrid means both grid are regular latlon grids
    !--- assume true and then set false if not
    latlongrid = .true.
    do jo = 1,nlat_o
    do io = 1,nlon_o
       if (abs( lat_o(io,jo) -  lat_o(1,jo)) > relerr .or. &
           abs(latn_o(io,jo) - latn_o(1,jo)) > relerr .or. &
           abs(lats_o(io,jo) - lats_o(1,jo)) > relerr) then
           latlongrid = .false.
       endif
       if (abs( lon_o(io,jo) -  lon_o(io,1)) > relerr .or. &
           abs(lone_o(io,jo) - lone_o(io,1)) > relerr .or. &
           abs(lonw_o(io,jo) - lonw_o(io,1)) > relerr) then
           latlongrid = .false.
       endif
    enddo
    enddo
    do ji = 1,nlat_i
    do ii = 1,nlon_i
       if (abs( lat_i(ii,ji) -  lat_i(1,ji)) > relerr .or. &
           abs(latn_i(ii,ji) - latn_i(1,ji)) > relerr .or. &
           abs(lats_i(ii,ji) - lats_i(1,ji)) > relerr) then
           latlongrid = .false.
       endif
       if (abs( lon_i(ii,ji) -  lon_i(ii,1)) > relerr .or. &
           abs(lone_i(ii,ji) - lone_i(ii,1)) > relerr .or. &
           abs(lonw_i(ii,ji) - lonw_i(ii,1)) > relerr) then
           latlongrid = .false.
       endif
    enddo
    enddo

    if (masterproc) write(6,*) 'setmapsFM overlapgrid,latlongrid = ',overlapgrid,latlongrid

    n_i2o = 0
    w_i2o = 0.0_r8
    do jo = 1, nlat_o
    do io = 1, nlon_o
    if (pftm_o(io,jo) >= 0) then       ! only real or fake points
       found  = .false.
       lat_o_local = min(max(lat_o(io,jo),-90.0_r8),90.0_r8)  !limit [-90,90]

       if (overlapgrid) then
          found = .true.
          if = io
          jf = jo
       elseif (latlongrid) then
          do ji = 1,nlat_i
             offset = n*doffset
             if ((ji == 1 .and. lat_o_local <= latn_i(1,ji) .and.   &
                                lat_o_local >= lats_i(1,ji)) .or.   &
                 (ji >  1 .and. lat_o_local <= latn_i(1,ji) .and.   &
                                lat_o_local >  lats_i(1,ji))) then
                if (found) then
                   write(6,*) 'gridmap_setmapsFM WARNING: found > 1 pt j ', &
                      io,jo,ji,jf,lon_o(io,jo),lat_o(io,jo),lat_o_local
                   call endrun()
                endif
                jf = ji
                found = .true.
             endif
          enddo  ! ji
          if (found) then     ! move on to i
             found = .false.
          else                ! stop
             write(6,*) 'gridmap_setmapsFM ERROR: pt not found, ', &
                io,jo,lon_o(io,jo),lat_o(io,jo),lat_o_local, '_o', &
                minval(lon_o),maxval(lon_o),           &
                minval(lat_o),maxval(lat_o),'_iwe',    &
                minval(lonw_i),maxval(lonw_i),         &
                minval(lone_i),maxval(lone_i),'_isn',  &
                minval(lats_i),maxval(lats_i),         &
                minval(latn_i),maxval(latn_i)
             call endrun()
          endif
          do ii = 1,nlon_i
          do n = -noffset,noffset
             offset = n*doffset
             if (lon_o(io,jo)+offset <= lone_i(ii,jf) .and.   &
                 lon_o(io,jo)+offset >= lonw_i(ii,jf)) then
                if (found) then
                   write(6,*) 'gridmap_setmapsFM WARNING: found > 1 pt i ', &
                      io,jo,lon_o(io,jo),lat_o(io,jo),lat_o_local
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
          do ji = 1,nlat_i
          do ii = 1,nlon_i
             if (lon_o(io,jo)+offset <= lone_i(ii,ji) .and.   &
                 lon_o(io,jo)+offset >= lonw_i(ii,ji) .and.   &
                 lat_o_local         <= latn_i(ii,ji) .and.   &
                 lat_o_local         >= lats_i(ii,ji)) then
                if (found) then
                   write(6,*) 'gridmap_setmapsFM WARNING: found > 1 pt', &
                      io,jo,lon_o(io,jo),lat_o(io,jo),lat_o_local
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

       if (found) then
          n_i2o(io,jo) = n_i2o(io,jo) + 1
          if (n_i2o(io,jo) > mx_i2o) then
            write(6,*) 'gridmap_setmapsFM ERROR: n_i2o > mx_i2o ', &
               io,jo,lon_o(io,jo),lat_o(io,jo),lat_o_local,mx_i2o
            call endrun()
          endif
          i_i2o(io,jo,n_i2o(io,jo)) = if
          j_i2o(io,jo,n_i2o(io,jo)) = jf
       else
          write(6,*) 'gridmap_setmapsFM ERROR: pt not found, ', &
             io,jo,lon_o(io,jo),lat_o(io,jo),lat_o_local, '_o', &
             minval(lon_o),maxval(lon_o),           &
             minval(lat_o),maxval(lat_o),'_iwe',    &
             minval(lonw_i),maxval(lonw_i),         &
             minval(lone_i),maxval(lone_i),'_isn',  &
             minval(lats_i),maxval(lats_i),         &
             minval(latn_i),maxval(latn_i)
          call endrun()
       endif

    endif
    enddo
    enddo

    !--- find aggregation overlap number (ncnta) for each ii,ji.
    !--- o2i is derived from i2o indices

    allocate(ncnta(nlon_i,nlat_i))
    ncnta = 0

    do jo = 1, nlat_o
    do io = 1, nlon_o
    do n = 1,n_i2o(io,jo)
       ii = i_i2o(io,jo,n)
       ji = j_i2o(io,jo,n)
       ncnta(ii,ji) = ncnta(ii,ji) + 1
    enddo
    enddo
    enddo

    !--- initialize gridmap_o2i ---

    mx_o2i = maxval(ncnta)
    lname = 'setmapsFM_l2a'
    if (present(name)) then
      lname = trim(name)//'_l2a'
    endif
    call gridmap_init(gridmap_o2i,domain_o,domain_i,mx_o2i,name=lname, &
       type=gridmap_typelocal)
    call gridmap_setptrs(gridmap_o2i,mx_ovr=mx_o2i,n_ovr=n_o2i,i_ovr=i_o2i, &
       j_ovr=j_o2i,w_ovr=w_o2i)

    !--- set *_o2i gridmap  ---

    n_o2i = 0
    w_o2i = 0.0_r8

    do jo = 1, nlat_o
    do io = 1, nlon_o
    do n = 1,n_i2o(io,jo)
       ii = i_i2o(io,jo,n)
       ji = j_i2o(io,jo,n)
       n_o2i(ii,ji) = n_o2i(ii,ji) + 1
       if (n_o2i(ii,ji) > mx_o2i .or. n_o2i(ii,ji) > ncnta(ii,ji)) then
          write(6,*) 'gridmap_setmapsFM ERROR: n_o2i > mx_o2i ', &
             ii,ji,lon_i(ii,ji),lat_i(ii,ji),n_o2i(ii,ji),mx_o2i,ncnta(ii,ji)
          call endrun()
       endif
       i_o2i(ii,ji,n_o2i(ii,ji)) = io
       j_o2i(ii,ji,n_o2i(ii,ji)) = jo
    enddo
    enddo
    enddo

    !--- remove fake land if there is at least one real land overlap point ---

    do ji = 1, nlat_i
    do ii = 1, nlon_i
       nold = n_o2i(ii,ji)
       found = .false.        ! found real land point
       do n = 1,nold
          if (pftm_o(i_o2i(ii,ji,n),j_o2i(ii,ji,n)) > 0 ) found = .true.
       enddo
       if (found) then        ! keep only real land points
          n_o2i(ii,ji) = 0
          do n = 1,nold
             if (pftm_o(i_o2i(ii,ji,n),j_o2i(ii,ji,n)) > 0 ) then
                n_o2i(ii,ji) = n_o2i(ii,ji) + 1
                i_o2i(ii,ji,n_o2i(ii,ji)) = i_o2i(ii,ji,n)
                j_o2i(ii,ji,n_o2i(ii,ji)) = j_o2i(ii,ji,n)
             endif
          enddo
       endif
    enddo
    enddo

    !--- set weights ---

    nara_o = 0.0_r8
    do ji = 1, nlat_i
    do ii = 1, nlon_i
       if (n_o2i(ii,ji) == 1) then
          w_o2i(ii,ji,1) = 1.0_r8
          io = i_o2i(ii,ji,1)
          jo = j_o2i(ii,ji,1)
          nara_o(io,jo) = area_i(ii,ji)
       else
          sum = 0.0_r8
          do n = 1,n_o2i(ii,ji)
             io = i_o2i(ii,ji,n)
             jo = j_o2i(ii,ji,n)
             sum = sum + area_o(io,jo)
             if (area_o(io,jo) <= 0.0_r8) then
                write(6,*) 'gridmap_setmapsFM ERROR: area_o <= 0., ', &
                  ii,ji,n,io,jo,lon_o(io,jo),lat_o(io,jo),area_o(io,jo)
                call endrun()          
             endif
          enddo
          do n = 1,n_o2i(ii,ji)
             io = i_o2i(ii,ji,n)
             jo = j_o2i(ii,ji,n)
             if (sum <= 0.0_r8) then
                write(6,*) 'gridmap_setmapsFM ERROR: sum <= 0., ', &
                   ii,ji,area_o(io,jo),sum
                call endrun()          
             endif
             w_o2i(ii,ji,n) = area_o(io,jo)/sum
             nara_o(io,jo) = (area_o(io,jo)/sum)*area_i(ii,ji)
          enddo
       endif
    enddo
    enddo

    do jo = 1, nlat_o
    do io = 1, nlon_o
       if (n_i2o(io,jo) == 1) then
          w_i2o(io,jo,1) = 1.0_r8
       else
          sum = 0.0_r8
          do n = 1,n_i2o(io,jo)
             ii = i_i2o(io,jo,n)
             ji = j_i2o(io,jo,n)
             sum = sum + area_i(ii,ji)
             if (area_i(ii,ji) <= 0.0_r8) then
                write(6,*) 'gridmap_setmapsFM ERROR: area_i <= 0., ', &
                   io,jo,n,ii,ji,lon_i(ii,ji),lat_i(ii,ji),area_i(ii,ji)
                call endrun()          
             endif
          enddo
          do n = 1,n_i2o(io,jo)
             ii = i_i2o(io,jo,n)
             ji = j_i2o(io,jo,n)
             if (sum <= 0.0_r8) then
                write(6,*) 'gridmap_setmapsFM ERROR: sum <= 0., ', &
                   ii,ji,area_o(io,jo),sum
                call endrun()          
             endif
             w_i2o(io,jo,n) = area_i(ii,ji)/sum
          enddo
       endif
    enddo
    enddo

    !--- check that valid fine grid points have coarse mapping gridpoints
    found = .false.
    do ji = 1, nlat_i
    do ii = 1, nlon_i
       if (mask_i(ii,ji) /= 0 .and. n_o2i(ii,ji) < 1) then
          write(6,*) 'gridmap_setmapsFM ERROR: invalid f->c index ', &
             ii,ji,mask_i(ii,ji),n_o2i(ii,ji)
          found = .true.
       endif
    enddo
    enddo
    if (found) call endrun()

    !--- clean up ---

    deallocate(ncnta)

    call gridmap_checkmap(gridmap_i2o)
    call gridmap_checkmap(gridmap_o2i)

end subroutine gridmap_setmapsFM

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
    real(r8), intent(in),optional,target :: fracin(:,:)
    real(r8), intent(in),optional,target :: fracout(:,:)
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
    real(r8),pointer :: area_i(:,:)  !input grid: cell area
    real(r8),pointer :: fland_i(:,:) !input grid: cell frac
    real(r8),pointer :: lone_i(:,:)  !input grid: longitude, E edge (degrees)
    real(r8),pointer :: lonw_i(:,:)  !input grid: longitude, W edge (degrees)
    real(r8),pointer :: latn_i(:,:)  !input grid: latitude, N edge (degrees)
    real(r8),pointer :: lats_i(:,:)  !input grid: latitude, S edge (degrees)
    integer          :: nlon_o       !output grid: max number of longitude pts
    integer          :: nlat_o       !output grid: number of latitude  points
    real(r8),pointer :: area_o(:,:)  !output grid: cell area
    real(r8),pointer :: fland_o(:,:) !output grid: cell frac
    real(r8),pointer :: lone_o(:,:)  !output grid: longitude, E edge  (degrees)
    real(r8),pointer :: lonw_o(:,:)  !output grid: longitude, W edge  (degrees)
    real(r8),pointer :: latn_o(:,:)  !output grid: latitude, N edge (degrees)
    real(r8),pointer :: lats_o(:,:)  !output grid: latitude, S edge (degrees)

    real(r8),allocatable :: fld_o(:,:) !output grid: dummy field
    real(r8),allocatable :: fld_i(:,:) !input grid: dummy field
    real(r8) :: sum_fldo               !global sum of dummy output field
    real(r8) :: sum_fldi               !global sum of dummy input field
    real(r8) :: relerr = 0.00001_r8    !relative error for error checks
    integer  :: mwts                   !max number of wts per cell in map
    integer  :: ii                     !input  grid longitude loop index
    integer  :: ji                     !input  grid latitude  loop index
    integer  :: io                     !output grid longitude loop index
    integer  :: jo                     !output grid latitude  loop index
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
    call gridmap_init(gridmap,domain_i,domain_o,mwts,name=lname,type=gridmap_typeglobal)

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

    dx_i = lone_i(nlon_i,1) - lonw_i(1,1)
    dx_o = lone_o(nlon_o,1) - lonw_o(1,1)

    if (latn_i(1,nlat_i) > latn_i(1,1)) then      !South to North grid
       dy_i = latn_i(1,nlat_i) - lats_i(1,1)
    else                                      !North to South grid
       dy_i = latn_i(1,1) - lats_i(1,nlat_i)
    end if
    if (latn_o(1,nlat_o) > latn_o(1,1)) then      !South to North grid
       dy_o = latn_o(1,nlat_o) - lats_o(1,1)
    else                                      !North to South grid
       dy_o = latn_o(1,1) - lats_o(1,nlat_o)
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
          fld_i(ii,ji) = ((ji-1)*nlon_i + ii) * fland_i(ii,ji)
          sum_fldi = sum_fldi + area_i(ii,ji)*fld_i(ii,ji)
       end do
    end do

    ! area-average output field from input field

    call areaave (fld_i , fld_o , gridmap)

    ! global sum of output field -- must multiply by fraction of output
    ! grid that is land as determined by input grid

    sum_fldo = 0._r8
    do jo = 1, nlat_o
       do io = 1, nlon_o
          sum_fldo = sum_fldo + area_o(io,jo)*fld_o(io,jo) * fland_o(io,jo)
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
    real(r8),optional ,intent(in) :: scale_i(:,:) !input scale field
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
    integer jo            !latitude index for output grid
    integer io            !longitude index for output grid
    integer ji            !latitude index for input grid
    integer ii            !longitude index for input grid
    integer n             !overlapping cell index
!------------------------------------------------------------------------
!dir$ inlinenever areaave

    call domain_setptrs(gridmap%domain_i,ni=nlon_i,nj=nlat_i)
    call domain_setptrs(gridmap%domain_o,ni=nlon_o,nj=nlat_o)
    call gridmap_setptrs(gridmap,mx_ovr=mx_ovr,n_ovr=n_ovr,i_ovr=i_ovr,j_ovr=j_ovr,w_ovr=w_ovr)

    if (trim(gridmap%type) /= trim(gridmap_typeglobal)) then
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
                if (present(scale_i)) then
                   fld_o(io,jo) = fld_o(io,jo) + &
                                  w_ovr(io,jo,n)*fld_i(ii,ji)*scale_i(ii,ji)
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

    real(r8),intent(in) :: fland_i(:,:)   ! input grid : mask (0, 1)
    real(r8),intent(in) :: fland_o(:,:)   ! output grid: fraction that is land
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
    integer :: io                   !output grid longitude loop index
    integer :: ii                   !input  grid longitude loop index
    integer :: jo                   !output grid latitude  loop index
    integer :: ji                   !input  grid latitude  loop index
    integer :: n                    !weights index loop
    real(r8) :: f_ovr               !sum of overlap weights
    real(r8) :: relerr = 0.00001_r8 !max error: sum overlap weights ne 1
    real(r8) :: dx_i                !input grid  longitudinal range
    real(r8) :: dy_i                !input grid  latitudinal  range
    real(r8) :: dx_o                !output grid longitudinal range
    real(r8) :: dy_o                !output grid latitudinal  range
    integer  :: nlon_i              !input size, i
    integer  :: nlat_i              !input size, j
    real(r8),pointer :: lone_i(:,:) !input grid: longitude, E edge (degrees)
    real(r8),pointer :: lonw_i(:,:) !input grid: longitude, W edge (degrees)
    real(r8),pointer :: latn_i(:,:) !input grid: latitude, N edge (degrees)
    real(r8),pointer :: lats_i(:,:) !input grid: latitude, S edge (degrees)
    integer  :: nlon_o              !output size, i
    integer  :: nlat_o              !output size, j
    real(r8),pointer :: lone_o(:,:) !output grid: longitude, E edge (degrees)
    real(r8),pointer :: lonw_o(:,:) !output grid: longitude, W edge (degrees)
    real(r8),pointer :: latn_o(:,:) !output grid: latitude, N edge (degrees)
    real(r8),pointer :: lats_o(:,:) !output grid: latitude, S edge (degrees)
    real(r8),pointer :: area_o(:,:) ! output grid: cell area
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

          ! find total land area of overlapping input cells

          f_ovr = 0._r8
          do n = 1, n_ovr(io,jo)
             ii = i_ovr(io,jo,n)
             ji = j_ovr(io,jo,n)
             f_ovr = f_ovr + w_ovr(io,jo,n)*fland_i(ii,ji)
          end do

          ! make sure area of overlap is less than or equal to output grid cell area

          if ((f_ovr-area_o(io,jo))/area_o(io,jo) > relerr) then
             write (6,*) 'AREAMAP error: area not conserved for lon,lat = ',io,jo
             write (6,'(a30,e20.10)') 'sum of overlap area = ',f_ovr
             write (6,'(a30,e20.10)') 'area of output grid = ',area_o(io,jo)
             call endrun
          end if

          ! make weights

          do n = 1, n_ovr(io,jo)
             ii = i_ovr(io,jo,n)
             ji = j_ovr(io,jo,n)
             if (f_ovr > 0._r8) then
                w_ovr(io,jo,n) = w_ovr(io,jo,n)*fland_i(ii,ji) / f_ovr
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

    dx_i = lone_i(nlon_i,1) - lonw_i(1,1)
    dx_o = lone_o(nlon_o,1) - lonw_o(1,1)

    if (latn_i(1,nlat_i) > latn_i(1,1)) then      !South to North grid
       dy_i = latn_i(1,nlat_i) - lats_i(1,1)
    else                                      !North to South grid
       dy_i = latn_i(1,1) - lats_i(1,nlat_i)
    end if
    if (latn_o(1,nlat_o) > latn_o(1,1)) then      !South to North grid
       dy_o = latn_o(1,nlat_o) - lats_o(1,1)
    else                                      !North to South grid
       dy_o = latn_o(1,1) - lats_o(1,nlat_o)
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
          f_ovr = 0._r8

          do n = 1, mx_ovr
             f_ovr = f_ovr + w_ovr(io,jo,n)
          end do

          ! error check only valid if output grid cell has land. non-land cells
          ! will have weights equal to zero

          if (fland_o(io,jo) > 0._r8) then
             if (abs(f_ovr-1._r8) > relerr) then
                write (6,*) 'AREAMAP error: area not conserved for lon,lat = ',io,jo
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
    real(r8), pointer :: lone_i(:,:)  !input grid : cell E edge lon (deg)
    real(r8), pointer :: lonw_i(:,:)  !input grid : cell W edge lon (deg)
    real(r8), pointer :: latn_i(:,:)  !input grid : cell N edge lat (deg)
    real(r8), pointer :: lats_i(:,:)  !input grid : cell S edge lat (deg)
    integer           :: nlon_o       !output grid: max lon points
    integer           :: nlat_o       !output grid: lat points
    real(r8), pointer :: lone_o(:,:)  !output grid: cell E edge lon (deg)
    real(r8), pointer :: lonw_o(:,:)  !output grid: cell W edge lon (deg)
    real(r8), pointer :: latn_o(:,:)  !output grid: cell N edge lat (deg)
    real(r8), pointer :: lats_o(:,:)  !output grid: cell S edge lat (deg)

    integer, parameter :: mx_ovr_ceiling = 100000 !emergency limit
    integer io             !output grid longitude loop index
    integer jo             !output grid latitude  loop index
    integer indexo         !output grid lat. index according to orientn
    integer ii             !input  grid longitude loop index
    integer ji             !input  grid latitude  loop index
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
    integer ,allocatable :: n_ovrl(:,:)    ! local copy of novr
    real(r8),allocatable :: lonw_il(:,:)  ! local copy of lonw_i with offset
    real(r8),allocatable :: lone_il(:,:)  ! local copy of lone_i with offset
    integer ier            ! error flag
!------------------------------------------------------------------------
 
    call domain_setptrs(domain_i,ni=nlon_i,nj=nlat_i, &
                        lats=lats_i,latn=latn_i,lonw=lonw_i,lone=lone_i)
    call domain_setptrs(domain_o,ni=nlon_o,nj=nlat_o, &
                        lats=lats_o,latn=latn_o,lonw=lonw_o,lone=lone_o)

    allocate(n_ovrl(nlon_o,nlat_o),lonw_il(nlon_i,nlat_i), &
                                  lone_il(nlon_i,nlat_i),stat=ier)
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

       if (lonw_i(1,1) < lonw_o(1,1)) then
          offset = (n*360)
       else
          offset = -(n*360)
       end if
       do ji = 1, nlat_i
       do ii = 1, nlon_i
          lonw_il(ii,ji) = lonw_i(ii,ji) + offset
          lone_il(ii,ji) = lone_i(ii,ji) + offset
       end do
       end do

       ! for all output grid cells-

       do jo = 1, nlat_o

       if (latn_o(1,nlat_o) > latn_o(1,1)) then
          indexo  = jo          !south to north at the center of cell
       else
          indexo  = nlat_o+1-jo !north to south at the center of cell
       end if

       do io = 1, nlon_o

          ! loop through all input grid cells to find overlap with output grid

          do ji = 1, nlat_i

          if (latn_i(1,nlat_i) > latn_i(1,1)) then
             indexi  = ji          !south to north at the center of cell
          else
             indexi  = nlat_i+1-ji !north to south at the center of cell
          end if

          ! lats overlap

          if ( lats_i(1,indexi)<latn_o(1,indexo) .and. &
               latn_i(1,indexi)>lats_o(1,indexo) ) then

          do ii = 1, nlon_i

             ! lons overlap

             if (lonw_il(ii,indexi)<lone_o(io,indexo) .and. &
                lone_il(ii,indexi)>lonw_o(io,indexo)) then

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

                   lone = min(lone_o(io,indexo),lone_il(ii,indexi))*deg2rad 
                   lonw = max(lonw_o(io,indexo),lonw_il(ii,indexi))*deg2rad 
                   dx = max(0.0_r8,(lone-lonw))
                   latn = min(latn_o(io,indexo),latn_i(ii,indexi))*deg2rad 
                   lats = max(lats_o(io,indexo),lats_i(ii,indexi))*deg2rad 
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
    real(r8), pointer :: lats(:,:)
    real(r8), pointer :: latn(:,:)
    real(r8), pointer :: lonw(:,:)
    real(r8), pointer :: lone(:,:)
    real(r8), pointer :: area(:,:)  
    integer i,j                 !indices
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
    do j = 1, nlat
       do i = 1, nlon
          global = global + area(i,j)
       end do
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
    real(r8), pointer :: lats(:,:)
    real(r8), pointer :: latn(:,:)
    real(r8), pointer :: lonw(:,:)
    real(r8), pointer :: lone(:,:)
    real(r8), pointer :: area(:,:)  
    integer i,j                 !indices
    real(r8) deg2rad            !pi/180
    real(r8) dx                 !cell width: E-W
    real(r8) dy                 !cell width: N-S
!------------------------------------------------------------------------

    ! Note: supports general lat/lon grids

    call domain_setptrs(domain,ni=nlon,nj=nlat,area=area, &
                        lats=lats,latn=latn,lonw=lonw,lone=lone)

    deg2rad = (SHR_CONST_PI) / 180._r8
    do j = 1, nlat
       do i = 1, nlon
          dx = (lone(i,j) - lonw(i,j)) * deg2rad
          dy = sin(latn(i,j)*deg2rad) - sin(lats(i,j)*deg2rad) 
          area(i,j) = dx*dy*re*re
       end do
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
    real(r8),pointer :: lonc(:,:) 
    real(r8),pointer :: latc(:,:) 
    real(r8),pointer :: lats(:,:)   
    real(r8),pointer :: latn(:,:)   
    real(r8),pointer :: lonw(:,:)   
    real(r8),pointer :: lone(:,:)   
    integer i,j             !indices
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
       lats(:,1)    = edges
       latn(:,nlat) = edgen
    elseif (latc(1,2) > latc(1,1)) then  ! South to North grid
       lats(:,1)    = edges
       latn(:,nlat) = edgen
       do j = 2, nlat
          lats(:,j) = (latc(1,j-1) + latc(1,j)) / 2._r8
          latn(:,j-1) = lats(:,j)
       end do
    else                                     ! North to South grid
       latn(:,1)    = edgen
       lats(:,nlat) = edges
       do j = 2, nlat
          latn(:,j) = (latc(1,j-1) + latc(1,j)) / 2._r8
          lats(:,j-1) = latn(:,j)
       end do
    end if

    ! Longitudes
    ! Western edge of first grid cell -- since grid starts with western
    ! edge on Dateline, lonw(1,j)=-180. This is the same as [edgew].

    do j = 1, nlat
       lonw(1,j)    = edgew
       lone(nlon,j) = edgee
       dx = (edgee - edgew) / nlon
       do i = 2, nlon
          lonw(i,j)   = lonw(1,j) + (i-1)*dx
          lone(i-1,j) = lonw(i,j)
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
    real(r8),pointer :: lonc(:,:) 
    real(r8),pointer :: latc(:,:) 
    real(r8),pointer :: lats(:,:)   
    real(r8),pointer :: latn(:,:)   
    real(r8),pointer :: lonw(:,:)   
    real(r8),pointer :: lone(:,:)   
    integer i,j             !indices
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
    lats(:,1)    = -90._r8
    latn(:,nlat) = 90._r8
    do j = 2, nlat   
       lats(:,j) = (latc(1,j-1) + latc(1,j)) / 2._r8
       latn(:,j-1) = lats(:,j)
    end do

    ! Longitudes

    if (lonc(1,1) >= 0._r8) then
       do j = 1, nlat
          dx = 360._r8/(nlon)
          lonw(1,j) = -dx/2._r8
          lone(nlon,j) = -dx/2._r8 + (nlon)*dx
          do i = 2, nlon
             lonw(i,j) = -dx/2._r8 + (i-1)*dx
             lone(i-1,j) = lonw(i,j)
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



















