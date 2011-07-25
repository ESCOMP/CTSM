
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
  use domainMod    , only : domain_type, domain_setptrs
  use shr_const_mod, only : SHR_CONST_PI, SHR_CONST_REARTH
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_sys_mod  , only : shr_sys_flush
  use nanMod       , only : nan, bigint
!
! !PUBLIC TYPES:
  implicit none
  private

  type gridmap_type
     private
     ! lower level in hierarchy
     character(len=32)         :: name
     character(len=16)         :: type               ! global, dst, src, etc
     type(domain_type),pointer :: domain_i           ! domain_i
     type(domain_type),pointer :: domain_o           ! domain_o
     integer                   :: mx_ovr             ! max num of overlapping cells
     integer          ,pointer :: n_ovr(:,:)         ! number of overlapping cells
     integer          ,pointer :: i_ovr(:,:,:)       ! i index of overlap input cell
     integer          ,pointer :: j_ovr(:,:,:)       ! j index of overlap input cell
     real(r8)         ,pointer :: a_ovr(:,:,:)       ! area of overlap input cell
     real(r8)         ,pointer :: w_ovr(:,:,:)       ! wt of overlap input cell
     real(r8)         ,pointer :: scale_pft_i(:,:)   ! PFT wt of overlap input cell
  end type gridmap_type
  public gridmap_type

  character(len=16),parameter,public :: gridmap_typelocal  = 'local'
  character(len=16),parameter,public :: gridmap_typeglobal = 'global'
  character(len=16),parameter,public :: gridmap_typedst    = 'dst'
  character(len=16),parameter,public :: gridmap_typesrc    = 'src'

  type (gridmap_type), public :: gridmap_d2l

!
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: gridmap_clean
  public :: gridmap_setptrs
  public :: areaini        ! area averaging initialization
  public :: areaave        ! area averaging of field from input to output grids
  public :: areaini_pft    ! area averaging initialization for plant function type average
  public :: areaave_pft    ! area averaging of field from input to output grids with pft
  interface areaave
     module procedure areaave
     module procedure areaave3D
     module procedure areaave4D
  end interface
  interface celledge
     module procedure celledge_regional
     module procedure celledge_global  
     module procedure celledge_global_new  
  end interface
  interface cellarea
     module procedure cellarea_regional
     module procedure cellarea_global
  end interface
  public :: celledge
  public :: cellarea

  public :: gridmap_checkmap
!
! !REVISION HISTORY:
! Created by Sam Levis
! Updated to clm2.1 data structures by Mariana Vertenstein
! 2005.11.01 Updated and cleaned by T Craig
!
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: gridmap_init
  private :: gridmap_init_pft   ! Initialize plant function type weigts
  private :: areaave_internal   ! area averaging of field from input to output grids
  private :: areamap            ! weights and indices for area of overlap between grids
  private :: areaovr            ! area of overlap between grid cells
  real(r8):: re = SHR_CONST_REARTH*0.001  ! radius of earth (km)
  logical :: masterproc=.true.
!EOP
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
!
! !LOCAL VARIABLES:
  character(len=*), parameter :: subName = "gridmap_init"
  integer ni,nj  ! size of domain_o
  integer ier    ! error flag
!EOP
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
           gridmap%j_ovr(ni,nj,mwts), gridmap%w_ovr(ni,nj,mwts), &
           gridmap%a_ovr(ni,nj,mwts), stat=ier)
  if (ier /= 0) then
     write(6,*) subName//' ERROR: allocate gridmap'
     stop
  endif

  gridmap%n_ovr = bigint
  gridmap%i_ovr = bigint
  gridmap%j_ovr = bigint
  gridmap%w_ovr = nan
  gridmap%a_ovr = nan

end subroutine gridmap_init

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gridmap_init_pft
!
! !INTERFACE:
  subroutine gridmap_init_pft(gridmap,ni,nj)
!
! !DESCRIPTION:
! This subroutine initializes the pft weights in the gridmap datatype
!
! !USES:
!
! !ARGUMENTS:
  implicit none
  type(gridmap_type), intent(inout) :: gridmap
  integer,            intent(in)    :: ni,nj                  ! size of pftin
!
! !REVISION HISTORY:
! 2007.01.29  E Kluzek  Creation.
!
!
! !LOCAL VARIABLES:
  character(len=*), parameter :: subName = "gridmap_init_pft"
  integer :: nio,njo          ! size of domain_o
  integer :: mwts             ! maximum number of overlap areas in output cell
  integer :: ier              ! error status
  logical, save :: first_time = .true.
!EOP
!------------------------------------------------------------------------------

  if ( first_time )then
     allocate( gridmap%scale_pft_i(ni,nj), stat=ier )
     if (ier /= 0) then
        write(6,*) subName//' ERROR: allocate gridmap scale_pft_i'
        stop
     endif
     nio = size(gridmap%w_ovr,1)
     njo = size(gridmap%w_ovr,2)
     call gridmap_setptrs( gridmap, mx_ovr=mwts )

     first_time = .false.
  end if
  gridmap%scale_pft_i(:,:) = nan
  gridmap%w_ovr(:,:,:)     = nan
end subroutine gridmap_init_pft

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gridmap_clean
!
! !INTERFACE:
  subroutine gridmap_clean(gridmap)
!
! !DESCRIPTION:
! This subroutine initializes the gridmap datatype
!
! !USES:
!
! !ARGUMENTS:
  implicit none
  type(gridmap_type), intent(inout)       :: gridmap
!
! !REVISION HISTORY:
! 2005.11.15  T Craig  Creation.
!
!
! !LOCAL VARIABLES:
  character(len=*), parameter :: subName = "gridmap_clean"
  integer ier    ! error flag
!EOP
!------------------------------------------------------------------------------

  nullify(gridmap%domain_i)
  nullify(gridmap%domain_o)
  gridmap%name = 'unset'
  gridmap%type = 'unset'
  gridmap%mx_ovr = bigint
  deallocate(gridmap%n_ovr, gridmap%i_ovr, &
             gridmap%j_ovr, stat=ier)
  if (ier /= 0) then
     write(6,*) SubName//' ERROR: deallocate gridmap'
     stop
  endif
  deallocate(gridmap%a_ovr, stat=ier)
  if (ier /= 0) then
     write(6,*) SubName//' ERROR: deallocate gridmap a_ovr'
     stop
  endif
  deallocate(gridmap%w_ovr, stat=ier)
  if (ier /= 0) then
     write(6,*) SubName//' ERROR: deallocate gridmap w_ovr'
     stop
  endif
  if ( associated( gridmap%scale_pft_i) )then
     deallocate(gridmap%scale_pft_i, stat=ier)
     if (ier /= 0) then
        write(6,*) SubName//' ERROR: deallocate gridmap scale_pft_i'
        stop
     endif
  end if

end subroutine gridmap_clean
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gridmap_setptrs
!
! !INTERFACE:
  subroutine gridmap_setptrs(gridmap,name,type,domain_i,domain_o, &
     mx_ovr,n_ovr,i_ovr,j_ovr,a_ovr,w_ovr,scale_pft_i)
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
    real(r8)         ,optional,pointer  :: a_ovr(:,:,:)
    real(r8)         ,optional,pointer  :: w_ovr(:,:,:)
    real(r8)         ,optional,pointer  :: scale_pft_i(:,:)
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
    character(len=*), parameter :: subName = "gridmap_setptrs"
!
!EOP
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
    if (present(a_ovr)) then
      a_ovr => gridmap%a_ovr
    endif
    if (present(w_ovr)) then
      w_ovr => gridmap%w_ovr
    endif
    if (present(scale_pft_i)) then
      if ( .not. associated(gridmap%scale_pft_i) )then
         write(6,*) subName//" ERROR:: scale_pft_i asked for but NOT allocated yet"
         stop
      end if
      scale_pft_i => gridmap%scale_pft_i
    endif

end subroutine gridmap_setptrs
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
!
! !LOCAL VARIABLES:
!EOP
    integer          :: nlon_i       !input  grid: max number of longitude pts
    integer          :: nlat_i       !input  grid: number of latitude  points
    integer          :: nlon_o       !output grid: max number of longitude pts
    integer          :: nlat_o       !output grid: number of latitude  points
    integer          :: mx_ovr       !max overlapping cells
    integer ,pointer :: n_ovr(:,:)   !lon index, overlapping input cell
    integer ,pointer :: i_ovr(:,:,:) !lon index, overlapping input cell
    integer ,pointer :: j_ovr(:,:,:) !lat index, overlapping input cell
    real(r8),pointer :: a_ovr(:,:,:) !overlap areas for input cells
    real(r8),pointer :: w_ovr(:,:,:) !overlap weights for input cells
    integer          :: i,j,n        !loop counters
    real(r8)         :: sum          !running sum
    real(r8)         :: rmin,rmax    !local min/max values
!
!------------------------------------------------------------------------

    !--- set pointers into domains ---
    call domain_setptrs(gridmap%domain_i,ni=nlon_i,nj=nlat_i)
    call domain_setptrs(gridmap%domain_o,ni=nlon_o,nj=nlat_o)
    call gridmap_setptrs(gridmap,mx_ovr=mx_ovr,n_ovr=n_ovr,i_ovr=i_ovr, &
       j_ovr=j_ovr,a_ovr=a_ovr,w_ovr=w_ovr)

    if (masterproc) then
       write(6,*) ' '
       write(6,*) 'gridmap_checkmap name         = ',trim(gridmap%name)
       write(6,*) 'gridmap_checkmap type         = ',trim(gridmap%type)
       write(6,*) 'gridmap_checkmap src grid     = ',nlon_i,nlat_i
       write(6,*) 'gridmap_checkmap dst grid     = ',nlon_o,nlat_o
       write(6,*) 'gridmap_checkmap mx_ovr        = ',mx_ovr
       write(6,*) 'gridmap_checkmap n_ovr min/max = ',minval(n_ovr),maxval(n_ovr)
       write(6,*) 'gridmap_checkmap i_ovr min/max = ',minval(i_ovr),maxval(i_ovr)
       write(6,*) 'gridmap_checkmap j_ovr min/max = ',minval(j_ovr),maxval(j_ovr)
       write(6,*) 'gridmap_checkmap a_ovr min/max = ',minval(a_ovr),maxval(a_ovr)
       write(6,*) 'gridmap_checkmap w_ovr min/max = ',minval(w_ovr),maxval(w_ovr)
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
    type(domain_type) ,intent(inout)     :: domain_o   ! output domain
    type(gridmap_type),intent(inout)     :: gridmap    ! gridmap
    real(r8), intent(in),optional,target :: fracin(:,:)
    real(r8), intent(in),optional,target :: fracout(:,:)
    character(len=*),intent(in),optional :: name
!
! !REVISION HISTORY:
! Created by Gordon Bonan
! 2005.11.20 Updated by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
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
    domain_o%frac = fland_o

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
       stop
    endif
    if (present(fracout)) then
       fland_o => fracout
    else
       write(6,*) 'areaini ERROR: fracout required'
       stop
    endif

    ! Dynamically allocate memory

    allocate (fld_o(nlon_o,nlat_o), fld_i(nlon_i,nlat_i), stat=ier)
    if (ier /= 0) then
       write (6,*) 'areaini(): allocation error'
       stop
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

    if (sum_fldi .ne. 0) then   ! NOT a null field
       if ( abs(sum_fldo/sum_fldi-1._r8) > relerr ) then
          write (6,*) 'AREAINI error: input field not conserved'
          write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
          write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
          stop
       end if
    end if

    deallocate (fld_o, fld_i)

  end subroutine areaini

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: areaini_pft
!
! !INTERFACE:
  subroutine areaini_pft( gridmap, ni, nj, pctpft_i, pft_indx )
!
! !DESCRIPTION:
! Initialize the Plant function type weights
!
! !USES:
    use mkvarctl, only: numpft
!
! !ARGUMENTS:
    type(gridmap_type),intent(inout) :: gridmap                  ! gridmap
    integer,           intent(in)    :: ni                       ! size of longitude of pft index
    integer,           intent(in)    :: nj                       ! size of latitude of pft index
    real(r8),          intent(in)    :: pctpft_i(ni,nj,0:numpft) ! plant function type %'s
    integer,           intent(in)    :: pft_indx                 ! PFT index
!
! !REVISION HISTORY:
! 2007.01.29 Created by Erik Kluek
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=*), parameter :: subName = "areaini_pft"
    integer  :: nlon_o                   !output  grid: max number of longitude pts
    integer  :: nlat_o                   !output  grid: number of latitude  points
    integer  :: nlon_i                   !input  grid: max number of longitude pts
    integer  :: nlat_i                   !input  grid: number of latitude  points
    integer  :: mwts                     !max number of wts per cell in map
    integer  :: ier                      !allocate error status
    integer  :: ii                       !input  grid longitude loop index
    integer  :: ji                       !input  grid latitude  loop index
    integer  :: io                       !output  grid longitude loop index
    integer  :: jo                       !output  grid latitude  loop index
    integer  :: n                        !weight index
    integer  :: p                        !plant function type loop index
    integer, pointer :: i_ovr(:,:,:)     !longitude index of overlap areas
    integer, pointer :: j_ovr(:,:,:)     !latitude index of overlap areas
    real(r8),pointer :: a_ovr_o(:,:,:)   !output grid: overlap areas
    real(r8),pointer :: w_ovr(:,:,:)     !output grid: overlap wts
    real(r8),pointer :: scale_pft_i(:,:) !input grid: PFT
    real(r8),pointer :: sumpft(:,:)      !input grid: sum of weights

    call domain_setptrs(gridmap%domain_i,ni=nlon_i,nj=nlat_i)
    if ( (ni /= nlon_i) .or. (nj /= nlat_i) )then
       write (6,*) subName//' error: size of input PFT does not match internal LAI grid'
       stop
    end if
    if ( (pft_indx < 0) .or. (pft_indx > numpft) )then
       write (6,*) subName//' error: pft_indx is out of range'
       stop
    end if
    call gridmap_init_pft(gridmap,ni,nj)
    call gridmap_setptrs(gridmap,mx_ovr=mwts, a_ovr=a_ovr_o,scale_pft_i=scale_pft_i, &
                         w_ovr=w_ovr, i_ovr=i_ovr,j_ovr=j_ovr )
    call domain_setptrs(gridmap%domain_o,ni=nlon_o,nj=nlat_o)
    allocate( sumpft(nlon_o,nlat_o), stat=ier )
    if (ier /= 0) then
       write(6,*) subName//' ERROR: allocate sumpft'
       stop
    endif
    do ji = 1, nlat_i
    do ii = 1, nlon_i
       scale_pft_i(ii,ji) = pctpft_i(ii,ji,pft_indx)
    end do
    end do
    sumpft(:,:) = 0.0_r8
    do n = 1, mwts
    do jo = 1, nlat_o
    do io = 1, nlon_o
       ii = i_ovr(io,jo,n)
       ji = j_ovr(io,jo,n)
       sumpft(io,jo) = sumpft(io,jo) + a_ovr_o(io,jo,n)*scale_pft_i(ii,ji)
#ifndef LINUX
       if ( scale_pft_i(ii,ji) == nan )then
          write (6,*) subName//' error: scale_pft_i == nan! at i, j=', ii,ji
          write (6,*) pctpft_i(ii,ji,pft_indx)
          stop
       end if
#endif
    end do
    end do
    end do
    ! Normalize weights by their sum (sum of weights includes a_ovr weights multiplied in)
    do n = 1, mwts
    do jo = 1, nlat_o
    do io = 1, nlon_o
       if ( sumpft(io,jo) > 0.0_r8 ) then
          w_ovr(io,jo,n) = a_ovr_o(io,jo,n)/sumpft(io,jo)
       else
          w_ovr(io,jo,n) = 0.0_r8
       end if
    end do
    end do
    end do
    if ( any(w_ovr== nan) )then
       write (6,*) subName//' error: w_ovr == nan!'
       stop
    end if
    deallocate( sumpft   )
    nullify( a_ovr_o     )
    nullify( i_ovr       )
    nullify( j_ovr       )
    nullify( w_ovr       )
    nullify( scale_pft_i )

  end subroutine areaini_pft

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: areaave_internal
!
! !INTERFACE:
  subroutine areaave_internal (fld_i , fld_o , gridmap, scale_i)
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
!
! !LOCAL VARIABLES:
!EOP
    character(len=*), parameter :: subName = "areaave_internal"
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

    call domain_setptrs(gridmap%domain_i,ni=nlon_i,nj=nlat_i)
    call domain_setptrs(gridmap%domain_o,ni=nlon_o,nj=nlat_o)
    call gridmap_setptrs(gridmap,mx_ovr=mx_ovr,n_ovr=n_ovr,i_ovr=i_ovr,j_ovr=j_ovr,w_ovr=w_ovr)

    if (trim(gridmap%type) /= trim(gridmap_typeglobal)) then
       write(6,*) subname, ' WARNING: gridmap type not global, ',gridmap%name,gridmap%type
    endif

    ! initialize field on output grid to zero everywhere

    do jo = 1, nlat_o
       do io = 1, nlon_o
          fld_o(io,jo) = 0._r8
       end do
    end do

    ! loop through overlapping cells on input grid to make area-average

    if (present(scale_i)) then
       do n = 1, mx_ovr
!$OMP PARALLEL DO PRIVATE (jo,io,ii,ji)
          do jo = 1, nlat_o
          do io = 1, nlon_o
             ii = i_ovr(io,jo,n)
             ji = j_ovr(io,jo,n)
             fld_o(io,jo) = fld_o(io,jo) + w_ovr(io,jo,n)*fld_i(ii,ji)*scale_i(ii,ji)
          end do
          end do
!$OMP END PARALLEL DO
       end do
    else
       do n = 1, mx_ovr
!$OMP PARALLEL DO PRIVATE (jo,io,ii,ji)
          do jo = 1, nlat_o
          do io = 1, nlon_o
             ii = i_ovr(io,jo,n)
             ji = j_ovr(io,jo,n)
             fld_o(io,jo) = fld_o(io,jo) + w_ovr(io,jo,n)*fld_i(ii,ji)
          end do
          end do
!$OMP END PARALLEL DO
       end do
    endif
    nullify( w_ovr )
    nullify( n_ovr )
    nullify( i_ovr )
    nullify( j_ovr )

    return
  end subroutine areaave_internal

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: areaave_internal3D
!
! !INTERFACE:
  subroutine areaave_internal3D (fld_i , fld_o , gridmap, scale_i)
!
! !DESCRIPTION:
! Mapping of field from input to output grids, 3d global fields
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8)          ,intent(in) :: fld_i(:,:,:)   !input grid : field
    real(r8)          ,intent(out):: fld_o(:,:,:)   !field for output grid
    type(gridmap_type),intent(in) :: gridmap      ! gridmap
    real(r8),optional ,intent(in) :: scale_i(:,:) !input scale field
!
! !REVISION HISTORY:
! Created by Erik Kluzek
!
!EOP
!
! LOCAL VARIABLES:
    character(len=*), parameter :: subName = "areaave_internal3D"
    integer  :: nlat_i    !input grid : number of latitude points
    integer  :: nlon_i    !input grid : max number longitude points
    integer  :: nlat_o    !output grid: number of latitude points
    integer  :: nlon_o    !output grid: max number of longitude points
    integer  :: nk        !input/output grid: max number of vertical levels
    integer          :: mx_ovr       !max overlapping cells
    integer ,pointer :: n_ovr(:,:)   !lon index, overlapping input cell
    integer ,pointer :: i_ovr(:,:,:) !lon index, overlapping input cell
    integer ,pointer :: j_ovr(:,:,:) !lat index, overlapping input cell
    real(r8),pointer :: w_ovr(:,:,:) !overlap weights for input cells
    integer jo            !latitude index for output grid
    integer k             !vertical level index for input or output grid
    integer io            !longitude index for output grid
    integer ji            !latitude index for input grid
    integer ii            !longitude index for input grid
    integer n             !overlapping cell index
!------------------------------------------------------------------------

    call domain_setptrs(gridmap%domain_i,ni=nlon_i,nj=nlat_i)
    call domain_setptrs(gridmap%domain_o,ni=nlon_o,nj=nlat_o)
    call gridmap_setptrs(gridmap,mx_ovr=mx_ovr,n_ovr=n_ovr,i_ovr=i_ovr,j_ovr=j_ovr,w_ovr=w_ovr)

    nk = size( fld_i, dim=3 )

    if (trim(gridmap%type) /= trim(gridmap_typeglobal)) then
       write(6,*) subname,' WARNING: gridmap type not global, ',gridmap%name,gridmap%type
    endif

    ! initialize field on output grid to zero everywhere

    do k  = 1, nk
    do jo = 1, nlat_o
    do io = 1, nlon_o
          fld_o(io,jo,k) = 0._r8
    end do
    end do
    end do

    ! loop through overlapping cells on input grid to make area-average

    if (present(scale_i)) then
       write(6,*) subname,' ERROR: scale input, but code does NOT handle it for a 3D field'
       stop
    end if
    do n = 1, mx_ovr
!$OMP PARALLEL DO PRIVATE (jo,io,ii,ji,k)
       do k  = 1, nk
       do jo = 1, nlat_o
       do io = 1, nlon_o
          ii = i_ovr(io,jo,n)
          ji = j_ovr(io,jo,n)
          fld_o(io,jo,k) = fld_o(io,jo,k) + w_ovr(io,jo,n)*fld_i(ii,ji,k)
       end do
       end do
       end do
!$OMP END PARALLEL DO
    end do

    nullify( w_ovr )
    nullify( n_ovr )
    nullify( i_ovr )
    nullify( j_ovr )

    return
  end subroutine areaave_internal3D

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: areaave_internal4D
!
! !INTERFACE:
  subroutine areaave_internal4D (fld_i , fld_o , gridmap)
!
! !DESCRIPTION:
! Mapping of field from input to output grids, 4d global fields
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8)          ,intent(in) :: fld_i(:,:,:,:)   !input grid : field
    real(r8)          ,intent(out):: fld_o(:,:,:,:)   !field for output grid
    type(gridmap_type),intent(in) :: gridmap          ! gridmap
!
! !REVISION HISTORY:
! Created by Erik Kluzek
!
!EOP
!
! LOCAL VARIABLES:
    character(len=*), parameter :: subName = "areaave_internal4D"
    integer  :: nlat_i    !input grid : number of latitude points
    integer  :: nlon_i    !input grid : max number longitude points
    integer  :: nlat_o    !output grid: number of latitude points
    integer  :: nlon_o    !output grid: max number of longitude points
    integer  :: nk        !input/output grid: max number of vertical levels
    integer  :: nl        !input/output grid: max number of second-vertical levels
    integer          :: mx_ovr       !max overlapping cells
    integer ,pointer :: n_ovr(:,:)   !lon index, overlapping input cell
    integer ,pointer :: i_ovr(:,:,:) !lon index, overlapping input cell
    integer ,pointer :: j_ovr(:,:,:) !lat index, overlapping input cell
    real(r8),pointer :: w_ovr(:,:,:) !overlap weights for input cells
    integer jo            !latitude index for output grid
    integer k             !vertical level index for input or output grid
    integer l             !second-vertical level index for input or output grid
    integer io            !longitude index for output grid
    integer ji            !latitude index for input grid
    integer ii            !longitude index for input grid
    integer n             !overlapping cell index
!------------------------------------------------------------------------

    call domain_setptrs(gridmap%domain_i,ni=nlon_i,nj=nlat_i)
    call domain_setptrs(gridmap%domain_o,ni=nlon_o,nj=nlat_o)
    call gridmap_setptrs(gridmap,mx_ovr=mx_ovr,n_ovr=n_ovr,i_ovr=i_ovr,j_ovr=j_ovr,w_ovr=w_ovr)

    nk = size( fld_i, dim=3 )
    nl = size( fld_i, dim=4 )

    if (trim(gridmap%type) /= trim(gridmap_typeglobal)) then
       write(6,*) subname,' WARNING: gridmap type not global, ',gridmap%name,gridmap%type
    endif

    ! initialize field on output grid to zero everywhere

    do l  = 1, nl
    do k  = 1, nk
    do jo = 1, nlat_o
    do io = 1, nlon_o
          fld_o(io,jo,k,l) = 0._r8
    end do
    end do
    end do
    end do

    ! loop through overlapping cells on input grid to make area-average

    do n = 1, mx_ovr
!$OMP PARALLEL DO PRIVATE (jo,io,ii,ji,k,l)
       do l  = 1, nl
       do k  = 1, nk
       do jo = 1, nlat_o
       do io = 1, nlon_o
          ii = i_ovr(io,jo,n)
          ji = j_ovr(io,jo,n)
          fld_o(io,jo,k,l) = fld_o(io,jo,k,l) + w_ovr(io,jo,n)*fld_i(ii,ji,k,l)
       end do
       end do
       end do
       end do
!$OMP END PARALLEL DO
    end do

    nullify( w_ovr )
    nullify( n_ovr )
    nullify( i_ovr )
    nullify( j_ovr )

    return
  end subroutine areaave_internal4D


!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: areaave
!
! !INTERFACE:
  subroutine areaave(fld_i , fld_o , gridmap)
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
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
!------------------------------------------------------------------------

    call areaave_internal (fld_i , fld_o , gridmap )

    return
  end subroutine areaave

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: areaave3D
!
! !INTERFACE:
  subroutine areaave3D(fld_i , fld_o , gridmap)
!
! !DESCRIPTION:
! Mapping of field from input to output grids, 2d global fields
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8)          ,intent(in) :: fld_i(:,:,:)   !input grid : field
    real(r8)          ,intent(out):: fld_o(:,:,:)   !field for output grid
    type(gridmap_type),intent(in) :: gridmap        ! gridmap
!
! !REVISION HISTORY:
! Created by Erik Kluzek
!
!EOP
!
! LOCAL VARIABLES:
!------------------------------------------------------------------------

    call areaave_internal3D (fld_i , fld_o , gridmap )

    return
  end subroutine areaave3D

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: areaave4D
!
! !INTERFACE:
  subroutine areaave4D(fld_i , fld_o , gridmap)
!
! !DESCRIPTION:
! Mapping of field from input to output grids, 4d global fields
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8)          ,intent(in) :: fld_i(:,:,:,:)   !input grid : field
    real(r8)          ,intent(out):: fld_o(:,:,:,:)   !field for output grid
    type(gridmap_type),intent(in) :: gridmap          ! gridmap
!
! !REVISION HISTORY:
! Created by Erik Kluzek
!
!EOP
!
! LOCAL VARIABLES:
!------------------------------------------------------------------------

    call areaave_internal4D (fld_i , fld_o , gridmap )

    return
  end subroutine areaave4D


!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: areaave_pft
!
! !INTERFACE:
  subroutine areaave_pft (fld_i , fld_o , gridmap)
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
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!
! !LOCAL VARIABLES:
    real(r8),pointer :: scale_pft_i(:,:)  !PFT %
!EOP
!------------------------------------------------------------------------

    call gridmap_setptrs(gridmap, scale_pft_i=scale_pft_i)
    call areaave_internal (fld_i , fld_o , gridmap, scale_i=scale_pft_i )
    nullify( scale_pft_i )

    return
  end subroutine areaave_pft

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
!
! !LOCAL VARIABLES:
!EOP
    integer          :: mx_ovr      ! max num input cells that overlap 
    integer ,pointer :: n_ovr(:,:)  ! number of overlapping input cells
    integer ,pointer :: i_ovr(:,:,:)! lon index, overlapping input cell
    integer ,pointer :: j_ovr(:,:,:)! lat index, overlapping input cell
    real(r8),pointer :: a_ovr(:,:,:)! overlap areas for input cells
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
       i_ovr=i_ovr,j_ovr=j_ovr,a_ovr=a_ovr,w_ovr=w_ovr)

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
             a_ovr(io,jo,n) = 0._r8
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
                  n_ovr=n_ovr, i_ovr=i_ovr, j_ovr=j_ovr, a_ovr=a_ovr  )

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

!$OMP PARALLEL DO PRIVATE (jo,io,f_ovr,n,ii,ji)
    do jo = 1, nlat_o
       do io = 1, nlon_o

          ! find total land area of overlapping input cells

          f_ovr = 0._r8
          do n = 1, n_ovr(io,jo)
             ii = i_ovr(io,jo,n)
             ji = j_ovr(io,jo,n)
             a_ovr(io,jo,n) = a_ovr(io,jo,n)*fland_i(ii,ji)
             f_ovr = f_ovr + a_ovr(io,jo,n)
          end do

          ! make sure area of overlap is less than or equal to output grid cell area

          if ((f_ovr-area_o(io,jo))/area_o(io,jo) > relerr) then
             write (6,*) 'AREAMAP error: area not conserved for lon,lat = ',io,jo
             write (6,'(a30,e20.10)') 'sum of overlap area = ',f_ovr
             write (6,'(a30,e20.10)') 'area of output grid = ',area_o(io,jo)
             stop
          end if

          ! make weights

          do n = 1, n_ovr(io,jo)
             ii = i_ovr(io,jo,n)
             ji = j_ovr(io,jo,n)
             if (f_ovr > 0._r8) then
                w_ovr(io,jo,n) = a_ovr(io,jo,n) / f_ovr
             else
                w_ovr(io,jo,n) = 0._r8
             end if
          end do

       end do
    end do
!$OMP END PARALLEL DO

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
                stop
             end if
          end if

       end do
    end do
    nullify( n_ovr, i_ovr, j_ovr, a_ovr, w_ovr, lone_i, lonw_i, latn_i, lats_i )
    nullify( lone_o, lonw_o, latn_o, lats_o, area_o )

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
                      mx_ovr , n_ovr  , i_ovr    , j_ovr , a_ovr )
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
    real(r8), intent(inout),optional :: a_ovr(:,:,:) !area of overlap pts
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
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
    integer  n             !overlapping cell index
    integer  noffsetl      !local, number of offsets to test, 0=none, default=1
    real(r8) offset        !value of offset used to shift x-grid 360 degrees
    integer ,allocatable :: n_ovrl(:,:)    ! local copy of novr
    real(r8),allocatable :: lonw_il(:,:)  ! local copy of lonw_i with offset
    real(r8),allocatable :: lone_il(:,:)  ! local copy of lone_i with offset
    logical W2E            ! Grid cells are west to east
    logical S2N            ! Grid cells are south to north
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
       stop
    end if

    if (latn_i(1,nlat_i) > latn_i(1,1)) then
       S2N = .true.             !south to north at the center of cell
    else
       S2N = .false.            !north to south at the center of cell
    end if

    if (lonw_i(1,1) < lonw_o(1,1)) then
       W2E = .true.             !west to east at the center of cell
    else
       W2E = .false.            !east to west at the center of cell
    end if

    deg2rad = SHR_CONST_PI / 180._r8
    noffsetl = 1
    if (present(noffset)) then
       noffsetl = noffset
    endif
    size3 = mx_ovr_ceiling
    if (present(a_ovr)) then
      size3 = size(a_ovr,3)
    endif
    n_ovrl(:,:) = 0

    do n = 0,noffsetl   ! loop through offsets

       if ( W2E ) then
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

       if ( S2N ) then
          indexo  = jo          !south to north at the center of cell
       else
          indexo  = nlat_o+1-jo !north to south at the center of cell
       end if

!$OMP PARALLEL DO PRIVATE (io,ii,ji,s2n,indexi,dx,dy,lone,lonw,latn,lats) SHARED (jo,indexo)
       do io = 1, nlon_o

          ! loop through all input grid cells to find overlap with output grid

          do ji = 1, nlat_i

          if ( S2N ) then
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
                   stop
                end if
                if (n_ovrl(io,indexo) > size3) then
                   write (6,*) 'AREAOVR error: n_ovr= ', &
                      n_ovrl(io,indexo),' exceeded size of arrays = ', &
                      size3,' for output lon,lat = ',io,indexo
                   stop
                end if

                if (present(i_ovr).and.present(j_ovr).and.present(a_ovr)) then
                   ! determine area of overlap

                   lone = min(lone_o(io,indexo),lone_il(ii,indexi))*deg2rad 
                   lonw = max(lonw_o(io,indexo),lonw_il(ii,indexi))*deg2rad 
                   dx = max(0.0_r8,(lone-lonw))
                   latn = min(latn_o(io,indexo),latn_i(ii,indexi))*deg2rad 
                   lats = max(lats_o(io,indexo),lats_i(ii,indexi))*deg2rad 
                   dy = max(0.0_r8,(sin(latn)-sin(lats)))

                   ! save lat, lon, area

                   i_ovr(io,indexo,n_ovrl(io,indexo)) = ii
                   j_ovr(io,indexo,n_ovrl(io,indexo)) = indexi
                   a_ovr(io,indexo,n_ovrl(io,indexo)) = dx*dy*re*re
                endif

             end if
          end do
          end if
          end do

       end do
!$OMP END PARALLEL DO
       end do

    enddo   ! offset loop

    if (present(n_ovr)) then
       n_ovr = n_ovrl
    endif
    if (present(mx_ovr)) then
       mx_ovr = maxval(n_ovrl)
    endif

    deallocate(n_ovrl,lonw_il,lone_il)
    nullify( lats_i, latn_i, lonw_i, lone_i, lats_o, latn_o, lonw_o, lone_o )

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
!
! !LOCAL VARIABLES:
!EOP
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
       stop
    end if
    nullify( lats, latn, lonw, lone, area )

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
!
! !LOCAL VARIABLES:
!EOP
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

    write(6,*) 'cellarea, using cellarea_global'

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
    nullify( lats, latn, lonw, lone, area )

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
! correspond to [longxy] (longitude at center of cell) and range from
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
!
! !LOCAL VARIABLES:
!EOP
    integer  :: nlon              
    integer  :: nlat              
    real(r8),pointer :: longxy(:,:) 
    real(r8),pointer :: latixy(:,:) 
    real(r8),pointer :: lats(:,:)   
    real(r8),pointer :: latn(:,:)   
    real(r8),pointer :: lonw(:,:)   
    real(r8),pointer :: lone(:,:)   
    integer i,j             !indices
    real(r8) dx             !cell width
!------------------------------------------------------------------------

    write(6,*) 'celledge, using celledge_regional'

    nlon = domain%ni
    nlat = domain%nj
    longxy => domain%longxy
    latixy => domain%latixy
    lats   => domain%lats
    latn   => domain%latn
    lonw   => domain%lonw
    lone   => domain%lone

    ! Latitudes
    ! Assumes lats are constant on an i line

    if (nlat == 1) then                      ! single latitude
       lats(:,1)    = edges
       latn(:,nlat) = edgen
    elseif (latixy(1,2) > latixy(1,1)) then  ! South to North grid
       lats(:,1)    = edges
       latn(:,nlat) = edgen
       do j = 2, nlat
          lats(:,j) = (latixy(1,j-1) + latixy(1,j)) / 2._r8
          latn(:,j-1) = lats(:,j)
       end do
    else                                     ! North to South grid
       latn(:,1)    = edgen
       lats(:,nlat) = edges
       do j = 2, nlat
          latn(:,j) = (latixy(1,j-1) + latixy(1,j)) / 2._r8
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
    nullify( longxy, latixy, lats, latn, lonw, lone )

    return

  end subroutine celledge_regional

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: celledge_global
!
! !INTERFACE:
  subroutine celledge_global (domain,type)
!
! !DESCRIPTION:
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout) :: domain
    character(len=*) ,intent(in)    :: type

!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 2005.11.20 Updated to domain datatype by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: nlon              
    integer  :: nlat              
    real(r8),pointer :: longxy(:,:) 
    real(r8),pointer :: latixy(:,:) 
    real(r8),pointer :: lats(:,:)   
    real(r8),pointer :: latn(:,:)   
    real(r8),pointer :: lonw(:,:)   
    real(r8),pointer :: lone(:,:)   
    integer i,j             !indices
    real(r8) dx             !cell width
!------------------------------------------------------------------------

    write(6,*) 'celledge, using celledge_global'

    nlon = domain%ni
    nlat = domain%nj
    longxy => domain%longxy
    latixy => domain%latixy
    lats   => domain%lats
    latn   => domain%latn
    lonw   => domain%lonw
    lone   => domain%lone

    ! Latitudes
    lats(:,1)    = -90._r8
    latn(:,nlat) = 90._r8
    do j = 2, nlat   
       lats(:,j) = (latixy(1,j-1) + latixy(1,j)) / 2._r8
       latn(:,j-1) = lats(:,j)
    end do

    ! Longitudes

    if (longxy(1,1) >= 0._r8) then
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
       stop
    endif
    nullify( longxy, latixy, lats, latn, lonw, lone )

    return
  end subroutine celledge_global

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: celledge_global_new
!
! !INTERFACE:
  subroutine celledge_global_new (domain)
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
! 2005.11.20 Updated to domain datatype by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: nlon              
    integer  :: nlat              
    real(r8),pointer :: longxy(:,:) 
    real(r8),pointer :: latixy(:,:) 
    real(r8),pointer :: lats(:,:)   
    real(r8),pointer :: latn(:,:)   
    real(r8),pointer :: lonw(:,:)   
    real(r8),pointer :: lone(:,:)   
    integer i,j             !indices
    real(r8) dx             !cell width
!------------------------------------------------------------------------

    write(6,*) 'celledge, using celledge_global_new'

    nlon = domain%ni
    nlat = domain%nj
    longxy => domain%longxy
    latixy => domain%latixy
    lats   => domain%lats
    latn   => domain%latn
    lonw   => domain%lonw
    lone   => domain%lone

    ! Latitudes
    lats(:,1)    = -90._r8
    latn(:,nlat) = 90._r8
    do j = 2, nlat   
    do i = 1, nlon
       lats(i,j) = (latixy(i,j-1) + latixy(i,j)) / 2._r8
       latn(i,j-1) = lats(i,j)
    enddo
    enddo

    ! Longitudes

    do j = 1, nlat
      lonw(1,j)    = 1.5_r8*longxy(1,j) - 0.5_r8*longxy(2,j)
      lone(nlon,j) = lonw(1,j) + 360._r8
    enddo

    do j = 1, nlat   
    do i = 2, nlon
      lonw(i,j) = 0.5_r8*longxy(i-1,j) + 0.5_r8*longxy(i,j)
      lone(i-1,j) = lonw(i,j)
    enddo
    enddo
    nullify( longxy, latixy, lats, latn, lonw, lone )

    return
  end subroutine celledge_global_new

!-----------------------------------------------------------------------

end module areaMod



















