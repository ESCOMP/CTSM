module areaMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: areaMod
!
! !DESCRIPTION:
! Area averaging routines. Thes routines are used for area-average 
! mapping of a field from one grid to another.
!
! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_const_mod, only : SHR_CONST_PI, SHR_CONST_REARTH
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  interface celledge
     module procedure celledge_create          !Southern and western edges of grid cells - regional grid
     module procedure celledge_input_regional  !Southern and western edges of grid cells - regional grid
     module procedure celledge_input_global    !Southern and western edges of grid cells - global grid
  end interface
  public :: areaini_point  ! area averaging initialization for single grid cell
  public :: mkmxovr        ! find maxinum numver of overlapping cells between input/output grid
!
! !REVISION HISTORY:
! Created by Sam Levis
! Updated to clm2.1 data structures by Mariana Vertenstein
!
!EOP
!
! PRIVATE MEMBER FUNCTIONS:
  private :: areamap_point ! weights and indices for area of overlap between grids for single grid cell
  private :: areaovr_point ! area of overlap between grid cells for single grid cell
  real(r8) :: re = SHR_CONST_REARTH*0.001  ! radius of earth (km)
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: celledge_create
!
! !INTERFACE:
  subroutine celledge_create(nlat   , nlon  , numlon , longxy ,  &
                             latixy , edgen , edgee  , edges  ,  &
                             edgew  , lats  , lonw, area   )
!
! !DESCRIPTION:
! Southern and western edges of grid cells of regular created grid 
! (can become global as special case)
! Latitudes -- southern/northern edges for each latitude strip.
! Grids must be oriented South to North, and West to East
! The southern and northern edges of latitude strip [j] are:
!        southern = lats(j  )
!        northern = lats(j+1)
! where [lats] must be dimensioned lats(lat+1)
! Longitudes for the western edge of the cells must increase continuously
! and span 360 degrees. Assume that grid starts at Dateline with western 
! edge on Dateline Western edges correspond to [longxy] 
! (longitude at center of cell) and range from
! -180 to 180 with negative longitudes west of Greenwich.
! Partial grids that do not span 360 degrees are allowed so long as they
! have the convention of Grid 1 with
!      western edge of grid: >= -180 and < 180
!      eastern edge of grid: > western edge  and <= 180
! [lonw] must be dimensioned lonw(lon+1,lat) because each latitude
! strip can have variable longitudinal resolution
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: nlat              ! dimension: number of latitude points
    integer , intent(in) :: nlon              ! dimension: number of longitude points
    integer , intent(in) :: numlon(nlat)      ! number of grid cells per latitude strip
    real(r8), intent(in) :: longxy(nlon,nlat) ! longitude at center of grid cell
    real(r8), intent(in) :: latixy(nlon,nlat) ! latitude at center of grid cell
    real(r8), intent(in) :: edgen             ! northern edge of grid (degrees)
    real(r8), intent(in) :: edgee             ! eastern edge of grid (degrees)
    real(r8), intent(in) :: edges             ! southern edge of grid (degrees)
    real(r8), intent(in) :: edgew             ! western edge of grid (degrees)
    real(r8), intent(out):: lats(nlat+1)      ! grid cell latitude, southern edge (degrees)
    real(r8), intent(out):: lonw(nlon+1,nlat) ! grid cell longitude, western edge (degrees)
    real(r8), intent(out), optional :: area(nlon,nlat) ! cell area (km**2)
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: i,j             ! indices
    real(r8) :: deg2rad         ! pi/180
    real(r8) :: global          ! summed area
    real(r8) :: dx              ! cell width: E-W
    real(r8) :: dy              ! cell width: N-S
    real(r8) :: error           ! true area for error check
    character(len=256) :: subname = 'celledge_create_regional'
!------------------------------------------------------------------------

    ! Latitudes

    if (nlat == 1) then
       lats(1) = edges
       lats(nlat+1) = edgen
    else
       if (latixy(1,2) > latixy(1,1)) then    !South to North grid
          lats(1) = edges
          lats(nlat+1) = edgen
       else                                   !North to South grid
          lats(1) = edgen
          lats(nlat+1) = edges
       end if
    end if
    do j = 2, nlat
       lats(j) = (latixy(1,j-1) + latixy(1,j)) / 2.
    end do

    ! Longitudes

    ! Western edge of first grid cell -- since grid starts with western
    ! edge on Dateline, lonw(1,j)=-180. This is the same as [edgew].
    ! Remaining grid cells. On a global grid lonw(nlon+1,j)=lonw(1,j)+360.
    ! This is the same as [edgee].  Set unused longitudes to non-valid number

    do j = 1, nlat
       dx = (edgee - edgew) / nlon
       lonw(1,j) = edgew
       do i = 2, nlon+1
          lonw(i,j) = lonw(1,j) + (i-1)*dx
       end do
    end do

    if (present(area)) then

       ! Note: only support S->N grid

       deg2rad = (SHR_CONST_PI) / 180.
       global = 0.

       do j = 1, nlat
          if (lats(j+1) < lats(j)) then 
             write(6,*)subname,' error: north->south latitudes not supported '
             call abort()
          end if
          do i = 1, nlon
             dx = (lonw(i+1,j) - lonw(i,j)) * deg2rad
             dy = sin(lats(j+1)*deg2rad) - sin(lats(j)*deg2rad)
             area(i,j) = dx*dy*re*re
             global = global + area(i,j)
          end do
       end do

       ! make sure total area from grid cells is same as area of grid
       ! as defined by its edges
       
       dx = (edgee - edgew) * deg2rad
       dy = sin(edgen*deg2rad) - sin(edges*deg2rad)
       error =  dx*dy*re*re
       
       if (abs(global-error)/error > 0.00001) then
          write (6,*)subname,' error: correct area is ',error, &
               ' but summed area of grid cells is ',global
          call abort()
       end if
       
    end if

  end subroutine celledge_create

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: celledge_input_global
!
! !INTERFACE:
  subroutine celledge_input_global (nlat, nlon, numlon, longxy, latixy, lats, lonw, area)
!
! !DESCRIPTION:
! Southern and western edges of grid cells - global grid
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
! cells must increase continuously and span 360 degrees.
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: nlat               ! dimension: number of latitude points
    integer , intent(in) :: nlon               ! dimension: number of longitude points
    integer , intent(in) :: numlon(nlat)       ! number of grid cells per latitude strip
    real(r8), intent(in) :: longxy(nlon,nlat)  ! longitude at center of grid cell
    real(r8), intent(in) :: latixy(nlon,nlat)  ! latitude at center of grid cell
    real(r8), intent(out):: lats(nlat+1)       ! grid cell latitude, southern edge (degrees)
    real(r8), intent(out):: lonw(nlon+1,nlat)  ! grid cell longitude, western edge (degrees)
    real(r8), intent(out), optional :: area(nlon,nlat) ! cell area (km**2)
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: i,j             !indices
    real(r8) :: deg2rad            !pi/180
    real(r8) :: dx                 !cell width: E-W
    real(r8) :: dy                 !cell width: N-S
    character(len=256) :: subname = 'celledge_input_global'
!------------------------------------------------------------------------

    ! Latitudes

    do j = 1, nlat+1              !southern edges
       if (j == 1) then           !south pole
          lats(j) = -90.
       else if (j == nlat+1) then !north pole
          lats(j) = 90.
       else                       !edge = average latitude
          lats(j) = (latixy(1,j-1) + latixy(1,j)) / 2.
       end if
    end do

    ! Longitudes

    if (longxy(1,1) >= 0.) then
       do j = 1, nlat
          dx = 360./(numlon(j))
          do i = 1, numlon(j)+1
             lonw(i,j) = -dx/2. + (i-1)*dx
          end do
          do i = numlon(j)+2, nlon
             lonw(i,j) = -999.
          end do
       end do
    else
       write(6,*)'global non-regional grids currently only supported ', &
            'for grids starting at greenwich and centered on Greenwich'
       call abort()
    endif

    if (present(area)) then

       ! Only support S->N grid. Note that cannot make sure total area from grid 
       ! cells is same as area of grid as defined by its edges as in regular created
       ! grid since the edges are not all at the same longitudes for every latitude
       ! for the possible case of a global reduced grid
       
       deg2rad = (SHR_CONST_PI) / 180.
       do j = 1, nlat
          if (lats(j+1) < lats(j)) then 
             write(6,*)subname,' error: north->south latitudes not supported '
             call abort()
          end if
          do i = 1, numlon(j)
             dx = (lonw(i+1,j) - lonw(i,j)) * deg2rad
             dy = sin(lats(j+1)*deg2rad) - sin(lats(j)*deg2rad)  
             area(i,j) = dx*dy*re*re
          end do
       end do

    end if

  end subroutine celledge_input_global

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: celledge_input_regional
!
! !INTERFACE:
  subroutine celledge_input_regional(nlat, nlon, longxy, latixy, lats, lonw, area)
!
! !DESCRIPTION:
! Southern and western edges of grid cells - regional grid
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
! cells must increase continuously and span 360 degrees.
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: nlat               ! dimension: number of latitude points
    integer , intent(in) :: nlon               ! dimension: number of longitude points
    real(r8), intent(in) :: longxy(nlon,nlat)  ! longitude at center of grid cell
    real(r8), intent(in) :: latixy(nlon,nlat)  ! latitude at center of grid cell
    real(r8), intent(out):: lats(nlat+1)       ! grid cell latitude, southern edge (degrees)
    real(r8), intent(out):: lonw(nlon+1,nlat)  ! grid cell longitude, western edge (degrees)
    real(r8), intent(out), optional :: area(nlon,nlat) ! cell area (km**2)
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: i,j                ! indices
    real(r8) :: deg2rad            ! pi/180
    real(r8) :: dx                 ! cell width: E-W
    real(r8) :: dy                 ! cell width: N-S
    character(len=256) :: subname = 'celledge_input_regional'
!------------------------------------------------------------------------

    ! southern edges of gridcell

    lats(1) = latixy(1,1) - (latixy(1,2) - latixy(1,1))/2.
    do j = 2, nlat
       lats(j) = (latixy(1,j-1) + latixy(1,j)) / 2.
    end do
    lats(nlat+1) = latixy(1,nlat) + (latixy(1,nlat) - latixy(1,nlat-1))/2.

    ! western edges of gridcell

    do j = 1,nlat
       lonw(1,j) = longxy(1,j) - (longxy(2,j) - longxy(1,j))/2.
       do i = 2,nlon
          lonw(i,j) = (longxy(i,j) + longxy(i-1,j)) / 2.
       end do
       lonw(nlon+1,j) = longxy(nlon,j) + (longxy(nlon,j) - longxy(nlon-1,j))/2.
    end do

    if (present(area)) then

       ! Only support S->N grid.
       
       deg2rad = (SHR_CONST_PI) / 180.
       do j = 1, nlat
          if (lats(j+1) < lats(j)) then 
             write(6,*)subname,' error: north->south latitudes not supported '
             call abort()
          end if
          do i = 1, nlon
             dx = (lonw(i+1,j) - lonw(i,j)) * deg2rad
             dy = sin(lats(j+1)*deg2rad) - sin(lats(j)*deg2rad)  
             area(i,j) = dx*dy*re*re
          end do
       end do

    end if

  end subroutine celledge_input_regional

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: areaini_point
!
! !INTERFACE:
  subroutine areaini_point (io      , jo          , nlon_i  , nlat_i  , &
                            numlon_i, lon_i   , lon_i_offset, lat_i   , &
                            area_i  , mask_i  , nlon_o  , nlat_o      , &
                            numlon_o, lon_o   , lat_o   , area_o      , &
                            fland_o , novr_i2o, iovr_i2o, jovr_i2o    , &
                            wovr_i2o, maxovr)
!
! !DESCRIPTION:
! area averaging initialization
! This subroutine is used for area-average mapping of a field from one
! grid to another.
!
!    areaini_point  - initializes indices and weights for area-averaging from
!                     input grid to output grid
!    areamap_point  - called by areaini_point: finds indices and weights
!    areaovr_point  - called by areamap_point: finds if cells overlap and area of overlap
!
! To map from one grid to another, must first call areaini to build
! the indices and weights (iovr_i2o, jovr_i2o, wovr_i2o). Then must
! call areaave to get new field on output grid.
!
! Not all grid cells on the input grid will be used in the area-averaging
! of a field to the output grid. Only input grid cells with [mask_i] = 1
! contribute to output grid cell average. If [mask_i] = 0, input grid cell
! does not contribute to output grid cell. This distinction is not usually
! required for atm -> land mapping, because all cells on the atm grid have
! data. But when going from land -> atm, only land grid cells have data.
! Non-land grid cells on surface grid do not have data. So if output grid cell
! overlaps with land and non-land cells (input grid), can only use land
! grid cells when computing area-average.
! o Input and output grids can be ANY resolution BUT:
!   a. Grid orientation -- Grids can be oriented south to north
!      (i.e. cell(lat+1) is north of cell(lat)) or from north to
!      south (i.e. cell(lat+1) is south of cell(lat)). Both grids must be
!      oriented from west to east, i.e., cell(lon+1) must be east of cell(lon)
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
!   c. Both grids can have variable number of longitude points for each
!      latitude strip. However, the western edge of the first point in each
!      latitude must be the same for all latitudes. Likewise, for the
!      eastern edge of the last point. That is, each latitude strip must span
!      the same longitudes, but the number of points to do this can be different
!   d. One grid can be a sub-set (i.e., smaller domain) than the other grid.
!      In this way, an atmospheric dataset for the entire globe can be
!      used in a simulation for a region 30N to 50N and 130W to 70W -- the
!      code will extract the appropriate data. The two grids do not have to
!      be the same resolution. Area-averaging will work for full => partial
!      grid but obviously will not work for partial => full grid.
! o Field values fld_i on an  input grid with dimensions nlon_i and nlat_i =>
!   field values fld_o on an output grid with dimensions nlon_o and nlat_o as
!
!   fld_o(io,jo) =
!   fld_i(i_ovr(io,jo,    1),j_ovr(io,jo,    1)) * w_ovr(io,jo,   1)
!                             ... + ... +
!   fld_i(i_ovr(io,jo,maxovr),j_ovr(io,jo,maxovr)) * w_ovr(io,jo,maxovr)
! o Error checks:
!   Overlap weights of input cells sum to 1 for each output cell.
!   Global sum of dummy field is conserved for input => output area-average.
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: io                     !output grid longitude index
    integer , intent(in)    :: jo                     !output grid latitude index
    integer , intent(in)    :: nlon_i                 !input  grid: max number of longitude points
    integer , intent(in)    :: nlat_i                 !input  grid: number of latitude  points
    integer , intent(in)    :: numlon_i(nlat_i)       !input  grid: number lon points at each lat
    real(r8), intent(in)    :: lon_i(nlon_i+1,nlat_i) !input grid: longitude, west edge (degrees)
    real(r8), intent(in)    :: lon_i_offset(nlon_i+1,nlat_i) !input grid : cell lons, west edge (deg)
    real(r8), intent(in)    :: lat_i(nlat_i+1)        !input grid: latitude, south edge (degrees)
    real(r8), intent(in)    :: area_i(nlon_i,nlat_i)  !input grid: cell area
    real(r8), intent(in)    :: mask_i(nlon_i,nlat_i)  !input  grid: mask (0, 1)
    integer , intent(in)    :: nlon_o                 !output grid: max number of longitude points
    integer , intent(in)    :: nlat_o                 !output grid: number of latitude  points
    integer , intent(in)    :: numlon_o(nlat_o)       !output grid: number lon points at each lat
    real(r8), intent(in)    :: lon_o(nlon_o+1,nlat_o) !output grid: longitude, west edge  (degrees)
    real(r8), intent(in)    :: lat_o(nlat_o+1)        !output grid: latitude, south edge (degrees)
    real(r8), intent(in)    :: area_o                 !output grid: cell area
    real(r8), intent(in)    :: fland_o                !output grid: fraction that is land
    integer , intent(out)   :: novr_i2o               !number of overlapping input cells
    integer , intent(in)    :: maxovr                 !maximum number of overlapping cells
    integer , intent(out)   :: iovr_i2o(maxovr)       !lon index of overlap input cell
    integer , intent(out)   :: jovr_i2o(maxovr)       !lat index of overlap input cell
    real(r8), intent(out)   :: wovr_i2o(maxovr)       !weight    of overlap input cell
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!EOP
!
! LOCAL VARIABLES:
    real(r8) :: relerr = 0.00001       !relative error for error checks
    integer  :: ii                     !input  grid longitude loop index
    integer  :: ji                     !input  grid latitude  loop index
    integer  :: n                      !overlap index
!------------------------------------------------------------------------

    ! Get indices and weights for mapping from input grid to output grid

    call areamap_point (io      , jo       , nlon_i  , nlat_i , nlon_o  , &
                        nlat_o  , numlon_i , lon_i   , lat_i  , mask_i  , &
                        lon_o   , lat_o    , fland_o , area_o , novr_i2o, &
                        iovr_i2o, jovr_i2o , wovr_i2o, lon_i_offset, maxovr)

  end subroutine areaini_point

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: areamap_point
!
! !INTERFACE:
  subroutine areamap_point(io     , jo       , nlon_i  , nlat_i , nlon_o , &
                           nlat_o , numlon_i , lon_i   , lat_i  , mask_i , &
                           lon_o  , lat_o    , fland_o , area_o , n_ovr  , &
                           i_ovr  , j_ovr    , w_ovr   , lon_i_offset, &
                           maxovr)
!
! !DESCRIPTION:
! weights and indices for area of overlap between grids
! Get indices and weights for area-averaging between input and output grids.
! For each output grid cell find:
!
!    o number of input grid cells that overlap with output grid cell (n_ovr)
!    o longitude index (1 <= i_ovr <= nlon_i) of the overlapping input grid cell
!    o latitude index  (1 <= j_ovr <= nlat_i) of the overlapping input grid cell
!    o fractional overlap of input grid cell (w_ovr)
!
! so that for
!
! field values fld_i on an  input grid with dimensions nlon_i and nlat_i
! field values fld_o on an output grid with dimensions nlon_o and nlat_o are
!
! fld_o(io,jo) =
! fld_i(i_ovr(io,jo,     1),j_ovr(io,jo,     1)) * w_ovr(io,jo,     1) +
!                             ... + ... +
! fld_i(i_ovr(io,jo,maxovr),j_ovr(io,jo,maxovr)) * w_ovr(io,jo,maxovr)
!
! Note: maxovr is some number greater than n_ovr. Weights of zero are
! used for the excess points
!
! !ARGUMENTS:
    implicit none
    integer ,intent(in)   :: io                     !output grid longitude index
    integer ,intent(in)   :: jo                     !output grid latitude index
    integer ,intent(in)   :: nlon_i                 !input grid : max number of long points
    integer ,intent(in)   :: nlat_i                 !input grid : number of latitude points
    integer ,intent(in)   :: nlon_o                 !output grid: max number of long points
    integer ,intent(in)   :: nlat_o                 !output grid: number of latitude points
    integer ,intent(in)   :: numlon_i(nlat_i)       !input grid : number long points for lat
    real(r8),intent(in)   :: lon_i(nlon_i+1,nlat_i) !input grid : cell lons, west edge (deg)
    real(r8),intent(in)   :: lon_i_offset(nlon_i+1,nlat_i) !input grid : cell lons, west edge (deg)
    real(r8),intent(in)   :: lat_i(nlat_i+1)        !input grid : cell lats, south edge (deg)
    real(r8),intent(in)   :: mask_i(nlon_i,nlat_i)  !input grid : mask (0, 1)
    real(r8),intent(in)   :: lon_o(nlon_o+1,nlat_o) !output grid: cell lons, west edge  (deg)
    real(r8),intent(in)   :: lat_o(nlat_o+1)        !output grid: cell lats, south edge (deg)
    real(r8),intent(in)   :: fland_o                !output grid: fraction that is land
    real(r8),intent(in)   :: area_o                 !output grid: cell area
    integer ,intent(out)  :: n_ovr                  !number of overlapping input cells
    integer ,intent(in)   :: maxovr                 !maximum number of overlapping cells
    integer ,intent(out)  :: i_ovr(maxovr)          !lon index, overlapping input cell
    integer ,intent(out)  :: j_ovr(maxovr)          !lat index, overlapping input cell
    real(r8),intent(out)  :: w_ovr(maxovr)          !overlap weights for input cells
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: ii                  !input  grid longitude loop index
    integer  :: ji                  !input  grid latitude  loop index
    integer  :: n                   !overlapping cell index
    real(r8) :: offset              !used to shift x-grid 360 degrees
    real(r8) :: f_ovr               !sum of overlap weights
    real(r8) :: relerr = 0.00001    !max error: sum overlap weights ne 1
    real(r8) :: dx_i                !input grid  longitudinal range
    real(r8) :: dy_i                !input grid  latitudinal  range
    real(r8) :: dx_o                !output grid longitudinal range
    real(r8) :: dy_o                !output grid latitudinal  range
!------------------------------------------------------------------------

    ! --------------------------------------------------------------------
    ! Initialize overlap weights on output grid to zero for maximum
    ! number of overlapping points. Set lat and lon indices of overlapping
    ! input cells to dummy values. Set number of overlapping cells to zero
    ! --------------------------------------------------------------------

    n_ovr           = 0
    i_ovr(1:maxovr) = 1
    j_ovr(1:maxovr) = 1
    w_ovr(1:maxovr) = 0.

    ! --------------------------------------------------------------------
    ! First pass to find cells that overlap and area of overlap
    ! --------------------------------------------------------------------

    call areaovr_point (io    , jo    , nlon_i , nlat_i , numlon_i , &
                        lon_i , lat_i , nlon_o , nlat_o , lon_o    , &
                        lat_o , n_ovr , i_ovr  , j_ovr  , w_ovr    , &
                        maxovr )

    ! --------------------------------------------------------------------
    ! Second pass to find cells that overlap and area of overlap with
    ! shifted grid
    ! --------------------------------------------------------------------

    call areaovr_point (io           , jo    , nlon_i , nlat_i , numlon_i , &
                        lon_i_offset , lat_i , nlon_o , nlat_o , lon_o    , &
                        lat_o        , n_ovr , i_ovr  , j_ovr  , w_ovr    , &
                        maxovr)

    ! --------------------------------------------------------------------
    ! Normalize areas of overlap to get fractional contribution of each
    ! overlapping grid cell (input grid) to grid cell average on output grid.
    ! Normally, do this by dividing area of overlap by area of output grid cell.
    ! But, only have data for land cells on input grid. So if output grid cell
    ! overlaps with land and non-land cells (input grid), do not have valid
    ! non-land data for area-average. Instead, weight by area of land using
    ! [mask_i], which has a value of one for land and zero for ocean. If
    ! [mask_i] = 1, input grid cell contributes to output grid cell average.
    ! If [mask_i] = 0, input grid cell does not contribute to output grid cell
    ! average.
    ! --------------------------------------------------------------------

    ! find total land area of overlapping input cells

    f_ovr = 0.
    do n = 1, n_ovr
       ii = i_ovr(n)
       ji = j_ovr(n)
       f_ovr = f_ovr + w_ovr(n)*mask_i(ii,ji)
    end do

    ! make sure area of overlap is less than or equal to output grid cell area

    if ((f_ovr-area_o)/area_o > relerr) then
       write (6,*) 'AREAMAP error: area not conserved for lon,lat = ',io,jo
       write (6,'(a30,e20.10)') 'sum of overlap area = ',f_ovr
       write (6,'(a30,e20.10)') 'area of output grid = ',area_o
       call abort()
    end if

    ! make weights

    do n = 1, n_ovr
       ii = i_ovr(n)
       ji = j_ovr(n)
       if (f_ovr > 0.) then
          w_ovr(n) = w_ovr(n)*mask_i(ii,ji) / f_ovr
       else
          w_ovr(n) = 0.
       end if
    end do

    ! --------------------------------------------------------------------
    ! Error check: overlap weights for input grid cells must sum to 1. This
    ! is always true if both grids span the same domain. However, if one
    ! grid is a subset of the other grid, this is only true when mapping
    ! from the full grid to the subset. When input grid covers a smaller
    ! domain than the output grid, this test is not valid.
    ! --------------------------------------------------------------------

    dx_i = lon_i(nlon_i+1,1) - lon_i(1,1)
    dx_o = lon_o(nlon_o+1,1) - lon_o(1,1)

    if (lat_i(nlat_i+1) > lat_i(1)) then      !South to North grid
       dy_i = lat_i(nlat_i+1) - lat_i(1)
    else                                      !North to South grid
       dy_i = lat_i(1) - lat_i(nlat_i+1)
    end if
    if (lat_o(nlat_o+1) > lat_o(1)) then      !South to North grid
       dy_o = lat_o(nlat_o+1) - lat_o(1)
    else                                      !North to South grid
       dy_o = lat_o(1) - lat_o(nlat_o+1)
    end if

    if (abs(dx_i-dx_o)>relerr .or. abs(dy_i-dy_o)>relerr) then
       if (dx_i<dx_o .or. dy_i<dy_o) then
          write (6,*) 'AREAMAP warning: area-average not valid for '
          write (6,*) '   input  grid of ',nlon_i,' x ',nlat_i
          write (6,*) '   output grid of ',nlon_o,' x ',nlat_o
          return
       end if
    end if

    ! error check only valid if output grid cell has land. non-land cells
    ! will have weights equal to zero

    f_ovr = 0.
    do n = 1, maxovr
       f_ovr = f_ovr + w_ovr(n)
    end do

    if ( (fland_o > 0.) .and. (abs(f_ovr-1.) > relerr)) then
       write (6,*) 'AREAMAP_POINT error: area not conserved for lon,lat = ',io,jo
       write (6,'(a30,e20.10)') 'sum of overlap weights = ',f_ovr
       call abort()
    end if

  end subroutine areamap_point

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: areaovr_point
!
! !INTERFACE:
  subroutine areaovr_point(io     , jo     , nlon_i , nlat_i , numlon_i , &
                           lon_i  , lat_i  , nlon_o , nlat_o , lon_o    , &
                           lat_o  , n_ovr  , i_ovr  , j_ovr , w_ovr, maxovr)
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
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: io                     !output grid lon index
    integer , intent(in) :: jo                     !output grid lat index
    integer , intent(in) :: nlon_i                 !input grid : max number of longitude points
    integer , intent(in) :: nlat_i                 !input grid : number of latitude points
    integer , intent(in) :: numlon_i(nlat_i)       !input grid : number of lon points for lat
    real(r8), intent(in) :: lon_i(nlon_i+1,nlat_i) !input grid : cell longitude, W edge (deg)
    real(r8), intent(in) :: lat_i(nlat_i+1)        !input grid : cell latitude, S edge (deg)
    integer , intent(in) :: nlon_o                 !output grid: max number of longitude points
    integer , intent(in) :: nlat_o                 !output grid: number of latitude points
    real(r8), intent(in) :: lon_o(nlon_o+1,nlat_o) !output grid: cell longitude, W edge (deg)
    real(r8), intent(in) :: lat_o(nlat_o+1)        !output grid: cell latitude, S edge (deg)
    integer , intent(inout) :: n_ovr               !number of overlapping input cells
    integer , intent(in)    :: maxovr              !maximum number of overlapping cells
    integer , intent(inout) :: i_ovr(maxovr)       !lon index, overlapping input cell
    integer , intent(inout) :: j_ovr(maxovr)       !lat index, overlapping input cell
    real(r8), intent(inout) :: w_ovr(maxovr)       !area of overlap for input cells
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!EOP
!
! LOCAL VARIABLES:
    integer indexo1        !output grid lat. index according to orientn
    integer indexo2        !output grid lat. index according to orientn
    integer indexo3        !output grid lat. index according to orientn
    integer ii             !input  grid longitude loop index
    integer ji             !input  grid latitude  loop index
    integer indexi1        !input grid lat. index according to orientn
    integer indexi2        !input grid lat. index according to orientn
    integer indexi3        !input grid lat. index according to orientn
    real(r8) lonw          !west longitudes of overlap
    real(r8) lone          !east longitudes of overlap
    real(r8) dx            !difference in longitudes
    real(r8) lats          !south latitudes of overlap
    real(r8) latn          !north latitudes of overlap
    real(r8) dy            !difference in latitudes
    real(r8) deg2rad       !pi/180
    real(r8) a_ovr         !area of overlap
!------------------------------------------------------------------------

    deg2rad = (SHR_CONST_PI) / 180.

    ! choose the right index according to the orientation of the data

    if (lat_o(nlat_o+1) > lat_o(1)) then
       indexo1 = jo+1        !south to north along the edges
       indexo2 = jo          !south to north along the edges
       indexo3 = jo          !south to north at the center of cell
    else
       indexo1 = nlat_o+1-jo !north to south along the edges
       indexo2 = nlat_o+2-jo !north to south along the edges
       indexo3 = nlat_o+1-jo !north to south at the center of cell
    end if

    ! loop through all input grid cells to find overlap with output grid

    do ji = 1, nlat_i

       ! choose the right index according to the orientation of the data

       if (lat_i(nlat_i+1) > lat_i(1)) then
          indexi1 = ji          !south to north along the edges
          indexi2 = ji+1        !south to north along the edges
          indexi3 = ji          !south to north at the center of cell
       else
          indexi1 = nlat_i+2-ji !north to south along the edges
          indexi2 = nlat_i+1-ji !north to south along the edges
          indexi3 = nlat_i+1-ji !north to south at the center of cell
       end if

       ! lat okay

       if ( lat_i(indexi1)<lat_o(indexo1) .and. &
            lat_i(indexi2)>lat_o(indexo2) ) then

          do ii = 1, numlon_i(indexi3)

             ! lon okay

             if (lon_i(ii,indexi3)<lon_o(io+1,indexo3) .and. &
                 lon_i(ii+1,indexi3)>lon_o(io,indexo3)) then

                ! increment number of overlapping cells. make sure 0 < n_ovr < maxovr

                n_ovr = n_ovr + 1
                if (n_ovr > maxovr) then
                   write (6,*)' AREAOVR_POINT error: n_ovr= ',n_ovr, &
                        ' exceeded parameter maxovr = ',maxovr, &
                        ' for output lon,lat = ',io,indexo3
                   write(6,*) ' increase parameter maxovr'
                   call abort()
                end if

                ! determine area of overlap

                lone = min(lon_o(io+1,indexo3),lon_i(ii+1,indexi3))*deg2rad !e edge
                lonw = max(lon_o(io  ,indexo3),lon_i(ii  ,indexi3))*deg2rad !w edge
                dx = max(0.0_r8,(lone-lonw))
                latn = min(lat_o(indexo1),lat_i(indexi2))*deg2rad !n edge
                lats = max(lat_o(indexo2),lat_i(indexi1))*deg2rad !s edge
                dy = max(0.0_r8,(sin(latn)-sin(lats)))
                a_ovr = dx*dy*re*re

                ! save lat, lon indices of overlapping cell and area of overlap

                i_ovr(n_ovr) = ii
                j_ovr(n_ovr) = indexi3
                w_ovr(n_ovr) = a_ovr

             end if   !end lon-okay if-block

          end do   !end input-lon loop

       end if   !end lat-okay if-block

    end do   !end input-lat loop

  end subroutine areaovr_point

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkmxovr
!
! !INTERFACE:
  subroutine mkmxovr (nlon_i, nlat_i, numlon_i, lon_i, lat_i, &
                      nlon_o, nlat_o, numlon_o, lon_o, lat_o, &
                      mxovr , n_ovr  )
!
! !DESCRIPTION:
! find maxinum numver of overlapping cells
! For each output grid cell: find overlapping input grid cells that
! that overlap with output grid cell. Cells overlap if:
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
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: nlon_i                 !input grid : max number of longitude points
    integer , intent(in) :: nlat_i                 !input grid : number of latitude points
    integer , intent(in) :: numlon_i(nlat_i)       !input grid : number of lon points for lat
    real(r8), intent(inout) :: lon_i(nlon_i+1,nlat_i) !input grid : cell longitude, W edge (deg)
    real(r8), intent(in) :: lat_i(nlat_i+1)        !input grid : cell latitude, S edge (deg)
    integer , intent(in) :: nlon_o                 !output grid: max number of longitude points
    integer , intent(in) :: nlat_o                 !output grid: number of latitude points
    integer , intent(in) :: numlon_o(nlat_o)       !output grid: number of lon points for lat
    real(r8), intent(in) :: lon_o(nlon_o+1,nlat_o) !output grid: cell longitude, W edge (deg)
    real(r8), intent(in) :: lat_o(nlat_o+1)        !output grid: cell latitude, S edge (deg)
    integer , intent(out):: n_ovr(nlon_o,nlat_o)   !number of overlapping input cells
    integer , intent(out):: mxovr                  !maximum number of overlapping input cells
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
!
    integer, parameter :: mxovr_ceiling = 100000 !very large value should only check for bad error
    integer :: ii          !input  grid longitude loop index
    integer :: ji          !input  grid latitude  loop index
    integer :: io          !output grid longitude loop index
    integer :: jo          !output grid latitude  loop index
    integer :: indexi1     !input  grid lat. index according to orientn
    integer :: indexi2     !input  grid lat. index according to orientn
    integer :: indexi3     !input  grid lat. index according to orientn
    integer :: indexo1     !output grid lat. index according to orientn
    integer :: indexo2     !output grid lat. index according to orientn
    integer :: indexo3     !output grid lat. index according to orientn
    real(r8) :: lonw       !west longitudes of overlap
    real(r8) :: lone       !east longitudes of overlap
    real(r8) :: dx         !difference in longitudes
    real(r8) :: lats       !south latitudes of overlap
    real(r8) :: latn       !north latitudes of overlap
    real(r8) :: dy         !difference in latitudes
    real(r8) :: deg2rad    !pi/180
    real(r8) :: offset     !shifted longitudinal offset
!-----------------------------------------------------------------------

    ! Set number of overlapping cells to zero and initialize mxovr and deg2rad

    mxovr = 0
    deg2rad = (SHR_CONST_PI) / 180.
    n_ovr(:,:) = 0

    ! loop through output grid cells
    ! choose the right index according to the orientation of the data
    do jo = 1, nlat_o
       if (lat_o(nlat_o+1) > lat_o(1)) then
          indexo1 = jo+1        !south to north along the edges
          indexo2 = jo          !south to north along the edges
          indexo3 = jo          !south to north at the center of cell
       else
          indexo1 = nlat_o+1-jo !north to south along the edges
          indexo2 = nlat_o+2-jo !north to south along the edges
          indexo3 = nlat_o+1-jo !north to south at the center of cell
       end if

       ! loop through all input grid cells to find overlap with output grid
       do io = 1, numlon_o(indexo3)
          do ji = 1, nlat_i
             ! choose the right index according to the orientation of the data
             if (lat_i(nlat_i+1) > lat_i(1)) then
                indexi1 = ji          !south to north along the edges
                indexi2 = ji+1        !south to north along the edges
                indexi3 = ji          !south to north at the center of cell
             else
                indexi1 = nlat_i+2-ji !north to south along the edges
                indexi2 = nlat_i+1-ji !north to south along the edges
                indexi3 = nlat_i+1-ji !north to south at the center of cell
             end if

             ! if lat and lon okay then increment number of overlapping cells
             ! make sure 0 < n_ovr < mxovr_ceiling
             if (lat_i(indexi1)<lat_o(indexo1) .and. lat_i(indexi2)>lat_o(indexo2)) then
                do ii = 1, numlon_i(indexi3)
                   if (lon_i(ii,indexi3)<lon_o(io+1,indexo3) .and. &
                       lon_i(ii+1,indexi3)>lon_o(io,indexo3)) then
                      n_ovr(io,indexo3) = n_ovr(io,indexo3) + 1
                      if (n_ovr(io,indexo3) > mxovr_ceiling) then
                         write (6,100) n_ovr(io,indexo3),mxovr_ceiling,io,indexo3
                         call abort()
                      end if
                      if (n_ovr(io,indexo3) > mxovr) then
                         mxovr = n_ovr(io,indexo3)
                      endif
                   end if
                end do
             end if
          end do
       end do
    end do

    ! Shift x-grid to locate periodic grid intersections. This
    ! assumes that all lon_i(1,j) have the same value for all
    ! latitudes j and that the same holds for lon_o(1,j)

    if (lon_i(1,1) < lon_o(1,1)) then
       offset = 360.0
    else
       offset = -360.0
    end if

    do ji = 1, nlat_i
       do ii = 1, numlon_i(ji) + 1
          lon_i(ii,ji) = lon_i(ii,ji) + offset
       end do
    end do

    ! loop through output grid cells
    ! choose the right index according to the orientation of the data
    do jo = 1, nlat_o
       if (lat_o(nlat_o+1) > lat_o(1)) then
          indexo1 = jo+1        !south to north along the edges
          indexo2 = jo          !south to north along the edges
          indexo3 = jo          !south to north at the center of cell
       else
          indexo1 = nlat_o+1-jo !north to south along the edges
          indexo2 = nlat_o+2-jo !north to south along the edges
          indexo3 = nlat_o+1-jo !north to south at the center of cell
       end if

       ! loop through all input grid cells to find overlap with output grid
       do io = 1, numlon_o(indexo3)
          do ji = 1, nlat_i
             ! choose the right index according to the orientation of the data
             if (lat_i(nlat_i+1) > lat_i(1)) then
                indexi1 = ji          !south to north along the edges
                indexi2 = ji+1        !south to north along the edges
                indexi3 = ji          !south to north at the center of cell
             else
                indexi1 = nlat_i+2-ji !north to south along the edges
                indexi2 = nlat_i+1-ji !north to south along the edges
                indexi3 = nlat_i+1-ji !north to south at the center of cell
             end if

             ! if lat and lon okay then increment number of overlapping cells
             ! make sure 0 < n_ovr < mxovr_ceiling
             if (lat_i(indexi1)<lat_o(indexo1) .and. lat_i(indexi2)>lat_o(indexo2)) then
                do ii = 1, numlon_i(indexi3)
                   if (lon_i(ii,indexi3)<lon_o(io+1,indexo3) .and. &
                        lon_i(ii+1,indexi3)>lon_o(io,indexo3)) then
                      n_ovr(io,indexo3) = n_ovr(io,indexo3) + 1
                      if (n_ovr(io,indexo3) > mxovr_ceiling) then
                         write (6,100) n_ovr(io,indexo3),mxovr_ceiling,io,indexo3
                         call abort()
                      end if
                      if (n_ovr(io,indexo3) > mxovr) then
                         mxovr = n_ovr(io,indexo3)
                      endif
                   end if
                end do
             end if
          end do
       end do
    end do

    ! restore x-grid (un-shift x-grid)
    do ji = 1, nlat_i
       do ii = 1, numlon_i(ji) + 1
          lon_i(ii,ji) = lon_i(ii,ji) - offset
       end do
    end do

100 format(' ','MKMXOVR error: n_ovr= ',i4,' exceeded mx_ceiling = ', &
         i4,' for output lon,lat = ',i4,i4)

  end subroutine mkmxovr

end module areaMod



















