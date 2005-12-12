module creategridMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkgridMod
!
! !DESCRIPTION:
! Routines to create land model grid
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use fileutils   , only : getfil
  use mkvarsur
  use areaMod
  use ncdio
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: creategrid    ! Generate land model grid.

! !PRIVATE MEMBER FUNCTIONS:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: creategrid
!
! !INTERFACE:
  subroutine creategrid(lsmlon, lsmlat, lsmedge, fnavyoro)
!
! !DESCRIPTION:
! Generate land model grid.
! Surface grid edges -- Grids do not have to be global. To allow this, grids
! must define the north, east, south, and west edges:
! namelist variables
!    o fnavyoro : 20 min navy orography dataset
!    o edgen (edge(1)) : northern edge of grid (degrees): >  -90 and <= 90
!    o edgee (edge(2)) : eastern edge of grid (degrees) : see following notes
!    o edges (edge(3)) : southern edge of grid (degrees): >= -90 and <  90
!    o edgew (edge(4)) : western edge of grid (degrees) : see following notes
! For partial grids, northern and southern edges are any latitude
! between 90 (North Pole) and -90 (South Pole). Western and eastern
! edges are any longitude between -180 and 180, with longitudes
! west of Greenwich negative. That is, western edge >= -180 and < 180;
! eastern edge > western edge and <= 180.
! For global grids, northern and southern edges are 90 (North Pole)
! and -90 (South Pole). The western edge of the longitude grid starts
! at the dateline if the grid is generated (the longitudes for each grid
! cell correspond with the edges (i.e., range from -180 to 180)).
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lsmlon
    integer , intent(in) :: lsmlat
    real(r8), intent(in) :: lsmedge(4)    
    character(len=*), intent(in) :: fnavyoro
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: locfn                !local file name

    integer  :: i,j,k,n                        !indices
    integer  :: ii,ji,io,jo                    !indices
    integer  :: ncid                           !netCDF file id
    integer  :: dimid                          !netCDF dimension id
    integer  :: varid                          !netCDF variable id
    integer  :: ier                            !error status

    integer  :: nlon_i                         !input number of longitudes
    integer  :: nlat_i                         !input number of latitudes
    real(r8) :: dx                             !land model cell width
    real(r8) :: dy                             !land model cell length
    real(r8) :: edge_i(4)                      !input grid: N,E,S,W edges (degrees)
    real(r8), allocatable :: latixy_i(:,:)     !input grid: latitude (degrees)
    real(r8), allocatable :: longxy_i(:,:)     !input grid: longitude (degrees)
    integer , allocatable :: numlon_i(:)       !input grid: number longitude points by lat
    real(r8), allocatable :: lon_i(:,:)        !input grid: longitude, west edge (degrees)
    real(r8), allocatable :: lon_i_offset(:,:) !input grid: longitude, west edge (degrees)
    real(r8), allocatable :: lat_i(:)          !input grid: latitude, south edge (degrees)
    real(r8), allocatable :: area_i(:,:)       !input grid: cell area
    real(r8), allocatable :: mask_i(:,:)       !input grid: mask (0, 1)
    real(r8), allocatable :: fland_i(:,:)      !input grid: fractional land

    integer, parameter :: maxovr = 100000
    real(r8) :: mask_o                         !output grid: mask (0, 1)
    integer  :: novr_i2o                       !number of overlapping input cells
    integer  :: iovr_i2o(maxovr)               !lon index of overlap input cell
    integer  :: jovr_i2o(maxovr)               !lat index of overlap input cell
    real(r8) :: wovr_i2o(maxovr)               !weight    of overlap input cell
    real(r8) :: offset                         !used to shift x-grid 360 degrees

    real(r8) :: fld_o(lsmlon,lsmlat)           !output grid: dummy field
    real(r8) :: fld_i                          !input grid: dummy field
    real(r8) :: sum_fldo                       !global sum of dummy output field
    real(r8) :: sum_fldi                       !global sum of dummy input field
    real(r8) :: relerr = 0.00001               !max error: sum overlap weights ne 1

    real(r8) :: flandmin = 0.001               !minimum land fraction for grid cell to be called land
    real(r8) :: edgen  
    real(r8) :: edgee  
    real(r8) :: edges  
    real(r8) :: edgew  
    character(len= 32) :: subname = 'create_grid'
!-----------------------------------------------------------------------

    ! Allocate memory for output grid (variables are in module mkvarsur)

    allocate(numlon(lsmlat)         , &                 
             latixy(lsmlon,lsmlat)  , &          
             longxy(lsmlon,lsmlat)  , &          
             landmask(lsmlon,lsmlat), &        
             landfrac(lsmlon,lsmlat), &
             area(lsmlon,lsmlat)    , &            
             lats(lsmlat+1)         , &               
             lonw(lsmlon+1,lsmlat), stat=ier)                 
    if (ier /= 0) call abort

    ! Determine output grid longitudes and latitudes in increments of dx and dy
    ! Global latitude grid goes from south pole to north pole
    ! Global longitude grid starts at Dateline with western edge on Dateline

    edgen  = lsmedge(1)
    edgee  = lsmedge(2)
    edges  = lsmedge(3)
    edgew  = lsmedge(4)
       
    numlon(:)   = lsmlon

    dx = (edgee - edgew) / lsmlon
    dy = (edgen - edges) / lsmlat
    do j = 1, lsmlat
       do i = 1, lsmlon
          longxy(i,j) = edgew + (2*i-1)*dx / 2.
          latixy(i,j) = edges + (2*j-1)*dy / 2
       end do
    end do

    ! Define edges and area of output land model grid cells

    call celledge (lsmlat     , lsmlon     , numlon     , longxy     ,&
                   latixy     , lsmedge(1) , lsmedge(2) , lsmedge(3) ,&
                   lsmedge(4) , lats       , lonw       , area)

    ! Read input navy orography data and obtain fractional land

    call getfil (fnavyoro, locfn, 0)
    call check_ret(nf_open(locfn, 0, ncid), subname)

    call check_ret(nf_inq_dimid  (ncid, 'lon', dimid), subname)
    call check_ret(nf_inq_dimlen (ncid, dimid, nlon_i), subname)

    call check_ret(nf_inq_dimid  (ncid, 'lat', dimid), subname)
    call check_ret(nf_inq_dimlen (ncid, dimid, nlat_i), subname)

    allocate (latixy_i(nlon_i,nlat_i), &
              longxy_i(nlon_i,nlat_i), &
              numlon_i(nlat_i), &
              lon_i(nlon_i+1,nlat_i), &
              lon_i_offset(nlon_i+1,nlat_i), &
              lat_i(nlat_i+1), &
              area_i(nlon_i,nlat_i), &
              mask_i(nlon_i,nlat_i), &
              fland_i(nlon_i,nlat_i), stat=ier)
    if (ier/=0) call abort

    call check_ret(nf_inq_varid (ncid, 'LATIXY', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, latixy_i), subname)

    call check_ret(nf_inq_varid (ncid, 'LONGXY', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, longxy_i), subname)

    call check_ret(nf_inq_varid (ncid, 'NUMLON', varid), subname)
    call check_ret(nf_get_var_int (ncid, varid, numlon_i), subname)

    call check_ret(nf_inq_varid (ncid, 'EDGEN', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, edge_i(1)), subname)

    call check_ret(nf_inq_varid (ncid, 'EDGEE', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, edge_i(2)), subname)

    call check_ret(nf_inq_varid (ncid, 'EDGES', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, edge_i(3)), subname)

    call check_ret(nf_inq_varid (ncid, 'EDGEW', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, edge_i(4)), subname)

    call check_ret(nf_inq_varid (ncid, 'LANDFRAC', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, fland_i), subname)

    call check_ret(nf_close(ncid), subname)

    ! determine maximim overlap for height resolution and model grids and
    ! and allocate dynamic memory for overlap arrays
    ! first determine input grid cell and cell areas

    numlon_i(:) = nlon_i

    call celledge (nlat_i    , nlon_i    , numlon_i  , longxy_i  , &
                   latixy_i  , edge_i(1) , edge_i(2) , edge_i(3) , &
                   edge_i(4) , lat_i     , lon_i     , area_i)

    do ji = 1, nlat_i
       do ii = 1, numlon_i(ji)
          mask_i(ii,ji) = 1.
       end do
    end do

    ! Shift x-grid to locate periodic grid intersections. This
    ! assumes that all lon_i(1,j) have the same value for all
    ! latitudes j and that the same holds for lon_o(1,j)

    if (lon_i(1,1) < lonw(1,1)) then
       offset = 360.0
    else
       offset = -360.0
    end if

    do ji = 1, nlat_i
       do ii = 1, numlon_i(ji) + 1
          lon_i_offset(ii,ji) = lon_i(ii,ji) + offset
       end do
    end do

  ! Process each cell on land model grid
  ! novr_i2o - number of input grid cells that overlap each land grid cell
  ! iovr_i2o - longitude index of overlapping input grid cell
  ! jovr_i2o - latitude  index of overlapping input grid cell
  ! wovr_i2o - fraction of land grid cell overlapped by input grid cell

!$OMP PARALLEL DO PRIVATE (io,jo,ii,ji,n,mask_o,novr_i2o,iovr_i2o,jovr_i2o,wovr_i2o,fld_i)
!CSD$ PARALLEL DO PRIVATE (io,jo,ii,ji,n,mask_o,novr_i2o,iovr_i2o,jovr_i2o,wovr_i2o,fld_i)
    do jo = 1, lsmlat
       do io = 1, numlon(jo)

        ! Determine areas of overlap and indices

          mask_o = 1.

          call areaini_point (io        , jo          , nlon_i  , nlat_i  , numlon_i , &
                              lon_i      , lon_i_offset, lat_i   , area_i  , mask_i   , &
                              lsmlon     , lsmlat      , numlon  , lonw    , lats     , &
                              area(io,jo), mask_o      , novr_i2o, iovr_i2o, jovr_i2o , &
                              wovr_i2o   , maxovr)

        ! Make area average

          landfrac(io,jo) = 0.
          do n = 1, novr_i2o   !overlap cell index
             ii = iovr_i2o(n)  !lon index (input grid) of overlap cell
             ji = jovr_i2o(n)  !lat index (input grid) of overlap cell
             landfrac(io,jo) = landfrac(io,jo) + fland_i(ii,ji) * wovr_i2o(n)
          end do
          if (landfrac(io,jo) > 100.000001_r8) then
             write (6,*) 'MKGRID error: fland = ',landfrac(io,jo), &
                  ' is greater than 100 for lon,lat = ',io,jo
             call abort
          end if

          ! Global sum of output field -- must multiply by fraction of
          ! output grid that is land as determined by input grid

          fld_o(io,jo) = 0.
          do n = 1, novr_i2o
             ii = iovr_i2o(n)
             ji = jovr_i2o(n)
             fld_i = ((ji-1)*nlon_i + ii)
             fld_o(io,jo) = fld_o(io,jo) + wovr_i2o(n) * fld_i
          end do

       end do  !end of output longitude loop
    end do     !end of output latitude  loop
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO

    ! -----------------------------------------------------------------
    ! Error check1
    ! Compare global sum fld_o to global sum fld_i.
    ! -----------------------------------------------------------------

    ! This check is true only if both grids span the same domain.
    ! To obtain global sum of input field must multiply by
    ! fraction of input grid that is land as determined by input grid

    sum_fldo = 0.
    do jo = 1,lsmlat
       do io = 1,numlon(jo)
          sum_fldo = sum_fldo + area(io,jo) * fld_o(io,jo)
       end do
    end do

    sum_fldi = 0.
    do ji = 1, nlat_i
       do ii = 1, numlon_i(ji)
          fld_i = ((ji-1)*nlon_i + ii)
          sum_fldi = sum_fldi + area_i(ii,ji) * fld_i
       end do
    end do

    if ( abs(edgen - edges) == 180. .and. &
         abs(edgee - edgew) == 360. ) then
       if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
          write (6,*) 'MKGRID error: input field not conserved'
          write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
          write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
          call abort
       end if
    end if

    ! Determine land mask

    where (landfrac(:,:) < flandmin)
       landmask(:,:) = 0     !ocean
    elsewhere
       landmask(:,:) = 1     !land
    endwhere

    ! Reset landfrac to zero where landmask has been set to zero

    where (landmask(:,:) == 0)
       landfrac(:,:) = 0
    endwhere

    ! deallocate dynamic memory

    deallocate (numlon_i)
    deallocate (latixy_i)
    deallocate (longxy_i)
    deallocate (lon_i)
    deallocate (lon_i_offset)
    deallocate (lat_i)
    deallocate (area_i)
    deallocate (mask_i)
    deallocate (fland_i)

    write (6,*) 'Successfully makde land grid data'
    write (6,*)

  end subroutine creategrid

end module creategridMod
