!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: mksoicol
!
! !INTERFACE:
subroutine mksoicol (lsmlon, lsmlat, fsoicol, ndiag, pctgla_o, soil_color_o, nsoicol)
!
! !DESCRIPTION:
! Make soil color classes for model grid from BATS T42 data
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use fileutils   , only : getfil
  use mkvarpar
  use mkvarsur
  use mkvarctl
  use areaMod 
  use ncdio
!
! !ARGUMENTS:
  implicit none
  integer , intent(in) :: lsmlon, lsmlat               ! clm grid resolution
  character(len=*), intent(in) :: fsoicol              ! input soicol dataset file name
  integer , intent(in) :: ndiag                        ! unit number for diagnostic output
  real(r8), intent(in) :: pctgla_o(lsmlon,lsmlat)      ! output grid: percent glacier
  integer , intent(out):: soil_color_o(lsmlon,lsmlat)  ! output grid: soil color classes
  integer , intent(out):: nsoicol                      ! number of soil colors 
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Gordon Bonan
!
!EOP
!
! !LOCAL VARIABLES:
  integer, parameter :: maxovr = 100000       ! 
  character(len=256) locfn                    ! local dataset file name
  integer :: nlon_i                           ! input grid : longitude points (read in)
  integer :: nlat_i                           ! input grid : latitude  points (read in)
  integer :: ncid,dimid,varid                 ! input netCDF id's
  integer :: ier                              ! error status
  integer :: ii                               ! longitude index for BATS grid
  integer :: io                               ! longitude index for model grid
  integer :: ji                               ! latitude  index for BATS grid
  integer :: jo                               ! latitude  index for model grid
  integer :: k                                ! temporary BATS or model soil color
  integer :: miss = 99999                     ! missing data indicator
  integer :: n                                ! loop index
  integer, parameter :: num=2                 ! get 1st and 2nd largest areas of overlap
  integer :: wsti(num)                        ! index to 1st and 2nd largest values in wst
  real(r8) ::edge_i(4)                        ! input grid: N,E,S,W edges (degrees)
  real(r8), allocatable :: latixy_i(:,:)      ! input grid: latitude (degrees)
  real(r8), allocatable :: longxy_i(:,:)      ! input grid: longitude (degrees)
  integer , allocatable :: numlon_i(:)        ! input grid: number longitude points by lat
  real(r8), allocatable :: lon_i(:,:)         ! input grid: longitude, west edge (degrees)
  real(r8), allocatable :: lon_i_offset(:,:)  ! input grid: longitude, west edge (degrees)
  real(r8), allocatable :: lat_i(:)           ! input grid: latitude, south edge (degrees)
  real(r8), allocatable :: area_i(:,:)        ! input grid: cell area
  real(r8), allocatable :: mask_i(:,:)        ! input grid: mask (0, 1)
  integer , allocatable :: soil_color_i(:,:)  ! input grid: BATS soil color
  real(r8), allocatable :: wst(:)             ! overlap weights, by surface type
  real(r8), allocatable :: gast_i(:)          ! input grid : global area, by surface type
  real(r8), allocatable :: gast_o(:)          ! output grid : global area, by surface type
  character(len=35), allocatable :: col(:)    ! name of each color
  real(r8) :: mask_o                          ! output grid: mask (0, 1)
  integer  :: novr_i2o                        ! number of overlapping input cells
  integer  :: iovr_i2o(maxovr)                ! lon index of overlap input cell
  integer  :: jovr_i2o(maxovr)                ! lat index of overlap input cell
  real(r8) :: wovr_i2o(maxovr)                ! weight    of overlap input cell
  real(r8) :: offset                          ! used to shift x-grid 360 degrees
  real(r8) :: fld_o(lsmlon,lsmlat)            ! output grid: dummy field
  real(r8) :: fld_i                           ! input grid: dummy field
  real(r8) :: sum_fldo                        ! global sum of dummy output field
  real(r8) :: sum_fldi                        ! global sum of dummy input field
  real(r8) :: relerr = 0.00001                ! max error: sum overlap weights ne 1
  integer  :: color                           ! 0: none; 1: some
  character(len=32) :: subname = 'mksoicol'   ! subroutine name
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make soil color classes .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input soil colors
  ! -----------------------------------------------------------------

  ! Obtain input grid info

  call getfil (fsoicol, locfn, 0)
  call check_ret(nf_open(locfn, 0, ncid), subname)

  call check_ret(nf_inq_dimid  (ncid, 'lon', dimid), subname)
  call check_ret(nf_inq_dimlen (ncid, dimid, nlon_i), subname)

  call check_ret(nf_inq_dimid  (ncid, 'lat', dimid), subname)
  call check_ret(nf_inq_dimlen (ncid, dimid, nlat_i), subname)

  allocate (latixy_i(nlon_i,nlat_i), longxy_i(nlon_i,nlat_i), numlon_i(nlat_i), &
            lon_i(nlon_i+1,nlat_i), lon_i_offset(nlon_i+1,nlat_i), lat_i(nlat_i+1), &
            area_i(nlon_i,nlat_i), mask_i(nlon_i,nlat_i), soil_color_i(nlon_i,nlat_i), &
            stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'LATIXY', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, latixy_i), subname)

  call check_ret(nf_inq_varid (ncid, 'LONGXY', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, longxy_i), subname)

  call check_ret(nf_inq_varid (ncid, 'EDGEN', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, edge_i(1)), subname)

  call check_ret(nf_inq_varid (ncid, 'EDGEE', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, edge_i(2)), subname)

  call check_ret(nf_inq_varid (ncid, 'EDGES', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, edge_i(3)), subname)

  call check_ret(nf_inq_varid (ncid, 'EDGEW', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, edge_i(4)), subname)

  ! Obtain input data

  call check_ret(nf_inq_varid (ncid, 'SOIL_COLOR', varid), subname)
  call check_ret(nf_get_var_int (ncid, varid, soil_color_i), subname)

  call check_ret(nf_close(ncid), subname)

  nsoicol = maxval(soil_color_i)
  write(6,*)'nsoicol = ',nsoicol

  allocate(wst(0:nsoicol), gast_i(0:nsoicol), gast_o(0:nsoicol), col(0:nsoicol))

  ! -----------------------------------------------------------------
  ! Define the model color classes: 0 to nsoicol
  ! -----------------------------------------------------------------

  if (nsoicol == 20) then
     col(0)  = 'no soil                            '
     col(1)  = 'class 1: light                     '
     col(2)  = 'class 2:                           '
     col(3)  = 'class 3:                           '
     col(4)  = 'class 4:                           '
     col(5)  = 'class 5:                           '
     col(6)  = 'class 6:                           '
     col(7)  = 'class 7:                           '
     col(8)  = 'class 8:                           '
     col(9)  = 'class 9:                           '
     col(10) = 'class 10:                          '
     col(11) = 'class 11:                          '
     col(12) = 'class 12:                          '
     col(13) = 'class 13:                          '
     col(14) = 'class 14:                          '
     col(15) = 'class 15:                          '
     col(16) = 'class 16:                          '
     col(17) = 'class 17:                          '
     col(18) = 'class 18:                          '
     col(19) = 'class 19:                          '
     col(20) = 'class 20: dark                     '
  else if (nsoicol == 8) then
     col(0) = 'no soil                            '
     col(1) = 'class 1: light                     '
     col(2) = 'class 2:                           '
     col(3) = 'class 3:                           '
     col(4) = 'class 4:                           '
     col(5) = 'class 5:                           '
     col(6) = 'class 6:                           '
     col(7) = 'class 7:                           '
     col(8) = 'class 8: dark                      '
  else
     write(6,*)'nsoicol value of ',nsoicol,' is not currently supported'
     call abort()
  end if

  ! -----------------------------------------------------------------
  ! Map data from input grid to land model grid. Get:
  ! -----------------------------------------------------------------

  ! Determine input grid cell and cell areas

  numlon_i(:) = nlon_i

  call celledge (nlat_i    , nlon_i    , numlon_i  , longxy_i  ,  &
                 latixy_i  , edge_i(1) , edge_i(2) , edge_i(3) ,  &
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

!$OMP PARALLEL DO PRIVATE (io,jo,ii,ji,n,k,mask_o,novr_i2o,iovr_i2o,jovr_i2o,wovr_i2o,fld_i, &
!$OMP & wst, wsti)
#if !defined (USE_OMP)
!CSD$ PARALLEL DO PRIVATE (io,jo,ii,ji,n,k,mask_o,novr_i2o,iovr_i2o,jovr_i2o,wovr_i2o,fld_i, &
!CSD$ & wst, wsti)
#endif
  do jo = 1, lsmlat
     do io = 1, numlon(jo)

        ! Determine areas of overlap and indices

        mask_o = 1.

        call areaini_point (io        , jo          , nlon_i  , nlat_i  , numlon_i , &
                           lon_i      , lon_i_offset, lat_i   , area_i  , mask_i   , &
                           lsmlon     , lsmlat      , numlon  , lonw    , lats     , &
                           area(io,jo), mask_o      , novr_i2o, iovr_i2o, jovr_i2o , &
                           wovr_i2o   , maxovr)

        ! Sum overlap weights by color class - make sure dominant non-zero soil
        ! color is used over land (e.g. if have 90% glacier and 10% soil use soil
        ! color not glacier color)

        color = 0
        do k = 0, nsoicol
           wst(k) = 0.
        end do
        do n = 1, novr_i2o         ! overlap cell index
           ii = iovr_i2o(n)        ! lon index (input grid) of overlap cell
           ji = jovr_i2o(n)        ! lat index (input grid) of overlap cell
           k = soil_color_i(ii,ji) ! color class (input grid)
           wst(k) = wst(k) + wovr_i2o(n)
           if (k>0 .and. wst(k)>0.) color = 1
        end do
        if (color == 1) wst(0) = 0.0

        ! Rank non-zero weights by color type. wsti(1) is the most extensive
        ! color type. wsti(2) is the second most extensive color type

        call mkrank (nsoicol, wst, miss, wsti, num)
        soil_color_o(io,jo) = wsti(1)

        ! If land but no color, set color to 15 (in older dataset generic soil color 4)

	if (nsoicol == 8) then
           if (landmask(io,jo)==1 .and. soil_color_o(io,jo)==0) soil_color_o(io,jo) = 4
        else if (nsoicol == 20) then
           if (landmask(io,jo)==1 .and. soil_color_o(io,jo)==0) soil_color_o(io,jo) = 15
        end if

        ! Set ocean colors to zero

        if (landmask(io,jo) == 0) soil_color_o(io,jo) = 0

        ! Set color for grid cells that are 100% glacier to zero. Otherwise,
        ! must have a soil color for the non-glacier portion of grid cell.

        if (abs(pctgla_o(io,jo)-100.)<1.e-06) soil_color_o(io,jo)=0

        ! Error checks

        if (soil_color_o(io,jo) < 0 .or. soil_color_o(io,jo) > nsoicol) then
           write (6,*) 'MKSOICOL error: land model soil color = ', &
                soil_color_o(io,jo),' is not valid for lon,lat = ',io,jo
           call abort()
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
#if !defined (USE_OMP)
!CSD$ END PARALLEL DO
#endif
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

  if ( mksrf_fgrid_global /= ' ') then
     if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
        write (6,*) 'MKSOILCOL error: input field not conserved'
        write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
        write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
        call abort()
     end if
  end if

  ! -----------------------------------------------------------------
  ! Error check2
  ! Compare global area of each soil color on input and output grids
  ! -----------------------------------------------------------------

  ! input grid

  gast_i(:) = 0.
  do ji = 1, nlat_i
     do ii = 1, nlon_i
        k = soil_color_i(ii,ji)
        gast_i(k) = gast_i(k) + area_i(ii,ji)
     end do
  end do

  ! output grid

  gast_o(:) = 0.
  do jo = 1, lsmlat
     do io = 1, numlon(jo)
        k = soil_color_o(io,jo)
        gast_o(k) = gast_o(k) + area(io,jo)
     end do
  end do

  ! area comparison

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',k=1,70)
  write (ndiag,*) 'Soil Color Output'
  write (ndiag,'(1x,70a1)') ('=',k=1,70)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,1001)
1001 format (1x,'soil color type',20x,' input grid area output grid area',/ &
             1x,33x,'     10**6 km**2','      10**6 km**2')
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,*)

  do k = 0, nsoicol
     write (ndiag,1002) col(k),gast_i(k)*1.e-6,gast_o(k)*1.e-6
1002 format (1x,a35,f16.3,f17.3)
  end do

  if (lsmlat > 1) then
     k = lsmlat/2
     write (ndiag,*)
     write (ndiag,*) 'For reference the area on the output grid of a cell near the equator is: '
     write (ndiag,'(f10.3,a14)')area(1,k)*1.e-06,' x 10**6 km**2'
     write (ndiag,*)
  endif

  write (6,*) 'Successfully made soil color classes'
  write (6,*)
  call shr_sys_flush(6)

  ! Deallocate dynamic memory

  deallocate (latixy_i,longxy_i, numlon_i, lon_i, lon_i_offset, &
              lat_i, area_i, mask_i, soil_color_i)

end subroutine mksoicol
