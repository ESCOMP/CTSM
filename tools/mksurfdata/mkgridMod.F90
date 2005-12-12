module mkgridMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkgridMod
!
! !DESCRIPTION:
! Routines to read in land model grid
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use fileutils, only : getfil
  use mkvarsur
  use mkvarctl
  use areaMod
  use ncdio
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: read_grid      ! Read land model grid.
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
! !IROUTINE: read_grid
!
! !INTERFACE:
  subroutine read_grid(lsmlon, lsmlat)
!
! !DESCRIPTION:
! Read land model grid.
! If namelist variable mksrf_fgrid is the empty string, then
! the corresponding file will be used to determine the land model grid
! Assume that file contains:
!   dimensions: lon,lat
!   variables : lon,lat or LONGXY,LATIXY (optionally numlon)  
!   variables : LANDFRAC (optionally LANDMASK)
!
! !ARGUMENTS:
    implicit none
    integer, intent(out) :: lsmlon, lsmlat
!
! !CALLED FROM:
! subroutine mkgrid_offline in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,k,n                 ! indices
    integer  :: ncid                    ! netCDF file id
    integer  :: dimid                   ! netCDF dimension id
    integer  :: varid                   ! netCDF variable id
    integer  :: ret                     ! netCDF return code
    real(r8), allocatable :: lon(:)     ! input longitude array (full grid)
    real(r8), allocatable :: lat(:)     ! input latitude array (full grid)
    character(len=256) :: locfn         ! local file name
    character(len= 32) :: subname = 'read_grid'
!-----------------------------------------------------------------------

    write (6,*) 'Attempting to read land grid data .....'
    write (6,'(72a1)') ("-",i=1,60)

    if (mksrf_fgrid_global /= ' ' .and. mksrf_fgrid_regional /= ' ') then 	
       write(6,*)' mksrf_fgrid_global  = ',mksrf_fgrid_global	 
       write(6,*)' mksrf_fgrid_regioanl= ',mksrf_fgrid_regional
       write(6,*)' must not specify both as non-blank namelist input if read in grid'
       stop
    end if

    ! Get file containing grid info

    if (mksrf_fgrid_global /= ' ') then

       call getfil (mksrf_fgrid_global, locfn, 0)
       call check_ret(nf_open(locfn, 0, ncid), subname)

    else if (mksrf_fgrid_regional /= ' ') then

       call getfil (mksrf_fgrid_regional, locfn, 0)
       call check_ret(nf_open(locfn, 0, ncid), subname)

    end if

    ! Set model grid edges (these will not be used)

    lsmedge(1) = spval
    lsmedge(2) = spval
    lsmedge(3) = spval
    lsmedge(4) = spval

    ! Read in longitude and latitude dimension info

    call check_ret(nf_inq_dimid  (ncid, 'lon', dimid), subname)
    call check_ret(nf_inq_dimlen (ncid, dimid, lsmlon), subname)

    call check_ret(nf_inq_dimid  (ncid, 'lat', dimid), subname)
    call check_ret(nf_inq_dimlen (ncid, dimid, lsmlat), subname)

    ! Allocate dynamic memory

    allocate(lon(lsmlon)            , &
             lat(lsmlat)            , &
             numlon(lsmlat)         , &                 
             latixy(lsmlon,lsmlat)  , &          
             longxy(lsmlon,lsmlat)  , &          
             landmask(lsmlon,lsmlat), &        
             landfrac(lsmlon,lsmlat), &
             area(lsmlon,lsmlat)    , &            
             lats(lsmlat+1)         , &               
             lonw(lsmlon+1,lsmlat), stat=ret)                 
    if (ret/=0) call abort

    numlon(:)     = 0
    latixy(:,:)   = spval
    longxy(:,:)   = spval
    landfrac(:,:) = spval
    landmask(:,:) = -999

    ! Determine grid longitudes.

    ret = nf_inq_varid (ncid, 'lon', varid)
    if (ret == NF_NOERR) then
       call check_ret(nf_get_var_double (ncid, varid, lon), subname)
       do j = 1,lsmlat
          do i = 1,lsmlon
             longxy(i,j) = lon(i)
          end do
       end do
    else 
       ret = nf_inq_varid (ncid, 'LONGXY', varid)
       if (ret == NF_NOERR) then
          call check_ret(nf_get_var_double (ncid, varid, longxy), subname)
       else
          write(6,*)subname,' error: either lon or LONGXY must be specified'
       end if
    end if

    ret = nf_inq_varid (ncid, 'numlon', varid)
    if (ret == NF_NOERR) then
       call check_ret(nf_get_var_double (ncid, varid, numlon), subname)
    else
       numlon(:) = lsmlon
    end if

    ! Determine grid latitudes.

    ret = nf_inq_varid (ncid, 'lat', varid)
    if (ret == NF_NOERR) then
       call check_ret(nf_get_var_double (ncid, varid, lat), subname)
       do j = 1,lsmlat
          do i =1,lsmlon
             latixy(i,j) = lat(j)
          end do
       end do
    else 
       ret = nf_inq_varid (ncid, 'LATIXY', varid)
       if (ret == NF_NOERR) then
          call check_ret(nf_get_var_double (ncid, varid, latixy), subname)
       else
          write(6,*)subname,' error: either lat or LATIXY must be specified'
       end if
    end if

    ! Get land fraction (this is required)

    call check_ret(nf_inq_varid (ncid, 'LANDFRAC' , varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, landfrac), subname)

    ! Determine land mask (this is optional - if not on dataset, set it using landfrac)

    ret = nf_inq_varid (ncid, 'LANDMASK', varid)
    if (ret == NF_NOERR) then
       call check_ret(nf_inq_varid (ncid, 'LANDMASK', varid), subname)
       call check_ret(nf_get_var_int(ncid, varid, landmask), subname)
    else
       do j = 1,lsmlat
          do i = 1,numlon(j)
             if (landfrac(i,j) > 0._r8) then
                landmask(i,j) = 1
             else
                landmask(i,j) = 0
             end if
          end do
       end do
    end if
    
    ! Define land grid edges and grid cell areas

    if (mksrf_fgrid_global /= ' ') then 	

       call celledge (lsmlat, lsmlon, numlon, longxy, latixy, lats, lonw, area)

    else if (mksrf_fgrid_regional /= ' ') then 	

       call celledge (lsmlat, lsmlon, longxy, latixy, lats, lonw, area)

    else

       write(6,*)' must specify either mksrf_fgrid_global or mksrf_fgrid_regional ',&
            ' to non-blank field if read in grid'
       stop

    end if

    write (6,'(72a1)') ("-",i=1,60)
    write (6,*) 'Successfully read land grid data'
    write (6,*)

    deallocate(lon,lat)

  end subroutine read_grid

end module mkgridMod
