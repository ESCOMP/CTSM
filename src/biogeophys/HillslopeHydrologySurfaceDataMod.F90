module HillslopeHydrologySurfaceDataMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate geomorphological quantities for hillslope columns 
  ! assuming independent profile and plan shapes.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use HillslopeHydrologyBaseMod, only : hillslope_geomorphology_type
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use spmdMod        , only : masterproc
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog
  use decompMod      , only : bounds_type
!scs
  use clm_varcon     , only : rpi
!scs
  implicit none
  private   
  save

  ! PRIVATE 
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  ! !PUBLIC TYPES:
  
  type, extends(hillslope_geomorphology_type), public :: &
       hillslope_geomorphology_surfacedata_type
     private
  
     ! variable declarations
     
   contains
     
     procedure :: Init
     procedure :: hcol_width
     procedure :: hcol_elevation
     procedure :: hcol_slope
     
  end type hillslope_geomorphology_surfacedata_type
  
  !-----------------------------------------------------------------------
  interface hillslope_geomorphology_surfacedata_type

  end interface hillslope_geomorphology_surfacedata_type


  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  function hcol_width(this,x,alpha,beta,hill_length,hill_width,hill_height) result(width)
    !
    ! !DESCRIPTION:
    ! Returns width of hillslope column.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(hillslope_geomorphology_surfacedata_type) , intent(in) :: this
    real(r8) :: width                   ! function result
    real(r8), intent(in) :: x           ! distance along hillslope
    real(r8), intent(in) :: alpha       ! profile curvature parameter
    real(r8), intent(in) :: beta        ! plan curvature parameter
    real(r8), intent(in) :: hill_length ! total hillslope length
    real(r8), intent(in) :: hill_width  ! total hillslope width
    real(r8), intent(in) :: hill_height ! total hillslope height
    !
    ! !LOCAL VARIABLES:
    real(r8) :: eps = 1.e-6_r8
    real(r8) :: x0,y0,yl                ! integration limits
    real(r8) :: slope, intercept
    character(len=*), parameter :: subname = 'hcol_width'
    !-----------------------------------------------------------------------

    ! width varies linearly with distance

    ! divergent
!    y0 = hill_width/2._r8
!    yl = 0.25_r8 * yl
    ! convergent
    yl = hill_width/2._r8
    y0 = 0.1_r8 * yl
!    y0 = 1.0_r8 * yl

    slope = (yl - y0)/(hill_length)
    intercept = y0
    width = slope*x + intercept

    ! hillslope width is twice integral [0:x]
    width = 2._r8 * width 

  end function hcol_width
             
  !-----------------------------------------------------------------------
  function hcol_elevation(this,xtop,xbottom,alpha,beta,hill_length,hill_width,hill_height) result(elev)
    !
    ! !DESCRIPTION:
    ! Returns mean elevation of column (relative to hillslope bottom).
    ! Area-weighted mean elevation is calculated by 
    ! numerically integrating using hcol_width function.    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(hillslope_geomorphology_surfacedata_type) , intent(in) :: this
    real(r8) :: elev                    ! function result
    real(r8), intent(in) :: xtop        ! upper integration limit
    real(r8), intent(in) :: xbottom     ! lower integration limit
    real(r8), intent(in) :: alpha       ! profile curvature parameter
    real(r8), intent(in) :: beta        ! plan curvature parameter
    real(r8), intent(in) :: hill_length ! total hillslope length
    real(r8), intent(in) :: hill_width  ! total hillslope width
    real(r8), intent(in) :: hill_height ! total hillslope height
    !
    ! !LOCAL VARIABLES:

    integer  :: n
    integer, parameter  :: ndiv = 100
    real(r8) :: x, y
    real(r8) :: dx
    real(r8) :: dA
    real(r8) :: area
    character(len=*), parameter :: subname = 'hcol_elevation'
    !-----------------------------------------------------------------------
    ! mean elevation of column relative to hillslope bottom
    ! elevation is considered 1-d

    dx = (xtop - xbottom)/real(ndiv)

    elev = 0._r8
    area = 0._r8
    do n = 0, ndiv-1
       x = xbottom + (n+0.5)*dx
       y = this%hcol_width(x,alpha,beta,hill_length,hill_width,hill_height)
       dA = dx * y
       area = area + dA
       elev = elev + dx * y * hill_height*(x/hill_length)**alpha
    enddo
    elev = elev / area

  end function hcol_elevation

  !-----------------------------------------------------------------------
  function hcol_slope(this,xtop,xbottom,alpha, hill_length, hill_height) result(slope)
    !
    ! !DESCRIPTION:
    ! Returns mean along-hillslope slope of hillslope column 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(hillslope_geomorphology_surfacedata_type) , intent(in) :: this
    real(r8) :: slope                   ! function result
    real(r8), intent(in) :: xtop        ! distance to upper edge of column
    real(r8), intent(in) :: xbottom     ! distance to lower edge of column
    real(r8), intent(in) :: alpha       ! hillslope profile curvature parameter
    real(r8), intent(in) :: hill_length ! total hillslope length
    real(r8), intent(in) :: hill_height ! total hillslope height
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'hcol_slope'
    !-----------------------------------------------------------------------

    ! mean along-hill slope of column
    slope = hill_height &
         * ((xtop/hill_length)**alpha &
         - (xbottom/hill_length)**alpha) &
         / (xtop - xbottom)

  end function hcol_slope

  !-----------------------------------------------------------------------

  subroutine Init(this,bounds,fsurdat)
    !
    ! !DESCRIPTION:
    ! Initialize hillslope geomorphology
    !
    ! !USES:
    use LandunitType    , only : lun                
    use GridcellType    , only : grc                
    use ColumnType      , only : col                
    use clm_varctl      , only : nhillslope, nmaxhillcol
    use clm_varcon      , only : zmin_bedrock, zisoi
    use clm_varpar      , only : nlevsoi
    use spmdMod         , only : masterproc
    use fileutils       , only : getfil
    use clm_varcon      , only : spval, ispval, grlnd 
    use landunit_varcon , only : istsoil
    use ncdio_pio

    !
    ! !ARGUMENTS:
    class(hillslope_geomorphology_surfacedata_type) , intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    character(len=*) , intent(in) :: fsurdat    ! surface data file name
    real(r8), pointer     :: ihillslope_in(:,:) ! read in - integer
    real(r8), pointer     :: fhillslope_in(:,:) ! read in - float
    integer,  allocatable :: pct_hillslope(:,:) ! percent of landunit occupied by hillslope
    integer,  allocatable :: hill_ndx(:,:)      ! hillslope index
    integer,  allocatable :: col_ndx(:,:)       ! column index
    integer,  allocatable :: col_dndx(:,:)      ! downhill column index
    real(r8), allocatable :: hill_slope(:,:)    ! hillslope slope  [m/m]
    real(r8), allocatable :: hill_area(:,:)     ! hillslope area   [m2]
    real(r8), allocatable :: hill_length(:,:)   ! hillslope length [m]
    real(r8), allocatable :: hill_width(:,:)    ! hillslope width  [m]
    real(r8), allocatable :: hill_height(:,:)   ! hillslope height [m]

    type(file_desc_t)     :: ncid                 ! netcdf id
    logical               :: readvar              ! check whether variable on file    
    character(len=256)    :: locfn                ! local filename
    integer               :: ierr                 ! error code
    integer               :: c, l, g, i, j, ci, nh       ! indices

    real(r8)              :: hillslope_area       ! total area of hillslope
    real(r8)              :: column_length        ! length of column [m]
    real(r8)              :: le_distance          ! distance of lower edge of column from bottom of hillslope
    real(r8)              :: ue_distance          ! distance of upper edge of column from bottom of hillslope
    integer               :: ctop, cbottom        ! hillslope top and bottom column indices

    character(len=*), parameter :: subname = 'Init'

    !-----------------------------------------------------------------------

    ! Open surface dataset to read in data below 

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    allocate(pct_hillslope(bounds%begl:bounds%endl,nhillslope),    &
         hill_ndx(bounds%begl:bounds%endl,nmaxhillcol),    &
         col_ndx(bounds%begl:bounds%endl,nmaxhillcol),    &
         col_dndx(bounds%begl:bounds%endl,nmaxhillcol),    &
         hill_slope(bounds%begl:bounds%endl,nmaxhillcol),    &
         hill_area(bounds%begl:bounds%endl,nmaxhillcol),   &
         hill_length(bounds%begl:bounds%endl,nmaxhillcol),   &
         hill_width(bounds%begl:bounds%endl,nmaxhillcol), &
         hill_height(bounds%begl:bounds%endl,nmaxhillcol), &
         stat=ierr)
       
       allocate(ihillslope_in(bounds%begg:bounds%endg,nhillslope))
       
       call ncd_io(ncid=ncid, varname='pct_hillslope', flag='read', data=ihillslope_in, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          if (masterproc) then
             call endrun( 'ERROR:: pct_hillslope not found on surface data set.'//errmsg(sourcefile, __LINE__) )
          end if
       end if
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          pct_hillslope(l,:) = ihillslope_in(g,:)
       enddo
       deallocate(ihillslope_in)
       
       allocate(ihillslope_in(bounds%begg:bounds%endg,nmaxhillcol))
       
       call ncd_io(ncid=ncid, varname='hillslope_index', flag='read', data=ihillslope_in, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          if (masterproc) then
             call endrun( 'ERROR:: hillslope_index not found on surface data set.'//errmsg(sourcefile, __LINE__) )
          end if
       end if
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          hill_ndx(l,:) = ihillslope_in(g,:)
       enddo

       call ncd_io(ncid=ncid, varname='column_index', flag='read', data=ihillslope_in, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          if (masterproc) then
             call endrun( 'ERROR:: column_index not found on surface data set.'//errmsg(sourcefile, __LINE__) )
          end if
       end if
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          col_ndx(l,:) = ihillslope_in(g,:)
       enddo

       call ncd_io(ncid=ncid, varname='downhill_column_index', flag='read', data=ihillslope_in, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          if (masterproc) then
             call endrun( 'ERROR:: downhill_column_index not found on surface data set.'//errmsg(sourcefile, __LINE__) )
          end if
       end if
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          col_dndx(l,:) = ihillslope_in(g,:)
       enddo
       deallocate(ihillslope_in)
       
       allocate(fhillslope_in(bounds%begg:bounds%endg,nmaxhillcol))
       call ncd_io(ncid=ncid, varname='h_slope', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          if (masterproc) then
             call endrun( 'ERROR:: h_slope not found on surface data set.'//errmsg(sourcefile, __LINE__) )
          end if
       end if
       
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          hill_slope(l,:) = fhillslope_in(g,:)
       enddo
       
       call ncd_io(ncid=ncid, varname='h_area', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          if (masterproc) then
             call endrun( 'ERROR:: h_area not found on surface data set.'//errmsg(sourcefile, __LINE__) )
          end if
       end if
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          hill_area(l,:) = fhillslope_in(g,:)
       enddo
       call ncd_io(ncid=ncid, varname='h_length', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          if (masterproc) then
             call endrun( 'ERROR:: h_length not found on surface data set.'//errmsg(sourcefile, __LINE__) )
          end if
       end if
       
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          hill_length(l,:) = fhillslope_in(g,:)
       enddo
       
       call ncd_io(ncid=ncid, varname='h_width', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          if (masterproc) then
             call endrun( 'ERROR:: h_width not found on surface data set.'//errmsg(sourcefile, __LINE__) )
          end if
       end if
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          hill_width(l,:) = fhillslope_in(g,:)
       enddo
       
       call ncd_io(ncid=ncid, varname='h_height', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          if (masterproc) then
             call endrun( 'ERROR:: h_height not found on surface data set.'//errmsg(sourcefile, __LINE__) )
          end if
       end if
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          hill_height(l,:) = fhillslope_in(g,:)
       enddo
       deallocate(fhillslope_in)
       
       !  Set hillslope hydrology column level variables
       !  This needs to match how columns set up in subgridMod
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          if(lun%itype(l) == istsoil) then
             ! lun%coli is the uppermost column in the hillslope, lun%colf is the lowermost
             
             ! map external column index to internal column index
             do c = lun%coli(l), lun%colf(l)
                ! ci should span [1:nhillcolumns(l)]
                ci = c-lun%coli(l)+1
                ! relative separation should be the same
                if (col_dndx(l,ci) <= -999) then
                   ! lowermost column of hillslope has no downstream neighbor
!scs                   col%cold(c) = col_dndx(l,ci)
                   col%cold(c) = ispval
                else
                   col%cold(c) = c + (col_dndx(l,ci) - col_ndx(l,ci))
                endif
             enddo


             do c = lun%coli(l), lun%colf(l)
                ci = c-lun%coli(l)+1

                col%hillslope_ndx(c) = hill_ndx(l,ci)

                ! Find uphill/downhill neighbors (must calculate 
                ! col%cold in a previous loop
!scs                col%colu(c) = -999 
                col%colu(c) = ispval
                do i = lun%coli(l), lun%colf(l)
                   if(i == col%cold(c)) then
                      col%colu(i) = c
                   endif
                enddo

! distance of lower edge of column from hillslope bottom
                col%hill_distance(c) = hill_length(l,ci)
! width of lower edge of column 
                col%hill_width(c) = hill_width(l,ci)
! mean elevation of column relative to gridcell mean elevation
                col%hill_elev(c) = hill_height(l,ci)
! mean along-hill slope of column
                col%hill_slope(c) = hill_slope(l,ci)
! area of column
                col%hill_area(c) = hill_area(l,ci)

             enddo

             ! Now that column areas are determined, column weights can be recalculated
             hillslope_area = 0._r8
             ! area weighted by pct_hillslope
             do c = lun%coli(l), lun%colf(l)
                nh = col%hillslope_ndx(c)
                hillslope_area = hillslope_area &
                     + col%hill_area(c)*(pct_hillslope(l,nh)*0.01_r8)
             enddo
             do c = lun%coli(l), lun%colf(l)
                nh = col%hillslope_ndx(c)
                col%wtlunit(c) = col%hill_area(c) &
                     * (pct_hillslope(l,nh)*0.01_r8)/hillslope_area      
             enddo
             
             !  Set column bedrock index
             !thin soil for non-riparian columns, thick for riparian
if(1==2) then
             do c =  lun%coli(l), lun%colf(l)
                if(col%cold(c) /= ispval) then 
                   do j = 1,nlevsoi
                      if(zisoi(j-1) > zmin_bedrock) then
                         if (zisoi(j-1) < 0.5_r8 .and. zisoi(j) >= 0.5_r8) then
                            col%nbedrock(c) = j
                         end if
                      end if
                   enddo
                else 
                   do j = 1,nlevsoi 
                      if(zisoi(j-1) > zmin_bedrock) then
                         if (zisoi(j-1) < 3.0_r8 .and. zisoi(j) >= 3.0_r8) then
                            col%nbedrock(c) = j
                         end if
                      end if
                   enddo
                endif
             end do
endif             
          endif ! end of istsoil
       enddo    ! end of loop over landunits
       
       deallocate(pct_hillslope,hill_ndx,col_ndx,col_dndx, &
            hill_slope,hill_area,hill_length, &
            hill_width,hill_height)

       call ncd_pio_closefile(ncid)

  end subroutine Init

end module HillslopeHydrologySurfaceDataMod
