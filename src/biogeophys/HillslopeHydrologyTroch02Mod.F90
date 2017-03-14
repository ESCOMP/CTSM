module HillslopeHydrologyTroch02Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate geomorphological quantities for hillslope columns 
  ! assuming Troch et al., AWR, 2002 profile and plan shapes.
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

  implicit none
  private   
  save

  ! PRIVATE 
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  ! !PUBLIC TYPES:
  
  type, extends(hillslope_geomorphology_type), public :: &
       hillslope_geomorphology_troch02_type
     private
  
     ! variable declarations
     
   contains
     
     procedure :: Init
     procedure :: hcol_width
     procedure :: hcol_elevation
     procedure :: hcol_slope
     
  end type hillslope_geomorphology_troch02_type
  
  !-----------------------------------------------------------------------
  interface hillslope_geomorphology_troch02_type

  end interface hillslope_geomorphology_troch02_type


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
    class(hillslope_geomorphology_troch02_type) , intent(in) :: this
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
    character(len=*), parameter :: subname = 'hcol_width'
    !-----------------------------------------------------------------------

    ! width function has special case for n = 2
    ! in this implementation, integration limits depend on sign of beta
    if (abs(alpha - 2._r8) < eps) then

       ! function blows up for x0=0; integration limits set by trial and error
       if(beta < 0._r8) then
          y0 = hill_width/2._r8
          yl = 0.1_r8
          x0=hill_length *(yl/y0)**(-hill_height/(beta*hill_length**2))
       else
          x0 = 0.2_r8
          y0 = (hill_width/2._r8)&
               *(x0/hill_length)**(beta*hill_length**2/hill_height)
       endif

       ! compiler does not like log(zero)
       if (x == 0._r8) then
          if (beta <  0._r8) then
             width = hill_width/2._r8
          else
             width = eps
          endif
       else
          width = y0*(x/x0)**(beta*hill_length**2/hill_height)
       endif
    else 
       ! alpha /= 2 case, x0 equals zero
       y0 = hill_width/2._r8
       if(beta > 0._r8) then
          y0 = y0 * exp(-(2._r8*beta*hill_length**2) &
               / (hill_height*(2._r8 - alpha)*alpha))
       endif
       ! compiler does not like zero to a negative power.
       if (x == 0._r8) then
          width = y0
       else
          width = y0*exp((((2._r8*beta*hill_length**2) &
               / (hill_height*(2._r8 - alpha)*alpha)) &
               * (x/hill_length)**(2._r8-alpha)))
       endif
    endif
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
    class(hillslope_geomorphology_troch02_type) , intent(in) :: this
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
    ! elevation is first integrated analytically in across-slope direction
    ! then summed in along-slope direction

    dx = (xtop - xbottom)/real(ndiv)

    elev = 0._r8
    area = 0._r8
    do n = 0, ndiv-1
       x = xbottom + (n+0.5)*dx
       y = this%hcol_width(x,alpha,beta,hill_length,hill_width,hill_height)
       dA = dx * y
       area = area + dA
       elev = elev + dx * (hill_height*y*(x/hill_length)**alpha &
            + beta*y**3/12._r8)
    enddo
    elev = elev / area

  end function hcol_elevation

  !-----------------------------------------------------------------------

  function hcol_slope(this, xtop,xbottom,alpha, hill_length, hill_height) result(slope)
    !
    ! !DESCRIPTION:
    ! Returns mean along-hillslope slope of hillslope column 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(hillslope_geomorphology_troch02_type) , intent(in) :: this
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
    use clm_varctl      , only : nhillslope
    use clm_varcon      , only : zmin_bedrock, zisoi
    use clm_varpar      , only : nlevsoi
    use spmdMod         , only : masterproc
    use fileutils       , only : getfil
    use clm_varcon      , only : spval, ispval, grlnd 
    use landunit_varcon , only : istsoil
    use ncdio_pio

    !
    ! !ARGUMENTS:
    class(hillslope_geomorphology_troch02_type) , intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    character(len=*) , intent(in) :: fsurdat      ! surface data file name
    real(r8), pointer     :: ihillslope_in(:,:)   ! read in - integer
    real(r8), pointer     :: fhillslope_in(:,:)   ! read in - float
    integer,  allocatable :: pct_hillslope(:,:)   ! percent of landunit occupied by hillslope
    real(r8), allocatable :: hill_alpha(:,:)      ! hillslope 'alpha' parameter
    real(r8), allocatable :: hill_beta(:,:)       ! hillslope 'beta' parameter
    real(r8), allocatable :: hill_length(:,:)     ! hillslope length [m]
    real(r8), allocatable :: hill_width(:,:)      ! hillslope width  [m]
    real(r8), allocatable :: hill_height(:,:)     ! hillslope height [m]

    type(file_desc_t)     :: ncid                 ! netcdf id
    logical               :: readvar              ! check whether variable on file    
    character(len=256)    :: locfn                ! local filename
    integer               :: ierr                 ! error code
    integer               :: c, l, g, j, nh       ! indices

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
         hill_alpha(bounds%begl:bounds%endl,nhillslope),    &
         hill_beta(bounds%begl:bounds%endl,nhillslope),   &
         hill_length(bounds%begl:bounds%endl,nhillslope),   &
         hill_width(bounds%begl:bounds%endl,nhillslope), &
         hill_height(bounds%begl:bounds%endl,nhillslope), &
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
       
       allocate(fhillslope_in(bounds%begg:bounds%endg,nhillslope))
       call ncd_io(ncid=ncid, varname='h_alpha', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          if (masterproc) then
             call endrun( 'ERROR:: h_alpha not found on surface data set.'//errmsg(sourcefile, __LINE__) )
          end if
       end if
       
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          hill_alpha(l,:) = fhillslope_in(g,:)
       enddo
       
       call ncd_io(ncid=ncid, varname='h_beta', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          if (masterproc) then
             call endrun( 'ERROR:: h_beta not found on surface data set.'//errmsg(sourcefile, __LINE__) )
          end if
       end if
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          hill_beta(l,:) = fhillslope_in(g,:)
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
             
             do c = lun%coli(l), lun%colf(l)
                nh = col%hillslope_ndx(c)
                ctop    = lun%coli(l)+(nh-1)*lun%ncolumns(l)/nhillslope
                cbottom = lun%coli(l)+(nh)*lun%ncolumns(l)/nhillslope - 1

                ! distance of lower edge of column from hillslope bottom
                col%hill_distance(c) = this%hcol_distance(c, &
                     ctop, cbottom, hill_length(l,nh))
                !              if (masterproc) write(iulog,*) 'hd: ',c,col%hill_distance(c)
                
                !distance of lower edge of column from hillslope bottom
                column_length = hill_length(l,nh)/(lun%ncolumns(l)/nhillslope)
                le_distance = col%hill_distance(c) - 0.5_r8*column_length
                ue_distance = col%hill_distance(c) + 0.5_r8*column_length

                ! width of lower edge of column from hillslope bottom
                col%hill_width(c) = this%hcol_width(le_distance, &
                     hill_alpha(l,nh), hill_beta(l,nh),hill_length(l,nh), &
                     hill_width(l,nh),hill_height(l,nh))
                !              if (masterproc) write(iulog,*) 'hw: ',c,col%hill_width(c)
                
                col%hill_area(c) = this%hcol_area(ue_distance, le_distance, &
                     hill_alpha(l,nh), hill_beta(l,nh),hill_length(l,nh), &
                     hill_width(l,nh),hill_height(l,nh))
                !              if (masterproc) write(iulog,*) 'ha: ',c,col%hill_area(c)
                
                ! mean elevation of column relative to mean gridcell elevation
                col%hill_elev(c) = this%hcol_elevation(ue_distance, le_distance, &
                     hill_alpha(l,nh), hill_beta(l,nh),hill_length(l,nh), &
                     hill_width(l,nh),hill_height(l,nh))
                !              if (masterproc) write(iulog,*) 'he: ',c,col%hill_elev(c)
                
                ! mean along-hill slope of column
                col%hill_slope(c) = this%hcol_slope(ue_distance, le_distance, &
                     hill_alpha(l,nh),hill_length(l,nh), hill_height(l,nh))
                !              if (masterproc) write(iulog,*) 'hs: ',c,col%hill_slope(c)
                
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
             
          endif ! end of istsoil
       enddo    ! end of loop over landunits
       
       deallocate(pct_hillslope,hill_alpha,hill_beta,hill_length, &
         hill_width,hill_height)

       call ncd_pio_closefile(ncid)

  end subroutine Init

end module HillslopeHydrologyTroch02Mod
