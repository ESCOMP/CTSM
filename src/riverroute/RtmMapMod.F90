module RtmMapMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: RtmMapMod
!
! !DESCRIPTION:
! Area averaging routines
! Thes routines are used for area-average mapping of a field from one
! grid to another.
!
! !USES:
  use domainMod    , only : latlon_type
  use clm_varcon   , only : re
  use clm_varctl   , only : iulog
  use shr_const_mod, only : SHR_CONST_PI
  use shr_kind_mod , only : r8 => shr_kind_r8
  use spmdMod      , only : masterproc
  use abortutils   , only : endrun
  use clm_mct_mod
!
! !PUBLIC TYPES:
  implicit none
  private

!
! !PUBLIC MEMBER FUNCTIONS:
  public :: map_setmapsAR
  public :: celledge
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
! !IROUTINE: map_setmapsAR
!
! !INTERFACE:
  subroutine map_setmapsAR (latlon_i, latlon_o, sMat, fracin, fracout)
!
! !DESCRIPTION:
! area averaging initialization
! This subroutine is used for area-average mapping of a field from one
! grid to another.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(latlon_type), intent(in)    :: latlon_i   ! input domain
    type(latlon_type), intent(in)    :: latlon_o   ! output domain
    type(mct_sMat)   , intent(inout) :: sMat       ! mct sparse matrix plux
    real(r8)         , intent(in),target :: fracin(:)
    real(r8)         , intent(in),target :: fracout(:)
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
    real(r8),pointer :: dx_i(:)      !input grid: dx length
    real(r8),pointer :: dy_i(:)      !input grid: dy length
    real(r8),pointer :: fland_i(:)   !input grid: cell frac
    real(r8),pointer :: lone_i(:)    !input grid: longitude, E edge (degrees)
    real(r8),pointer :: lonw_i(:)    !input grid: longitude, W edge (degrees)
    real(r8),pointer :: latn_i(:)    !input grid: latitude, N edge (degrees)
    real(r8),pointer :: lats_i(:)    !input grid: latitude, S edge (degrees)
    integer          :: nlon_o       !output grid: max number of longitude pts
    integer          :: nlat_o       !output grid: number of latitude  points
    real(r8),pointer :: dx_o(:)      !output grid: dx length
    real(r8),pointer :: dy_o(:)      !output grid: dy length
    real(r8),pointer :: fland_o(:)   !output grid: cell frac
    real(r8),pointer :: lone_o(:)    !output grid: longitude, E edge  (degrees)
    real(r8),pointer :: lonw_o(:)    !output grid: longitude, W edge  (degrees)
    real(r8),pointer :: latn_o(:)    !output grid: latitude, N edge (degrees)
    real(r8),pointer :: lats_o(:)    !output grid: latitude, S edge (degrees
    integer          :: nwts         !local num of weights
    integer          :: n,ns,nw,k    !loop counters
    integer nb,na          !grid sizes
    integer i,j            !loop indices
    integer io,jo,no       !output grid loop index
    integer ii,ji,ni       !input  grid loop index
    integer nji, njilap    !counter and number of j overlaps
    integer,pointer  :: jilap(:) ! j overlap indices
    integer,pointer  :: iilaps(:),iilape(:)  ! i start/end overlap indices
    real(r8) lonw          !west longitudes of overlap
    real(r8) lone          !east longitudes of overlap
    real(r8) dx            !difference in longitudes
    real(r8) lats          !south latitudes of overlap
    real(r8) latn          !north latitudes of overlap
    real(r8) dy            !difference in latitudes
    real(r8) deg2rad       !pi/180
    real(r8) a_ovr         !area of overlap
    integer  noffsetl      !local, number of offsets to test, 0=none, default=1
    real(r8) offset        !value of offset used to shift x-grid 360 degrees
    real(r8) :: sum        !local sum of wts
    integer  :: igrow,igcol,iwgt   ! mct smat field indices
    integer  :: nr         !mct field index
    real(r8) :: wt         !mct weight

    real(r8),pointer :: fld_i(:), fld_o(:)    !dummy fields for testing
    real(r8) :: sum_fldi               !global sum of dummy input field
    real(r8) :: sum_fldo               !global sum of dummy output field
    real(r8) :: sum_areo               !global sum of dummy output field
    real(r8) :: sum_arei               !global sum of dummy input field
    real(r8) :: relerr = 0.001_r8      !relative error for error checks
    real(r8) :: sdx_i                   !input grid  longitudinal range
    real(r8) :: sdy_i                   !input grid  latitudinal  range
    real(r8) :: sdx_o                   !output grid longitudinal range
    real(r8) :: sdy_o                   !output grid latitudinal  range
    integer  :: ier                    !error status
!------------------------------------------------------------------------

    nlon_i =  latlon_i%ni
    nlat_i =  latlon_i%nj
    latn_i => latlon_i%latn
    lats_i => latlon_i%lats
    lonw_i => latlon_i%lonw
    lone_i => latlon_i%lone
    fland_i => fracin

    nlon_o =  latlon_o%ni
    nlat_o =  latlon_o%nj
    latn_o => latlon_o%latn
    lats_o => latlon_o%lats
    lonw_o => latlon_o%lonw
    lone_o => latlon_o%lone
    fland_o => fracout

    ! Dynamically allocate memory

    na = nlon_i*nlat_i
    nb = nlon_o*nlat_o
    noffsetl = 1

    ! Get indices and weights for mapping from input grid to output grid

    allocate(jilap(nlat_i))
    allocate(iilaps(nlon_o),iilape(nlon_o))

    ! Find start/end indices for i overlap to save time later, set default never
    iilaps = nlon_i
    iilape = 1
    do io = 1, nlon_o
       !------ offset -------
       do n = 0,noffsetl   ! loop through offsets
          if (lonw_i(1) < lonw_o(1)) then
             offset = (n*360)
          else
             offset = -(n*360)
          end if

          do ii = 1, nlon_i
             !--- lons overlap ---
             if (lonw_i(ii)+offset < lone_o(io) .and. &
                 lone_i(ii)+offset > lonw_o(io)) then
                 if (ii < iilaps(io)) iilaps(io) = ii
                 if (ii > iilape(io)) iilape(io) = ii
             endif
          enddo
       enddo
    enddo

    do k = 1,2    ! k = 1 compute num of wts, k = 2 set wts
       if (k == 2) then
          call mct_sMat_init(sMat, nb, na, nwts)
          igrow = mct_sMat_indexIA(sMat,'grow')
          igcol = mct_sMat_indexIA(sMat,'gcol')
          iwgt  = mct_sMat_indexRA(sMat,'weight')
       endif
       nw = 0
       deg2rad = SHR_CONST_PI / 180._r8

       !------ output grid -------
       do jo = 1, nlat_o
          njilap = 0
          do ji = 1, nlat_i
             !--- lats overlap ---
             if ( lats_i(ji)<latn_o(jo) .and. &
                  latn_i(ji)>lats_o(jo)) then
                njilap = njilap + 1
                jilap(njilap) = ji
             endif
          enddo

       do io = 1, nlon_o
          ns = nw + 1

          !------ input grid -------
          sum = 0._r8
          do nji = 1, njilap
             ji = jilap(nji)

             !------ offset -------
             do n = 0,noffsetl   ! loop through offsets
                if (lonw_i(1) < lonw_o(1)) then
                   offset = (n*360)
                else
                   offset = -(n*360)
                end if

!             do ii = 1, nlon_i
              do ii = iilaps(io),iilape(io)
                ni = (ji-1)*nlon_i + ii
                no = (jo-1)*nlon_o + io
                !--- lons overlap ---
                if (k == 1) then
                   if (lonw_i(ii)+offset < lone_o(io) .and. &
                       lone_i(ii)+offset > lonw_o(io)) then
                       
                      !------- found overlap ------
                      if (fland_i(ni) > 0._r8 .and. fland_o(no) > 0._r8 ) then
                         nw = nw + 1
                      endif
                   endif
                elseif (k == 2) then
                   if (lonw_i(ii)+offset < lone_o(io) .and. &
                       lone_i(ii)+offset > lonw_o(io)) then

                      ! determine area of overlap
                      lone = min(lone_o(io),lone_i(ii)+offset)*deg2rad 
                      lonw = max(lonw_o(io),lonw_i(ii)+offset)*deg2rad 
                      dx = max(0.0_r8,(lone-lonw))
                      latn = min(latn_o(jo),latn_i(ji))*deg2rad 
                      lats = max(lats_o(jo),lats_i(ji))*deg2rad 
                      dy = max(0.0_r8,(sin(latn)-sin(lats)))
                      a_ovr = dx*dy*re*re

                      sum = sum + a_ovr

                      !------- found overlap ------
                      if (fland_i(ni) > 0._r8 .and. fland_o(no) > 0._r8 ) then
                         nw = nw + 1
                         ! make sure nw <= nwts
                         if (nw > nwts) then
                            write(iulog,*) 'AREAOVR error: nw= ', &
                               nw,' exceeded nwts ceiling = ', &
                               nwts,' for output lon,lat = ',io,jo
                            call endrun
                         end if

                         ! save cell indices, area
                         sMat%data%iAttr(igrow,nw) = no
                         sMat%data%iAttr(igcol,nw) = ni
                         sMat%data%rAttr(iwgt ,nw) = a_ovr
                      endif
                   endif
                end if   ! found overlap lon
             end do   ! ii
             enddo    ! offset loop
          end do   ! ji

          !--- normalize ---
          if (k == 2) then
             do n = ns,nw
                if (sum > 0._r8) then
                   sMat%data%rAttr(iwgt,n) = sMat%data%rAttr(iwgt,n) / sum
                else 
                   sMat%data%rAttr(iwgt,n) = 0._r8
                endif
             enddo
           endif

       end do  ! io
       end do  ! jo

       nwts = nw

    enddo   ! k loop

    deallocate(jilap)
    deallocate(iilaps,iilape)


    ! Error check: global sum fld_o = global sum fld_i.
    ! This true only if both grids span the same domain.

    sdx_i = lone_i(nlon_i) - lonw_i(1)
    sdx_o = lone_o(nlon_o) - lonw_o(1)

    sdy_i = max(abs(latn_i(nlat_i)-lats_i(1)),abs(lats_i(nlat_i)-latn_i(1)))
    sdy_o = max(abs(latn_o(nlat_o)-lats_o(1)),abs(lats_o(nlat_o)-latn_o(1)))

    if (abs(sdx_i-sdx_o)>relerr .or. abs(sdy_i-sdy_o)>relerr) then
       if (masterproc) then
          write(iulog,*) 'MAP_SETMAPSAR warning: lats/lons not overlapping'
          write(iulog,*) '   input  grid of ',nlon_i,' x ',nlat_i,' dx,dy= ',sdx_i,sdy_i
          write(iulog,*) '   output grid of ',nlon_o,' x ',nlat_o,' dx,dy= ',sdx_o,sdy_o
       endif
       return
    end if


    ! Compute dx,dy for checking areas later

    allocate(dx_i(nlon_i),dy_i(nlat_i),dx_o(nlon_o),dy_o(nlat_o))

    do i = 1,nlon_i
       dx_i(i) = (lone_i(i) - lonw_i(i))*deg2rad
    enddo
    do j = 1,nlat_i
       dy_i(j) = sin(latn_i(j)*deg2rad) - sin(lats_i(j)*deg2rad)
    enddo

    do i = 1,nlon_o
       dx_o(i) = (lone_o(i) - lonw_o(i))*deg2rad
    enddo
    do j = 1,nlat_o
       dy_o(j) = sin(latn_o(j)*deg2rad) - sin(lats_o(j)*deg2rad)
    enddo

    ! check for conservation

    allocate(fld_i(na),fld_o(nb))

    sum_fldi = 0._r8
    sum_arei = 0._r8
    do j = 1,nlat_i
    do i = 1,nlon_i
       ni = (j-1)*nlon_i + i
       fld_i(ni) = ni
       sum_fldi = sum_fldi + fld_i(ni)*dx_i(i)*dy_i(j)*fland_i(ni)
       sum_arei = sum_arei + dx_i(i)*dy_i(j)*fland_i(ni)
    end do
    end do

    fld_o = 0._r8
    do n = 1,mct_sMat_lsize(sMat)
       no = sMat%data%iAttr(igrow,n)
       ni = sMat%data%iAttr(igcol,n)
       wt = sMat%data%rAttr(iwgt ,n)
       fld_o(no) = fld_o(no) + wt*fld_i(ni)
    enddo

    sum_fldo = 0._r8
    sum_areo = 0._r8
    do j = 1,nlat_o
    do i = 1,nlon_o
       no = (j-1)*nlon_o + i
       sum_fldo = sum_fldo + dx_o(i)*dy_o(j)*fld_o(no)
       sum_areo = sum_areo + dx_o(i)*dy_o(j)*fland_o(no)
    end do
    end do

    if ( abs(sum_fldo/sum_fldi-1._r8) > relerr ) then
       if ( abs(sum_arei-sum_areo)/(sum_arei+sum_areo) > relerr ) then
          if (masterproc) then
             write(iulog,*) 'MAP_SETMAPSAR: no conservation check done, mask incompatable'
             write(iulog,'(a30,e20.10)') 'global sum input  field = ',sum_arei
             write(iulog,'(a30,e20.10)') 'global sum output field = ',sum_areo
          endif
       else
          if (masterproc) then
             write(iulog,*) 'MAP_SETMAPSAR error: conservation check fail'
             write(iulog,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
             write(iulog,'(a30,e20.10)') 'global sum output field = ',sum_fldo
           endif
       endif
    else
       if (masterproc) then
          write(iulog,*) 'MAP_SETMAPSAR: conservation check passed'
          write(iulog,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
          write(iulog,'(a30,e20.10)') 'global sum output field = ',sum_fldo
       endif
    endif

    deallocate(dx_i,dy_i,dx_o,dy_o)
    deallocate(fld_i,fld_o)

  end subroutine map_setmapsAR

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: celledge
!
! !INTERFACE:
  subroutine celledge(latlon, edgen, edgee, edges, edgew)
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
    type(latlon_type),intent(inout) :: latlon
    real(r8), intent(in) :: edgen             !northern edge of grid (degrees)
    real(r8), intent(in) :: edgee             !eastern edge of grid (degrees)
    real(r8), intent(in) :: edges             !southern edge of grid (degrees)
    real(r8), intent(in) :: edgew             !western edge of grid (degrees)
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 2005.11.20 Updated to latlon datatype by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: i,j     
    real(r8) :: dx      
    integer  :: nlon              
    integer  :: nlat              
    real(r8),pointer :: latc(:) 
    real(r8),pointer :: lats(:)   
    real(r8),pointer :: latn(:)   
    real(r8),pointer :: lonw(:)   
    real(r8),pointer :: lone(:)   
!------------------------------------------------------------------------

    nlon = latlon%ni
    nlat = latlon%nj
    latc => latlon%latc
    lats => latlon%lats
    latn => latlon%latn
    lonw => latlon%lonw
    lone => latlon%lone

    ! Latitudes
    ! Assumes lats are constant on an i line

    if (nlat == 1) then                      ! single latitude
       lats(:)    = edges
       latn(:)    = edgen
    elseif (latc(2) > latc(1)) then  ! South to North grid
       lats(:) = edges
       latn(:) = edgen
       do j = 2, nlat
          lats(j) = (latc(j-1) + latc(j)) / 2._r8
          latn(j-1) = lats(j)
       end do
    else                                     ! North to South grid
       lats(:) = edges
       latn(:) = edgen
       do j = 2, nlat
          latn(j)    = (latc(j-1) + latc(j)) / 2._r8
          lats(j-1) = latn(j)
       end do
    end if

    ! Longitudes
    ! Western edge of first grid cell -- since grid starts with western
    ! edge on Dateline, lonw(1,j)=-180. This is the same as [edgew].

    lonw(:) = edgew
    lone(:) = edgee
    dx = (edgee - edgew) / nlon
    do i = 2, nlon
       lonw(i)    = lonw(i) + (i-1)*dx
       lone(i-1) = lonw(i)
    end do

    latlon%regional = .true.
    latlon%edges(1) = edgen
    latlon%edges(2) = edgee
    latlon%edges(3) = edges
    latlon%edges(4) = edgew

  end subroutine celledge

end module RtmMapMod



















