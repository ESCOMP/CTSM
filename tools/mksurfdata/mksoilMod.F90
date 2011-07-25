module mksoilMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mksoilMod
!
! !DESCRIPTION:
! Make soil data (texture, color and organic)
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!-----------------------------------------------------------------------
!!USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  implicit none

  SAVE
  private           ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public mksoilInit    ! Soil Initialization

  public mksoilAtt     ! Add attributes to output file

  public mksoitex      ! Set soil texture
  public mkorganic     ! Set organic soil
  public mksoicol      ! Set soil color
  public mkfmax        ! Make percent fmax
!
! !PUBLIC DATA MEMBERS:
!
  real(r8), public, parameter :: unset = -999.99_r8 ! Flag to signify soil texture override not set
  real(r8), public            :: soil_sand = unset  ! soil texture sand % to override with
  real(r8), public            :: soil_clay = unset  ! soil texture clay % to override with
  real(r8), public            :: soil_fmax = unset  ! soil max saturation frac to override with
  integer, parameter :: unsetcol   = -999      ! flag to indicate soil color NOT set
  integer, public    :: soil_color = unsetcol  ! soil color to override with
!
! !PRIVATE DATA MEMBERS:
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: mkrank
  private :: mksoitexInit  ! Soil texture Initialization
  private :: mksoicolInit  ! Soil color Initialization
  private :: mksoifmaxInit ! Soil fmax Initialization

!EOP
!===============================================================
contains
!===============================================================

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksoilInit
!
! !INTERFACE:
subroutine mksoilInit( )
!
! !DESCRIPTION:
! Initialize the different soil types
! !USES:
!
! !ARGUMENTS:
  implicit none
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
  character(len=32) :: subname = 'mksoilInit'
!-----------------------------------------------------------------------
  call mksoitexInit()
  call mksoicolInit()
  call mksoifmaxInit()

end subroutine mksoilInit

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksoitexInit
!
! !INTERFACE:
subroutine mksoitexInit( )
!
! !DESCRIPTION:
! Initialize of make soil texture
! !USES:
!
! !ARGUMENTS:
  implicit none
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8) :: sumtex
  character(len=32) :: subname = 'mksoitexInit'
!-----------------------------------------------------------------------
    if ( soil_clay /= unset )then
       write(6,*) 'Replace soil clay % for all points with: ', soil_clay
       if ( soil_sand == unset )then
           write (6,*) subname//':error: soil_clay set, but NOT soil_sand'
           call abort()
       end if
    end if
    if ( soil_sand /= unset )then
       write(6,*) 'Replace soil sand % for all points with: ', soil_sand
       if ( soil_clay == unset )then
           write (6,*) subname//':error: soil_sand set, but NOT soil_clay'
           call abort()
       end if
       sumtex = soil_sand + soil_clay
       if ( sumtex < 0.0_r8 .or. sumtex > 100.0_r8 )then
           write (6,*) subname//':error: soil_sand and soil_clay out of bounds: sand, clay = ', &
                       soil_sand, soil_clay
           call abort()
       end if
    end if

end subroutine mksoitexInit

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksoitex
!
! !INTERFACE:
subroutine mksoitex(lsmlon, lsmlat, fname, ndiag, pctglac_o, sand_o, clay_o)
!
! !DESCRIPTION:
! make %sand and %clay from IGBP soil data, which includes
! igbp soil 'mapunits' and their corresponding textures
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use domainMod   , only : domain_type,domain_clean,domain_setptrs
  use creategridMod, only : read_domain
  use mkvarpar	
  use mkvarsur    , only : ldomain
  use mkvarctl    
  use areaMod     , only : gridmap_type,gridmap_clean,gridmap_setptrs
  use areaMod     , only : areaini,areaave
  use ncdio
!
! !ARGUMENTS:
  implicit none
  integer , intent(in) :: lsmlon, lsmlat          ! clm grid resolution
  character(len=*), intent(in) :: fname           ! input dataset file name
  integer , intent(in) :: ndiag                   ! unit number for diag out
  real(r8), intent(in) :: pctglac_o(lsmlon,lsmlat)      ! % glac (output grid)
  real(r8), intent(out):: sand_o(lsmlon,lsmlat,nlevsoi) ! % sand (output grid)
  real(r8), intent(out):: clay_o(lsmlon,lsmlat,nlevsoi) ! % clay (output grid)
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
  integer  :: nlon_i                        ! input grid : lon points
  integer  :: nlat_i                        ! input grid : lat points

  type(domain_type)     :: tdomain          ! local domain

  type(gridmap_type)    :: tgridmap         ! local gridmap
  integer          :: mxovr                 ! total num of overlapping cells
  integer ,pointer :: novr(:,:)             ! number of overlapping input cells
  integer ,pointer :: iovr(:,:,:)           ! lon index of overlap input cell
  integer ,pointer :: jovr(:,:,:)           ! lat index of overlap input cell
  real(r8),pointer :: wovr(:,:,:)           ! weight    of overlap input cell

  integer, parameter :: nlsm=4              ! number of soil textures 
  character(len=38)  :: soil(0:nlsm)        ! name of each soil texture
  character(len=38)  :: typ                 ! soil texture based on ...
  integer  :: nlay                          ! number of soil layers
  integer  :: mapunitmax                    ! max value of igbp soil mapunits
  integer  :: mapunittemp                   ! temporary igbp soil mapunit

  integer, parameter :: num=2               ! set soil mapunit number
  integer  :: wsti(num)                     ! index to 1st and 2nd largest wst
  real(r8),allocatable :: wst(:)            ! overlap weights, by soil mapunit
  real(r8) :: gast_i(0:nlsm)                ! global area, by texture type
  real(r8) :: gast_o(0:nlsm)                ! global area, by texture type
  real(r8), allocatable :: sand_i(:,:)      ! input grid: percent sand
  real(r8), allocatable :: clay_i(:,:)      ! input grid: percent clay
  real(r8), allocatable :: mask_i(:,:)      ! input grid: mask (0, 1)
  real(r8), allocatable :: mask_o(:,:)      ! output grid: mask (0, 1)
  real(r8), allocatable :: fld_i(:,:)       ! input grid: dummy field
  real(r8), allocatable :: fld_o(:,:)       ! output grid: dummy field
  real(r8), allocatable :: mapunit_i(:,:)   ! input grid: igbp soil mapunits
  real(r8), allocatable :: mapunit_o(:,:)   ! output grid: igbp soil mapunits
  real(r8) :: sum_fldi                      ! global sum of dummy input fld
  real(r8) :: sum_fldo                      ! global sum of dummy output fld

  integer  :: ii,ji                         ! indices
  integer  :: io,jo                         ! indices
  integer  :: k,l,n,m                       ! indices
  integer  :: ncid,dimid,varid              ! input netCDF id's
  integer  :: ier                           ! error status
  integer  :: miss = 99999                  ! missing data indicator
  real(r8) :: relerr = 0.00001              ! max error: sum overlap wts ne 1
  integer  :: t1,t2,t3,t4                   ! timers
  character(len=32) :: subname = 'mksoitex'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make %sand and %clay .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Define the model surface types: 0 to nlsm
  ! -----------------------------------------------------------------

  soil(0) = 'no soil: ocean, glacier, lake, no data'
  soil(1) = 'clays                                 '
  soil(2) = 'sands                                 '
  soil(3) = 'loams                                 '
  soil(4) = 'silts                                 '

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call read_domain(tdomain,fname)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(fname, 0, ncid), subname)

  call check_ret(nf_inq_dimid  (ncid, 'number_of_layers', dimid), subname)
  call check_ret(nf_inq_dimlen (ncid, dimid, nlay), subname)

  call check_ret(nf_inq_dimid  (ncid, 'max_value_mapunit', dimid), subname)
  call check_ret(nf_inq_dimlen (ncid, dimid, mapunitmax), subname)

  allocate(sand_i(mapunitmax,nlay), clay_i(mapunitmax,nlay), stat=ier)
  if (ier/=0) call abort()
  allocate(mapunit_i(nlon_i,nlat_i),mapunit_o(lsmlon,lsmlat), stat=ier)
  if (ier/=0) call abort()
  allocate(wst(0:mapunitmax), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'MAPUNITS', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, mapunit_i), subname)

  call check_ret(nf_inq_varid (ncid, 'PCT_SAND', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, sand_i), subname)

  call check_ret(nf_inq_varid (ncid, 'PCT_CLAY', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, clay_i), subname)

  call check_ret(nf_close(ncid), subname)

  ! Compute local fields _o

  allocate(mask_i(nlon_i,nlat_i),mask_o(lsmlon,lsmlat))
  allocate( fld_i(nlon_i,nlat_i), fld_o(lsmlon,lsmlat))

  mask_i = 1.0_r8
  mask_o = 1.0_r8
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  call gridmap_setptrs(tgridmap,mx_ovr=mxovr,n_ovr=novr,i_ovr=iovr,j_ovr=jovr,w_ovr=wovr)

  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)
     do k = 0, mapunitmax
        wst(k) = 0
     enddo
     do n = 1, novr(io,jo)
        ii = iovr(io,jo,n)
        ji = jovr(io,jo,n)
        k = mapunit_i(ii,ji) * tdomain%mask(ii,ji)
        wst(k) = wst(k) + wovr(io,jo,n)
     enddo

     ! Rank non-zero weights by soil mapunit.
     ! wsti(1) is the most extensive mapunit.
     ! wsti(2) is the second most extensive mapunit.

     call mkrank (mapunitmax, wst, miss, wsti, num)

     ! Set soil texture as follows:
     ! If land grid cell is ocean or 100% glacier: cell has no soil
     ! Otherwise, grid cell needs soil:
     !   a. Use dominant igbp soil mapunit based on area of overlap unless
     !     'no data' is dominant
     !   b. In this case use second most dominant mapunit if it has data
     !   c. If this has no data or if there isn't a second most dominant
     !      mapunit, use loam for soil texture

     if (abs(pctglac_o(io,jo)-100.) < 1.e-06) then    !---glacier
        mapunit_o(io,jo) = 0.
        do l = 1, nlay
           sand_o(io,jo,l) = 0.
           clay_o(io,jo,l) = 0.
        end do
     else                                             !---need soil

                                                      !---if soil texture is input
        if ( soil_sand /= unset .and. soil_clay /= unset ) then  
           do l = 1, nlay
              sand_o(io,jo,l) = soil_sand
              clay_o(io,jo,l) = soil_clay
           end do
        else if (wsti(1) /= 0) then           !---not 'no data'
           mapunit_o(io,jo) = wsti(1)
           do l = 1, nlay
              sand_o(io,jo,l) = sand_i(wsti(1),l)
              clay_o(io,jo,l) = clay_i(wsti(1),l)
           end do
        else                                  !---if (wsti(1) == 0) then
           if (wsti(2) == 0 .or. wsti(2) == miss) then     !---no data
              mapunit_o(io,jo) = wsti(2)
              do l = 1, nlay
                 sand_o(io,jo,l) = 43.        !---use loam
                 clay_o(io,jo,l) = 18.
              end do
           else                               !---if (wsti(2) /= 0 and /= miss)
              mapunit_o(io,jo) = wsti(2)
              do l = 1, nlay
                 sand_o(io,jo,l) = sand_i(wsti(2),l)
                 clay_o(io,jo,l) = clay_i(wsti(2),l)
              end do
           end if       !---end of wsti(2) if-block
        end if          !---end of wsti(1) if-block
     end if             !---end of land/ocean if-block
  enddo
  enddo

  ! Global sum of output field 

  sum_fldi = 0.0_r8
  do ji = 1, nlat_i
  do ii = 1, nlon_i
    fld_i(ii,ji) = ((ji-1)*nlon_i + ii) * tdomain%mask(ii,ji)
    sum_fldi = sum_fldi + tdomain%area(ii,ji) * fld_i(ii,ji)
  enddo
  enddo

  call areaave(fld_i,fld_o,tgridmap)

  sum_fldo = 0.
  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)
     fld_o(io,jo) = fld_o(io,jo)*mask_o(io,jo)
     sum_fldo = sum_fldo + ldomain%area(io,jo) * fld_o(io,jo)
  end do
  end do

  ! -----------------------------------------------------------------
  ! Error check1
  ! Compare global sum fld_o to global sum fld_i.
  ! -----------------------------------------------------------------

  if ( trim(mksrf_gridtype) == 'global') then
     if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
        write (6,*) 'MKSOITEX error: input field not conserved'
        write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
        write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
        stop
     end if
  end if

  ! -----------------------------------------------------------------
  ! Error check2
  ! Compare global area of each soil type on input and output grids
  ! -----------------------------------------------------------------

  mask_i = float(tdomain%mask(:,:))
  call areaave(mask_i,mask_o,tgridmap)
  ldomain%frac = mask_o

  call gridmap_clean(tgridmap)
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  ! input grid: global areas by texture class

  gast_i(:) = 0.
  do l = 1, nlay
     do ji = 1, nlat_i
        do ii = 1, nlon_i
           mapunittemp = nint(mapunit_i(ii,ji))
           if (mapunittemp==0) then
              typ = 'no soil: ocean, glacier, lake, no data'
           else if (clay_i(mapunittemp,l) >= 40.) then
              typ = 'clays'
           else if (sand_i(mapunittemp,l) >= 50.) then
              typ = 'sands'
           else if (clay_i(mapunittemp,l)+sand_i(mapunittemp,l) < 50.) then
              if (tdomain%mask(ii,ji) /= 0.) then
                 typ = 'silts'
              else            !if (tdomain%mask(ii,ji) == 0.) then no data
                 typ = 'no soil: ocean, glacier, lake, no data'
              end if
           else
              typ = 'loams'
           end if
           do m = 0, nlsm
              if (typ == soil(m)) go to 101
           end do
           write (6,*) 'MKSOITEX error: sand = ',sand_i(mapunittemp,l), &
             ' clay = ',clay_i(mapunittemp,l), &
             ' not assigned to soil type for input grid lon,lat,layer = ',ii,ji,l
           call abort()
101        continue
           gast_i(m) = gast_i(m) + tdomain%area(ii,ji) * tdomain%frac(ii,ji)
        end do
     end do
  end do

  ! output grid: global areas by texture class

  gast_o(:) = 0.
  do l = 1, nlay
     do jo = 1, ldomain%nj
        do io = 1, ldomain%numlon(jo)
           if (clay_o(io,jo,l)==0. .and. sand_o(io,jo,l)==0.) then
              typ = 'no soil: ocean, glacier, lake, no data'
           else if (clay_o(io,jo,l) >= 40.) then
              typ = 'clays'
           else if (sand_o(io,jo,l) >= 50.) then
              typ = 'sands'
           else if (clay_o(io,jo,l)+sand_o(io,jo,l) < 50.) then
              typ = 'silts'
           else
              typ = 'loams'
           end if
           do m = 0, nlsm
              if (typ == soil(m)) go to 102
           end do
           write (6,*) 'MKSOITEX error: sand = ',sand_o(io,jo,l), &
             ' clay = ',clay_o(io,jo,l), &
             ' not assigned to soil type for output grid lon,lat,layer = ',io,jo,l
           call abort()
102        continue
           gast_o(m) = gast_o(m) + ldomain%area(io,jo) * ldomain%frac(io,jo)
        end do
     end do
  end do

  ! Diagnostic output

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',l=1,70)
  write (ndiag,*) 'Soil Texture Output'
  write (ndiag,'(1x,70a1)') ('=',l=1,70)
  write (ndiag,*)

  write (ndiag,*) 'The following table of soil texture classes is for comparison only.'
  write (ndiag,*) 'The actual data is continuous %sand, %silt and %clay not textural classes'
  write (ndiag,*)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',l=1,70)
  write (ndiag,1001)
1001 format (1x,'soil texture class',17x,' input grid area output grid area',/ &
             1x,33x,'     10**6 km**2','      10**6 km**2')
  write (ndiag,'(1x,70a1)') ('.',l=1,70)
  write (ndiag,*)

  do l = 0, nlsm
     write (ndiag,1002) soil(l),gast_i(l)*1.e-6,gast_o(l)*1.e-6
1002 format (1x,a38,f16.3,f17.3)
  end do

  if (lsmlat > 1) then
     k = lsmlat/2
     write (ndiag,*)
     write (ndiag,*) 'For reference the area on the output grid of a cell near the equator is: '
     write (ndiag,'(f10.3,a14)')ldomain%area(1,k)*1.e-06,' x 10**6 km**2'
     write (ndiag,*)
  endif

  write (6,*) 'Successfully made %sand and %clay'
  write (6,*)
  call shr_sys_flush(6)

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (sand_i,clay_i,mapunit_i,mapunit_o,wst)
  deallocate (mask_i,mask_o)
  deallocate ( fld_i, fld_o)

end subroutine mksoitex

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksoicolInit
!
! !INTERFACE:
subroutine mksoicolInit( )
!
! !DESCRIPTION:
! Initialize of make soil color
! !USES:
!
! !ARGUMENTS:
  implicit none
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8) :: sumtex
  character(len=32) :: subname = 'mksoicolInit'
!-----------------------------------------------------------------------

  ! Error check soil_color if it is set
  if ( soil_color /= unsetcol )then
     if ( soil_color < 0 .or. soil_color > 20 )then
        write(6,*)'soil_color is out of range = ', soil_color
        call abort()
     end if
     write(6,*) 'Replace soil color for all points with: ', soil_color
  end if
end subroutine mksoicolInit


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksoicol
!
! !INTERFACE:
subroutine mksoicol(lsmlon, lsmlat, fname, ndiag, pctglac_o, soil_color_o, nsoicol)
!
! !DESCRIPTION:
! make %sand and %clay from IGBP soil data, which includes
! igbp soil 'mapunits' and their corresponding textures
!
! !USES:
  use domainMod   , only : domain_type,domain_clean,domain_setptrs
  use creategridMod, only : read_domain
  use mkvarpar	
  use mkvarsur    , only : ldomain
  use mkvarctl    
  use areaMod     , only : gridmap_type,gridmap_clean,gridmap_setptrs
  use areaMod     , only : areaini,areaave
  use ncdio
!
! !ARGUMENTS:
  implicit none
  integer , intent(in) :: lsmlon, lsmlat          ! clm grid resolution
  character(len=*), intent(in) :: fname           ! input dataset file name
  integer , intent(in) :: ndiag                   ! unit number for diag out
  real(r8), intent(in) :: pctglac_o(lsmlon,lsmlat)     ! % glac (output grid)
  integer , intent(out):: soil_color_o(lsmlon,lsmlat)  ! soil color classes
  integer , intent(out):: nsoicol                      ! number of soil colors 
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
  integer  :: nlon_i                        ! input grid : lon points
  integer  :: nlat_i                        ! input grid : lat points

  type(domain_type)     :: tdomain          ! local domain

  type(gridmap_type)    :: tgridmap         ! local gridmap
  integer          :: mxovr                 ! total num of overlapping cells
  integer ,pointer :: novr(:,:)             ! number of overlapping input cells
  integer ,pointer :: iovr(:,:,:)           ! lon index of overlap input cell
  integer ,pointer :: jovr(:,:,:)           ! lat index of overlap input cell
  real(r8),pointer :: wovr(:,:,:)           ! weight    of overlap input cell


  integer, parameter :: num=2               ! set soil mapunit number
  integer  :: wsti(num)                     ! index to 1st and 2nd largest wst
  real(r8), allocatable :: wst(:)           ! overlap weights, by surface type
  real(r8), allocatable :: gast_i(:)        ! global area, by surface type
  real(r8), allocatable :: gast_o(:)        ! global area, by surface type
  integer , allocatable :: soil_color_i(:,:)! input grid: BATS soil color
  real(r8), allocatable :: mask_i(:,:)      ! input grid: mask (0, 1)
  real(r8), allocatable :: mask_o(:,:)      ! output grid: mask (0, 1)
  real(r8), allocatable :: fld_i(:,:)       ! input grid: dummy field
  real(r8), allocatable :: fld_o(:,:)       ! output grid: dummy field
  real(r8) :: sum_fldi                      ! global sum of dummy input fld
  real(r8) :: sum_fldo                      ! global sum of dummy output fld
  character(len=35), allocatable :: col(:)  ! name of each color
  integer  :: color                         ! 0: none; 1: some

  integer  :: ii,ji                         ! indices
  integer  :: io,jo                         ! indices
  integer  :: k,l,n,m                       ! indices
  integer  :: ncid,dimid,varid              ! input netCDF id's
  integer  :: ier                           ! error status
  integer  :: miss = 99999                  ! missing data indicator
  real(r8) :: relerr = 0.00001              ! max error: sum overlap wts ne 1
  integer  :: t1,t2,t3,t4                   ! timers
  character(len=32) :: subname = 'mksoicol'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make soil color classes .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call read_domain(tdomain,fname)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(fname, 0, ncid), subname)

  allocate(soil_color_i(nlon_i,nlat_i), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'SOIL_COLOR', varid), subname)
  call check_ret(nf_get_var_int (ncid, varid, soil_color_i), subname)

  call check_ret(nf_close(ncid), subname)

  nsoicol = maxval(soil_color_i)
  write(6,*)'nsoicol = ',nsoicol

  allocate(gast_i(0:nsoicol),gast_o(0:nsoicol),col(0:nsoicol))

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

  ! Error check soil_color if it is set
  if ( soil_color /= unsetcol )then
     if ( soil_color > nsoicol )then
        write(6,*)'soil_color is out of range = ', soil_color
        call abort()
     end if
  end if

  allocate(mask_i(nlon_i,nlat_i),mask_o(lsmlon,lsmlat))
  allocate( fld_i(nlon_i,nlat_i), fld_o(lsmlon,lsmlat))

  mask_i = 1.0_r8
  mask_o = 1.0_r8
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  if ( soil_color /= unsetcol )then
     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        soil_color_o(io,jo) = soil_color
     end do
     end do
  else
     allocate(wst(0:nsoicol))

     ! Compute local fields _o

     call gridmap_setptrs(tgridmap,mx_ovr=mxovr,n_ovr=novr,i_ovr=iovr,j_ovr=jovr,w_ovr=wovr)

     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        color = 0
        do k = 0, nsoicol
           wst(k) = 0
        enddo
        do n = 1, novr(io,jo)
           ii = iovr(io,jo,n)
           ji = jovr(io,jo,n)
           k = soil_color_i(ii,ji) * tdomain%mask(ii,ji)
           wst(k) = wst(k) + wovr(io,jo,n)
           if (k>0 .and. wst(k)>0.) color = 1
        enddo
        if (color == 1) wst(0) = 0.0

        ! Rank non-zero weights by color type. wsti(1) is the most extensive
        ! color type. wsti(2) is the second most extensive color type

        call mkrank (nsoicol, wst, miss, wsti, num)
        soil_color_o(io,jo) = wsti(1)
        ! Set everything to input soil_color if it's set


        ! If land but no color, set color to 15 (in older dataset generic 
        ! soil color 4)

        if (nsoicol == 8) then
           if (soil_color_o(io,jo)==0) &
              soil_color_o(io,jo) = 4
        else if (nsoicol == 20) then
           if (soil_color_o(io,jo)==0) &
              soil_color_o(io,jo) = 15
        end if

        ! Set color for grid cells that are 100% glacier to zero. Otherwise,
        ! must have a soil color for the non-glacier portion of grid cell.

        if (abs(pctglac_o(io,jo)-100.)<1.e-06) soil_color_o(io,jo)=0

        ! Error checks

        if (soil_color_o(io,jo) < 0 .or. soil_color_o(io,jo) > nsoicol) then
           write (6,*) 'MKSOICOL error: land model soil color = ', &
                soil_color_o(io,jo),' is not valid for lon,lat = ',io,jo
           call abort()
        end if

     enddo
     enddo
     deallocate (wst)

  end if

  ! Global sum of output field 

  sum_fldi = 0.0_r8
  do ji = 1, nlat_i
  do ii = 1, nlon_i
    fld_i(ii,ji) = ((ji-1)*nlon_i + ii) * tdomain%mask(ii,ji)
    sum_fldi = sum_fldi + tdomain%area(ii,ji) * fld_i(ii,ji)
  enddo
  enddo

  call areaave(fld_i,fld_o,tgridmap)

  sum_fldo = 0.
  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)
     fld_o(io,jo) = fld_o(io,jo)*mask_o(io,jo)
     sum_fldo = sum_fldo + ldomain%area(io,jo) * fld_o(io,jo)
  end do
  end do

  ! -----------------------------------------------------------------
  ! Error check1
  ! Compare global sum fld_o to global sum fld_i.
  ! -----------------------------------------------------------------

  if ( trim(mksrf_gridtype) == 'global') then
     if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
        write (6,*) 'MKSOICOL error: input field not conserved'
        write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
        write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
        stop
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
        gast_i(k) = gast_i(k) + tdomain%area(ii,ji) * tdomain%frac(ii,ji)
     end do
  end do

  ! output grid

  mask_i = float(tdomain%mask(:,:))
  call areaave(mask_i,mask_o,tgridmap)
  ldomain%frac = mask_o

  call gridmap_clean(tgridmap)
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  gast_o(:) = 0.
  do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        k = soil_color_o(io,jo)
        gast_o(k) = gast_o(k) + ldomain%area(io,jo) * ldomain%frac(io,jo)
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
     write (ndiag,'(f10.3,a14)')ldomain%area(1,k)*1.e-06,' x 10**6 km**2'
     write (ndiag,*)
  endif

  write (6,*) 'Successfully made soil color classes'
  write (6,*)
  call shr_sys_flush(6)

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (soil_color_i,gast_i,gast_o,col)
  deallocate (mask_i,mask_o)
  deallocate ( fld_i, fld_o)

end subroutine mksoicol

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkorganic
!
! !INTERFACE:
subroutine mkorganic(lsmlon, lsmlat, fname, ndiag, organic_o)
!
! !DESCRIPTION:
! make organic matter dataset
!
! !USES:
  use domainMod   , only : domain_type,domain_clean,domain_setptrs
  use creategridMod, only : read_domain
  use mkvarpar	
  use mkvarsur    , only : ldomain
  use mkvarctl    
  use areaMod     , only : areaini,areaave,gridmap_type,gridmap_clean  
  use ncdio
!
! !ARGUMENTS:
  implicit none
  integer , intent(in) :: lsmlon, lsmlat          ! clm grid resolution
  character(len=*), intent(in) :: fname           ! input dataset file name
  integer , intent(in) :: ndiag                   ! unit number for diag out
  real(r8), intent(out):: organic_o(lsmlon,lsmlat,nlevsoi)    ! output grid:
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! 
! Author: David Lawrence
!
!
! !LOCAL VARIABLES:
!EOP
  integer  :: nlon_i                          ! input grid : lon points
  integer  :: nlat_i                          ! input grid : lat points

  type(domain_type)     :: tdomain            ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap

  real(r8), allocatable :: organic_i(:,:,:)  ! input grid: total column organic matter
  real(r8), allocatable :: mask_i(:,:)        ! input grid: mask (0, 1)
  real(r8), allocatable :: mask_o(:,:)        ! output grid: mask (0, 1)
  real(r8), allocatable :: omlev_i(:,:)       ! input grid: organic matter on lev
  real(r8), allocatable :: omlev_o(:,:)       ! output grid: organic matter on lev
  real(r8) :: sum_fldi                        ! global sum of dummy input fld
  real(r8) :: sum_fldo                        ! global sum of dummy output fld
  real(r8) :: gomlev_i                         ! input  grid: global organic on lev
  real(r8) :: garea_i                         ! input  grid: global area
  real(r8) :: gomlev_o                         ! output grid: global organic on lev
  real(r8) :: garea_o                         ! output grid: global area

  integer  :: ii,ji                           ! indices
  integer  :: io,jo                           ! indices
  integer  :: k,n,m                           ! indices
  integer  :: lev                             ! level index
  integer   :: nlay                            ! number of soil layers
  integer  :: ncid,dimid,varid                ! input netCDF id's
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001                ! max error: sum overlap wts ne 1
  character(len=32) :: subname = 'mkorganic'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make organic matter dataset .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call read_domain(tdomain,fname)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(fname, 0, ncid), subname)

  call check_ret(nf_inq_dimid  (ncid, 'number_of_layers', dimid), subname)
  call check_ret(nf_inq_dimlen (ncid, dimid, nlay), subname)

  allocate(organic_i(nlon_i,nlat_i,nlay), omlev_i(nlon_i,nlat_i), &
           omlev_o(lsmlon,lsmlat),stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'ORGANIC', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, organic_i), subname)

  call check_ret(nf_close(ncid), subname)

  ! Compute local fields _o

  allocate(mask_i(nlon_i,nlat_i),mask_o(lsmlon,lsmlat))

  mask_i = 1.0_r8
  mask_o = 1.0_r8

  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  mask_i = float(tdomain%mask(:,:))
  call areaave(mask_i,mask_o,tgridmap)

  call gridmap_clean(tgridmap)
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  ldomain%frac = mask_o

  ! Area-average percent cover on input grid to output grid 
  ! and correct according to land landmask
  ! Note that percent cover is in terms of total grid area.

  call areaave(organic_i,organic_o,tgridmap)

  do lev = 1,nlevsoi

     ! Check for conservation

     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        if ((organic_o(io,jo,lev)) > 130.000001_r8) then
           write (6,*) 'MKORGANIC error: organic = ',organic_o(io,jo,lev), &
                ' greater than 130.000001 for column, row = ',io,jo
           call shr_sys_flush(6)
           stop
        end if
     enddo
     enddo

!    ! Diagnostic output

!    write (ndiag,*)
!    write (ndiag,'(1x,70a1)') ('.',k=1,70)
!    write (ndiag,2001)
!2001 format (1x,'surface type   input grid area  output grid area'/ &
!            1x,'                 10**6 km**2      10**6 km**2   ')
!    write (ndiag,'(1x,70a1)') ('.',k=1,70)
!    write (ndiag,*)
!    write (ndiag,2002) gomlev_i*1.e-06,gomlev_o*1.e-06
!    write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
!2002 format (1x,'organic    ',f14.3,f17.3)
!2004 format (1x,'all surface ',f14.3,f17.3)

!    if (lsmlat > 1) then
!       k = lsmlat/2
!       write (ndiag,*)
!       write (ndiag,*) 'For reference the area on the output grid of a cell near the equator is: '
!       write (ndiag,'(f10.3,a14)')ldomain%area(1,k)*1.e-06,' x 10**6 km**2'
!       write (ndiag,*)
!    endif
     call shr_sys_flush(ndiag)

     write (6,*) 'Successfully made organic matter, level = ', lev
     write (6,*)
     call shr_sys_flush(6)

  end do   ! lev

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (organic_i, omlev_i, omlev_o)
  deallocate (mask_i,mask_o)

end subroutine mkorganic

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: mkrank
!
! !INTERFACE:
subroutine mkrank (n, a, miss, iv, num)
!
! !DESCRIPTION:
! Return indices of largest [num] values in array [a]. Private method
! only used for soil color and soil texture.
!
! !USES:
!
! !ARGUMENTS:
  implicit none
  integer , intent(in) :: n        !array length
  real(r8), intent(in) :: a(0:n)   !array to be ranked
  integer , intent(in) :: miss     !missing data value
  integer , intent(in) :: num      !number of largest values requested
  integer , intent(out):: iv(num)  !index to [num] largest values in array [a]
!
! !CALLED FROM:
! subroutine mkpft
! subroutine mksoicol
! subroutine mksoitex
!
! !REVISION HISTORY:
! Author: Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8) a_max       !maximum value in array
  integer i            !array index
  real(r8) delmax      !tolerance for finding if larger value
  integer m            !do loop index
  integer k            !do loop index
  logical exclude      !true if data value has already been chosen
!-----------------------------------------------------------------------

  delmax = 1.e-06

  ! Find index of largest non-zero number

  iv(1) = miss
  a_max = -9999.

  do i = 0, n
     if (a(i)>0. .and. (a(i)-a_max)>delmax) then
        a_max = a(i)
        iv(1)  = i
     end if
  end do

  ! iv(1) = miss indicates no values > 0. this is an error

  if (iv(1) == miss) then
     write (6,*) 'MKRANK error: iv(1) = missing'
     call abort()
  end if

  ! Find indices of the next [num]-1 largest non-zero number.
  ! iv(m) = miss if there are no more values > 0

  do m = 2, num
     iv(m) = miss
     a_max = -9999.
     do i = 0, n

        ! exclude if data value has already been chosen

        exclude = .false.
        do k = 1, m-1
           if (i == iv(k)) exclude = .true.
        end do

        ! if not already chosen, see if it is the largest of
        ! the remaining values

        if (.not. exclude) then
           if (a(i)>0. .and. (a(i)-a_max)>delmax) then
              a_max = a(i)
              iv(m)  = i
           end if
        end if
     end do
  end do

  return
end subroutine mkrank

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksoifmaxInit
!
! !INTERFACE:
subroutine mksoifmaxInit( )
!
! !DESCRIPTION:
! Initialize of make soil fmax
! !USES:
!
! !ARGUMENTS:
  implicit none
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8) :: sumtex
  character(len=32) :: subname = 'mksoifmaxInit'
!-----------------------------------------------------------------------

  ! Error check soil_fmax if it is set
  if ( soil_fmax /= unset )then
     if ( soil_fmax < 0.0 .or. soil_fmax > 1.0 )then
        write(6,*)'soil_fmax is out of range = ', soil_fmax
        call abort()
     end if
     write(6,*) 'Replace soil fmax for all points with: ', soil_fmax
  end if
end subroutine mksoifmaxInit

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkfmax
!
! !INTERFACE:
subroutine mkfmax(lsmlon, lsmlat, fname, ndiag, fmax_o)
!
! !DESCRIPTION:
! make percent fmax
!
! !USES:
  use domainMod   , only : domain_type,domain_clean,domain_setptrs
  use creategridMod, only : read_domain
  use mkvarpar	
  use mkvarsur    , only : ldomain
  use mkvarctl    
  use areaMod     , only : areaini,areaave,gridmap_type,gridmap_clean  
  use ncdio
!
! !ARGUMENTS:
  implicit none
  integer , intent(in) :: lsmlon, lsmlat          ! clm grid resolution
  character(len=*), intent(in) :: fname           ! input dataset file name
  integer , intent(in) :: ndiag                   ! unit number for diag out
  real(r8), intent(out):: fmax_o(lsmlon,lsmlat)    ! output grid: %fmax
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Revised: Nan Rosenbloom - used mkglacier.F90 as template.
! Original Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
  integer  :: nlon_i                          ! input grid : lon points
  integer  :: nlat_i                          ! input grid : lat points

  type(domain_type)     :: tdomain            ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap

  real(r8), allocatable :: fmax_i(:,:)        ! input grid: percent fmax
  real(r8), allocatable :: mask_i(:,:)        ! input grid: mask (0, 1)
  real(r8), allocatable :: mask_o(:,:)        ! output grid: mask (0, 1)
  real(r8), allocatable :: fld_i(:,:)         ! input grid: dummy field
  real(r8), allocatable :: fld_o(:,:)         ! output grid: dummy field
  real(r8) :: sum_fldi                        ! global sum of dummy input fld
  real(r8) :: sum_fldo                        ! global sum of dummy output fld
  real(r8) :: gfmax_i                         ! input  grid: global fmax
  real(r8) :: garea_i                         ! input  grid: global area
  real(r8) :: gfmax_o                         ! output grid: global fmax
  real(r8) :: garea_o                         ! output grid: global area

  integer  :: ii,ji                           ! indices
  integer  :: io,jo                           ! indices
  integer  :: k,n,m                           ! indices
  integer  :: ncid,dimid,varid                ! input netCDF id's
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001                ! max error: sum overlap wts ne 1
  character(len=32) :: subname = 'mkfmax'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make %fmax .....'
  call shr_sys_flush(6)

  allocate( fld_o(lsmlon,lsmlat))

  if ( soil_fmax  /= unset )then
     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        fmax_o(io,jo) = soil_fmax
     end do
     end do
  else
     ! -----------------------------------------------------------------
     ! Read input file
     ! -----------------------------------------------------------------

     ! Obtain input grid info, read local fields

     call read_domain(tdomain,fname)
     call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

     call check_ret(nf_open(fname, 0, ncid), subname)

     allocate(fmax_i(nlon_i,nlat_i), stat=ier)
     if (ier/=0) call abort()

     call check_ret(nf_inq_varid (ncid, 'FMAX', varid), subname)
     call check_ret(nf_get_var_double (ncid, varid, fmax_i), subname)

     call check_ret(nf_close(ncid), subname)

     ! Compute local fields _o

     allocate(mask_i(nlon_i,nlat_i),mask_o(lsmlon,lsmlat))
     allocate( fld_i(nlon_i,nlat_i) )

     mask_i = 1.0_r8
     mask_o = 1.0_r8
     call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

     mask_i = float(tdomain%mask(:,:))
     call areaave(mask_i,mask_o,tgridmap)

     call gridmap_clean(tgridmap)
     call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)
     ldomain%frac = mask_o

     ! Area-average percent cover on input grid to output grid 
     ! and correct according to land landmask
     ! Note that percent cover is in terms of total grid area.

     call areaave(fmax_i,fmax_o,tgridmap)

     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        if (fmax_o(io,jo) == 0.0) then
           fmax_o(io,jo) = .365783  ! fmax_o(io,jo) = globalAvg
        end if
        if (fmax_o(io,jo) == -999.99) then
            fmax_o(io,jo) = .365783 ! fmax_o(io,jo) = globalAvg
        end if
     enddo
     enddo

     ! Check for conservation

     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        if ((fmax_o(io,jo)) > 1.000001_r8) then
           write (6,*) 'MKFMAX error: fmax = ',fmax_o(io,jo), &
                   ' greater than 1.000001 for column, row = ',io,jo
           call shr_sys_flush(6)
           stop
        end if
     enddo
     enddo

     ! Global sum of output field -- must multiply by fraction of
     ! output grid that is land as determined by input grid

     sum_fldi = 0.0_r8
     do ji = 1, nlat_i
     do ii = 1, nlon_i
       fld_i(ii,ji) = ((ji-1)*nlon_i + ii) * tdomain%mask(ii,ji)
       sum_fldi = sum_fldi + tdomain%area(ii,ji) * fld_i(ii,ji)
     enddo
     enddo

     call areaave(fld_i,fld_o,tgridmap)

     sum_fldo = 0.
     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        fld_o(io,jo) = fld_o(io,jo)*mask_o(io,jo)
        sum_fldo = sum_fldo + ldomain%area(io,jo) * fld_o(io,jo)
     end do
     end do

     ! -----------------------------------------------------------------
     ! Error check1
     ! Compare global sum fld_o to global sum fld_i.
     ! -----------------------------------------------------------------

     if ( trim(mksrf_gridtype) == 'global') then
        if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
           write (6,*) 'MKFMAX error: input field not conserved'
           write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
           write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
           stop
        end if
     end if

     ! -----------------------------------------------------------------
     ! Error check2
     ! Compare global areas on input and output grids
     ! -----------------------------------------------------------------

     ! Input grid

     gfmax_i = 0.
     garea_i = 0.

     do ji = 1, nlat_i
     do ii = 1, nlon_i
        garea_i = garea_i + tdomain%area(ii,ji)
        gfmax_i = gfmax_i + fmax_i(ii,ji)*(tdomain%area(ii,ji)/100.) * &
                                           tdomain%frac(ii,ji)
     end do
     end do

     ! Output grid

     gfmax_o = 0.
     garea_o = 0.

     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        garea_o = garea_o + ldomain%area(io,jo)
        gfmax_o = gfmax_o + fmax_o(io,jo)*(ldomain%area(io,jo)/100.) * &
                                           ldomain%frac(io,jo)
        if ( (ldomain%frac(io,jo) < 0.0) .or. (ldomain%frac(io,jo) > 1.0001) )then
           write(6,*) "ERROR:: frac out of range: ", ldomain%frac(io,jo), io, jo
           stop
        end if
     end do
     end do

     ! Diagnostic output

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('=',k=1,70)
     write (ndiag,*) 'Maximum Fractional Saturated Area Output'
     write (ndiag,'(1x,70a1)') ('=',k=1,70)

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/ &
             1x,'                 10**6 km**2      10**6 km**2   ')
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,*)
     write (ndiag,2002) gfmax_i*1.e-06,gfmax_o*1.e-06
     write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2002 format (1x,'fmax    ',f14.3,f17.3)
2004 format (1x,'all surface ',f14.3,f17.3)

     if (lsmlat > 1) then
        k = lsmlat/2
        write (ndiag,*)
        write (ndiag,*) 'For reference the area on the output grid of a cell near the equator is: '
        write (ndiag,'(f10.3,a14)')ldomain%area(1,k)*1.e-06,' x 10**6 km**2'
        write (ndiag,*)
     endif
     call shr_sys_flush(ndiag)

     ! Deallocate dynamic memory

     call domain_clean(tdomain)
     call gridmap_clean(tgridmap)
     deallocate (fmax_i)
     deallocate (mask_i,mask_o)
     deallocate ( fld_i )

  end if

  write (6,*) 'Successfully made %fmax'
  write (6,*)
  call shr_sys_flush(6)

  deallocate ( fld_o)

end subroutine mkfmax

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkfmax
!
! !INTERFACE:
subroutine mksoilAtt( ncid, dynlanduse, xtype )
!
! !DESCRIPTION:
! add atttributes to output file regarding the soil module
!
! !USES:
  use ncdio       , only : check_ret, ncd_defvar
  use mkvarpar    , only : nlevsoi
  use mkvarctl    , only : mksrf_fsoitex, mksrf_fsoicol, mksrf_firrig, &
                           mksrf_forganic, mksrf_fmax
  use fileutils   , only : get_filename
! !ARGUMENTS:
  implicit none
  include 'netcdf.inc'
  integer, intent(in) :: ncid ! NetCDF file ID to write out to
  logical, intent(in) :: dynlanduse   ! if dynamic land-use file
  integer, intent(in) :: xtype        ! external type to output real data as
!
! !CALLED FROM:
! subroutine mkfile in module mkfileMod
!
! !REVISION HISTORY:
! Original Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
  integer :: dimid                ! temporary
  character(len=256) :: str       ! global attribute string
  character(len=32) :: subname = 'mksoilAtt'
!-----------------------------------------------------------------------

  if (.not. dynlanduse) then
  ! Define dimensions unique to soil

     call check_ret(nf_def_dim (ncid, 'nlevsoi',  &
                                       nlevsoi    , dimid), subname)

     ! Add global attributes to file

     if ( soil_clay /= unset .and. soil_sand /= unset )then
        str = 'TRUE'
        call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
          'soil_clay_override', len_trim(str), trim(str)), subname)
        str = 'TRUE'
        call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
          'soil_sand_override', len_trim(str), trim(str)), subname)
     else
        str = get_filename(mksrf_fsoitex)
        call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Soil_texture_raw_data_file_name', len_trim(str), trim(str)), subname)
     end if
     if ( soil_color /= unsetcol )then
        str = 'TRUE'
        call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
          'soil_color_override', len_trim(str), trim(str)), subname)
     else
        str = get_filename(mksrf_fsoicol)
        call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Soil_color_raw_data_file_name', len_trim(str), trim(str)), subname)
     end if
     if ( soil_fmax /= unset )then
        str = 'TRUE'
        call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
          'soil_fmax_override', len_trim(str), trim(str)), subname)
     else
        str = get_filename(mksrf_fmax)
        call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Fmax_raw_data_file_name', len_trim(str), trim(str)), subname)
     end if

     str = get_filename(mksrf_forganic)
     call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Organic_matter_raw_data_file_name', len_trim(str), trim(str)), subname)

     call ncd_defvar(ncid=ncid, varname='mxsoil_color', xtype=nf_int, &
            long_name='maximum numbers of soil colors', units='unitless')

     call ncd_defvar(ncid=ncid, varname='SOIL_COLOR', xtype=nf_int, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='soil color', units='unitless')

     call ncd_defvar(ncid=ncid, varname='PCT_SAND', xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevsoi', &
            long_name='percent sand', units='unitless')

     call ncd_defvar(ncid=ncid, varname='PCT_CLAY', xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevsoi', &
            long_name='percent clay', units='unitless')

     call ncd_defvar(ncid=ncid, varname='ORGANIC', xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevsoi', &
            long_name='organic matter density at soil levels', &
            units='kg/m3 (assumed carbon content 0.58 gC per gOM)')

     call ncd_defvar(ncid=ncid, varname='FMAX', xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='maximum fractional saturated area', units='unitless')

  end if

end subroutine mksoilAtt

!-----------------------------------------------------------------------

end module mksoilMod
