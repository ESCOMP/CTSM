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
  use shr_sys_mod , only : shr_sys_flush
  use fileutils   , only : getfil
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
  character(len=256) locfn                  ! local dataset file name
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

  call getfil (fname, locfn, 0)

  call read_domain(tdomain,locfn)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(locfn, 0, ncid), subname)

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
     else                                                  !---need soil
        if (wsti(1) /= 0) then                !---not 'no data'
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

