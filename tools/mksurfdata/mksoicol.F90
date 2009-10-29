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
  character(len=256) locfn                  ! local dataset file name
  character(len=32) :: subname = 'mksoicol'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make soil color classes .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call getfil (fname, locfn, 0)

  call read_domain(tdomain,locfn)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(locfn, 0, ncid), subname)

  allocate(soil_color_i(nlon_i,nlat_i), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'SOIL_COLOR', varid), subname)
  call check_ret(nf_get_var_int (ncid, varid, soil_color_i), subname)

  call check_ret(nf_close(ncid), subname)

  nsoicol = maxval(soil_color_i)
  write(6,*)'nsoicol = ',nsoicol

  allocate(wst(0:nsoicol),gast_i(0:nsoicol),gast_o(0:nsoicol),col(0:nsoicol))

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

  ! Compute local fields _o

  allocate(mask_i(nlon_i,nlat_i),mask_o(lsmlon,lsmlat))
  allocate( fld_i(nlon_i,nlat_i), fld_o(lsmlon,lsmlat))

  mask_i = 1.0_r8
  mask_o = 1.0_r8
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

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
  deallocate (soil_color_i,wst,gast_i,gast_o,col)
  deallocate (mask_i,mask_o)
  deallocate ( fld_i, fld_o)

end subroutine mksoicol

