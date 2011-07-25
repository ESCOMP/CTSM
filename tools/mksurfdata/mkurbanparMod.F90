module mkurbanparMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkurbanpar
!
! !DESCRIPTION:
! Make Urban Parameter data
!
! !REVISION HISTORY:
! Author: Keith Oleson
!
!-----------------------------------------------------------------------
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use areaMod     , only : gridmap_type

  implicit none

  private

! !PUBLIC MEMBER FUNCTIONS:
  public  :: mkelev        ! Get elevation to reduce urban for high elevation areas
  public  :: mkurban       ! Get the urban percentage
  public  :: mkurbanpar    ! Make the urban parameters

!EOP
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkurban
!
! !INTERFACE:
subroutine mkurban(lsmlon, lsmlat, fname, ndiag, zero_out, urbn_o)
!
! !DESCRIPTION:
! make percent urban
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
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
  logical , intent(in) :: zero_out                ! if should zero urban out
  real(r8), intent(out):: urbn_o(lsmlon,lsmlat)    ! output grid: %urban
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
  integer  :: nlon_i                          ! input grid : lon points
  integer  :: nlat_i                          ! input grid : lat points

  type(domain_type)     :: tdomain            ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap

  real(r8), allocatable :: urbn_i(:,:)        ! input grid: percent urbn
  real(r8), allocatable :: mask_i(:,:)        ! input grid: mask (0, 1)
  real(r8), allocatable :: mask_o(:,:)        ! output grid: mask (0, 1)
  real(r8), allocatable :: fld_i(:,:)         ! input grid: dummy field
  real(r8), allocatable :: fld_o(:,:)         ! output grid: dummy field
  real(r8) :: sum_fldi                        ! global sum of dummy input fld
  real(r8) :: sum_fldo                        ! global sum of dummy output fld
  real(r8) :: gurbn_i                          ! input  grid: global urbn
  real(r8) :: garea_i                         ! input  grid: global area
  real(r8) :: gurbn_o                          ! output grid: global urbn
  real(r8) :: garea_o                         ! output grid: global area

  integer  :: ii,ji                           ! indices
  integer  :: io,jo                           ! indices
  integer  :: k,n,m                           ! indices
  integer  :: ncid,dimid,varid                ! input netCDF id's
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001_r8             ! max error: sum overlap wts ne 1
  character(len=32) :: subname = 'mkurban'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make %urban .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call read_domain(tdomain,fname)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(fname, 0, ncid), subname)

  allocate(urbn_i(nlon_i,nlat_i), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'PCT_URBAN', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, urbn_i), subname)

  call check_ret(nf_close(ncid), subname)

  ! Compute local fields _o

  allocate(mask_i(nlon_i,nlat_i),mask_o(lsmlon,lsmlat))
  allocate( fld_i(nlon_i,nlat_i), fld_o(lsmlon,lsmlat))

  mask_i = 1.0_r8
  mask_o = 1.0_r8
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  if ( all_urban .or. zero_out )then
     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
           if (all_urban )then
              urbn_o(io,jo) = 100.00_r8
           else 
              urbn_o(io,jo) =   0.00_r8
           end if
     enddo
     enddo
  else
     mask_i = float(tdomain%mask(:,:))
     call areaave(mask_i,mask_o,tgridmap)

     call gridmap_clean(tgridmap)
     call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

     ! Area-average percent cover on input grid to output grid 
     ! and correct according to land landmask
     ! Note that percent cover is in terms of total grid area.

     call areaave(urbn_i,urbn_o,tgridmap)

     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
           if (urbn_o(io,jo) < 0.1_r8) urbn_o(io,jo) = 0._r8
     enddo
     enddo

  end if

  ! Check for conservation

  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)
     if ((urbn_o(io,jo)) > 100.000001_r8) then
        write (6,*) 'MKURBAN error: urban = ',urbn_o(io,jo), &
                ' greater than 100.000001 for column, row = ',io,jo
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

  sum_fldo = 0._r8
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

  if ( .not. all_urban .and. trim(mksrf_gridtype) == 'global') then
     if ( abs(sum_fldo/sum_fldi-1._r8) > relerr ) then
        write (6,*) 'MKURBAN error: input field not conserved'
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

  gurbn_i = 0._r8
  garea_i = 0._r8

  do ji = 1, nlat_i
  do ii = 1, nlon_i
     garea_i = garea_i + tdomain%area(ii,ji)
     gurbn_i = gurbn_i + urbn_i(ii,ji)*(tdomain%area(ii,ji)/100._r8) * &
                         tdomain%frac(ii,ji)
  end do
  end do

  ! Output grid

  gurbn_o = 0._r8
  garea_o = 0._r8

  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)
     garea_o = garea_o + ldomain%area(io,jo)
     gurbn_o = gurbn_o + urbn_o(io,jo)*(ldomain%area(io,jo)/100._r8) * &
                         ldomain%frac(io,jo)
  end do
  end do

  ! Diagnostic output

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',k=1,70)
  write (ndiag,*) 'Urban Output'
  write (ndiag,'(1x,70a1)') ('=',k=1,70)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/ &
             1x,'                 10**6 km**2      10**6 km**2   ')
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,*)
! write (ndiag,2002) gurbn_i*1.e-06,gurbn_o*1.e-06
  write (ndiag,2003) gurbn_i*1.e-06,gurbn_o*1.e-06
  write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2002 format (1x,'urban       ',f14.3,f17.3)
2003 format (1x,'urban       ',f14.3,f22.8)
2004 format (1x,'all surface ',f14.3,f17.3)

  if (lsmlat > 1) then
     k = lsmlat/2
     write (ndiag,*)
     write (ndiag,*) 'For reference the area on the output grid of a cell near the equator is: '
     write (ndiag,'(f10.3,a14)')ldomain%area(1,k)*1.e-06,' x 10**6 km**2'
     write (ndiag,*)
  endif
  call shr_sys_flush(ndiag)

  write (6,*) 'Successfully made %urban'
  write (6,*)
  call shr_sys_flush(6)

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (urbn_i)
  deallocate (mask_i,mask_o)
  deallocate ( fld_i, fld_o)

end subroutine mkurban

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkurbanpar
!
! !INTERFACE:
subroutine mkurbanpar(lsmlon, lsmlat, fname, ndiag, ncido)
!
! !DESCRIPTION:
! Make Urban Parameter data
!
! !USES:
  use domainMod   , only : domain_type,domain_clean,domain_setptrs
  use creategridMod, only : read_domain
  use mkvarpar	  , only : nlevurb, numsolar, numrad
  use mkvarsur    , only : ldomain
  use mkvarctl    
  use areaMod     , only : areaini,areaave,gridmap_clean,gridmap_setptrs
  use ncdio
!
! !ARGUMENTS:
  implicit none
  integer , intent(in) :: lsmlon, lsmlat          ! clm grid resolution
  character(len=256), intent(in) :: fname         ! input dataset file name
  integer , intent(in) :: ndiag                   ! unit number for diag out
  integer , intent(in) :: ncido                   ! output netcdf file id
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
  integer  :: nlon_i                          ! input grid : lon points
  integer  :: nlat_i                          ! input grid : lat points

  type(domain_type)     :: tdomain            ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap

  real(r8), allocatable :: canyon_hwr_o(:,:)  ! canyon height to width ratio out
  real(r8), allocatable :: canyon_hwr_i(:,:)  ! canyon_height to width ratio in
  real(r8), allocatable :: em_improad_o(:,:)  ! emissivity of impervious road out
  real(r8), allocatable :: em_improad_i(:,:)  ! emissivity of impervious road in
  real(r8), allocatable :: em_perroad_o(:,:)  ! emissivity of pervious road out
  real(r8), allocatable :: em_perroad_i(:,:)  ! emissivity of pervious road in
  real(r8), allocatable :: em_roof_o(:,:)     ! emissivity of roof out
  real(r8), allocatable :: em_roof_i(:,:)     ! emissivity of roof in
  real(r8), allocatable :: em_wall_o(:,:)     ! emissivity of wall out
  real(r8), allocatable :: em_wall_i(:,:)     ! emissivity of wall in
  real(r8), allocatable :: ht_roof_o(:,:)     ! height of roof out
  real(r8), allocatable :: ht_roof_i(:,:)     ! height of roof in
  real(r8), allocatable :: thick_roof_o(:,:)  ! thickness of roof out
  real(r8), allocatable :: thick_roof_i(:,:)  ! thickness of roof in
  real(r8), allocatable :: thick_wall_o(:,:)  ! thickness of wall out
  real(r8), allocatable :: thick_wall_i(:,:)  ! thickness of wall in
  real(r8), allocatable :: t_building_max_o(:,:)  ! maximum interior building temperature out
  real(r8), allocatable :: t_building_max_i(:,:)  ! maximum interior building temperature in
  real(r8), allocatable :: t_building_min_o(:,:)  ! minimum interior building temperature out
  real(r8), allocatable :: t_building_min_i(:,:)  ! minimum interior building temperature in
  real(r8), allocatable :: wind_hgt_canyon_o(:,:) ! height of wind in canyon out
  real(r8), allocatable :: wind_hgt_canyon_i(:,:) ! height of wind in canyon in
  real(r8), allocatable :: wtlunit_roof_o(:,:)    ! fraction of roof out
  real(r8), allocatable :: wtlunit_roof_i(:,:)    ! fraction of roof in
  real(r8), allocatable :: wtroad_perv_o(:,:)     ! fraction of pervious road out
  real(r8), allocatable :: wtroad_perv_i(:,:)     ! fraction of pervious road in
  real(r8), allocatable :: alb_improad_o(:,:,:,:)  ! albedo of impervious road out (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_improad_i(:,:,:,:)  ! albedo of impervious road in (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_perroad_o(:,:,:,:)  ! albedo of pervious road out (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_perroad_i(:,:,:,:)  ! albedo of pervious road in (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_roof_o(:,:,:,:)  ! albedo of roof out (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_roof_i(:,:,:,:)  ! albedo of roof in (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_wall_o(:,:,:,:)  ! albedo of wall out (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_wall_i(:,:,:,:)  ! albedo of wall in (lon,lat,numrad,numsolar)
  real(r8), allocatable :: tk_roof_o(:,:,:)     ! thermal conductivity of roof out (lon,lat,nlevurb)
  real(r8), allocatable :: tk_roof_i(:,:,:)     ! thermal conductivity of roof in (lon,lat,nlevurb)
  real(r8), allocatable :: tk_wall_o(:,:,:)     ! thermal conductivity of wall out (lon,lat,nlevurb)
  real(r8), allocatable :: tk_wall_i(:,:,:)     ! thermal conductivity of wall in (lon,lat,nlevurb)
  real(r8), allocatable :: tk_improad_o(:,:,:)   ! thermal conductivity of impervious road out (lon,lat,nlevurb)
  real(r8), allocatable :: tk_improad_i(:,:,:)   ! thermal conductivity of impervious road in (lon,lat,nlevurb)
  real(r8), allocatable :: cv_roof_o(:,:,:)     ! volumetric heat capacity of roof out (lon,lat,nlevurb)
  real(r8), allocatable :: cv_roof_i(:,:,:)     ! volumetric heat capacity of roof in (lon,lat,nlevurb)
  real(r8), allocatable :: cv_wall_o(:,:,:)     ! volumetric heat capacity of wall out (lon,lat,nlevurb)
  real(r8), allocatable :: cv_wall_i(:,:,:)     ! volumetric heat capacity of wall in (lon,lat,nlevurb)
  real(r8), allocatable :: cv_improad_o(:,:,:)   ! volumetric heat capacity of impervious road out (lon,lat,nlevurb)
  real(r8), allocatable :: cv_improad_i(:,:,:)   ! volumetric heat capacity of impervious road in (lon,lat,nlevurb)
  integer,  allocatable :: nlev_improad_o(:,:)! number of impervious road layers out
  real(r8), allocatable :: mask_i(:,:)        ! input grid: mask (0, 1)
  real(r8), allocatable :: mask_o(:,:)        ! output grid: mask (0, 1)
  integer , allocatable :: num_good_o(:,:)    ! number of non-missing inputs within output gridcell
  integer , allocatable :: num_miss_o(:,:)    ! number of missing inputs within output gridcell
  integer  :: nlevurb_i                       ! input  grid: number of urban vertical levels
  integer  :: numsolar_i                      ! input  grid: number of solar type (DIR/DIF)
  integer  :: numrad_i                        ! input  grid: number of solar bands (VIS/NIR)
  integer          :: mxovr                 ! total num of overlapping cells
  integer ,pointer :: novr(:,:)             ! number of overlapping input cells
  integer ,pointer :: iovr(:,:,:)           ! lon index of overlap input cell
  integer ,pointer :: jovr(:,:,:)           ! lat index of overlap input cell
  real(r8),pointer :: wovr(:,:,:)           ! weight    of overlap input cell

  integer  :: ii,ji                           ! indices
  integer  :: io,jo                           ! indices
  integer  :: k,l,n,m                         ! indices
  integer  :: nsolar,nrad,nurb                ! indices
  integer  :: ncidi,dimid,varid               ! input netCDF id's
  integer  :: numlev                          ! number of valid impervious road layers
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001                ! max error: sum overlap wts ne 1
  real(r8) :: wt                              ! temporary weight of overlap input cell
  character(len=256) :: name                  ! name of attribute
  character(len=256) :: unit                  ! units of attribute
  character(len= 32) :: subname = 'mkurbanpar'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make Urban Parameters .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call read_domain(tdomain,fname,readmask=.true.)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(fname, 0, ncidi), subname)

  call check_ret(nf_inq_dimid(ncidi, 'nlevurb', dimid), subname)
  call check_ret(nf_inq_dimlen(ncidi, dimid, nlevurb_i), subname)
  if (nlevurb_i /= nlevurb) then
     write(6,*)'MKURBANPAR: parameter nlevurb= ',nlevurb, &
          'does not equal input dataset nlevurb= ',nlevurb_i
     stop
  endif
  call check_ret(nf_inq_dimid(ncidi, 'numsolar', dimid), subname)
  call check_ret(nf_inq_dimlen(ncidi, dimid, numsolar_i), subname)
  if (numsolar_i /= numsolar) then
     write(6,*)'MKURBANPAR: parameter numsolar= ',numsolar, &
          'does not equal input dataset numsolar= ',numsolar_i
     stop
  endif
  call check_ret(nf_inq_dimid(ncidi, 'numrad', dimid), subname)
  call check_ret(nf_inq_dimlen(ncidi, dimid, numrad_i), subname)
  if (numrad_i /= numrad) then
     write(6,*)'MKURBANPAR: parameter numrad= ',numrad, &
          'does not equal input dataset numrad= ',numrad_i
     stop
  endif

  ! Allocation

  allocate(canyon_hwr_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(em_improad_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(em_perroad_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(em_roof_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(em_wall_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(ht_roof_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(thick_roof_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(thick_wall_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(t_building_max_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(t_building_min_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(wind_hgt_canyon_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(wtlunit_roof_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(wtroad_perv_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_improad_i(nlon_i,nlat_i,numrad_i,numsolar_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_perroad_i(nlon_i,nlat_i,numrad_i,numsolar_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_roof_i(nlon_i,nlat_i,numrad_i,numsolar_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_wall_i(nlon_i,nlat_i,numrad_i,numsolar_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if

  allocate(canyon_hwr_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(em_improad_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(em_perroad_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(em_roof_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(em_wall_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(ht_roof_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(thick_roof_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(thick_wall_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(t_building_max_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(t_building_min_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(wind_hgt_canyon_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(wtlunit_roof_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(wtroad_perv_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_improad_o(lsmlon,lsmlat,numrad,numsolar), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_perroad_o(lsmlon,lsmlat,numrad,numsolar), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_roof_o(lsmlon,lsmlat,numrad,numsolar), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_wall_o(lsmlon,lsmlat,numrad,numsolar), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(mask_i(nlon_i,nlat_i),mask_o(lsmlon,lsmlat), stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if

  ! Compute local fields _o
  ! Area average and then deallocate input data

  mask_i = 1.0_r8
  mask_o = 1.0_r8
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  mask_i = float(tdomain%mask(:,:))

  call areaave(mask_i,mask_o,tgridmap)

  call gridmap_clean(tgridmap)
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  call check_ret(nf_inq_varid (ncidi, 'CANYON_HWR', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, canyon_hwr_i), subname)
  call areaave(canyon_hwr_i,canyon_hwr_o,tgridmap)
  deallocate (canyon_hwr_i)

  call check_ret(nf_inq_varid (ncidi, 'EM_IMPROAD', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, em_improad_i), subname)
  call areaave(em_improad_i,em_improad_o,tgridmap)
  deallocate (em_improad_i)

  call check_ret(nf_inq_varid (ncidi, 'EM_PERROAD', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, em_perroad_i), subname)
  call areaave(em_perroad_i,em_perroad_o,tgridmap)
  deallocate (em_perroad_i)

  call check_ret(nf_inq_varid (ncidi, 'EM_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, em_roof_i), subname)
  call areaave(em_roof_i,em_roof_o,tgridmap)
  deallocate (em_roof_i)

  call check_ret(nf_inq_varid (ncidi, 'EM_WALL', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, em_wall_i), subname)
  call areaave(em_wall_i,em_wall_o,tgridmap)
  deallocate (em_wall_i)

  call check_ret(nf_inq_varid (ncidi, 'HT_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, ht_roof_i), subname)
  call areaave(ht_roof_i,ht_roof_o,tgridmap)
  deallocate (ht_roof_i)

  call check_ret(nf_inq_varid (ncidi, 'THICK_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, thick_roof_i), subname)
  call areaave(thick_roof_i,thick_roof_o,tgridmap)
  deallocate (thick_roof_i)

  call check_ret(nf_inq_varid (ncidi, 'THICK_WALL', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, thick_wall_i), subname)
  call areaave(thick_wall_i,thick_wall_o,tgridmap)
  deallocate (thick_wall_i)

  call check_ret(nf_inq_varid (ncidi, 'T_BUILDING_MAX', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, t_building_max_i), subname)
  call areaave(t_building_max_i,t_building_max_o,tgridmap)
  deallocate (t_building_max_i)

  call check_ret(nf_inq_varid (ncidi, 'T_BUILDING_MIN', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, t_building_min_i), subname)
  call areaave(t_building_min_i,t_building_min_o,tgridmap)
  deallocate (t_building_min_i)

  call check_ret(nf_inq_varid (ncidi, 'WIND_HGT_CANYON', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, wind_hgt_canyon_i), subname)
  call areaave(wind_hgt_canyon_i,wind_hgt_canyon_o,tgridmap)
  deallocate (wind_hgt_canyon_i)

  call check_ret(nf_inq_varid (ncidi, 'WTLUNIT_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, wtlunit_roof_i), subname)
  call areaave(wtlunit_roof_i,wtlunit_roof_o,tgridmap)
  deallocate (wtlunit_roof_i)

  call check_ret(nf_inq_varid (ncidi, 'WTROAD_PERV', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, wtroad_perv_i), subname)
  call areaave(wtroad_perv_i,wtroad_perv_o,tgridmap)
  deallocate (wtroad_perv_i)

  call check_ret(nf_inq_varid (ncidi, 'ALB_IMPROAD', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, alb_improad_i), subname)
  call areaave(alb_improad_i,alb_improad_o,tgridmap)
  deallocate (alb_improad_i)

  call check_ret(nf_inq_varid (ncidi, 'ALB_PERROAD', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, alb_perroad_i), subname)
  call areaave(alb_perroad_i,alb_perroad_o,tgridmap)
  deallocate (alb_perroad_i)

  call check_ret(nf_inq_varid (ncidi, 'ALB_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, alb_roof_i), subname)
  call areaave(alb_roof_i,alb_roof_o,tgridmap)
  deallocate (alb_roof_i)

  call check_ret(nf_inq_varid (ncidi, 'ALB_WALL', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, alb_wall_i), subname)
  call areaave(alb_wall_i,alb_wall_o,tgridmap)
  deallocate (alb_wall_i)

  ! Now write output data to the file and then deallocate
  call check_ret(nf_inq_varid(ncido, 'CANYON_HWR', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, canyon_hwr_o), subname)
  deallocate (canyon_hwr_o)

  call check_ret(nf_inq_varid(ncido, 'EM_IMPROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, em_improad_o), subname)
  deallocate (em_improad_o)

  call check_ret(nf_inq_varid(ncido, 'EM_PERROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, em_perroad_o), subname)
  deallocate (em_perroad_o)

  call check_ret(nf_inq_varid(ncido, 'EM_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, em_roof_o), subname)
  deallocate (em_roof_o)

  call check_ret(nf_inq_varid(ncido, 'EM_WALL', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, em_wall_o), subname)
  deallocate (em_wall_o)

  call check_ret(nf_inq_varid(ncido, 'HT_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, ht_roof_o), subname)
  deallocate (ht_roof_o)

  call check_ret(nf_inq_varid(ncido, 'THICK_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, thick_roof_o), subname)
  deallocate (thick_roof_o)

  call check_ret(nf_inq_varid(ncido, 'THICK_WALL', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, thick_wall_o), subname)
  deallocate (thick_wall_o)

  call check_ret(nf_inq_varid(ncido, 'T_BUILDING_MAX', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, t_building_max_o), subname)
  deallocate (t_building_max_o)

  call check_ret(nf_inq_varid(ncido, 'T_BUILDING_MIN', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, t_building_min_o), subname)
  deallocate (t_building_min_o)

  call check_ret(nf_inq_varid(ncido, 'WIND_HGT_CANYON', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, wind_hgt_canyon_o), subname)
  deallocate (wind_hgt_canyon_o)

  call check_ret(nf_inq_varid(ncido, 'WTLUNIT_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, wtlunit_roof_o), subname)
  deallocate (wtlunit_roof_o)

  call check_ret(nf_inq_varid(ncido, 'WTROAD_PERV', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, wtroad_perv_o), subname)
  deallocate (wtroad_perv_o)

  call check_ret(nf_inq_varid(ncido, 'ALB_IMPROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, alb_improad_o), subname)
  deallocate (alb_improad_o)

  call check_ret(nf_inq_varid(ncido, 'ALB_PERROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, alb_perroad_o), subname)
  deallocate (alb_perroad_o)

  call check_ret(nf_inq_varid(ncido, 'ALB_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, alb_roof_o), subname)
  deallocate (alb_roof_o)

  call check_ret(nf_inq_varid(ncido, 'ALB_WALL', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, alb_wall_o), subname)
  deallocate (alb_wall_o)

  !
  ! 3D nlevurb fields
  !
  ! First allocate data
  allocate(cv_improad_i(nlon_i,nlat_i,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if

  allocate(tk_roof_i(nlon_i,nlat_i,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(tk_wall_i(nlon_i,nlat_i,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(tk_improad_i(nlon_i,nlat_i,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(cv_roof_i(nlon_i,nlat_i,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(cv_wall_i(nlon_i,nlat_i,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if

  allocate(tk_roof_o(lsmlon,lsmlat,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(tk_wall_o(lsmlon,lsmlat,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(tk_improad_o(lsmlon,lsmlat,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(cv_roof_o(lsmlon,lsmlat,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(cv_wall_o(lsmlon,lsmlat,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(cv_improad_o(lsmlon,lsmlat,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if

  ! Next read in input data
  call check_ret(nf_inq_varid (ncidi, 'TK_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, tk_roof_i), subname)

  call check_ret(nf_inq_varid (ncidi, 'TK_WALL', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, tk_wall_i), subname)

  call check_ret(nf_inq_varid (ncidi, 'CV_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, cv_roof_i), subname)

  call check_ret(nf_inq_varid (ncidi, 'CV_WALL', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, cv_wall_i), subname)

  ! Do the areaaveraging and then deallocate data
  call areaave(tk_roof_i,tk_roof_o,tgridmap)
  deallocate (tk_roof_i)

  call areaave(tk_wall_i,tk_wall_o,tgridmap)
  deallocate (tk_wall_i)

  call areaave(cv_roof_i,cv_roof_o,tgridmap)
  deallocate (cv_roof_i)

  call areaave(cv_wall_i,cv_wall_o,tgridmap)
  deallocate (cv_wall_i)

  call check_ret(nf_inq_varid(ncido, 'TK_WALL', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, tk_wall_o), subname)
  deallocate (tk_wall_o)

  call check_ret(nf_inq_varid(ncido, 'TK_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, tk_roof_o), subname)
  deallocate (tk_roof_o)

  call check_ret(nf_inq_varid(ncido, 'CV_WALL', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, cv_wall_o), subname)
  deallocate (cv_wall_o)

  call check_ret(nf_inq_varid(ncido, 'CV_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, cv_roof_o), subname)
  deallocate (cv_roof_o)

  ! Get fields from input file
  call check_ret(nf_inq_varid (ncidi, 'CV_IMPROAD', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, cv_improad_i), subname)

  call check_ret(nf_inq_varid (ncidi, 'TK_IMPROAD', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, tk_improad_i), subname)

  ! Impervious road thermal conductivity and heat capacity need to be
  ! handled differently because of varying levels of data.

  allocate(num_good_o(lsmlon,lsmlat),num_miss_o(lsmlon,lsmlat), stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if

  do nurb = 1,nlevurb
      write(6,*)'nlevurb: ',nurb
      mask_i = 1.0_r8
      mask_o = 1.0_r8
      call gridmap_clean(tgridmap)
      call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)
      call gridmap_setptrs(tgridmap,mx_ovr=mxovr,n_ovr=novr,i_ovr=iovr,j_ovr=jovr,w_ovr=wovr)
      ! Create mask for input data from missing values
      do ji = 1, nlat_i
      do ii = 1, nlon_i
        if (tk_improad_i(ii,ji,nurb) .eq. -999.) then
           mask_i(ii,ji) = 0.
        else
           mask_i(ii,ji) = 1.
        end if
      end do
      end do
      num_good_o = 0
      num_miss_o = 0
      do jo = 1, ldomain%nj
         do io = 1, ldomain%numlon(jo)
            mask_o(io,jo) = 0.
            do n = 1, novr(io,jo) !overlap cell index
               ii = iovr(io,jo,n) !lon index (input grid) of overlap cell
               ji = jovr(io,jo,n) !lat index (input grid) of overlap cell
               mask_o(io,jo) = mask_o(io,jo) + mask_i(ii,ji) * wovr(io,jo,n)
               if (mask_i(ii,ji) .eq. 1.) then
                 num_good_o(io,jo) = num_good_o(io,jo) + 1
               else
                 num_miss_o(io,jo) = num_miss_o(io,jo) + 1
               end if
            end do
         end do
      end do
      call gridmap_clean(tgridmap)
      call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)
      call gridmap_setptrs(tgridmap,mx_ovr=mxovr,n_ovr=novr,i_ovr=iovr,j_ovr=jovr,w_ovr=wovr)
      do jo = 1, ldomain%nj
         do io = 1, ldomain%numlon(jo)
            tk_improad_o(io,jo,nurb) = 0.
            cv_improad_o(io,jo,nurb) = 0.
            do n = 1, novr(io,jo)  !overlap cell index
               ii = iovr(io,jo,n)  !lon index (input grid) of overlap cell
               ji = jovr(io,jo,n)  !lat index (input grid) of overlap cell
               wt = wovr(io,jo,n)  !overlap weight
               if (num_good_o(io,jo) .gt. 0) then
                 tk_improad_o(io,jo,nurb) = tk_improad_o(io,jo,nurb) + tk_improad_i(ii,ji,nurb) * wt
                 cv_improad_o(io,jo,nurb) = cv_improad_o(io,jo,nurb) + cv_improad_i(ii,ji,nurb) * wt
               end if
            end do
         end do
      end do
  end do
  deallocate (num_miss_o)
  deallocate (num_good_o)
  deallocate (cv_improad_i)
  deallocate (tk_improad_i)


  ! Deallocate dynamic memory needed for regridding

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)

  ! Output fields to file
  call check_ret(nf_inq_varid(ncido, 'TK_IMPROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, tk_improad_o), subname)

  call check_ret(nf_inq_varid(ncido, 'CV_IMPROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, cv_improad_o), subname)

  allocate(nlev_improad_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  nlev_improad_o(:,:)  = 0
  do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
       numlev = 0
       do nurb = 1,nlevurb
         if (tk_improad_o(io,jo,nurb) .gt. 0. .and. cv_improad_o(io,jo,nurb) .gt. 0.) then
           numlev = numlev+1 
         end if
       end do
       nlev_improad_o(io,jo) = numlev
     end do
  end do

  call check_ret(nf_inq_varid(ncido, 'NLEV_IMPROAD', varid), subname)
  call check_ret(nf_put_var_int(ncido, varid, nlev_improad_o), subname)
  ! Deallocate dynamic memory
  deallocate (nlev_improad_o)
  deallocate (cv_improad_o)
  deallocate (tk_improad_o)
  deallocate (mask_i,mask_o)
     
  call check_ret(nf_sync(ncido), subname)

  write (6,*) 'Successfully made Urban Parameters'
  write (6,*)
  call shr_sys_flush(6)

  call check_ret(nf_close(ncidi), subname)

end subroutine mkurbanpar

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkelev
!
! !INTERFACE:
subroutine mkelev(lsmlon, lsmlat, fname1, fname2, ndiag, elev_o, ncido)
!
! !DESCRIPTION:
! Make elevation data
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
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
  character(len=256), intent(in) :: fname1        ! input dataset file name
  character(len=256), intent(in) :: fname2        ! input dataset file name
  integer , intent(in) :: ndiag                   ! unit number for diag out
  integer , intent(in) :: ncido                   ! output netcdf file id
  real(r8), intent(out) :: elev_o(lsmlon,lsmlat)  !
!
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Keith Oleson
!
!
! !LOCAL VARIABLES:
!EOP
  integer  :: nlon_i                          ! input grid : lon points
  integer  :: nlat_i                          ! input grid : lat points

  type(domain_type)     :: tdomain            ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap

  real(r8), allocatable :: elev_i(:,:)  ! canyon_height to width ratio in
  real(r8), allocatable :: mask_i(:,:)        ! input grid: mask (0, 1)
  real(r8), allocatable :: mask_o(:,:)        ! output grid: mask (0, 1)

  integer  :: ii,ji                           ! indices
  integer  :: io,jo                           ! indices
  integer  :: k,l,n,m                         ! indices
  integer  :: ncidi,dimid,varid               ! input netCDF id's
  integer  :: numlev                          ! number of valid impervious road layers
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001                ! max error: sum overlap wts ne 1
  real(r8) :: wt                              ! temporary weight of overlap input cell
  character(len=256) :: name                  ! name of attribute
  character(len=256) :: unit                  ! units of attribute
  character(len= 32) :: subname = 'mkelev'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make elevation .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call read_domain(tdomain,fname1)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(fname1, 0, ncidi), subname)

  ! Allocation

  allocate(elev_i(nlon_i,nlat_i), stat=ier)
  if (ier /= 0) then
     write(6,*)'mkelev allocation error'; call abort()
  end if

  allocate(mask_i(nlon_i,nlat_i),mask_o(lsmlon,lsmlat), stat=ier)
  if (ier /= 0) then
     write(6,*)'mkelev allocation error'; call abort()
  end if

  call check_ret(nf_inq_varid (ncidi, 'TOPO_ICE', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, elev_i), subname)
  call check_ret(nf_close(ncidi), subname)

  ! Compute local fields _o

  mask_i = 1.0_r8
  mask_o = 1.0_r8
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call read_domain(tdomain,fname2)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(fname2, 0, ncidi), subname)
  call check_ret(nf_inq_varid (ncidi, 'LANDMASK', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, mask_i), subname)
  call check_ret(nf_close(ncidi), subname)

  tdomain%mask(:,:) = mask_i

  call areaave(mask_i,mask_o,tgridmap)

  call gridmap_clean(tgridmap)
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  elev_o(:,:)  = 0.
  call areaave(elev_i,elev_o,tgridmap)

  write (6,*) 'Successfully made elevation' 
  write (6,*)
  call shr_sys_flush(6)

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (elev_i)
  deallocate (mask_i,mask_o)

end subroutine mkelev

!-----------------------------------------------------------------------

end module mkurbanparMod
