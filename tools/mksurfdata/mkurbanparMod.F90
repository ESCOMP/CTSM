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
!EOP
!-----------------------------------------------------------------------
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use areaMod     , only : gridmap_type

  implicit none

  private

  public  :: mkurbanpar

contains

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
  use fileutils   , only : getfil
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
!EOP
!
! !LOCAL VARIABLES:
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
  real(r8), allocatable :: alb_improad_sdiv_o(:,:) ! albedo of impervious road out on sdiv
  real(r8), allocatable :: alb_improad_sdiv_i(:,:) ! albedo of impervious road in on sdiv
  real(r8), allocatable :: alb_perroad_o(:,:,:,:)  ! albedo of pervious road out (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_perroad_i(:,:,:,:)  ! albedo of pervious road in (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_perroad_sdiv_o(:,:) ! albedo of pervious road out on sdiv
  real(r8), allocatable :: alb_perroad_sdiv_i(:,:) ! albedo of pervious road in on sdiv
  real(r8), allocatable :: alb_roof_o(:,:,:,:)  ! albedo of roof out (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_roof_i(:,:,:,:)  ! albedo of roof in (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_roof_sdiv_o(:,:) ! albedo of roof out on sdiv
  real(r8), allocatable :: alb_roof_sdiv_i(:,:) ! albedo of roof in on sdiv
  real(r8), allocatable :: alb_wall_o(:,:,:,:)  ! albedo of wall out (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_wall_i(:,:,:,:)  ! albedo of wall in (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_wall_sdiv_o(:,:) ! albedo of wall out on sdiv
  real(r8), allocatable :: alb_wall_sdiv_i(:,:) ! albedo of wall in on sdiv
  real(r8), allocatable :: tk_roof_o(:,:,:)     ! thermal conductivity of roof out (lon,lat,nlevurb)
  real(r8), allocatable :: tk_roof_i(:,:,:)     ! thermal conductivity of roof in (lon,lat,nlevurb)
  real(r8), allocatable :: tk_roof_lev_o(:,:)   ! thermal conductivity of roof out on lev
  real(r8), allocatable :: tk_roof_lev_i(:,:)   ! thermal conductivity of roof in on lev
  real(r8), allocatable :: tk_wall_o(:,:,:)     ! thermal conductivity of wall out (lon,lat,nlevurb)
  real(r8), allocatable :: tk_wall_i(:,:,:)     ! thermal conductivity of wall in (lon,lat,nlevurb)
  real(r8), allocatable :: tk_wall_lev_o(:,:)   ! thermal conductivity of wall out on lev
  real(r8), allocatable :: tk_wall_lev_i(:,:)   ! thermal conductivity of wall in on lev
  real(r8), allocatable :: tk_improad_o(:,:,:)   ! thermal conductivity of impervious road out (lon,lat,nlevurb)
  real(r8), allocatable :: tk_improad_i(:,:,:)   ! thermal conductivity of impervious road in (lon,lat,nlevurb)
  real(r8), allocatable :: tk_improad_lev_o(:,:) ! thermal conductivity of impervious road out on lev
  real(r8), allocatable :: tk_improad_lev_i(:,:) ! thermal conductivity of impervious road in on lev
  real(r8), allocatable :: cv_roof_o(:,:,:)     ! volumetric heat capacity of roof out (lon,lat,nlevurb)
  real(r8), allocatable :: cv_roof_i(:,:,:)     ! volumetric heat capacity of roof in (lon,lat,nlevurb)
  real(r8), allocatable :: cv_roof_lev_o(:,:)   ! volumetric heat capacity of roof out on lev
  real(r8), allocatable :: cv_roof_lev_i(:,:)   ! volumetric heat capactiy of roof in on lev
  real(r8), allocatable :: cv_wall_o(:,:,:)     ! volumetric heat capacity of wall out (lon,lat,nlevurb)
  real(r8), allocatable :: cv_wall_i(:,:,:)     ! volumetric heat capacity of wall in (lon,lat,nlevurb)
  real(r8), allocatable :: cv_wall_lev_o(:,:)   ! volumetric heat capacity of wall out on lev
  real(r8), allocatable :: cv_wall_lev_i(:,:)   ! volumetric heat capacity of wall in on lev
  real(r8), allocatable :: cv_improad_o(:,:,:)   ! volumetric heat capacity of impervious road out (lon,lat,nlevurb)
  real(r8), allocatable :: cv_improad_i(:,:,:)   ! volumetric heat capacity of impervious road in (lon,lat,nlevurb)
  real(r8), allocatable :: cv_improad_lev_o(:,:) ! volumetric heat capacity of impervious road out on lev
  real(r8), allocatable :: cv_improad_lev_i(:,:) ! volumetric heat capacity of impervious road in on lev
  integer,  allocatable :: nlev_improad_o(:,:)! number of impervious road layers out
  integer,  allocatable :: nlev_improad_i(:,:)! number of impervious road layers in
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
  character(len=256) :: locfn                 ! local dataset file name
  character(len= 32) :: subname = 'mkurbanpar'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make Urban Parameters .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call getfil (fname, locfn, 0)

  call read_domain(tdomain,locfn,readmask=.true.)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(locfn, 0, ncidi), subname)

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
  allocate(nlev_improad_i(nlon_i,nlat_i), &
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
  allocate(alb_improad_sdiv_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_perroad_i(nlon_i,nlat_i,numrad_i,numsolar_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_perroad_sdiv_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_roof_i(nlon_i,nlat_i,numrad_i,numsolar_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_roof_sdiv_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_wall_i(nlon_i,nlat_i,numrad_i,numsolar_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_wall_sdiv_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(tk_roof_i(nlon_i,nlat_i,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(tk_roof_lev_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(tk_wall_i(nlon_i,nlat_i,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(tk_wall_lev_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(tk_improad_i(nlon_i,nlat_i,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(tk_improad_lev_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(cv_roof_i(nlon_i,nlat_i,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(cv_roof_lev_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(cv_wall_i(nlon_i,nlat_i,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(cv_wall_lev_i(nlon_i,nlat_i), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(cv_improad_i(nlon_i,nlat_i,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(cv_improad_lev_i(nlon_i,nlat_i), &
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
  allocate(nlev_improad_o(lsmlon,lsmlat), &
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
  allocate(alb_improad_sdiv_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_perroad_o(lsmlon,lsmlat,numrad,numsolar), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_perroad_sdiv_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_roof_o(lsmlon,lsmlat,numrad,numsolar), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_roof_sdiv_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_wall_o(lsmlon,lsmlat,numrad,numsolar), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(alb_wall_sdiv_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(tk_roof_o(lsmlon,lsmlat,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(tk_roof_lev_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(tk_wall_o(lsmlon,lsmlat,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(tk_wall_lev_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(tk_improad_o(lsmlon,lsmlat,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(tk_improad_lev_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(cv_roof_o(lsmlon,lsmlat,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(cv_roof_lev_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(cv_wall_o(lsmlon,lsmlat,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(cv_wall_lev_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(cv_improad_o(lsmlon,lsmlat,nlevurb), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(cv_improad_lev_o(lsmlon,lsmlat), &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if

  allocate(mask_i(nlon_i,nlat_i),mask_o(lsmlon,lsmlat), stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
  allocate(num_good_o(lsmlon,lsmlat),num_miss_o(lsmlon,lsmlat), stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if

  ! Compute local fields _o

  mask_i = 1.0_r8
  mask_o = 1.0_r8
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  mask_i = float(tdomain%mask(:,:))

  call areaave(mask_i,mask_o,tgridmap)

  call gridmap_clean(tgridmap)
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  call check_ret(nf_inq_varid (ncidi, 'CANYON_HWR', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, canyon_hwr_i), subname)
  canyon_hwr_o(:,:)  = 0.
  call areaave(canyon_hwr_i,canyon_hwr_o,tgridmap)

  call check_ret(nf_inq_varid (ncidi, 'EM_IMPROAD', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, em_improad_i), subname)
  em_improad_o(:,:)  = 0.
  call areaave(em_improad_i,em_improad_o,tgridmap)

  call check_ret(nf_inq_varid (ncidi, 'EM_PERROAD', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, em_perroad_i), subname)
  em_perroad_o(:,:)  = 0.
  call areaave(em_perroad_i,em_perroad_o,tgridmap)

  call check_ret(nf_inq_varid (ncidi, 'EM_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, em_roof_i), subname)
  em_roof_o(:,:)  = 0.
  call areaave(em_roof_i,em_roof_o,tgridmap)

  call check_ret(nf_inq_varid (ncidi, 'EM_WALL', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, em_wall_i), subname)
  em_wall_o(:,:)  = 0.
  call areaave(em_wall_i,em_wall_o,tgridmap)

  call check_ret(nf_inq_varid (ncidi, 'HT_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, ht_roof_i), subname)
  ht_roof_o(:,:)  = 0.
  call areaave(ht_roof_i,ht_roof_o,tgridmap)

  call check_ret(nf_inq_varid (ncidi, 'THICK_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, thick_roof_i), subname)
  thick_roof_o(:,:)  = 0.
  call areaave(thick_roof_i,thick_roof_o,tgridmap)

  call check_ret(nf_inq_varid (ncidi, 'THICK_WALL', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, thick_wall_i), subname)
  thick_wall_o(:,:)  = 0.
  call areaave(thick_wall_i,thick_wall_o,tgridmap)

  call check_ret(nf_inq_varid (ncidi, 'T_BUILDING_MAX', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, t_building_max_i), subname)
  t_building_max_o(:,:)  = 0.
  call areaave(t_building_max_i,t_building_max_o,tgridmap)

  call check_ret(nf_inq_varid (ncidi, 'T_BUILDING_MIN', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, t_building_min_i), subname)
  t_building_min_o(:,:)  = 0.
  call areaave(t_building_min_i,t_building_min_o,tgridmap)

  call check_ret(nf_inq_varid (ncidi, 'WIND_HGT_CANYON', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, wind_hgt_canyon_i), subname)
  wind_hgt_canyon_o(:,:)  = 0.
  call areaave(wind_hgt_canyon_i,wind_hgt_canyon_o,tgridmap)

  call check_ret(nf_inq_varid (ncidi, 'WTLUNIT_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, wtlunit_roof_i), subname)
  wtlunit_roof_o(:,:)  = 0.
  call areaave(wtlunit_roof_i,wtlunit_roof_o,tgridmap)

  call check_ret(nf_inq_varid (ncidi, 'WTROAD_PERV', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, wtroad_perv_i), subname)
  wtroad_perv_o(:,:)  = 0.
  call areaave(wtroad_perv_i,wtroad_perv_o,tgridmap)

  call check_ret(nf_inq_varid (ncidi, 'ALB_IMPROAD', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, alb_improad_i), subname)
  call check_ret(nf_inq_varid (ncidi, 'ALB_PERROAD', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, alb_perroad_i), subname)
  call check_ret(nf_inq_varid (ncidi, 'ALB_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, alb_roof_i), subname)
  call check_ret(nf_inq_varid (ncidi, 'ALB_WALL', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, alb_wall_i), subname)
  alb_improad_o(:,:,:,:)  = 0.
  alb_perroad_o(:,:,:,:)  = 0.
  alb_roof_o(:,:,:,:)  = 0.
  alb_wall_o(:,:,:,:)  = 0.
  do nsolar = 1,numsolar
    do nrad = 1,numrad
      do ji = 1, nlat_i
      do ii = 1, nlon_i
        alb_improad_sdiv_i(ii,ji) = alb_improad_i(ii,ji,nrad,nsolar)
        alb_perroad_sdiv_i(ii,ji) = alb_perroad_i(ii,ji,nrad,nsolar)
        alb_roof_sdiv_i(ii,ji)    = alb_roof_i(ii,ji,nrad,nsolar)
        alb_wall_sdiv_i(ii,ji)    = alb_wall_i(ii,ji,nrad,nsolar)
      end do
      end do 
      call areaave(alb_improad_sdiv_i,alb_improad_sdiv_o,tgridmap)
      call areaave(alb_perroad_sdiv_i,alb_perroad_sdiv_o,tgridmap)
      call areaave(alb_roof_sdiv_i,alb_roof_sdiv_o,tgridmap)
      call areaave(alb_wall_sdiv_i,alb_wall_sdiv_o,tgridmap)
      do jo = 1, ldomain%nj
      do io = 1, ldomain%numlon(jo)
        alb_improad_o(io,jo,nrad,nsolar) = alb_improad_sdiv_o(io,jo)
        alb_perroad_o(io,jo,nrad,nsolar) = alb_perroad_sdiv_o(io,jo)
        alb_roof_o(io,jo,nrad,nsolar)    = alb_roof_sdiv_o(io,jo)
        alb_wall_o(io,jo,nrad,nsolar)    = alb_wall_sdiv_o(io,jo)
      end do
      end do
    end do
  end do 

  call check_ret(nf_inq_varid (ncidi, 'TK_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, tk_roof_i), subname)
  call check_ret(nf_inq_varid (ncidi, 'TK_WALL', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, tk_wall_i), subname)
  call check_ret(nf_inq_varid (ncidi, 'CV_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, cv_roof_i), subname)
  call check_ret(nf_inq_varid (ncidi, 'CV_WALL', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, cv_wall_i), subname)
  call check_ret(nf_inq_varid (ncidi, 'CV_IMPROAD', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, cv_improad_i), subname)
  tk_roof_o(:,:,:)  = 0.
  tk_wall_o(:,:,:)  = 0.
  cv_roof_o(:,:,:)  = 0.
  cv_wall_o(:,:,:)  = 0.
  do nurb = 1,nlevurb
      do ji = 1, nlat_i
      do ii = 1, nlon_i
        tk_roof_lev_i(ii,ji)    = tk_roof_i(ii,ji,nurb)
        tk_wall_lev_i(ii,ji)    = tk_wall_i(ii,ji,nurb)
        cv_roof_lev_i(ii,ji)    = cv_roof_i(ii,ji,nurb)
        cv_wall_lev_i(ii,ji)    = cv_wall_i(ii,ji,nurb)
      end do
      end do 
      call areaave(tk_roof_lev_i,tk_roof_lev_o,tgridmap)
      call areaave(tk_wall_lev_i,tk_wall_lev_o,tgridmap)
      call areaave(cv_roof_lev_i,cv_roof_lev_o,tgridmap)
      call areaave(cv_wall_lev_i,cv_wall_lev_o,tgridmap)
      do jo = 1, ldomain%nj
      do io = 1, ldomain%numlon(jo)
        tk_roof_o(io,jo,nurb)    = tk_roof_lev_o(io,jo)
        tk_wall_o(io,jo,nurb)    = tk_wall_lev_o(io,jo)
        cv_roof_o(io,jo,nurb)    = cv_roof_lev_o(io,jo)
        cv_wall_o(io,jo,nurb)    = cv_wall_lev_o(io,jo)
      end do
      end do
  end do 

! Impervious road thermal conductivity and heat capacity need to be
! handled differently because of varying levels of data.

  call check_ret(nf_inq_varid (ncidi, 'TK_IMPROAD', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, tk_improad_i), subname)
  call check_ret(nf_inq_varid (ncidi, 'CV_IMPROAD', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, cv_improad_i), subname)
  tk_improad_o(:,:,:)  = 0.
  cv_improad_o(:,:,:)  = 0.
  do nurb = 1,nlevurb
      write(6,*)'nlevurb: ',nurb
      mask_i = 1.0_r8
      mask_o = 1.0_r8
      call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)
      call gridmap_setptrs(tgridmap,mx_ovr=mxovr,n_ovr=novr,i_ovr=iovr,j_ovr=jovr,w_ovr=wovr)
      ! Create mask for input data from missing values
      do ji = 1, nlat_i
      do ii = 1, nlon_i
        tk_improad_lev_i(ii,ji)    = tk_improad_i(ii,ji,nurb)
        cv_improad_lev_i(ii,ji)    = cv_improad_i(ii,ji,nurb)
        if (tk_improad_lev_i(ii,ji) .eq. -999.) then
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
      call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)
      call gridmap_setptrs(tgridmap,mx_ovr=mxovr,n_ovr=novr,i_ovr=iovr,j_ovr=jovr,w_ovr=wovr)
      do jo = 1, ldomain%nj
         do io = 1, ldomain%numlon(jo)
            tk_improad_lev_o(io,jo) = 0.
            cv_improad_lev_o(io,jo) = 0.
            do n = 1, novr(io,jo)  !overlap cell index
               ii = iovr(io,jo,n)  !lon index (input grid) of overlap cell
               ji = jovr(io,jo,n)  !lat index (input grid) of overlap cell
               wt = wovr(io,jo,n)  !overlap weight
               if (num_good_o(io,jo) .gt. 0) then
                 tk_improad_lev_o(io,jo) = tk_improad_lev_o(io,jo) + tk_improad_lev_i(ii,ji) * wt
                 cv_improad_lev_o(io,jo) = cv_improad_lev_o(io,jo) + cv_improad_lev_i(ii,ji) * wt
               end if
            end do
            tk_improad_o(io,jo,nurb) = tk_improad_lev_o(io,jo)
            cv_improad_o(io,jo,nurb) = cv_improad_lev_o(io,jo)
         end do
      end do
  end do

  ! -----------------------------------------------------------------
  ! Output model resolution Urban Parameter data
  ! -----------------------------------------------------------------

  call check_ret(nf_inq_varid(ncido, 'CANYON_HWR', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, canyon_hwr_o), subname)

  call check_ret(nf_inq_varid(ncido, 'EM_IMPROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, em_improad_o), subname)

  call check_ret(nf_inq_varid(ncido, 'EM_PERROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, em_perroad_o), subname)

  call check_ret(nf_inq_varid(ncido, 'EM_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, em_roof_o), subname)

  call check_ret(nf_inq_varid(ncido, 'EM_WALL', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, em_wall_o), subname)

  call check_ret(nf_inq_varid(ncido, 'HT_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, ht_roof_o), subname)

  call check_ret(nf_inq_varid(ncido, 'THICK_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, thick_roof_o), subname)

  call check_ret(nf_inq_varid(ncido, 'THICK_WALL', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, thick_wall_o), subname)

  call check_ret(nf_inq_varid(ncido, 'T_BUILDING_MAX', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, t_building_max_o), subname)

  call check_ret(nf_inq_varid(ncido, 'T_BUILDING_MIN', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, t_building_min_o), subname)

  call check_ret(nf_inq_varid(ncido, 'WIND_HGT_CANYON', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, wind_hgt_canyon_o), subname)

  call check_ret(nf_inq_varid(ncido, 'WTLUNIT_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, wtlunit_roof_o), subname)

  call check_ret(nf_inq_varid(ncido, 'WTROAD_PERV', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, wtroad_perv_o), subname)

  call check_ret(nf_inq_varid(ncido, 'ALB_IMPROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, alb_improad_o), subname)

  call check_ret(nf_inq_varid(ncido, 'ALB_PERROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, alb_perroad_o), subname)

  call check_ret(nf_inq_varid(ncido, 'ALB_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, alb_roof_o), subname)

  call check_ret(nf_inq_varid(ncido, 'ALB_WALL', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, alb_wall_o), subname)

  call check_ret(nf_inq_varid(ncido, 'TK_WALL', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, tk_wall_o), subname)

  call check_ret(nf_inq_varid(ncido, 'TK_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, tk_roof_o), subname)

  call check_ret(nf_inq_varid(ncido, 'TK_IMPROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, tk_improad_o), subname)

  call check_ret(nf_inq_varid(ncido, 'CV_WALL', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, cv_wall_o), subname)

  call check_ret(nf_inq_varid(ncido, 'CV_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, cv_roof_o), subname)

  call check_ret(nf_inq_varid(ncido, 'CV_IMPROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, cv_improad_o), subname)

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
     
  call check_ret(nf_sync(ncido), subname)

  write (6,*) 'Successfully made Urban Parameters'
  write (6,*)
  call shr_sys_flush(6)

  call check_ret(nf_close(ncidi), subname)

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (canyon_hwr_o)
  deallocate (canyon_hwr_i)
  deallocate (em_improad_o)
  deallocate (em_improad_i)
  deallocate (em_perroad_o)
  deallocate (em_perroad_i)
  deallocate (em_roof_o)
  deallocate (em_roof_i)
  deallocate (em_wall_o)
  deallocate (em_wall_i)
  deallocate (ht_roof_o)
  deallocate (ht_roof_i)
  deallocate (nlev_improad_o)
  deallocate (nlev_improad_i)
  deallocate (thick_roof_o)
  deallocate (thick_roof_i)
  deallocate (thick_wall_o)
  deallocate (thick_wall_i)
  deallocate (t_building_max_o)
  deallocate (t_building_max_i)
  deallocate (t_building_min_o)
  deallocate (t_building_min_i)
  deallocate (wind_hgt_canyon_o)
  deallocate (wind_hgt_canyon_i)
  deallocate (wtlunit_roof_o)
  deallocate (wtlunit_roof_i)
  deallocate (wtroad_perv_o)
  deallocate (wtroad_perv_i)
  deallocate (alb_improad_o)
  deallocate (alb_improad_i)
  deallocate (alb_improad_sdiv_o)
  deallocate (alb_improad_sdiv_i)
  deallocate (alb_perroad_o)
  deallocate (alb_perroad_i)
  deallocate (alb_perroad_sdiv_o)
  deallocate (alb_perroad_sdiv_i)
  deallocate (alb_roof_o)
  deallocate (alb_roof_i)
  deallocate (alb_roof_sdiv_o)
  deallocate (alb_roof_sdiv_i)
  deallocate (alb_wall_o)
  deallocate (alb_wall_i)
  deallocate (alb_wall_sdiv_o)
  deallocate (alb_wall_sdiv_i)
  deallocate (tk_roof_o)
  deallocate (tk_roof_i)
  deallocate (tk_roof_lev_o)
  deallocate (tk_roof_lev_i)
  deallocate (tk_wall_o)
  deallocate (tk_wall_i)
  deallocate (tk_wall_lev_o)
  deallocate (tk_wall_lev_i)
  deallocate (tk_improad_o)
  deallocate (tk_improad_i)
  deallocate (tk_improad_lev_o)
  deallocate (tk_improad_lev_i)
  deallocate (cv_roof_o)
  deallocate (cv_roof_i)
  deallocate (cv_roof_lev_o)
  deallocate (cv_roof_lev_i)
  deallocate (cv_wall_o)
  deallocate (cv_wall_i)
  deallocate (cv_wall_lev_o)
  deallocate (cv_wall_lev_i)
  deallocate (cv_improad_o)
  deallocate (cv_improad_i)
  deallocate (cv_improad_lev_o)
  deallocate (cv_improad_lev_i)
  deallocate (mask_i,mask_o)

end subroutine mkurbanpar

end module mkurbanparMod
