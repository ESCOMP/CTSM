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
  use fileutils   , only : getfil
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
  character(len=256) :: locfn                 ! local dataset file name
  character(len= 32) :: subname = 'mkelev'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make elevation .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call getfil (fname1, locfn, 0)

  call read_domain(tdomain,locfn)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(locfn, 0, ncidi), subname)

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

  call getfil (fname2, locfn, 0)

  call read_domain(tdomain,locfn)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(locfn, 0, ncidi), subname)
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
