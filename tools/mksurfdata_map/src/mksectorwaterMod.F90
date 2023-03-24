module mksectorwaterMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mksectorwater
!
! !DESCRIPTION:
! Make sectoral withdrawal/consumption data
!
! !REVISION HISTORY:
! Author: Ioan Sabin Taranu
!
!EOP
!-----------------------------------------------------------------------
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use mkdomainMod , only : domain_checksame
  use mkvarctl   
  

  implicit none

  private

  public  :: mksectorwater

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksectorwater
!
! !INTERFACE:
subroutine mksectorwater(ldomain, mapfname, datfname, ndiag, ncido)
!
! !DESCRIPTION:
! Make sectoral water withdrawal/consumption data
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkvarpar    , only : re
  use mkvarctl    
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname     ! input mapping file name
  character(len=*)  , intent(in) :: datfname     ! input data file name
  integer           , intent(in) :: ndiag        ! unit number for diag out
  integer           , intent(in) :: ncido        ! output netcdf file id

!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)    :: tgridmap
  type(domain_type)     :: tdomain              ! local domain
  real(r8), allocatable :: mdom_withd_o(:)      ! monthly withdrawal domestic output
  real(r8), allocatable :: mliv_withd_o(:)      ! monthly withdrawal livestock output
  real(r8), allocatable :: melec_withd_o(:)     ! monthly withdrawal thermoelectric output
  real(r8), allocatable :: mmfc_withd_o(:)      ! monthly withdrawal manufacturing output 
  real(r8), allocatable :: mmin_withd_o(:)      ! monthly withdrawal mining output

  real(r8), allocatable :: mdom_withd_i(:)      ! monthly withdrawal domestic input
  real(r8), allocatable :: mliv_withd_i(:)      ! monthly withdrawal livestock input
  real(r8), allocatable :: melec_withd_i(:)     ! monthly withdrawal thermoelectric input
  real(r8), allocatable :: mmfc_withd_i(:)      ! monthly withdrawal manufacturing input 
  real(r8), allocatable :: mmin_withd_i(:)      ! monthly withdrawal mining input

  real(r8), allocatable :: mdom_cons_o(:)       ! monthly consumption domestic output
  real(r8), allocatable :: mliv_cons_o(:)       ! monthly consumption livestock output
  real(r8), allocatable :: melec_cons_o(:)      ! monthly consumption thermoelectric output
  real(r8), allocatable :: mmfc_cons_o(:)       ! monthly consumption manufacturing output 
  real(r8), allocatable :: mmin_cons_o(:)       ! monthly consumption mining output

  real(r8), allocatable :: mdom_cons_i(:)       ! monthly consumption domestic input
  real(r8), allocatable :: mliv_cons_i(:)       ! monthly consumption livestock input
  real(r8), allocatable :: melec_cons_i(:)      ! monthly consumption thermoelectric input
  real(r8), allocatable :: mmfc_cons_i(:)       ! monthly consumption manufacturing input 
  real(r8), allocatable :: mmin_cons_i(:)       ! monthly consumption mining input

  real(r8) :: gdom_withd_o              ! output grid
  real(r8) :: gliv_withd_o              ! output grid 
  real(r8) :: gelec_withd_o             ! output grid 
  real(r8) :: gmfc_withd_o              ! output grid 
  real(r8) :: gmin_withd_o              ! output grid

  real(r8) :: gdom_cons_o               ! output grid
  real(r8) :: gliv_cons_o               ! output grid 
  real(r8) :: gelec_cons_o              ! output grid 
  real(r8) :: gmfc_cons_o               ! output grid 
  real(r8) :: gmin_cons_o               ! output grid

  real(r8) :: gdom_withd_i              ! output grid
  real(r8) :: gliv_withd_i              ! output grid 
  real(r8) :: gelec_withd_i             ! output grid 
  real(r8) :: gmfc_withd_i              ! output grid 
  real(r8) :: gmin_withd_i              ! output grid

  real(r8) :: gdom_cons_i               ! output grid
  real(r8) :: gliv_cons_i               ! output grid 
  real(r8) :: gelec_cons_i              ! output grid 
  real(r8) :: gmfc_cons_i               ! output grid 
  real(r8) :: gmin_cons_i               ! output grid


  real(r8), allocatable :: frac_dst(:)      ! output fractions: same as frac_dst
  real(r8) :: garea_i                       ! input  grid: global area
  real(r8) :: garea_o                       ! output grid: global area
  integer  :: mwts                          ! number of weights
  integer  :: ni,no,ns_i,ns_o               ! indices
  integer  :: k,l,n,m                       ! indices
  integer  :: ncidi,dimid,varid             ! input netCDF id's
  integer  :: ndimsi,ndimso                 ! netCDF dimension sizes 
  integer  :: dimids(4)                     ! netCDF dimension ids
  integer  :: bego(4),leno(4)               ! netCDF bounds
  integer  :: begi(4),leni(4)               ! netCDF bounds 
  integer  :: ntim                          ! number of input time samples
  real(r8), allocatable :: mask_r8(:)  ! float of tdomain%mask
  integer  :: ier                           ! error status
  real(r8) :: relerr = 0.00001              ! max error: sum overlap wts ne 1
  character(len=256) :: name                ! name of attribute
  character(len=256) :: unit                ! units of attribute
  character(len= 32) :: subname = 'mksectorwater'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make sectoral water withdrawal/consumption data .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  ns_o = ldomain%ns
  PRINT *, "ns_o:", ns_o

  call domain_read(tdomain,datfname)
  ns_i = tdomain%ns
  PRINT *, "ns_i:", ns_i

  write (6,*) 'Open sectorwater file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncidi), subname)
  call check_ret(nf_inq_dimid(ncidi, 'time', dimid), subname)
  call check_ret(nf_inq_dimlen(ncidi, dimid, ntim), subname)

  
  if (ntim /= 12) then
     write(6,*)'MKsectorwater: must have 12 time samples on input data'
     call abort()
  endif

  ! NOTE - close data set at bottom of routine

  ! Dynamic allocation of variables

  allocate(mdom_withd_i(ns_i),  &
           mliv_withd_i(ns_i),  &
           melec_withd_i(ns_i), &
           mmfc_withd_i(ns_i),  &
           mmin_withd_i(ns_i),  &
           mdom_withd_o(ns_o),  &
           mliv_withd_o(ns_o),  &
           melec_withd_o(ns_o), &
           mmfc_withd_o(ns_o),  &
           mmin_withd_o(ns_o),  &
           mdom_cons_i(ns_i),  &
           mliv_cons_i(ns_i),  &
           melec_cons_i(ns_i), &
           mmfc_cons_i(ns_i),  &
           mmin_cons_i(ns_i),  &
           mdom_cons_o(ns_o),  &
           mliv_cons_o(ns_o),  &
           melec_cons_o(ns_o), &
           mmfc_cons_o(ns_o),  &
           mmin_cons_o(ns_o),  &
           frac_dst(ns_o),     &
           stat=ier)
  if (ier /= 0) then
     write(6,*)'mksectorwater allocation error'; call abort()
  end if

  ! Determine mapping weights and map

  call gridmap_mapread(tgridmap, mapfname)

  ! Error checks for domain and map consistencies

  allocate(mask_r8(tdomain%ns), stat=ier)

  call domain_checksame( tdomain, ldomain, tgridmap )

  ! Determine number of dimensions in input by querying withd_dom

  call check_ret(nf_inq_varid(ncidi, 'withd_dom', varid), subname)
  call check_ret(nf_inq_vardimid(ncidi, varid, dimids), subname)
  call check_ret(nf_inq_varndims(ncidi, varid, ndimsi), subname)
  if (ndimsi ==4) then
     begi(1) = 1
     begi(2) = 1
     begi(3) = 1
     leni(4) = 1
     call check_ret(nf_inq_dimlen(ncidi, dimids(1), leni(1)), subname)
     call check_ret(nf_inq_dimlen(ncidi, dimids(2), leni(2)), subname)
     call check_ret(nf_inq_dimlen(ncidi, dimids(3), leni(3)), subname)
  else if (ndimsi== 3) then
     begi(1) = 1
     begi(2) = 1
     leni(3) = 1
     call check_ret(nf_inq_dimlen(ncidi, dimids(1), leni(1)), subname)
     call check_ret(nf_inq_dimlen(ncidi, dimids(2), leni(2)), subname)
  end if
  PRINT *, 'ndimsi:  ', ndimsi
  PRINT *, 'leni(1): ', leni(1)
  PRINT *, 'leni(2): ', leni(2) 

  ! Determine number of dimensions in output by querying withd_dom

  call check_ret(nf_inq_varid(ncido, 'withd_dom', varid), subname)
  call check_ret(nf_inq_varndims(ncido, varid, ndimso), subname)
  call check_ret(nf_inq_vardimid(ncido, varid, dimids), subname)
  if (ndimso ==4) then
     bego(1) = 1
     bego(2) = 1
     bego(3) = 1
     leno(4) = 1
     call check_ret(nf_inq_dimlen(ncido, dimids(1), leno(1)), subname)
     call check_ret(nf_inq_dimlen(ncido, dimids(2), leno(2)), subname)
     call check_ret(nf_inq_dimlen(ncido, dimids(3), leno(3)), subname)
  else if (ndimso== 3) then
     bego(1) = 1
     bego(2) = 1
     leno(3) = 1
     call check_ret(nf_inq_dimlen(ncido, dimids(1), leno(1)), subname)
     call check_ret(nf_inq_dimlen(ncido, dimids(2), leno(2)), subname)
  end if

  PRINT *, 'ndimso:  ', ndimso
  PRINT *, 'leno(1):  ', leno(1)
  PRINT *, 'leno(2):  ', leno(2)


  ! Loop over months 

  do m = 1, ntim

     if (ndimsi == 4) begi(4)=m
     if (ndimsi == 3) begi(3)=m
     
     call check_ret(nf_inq_varid (ncidi, 'withd_dom', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, begi(1:ndimsi), leni(1:ndimsi), &
          mdom_withd_i), subname)
     !call writenetcdffile(mdom_withd_i)
     !PRINT *, 'nf_get_vara_double is working for withd_dom'

     call check_ret(nf_inq_varid (ncidi, 'cons_dom', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, begi(1:ndimsi), leni(1:ndimsi), &
          mdom_cons_i), subname)
     !PRINT *, 'nf_get_vara_double is working for cons_dom'

     call check_ret(nf_inq_varid (ncidi, 'withd_liv', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, begi(1:ndimsi), leni(1:ndimsi), &
          mliv_withd_i), subname)
     !PRINT *, 'nf_get_vara_double is working for withd_liv'

     call check_ret(nf_inq_varid (ncidi, 'cons_liv', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, begi(1:ndimsi), leni(1:ndimsi), &
          mliv_cons_i), subname)
     !PRINT *, 'nf_get_vara_double is working for cons_liv'

     call check_ret(nf_inq_varid (ncidi, 'withd_elec', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, begi(1:ndimsi), leni(1:ndimsi), &
          melec_withd_i), subname)
     !PRINT *, 'nf_get_vara_double is working for withd_elec'

     call check_ret(nf_inq_varid (ncidi, 'cons_elec', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, begi(1:ndimsi), leni(1:ndimsi), &
          melec_cons_i), subname)
     !PRINT *, 'nf_get_vara_double is working for cons_elec'

     call check_ret(nf_inq_varid (ncidi, 'withd_mfg', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, begi(1:ndimsi), leni(1:ndimsi), &
          mmfc_withd_i), subname)
     !PRINT *, 'nf_get_vara_double is working for withd_mfc'

     call check_ret(nf_inq_varid (ncidi, 'cons_mfg', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, begi(1:ndimsi), leni(1:ndimsi), &
          mmfc_cons_i), subname)
     !PRINT *, 'nf_get_vara_double is working for cons_mfc'

     call check_ret(nf_inq_varid (ncidi, 'withd_min', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, begi(1:ndimsi), leni(1:ndimsi), &
          mmin_withd_i), subname)
     !PRINT *, 'nf_get_vara_double is working for withd_min'

     call check_ret(nf_inq_varid (ncidi, 'cons_min', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, begi(1:ndimsi), leni(1:ndimsi), &
          mmin_cons_i), subname)
     !PRINT *, 'nf_get_vara_double is working for cons_min'

     !PRINT *, "Till now everything seems to work well"

     mdom_withd_o(:)  = 0.
     mliv_withd_o(:)  = 0.
     melec_withd_o(:) = 0.
     mmfc_withd_o(:) = 0.
     mmin_withd_o(:) = 0.

     mdom_cons_o(:)  = 0.
     mliv_cons_o(:)  = 0.
     melec_cons_o(:) = 0.
     mmfc_cons_o(:) = 0.
     mmin_cons_o(:) = 0.

     ! PRINT *, "The output variables are initialized to 0"
     ! Obtain frac_dst
     !allocate(frac_dst(ldomain%ns), stat=ier)
     call gridmap_calc_frac_dst(tgridmap, tdomain%mask, frac_dst)
     mask_r8 = tdomain%mask
     !call writenetcdffile(mask_r8)
     call gridmap_check( tgridmap, mask_r8, frac_dst, subname )
        
     ! Do the mapping
     call gridmap_areaave_srcmask(tgridmap, mdom_withd_i , mdom_withd_o, nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)
     ! PRINT *,"mdom_withd"
     call gridmap_areaave_srcmask(tgridmap, mliv_withd_i(:) , mliv_withd_o(:) , nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)
     ! PRINT *,"mliv_withd"
     call gridmap_areaave_srcmask(tgridmap, melec_withd_i(:) , melec_withd_o(:) , nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)
     ! PRINT *,"melec_withd"
     call gridmap_areaave_srcmask(tgridmap, mmfc_withd_i(:) , mmfc_withd_o(:) , nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)
     ! PRINT *,"mmfc_withd"
     call gridmap_areaave_srcmask(tgridmap, mmin_withd_i(:) , mmin_withd_o(:) , nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)
     ! PRINT *,"mmin_withd"


     call gridmap_areaave_srcmask(tgridmap, mdom_cons_i(:) , mdom_cons_o(:) , nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)
     !PRINT *,"mdom_cons"
     call gridmap_areaave_srcmask(tgridmap, mliv_cons_i(:) , mliv_cons_o(:) , nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)
     !PRINT *,"mliv_cons"
     call gridmap_areaave_srcmask(tgridmap, melec_cons_i(:) , melec_cons_o(:) , nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)
     !PRINT *,"melec_cons"
     call gridmap_areaave_srcmask(tgridmap, mmfc_cons_i(:) , mmfc_cons_o(:) , nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)
     !PRINT *,"mmfc_cons"
     call gridmap_areaave_srcmask(tgridmap, mmin_cons_i(:) , mmin_cons_o(:) , nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)
     !PRINT *,"mmin_cons"


     ! PRINT *, "regridding is done correctly"
     
     ! -----------------------------------------------------------------
     ! Output model resolution sectoral Water withdrawal/consumption data
     ! -----------------------------------------------------------------

     ! Now write out all variables

     if (ndimso == 4) bego(4)=m
     if (ndimso == 3) bego(3)=m

     call check_ret(nf_inq_varid(ncido, 'withd_dom', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, bego, leno, mdom_withd_o), subname)
     ! PRINT *, "withd_dom saved to output file"

     call check_ret(nf_inq_varid(ncido, 'cons_dom', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, bego, leno, mdom_cons_o), subname)
     ! PRINT *, "cons_dom saved to output file"

     call check_ret(nf_inq_varid(ncido, 'withd_liv', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, bego, leno, mliv_withd_o), subname)
     ! PRINT *, "withd_liv saved to output file"

     call check_ret(nf_inq_varid(ncido, 'cons_liv', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, bego, leno, mliv_cons_o), subname)
     ! PRINT *, "cons_dom saved to output file"

     call check_ret(nf_inq_varid(ncido, 'withd_elec', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, bego, leno, melec_withd_o), subname)
     ! PRINT *, "withd_elec saved to output file"

     call check_ret(nf_inq_varid(ncido, 'cons_elec', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, bego, leno, melec_cons_o), subname)
     ! PRINT *, "cons_elec saved to output file"

     call check_ret(nf_inq_varid(ncido, 'withd_mfc', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, bego, leno, mmfc_withd_o), subname)
     ! PRINT *, "withd_mfc saved to output file"

     call check_ret(nf_inq_varid(ncido, 'cons_mfc', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, bego, leno, mmfc_cons_o), subname)
     ! PRINT *, "cons_mfc saved to output file"

     call check_ret(nf_inq_varid(ncido, 'withd_min', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, bego, leno, mmin_withd_o), subname)
     ! PRINT *, "withd_min saved to output file"

     call check_ret(nf_inq_varid(ncido, 'cons_min', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, bego, leno, mmin_cons_o), subname)
     ! PRINT *, "cons_min saved to output file"

     call check_ret(nf_inq_varid(ncido, 'time', varid), subname)
     call check_ret(nf_put_vara_int(ncido, varid, bego(ndimso), leno(ndimso), m), subname)
     ! PRINT *, "time var saved to output file"

     call check_ret(nf_sync(ncido), subname)


     ! -----------------------------------------------------------------
     ! Error check2
     ! Compare global areas on input and output grids
     ! -----------------------------------------------------------------

     ! Input grid global area

     garea_i    = 0.
     do ni = 1,ns_i
        garea_i = garea_i + tgridmap%area_dst(ni)
     end do

     gdom_withd_i  = 0.
     gliv_withd_i  = 0.
     gelec_withd_i = 0.
     gmfc_withd_i  = 0.
     gmin_withd_i  = 0.

     gdom_cons_i  = 0.
     gliv_cons_i  = 0.
     gelec_cons_i = 0.
     gmfc_cons_i  = 0.
     gmin_cons_i  = 0.

     do ni = 1,ns_i
         gdom_withd_i  = gdom_withd_i + mdom_withd_i(ni)*tgridmap%area_dst(ni)* &
             frac_dst(ni)*re**2
         gliv_withd_i  = gliv_withd_i + mliv_withd_i(ni)*tgridmap%area_dst(ni)* &
             frac_dst(ni)*re**2
         gelec_withd_i  = gelec_withd_i + melec_withd_i(ni)*tgridmap%area_dst(ni)* &
             frac_dst(ni)*re**2
         gmfc_withd_i  = gmfc_withd_i + mmfc_withd_i(ni)*tgridmap%area_dst(ni)* &
             frac_dst(ni)*re**2
         gmin_withd_i  = gmin_withd_i + mmin_withd_i(ni)*tgridmap%area_dst(ni)* &
             frac_dst(ni)*re**2

         gdom_cons_i  = gdom_cons_i + mdom_cons_i(ni)*tgridmap%area_dst(ni)* &
             frac_dst(ni)*re**2
         gliv_cons_i  = gliv_cons_i + mliv_cons_i(ni)*tgridmap%area_dst(ni)* &
             frac_dst(ni)*re**2
         gelec_cons_i  = gelec_cons_i + melec_cons_i(ni)*tgridmap%area_dst(ni)* &
             frac_dst(ni)*re**2
         gmfc_cons_i  = gmfc_cons_i + mmfc_cons_i(ni)*tgridmap%area_dst(ni)* &
             frac_dst(ni)*re**2
         gmin_cons_i  = gmin_cons_i + mmin_cons_i(ni)*tgridmap%area_dst(ni)* &
             frac_dst(ni)*re**2

     end do

     ! Output grid global area

     garea_o    = 0.
     do no = 1,ns_o
        garea_o = garea_o + tgridmap%area_dst(no)
     end do

     gdom_withd_o  = 0.
     gliv_withd_o  = 0.
     gelec_withd_o = 0.
     gmfc_withd_o  = 0.
     gmin_withd_o  = 0.

     gdom_cons_o  = 0.
     gliv_cons_o  = 0.
     gelec_cons_o = 0.
     gmfc_cons_o  = 0.
     gmin_cons_o  = 0.

     do no = 1,ns_o
         gdom_withd_o  = gdom_withd_o + mdom_withd_o(no)*tgridmap%area_dst(no)* &
             frac_dst(no)*re**2
         gliv_withd_o  = gliv_withd_o + mliv_withd_o(no)*tgridmap%area_dst(no)* &
             frac_dst(no)*re**2
         gelec_withd_o  = gelec_withd_o + melec_withd_o(no)*tgridmap%area_dst(no)* &
             frac_dst(no)*re**2
         gmfc_withd_o  = gmfc_withd_o + mmfc_withd_o(no)*tgridmap%area_dst(no)* &
             frac_dst(no)*re**2
         gmin_withd_o  = gmin_withd_o + mmin_withd_o(no)*tgridmap%area_dst(no)* &
             frac_dst(no)*re**2

         gdom_cons_o  = gdom_cons_o + mdom_cons_o(no)*tgridmap%area_dst(no)* &
             frac_dst(no)*re**2
         gliv_cons_o  = gliv_cons_o + mliv_cons_o(no)*tgridmap%area_dst(no)* &
             frac_dst(no)*re**2
         gelec_cons_o  = gelec_cons_o + melec_cons_o(no)*tgridmap%area_dst(no)* &
             frac_dst(no)*re**2
         gmfc_cons_o  = gmfc_cons_o + mmfc_cons_o(no)*tgridmap%area_dst(no)* &
             frac_dst(no)*re**2
         gmin_cons_o  = gmin_cons_o + mmin_cons_o(no)*tgridmap%area_dst(no)* &
             frac_dst(no)*re**2

     end do

     write (6,*) 'Successfully made sector water withdrawal and consumption for month ', m
     call shr_sys_flush(6)

  enddo
  write (6,*)

  ! Close input file
  call check_ret(nf_close(ncidi), subname)

  ! Deallocate dynamic memory
  deallocate(mdom_withd_i)
  deallocate(mliv_withd_i)
  deallocate(melec_withd_i)
  deallocate(mmfc_withd_i)
  deallocate(mmin_withd_i)

  deallocate(mdom_cons_i)
  deallocate(mliv_cons_i)
  deallocate(melec_cons_i)
  deallocate(mmfc_cons_i)
  deallocate(mmin_cons_i)

  deallocate(mdom_withd_o)
  deallocate(mliv_withd_o)
  deallocate(melec_withd_o)
  deallocate(mmfc_withd_o)
  deallocate(mmin_withd_o)

  deallocate(mdom_cons_o)
  deallocate(mliv_cons_o)
  deallocate(melec_cons_o)
  deallocate(mmfc_cons_o)
  deallocate(mmin_cons_o)
  
  deallocate(frac_dst)

  call gridmap_clean(tgridmap)
  call domain_clean(tdomain) 

end subroutine mksectorwater
end module mksectorwaterMod
