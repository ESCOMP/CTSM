module mkirrigation_methodMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkirrigation_method
!
! !DESCRIPTION:
! Make irrigation method data
!
! !REVISION HISTORY:
! Author: Yi Yao
!
!EOP
!-----------------------------------------------------------------------
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use mkdomainMod , only : domain_checksame
  use mkvarctl    

  implicit none

  private

  public  :: mkirrigation_method
  

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkirrigation_method
!
! !INTERFACE:
subroutine mkirrigation_method(ldomain, mapfname, datfname, ndiag, ncido)
!
! !DESCRIPTION:
! Make irrigation_method data
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
! !REVISION HISTORY:
! Author: 
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)    :: tgridmap
  type(domain_type)    :: tdomain          ! local domain
  integer  :: num_cft_i                      ! number of plant types on input
  

  real(r8), allocatable :: flood_fraction_o(:,:)      ! flood fraction
  real(r8), allocatable :: sprinkler_fraction_o(:,:)      ! sprinkler fraction
  real(r8), allocatable :: drip_fraction_o(:,:)     ! drip fraction
  
  integer, allocatable :: irrig_method_o(:,:)     ! irrig method
  
  
  
  real(r8), allocatable :: flood_fraction_i(:,:)      ! flood fraction in
  real(r8), allocatable :: sprinkler_fraction_i(:,:)      ! sprinkler fraction in
  real(r8), allocatable :: drip_fraction_i(:,:)     ! drip fraction in
  
  real(r8), allocatable :: frac_dst(:)        ! output fractions: same as frac_dst
  
  real(r8), allocatable :: mask_src(:)      ! input grid: mask (0, 1)

  integer  :: ni,no,ns_i,ns_o               ! indices
  integer  :: k,l,n,m                       ! indices
  integer  :: ncidi,dimid,varid             ! input netCDF id's
  integer  :: ndimsi,ndimso                 ! netCDF dimension sizes
  integer  :: num_cft = 64  
  integer  :: dimids(4)                     ! netCDF dimension ids
  integer  :: bego(3),leno(3)               ! netCDF bounds
  integer  :: begi(4),leni(4)               ! netCDF bounds 
  integer  :: nmethod                       ! number of input time samples
  integer  :: ier                           ! error status
  real(r8) :: relerr = 0.00001              ! max error: sum overlap wts ne 1
  character(len=256) :: name                ! name of attribute
  character(len=256) :: unit                ! units of attribute
  character(len= 32) :: subname = 'mkirrigation_method'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make irrigation_method .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  ns_o = ldomain%ns

  call domain_read(tdomain,datfname)
  ns_i = tdomain%ns

  write (6,*) 'Open irrigation method file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncidi), subname)		
  call check_ret(nf_inq_dimid(ncidi, 'cft', dimid), subname)	
  call check_ret(nf_inq_dimlen(ncidi, dimid, num_cft_i), subname)
  write(6,*) 'num_cft_i = ', num_cft_i
  
  call check_ret(nf_inq_dimid(ncidi, 'irrig_method', dimid), subname)
  call check_ret(nf_inq_dimlen(ncidi, dimid, nmethod), subname)
  
  if (num_cft_i /= num_cft) then
     write(6,*) 'WARNING: ' // trim(subname) // '(): parameter num_cft = ', num_cft, &
          'does not equal input dataset num_cft = ', num_cft_i
     write(6,*)'This inconsistency will stop the program.'
     stop
  endif
  if (nmethod /= 3) then	
     write(6,*)'MKIRRIGATION_METHOD: must have 3 method samples on input data'
     call abort()
  endif

  ! NOTE - close data set at bottom of routine

  ! Dynamic allocation of variables

  allocate(flood_fraction_i(ns_i,num_cft_i),  &
           sprinkler_fraction_i(ns_i,num_cft_i),  &
           drip_fraction_i(ns_i,num_cft_i), &
		   
           mask_src(ns_i),         &
		   frac_dst(ns_o),           &
           flood_fraction_o(ns_o, num_cft),  &
           sprinkler_fraction_o(ns_o, num_cft),  &
           drip_fraction_o(ns_o, num_cft), &
		   
           irrig_method_o(ns_o,num_cft), &
		   stat=ier )
  if (ier /= 0) then
     write(6,*)'mkirrigation_method allocation error'; call abort()
  end if
  ! Determine mapping weights and map

  call gridmap_mapread(tgridmap, mapfname)

  ! Error checks for domain and map consistencies
  
  call domain_checksame( tdomain, ldomain, tgridmap )					

  ! Determine number of dimensions in input by IRRIGATION_METHOD_FRACTION

  call check_ret(nf_inq_varid(ncidi, 'IRRIGATION_METHOD_FRACTION', varid), subname)
  call check_ret(nf_inq_vardimid(ncidi, varid, dimids), subname)
  call check_ret(nf_inq_varndims(ncidi, varid, ndimsi), subname)
  
  begi(1) = 1
  begi(2) = 1
  begi(3) = 1
  leni(4) = 1
  
  call check_ret(nf_inq_dimlen(ncidi, dimids(1), leni(1)), subname)
  call check_ret(nf_inq_dimlen(ncidi, dimids(2), leni(2)), subname)
  call check_ret(nf_inq_dimlen(ncidi, dimids(3), leni(3)), subname)
  

  ! Determine number of dimensions in output by querying irrigation_method

  call check_ret(nf_inq_varid(ncido, 'irrigation_method', varid), subname)
  call check_ret(nf_inq_varndims(ncido, varid, ndimso), subname)
  call check_ret(nf_inq_vardimid(ncido, varid, dimids), subname)
  
  
  bego(1) = 1
  bego(2) = 1
  !leno(3) = 1
  bego(3) = 1
  call check_ret(nf_inq_dimlen(ncido, dimids(1), leno(1)), subname)
  call check_ret(nf_inq_dimlen(ncido, dimids(2), leno(2)), subname)
  call check_ret(nf_inq_dimlen(ncido, dimids(3), leno(3)), subname)
  
  
  
  
  begi(4) = 1
  call check_ret(nf_inq_varid (ncidi, 'IRRIGATION_METHOD_FRACTION', varid), subname)
  call check_ret(nf_get_vara_double(ncidi, varid, begi(1:ndimsi), leni(1:ndimsi), &
       flood_fraction_i), subname)
	   
  begi(4) = 2
  call check_ret(nf_inq_varid(ncidi, 'IRRIGATION_METHOD_FRACTION', varid), subname)
  call check_ret(nf_get_vara_double(ncidi, varid, begi(1:ndimsi), leni(1:ndimsi), &
       sprinkler_fraction_i), subname)
	   
  begi(4) = 3
  call check_ret(nf_inq_varid(ncidi, 'IRRIGATION_METHOD_FRACTION', varid), subname)
  call check_ret(nf_get_vara_double(ncidi, varid, begi(1:ndimsi), leni(1:ndimsi), &
       drip_fraction_i), subname)
  


  
  flood_fraction_o(:,:)  = 0.
  sprinkler_fraction_o(:,:)  = 0.
  drip_fraction_o(:,:) = 0.
  
  
  irrig_method_o(:,:) = 0.
  
  call gridmap_calc_frac_dst(tgridmap, tdomain%mask, frac_dst)
  ! Loop over cft types to do mapping

  do l = 1, num_cft_i
     mask_src(:) = 1._r8 
     call gridmap_areaave_srcmask(tgridmap, flood_fraction_i(:,l) , flood_fraction_o(:,l) , nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)	
     call gridmap_areaave_srcmask(tgridmap, sprinkler_fraction_i(:,l) , sprinkler_fraction_o(:,l) , nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)
     call gridmap_areaave_srcmask(tgridmap, drip_fraction_i(:,l), drip_fraction_o(:,l), nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)
  end do
  
  
  ! Determine the irrigation method based on the dominant irrigation method
  do k = 1, ns_o
	 do l = 1, num_cft_i
		if (flood_fraction_o(k, l) > sprinkler_fraction_o(k, l) .and. flood_fraction_o(k, l) > drip_fraction_o(k, l)) then
			irrig_method_o(k, l) = 3
		else if (sprinkler_fraction_o(k, l) > flood_fraction_o(k, l) .and. sprinkler_fraction_o(k, l) > drip_fraction_o(k, l)) then
			irrig_method_o(k, l) = 2
	    else if (drip_fraction_o(k, l) > flood_fraction_o(k, l) .and. drip_fraction_o(k, l) > sprinkler_fraction_o(k, l)) then
   			irrig_method_o(k, l) = 1
		else
			irrig_method_o(k, l) = 0
		end if
	 end do
  end do

  call check_ret(nf_inq_varid(ncido, 'irrigation_method', varid), subname)
  !call check_ret(nf_put_vara_int(ncido, varid, bego, leno, irrig_method_o), subname)
  call check_ret(nf_put_var_int(ncido, varid, irrig_method_o), subname)

  
  call check_ret(nf_close(ncidi), subname)


  ! Deallocate dynamic memory
  deallocate(flood_fraction_i)
  deallocate(sprinkler_fraction_i)
  deallocate(drip_fraction_i)
  
  deallocate(mask_src)
  deallocate(flood_fraction_o)
  deallocate(sprinkler_fraction_o)
  deallocate(drip_fraction_o)
  deallocate(irrig_method_o)
  

  call gridmap_clean(tgridmap)
  call domain_clean(tdomain) 

end subroutine mkirrigation_method

  
end module mkirrigation_methodMod
