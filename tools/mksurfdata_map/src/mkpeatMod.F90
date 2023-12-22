module mkpeatMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkpeatMod
!
! !DESCRIPTION:
! make fraction peat from input peat data
!
! !REVISION HISTORY:
! Author: Sam Levis and Bill Sacks
!
!-----------------------------------------------------------------------
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use mkdomainMod , only : domain_checksame

  implicit none

  private

! !PUBLIC MEMBER FUNCTIONS:
  public mkpeat           ! regrid peat data
!
!EOP
!===============================================================
contains
!===============================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkpeat
!
! !INTERFACE:
subroutine mkpeat(ldomain, mapfname, datfname, ndiag, peat_o)
!
! !DESCRIPTION:
! make peat
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkncdio
  use mkdiagnosticsMod, only : output_diagnostics_area
  use mkchecksMod, only : min_bad, max_bad
!
! !ARGUMENTS:
  
  implicit none
  type(domain_type) , intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname  ! input mapping file name
  character(len=*)  , intent(in) :: datfname  ! input data file name
  integer           , intent(in) :: ndiag     ! unit number for diag out
  real(r8)          , intent(out):: peat_o(:) ! output grid: fraction peat
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Sam Levis and Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)    :: tgridmap
  type(domain_type)     :: tdomain            ! local domain
  real(r8), allocatable :: data_i(:)          ! data on input grid
  real(r8), allocatable :: frac_dst(:)        ! output fractions
  real(r8), allocatable :: mask_r8(:)  ! float of tdomain%mask
  integer  :: ncid,varid                      ! input netCDF id's
  integer  :: ier                             ! error status
 
  real(r8), parameter :: min_valid = 0._r8          ! minimum valid value
  real(r8), parameter :: max_valid = 100.000001_r8  ! maximum valid value
  character(len=32) :: subname = 'mkpeat'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make peat .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read domain and mapping information, check for consistency
  ! -----------------------------------------------------------------

  call domain_read( tdomain, datfname )

  call gridmap_mapread( tgridmap, mapfname )

  ! Obtain frac_dst
  allocate(frac_dst(ldomain%ns), stat=ier)
  if (ier/=0) call abort()
  call gridmap_calc_frac_dst(tgridmap, tdomain%mask, frac_dst)

  allocate(mask_r8(tdomain%ns), stat=ier)
  if (ier/=0) call abort()
  mask_r8 = tdomain%mask
  call gridmap_check( tgridmap, mask_r8, frac_dst, subname )

  call domain_checksame( tdomain, ldomain, tgridmap )

  ! -----------------------------------------------------------------
  ! Open input file, allocate memory for input data
  ! -----------------------------------------------------------------

  write(6,*)'Open peat file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)

  allocate(data_i(tdomain%ns), stat=ier)
  if (ier/=0) call abort()

  ! -----------------------------------------------------------------
  ! Regrid peat
  ! -----------------------------------------------------------------

  call check_ret(nf_inq_varid (ncid, 'peatf', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
  call gridmap_areaave_srcmask(tgridmap, data_i, peat_o, nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)

  ! Check validity of output data
  if (min_bad(peat_o, min_valid, 'peat') .or. &
      max_bad(peat_o, max_valid, 'peat')) then
     call abort()
  end if
  
  call output_diagnostics_area(data_i, peat_o, tgridmap, "Peat", percent=.false., ndiag=ndiag, mask_src=tdomain%mask, frac_dst=frac_dst)
  
  ! -----------------------------------------------------------------
  ! Close files and deallocate dynamic memory
  ! -----------------------------------------------------------------

  call check_ret(nf_close(ncid), subname)
  call domain_clean(tdomain) 
  call gridmap_clean(tgridmap)
  deallocate (data_i)
  deallocate (frac_dst)
  deallocate (mask_r8)

  write (6,*) 'Successfully made peat'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mkpeat


end module mkpeatMod
