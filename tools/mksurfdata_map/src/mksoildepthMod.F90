module mksoildepthMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mksoildepthMod
!
! !DESCRIPTION:
! make fraction soildepth from input soildepth data
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
  public mksoildepth           ! regrid soildepth data
!
!EOP
!===============================================================
contains
!===============================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksoildepth
!
! !INTERFACE:
subroutine mksoildepth(ldomain, mapfname, datfname, ndiag, soildepth_o)
!
! !DESCRIPTION:
! make soildepth
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
  real(r8)          , intent(out):: soildepth_o(:) ! output grid: fraction soildepth
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
  integer  :: ncid,varid                      ! input netCDF id's
  integer  :: ier                             ! error status
 
  real(r8), parameter :: min_valid = 0._r8          ! minimum valid value
  real(r8), parameter :: max_valid = 100.000001_r8  ! maximum valid value
  character(len=32) :: subname = 'mksoildepth'
  character(len=32) :: varname 
  integer  :: varnum
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make soildepth .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read domain and mapping information, check for consistency
  ! -----------------------------------------------------------------

  call domain_read( tdomain, datfname )

  call gridmap_mapread( tgridmap, mapfname )
  call gridmap_check( tgridmap, subname )

  call domain_checksame( tdomain, ldomain, tgridmap )

  ! -----------------------------------------------------------------
  ! Open input file, allocate memory for input data
  ! -----------------------------------------------------------------

  write(6,*)'Open soildepth file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)

  allocate(data_i(tdomain%ns), stat=ier)
  if (ier/=0) call abort()

  ! -----------------------------------------------------------------
  ! Regrid soildepth
  ! -----------------------------------------------------------------

  varnum = 1
  select case (varnum)
  case(1)
     varname = 'Avg_Depth_Median'
  case(2)
     varname = 'Avg_Depth_Mean'
  case(3)
     varname = 'Upland_Valley_Depth_Median'
  case(4)
     varname = 'Upland_Valley_Depth_Mean'
  case(5)
     varname = 'Upland_Hillslope_Depth_Median'
  case(6)
     varname = 'Upland_Hillslope_Depth_Mean'
  case(7)
     varname = 'Lowland_Depth_Mean'
  case(8)
     varname = 'Lowland_Depth_Mean'
  end select

!  call check_ret(nf_inq_varid (ncid, 'Avg_Depth_Median', varid), subname)
  call check_ret(nf_inq_varid (ncid, varname, varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
  call gridmap_areaave(tgridmap, data_i, soildepth_o, nodata=0._r8)

  ! Check validity of output data
  if (min_bad(soildepth_o, min_valid, 'soildepth') .or. &
      max_bad(soildepth_o, max_valid, 'soildepth')) then
     stop
  end if
  
  call output_diagnostics_area(data_i, soildepth_o, tgridmap, "Soildepth", percent=.false., ndiag=ndiag)
  
  ! -----------------------------------------------------------------
  ! Close files and deallocate dynamic memory
  ! -----------------------------------------------------------------

  call check_ret(nf_close(ncid), subname)
  call domain_clean(tdomain) 
  call gridmap_clean(tgridmap)
  deallocate (data_i)

  write (6,*) 'Successfully made soildepth'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mksoildepth


end module mksoildepthMod
