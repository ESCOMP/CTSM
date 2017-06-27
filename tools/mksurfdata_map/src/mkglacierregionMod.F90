module mkglacierregionMod

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: mkglacierregionMod
  !
  ! !DESCRIPTION:
  ! make glacier region ID
  !
  ! !REVISION HISTORY:
  ! Author: Bill Sacks
  !
  !-----------------------------------------------------------------------
  !
  ! !USES:
  use shr_sys_mod , only : shr_sys_flush
  implicit none

  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public mkglacierregion  ! make glacier region ID
  !
  !EOP

contains

  !-----------------------------------------------------------------------
  subroutine mkglacierregion(ldomain, mapfname, datfname, ndiag, &
       glacier_region_o)
    !
    ! !DESCRIPTION:
    ! Make glacier region ID
    !
    ! Regridding is done by finding the max index that overlaps each destination cell,
    ! without regard to the weight of overlap or dominance of each overlapping index.
    !
    ! !USES:
    use mkdomainMod, only : domain_type, domain_clean, domain_read, domain_checksame
    use mkgridmapMod
    use mkncdio
    use mkindexmapMod, only : get_max_indices
    use mkdiagnosticsMod, only : output_diagnostics_index
    use mkchecksMod, only : min_bad
    !
    ! !ARGUMENTS:
    type(domain_type), intent(in) :: ldomain
    character(len=*) , intent(in)  :: mapfname            ! input mapping file name
    character(len=*) , intent(in)  :: datfname            ! input data file name
    integer          , intent(in)  :: ndiag               ! unit number for diag out
    integer          , intent(out) :: glacier_region_o(:) ! glacier region
    !
    ! !LOCAL VARIABLES:
    type(gridmap_type)   :: tgridmap
    type(domain_type)    :: tdomain             ! local domain
    integer, allocatable :: glacier_region_i(:) ! glacier region on input grid
    integer              :: ncid,varid          ! input netCDF id's
    integer              :: ier                 ! error status
    integer              :: max_region          ! max region ID

    character(len=*), parameter :: subname = 'mkglacierregion'
    !-----------------------------------------------------------------------

    write (6,*) 'Attempting to make glacier region .....'
    call shr_sys_flush(6)

    ! ------------------------------------------------------------------------
    ! Read domain and mapping information, check for consistency
    ! ------------------------------------------------------------------------

    call domain_read(tdomain, datfname)

    call gridmap_mapread(tgridmap, mapfname)
    call gridmap_check(tgridmap, subname)

    call domain_checksame(tdomain, ldomain, tgridmap)

    ! ------------------------------------------------------------------------
    ! Open input file, allocate memory for input data
    ! ------------------------------------------------------------------------

    write (6,*) 'Open glacier region raw data file: ', trim(datfname)
    call check_ret(nf_open(datfname, 0, ncid), subname)

    allocate(glacier_region_i(tdomain%ns), stat=ier)
    if (ier/=0) call abort()

    ! ------------------------------------------------------------------------
    ! Regrid glacier_region
    ! ------------------------------------------------------------------------

    call check_ret(nf_inq_varid(ncid, 'GLACIER_REGION', varid), subname)
    call check_ret(nf_get_var_int(ncid, varid, glacier_region_i), subname)
    if (min_bad(glacier_region_i, 0, 'GLACIER_REGION')) then
       stop
    end if

    call get_max_indices( &
         gridmap = tgridmap, &
         src_array = glacier_region_i, &
         dst_array = glacier_region_o, &
         nodata = 0)

    max_region = maxval(glacier_region_i)
    call output_diagnostics_index(glacier_region_i, glacier_region_o, tgridmap, &
         'Glacier Region ID', 0, max_region, ndiag)

    ! ------------------------------------------------------------------------
    ! Deallocate dynamic memory & other clean up
    ! ------------------------------------------------------------------------

    call check_ret(nf_close(ncid), subname)
    call domain_clean(tdomain)
    call gridmap_clean(tgridmap)
    deallocate(glacier_region_i)

    write (6,*) 'Successfully made glacier region'
    write (6,*)
    call shr_sys_flush(6)

  end subroutine mkglacierregion

end module mkglacierregionMod
