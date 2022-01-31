module mkvocefMod

  !-----------------------------------------------------------------------
  ! Make VOC percentage emissions for surface dataset
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use mkvarpar	
  use mkvarctl    
  use mkncdio

  implicit none
  private

  public :: mkvocef  ! Get the percentage emissions for VOC for different land cover types

contains

  !-----------------------------------------------------------------------
  subroutine mkvocef(ldomain, mapfname, datfname, ndiag, &
       ef_btr_o, ef_fet_o, ef_fdt_o, ef_shr_o, ef_grs_o, ef_crp_o)

    ! make volatile organic coumpunds (VOC) emission factors.

    ! input/output variables
    type(domain_type) , intent(in) :: ldomain
    character(len=*)  , intent(in) :: mapfname    ! input mapping file name
    character(len=*)  , intent(in) :: datfname    ! input data file name
    integer           , intent(in) :: ndiag       ! unit number for diagnostic output
    real(r8)          , intent(out):: ef_btr_o(:) ! output grid: EFs for broadleaf trees
    real(r8)          , intent(out):: ef_fet_o(:) ! output grid: EFs for fineleaf evergreen
    real(r8)          , intent(out):: ef_fdt_o(:) ! output grid: EFs for fineleaf deciduous
    real(r8)          , intent(out):: ef_shr_o(:) ! output grid: EFs for shrubs
    real(r8)          , intent(out):: ef_grs_o(:) ! output grid: EFs for grasses
    real(r8)          , intent(out):: ef_crp_o(:) ! output grid: EFs for crops
    !
    ! local variables:
    real(r8), allocatable :: ef_btr_i(:)      ! input grid: EFs for broadleaf trees
    real(r8), allocatable :: ef_fet_i(:)      ! input grid: EFs for fineleaf evergreen
    real(r8), allocatable :: ef_fdt_i(:)      ! input grid: EFs for fineleaf deciduous
    real(r8), allocatable :: ef_shr_i(:)      ! input grid: EFs for shrubs
    real(r8), allocatable :: ef_grs_i(:)      ! input grid: EFs for grasses
    real(r8), allocatable :: ef_crp_i(:)      ! input grid: EFs for crops
    real(r8), allocatable :: frac_dst(:)      ! output fractions
    real(r8), allocatable :: mask_r8(:)  ! float of tdomain%mask
    real(r8) :: sum_fldo                      ! global sum of dummy input fld
    real(r8) :: sum_fldi                      ! global sum of dummy input fld
    integer  :: k,n,no,ni,ns_o,ns_i           ! indices
    integer  :: ncid,dimid,varid              ! input netCDF id's
    integer  :: ier                           ! error status
    real(r8) :: relerr = 0.00001_r8           ! max error: sum overlap wts ne 1
    character(len=32) :: subname = 'mkvocef'
    !-----------------------------------------------------------------------

    write (6,*) 'Attempting to make VOC emission factors .....'

    ns_o = ldomain%ns

    ! -----------------------------------------------------------------
    ! Read input Emission Factors
    ! -----------------------------------------------------------------

    ! Obtain input grid info, read local fields

    call domain_read(tdomain,datfname)
    ns_i = tdomain%ns
    allocate(ef_btr_i(ns_i), ef_fet_i(ns_i), ef_fdt_i(ns_i), &
             ef_shr_i(ns_i), ef_grs_i(ns_i), ef_crp_i(ns_i), &
             frac_dst(ns_o), stat=ier)
    if (ier/=0) call shr_sys_abort()

    write (6,*) 'Open VOC file: ', trim(datfname)
    call check_ret(nf_open(datfname, 0, ncid), subname)
    call check_ret(nf_inq_varid (ncid, 'ef_btr', varid), subname)
    call check_ret(nf_get_var_double(ncid, varid, ef_btr_i), subname)
    call check_ret(nf_inq_varid (ncid, 'ef_fet', varid), subname)
    call check_ret(nf_get_var_double(ncid, varid, ef_fet_i), subname)
    call check_ret(nf_inq_varid (ncid, 'ef_fdt', varid), subname)
    call check_ret(nf_get_var_double(ncid, varid, ef_fdt_i), subname)
    call check_ret(nf_inq_varid (ncid, 'ef_shr', varid), subname)
    call check_ret(nf_get_var_double(ncid, varid, ef_shr_i), subname)
    call check_ret(nf_inq_varid (ncid, 'ef_grs', varid), subname)
    call check_ret(nf_get_var_double(ncid, varid, ef_grs_i), subname)
    call check_ret(nf_inq_varid (ncid, 'ef_crp', varid), subname)
    call check_ret(nf_get_var_double(ncid, varid, ef_crp_i), subname)
    call check_ret(nf_close(ncid), subname)

    ! Area-average percent cover on input grid to output grid 
    ! and correct according to land landmask
    ! Note that percent cover is in terms of total grid area.

    call gridmap_mapread(tgridmap, mapfname )

    ! Error checks for domain and map consistencies

    call domain_checksame( tdomain, ldomain, tgridmap )

    ! Obtain frac_dst
    call gridmap_calc_frac_dst(tgridmap, tdomain%mask, frac_dst)

    ! Do mapping from input to output grid

    call gridmap_areaave_srcmask(tgridmap, ef_btr_i, ef_btr_o, nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)
    call gridmap_areaave_srcmask(tgridmap, ef_fet_i, ef_fet_o, nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)
    call gridmap_areaave_srcmask(tgridmap, ef_fdt_i, ef_fdt_o, nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)
    call gridmap_areaave_srcmask(tgridmap, ef_shr_i, ef_shr_o, nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)
    call gridmap_areaave_srcmask(tgridmap, ef_grs_i, ef_grs_o, nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)
    call gridmap_areaave_srcmask(tgridmap, ef_crp_i, ef_crp_o, nodata=0._r8, mask_src=tdomain%mask, frac_dst=frac_dst)

    ! Check for conservation

    do no = 1, ns_o
       if ( ef_btr_o(no) < 0._r8 ) then
          write (6,*) 'MKVOCEF error: EF btr = ',ef_btr_o(no), ' is negative for no = ',no
          call shr_sys_abort()
       end if
       if ( ef_fet_o(no) < 0._r8 ) then
          write (6,*) 'MKVOCEF error: EF fet = ',ef_fet_o(no), ' is negative for no = ',no
          call shr_sys_abort()
       end if
       if ( ef_fdt_o(no) < 0._r8 ) then
          write (6,*) 'MKVOCEF error: EF fdt = ',ef_fdt_o(no), ' is negative for no = ',no
          call shr_sys_abort()
       end if
       if ( ef_shr_o(no) < 0._r8 ) then
          write (6,*) 'MKVOCEF error: EF shr = ',ef_shr_o(no), ' is negative for no = ',no
          call shr_sys_abort()
       end if
       if ( ef_grs_o(no) < 0._r8 ) then
          write (6,*) 'MKVOCEF error: EF grs = ',ef_grs_o(no), ' is negative for no = ',no
          call shr_sys_abort()
       end if
       if ( ef_crp_o(no) < 0._r8 ) then
          write (6,*) 'MKVOCEF error: EF crp = ',ef_crp_o(no), ' is negative for no = ',no
          call shr_sys_abort()
       end if
    enddo

    ! -----------------------------------------------------------------
    ! Error check1
    ! Compare global sum fld_o to global sum fld_i.
    ! -----------------------------------------------------------------

    ! Global sum of output field -- must multiply by fraction of
    ! output grid that is land as determined by input grid

    allocate(mask_r8(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    mask_r8 = tdomain%mask
    call gridmap_check( tgridmap, mask_r8, frac_dst, subname )

    write (6,*) 'Successfully made VOC Emission Factors'
    write (6,*)

    ! Deallocate dynamic memory

    deallocate ( ef_btr_i, ef_fet_i, ef_fdt_i, &
                 ef_shr_i, ef_grs_i, ef_crp_i, frac_dst, mask_r8 )

  end subroutine mkvocef

end module mkvocefMod
