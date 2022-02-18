module mkVICparamsMod

  !-----------------------------------------------------------------------
  ! make parameters for VIC
  !-----------------------------------------------------------------------
  !
  use ESMF
  use pio
  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_sys_mod      , only : shr_sys_abort
  use mkpioMod         , only : mkpio_get_rawdata
  use mkpioMod         , only : mkpio_iodesc_rawdata, pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod        , only : regrid_rawdata, create_routehandle_r8, get_meshareas
  use mkutilsMod       , only : chkerr
  use mkvarctl         , only : root_task, ndiag, mpicom, MPI_INTEGER, MPI_MAX
  use mkvarctl         , only : soil_color_override, unsetcol
  use mkvarpar         , only : re
  use mkchecksMod      , only : min_bad
  use mkdiagnosticsMod , only : output_diagnostics_continuous

  implicit none
  private

  public :: mkVICparams            ! make VIC parameters

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mkVICparams(file_mesh_i, file_data_i, mesh_o, binfl_o, ws_o, dsmax_o, ds_o, rc)
    !
    ! make VIC parameters
    !
    ! input/output variables
    character(len=*)  , intent(in)  :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)  :: file_data_i ! input data file name
    type(ESMF_Mesh)   , intent(in)  :: mesh_o      ! model mesho
    real(r8)         , intent(out)  :: binfl_o(:)  ! output VIC b parameter for the Variable Infiltration Capacity Curve (unitless)
    real(r8)         , intent(out)  :: ws_o(:)     ! output VIC Ws parameter for the ARNO curve (unitless)
    real(r8)         , intent(out)  :: dsmax_o(:)  ! output VIC Dsmax parameter for the ARNO curve (mm/day)
    real(r8)         , intent(out)  :: ds_o(:)     ! output VIC Ds parameter for the ARNO curve (unitless)
    integer           , intent(out) :: rc
    !
    ! !LOCAL VARIABLES:
    real(r8), allocatable :: data_i(:)   ! data on input grid
    real(r8), allocatable :: frac_dst(:) ! output fractions
    real(r8), allocatable :: mask_r8(:)  ! float of tdomain%mask
    integer               :: ncid,varid  ! input netCDF id's
    integer               :: ier         ! error status
    real(r8), parameter   :: min_valid_binfl = 0._r8
    real(r8), parameter   :: min_valid_ws    = 0._r8
    real(r8), parameter   :: min_valid_dsmax = 0._r8
    real(r8), parameter   :: min_valid_ds    = 0._r8
    character(len=*), parameter :: subname = 'mkVICparams'
    !-----------------------------------------------------------------------

    write (6,*) 'Attempting to make VIC parameters.....'

    ! -----------------------------------------------------------------
    ! Read domain and mapping information, check for consistency
    ! -----------------------------------------------------------------

    call domain_read(tdomain,datfname)

    call gridmap_mapread(tgridmap, mapfname )

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

    write(6,*)'Open VIC parameter file: ', trim(datfname)
    call check_ret(nf_open(datfname, 0, ncid), subname)

    allocate(data_i(tdomain%ns), stat=ier)
    if (ier/=0) call abort()

    ! -----------------------------------------------------------------
    ! Regrid binfl
    ! -----------------------------------------------------------------

    call check_ret(nf_inq_varid (ncid, 'binfl', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
    call gridmap_areaave_srcmask(tgridmap, data_i, binfl_o, nodata=0.1_r8, mask_src=tdomain%mask, frac_dst=frac_dst)

    ! Check validity of output data
    if (min_bad(binfl_o, min_valid_binfl, 'binfl')) then
       call abort()
    end if

    call output_diagnostics_continuous(data_i, binfl_o, tgridmap, "VIC b parameter", "unitless", ndiag, tdomain%mask, frac_dst)

    ! -----------------------------------------------------------------
    ! Regrid Ws
    ! -----------------------------------------------------------------

    call check_ret(nf_inq_varid (ncid, 'Ws', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
    call gridmap_areaave_srcmask(tgridmap, data_i, ws_o, nodata=0.75_r8, mask_src=tdomain%mask, frac_dst=frac_dst)

    ! Check validity of output data
    if (min_bad(ws_o, min_valid_ws, 'Ws')) then
       call abort()
    end if

    call output_diagnostics_continuous(data_i, ws_o, tgridmap, "VIC Ws parameter", "unitless", ndiag, tdomain%mask, frac_dst)

    ! -----------------------------------------------------------------
    ! Regrid Dsmax
    ! -----------------------------------------------------------------

    call check_ret(nf_inq_varid (ncid, 'Dsmax', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
    call gridmap_areaave_srcmask(tgridmap, data_i, dsmax_o, nodata=10._r8, mask_src=tdomain%mask, frac_dst=frac_dst)

    ! Check validity of output data
    if (min_bad(dsmax_o, min_valid_dsmax, 'Dsmax')) then
       call abort()
    end if

    call output_diagnostics_continuous(data_i, dsmax_o, tgridmap, "VIC Dsmax parameter", "mm/day", ndiag, tdomain%mask, frac_dst)

    ! -----------------------------------------------------------------
    ! Regrid Ds
    ! -----------------------------------------------------------------

    call check_ret(nf_inq_varid (ncid, 'Ds', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
    call gridmap_areaave_srcmask(tgridmap, data_i, ds_o, nodata=0.1_r8, mask_src=tdomain%mask, frac_dst=frac_dst)

    ! Check validity of output data
    if (min_bad(ds_o, min_valid_ds, 'Ds')) then
       call abort()
    end if

    call output_diagnostics_continuous(data_i, ds_o, tgridmap, "VIC Ds parameter", "unitless", ndiag, tdomain%mask, frac_dst)

    ! -----------------------------------------------------------------
    ! Close files and deallocate dynamic memory
    ! -----------------------------------------------------------------

    call check_ret(nf_close(ncid), subname)
    call domain_clean(tdomain) 
    call gridmap_clean(tgridmap)
    deallocate (data_i)
    deallocate (frac_dst)
    deallocate (mask_r8)

    write (6,*) 'Successfully made VIC parameters'
    write (6,*)

  end subroutine mkVICparams

end module mkVICparamsMod
