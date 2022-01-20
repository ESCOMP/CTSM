module mkfileMod

  implicit none
  private

  public :: mkfile

!=================================================================================
contains
!=================================================================================

#ifdef TODO
  subroutine mkfile(vm, nx, ny, fname, harvdata, dynlanduse, pioid)
#else
  subroutine mkfile(vm, nx, ny, fname, dynlanduse, pioid)
#endif

    use ESMF
    use pio
    use shr_kind_mod , only : r8 => shr_kind_r8
    use shr_sys_mod  , only : shr_sys_getenv
    use fileutils    , only : get_filename
    use mkvarpar     , only : nlevsoi, numrad, numstdpft
#ifdef TODO
    use mkurbanparMod, only : numurbl, nlevurb
    use mkglcmecMod  , only : nglcec
    use mkpftMod     , only : mkpftAtt
    use mksoilMod    , only : mksoilAtt
    use mkharvestMod , only : mkharvest_fieldname, mkharvest_numtypes, mkharvest_longname
    use mkharvestMod , only : mkharvest_units, harvestDataType
#endif
    use mkpioMod ! TODO: add only
    use mkvarctl 

    implicit none

    ! input/output variables
    type(ESMF_VM)    , intent(inout) :: vm
    integer          , intent(in)    :: nx
    integer          , intent(in)    :: ny
    character(len=*) , intent(in)    :: fname
    logical          , intent(in)    :: dynlanduse
#ifdef TODO    
    type(harvestDataType) , intent(in) :: harvdata
#endif

    ! local variables
    type(file_desc_t)    :: pioid
    integer              :: j                  ! index
    integer              :: dimid              ! temporary
    integer              :: values(8)          ! temporary
    character(len=256)   :: str                ! global attribute string
    character(len=256)   :: name               ! name of attribute
    character(len=256)   :: unit               ! units of attribute
    character(len= 18)   :: datetime           ! temporary
    character(len=  8)   :: date               ! temporary
    character(len= 10)   :: time               ! temporary
    character(len=  5)   :: zone               ! temporary
    integer              :: ier                ! error status
    integer              :: omode              ! netCDF output mode
    integer              :: xtype              ! external type
    integer, allocatable :: ind1D(:)           ! Indices of 1D harvest variables
    integer, allocatable :: ind2D(:)           ! Indices of 2D harvest variables
    integer              :: rcode
    character(len=32)    :: subname = 'mkfile' ! subroutine name
    !-----------------------------------------------------------------------

    !---------------------------
    ! Create and open file
    !---------------------------

    ! TODO: how to translate the following into the pio call
    !call check_ret(nf_create(trim(fname), ior(nf_clobber,nf_64bit_offset), ncid), subname)
    call mkpio_wopen(pioid, trim(fname), vm, clobber=.true.)
    ! TODO: what about setting no fill values?

    !---------------------------
    ! Define dimensions.
    !---------------------------

    if (outnc_1d) then
       rcode = pio_def_dim(pioid, 'gridcell', nx, dimid)
    else
       rcode = pio_def_dim(pioid, 'lsmlon', nx, dimid)
       rcode = pio_def_dim(pioid, 'lsmlat', ny, dimid)
    end if

    if (.not. dynlanduse) then
#ifdef TODO
       rcode = pio_def_dim(pioid, 'nglcec', nglcec, dimid)
       rcode = pio_def_dim(pioid, 'nglcecp1', nglcec+1, dimid)
#endif
    else
#ifdef TODO
       rcode = pio_def_dim(pioid, 'numurbl', numurbl, dimid)
       rcode = pio_def_dim(pioid, 'nlevurb', nlevurb, dimid)
#endif
       rcode = pio_def_dim(pioid, 'numrad' , numrad, dimid)
       rcode = pio_def_dim(pioid, 'nchar'  , 256, dimid)
    end if

    !---------------------------
    ! Set global attributes.
    !---------------------------

    str = 'NCAR-CSM'
    rcode = pio_put_att(pioid, pio_global, "Conventions", trim(str))

    call date_and_time (date, time, zone, values)
    datetime(1:8) =        date(5:6) // '-' // date(7:8) // '-' // date(3:4)
    datetime(9:)  = ' ' // time(1:2) // ':' // time(3:4) // ':' // time(5:6) // ' '
    str = 'created on: ' // datetime
    rcode = pio_put_att (pioid, pio_global, 'History_Log', trim(str))

    call shr_sys_getenv ('LOGNAME', str, ier)
    rcode = pio_put_att (pioid, pio_global, 'Logname', trim(str))

    call shr_sys_getenv ('HOST', str, ier)
    rcode = pio_put_att (pioid, pio_global, 'Host', trim(str))

    str = 'Community Land Model: CLM5'
    rcode = pio_put_att (pioid, pio_global, 'Source', trim(str))
    !rcode = pio_put_att (pioid, pio_global, 'Version', trim(gitdescribe), trim(gitdescribe))
#ifdef OPT
    str = 'TRUE'
#else
    str = 'FALSE'
#endif
    rcode = pio_put_att (pioid, pio_global, 'Compiler_Optimized', trim(str))

    if ( all_urban )then
       str = 'TRUE'
       rcode = pio_put_att(pioid, pio_global, 'all_urban', trim(str))
    end if
    if ( no_inlandwet )then
       str = 'TRUE'
       rcode = pio_put_att(pioid, pio_global, 'no_inlandwet', trim(str))
    end if
    !rcode = pio_put_att_int(pioid, pio_global, 'nglcec', nglcec)

    ! Raw data file names
    str = get_filename(mksrf_fgrid_mesh)
    rcode = pio_put_att(pioid, pio_global, 'Input_grid_dataset', trim(str))
    str = trim(mksrf_gridtype)
    rcode = pio_put_att(pioid, pio_global, 'Input_gridtype', trim(str))

    if (.not. dynlanduse) then
       str = get_filename(mksrf_fvocef)
       rcode = pio_put_att(pioid, pio_global, 'VOC_EF_raw_data_file_name', trim(str))
    end if
    str = get_filename(mksrf_flakwat)
    rcode = pio_put_att(pioid, pio_global, 'Inland_lake_raw_data_file_name', trim(str))
    str = get_filename(mksrf_fwetlnd)
    rcode = pio_put_att(pioid, pio_global, 'Inland_wetland_raw_data_file_name', trim(str))
    str = get_filename(mksrf_fglacier)
    rcode = pio_put_att(pioid, pio_global, 'Glacier_raw_data_file_name', trim(str))
    str = get_filename(mksrf_fglacierregion)
    rcode = pio_put_att(pioid, pio_global, 'Glacier_region_raw_data_file_name', trim(str))
    str = get_filename(mksrf_furbtopo)
    rcode = pio_put_att(pioid, pio_global, 'Urban_Topography_raw_data_file_name', trim(str))
    str = get_filename(mksrf_furban)
    rcode = pio_put_att(pioid, pio_global, 'Urban_raw_data_file_name', trim(str))
    if (.not. dynlanduse .and. (numpft == numstdpft) ) then
       str = get_filename(mksrf_flai)
       rcode = pio_put_att(pioid, pio_global, 'Lai_raw_data_file_name', trim(str))
    end if
    str = get_filename(mksrf_fabm)
    rcode = pio_put_att(pioid, pio_global, 'agfirepkmon_raw_data_file_name', trim(str))
    str = get_filename(mksrf_fgdp)
    rcode = pio_put_att(pioid, pio_global, 'gdp_raw_data_file_name', trim(str))
    str = get_filename(mksrf_fpeat)
    rcode = pio_put_att(pioid, pio_global, 'peatland_raw_data_file_name', trim(str))
    str = get_filename(mksrf_fsoildepth)
    rcode = pio_put_att(pioid, pio_global, 'soildepth_raw_data_file_name', trim(str))
    str = get_filename(mksrf_ftopostats)
    rcode = pio_put_att(pioid, pio_global, 'topography_stats_raw_data_file_name', trim(str))
    if ( outnc_vic )then
       str = get_filename(mksrf_fvic)
       rcode = pio_put_att(pioid, pio_global,    'vic_raw_data_file_name', trim(str))
    end if

    ! Mesh file names
    str = get_filename(mksrf_fvegtyp_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_pft_file_name', trim(str))
    str = get_filename(mksrf_flakwat_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_lakwat_file', trim(str))
    str = get_filename(mksrf_fwetlnd_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_wetlnd_file', trim(str))
    str = get_filename(mksrf_fglacier_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_glacier_file', trim(str))
    str = get_filename(mksrf_fglacierregion_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_glacier_region_file', trim(str))
    str = get_filename(mksrf_fsoitex_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_soil_texture_file', trim(str))
    str = get_filename(mksrf_fsoicol_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_soil_color_file', trim(str))
    str = get_filename(mksrf_forganic_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_soil_organic_file', trim(str))
    str = get_filename(mksrf_furban_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_urban_file', trim(str))
    str = get_filename(mksrf_fmax_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_fmax_file', trim(str))
    str = get_filename(mksrf_fvocef_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_VOC_EF_file', trim(str))
    str = get_filename(mksrf_fhrvtyp_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_harvest_file', trim(str))
    if ( numpft == numstdpft )then
       str = get_filename(mksrf_flai_mesh)
       rcode = pio_put_att(pioid, pio_global, 'mesh_lai_sai_file', trim(str))
    end if
    str = get_filename(mksrf_furbtopo_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_urban_topography_file', trim(str))
    str = get_filename(mksrf_fabm_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_agfirepkmon_file', trim(str))
    str = get_filename(mksrf_fgdp_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_gdp_file', trim(str))
    str = get_filename(mksrf_fpeat_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_peatland_file', trim(str))
    str = get_filename(mksrf_fsoildepth_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_soildepth_file', trim(str))
    str = get_filename(mksrf_ftopostats_mesh)
    rcode = pio_put_att(pioid, pio_global, 'mesh_topography_stats_file', trim(str))
    if ( outnc_vic )then
       str = get_filename(mksrf_fvic_mesh)
       rcode = pio_put_att(pioid, pio_global, 'mesh_vic_file', trim(str))
    end if

    ! ----------------------------------------------------------------------
    ! Define variables
    ! ----------------------------------------------------------------------

    if ( .not. outnc_double ) then
       xtype = PIO_REAL
    else
       xtype = PIO_DOUBLE
    end if

    !call mksoilAtt( pioid, dynlanduse, xtype )

    !call mkpftAtt( pioid, dynlanduse, xtype )

    call mkpio_def_spatial_var(pioid, varname='AREA' , xtype=PIO_DOUBLE, &
         long_name='area', units='km^2')

    call mkpio_def_spatial_var(pioid, varname='LONGXY', xtype=PIO_DOUBLE, &
         long_name='longitude', units='degrees east')

    call mkpio_def_spatial_var(pioid, varname='LATIXY', xtype=PIO_double, &
         long_name='latitude', units='degrees north')

    if (.not. dynlanduse) then
       call mkpio_def_spatial_var(pioid, varname='EF1_BTR', xtype=xtype, &
            long_name='EF btr (isoprene)', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='EF1_FET', xtype=xtype, &
            long_name='EF fet (isoprene)', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='EF1_FDT', xtype=xtype, &
            long_name='EF fdt (isoprene)', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='EF1_SHR', xtype=xtype, &
            long_name='EF shr (isoprene)', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='EF1_GRS', xtype=xtype, &
            long_name='EF grs (isoprene)', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='EF1_CRP', xtype=xtype, &
            long_name='EF crp (isoprene)', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='CANYON_HWR', xtype=xtype, &
            lev1name='numurbl', &
            long_name='canyon height to width ratio', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='EM_IMPROAD', xtype=xtype, &
            lev1name='numurbl', &
            long_name='emissivity of impervious road', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='EM_PERROAD', xtype=xtype, &
            lev1name='numurbl', &
            long_name='emissivity of pervious road', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='EM_ROOF', xtype=xtype, &
            lev1name='numurbl', &
            long_name='emissivity of roof', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='EM_WALL', xtype=xtype, &
            lev1name='numurbl', &
            long_name='emissivity of wall', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='HT_ROOF', xtype=xtype, &
            lev1name='numurbl', &
            long_name='height of roof', units='meters')

       call mkpio_def_spatial_var(pioid, varname='THICK_ROOF', xtype=xtype, &
            lev1name='numurbl', &
            long_name='thickness of roof', units='meters')

       call mkpio_def_spatial_var(pioid, varname='THICK_WALL', xtype=xtype, &
            lev1name='numurbl', &
            long_name='thickness of wall', units='meters')

       call mkpio_def_spatial_var(pioid, varname='T_BUILDING_MIN', xtype=xtype, &
            lev1name='numurbl', &
            long_name='minimum interior building temperature', units='K')

       call mkpio_def_spatial_var(pioid, varname='WIND_HGT_CANYON', xtype=xtype, &
            lev1name='numurbl', &
            long_name='height of wind in canyon', units='meters')

       call mkpio_def_spatial_var(pioid, varname='WTLUNIT_ROOF', xtype=xtype, &
            lev1name='numurbl', &
            long_name='fraction of roof', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='WTROAD_PERV', xtype=xtype, &
            lev1name='numurbl', &
            long_name='fraction of pervious road', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='ALB_IMPROAD_DIR', xtype=xtype, &
            lev1name='numurbl', lev2name='numrad', &
            long_name='direct albedo of impervious road', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='ALB_IMPROAD_DIF', xtype=xtype, &
            lev1name='numurbl', lev2name='numrad', &
            long_name='diffuse albedo of impervious road', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='ALB_PERROAD_DIR', xtype=xtype, &
            lev1name='numurbl', lev2name='numrad', &
            long_name='direct albedo of pervious road', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='ALB_PERROAD_DIF', xtype=xtype, &
            lev1name='numurbl', lev2name='numrad', &
            long_name='diffuse albedo of pervious road', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='ALB_ROOF_DIR', xtype=xtype, &
            lev1name='numurbl', lev2name='numrad', &
            long_name='direct albedo of roof', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='ALB_ROOF_DIF', xtype=xtype, &
            lev1name='numurbl', lev2name='numrad', &
            long_name='diffuse albedo of roof', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='ALB_WALL_DIR', xtype=xtype, &
            lev1name='numurbl', lev2name='numrad', &
            long_name='direct albedo of wall', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='ALB_WALL_DIF', xtype=xtype, &
            lev1name='numurbl', lev2name='numrad', &
            long_name='diffuse albedo of wall', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='TK_ROOF', xtype=xtype, &
            lev1name='numurbl', lev2name='nlevurb', &
            long_name='thermal conductivity of roof', units='W/m*K')

       call mkpio_def_spatial_var(pioid, varname='TK_WALL', xtype=xtype, &
            lev1name='numurbl', lev2name='nlevurb', &
            long_name='thermal conductivity of wall', units='W/m*K')

       call mkpio_def_spatial_var(pioid, varname='TK_IMPROAD', xtype=xtype, &
            lev1name='numurbl', lev2name='nlevurb', &
            long_name='thermal conductivity of impervious road', units='W/m*K')

       call mkpio_def_spatial_var(pioid, varname='CV_ROOF', xtype=xtype, &
            lev1name='numurbl', lev2name='nlevurb', &
            long_name='volumetric heat capacity of roof', units='J/m^3*K')

       call mkpio_def_spatial_var(pioid, varname='CV_WALL', xtype=xtype, &
            lev1name='numurbl', lev2name='nlevurb', &
            long_name='volumetric heat capacity of wall', units='J/m^3*K')

       call mkpio_def_spatial_var(pioid, varname='CV_IMPROAD', xtype=xtype, &
            lev1name='numurbl', lev2name='nlevurb', &
            long_name='volumetric heat capacity of impervious road', units='J/m^3*K')

       call mkpio_def_spatial_var(pioid, varname='NLEV_IMPROAD', xtype=PIO_INT, &
            lev1name='numurbl', &
            long_name='number of impervious road layers', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='peatf', xtype=xtype, &
            long_name='peatland fraction', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='zbedrock', xtype=xtype, &
            long_name='soil depth', units='m')

       call mkpio_def_spatial_var(pioid, varname='abm', xtype=PIO_INT, &
            long_name='agricultural fire peak month', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='gdp', xtype=xtype, &
            long_name='gdp', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='SLOPE', xtype=xtype, &
            long_name='mean topographic slope', units='degrees')

       call mkpio_def_spatial_var(pioid, varname='STD_ELEV', xtype=xtype, &
            long_name='standard deviation of elevation', units='m')

       if ( outnc_vic )then
          call mkpio_def_spatial_var(pioid, varname='binfl', xtype=xtype, &
               long_name='VIC b parameter for the Variable Infiltration Capacity Curve', units='unitless')

          call mkpio_def_spatial_var(pioid, varname='Ws', xtype=xtype, &
               long_name='VIC Ws parameter for the ARNO curve', units='unitless')

          call mkpio_def_spatial_var(pioid, varname='Dsmax', xtype=xtype, &
               long_name='VIC Dsmax parameter for the ARNO curve', units='mm/day')

          call mkpio_def_spatial_var(pioid, varname='Ds', xtype=xtype, &
               long_name='VIC Ds parameter for the ARNO curve', units='unitless')

       end if
       call mkpio_def_spatial_var(pioid, varname='LAKEDEPTH', xtype=xtype, &
            long_name='lake depth', units='m')

       call mkpio_def_spatial_var(pioid, varname='PCT_WETLAND', xtype=xtype, &
            long_name='percent wetland', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='PCT_LAKE', xtype=xtype, &
            long_name='percent lake', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='PCT_GLACIER', xtype=xtype, &
            long_name='percent glacier', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='GLACIER_REGION', xtype=PIO_INT, &
            long_name='glacier region ID', units='unitless')

       call mkpio_defvar(pioid, varname='GLC_MEC', xtype=xtype, &
            dim1name='nglcecp1', long_name='Glacier elevation class', units='m')

       call mkpio_def_spatial_var(pioid, varname='PCT_GLC_MEC', xtype=xtype, &
            lev1name='nglcec', &
            long_name='percent glacier for each glacier elevation class (% of landunit)', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='TOPO_GLC_MEC', xtype=xtype, &
            lev1name='nglcec', &
            long_name='mean elevation on glacier elevation classes', units='m')

       if ( outnc_3dglc ) then
          call mkpio_def_spatial_var(pioid, varname='PCT_GLC_MEC_GIC', xtype=xtype, &
               lev1name='nglcec', &
               long_name='percent smaller glaciers and ice caps for each glacier elevation class (% of landunit)', units='unitless')

          call mkpio_def_spatial_var(pioid, varname='PCT_GLC_MEC_ICESHEET', xtype=xtype, &
               lev1name='nglcec', &
               long_name='percent ice sheet for each glacier elevation class (% of landunit)', units='unitless')

          call mkpio_def_spatial_var(pioid, varname='PCT_GLC_GIC', xtype=xtype, &
               long_name='percent ice caps/glaciers (% of landunit)', units='unitless')

          call mkpio_def_spatial_var(pioid, varname='PCT_GLC_ICESHEET', xtype=xtype, &
               long_name='percent ice sheet (% of landunit)', units='unitless')

       end if

       if ( outnc_3dglc ) then
          call mkpio_def_spatial_var(pioid, varname='PCT_GLC_MEC_GIC', xtype=xtype, &
               lev1name='nglcec', &
               long_name='percent smaller glaciers and ice caps for each glacier elevation class (% of landunit)', units='unitless')

          call mkpio_def_spatial_var(pioid, varname='PCT_GLC_MEC_ICESHEET', xtype=xtype, &
               lev1name='nglcec', &
               long_name='percent ice sheet for each glacier elevation class (% of landunit)', units='unitless')

          call mkpio_def_spatial_var(pioid, varname='PCT_GLC_GIC', xtype=xtype, &
               long_name='percent ice caps/glaciers (% of landunit)', units='unitless')

          call mkpio_def_spatial_var(pioid, varname='PCT_GLC_ICESHEET', xtype=xtype, &
               long_name='percent ice sheet (% of landunit)', units='unitless')

       end if

       call mkpio_def_spatial_var(pioid, varname='PCT_URBAN', xtype=xtype, &
            lev1name='numurbl', &
            long_name='percent urban for each density type', units='unitless')

       call mkpio_def_spatial_var(pioid, varname='URBAN_REGION_ID', xtype=PIO_INT, &
            long_name='urban region ID', units='unitless')

       ! call harvdata%getFieldsIdx( ind1D, ind2D )
       ! do j = 1, harvdata%num1Dfields()
       !    call mkpio_def_spatial_var(pioid, varname=mkharvest_fieldname(ind1D(j),constant=.true.), xtype=xtype, &
       !         long_name=mkharvest_longname(ind1D(j)), units=mkharvest_units(ind1D(j)) )
       ! end do
       ! do j = 1, harvdata%num2Dfields()
       !    call mkpio_def_spatial_var(pioid, varname=mkharvest_fieldname(ind2D(j),constant=.true.), xtype=xtype, &
       !         lev1name=harvdata%getFieldsDim(ind2D(j)), &
       !         long_name=mkharvest_longname(ind2D(j)), units=mkharvest_units(ind2D(j)) )
       ! end do
       ! deallocate(ind1D, ind2D)

    else

       ! call harvdata%getFieldsIdx( ind1D, ind2D )
       ! do j = 1, harvdata%num1Dfields()
       !    call mkpio_def_spatial_var(pioid, varname=mkharvest_fieldname(ind1D(j),constant=.false.), xtype=xtype, &
       !         lev1name='time', &
       !         long_name=mkharvest_longname(ind1D(j)), units=mkharvest_units(ind1D(j)) )
       ! end do
       ! do j = 1, harvdata%num2Dfields()
       !    call mkpio_def_spatial_var(pioid, varname=mkharvest_fieldname(ind2D(j),constant=.false.), xtype=xtype, &
       !         lev1name=harvdata%getFieldsDim(ind2D(j)), lev2name="time", &
       !         long_name=mkharvest_longname(ind2D(j)), units=mkharvest_units(ind2D(j)) )
       ! end do
       ! deallocate(ind1D, ind2D)

    end if  ! .not. dynlanduse

    ! End of define mode

    call mkpio_enddef(pioid)

    call pio_closefile(pioid)

  end subroutine mkfile

end module mkfileMod
