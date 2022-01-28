module mkfileMod

  use ESMF
  use pio
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_sys_mod  , only : shr_sys_getenv, shr_sys_abort
  use mkutilsMod   , only : get_filename, chkerr
  use mkvarpar     , only : nlevsoi, numrad, numstdpft
  use mkurbanparMod, only : numurbl, nlevurb
#ifdef TODO
  use mkglcmecMod  , only : nglcec
  use mkpftMod     , only : mkpftAtt
  use mksoilMod    , only : mksoilAtt
  use mkharvestMod , only : mkharvest_fieldname, mkharvest_numtypes, mkharvest_longname
  use mkharvestMod , only : mkharvest_units, harvestDataType
#endif
  use mkpioMod ! TODO: add only
  use mkvarctl

  implicit none
  private

  interface mkfile_output
     module procedure mkfile_output_int1d
     module procedure mkfile_output_int2d
     module procedure mkfile_output_real1d
     module procedure mkfile_output_real2d
  end interface mkfile_output

  public :: mkfile_fsurdat

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mkfile_fsurdat(nx, ny, mesh_o, dynlanduse, &
       pctlak, pctwet, lakedepth, organic, urban_classes, urban_region, soil_color, nsoilcol)

    ! input/output variables
    integer          , intent(in) :: nx
    integer          , intent(in) :: ny
    logical          , intent(in) :: dynlanduse
    type(ESMF_Mesh)  , intent(in) :: mesh_o
    real(r8), pointer, intent(in) :: pctlak(:)               ! percent of grid cell that is lake
    real(r8), pointer, intent(in) :: pctwet(:)               ! percent of grid cell that is wetland
    real(r8), pointer, intent(in) :: lakedepth(:)            ! lake depth (m)
    real(r8), pointer, intent(in) :: organic(:,:)            ! organic
    real(r8), pointer, intent(in) :: urban_classes(:,:)      ! percent cover of each urban class, as % of total urban area
    integer , pointer, intent(in) :: urban_region(:)         ! urban region ID
    integer , pointer, intent(in) :: soil_color(:)
    integer          , intent(in) :: nsoilcol 
#ifdef TODO
    type(harvestDataType) , intent(in) :: harvdata
#endif

    ! local variables
    type(file_desc_t)    :: pioid
    character(len=256)   :: varname
    character(len=256)   :: longname
    character(len=256)   :: units
    integer              :: xtype              ! external type
    integer, allocatable :: ind1D(:)           ! Indices of 1D harvest variables
    integer, allocatable :: ind2D(:)           ! Indices of 2D harvest variables
    integer              :: rcode
    integer              :: rc
    integer              :: n, i
    logical              :: define_mode
    character(len=256)   :: lev1name
    character(len=*), parameter :: subname=' (mkfile_fsurdat) '
    !-----------------------------------------------------------------------

    !---------------------------
    ! Create and open file
    !---------------------------

    ! TODO: what about setting no fill values?
     call mkpio_wopen(trim(fsurdat), clobber=.true., pioid=pioid)

    ! ----------------------------------------------------------------------
    ! Define dimensions and global attributes
    ! ----------------------------------------------------------------------

     call mkfile_define_dims(pioid, nx, ny, dynlanduse)
     call mkfile_define_atts(pioid, dynlanduse)

    ! ----------------------------------------------------------------------
    ! Define and outut variables
    ! ----------------------------------------------------------------------

    call ESMF_LogWrite(subname//'defining variables', ESMF_LOGMSG_INFO)

    if ( outnc_double ) then
       xtype = PIO_DOUBLE
    else
       xtype = PIO_REAL
    end if

    do n = 1,2

       define_mode = (n == 1)
       if (.not. define_mode) then
          rcode = pio_enddef(pioid)
       end if

       if (.not. dynlanduse) then

          if (root_task)  write(ndiag, '(a)') trim(subname)//"writing out pct_lake"
          call mkfile_output(pioid, define_mode, mesh_o, xtype, 'PCT_LAKE', 'percent_lake', &
               'unitless', pctlak, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          if (root_task)  write(ndiag, '(a)') trim(subname)//"writing out pct_wetland"
          call mkfile_output(pioid,define_mode, mesh_o, xtype, 'PCT_WETLAND', 'percent_wetland', &
               'unitless', pctwet, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          if (root_task)  write(ndiag, '(a)') trim(subname)//"writing out lakedepth"
          call mkfile_output(pioid, define_mode, mesh_o, xtype, 'LAKEDEPTH', 'lake_depth', &
               'm', lakedepth, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          if (root_task)  write(ndiag, '(a)') trim(subname)//"writing out organic"
          call mkfile_output(pioid, define_mode, mesh_o, xtype, 'ORGANIC', 'organic matter density at soil levels', &
               'kg/m3 (assumed carbon content 0.58 gC per gOM)', organic, lev1name='nlevsoi', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          if (root_task)  write(ndiag, '(a)') trim(subname)//"writing out soil color"
          call mkfile_output(pioid, define_mode, mesh_o, 'SOIL_COLOR', 'soil color', &
               'unitless', soil_color,  rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          if (root_task)  write(ndiag, '(a)') trim(subname)//"writing out urban region id"
          call mkfile_output(pioid, define_mode, mesh_o, 'URBAN_REGION_ID', 'urban region ID', &
               'unitless', urban_region, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          if (root_task)  write(ndiag, '(a)') trim(subname)//"writing out percnt urban"
          call mkfile_output(pioid, define_mode, mesh_o, xtype, 'PCT_URBAN', 'percent urban for each density type', &
               'unitless', urban_classes, lev1name='numurbl', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! ----------------------------------------------------------------------
          ! Make Urban Parameters from raw input data and write to surface dataset 
          ! Write to netcdf file is done inside mkurbanpar routine
          ! ----------------------------------------------------------------------

          ! call mkurbanpar(datfname=mksrf_furban, ncido=ncid, region_o=urban_region, &
          !      urbn_classes_gcell_o=urbn_classes_g, &
          !      urban_skip_abort_on_invalid_data_check=urban_skip_abort_on_invalid_data_check)

       end if
    end do

    ! Close surface dataset
    call pio_closefile(pioid)

    if (root_task) then
       write (ndiag,'(72a1)') ("-",i=1,60)
       write (ndiag,'(a)')' land model surface data set successfully created for '
    end if

  end subroutine mkfile_fsurdat

!=================================================================================
  subroutine mkfile_define_dims(pioid, nx, ny, dynlanduse)

    ! Define dimensions.

    ! input/output variables
    type(file_desc_t) , intent(in) :: pioid
    integer           , intent(in) :: nx, ny
    logical           , intent(in) :: dynlanduse

    ! local variables
    integer :: dimid              ! temporary
    integer :: rcode
    character(len=*), parameter :: subname = 'mkfile_define_dims'
    !-----------------------------------------------------------------------

    call ESMF_LogWrite(subname//' defining dimensions', ESMF_LOGMSG_INFO)

    if (outnc_1d) then
       rcode = pio_def_dim(pioid, 'gridcell', nx, dimid)
    else
       rcode = pio_def_dim(pioid, 'lsmlon', nx, dimid)
       rcode = pio_def_dim(pioid, 'lsmlat', ny, dimid)
    end if
    rcode = pio_def_dim(pioid, 'nlevsoi', nlevsoi       , dimid)
    rcode = pio_def_dim(pioid, 'nlevurb', nlevurb       , dimid)
    rcode = pio_def_dim(pioid, 'numurbl', numurbl       , dimid)

#ifdef TODO
    rcode = pio_def_dim(pioid, 'time'   , PIO_UNLIMITED , dimid)
    rcode = pio_def_dim(pioid, 'nchar'  , 256           , dimid)
    rcode = pio_def_dim(pioid, 'numrad' , numrad        , dimid)
    if (.not. dynlanduse) then
       rcode = pio_def_dim(pioid, 'nglcec', nglcec, dimid)
       rcode = pio_def_dim(pioid, 'nglcecp1', nglcec+1, dimid)
    end if
#endif

  end subroutine mkfile_define_dims

  !=================================================================================
  subroutine mkfile_define_atts(pioid, dynlanduse)

    ! input/output variables
    type(file_desc_t) , intent(in) :: pioid
    logical           , intent(in) :: dynlanduse

    ! local variables
    integer              :: values(8)          ! temporary
    character(len=256)   :: str                ! global attribute string
    character(len=256)   :: name               ! name of attribute
    character(len=256)   :: unit               ! units of attribute
    character(len= 18)   :: datetime           ! temporary
    character(len=  8)   :: date               ! temporary
    character(len= 10)   :: time               ! temporary
    character(len=  5)   :: zone               ! temporary
    integer              :: rcode
    character(len=*), parameter :: subname = 'mkfile_define_atts'
    !-----------------------------------------------------------------------

    !---------------------------
    ! Set global attributes.
    !---------------------------

    call ESMF_LogWrite(subname//'setting global attributes', ESMF_LOGMSG_INFO)

    str = 'NCAR-CSM'
    rcode = pio_put_att(pioid, pio_global, "Conventions", trim(str))

    call date_and_time (date, time, zone, values)
    datetime(1:8) =        date(5:6) // '-' // date(7:8) // '-' // date(3:4)
    datetime(9:)  = ' ' // time(1:2) // ':' // time(3:4) // ':' // time(5:6) // ' '
    str = 'created on: ' // datetime
    rcode = pio_put_att (pioid, pio_global, 'History_Log', trim(str))

#ifdef TODO
    call shr_sys_getenv ('LOGNAME', str, ier)
    rcode = pio_put_att (pioid, pio_global, 'Logname', trim(str))

    call shr_sys_getenv ('HOST', str, ier)
    rcode = pio_put_att (pioid, pio_global, 'Host', trim(str))
#endif

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

  end subroutine mkfile_define_atts

  !=================================================================================
  subroutine mkfile_output_int1d(pioid, define_mode, mesh, varname, longname, units, ipointer, rc)

    ! input/output variables
    type(file_desc_t) , intent(inout)  :: pioid
    logical           , intent(in)  :: define_mode
    type(ESMF_Mesh)   , intent(in)  :: mesh
    character(len=*)  , intent(in)  :: varname
    character(len=*)  , intent(in)  :: longname
    character(len=*)  , intent(in)  :: units
    integer           , pointer     :: ipointer(:)
    integer           , intent(out) :: rc

    ! local variables
    type(io_desc_t)      :: pio_iodesc
    type(var_desc_t)     :: pio_varid
    integer              :: rcode
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (define_mode) then
       call mkpio_def_spatial_var(pioid, trim(varname), PIO_INT, trim(longname), trim(units))
    else
       call mkpio_iodesc_output(pioid, mesh, trim(varname), pio_iodesc, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
       rcode = pio_inq_varid(pioid, trim(varname), pio_varid)
       call pio_write_darray(pioid, pio_varid, pio_iodesc, ipointer, rcode)
       call pio_freedecomp(pioid, pio_iodesc)
    end if
  end subroutine mkfile_output_int1d

  !=================================================================================
  subroutine mkfile_output_int2d(pioid, define_mode, mesh, varname, longname, units, ipointer, rc)

    ! input/output variables
    type(file_desc_t) , intent(inout)  :: pioid
    logical           , intent(in)  :: define_mode
    type(ESMF_Mesh)   , intent(in)  :: mesh
    character(len=*)  , intent(in)  :: varname
    character(len=*)  , intent(in)  :: longname
    character(len=*)  , intent(in)  :: units
    integer           , pointer     :: ipointer(:,:)
    integer           , intent(out) :: rc

    ! local variables
    type(io_desc_t)      :: pio_iodesc
    type(var_desc_t)     :: pio_varid
    integer              :: rcode
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (define_mode) then
       call mkpio_def_spatial_var(pioid, trim(varname), PIO_INT, trim(longname), trim(units))
    else
       call mkpio_iodesc_output(pioid, mesh, trim(varname), pio_iodesc, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
       rcode = pio_inq_varid(pioid, trim(varname), pio_varid)
       call pio_write_darray(pioid, pio_varid, pio_iodesc, ipointer, rcode)
       call pio_freedecomp(pioid, pio_iodesc)
    end if
  end subroutine mkfile_output_int2d

  !=================================================================================
  subroutine mkfile_output_real1d(pioid, define_mode, mesh, xtype, varname, longname, units, rpointer, rc)

    ! input/output variables
    type(file_desc_t), intent(inout) :: pioid
    logical         , intent(in)  :: define_mode
    type(ESMF_Mesh) , intent(in) :: mesh
    integer         , intent(in)  :: xtype
    character(len=*), intent(in)  :: varname
    character(len=*), intent(in)  :: longname
    character(len=*), intent(in)  :: units
    real(r8)        , pointer     :: rpointer(:)
    integer         , intent(out) :: rc

    ! local variables
    type(io_desc_t)      :: pio_iodesc
    type(var_desc_t)     :: pio_varid
    integer              :: rcode
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (define_mode) then
       call mkpio_def_spatial_var(pioid, trim(varname), xtype, trim(longname), trim(units))
    else
       call mkpio_iodesc_output(pioid, mesh, trim(varname), pio_iodesc, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
       rcode = pio_inq_varid(pioid, trim(varname), pio_varid)
       call pio_write_darray(pioid, pio_varid, pio_iodesc, rpointer, rcode)
       call pio_freedecomp(pioid, pio_iodesc)
    end if
  end subroutine mkfile_output_real1d

  !=================================================================================
  subroutine mkfile_output_real2d(pioid, define_mode, mesh, xtype, varname, longname, units, rpointer, lev1name, rc)

    ! input/output variables
    type(file_desc_t), intent(inout) :: pioid
    logical         , intent(in) :: define_mode
    type(ESMF_Mesh) , intent(in) :: mesh
    integer         , intent(in) :: xtype
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: units
    character(len=*), intent(in) :: lev1name
    real(r8)        , pointer    :: rpointer(:,:)
    integer         , intent(out) :: rc

    ! local variables
    type(io_desc_t)      :: pio_iodesc
    type(var_desc_t)     :: pio_varid
    integer              :: rcode
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (define_mode) then
       call mkpio_def_spatial_var(pioid, trim(varname), xtype, trim(lev1name), trim(longname), trim(units))
    else
       call mkpio_iodesc_output(pioid, mesh, trim(varname), pio_iodesc, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
       rcode = pio_inq_varid(pioid, trim(varname), pio_varid)
       call pio_write_darray(pioid, pio_varid, pio_iodesc, rpointer, rcode)
       call pio_freedecomp(pioid, pio_iodesc)
    end if
  end subroutine mkfile_output_real2d

end module mkfileMod
