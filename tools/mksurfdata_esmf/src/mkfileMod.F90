module mkfileMod

  use ESMF
  use pio
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_sys_mod       , only : shr_sys_getenv, shr_sys_abort
  use mkutilsMod        , only : get_filename, chkerr
  use mkvarpar          , only : nlevsoi, numrad, numstdpft
  use mkurbanparMod     , only : numurbl, nlevurb, mkurbanpar
  use mkglcmecMod       , only : nglcec
  use mklaiMod          , only : mklai          
  use mkpftConstantsMod , only : natpft_lb, natpft_ub, cft_lb, cft_ub, num_cft, num_natpft
  use mkpctPftTypeMod   , only : pct_pft_type, get_pct_p2l_array, get_pct_l2g_array, update_max_array
  !use mkharvestMod      , only : mkharvest_fieldname, mkharvest_numtypes, mkharvest_longname, mkharvest_units, harvestDataType
  use mkpioMod ! TODO: add only
  use mkvarctl

  implicit none
  private

  public :: mkfile_fdyndat
  public :: mkfile_define_dims
  public :: mkfile_define_atts
  public :: mkfile_define_vars
  public :: mkfile_output

  interface mkfile_output
     module procedure mkfile_output_int1d
     module procedure mkfile_output_int2d
     module procedure mkfile_output_real1d
     module procedure mkfile_output_real2d
  end interface mkfile_output

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mkfile_fdyndat(nx, ny, mesh_o, ns_o, lon, lat, rc) 

    ! input/output variables
    integer         , intent(in) :: nx
    integer         , intent(in) :: ny
    type(ESMF_Mesh) , intent(in) :: mesh_o
    integer         , intent(in) :: ns_o
    real(r8)        , intent(in) :: lon(:) 
    real(r8)        , intent(in) :: lat(:) 
    integer         , intent(out) :: rc

    ! local variables
    type(file_desc_t)     :: pioid_o
    type(var_desc_t)      :: pio_varid_o
    type(file_desc_t)     :: pioid_i
    type(var_desc_t)      :: pio_varid_i
    real(r8), allocatable :: pctlnd_pft_dyn(:) ! PFT data: % of gridcell for dyn landuse PFTs
    character(len=256)    :: varname
    character(len=256)    :: longname
    character(len=256)    :: units
    integer               :: xtype             ! external type
    integer               :: dimid
    logical               :: define_mode
    character(len=256)    :: lev1name
    integer, allocatable  :: ind1D(:)          ! Indices of 1D harvest variables
    integer, allocatable  :: ind2D(:)          ! Indices of 2D harvest variables
    integer               :: n, i              ! Indices
    integer               :: rcode             ! Error status
    character(len=*), parameter :: subname=' (mkfile_fdyndat) '
    ! ------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! if (root_task) then
    !    write(ndiag,'(a)')'creating dynamic land use dataset'
    ! end if

    ! allocate(pctlnd_pft_dyn(ns_o))
    ! !call mkharvest_init( lsize_o, spval, harvdata, mksrf_fhrvtyp )

    ! ! open output file
    ! call mkpio_wopen(trim(fdyndat), clobber=.true., pioid=pioid_o)

    ! ! Define dimensions and global attributes
    ! call mkfile_define_dims(pioid, nx, ny, dynlanduse=.true.)
    ! call mkfile_define_atts(pioid, dynlanduse = .true.)

    ! if ( outnc_double ) then
    !    xtype = PIO_DOUBLE
    ! else
    !    xtype = PIO_REAL
    ! end if

    ! do n = 1,2

    !    define_mode = (n == 1)
    !    if (.not. define_mode) then
    !       rcode = pio_enddef(pioid_o)
    !    end if

    !    ! --------------------------------
    !    ! Write out model grid
    !    ! --------------------------------
    !    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out model grid"

    !    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out LONGXY"
    !    call mkfile_output(pioid_o,  mesh_o,  'LONGXY', 'model longitudes', &
    !         'degrees', lat, rc=rc) 
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out LATIXY"
    !    call mkfile_output(pioid_o,  mesh_o,  'LATIXY', 'model latitudes', &
    !         'degrees', lon, rc=rc) 
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !    !rcode = pio_open(trim(fdyndat), nf_write, pioid_o))
    !    !rcode = pio_set_fill (pioid_o, nf_nofill, omode))

    !    ! --------------------------------
    !    ! Write out non-spatial variables
    !    ! --------------------------------
    !    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out non-spatial variables"

    !    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out natpft"
    !    if (define_mode) then
    !       rcode = pio_def_var(pioid_o, 'natpft', PIO_INT, pio_varid)
    !    else
    !       rcode = pio_inq_varid(pioid_o, 'natpft', pio_varid)
    !       rcode = pio_put_var_int(pioid_o, pio_varid, (/(n,n=natpft_lb,natpft_ub)/))
    !    end if

    !    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out cft"
    !    if (define_mode) then
    !       rcode = pio_def_var(pioid_o, 'cft', PIO_INT, pio_varid)
    !    else
    !       rcode = pio_inq_varid(pioid_o, 'cft', pio_varid)
    !       rcode = pio_put_var_int(pioid_o, pio_varid, (/(n,n=cft_lb,cft_ub)/))
    !    end if

    !    ! --------------------------------
    !    ! Write out spatial variables
    !    ! --------------------------------

    !    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing land mask from pft dataset"
    !    call mkfile_output(pioid_o,  mesh_o, 'PFTDATA_MASK', &
    !         'land mask from pft dataset, indicative of real/fake points', 'unitless', pftdata_mask, rc=rc)
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out soil depth"
    !    call mkfile_output(pioid_o,  mesh_o,  'LANDFRAC_PFT', &
    !         'land fraction from pft dataset', 'unitless', landfrac_pft, rc=rc)
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! end do

    ! ! -----------------------------------------
    ! ! Read in each dynamic pft landuse dataset
    ! ! -----------------------------------------

    ! ! Open txt file
    ! open (newunit=nfdyn, file=trim(mksrf_fdynuse), form='formatted', iostat=ier)
    ! if (ioe /= 0) then
    !    call shr_sys_abort(subname//" failed to open file "//trim(mksrf_fdynuse))
    ! end if

    ! pctnatpft_max = pctnatpft
    ! pctcft_max = pctcft

    ! ntim = 0
    ! do

    !    ! Determine file name
    !    read(nfdyn, '(A195,1x,I4)', iostat=ier) string, year
    !    if (ier /= 0) exit

    !    ! If pft fraction override is set, than intrepret string as PFT and harvesting override values
    !    ! Otherwise intrepret string as a filename with PFT and harvesting values in it
    !    if ( all_veg )then
    !       fname = ' '
    !       fhrvname  = ' '
    !       call mkpft_parse_oride(string)
    !       call mkharvest_parse_oride(string)
    !       write(6, '(a, i4, a)') 'PFT and harvesting values for year ', year, ' :'
    !       write(6, '(a, a)') '    ', trim(string)
    !    else
    !       fname = string
    !       read(nfdyn, '(A195,1x,I4)', iostat=ier) fhrvname, year2
    !       if (root_task) then
    !          write(ndiat,'(a,i8,a)')' input pft dynamic dataset for year ', year,' is : '//trim(fname)
    !       end if
    !       if ( year2 /= year ) then
    !          write(6,*) subname, ' error: year for harvest not equal to year for PFT files'
    !          call shr_sys_abort()
    !       end if
    !    end if
    !    ntim = ntim + 1

    !    ! Create pctpft data at model resolution from file fname
    !    call mkpft( mksrf_fvegtyp_mesh, fname, mesh_model, &
    !         pctlnd_o=pctlnd_pft, pctnatpft_o=pctnatpft, pctcft_o=pctcft, rc=rc)
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkdomain')

    !    ! Create harvesting data at model resolution
    !    call mkharvest(  mksrf_fhrvtyp, mesh_model, harvdata=harvdata, rc )
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkdomain')

    !    ! Consistency check on input land fraction
    !    do n = 1,lsize_o
    !       if (pctlnd_pft_dyn(n) /= pctlnd_pft(n)) then
    !          write(6,*) subname,' error: pctlnd_pft for dynamics data = ',&
    !               pctlnd_pft_dyn(n), ' not equal to pctlnd_pft for surface data = ',&
    !               pctlnd_pft(n),' at n= ',n
    !          if ( trim(fname) == ' ' )then
    !             write(6,*) ' PFT string = ', string
    !          else
    !             write(6,*) ' PFT file = ', fname
    !          end if
    !          call shr_sys_abort()
    !       end if
    !    end do

    !    call change_landuse( dynpft=.true.)

    !    call normalizencheck_landuse(ldomain)

    !    call update_max_array(pctnatpft_max,pctnatpft)
    !    call update_max_array(pctcft_max,pctcft)

    !    ! Output time-varying data for current year

    !    rcode = pio_inq_varid(pioid, 'PCT_NAT_PFT', pio_varid)
    !    call mkpio_put_time_slice(pioid, pio_varid, ntim, get_pct_p2l_array(pctnatpft))

    !    rcode = pio_inq_varid(pioid, 'PCT_CROP', pio_varid)
    !    call mkpio_put_time_slice(pioid, pio_varid, ntim, get_pct_l2g_array(pctcft))

    !    if (num_cft > 0) then
    !       rcode = pio_inq_varid(pioid, 'PCT_CFT', pio_varid)
    !       call mkpio_put_time_slice(pioid, pio_varid, ntim, get_pct_p2l_array(pctcft))
    !    end if

    !    call harvdata%getFieldsIdx( harvind1D, harvind2D)
    !    do k = 1, harvdata%num1Dfields()
    !       rcode = pio_inq_varid(pioid, trim(mkharvest_fieldname(harvind1D(k),constant=.false.)), pio_varid)
    !       harvest1D => harvdata%get1DFieldPtr( harvind1D(k), output=.true. )
    !       call mkpio_put_time_slice(pioid, pio_varid, ntim, harvest1D)
    !    end do
    !    do k = 1, harvdata%num2Dfields()
    !       rcode = pio_inq_varid(pioid, trim(mkharvest_fieldname(harvind2D(k),constant=.false.)), pio_varid)
    !       harvest2D => harvdata%get2DFieldPtr( harvind2D(k), output=.true. )
    !       call mkpio_put_time_slice(pioid, pio_varid, ntim, harvest2D)
    !    end do
    !    deallocate( harvind1D, harvind2D )

    !    rcode = pio_inq_varid(pioid, 'YEAR', pio_varid)
    !    rcode = pio_put_vara_int(pioid, pio_varid, ntim, 1, year)

    !    rcode = pio_inq_varid(pioid, 'time', pio_varid)
    !    rcode = pio_put_vara_int(pioid, pio_varid, ntim, 1, year)

    !    rcode = pio_inq_varid(pioid, 'input_pftdata_filename', pio_varid)
    !    rcode = pio_put_vara_text(pioid, pio_varid, (/ 1, ntim /), (/ len_trim(string), 1 /), trim(string))

    ! end do   ! end of read loop

    ! rcode = pio_inq_varid(pioid, 'PCT_NAT_PFT_MAX', pio_varid)
    ! rcode = pio_put_var_double(pioid, pio_varid, get_pct_p2l_array(pctnatpft_max))

    ! rcode = pio_inq_varid(pioid, 'PCT_CROP_MAX', pio_varid)
    ! rcode = pio_put_var_double(pioid, pio_varid, get_pct_l2g_array(pctcft_max))

    ! if (num_cft > 0) then
    !    rcode = pio_inq_varid(pioid, 'PCT_CFT_MAX', pio_varid)
    !    rcode = pio_put_var_double(pioid, pio_varid, get_pct_p2l_array(pctcft_max))
    ! end if

    ! rcode = pio_closefile(pioid)

  end subroutine mkfile_fdyndat

  !=================================================================================
  subroutine mkfile_define_dims(pioid, nx, ny, dynlanduse)
    !
    ! Define dimensions.
    !
    ! input/output variables
    type(file_desc_t) , intent(in) :: pioid
    integer           , intent(in) :: nx, ny
    logical           , intent(in) :: dynlanduse

    ! local variables
    integer :: dimid  ! temporary
    integer :: rcode  ! error status
    integer :: pftsize
    integer :: natpftsize
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
    rcode = pio_def_dim(pioid, 'numrad' , numrad, dimid)
    if (.not. dynlanduse) then
       rcode = pio_def_dim(pioid, 'nglcec'  , nglcec  , dimid)
       rcode = pio_def_dim(pioid, 'nglcecp1', nglcec+1, dimid)
    end if
    rcode = pio_def_dim(pioid, 'time'   , PIO_UNLIMITED , dimid)
    rcode = pio_def_dim(pioid, 'nchar'  , 256           , dimid)

    if (.not. dynlanduse) then
       pftsize = numpft + 1
       rcode = pio_def_dim(pioid, 'lsmpft' , pftsize, dimid)
    end if
    natpftsize = num_natpft + 1
    rcode = pio_def_dim (pioid, 'natpft', natpftsize, dimid)

    ! zero-size dimensions can cause problems, so we only include the
    ! cft dimension if num_cft > 0 Note that this implies that we can
    ! only include PCT_CFT on the dataset if num_cft > 0
    if (num_cft > 0) then
       rcode = pio_def_dim (pioid, 'cft', num_cft, dimid)
    end if

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

    if (.not. dynlanduse) then
       if ( soil_color_override /= unsetcol )then
          rcode = pio_put_att(pioid, pio_global, 'soil_color_override', 'TRUE')
       else
          str = get_filename(mksrf_fsoicol)
          rcode = pio_put_att(pioid, pio_global, 'soil_color_raw_data_file_name', trim(str))
       end if
       if ( soil_clay_override /= unsetsoil .and. soil_sand_override /= unsetsoil) then
          rcode = pio_put_att(pioid, pio_global, 'soil_clay_override', 'TRUE')
       else
          str = get_filename(mksrf_fsoitex)
          rcode = pio_put_att(pioid, pio_global, 'soil_texture_raw_data_file_name', trim(str))
       end if
       if ( soil_fmax_override /= unsetsoil )then
          rcode = pio_put_att(pioid, pio_global, 'soil_fmax_override', 'TRUE')
       else
          str = get_filename(mksrf_fmax)
          rcode = pio_put_att(pioid, pio_global, 'fmax_raw_data_file_name', trim(str))
       end if
       str = get_filename(mksrf_forganic)
       rcode = pio_put_att(pioid, pio_global, 'organic_matter_raw_data_file_name', trim(str))
    end if

  end subroutine mkfile_define_atts

  !=================================================================================
  subroutine mkfile_define_vars(pioid, dynlanduse)

    ! Define fsurdat variables

    ! input/output variables
    type(file_desc_t) , intent(in) :: pioid
    logical           , intent(in) :: dynlanduse

    ! local variables
    integer :: xtype    ! external type
    character(len=*), parameter :: subname = 'mkfile_define_vars'
    !-----------------------------------------------------------------------

    if ( outnc_double ) then
       xtype = PIO_DOUBLE
    else
       xtype = PIO_REAL
    end if

    call mkpio_def_spatial_var(pioid=pioid, varname='LONGXY', xtype=xtype, &
         long_name='longitude', units='degrees east')

    call mkpio_def_spatial_var(pioid=pioid, varname='LATIXY', xtype=xtype, &
         long_name='latitude', units='degrees north')

    if (.not. dynlanduse) then

       call mkpio_defvar(pioid=pioid, varname='mxsoil_color', xtype=PIO_INT, &
          long_name='maximum numbers of soil colors', units='unitless')
     
       call mkpio_def_spatial_var(pioid=pioid, varname='SOIL_COLOR', xtype=PIO_INT, &
            long_name='soil color', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='PCT_SAND', xtype=xtype, &
            lev1name='nlevsoi', &
            long_name='percent sand', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='PCT_CLAY', xtype=xtype, &
            lev1name='nlevsoi', &
            long_name='percent clay', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='ORGANIC', xtype=xtype, &
            lev1name='nlevsoi', &
            long_name='organic matter density at soil levels', &
            units='kg/m3 (assumed carbon content 0.58 gC per gOM)')

       call mkpio_def_spatial_var(pioid=pioid, varname='FMAX', xtype=xtype, &
            long_name='maximum fractional saturated area', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='EF1_BTR', xtype=xtype, &
            long_name='EF btr (isoprene)', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='EF1_FET', xtype=xtype, &
            long_name='EF fet (isoprene)', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='EF1_FDT', xtype=xtype, &
            long_name='EF fdt (isoprene)', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='EF1_SHR', xtype=xtype, &
            long_name='EF shr (isoprene)', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='EF1_GRS', xtype=xtype, &
            long_name='EF grs (isoprene)', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='EF1_CRP', xtype=xtype, &
            long_name='EF crp (isoprene)', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='CANYON_HWR', xtype=xtype, &
            lev1name='numurbl', &
            long_name='canyon height to width ratio', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='EM_IMPROAD', xtype=xtype, &
            lev1name='numurbl', &
            long_name='emissivity of impervious road', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='EM_PERROAD', xtype=xtype, &
            lev1name='numurbl', &
            long_name='emissivity of pervious road', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='EM_ROOF', xtype=xtype, &
            lev1name='numurbl', &
            long_name='emissivity of roof', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='EM_WALL', xtype=xtype, &
            lev1name='numurbl', &
            long_name='emissivity of wall', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='HT_ROOF', xtype=xtype, &
            lev1name='numurbl', &
            long_name='height of roof', units='meters')

       call mkpio_def_spatial_var(pioid=pioid, varname='THICK_ROOF', xtype=xtype, &
            lev1name='numurbl', &
            long_name='thickness of roof', units='meters')

       call mkpio_def_spatial_var(pioid=pioid, varname='THICK_WALL', xtype=xtype, &
            lev1name='numurbl', &
            long_name='thickness of wall', units='meters')

       call mkpio_def_spatial_var(pioid=pioid, varname='T_BUILDING_MIN', xtype=xtype, &
            lev1name='numurbl', &
            long_name='minimum interior building temperature', units='K')

       call mkpio_def_spatial_var(pioid=pioid, varname='WIND_HGT_CANYON', xtype=xtype, &
            lev1name='numurbl', &
            long_name='height of wind in canyon', units='meters')

       call mkpio_def_spatial_var(pioid=pioid, varname='WTLUNIT_ROOF', xtype=xtype, &
            lev1name='numurbl', &
            long_name='fraction of roof', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='WTROAD_PERV', xtype=xtype, &
            lev1name='numurbl', &
            long_name='fraction of pervious road', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='ALB_IMPROAD_DIR', xtype=xtype, &
            lev1name='numurbl', lev2name='numrad', &
            long_name='direct albedo of impervious road', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='ALB_IMPROAD_DIF', xtype=xtype, &
            lev1name='numurbl', lev2name='numrad', &
            long_name='diffuse albedo of impervious road', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='ALB_PERROAD_DIR', xtype=xtype, &
            lev1name='numurbl', lev2name='numrad', &
            long_name='direct albedo of pervious road', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='ALB_PERROAD_DIF', xtype=xtype, &
            lev1name='numurbl', lev2name='numrad', &
            long_name='diffuse albedo of pervious road', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='ALB_ROOF_DIR', xtype=xtype, &
            lev1name='numurbl', lev2name='numrad', &
            long_name='direct albedo of roof', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='ALB_ROOF_DIF', xtype=xtype, &
            lev1name='numurbl', lev2name='numrad', &
            long_name='diffuse albedo of roof', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='ALB_WALL_DIR', xtype=xtype, &
            lev1name='numurbl', lev2name='numrad', &
            long_name='direct albedo of wall', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='ALB_WALL_DIF', xtype=xtype, &
            lev1name='numurbl', lev2name='numrad', &
            long_name='diffuse albedo of wall', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='TK_ROOF', xtype=xtype, &
            lev1name='numurbl', lev2name='nlevurb', &
            long_name='thermal conductivity of roof', units='W/m*K')

       call mkpio_def_spatial_var(pioid=pioid, varname='TK_WALL', xtype=xtype, &
            lev1name='numurbl', lev2name='nlevurb', &
            long_name='thermal conductivity of wall', units='W/m*K')

       call mkpio_def_spatial_var(pioid=pioid, varname='TK_IMPROAD', xtype=xtype, &
            lev1name='numurbl', lev2name='nlevurb', &
            long_name='thermal conductivity of impervious road', units='W/m*K')

       call mkpio_def_spatial_var(pioid=pioid, varname='CV_ROOF', xtype=xtype, &
            lev1name='numurbl', lev2name='nlevurb', &
            long_name='volumetric heat capacity of roof', units='J/m^3*K')

       call mkpio_def_spatial_var(pioid=pioid, varname='CV_WALL', xtype=xtype, &
            lev1name='numurbl', lev2name='nlevurb', &
            long_name='volumetric heat capacity of wall', units='J/m^3*K')

       call mkpio_def_spatial_var(pioid=pioid, varname='CV_IMPROAD', xtype=xtype, &
            lev1name='numurbl', lev2name='nlevurb', &
            long_name='volumetric heat capacity of impervious road', units='J/m^3*K')

       call mkpio_def_spatial_var(pioid=pioid, varname='NLEV_IMPROAD', xtype=PIO_INT, &
            lev1name='numurbl', &
            long_name='number of impervious road layers', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='peatf', xtype=xtype, &
            long_name='peatland fraction', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='zbedrock', xtype=xtype, &
            long_name='soil depth', units='m')

       call mkpio_def_spatial_var(pioid=pioid, varname='abm', xtype=PIO_INT, &
            long_name='agricultural fire peak month', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='gdp', xtype=xtype, &
            long_name='gdp', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='SLOPE', xtype=xtype, &
            long_name='mean topographic slope', units='degrees')

       call mkpio_def_spatial_var(pioid=pioid, varname='STD_ELEV', xtype=xtype, &
            long_name='standard deviation of elevation', units='m')

       if ( outnc_vic )then
          call mkpio_def_spatial_var(pioid=pioid, varname='binfl', xtype=xtype, &
               long_name='VIC b parameter for the Variable Infiltration Capacity Curve', units='unitless')

          call mkpio_def_spatial_var(pioid=pioid, varname='Ws', xtype=xtype, &
               long_name='VIC Ws parameter for the ARNO curve', units='unitless')

          call mkpio_def_spatial_var(pioid=pioid, varname='Dsmax', xtype=xtype, &
               long_name='VIC Dsmax parameter for the ARNO curve', units='mm/day')

          call mkpio_def_spatial_var(pioid=pioid, varname='Ds', xtype=xtype, &
               long_name='VIC Ds parameter for the ARNO curve', units='unitless')

       end if
       call mkpio_def_spatial_var(pioid=pioid, varname='LAKEDEPTH', xtype=xtype, &
            long_name='lake depth', units='m')

       call mkpio_def_spatial_var(pioid=pioid, varname='PCT_WETLAND', xtype=xtype, &
            long_name='percent wetland', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='PCT_LAKE', xtype=xtype, &
            long_name='percent lake', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='PCT_GLACIER', xtype=xtype, &
            long_name='percent glacier', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='GLACIER_REGION', xtype=PIO_INT, &
            long_name='glacier region ID', units='unitless')

       call mkpio_defvar(pioid=pioid, varname='GLC_MEC', xtype=xtype, &
            dim1name='nglcecp1', long_name='Glacier elevation class', units='m')

       call mkpio_def_spatial_var(pioid=pioid, varname='PCT_GLC_MEC', xtype=xtype, &
            lev1name='nglcec', &
            long_name='percent glacier for each glacier elevation class (% of landunit)', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='TOPO_GLC_MEC', xtype=xtype, &
            lev1name='nglcec', &
            long_name='mean elevation on glacier elevation classes', units='m')

       if ( outnc_3dglc ) then
          call mkpio_def_spatial_var(pioid=pioid, varname='PCT_GLC_MEC_GIC', xtype=xtype, &
               lev1name='nglcec', &
               long_name='percent smaller glaciers and ice caps for each glacier elevation class (% of landunit)', units='unitless')

          call mkpio_def_spatial_var(pioid=pioid, varname='PCT_GLC_MEC_ICESHEET', xtype=xtype, &
               lev1name='nglcec', &
               long_name='percent ice sheet for each glacier elevation class (% of landunit)', units='unitless')

          call mkpio_def_spatial_var(pioid=pioid, varname='PCT_GLC_GIC', xtype=xtype, &
               long_name='percent ice caps/glaciers (% of landunit)', units='unitless')

          call mkpio_def_spatial_var(pioid=pioid, varname='PCT_GLC_ICESHEET', xtype=xtype, &
               long_name='percent ice sheet (% of landunit)', units='unitless')

       end if

       if ( outnc_3dglc ) then
          call mkpio_def_spatial_var(pioid=pioid, varname='PCT_GLC_MEC_GIC', xtype=xtype, &
               lev1name='nglcec', &
               long_name='percent smaller glaciers and ice caps for each glacier elevation class (% of landunit)', &
               units='unitless')

          call mkpio_def_spatial_var(pioid=pioid, varname='PCT_GLC_MEC_ICESHEET', xtype=xtype, &
               lev1name='nglcec', &
               long_name='percent ice sheet for each glacier elevation class (% of landunit)', units='unitless')

          call mkpio_def_spatial_var(pioid=pioid, varname='PCT_GLC_GIC', xtype=xtype, &
               long_name='percent ice caps/glaciers (% of landunit)', units='unitless')

          call mkpio_def_spatial_var(pioid=pioid, varname='PCT_GLC_ICESHEET', xtype=xtype, &
               long_name='percent ice sheet (% of landunit)', units='unitless')

       end if

       call mkpio_def_spatial_var(pioid=pioid, varname='PCT_URBAN', xtype=xtype, &
            lev1name='numurbl', &
            long_name='percent urban for each density type', units='unitless')

       call mkpio_def_spatial_var(pioid=pioid, varname='URBAN_REGION_ID', xtype=PIO_INT, &
            long_name='urban region ID', units='unitless')
    end if

    if (.not. dynlanduse) then
       call mkpio_def_spatial_var(pioid=pioid, varname='CONST_HARVEST_VH1', xtype=xtype, &
            long_name = "harvest from primary forest", units = "gC/m2/yr")
       call mkpio_def_spatial_var(pioid=pioid, varname='CONST_HARVEST_VH2', xtype=xtype, &
            long_name = "harvest from primary non-forest", units = "gC/m2/yr")
       call mkpio_def_spatial_var(pioid=pioid, varname='CONST_HARVEST_SH1', xtype=xtype, &
            long_name = "harvest from secondary mature-forest", units = "gC/m2/yr")
       call mkpio_def_spatial_var(pioid=pioid, varname='CONST_HARVEST_SH2', xtype=xtype, &
            long_name = "harvest from secondary young-forest", units = "gC/m2/yr")
       call mkpio_def_spatial_var(pioid=pioid, varname='CONST_HARVEST_SH3', xtype=xtype, &
            long_name = "harvest from secondary non-forest", units = "gC/m2/yr")
       call mkpio_def_spatial_var(pioid=pioid, varname='CONST_GRAZING', xtype=xtype, &
            long_name = "grazing of herbacous pfts", units = "gC/m2/yr")
       call mkpio_def_spatial_var(pioid=pioid, varname='CONST_FERTNITRO_CFT', xtype=xtype, &
            lev1name = 'cft', &
            long_name = "nitrogen fertilizer for each crop", units = "gN/m2/yr")
       call mkpio_def_spatial_var(pioid=pioid, varname='UNREPRESENTED_PFT_LULCC', xtype=xtype,&
            lev1name = 'natpft', &
            long_name = "unrepresented PFT gross LULCC transitions", units = "unitless")
       call mkpio_def_spatial_var(pioid=pioid, varname='UNREPRESENTED_CFT_LULCC', xtype=xtype, &
            lev1name = 'cft', &
            long_name = "unrepresented crop gross LULCC transitions", units = "unitless")
    else
       call mkpio_def_spatial_var(pioid=pioid, varname='HARVEST_VH1', xtype=xtype, &
            long_name = "harvest from primary forest", units = "gC/m2/yr")
       call mkpio_def_spatial_var(pioid=pioid, varname='HARVEST_VH2', xtype=xtype, &
            long_name = "harvest from primary non-forest", units = "gC/m2/yr")
       call mkpio_def_spatial_var(pioid=pioid, varname='HARVEST_SH1', xtype=xtype, &
            long_name = "harvest from secondary mature-forest", units = "gC/m2/yr")
       call mkpio_def_spatial_var(pioid=pioid, varname='HARVEST_SH2', xtype=xtype, &
            long_name = "harvest from secondary young-forest", units = "gC/m2/yr")
       call mkpio_def_spatial_var(pioid=pioid, varname='HARVEST_SH3', xtype=xtype, &
            long_name = "harvest from secondary non-forest", units = "gC/m2/yr")
       call mkpio_def_spatial_var(pioid=pioid, varname='GRAZING', xtype=xtype, &
            long_name = "grazing of herbacous pfts", units = "gC/m2/yr")
       call mkpio_def_spatial_var(pioid=pioid, varname='FERTNITRO_CFT', xtype=xtype, &
            lev1name = 'cft', &
            long_name = "nitrogen fertilizer for each crop", units = "gN/m2/yr")
       call mkpio_def_spatial_var(pioid=pioid, varname='UNREPRESENTED_PFT_LULCC', xtype=xtype, &
            lev1name = 'natpft', &
            long_name = "unrepresented PFT gross LULCC transitions", units = "unitless")
       call mkpio_def_spatial_var(pioid=pioid, varname='UNREPRESENTED_CFT_LULCC', xtype=xtype, &
            lev1name = 'cft', &
            long_name = "unrepresented crop gross LULCC transitions", units = "unitless")
    end if  ! .not. dynlanduse

    ! Coordinate variable for indices of natural PFTs
    call mkpio_defvar(pioid=pioid, varname='natpft', xtype=PIO_INT, &
         dim1name='natpft', long_name='indices of natural PFTs', units='index')

    ! Coordinate variable for indices of CFTs
    if (num_cft > 0) then
       call mkpio_defvar(pioid=pioid, varname='cft', xtype=PIO_INT, &
            dim1name='cft', long_name='indices of CFTs', units='index')
    end if

    call mkpio_def_spatial_var(pioid=pioid, varname='LANDFRAC_PFT', xtype=xtype, &
         long_name='land fraction from pft dataset', units='unitless')

    call mkpio_def_spatial_var(pioid=pioid, varname='PFTDATA_MASK', xtype=PIO_INT, &
         long_name='land mask from pft dataset, indicative of real/fake points', units='unitless')

    if (.not. dynlanduse) then
       call mkpio_def_spatial_var(pioid=pioid, varname='PCT_NATVEG', xtype=xtype, &
            long_name='total percent natural vegetation landunit', units='unitless')
    end  if

    if (.not. dynlanduse) then
       call mkpio_def_spatial_var(pioid=pioid, varname='PCT_CROP', xtype=xtype, &
            long_name='total percent crop landunit', units='unitless')
    else
       call mkpio_def_spatial_var(pioid=pioid, varname='PCT_CROP', xtype=xtype, &
            lev1name='time', &
            long_name='total percent crop landunit', units='unitless')
       call mkpio_def_spatial_var(pioid=pioid, varname='PCT_CROP_MAX', xtype=xtype, &
            long_name='maximum total percent crop landunit during time period', units='unitless')
    end if

    if (.not. dynlanduse) then
       call mkpio_def_spatial_var(pioid=pioid, varname='PCT_NAT_PFT', xtype=xtype, &
            lev1name='natpft', &
            long_name='percent plant functional type on the natural veg landunit (% of landunit)', units='unitless')
    else
       call mkpio_def_spatial_var(pioid=pioid, varname='PCT_NAT_PFT', xtype=xtype, &
            lev1name='natpft', lev2name='time', &
            long_name='percent plant functional type on the natural veg landunit (% of landunit)', units='unitless')
       call mkpio_def_spatial_var(pioid=pioid, varname='PCT_NAT_PFT_MAX', xtype=xtype, &
            lev1name='natpft', &
            long_name='maximum percent plant functional type during time period (% of landunit)', units='unitless')
    end if

    if (num_cft > 0) then
       if (.not. dynlanduse) then
          call mkpio_def_spatial_var(pioid=pioid, varname='PCT_CFT', xtype=xtype, &
               lev1name='cft', &
               long_name='percent crop functional type on the crop landunit (% of landunit)', units='unitless')
       else
          call mkpio_def_spatial_var(pioid=pioid, varname='PCT_CFT', xtype=xtype, &
               lev1name='cft', lev2name='time', &
               long_name='percent crop functional type on the crop landunit (% of landunit)', units='unitless')
          call mkpio_def_spatial_var(pioid=pioid, varname='PCT_CFT_MAX', xtype=xtype, &
               lev1name='cft', &
               long_name='maximum percent crop functional type during time period (% of landunit)', units='unitless')
       end if
    end if

    if (.not. dynlanduse) then
       call mkpio_def_spatial_var(pioid=pioid, varname='MONTHLY_LAI', xtype=xtype,  &
            lev1name='lsmpft', lev2name='time', &
            long_name='monthly leaf area index', units='unitless')
       call mkpio_def_spatial_var(pioid=pioid, varname='MONTHLY_SAI', xtype=xtype,  &
            lev1name='lsmpft', lev2name='time', &
            long_name='monthly stem area index', units='unitless')
       call mkpio_def_spatial_var(pioid=pioid, varname='MONTHLY_HEIGHT_TOP', xtype=xtype,  &
            lev1name='lsmpft', lev2name='time', &
            long_name='monthly height top', units='meters')
       call mkpio_def_spatial_var(pioid=pioid, varname='MONTHLY_HEIGHT_BOT', xtype=xtype,  &
            lev1name='lsmpft', lev2name='time', &
            long_name='monthly height bottom', units='meters')
    end if

    if (dynlanduse) then
       call mkpio_defvar(pioid=pioid, varname='YEAR', xtype=PIO_INT,  &
            dim1name='time', &
            long_name='Year of PFT data', units='unitless')
       call mkpio_defvar(pioid=pioid, varname='time', xtype=PIO_INT,  &
            dim1name='time', &
            long_name='year', units='unitless')
       call mkpio_defvar(pioid=pioid, varname='input_pftdata_filename', xtype=PIO_CHAR,  &
            dim1name='nchar', dim2name='time',  &
            long_name='Input filepath for PFT values for this year', units='unitless')
    else
       call mkpio_defvar(pioid=pioid, varname='time', xtype=PIO_INT,  &
            dim1name='time', &
            long_name='Calendar month', units='month')
    end if

  end subroutine mkfile_define_vars

  !=================================================================================
  subroutine mkfile_output_int1d(pioid, mesh, varname, data, rc)

    ! input/output variables
    type(file_desc_t) , intent(inout) :: pioid
    type(ESMF_Mesh)   , intent(in)    :: mesh
    character(len=*)  , intent(in)    :: varname
    integer           , intent(in)    :: data(:)
    integer           , intent(out)   :: rc

    ! local variables
    type(io_desc_t)      :: pio_iodesc
    type(var_desc_t)     :: pio_varid
    integer              :: rcode
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call mkpio_iodesc_output(pioid, mesh, trim(varname), pio_iodesc, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
    rcode = pio_inq_varid(pioid, trim(varname), pio_varid)
    call pio_write_darray(pioid, pio_varid, pio_iodesc, data, rcode)
    call pio_freedecomp(pioid, pio_iodesc)

  end subroutine mkfile_output_int1d

  !=================================================================================
  subroutine mkfile_output_int2d(pioid, mesh, varname, data, rc)

    ! input/output variables
    type(file_desc_t) , intent(inout) :: pioid
    type(ESMF_Mesh)   , intent(in)    :: mesh
    character(len=*)  , intent(in)    :: varname
    integer           , intent(in)    :: data(:,:)
    integer           , intent(out)   :: rc

    ! local variables
    type(io_desc_t)      :: pio_iodesc
    type(var_desc_t)     :: pio_varid
    integer              :: rcode
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call mkpio_iodesc_output(pioid, mesh, trim(varname), pio_iodesc, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
    rcode = pio_inq_varid(pioid, trim(varname), pio_varid)
    call pio_write_darray(pioid, pio_varid, pio_iodesc, data, rcode)
    call pio_freedecomp(pioid, pio_iodesc)

  end subroutine mkfile_output_int2d

  !=================================================================================
  subroutine mkfile_output_real1d(pioid, mesh, varname, data, rc)

    ! input/output variables
    type(file_desc_t) , intent(inout) :: pioid
    type(ESMF_Mesh)   , intent(in)    :: mesh
    character(len=*)  , intent(in)    :: varname
    real(r8)          , intent(in)    :: data(:)
    integer           , intent(out)   :: rc

    ! local variables
    type(io_desc_t)      :: pio_iodesc
    type(var_desc_t)     :: pio_varid
    integer              :: rcode
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call mkpio_iodesc_output(pioid, mesh, trim(varname), pio_iodesc, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
    rcode = pio_inq_varid(pioid, trim(varname), pio_varid)
    call pio_write_darray(pioid, pio_varid, pio_iodesc, data, rcode)
    call pio_freedecomp(pioid, pio_iodesc)

  end subroutine mkfile_output_real1d

  !=================================================================================
  subroutine mkfile_output_real2d(pioid, mesh, varname, data, lev1name, rc)

    ! input/output variables
    type(file_desc_t) , intent(inout) :: pioid
    type(ESMF_Mesh)   , intent(in)    :: mesh
    character(len=*)  , intent(in)    :: varname
    character(len=*)  , intent(in)    :: lev1name
    real(r8)          , intent(in)    :: data(:,:)
    integer           , intent(out)   :: rc

    ! local variables
    type(io_desc_t)      :: pio_iodesc
    type(var_desc_t)     :: pio_varid
    integer              :: rcode
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call mkpio_iodesc_output(pioid, mesh, trim(varname), pio_iodesc, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
    rcode = pio_inq_varid(pioid, trim(varname), pio_varid)
    call pio_write_darray(pioid, pio_varid, pio_iodesc, data, rcode)
    call pio_freedecomp(pioid, pio_iodesc)
  end subroutine mkfile_output_real2d

end module mkfileMod
