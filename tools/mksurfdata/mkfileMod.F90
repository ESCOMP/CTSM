module mkfileMod

contains

!-----------------------------------------------------------------------
  subroutine mkfile(lsmlon, lsmlat, fname, dynlanduse)

    use shr_kind_mod, only : r8 => shr_kind_r8
    use shr_sys_mod , only : shr_sys_getenv
    use fileutils   , only : get_filename
    use mkvarpar    , only : nlevsoi, numpft
    use mkvarctl
    use ncdio

    implicit none
    integer, intent(in) :: lsmlon, lsmlat
    character(len=*),intent(in) :: fname
    logical, intent(in) :: dynlanduse	

    integer :: ncid
    integer :: j                    ! index
    integer :: pftsize              ! size of lsmpft dimension
    integer :: dimid                ! temporary
    integer :: values(8)            ! temporary
    character(len=256) :: str       ! global attribute string
    character(len=256) :: name      ! name of attribute
    character(len=256) :: unit      ! units of attribute
    character(len= 18) :: datetime  ! temporary
    character(len=  8) :: date      ! temporary
    character(len= 10) :: time      ! temporary
    character(len=  5) :: zone      ! temporary
    integer            :: ier       ! error status
    integer            :: omode     ! netCDF output mode
    character(len=32) :: subname = 'mkfile'  ! subroutine name
!-----------------------------------------------------------------------

    call check_ret(nf_create(trim(fname), nf_clobber, ncid), subname)
    call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

    ! Define dimensions.

    pftsize = numpft + 1
    call check_ret(nf_def_dim (ncid, 'lsmlon' , lsmlon      , dimid), subname)
    call check_ret(nf_def_dim (ncid, 'lsmlat' , lsmlat      , dimid), subname)
    if (.not. dynlanduse) then
       call check_ret(nf_def_dim (ncid, 'nlevsoi', nlevsoi     , dimid), subname)
    end if
    call check_ret(nf_def_dim (ncid, 'lsmpft' , pftsize     , dimid), subname)
    call check_ret(nf_def_dim (ncid, 'time'   , nf_unlimited, dimid), subname)
    call check_ret(nf_def_dim (ncid, 'nchar'  , 128         , dimid), subname)

    ! Create global attributes.

    str = 'NCAR-CSM'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Conventions', len_trim(str), trim(str)), subname)

    call date_and_time (date, time, zone, values)
    datetime(1:8) =        date(5:6) // '-' // date(7:8) // '-' // date(3:4)
    datetime(9:)  = ' ' // time(1:2) // ':' // time(3:4) // ':' // time(5:6) // ' '
    str = 'created on: ' // datetime
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'History_Log', len_trim(str), trim(str)), subname)

    call shr_sys_getenv ('LOGNAME', str, ier)
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Logname', len_trim(str), trim(str)), subname)

    call shr_sys_getenv ('HOST', str, ier)
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Host', len_trim(str), trim(str)), subname)

    str = 'Community Land Model: CLM3'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Source', len_trim(str), trim(str)), subname)

    str = '$Name: clm3_expa_48_brnchT_fmesh13 $'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Version', len_trim(str), trim(str)), subname)

    str = '$Id: mkfileMod.F90,v 1.1.2.4.2.1 2005/12/22 16:25:37 tcraig Exp $'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Revision_Id', len_trim(str), trim(str)), subname)

    if (mksrf_fgrid_global /= ' ') then
       str = mksrf_fgrid_global
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Input_global_grid_dataset', len_trim(str), trim(str)), subname)
    else if (mksrf_fgrid_regional /= ' ') then
       str = mksrf_fgrid_regional
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Input_regional_grid_dataset', len_trim(str), trim(str)), subname)
    endif

    str = get_filename(mksrf_fvegtyp)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Vegetation_type_raw_data_filename', len_trim(str), trim(str)), subname)

    if (.not. dynlanduse) then
       str = get_filename(mksrf_fsoitex)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Soil_texture_raw_data_file_name', len_trim(str), trim(str)), subname)

       str = get_filename(mksrf_fsoicol)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Soil_color_raw_data_file_name', len_trim(str), trim(str)), subname)
    end if
       
    str = get_filename(mksrf_flanwat)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Inland_water_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_fglacier)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Glacier_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_furban)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Urban_raw_data_file_name', len_trim(str), trim(str)), subname)

    if (.not. dynlanduse) then
       str = get_filename(mksrf_flai)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Lai_raw_data_file_name', len_trim(str), trim(str)), subname)
    end if

    ! ----------------------------------------------------------------------
    ! Define variables
    ! ----------------------------------------------------------------------

    call ncd_defvar(ncid=ncid, varname='EDGEN', xtype=nf_double, &
         long_name='northern edge of surface grid', units='degrees north')
    
    call ncd_defvar(ncid=ncid, varname='EDGEE', xtype=nf_double, &
         long_name='eastern edge of surface grid', units='degrees east')
    
    call ncd_defvar(ncid=ncid, varname='EDGES', xtype=nf_double, &
         long_name='southern edge of surface grid', units='degrees north')
    
    call ncd_defvar(ncid=ncid, varname='EDGEW', xtype=nf_double, &
         long_name='western edge of surface grid', units='degrees east')

    call ncd_defvar(ncid=ncid, varname='LATN' , xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='latitude of north edge', units='degrees north')

    call ncd_defvar(ncid=ncid, varname='LONE' , xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='longitude of east edge', units='degrees east')

    call ncd_defvar(ncid=ncid, varname='LATS' , xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='latitude of south edge', units='degrees north')

    call ncd_defvar(ncid=ncid, varname='LONW' , xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='longitude of west edge', units='degrees east')

    call ncd_defvar(ncid=ncid, varname='AREA' , xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='area', units='km^2')

    call ncd_defvar(ncid=ncid, varname='NUMLON', xtype=nf_int, &
         dim1name='lsmlat', long_name='number of longitudes for each latitude', units='unitless')

    call ncd_defvar(ncid=ncid, varname='LONGXY', xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='longitude', units='degrees east')

    call ncd_defvar(ncid=ncid, varname='LATIXY', xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='latitude', units='degrees north')

    call ncd_defvar(ncid=ncid, varname='LANDMASK', xtype=nf_int, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='land/ocean mask', units='0=ocean and 1=land')

    call ncd_defvar(ncid=ncid, varname='LANDFRAC', xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='land fraction', units='unitless')

    call ncd_defvar(ncid=ncid, varname='LANDFRAC_PFT', xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='land fraction from pft dataset', units='unitless')

    call ncd_defvar(ncid=ncid, varname='PFTDATA_MASK', xtype=nf_int, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='land mask from pft dataset, indicative of real/fake points', units='unitless')

    if (.not. dynlanduse) then
       call ncd_defvar(ncid=ncid, varname='mxsoil_color', xtype=nf_int, &
            long_name='maximum numbers of soil colors', units='unitless')

       call ncd_defvar(ncid=ncid, varname='SOIL_COLOR', xtype=nf_int, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='soil color', units='unitless')

       call ncd_defvar(ncid=ncid, varname='PCT_SAND', xtype=nf_float, &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevsoi', &
            long_name='percent sand', units='unitless')
       
       call ncd_defvar(ncid=ncid, varname='PCT_CLAY', xtype=nf_float, &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevsoi', &
            long_name='percent clay', units='unitless')
    endif

    call ncd_defvar(ncid=ncid, varname='PCT_WETLAND', xtype=nf_float, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='percent wetland', units='unitless')

    call ncd_defvar(ncid=ncid, varname='PCT_LAKE', xtype=nf_float, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='percent lake', units='unitless')

    call ncd_defvar(ncid=ncid, varname='PCT_GLACIER', xtype=nf_float, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='percent glacier', units='unitless')

    call ncd_defvar(ncid=ncid, varname='PCT_URBAN', xtype=nf_float, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='percent urban', units='unitless')

    if (.not. dynlanduse) then
       call ncd_defvar(ncid=ncid, varname='PCT_PFT', xtype=nf_float, &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', &
            long_name='percent plant functional type of gridcell', units='unitless')
    else
       call ncd_defvar(ncid=ncid, varname='PCT_PFT', xtype=nf_float, &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
            long_name='percent plant functional type of gridcell', units='unitless')
    end if

    if (.not. dynlanduse) then
       call ncd_defvar(ncid=ncid, varname='MONTHLY_LAI', xtype=nf_float,  &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
            long_name='monthly leaf area index', units='unitless')
       
       call ncd_defvar(ncid=ncid, varname='MONTHLY_SAI', xtype=nf_float,  &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
            long_name='monthly stem area index', units='unitless')
       
       call ncd_defvar(ncid=ncid, varname='MONTHLY_HEIGHT_TOP', xtype=nf_float,  &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
            long_name='monthly height top', units='meters')
       
       call ncd_defvar(ncid=ncid, varname='MONTHLY_HEIGHT_BOT', xtype=nf_float,  &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
            long_name='monthly height bottom', units='meters')
    end if
       
    if (dynlanduse) then
       call ncd_defvar(ncid=ncid, varname='YEAR', xtype=nf_int,  &
            dim1name='time', &
            long_name='year of PFT data', units='unitless')
    end if

    ! End of define mode

    call check_ret(nf_enddef(ncid), subname)
    call check_ret(nf_close(ncid), subname)

  end subroutine mkfile

end module mkfileMod
