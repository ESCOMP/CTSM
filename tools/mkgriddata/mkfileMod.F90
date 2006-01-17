module mkfileMod

contains

!-----------------------------------------------------------------------
  subroutine mkfile(lsmlon, lsmlat, fname, finfo)

    use shr_kind_mod, only : r8 => shr_kind_r8
    use shr_sys_mod , only : shr_sys_getenv
    use fileutils   , only : get_filename
    use mkvarctl
    use ncdio

    implicit none
    integer, intent(in) :: lsmlon, lsmlat
    character(len=*),intent(in) :: fname
    character(len=*),intent(in) :: finfo

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

    call check_ret(nf_def_dim (ncid, 'lsmlon'    , lsmlon      , dimid), subname)
    call check_ret(nf_def_dim (ncid, 'lsmlat'    , lsmlat      , dimid), subname)
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

    str = finfo
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Input_Filename', len_trim(str), trim(str)), subname)

    str = 'Community Land Model: CLM3'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Source', len_trim(str), trim(str)), subname)

    str = '$Name: clm3_expa_48_brnchT_fmesh13 $'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Version', len_trim(str), trim(str)), subname)

    str = '$Id: mkfileMod.F90,v 1.1.2.1.2.1 2005/12/22 16:25:18 tcraig Exp $'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Revision_Id', len_trim(str), trim(str)), subname)

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

    call ncd_defvar(ncid=ncid, varname='NUMLON', xtype=nf_int, &
         dim1name='lsmlat', &
         long_name='number of grid cells at each latitude', units='unitless')

    call ncd_defvar(ncid=ncid, varname='AREA' , xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='area', units='km^2')

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

    ! End of define mode

    call check_ret(nf_enddef(ncid), subname)
    call check_ret(nf_close(ncid), subname)

  end subroutine mkfile

end module mkfileMod
