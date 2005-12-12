module mkfileMod

contains

!-----------------------------------------------------------------------
  subroutine mkfile(lsmlon, lsmlat, ncid)

    use shr_kind_mod, only : r8 => shr_kind_r8
    use shr_sys_mod , only : shr_sys_getenv
    use fileutils   , only : get_filename
    use mkvarctl
    use ncdio

    implicit none
    integer, intent(in) :: lsmlon, lsmlat
    integer, intent(in) :: ncid

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
    character(len=32) :: subname = 'mkfile'  ! subroutine name
!-----------------------------------------------------------------------

    ! Define dimensions.

    call check_ret(nf_def_dim (ncid, 'lsmlon' , lsmlon      , dimid), subname)
    call check_ret(nf_def_dim (ncid, 'lsmlat' , lsmlat      , dimid), subname)
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

    str = '$Name$'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Version', len_trim(str), trim(str)), subname)

    str = '$Id$'
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

  end subroutine mkfile

end module mkfileMod
