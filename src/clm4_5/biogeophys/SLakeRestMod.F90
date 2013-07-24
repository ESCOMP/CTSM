module SLakeRestMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Reads from or writes restart data
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  ! save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SLakeRest
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SLakeRest( ncid, flag )
    !
    ! !DESCRIPTION:
    ! Read/Write biogeophysics information to/from restart file.
    !
    ! !USES:
    use clmtype
    use ncdio_pio
    use clm_time_manager , only : is_restart
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid ! netcdf id
    character(len=*), intent(in) :: flag     ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: c,l,g,j      ! indices
    logical :: readvar      ! determine if variable is on initial file
    character(len=128) :: varname         ! temporary
    !-----------------------------------------------------------------------

    ! Note t_lake is already in BiogeophysRest.

    ! column water state variable - lake_icefrac

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='LAKE_ICEFRAC', xtype=ncd_double, &
            dim1name='column', dim2name='levlak', switchdim=.true., &
            long_name='lake layer ice fraction', units='kg/kg')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='LAKE_ICEFRAC', data=cws%lake_icefrac, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column physical state variable - savedtke1

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='SAVEDTKE1', xtype=ncd_double,  &
            dim1name='column', &
            long_name='top lake layer eddy conductivity', units='W/(m K)')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='SAVEDTKE1', data=cps%savedtke1, &
            dim1name='column', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column physical state variable - ust_lake

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='USTLAKE', xtype=ncd_double,  &
            dim1name='column', &
            long_name='friction velocity for lakes', units='m/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='USTLAKE', data=cps%ust_lake, &
            dim1name='column', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column physical state variable - z0mg

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='Z0MG', xtype=ncd_double,  &
            dim1name='column', &
            long_name='ground momentum roughness length', units='m')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='Z0MG', data=cps%z0mg, &
            dim1name='column', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

  end subroutine SLakeRest

end module SLakeRestMod
