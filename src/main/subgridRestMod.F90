module subgridRestMod

#include "shr_assert.h"

  !------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use glc_elevclass_mod  , only : glc_get_num_elevation_classes, glc_get_elevclass_bounds
  use abortutils         , only : endrun
  use decompMod          , only : bounds_type, bounds_level_proc, gindex_global, get_global_index_array
  use decompMod          , only : subgrid_level_gridcell, subgrid_level_landunit, subgrid_level_column, subgrid_level_patch
  use domainMod          , only : ldomain
  use clm_time_manager   , only : get_curr_date
  use clm_varpar         , only : nlevsno, nlevmaxurbgrnd
  use pio                , only : file_desc_t
  use ncdio_pio          , only : ncd_int, ncd_double
  use GridcellType       , only : grc
  use LandunitType       , only : lun                
  use ColumnType         , only : col                
  use PatchType          , only : patch                
  use restUtilMod
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: subgridRestWrite              ! handle restart writes of subgrid variables
  public :: subgridRestRead               ! handle restart reads of subgrid variables
  public :: subgridRest_check_consistency ! check consistency of variables read by subgridRest
  public :: subgridRest_read_cleanup      ! do cleanup of variables allocated when reading the restart file; should be called after subgridRest and subgridRest_check_consistency are complete

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: subgridRest_write_only     ! handle restart of subgrid variables that only need to be written, not read
  private :: subgridRest_write_and_read ! handle restart of subgrid variables that need to be read as well as written
  private :: save_old_weights

  ! !PRIVATE TYPES:
  real(r8), allocatable :: pft_wtlunit_before_rest_read(:)  ! patch%wtlunit weights - saved values from before the restart read

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine subgridRestWrite(bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Handle restart writes (and defines) of subgrid variables
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds ! bounds
    type(file_desc_t), intent(inout) :: ncid   ! netCDF dataset id
    character(len=*) , intent(in)    :: flag   ! flag to determine if define or write data
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'subgridRestWrite'
    !-----------------------------------------------------------------------

    call subgridRest_write_only(bounds, ncid, flag)
    call subgridRest_write_and_read(bounds, ncid, flag)

  end subroutine subgridRestWrite


  !------------------------------------------------------------------------
  subroutine subgridRestRead(bounds, ncid)
    !
    ! !DESCRIPTION:
    ! Handle restart reads of subgrid variables
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds ! bounds
    type(file_desc_t), intent(inout) :: ncid   ! netCDF dataset id
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname='subgridRestRead' ! subroutine name
    !------------------------------------------------------------------------

    call subgridRest_write_and_read(bounds, ncid, 'read')

  end subroutine subgridRestRead

  !-----------------------------------------------------------------------
  subroutine subgridRest_write_only(bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Handle restart for variables that only need to be written, not read. This applies
    ! to variables that are time-constant and are only put on the restart file for the
    ! sake of having some additional metadata there.
    !
    ! Note that 'active' flags appear in this routine: they don't need to be read because
    ! they can be computed using other info on the restart file (particularly subgrid
    ! weights).
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds ! bounds
    type(file_desc_t), intent(inout) :: ncid   ! netCDF dataset id
    character(len=*) , intent(in)    :: flag   ! flag to determine if define, write or read data
    !
    ! !LOCAL VARIABLES:
    integer :: g,l,c,p,i             ! indices
    logical :: readvar               ! temporary
    real(r8), pointer :: rgarr(:)    ! temporary
    real(r8), pointer :: rlarr(:)    ! temporary
    real(r8), pointer :: rcarr(:)    ! temporary
    real(r8), pointer :: rparr(:)    ! temporary
    integer , pointer :: igarr(:)    ! temporary
    integer , pointer :: ilarr(:)    ! temporary
    integer , pointer :: icarr(:)    ! temporary
    integer , pointer :: iparr(:)    ! temporary
    integer           :: gindex      ! global index

    real(r8), pointer :: elevclass_bounds(:)

    real(r8), pointer :: temp2d_r(:,:) ! temporary for multi-level variables
    integer , pointer :: temp2d_i(:,:) ! temporary for multi-level variables

    character(len=*), parameter :: subname = 'subgridRest_write_only'
    !-----------------------------------------------------------------------
    
    !------------------------------------------------------------------
    ! Write gridcell info
    !------------------------------------------------------------------

    allocate(rgarr(bounds%begg:bounds%endg), igarr(bounds%begg:bounds%endg))

    call restartvar(ncid=ncid, flag=flag, varname='grid1d_lon', xtype=ncd_double, &
         dim1name='gridcell',                                          &
         long_name='gridcell longitude', units='degrees_east',         &
         interpinic_flag='skip', readvar=readvar, data=grc%londeg)

    call restartvar(ncid=ncid, flag=flag, varname='grid1d_lat', xtype=ncd_double, &
         dim1name='gridcell',                                          &
         long_name='gridcell latitude', units='degrees_north',         &
         interpinic_flag='skip', readvar=readvar, data=grc%latdeg)

    do g=bounds%begg,bounds%endg
       gindex = gindex_global(g-bounds%begg+1)
       igarr(g)= mod(gindex-1,ldomain%ni) + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='grid1d_ixy', xtype=ncd_int,    &
         dim1name='gridcell',                                          &
         long_name='2d longitude index of corresponding gridcell',     &
         interpinic_flag='skip', readvar=readvar, data=igarr)

    do g=bounds%begg,bounds%endg
       gindex = gindex_global(g-bounds%begg+1)
       igarr(g)= (gindex - 1)/ldomain%ni + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='grid1d_jxy', xtype=ncd_int,    &
         dim1name='gridcell',                                          &
         long_name='2d latitude index of corresponding gridcell',      &
         interpinic_flag='skip', readvar=readvar, data=igarr)

    deallocate(rgarr,igarr)

    !------------------------------------------------------------------
    ! Write landunit info
    !------------------------------------------------------------------

    allocate(rlarr(bounds%begl:bounds%endl), ilarr(bounds%begl:bounds%endl))

    do l=bounds%begl,bounds%endl
       rlarr(l) = grc%londeg(lun%gridcell(l))
    enddo

    call restartvar(ncid=ncid, flag=flag, varname='land1d_lon', xtype=ncd_double,  &
         dim1name='landunit',                                                      &
         long_name='landunit longitude', units='degrees_east',                     &
         interpinic_flag='skip', readvar=readvar, data=rlarr)
    
    do l=bounds%begl,bounds%endl
       rlarr(l) = grc%latdeg(lun%gridcell(l))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='land1d_lat', xtype=ncd_double,  &
         dim1name='landunit',                                                      &
         long_name='landunit latitude', units='degrees_north',                     &
         interpinic_flag='skip', readvar=readvar, data=rlarr)

    do l=bounds%begl,bounds%endl
       gindex = gindex_global(lun%gridcell(l)-bounds%begg+1)
       ilarr(l) = mod(gindex-1,ldomain%ni) + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='land1d_ixy', xtype=ncd_int,     &
         dim1name='landunit',                                                      &
         long_name='2d longitude index of corresponding landunit',                 &
         interpinic_flag='skip', readvar=readvar, data=ilarr)

    do l=bounds%begl,bounds%endl
       gindex = gindex_global(lun%gridcell(l)-bounds%begg+1)
       ilarr(l) = (gindex-1)/ldomain%ni + 1
    end do
    call restartvar(ncid=ncid, flag=flag, varname='land1d_jxy', xtype=ncd_int,     &
         dim1name='landunit',                                                      &
         long_name='2d latitude index of corresponding landunit',                  &
         interpinic_flag='skip', readvar=readvar, data=ilarr)

    ilarr = get_global_index_array(lun%gridcell(bounds%begl:bounds%endl), bounds%begl, bounds%endl, &
         subgrid_level=subgrid_level_gridcell)
    call restartvar(ncid=ncid, flag=flag, varname='land1d_gridcell_index', xtype=ncd_int, &
         dim1name='landunit',                                                             &
         long_name='gridcell index of corresponding landunit',                            &
         interpinic_flag='skip', readvar=readvar, data=ilarr)

    call restartvar(ncid=ncid, flag=flag, varname='land1d_ityplun', xtype=ncd_int, &
         dim1name='landunit',                                                      &
         long_name='landunit type (see global attributes)', units=' ',             &
         interpinic_flag='skip', readvar=readvar, data=lun%itype)

    do l=bounds%begl,bounds%endl
       if (lun%active(l)) then
          ilarr(l) = 1
       else
          ilarr(l) = 0
       end if
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='land1d_active', xtype=ncd_int,  &
         dim1name='landunit',                                                      &
         long_name='landunit active flag (1=active, 0=inactive)',                  &
         interpinic_flag='skip', readvar=readvar, data=ilarr)

    deallocate(rlarr, ilarr)

    !------------------------------------------------------------------
    ! Write column info
    !------------------------------------------------------------------

    allocate(rcarr(bounds%begc:bounds%endc), icarr(bounds%begc:bounds%endc))

    do c= bounds%begc, bounds%endc
       rcarr(c) = grc%londeg(col%gridcell(c))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_lon', xtype=ncd_double,   &
         dim1name='column',                                                         &
         long_name='column longitude', units='degrees_east',                        &
         interpinic_flag='skip', readvar=readvar, data=rcarr)

    do c= bounds%begc, bounds%endc
       rcarr(c) = grc%latdeg(col%gridcell(c))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_lat', xtype=ncd_double,   &
         dim1name='column',                                                         &
         long_name='column latitude', units='degrees_north',                        &
         interpinic_flag='skip', readvar=readvar, data=rcarr)

    do c= bounds%begc, bounds%endc
       gindex = gindex_global(col%gridcell(c)-bounds%begg+1)
       icarr(c) = mod(gindex-1,ldomain%ni) + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_ixy', xtype=ncd_int,      &
         dim1name='column',                                                         &
         long_name='2d longitude index of corresponding column', units=' ',         &
         interpinic_flag='skip', readvar=readvar, data=icarr)

    do c= bounds%begc, bounds%endc
       gindex = gindex_global(col%gridcell(c)-bounds%begg+1)
       icarr(c) = (gindex-1)/ldomain%ni + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_jxy', xtype=ncd_int,      &
         dim1name='column',                                                         &
         long_name='2d latitude index of corresponding column', units=' ',          &
         interpinic_flag='skip', readvar=readvar, data=icarr)

    icarr = get_global_index_array(col%gridcell(bounds%begc:bounds%endc), bounds%begc, bounds%endc, &
         subgrid_level=subgrid_level_gridcell)
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_gridcell_index', xtype=ncd_int, &
         dim1name='column',                                                               &
         long_name='gridcell index of corresponding column',                              &
         interpinic_flag='skip', readvar=readvar, data=icarr)

    icarr = get_global_index_array(col%landunit(bounds%begc:bounds%endc), bounds%begc, bounds%endc, &
         subgrid_level=subgrid_level_landunit)
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_landunit_index', xtype=ncd_int, &
         dim1name='column',                                                               &
         long_name='landunit index of corresponding column',                              &
         interpinic_flag='skip', readvar=readvar, data=icarr)

    do c= bounds%begc, bounds%endc
       icarr(c) = lun%itype(col%landunit(c))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_ityplun', xtype=ncd_int,  &
         dim1name='column',                                                         &
         long_name='column landunit type (see global attributes)', units=' ',       &
         interpinic_flag='skip', readvar=readvar, data=icarr)

    call restartvar(ncid=ncid, flag=flag, varname='cols1d_ityp', xtype=ncd_int,     &
         dim1name='column',                                                         &
         long_name='column type (see global attributes)', units=' ',                &
         interpinic_flag='skip', readvar=readvar, data=col%itype)

    do c=bounds%begc,bounds%endc
       if (col%active(c)) then 
          icarr(c) = 1
       else
          icarr(c) = 0
       end if
    end do
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_active', xtype=ncd_int,   &
         dim1name='column',                                                         &
         long_name='column active flag (1=active, 0=inactive)', units=' ',          &
         interpinic_flag='skip', readvar=readvar, data=icarr)

    call restartvar(ncid=ncid, flag=flag, varname='LEVGRND_CLASS', xtype=ncd_int,   &
         dim1name='column', dim2name='levmaxurbgrnd', switchdim=.true.,                   &
         long_name='class in which each layer falls', units=' ',                    &
         scale_by_thickness=.false., &
         interpinic_flag='skip', readvar=readvar, data=col%levgrnd_class)

    allocate(temp2d_r(bounds%begc:bounds%endc, 1:nlevmaxurbgrnd))
    temp2d_r(bounds%begc:bounds%endc, 1:nlevmaxurbgrnd) = col%z(bounds%begc:bounds%endc, 1:nlevmaxurbgrnd)
    call restartvar(ncid=ncid, flag=flag, varname='COL_Z', xtype=ncd_double,  & 
         dim1name='column', dim2name='levmaxurbgrnd', switchdim=.true., &
         long_name='layer depth, excluding snow layers', units='m', &
         scale_by_thickness=.false., &
         interpinic_flag='skip', readvar=readvar, data=temp2d_r)
    deallocate(temp2d_r)

    allocate(temp2d_r(bounds%begc:bounds%endc, 1:nlevmaxurbgrnd))
    temp2d_r(bounds%begc:bounds%endc, 1:nlevmaxurbgrnd) = col%dz(bounds%begc:bounds%endc, 1:nlevmaxurbgrnd)
    call restartvar(ncid=ncid, flag=flag, varname='DZSOI', xtype=ncd_double,  &
         dim1name='column', dim2name='levmaxurbgrnd', switchdim=.true., &
         long_name='soil layer thickness', units='m', &
         scale_by_thickness=.false., &
         interpinic_flag='skip', readvar=readvar, data=temp2d_r)
    deallocate(temp2d_r)

    deallocate(rcarr, icarr)

    !------------------------------------------------------------------
    ! Write patch info
    !------------------------------------------------------------------

    allocate(rparr(bounds%begp:bounds%endp), iparr(bounds%begp:bounds%endp))

    do p=bounds%begp,bounds%endp
       rparr(p) = grc%londeg(patch%gridcell(p))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_lon', xtype=ncd_double, &
         dim1name='pft',                                                          &
         long_name='pft longitude', units='degrees_east',                         &
         interpinic_flag='skip', readvar=readvar, data=rparr)

    do p=bounds%begp,bounds%endp
       rparr(p) = grc%latdeg(patch%gridcell(p))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_lat', xtype=ncd_double, &
         dim1name='pft',                                                          &
         long_name='pft latitude', units='degrees_north',                         &
         interpinic_flag='skip', readvar=readvar, data=rparr)

    do p=bounds%begp,bounds%endp
       gindex = gindex_global(patch%gridcell(p)-bounds%begg+1)
       iparr(p) = mod(gindex-1,ldomain%ni) + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_ixy', xtype=ncd_int, &
         dim1name='pft',                                                       &
         long_name='2d longitude index of corresponding pft', units='',        &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    do p=bounds%begp,bounds%endp
       gindex = gindex_global(patch%gridcell(p)-bounds%begg+1)
       iparr(p) = (gindex-1)/ldomain%ni + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_jxy', xtype=ncd_int, &
         dim1name='pft',                                                       &
         long_name='2d latitude index of corresponding pft', units='',         &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    iparr = get_global_index_array(patch%gridcell(bounds%begp:bounds%endp), bounds%begp, bounds%endp, &
         subgrid_level=subgrid_level_gridcell)
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_gridcell_index', xtype=ncd_int, &
         dim1name='pft',                                                                  &
         long_name='gridcell index of corresponding pft',                                 &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    iparr = get_global_index_array(patch%landunit(bounds%begp:bounds%endp), bounds%begp, bounds%endp, &
         subgrid_level=subgrid_level_landunit)
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_landunit_index', xtype=ncd_int, &
         dim1name='pft',                                                                  &
         long_name='landunit index of corresponding pft',                                 &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    iparr = get_global_index_array(patch%column(bounds%begp:bounds%endp), bounds%begp, bounds%endp, &
         subgrid_level=subgrid_level_column)
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_column_index', xtype=ncd_int,   &
         dim1name='pft',                                                                  &
         long_name='column index of corresponding pft',                                   &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_itypveg', xtype=ncd_int,  &
         dim1name='pft',                                                            &
         long_name='pft vegetation type', units='',                                 &
         interpinic_flag='skip', readvar=readvar, data=patch%itype)

    do p=bounds%begp,bounds%endp
       iparr(p) = col%itype(patch%column(p))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_itypcol', xtype=ncd_int, &
         dim1name='pft',                                                           &
         long_name='pft column type (see global attributes)', units='',          &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    do p=bounds%begp,bounds%endp
       iparr(p) = lun%itype(patch%landunit(p))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_ityplun', xtype=ncd_int, &
         dim1name='pft',                                                           &
         long_name='pft landunit type (see global attributes)', units='',          &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    do p=bounds%begp,bounds%endp
       if (patch%active(p)) then
          iparr(p) = 1
       else
          iparr(p) = 0
       end if
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_active', xtype=ncd_int, &
         dim1name='pft',                                                          &
         long_name='pft active flag (1=active, 0=inactive)', units='',            &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    allocate(temp2d_i(bounds%begp:bounds%endp, 1:nlevmaxurbgrnd))
    do p=bounds%begp,bounds%endp
       c = patch%column(p)
       temp2d_i(p, 1:nlevmaxurbgrnd) = col%levgrnd_class(c, 1:nlevmaxurbgrnd)
    end do
    call restartvar(ncid=ncid, flag=flag, varname='LEVGRND_CLASS_p', xtype=ncd_int, &
         dim1name='pft', dim2name='levmaxurbgrnd', switchdim=.true., &
         long_name='class in which each layer falls, patch-level', units=' ', &
         scale_by_thickness=.false., &
         interpinic_flag='skip', readvar=readvar, data=temp2d_i)
    deallocate(temp2d_i)

    allocate(temp2d_r(bounds%begp:bounds%endp, 1:nlevmaxurbgrnd))
    do p=bounds%begp,bounds%endp
       c = patch%column(p)
       temp2d_r(p, 1:nlevmaxurbgrnd) = col%z(c, 1:nlevmaxurbgrnd)
    end do
    call restartvar(ncid=ncid, flag=flag, varname='COL_Z_p', xtype=ncd_double, &
         dim1name='pft', dim2name='levmaxurbgrnd', switchdim=.true., &
         long_name='layer depth, excluding snow layers, patch-level', units='m', &
         scale_by_thickness=.false., &
         interpinic_flag='skip', readvar=readvar, data=temp2d_r)
    deallocate(temp2d_r)

    allocate(temp2d_r(bounds%begp:bounds%endp, 1:nlevmaxurbgrnd))
    do p=bounds%begp,bounds%endp
       c = patch%column(p)
       temp2d_r(p, 1:nlevmaxurbgrnd) = col%dz(c, 1:nlevmaxurbgrnd)
    end do
    call restartvar(ncid=ncid, flag=flag, varname='DZSOI_p', xtype=ncd_double, &
         dim1name='pft', dim2name='levmaxurbgrnd', switchdim=.true., &
         long_name='soil layer thickness, patch-level', units='m', &
         scale_by_thickness=.false., &
         interpinic_flag='skip', readvar=readvar, data=temp2d_r)
    deallocate(temp2d_r)

    deallocate(rparr, iparr)

    ! ------------------------------------------------------------------------
    ! Write other subgrid-related metadata
    ! ------------------------------------------------------------------------

    allocate(elevclass_bounds(0:glc_get_num_elevation_classes()))
    elevclass_bounds = glc_get_elevclass_bounds()
    call restartvar(ncid=ncid, flag=flag, varname='glc_elevclass_bounds', xtype=ncd_double, &
         dim1name='glc_nec1', is_spatial = .false., &
         long_name='glacier elevation class bounds', units='m', &
         interpinic_flag='skip', readvar=readvar, data=elevclass_bounds)
    deallocate(elevclass_bounds)

  end subroutine subgridRest_write_only

  !-----------------------------------------------------------------------
  subroutine subgridRest_write_and_read(bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds ! bounds
    type(file_desc_t), intent(inout) :: ncid   ! netCDF dataset id
    character(len=*) , intent(in)    :: flag   ! flag to determine if define, write or read data
    !
    ! !LOCAL VARIABLES:
    logical :: readvar              ! temporary
    real(r8), pointer :: temp2d(:,:) ! temporary for sno column variables
    
    character(len=*), parameter :: subname = 'subgridRest_write_and_read'
    !-----------------------------------------------------------------------
    
    if (flag == 'read') then
       call save_old_weights(bounds)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='land1d_wtxy', xtype=ncd_double, &
         dim1name='landunit',                                                      &
         long_name='landunit weight relative to corresponding gridcell',           &
         interpinic_flag='area', readvar=readvar, data=lun%wtgcell)

    call restartvar(ncid=ncid, flag=flag, varname='cols1d_wtxy', xtype=ncd_double,  &
         dim1name='column',                                                         &
         long_name='column weight relative to corresponding gridcell', units=' ',   &
         interpinic_flag='area', readvar=readvar, data=col%wtgcell)

    call restartvar(ncid=ncid, flag=flag, varname='cols1d_wtlnd', xtype=ncd_double, &
         dim1name='column',                                                         &
         long_name='column weight relative to corresponding landunit', units=' ',   &
         interpinic_flag='area', readvar=readvar, data=col%wtlunit)

    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_wtxy', xtype=ncd_double,  &
         dim1name='pft',                                                            &
         long_name='pft weight relative to corresponding gridcell', units='',       &  
         interpinic_flag='area', readvar=readvar, data=patch%wtgcell)

    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_wtlnd', xtype=ncd_double, &
         dim1name='pft',                                                            &
         long_name='pft weight relative to corresponding landunit', units='',       & 
         interpinic_flag='area', readvar=readvar, data=patch%wtlunit)

    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_wtcol', xtype=ncd_double, &
         dim1name='pft',                                                            &
         long_name='pft weight relative to corresponding column', units='',         &
         interpinic_flag='area', readvar=readvar, data=patch%wtcol)

    ! Snow column variables

    call restartvar(ncid=ncid, flag=flag, varname='SNLSNO', xtype=ncd_int,  & 
         dim1name='column', &
         long_name='negative number of snow layers', units='unitless', &
         interpinic_flag='interp', readvar=readvar, data=col%snl)

    allocate(temp2d(bounds%begc:bounds%endc,-nlevsno+1:0))
    if (flag == 'write') then
       temp2d(bounds%begc:bounds%endc,-nlevsno+1:0) = col%dz(bounds%begc:bounds%endc,-nlevsno+1:0)
    end if
    call restartvar(ncid=ncid, flag=flag, varname='DZSNO', xtype=ncd_double,  & 
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer thickness', units='m', &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=temp2d)
    if (flag == 'read') then
       col%dz(bounds%begc:bounds%endc,-nlevsno+1:0) = temp2d(bounds%begc:bounds%endc,-nlevsno+1:0) 
    end if
    deallocate(temp2d)

    allocate(temp2d(bounds%begc:bounds%endc,-nlevsno+1:0))
    if (flag == 'write') then
       temp2d(bounds%begc:bounds%endc,-nlevsno+1:0) = col%z(bounds%begc:bounds%endc,-nlevsno+1:0)
    end if
    call restartvar(ncid=ncid, flag=flag, varname='ZSNO', xtype=ncd_double,  & 
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer depth', units='m', &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=temp2d)
    if (flag == 'read') then
       col%z(bounds%begc:bounds%endc,-nlevsno+1:0) = temp2d(bounds%begc:bounds%endc,-nlevsno+1:0) 
    end if
    deallocate(temp2d)

    allocate(temp2d(bounds%begc:bounds%endc,-nlevsno:-1))
    if (flag == 'write') then
       temp2d(bounds%begc:bounds%endc,-nlevsno:-1) = col%zi(bounds%begc:bounds%endc,-nlevsno:-1)
    end if
    call restartvar(ncid=ncid, flag=flag, varname='ZISNO', xtype=ncd_double,  & 
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno, upperb2=-1, &
         long_name='snow interface depth at the top of the given layer', units='m', &
         scale_by_thickness=.false., &
         interpinic_flag='interp', readvar=readvar, data=temp2d)
    if (flag == 'read') then
       col%zi(bounds%begc:bounds%endc,-nlevsno:-1) = temp2d(bounds%begc:bounds%endc,-nlevsno:-1) 
    end if
    deallocate(temp2d)

  end subroutine subgridRest_write_and_read

  !-----------------------------------------------------------------------
  subroutine save_old_weights(bounds)
    !
    ! !DESCRIPTION:
    ! Save old weights, from before the restart read, for later consistency checks.
    !
    ! !USES:
    type(bounds_type), intent(in)    :: bounds ! bounds (expected to be proc-level)
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'save_old_weights'
    !-----------------------------------------------------------------------
    
    SHR_ASSERT(bounds%level == bounds_level_proc, subname//' ERROR: expect proc-level bounds')

    allocate(pft_wtlunit_before_rest_read(bounds%begp:bounds%endp))
    pft_wtlunit_before_rest_read(bounds%begp:bounds%endp) = patch%wtlunit(bounds%begp:bounds%endp)

  end subroutine save_old_weights


  !-----------------------------------------------------------------------
  subroutine subgridRest_check_consistency(bounds)
    !
    ! !DESCRIPTION:
    ! Check consistency of variables read by subgridRest.
    !
    ! This should be called AFTER subgridRest is called to read the restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds ! bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'subgridRest_check_consistency'
    !-----------------------------------------------------------------------

    if (do_check_weights()) then
       call check_weights(bounds)
    end if

  contains

    !-----------------------------------------------------------------------
    logical function do_check_weights()
      !
      ! !DESCRIPTION:
      ! Return true if we should check weights
      !
      ! !USES:
      use clm_varctl, only : nsrest, nsrContinue, use_cndv, use_fates
      use dynSubgridControlMod, only : get_do_transient_pfts
      !
      ! !ARGUMENTS:
      !
      ! !LOCAL VARIABLES:
      
      character(len=*), parameter :: subname = 'do_check_weights'
      !-----------------------------------------------------------------------
      
      if (get_do_transient_pfts()) then
         ! Don't check weights for a transient PATCH case, because it's harder to come up with the
         ! correct weights to check against
         do_check_weights = .false.
      else if (nsrest == nsrContinue) then
         ! Don't check weights for a restart run
         !
         ! WJS (3-25-14): I'm not sure why we don't do the check in this case, but I'm
         ! maintaining the logic that used to be in BiogeophysRestMod regarding these
         ! weight checks
         do_check_weights = .false.
      else if (use_cndv) then
         ! Don't check weights for a cndv case, because the weights will almost certainly
         ! differ from the surface dataset in this case
         do_check_weights = .false.
      else if (use_fates) then
         ! Don't check weights for a fates case, because the weights will almost certainly
         ! differ from the surface dataset in this case
         do_check_weights = .false.
      else
         do_check_weights = .true.
      end if

    end function do_check_weights

    !-----------------------------------------------------------------------
    subroutine check_weights(bounds)
      !
      ! !DESCRIPTION:
      ! Make sure that patch weights on the landunit agree with the weights read from the
      ! surface dataset, for the natural veg landunit.
      !
      ! Note that we do NOT do a more general check of all subgrid weights, because it's
      ! possible that some other subgrid weights have changed relative to the surface
      ! dataset, e.g., due to dynamic landunits. It would probably be possible to do more
      ! checking than is done here, but the check here should be sufficient to catch major
      ! inconsistencies between the restart file and the surface dataset.
      !
      ! !USES:
      use landunit_varcon, only : istsoil
      use clm_varctl, only : iulog
      !
      ! !ARGUMENTS:
      type(bounds_type), intent(in)    :: bounds ! bounds
      !
      ! !LOCAL VARIABLES:
      integer  :: p, l ! indices
      real(r8) :: diff ! difference in weights

      real(r8), parameter :: tol = 5.e-3  ! tolerance for checking weights
      
      character(len=*), parameter :: subname = 'check_weights'
      !-----------------------------------------------------------------------
      
      do p = bounds%begp, bounds%endp
         l = patch%landunit(p)
         if (lun%itype(l) == istsoil) then
            diff = abs(patch%wtlunit(p) - pft_wtlunit_before_rest_read(p))
            if (diff > tol .and. patch%wtgcell(p) > 1.0e-16_r8) then
               write(iulog,*) 'ERROR: PATCH weights are SIGNIFICANTLY different between :'
               write(iulog,*) 'the restart (finidat) file : ', patch%wtlunit(p)
               write(iulog,*) 'and the surface dataset (fsurdat): ', pft_wtlunit_before_rest_read(p)
               write(iulog,*) 'weight gridcell: ', patch%wtgcell(p)
               write(iulog,*)
               write(iulog,*) 'Maximum allowed difference: ', tol
               write(iulog,*) 'Difference found: ', diff
               write(iulog,*) 'This match is a requirement for non-transient runs'
               write(iulog,*)
               write(iulog,*) 'Possible solutions to this problem:'
               write(iulog,*) '(1) Make sure you are using the intended finidat and fsurdat files'
               write(iulog,*) '(2) If you are running a present-day simulation, then make sure that your'
               write(iulog,*) '    initial conditions file is from the END of a 20th century transient run'
               write(iulog,*) '(3) If you are confident that you are using the correct finidat and fsurdat files,'
               write(iulog,*) '    yet are still experiencing this error, then you can bypass this check by setting:'
               write(iulog,*) '      check_finidat_pct_consistency = .false.'
               write(iulog,*) '    in user_nl_clm'
               write(iulog,*) '    In this case, CLM will take the weights from the initial conditions file.'
               write(iulog,*) ' '
               call endrun(subgrid_index=p, subgrid_level=subgrid_level_patch, msg=errMsg(sourcefile, __LINE__))
            end if
         end if
      end do

    end subroutine check_weights

  end subroutine subgridRest_check_consistency


  !-----------------------------------------------------------------------
  subroutine subgridRest_read_cleanup
    !
    ! !DESCRIPTION:
    ! Do cleanup of variables allocated when reading the restart file
    !
    ! Should be called after subgridRest and subgridRest_check_consistency are complete.
    ! Note that this must be called after subgridRest is called to read the restart file,
    ! in order to avoid a memory leak.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'subgridRest_read_cleanup'
    !-----------------------------------------------------------------------
    
    deallocate(pft_wtlunit_before_rest_read)

  end subroutine subgridRest_read_cleanup


end module subgridRestMod
