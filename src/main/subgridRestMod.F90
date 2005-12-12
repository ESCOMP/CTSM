module subgridRestMod

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subgridRest
!
! !INTERFACE:
  subroutine subgridRest( ncid, flag )

    use shr_kind_mod    , only : r8 => shr_kind_r8
    use clmtype
    use ncdio           
    use decompMod       , only : get_proc_bounds, get_proc_global, map_dc2sn
    use initGridCellsMod, only : get_sn_land1d, get_sn_cols1d, get_sn_pfts1d
    use time_manager    , only : get_curr_date
    use spmdMod         , only : masterproc
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid           ! netCDF dataset id
    character(len=*), intent(in) :: flag  ! flag to determine if define, write or read data
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: g,l,c,p,j,i         ! indices
    integer :: yr                  ! current year (0 -> ...)
    integer :: mon                 ! current month (1 -> 12)
    integer :: day                 ! current day (1 -> 31)
    integer :: mcsec               ! seconds of current date
    integer :: mcdate              ! current date
    logical :: readvar             ! determine if variable is on initial file
    integer :: begp, endp          ! per-proc beginning and ending pft indices
    integer :: begc, endc          ! per-proc beginning and ending column indices
    integer :: begl, endl          ! per-proc beginning and ending landunit indices
    integer :: begg, endg          ! per-proc gridcell ending gridcell indices
    integer :: numg                ! total number of gridcells across all processors
    integer :: numl                ! total number of landunits across all processors
    integer :: numc                ! total number of columns across all processors
    integer :: nump                ! total number of pfts across all processors
    integer :: ier                 ! error status
    integer, pointer :: iltemp(:)  ! temporary
    integer, pointer :: ictemp(:)  ! temporary
    integer, pointer :: iptemp(:)  ! temporary
    integer, pointer :: ictype(:)  ! temporary
    integer, pointer :: iptype(:)  ! temporary
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
    character(len=32) :: subname='SubgridRest' ! subroutine name
!------------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Get relevant sizes

    call get_proc_global(numg, numl, numc, nump)
    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! Allocate dynamic memory

    if (flag == 'write') then
       allocate(iltemp(numl), ictemp(numc), iptemp(nump), ictype(numc), iptype(nump), stat=ier)
       if (ier /= 0) then
          write(6,*)'allocation error from inicfile_fields'; call endrun()
       end if
    end if

    ! Write output data (first write current date and seconds of current date)

    if (flag == 'define') then
       if (masterproc) then
          call ncd_defvar(ncid=ncid, varname='mcdate', xtype=nf_int, &
               long_name='current date as 8 digit integer (YYYYMMDD)')
          call ncd_defvar(ncid=ncid, varname='mcsec', xtype=nf_int,  &
               long_name='current seconds of current date', units='s')
       else if (flag == 'write') then
          call get_curr_date (yr, mon, day, mcsec)
          mcdate = yr*10000 + mon*100 + day
          call ncd_ioglobal(varname='mcdate', data=mcdate, ncid=ncid, flag=flag, readvar=readvar)
          call ncd_ioglobal(varname='mcsec' , data=mcsec , ncid=ncid, flag=flag, readvar=readvar)
       end if
    end if

    ! Write gridcell info

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='grid1d_lon', xtype=nf_double,  &
            dim1name='gridcell', long_name='gridcell longitude', units='degrees_east')
       call ncd_defvar(ncid=ncid, varname='grid1d_lat', xtype=nf_double,  &
            dim1name='gridcell', long_name='gridcell latitude', units='degrees_north')
       call ncd_defvar(ncid=ncid, varname='grid1d_ixy', xtype=nf_int,  &
            dim1name='gridcell', long_name='2d longitude index of corresponding gridcell')
       call ncd_defvar(ncid=ncid, varname='grid1d_jxy', xtype=nf_int,  &
            dim1name='gridcell', long_name='2d latitude index of corresponding gridcell')
    else if (flag == 'write') then
       call ncd_iolocal(varname='grid1d_lon', data=gptr%londeg, dim1name='gridcell', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='grid1d_lat', data=gptr%latdeg, dim1name='gridcell', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='grid1d_ixy', data=gptr%ixy   , dim1name='gridcell', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='grid1d_jxy', data=gptr%jxy   , dim1name='gridcell', ncid=ncid, flag=flag)
    end if

    ! Write landunit info

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='land1d_lon', xtype=nf_double,  &
            dim1name='landunit', long_name='landunit longitude', units='degrees_east')
       call ncd_defvar(ncid=ncid, varname='land1d_lat', xtype=nf_double,  &
            dim1name='landunit', long_name='landunit latitude', units='degrees_north')
       call ncd_defvar(ncid=ncid, varname='land1d_ixy', xtype=nf_int,  &
            dim1name='landunit', long_name='2d longitude index of corresponding landunit')
       call ncd_defvar(ncid=ncid, varname='land1d_jxy', xtype=nf_int,  &
            dim1name='landunit', long_name='2d latitude index of corresponding landunit')
       call ncd_defvar(ncid=ncid, varname='land1d_gi', xtype=nf_int,  &
            dim1name='landunit', long_name='1d grid index of corresponding landunit')
       call ncd_defvar(ncid=ncid, varname='land1d_wtxy', xtype=nf_double,  &
            dim1name='landunit', long_name='landunit weight relative to corresponding gridcell')
       call ncd_defvar(ncid=ncid, varname='land1d_ityplun', xtype=nf_int,  &
            dim1name='landunit', long_name='landunit type (vegetated,urban,lake,wetland or glacier)')
    else if (flag == 'write') then
       call ncd_iolocal(varname='land1d_lon'    , data=lptr%londeg , dim1name='landunit', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='land1d_lat'    , data=lptr%latdeg , dim1name='landunit', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='land1d_ixy'    , data=lptr%ixy    , dim1name='landunit', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='land1d_jxy'    , data=lptr%jxy    , dim1name='landunit', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='land1d_wtxy'   , data=lptr%wtgcell, dim1name='landunit', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='land1d_ityplun', data=lptr%itype  , dim1name='landunit', ncid=ncid, flag=flag)
       if (masterproc) then
          call get_sn_land1d(iltemp, type1d='gridcell',numl=numl)
          call ncd_ioglobal(varname='land1d_gi', data=iltemp, ncid=ncid, flag=flag)
       end if
    end if

    ! Write column info

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cols1d_lon', xtype=nf_double,  &
            dim1name='column', long_name='column longitude', units='degrees_east')
       call ncd_defvar(ncid=ncid, varname='cols1d_lat', xtype=nf_double,  &
            dim1name='column', long_name='column latitude', units='degrees_north')
       call ncd_defvar(ncid=ncid, varname='cols1d_ixy', xtype=nf_int,   &
            dim1name='column', long_name='2d longitude index of corresponding column')
       call ncd_defvar(ncid=ncid, varname='cols1d_jxy', xtype=nf_int,   &
            dim1name='column', long_name='2d latitude index of corresponding column')
       call ncd_defvar(ncid=ncid, varname='cols1d_gi', xtype=nf_int,   &
            dim1name='column', long_name='1d grid index of corresponding column')
       call ncd_defvar(ncid=ncid, varname='cols1d_li', xtype=nf_int,   &
            dim1name='column', long_name='1d landunit index of corresponding column')
       call ncd_defvar(ncid=ncid, varname='cols1d_wtxy', xtype=nf_double,   &
            dim1name='column', long_name='column weight relative to corresponding gridcell')
       call ncd_defvar(ncid=ncid, varname='cols1d_wtlnd', xtype=nf_double,   &
            dim1name='column', long_name='column weight relative to corresponding landunit')
       call ncd_defvar(ncid=ncid, varname='cols1d_ityplun', xtype=nf_int,   &
            dim1name='column', long_name='column landunit type (vegetated,urban,lake,wetland or glacier)')
    else if (flag == 'write') then
       call ncd_iolocal(varname='cols1d_lon'  , data=cptr%londeg , dim1name='column', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='cols1d_lat'  , data=cptr%latdeg , dim1name='column', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='cols1d_ixy'  , data=cptr%ixy    , dim1name='column', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='cols1d_jxy'  , data=cptr%jxy    , dim1name='column', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='cols1d_wtxy' , data=cptr%wtgcell, dim1name='column', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='cols1d_wtlnd', data=cptr%wtlunit, dim1name='column', ncid=ncid, flag=flag)
       if (masterproc) then
          call get_sn_cols1d(ictemp, type1d='gridcell',numc=numc)
          call ncd_ioglobal(varname='cols1d_gi', data=ictemp, ncid=ncid, flag=flag)
          call get_sn_cols1d(ictemp, type1d='landunit',numc=numc)
          call ncd_ioglobal(varname='cols1d_li', data=ictemp, ncid=ncid, flag=flag)
          do c = 1,numc
             l = cptr%landunit(c)
             ictype(c) = lptr%itype(l)
          end do
          call map_dc2sn(ictype, ictemp, type1d=namec)
          call ncd_ioglobal(varname='cols1d_ityplun', data=ictemp, ncid=ncid, flag=flag)
       end if
    end if

    ! Write pft info

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pfts1d_lon', xtype=nf_double,  &
            dim1name='pft', long_name='pft longitude', units='degrees_east')
       call ncd_defvar(ncid=ncid, varname='pfts1d_lat', xtype=nf_double,  &
            dim1name='pft', long_name='pft latitude', units='degrees_north')
       call ncd_defvar(ncid=ncid, varname='pfts1d_ixy', xtype=nf_int,  &
            dim1name='pft', long_name='2d longitude index of corresponding pft')
       call ncd_defvar(ncid=ncid, varname='pfts1d_jxy', xtype=nf_int,  &
            dim1name='pft', long_name='2d latitude index of corresponding pft')
       call ncd_defvar(ncid=ncid, varname='pfts1d_gi', xtype=nf_int,  &
            dim1name='pft', long_name='1d grid index of corresponding pft')
       call ncd_defvar(ncid=ncid, varname='pfts1d_li', xtype=nf_int,  &
            dim1name='pft', long_name='1d landunit index of corresponding pft')
       call ncd_defvar(ncid=ncid, varname='pfts1d_ci', xtype=nf_int,  &
            dim1name='pft', long_name='1d column index of corresponding pft')
       call ncd_defvar(ncid=ncid, varname='pfts1d_wtxy', xtype=nf_double,  &
            dim1name='pft', long_name='pft weight relative to corresponding gridcell')
       call ncd_defvar(ncid=ncid, varname='pfts1d_wtlnd', xtype=nf_double,  &
            dim1name='pft', long_name='pft weight relative to corresponding landunit')
       call ncd_defvar(ncid=ncid, varname='pfts1d_wtcol', xtype=nf_double,  &
            dim1name='pft', long_name='pft weight relative to corresponding column')
       call ncd_defvar(ncid=ncid, varname='pfts1d_itypveg', xtype=nf_int,  &
            dim1name='pft', long_name='pft vegetation type')
       call ncd_defvar(ncid=ncid, varname='pfts1d_ityplun', xtype=nf_int,  &
            dim1name='pft', long_name='pft landunit type (vegetated,urban,lake,wetland or glacier)')
    else if (flag == 'write') then
       call ncd_iolocal(varname='pfts1d_lon'    , data=pptr%londeg , dim1name='pft', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='pfts1d_lat'    , data=pptr%latdeg , dim1name='pft', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='pfts1d_ixy'    , data=pptr%ixy    , dim1name='pft', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='pfts1d_jxy'    , data=pptr%jxy    , dim1name='pft', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='pfts1d_wtxy'   , data=pptr%wtgcell, dim1name='pft', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='pfts1d_wtlnd'  , data=pptr%wtlunit, dim1name='pft', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='pfts1d_wtcol'  , data=pptr%wtcol  , dim1name='pft', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='pfts1d_itypveg',data=pptr%itype   , dim1name='pft', ncid=ncid, flag=flag)
       if (masterproc) then
          call get_sn_pfts1d(iptemp, type1d='gridcell', nump=nump)
          call ncd_ioglobal(varname='pfts1d_gi', data=iptemp, ncid=ncid, flag=flag)
          call get_sn_pfts1d(iptemp, type1d='landunit', nump=nump)
          call ncd_ioglobal(varname='pfts1d_li', data=iptemp, ncid=ncid, flag=flag)
          call get_sn_pfts1d(iptemp, type1d='column', nump=nump)
          call ncd_ioglobal(varname='pfts1d_ci', data=iptemp, ncid=ncid, flag=flag)
          do p = 1,nump
             l = pptr%landunit(p)
             iptype(p) = lptr%itype(l)
          end do
          call map_dc2sn(iptype, iptemp, type1d=namep)
          call ncd_ioglobal(varname='pfts1d_ityplun', data=iptemp, ncid=ncid, flag=flag)
       end if
    end if

    if (flag == 'write') then
       deallocate(iltemp, ictemp, iptemp, ictype, iptype)
    end if

  end subroutine subgridRest

end module subgridRestMod
