module subgridRestMod

  ! !PUBLIC TYPES:
  implicit none
  save
  public

contains

  !------------------------------------------------------------------------
  subroutine subgridRest( bounds, ncid, flag )

    use shr_kind_mod     , only : r8 => shr_kind_r8
    use decompMod        , only : bounds_type, ldecomp
    use domainMod        , only : ldomain
    use clm_time_manager , only : get_curr_date
    use pio              , only : file_desc_t
    use ncdio_pio        , only : ncd_int, ncd_double
    use restUtilMod
    use clmtype
   !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in)    :: bounds ! bounds
    type(file_desc_t), intent(inout) :: ncid   ! netCDF dataset id
    character(len=*) , intent(in)    :: flag   ! flag to determine if define, write or read data
    !
    ! !LOCAL VARIABLES:
    integer :: g,l,c,p,i           ! indices
    integer :: yr                  ! current year (0 -> ...)
    integer :: mon                 ! current month (1 -> 12)
    integer :: day                 ! current day (1 -> 31)
    integer :: mcsec               ! seconds of current date
    integer :: mcdate              ! current date
    logical :: readvar             ! temporary
    real(r8),pointer :: rgarr(:)   ! temporary
    real(r8),pointer :: rlarr(:)   ! temporary
    real(r8),pointer :: rcarr(:)   ! temporary
    real(r8),pointer :: rparr(:)   ! temporary
    integer ,pointer :: igarr(:)   ! temporary
    integer ,pointer :: ilarr(:)   ! temporary
    integer ,pointer :: icarr(:)   ! temporary
    integer ,pointer :: iparr(:)   ! temporary
    character(len=32) :: subname='SubgridRest' ! subroutine name
    !------------------------------------------------------------------------

    ! Below argument readvar is ONLY needed for API consistency - it is not used

    ! Write output data (first write current date and seconds of current date)
    call get_curr_date (yr, mon, day, mcsec)
    mcdate = yr*10000 + mon*100 + day

    call restartvar(ncid=ncid, flag=flag, varname='mcdate', xtype=ncd_int,  &
         long_name='current date as 8 digit integer (YYYYMMDD)', &
         interpinic_flag='skip', readvar=readvar, data=mcdate)

    call restartvar(ncid=ncid, flag=flag, varname='mcsec', xtype=ncd_int,   &
         long_name='current seconds of current date', units='s', &
         interpinic_flag='skip', readvar=readvar, data=mcsec)

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
       igarr(g)= mod(ldecomp%gdc2glo(g)-1,ldomain%ni) + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='grid1d_ixy', xtype=ncd_int,    &
         dim1name='gridcell',                                          &
         long_name='2d longitude index of corresponding gridcell',     &
         interpinic_flag='skip', readvar=readvar, data=igarr)

    do g=bounds%begg,bounds%endg
       igarr(g)= (ldecomp%gdc2glo(g) - 1)/ldomain%ni + 1
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
       ilarr(l) = mod(ldecomp%gdc2glo(lun%gridcell(l))-1,ldomain%ni) + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='land1d_ixy', xtype=ncd_int,     &
         dim1name='landunit',                                                      &
         long_name='2d longitude index of corresponding landunit',                 &
         interpinic_flag='skip', readvar=readvar, data=ilarr)

    do l=bounds%begl,bounds%endl
       ilarr(l) = (ldecomp%gdc2glo(lun%gridcell(l))-1)/ldomain%ni + 1
    end do
    call restartvar(ncid=ncid, flag=flag, varname='land1d_jxy', xtype=ncd_int,     &
         dim1name='landunit',                                                      &
         long_name='2d latitude index of corresponding landunit',                  &
         interpinic_flag='skip', readvar=readvar, data=ilarr)

    call restartvar(ncid=ncid, flag=flag, varname='land1d_wtxy', xtype=ncd_double, &
         dim1name='landunit',                                                      &
         long_name='landunit weight relative to corresponding gridcell',           &
         interpinic_flag='skip', readvar=readvar, data=lun%wtgcell)

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
       icarr(c) = mod(ldecomp%gdc2glo(col%gridcell(c))-1,ldomain%ni) + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_ixy', xtype=ncd_int,      &
         dim1name='column',                                                         &
         long_name='2d longitude index of corresponding column', units=' ',         &
         interpinic_flag='skip', readvar=readvar, data=icarr)

    do c= bounds%begc, bounds%endc
       icarr(c) = (ldecomp%gdc2glo(col%gridcell(c))-1)/ldomain%ni + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_jxy', xtype=ncd_int,      &
         dim1name='column',                                                         &
         long_name='2d latitude index of corresponding column', units=' ',          &
         interpinic_flag='skip', readvar=readvar, data=icarr)

    call restartvar(ncid=ncid, flag=flag, varname='cols1d_wtxy', xtype=ncd_double,  &
         dim1name='column',                                                         &
         long_name='column weight relative to corresponding gridcell', units=' ',   &
         interpinic_flag='skip', readvar=readvar, data=col%wtgcell)

    call restartvar(ncid=ncid, flag=flag, varname='cols1d_wtlnd', xtype=ncd_double, &
         dim1name='column',                                                         &
         long_name='column weight relative to corresponding landunit', units=' ',   &
         interpinic_flag='skip', readvar=readvar, data=col%wtlunit)

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

    call restartvar(ncid=ncid, flag=flag, varname='cols1d_topoglc', xtype=ncd_double,   &
         dim1name='column',                                                             &
         long_name='mean elevation on glacier elevation classes', units='m',            &
         interpinic_flag='skip', readvar=readvar, data=cps%glc_topo)

    deallocate(rcarr, icarr)

    !------------------------------------------------------------------
    ! Write pft info
    !------------------------------------------------------------------

    allocate(rparr(bounds%begp:bounds%endp), iparr(bounds%begp:bounds%endp))

    do p=bounds%begp,bounds%endp
       rparr(p) = grc%londeg(pft%gridcell(p))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_lon', xtype=ncd_double, &
         dim1name='pft',                                                          &
         long_name='pft longitude', units='degrees_east',                         &
         interpinic_flag='skip', readvar=readvar, data=rparr)

    do p=bounds%begp,bounds%endp
       rparr(p) = grc%latdeg(pft%gridcell(p))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_lat', xtype=ncd_double, &
         dim1name='pft',                                                          &
         long_name='pft latitude', units='degrees_north',                         &
         interpinic_flag='skip', readvar=readvar, data=rparr)

    do p=bounds%begp,bounds%endp
       iparr(p) = mod(ldecomp%gdc2glo(pft%gridcell(p))-1,ldomain%ni) + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_ixy', xtype=ncd_int, &
         dim1name='pft',                                                       &
         long_name='2d longitude index of corresponding pft', units='',        &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    do p=bounds%begp,bounds%endp
       iparr(p) = (ldecomp%gdc2glo(pft%gridcell(p))-1)/ldomain%ni + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_jxy', xtype=ncd_int, &
         dim1name='pft',                                                       &
         long_name='2d latitude index of corresponding pft', units='',         &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_wtxy', xtype=ncd_double,  &
         dim1name='pft',                                                            &
         long_name='pft weight relative to corresponding gridcell', units='',       &  
         interpinic_flag='skip', readvar=readvar, data=pft%wtgcell)

    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_wtlnd', xtype=ncd_double, &
         dim1name='pft',                                                            &
         long_name='pft weight relative to corresponding landunit', units='',       & 
         interpinic_flag='skip', readvar=readvar, data=pft%wtlunit)

    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_wtcol', xtype=ncd_double, &
         dim1name='pft',                                                            &
         long_name='pft weight relative to corresponding column', units='',         &
         interpinic_flag='skip', readvar=readvar, data=pft%wtcol)

    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_itypveg', xtype=ncd_int,  &
         dim1name='pft',                                                            &
         long_name='pft vegetation type', units='',                                 &
         interpinic_flag='skip', readvar=readvar, data=pft%itype)

    do p=bounds%begp,bounds%endp
       iparr(p) = col%itype(pft%column(p))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_itypcol', xtype=ncd_int, &
         dim1name='pft',                                                           &
         long_name='pft column type (see global attributes)', units='',          &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    do p=bounds%begp,bounds%endp
       iparr(p) = lun%itype(pft%landunit(p))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_ityplun', xtype=ncd_int, &
         dim1name='pft',                                                           &
         long_name='pft landunit type (see global attributes)', units='',          &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    do p=bounds%begp,bounds%endp
       if (pft%active(p)) then
          iparr(p) = 1
       else
          iparr(p) = 0
       end if
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_active', xtype=ncd_int, &
         dim1name='pft',                                                          &
         long_name='pft active flag (1=active, 0=inactive)', units='',            &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    do p=bounds%begp,bounds%endp
       c = pft%column(p)
       rparr(p) = cps%glc_topo(c)
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_topoglc', xtype=ncd_double,   &
         dim1name='column',                                                             &
         long_name='mean elevation on glacier elevation classes', units='m',            &
         interpinic_flag='skip', readvar=readvar, data=rparr)

    deallocate(rparr, iparr)

  end subroutine subgridRest

end module subgridRestMod
