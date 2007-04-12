program mkgriddata

!-----------------------------------------------------------------------
!BOP
!
! !PROGRAM: mgriddata
!
! !DESCRIPTION:
! Creates land model grid dataset from original "raw" topo file
! or ccsm or cam grid file.
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use shr_sys_mod  , only : shr_sys_getenv
    use fileutils    , only : getfil, putfil, opnfil, getavu, get_filename
    use creategridMod, only : creategrid, write_domain, mkfile, settopo
    use mkvarctl
    use mkvarsur
    use areaMod
    use ncdio
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Authors: Gordon Bonan, Sam Levis and Mariana Vertenstein
! 2005.12.15 T Craig Updated
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: lsmlon, lsmlat       ! clm grid resolution
    integer  :: i,j,k,m              ! indices
    integer  :: ier                  ! error status
    character(len= 9) :: resol       ! resolution for file name
    character(len=64) :: fgriddat    ! output filename
    character(len=64) :: ffracdat    ! output filename
    character(len=64) :: ftopodat    ! output filename
    character(len=256):: fileinfo    ! output filename
    character(len=32) :: subname = 'mkgriddata'  ! program name

    namelist /clmexp/    &
         mksrf_fnavyoro, &
         mksrf_frawtopo, &
         mksrf_fcamfile, &
         mksrf_fclmgrid, &
         mksrf_fccsmdom, &
         mksrf_fcamtopo, &
         mksrf_lsmlon  , &
         mksrf_lsmlat  , &
         mksrf_edgen   , &
         mksrf_edgee   , &
         mksrf_edges   , &
         mksrf_edgew   
!-----------------------------------------------------------------------

    !--- read namelist ---
    mksrf_lsmlon = -1
    mksrf_lsmlat = -1
    mksrf_edgen  = -999.
    mksrf_edges  = -999.
    mksrf_edgee  = -999.
    mksrf_edgew  = -999.

    write(6,*) 'Attempting to initialize control settings .....'

    read(5, clmexp, iostat=ier)
    if (ier /= 0) then
       write(6,*)'error: namelist input resulted in error code ',ier
       call abort()
    endif

    !--- set ldomain ---

    write (6,*) 'Attempting to create grid dataset .....'

    fileinfo  = ' '
    if     (mksrf_fccsmdom /= ' ')  then
       write(6,*) 'Setting grid from ccsm domain file ',trim(mksrf_fccsmdom)
       area_units = 1
       area_valid = .false.
       call creategrid(mksrf_fccsmdom,'external')
    elseif (mksrf_fcamfile /= ' ') then
       write(6,*) 'Setting grid from cam grid file ',trim(mksrf_fcamfile)
       call creategrid(mksrf_fcamfile,'external')
    elseif (mksrf_fcamtopo /= ' ') then
       write(6,*) 'Setting grid from cam topo file ',trim(mksrf_fcamtopo)
       call creategrid(mksrf_fcamtopo,'external')
    elseif (mksrf_fclmgrid /= ' ') then
       write(6,*) 'Setting grid from clm grid file ',trim(mksrf_fclmgrid)
       call creategrid(mksrf_fclmgrid,'external')
       fileinfo = mksrf_fclmgrid
    elseif (mksrf_fnavyoro /= ' ') then
       write(6,*) 'Setting grid from navy oro file ',trim(mksrf_fnavyoro)
       call creategrid(mksrf_fnavyoro,'internal')
       fileinfo = mksrf_fnavyoro
    else
       write(6,*) 'MKGRID namelist error, must specify file'
    endif

    if (mksrf_frawtopo /= ' ') then
       call settopo(mksrf_frawtopo)
    endif

    !--- final comments ---

    lsmlon = ldomain%ni
    lsmlat = ldomain%nj

    write (resol,'(i4.4,"x",i4.4)') lsmlat,lsmlon
    fgriddat = './griddata_'//trim(resol)//'.nc'
    ffracdat = './fracdata_'//trim(resol)//'.nc'
    ftopodat = './topodata_'//trim(resol)//'.nc'

    call mkfile(lsmlon, lsmlat, fgriddat, fileinfo, itype=1)
    call mkfile(lsmlon, lsmlat, ffracdat, fileinfo, itype=2)
    call mkfile(lsmlon, lsmlat, ftopodat, fileinfo, itype=3)

    call write_domain(ldomain,fgriddat, itype=1)
    call write_domain(ldomain,ffracdat, itype=2)
    call write_domain(ldomain,ftopodat, itype=3)

    write (6,'(72a1)') ("-",i=1,60)
    write (6,'(a46,f5.1,a4,f5.1,a5)') 'land model grid data set successfully created for ', &
         360./lsmlon,' by ',180./lsmlat,' grid'

end program mkgriddata


