program mkgriddata

!-----------------------------------------------------------------------
!BOP
!
! !PROGRAM: mgriddata
!
! !DESCRIPTION:
! Creates land model grid dataset from original "raw" topo file
! or cesm or cam grid file.
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use shr_sys_mod  , only : shr_sys_getenv
    use creategridMod, only : creategrid, write_domain, mkfile, settopo
    use domainMod    , only : domain_type
    use mkvarctl
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
!
! !LOCAL VARIABLES:
!EOP
    integer  :: lsmlon, lsmlat              ! clm grid resolution
    integer  :: i,j,k,m                     ! indices
    integer  :: nfile                       ! number of files read in
    integer  :: nused                       ! number of files used
    integer  :: ier                         ! error status
    character(len= 9) :: resol              ! resolution for file name
    character(len=64) :: fgriddat           ! output filename
    character(len=64) :: ffracdat           ! output filename
    character(len=64) :: ftopodat           ! output filename
    character(len=256):: fileinfo     = ' ' ! output filename
    character(len=256):: Topofileinfo = ' ' ! output filename
    character(len=256):: Fracfileinfo = ' ' ! output filename
    logical           :: grid               ! grid is being set
    logical           :: writeTopo          ! flag to write topo file out or not
    logical           :: writeLFrc          ! flag to write land-frac file out or not
    logical           :: writeTopo1         ! tmp flag to write topo file out or not
    logical           :: writeLFrc1         ! tmp flag to write land-frac file out or not
    type(domain_type) :: ldomain            ! local domain
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
    writeTopo = .false.
    writeLFrc = .false.
    grid      = .false.
    nfile     = 0

    if (mksrf_fcamfile /= ' ') then
       write(6,*) 'Setting grid from cam grid file ',trim(mksrf_fcamfile)
       call creategrid(mksrf_fcamfile,'mksrf_fcamfile','external',grid,ldomain,writeTopo1, writeLFrc1)
       call updateFileInfo( mksrf_fcamfile )
    end if

    if (mksrf_fcamtopo /= ' ') then
       write(6,*) 'Setting grid from cam topo file ',trim(mksrf_fcamtopo)
       call creategrid(mksrf_fcamtopo,'mksrf_fcamtopo','external',grid,ldomain,writeTopo1, writeLFrc1)
       call updateFileInfo( mksrf_fcamtopo )
    end if

    if (mksrf_fclmgrid /= ' ') then
       write(6,*) 'Setting grid from clm grid file ',trim(mksrf_fclmgrid)
       call creategrid(mksrf_fclmgrid,'mksrf_fclmgrid','external',grid,ldomain,writeTopo1, writeLFrc1)
       call updateFileInfo( mksrf_fclmgrid )
    end if

    if (mksrf_fnavyoro /= ' ') then
       write(6,*) 'Setting grid from navy oro file ',trim(mksrf_fnavyoro)
       call creategrid(mksrf_fnavyoro,'mksrf_fnavyoro','internal',grid,ldomain,writeTopo1, writeLFrc1)
       call updateFileInfo( mksrf_fnavyoro )
    end if

    if     (mksrf_fccsmdom /= ' ')  then
       write(6,*) 'Setting grid from cesm domain file ',trim(mksrf_fccsmdom)
       area_units = 1
       area_valid = .false.
       call creategrid(mksrf_fccsmdom,'mksrf_fccsmdom','external',grid,ldomain,writeTopo1, writeLFrc1)
       call updateFileInfo( mksrf_fccsmdom )
       area_units = 0
       area_valid = .true.
    end if

    if ( grid .and. mksrf_frawtopo /= ' ') then
       write(6,*) 'Setting topo from raw topo file ',trim(mksrf_frawtopo)
       writeTopo1 = .true.
       writeLFrc1 = .false.
       call settopo(mksrf_frawtopo,ldomain)
       call updateFileInfo( mksrf_frawtopo )
    else if ( mksrf_frawtopo /= ' ') then
       write(6,*) 'MKGRID namelist error, mksrf_frawtopo used without a grid file'
       call abort()
    endif

    if ( nfile == 0 )then
       write(6,*) 'MKGRID namelist error, must specify at least one grid file'
       call abort()
    else if ( nfile > 1 )then
       write(6,*) 'More than one file selected -- files used with this precedence:'
       write(6,*) ' first cam files'
       write(6,*) ' second clm files'
       write(6,*) ' third navy oro'
       write(6,*) ' fourth cesm-domain'
       write(6,*) ' last rawtopo'
       if ( nfile > nused )then
          write(6,*) 'More files than can be used'
          call abort()
       end if
    endif


    !--- final comments ---

    lsmlon = ldomain%ni
    lsmlat = ldomain%nj

    write (resol,'(i4.4,"x",i4.4)') lsmlat,lsmlon
    fgriddat = './griddata_'//trim(resol)//'.nc'
    ffracdat = './fracdata_'//trim(resol)//'.nc'
    ftopodat = './topodata_'//trim(resol)//'.nc'

    write(6,*) 'Write land grid file:', trim(fgriddat)
    call mkfile(lsmlon, lsmlat, fgriddat, fileinfo, itype=1)
    if ( writeLfrc )then
        write(6,*) 'Write land frac file:', trim(ffracdat)
        call mkfile(lsmlon, lsmlat, ffracdat, Fracfileinfo, itype=2)
    end if
    if ( writeTopo )then
        write(6,*) 'Write land topo file:', trim(ftopodat)
        call mkfile(lsmlon, lsmlat, ftopodat, Topofileinfo, itype=3)
    end if

    call write_domain(ldomain,fgriddat, itype=1)
    if ( writeLfrc ) call write_domain(ldomain,ffracdat, itype=2)
    if ( writeTopo ) call write_domain(ldomain,ftopodat, itype=3)

    write (6,'(72a1)') ("-",i=1,60)
    write (6,'(a46,f5.1,a4,f5.1,a5)') 'Successfully created land model grid data set for ', &
         360./lsmlon,' by ',180./lsmlat,' grid'

contains

subroutine updateFileInfo( filename )
!
! Update fileinfo information when a file is being processed
!
  character(len=*), intent(IN) :: filename

  if ( .not. grid )then
     fileinfo = filename
  end if

  grid  = .true.
  nfile = nfile + 1
  if ( writeTopo1 )then
     Topofileinfo = filename
     writeTopo    = .true.
  end if
  if ( writeLFrc1 )then
     Fracfileinfo = filename
     writeLFrc    = .true.
  end if

  nused = 1
  if ( writeTopo .and. trim(fileinfo) /= trim(topofileinfo) ) nused = nused + 1
  if ( writeLFrc .and. trim(fileinfo) /= trim(fracfileinfo) ) nused = nused + 1
  if ( nused == 3 .and. trim(topofileinfo) == trim(fracfileinfo) ) nused = 2

end subroutine updateFileInfo

end program mkgriddata
