program mkgriddata

!-----------------------------------------------------------------------
!BOP
!
! !PROGRAM: mgriddata
!
! !DESCRIPTION:
! Creates land model surface dataset from original "raw" data files.
! Surface dataset contains model grid, pfts, inland water, glacier,
! soil texture, soil color, LAI and SAI and urban fraction.
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use shr_sys_mod  , only : shr_sys_getenv
    use fileutils    , only : getfil, putfil, opnfil, getavu, get_filename
    use creategridMod, only : creategrid
    use mkfileMod    , only : mkfile
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
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: lsmlon, lsmlat       ! clm grid resolution
    integer  :: i,j,k,m              ! indices
    integer  :: ier                  ! error status
    integer  :: ncid                 ! netCDF id
    integer  :: omode                ! netCDF output mode
    integer  :: varid                ! netCDF variable id
    real(r8) :: lsmedge(4)           ! North,East,South,West edges of grid (deg)
    character(len=256) :: fgriddat   ! grid data file name
    character(len=256) :: loc_fn     ! local file name
    character(len=  7) :: resol      ! resolution for file name
    character(len=32) :: subname = 'mkgriddata'  ! program name

    namelist /clmexp/    &
         mksrf_fnavyoro, &
         mksrf_lsmlon  , &
         mksrf_lsmlat  , &
         mksrf_edgen   , &
         mksrf_edgee   , &
         mksrf_edges   , &
         mksrf_edgew   
!-----------------------------------------------------------------------

    write(6,*) 'Attempting to initialize control settings .....'

    read(5, clmexp, iostat=ier)
    if (ier /= 0) then
       write(6,*)'error: namelist input resulted in error code ',ier
       call abort()
    endif

    if (mksrf_lsmlon==0 .and. mksrf_lsmlat==0) then
       write(6,*)'must specify mksrf_lsmlon and mksrf_lsmlat if are creating grid'
       stop
    else
       lsmlon     = mksrf_lsmlon
       lsmlat     = mksrf_lsmlat
       lsmedge(1) = mksrf_edgen
       lsmedge(2) = mksrf_edgee
       lsmedge(3) = mksrf_edges
       lsmedge(4) = mksrf_edgew 
    end if

    write (6,*) 'Attempting to create surface boundary data .....'
    write (6,'(72a1)') ("-",i=1,60)

    ! Create grid data 

    call creategrid(mksrf_lsmlon, mksrf_lsmlat, lsmedge, mksrf_fnavyoro)

    ! Create netCDF surface dataset.  

    write (resol,'(i3.3,"x",i3.3)') lsmlon,lsmlat
    fgriddat = './grid-data.'//trim(resol)//'.nc'
    call check_ret(nf_create(trim(fgriddat), nf_clobber, ncid), subname)

    ! File will be in define mode. Set fill mode to "no fill" to optimize performance

    call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

    ! Define dimensions and global attributes

    call mkfile(lsmlon, lsmlat, ncid)

    ! Write fields other than lai, sai, and heights to netcdf surface dataset

    call ncd_ioglobal(varname='EDGEN'   , data=lsmedge(1), ncid=ncid, flag='write')
    call ncd_ioglobal(varname='EDGEE'   , data=lsmedge(2), ncid=ncid, flag='write')
    call ncd_ioglobal(varname='EDGES'   , data=lsmedge(3), ncid=ncid, flag='write')
    call ncd_ioglobal(varname='EDGEW'   , data=lsmedge(4), ncid=ncid, flag='write')
    call ncd_ioglobal(varname='LONGXY'  , data=longxy    , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='LATIXY'  , data=latixy    , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='LANDMASK', data=landmask  , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='LANDFRAC', data=landfrac  , ncid=ncid, flag='write')

    ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

    call check_ret(nf_sync(ncid), subname)

    ! Close grid data dataset

    call check_ret(nf_close(ncid), subname)

    write (6,'(72a1)') ("-",i=1,60)
    write (6,'(a46,f5.1,a4,f5.1,a5)') 'land model surface data set successfully created for ', &
         360./lsmlon,' by ',180./lsmlat,' grid'

end program mkgriddata
