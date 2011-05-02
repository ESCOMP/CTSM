module surfrdMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: surfrdMod
!
! !DESCRIPTION:
! Contains methods for reading in surface data file and determining
! two-dimensional subgrid weights as well as writing out new surface
! dataset. When reading in the surface dataset, determines array
! which sets the PFT for each of the [maxpatch] patches and
! array which sets the relative abundance of the PFT.
! Also fills in the PFTs for vegetated portion of each grid cell.
! Fractional areas for these points pertain to "vegetated"
! area not to total grid area. Need to adjust them for fraction of grid
! that is vegetated. Also fills in urban, lake, wetland, and glacier patches.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use clm_varpar  , only : lsmlon, lsmlat, nlevsoi, numpft, &
                           maxpatch_pft, numcft, maxpatch, &
                           npatch_urban, npatch_lake, npatch_wet, npatch_glacier, &
                           maxpatch_urb, npatch_glacier_mec
  use clm_varctl  , only : create_glacier_mec_landunit, &
                           iulog, scmlat, scmlon, single_column
  use clm_varsur  , only : wtxy, vegxy, topoxy, pctspec
  use decompMod   , only : get_proc_bounds,gsMap_lnd_gdc2glo,ldecomp
  use clmtype
  use spmdMod                         
  use ncdio_pio
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: surfrd_get_latlon  ! Read surface dataset into domain (before domain decomp)
  public :: surfrd_get_grid    ! Read surface dataset into domain (after domain decomp)
  public :: surfrd_get_frac    ! Read land fraction into domain
  public :: surfrd_get_topo    ! Read topography into domain
  public :: surfrd             ! Read surface dataset and determine subgrid weights
!
! !PUBLIC DATA MEMBERS:
  logical, public :: crop_prog = .false. ! If prognostic crops is turned on
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! Updated by T Craig
!
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: surfrd_wtxy_special
  private :: surfrd_wtxy_veg_rank
  private :: surfrd_wtxy_veg_all
  private :: surfrd_wtxy_veg_dgvm
  private :: surfrd_mkrank
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_get_latlon
!
! !INTERFACE:
  subroutine surfrd_get_latlon(latlon,filename,mask,mfilename,pftmflag)
!
! !DESCRIPTION:
! Read the surface dataset grid related information:
! This is the first routine called by clm_initialize and no domain decomposition
! has been set yet
!
! !USES:
    use clm_varcon, only : spval
    use domainMod , only : latlon_type, latlon_init
    use fileutils , only : getfil
!
! !ARGUMENTS:
    implicit none
    type(latlon_type)         ,intent(inout) :: latlon   ! domain to init
    character(len=*)          ,intent(in)    :: filename ! grid filename
    integer,pointer  ,optional               :: mask(:)
    character(len=*) ,optional,intent(in)    :: mfilename ! grid filename
    logical          ,optional,intent(in)    :: pftmflag   ! is mask pft mask?
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: dimid,varid         ! netCDF id's
    integer :: ni,nj               ! size of grid on file
    integer :: n                   ! index
    integer :: ier                 ! error status
    type(file_desc_t) :: ncid      ! netcdf id
    type(file_desc_t) :: ncidm     ! netcdf id
    type(var_desc_t)  :: vardesc   ! variable descriptor
    character(len=256)  :: locfn   ! local file name
    real(r8),pointer :: rdata(:,:) ! temporary data
    integer ,pointer :: idata(:,:)
    logical :: NSEWset             ! true if lat/lon NSEW read from grid file
    logical :: EDGEset             ! true if EDGE NSEW read from grid file
    logical :: lpftmflag           ! is mask a pft mask, local copy
    logical :: readvar             ! read variable in or not
    logical :: dimexists           ! if dimension exists or not
    character(len=32) :: subname = 'surfrd_get_latlon'     ! subroutine name
!DEBUG
    integer :: i,j	
!DEBUG
!-----------------------------------------------------------------------

    NSEWset = .false.
    EDGEset = .false.
    lpftmflag = .false.

    ni = 0
    nj = 0

    if (present(pftmflag)) then
       lpftmflag = pftmflag
    endif

    if (masterproc) then
       if (filename == ' ') then
          write(iulog,*) trim(subname),' ERROR: filename must be specified '
          call endrun()
       endif
    end if

    call getfil( filename, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    if (single_column) then
       ni = lsmlon
       nj = lsmlat
    else
! Set PIO Mark
       call ncd_inqdid(ncid,'lon',dimid,dimexists)
       if (dimexists) call ncd_inqdlen(ncid,dimid,ni)
       call ncd_inqdid(ncid,'lat',dimid,dimexists)
       if (dimexists) call ncd_inqdlen(ncid,dimid,nj)
       call ncd_inqdid(ncid,'lsmlon',dimid,dimexists)
       if (dimexists) call ncd_inqdlen(ncid,dimid,ni)
       call ncd_inqdid(ncid,'lsmlat',dimid,dimexists)
       if (dimexists) call ncd_inqdlen(ncid,dimid,nj)
       if (ni == 0 .or. nj == 0) then
          write(iulog,*) trim(subname),' ERROR: ni or nj not set',ni,nj
          call endrun()
       end if
    endif
    
    call latlon_init(latlon,ni,nj)
    nullify(idata)
    if (present(mask)) then
       allocate(mask(ni*nj))
       allocate(idata(ni,nj))	
       mask = 1
       idata = 1	
    endif

    allocate(rdata(ni,nj))

    call ncd_io(ncid=ncid, varname='LONGXY', data=rdata, flag='read', readvar=readvar)
    if ( .not. readvar) call endrun( trim(subname)//' ERROR: LONGXY NOT on file' )
    latlon%lonc(:) = rdata(:,1)

    call ncd_io(ncid=ncid, varname='LATIXY', data=rdata, flag='read',readvar=readvar)
    if ( .not. readvar) call endrun( trim(subname)//' ERROR: LONGXY NOT on file' )
    latlon%latc(:) = rdata(1,:)

    if (single_column) then
       EDGEset = .true.
       latlon%edges(1) = 90.0_r8
       latlon%edges(2) = latlon%lonc(1) + 180._r8
       latlon%edges(3) = -90.0_r8
       latlon%edges(4) = latlon%lonc(1) - 180._r8
    else
       latlon%edges(:) = spval
       call ncd_inqvid(ncid,'EDGEN',varid,vardesc,readvar)
       if (readvar)then
          EDGEset = .true.
          call ncd_io(ncid=ncid,varname='EDGEN',data=latlon%edges(1),flag='read',readvar=readvar)
          call ncd_io(ncid=ncid,varname='EDGEE',data=latlon%edges(2),flag='read',readvar=readvar)
          call ncd_io(ncid=ncid,varname='EDGES',data=latlon%edges(3),flag='read',readvar=readvar)
          call ncd_io(ncid=ncid,varname='EDGEW',data=latlon%edges(4),flag='read',readvar=readvar)
          if (maxval(latlon%edges) > 1.0e35) EDGEset = .false. ! read garbage
       endif
    endif

    call ncd_inqvid(ncid,'LATN',varid,vardesc,readvar)
    if (readvar)then
       NSEWset = .true.
       call ncd_io(ncid=ncid, varname='LATN',data=rdata,flag='read')
       latlon%latn(:) = rdata(1,:)
       call ncd_io(ncid=ncid, varname='LONE',data=rdata,flag='read')
       latlon%lone(:) = rdata(:,1)
       call ncd_io(ncid=ncid, varname='LATS',data=rdata,flag='read')
       latlon%lats(:) = rdata(1,:)
       call ncd_io(ncid=ncid, varname='LONW',data=rdata,flag='read')
       latlon%lonw(:) = rdata(:,1)
    endif
    
    if (present(mask)) then
       if (present(mfilename)) then
          if (mfilename == ' ') then
             write(iulog,*) trim(subname),' ERROR: mfilename must be specified '
             call endrun()
          endif
          call getfil( mfilename, locfn, 0 )
          call ncd_pio_openfile (ncidm, trim(locfn), 0)
       else
          ncidm = ncid
       endif
       
       mask = 1
       call ncd_inqvid(ncidm,'LANDMASK',varid,vardesc,readvar)
       if (readvar)then
	  ! ASSUME that LANDMASK is 2d here !TODO? talk to jim about this
          call ncd_io(ncid=ncidm, varname='LANDMASK', data=idata, flag='read')
          ! mask = reshape(idata, (/1/)) !TODO? why does this not work?
	  do j = 1,nj
          do i = 1,ni
             n = (j-1)*ni + i	
             mask(n) = idata(i,j)
          enddo
          enddo
       endif
       
       !--- if this is a pft mask, then modify and look for pftdata_mask array on dataset ---
       if (lpftmflag) then
          do n = 1,ni*nj
             if (mask(n) <= 0) mask(n) = -1
          enddo
          call ncd_inqvid(ncidm,'PFTDATA_MASK',varid,vardesc,readvar)
          if (readvar)then
             call ncd_io(ncid=ncidm, varname='PFTDATA_MASK', data=idata, flag='read')
             do j = 1,nj
             do i = 1,ni
                n = (j-1)*ni + i	
                mask(n) = idata(i,j)
             enddo
             enddo
          endif
       endif
       
       if (present(mfilename)) then
          call ncd_pio_closefile(ncidm)
       endif
    endif

    deallocate(rdata)
    if(associated(idata)) deallocate(idata)
        
    call ncd_pio_closefile(ncid)

  end subroutine surfrd_get_latlon

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_get_grid
!
! !INTERFACE:
  subroutine surfrd_get_grid(domain,filename,beg,end,clmlevel)
!
! !DESCRIPTION:
! This is called after the domain decomposition has been created
! Read the surface dataset grid related information:
! o real edges of grid
! o integer  number of longitudes per latitude
! o real latitude  of grid cell (degrees)
! o real longitude of grid cell (degrees)
! If grid is read in from dataset, grid is assumed to be global but
! does not have to be regular.
! If grid is generated by model, grid does not have to be global but 
! must then define the north, east, south, and west edges:
!
! o edges(1)    = northern edge of grid (degrees): >  -90 and <= 90
! o edges(2)    = eastern edge of grid (degrees) : see following notes
! o edges(3)    = southern edge of grid (degrees): >= -90 and <  90
! o edges(4)    = western edge of grid (degrees) : see following notes
!
!   For partial grids, northern and southern edges are any latitude
!   between 90 (North Pole) and -90 (South Pole). Western and eastern
!   edges are any longitude between -180 and 180, with longitudes
!   west of Greenwich negative. That is, western edge >= -180 and < 180;
!   eastern edge > western edge and <= 180.
!
!   For global grids, northern and southern edges are 90 (North Pole)
!   and -90 (South Pole). The western and eastern edges depend on
!   whether the grid starts at Dateline or Greenwich. Regardless,
!   these edges must span 360 degrees. Examples:
!
!                              West edge    East edge
!                            --------------------------------------------------
!  (1) Dateline            :        -180 to 180       (negative W of Greenwich)
!  (2) Greenwich (centered):    0 - dx/2 to 360 - dx/2
!
!    Grid 2 is the grid for cam and cesm mode since the NCAR CAM
!    starts at Greenwich, centered on Greenwich
!
! !USES:
    use clm_varcon, only : spval
    use domainMod , only : domain_type,domain_init
    use fileutils , only : getfil
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout) :: domain   ! domain to init
    character(len=*) ,intent(in)    :: filename ! grid filename
    integer          ,intent(in)    :: beg      ! local beg index
    integer          ,intent(in)    :: end      ! local end index
    character(len=*) ,intent(in)    :: clmlevel ! type of grid
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    type(file_desc_t) :: ncid               ! netcdf id
    integer :: ni,nj                        ! size of grid on file
    integer :: dimid,varid                  ! netCDF id's
    integer :: ier,ret                      ! error status
    logical :: readvar 
    character(len=256):: locfn              ! local file name
    character(len=32) :: subname = 'surfrd_get_grid'     ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then
       if (filename == ' ') then
          write(iulog,*) trim(subname),' ERROR: filename must be specified '
          call endrun()
       endif
    end if

    call getfil( filename, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    if (single_column) then
       ni = lsmlon
       nj = lsmlat
    else
       call ncd_inqdid (ncid, 'lsmlon', dimid)
       call ncd_inqdlen(ncid, dimid, ni)
       call ncd_inqdid (ncid, 'lsmlat', dimid)
       call ncd_inqdlen(ncid, dimid, nj)
    endif

    call domain_init(domain,ni,nj,beg,end,clmlevel)

    call ncd_io(ncid=ncid, varname= 'AREA', flag='read', data=domain%area, &
         dim1name=clmlevel, readvar=readvar)
    if (readvar) domain%areaset = .true.

    call ncd_io(ncid=ncid, varname= 'LONGXY', flag='read', data=domain%lonc, &
         dim1name=clmlevel, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: LONGXY NOT on file' )

    call ncd_io(ncid=ncid, varname= 'LATIXY', flag='read', data=domain%latc, &
         dim1name=clmlevel, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: LATIXY NOT on file' )

    ! set mask to 1 everywhere by default, override if LANDMASK exists
    ! if landmask exists, use it to set pftm (for older datasets)
    ! pftm should be overwritten below for newer datasets

    domain%mask = 1
    call ncd_io(ncid=ncid, varname= 'LANDMASK', flag='read', data=domain%mask, &
         dim1name=clmlevel, readvar=readvar)
    domain%pftm = domain%mask
    where (domain%mask <= 0)
       domain%pftm = -1
    endwhere

    call ncd_io(ncid=ncid, varname= 'PFTDATA_MASK', flag='read', data=domain%pftm, &
         dim1name=clmlevel, readvar=readvar)

    call ncd_pio_closefile(ncid)

  end subroutine surfrd_get_grid

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd
!
! !INTERFACE:
!  subroutine surfrd(vegxy, wtxy, lfsurdat, domain)
  subroutine surfrd(lfsurdat, domain)
!
! !DESCRIPTION:
! Read the surface dataset and create subgrid weights.
! The model's surface dataset recognizes 6 basic land cover types within a grid
! cell: lake, wetland, urban, glacier, glacier_mec and vegetated. The vegetated
! portion of the grid cell is comprised of up to [maxpatch_pft] PFTs. These
! subgrid patches are read in explicitly for each grid cell. This is in
! contrast to LSMv1, where the PFTs were built implicitly from biome types.
!    o real edges of grid
!    o integer  number of longitudes per latitude
!    o real latitude  of grid cell (degrees)
!    o real longitude of grid cell (degrees)
!    o integer surface type: 0 = ocean or 1 = land
!    o integer soil color (1 to 20) for use with soil albedos
!    o real soil texture, %sand, for thermal and hydraulic properties
!    o real soil texture, %clay, for thermal and hydraulic properties
!    o real % of cell covered by lake    for use as subgrid patch
!    o real % of cell covered by wetland for use as subgrid patch
!    o real % of cell that is urban      for use as subgrid patch
!    o real % of cell that is glacier    for use as subgrid patch
!    o real % of cell that is glacier_mec for use as subgrid patch
!    o integer PFTs
!    o real % abundance PFTs (as a percent of vegetated area)
!
! !USES:
    use clm_varctl  , only : allocate_all_vegpfts, create_crop_landunit
    use pftvarcon   , only : noveg
    use fileutils   , only : getfil
    use domainMod   , only : domain_type
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: lfsurdat               ! surf filename
    type(domain_type),intent(in) :: domain ! domain associated with wtxy
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein, Sam Levis and Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=256):: locfn                          ! local file name
    type(file_desc_t) :: ncid                           ! netcdf id
    integer           :: begg,endg   
    character(len=32) :: subname = 'surfrd'             ! subroutine name
!-----------------------------------------------------------------------

    call get_proc_bounds(begg,endg)
    allocate(pctspec(begg:endg))

    vegxy(:,:) = noveg
    wtxy(:,:)  = 0._r8
    pctspec(:) = 0._r8
    if (allocated(topoxy)) topoxy(:,:) = 0._r8

    if (masterproc) then
       write(iulog,*) 'Attempting to read surface boundary data .....'
       if (lfsurdat == ' ') then
          call endrun( trim(subname)//' ERROR: lfsurdat must be specified' )
       endif
    endif

    ! Read surface data

    call getfil( lfsurdat, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Obtain surface dataset special landunit info

    call surfrd_wtxy_special(ncid, domain)

    ! Obtain surface dataset vegetated landunit info

#if (defined CNDV)
    if (create_crop_landunit) then ! CNDV means allocate_all_vegpfts = .true.
       call surfrd_wtxy_veg_all(ncid, domain)
    end if
    call surfrd_wtxy_veg_dgvm(domain)
#else
    if (allocate_all_vegpfts) then
       call surfrd_wtxy_veg_all(ncid, domain)
    else
       call surfrd_wtxy_veg_rank(ncid, domain)
    end if
#endif

#ifdef CROP
    if ( .not. crop_prog )then
       call endrun( trim(subname)//' ERROR: surface dataset MUST have '// & 
                    'prognostic crop on it if CROP #ifdef is set'//       &
                    ' -- use a different surface dataset' )
    end if
#else
    if ( crop_prog )then
       call endrun( trim(subname)//' ERROR: surface dataset can NOT have '// & 
                    'prognostic crop on it if CROP #ifdef is NOT set'//      &
                    ' -- use a different surface dataset' )
    end if
#endif

    call ncd_pio_closefile(ncid)
    if ( masterproc )then
       write(iulog,*) 'Successfully read surface boundary data'
       write(iulog,*)
    end if

  end subroutine surfrd

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_get_frac
!
! !INTERFACE:
  subroutine surfrd_get_frac(domain,filename,glcfilename)
!
! !DESCRIPTION:
! Read the landfrac dataset grid related information:
! Assume domain has already been initialized and read
!
! !USES:
    use domainMod , only : domain_type
    use fileutils , only : getfil
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout) :: domain   ! domain to init
    character(len=*) ,intent(in)    :: filename ! grid filename
    character(len=*) ,optional, intent(in) :: glcfilename ! glc mask filename
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    type(file_desc_t) :: ncid      ! netcdf file id
    type(file_desc_t) :: ncidg     ! netCDF id for glcmask
    integer :: n                   ! indices
    integer :: ni,nj,ns            ! size of grid on file
    integer :: dimid,varid         ! netCDF id's
    integer :: ier                 ! error status
    real(r8):: eps = 1.0e-12_r8    ! lat/lon error tolerance
    integer :: beg,end             ! local beg,end indices
    character(len=8)   :: clmlevel         ! grid type
    logical            :: readvar          ! is variable in file
    real(r8),pointer   :: lonc(:),latc(:)  ! local lat/lon
    character(len=256) :: locfn            ! local file name
    character(len=32) :: subname = 'surfrd_get_frac'     ! subroutine name
!-----------------------------------------------------------------------

    if (filename == ' ') then
       write(iulog,*) trim(subname),' ERROR: filename must be specified '
       call endrun()
    endif
    
    call getfil( filename, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    if (single_column) then
       ni = lsmlon
       nj = lsmlat
    else
       call ncd_inqdid (ncid, 'lsmlon', dimid)
       call ncd_inqdlen(ncid, dimid, ni)
       call ncd_inqdid (ncid, 'lsmlat', dimid)
       call ncd_inqdlen(ncid, dimid, nj)
    endif
    
    ns = ni*nj
    if (domain%ni /= ni .or. domain%nj /= nj .or. domain%ns /= ns) then
       write(iulog,*) trim(subname),' ERROR: landfrac file mismatch ni,nj',&
            domain%ni,ni,domain%nj,nj,domain%ns,ns
       call endrun()
    endif
    
    ni  = domain%ni
    nj  = domain%nj
    beg = domain%nbeg
    end = domain%nend
    clmlevel = domain%clmlevel

    allocate(latc(beg:end),lonc(beg:end))

    call ncd_io(ncid=ncid, varname='LONGXY', flag='read', data=lonc, &
         dim1name=clmlevel, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: LONGXY NOT on fracdata file' )

    call ncd_io(ncid=ncid, varname='LATIXY', flag='read', data=latc, &
         dim1name=clmlevel, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: LATIXY NOT on fracdata file' )

    do n = beg,end
       if (abs(latc(n)-domain%latc(n)) > eps .or. &
           abs(lonc(n)-domain%lonc(n)) > eps) then
          write(iulog,*) trim(subname),' ERROR: landfrac file mismatch lat,lon',&
               latc(n),domain%latc(n),lonc(n),domain%lonc(n),eps
          call endrun()
       endif
    enddo
          
    call ncd_io(ncid=ncid, varname='LANDMASK', flag='read', data=domain%mask, &
         dim1name=clmlevel, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: LANDMASK NOT on fracdata file' )

    call ncd_io(ncid=ncid, varname='LANDFRAC', flag='read', data=domain%frac, &
         dim1name=clmlevel, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: LANDFRAC NOT on fracdata file' )

    deallocate(latc,lonc)

    call ncd_pio_closefile(ncid)

    if (present(glcfilename)) then
       if (masterproc) then
          if (glcfilename == ' ') then
             write(iulog,*) trim(subname),' ERROR: glc filename must be specified '
             call endrun()
          endif
       end if
       call getfil( glcfilename, locfn, 0 )
       call ncd_pio_openfile (ncidg, trim(locfn), 0)

       domain%glcmask = 0
       call ncd_io(ncid=ncidg, varname='GLCMASK', flag='read', data=domain%glcmask, &
            dim1name=clmlevel, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: GLCMASK NOT in file' )

       call ncd_pio_closefile(ncidg)
    endif   ! present(glcfilename)

  end subroutine surfrd_get_frac

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_get_topo
!
! !INTERFACE:
  subroutine surfrd_get_topo(domain,filename)
!
! !DESCRIPTION:
! Read the topo dataset grid related information:
! Assume domain has already been initialized and read
!
! !USES:
    use domainMod , only : domain_type
    use fileutils , only : getfil
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout) :: domain   ! domain to init
    character(len=*) ,intent(in)    :: filename ! grid filename
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    type(file_desc_t) :: ncid      ! netcdf file id
    integer :: n                   ! indices
    integer :: ni,nj,ns            ! size of grid on file
    integer :: dimid,varid         ! netCDF id's
    integer :: ier                 ! error status
    real(r8):: eps = 1.0e-12_r8    ! lat/lon error tolerance
    integer :: beg,end             ! local beg,end indices
    character(len=8)    :: clmlevel   ! grid type
    real(r8),pointer    :: lonc(:),latc(:)  ! local lat/lon
    character(len=256)  :: locfn   ! local file name
    logical :: readvar             ! is variable on file
    character(len=32) :: subname = 'surfrd_get_topo'     ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then
       if (filename == ' ') then
          write(iulog,*) trim(subname),' ERROR: filename must be specified '
          call endrun()
       endif
    end if

    call getfil( filename, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    if (single_column) then
       ni = lsmlon
       nj = lsmlat
    else
       call ncd_inqdid (ncid, 'lsmlon', dimid)
       call ncd_inqdlen(ncid, dimid, ni)
       call ncd_inqdid (ncid, 'lsmlat', dimid)
       call ncd_inqdlen(ncid, dimid, nj)
    endif
    
    ns = ni*nj
    if (domain%ni /= ni .or. domain%nj /= nj .or. domain%ns /= ns) then
       write(iulog,*) trim(subname),' ERROR: topo file mismatch ni,nj',domain%ni,ni,domain%nj,nj,domain%ns,ns
       call endrun()
    endif
    
    ni  = domain%ni
    nj  = domain%nj
    beg = domain%nbeg
    end = domain%nend
    clmlevel = domain%clmlevel

    allocate(latc(beg:end),lonc(beg:end))

    call ncd_io(ncid=ncid, varname='LONGXY', flag='read', data=lonc, &
         dim1name=clmlevel, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: LONGXY  NOT on topodata file' )

    call ncd_io(ncid=ncid, varname='LATIXY', flag='read', data=latc, &
         dim1name=clmlevel, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: LONGXY  NOT on topodata file' )

    do n = beg,end
       if (abs(latc(n)-domain%latc(n)) > eps .or. &
           abs(lonc(n)-domain%lonc(n)) > eps) then
          write(iulog,*) trim(subname),' ERROR: topo file mismatch lat,lon',latc(n),&
               domain%latc(n),lonc(n),domain%lonc(n),eps
          call endrun()
       endif
    enddo

    call ncd_io(ncid=ncid, varname='TOPO', flag='read', data=domain%topo, &
         dim1name=clmlevel, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: LONGXY  NOT on topodata file' )

    deallocate(latc,lonc)

    call ncd_pio_closefile(ncid)

  end subroutine surfrd_get_topo

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_wtxy_special 
!
! !INTERFACE:
  subroutine surfrd_wtxy_special(ncid, domain)
!
! !DESCRIPTION:
! Determine weight with respect to gridcell of all special "pfts" as well
! as soil color and percent sand and clay
!
! !USES:
    use pftvarcon     , only : noveg
    use UrbanInputMod , only : urbinp
    use domainMod     , only : domain_type
    use clm_varctl    , only : glc_nec, glc_topomax
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    type(domain_type), intent(in) :: domain ! domain associated with wtxy
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein, Sam Levis and Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: n,nl,ns,nurb,g             ! indices
    integer  :: begg,endg                  ! gcell beg/end
    integer  :: dimid,varid                ! netCDF id's
    real(r8) :: nlevsoidata(nlevsoi)
    logical  :: found                      ! temporary for error check
    integer  :: nindx                      ! temporary for error check
    integer  :: ier                        ! error status
    integer  :: nlev                       ! level
    integer  :: npatch
    logical  :: readvar
    real(r8),pointer :: pctgla(:)      ! percent of grid cell is glacier
    real(r8),pointer :: pctlak(:)      ! percent of grid cell is lake
    real(r8),pointer :: pctwet(:)      ! percent of grid cell is wetland
    real(r8),pointer :: pcturb(:)      ! percent of grid cell is urbanized
    real(r8),pointer :: pctglc_mec(:,:)   ! percent of grid cell is glacier_mec (in each elev class)
    real(r8),pointer :: pctglc_mec_tot(:) ! percent of grid cell is glacier (sum over classes)
    real(r8),pointer :: topoglc_mec(:,:)  ! surface elevation in each elev class
    character(len=32) :: subname = 'surfrd_wtxy_special'  ! subroutine name
    real(r8) closelat,closelon
!!-----------------------------------------------------------------------

    ns = domain%ns
    call get_proc_bounds(begg,endg)

    allocate(pctgla(begg:endg),pctlak(begg:endg))
    allocate(pctwet(begg:endg),pcturb(begg:endg))
    if (create_glacier_mec_landunit) then
       allocate(pctglc_mec(begg:endg,glc_nec))
       allocate(pctglc_mec_tot(begg:endg))
       allocate(topoglc_mec(begg:endg,glc_nec))
    endif

    call check_dim(ncid, 'nlevsoi', nlevsoi)

       ! Obtain non-grid surface properties of surface dataset other than percent pft

    call ncd_io(ncid=ncid, varname='PCT_WETLAND', flag='read', data=pctwet, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_WETLAND  NOT on surfdata file' )

    call ncd_io(ncid=ncid, varname='PCT_LAKE'   , flag='read', data=pctlak, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_LAKE NOT on surfdata file' )

    call ncd_io(ncid=ncid, varname='PCT_GLACIER', flag='read', data=pctgla, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_GLACIER NOT on surfdata file' )

    call ncd_io(ncid=ncid, varname='PCT_URBAN'  , flag='read', data=pcturb, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_URBAN NOT on surfdata file' )

    if (create_glacier_mec_landunit) then          ! call ncd_io_gs_int2d

       call check_dim(ncid, 'nglcec',   glc_nec   )
       call check_dim(ncid, 'nglcecp1', glc_nec+1 )

       call ncd_io(ncid=ncid, varname='GLC_MEC', flag='read', data=glc_topomax, &
            readvar=readvar)
       if ( .not. readvar) call endrun( trim(subname)//'ERROR: GLC_MEC was NOT on the input surfdata file' )

       call ncd_io(ncid=ncid, varname='PCT_GLC_MEC', flag='read', data=pctglc_mec, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_GLC_MEC NOT on surfdata file' )

       call ncd_io(ncid=ncid, varname='TOPO_GLC_MEC',  flag='read', data=topoglc_mec, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: TOPO_GLC_MEC NOT on surfdata file' )

       pctglc_mec_tot(:) = 0._r8
       do n = 1, glc_nec
          do nl = begg,endg
             pctglc_mec_tot(nl) = pctglc_mec_tot(nl) + pctglc_mec(nl,n)
          enddo
       enddo

       ! Make sure sum of pctglc_mec = pctgla, then zero out pctgla
       ! (assumes glc_mec values are double precision)

       do nl = begg,endg
          if (abs(pctgla(nl) - pctglc_mec_tot(nl)) > 1.0e-11) then
             write(iulog,*) ' '
             write(iulog,*) 'surfrd error: pctgla not equal to sum of pctglc_mec for nl=', nl
             write(iulog,*) 'pctgla =', pctgla(nl)
             write(iulog,*) 'pctglc_mec_tot =', pctglc_mec_tot(nl)
             call endrun()
          endif
          pctgla(nl) = 0._r8
       enddo

       ! If pctglc_mec_tot is very close to 100%, round to 100%

       do nl = begg,endg
          if (abs(pctglc_mec_tot(nl) - 100._r8) < 1.0e-8) then
             pctglc_mec(nl,:) = pctglc_mec(nl,:) * 100._r8 / pctglc_mec_tot(nl)
             pctglc_mec_tot(nl) = 100._r8
          endif
       enddo

       pctspec = pctwet + pctlak + pcturb + pctglc_mec_tot

       if ( masterproc ) write(iulog,*) '   elevation limits = ', glc_topomax

    else

       pctspec = pctwet + pctlak + pcturb + pctgla
 
    endif

    ! Error check: glacier, lake, wetland, urban sum must be less than 100

    found = .false.
    do nl = begg,endg
       if (pctspec(nl) > 100._r8+1.e-04_r8) then
          found = .true.
          nindx = nl
          exit
       end if
       if (found) exit
    end do
    if ( found ) then
       write(iulog,*)'surfrd error: PFT cover>100 for nl=',nindx
       call endrun()
    end if

    ! Determine veg and wtxy for special landunits

    do nl = begg,endg

       vegxy(nl,npatch_lake)   = noveg
       wtxy(nl,npatch_lake)    = pctlak(nl)/100._r8

       vegxy(nl,npatch_wet)    = noveg
       wtxy(nl,npatch_wet)     = pctwet(nl)/100._r8

       vegxy(nl,npatch_glacier)= noveg
       wtxy(nl,npatch_glacier) = pctgla(nl)/100._r8

       ! Initialize urban weights

       do nurb = npatch_urban, npatch_lake-1 
          vegxy(nl,nurb) = noveg
          wtxy(nl,nurb)  = pcturb(nl) / 100._r8
       end do
       if ( pcturb(nl) > 0.0_r8 )then
          wtxy(nl,npatch_urban)   = wtxy(nl,npatch_urban)*urbinp%wtlunit_roof(nl)
          wtxy(nl,npatch_urban+1) = wtxy(nl,npatch_urban+1)*(1 - urbinp%wtlunit_roof(nl))/3
          wtxy(nl,npatch_urban+2) = wtxy(nl,npatch_urban+2)*(1 - urbinp%wtlunit_roof(nl))/3
          wtxy(nl,npatch_urban+3) = wtxy(nl,npatch_urban+3)*(1 - urbinp%wtlunit_roof(nl))/3 * (1.-urbinp%wtroad_perv(nl))
          wtxy(nl,npatch_urban+4) = wtxy(nl,npatch_urban+4)*(1 - urbinp%wtlunit_roof(nl))/3 * urbinp%wtroad_perv(nl)
       end if

    end do

    ! Check to make sure we have valid urban data for each urban patch

    found = .false.
    do nl = begg,endg
       if ( pcturb(nl) > 0.0_r8 )then
         if (urbinp%canyon_hwr(nl)            .le. 0._r8 .or. &
             urbinp%em_improad(nl)            .le. 0._r8 .or. &
             urbinp%em_perroad(nl)            .le. 0._r8 .or. &
             urbinp%em_roof(nl)               .le. 0._r8 .or. &
             urbinp%em_wall(nl)               .le. 0._r8 .or. &
             urbinp%ht_roof(nl)               .le. 0._r8 .or. &
             urbinp%thick_roof(nl)            .le. 0._r8 .or. &
             urbinp%thick_wall(nl)            .le. 0._r8 .or. &
             urbinp%t_building_max(nl)        .le. 0._r8 .or. &
             urbinp%t_building_min(nl)        .le. 0._r8 .or. &
             urbinp%wind_hgt_canyon(nl)       .le. 0._r8 .or. &
             urbinp%wtlunit_roof(nl)          .le. 0._r8 .or. &
             urbinp%wtroad_perv(nl)           .le. 0._r8 .or. &
             any(urbinp%alb_improad_dir(nl,:) .le. 0._r8) .or. &
             any(urbinp%alb_improad_dif(nl,:) .le. 0._r8) .or. &
             any(urbinp%alb_perroad_dir(nl,:) .le. 0._r8) .or. &
             any(urbinp%alb_perroad_dif(nl,:) .le. 0._r8) .or. &
             any(urbinp%alb_roof_dir(nl,:)    .le. 0._r8) .or. &
             any(urbinp%alb_roof_dif(nl,:)    .le. 0._r8) .or. &
             any(urbinp%alb_wall_dir(nl,:)    .le. 0._r8) .or. &
             any(urbinp%alb_wall_dif(nl,:)    .le. 0._r8) .or. &
             any(urbinp%tk_roof(nl,:)         .le. 0._r8) .or. &
             any(urbinp%tk_wall(nl,:)         .le. 0._r8) .or. &
             any(urbinp%cv_roof(nl,:)         .le. 0._r8) .or. &
             any(urbinp%cv_wall(nl,:)         .le. 0._r8)) then
            found = .true.
            nindx = nl
            exit
         else
            if (urbinp%nlev_improad(nl) .gt. 0) then
               nlev = urbinp%nlev_improad(nl)
               if (any(urbinp%tk_improad(nl,1:nlev) .le. 0._r8) .or. &
                   any(urbinp%cv_improad(nl,1:nlev) .le. 0._r8)) then
                  found = .true.
                  nindx = nl
                  exit
               end if
            end if
         end if
         if (found) exit
       end if
    end do
    if ( found ) then
       write(iulog,*)'surfrd error: no valid urban data for nl=',nindx
       write(iulog,*)'canyon_hwr: ',urbinp%canyon_hwr(nindx)
       write(iulog,*)'em_improad: ',urbinp%em_improad(nindx)
       write(iulog,*)'em_perroad: ',urbinp%em_perroad(nindx)
       write(iulog,*)'em_roof: ',urbinp%em_roof(nindx)
       write(iulog,*)'em_wall: ',urbinp%em_wall(nindx)
       write(iulog,*)'ht_roof: ',urbinp%ht_roof(nindx)
       write(iulog,*)'thick_roof: ',urbinp%thick_roof(nindx)
       write(iulog,*)'thick_wall: ',urbinp%thick_wall(nindx)
       write(iulog,*)'t_building_max: ',urbinp%t_building_max(nindx)
       write(iulog,*)'t_building_min: ',urbinp%t_building_min(nindx)
       write(iulog,*)'wind_hgt_canyon: ',urbinp%wind_hgt_canyon(nindx)
       write(iulog,*)'wtlunit_roof: ',urbinp%wtlunit_roof(nindx)
       write(iulog,*)'wtroad_perv: ',urbinp%wtroad_perv(nindx)
       write(iulog,*)'alb_improad_dir: ',urbinp%alb_improad_dir(nindx,:)
       write(iulog,*)'alb_improad_dif: ',urbinp%alb_improad_dif(nindx,:)
       write(iulog,*)'alb_perroad_dir: ',urbinp%alb_perroad_dir(nindx,:)
       write(iulog,*)'alb_perroad_dif: ',urbinp%alb_perroad_dif(nindx,:)
       write(iulog,*)'alb_roof_dir: ',urbinp%alb_roof_dir(nindx,:)
       write(iulog,*)'alb_roof_dif: ',urbinp%alb_roof_dif(nindx,:)
       write(iulog,*)'alb_wall_dir: ',urbinp%alb_wall_dir(nindx,:)
       write(iulog,*)'alb_wall_dif: ',urbinp%alb_wall_dif(nindx,:)
       write(iulog,*)'tk_roof: ',urbinp%tk_roof(nindx,:)
       write(iulog,*)'tk_wall: ',urbinp%tk_wall(nindx,:)
       write(iulog,*)'cv_roof: ',urbinp%cv_roof(nindx,:)
       write(iulog,*)'cv_wall: ',urbinp%cv_wall(nindx,:)
       if (urbinp%nlev_improad(nindx) .gt. 0) then
          nlev = urbinp%nlev_improad(nindx)
          write(iulog,*)'tk_improad: ',urbinp%tk_improad(nindx,1:nlev)
          write(iulog,*)'cv_improad: ',urbinp%cv_improad(nindx,1:nlev)
       end if
       call endrun()
    end if

    ! Initialize glacier_mec weights

    if (create_glacier_mec_landunit) then

       do n = 1, glc_nec
          npatch = npatch_glacier_mec - glc_nec + n

          do nl = begg, endg
             vegxy (nl,npatch) = noveg
             wtxy  (nl,npatch) = pctglc_mec(nl,n)/100._r8
             topoxy(nl,npatch) = topoglc_mec(nl,n)

          enddo   ! nl
       enddo      ! glc_nec

       deallocate(pctglc_mec, pctglc_mec_tot, topoglc_mec)

   endif    ! create_glacier_mec_landunit

   deallocate(pctgla,pctlak,pctwet,pcturb)

  end subroutine surfrd_wtxy_special

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_wtxy_veg_rank
!
! !INTERFACE:
  subroutine surfrd_wtxy_veg_rank(ncid, domain)
!
! !DESCRIPTION:
! Determine wtxy and veg arrays for non-dynamic landuse mode
!
! !USES:
    use clm_varctl, only : create_crop_landunit
    use pftvarcon , only : crop, noveg, nirrig
    use domainMod , only : domain_type
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! netcdf id
    type(domain_type),intent(in) :: domain ! domain associated with wtxy
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein, Sam Levis and Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: k,m,k1,k2,n,nl,ns               ! indices
    integer  :: begg,endg                       ! beg/end gcell index
    integer  :: dimid,varid                     ! netCDF id's
    integer  :: start(3),count(3)               ! netCDF start/count arrays
    integer  :: cropcount                       ! temporary counter
    real(r8),allocatable :: sumvec(:)           ! temporary vector sum
    logical  :: found                           ! temporary for error check
    integer  :: nindx                           ! temporary for error check
    integer  :: miss = 99999                    ! missing data indicator
    real(r8) :: wst(0:numpft)                   ! as pft at specific i, j
    integer ,allocatable :: wsti(:)             ! ranked indices largest wst values
    real(r8) :: wst_sum                         ! sum of %pft
    real(r8) :: sumpct                          ! sum of %pft over maxpatch_pft
    real(r8) :: diff                            ! the difference (wst_sum - sumpct)
    real(r8) :: rmax                            ! maximum patch cover
    integer ,allocatable :: pft(:,:)            ! PFT
    integer ,allocatable :: cft(:,:)            ! CFT
    real(r8),allocatable :: pctcft_lunit(:,:)   ! % of crop lunit area for CFTs
    real(r8),allocatable :: pctpft_lunit(:,:)   ! % of vegetated lunit area PFTs
    integer  :: ier                             ! error status
    real(r8),allocatable :: pctpft(:,:)         ! percent of vegetated gridcell area for PFTs
    real(r8),pointer :: arrayl(:,:)             ! local array
    integer ,pointer :: irrayg(:)               ! global array
    real(r8),allocatable :: rmaxpatchdata(:)
    integer ,allocatable :: imaxpatchdata(:)
    real(r8) :: numpftp1data(0:numpft)         
    character(len=32) :: subname = 'surfrd_wtxy_veg_rank'  ! subroutine name
!-----------------------------------------------------------------------

    ns = domain%ns
    call get_proc_bounds(begg,endg)

    allocate(sumvec(begg:endg))
    allocate(cft(begg:endg,numcft))
    allocate(pft(begg:endg,maxpatch_pft))
    allocate(pctcft_lunit(begg:endg,numcft))
    allocate(pctpft_lunit(begg:endg,maxpatch_pft))
    allocate(pctpft(begg:endg,0:numpft))
    allocate(wsti(maxpatch_pft))

    call check_dim(ncid, 'lsmpft', numpft+1)
    allocate(arrayl(begg:endg,0:numpft))
    call ncd_io(ncid=ncid, varname='PCT_PFT', flag='read', data=arrayl, dim1name=grlnd)
    pctpft(begg:endg,0:numpft) = arrayl(begg:endg,0:numpft)
    deallocate(arrayl)

    ! 1. pctpft must go back to %vegetated landunit instead of %gridcell
    ! 2. pctpft bare = 100 when landmask = 1 and 100% special landunit
    ! NB: (1) and (2) do not apply to crops.
    ! For now keep all cfts instead of most dominant cfts

    do nl = begg,endg

       cft(nl,:) = 0
       pctcft_lunit(nl,:) = 0._r8
       cropcount = 0

       if (pctspec(nl) < 100._r8) then

          do m = 0, numpft
             if (create_crop_landunit) then
                ! Separate crop landunit is to be created

                if (crop(m) == 1._r8 .and. pctpft(nl,m) > 0._r8) then
                   cropcount = cropcount + 1
                   if (cropcount > numcft) then
                      write(iulog,*) 'ERROR surfrdMod: cropcount>numcft'
                      call endrun()
                   end if
                   cft(nl,cropcount) = m
                   pctcft_lunit(nl,cropcount) = pctpft(nl,m) * 100._r8/(100._r8-pctspec(nl))
                   pctpft(nl,m) = 0.0_r8
                else if (crop(m) == 0._r8) then
                   pctpft(nl,m) = pctpft(nl,m) * 100._r8/(100._r8-pctspec(nl))
                end if

             else
                ! Separate crop landunit is not created

                pctpft(nl,m) = pctpft(nl,m) * 100._r8/(100._r8-pctspec(nl))
                if (m == nirrig) then
                   if (pctpft(nl,m) > 0._r8) then
                      write(iulog,*) 'ERROR surfrdMod: irrigated crop pft requires create_crop_landunit=.true. to be irrigated'
                      call endrun()
                   end if
                end if
             end if
          end do


       else if (pctspec(nl) == 100._r8) then

          pctpft(nl,0)        = 100._r8
          pctpft(nl,1:numpft) =   0._r8

       else

          write(iulog,*)subname, 'error: pcturb+pctgla+pctlak+pctwet = ',pctspec(nl), &
               ' must be less than or equal to 100'
          call endrun()

       end if
    end do

    ! Find pft and pct arrays 
    ! Save percent cover by PFT [wst] and total percent cover [wst_sum]

    do nl = begg,endg

       wst_sum = 0._r8
       sumpct = 0
       do m = 0, numpft
          wst(m) = pctpft(nl,m)
          wst_sum = wst_sum + pctpft(nl,m)
       end do

       if (domain%pftm(nl) >= 0) then

          ! Rank [wst] in ascendg order to obtain the top [maxpatch_pft] PFTs

          call surfrd_mkrank (numpft, wst, miss, wsti, maxpatch_pft)

          ! Fill in [pft] and [pctpft] with data for top [maxpatch_pft] PFTs.
          ! If land model grid cell is ocean, set to no PFTs.
          ! If land model grid cell is land then:
          !  1. If [pctlnd_o] = 0, there is no PFT data from the input grid.
          !     Since need land data, use bare ground.
          !  2. If [pctlnd_o] > 0, there is PFT data from the input grid but:
          !     a. use the chosen PFT so long as it is not a missing value
          !     b. missing value means no more PFTs with cover > 0
                   
          do m = 1, maxpatch_pft
             if (wsti(m) /=  miss) then
                pft(nl,m) = wsti(m)
                pctpft_lunit(nl,m) = wst(wsti(m))
             else
                pft(nl,m) = noveg
                pctpft_lunit(nl,m) = 0._r8
             end if
             sumpct = sumpct + pctpft_lunit(nl,m)
          end do

       else                               ! model grid wants ocean

          do m = 1, maxpatch_pft
             pft(nl,m) = 0
             pctpft_lunit(nl,m) = 0._r8
          end do

       end if

       ! Correct for the case of more than [maxpatch_pft] PFTs present
                
       if (sumpct < wst_sum) then
          diff  = wst_sum - sumpct
          sumpct = 0._r8
          do m = 1, maxpatch_pft
             pctpft_lunit(nl,m) = pctpft_lunit(nl,m) + diff/maxpatch_pft
             sumpct = sumpct + pctpft_lunit(nl,m)
          end do
       end if

       ! Error check: make sure have a valid PFT

       do m = 1,maxpatch_pft
          if (pft(nl,m) < 0 .or. pft(nl,m) > numpft) then
             write(iulog,*)'surfrd error: invalid PFT at gridcell nl=',nl,pft(nl,m)
             call endrun()
          end if
       end do

       ! As done in mksrfdatMod.F90 for other percentages, truncate pctpft to
       ! ensure that weight relative to landunit is not nonzero
       ! (i.e. a very small number such as 1e-16) where it really should be zero
       ! The following if-block is here to preserve roundoff level differences
       ! between the call to surfrd_wtxy_veg_all and surfrd_wtxy_veg_rank

       if (maxpatch_pft < numpft+1) then
          do m=1,maxpatch_pft
             pctpft_lunit(nl,m) = float(nint(pctpft_lunit(nl,m)))
          end do
          do m=1,numcft
             pctcft_lunit(nl,m) = float(nint(pctcft_lunit(nl,m)))
          end do
       end if
                   
       ! Make sure sum of PFT cover == 100 for land points. If not,
       ! subtract excess from most dominant PFT.

       rmax = -9999._r8
       k1 = -9999
       k2 = -9999
       sumpct = 0._r8
       do m = 1, maxpatch_pft
          sumpct = sumpct + pctpft_lunit(nl,m)
          if (pctpft_lunit(nl,m) > rmax) then
             k1 = m
             rmax = pctpft_lunit(nl,m)
          end if
       end do
       do m = 1, numcft
          sumpct = sumpct + pctcft_lunit(nl,m)
          if (pctcft_lunit(nl,m) > rmax) then
             k2 = m
             rmax = pctcft_lunit(nl,m)
          end if
       end do
       if (k1 == -9999 .and. k2 == -9999) then
          write(iulog,*)'surfrd error: largest PFT patch not found'
          call endrun()
       else if (domain%pftm(nl) >= 0) then
          if (sumpct < 95 .or. sumpct > 105._r8) then
             write(iulog,*)'surfrd error: sum of PFT cover =',sumpct,' at nl=',nl
             call endrun()
          else if (sumpct /= 100._r8 .and. k2 /= -9999) then
             pctcft_lunit(nl,k2) = pctcft_lunit(nl,k2) - (sumpct-100._r8)
          else if (sumpct /= 100._r8) then
             pctpft_lunit(nl,k1) = pctpft_lunit(nl,k1) - (sumpct-100._r8)
          end if
       end if

       ! Error check: make sure PFTs sum to 100% cover

       sumpct = 0._r8
       do m = 1, maxpatch_pft
          sumpct = sumpct + pctpft_lunit(nl,m)
       end do
       do m = 1, numcft
          sumpct = sumpct + pctcft_lunit(nl,m)
       end do
       if (domain%pftm(nl) >= 0) then
          if (abs(sumpct - 100._r8) > 0.000001_r8) then
             write(iulog,*)'surfrdMod error: sum(pct) over maxpatch_pft is not = 100.'
             write(iulog,*)sumpct, nl
             call endrun()
          end if
          if (sumpct < -0.000001_r8) then
             write(iulog,*)'surfrdMod error: sum(pct) over maxpatch_pft is < 0.'
             write(iulog,*)sumpct, nl
             call endrun()
          end if
       end if

    end do   ! end of latitude loop

    ! Determine array [veg], which sets the PFT type for each of the [maxpatch]
    ! patches and array [wtxy], which sets the relative abundance of the PFT.
    ! Fill in PFTs for vegetated portion of grid cell. Fractional areas for
    ! these points [pctpft] pertain to "vegetated" area not to total grid area.
    ! So need to adjust them for fraction of grid that is vegetated.
    ! Next, fill in urban, lake, wetland, and glacier patches.

    do nl = begg,endg
       if (domain%pftm(nl) >= 0) then

          ! Naturally vegetated landunit

          do m = 1, maxpatch_pft
             vegxy(nl,m) = pft(nl,m)
             wtxy(nl,m) = pctpft_lunit(nl,m) * (100._r8-pctspec(nl))/10000._r8
          end do

          ! Crop landunit

          if (create_crop_landunit) then
             do m = 1,numcft
                vegxy(nl,npatch_glacier_mec+m) = cft(nl,m)
                wtxy(nl,npatch_glacier_mec+m) = pctcft_lunit(nl,m) * (100._r8-pctspec(nl))/10000._r8
            end do
          end if

       end if
    end do

    ! Error check

    found = .false.
    sumvec(:) = abs(sum(wtxy,dim=2)-1._r8)
    do nl = begg,endg
       if (sumvec(nl) > 1.e-06_r8 .and. domain%pftm(nl)>=0) then
          found = .true.
          nindx = nl
          exit
       endif
    end do
    if ( found ) then
       write(iulog,*)'surfrd error: WTXY > 1 occurs at nl= ',nindx; call endrun()
    end if

    deallocate(sumvec,cft,pft)
    deallocate(pctcft_lunit,pctpft_lunit,pctpft)
    deallocate(wsti)

  end subroutine surfrd_wtxy_veg_rank

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_wtxy_veg_all
!
! !INTERFACE:
  subroutine surfrd_wtxy_veg_all(ncid, domain)
!
! !DESCRIPTION:
! Determine wtxy and veg arrays for non-dynamic landuse mode
!
! !USES:
    use domainMod   , only : domain_type
    use clm_varctl  , only : create_crop_landunit
    use pftvarcon   , only : nirrig, npcropmin
    use spmdMod     , only : mpicom, MPI_LOGICAL, MPI_LOR
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! netcdf id
    type(domain_type),intent(in) :: domain ! domain associated with wtxy
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein, Sam Levis and Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: m,mp7,mp8,mp11,n,nl,ns         ! indices
    integer  :: begg,endg                      ! beg/end gcell index
    integer  :: dimid,varid                    ! netCDF id's
    integer  :: start(3),count(3)              ! netcdf start/count arrays
    integer  :: ier                            ! error status	
    logical  :: readvar                        ! is variable on dataset
    real(r8) :: sumpct                         ! sum of %pft over maxpatch_pft
    real(r8),allocatable :: pctpft(:,:)        ! percent of vegetated gridcell area for PFTs
    real(r8),pointer :: arrayl(:,:)            ! local array
    real(r8) :: numpftp1data(0:numpft)         
    logical  :: crop = .false.                 ! if crop data on this section of file
    character(len=32) :: subname = 'surfrd_wtxy_veg_all'  ! subroutine name
!-----------------------------------------------------------------------

    ns = domain%ns
    call get_proc_bounds(begg,endg)
    allocate(pctpft(begg:endg,0:numpft))

    call check_dim(ncid, 'lsmpft', numpft+1)

    allocate(arrayl(begg:endg,0:numpft))
    call ncd_io(ncid=ncid, varname='PCT_PFT', flag='read', data=arrayl, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_PFT NOT on surfdata file' )
    pctpft(begg:endg,0:numpft) = arrayl(begg:endg,0:numpft)
    deallocate(arrayl)

    do nl = begg,endg
       if (domain%pftm(nl) >= 0) then

          ! Error check: make sure PFTs sum to 100% cover for vegetated landunit 
          ! (convert pctpft from percent with respect to gridcel to percent with 
          ! respect to vegetated landunit)

          if (pctspec(nl) < 100._r8) then
             sumpct = 0._r8
             do m = 0,numpft
                sumpct = sumpct + pctpft(nl,m) * 100._r8/(100._r8-pctspec(nl))
                if (m == nirrig .and. .not. create_crop_landunit) then
                   if (pctpft(nl,m) > 0._r8) then
                      call endrun( trim(subname)//' ERROR surfrdMod: irrigated crop'// &
                                   ' PFT requires create_crop_landunit=.true.' )
                   end if
                end if
             end do
             if (abs(sumpct - 100._r8) > 0.1e-4_r8) then
                write(iulog,*) trim(subname)//' ERROR: sum(pct) over numpft+1 is not = 100.'
                write(iulog,*) sumpct, sumpct-100._r8, nl
                call endrun()
             end if
             if (sumpct < -0.000001_r8) then
                write(iulog,*) trim(subname)//' ERROR: sum(pct) over numpft+1 is < 0.'
                write(iulog,*) sumpct, nl
                call endrun()
             end if
             do m = npcropmin, numpft
                if ( pctpft(nl,m) > 0.0_r8 ) crop = .true.
             end do
          end if

          ! Set weight of each pft wrt gridcell (note that maxpatch_pft = numpft+1 here)

          do m = 1,numpft+1
             vegxy(nl,m)  = m - 1 ! 0 (bare ground) to numpft
             wtxy(nl,m) = pctpft(nl,m-1) / 100._r8
          end do

       end if
    end do
    call mpi_allreduce(crop,crop_prog,1,MPI_LOGICAL,MPI_LOR,mpicom,ier)
    if (ier /= 0) then
       write(iulog,*) trim(subname)//' mpi_allreduce error = ',ier
       call endrun( trim(subname)//' ERROR: error in reduce of crop_prog' )
    endif
    if (crop_prog .and. .not. create_crop_landunit) then
       call endrun( trim(subname)//' ERROR: prognostic crop '// &
                    'PFTs requires create_crop_landunit=.true.' )
    end if

    deallocate(pctpft)

  end subroutine surfrd_wtxy_veg_all

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_wtxy_veg_dgvm
!
! !INTERFACE:
  subroutine surfrd_wtxy_veg_dgvm(domain)
!
! !DESCRIPTION:
! Determine wtxy and vegxy for CNDV mode.
!
! !USES:
    use domainMod, only : domain_type
    use pftvarcon, only : noveg, crop
    use clm_varctl, only : create_crop_landunit
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(in) :: domain ! domain associated with wtxy
!
! !CALLED FROM:
! subroutine surfrd in this module
!
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/04
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: m,nl         ! indices
    integer  :: begg,endg   ! beg/end gcell index
!-----------------------------------------------------------------------

    call get_proc_bounds(begg,endg)
    do nl = begg,endg
       do m = 1, maxpatch_pft ! CNDV means allocate_all_vegpfts = .true.
          if (create_crop_landunit) then ! been through surfrd_wtxy_veg_all
             if (crop(m-1) == 0) then    ! so update natural vegetation only
                wtxy(nl,m) = 0._r8       ! crops should have values >= 0.
             end if
          else                   ! not been through surfrd_wtxy_veg_all
             wtxy(nl,m) = 0._r8  ! so update all vegetation
             vegxy(nl,m) = m - 1 ! 0 (bare ground) to maxpatch_pft-1 (= 16)
          end if
       end do
       ! bare ground weights
       wtxy(nl,noveg+1) = max(0._r8, 1._r8 - sum(wtxy(nl,:)))
       if (abs(sum(wtxy(nl,:)) - 1._r8) > 1.e-5_r8) then
          write(iulog,*) 'all wtxy =', wtxy(nl,:)
          call endrun()
       end if ! error check
    end do

  end subroutine surfrd_wtxy_veg_dgvm
   
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: surfrd_mkrank
!
! !INTERFACE:
  subroutine surfrd_mkrank (n, a, miss, iv, num)
!
! !DESCRIPTION:
! Return indices of largest [num] values in array [a]
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use abortutils, only : endrun
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: n        ! array length
    real(r8), intent(in) :: a(0:n)   ! array to be ranked
    integer , intent(in) :: miss     ! missing data value
    integer , intent(in) :: num      ! number of largest values requested
    integer , intent(out):: iv(num)  ! index to [num] largest values in array [a]
!
! !CALLED FROM:
! ! subroutine surfrd_wtxy_veg_rank in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
    real(r8) :: a_max       ! maximum value in array
    integer  :: i           ! array index
    real(r8) :: delmax      ! tolerance for finding if larger value
    integer  :: m           ! do loop index
    integer  :: k           ! do loop index
    logical  :: exclude     ! true if data value has already been chosen
!-----------------------------------------------------------------------

    delmax = 1.e-06_r8

    ! Find index of largest non-zero number
    
    iv(1) = miss
    a_max = -9999._r8

    do i = 0, n
       if (a(i)>0._r8 .and. (a(i)-a_max)>delmax) then
          a_max = a(i)
          iv(1)  = i
       end if
    end do

    ! iv(1) = miss indicates no values > 0. this is an error

    if (iv(1) == miss) then
       write(iulog,*) 'surfrd_mkrank error: iv(1) = missing'
       call endrun
    end if

    ! Find indices of the next [num]-1 largest non-zero number.
    ! iv(m) = miss if there are no more values > 0

    do m = 2, num
       iv(m) = miss
       a_max = -9999._r8
       do i = 0, n

          ! exclude if data value has already been chosen

          exclude = .false.
          do k = 1, m-1
             if (i == iv(k)) exclude = .true.
          end do

          ! if not already chosen, see if it is the largest of
          ! the remaining values

          if (.not. exclude) then
             if (a(i)>0._r8 .and. (a(i)-a_max)>delmax) then
                a_max = a(i)
                iv(m)  = i
             end if
          end if
       end do
    end do

  end subroutine surfrd_mkrank

end module surfrdMod
