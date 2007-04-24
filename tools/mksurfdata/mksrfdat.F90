program mksrfdat

!-----------------------------------------------------------------------
!BOP
!
! !PROGRAM: mksrfdat
!
! !DESCRIPTION:
! Creates land model surface dataset from original "raw" data files.
! Surface dataset contains model grid, pfts, inland water, glacier,
! soil texture, soil color, LAI and SAI and urban fraction.
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use shr_sys_mod , only : shr_sys_getenv
    use shr_timer_mod
    use fileutils   , only : getfil, putfil, opnfil, getavu, get_filename
    use mklaiMod    , only : mklai
    use mkpftMod    , only : mkpft
    use creategridMod, only : read_domain,write_domain
    use domainMod   , only : domain_setptrs, domain_type, domain_init
    use mkfileMod   , only : mkfile
    use mkvarpar    , only : numpft, nlevsoi
    use mkvarsur
    use mkvarctl
    use areaMod
    use ncdio
    use nanMod
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Authors: Gordon Bonan, Sam Levis and Mariana Vertenstein
! Revised: Nan Rosenbloom to add fmax processing.
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: lsmlon, lsmlat       ! clm grid resolution
    integer  :: nsoicol              ! number of model color classes
    integer  :: i,j,k,m,n            ! indices
    integer  :: ni, nj               ! size of pft index
    integer  :: ier                  ! error status
    integer  :: ndiag,nfdyn          ! unit numbers
    integer  :: ncid                 ! netCDF id
    integer  :: omode                ! netCDF output mode
    integer  :: varid                ! netCDF variable id
    integer  :: beg4d(4),len4d(4)    ! netCDF variable edges
    integer  :: ret                  ! netCDF return status
    integer  :: ntim                 ! time sample for dynamic land use
    integer  :: year                 ! year for dynamic land use
    real(r8) :: sum                  ! sum for error check
    real(r8) :: rmax                 ! maximum patch cover
    character(len=256) :: fgrddat    ! grid data file
    character(len=256) :: fsurdat    ! surface data file name
    character(len=256) :: fdyndat    ! dynamic landuse data file name
    character(len=256) :: fname      ! generic filename
    character(len=256) :: loc_fn     ! local file name
    character(len=  9) :: resol      ! resolution for file name
    integer  :: t1                   ! timer
    real(r8),parameter :: p5  = 0.5_r8   ! constant
    real(r8),parameter :: p25 = 0.25_r8  ! constant

    real(r8), allocatable  :: landfrac_pft(:,:)    ! PFT data: % land per gridcell
    real(r8), allocatable  :: pctlnd_pft(:,:)      ! PFT data: % of gridcell for PFTs
    real(r8), allocatable  :: pctlnd_pft_dyn(:,:)  ! PFT data: % of gridcell for dyn landuse PFTs
    integer , allocatable  :: pftdata_mask(:,:)    ! mask indicating real or fake land type
    real(r8), allocatable  :: pctpft(:,:,:)        ! PFT data: land fraction per gridcell
    real(r8), pointer      :: pctpft_i(:,:,:)      ! PFT data: % fraction on input grid
    real(r8), allocatable  :: pctgla(:,:)          ! percent of grid cell that is glacier  
    real(r8), allocatable  :: pctlak(:,:)          ! percent of grid cell that is lake     
    real(r8), allocatable  :: pctwet(:,:)          ! percent of grid cell that is wetland  
    real(r8), allocatable  :: pcturb(:,:)          ! percent of grid cell that is urbanized
    real(r8), allocatable  :: fmax(:,:)            ! fractional saturated area
    integer , allocatable  :: soic2d(:,:)          ! soil color                            
    real(r8), allocatable  :: sand3d(:,:,:)        ! soil texture: percent sand            
    real(r8), allocatable  :: clay3d(:,:,:)        ! soil texture: percent clay            

    logical  :: do_double_res_too     ! write matching double resolution dataset
    character(len=256) :: fdfile      ! double surface data file name
    integer  :: dbllon, dbllat        ! double clm grid resolution
    type(domain_type) :: dbldomain    ! double domain
    integer  :: di,dj                 ! double i,j
    integer  :: ncidi                 ! netCDF id
    integer  :: dnsoicol              ! number of model color classes
    real(r8), allocatable  :: dlandfrac_pft(:,:)    ! PFT data: % land per gridcell
    real(r8), allocatable  :: dpctlnd_pft(:,:)      ! PFT data: % of gridcell for PFTs
    real(r8), allocatable  :: dpctlnd_pft_dyn(:,:)  ! PFT data: % of gridcell for dyn landuse PFTs
    integer , allocatable  :: dpftdata_mask(:,:)    ! mask indicating real or fake land type
    real(r8), allocatable  :: dpctpft(:,:,:)        ! PFT data: land fraction per gridcell
    real(r8), allocatable  :: dpctgla(:,:)          ! percent of grid cell that is glacier  
    real(r8), allocatable  :: dpctlak(:,:)          ! percent of grid cell that is lake     
    real(r8), allocatable  :: dpctwet(:,:)          ! percent of grid cell that is wetland  
    real(r8), allocatable  :: dpcturb(:,:)          ! percent of grid cell that is urbanized
    real(r8), allocatable  :: dfmax(:,:)            ! fractional saturated area
    integer , allocatable  :: dsoic2d(:,:)          ! soil color                            
    real(r8), allocatable  :: dsand3d(:,:,:)        ! soil texture: percent sand            
    real(r8), allocatable  :: dclay3d(:,:,:)        ! soil texture: percent clay            
    real(r8), allocatable  :: mlai_i(:,:,:)         ! baseline monthly lai
    real(r8), allocatable  :: msai_i(:,:,:)         ! baseline monthly sai
    real(r8), allocatable  :: mhgtt_i(:,:,:)        ! baseline monthly height (top)
    real(r8), allocatable  :: mhgtb_i(:,:,:)        ! baseline monthly height (bottom)
    real(r8), allocatable  :: mlai_o(:,:,:)         ! double monthly lai
    real(r8), allocatable  :: msai_o(:,:,:)         ! double monthly sai
    real(r8), allocatable  :: mhgtt_o(:,:,:)        ! double monthly height (top)
    real(r8), allocatable  :: mhgtb_o(:,:,:)        ! double monthly height (bottom)

    character(len=32) :: subname = 'mksrfdat'  ! program name

    namelist /clmexp/              &
	 mksrf_fgrid,              &	
	 mksrf_gridtype,           &	
         mksrf_fvegtyp,            &
	 mksrf_fsoitex,            &
         mksrf_fsoicol,            &
         mksrf_flanwat,            &
         mksrf_fglacier,           &
         mksrf_fmax,               &
         mksrf_furban,             &
         mksrf_flai,               &
         mksrf_fdynuse,            &
         do_double_res_too      
!-----------------------------------------------------------------------

    ! ======================================================================
    ! Read input namelist
    ! ======================================
    ! Must specify settings for either:
    ! ======================================
    !	 mksrf_fgrid
    ! ======================================
    ! Must specify settings for:
    ! ======================================
    !    mksrf_fvegtyp
    !	 mksrf_fsoitex
    !    mksrf_fsoicol
    !    mksrf_flanwat
    !    mksrf_fmax
    !    mksrf_fglacier
    !    mksrf_furban
    !    mksrf_flai
    ! ======================================
    ! Optionally specify setting for:
    ! ======================================
    !    mksrf_gridtype
    !    mksrf_fdynuse
    !    do_double_res_too
    ! ======================================================================

    call shr_timer_init()
    call shr_timer_get(t1,'accumulating timer')
    call shr_timer_start(t1)

    write(6,*) 'Attempting to initialize control settings .....'

    do_double_res_too = .false.
    mksrf_gridtype = 'global'
    read(5, clmexp, iostat=ier)
    if (ier /= 0) then
       write(6,*)'error: namelist input resulted in error code ',ier
       call abort()
    endif

    write (6,*) 'Attempting to create surface boundary data .....'
    write (6,'(72a1)') ("-",i=1,60)

    ! ----------------------------------------------------------------------
    ! Error check namelist input
    ! ----------------------------------------------------------------------
    
    if (mksrf_fgrid /= ' ')then
       fgrddat = mksrf_fgrid
       write(6,*)'mksrf_fgrid = ',mksrf_fgrid
    else
       write (6,*)'must specify mksrf_fgrid'
       stop
    endif

    if (trim(mksrf_gridtype) == 'global' .or. &
        trim(mksrf_gridtype) == 'regional') then
       write(6,*)'mksrf_gridtype = ',trim(mksrf_gridtype)
    else
       write(6,*)'mksrf_gridtype = ',trim(mksrf_gridtype)
       write (6,*)'illegal mksrf_gridtype, must be global or regional '
       stop
    endif

    ! ----------------------------------------------------------------------
    ! Interpolate input dataset to model resolution
    ! ----------------------------------------------------------------------
    
    ! Determine land model grid, fractional land and land mask
    
    call read_domain(ldomain,fgrddat)
    ! --- invalidate mask and frac for ldomain ---
    ldomain%mask = bigint
    ldomain%frac = nan

    call domain_setptrs(ldomain,lsmlon,lsmlat)
    write(6,*)'lsmlon= ',lsmlon,' lsmlat= ',lsmlat

    ! Allocate and initialize dynamic memory

    allocate ( landfrac_pft(lsmlon,lsmlat)    , &
               pctlnd_pft(lsmlon,lsmlat)      , & 
               pctlnd_pft_dyn(lsmlon,lsmlat)  , & 
               pftdata_mask(lsmlon,lsmlat)    , & 
               pctpft(lsmlon,lsmlat,0:numpft) , & 
               pctgla(lsmlon,lsmlat)          , & 
               pctlak(lsmlon,lsmlat)          , & 
               pctwet(lsmlon,lsmlat)          , & 
               pcturb(lsmlon,lsmlat)          , & 
               fmax(lsmlon,lsmlat)            , & 
               sand3d(lsmlon,lsmlat,nlevsoi)  , & 
               clay3d(lsmlon,lsmlat,nlevsoi)  , & 
               soic2d(lsmlon,lsmlat))

    landfrac_pft(:,:) = spval 
    pctlnd_pft(:,:)   = spval
    pftdata_mask(:,:) = -999
    pctpft(:,:,:)     = spval
    pctgla(:,:)       = spval
    pctlak(:,:)       = spval
    pctwet(:,:)       = spval
    pcturb(:,:)       = spval
    fmax(:,:)         = spval
    sand3d(:,:,:)     = spval
    clay3d(:,:,:)     = spval
    soic2d(:,:)       = -999

    write(6,*) ' timer_a2 init-----'
    call shr_timer_print(t1)

    ! ----------------------------------------------------------------------
    ! Open diagnostic output log file
    ! ----------------------------------------------------------------------
    
    write (resol,'(i4.4,"x",i4.4)') lsmlat,lsmlon
    loc_fn= './surfdata_'//trim(resol)//'.log'
    ndiag = getavu()
    call opnfil (loc_fn, ndiag, 'f')
    
    if (mksrf_fgrid /= ' ')then
       write (ndiag,*)'using fractional land data from file= ', &
            trim(mksrf_fgrid),' to create the surface dataset'
    endif

    if (trim(mksrf_gridtype) == 'global' .or. &
        trim(mksrf_gridtype) == 'regional') then
       write(6,*)'mksrf_gridtype = ',trim(mksrf_gridtype)
    endif

    write (ndiag,*) 'PFTs from:         ',trim(mksrf_fvegtyp)
    write (ndiag,*) 'fmax from:         ',trim(mksrf_fmax)
    write (ndiag,*) 'glaciers from:     ',trim(mksrf_fglacier)
    write (ndiag,*) 'urban from:        ',trim(mksrf_furban)
    write (ndiag,*) 'inland water from: ',trim(mksrf_flanwat)
    write (ndiag,*) 'soil texture from: ',trim(mksrf_fsoitex)
    write (ndiag,*) 'soil color from:   ',trim(mksrf_fsoicol)

    write(6,*) ' timer_a1 init-----'
    call shr_timer_print(t1)


    ! Make PFTs [pctpft] from dataset [fvegtyp] (1/2 degree PFT data)

    call mkpft(lsmlon, lsmlat, mksrf_fvegtyp, ndiag, pctlnd_pft, pctpft, pctpft_i)

    write(6,*) ' timer_b mkpft-----'
    call shr_timer_print(t1)

    ! Make inland water [pctlak, pctwet] from Cogley's one degree data [flanwat]

    call mklanwat (lsmlon, lsmlat, mksrf_flanwat, ndiag, pctlak, pctwet)

    write(6,*) ' timer_c mklanwat-----'
    call shr_timer_print(t1)

    ! Make fmax [fmax] from [fmax] dataset

    call mkfmax (lsmlon, lsmlat, mksrf_fmax, ndiag, fmax)

    write(6,*) ' timer_d mkfmax-----'
    call shr_timer_print(t1)

    ! Make glacier fraction [pctgla] from [fglacier] dataset

    call mkglacier (lsmlon, lsmlat, mksrf_fglacier, ndiag, pctgla)

    write(6,*) ' timer_d mkglacier-----'
    call shr_timer_print(t1)

    ! Make soil texture [sand3d, clay3d] from IGBP 5 minute data [fsoitex]

    call mksoitex (lsmlon, lsmlat, mksrf_fsoitex, ndiag, pctgla, sand3d, clay3d)

    write(6,*) ' timer_e mksoitex-----'
    call shr_timer_print(t1)

    ! Make soil color classes [soic2d] from BATS T42 data [fsoicol]

    call mksoicol (lsmlon, lsmlat, mksrf_fsoicol, ndiag, pctgla, soic2d, nsoicol)

    write(6,*) ' timer_f mksoicol-----'
    call shr_timer_print(t1)

    ! Make urban fraction [pcturb] from [furban] dataset

    call mkurban (lsmlon, lsmlat, mksrf_furban, ndiag, pcturb)

    write(6,*) ' timer_g mkurban-----'
    call shr_timer_print(t1)

    ! Set land values on Ross ice shelf to glacier

    do j = 1,ldomain%nj
       do i = 1,ldomain%numlon(j)
          if (ldomain%latixy(i,j) < -79.) then
             soic2d(i,j) = 0
             pctlak(i,j) = 0.
             pctwet(i,j) = 0.
             pcturb(i,j) = 0.
             pctgla(i,j) = 100.
             pctpft(i,j,:)  = 0.
             sand3d(i,j,1:nlevsoi) = 0.
             clay3d(i,j,1:nlevsoi) = 0.
          end if
       end do
    end do

    ! Assume wetland and/or lake when dataset landmask implies ocean 
    ! (assume medium soil color (15) and loamy texture).
    ! Also set pftdata_mask here

    do j = 1,ldomain%nj
       do i = 1,ldomain%numlon(j)
          if (pctlnd_pft(i,j) < 1.e-6) then
             pftdata_mask(i,j) = 0
             soic2d(i,j) = 15
             pctwet(i,j) = 100. - pctlak(i,j)
             pcturb(i,j) = 0.
             pctgla(i,j) = 0.
             pctpft(i,j,:) = 0.
             sand3d(i,j,1:nlevsoi) = 43.
             clay3d(i,j,1:nlevsoi) = 18.
          else
             pftdata_mask(i,j) = 1
          end if
       end do
    end do

    ! If have pole points on grid - set south pole to glacier
    ! north pole is as assumed as non-land

    if (abs((ldomain%latixy(1,lsmlat) - 90.)) < 1.e-6) then
       write(6,*)'MKSRFDAT: grid has pole_points'
       do i = 1,ldomain%numlon(1)
          soic2d(i,1)   = 0
          pctlak(i,1)   = 0.
          pctwet(i,1)   = 0.
          pcturb(i,1)   = 0.
          sand3d(i,1,:) = 0.
          clay3d(i,1,:) = 0.
          pctgla(i,1)   = 100.
          pctpft(i,1,:) = 0.
       end do
    end if

    ! Truncate all percentage fields on output grid. This is needed to
    ! insure that wt is not nonzero (i.e. a very small number such as
    ! 1e-16) where it really should be zero

    do j = 1,ldomain%nj
       do i = 1,ldomain%numlon(j)
          do k = 1,nlevsoi
             sand3d(i,j,k) = float(nint(sand3d(i,j,k)))
             clay3d(i,j,k) = float(nint(clay3d(i,j,k)))
          end do
          pctlak(i,j) = float(nint(pctlak(i,j)))
          pctwet(i,j) = float(nint(pctwet(i,j)))
          pcturb(i,j) = float(nint(pcturb(i,j)))
          pctgla(i,j) = float(nint(pctgla(i,j)))
       end do
    end do

    ! Make sure sum of land cover types does not exceed 100. If it does,
    ! subtract excess from most dominant land cover.

    do j = 1,ldomain%nj
       do i = 1,ldomain%numlon(j)

          sum = pctlak(i,j)+pctwet(i,j)+pcturb(i,j)+pctgla(i,j)
          if (sum > 120.) then
             write (6,*) 'MKSRFDAT error: sum of pctlak, pctwet,', &
                  'pcturb and pctgla is greater than 120%'
             write (6,*)'i,j,pctlak,pctwet,pcturb,pctgla= ', &
                  i,j,pctlak(i,j),pctwet(i,j),pcturb(i,j),pctgla(i,j)
             call abort()
          else if (sum > 100._r8) then
             pctlak(i,j) = pctlak(i,j) * 100._r8/sum
             pctwet(i,j) = pctwet(i,j) * 100._r8/sum
             pcturb(i,j) = pcturb(i,j) * 100._r8/sum
             pctgla(i,j) = pctgla(i,j) * 100._r8/sum
          end if

          ! Normalize pctpft to be the remainder of [100 - (special landunits)]

          sum = pctlak(i,j)+pctwet(i,j)+pcturb(i,j)+pctgla(i,j)
          do m = 0, numpft
             pctpft(i,j,m) = 0.01_r8 * pctpft(i,j,m) * (100._r8 - sum)
          end do

          sum = pctlak(i,j)+pctwet(i,j)+pcturb(i,j)+pctgla(i,j)
          do m = 0,numpft
             sum = sum + pctpft(i,j,m)
          end do

          if (sum < 90._r8) then
             write (6,*) 'MKSRFDAT error: sum of pctlak, pctwet,', &
                  'pcturb, pctgla and pctpft is less than 90'
             write (6,*)'i,j,pctlak,pctwet,pcturb,pctgla,pctpft= ', &
                  i,j,pctlak(i,j),pctwet(i,j),pcturb(i,j),pctgla(i,j),&
                  pctpft(i,j,:)
             call abort()
          else if (sum > 100._r8 + 1.e-6) then
             write (6,*) 'MKSRFDAT error: sum of pctlak, pctwet,', &
                  'pcturb, pctgla and pctpft is greater than 100'
             write (6,*)'i,j,pctlak,pctwet,pcturb,pctgla,pctpft= ', &
                  i,j,pctlak(i,j),pctwet(i,j),pcturb(i,j),pctgla(i,j),&
                  pctpft(i,j,:)
             call abort()
          else
             pctlak(i,j)   = pctlak(i,j) * 100./sum
             pctwet(i,j)   = pctwet(i,j) * 100./sum
             pcturb(i,j)   = pcturb(i,j) * 100./sum
             pctgla(i,j)   = pctgla(i,j) * 100./sum
             pctpft(i,j,:) = pctpft(i,j,:) * 100./sum
          end if

       end do
    end do

    do k = 0,numpft
       sum = 0.
       do j = 1,ldomain%nj
          do i = 1,ldomain%numlon(j)
             sum = sum + pctpft(i,j,k)
          enddo
       enddo
       write(6,*) 'sum over domain of pft ',k,sum
    enddo

    ! Determine fractional land from pft dataset

    do j = 1,lsmlat
       do i = 1,lsmlon
          landfrac_pft(i,j) = pctlnd_pft(i,j)/100.
       end do
    end do

    write(6,*) ' timer_h final-----'
    call shr_timer_print(t1)

    ! ----------------------------------------------------------------------
    ! Create surface dataset
    ! ----------------------------------------------------------------------

    ! Create netCDF surface dataset.  

    fsurdat = './surfdata_'//trim(resol)//'.nc'

    call mkfile(ldomain%ni, ldomain%nj, fsurdat, dynlanduse = .false.)
    call write_domain(ldomain,fsurdat)

    call check_ret(nf_open(trim(fsurdat), nf_write, ncid), subname)
    call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

    ! Write fields other than lai, sai, and heights to netcdf surface dataset

    call ncd_ioglobal(varname='PFTDATA_MASK', data=pftdata_mask, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='LANDFRAC_PFT', data=landfrac_pft, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='mxsoil_color', data=nsoicol     , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='SOIL_COLOR'  , data=soic2d      , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='PCT_SAND'    , data=sand3d      , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='PCT_CLAY'    , data=clay3d      , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='PCT_WETLAND' , data=pctwet      , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='PCT_LAKE'    , data=pctlak      , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='PCT_GLACIER' , data=pctgla      , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='PCT_URBAN'   , data=pcturb      , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='PCT_PFT'     , data=pctpft      , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='FMAX'        , data=fmax        , ncid=ncid, flag='write')

    ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

    call check_ret(nf_sync(ncid), subname)

    write(6,*) ' timer_i writesurf-----'
    call shr_timer_print(t1)

    ! ----------------------------------------------------------------------
    ! Make LAI and SAI from 1/2 degree data and write to surface dataset 
    ! Write to netcdf file is done inside mklai routine
    ! ----------------------------------------------------------------------

    if ( .not. associated(pctpft_i) )then
       write(6,*)'error: pctpft_i is not allocated at this point'
       call abort()
    end if
    ni = size(pctpft_i,1)
    nj = size(pctpft_i,2)
    call mklai(lsmlon, lsmlat, mksrf_flai, ndiag, ncid, ni, nj, pctpft_i)
    deallocate( pctpft_i )

    write(6,*) ' timer_j mklai-----'
    call shr_timer_print(t1)

    ! Close surface dataset

    call check_ret(nf_close(ncid), subname)

    write (6,'(72a1)') ("-",i=1,60)
    write (6,'(a52,f5.1,a4,f5.1,a5)') 'land model surface data set successfully created for ', &
         360./lsmlon,' by ',180./lsmlat,' grid'

    ! ----------------------------------------------------------------------
    ! Create dynamic land use dataset if appropriate
    ! ----------------------------------------------------------------------

    if (mksrf_fdynuse /= ' ') then

       write (resol,'(i4.4,"x",i4.4)') lsmlon,lsmlat
       fdyndat = './surface-data.dynpft.'//trim(resol)//'.nc'

       ! Define dimensions and global attributes

       call mkfile(lsmlon, lsmlat, fdyndat, dynlanduse = .true.)
       call write_domain(ldomain,fdyndat)

       ! Write fields other pft to dynamic land use dataset

       call check_ret(nf_open(trim(fdyndat), nf_write, ncid), subname)
       call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

       call ncd_ioglobal(varname='PFTDATA_MASK', data=pftdata_mask, ncid=ncid, flag='write')
       call ncd_ioglobal(varname='LANDFRAC_PFT', data=landfrac_pft, ncid=ncid, flag='write')
       call ncd_ioglobal(varname='PCT_WETLAND' , data=pctwet      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='PCT_LAKE'    , data=pctlak      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='PCT_GLACIER' , data=pctgla      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='PCT_URBAN'   , data=pcturb      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='FMAX'        , data=fmax        , ncid=ncid, flag='write')

       ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

       call check_ret(nf_sync(ncid), subname)

       ! Read in each dynamic pft landuse dataset

       call getfil (mksrf_fdynuse, loc_fn, 0)
       nfdyn = getavu(); call opnfil (loc_fn, nfdyn, 'f')

       ntim = 0
       do 
          ! Read input pft data filename

          read(nfdyn, *, iostat=ier) fname, year
	  write(6,*)'input pft dynamic dataset is ',fname,' year is ',year
          if (ier /= 0) exit
          ntim = ntim + 1

          ! Create pctpft data at model resolution

          call mkpft(lsmlon, lsmlat, fname, ndiag, pctlnd_pft_dyn, pctpft)

          ! If have pole points, set south pole to glacier (north pole is as assumed as non-land)

          if (abs((ldomain%latixy(1,lsmlat) - 90.)) < 1.e-6) then
             write(6,*)'MKSRFDAT: grid has pole_points'
             do i = 1,ldomain%numlon(1)
                pctpft(i,1,0)        = 100.
                pctpft(i,1,1:numpft) = 0.
             end do
          end if

          ! Adjust pctpft before output

          do j = 1,lsmlat
             do i = 1,lsmlon

                ! Consistency check on input land fraction

                if (pctlnd_pft_dyn(i,j) /= pctlnd_pft(i,j)) then
                   write(6,*) subname,' error: pctlnd_pft for dynamics data = ',&
                        pctlnd_pft_dyn(i,j), ' not equal to pctlnd_pft for surface data = ',&
                        pctlnd_pft(i,j),' at i,j= ',i,j,' and filename = ',trim(fname)
                   call abort()
                end if

                ! Set land values on Ross ice shelf to glacier

                if (ldomain%latixy(i,j) < -79.) then
                   pctpft(i,j,0) = 100.
                   pctpft(i,j,1:numpft)  = 0.
                end if

                ! Assume nonvegetated wetland when input landmask says land and 
                ! pft landmask says ocean 

                if (pctlnd_pft(i,j)==0.) then
                   pctpft(i,j,0) = 100.
                   pctpft(i,j,1:numpft) = 0.
                   pftdata_mask(i,j) = 0
                else
                   pftdata_mask(i,j) = 1
                end if

                ! Normalize pctpft to be the remainder of [100 - (special landunits)]

                sum = pctlak(i,j)+pctwet(i,j)+pcturb(i,j)+pctgla(i,j)
                do m = 0, numpft
                   pctpft(i,j,m) = 0.01_r8 * pctpft(i,j,m) * (100._r8 - sum)
                end do

             end do
          end do

          ! Output pctpft data for current year

          beg4d(1) = 1     ;  len4d(1) = lsmlon
          beg4d(2) = 1     ;  len4d(2) = lsmlat
          beg4d(3) = 1     ;  len4d(3) = numpft+1
          beg4d(4) = ntim  ;  len4d(4) = 1

          call check_ret(nf_inq_varid(ncid, 'PCT_PFT', varid), subname)
          call check_ret(nf_put_vara_double(ncid, varid, beg4d, len4d, pctpft), subname)

          call check_ret(nf_inq_varid(ncid, 'YEAR', varid), subname)
          call check_ret(nf_put_vara_int(ncid, varid, beg4d(4), len4d(4), year), subname)

	  ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

	  call check_ret(nf_sync(ncid), subname)

       end do   ! end of read loop

       call check_ret(nf_close(ncid), subname)

    end if   ! end of if-create dynamic landust dataset   

    write(6,*) ' timer_k writedyn-----'
    call shr_timer_print(t1)

    ! ----------------------------------------------------------------------
    ! Close diagnostic dataset
    ! ----------------------------------------------------------------------

    close (ndiag)
    write (6,*)
    write (6,*) 'Surface data output file = ',trim(fsurdat)
    write (6,*) '   This file contains the land model surface data'
    write (6,*) 'Diagnostic log file      = ',trim(loc_fn)
    write (6,*) '   See this file for a summary of the dataset'
    write (6,*)

    write(6,*) ' timer_v double_res-----'
    call shr_timer_print(t1)

    if (do_double_res_too) then
       dbllon = 2*lsmlon
       dbllat = 2*lsmlat
       call domain_init(dbldomain,dbllon,dbllat)

       allocate ( dlandfrac_pft(dbllon,dbllat)    , &
                  dpctlnd_pft(dbllon,dbllat)      , & 
                  dpctlnd_pft_dyn(dbllon,dbllat)  , & 
                  dpftdata_mask(dbllon,dbllat)    , & 
                  dpctpft(dbllon,dbllat,0:numpft) , & 
                  dpctgla(dbllon,dbllat)          , & 
                  dpctlak(dbllon,dbllat)          , & 
                  dpctwet(dbllon,dbllat)          , & 
                  dpcturb(dbllon,dbllat)          , & 
                  dfmax(dbllon,dbllat)            , & 
                  dsand3d(dbllon,dbllat,nlevsoi)  , & 
                  dclay3d(dbllon,dbllat,nlevsoi)  , & 
                  dsoic2d(dbllon,dbllat))
       allocate ( mlai_i(lsmlon,lsmlat,0:numpft)  , &
                  msai_i(lsmlon,lsmlat,0:numpft)  , &
                  mhgtt_i(lsmlon,lsmlat,0:numpft) , &
                  mhgtb_i(lsmlon,lsmlat,0:numpft))
       allocate ( mlai_o(dbllon,dbllat,0:numpft)  , &
                  msai_o(dbllon,dbllat,0:numpft)  , &
                  mhgtt_o(dbllon,dbllat,0:numpft) , &
                  mhgtb_o(dbllon,dbllat,0:numpft))

       dbldomain%ni = dbllon
       dbldomain%nj = dbllat
       dbldomain%edgen = ldomain%edgen
       dbldomain%edgee = ldomain%edgee
       dbldomain%edges = ldomain%edges
       dbldomain%edgew = ldomain%edgew
       dbldomain%numlon = dbllon

       dnsoicol = nsoicol

       do j = 1,lsmlat
       do i = 1,lsmlon
         di = 2*(i-1)+1
         dj = 2*(j-1)+1
         dbldomain%mask  (di:di+1,dj:dj+1) = ldomain%mask(i,j)
         dbldomain%frac  (di:di+1,dj:dj+1) = ldomain%frac(i,j)
         dbldomain%latixy(di:di+1,dj     ) = ldomain%latixy(i,j)*p5 + ldomain%lats(i,j)*p5
         dbldomain%latixy(di:di+1,   dj+1) = ldomain%latixy(i,j)*p5 + ldomain%latn(i,j)*p5
         dbldomain%longxy(di     ,dj:dj+1) = ldomain%longxy(i,j)*p5 + ldomain%lonw(i,j)*p5
         dbldomain%longxy(   di+1,dj:dj+1) = ldomain%longxy(i,j)*p5 + ldomain%lone(i,j)*p5
         dbldomain%area  (di:di+1,dj:dj+1) = ldomain%area(i,j)*p25
         dbldomain%lats  (di:di+1,dj     ) = ldomain%lats(i,j)
         dbldomain%lats  (di:di+1,   dj+1) = ldomain%latixy(i,j)
         dbldomain%latn  (di:di+1,dj     ) = ldomain%latixy(i,j)
         dbldomain%latn  (di:di+1,   dj+1) = ldomain%latn(i,j)
         dbldomain%lonw  (di     ,dj:dj+1) = ldomain%lonw(i,j)
         dbldomain%lonw  (   di+1,dj:dj+1) = ldomain%longxy(i,j)
         dbldomain%lone  (di     ,dj:dj+1) = ldomain%longxy(i,j)
         dbldomain%lone  (   di+1,dj:dj+1) = ldomain%lone(i,j)
         dlandfrac_pft   (di:di+1,dj:dj+1) = landfrac_pft(i,j)
         dpctlnd_pft     (di:di+1,dj:dj+1) = pctlnd_pft(i,j)
         dpctlnd_pft_dyn (di:di+1,dj:dj+1) = pctlnd_pft_dyn(i,j)
         dpftdata_mask   (di:di+1,dj:dj+1) = pftdata_mask(i,j)
         do n = 0,numpft
            dpctpft      (di:di+1,dj:dj+1,n) = pctpft(i,j,n)
         enddo
         dpctgla         (di:di+1,dj:dj+1) = pctgla(i,j)
         dpctlak         (di:di+1,dj:dj+1) = pctlak(i,j)
         dpctwet         (di:di+1,dj:dj+1) = pctwet(i,j)
         dpcturb         (di:di+1,dj:dj+1) = pcturb(i,j)
         dfmax           (di:di+1,dj:dj+1) = fmax(i,j)
         do n = 1,nlevsoi
            dsand3d      (di:di+1,dj:dj+1,n) = sand3d(i,j,n)
            dclay3d      (di:di+1,dj:dj+1,n) = clay3d(i,j,n)
         enddo
         dsoic2d         (di:di+1,dj:dj+1) = soic2d(i,j)
       enddo
       enddo

       write (resol,'(i4.4,"x",i4.4)') dbllat,dbllon
       fdfile = './surfdata_'//trim(resol)//'.double.nc'

       call mkfile(dbldomain%ni, dbldomain%nj, fdfile, dynlanduse = .false.)
       call write_domain(dbldomain,fdfile)

       call check_ret(nf_open(trim(fdfile), nf_write, ncid), subname)
       call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

       ! Write fields other than lai, sai, and heights to netcdf surface dataset

       call ncd_ioglobal(varname='PFTDATA_MASK', data=dpftdata_mask, ncid=ncid, flag='write')
       call ncd_ioglobal(varname='LANDFRAC_PFT', data=dlandfrac_pft, ncid=ncid, flag='write')
       call ncd_ioglobal(varname='mxsoil_color', data=dnsoicol     , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='SOIL_COLOR'  , data=dsoic2d      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='PCT_SAND'    , data=dsand3d      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='PCT_CLAY'    , data=dclay3d      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='PCT_WETLAND' , data=dpctwet      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='PCT_LAKE'    , data=dpctlak      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='PCT_GLACIER' , data=dpctgla      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='PCT_URBAN'   , data=dpcturb      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='PCT_PFT'     , data=dpctpft      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='FMAX'        , data=dfmax        , ncid=ncid, flag='write')

       call check_ret(nf_open(fsurdat, 0, ncidi), subname)
       do m = 1,12
          beg4d(1) = 1  ;  len4d(1) = lsmlon
          beg4d(2) = 1  ;  len4d(2) = lsmlat
          beg4d(3) = 1  ;  len4d(3) = numpft+1
          beg4d(4) = m  ;  len4d(4) = 1

          call check_ret(nf_inq_varid(ncidi, 'MONTHLY_LAI', varid), subname)
          call check_ret(nf_get_vara_double(ncidi, varid, beg4d, len4d, mlai_i), subname)
          call check_ret(nf_inq_varid(ncidi, 'MONTHLY_SAI', varid), subname)
          call check_ret(nf_get_vara_double(ncidi, varid, beg4d, len4d, msai_i), subname)
          call check_ret(nf_inq_varid(ncidi, 'MONTHLY_HEIGHT_TOP', varid), subname)
          call check_ret(nf_get_vara_double(ncidi, varid, beg4d, len4d, mhgtt_i), subname)
          call check_ret(nf_inq_varid(ncidi, 'MONTHLY_HEIGHT_BOT', varid), subname)
          call check_ret(nf_get_vara_double(ncidi, varid, beg4d, len4d, mhgtb_i), subname)

          do j = 1,lsmlat
          do i = 1,lsmlon
             di = 2*(i-1)+1
             dj = 2*(j-1)+1
             do n = 0,numpft
                mlai_o (di:di+1,dj:dj+1,n) = mlai_i (i,j,n)
                msai_o (di:di+1,dj:dj+1,n) = msai_i (i,j,n)
                mhgtt_o(di:di+1,dj:dj+1,n) = mhgtt_i(i,j,n)
                mhgtb_o(di:di+1,dj:dj+1,n) = mhgtb_i(i,j,n)
             enddo
          enddo
          enddo

          beg4d(1) = 1  ;  len4d(1) = dbllon
          beg4d(2) = 1  ;  len4d(2) = dbllat
          beg4d(3) = 1  ;  len4d(3) = numpft+1
          beg4d(4) = m  ;  len4d(4) = 1

          call check_ret(nf_inq_varid(ncid, 'MONTHLY_LAI', varid), subname)
          call check_ret(nf_put_vara_double(ncid, varid, beg4d, len4d, mlai_o), subname)
          call check_ret(nf_inq_varid(ncid, 'MONTHLY_SAI', varid), subname)
          call check_ret(nf_put_vara_double(ncid, varid, beg4d, len4d, msai_o), subname)
          call check_ret(nf_inq_varid(ncid, 'MONTHLY_HEIGHT_TOP', varid), subname)
          call check_ret(nf_put_vara_double(ncid, varid, beg4d, len4d, mhgtt_o), subname)
          call check_ret(nf_inq_varid(ncid, 'MONTHLY_HEIGHT_BOT', varid), subname)
          call check_ret(nf_put_vara_double(ncid, varid, beg4d, len4d, mhgtb_o), subname)
       enddo

       call check_ret(nf_sync(ncid), subname)
       ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

       call check_ret(nf_sync(ncid), subname)

       call check_ret(nf_close(ncid), subname)
 
       write (6,'(72a1)') ("-",i=1,60)
       write (6, *) 'Generated file ',trim(fdfile)
       write (6,'(a52,f5.1,a4,f5.1,a5)') 'land model surface data set successfully created for ', &
          360./dbllon,' by ',180./dbllat,' double grid'
    endif

    write(6,*) ' timer_z end-----'
    call shr_timer_print(t1)

end program mksrfdat
