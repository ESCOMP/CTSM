program mksrfdat

!-----------------------------------------------------------------------
!BOP
!
! !PROGRAM: mksrfdat
!
! !DESCRIPTION:
! Creates land model surface dataset from original "raw" data files.
! Surface dataset contains model grid, pfts, inland water, glacier,
! soil texture, soil color, LAI and SAI, urban fraction, and urban
! parameters.
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8, r4 => shr_kind_r4
    use shr_sys_mod , only : shr_sys_getenv
    use shr_timer_mod
    use fileutils   , only : opnfil, getavu, get_filename
    use mklaiMod    , only : mklai
    use mkpftMod    , only : pft_idx, pft_frc, mkpft, mkpftInit, mkpft_parse_oride, &
                             mkirrig
    use mksoilMod   , only : soil_sand, soil_clay, mksoitex, mksoilInit, &
                             soil_color, mksoicol, mkorganic, soil_fmax, mkfmax
    use mkvocefMod  , only : mkvocef
    use mklanwatMod , only : mklanwat
    use mkglcmecMod , only : nglcec, mkglcmec, mkglcmecInit, mkglacier
    use mkharvestMod, only : mkharvest, mkharvest_init, mkharvest_fieldname, &
                             mkharvest_numtypes, mkharvest_parse_oride
    use mkurbanparMod, only : mkurban, mkurbanpar, mkelev
    use creategridMod, only : read_domain,write_domain
    use domainMod   , only : domain_setptrs, domain_type, domain_init
    use mkfileMod   , only : mkfile
    use mkvarpar    , only : nlevsoi, elev_thresh
    use mkvarsur    , only : spval, ldomain
    use mkvarctl
    use areaMod
    use ncdio       , only : check_ret, ncd_ioglobal
    use nanMod
!
! !ARGUMENTS:
    implicit none

    include 'netcdf.inc'
!
! !REVISION HISTORY:
! Authors: Gordon Bonan, Sam Levis and Mariana Vertenstein
! Revised: Nan Rosenbloom to add fmax processing.
! 3/18/08: David Lawrence added organic matter processing
! 1/22/09: Keith Oleson added urban parameter processing
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: lsmlon, lsmlat              ! clm grid resolution
    integer  :: nsoicol                     ! number of model color classes
    integer  :: i,j,k,m,n                   ! indices
    integer  :: ni, nj                      ! size of pft index
    integer  :: ier                         ! error status
    integer  :: ndiag,nfdyn                 ! unit numbers
    integer  :: ncid                        ! netCDF id
    integer  :: omode                       ! netCDF output mode
    integer  :: varid                       ! netCDF variable id
    integer  :: beg4d(4),len4d(4)           ! netCDF variable edges 4Dimensional
    integer  :: beg3d(3),len3d(3)           ! netCDF variable edges 3Dimensional
    integer  :: ret                         ! netCDF return status
    integer  :: ntim                        ! time sample for dynamic land use
    integer  :: year                        ! year for dynamic land use
    logical  :: all_veg                     ! if gridcell will be 100% vegetated land-cover
    real(r8) :: suma                        ! sum for error check
    character(len=256) :: fgrddat           ! grid data file
    character(len=256) :: fsurdat           ! surface data file name
    character(len=256) :: fsurlog           ! surface log file name
    character(len=256) :: fdyndat           ! dynamic landuse data file name
    character(len=256) :: fname             ! generic filename
    character(len=256) :: string            ! string read in
    character(len= 32) :: resol             ! resolution for file name
    integer  :: t1                          ! timer
    real(r8),parameter :: p5  = 0.5_r8      ! constant
    real(r8),parameter :: p25 = 0.25_r8     ! constant

    real(r8), allocatable  :: landfrac_pft(:,:)    ! PFT data: % land per gridcell
    real(r8), allocatable  :: pctlnd_pft(:,:)      ! PFT data: % of gridcell for PFTs
    real(r8), allocatable  :: pctlnd_pft_dyn(:,:)  ! PFT data: % of gridcell for dyn landuse PFTs
    integer , allocatable  :: pftdata_mask(:,:)    ! mask indicating real or fake land type
    real(r8), allocatable  :: pctpft(:,:,:)        ! PFT data: land fraction per gridcell
    real(r8), pointer      :: harvest(:,:,:)       ! harvest data: normalized harvesting
    real(r8), pointer      :: pctpft_i(:,:,:)      ! PFT data: % fraction on input grid
    real(r8), allocatable  :: pctgla(:,:)          ! percent of grid cell that is glacier  
    real(r8), allocatable  :: pctglcmec(:,:,:)     ! glacier_mec pct coverage in each gridcell and class
    real(r8), allocatable  :: topoglcmec(:,:,:)    ! glacier_mec sfc elevation in each gridcell and class
    real(r8), allocatable  :: thckglcmec(:,:,:)    ! glacier_mec ice sheet thcknss in each gridcell and class
    real(r8), allocatable  :: elevclass(:)         ! glacier_mec elevation classes
    real(r8), allocatable  :: pctlak(:,:)          ! percent of grid cell that is lake     
    real(r8), allocatable  :: pctwet(:,:)          ! percent of grid cell that is wetland  
    real(r8), allocatable  :: pctirr(:,:)          ! percent of grid cell that is irrigated  
    real(r8), allocatable  :: pcturb(:,:)          ! percent of grid cell that is urbanized
    real(r8), allocatable  :: elev(:,:)            ! elevation (m)
    real(r8), allocatable  :: fmax(:,:)            ! fractional saturated area
    integer , allocatable  :: soic2d(:,:)          ! soil color                            
    real(r8), allocatable  :: sand3d(:,:,:)        ! soil texture: percent sand            
    real(r8), allocatable  :: clay3d(:,:,:)        ! soil texture: percent clay            
    real(r8), allocatable  :: ef1_btr(:,:)         ! Isoprene emission factor for broadleaf
    real(r8), allocatable  :: ef1_fet(:,:)         ! Isoprene emission factor for fine/everg
    real(r8), allocatable  :: ef1_fdt(:,:)         ! Isoprene emission factor for fine/dec
    real(r8), allocatable  :: ef1_shr(:,:)         ! Isoprene emission factor for shrubs
    real(r8), allocatable  :: ef1_grs(:,:)         ! Isoprene emission factor for grasses
    real(r8), allocatable  :: ef1_crp(:,:)         ! Isoprene emission factor for crops
    real(r8), allocatable  :: organic3d(:,:,:)     ! organic matter density (kg/m3)            

    character(len=32) :: subname = 'mksrfdat'  ! program name

    namelist /clmexp/              &
	 mksrf_fgrid,              &	
	 mksrf_gridnm,             &	
	 mksrf_gridtype,           &	
         mksrf_fvegtyp,            &
	 mksrf_fsoitex,            &
         mksrf_forganic,           &
         mksrf_fsoicol,            &
         mksrf_fvocef,             &
         mksrf_flanwat,            &
         mksrf_fglacier,           &
         mksrf_ftopo,              &
         mksrf_ffrac,              &
         mksrf_fmax,               &
         mksrf_furban,             &
         mksrf_flai,               &
         mksrf_firrig,             &
         mksrf_fdynuse,            &
         outnc_large_files,        &
         outnc_double,             &
         numpft,                   &
         nglcec,                   &
         soil_color,               &
         soil_sand,                &
         soil_clay,                &
         soil_fmax,                &
         pft_idx,                  &
         pft_frc,                  &
         all_urban
!-----------------------------------------------------------------------

    ! ======================================================================
    ! Read input namelist
    ! ======================================
    ! Must specify settings for the output grid:
    ! ======================================
    !    mksrf_fgrid -- Grid dataset
    ! ======================================
    ! Must specify settings for input high resolution datafiles
    ! ======================================
    !    mksrf_ffrac ---- land fraction and land mask dataset
    !    mksrf_fglacier - Glacier dataset
    !    mksrf_flai ----- Leaf Area Index dataset
    !    mksrf_flanwat -- Land water dataset
    !    mksrf_forganic - Organic soil carbon dataset
    !    mksrf_fmax ----- Max fractional saturated area dataset
    !    mksrf_fsoicol -- Soil color dataset
    !    mksrf_fsoitex -- Soil texture dataset
    !    mksrf_ftopo ---- Topography dataset (for glacier multiple elevation classes)
    !                     (and for limiting urban areas)
    !    mksrf_furban --- Urban dataset
    !    mksrf_fvegtyp -- PFT vegetation type dataset
    !    mksrf_fvocef  -- Volatile Organic Compund Emission Factor dataset
    !    mksrf_fdynuse -- ASCII text file that lists each year of pft files to use
    ! ======================================
    ! Optionally specify setting for:
    ! ======================================
    !    mksrf_firrig ------ Irrigation dataset
    !    mksrf_gridtype ---- Type of grid (default is 'global')
    !    mksrf_gridnm ------ Name of output grid resolution
    !    outnc_double ------ If output should be in double precision
    !    outnc_large_files - If output should be in NetCDF large file format
    !    nglcec ------------ If you want to change the number of Glacier elevation classes
    ! ======================================
    ! Optional settings to change values for entire area
    ! ======================================
    !    all_urban --------- If entire area is urban
    !    soil_color -------- If you want to change the soil_color to this value everywhere
    !    soil_clay --------- If you want to change the soil_clay % to this value everywhere
    !    soil_fmax --------- If you want to change the soil_fmax  to this value everywhere
    !    soil_sand --------- If you want to change the soil_sand % to this value everywhere
    !    pft_idx ----------- If you want to change to 100% veg covered with given PFT indices
    !    pft_frc ----------- Fractions that correspond to the pft_idx above
    ! ==================
    !    numpft            (if different than default of 16)
    ! ======================================================================

    call shr_timer_init()
    call shr_timer_get(t1,'accumulating timer')
    call shr_timer_start(t1)

    write(6,*) 'Attempting to initialize control settings .....'

    mksrf_gridtype    = 'global'
    outnc_large_files = .false.
    outnc_double      = .true.
    all_urban         = .false.
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
       call abort()
    endif

    if (trim(mksrf_gridtype) == 'global' .or. &
        trim(mksrf_gridtype) == 'regional') then
       write(6,*)'mksrf_gridtype = ',trim(mksrf_gridtype)
    else
       write(6,*)'mksrf_gridtype = ',trim(mksrf_gridtype)
       write (6,*)'illegal mksrf_gridtype, must be global or regional '
       call abort()
    endif
    if ( outnc_large_files )then
       write(6,*)'Output file in NetCDF 64-bit large_files format'
    end if
    if ( outnc_double )then
       write(6,*)'Output ALL data in file as 64-bit'
    end if
    if ( all_urban )then
       write(6,*) 'Output ALL data in file as 100% urban'
    end if
    !
    ! Call module initialization routines
    !
    call mksoilInit( )
    call mkpftInit( all_urban, all_veg )
    allocate ( elevclass(nglcec+1) )
    call mkglcmecInit (elevclass)

    if ( all_veg )then
       write(6,*) 'Output ALL data in file as 100% vegetated'
    end if

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

    allocate ( landfrac_pft(lsmlon,lsmlat)      , &
               pctlnd_pft(lsmlon,lsmlat)        , & 
               pftdata_mask(lsmlon,lsmlat)      , & 
               pctpft(lsmlon,lsmlat,0:numpft)   , & 
               pctgla(lsmlon,lsmlat)            , & 
               pctlak(lsmlon,lsmlat)            , & 
               pctwet(lsmlon,lsmlat)            , & 
               pcturb(lsmlon,lsmlat)            , & 
               sand3d(lsmlon,lsmlat,nlevsoi)    , & 
               clay3d(lsmlon,lsmlat,nlevsoi)    , & 
               soic2d(lsmlon,lsmlat)              )
    landfrac_pft(:,:) = spval 
    pctlnd_pft(:,:)   = spval
    pftdata_mask(:,:) = -999
    pctpft(:,:,:)     = spval
    pctgla(:,:)       = spval
    pctlak(:,:)       = spval
    pctwet(:,:)       = spval
    pcturb(:,:)       = spval
    sand3d(:,:,:)     = spval
    clay3d(:,:,:)     = spval
    soic2d(:,:)       = -999

    write(6,*) ' timer_a2 init-----'
    call shr_timer_print(t1)

    ! ----------------------------------------------------------------------
    ! Open diagnostic output log file
    ! ----------------------------------------------------------------------
    
    if ( mksrf_gridnm == ' ' )then
       write (resol,'(i4.4,"x",i4.4)') lsmlat,lsmlon
    else
       resol = trim(mksrf_gridnm)
    end if
    fsurlog = './surfdata_'//trim(resol)//'.log'
    ndiag = getavu()
    call opnfil (fsurlog, ndiag, 'f')
    
    if (mksrf_fgrid /= ' ')then
       write (ndiag,*)'using fractional land data from file= ', &
            trim(mksrf_fgrid),' to create the surface dataset'
    endif

    if (trim(mksrf_gridtype) == 'global' .or. &
        trim(mksrf_gridtype) == 'regional') then
       write(6,*)'mksrf_gridtype = ',trim(mksrf_gridtype)
    endif

    write (ndiag,*) 'PFTs from:                 ',trim(mksrf_fvegtyp)
    write (ndiag,*) 'fmax from:                 ',trim(mksrf_fmax)
    write (ndiag,*) 'glaciers from:             ',trim(mksrf_fglacier)
    write (ndiag,*) 'topography from:           ',trim(mksrf_ftopo)
    write (ndiag,*) '           with:           ', nglcec, ' glacier elevation classes'
    write (ndiag,*) 'fracdata from:             ',trim(mksrf_ffrac)
    write (ndiag,*) 'urban from:                ',trim(mksrf_furban)
    write (ndiag,*) 'inland water from:         ',trim(mksrf_flanwat)
    write (ndiag,*) 'soil texture from:         ',trim(mksrf_fsoitex)
    write (ndiag,*) 'soil organic from:         ',trim(mksrf_forganic)
    write (ndiag,*) 'soil color from:           ',trim(mksrf_fsoicol)
    write (ndiag,*) 'VOC emission factors from: ',trim(mksrf_fvocef)
    if (mksrf_firrig /= ' ') then
       write (ndiag,*) &
                    'irrigated area from:       ',trim(mksrf_firrig)
    endif
    if (mksrf_fdynuse /= ' ') then
       write(6,*)'mksrf_fdynuse = ',trim(mksrf_fdynuse)
    else
       write (6,*) subname, ' error: fdynuse file is required'
       call abort()
    end if

    write(6,*) ' timer_a1 init-----'
    call shr_timer_print(t1)

    ! Make irrigated area fraction [pctirr] from [firrig] dataset if requested in namelist

    allocate ( pctirr(lsmlon,lsmlat) )
    pctirr(:,:)       = spval

    if (mksrf_firrig /= ' ') then

       call mkirrig (lsmlon, lsmlat, mksrf_firrig, ndiag, pctirr)

       write(6,*) ' timer_d mkirrig-----'
       call shr_timer_print(t1)

    endif

    ! Make PFTs [pctpft] from dataset [fvegtyp] (1/2 degree PFT data)

    nullify(pctpft_i)
    call mkpft(lsmlon, lsmlat, mksrf_fvegtyp, mksrf_firrig, ndiag, &
               pctlnd_pft, pctirr, pctpft, pctpft_i)

    write(6,*) ' timer_b mkpft-----'
    call shr_timer_print(t1)

    ! Make inland water [pctlak, pctwet] from Cogley's one degree data [flanwat]

    call mklanwat (lsmlon, lsmlat, mksrf_flanwat, ndiag, all_urban.or.all_veg, &
                   pctlak, pctwet)

    write(6,*) ' timer_c mklanwat-----'
    call shr_timer_print(t1)

    ! Make glacier fraction [pctgla] from [fglacier] dataset

    call mkglacier (lsmlon, lsmlat, mksrf_fglacier, ndiag, all_urban.or.all_veg, &
                    pctgla)

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

    call mkurban (lsmlon, lsmlat, mksrf_furban, ndiag, all_veg, pcturb)

    write(6,*) ' timer_g mkurban-----'
    call shr_timer_print(t1)

    ! Make elevation [elev] from [ftopo, ffrac] dataset
    ! Used only to screen pcturb

    allocate(  elev(lsmlon,lsmlat)  )
    elev(:,:)         = spval
    call mkelev (lsmlon, lsmlat, mksrf_ftopo, mksrf_ffrac, ndiag, elev, ncid)
    write(6,*) ' timer_f mkelev-----'
    call shr_timer_print(t1)

    ! Screen pcturb by elevation threshold from elev dataset
    if ( .not. all_urban )then
       where (elev .gt. elev_thresh)
         pcturb = 0._r8
       end where
    end if
    deallocate( elev )

    ! Make fmax [fmax] from [fmax] dataset
    allocate ( fmax(lsmlon,lsmlat) )
    fmax(:,:)         = spval

    call mkfmax (lsmlon, lsmlat, mksrf_fmax, ndiag, fmax)

    write(6,*) ' timer_d mkfmax-----'
    call shr_timer_print(t1)

    ! Make organic matter density [organic3d] from Global Soil Data Task [forganic]
    allocate ( organic3d(lsmlon,lsmlat,nlevsoi) )
    organic3d(:,:,:)  = spval

    call mkorganic (lsmlon, lsmlat, mksrf_forganic, ndiag, organic3d)

    write(6,*) ' timer_f mkorganic-----'
    call shr_timer_print(t1)

    ! Make VOC emission factors for isoprene [ef1_btr,ef1_fet,ef1_fdt,ef1_shr,ef1_grs,ef1_crp] 
    ! from MEGAN data [fvocef] (heald, 04/06)
    allocate ( ef1_btr(lsmlon,lsmlat)         , & 
               ef1_fet(lsmlon,lsmlat)         , & 
               ef1_fdt(lsmlon,lsmlat)         , & 
               ef1_shr(lsmlon,lsmlat)         , & 
               ef1_grs(lsmlon,lsmlat)         , & 
               ef1_crp(lsmlon,lsmlat) )
    ef1_btr(:,:)      = 0._r8
    ef1_fet(:,:)      = 0._r8
    ef1_fdt(:,:)      = 0._r8
    ef1_shr(:,:)      = 0._r8
    ef1_grs(:,:)      = 0._r8
    ef1_crp(:,:)      = 0._r8
    
    call mkvocef (lsmlon, lsmlat, mksrf_fvocef, ndiag, ef1_btr,ef1_fet,ef1_fdt,  &
                  ef1_shr,ef1_grs,ef1_crp)

    write(6,*) ' timer_f mkvocef-----'
    call shr_timer_print(t1)

    call change_landuse( lsmlon, lsmlat, dynpft=.false. )

    do j = 1,lsmlat
       do i = 1,lsmlon

          ! Assume wetland and/or lake when dataset landmask implies ocean 
          ! (assume medium soil color (15) and loamy texture).
          ! Also set pftdata_mask here

          if (pctlnd_pft(i,j) < 1.e-6_r8) then
             pftdata_mask(i,j) = 0

             soic2d(i,j)      = 15
             pctwet(i,j)      = 100._r8 - pctlak(i,j)
             pcturb(i,j)      = 0._r8
             if (mksrf_firrig /= ' ') &
             pctirr(i,j)      = 0._r8
             pctgla(i,j)      = 0._r8
             pctpft(i,j,:)    = 0._r8
             sand3d(i,j,:)    = 43._r8
             clay3d(i,j,:)    = 18._r8
             organic3d(i,j,:) = 0._r8
          else
             pftdata_mask(i,j) = 1
          end if

          ! Truncate all percentage fields on output grid. This is needed to
          ! insure that wt is not nonzero (i.e. a very small number such as
          ! 1e-16) where it really should be zero

          do k = 1,nlevsoi
             sand3d(i,j,k) = float(nint(sand3d(i,j,k)))
             clay3d(i,j,k) = float(nint(clay3d(i,j,k)))
          end do
          pctlak(i,j) = float(nint(pctlak(i,j)))
          pctwet(i,j) = float(nint(pctwet(i,j)))
          pctgla(i,j) = float(nint(pctgla(i,j)))
          if (mksrf_firrig /= ' ') &
          pctirr(i,j) = float(nint(pctirr(i,j)))

          ! Make sure sum of land cover types does not exceed 100. If it does,
          ! subtract excess from most dominant land cover.

          suma = pctlak(i,j) + pctwet(i,j) + pcturb(i,j) + pctgla(i,j)
          if (suma > 130._r4) then
             write (6,*) subname, ' error: sum of pctlak, pctwet,', &
                  'pcturb and pctgla is greater than 130%'
             write (6,*)'i,j,pctlak,pctwet,pcturb,pctgla= ', &
                  i,j,pctlak(i,j),pctwet(i,j),pcturb(i,j),pctgla(i,j)
             call abort()
          else if (suma > 100._r4) then
             pctlak(i,j) = pctlak(i,j) * 100._r8/suma
             pctwet(i,j) = pctwet(i,j) * 100._r8/suma
             pcturb(i,j) = pcturb(i,j) * 100._r8/suma
             pctgla(i,j) = pctgla(i,j) * 100._r8/suma
          end if

       end do
    end do

    call normalizencheck_landuse( lsmlon, lsmlat )

    ! Write out sum of PFT's

    do k = 0,numpft
       suma = 0._r8
       do j = 1,ldomain%nj
          do i = 1,ldomain%numlon(j)
             suma = suma + pctpft(i,j,k)
          enddo
       enddo
       write(6,*) 'sum over domain of pft ',k,suma
    enddo

    if ( nglcec > 0 )then
       ! Make glacier multiple elevation classes [pctglcmec,topoglcmec] from [fglacier,ftopo] dataset
       ! This call needs to occur after pctgla has been adjusted for the final time
       allocate ( pctglcmec(lsmlon,lsmlat,nglcec),  &
                  topoglcmec(lsmlon,lsmlat,nglcec), &
                  thckglcmec(lsmlon,lsmlat,nglcec) )

       pctglcmec(:,:,:)  = spval
       topoglcmec(:,:,:) = spval
       thckglcmec(:,:,:) = spval
       call mkglcmec (lsmlon, lsmlat, mksrf_ftopo, mksrf_ffrac, mksrf_fglacier, ndiag, pctgla, &
                      pctglcmec, topoglcmec, thckglcmec )

    end if

    write(6,*) ' timer_d mkglcmec-----'
    call shr_timer_print(t1)

    ! Determine fractional land from pft dataset

    do j = 1,lsmlat
       do i = 1,lsmlon
          landfrac_pft(i,j) = pctlnd_pft(i,j)/100._r8
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

    ! Write fields other than lai, sai, heights, and urban parameters to netcdf surface dataset

    call ncd_ioglobal(varname='PFTDATA_MASK', data=pftdata_mask, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='LANDFRAC_PFT', data=landfrac_pft, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='mxsoil_color', data=nsoicol     , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='SOIL_COLOR'  , data=soic2d      , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='PCT_SAND'    , data=sand3d      , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='PCT_CLAY'    , data=clay3d      , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='PCT_WETLAND' , data=pctwet      , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='PCT_LAKE'    , data=pctlak      , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='PCT_GLACIER' , data=pctgla      , ncid=ncid, flag='write')
    if ( nglcec > 0 )then
       call ncd_ioglobal(varname='PCT_GLC_MEC' , data=pctglcmec   , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='GLC_MEC'     , data=elevclass   , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='TOPO_GLC_MEC', data=topoglcmec  , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='THCK_GLC_MEC', data=thckglcmec  , ncid=ncid, flag='write')
    end if
    call ncd_ioglobal(varname='PCT_URBAN'   , data=pcturb      , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='PCT_PFT'     , data=pctpft      , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='FMAX'        , data=fmax        , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='EF1_BTR'     , data=ef1_btr     , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='EF1_FET'     , data=ef1_fet     , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='EF1_FDT'     , data=ef1_fdt     , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='EF1_SHR'     , data=ef1_shr     , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='EF1_GRS'     , data=ef1_grs     , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='EF1_CRP'     , data=ef1_crp     , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='ORGANIC'     , data=organic3d   , ncid=ncid, flag='write')
    if (mksrf_firrig /= ' ') then
       call ncd_ioglobal(varname='PCT_IRRIG', data=pctirr      , ncid=ncid, flag='write')
    endif
    ! Deallocate arrays NOT needed for dynamic-pft section of code
    deallocate ( organic3d )
    deallocate ( ef1_btr, ef1_fet, ef1_fdt, ef1_shr, ef1_grs, ef1_crp )
    if ( nglcec > 0 ) deallocate ( pctglcmec, topoglcmec, thckglcmec )
    deallocate ( elevclass )
    deallocate ( fmax )
    deallocate ( sand3d, clay3d )
    deallocate ( soic2d )

    ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

    call check_ret(nf_sync(ncid), subname)

    write(6,*) ' timer_i writesurf-----'
    call shr_timer_print(t1)

    ! ----------------------------------------------------------------------
    ! Make Urban Parameters from 1/2 degree data and write to surface dataset 
    ! Write to netcdf file is done inside mkurbanpar routine
    ! Only call this routine if pcturb is greater than zero somewhere.  Raw urban
    ! datasets will have no associated parameter fields if there is no urban 
    ! (e.g., mksrf_urban.060929.nc).
    ! ----------------------------------------------------------------------

    if (any(pcturb > 0._r8)) then
      call mkurbanpar(lsmlon, lsmlat, mksrf_furban, ndiag, ncid)
    else
      write(6,*) 'PCT_URBAN is zero everywhere, no urban parameter fields will be created'
    end if

    write(6,*) ' timer_j mkurbanpar-----'
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
    call mklai(lsmlon, lsmlat, mksrf_flai, mksrf_firrig, &
               ndiag, ncid, ni, nj, pctpft_i)
    deallocate( pctpft_i )

    write(6,*) ' timer_k mklai-----'
    call shr_timer_print(t1)

    ! Close surface dataset

    call check_ret(nf_close(ncid), subname)

    write (6,'(72a1)') ("-",i=1,60)
    write (6,'(a,f5.1,a4,f5.1,a5)') 'land model surface data set successfully created for ', &
         360._r8/lsmlon,' by ',180._r8/lsmlat,' grid'

    ! ----------------------------------------------------------------------
    ! Create dynamic land use dataset
    ! ----------------------------------------------------------------------

    allocate( pctlnd_pft_dyn(lsmlon,lsmlat) )
    call mkharvest_init( lsmlon, lsmlat, spval, harvest, mksrf_fvegtyp )

    fdyndat = './surfdata.pftdyn_'//trim(resol)//'.nc'

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

    ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

    call check_ret(nf_sync(ncid), subname)

    ! Read in each dynamic pft landuse dataset

    nfdyn = getavu(); call opnfil (mksrf_fdynuse, nfdyn, 'f')

    ntim = 0
    do 
       ! Read input pft data

       read(nfdyn, '(A125,1x,I4)', iostat=ier) string, year
       if (ier /= 0) exit
       !
       ! If pft fraction override is set, than intrepret string as PFT and harvesting override values
       !
       if ( any(pft_frc > 0.0_r8 ) )then
          fname = ' '
          call mkpft_parse_oride(     string )
          call mkharvest_parse_oride( string )
          write(6,*)'PFT and harvesting values are ',trim(string),' year is ',year
          !
          ! Otherwise intrepret string as a filename with PFT and harvesting values in it
          !
       else
          fname = string
          write(6,*)'input pft dynamic dataset is  ',trim(fname),' year is ',year
       end if
       ntim = ntim + 1

       ! Create pctpft data at model resolution

       call mkpft(lsmlon, lsmlat, fname, mksrf_firrig, ndiag, &
                  pctlnd_pft_dyn, pctirr, pctpft)

       ! Create harvesting data at model resolution

       call mkharvest( lsmlon, lsmlat, fname, ndiag, harvest )

       ! Consistency check on input land fraction

       do j = 1,lsmlat
          do i = 1,lsmlon

             if (pctlnd_pft_dyn(i,j) /= pctlnd_pft(i,j)) then
                write(6,*) subname,' error: pctlnd_pft for dynamics data = ',&
                     pctlnd_pft_dyn(i,j), ' not equal to pctlnd_pft for surface data = ',&
                     pctlnd_pft(i,j),' at i,j= ',i,j
                if ( trim(fname) == ' ' )then
                   write(6,*) ' PFT string = ', string
                else
                   write(6,*) ' PFT file = ', fname
                end if
                call abort()
             end if

          end do
       end do

       call change_landuse( lsmlon, lsmlat, dynpft=.true. )

       call normalizencheck_landuse( lsmlon, lsmlat )

       ! Output pctpft data for current year

       beg4d(1) = 1     ;  len4d(1) = lsmlon
       beg4d(2) = 1     ;  len4d(2) = lsmlat
       beg4d(3) = 1     ;  len4d(3) = numpft+1
       beg4d(4) = ntim  ;  len4d(4) = 1

       call check_ret(nf_inq_varid(ncid, 'PCT_PFT', varid), subname)
       call check_ret(nf_put_vara_double(ncid, varid, beg4d, len4d, pctpft), subname)

    
       beg3d(1) = 1     ;  len3d(1) = lsmlon
       beg3d(2) = 1     ;  len3d(2) = lsmlat
       beg3d(3) = ntim  ;  len3d(3) = 1
       do k = 1, mkharvest_numtypes()
          call check_ret(nf_inq_varid(ncid, trim(mkharvest_fieldname(k)), varid), subname)
          call check_ret(nf_put_vara_double(ncid, varid, beg3d, len3d, harvest(:,:,k) ), subname)
       end do

       call check_ret(nf_inq_varid(ncid, 'YEAR', varid), subname)
       call check_ret(nf_put_vara_int(ncid, varid, ntim, 1, year), subname)

       call check_ret(nf_inq_varid(ncid, 'time', varid), subname)
       call check_ret(nf_put_vara_int(ncid, varid, ntim, 1, year), subname)

       call check_ret(nf_inq_varid(ncid, 'input_pftdata_filename', varid), subname)
       call check_ret(nf_put_vara_text(ncid, varid, (/ 1, ntim /), (/ len_trim(string), 1 /), trim(string) ), subname)

       ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

       call check_ret(nf_sync(ncid), subname)

    end do   ! end of read loop

    call check_ret(nf_close(ncid), subname)

    write(6,*) ' timer_l writedyn-----'
    call shr_timer_print(t1)

    ! ----------------------------------------------------------------------
    ! Close diagnostic dataset
    ! ----------------------------------------------------------------------

    close (ndiag)
    write (6,*)
    write (6,*) 'Surface data output file = ',trim(fsurdat)
    write (6,*) '   This file contains the land model surface data'
    write (6,*) 'Diagnostic log file      = ',trim(fsurlog)
    write (6,*) '   See this file for a summary of the dataset'
    write (6,*)

    write(6,*) ' timer_z end-----'
    call shr_timer_print(t1)
    write(6,*) 'Successfully created surface dataset'

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: change_landuse
!
! !INTERFACE:
subroutine change_landuse( lsmlon, lsmlat, dynpft )
!
! !DESCRIPTION:
!
! Do landuse changes such as for the poles, Ross ice-shelf and etc.
!
! !USES:
    use mkvarsur    , only : ldomain
    implicit none
!
! !ARGUMENTS:
    integer, intent(in)  :: lsmlon, lsmlat  ! clm grid resolution
    logical, intent(in)  :: dynpft          ! if part of the dynpft section of code

!
! !REVISION HISTORY:
! 9/10/09: Erik Kluzek spin off subroutine from original embedded code
!
!EOP
!
! !LOCAL VARIABLES:
    logical  :: first_time                  ! flag if this is the first pass through or not
    integer ,parameter :: bdtroptree = 6    ! Index for broadleaf decidious tropical tree
    integer ,parameter :: bdtemptree = 7    ! Index for broadleaf decidious temperate tree
    integer ,parameter :: bdtempshrub = 10  ! Index for broadleaf decidious temperate shrub
    real(r8),parameter :: troplat = 23.5_r8 ! Latitude to define as tropical
    real(r8),parameter :: rosslat = -79._r8 ! Latitude to define as Ross ice-shelf
    integer  :: i,j                         ! indices
    character(len=32) :: subname = 'change_landuse'  ! subroutine name
!-----------------------------------------------------------------------

    first_time = .true.
    do j = 1,lsmlat
       do i = 1,lsmlon

          ! Set pfts 7 and 10 to 6 in the tropics to avoid lais > 1000
          ! Using P. Thornton's method found in surfrdMod.F90 in clm3.5

          if (abs(ldomain%latixy(i,j))<troplat .and. pctpft(i,j,bdtemptree)>0._r8) then
             pctpft(i,j,bdtroptree) = pctpft(i,j,bdtroptree) + pctpft(i,j,bdtemptree)
             pctpft(i,j,bdtemptree) = 0._r8
             if ( first_time ) write (6,*) subname, ' Warning: all wgt of pft ', bdtemptree, ' now added to pft ', bdtroptree
          end if
          if (abs(ldomain%latixy(i,j))<troplat .and. pctpft(i,j,bdtempshrub)>0._r8) then
             pctpft(i,j,bdtroptree) = pctpft(i,j,bdtroptree) + pctpft(i,j,bdtempshrub)
             pctpft(i,j,bdtempshrub) = 0._r8
             if ( first_time ) write (6,*) subname, ' Warning: all wgt of pft ', bdtempshrub, ' now added to pft ', bdtroptree
          end if
          first_time = .false.

          ! Set land values on Ross ice shelf to glacier

          if (ldomain%latixy(i,j) < rosslat) then

             pctpft(i,j,:)  = 0._r8

             pctlak(i,j)      = 0._r8
             pctwet(i,j)      = 0._r8
             pcturb(i,j)      = 0._r8

             pctgla(i,j)      = 100._r8

          end if

       end do
    end do
    ! Set land values on Ross ice shelf to glacier -- non-dynamic-PFT part of the code
    if ( .not. dynpft )then
       do j = 1,lsmlat
          do i = 1,lsmlon
             if (ldomain%latixy(i,j) < rosslat) then
                ef1_btr(i,j)     = 0._r8
                ef1_fet(i,j)     = 0._r8
                ef1_fdt(i,j)     = 0._r8
                ef1_shr(i,j)     = 0._r8
                ef1_grs(i,j)     = 0._r8
                ef1_crp(i,j)     = 0._r8
                organic3d(i,j,:) = 0._r8
                soic2d(i,j)      = 0
                if (mksrf_firrig /= ' ') &
                pctirr(i,j)      = 0._r8
                sand3d(i,j,:)    = 0._r8
                clay3d(i,j,:)    = 0._r8
             end if
          end do
       end do
    end if

    ! If have pole points on grid - set south pole to glacier
    ! north pole is assumed as non-land

    if (abs((ldomain%latixy(1,lsmlat) - 90._r8)) < 1.e-6_r8) then
       write(6,*) subname, ': grid has pole_points'
       do i = 1,lsmlon
          pctpft(i,1,:)    = 0._r8

          pctlak(i,1)      = 0._r8
          pctwet(i,1)      = 0._r8
          pcturb(i,1)      = 0._r8

          pctgla(i,1)      = 100._r8
       end do
       if ( .not. dynpft )then
          do i = 1,lsmlon
             organic3d(i,1,:) = 0._r8
             ef1_btr(i,1)     = 0._r8
             ef1_fet(i,1)     = 0._r8
             ef1_fdt(i,1)     = 0._r8
             ef1_shr(i,1)     = 0._r8
             ef1_grs(i,1)     = 0._r8
             ef1_crp(i,1)     = 0._r8
             soic2d(i,1)      = 0
             if (mksrf_firrig /= ' ') &
             pctirr(i,1)      = 0._r8
             sand3d(i,1,:)    = 0._r8
             clay3d(i,1,:)    = 0._r8
          end do
       end if
    end if
!-----------------------------------------------------------------------
end subroutine change_landuse

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: normalizencheck_landuse
!
! !INTERFACE:
subroutine normalizencheck_landuse( lsmlon, lsmlat )
!
! !DESCRIPTION:
!
! Normalize land use and make sure things add up to 100% as well as
! checking that things are as they should be.
!
! !USES:
    use mkvarsur    , only : ldomain
    implicit none
! !ARGUMENTS:
    integer, intent(in)  :: lsmlon, lsmlat  ! clm grid resolution

!
! !REVISION HISTORY:
! 9/10/09: Erik Kluzek spin off subroutine from original embedded code
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,m,k                     ! indices
    real(r8) :: suma                        ! sum for error check
    real(r8) :: bare_urb_diff               ! difference between bare soil and urban %
    real(r8) :: pcturb_excess               ! excess urban % not accounted for by bare soil
    real(r8) :: sumpft                      ! sum of non-baresoil pfts
    real(r8) :: sum8, sum8a                 ! sum for error check
    real(r4) :: sum4a                       ! sum for error check
    character(len=32) :: subname = 'normalizencheck_landuse'  ! subroutine name
!-----------------------------------------------------------------------

    do j = 1,lsmlat
       do i = 1,lsmlon

          if (pcturb(i,j) .gt. 0._r8) then
          ! Replace bare soil preferentially with urban

             suma = pctlak(i,j)+pctwet(i,j)+pctgla(i,j)
             bare_urb_diff = 0.01_r8 * pctpft(i,j,0) * (100._r8 - suma) - pcturb(i,j)
             pctpft(i,j,0) = max(0._r8,bare_urb_diff)
             pcturb_excess = abs(min(0._r8,bare_urb_diff))

             ! Normalize pctpft to be the remainder of [100 - (special landunits)]
             ! including any urban not accounted for by bare soil above

             sumpft = sum(pctpft(i,j,1:numpft))
             if (sumpft > 0._r8) then
                suma = pctlak(i,j)+pctwet(i,j)+pctgla(i,j)
                do m = 1, numpft
                   pctpft(i,j,m) = 0.01_r8 * pctpft(i,j,m) * (100._r8 - suma) - &
                                   pcturb_excess*pctpft(i,j,m)/sumpft
                end do
             end if
          else

             ! Normalize pctpft to be the remainder of [100 - (special landunits)]

             suma = pctlak(i,j) + pctwet(i,j) + pcturb(i,j) + pctgla(i,j)
             do m = 0, numpft
                pctpft(i,j,m) = 0.01_r8 * pctpft(i,j,m) * (100._r8 - suma)
             end do
          end if

          suma = pctlak(i,j) + pctwet(i,j) + pcturb(i,j) + pctgla(i,j)
          do m = 0,numpft
             suma = suma + pctpft(i,j,m)
          end do

          if (suma < 90._r8) then
             write (6,*) subname, ' error: sum of pctlak, pctwet,', &
                  'pcturb, pctgla and pctpft is less than 90'
             write (6,*)'i,j,pctlak,pctwet,pcturb,pctgla,pctpft= ', &
                  i,j,pctlak(i,j),pctwet(i,j),pcturb(i,j),pctgla(i,j),&
                  pctpft(i,j,:)
             call abort()
          else if (suma > 100._r8 + 1.e-4_r8) then
             write (6,*) subname, ' error: sum of pctlak, pctwet,', &
                  'pcturb, pctgla and pctpft is greater than 100'
             write (6,*)'i,j,pctlak,pctwet,pcturb,pctgla,pctpft,sum= ', &
                  i,j,pctlak(i,j),pctwet(i,j),pcturb(i,j),pctgla(i,j),&
                  pctpft(i,j,:), suma
             call abort()
          else
             pctlak(i,j)   = pctlak(i,j)   * 100._r8/suma
             pctwet(i,j)   = pctwet(i,j)   * 100._r8/suma
             pcturb(i,j)   = pcturb(i,j)   * 100._r8/suma
             pctgla(i,j)   = pctgla(i,j)   * 100._r8/suma
             pctpft(i,j,:) = pctpft(i,j,:) * 100._r8/suma
          end if

          ! Roundoff error fix
          suma = pctlak(i,j) + pctwet(i,j) + pcturb(i,j) + pctgla(i,j)
          if (suma < 100._r8 .and. suma > (100._r8 - 100._r8*epsilon(suma))) then
             write (6,*) 'Special land units near 100%, but not quite for i,j,suma =',i,j,suma
             write (6,*) 'Adjusting special land units to 100%'
             if (pctlak(i,j) >= 25._r8) then
                pctlak(i,j) = 100._r8 - (pctwet(i,j) + pcturb(i,j) + pctgla(i,j))
             else if (pctwet(i,j) >= 25._r8) then
                pctwet(i,j) = 100._r8 - (pctlak(i,j) + pcturb(i,j) + pctgla(i,j))
             else if (pcturb(i,j) >= 25._r8) then
                pcturb(i,j) = 100._r8 - (pctlak(i,j) + pctwet(i,j) + pctgla(i,j))
             else if (pctgla(i,j) >= 25._r8) then
                pctgla(i,j) = 100._r8 - (pctlak(i,j) + pctwet(i,j) + pcturb(i,j))
             else
                write (6,*) subname, 'Error: sum of special land units nearly 100% but none is >= 25% at ', &
                   'i,j,pctlak(i,j),pctwet(i,j),pcturb(i,j),pctgla(i,j),pctpft(i,j,:),suma = ', &
                   i,j,pctlak(i,j),pctwet(i,j),pcturb(i,j),pctgla(i,j),pctpft(i,j,:),suma
                call abort()
             end if
             pctpft(i,j,:) = 0._r8
          end if

          suma = pctlak(i,j) + pctwet(i,j) + pcturb(i,j) + pctgla(i,j)
          if (suma < 100._r8-epsilon(suma) .and. suma > (100._r8 - 4._r8*epsilon(suma))) then
             write (6,*) subname, 'i,j,pctlak,pctwet,pcturb,pctgla,pctpft= ', &
                  i,j,pctlak(i,j),pctwet(i,j),pcturb(i,j),pctgla(i,j),&
                  pctpft(i,j,:)
             call abort()
          end if
          do m = 0,numpft
             suma = suma + pctpft(i,j,m)
          end do
          if ( abs(suma-100._r8) > 1.e-10_r8) then
             write (6,*) subname, ' error: sum of pctlak, pctwet,', &
                  'pcturb, pctgla and pctpft is NOT equal to 100'
             write (6,*)'i,j,pctlak,pctwet,pcturb,pctgla,pctpft,sum= ', &
                  i,j,pctlak(i,j),pctwet(i,j),pcturb(i,j),pctgla(i,j),&
                  pctpft(i,j,:), sum8
             call abort()
          end if

       end do
    end do

    ! Check that when pctpft identically zero, sum of special landunits is identically 100%

    if ( .not. outnc_double )then
       do j = 1,ldomain%nj
          do i = 1,ldomain%numlon(j)
            sum8  =        real(pctlak(i,j),r4)
            sum8  = sum8 + real(pctwet(i,j),r4)
            sum8  = sum8 + real(pcturb(i,j),r4)
            sum8  = sum8 + real(pctgla(i,j),r4)
            sum4a = 0.0_r4
            do k = 0,numpft
               sum4a = sum4a + real(pctpft(i,j,k),r4)
            end do
            if ( sum4a==0.0_r4 .and. sum8 < 100._r4-2._r4*epsilon(sum4a) )then
               write (6,*) subname, ' error: sum of pctlak, pctwet,', &
                    'pcturb, pctgla is < 100% when pctpft==0 sum = ', sum8
               write (6,*)'i,j,pctlak,pctwet,pcturb,pctgla= ', &
                    i,j,pctlak(i,j),pctwet(i,j),pcturb(i,j),pctgla(i,j), pctpft(i,j,:)
               call abort()
            end if
          enddo
       enddo
    else
       do j = 1,ldomain%nj
          do i = 1,ldomain%numlon(j)
            sum8  =        pctlak(i,j)
            sum8  = sum8 + pctwet(i,j)
            sum8  = sum8 + pcturb(i,j)
            sum8  = sum8 + pctgla(i,j)
            sum8a = 0._r8
            do k = 0,numpft
               sum8a = sum8a + pctpft(i,j,k)
            end do
            if ( sum8a==0._r8 .and. sum8 < (100._r8-4._r8*epsilon(sum8)) )then
               write (6,*) subname, ' error: sum of pctlak, pctwet,', &
                    'pcturb, pctgla is < 100% when pctpft==0 sum = ', sum8
               write (6,*) 'Total error, error/epsilon = ',100._r8-sum8, ((100._r8-sum8)/epsilon(sum8))
               write (6,*)'i,j,pctlak,pctwet,pcturb,pctgla,epsilon= ', &
                    i,j,pctlak(i,j),pctwet(i,j),pcturb(i,j),pctgla(i,j), pctpft(i,j,:), epsilon(sum8)
               call abort()
            end if
          enddo
       enddo
    end if

!-----------------------------------------------------------------------
end subroutine normalizencheck_landuse

end program mksrfdat
