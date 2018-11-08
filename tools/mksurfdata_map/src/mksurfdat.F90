program mksurfdat

!-----------------------------------------------------------------------
!BOP
!
! !PROGRAM: mksurfdat
!
! !DESCRIPTION:
! Creates land model surface dataset from original "raw" data files.
! Surface dataset contains model grid, pfts, inland water, glacier,
! soil texture, soil color, LAI and SAI, urban fraction, and urban
! parameters.
!
! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8, r4 => shr_kind_r4
    use fileutils          , only : opnfil, getavu
    use mklaiMod           , only : mklai
    use mkpctPftTypeMod    , only : pct_pft_type, get_pct_p2l_array, get_pct_l2g_array, update_max_array
    use mkpftConstantsMod  , only : natpft_lb, natpft_ub, cft_lb, cft_ub, num_cft
    use mkpftMod           , only : pft_idx, pft_frc, mkpft, mkpftInit, mkpft_parse_oride
    use mksoilMod          , only : soil_sand, soil_clay, mksoiltex, mksoilInit, &
                                    soil_color, mksoilcol, mkorganic, &
                                    soil_fmax, mkfmax
    use mkvocefMod         , only : mkvocef
    use mklanwatMod        , only : mklakwat, mkwetlnd, mklakparams
    use mkglacierregionMod , only : mkglacierregion
    use mkglcmecMod        , only : nglcec, mkglcmec, mkglcmecInit, mkglacier
    use mkharvestMod       , only : mkharvest, mkharvest_init, mkharvest_fieldname
    use mkharvestMod       , only : mkharvest_numtypes, mkharvest_parse_oride
    use mkharvestMod       , only : harvestDataType
    use mkurbanparCommonMod, only : mkelev
    use mkurbanparMod      , only : mkurbanInit, mkurban, mkurbanpar, numurbl
    use mkutilsMod         , only : normalize_classes_by_gcell
    use mkfileMod          , only : mkfile
    use mkvarpar           , only : nlevsoi, elev_thresh
    use mkvarctl
    use nanMod             , only : nan, bigint
    use mkncdio            , only : check_ret, ncd_put_time_slice
    use mkdomainMod        , only : domain_type, domain_read_map, domain_read, &
                                    domain_write
    use mkgdpMod           , only : mkgdp
    use mkpeatMod          , only : mkpeat
    use mksoildepthMod          , only : mksoildepth
    use mkagfirepkmonthMod , only : mkagfirepkmon
    use mktopostatsMod     , only : mktopostats
    use mkVICparamsMod     , only : mkVICparams
    use mkCH4inversionMod  , only : mkCH4inversion
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
! 2/11/13: Sam Levis added abm, peat, and gdp processing for new fire model
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: nsoicol                     ! number of model color classes
    integer  :: k,m,n                       ! indices
    integer  :: ni,nj,ns_o                  ! indices
    integer  :: ier                         ! error status
    integer  :: ndiag,nfdyn                 ! unit numbers
    integer  :: ncid                        ! netCDF id
    integer  :: omode                       ! netCDF output mode
    integer  :: varid                       ! netCDF variable id
    integer  :: ret                         ! netCDF return status
    integer  :: ntim                        ! time sample for dynamic land use
    integer  :: year                        ! year for dynamic land use
    integer  :: year2                       ! year for dynamic land use for harvest file
    logical  :: all_veg                     ! if gridcell will be 100% vegetated land-cover
    real(r8) :: suma                        ! sum for error check
    character(len=256) :: fgrddat           ! grid data file
    character(len=256) :: fsurdat           ! output surface data file name (if blank, do not output a surface dataset)
    character(len=256) :: fsurlog           ! output surface log file name
    character(len=256) :: fdyndat           ! dynamic landuse data file name
    character(len=256) :: fname             ! generic filename
    character(len=256) :: fhrvname          ! generic harvest filename
    character(len=256) :: string            ! string read in
    integer  :: t1                          ! timer
    real(r8),parameter :: p5  = 0.5_r8      ! constant
    real(r8),parameter :: p25 = 0.25_r8     ! constant

    real(r8), allocatable  :: landfrac_pft(:)    ! PFT data: % land per gridcell
    real(r8), allocatable  :: pctlnd_pft(:)      ! PFT data: % of gridcell for PFTs
    real(r8), allocatable  :: pctlnd_pft_dyn(:)  ! PFT data: % of gridcell for dyn landuse PFTs
    integer , allocatable  :: pftdata_mask(:)    ! mask indicating real or fake land type
    type(pct_pft_type), allocatable :: pctnatpft(:)     ! % of grid cell that is nat veg, and breakdown into PFTs
    type(pct_pft_type), allocatable :: pctnatpft_max(:) ! % of grid cell maximum PFTs of the time series
    type(pct_pft_type), allocatable :: pctcft(:)        ! % of grid cell that is crop, and breakdown into CFTs
    type(pct_pft_type), allocatable :: pctcft_max(:)    ! % of grid cell maximum CFTs of the time series
    type(pct_pft_type), allocatable :: pctcft_saved(:)  ! version of pctcft saved from the initial call to mkpft
    real(r8), pointer      :: harvest1D(:)       ! harvest 1D data: normalized harvesting
    real(r8), pointer      :: harvest2D(:,:)     ! harvest 1D data: normalized harvesting
    real(r8), allocatable  :: pctgla(:)          ! percent of grid cell that is glacier  
    real(r8), allocatable  :: pctglc_gic(:)      ! percent of grid cell that is gic (% of glc landunit)
    real(r8), allocatable  :: pctglc_icesheet(:) ! percent of grid cell that is ice sheet (% of glc landunit)
    real(r8), allocatable  :: pctglcmec(:,:)     ! glacier_mec pct coverage in each class (% of landunit)
    real(r8), allocatable  :: topoglcmec(:,:)    ! glacier_mec sfc elevation in each gridcell and class
    real(r8), allocatable  :: pctglcmec_gic(:,:) ! GIC pct coverage in each class (% of landunit)
    real(r8), allocatable  :: pctglcmec_icesheet(:,:) ! icesheet pct coverage in each class (% of landunit)
    real(r8), allocatable  :: elevclass(:)       ! glacier_mec elevation classes
    integer,  allocatable  :: glacier_region(:)  ! glacier region ID
    real(r8), allocatable  :: pctlak(:)          ! percent of grid cell that is lake
    real(r8), allocatable  :: pctwet(:)          ! percent of grid cell that is wetland  
    real(r8), allocatable  :: pcturb(:)          ! percent of grid cell that is urbanized (total across all urban classes)
    real(r8), allocatable  :: urbn_classes(:,:)  ! percent cover of each urban class, as % of total urban area
    real(r8), allocatable  :: urbn_classes_g(:,:)! percent cover of each urban class, as % of grid cell
    real(r8), allocatable  :: elev(:)            ! glc elevation (m)
    real(r8), allocatable  :: fmax(:)            ! fractional saturated area
    integer , allocatable  :: soicol(:)          ! soil color                            
    real(r8), allocatable  :: pctsand(:,:)       ! soil texture: percent sand            
    real(r8), allocatable  :: pctclay(:,:)       ! soil texture: percent clay            
    real(r8), allocatable  :: ef1_btr(:)         ! Isoprene emission factor for broadleaf
    real(r8), allocatable  :: ef1_fet(:)         ! Isoprene emission factor for fine/everg
    real(r8), allocatable  :: ef1_fdt(:)         ! Isoprene emission factor for fine/dec
    real(r8), allocatable  :: ef1_shr(:)         ! Isoprene emission factor for shrubs
    real(r8), allocatable  :: ef1_grs(:)         ! Isoprene emission factor for grasses
    real(r8), allocatable  :: ef1_crp(:)         ! Isoprene emission factor for crops
    real(r8), allocatable  :: organic(:,:)       ! organic matter density (kg/m3)            
    real(r8), allocatable  :: gdp(:)             ! GDP (x1000 1995 US$/capita)
    real(r8), allocatable  :: fpeat(:)           ! peatland fraction of gridcell
    real(r8), allocatable  :: soildepth(:)       ! soil depth (m)
    integer , allocatable  :: agfirepkmon(:)     ! agricultural fire peak month
    integer , allocatable  :: urban_region(:)    ! urban region ID
    real(r8), allocatable  :: topo_stddev(:)     ! standard deviation of elevation (m)
    real(r8), allocatable  :: slope(:)           ! topographic slope (degrees)
    real(r8), allocatable  :: vic_binfl(:)       ! VIC b parameter (unitless)
    real(r8), allocatable  :: vic_ws(:)          ! VIC Ws parameter (unitless)
    real(r8), allocatable  :: vic_dsmax(:)       ! VIC Dsmax parameter (mm/day)
    real(r8), allocatable  :: vic_ds(:)          ! VIC Ds parameter (unitless)
    real(r8), allocatable  :: lakedepth(:)       ! lake depth (m)
    real(r8), allocatable  :: f0(:)              ! max fractional inundated area (unitless)
    real(r8), allocatable  :: p3(:)              ! coefficient for qflx_surf_lag for finundated (s/mm)
    real(r8), allocatable  :: zwt0(:)            ! decay factor for finundated (m)

    real(r8) :: std_elev = -999.99_r8            ! Standard deviation of elevation (m) to use for entire grid

    integer, allocatable :: harvind1D(:)         ! Indices of 1D harvest fields
    integer, allocatable :: harvind2D(:)         ! Indices of 2D harvest fields

    ! NOTE(bja, 2015-01) added to work around a ?bug? causing 1x1_urbanc_alpha to abort. See
    !/glade/p/cesm/cseg/inputdata/lnd/clm2/surfdata_map/README_c141219
    logical :: urban_skip_abort_on_invalid_data_check

    type(domain_type) :: ldomain

    character(len=32) :: subname = 'mksrfdat'  ! program name
    type(harvestDataType) :: harvdata

    namelist /clmexp/              &
	 mksrf_fgrid,              &	
	 mksrf_gridtype,           &	
         mksrf_fvegtyp,            &
         mksrf_fhrvtyp,            &
	 mksrf_fsoitex,            &
         mksrf_forganic,           &
         mksrf_fsoicol,            &
         mksrf_fvocef,             &
         mksrf_flakwat,            &
         mksrf_fwetlnd,            &
         mksrf_fglacier,           &
         mksrf_fglacierregion,     &
         mksrf_furbtopo,           &
         mksrf_fmax,               &
         mksrf_furban,             &
         mksrf_flai,               &
         mksrf_fdynuse,            &
         mksrf_fgdp,               &
         mksrf_fpeat,              &
         mksrf_fsoildepth,         &
         mksrf_fabm,               &
         mksrf_ftopostats,         &
         mksrf_fvic,               &
         mksrf_fch4,               &
         nglcec,                   &
         numpft,                   &
         soil_color,               &
         soil_sand,                &
         soil_fmax,                &
         soil_clay,                &
         pft_idx,                  &
         pft_frc,                  &
         all_urban,                &
         no_inlandwet,             &
         map_fpft,                 &
         map_flakwat,              &
         map_fwetlnd,              &
         map_fglacier,             &
         map_fglacierregion,       &
         map_fsoitex,              &
         map_fsoicol,              &
         map_furban,               &
         map_furbtopo,             &
         map_fmax,                 &
         map_forganic,             &
         map_fvocef,               &
         map_flai,                 &
         map_fharvest,             &
         map_fgdp,                 &
         map_fpeat,                &
         map_fsoildepth,           &
         map_fabm,                 &
         map_ftopostats,           &
         map_fvic,                 &
         map_fch4,                 &
         gitdescribe,              &
         outnc_large_files,        &
         outnc_double,             &
         outnc_dims,               &
         fsurdat,                  &
         fdyndat,                  &   
         fsurlog,                  &
         std_elev,                 &
         urban_skip_abort_on_invalid_data_check

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
    !    mksrf_fglacier - Glacier dataset
    !    mksrf_fglacierregion - Glacier region ID dataset
    !    mksrf_flai ----- Leaf Area Index dataset
    !    mksrf_flakwat -- Lake water dataset
    !    mksrf_fwetlnd -- Wetland water dataset
    !    mksrf_forganic - Organic soil carbon dataset
    !    mksrf_fmax ----- Max fractional saturated area dataset
    !    mksrf_fsoicol -- Soil color dataset
    !    mksrf_fsoitex -- Soil texture dataset
    !    mksrf_furbtopo-- Topography dataset (for limiting urban areas)
    !    mksrf_furban --- Urban dataset
    !    mksrf_fvegtyp -- PFT vegetation type dataset
    !    mksrf_fhrvtyp -- harvest type dataset
    !    mksrf_fvocef  -- Volatile Organic Compund Emission Factor dataset
    !    mksrf_fgdp ----- GDP dataset
    !    mksrf_fpeat ---- Peatland dataset
    !    mksrf_fsoildepth Soil depth dataset
    !    mksrf_fabm ----- Agricultural fire peak month dataset
    !    mksrf_ftopostats Topography statistics dataset
    !    mksrf_fvic ----- VIC parameters dataset
    !    mksrf_fch4 ----- inversion-derived CH4 parameters dataset
    ! ======================================
    ! Must specify mapping file for the different datafiles above
    ! ======================================
    !    map_fpft -------- Mapping for mksrf_fvegtyp
    !    map_flakwat ----- Mapping for mksrf_flakwat
    !    map_fwetlnd ----- Mapping for mksrf_fwetlnd
    !    map_fglacier ---- Mapping for mksrf_fglacier
    !    map_fglacierregion - Mapping for mksrf_fglacierregion
    !    map_fsoitex ----- Mapping for mksrf_fsoitex
    !    map_fsoicol ----- Mapping for mksrf_fsoicol
    !    map_furban ------ Mapping for mksrf_furban
    !    map_furbtopo ---- Mapping for mksrf_furbtopo
    !    map_fmax -------- Mapping for mksrf_fmax
    !    map_forganic ---- Mapping for mksrf_forganic
    !    map_fvocef ------ Mapping for mksrf_fvocef
    !    map_flai -------- Mapping for mksrf_flai
    !    map_fharvest ---- Mapping for mksrf_flai harvesting
    !    map_fgdp -------- Mapping for mksrf_fgdp
    !    map_fpeat ------- Mapping for mksrf_fpeat
    !    map_fsoildepth -- Mapping for mksrf_fsoildepth
    !    map_fabm -------- Mapping for mksrf_fabm
    !    map_ftopostats -- Mapping for mksrf_ftopostats
    !    map_fvic -------- Mapping for mksrf_fvic
    !    map_fch4 -------- Mapping for mksrf_fch4
    ! ======================================
    ! Optionally specify setting for:
    ! ======================================
    !    mksrf_fdynuse ----- ASCII text file that lists each year of pft files to use
    !    mksrf_gridtype ---- Type of grid (default is 'global')
    !    outnc_double ------ If output should be in double precision
    !    outnc_large_files - If output should be in NetCDF large file format
    !    nglcec ------------ If you want to change the number of Glacier elevation classes
    !    gitdescribe ------- Description of this version from git
    ! ======================================
    ! Optional settings to change values for entire area
    ! ======================================
    !    all_urban --------- If entire area is urban
    !    no_inlandwet ------ If wetland should be set to 0% over land
    !    soil_color -------- If you want to change the soil_color to this value everywhere
    !    soil_clay --------- If you want to change the soil_clay % to this value everywhere
    !    soil_fmax --------- If you want to change the soil_fmax  to this value everywhere
    !    soil_sand --------- If you want to change the soil_sand % to this value everywhere
    !    pft_idx ----------- If you want to change to 100% veg covered with given PFT indices
    !    pft_frc ----------- Fractions that correspond to the pft_idx above
    ! ==================
    !    numpft            (if different than default of 16)
    ! ======================================
    ! Optional settings to work around urban bug?
    ! ======================================
    !    urban_skip_abort_on_invalid_data_check
    ! ======================================================================

    write(6,*) 'Attempting to initialize control settings .....'

    mksrf_gridtype    = 'global'
    outnc_large_files = .false.
    outnc_double      = .true.
    all_urban         = .false.
    no_inlandwet      = .true.

    ! default value for bug work around
    urban_skip_abort_on_invalid_data_check = .false.

    read(5, clmexp, iostat=ier)
    if (ier /= 0) then
       write(6,*)'error: namelist input resulted in error code ',ier
       call abort()
    endif

    write (6,*) 'Attempting to create surface boundary data .....'
    write (6,'(72a1)') ("-",n=1,60)

    ! ----------------------------------------------------------------------
    ! Error check namelist input
    ! ----------------------------------------------------------------------
    
    if (urban_skip_abort_on_invalid_data_check) then
       write(6, *) "WARNING: aborting on invalid data check in urban has been disabled!"
       write(6, *) "WARNING: urban data may be invalid!"
    end if
    
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
    if ( no_inlandwet )then
       write(6,*) 'Set wetland to 0% over land'
    end if
    if (nglcec <= 0) then
       write(6,*) 'nglcec must be at least 1'
       call abort()
    end if

    !
    ! Call module initialization routines
    !
    call mksoilInit( )
    call mkpftInit( all_urban, all_veg )
    allocate ( elevclass(nglcec+1) )
    call mkglcmecInit (elevclass)
    call mkurbanInit (mksrf_furban)

    if ( all_veg )then
       write(6,*) 'Output ALL data in file as 100% vegetated'
    end if

    ! ----------------------------------------------------------------------
    ! Determine land model grid, fractional land and land mask
    ! ----------------------------------------------------------------------
    
    write(6,*)'calling domain_read'
    if ( .not. domain_read_map(ldomain, fgrddat) )then
        call domain_read(ldomain, fgrddat)
    end if
    write(6,*)'finished domain_read'
    
    ! Invalidate mask and frac for ldomain 

    !ldomain%mask = bigint
    !ldomain%frac = nan

    ! Determine if will have 1d output

    if (ldomain%ni /= -9999 .and. ldomain%nj /= -9999) then
       write(6,*)'fsurdat is 2d lat/lon grid'
       write(6,*)'nlon= ',ldomain%ni,' nlat= ',ldomain%nj
       if (outnc_dims == 1) then
          write(6,*)' writing output file in 1d gridcell format'
       end if
    else
       write(6,*)'fsurdat is 1d gridcell grid'
       outnc_dims = 1
    end if

    outnc_1d = .false.
    if ((ldomain%ni == -9999 .and. ldomain%nj == -9999) .or. outnc_dims==1) then
       outnc_1d = .true.
       write(6,*)'output file will be 1d'
    end if

    ! ----------------------------------------------------------------------
    ! Allocate and initialize dynamic memory
    ! ----------------------------------------------------------------------

    ns_o = ldomain%ns
    allocate ( landfrac_pft(ns_o)                 , &
               pctlnd_pft(ns_o)                   , & 
               pftdata_mask(ns_o)                 , & 
               pctnatpft(ns_o)                    , &
               pctnatpft_max(ns_o)                , &
               pctcft(ns_o)                       , &
               pctcft_max(ns_o)                   , &
               pctcft_saved(ns_o)                 , &
               pctgla(ns_o)                       , & 
               pctlak(ns_o)                       , & 
               pctwet(ns_o)                       , & 
               pcturb(ns_o)                       , &
               urban_region(ns_o)                 , &
               urbn_classes(ns_o,numurbl)         , &
               urbn_classes_g(ns_o,numurbl)       , &
               pctsand(ns_o,nlevsoi)              , & 
               pctclay(ns_o,nlevsoi)              , & 
               soicol(ns_o)                       , & 
               gdp(ns_o)                          , & 
               fpeat(ns_o)                        , & 
               soildepth(ns_o)                    , & 
               agfirepkmon(ns_o)                  , & 
               topo_stddev(ns_o)                  , &
               slope(ns_o)                        , &
               vic_binfl(ns_o)                    , &
               vic_ws(ns_o)                       , &
               vic_dsmax(ns_o)                    , &
               vic_ds(ns_o)                       , &
               lakedepth(ns_o)                    , &
               f0(ns_o)                           , &
               p3(ns_o)                           , &
               zwt0(ns_o)                         , &
               glacier_region(ns_o)               )       
    landfrac_pft(:)       = spval 
    pctlnd_pft(:)         = spval
    pftdata_mask(:)       = -999
    pctgla(:)             = spval
    pctlak(:)             = spval
    pctwet(:)             = spval
    pcturb(:)             = spval
    urban_region(:)       = -999
    urbn_classes(:,:)     = spval
    urbn_classes_g(:,:)   = spval
    pctsand(:,:)          = spval
    pctclay(:,:)          = spval
    soicol(:)             = -999
    gdp(:)                = spval
    fpeat(:)              = spval
    soildepth(:)          = spval
    agfirepkmon(:)        = -999
    topo_stddev(:)        = spval
    slope(:)              = spval
    vic_binfl(:)          = spval
    vic_ws(:)             = spval
    vic_dsmax(:)          = spval
    vic_ds(:)             = spval
    lakedepth(:)          = spval
    f0(:)                 = spval
    p3(:)                 = spval
    zwt0(:)               = spval
    glacier_region(:)     = -999

    ! ----------------------------------------------------------------------
    ! Open diagnostic output log file
    ! ----------------------------------------------------------------------
    
    if (fsurlog == ' ') then
       write(6,*)' must specify fsurlog in namelist'
       stop
    else
       ndiag = getavu(); call opnfil (fsurlog, ndiag, 'f')
    end if
    
    if (urban_skip_abort_on_invalid_data_check) then
       write(ndiag, *) "WARNING: aborting on invalid data check in urban has been disabled!"
       write(ndiag, *) "WARNING: urban data may be invalid!"
    end if
    
    if (mksrf_fgrid /= ' ')then
       write (ndiag,*)'using fractional land data from file= ', &
            trim(mksrf_fgrid),' to create the surface dataset'
    endif

    if (trim(mksrf_gridtype) == 'global' .or. &
        trim(mksrf_gridtype) == 'regional') then
       write(6,*)'mksrf_gridtype = ',trim(mksrf_gridtype)
    endif

    write(ndiag,*) 'PFTs from:                   ',trim(mksrf_fvegtyp)
    write(ndiag,*) 'harvest from:                ',trim(mksrf_fhrvtyp)
    write(ndiag,*) 'fmax from:                   ',trim(mksrf_fmax)
    write(ndiag,*) 'glaciers from:               ',trim(mksrf_fglacier)
    write(ndiag,*) '           with:             ', nglcec, ' glacier elevation classes'
    write(ndiag,*) 'glacier region ID from:      ',trim(mksrf_fglacierregion)
    write(ndiag,*) 'urban topography from:       ',trim(mksrf_furbtopo)
    write(ndiag,*) 'urban from:                  ',trim(mksrf_furban)
    write(ndiag,*) 'inland lake from:            ',trim(mksrf_flakwat)
    write(ndiag,*) 'inland wetland from:         ',trim(mksrf_fwetlnd)
    write(ndiag,*) 'soil texture from:           ',trim(mksrf_fsoitex)
    write(ndiag,*) 'soil organic from:           ',trim(mksrf_forganic)
    write(ndiag,*) 'soil color from:             ',trim(mksrf_fsoicol)
    write(ndiag,*) 'VOC emission factors from:   ',trim(mksrf_fvocef)
    write(ndiag,*) 'gdp from:                    ',trim(mksrf_fgdp)
    write(ndiag,*) 'peat from:                   ',trim(mksrf_fpeat)
    write(ndiag,*) 'soil depth from:             ',trim(mksrf_fsoildepth)
    write(ndiag,*) 'abm from:                    ',trim(mksrf_fabm)
    write(ndiag,*) 'topography statistics from:  ',trim(mksrf_ftopostats)
    write(ndiag,*) 'VIC parameters from:         ',trim(mksrf_fvic)
    write(ndiag,*) 'CH4 parameters from:         ',trim(mksrf_fch4)
    write(ndiag,*)' mapping for pft              ',trim(map_fpft)
    write(ndiag,*)' mapping for lake water       ',trim(map_flakwat)
    write(ndiag,*)' mapping for wetland          ',trim(map_fwetlnd)
    write(ndiag,*)' mapping for glacier          ',trim(map_fglacier)
    write(ndiag,*)' mapping for glacier region   ',trim(map_fglacierregion)
    write(ndiag,*)' mapping for soil texture     ',trim(map_fsoitex)
    write(ndiag,*)' mapping for soil color       ',trim(map_fsoicol)
    write(ndiag,*)' mapping for soil organic     ',trim(map_forganic)
    write(ndiag,*)' mapping for urban            ',trim(map_furban)
    write(ndiag,*)' mapping for fmax             ',trim(map_fmax)
    write(ndiag,*)' mapping for VOC pct emis     ',trim(map_fvocef)
    write(ndiag,*)' mapping for harvest          ',trim(map_fharvest)
    write(ndiag,*)' mapping for lai/sai          ',trim(map_flai)
    write(ndiag,*)' mapping for urb topography   ',trim(map_furbtopo)
    write(ndiag,*)' mapping for GDP              ',trim(map_fgdp)
    write(ndiag,*)' mapping for peatlands        ',trim(map_fpeat)
    write(ndiag,*)' mapping for soil depth       ',trim(map_fsoildepth)
    write(ndiag,*)' mapping for ag fire pk month ',trim(map_fabm)
    write(ndiag,*)' mapping for topography stats ',trim(map_ftopostats)
    write(ndiag,*)' mapping for VIC parameters   ',trim(map_fvic)
    write(ndiag,*)' mapping for CH4 parameters   ',trim(map_fch4)

    if (mksrf_fdynuse /= ' ') then
       write(6,*)'mksrf_fdynuse = ',trim(mksrf_fdynuse)
    end if

    ! ----------------------------------------------------------------------
    ! Make surface dataset fields
    ! ----------------------------------------------------------------------

    ! Make PFTs [pctnatpft, pctcft] from dataset [fvegtyp]

    call mkpft(ldomain, mapfname=map_fpft, fpft=mksrf_fvegtyp, &
         ndiag=ndiag, allow_no_crops=.false., &
         pctlnd_o=pctlnd_pft, pctnatpft_o=pctnatpft, pctcft_o=pctcft)

    ! Create harvesting data at model resolution
    call mkharvest_init( ns_o, spval, harvdata, mksrf_fhrvtyp )
    if ( .not. any(pft_frc > 0.0_r8 ) )then

       call mkharvest( ldomain, mapfname=map_fharvest, datfname=mksrf_fhrvtyp, &
                       ndiag=ndiag, harvdata=harvdata )
    end if

    ! Save the version of pctcft before any corrections are made. In particular, we want
    ! to save the version before remove_small_cover is called.
    pctcft_saved = pctcft

    ! Make inland water [pctlak, pctwet] [flakwat] [fwetlnd]

    call mklakwat (ldomain, mapfname=map_flakwat, datfname=mksrf_flakwat, &
         ndiag=ndiag, zero_out=all_urban.or.all_veg, lake_o=pctlak)

    call mkwetlnd (ldomain, mapfname=map_fwetlnd, datfname=mksrf_fwetlnd, &
         ndiag=ndiag, zero_out=all_urban.or.all_veg.or.no_inlandwet, swmp_o=pctwet)

    ! Make glacier fraction [pctgla] from [fglacier] dataset

    call mkglacier (ldomain, mapfname=map_fglacier, datfname=mksrf_fglacier, &
         ndiag=ndiag, zero_out=all_urban.or.all_veg, glac_o=pctgla)

    ! Make glacier region ID [glacier_region] from [fglacierregion] dataset

    call mkglacierregion (ldomain, mapfname=map_fglacierregion, &
         datfname=mksrf_fglacierregion, ndiag=ndiag, &
         glacier_region_o = glacier_region)

    ! Make soil texture [pctsand, pctclay]  [fsoitex]

    call mksoiltex (ldomain, mapfname=map_fsoitex, datfname=mksrf_fsoitex, &
         ndiag=ndiag, sand_o=pctsand, clay_o=pctclay)
    ! Make soil color classes [soicol] [fsoicol]

    call mksoilcol (ldomain, mapfname=map_fsoicol, datfname=mksrf_fsoicol, &
         ndiag=ndiag, soil_color_o=soicol, nsoicol=nsoicol)

    ! Make fmax [fmax] from [fmax] dataset

    allocate(fmax(ns_o))
    fmax(:) = spval
    call mkfmax (ldomain, mapfname=map_fmax, datfname=mksrf_fmax, &
         ndiag=ndiag, fmax_o=fmax)

    ! Make GDP data [gdp] from [gdp]

    call mkgdp (ldomain, mapfname=map_fgdp, datfname=mksrf_fgdp, &
         ndiag=ndiag, gdp_o=gdp)

    ! Make peat data [fpeat] from [peatf]

    call mkpeat (ldomain, mapfname=map_fpeat, datfname=mksrf_fpeat, &
         ndiag=ndiag, peat_o=fpeat)

    ! Make soil depth data [soildepth] from [soildepthf]

    call mksoildepth (ldomain, mapfname=map_fsoildepth, datfname=mksrf_fsoildepth, &
         ndiag=ndiag, soildepth_o=soildepth)

    ! Make agricultural fire peak month data [abm] from [abm]

    call mkagfirepkmon (ldomain, mapfname=map_fabm, datfname=mksrf_fabm, &
         ndiag=ndiag, agfirepkmon_o=agfirepkmon)

    ! Make urban fraction [pcturb] from [furban] dataset

    call mkurban (ldomain, mapfname=map_furban, datfname=mksrf_furban, &
         ndiag=ndiag, zero_out=all_veg, urbn_o=pcturb, urbn_classes_o=urbn_classes, &
         region_o=urban_region)

    ! Make elevation [elev] from [ftopo, ffrac] dataset
    ! Used only to screen pcturb  
    ! Screen pcturb by elevation threshold from elev dataset

    if ( .not. all_urban .and. .not. all_veg )then
       allocate(elev(ns_o))
       elev(:) = spval
       ! NOTE(wjs, 2016-01-15) This uses the 'TOPO_ICE' variable for historical reasons
       ! (this same dataset used to be used for glacier-related purposes as well).
       ! TODO(wjs, 2016-01-15) A better solution for this urban screening would probably
       ! be to modify the raw urban data; in that case, I believe we could remove
       ! furbtopo.
       call mkelev (ldomain, mapfname=map_furbtopo, datfname=mksrf_furbtopo, &
         varname='TOPO_ICE', ndiag=ndiag, elev_o=elev)

       where (elev .gt. elev_thresh)
         pcturb = 0._r8
       end where
       deallocate(elev)
    end if
    
    ! Compute topography statistics [topo_stddev, slope] from [ftopostats]
    call mktopostats (ldomain, mapfname=map_ftopostats, datfname=mksrf_ftopostats, &
         ndiag=ndiag, topo_stddev_o=topo_stddev, slope_o=slope, std_elev=std_elev)

    ! Make VIC parameters [binfl, ws, dsmax, ds] from [fvic]
    call mkVICparams (ldomain, mapfname=map_fvic, datfname=mksrf_fvic, ndiag=ndiag, &
         binfl_o=vic_binfl, ws_o=vic_ws, dsmax_o=vic_dsmax, ds_o=vic_ds)

    ! Make lake depth [lakedepth] from [flakwat]
    call mklakparams (ldomain, mapfname=map_flakwat, datfname=mksrf_flakwat, ndiag=ndiag, &
         lakedepth_o=lakedepth)

    ! Make inversion-derived CH4 parameters [f0, p3, zwt0] from [fch4]
    call mkCH4inversion (ldomain, mapfname=map_fch4, datfname=mksrf_fch4, ndiag=ndiag, &
         f0_o=f0, p3_o=p3, zwt0_o=zwt0)

    ! Make organic matter density [organic] [forganic]
    allocate (organic(ns_o,nlevsoi))
    organic(:,:) = spval
    call mkorganic (ldomain, mapfname=map_forganic, datfname=mksrf_forganic, &
         ndiag=ndiag, organic_o=organic)

    ! Make VOC emission factors for isoprene &
    ! [ef1_btr,ef1_fet,ef1_fdt,ef1_shr,ef1_grs,ef1_crp]

    allocate ( ef1_btr(ns_o) , & 
               ef1_fet(ns_o) , & 
               ef1_fdt(ns_o) , & 
               ef1_shr(ns_o) , & 
               ef1_grs(ns_o) , & 
               ef1_crp(ns_o) )
    ef1_btr(:) = 0._r8
    ef1_fet(:) = 0._r8
    ef1_fdt(:) = 0._r8
    ef1_shr(:) = 0._r8
    ef1_grs(:) = 0._r8
    ef1_crp(:) = 0._r8
    
    call mkvocef (ldomain, mapfname=map_fvocef, datfname=mksrf_fvocef, ndiag=ndiag, &
         ef_btr_o=ef1_btr, ef_fet_o=ef1_fet, ef_fdt_o=ef1_fdt,  &
         ef_shr_o=ef1_shr, ef_grs_o=ef1_grs, ef_crp_o=ef1_crp)

    ! Do landuse changes such as for the poles, etc.

    call change_landuse( ldomain, dynpft=.false. )

    do n = 1,ns_o

       ! Assume wetland and/or lake when dataset landmask implies ocean 
       ! (assume medium soil color (15) and loamy texture).
       ! Also set pftdata_mask here

       if (pctlnd_pft(n) < 1.e-6_r8) then
          pftdata_mask(n)  = 0
          soicol(n)        = 15
          pctwet(n)        = 100._r8 - pctlak(n)
          pcturb(n)        = 0._r8
          pctgla(n)        = 0._r8
          call pctnatpft(n)%set_pct_l2g(0._r8)
          call pctcft(n)%set_pct_l2g(0._r8)
          pctsand(n,:)     = 43._r8
          pctclay(n,:)     = 18._r8
          organic(n,:)   = 0._r8
       else
          pftdata_mask(n) = 1
       end if

       ! Truncate all percentage fields on output grid. This is needed to
       ! insure that wt is zero (not a very small number such as
       ! 1e-16) where it really should be zero
       
       do k = 1,nlevsoi
          pctsand(n,k) = float(nint(pctsand(n,k)))
          pctclay(n,k) = float(nint(pctclay(n,k)))
       end do
       pctlak(n) = float(nint(pctlak(n)))
       pctwet(n) = float(nint(pctwet(n)))
       pctgla(n) = float(nint(pctgla(n)))
       
       ! Make sure sum of land cover types does not exceed 100. If it does,
       ! subtract excess from most dominant land cover.
       
       suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n)
       if (suma > 250._r4) then
          write (6,*) subname, ' error: sum of pctlak, pctwet,', &
               'pcturb and pctgla is greater than 250%'
          write (6,*)'n,pctlak,pctwet,pcturb,pctgla= ', &
               n,pctlak(n),pctwet(n),pcturb(n),pctgla(n)
          call abort()
       else if (suma > 100._r4) then
          pctlak(n) = pctlak(n) * 100._r8/suma
          pctwet(n) = pctwet(n) * 100._r8/suma
          pcturb(n) = pcturb(n) * 100._r8/suma
          pctgla(n) = pctgla(n) * 100._r8/suma
       end if
       
    end do

    call normalizencheck_landuse(ldomain)

    ! Write out sum of PFT's

    do k = natpft_lb,natpft_ub
       suma = 0._r8
       do n = 1,ns_o
          suma = suma + pctnatpft(n)%get_one_pct_p2g(k)
       enddo
       write(6,*) 'sum over domain of pft ',k,suma
    enddo
    write(6,*)

    do k = cft_lb,cft_ub
       suma = 0._r8
       do n = 1,ns_o
          suma = suma + pctcft(n)%get_one_pct_p2g(k)
       enddo
       write(6,*) 'sum over domain of cft ',k,suma
    enddo
    write(6,*)

    ! Make final values of percent urban by class
    ! This call needs to occur after all corrections are made to pcturb

    call normalize_classes_by_gcell(urbn_classes, pcturb, urbn_classes_g)


    ! Make glacier multiple elevation classes [pctglcmec,topoglcmec] from [fglacier] dataset
    ! This call needs to occur after pctgla has been adjusted for the final time

    allocate (pctglcmec(ns_o,nglcec), &
         topoglcmec(ns_o,nglcec), &
         pctglcmec_gic(ns_o,nglcec), &
         pctglcmec_icesheet(ns_o,nglcec))
    allocate (pctglc_gic(ns_o))
    allocate (pctglc_icesheet(ns_o))

    pctglcmec(:,:)          = spval
    topoglcmec(:,:)         = spval
    pctglcmec_gic(:,:)      = spval
    pctglcmec_icesheet(:,:) = spval
    pctglc_gic(:)           = spval
    pctglc_icesheet(:)      = spval

    call mkglcmec (ldomain, mapfname=map_fglacier, &
         datfname_fglacier=mksrf_fglacier, ndiag=ndiag, &
         pctglcmec_o=pctglcmec, topoglcmec_o=topoglcmec, &
         pctglcmec_gic_o=pctglcmec_gic, pctglcmec_icesheet_o=pctglcmec_icesheet, &
         pctglc_gic_o=pctglc_gic, pctglc_icesheet_o=pctglc_icesheet)

    ! Determine fractional land from pft dataset

    do n = 1,ns_o
       landfrac_pft(n) = pctlnd_pft(n)/100._r8
    end do

    ! ----------------------------------------------------------------------
    ! Create surface dataset
    ! ----------------------------------------------------------------------

    ! Create netCDF surface dataset.  

    ! If fsurdat is blank, then we do not write a surface dataset - but we may still
    ! write a dynamic landuse file. This is useful if we are creating many datasets at
    ! once, and don't want duplicate surface datasets.
    !
    ! TODO(wjs, 2016-01-26) Ideally, we would also avoid doing the processing of
    ! variables that are just needed by the surface dataset (not by the dynamic landuse
    ! file). However, this would require some analysis of the above code, to determine
    ! which processing is needed (directly or indirectly) in order to create a dynamic
    ! landuse file.

    if (fsurdat /= ' ') then

       call mkfile(ldomain, trim(fsurdat), harvdata, dynlanduse = .false.)

       call domain_write(ldomain, fsurdat)

       call check_ret(nf_open(trim(fsurdat), nf_write, ncid), subname)
       call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

       ! Write fields OTHER THAN lai, sai, heights, and urban parameters to netcdf surface dataset

       call check_ret(nf_inq_varid(ncid, 'natpft', varid), subname)
       call check_ret(nf_put_var_int(ncid, varid, (/(n,n=natpft_lb,natpft_ub)/)), subname)

       if (num_cft > 0) then
          call check_ret(nf_inq_varid(ncid, 'cft', varid), subname)
          call check_ret(nf_put_var_int(ncid, varid, (/(n,n=cft_lb,cft_ub)/)), subname)
       end if

       call check_ret(nf_inq_varid(ncid, 'PFTDATA_MASK', varid), subname)
       call check_ret(nf_put_var_int(ncid, varid, pftdata_mask), subname)

       call check_ret(nf_inq_varid(ncid, 'LANDFRAC_PFT', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, landfrac_pft), subname)

       call check_ret(nf_inq_varid(ncid, 'mxsoil_color', varid), subname)
       call check_ret(nf_put_var_int(ncid, varid, nsoicol), subname)

       call check_ret(nf_inq_varid(ncid, 'SOIL_COLOR', varid), subname)
       call check_ret(nf_put_var_int(ncid, varid, soicol), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_SAND', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, pctsand), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_CLAY', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, pctclay), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_WETLAND', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, pctwet), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_LAKE', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, pctlak), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_GLACIER', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, pctgla), subname)

       call check_ret(nf_inq_varid(ncid, 'GLACIER_REGION', varid), subname)
       call check_ret(nf_put_var_int(ncid, varid, glacier_region), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_GLC_MEC', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, pctglcmec), subname)

       call check_ret(nf_inq_varid(ncid, 'GLC_MEC', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, elevclass), subname)

       call check_ret(nf_inq_varid(ncid, 'TOPO_GLC_MEC', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, topoglcmec), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_GLC_MEC_GIC', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, pctglcmec_gic), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_GLC_MEC_ICESHEET', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, pctglcmec_icesheet), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_GLC_GIC', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, pctglc_gic), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_GLC_ICESHEET', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, pctglc_icesheet), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_URBAN', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, urbn_classes_g), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_NATVEG', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, get_pct_l2g_array(pctnatpft)), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_CROP', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, get_pct_l2g_array(pctcft)), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_NAT_PFT', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, get_pct_p2l_array(pctnatpft)), subname)

       if (num_cft > 0) then
          call check_ret(nf_inq_varid(ncid, 'PCT_CFT', varid), subname)
          call check_ret(nf_put_var_double(ncid, varid, get_pct_p2l_array(pctcft)), subname)
       end if

       call harvdata%getFieldsIdx( harvind1D, harvind2D )
       do k = 1, harvdata%num1Dfields()
          call check_ret(nf_inq_varid(ncid, trim(mkharvest_fieldname(harvind1D(k),constant=.true.)), varid), subname)
          harvest1D => harvdata%get1DFieldPtr( harvind1D(k), output=.true. )
          call check_ret(nf_put_var_double(ncid, varid, harvest1D), subname)
       end do
       do k = 1, harvdata%num2Dfields()
          call check_ret(nf_inq_varid(ncid, trim(mkharvest_fieldname(harvind2D(k),constant=.true.)), varid), subname)
          harvest2D => harvdata%get2DFieldPtr( harvind2D(k), output=.true. )
          call check_ret(nf_put_var_double(ncid, varid, harvest2D), subname)
       end do
       deallocate( harvind1D, harvind2D )

       call check_ret(nf_inq_varid(ncid, 'FMAX', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, fmax), subname)

       call check_ret(nf_inq_varid(ncid, 'gdp', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, gdp), subname)

       call check_ret(nf_inq_varid(ncid, 'peatf', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, fpeat), subname)


       !    call check_ret(nf_inq_varid(ncid, 'Avg_Depth_Median', varid), subname)
       call check_ret(nf_inq_varid(ncid, 'zbedrock', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, soildepth), subname)

       call check_ret(nf_inq_varid(ncid, 'abm', varid), subname)
       call check_ret(nf_put_var_int(ncid, varid, agfirepkmon), subname)

       call check_ret(nf_inq_varid(ncid, 'SLOPE', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, slope), subname)

       call check_ret(nf_inq_varid(ncid, 'STD_ELEV', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, topo_stddev), subname)

       call check_ret(nf_inq_varid(ncid, 'binfl', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, vic_binfl), subname)

       call check_ret(nf_inq_varid(ncid, 'Ws', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, vic_ws), subname)

       call check_ret(nf_inq_varid(ncid, 'Dsmax', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, vic_dsmax), subname)

       call check_ret(nf_inq_varid(ncid, 'Ds', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, vic_ds), subname)

       call check_ret(nf_inq_varid(ncid, 'LAKEDEPTH', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, lakedepth), subname)

       call check_ret(nf_inq_varid(ncid, 'F0', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, f0), subname)

       call check_ret(nf_inq_varid(ncid, 'P3', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, p3), subname)

       call check_ret(nf_inq_varid(ncid, 'ZWT0', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, zwt0), subname)

       call check_ret(nf_inq_varid(ncid, 'EF1_BTR', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, ef1_btr), subname)

       call check_ret(nf_inq_varid(ncid, 'EF1_FET', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, ef1_fet), subname)

       call check_ret(nf_inq_varid(ncid, 'EF1_FDT', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, ef1_fdt), subname)

       call check_ret(nf_inq_varid(ncid, 'EF1_SHR', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, ef1_shr), subname)

       call check_ret(nf_inq_varid(ncid, 'EF1_GRS', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, ef1_grs), subname)

       call check_ret(nf_inq_varid(ncid, 'EF1_CRP', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, ef1_crp), subname)

       call check_ret(nf_inq_varid(ncid, 'ORGANIC', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, organic), subname)

       call check_ret(nf_inq_varid(ncid, 'URBAN_REGION_ID', varid), subname)
       call check_ret(nf_put_var_int(ncid, varid, urban_region), subname)

       ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

       call check_ret(nf_sync(ncid), subname)

       ! ----------------------------------------------------------------------
       ! Make Urban Parameters from raw input data and write to surface dataset 
       ! Write to netcdf file is done inside mkurbanpar routine
       ! ----------------------------------------------------------------------

       write(6,*)'calling mkurbanpar'
       call mkurbanpar(datfname=mksrf_furban, ncido=ncid, region_o=urban_region, &
            urbn_classes_gcell_o=urbn_classes_g, &
            urban_skip_abort_on_invalid_data_check=urban_skip_abort_on_invalid_data_check)

       ! ----------------------------------------------------------------------
       ! Make LAI and SAI from 1/2 degree data and write to surface dataset 
       ! Write to netcdf file is done inside mklai routine
       ! ----------------------------------------------------------------------

       write(6,*)'calling mklai'
       call mklai(ldomain, mapfname=map_flai, datfname=mksrf_flai, &
            ndiag=ndiag, ncido=ncid )

       ! Close surface dataset

       call check_ret(nf_close(ncid), subname)

       write (6,'(72a1)') ("-",n=1,60)
       write (6,*)' land model surface data set successfully created for ', &
            'grid of size ',ns_o

    else  ! fsurdat == ' '

       write (6,*) 'fsurdat is blank: skipping writing surface dataset'

    end if  ! if (fsurdat /= ' ')

    ! Deallocate arrays NOT needed for dynamic-pft section of code

    deallocate ( organic )
    deallocate ( ef1_btr, ef1_fet, ef1_fdt, ef1_shr, ef1_grs, ef1_crp )
    deallocate ( pctglcmec, topoglcmec)
    deallocate ( pctglc_gic, pctglc_icesheet)
    deallocate ( elevclass )
    deallocate ( fmax )
    deallocate ( pctsand, pctclay )
    deallocate ( soicol )
    deallocate ( gdp, fpeat, agfirepkmon )
    deallocate ( soildepth )
    deallocate ( topo_stddev, slope )
    deallocate ( vic_binfl, vic_ws, vic_dsmax, vic_ds )
    deallocate ( lakedepth )
    deallocate ( f0, p3, zwt0 )
    deallocate ( glacier_region )

    call harvdata%clean()

    ! ----------------------------------------------------------------------
    ! Create dynamic land use dataset if appropriate
    ! ----------------------------------------------------------------------

    if (mksrf_fdynuse /= ' ') then

       write(6,*)'creating dynamic land use dataset'

       allocate(pctlnd_pft_dyn(ns_o))
       call mkharvest_init( ns_o, spval, harvdata, mksrf_fhrvtyp )

       if (fdyndat == ' ') then
          write(6,*)' must specify fdyndat in namelist if mksrf_fdynuse is not blank'
          stop
       end if

       ! Define dimensions and global attributes

       call mkfile(ldomain, fdyndat, harvdata, dynlanduse=.true.)

       ! Write fields other pft to dynamic land use dataset

       call domain_write(ldomain, fdyndat)

       call check_ret(nf_open(trim(fdyndat), nf_write, ncid), subname)
       call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

       call check_ret(nf_inq_varid(ncid, 'natpft', varid), subname)
       call check_ret(nf_put_var_int(ncid, varid, (/(n,n=natpft_lb,natpft_ub)/)), subname)

       if (num_cft > 0) then
          call check_ret(nf_inq_varid(ncid, 'cft', varid), subname)
          call check_ret(nf_put_var_int(ncid, varid, (/(n,n=cft_lb,cft_ub)/)), subname)
       end if

       call check_ret(nf_inq_varid(ncid, 'PFTDATA_MASK', varid), subname)
       call check_ret(nf_put_var_int(ncid, varid, pftdata_mask), subname)

       call check_ret(nf_inq_varid(ncid, 'LANDFRAC_PFT', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, landfrac_pft), subname)

       ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

       call check_ret(nf_sync(ncid), subname)

       ! Read in each dynamic pft landuse dataset

       nfdyn = getavu(); call opnfil (mksrf_fdynuse, nfdyn, 'f')

       pctnatpft_max = pctnatpft
       pctcft_max = pctcft

       ntim = 0
       do 
          ! Read input pft data

          read(nfdyn, '(A195,1x,I4)', iostat=ier) string, year
          if (ier /= 0) exit
          !
          ! If pft fraction override is set, than intrepret string as PFT and harvesting override values
          !
          if ( any(pft_frc > 0.0_r8 ) )then
             fname = ' '
             fhrvname  = ' '
             call mkpft_parse_oride(string)
             call mkharvest_parse_oride(string)
	     write(6, '(a, i4, a)') 'PFT and harvesting values for year ', year, ' :'
             write(6, '(a, a)') '    ', trim(string)
          !
          ! Otherwise intrepret string as a filename with PFT and harvesting values in it
          !
          else
             fname = string
	     write(6,*)'input pft dynamic dataset for year ', year, ' is : ', trim(fname)
             read(nfdyn, '(A195,1x,I4)', iostat=ier) fhrvname, year2
             if ( year2 /= year ) then
                write(6,*) subname, ' error: year for harvest not equal to year for PFT files'
                call abort()
             end if
          end if
          ntim = ntim + 1

          ! Create pctpft data at model resolution
          
          call mkpft(ldomain, mapfname=map_fpft, fpft=fname, &
               ndiag=ndiag, allow_no_crops=.false., &
               pctlnd_o=pctlnd_pft_dyn, pctnatpft_o=pctnatpft, pctcft_o=pctcft, &
               pctcft_o_saved=pctcft_saved)

          ! Create harvesting data at model resolution

          call mkharvest( ldomain, mapfname=map_fharvest, datfname=fhrvname, &
               ndiag=ndiag, harvdata=harvdata )

          ! Consistency check on input land fraction

          do n = 1,ns_o
             if (pctlnd_pft_dyn(n) /= pctlnd_pft(n)) then
                write(6,*) subname,' error: pctlnd_pft for dynamics data = ',&
                     pctlnd_pft_dyn(n), ' not equal to pctlnd_pft for surface data = ',&
                     pctlnd_pft(n),' at n= ',n
                if ( trim(fname) == ' ' )then
                   write(6,*) ' PFT string = ', string
                else
                   write(6,*) ' PFT file = ', fname
                end if
                call abort()
             end if
          end do

          call change_landuse(ldomain, dynpft=.true.)

          call normalizencheck_landuse(ldomain)
	  
          call update_max_array(pctnatpft_max,pctnatpft)
          call update_max_array(pctcft_max,pctcft)

          ! Output time-varying data for current year

          call check_ret(nf_inq_varid(ncid, 'PCT_NAT_PFT', varid), subname)
          call ncd_put_time_slice(ncid, varid, ntim, get_pct_p2l_array(pctnatpft))

          call check_ret(nf_inq_varid(ncid, 'PCT_CROP', varid), subname)
          call ncd_put_time_slice(ncid, varid, ntim, get_pct_l2g_array(pctcft))

          if (num_cft > 0) then
             call check_ret(nf_inq_varid(ncid, 'PCT_CFT', varid), subname)
             call ncd_put_time_slice(ncid, varid, ntim, get_pct_p2l_array(pctcft))
          end if

          call harvdata%getFieldsIdx( harvind1D, harvind2D )
          do k = 1, harvdata%num1Dfields()
             call check_ret(nf_inq_varid(ncid, trim(mkharvest_fieldname(harvind1D(k),constant=.false.)), varid), subname)
             harvest1D => harvdata%get1DFieldPtr( harvind1D(k), output=.true. )
             call ncd_put_time_slice(ncid, varid, ntim, harvest1D)
          end do
          do k = 1, harvdata%num2Dfields()
             call check_ret(nf_inq_varid(ncid, trim(mkharvest_fieldname(harvind2D(k),constant=.false.)), varid), subname)
             harvest2D => harvdata%get2DFieldPtr( harvind2D(k), output=.true. )
             call ncd_put_time_slice(ncid, varid, ntim, harvest2D)
          end do
          deallocate( harvind1D, harvind2D )

          call check_ret(nf_inq_varid(ncid, 'YEAR', varid), subname)
          call check_ret(nf_put_vara_int(ncid, varid, ntim, 1, year), subname)

          call check_ret(nf_inq_varid(ncid, 'time', varid), subname)
          call check_ret(nf_put_vara_int(ncid, varid, ntim, 1, year), subname)

          call check_ret(nf_inq_varid(ncid, 'input_pftdata_filename', varid), subname)
          call check_ret(nf_put_vara_text(ncid, varid, (/ 1, ntim /), (/ len_trim(string), 1 /), trim(string) ), subname)

	  ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

	  call check_ret(nf_sync(ncid), subname)

       end do   ! end of read loop

       call check_ret(nf_inq_varid(ncid, 'PCT_NAT_PFT_MAX', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, get_pct_p2l_array(pctnatpft_max)), subname)

       call check_ret(nf_inq_varid(ncid, 'PCT_CROP_MAX', varid), subname)
       call check_ret(nf_put_var_double(ncid, varid, get_pct_l2g_array(pctcft_max)), subname)

       if (num_cft > 0) then
          call check_ret(nf_inq_varid(ncid, 'PCT_CFT_MAX', varid), subname)
          call check_ret(nf_put_var_double(ncid, varid, get_pct_p2l_array(pctcft_max)), subname)
       end if

       call check_ret(nf_close(ncid), subname)

    end if   ! end of if-create dynamic landust dataset   

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

    write (6,*) 'Successfully created surface dataset'

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: change_landuse
!
! !INTERFACE:
subroutine change_landuse( ldomain, dynpft )
!
! !DESCRIPTION:
!
! Do landuse changes such as for the poles, etc.
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    type(domain_type)   :: ldomain
    logical, intent(in)  :: dynpft   ! if part of the dynpft section of code

!
! !REVISION HISTORY:
! 9/10/09: Erik Kluzek spin off subroutine from original embedded code
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: n,ns_o                       ! indices
    character(len=32) :: subname = 'change_landuse'  ! subroutine name
!-----------------------------------------------------------------------

    ns_o = ldomain%ns
    do n = 1,ns_o

       ! If have pole points on grid - set south pole to glacier
       ! north pole is assumed as non-land
       
       if (abs((ldomain%latc(n) - 90._r8)) < 1.e-6_r8) then
          pctlak(n)   = 0._r8
          pctwet(n)   = 0._r8
          pcturb(n)   = 0._r8
          pctgla(n)   = 100._r8
          call pctnatpft(n)%set_pct_l2g(0._r8)
          call pctcft(n)%set_pct_l2g(0._r8)
          if ( .not. dynpft )then
             organic(n,:)   = 0._r8
             ef1_btr(n)     = 0._r8
             ef1_fet(n)     = 0._r8
             ef1_fdt(n)     = 0._r8
             ef1_shr(n)     = 0._r8
             ef1_grs(n)     = 0._r8
             ef1_crp(n)     = 0._r8
          end if
       end if

    end do

end subroutine change_landuse

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: normalizencheck_landuse
!
! !INTERFACE:
subroutine normalizencheck_landuse(ldomain)
!
! !DESCRIPTION:
!
! Normalize land use and make sure things add up to 100% as well as
! checking that things are as they should be.
!
! Precondition: pctlak + pctwet + pcturb + pctgla <= 100 (within roundoff)
!
! !USES:
    use mkpftConstantsMod , only : baregroundindex
    use mkpftUtilsMod     , only : adjust_total_veg_area
    implicit none
! !ARGUMENTS:
    type(domain_type)   :: ldomain
!
! !REVISION HISTORY:
! 9/10/09: Erik Kluzek spin off subroutine from original embedded code
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m,k,n,ns_o                  ! indices
    integer  :: nsmall                      ! number of small PFT values for a single check
    integer  :: nsmall_tot                  ! total number of small PFT values in all grid cells
    real(r8) :: suma                        ! sum for error check
    real(r8) :: suma2                       ! another sum for error check
    real(r8) :: new_total_veg_pct           ! new % veg (% of grid cell, total of natural veg & crop)
    real(r8) :: bare_pct_p2g                ! % of bare soil, as % of grid cell
    real(r8) :: bare_urb_diff               ! difference between bare soil and urban %
    real(r8) :: pcturb_excess               ! excess urban % not accounted for by bare soil
    real(r8) :: sum8, sum8a                 ! sum for error check
    real(r4) :: sum4a                       ! sum for error check
    real(r8), parameter :: tol_loose = 1.e-4_r8               ! tolerance for some 'loose' error checks
    real(r8), parameter :: toosmallPFT = 1.e-10_r8            ! tolerance for PFT's to ignore
    character(len=32) :: subname = 'normalizencheck_landuse'  ! subroutine name
!-----------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! Normalize vegetated area so that vegetated + special area is 100%
    ! ------------------------------------------------------------------------

    ns_o = ldomain%ns
    do n = 1,ns_o

       ! Check preconditions

       suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n)
       if (suma > (100._r8 + tol_loose)) then
          write(6,*) subname, ' ERROR: pctlak + pctwet + pcturb + pctgla must be'
          write(6,*) '<= 100% before calling this subroutine'
          write(6,*) 'n, pctlak, pctwet, pcturb, pctgla = ', &
               n, pctlak(n), pctwet(n), pcturb(n), pctgla(n)
          call abort()
       end if

       ! First normalize vegetated (natural veg + crop) cover so that the total of
       ! (vegetated + (special excluding urban)) is 100%. We'll deal with urban later.
       !
       ! Note that, in practice, the total area of natural veg + crop is typically 100%
       ! going into this routine. However, the following code does NOT rely on this, and
       ! will work properly regardless of the initial area of natural veg + crop (even if
       ! that initial area is 0%).
       
       suma = pctlak(n)+pctwet(n)+pctgla(n)
       new_total_veg_pct = 100._r8 - suma
       ! correct for rounding error:
       new_total_veg_pct = max(new_total_veg_pct, 0._r8)

       call adjust_total_veg_area(new_total_veg_pct, pctnatpft=pctnatpft(n), pctcft=pctcft(n))

       ! Make sure we did the above rescaling correctly

       suma = suma + pctnatpft(n)%get_pct_l2g() + pctcft(n)%get_pct_l2g()
       if (abs(suma - 100._r8) > tol_loose) then
          write(6,*) subname, ' ERROR in rescaling veg based on (special excluding urban'
          write(6,*) 'suma = ', suma
          call abort()
       end if

       ! Now decrease the vegetated area to account for urban area. Urban needs to be
       ! handled specially because we replace bare soil preferentially with urban, rather
       ! than rescaling all PFTs equally.

       if (pcturb(n) > 0._r8) then
          
          ! Replace bare soil preferentially with urban
          bare_pct_p2g = pctnatpft(n)%get_one_pct_p2g(baregroundindex)
          bare_urb_diff = bare_pct_p2g - pcturb(n)
          bare_pct_p2g = max(0._r8, bare_urb_diff)
          call pctnatpft(n)%set_one_pct_p2g(baregroundindex, bare_pct_p2g)
          pcturb_excess = abs(min(0._r8,bare_urb_diff))
          
          ! For any urban not accounted for by bare soil, replace other PFTs
          ! proportionally
          if (pcturb_excess > 0._r8) then
             ! Note that, in this case, we will have already reduced bare ground to 0%

             new_total_veg_pct = pctnatpft(n)%get_pct_l2g() + pctcft(n)%get_pct_l2g() - pcturb_excess
             if (new_total_veg_pct < 0._r8) then
                if (abs(new_total_veg_pct) < tol_loose) then
                   ! only slightly less than 0; correct it
                   new_total_veg_pct = 0._r8
                else
                   write(6,*) subname, ' ERROR: trying to replace veg with urban,'
                   write(6,*) 'but pcturb_excess exceeds current vegetation percent'
                   call abort()
                end if
             end if

             call adjust_total_veg_area(new_total_veg_pct, pctnatpft=pctnatpft(n), pctcft=pctcft(n))
          end if

       end if ! pcturb(n) > 0

       ! Confirm that we have done the rescaling correctly: now the sum of all landunits
       ! should be 100%
       suma = pctlak(n)+pctwet(n)+pctgla(n)+pcturb(n)
       suma = suma + pctnatpft(n)%get_pct_l2g() + pctcft(n)%get_pct_l2g()
       if (abs(suma - 100._r8) > tol_loose) then
          write(6,*) subname, ' ERROR: landunits do not sum to 100%'
          write(6,*) 'n, suma, pctlak, pctwet, pctgla, pcturb, pctnatveg, pctcrop = '
          write(6,*) n, suma, pctlak(n), pctwet(n), pctgla(n), pcturb(n), &
               pctnatpft(n)%get_pct_l2g(), pctcft(n)%get_pct_l2g()
          call abort()
       end if
          
    end do

    ! ------------------------------------------------------------------------
    ! Do other corrections and error checks
    ! ------------------------------------------------------------------------

    nsmall_tot = 0

    do n = 1,ns_o
       
       ! If the coverage of any PFT or CFT is too small at the gridcell level, set its
       ! % cover to 0, then renormalize everything else as needed
       call pctnatpft(n)%remove_small_cover(toosmallPFT, nsmall)
       nsmall_tot = nsmall_tot + nsmall
       call pctcft(n)%remove_small_cover(toosmallPFT, nsmall)
       nsmall_tot = nsmall_tot + nsmall

       suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n)
       suma = suma + pctnatpft(n)%get_pct_l2g() + pctcft(n)%get_pct_l2g()
       if ( abs(suma - 100.0_r8) > 2.0*epsilon(suma) )then
          pctlak(n)    = pctlak(n)    * 100._r8/suma
          pctwet(n)    = pctwet(n)    * 100._r8/suma
          pcturb(n)    = pcturb(n)    * 100._r8/suma
          pctgla(n)    = pctgla(n)    * 100._r8/suma
          call pctnatpft(n)%set_pct_l2g(pctnatpft(n)%get_pct_l2g() * 100._r8/suma)
          call pctcft(n)%set_pct_l2g(pctcft(n)%get_pct_l2g() * 100._r8/suma)
       end if
       
       ! Roundoff error fix
       suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n)
       suma2 = pctnatpft(n)%get_pct_l2g() + pctcft(n)%get_pct_l2g()
       if ( (suma < 100._r8 .and. suma > (100._r8 - 1.e-6_r8)) .or. &
            (suma2 > 0.0_r8 .and. suma2 <  1.e-6_r8) ) then
          write (6,*) 'Special land units near 100%, but not quite for n,suma =',n,suma
          write (6,*) 'Adjusting special land units to 100%'
          if (pctlak(n) >= 25._r8) then
             pctlak(n) = 100._r8 - (pctwet(n) + pcturb(n) + pctgla(n))
          else if (pctwet(n) >= 25._r8) then
             pctwet(n) = 100._r8 - (pctlak(n) + pcturb(n) + pctgla(n))
          else if (pcturb(n) >= 25._r8) then
             pcturb(n) = 100._r8 - (pctlak(n) + pctwet(n) + pctgla(n))
          else if (pctgla(n) >= 25._r8) then
             pctgla(n) = 100._r8 - (pctlak(n) + pctwet(n) + pcturb(n))
          else
             write (6,*) subname, 'Error: sum of special land units nearly 100% but none is >= 25% at ', &
                  'n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),pctnatveg(n),pctcrop(n),suma = ', &
                  n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),&
                  pctnatpft(n)%get_pct_l2g(),pctcft(n)%get_pct_l2g(),suma
             call abort()
          end if
          call pctnatpft(n)%set_pct_l2g(0._r8)
          call pctcft(n)%set_pct_l2g(0._r8)
       end if
       if ( any(pctnatpft(n)%get_pct_p2g() > 0.0_r8 .and. pctnatpft(n)%get_pct_p2g() < toosmallPFT ) .or. &
            any(pctcft(n)%get_pct_p2g()    > 0.0_r8 .and. pctcft(n)%get_pct_p2g()    < toosmallPFT )) then
          write (6,*) 'pctnatpft or pctcft is small at n=', n
          write (6,*) 'pctnatpft%pct_p2l = ', pctnatpft(n)%get_pct_p2l()
          write (6,*) 'pctcft%pct_p2l = ', pctcft(n)%get_pct_p2l()
          write (6,*) 'pctnatpft%pct_l2g = ', pctnatpft(n)%get_pct_l2g()
          write (6,*) 'pctcft%pct_l2g = ', pctcft(n)%get_pct_l2g()
          call abort()
       end if
       
       suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n)
       if (suma < 100._r8-epsilon(suma) .and. suma > (100._r8 - 4._r8*epsilon(suma))) then
          write (6,*) subname, 'n,pctlak,pctwet,pcturb,pctgla,pctnatveg,pctcrop= ', &
               n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),&
               pctnatpft(n)%get_pct_l2g(), pctcft(n)%get_pct_l2g()
          call abort()
       end if
       suma = suma + pctnatpft(n)%get_pct_l2g() + pctcft(n)%get_pct_l2g()
       if ( abs(suma-100._r8) > 1.e-10_r8) then
          write (6,*) subname, ' error: sum of pctlak, pctwet,', &
               'pcturb, pctgla, pctnatveg and pctcrop is NOT equal to 100'
          write (6,*)'n,pctlak,pctwet,pcturb,pctgla,pctnatveg,pctcrop,sum= ', &
               n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),&
               pctnatpft(n)%get_pct_l2g(),pctcft(n)%get_pct_l2g(), suma
          call abort()
       end if
       
    end do

    ! Check that when pctnatveg+pctcrop identically zero, sum of special landunits is identically 100%

    if ( .not. outnc_double )then
       do n = 1,ns_o
          sum8  =         real(pctlak(n),r4)
          sum8  = sum8  + real(pctwet(n),r4)
          sum8  = sum8  + real(pcturb(n),r4)
          sum8  = sum8  + real(pctgla(n),r4)
          sum4a =         real(pctnatpft(n)%get_pct_l2g(),r4)
          sum4a = sum4a + real(pctcft(n)%get_pct_l2g(),r4)
          if ( sum4a==0.0_r4 .and. sum8 < 100._r4-2._r4*epsilon(sum4a) )then
             write (6,*) subname, ' error: sum of pctlak, pctwet,', &
                  'pcturb, pctgla is < 100% when pctnatveg+pctcrop==0 sum = ', sum8
             write (6,*)'n,pctlak,pctwet,pcturb,pctgla,pctnatveg,pctcrop= ', &
                  n,pctlak(n),pctwet(n),pcturb(n),pctgla(n), &
                  pctnatpft(n)%get_pct_l2g(),pctcft(n)%get_pct_l2g()
             call abort()
          end if
       end do
    else
       do n = 1,ns_o
          sum8  =         pctlak(n)
          sum8  = sum8  + pctwet(n)
          sum8  = sum8  + pcturb(n)
          sum8  = sum8  + pctgla(n)
          sum8a =         pctnatpft(n)%get_pct_l2g()
          sum8a = sum8a + pctcft(n)%get_pct_l2g()
          if ( sum8a==0._r8 .and. sum8 < (100._r8-4._r8*epsilon(sum8)) )then
             write (6,*) subname, ' error: sum of pctlak, pctwet,', &
                  'pcturb, pctgla is < 100% when pctnatveg+pctcrop==0 sum = ', sum8
             write (6,*) 'Total error, error/epsilon = ',100._r8-sum8, ((100._r8-sum8)/epsilon(sum8))
             write (6,*)'n,pctlak,pctwet,pcturb,pctgla,pctnatveg,pctcrop,epsilon= ', &
                  n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),&
                  pctnatpft(n)%get_pct_l2g(),pctcft(n)%get_pct_l2g(), epsilon(sum8)
             call abort()
          end if
       end do
    end if

    ! Make sure that there is no vegetation outside the pft mask
    do n = 1,ns_o
       if (pftdata_mask(n) == 0 .and. (pctnatpft(n)%get_pct_l2g() > 0 .or. pctcft(n)%get_pct_l2g() > 0)) then
          write (6,*)'vegetation found outside the pft mask at n=',n
          write (6,*)'pctnatveg,pctcrop=', pctnatpft(n)%get_pct_l2g(), pctcft(n)%get_pct_l2g()
          call abort()
       end if
    end do

    ! Make sure that sums at the landunit level all add to 100%
    ! (Note that we don't check pctglcmec here, because it isn't computed at the point
    ! that this subroutine is called -- but the check of sum(pctglcmec) is done in
    ! mkglcmecMod)
    ! (Also note that we don't need to check pctnatpft or pctcft, because a similar check
    ! is done internally by the pct_pft_type routines.)
    do n = 1,ns_o
       if (abs(sum(urbn_classes(n,:)) - 100._r8) > 1.e-12_r8) then
          write(6,*) 'sum(urbn_classes(n,:)) != 100: ', n, sum(urbn_classes(n,:))
          call abort()
       end if
    end do

    if ( nsmall_tot > 0 )then
       write (6,*)'number of small pft = ', nsmall_tot
    end if

end subroutine normalizencheck_landuse

end program mksurfdat
