program mksurfdata

  !-----------------------------------------------------------------------
  ! mksurfdata creates land model surface dataset from original "raw"
  ! data files.  Surface dataset contains model grid, pfts, inland
  ! water, glacier, soil texture, soil color, LAI and SAI, urban
  ! fraction, and urban parameters.
  ! -----------------------------------------------------------------------

  ! ======================================================================
  ! Summary of  input namelist
  ! ======================================
  ! Must specify settings for the output grid:
  ! ======================================
  !    mksrf_fgrid_mesh -- mesh for output grid
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
  ! ======================================
  ! Must specify meshes for input high resolutin datasets
  ! ======================================
  !    mksrf_fpft_mesh --------    Mesh for mksrf_fvegtyp
  !    mksrf_flakwat_mesh -----    Mesh for mksrf_flakwat
  !    mksrf_fwetlnd_mesh -----    Mesh for mksrf_fwetlnd
  !    mksrf_fglacier_mesh ----    Mesh for mksrf_fglacier
  !    mksrf_fglacierregion_mesh   Mesh for mksrf_fglacierregion
  !    mksrf_fsoitex_mesh -----    Mesh for mksrf_fsoitex
  !    mksrf_fsoicol_mesh -----    Mesh for mksrf_fsoicol
  !    mksrf_furban_mesh ------    Mesh for mksrf_furban
  !    mksrf_furbtopo_mesh ----    Mesh for mksrf_furbtopo
  !    mksrf_fmax_mesh --------    Mesh for mksrf_fmax
  !    mksrf_forganic_mesh ----    Mesh for mksrf_forganic
  !    mksrf_fvocef_mesh ------    Mesh for mksrf_fvocef
  !    mksrf_flai_mesh --------    Mesh for mksrf_flai
  !    mksrf_fhrv_mesh --------    Mesh for mksrf_flai harvesting
  !    mksrf_fgdp_mesh --------    Mesh for mksrf_fgdp
  !    mksrf_fpeat_mesh -------    Mesh for mksrf_fpeat
  !    mksrf_fsoildepth_mesh --    Mesh for mksrf_fsoildepth
  !    mksrf_fabm_mesh --------    Mesh for mksrf_fabm
  !    mksrf_ftopostats_mesh --    Mesh for mksrf_ftopostats
  !    mksrf_fvic_mesh --------    Mesh for mksrf_fvic
  ! ======================================
  ! Optionally specify setting for:
  ! ======================================
  !    mksrf_fdynuse ----- ASCII text file that lists each year of pft files to use
  !    mksrf_gridtype ---- Type of grid (default is 'global')
  !    outnc_double ------ If output should be in double precision
  !    outnc_large_files - If output should be in NetCDF large file format
  !    outnc_vic --------- Output fields needed for VIC
  !    outnc_3dglc ------- Output 3D glacier fields (normally only needed for comparasion)
  !    nglcec ------------ If you want to change the number of Glacier elevation classes
  !    gitdescribe ------- Description of this version from git
  ! ======================================
  ! Optional settings to change values for entire area
  ! ======================================
  !    all_urban --------- If entire area is urban
  !    all_veg ----------- If entire area is to be vegetated (pft_idx and pft_frc then required)
  !    no_inlandwet ------ If wetland should be set to 0% over land
  !    pft_idx ----------- If you want to change to 100% veg covered with given PFT indices
  !    pft_frc ----------- Fractions that correspond to the pft_idx above
  ! ==================
  !    numpft            (if different than default of 16)
  ! ======================================
  ! Optional settings to work around urban bug?
  ! ======================================
  !    urban_skip_abort_on_invalid_data_check
  ! ======================================================================

  ! !USES:
  use ESMF
  use pio
  use shr_kind_mod       , only : r8 => shr_kind_r8, r4 => shr_kind_r4, cs => shr_kind_cs
  use shr_sys_mod        , only : shr_sys_abort

#ifdef TODO  
  use mklaiMod           , only : mklai
  use mkpctPftTypeMod    , only : pct_pft_type, get_pct_p2l_array, get_pct_l2g_array, update_max_array
  use mkpftConstantsMod  , only : natpft_lb, natpft_ub, cft_lb, cft_ub, num_cft
  use mkpftMod           , only : pft_idx, pft_frc, mkpft, mkpftInit, mkpft_parse_oride
  use mksoilMod          , only : soil_sand, soil_clay, mksoiltex, mksoilInit
  use mksoilMod          , only : soil_color, mksoilcol, mkorganic
  use mksoilMod          , only : soil_fmax, mkfmax
  use mkvocefMod         , only : mkvocef
  use mkglacierregionMod , only : mkglacierregion
  use mkglcmecMod        , only : nglcec, mkglcmec, mkglcmecInit, mkglacier
  use mkharvestMod       , only : mkharvest, mkharvest_init, mkharvest_fieldname
  use mkharvestMod       , only : mkharvest_numtypes, mkharvest_parse_oride
  use mkharvestMod       , only : harvestDataType
  use mkurbanparCommonMod, only : mkelev
  use mkurbanparMod      , only : mkurbanInit, mkurban, mkurbanpar, numurbl
  use mkgdpMod           , only : mkgdp
  use mkpeatMod          , only : mkpeat
  use mksoildepthMod     , only : mksoildepth
  use mkagfirepkmonthMod , only : mkagfirepkmon
  use mktopostatsMod     , only : mktopostats
  use mkVICparamsMod     , only : mkVICparams
#endif
  use mklanwatMod        , only : mklakwat
  use mkutilsMod         , only : normalize_classes_by_gcell, chkerr
  use mkfileMod          , only : mkfile
  use mkvarpar           , only : nlevsoi, elev_thresh, numstdpft
  use nanMod             , only : nan, bigint
  use mkpioMod           , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkpioMod           , only : mkpio_put_time_slice, mkpio_iodesc_output
  use mkvarctl

  implicit none

  ! local variables
  type(ESMF_Mesh)               :: mesh_model 
  integer                       :: nsoicol                 ! number of model color classes
  integer                       :: k,m,n                   ! indices
  integer                       :: ni,nj,ns_o              ! indices
  integer                       :: lsize_o 
  integer                       :: ier                     ! error status
  integer                       :: nfdyn                   ! unit numbers
  integer                       :: omode                   ! netCDF output mode
  integer                       :: ret                     ! netCDF return status
  integer                       :: ntim                    ! time sample for dynamic land use
  integer                       :: year                    ! year for dynamic land use
  integer                       :: year2                   ! year for dynamic land use for harvest file
  real(r8)                      :: suma                    ! sum for error check
  character(len=256)            :: string                  ! string read in
  integer                       :: rcode

  ! pio variables
  type(file_desc_t)             :: pioid
  type(var_desc_t)              :: pio_varid
  type(io_desc_t)               :: pio_iodesc 

  ! data arrays
  real(r8), pointer             :: landfrac_pft(:)         ! PFT data: % land per gridcell
  real(r8), pointer             :: pctlnd_pft(:)           ! PFT data: % of gridcell for PFTs
  real(r8), pointer             :: pctlnd_pft_dyn(:)       ! PFT data: % of gridcell for dyn landuse PFTs
  integer , pointer             :: pftdata_mask(:)         ! mask indicating real or fake land type
#ifdef TODO
  type(pct_pft_type), pointer   :: pctnatpft(:)            ! % of grid cell that is nat veg, and breakdown into PFTs
  type(pct_pft_type), pointer   :: pctnatpft_max(:)        ! % of grid cell maximum PFTs of the time series
  type(pct_pft_type), pointer   :: pctcft(:)               ! % of grid cell that is crop, and breakdown into CFTs
  type(pct_pft_type), pointer   :: pctcft_max(:)           ! % of grid cell maximum CFTs of the time series
#endif
  real(r8)                      :: harvest_initval         ! initial value for harvest variables
  real(r8), pointer             :: harvest1D(:)            ! harvest 1D data: normalized harvesting
  real(r8), pointer             :: harvest2D(:,:)          ! harvest 1D data: normalized harvesting
  real(r8), pointer             :: pctgla(:)               ! percent of grid cell that is glacier
  real(r8), pointer             :: pctglc_gic(:)           ! percent of grid cell that is gic (% of glc landunit)
  real(r8), pointer             :: pctglc_icesheet(:)      ! percent of grid cell that is ice sheet (% of glc landunit)
  real(r8), pointer             :: pctglcmec(:,:)          ! glacier_mec pct coverage in each class (% of landunit)
  real(r8), pointer             :: topoglcmec(:,:)         ! glacier_mec sfc elevation in each gridcell and class
  real(r8), pointer             :: pctglcmec_gic(:,:)      ! GIC pct coverage in each class (% of landunit)
  real(r8), pointer             :: pctglcmec_icesheet(:,:) ! icesheet pct coverage in each class (% of landunit)
  real(r8), pointer             :: elevclass(:)            ! glacier_mec elevation classes
  integer,  pointer             :: glacier_region(:)       ! glacier region ID
  real(r8), pointer             :: pctlak(:)               ! percent of grid cell that is lake
  real(r8), pointer             :: pctwet(:)               ! percent of grid cell that is wetland
  real(r8), pointer             :: pcturb(:)               ! percent of grid cell that is urbanized (total across all urban classes)
  real(r8), pointer             :: urbn_classes(:,:)       ! percent cover of each urban class, as % of total urban area
  real(r8), pointer             :: urbn_classes_g(:,:)     ! percent cover of each urban class, as % of grid cell
  real(r8), pointer             :: elev(:)                 ! glc elevation (m)
  real(r8), pointer             :: fmax(:)                 ! fractional saturated area
  integer , pointer             :: soicol(:)               ! soil color
  real(r8), pointer             :: pctsand(:,:)            ! soil texture: percent sand
  real(r8), pointer             :: pctclay(:,:)            ! soil texture: percent clay
  real(r8), pointer             :: ef1_btr(:)              ! Isoprene emission factor for broadleaf
  real(r8), pointer             :: ef1_fet(:)              ! Isoprene emission factor for fine/everg
  real(r8), pointer             :: ef1_fdt(:)              ! Isoprene emission factor for fine/dec
  real(r8), pointer             :: ef1_shr(:)              ! Isoprene emission factor for shrubs
  real(r8), pointer             :: ef1_grs(:)              ! Isoprene emission factor for grasses
  real(r8), pointer             :: ef1_crp(:)              ! Isoprene emission factor for crops
  real(r8), pointer             :: organic(:,:)            ! organic matter density (kg/m3)
  real(r8), pointer             :: gdp(:)                  ! GDP (x1000 1995 US$/capita)
  real(r8), pointer             :: fpeat(:)                ! peatland fraction of gridcell
  real(r8), pointer             :: soildepth(:)            ! soil depth (m)
  integer , pointer             :: agfirepkmon(:)          ! agricultural fire peak month
  integer , pointer             :: urban_region(:)         ! urban region ID
  real(r8), pointer             :: topo_stddev(:)          ! standard deviation of elevation (m)
  real(r8), pointer             :: slope(:)                ! topographic slope (degrees)
  real(r8), pointer             :: vic_binfl(:)            ! VIC b

  ! parameters (unitless)
  real(r8), pointer             :: vic_ws(:)               ! VIC Ws parameter (unitless)
  real(r8), pointer             :: vic_dsmax(:)            ! VIC Dsmax parameter (mm/day)
  real(r8), pointer             :: vic_ds(:)               ! VIC Ds parameter (unitless)
  real(r8), pointer             :: lakedepth(:)            ! lake depth (m)
  integer , pointer             :: harvind1D(:)            ! Indices of 1D harvest fields
  integer , pointer             :: harvind2D(:)            ! Indices of 2D harvest fields
  logical                       :: zero_out_lake           ! if should zero glacier out
  logical                       :: zero_out_wetland        ! if should zero glacier out
#ifdef TODO
  type(harvestDataType)         :: harvdata
#endif

  ! esmf variables
  integer                       :: rc
  type(ESMF_LogKind_Flag)       :: logkindflag
  type(ESMF_VM)                 :: vm
  character(len=CS)             :: varname
  logical                       :: create_esmf_pet_files = .true.
  character(len=32)             :: subname = 'mksrfdata'    ! program name

  ! NOTE(bja, 2015-01) added to work around a ?bug? causing 1x1_urbanc_alpha to abort. See
  !/glade/p/cesm/cseg/inputdata/lnd/clm2/surfdata_map/README_c141219
  logical :: urban_skip_abort_on_invalid_data_check

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__
  ! ------------------------------------------------------------

  ! ------------------
  ! Initialize MPI
  ! ------------------

  call MPI_init(rc)
  mpicom = mpi_comm_world

  ! ------------------
  ! Initialize ESMF and get mpicom from ESMF
  ! ------------------

  if (create_esmf_pet_files) then
     logkindflag = ESMF_LOGKIND_MULTI
  else
     logkindflag = ESMF_LOGKIND_MULTI_ON_ERROR
  end if
  call ESMF_Initialize(mpiCommunicator=MPICOM, logkindflag=logkindflag, logappendflag=.false., &
       ioUnitLBound=5001, ioUnitUBound=5101, vm=vm, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
  call ESMF_VMGetGlobal(vm, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
  call ESMF_VMGet(vm, mpicommunicator=mpicom, localPet=iam, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
  call ESMF_LogSet(flush=.true.)
  call ESMF_LogWrite("mksurfdata starting", ESMF_LOGMSG_INFO)

  ! ------------------
  ! Initialize PIO - first phase
  ! ------------------

  ! the following returns pio_iosystem
  ! call pio_init(iam, mpicom, 1, 0, 36, PIO_REARR_SUBSET, pio_iosystem)
  ! TODO: generalize iotype

  call pio_init(iam, mpicom, 4, 0, 36, PIO_REARR_SUBSET, pio_iosystem)

  ! call pio_init(comp_rank=iam, &
  !      comp_comm=mpicom, &
  !      num_iotasks = 8, &
  !      num_aggregator = 0, &
  !      stride = 36, &
  !      rearr = PIO_REARR_SUBSET, &
  !      iosystem = pio_iosystem)

  pio_iotype   =  PIO_IOTYPE_PNETCDF
  pio_ioformat =  PIO_64BIT_DATA

  call ESMF_LogWrite("finished initializing PIO", ESMF_LOGMSG_INFO)

  ! ------------------
  ! Read input namelist on root_task and broadcast to all pes
  ! ------------------

  ! root_task is a module variable in mkvarctl
  root_task = (iam == 0)
  call read_namelist_input()

  ! ------------------
  ! open output ndiag file
  ! ------------------

  if (fsurlog == ' ') then
     call shr_sys_abort(' ERROR: must specify fsurlog in namelist')
  end if
  if (root_task) then
     open (newunit=ndiag, file=trim(fsurlog), iostat=ier)
     if (ier /= 0) then
        call shr_sys_abort(' failed to open ndiag file '//trim(fsurlog))
     end if
     write (ndiag,*) 'Attempting to create surface boundary data .....'
     write (ndiag,'(72a1)') ("-",n=1,60)
  else
     ndiag = 6
  end if

  ! ------------------
  ! Write out namelist input to ndiag
  ! ------------------

  call check_namelist_input()
  call write_namelist_input()

  ! ----------------------------------------------------------------------
  ! Call module initialization routines
  ! ----------------------------------------------------------------------

#ifdef TODO
  call mksoilInit( )
  call mkpftInit( zero_out_l=all_urban, all_veg_l=all_veg )
  allocate ( elevclass(nglcec+1) )
  call mkglcmecInit (elevclass)
  call mkurbanInit (mksrf_furban)
#endif 
  ! ----------------------------------------------------------------------
  ! Allocate and initialize dynamic memory for local variables
  ! ----------------------------------------------------------------------

  ! Read in model mesh to determine the number of local points
  mesh_model = ESMF_MeshCreate(filename=trim(mksrf_fgrid_mesh), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

  ! Get the number of local destination points on my processor (lsize_o)
  call ESMF_MeshGet(mesh_model, numOwnedElements=lsize_o, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

  allocate ( landfrac_pft(lsize_o))           ; landfrac_pft(:)     = spval
  allocate ( pctlnd_pft(lsize_o))             ; pctlnd_pft(:)       = spval
  allocate ( pftdata_mask(lsize_o))           ; pftdata_mask(:)     = -999
  !allocate ( pctnatpft(lsize_o))              ;
  !allocate ( pctnatpft_max(lsize_o))          ;
  !allocate ( pctcft(lsize_o))                 ;
  !allocate ( pctcft_max(lsize_o))             ;
  allocate ( pctgla(lsize_o))                 ; pctgla(:)           = spval
  allocate ( pctlak(lsize_o))                 ; pctlak(:)           = spval
  allocate ( pctwet(lsize_o))                 ; pctwet(:)           = spval
  allocate ( pcturb(lsize_o))                 ; pcturb(:)           = spval
  allocate ( urban_region(lsize_o))           ; urban_region(:)     = -999
  !allocate ( urbn_classes(lsize_o,numurbl))   ; urbn_classes(:,:)   = spval
  !allocate ( urbn_classes_g(lsize_o,numurbl)) ; urbn_classes_g(:,:) = spval
  allocate ( pctsand(lsize_o,nlevsoi))        ; pctsand(:,:)        = spval
  allocate ( pctclay(lsize_o,nlevsoi))        ; pctclay(:,:)        = spval
  allocate ( soicol(lsize_o))                 ; soicol(:)           = -999
  allocate ( gdp(lsize_o))                    ; gdp(:)              = spval
  allocate ( fpeat(lsize_o))                  ; fpeat(:)            = spval
  allocate ( soildepth(lsize_o))              ; soildepth(:)        = spval
  allocate ( agfirepkmon(lsize_o))            ; agfirepkmon(:)      = -999
  allocate ( topo_stddev(lsize_o))            ; topo_stddev(:)      = spval
  allocate ( slope(lsize_o))                  ; slope(:)            = spval
  allocate ( vic_binfl(lsize_o))              ; vic_binfl(:)        = spval
  allocate ( vic_ws(lsize_o))                 ; vic_ws(:)           = spval
  allocate ( vic_dsmax(lsize_o))              ; vic_dsmax(:)        = spval
  allocate ( vic_ds(lsize_o))                 ; vic_ds(:)           = spval
  allocate ( lakedepth(lsize_o))              ; lakedepth(:)        = spval
  allocate ( glacier_region(lsize_o))         ; glacier_region(:)   = -999

  ! ----------------------------------------------------------------------
  ! Read in and interpolate surface data set fields
  ! ----------------------------------------------------------------------

  ! ! Make PFTs [pctnatpft, pctcft] from dataset [fvegtyp]
  ! call mkpft( mapfname=map_fpft, fpft=mksrf_fvegtyp, ndiag=ndiag, pctlnd_o=pctlnd_pft, pctnatpft_o=pctnatpft, pctcft_o=pctcft)

  ! ! Create harvesting data at model resolution
  ! if (all_veg) then
  !    ! In this case, we don't call mkharvest, so we want the harvest variables to be
  !    ! initialized reasonably.
  !    harvest_initval = 0._r8
  ! else
  !    harvest_initval = spval
  ! end if
  ! call mkharvest_init( lsize_o, harvest_initval, harvdata, mksrf_fhrvtyp )
  ! if ( .not. all_veg )then
  !    call mkharvest(  mapfname=map_fharvest, datfname=mksrf_fhrvtyp, ndiag=ndiag, harvdata=harvdata )
  ! end if

  ! Make inland water [pctlak, pctwet] [flakwat] [fwetlnd]
  zero_out_lake = all_urban .or. all_veg
  zero_out_wetland = all_urban .or. all_veg .or. no_inlandwet
  call mklakwat(mksrf_flakwat_mesh, mksrf_flakwat, mesh_model, &
       zero_out_lake, zero_out_wetland, pctlak, pctwet, lakedepth, rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mklatwat')

  ! ! Make glacier fraction [pctgla] from [fglacier] dataset
  ! call mkglacier ( mapfname=map_fglacier, datfname=mksrf_fglacier, ndiag=ndiag, zero_out=all_urban.or.all_veg, glac_o=pctgla)

  ! ! Make glacier region ID [glacier_region] from [fglacierregion] dataset
  ! call mkglacierregion ( mapfname=map_fglacierregion, &
  !      datfname=mksrf_fglacierregion, ndiag=ndiag, glacier_region_o = glacier_region)

  ! ! Make soil texture [pctsand, pctclay]  [fsoitex]
  ! call mksoiltex ( mapfname=map_fsoitex, datfname=mksrf_fsoitex, ndiag=ndiag, sand_o=pctsand, clay_o=pctclay)

  ! ! Make soil color classes [soicol] [fsoicol]
  ! call mksoicol ( mapfname=map_fsoicol, datfname=mksrf_fsoicol, ndiag=ndiag, soil_color_o=soicol, nsoicol=nsoicol)

  ! ! Make fmax [fmax] from [fmax] dataset
  ! allocate(fmax(lsize_o)); fmax(:) = spval
  ! call mkfmax ( mapfname=map_fmax, datfname=mksrf_fmax, ndiag=ndiag, fmax_o=fmax)

  ! ! Make GDP data [gdp] from [gdp]
  ! call mkgdp ( mapfname=map_fgdp, datfname=mksrf_fgdp, ndiag=ndiag, gdp_o=gdp)

  ! ! Make peat data [fpeat] from [peatf]
  ! call mkpeat ( mapfname=map_fpeat, datfname=mksrf_fpeat, ndiag=ndiag, peat_o=fpeat)

  ! ! Make soil depth data [soildepth] from [soildepthf]
  ! call mksoildepth ( mapfname=map_fsoildepth, datfname=mksrf_fsoildepth, ndiag=ndiag, soildepth_o=soildepth)

  ! ! Make agricultural fire peak month data [abm] from [abm]
  ! call mkagfirepkmon ( mapfname=map_fabm, datfname=mksrf_fabm, ndiag=ndiag, agfirepkmon_o=agfirepkmon)

  ! ! Make urban fraction [pcturb] from [furban] dataset
  ! call mkurban ( mapfname=map_furban, datfname=mksrf_furban, &
  !      ndiag=ndiag, zero_out=all_veg, urbn_o=pcturb, urbn_classes_o=urbn_classes, region_o=urban_region)

  ! ! Make elevation [elev] from [ftopo, ffrac] dataset
  ! ! Used only to screen pcturb, screen pcturb by elevation threshold from elev dataset
  ! if ( .not. all_urban .and. .not. all_veg )then
  !    allocate(elev(lsize_o))
  !    elev(:) = spval
  !    ! NOTE(wjs, 2016-01-15) This uses the 'TOPO_ICE' variable for historical reasons
  !    ! (this same dataset used to be used for glacier-related purposes as well).
  !    ! TODO(wjs, 2016-01-15) A better solution for this urban screening would probably
  !    ! be to modify the raw urban data; in that case, I believe we could remove furbtopo.
  !    call mkelev ( mapfname=map_furbtopo, datfname=mksrf_furbtopo, varname='TOPO_ICE', ndiag=ndiag, elev_o=elev)
  !    where (elev .gt. elev_thresh)
  !       pcturb = 0._r8
  !    end where
  !    deallocate(elev)
  ! end if

  ! ! Compute topography statistics [topo_stddev, slope] from [ftopostats]
  ! call mktopostats ( mapfname=map_ftopostats, datfname=mksrf_ftopostats, &
  !      ndiag=ndiag, topo_stddev_o=topo_stddev, slope_o=slope, std_elev=std_elev)

  ! ! Make VIC parameters [binfl, ws, dsmax, ds] from [fvic]
  ! if ( outnc_vic )then
  !    call mkVICparams ( mapfname=map_fvic, datfname=mksrf_fvic, ndiag=ndiag, &
  !         binfl_o=vic_binfl, ws_o=vic_ws, dsmax_o=vic_dsmax, ds_o=vic_ds)
  ! end if

  ! ! Make organic matter density [organic] [forganic]
  ! allocate (organic(lsize_o,nlevsoi))
  ! organic(:,:) = spval
  ! call mkorganic ( mapfname=map_forganic, datfname=mksrf_forganic, ndiag=ndiag, organic_o=organic)

  ! ! Make VOC emission factors for isoprene &
  ! ! [ef1_btr,ef1_fet,ef1_fdt,ef1_shr,ef1_grs,ef1_crp]
  ! allocate ( ef1_btr(lsize_o) , &
  !      ef1_fet(lsize_o) , &
  !      ef1_fdt(lsize_o) , &
  !      ef1_shr(lsize_o) , &
  !      ef1_grs(lsize_o) , &
  !      ef1_crp(lsize_o) )
  ! ef1_btr(:) = 0._r8
  ! ef1_fet(:) = 0._r8
  ! ef1_fdt(:) = 0._r8
  ! ef1_shr(:) = 0._r8
  ! ef1_grs(:) = 0._r8
  ! ef1_crp(:) = 0._r8

  ! call mkvocef ( mapfname=map_fvocef, datfname=mksrf_fvocef, ndiag=ndiag, &
  !      ef_btr_o=ef1_btr, ef_fet_o=ef1_fet, ef_fdt_o=ef1_fdt,  &
  !      ef_shr_o=ef1_shr, ef_grs_o=ef1_grs, ef_crp_o=ef1_crp)

  ! ! Do landuse changes such as for the poles, etc.
  ! call change_landuse(  dynpft=.false. )

  ! ----------------------------------------------------------------------
  ! Modify interpolated fields based on additional constrants
  ! ----------------------------------------------------------------------

#ifdef TODO
  do n = 1,lsize_o
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

     ! Assume wetland, glacier and/or lake when dataset landmask implies ocean
     ! (assume medium soil color (15) and loamy texture).
     ! Also set pftdata_mask here

     if (pctlnd_pft(n) < 1.e-6_r8) then
        pftdata_mask(n)  = 0
        soicol(n)        = 15
        if (pctgla(n) < 1.e-6_r8) then
           pctwet(n)    = 100._r8 - pctlak(n)
           pctgla(n)    = 0._r8
        else
           pctwet(n)    = max(100._r8 - pctgla(n) - pctlak(n), 0.0_r8)
        end if
        pcturb(n)        = 0._r8
        !call pctnatpft(n)%set_pct_l2g(0._r8)
        !call pctcft(n)%set_pct_l2g(0._r8)
        pctsand(n,:)     = 43._r8
        pctclay(n,:)     = 18._r8
        organic(n,:)   = 0._r8
     else
        pftdata_mask(n) = 1
     end if

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

  call normalizencheck_landuse()

  !Write out sum of PFT's

  do k = natpft_lb,natpft_ub
     suma = 0._r8
     do n = 1,lsize_o
        suma = suma + pctnatpft(n)%get_one_pct_p2g(k)
     enddo
     write(6,*) 'sum over domain of pft ',k,suma
  enddo
  write(6,*)

  do k = cft_lb,cft_ub
     suma = 0._r8
     do n = 1,lsize_o
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

  allocate (pctglcmec(lsize_o,nglcec), &
       topoglcmec(lsize_o,nglcec) )
  if ( outnc_3dglc )then
     allocate( &
          pctglcmec_gic(lsize_o,nglcec), &
          pctglcmec_icesheet(lsize_o,nglcec))
     allocate (pctglc_gic(lsize_o))
     allocate (pctglc_icesheet(lsize_o))
  end if

  pctglcmec(:,:)          = spval
  topoglcmec(:,:)         = spval

  if ( outnc_3dglc )then
     call mkglcmec ( mapfname=map_fglacier, &
          datfname_fglacier=mksrf_fglacier, ndiag=ndiag, &
          pctglcmec_o=pctglcmec, topoglcmec_o=topoglcmec, &
          pctglcmec_gic_o=pctglcmec_gic, pctglcmec_icesheet_o=pctglcmec_icesheet, &
          pctglc_gic_o=pctglc_gic, pctglc_icesheet_o=pctglc_icesheet)
  else
     call mkglcmec ( mapfname=map_fglacier, &
          datfname_fglacier=mksrf_fglacier, ndiag=ndiag, &
          pctglcmec_o=pctglcmec, topoglcmec_o=topoglcmec )
  end if

  ! Determine fractional land from pft dataset

  do n = 1,lsize_o
     landfrac_pft(n) = pctlnd_pft(n)/100._r8
  end do
#endif

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

#ifdef TODO
     ! The following variables need to be output as int
     ! nf_put_var_int(ncid, varid, (/(n,n=natpft_lb,natpft_ub)/))
     ! nf_put_var_int(ncid, varid, (/(n,n=cft_lb,cft_ub)/))
     ! nf_put_var_int(ncid, varid, pftdata_mask) ! derived in mksurfdata
     ! nf_put_var_int(ncid, varid, nsoicol)      ! set in mksoicol
     ! nf_put_var_int(ncid, varid, soicol)       ! set in mksoicol
     ! nf_put_var_int(ncid, varid, glacier_region)
     ! nf_put_var_int(ncid, varid, agfirepkmon)
     ! nf_put_var_int(ncid, varid, urban_region)
#endif

     !call mkfile( vm, nx, ny, trim(fsurdat), harvdata, dynlanduse = .false., pioid=pioid)
     call mkfile( vm, mksrf_fgrid_mesh_nx, mksrf_fgrid_mesh_ny, trim(fsurdat), dynlanduse=.false., pioid=pioid)

#ifdef TODO
     ! write area, longxy and latixy
     ! TODO: fill this in
     ! call check_ret(nf_open(trim(fname), nf_write, ncid), subname)
     ! File will be in define mode. Set fill mode to "no fill" to optimize performance
     ! call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)
     ! call check_ret(nf_inq_varid(ncid, 'AREA', varid), subname)
     ! call check_ret(nf_put_var_double(ncid, varid, domain%area), subname)
     ! call check_ret(nf_inq_varid(ncid, 'LONGXY', varid), subname)
     ! call check_ret(nf_put_var_double(ncid, varid, domain%lonc), subname)
     ! call check_ret(nf_inq_varid(ncid, 'LATIXY', varid), subname)
     ! call check_ret(nf_put_var_double(ncid, varid, domain%latc), subname)

     ! Synchronize the disk copy of a netCDF dataset with in-memory buffers
     !call check_ret(nf_sync(ncid), subname)

     ! Close grid data dataset
     !call check_ret(nf_close(ncid), subname)

     ! Write fields OTHER THAN lai, sai, heights, and urban parameters to netcdf surface dataset

     ! rcode = pio_inq_varid(pioid, 'natpft', pio_varid)
     ! rcode = pio_put_var_int(pioid, pio_varid, (/(n,n=natpft_lb,natpft_ub)/)))

     ! if (num_cft > 0) then
     !    rcode = pio_inq_varid(pioid, 'cft', pio_varid))
     !    rcode = pio_put_var_int(pioid, pio_varid, (/(n,n=cft_lb,cft_ub)/)))
     ! end if

     varname = 'PFTDATA_MASK'
     call mkpio_iodesc_output(pioid, mesh_model, trim(varname), pio_iodesc, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
     call pio_write_darray(pioid,              pio_varid,       pio_iodesc, pftdata_mask, rcode)
    !call pio_write_darray(io_file(lfile_ind), varid, iodesc, fldptr2(:,n), rcode, fillval=lfillvalue)
     call pio_freedecomp(pioid, pio_iodesc)

     varname = 'LANDFRAC_PFT'
     call mkpio_iodesc_output(pioid, mesh_model, trim(varname), pio_iodesc, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
     call pio_write_darray(pioid, pio_varid, pio_iodesc, landfrac_pft, rcode)
     call pio_freedecomp(pioid, pio_iodesc)

     varname = 'mxsoil_color'
     ! TODO: output this as an integer - not calling darray
     !call pio_write_darray(pioid, pio_varid, pio_iodesc, nsoicol, rcode)

     varname = 'SOIL_COLOR'
     call mkpio_iodesc_output(pioid, mesh_model, trim(varname), pio_iodesc, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
     call pio_write_darray(pioid, pio_varid, pio_iodesc, soicol, rcode)
     call pio_freedecomp(pioid, pio_iodesc)

     varname = 'PCT_SAND'
     call mkpio_iodesc_output(pioid, mesh_model, trim(varname), pio_iodesc, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
     call pio_write_darray(pioid, pio_varid, pio_iodesc, pctsand, rcode)
     call pio_freedecomp(pioid, pio_iodesc)

     varname = 'PCT_CLAY'
     call mkpio_iodesc_output(pioid, mesh_model, trim(varname), pio_iodesc, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
     call pio_write_darray(pioid, pio_varid, pio_iodesc, pctclay, rcode)
     call pio_freedecomp(pioid, pio_iodesc)
#endif

     varname = 'PCT_WETLAND'
     call mkpio_iodesc_output(pioid, mesh_model, trim(varname), pio_iodesc, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
     call pio_write_darray(pioid, pio_varid, pio_iodesc, pctwet, rcode)
     call pio_freedecomp(pioid, pio_iodesc)

     varname = 'PCT_LAKE'
     call mkpio_iodesc_output(pioid, mesh_model, trim(varname), pio_iodesc, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
     call pio_write_darray(pioid, pio_varid, pio_iodesc, pctlak, rcode)
     call pio_freedecomp(pioid, pio_iodesc)

#ifdef TODO
     varname = 'PCT_GLACIER'
     call mkpio_iodesc_output(pioid, mesh_model, trim(varname), pio_iodesc, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
     call pio_write_darray(pioid, pio_varid, pio_iodesc, pctgla, rcode)
     call pio_freedecomp(pioid, pio_iodesc)

     varname = 'GLACIER_REGION'
     call mkpio_iodesc_output(pioid, mesh_model, trim(varname), pio_iodesc, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
     call pio_write_darray(pioid, pio_varid, pio_iodesc, glacier_region, rcode)
     call pio_freedecomp(pioid, pio_iodesc)

     varname = 'PCT_GLC_MEC'
     call mkpio_iodesc_output(pioid, mesh_model, trim(varname), pio_iodesc, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
     call pio_write_darray(pioid, pio_varid, pio_iodesc, pctglcmec, rcode)
     call pio_freedecomp(pioid, pio_iodesc)

     varname = 'GLC_MEC'
     call mkpio_iodesc_output(pioid, mesh_model, trim(varname), pio_iodesc, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//trim(varname))
     call pio_write_darray(pioid, pio_varid, pio_iodesc, elevclass, rcode)
     call pio_freedecomp(pioid, pio_iodesc)

     if ( outnc_3dglc )then
        !'TOPO_GLC_MEC', topoglcmec
        ! 'PCT_GLC_MEC_ICESHEET', pctglcmec_icesheet
        ! 'PCT_GLC_GIC', pctglc_gic
        ! 'PCT_GLC_ICESHEET', pctglc_icesheet
     end if

     !'PCT_URBAN', ut_urbn_classes_g
     !'PCT_NATVEG', get_pct_l2g_array(pctnatpft)

     ! rcode = pio_inq_varid(pioid, 'PCT_CROP', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, get_pct_l2g_array(pctcft)))

     ! rcode = pio_inq_varid(pioid, 'PCT_NAT_PFT', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, get_pct_p2l_array(pctnatpft)))

     ! if (num_cft > 0) then
     !    rcode = pio_inq_varid(pioid, 'PCT_CFT', pio_varid))
     !    rcode = pio_put_var_double(pioid, pio_varid, get_pct_p2l_array(pctcft)))
     ! end if

     ! call harvdata%getFieldsIdx( harvind1D, harvind2D )
     ! do k = 1, harvdata%num1Dfields()
     !    rcode = pio_inq_varid(pioid, trim(mkharvest_fieldname(harvind1D(k),constant=.true.)), pio_varid))
     !    harvest1D => harvdata%get1DFieldPtr( harvind1D(k), output=.true. )
     !    rcode = pio_put_var_double(pioid, pio_varid, harvest1D))
     ! end do
     ! do k = 1, harvdata%num2Dfields()
     !    rcode = pio_inq_varid(pioid, trim(mkharvest_fieldname(harvind2D(k),constant=.true.)), pio_varid))
     !    harvest2D => harvdata%get2DFieldPtr( harvind2D(k), output=.true. )
     !    rcode = pio_put_var_double(pioid, pio_varid, harvest2D))
     ! end do
     ! deallocate( harvind1D, harvind2D )

     ! rcode = pio_inq_varid(pioid, 'FMAX', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, fmax))

     ! rcode = pio_inq_varid(pioid, 'gdp', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, gdp))

     ! rcode = pio_inq_varid(pioid, 'peatf', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, fpeat))


     ! !    rcode = pio_inq_varid(pioid, 'Avg_Depth_Median', pio_varid))
     ! rcode = pio_inq_pio_varid(pioid, 'zbedrock', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, soildepth))

     ! rcode = pio_inq_varid(pioid, 'abm', pio_varid))
     ! rcode = pio_put_var_int(pioid, pio_varid, agfirepkmon))

     ! rcode = pio_inq_pio_varid(pioid, 'SLOPE', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, slope))

     ! rcode = pio_inq_varid(pioid, 'STD_ELEV', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, topo_stddev))

     ! if ( outnc_vic )then
     !    rcode = pio_inq_varid(pioid, 'binfl', pio_varid))
     !    rcode = pio_put_var_double(pioid, pio_varid, vic_binfl))

     !    rcode = pio_inq_varid(pioid, 'Ws', pio_varid))
     !    rcode = pio_put_var_double(pioid, pio_varid, vic_ws))

     !    rcode = pio_inq_varid(pioid, 'Dsmax', pio_varid))
     !    rcode = pio_put_var_double(pioid, pio_varid, vic_dsmax))

     !    rcode = pio_inq_varid(pioid, 'Ds', pio_varid))
     !    rcode = pio_put_var_double(pioid, pio_varid, vic_ds))
     ! end if

     ! rcode = pio_inq_varid(pioid, 'LAKEDEPTH', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, lakedepth))

     ! rcode = pio_inq_varid(pioid, 'EF1_BTR', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, ef1_btr))

     ! rcode = pio_inq_varid(pioid, 'EF1_FET', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, ef1_fet))

     ! rcode = pio_inq_varid(pioid, 'EF1_FDT', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, ef1_fdt))

     ! rcode = pio_inq_varid(pioid, 'EF1_SHR', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, ef1_shr))

     ! rcode = pio_inq_varid(pioid, 'EF1_GRS', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, ef1_grs))

     ! rcode = pio_inq_varid(pioid, 'EF1_CRP', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, ef1_crp))

     ! rcode = pio_inq_varid(pioid, 'ORGANIC', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, organic))

     ! rcode = pio_inq_varid(pioid, 'URBAN_REGION_ID', pio_varid))
     ! rcode = pio_put_var_int(pioid, pio_varid, urban_region))

     ! ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

     ! rcode = pio_sync(pioid))

     ! ----------------------------------------------------------------------
     ! Make Urban Parameters from raw input data and write to surface dataset
     ! Write to netcdf file is done inside mkurbanpar routine
     ! ----------------------------------------------------------------------

     ! write(6,*)'calling mkurbanpar'
     ! call mkurbanpar(datfname=mksrf_furban, pioido=pioid, region_o=urban_region, &
     !      urbn_classes_gcell_o=urbn_classes_g, &
     !      urban_skip_abort_on_invalid_data_check=urban_skip_abort_on_invalid_data_check)

     ! ----------------------------------------------------------------------
     ! Make LAI and SAI from 1/2 degree data and write to surface dataset
     ! Write to netcdf file is done inside mklai routine
     ! ----------------------------------------------------------------------

     ! write(6,*)'calling mklai'
     ! call mklai( mapfname=map_flai, datfname=mksrf_flai, ndiag=ndiag, pioido=pioid )
#endif

     ! Close surface dataset
     call pio_closefile(pioid)

     write (6,'(72a1)') ("-",n=1,60)
     write (6,*)' land model surface data set successfully created for ', &
          'grid of size ',lsize_o

  else  ! fsurdat == ' '

     write (6,*) 'fsurdat is blank: skipping writing surface dataset'

  end if  ! if (fsurdat /= ' ')

  ! Deallocate arrays NOT needed for dynamic-pft section of code

  deallocate ( organic )
  deallocate ( ef1_btr, ef1_fet, ef1_fdt, ef1_shr, ef1_grs, ef1_crp )
  deallocate ( pctglcmec, topoglcmec)
  if ( outnc_3dglc ) deallocate ( pctglc_gic, pctglc_icesheet)
  deallocate ( elevclass )
  deallocate ( fmax )
  deallocate ( pctsand, pctclay )
  deallocate ( soicol )
  deallocate ( gdp, fpeat, agfirepkmon )
  deallocate ( soildepth )
  deallocate ( topo_stddev, slope )
  deallocate ( vic_binfl, vic_ws, vic_dsmax, vic_ds )
  deallocate ( lakedepth )
  deallocate ( glacier_region )

#ifdef TODO
  call harvdata%clean()
#endif

  ! ----------------------------------------------------------------------
  ! Create dynamic land use dataset if appropriate
  ! ----------------------------------------------------------------------

  if (mksrf_fdynuse /= ' ') then

     ! write(6,*)'creating dynamic land use dataset'

     ! allocate(pctlnd_pft_dyn(lsize_o))
     ! call mkharvest_init( lsize_o, spval, harvdata, mksrf_fhrvtyp )

     ! if (fdyndat == ' ') then
     !    write(6,*)' must specify fdyndat in namelist if mksrf_fdynuse is not blank'
     !    call abort()
     ! end if

     ! ! Define dimensions and global attributes

     ! call mkfile( fdyndat, harvdata, dynlanduse=.true., pioid)

     ! ! Write fields other pft to dynamic land use dataset

     ! call domain_write( fdyndat)

     ! rcode = pio_open(trim(fdyndat), nf_write, pioid))
     ! rcode = pio_set_fill (pioid, nf_nofill, omode))

     ! rcode = pio_inq_varid(pioid, 'natpft', pio_varid))
     ! rcode = pio_put_var_int(pioid, pio_varid, (/(n,n=natpft_lb,natpft_ub)/))

     ! if (num_cft > 0) then
     !    rcode = pio_inq_varid(pioid, 'cft', pio_varid))
     !    rcode = pio_put_var_int(pioid, pio_varid, (/(n,n=cft_lb,cft_ub)/)
     ! end if

     ! rcode = pio_inq_varid(pioid, 'PFTDATA_MASK', pio_varid))
     ! rcode = pio_put_var_int(pioid, pio_varid, pftdata_mask))

     ! rcode = pio_inq_varid(pioid, 'LANDFRAC_PFT', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, landfrac_pft))

     ! ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

     ! rcode = pio_sync(pioid))

     ! ! Read in each dynamic pft landuse dataset

     ! call opnfil (newunit=nfdyn, mksrf_fdynuse, nfdyn, 'f')

     ! pctnatpft_max = pctnatpft
     ! pctcft_max = pctcft

     ! ntim = 0
     ! do
     !    ! Read input pft data

     !    read(nfdyn, '(A195,1x,I4)', iostat=ier) string, year
     !    if (ier /= 0) exit
     !    !
     !    ! If pft fraction override is set, than intrepret string as PFT and harvesting override values
     !    !
     !    if ( all_veg )then
     !       fname = ' '
     !       fhrvname  = ' '
     !       call mkpft_parse_oride(string)
     !       call mkharvest_parse_oride(string)
     !       write(6, '(a, i4, a)') 'PFT and harvesting values for year ', year, ' :'
     !       write(6, '(a, a)') '    ', trim(string)
     !       !
     !       ! Otherwise intrepret string as a filename with PFT and harvesting values in it
     !       !
     !    else
     !       fname = string
     !       write(6,*)'input pft dynamic dataset for year ', year, ' is : ', trim(fname)
     !       read(nfdyn, '(A195,1x,I4)', iostat=ier) fhrvname, year2
     !       if ( year2 /= year ) then
     !          write(6,*) subname, ' error: year for harvest not equal to year for PFT files'
     !          call abort()
     !       end if
     !    end if
     !    ntim = ntim + 1

     !    ! Create pctpft data at model resolution

     !    call mkpft( mapfname=map_fpft, fpft=fname, &
     !         ndiag=ndiag, pctlnd_o=pctlnd_pft_dyn, pctnatpft_o=pctnatpft, pctcft_o=pctcft )

     !    ! Create harvesting data at model resolution

     !    call mkharvest(  mapfname=map_fharvest, datfname=fhrvname, &
     !         ndiag=ndiag, harvdata=harvdata )

     !    ! Consistency check on input land fraction

     !    do n = 1,lsize_o
     !       if (pctlnd_pft_dyn(n) /= pctlnd_pft(n)) then
     !          write(6,*) subname,' error: pctlnd_pft for dynamics data = ',&
     !               pctlnd_pft_dyn(n), ' not equal to pctlnd_pft for surface data = ',&
     !               pctlnd_pft(n),' at n= ',n
     !          if ( trim(fname) == ' ' )then
     !             write(6,*) ' PFT string = ', string
     !          else
     !             write(6,*) ' PFT file = ', fname
     !          end if
     !          call abort()
     !       end if
     !    end do

     !    call change_landuse( dynpft=.true.)

     !    call normalizencheck_landuse(ldomain)

     !    call update_max_array(pctnatpft_max,pctnatpft)
     !    call update_max_array(pctcft_max,pctcft)

     !    ! Output time-varying data for current year

     !    rcode = pio_inq_varid(pioid, 'PCT_NAT_PFT', pio_varid))
     !    call mkpio_put_time_slice(pioid, pio_varid, ntim, get_pct_p2l_array(pctnatpft))

     !    rcode = pio_inq_varid(pioid, 'PCT_CROP', pio_varid))
     !    call mkpio_put_time_slice(pioid, pio_varid, ntim, get_pct_l2g_array(pctcft))

     !    if (num_cft > 0) then
     !       rcode = pio_inq_varid(pioid, 'PCT_CFT', pio_varid))
     !       call mkpio_put_time_slice(pioid, pio_varid, ntim, get_pct_p2l_array(pctcft))
     !    end if

     !    call harvdata%getFieldsIdx( harvind1D, harvind2D )
     !    do k = 1, harvdata%num1Dfields()
     !       rcode = pio_inq_varid(pioid, trim(mkharvest_fieldname(harvind1D(k),constant=.false.)), pio_varid))
     !       harvest1D => harvdata%get1DFieldPtr( harvind1D(k), output=.true. )
     !       call mkpio_put_time_slice(pioid, pio_varid, ntim, harvest1D)
     !    end do
     !    do k = 1, harvdata%num2Dfields()
     !       rcode = pio_inq_varid(pioid, trim(mkharvest_fieldname(harvind2D(k),constant=.false.)), pio_varid))
     !       harvest2D => harvdata%get2DFieldPtr( harvind2D(k), output=.true. )
     !       call mkpio_put_time_slice(pioid, pio_varid, ntim, harvest2D)
     !    end do
     !    deallocate( harvind1D, harvind2D )

     !    rcode = pio_inq_varid(pioid, 'YEAR', pio_varid))
     !    rcode = pio_put_vara_int(pioid, pio_varid, ntim, 1, year))

     !    rcode = pio_inq_varid(pioid, 'time', pio_varid))
     !    rcode = pio_put_vara_int(pioid, pio_varid, ntim, 1, year))

     !    rcode = pio_inq_varid(pioid, 'input_pftdata_filename', pio_varid))
     !    rcode = pio_put_vara_text(pioid, pio_varid, (/ 1, ntim /), (/ len_trim(string), 1 /), trim(string) ))

     !    ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

     !    rcode = pio_sync(pioid))

     ! end do   ! end of read loop

     ! rcode = pio_inq_varid(pioid, 'PCT_NAT_PFT_MAX', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, get_pct_p2l_array(pctnatpft_max)))

     ! rcode = pio_inq_varid(pioid, 'PCT_CROP_MAX', pio_varid))
     ! rcode = pio_put_var_double(pioid, pio_varid, get_pct_l2g_array(pctcft_max)))

     ! if (num_cft > 0) then
     !    rcode = pio_inq_varid(pioid, 'PCT_CFT_MAX', pio_varid))
     !    rcode = pio_put_var_double(pioid, pio_varid, get_pct_p2l_array(pctcft_max)))
     ! end if

     ! rcode = pio_closefile(pioid))

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

#ifdef TODO
  subroutine change_landuse(  dynpft, lsize_o, mesh_o)

    ! Do landuse changes such as for the poles, etc.
    ! If have pole points on grid - set south pole to glacier
    ! north pole is assumed as non-land

    ! input/output variables
    logical, intent(in)  :: dynpft   ! if part of the dynpft section of code

    ! !LOCAL VARIABLES:
    integer  :: n,lsize_o                       ! indices
    character(len=32) :: subname = 'change_landuse'  ! subroutine name
    !-----------------------------------------------------------------------

    ! 

    do n = 1,lsize_o
       ! TODO: define latc
       if (abs(latc(n) - 90._r8) < 1.e-6_r8) then
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
#endif

#ifdef TODO
  !-----------------------------------------------------------------------
  subroutine normalizencheck_landuse(ldomain)
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
       if ( pctlak(n) < 0.0_r8 )then
          write(6,*) subname, ' ERROR: pctlak is negative!'
          write(6,*) 'n, pctlak = ', n, pctlak(n)
          call abort()
       end if
       if ( pctwet(n) < 0.0_r8 )then
          write(6,*) subname, ' ERROR: pctwet is negative!'
          write(6,*) 'n, pctwet = ', n, pctwet(n)
          call abort()
       end if
       if ( pcturb(n) < 0.0_r8 )then
          write(6,*) subname, ' ERROR: pcturb is negative!'
          write(6,*) 'n, pcturb = ', n, pcturb(n)
          call abort()
       end if
       if ( pctgla(n) < 0.0_r8 )then
          write(6,*) subname, ' ERROR: pctgla is negative!'
          write(6,*) 'n, pctgla = ', n, pctgla(n)
          call abort()
       end if

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
#endif

end program mksurfdata
