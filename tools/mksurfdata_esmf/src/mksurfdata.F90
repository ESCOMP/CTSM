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
  ! Must specify settings for input high resolution datafiles and corresponding meshes
  ! ======================================
  !    mksrf_fglacier            - Glacier dataset
  !    mksrf_fglacier_mesh       - Mesh for mksrf_fglacier
  !    mksrf_fglacierregion      - Glacier region ID dataset
  !    mksrf_fglacierregion_mesh - Mesh for mksrf_fglacierregion
  !    mksrf_flai                - Leaf Area Index dataset
  !    mksrf_flai_mesh           - Mesh for mksrf_flai
  !    mksrf_flakwat             - Lake water dataset
  !    mksrf_flakwat_mesh        - Mesh for mksrf_flakwat
  !    mksrf_fwetlnd             - Wetland water dataset
  !    mksrf_fwetlnd_mesh        - Mesh for mksrf_fwetlnd
  !    mksrf_forganic            - Organic soil carbon dataset
  !    mksrf_forganic_mesh       - Mesh for mksrf_forganic
  !    mksrf_fmax                - Max fractional saturated area dataset
  !    mksrf_fmax_mesh           - Mesh for mksrf_fmax
  !    mksrf_fsoicol             - Soil color dataset
  !    mksrf_fsoicol_mesh        - Mesh for mksrf_fsoicol
  !    mksrf_fsoitex             - Soil texture dataset
  !    mksrf_fsoitex_mesh        - Mesh for mksrf_fsoitex
  !    mksrf_furbtopo            - Topography dataset (for limiting urban areas)
  !    mksrf_furbtopo_mesh       - Mesh for mksrf_furbtopo
  !    mksrf_furban              - Urban dataset
  !    mksrf_furban_mesh         - Mesh for mksrf_furban
  !    mksrf_fvegtyp             - PFT vegetation type dataset
  !    mksrf_fpft_mesh           - Mesh for mksrf_fvegtyp
  !    mksrf_fhrvtyp             - harvest type dataset
  !    mksrf_fhrvtyp_mesh        - Mesh for mksrf_flai harvesting
  !    mksrf_fvocef              - Volatile Organic Compund Emission Factor dataset
  !    mksrf_fvocef_mesh         - Mesh for mksrf_fvocef
  !    mksrf_fgdp                - GDP dataset
  !    mksrf_fgdp_mesh           - Mesh for mksrf_fgdp
  !    mksrf_fpeat               - Peatland dataset
  !    mksrf_fpeat_mesh          - Mesh for mksrf_fpeat
  !    mksrf_fsoildepth          - Soil depth dataset
  !    mksrf_fsoildepth_mesh     - Mesh for mksrf_fsoildepth
  !    mksrf_fabm                - Agricultural fire peak month dataset
  !    mksrf_fabm_mesh           - Mesh for mksrf_fabm
  !    mksrf_ftopostats          - Topography statistics dataset
  !    mksrf_ftopostats_mesh     - Mesh for mksrf_ftopostats
  !    mksrf_fvic                - VIC parameters dataset
  !    mksrf_fvic_mesh           - Mesh for mksrf_fvic
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
  !    numpft ------------ Iif different than default of 16
  !    urban_skip_abort_on_invalid_data_check--- work around urban bug
  ! ======================================
  ! Note: the folloiwng Optional settings have been REMOVED -
  !  instead should now use tools subset_data and modify_fsurdat
  ! ======================================
  !    all_veg ----------- If entire area is to be vegetated (pft_idx and pft_frc then required)
  !    all_urban --------- If entire area is urban
  !    no_inlandwet ------ If wetland should be set to 0% over land
  !    soil_clay --------- If you want to change the soil_clay % to this value everywhere
  !    soil_fmax --------- If you want to change the soil_fmax  to this value everywhere
  !    soil_sand --------- If you want to change the soil_sand % to this value everywhere
  !    pft_idx ----------- If you want to change to 100% veg covered with given PFT indices
  !    pft_frc ----------- Fractions that correspond to the pft_idx above
  ! ======================================================================

  use ESMF
  use pio
  use shr_kind_mod       , only : r8 => shr_kind_r8, r4 => shr_kind_r4, cs => shr_kind_cs, cl => shr_kind_cl
  use shr_sys_mod        , only : shr_sys_abort
  use mkVICparamsMod     , only : mkVICparams
  use mktopostatsMod     , only : mktopostats
  use mkpftMod           , only : mkpft, mkpftInit
  use mkpctPftTypeMod    , only : pct_pft_type, get_pct_p2l_array, get_pct_l2g_array, update_max_array
  use mkpftConstantsMod  , only : natpft_lb, natpft_ub, cft_lb, cft_ub, num_cft, num_natpft
  use mkdomainMod        , only : mkdomain
  use mkharvestMod       , only : mkharvest
  use mkgdpMod           , only : mkgdp
  use mkagfirepkmonthMod , only : mkagfirepkmon
  use mklaiMod           , only : mklai
  use mkpeatMod          , only : mkpeat
  use mkvocefMod         , only : mkvocef
  use mkglcmecMod        , only : mkglcmecInit, mkglcmec, mkglacier
  use mkglacierregionMod , only : mkglacierregion
  use mksoiltexMod       , only : mksoiltex
  use mksoilfmaxMod      , only : mksoilfmax
  use mksoildepthMod     , only : mksoildepth
  use mksoilcolMod       , only : mksoilcol
  use mkurbanparMod      , only : mkurbanInit, mkurban, mkurbanpar, mkurban_topo, numurbl
  use mklanwatMod        , only : mklakwat
  use mkorganicMod       , only : mkorganic
  use mkutilsMod         , only : normalize_classes_by_gcell, chkerr
  use mkfileMod          , only : mkfile_define_dims, mkfile_define_atts, mkfile_define_vars
  use mkfileMod          , only : mkfile_output
  use mkvarpar           , only : nlevsoi, elev_thresh, numstdpft
  use nanMod             , only : nan, bigint
  use mkpioMod           , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkpioMod           , only : mkpio_put_time_slice, mkpio_iodesc_output, mkpio_wopen
  use mkinputMod
  use mkvarctl

  implicit none

#include <mpif.h>

  ! indices
  integer                         :: k,n                     ! indices
  integer                         :: lsize_o

  ! error status
  integer                         :: ier,rcode               ! error status

  ! dynamic land use
  integer                         :: nfdyn                   ! unit numbers
  integer                         :: ntim                    ! time sample for dynamic land use
  integer                         :: year                    ! year for dynamic land use
  integer                         :: year2                   ! year for dynamic land use for harvest file
  real(r8)                        :: suma                    ! sum for error check

  ! model grid
  real(r8), allocatable           :: lon(:)
  real(r8), allocatable           :: lat(:)

  ! pct vegetation data
  real(r8), allocatable           :: landfrac_pft(:)         ! PFT data: % land per gridcell
  real(r8), allocatable           :: pctlnd_pft(:)           ! PFT data: % of gridcell for PFTs
  integer , allocatable           :: pftdata_mask(:)         ! mask indicating real or fake land type
  type(pct_pft_type), allocatable :: pctnatpft(:)            ! % of grid cell that is nat veg, and breakdown into PFTs
  type(pct_pft_type), allocatable :: pctcft(:)               ! % of grid cell that is crop, and breakdown into CFTs

  ! dynamic land use
  real(r8)           , allocatable :: pctlnd_pft_dyn(:)       ! PFT data: % of gridcell for dyn landuse PFTs
  type(pct_pft_type) , allocatable :: pctnatpft_max(:)        ! % of grid cell maximum PFTs of the time series
  type(pct_pft_type) , allocatable :: pctcft_max(:)           ! % of grid cell maximum CFTs of the time series
  real(r8)           , allocatable :: pctnatveg(:)
  real(r8)           , allocatable :: pctcrop(:)
  real(r8)           , allocatable :: pct_nat_pft(:,:)
  real(r8)           , allocatable :: pct_cft(:,:)
  logical                          :: end_of_fdynloop

  ! inland water data, glacier data and urban data
  real(r8), allocatable           :: pctlak(:)               ! percent of grid cell that is lake
  real(r8), allocatable           :: pctwet(:)               ! percent of grid cell that is wetland
  real(r8), allocatable           :: pctgla(:)               ! percent of grid cell that is glacier
  integer , allocatable           :: urban_region(:)         ! urban region ID
  real(r8), allocatable           :: pcturb(:)               ! percent of grid cell that is urbanized (total across all urban classes)
  real(r8), allocatable           :: urban_classes(:,:)      ! percent cover of each urban class, as % of total urban area
  real(r8), allocatable           :: urban_classes_g(:,:)    ! percent cover of each urban class, as % of grid cell
  real(r8), allocatable           :: elev(:)                 ! glc elevation (m)

  ! pio/esmf variables
  type(file_desc_t)               :: pioid
  type(var_desc_t)                :: pio_varid
  type(io_desc_t)                 :: pio_iodesc
  integer                         :: petcount
  integer                         :: stride
  type(ESMF_Mesh)                 :: mesh_model
  type(ESMF_Field)                :: field_model
  type(ESMF_LogKind_Flag)         :: logkindflag
  type(ESMF_VM)                   :: vm
  integer                         :: rc
  logical                         :: create_esmf_pet_files = .true.

  ! character variables
  character(len=CL)               :: string                  ! string read in
  character(len=CL)               :: fname
  character(len=*), parameter     :: subname = 'mksrfdata'   ! program name

  character(len=*), parameter :: u_FILE_u = &
       __FILE__
  ! ------------------------------------------------------------

  ! ======================================================================
  ! Initialize MPI
  ! ======================================================================

  call MPI_init(rc)
  mpicom = mpi_comm_world

  ! ======================================================================
  ! Initialize ESMF and get mpicom from ESMF
  ! ======================================================================

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
  call ESMF_VMGet(vm, mpicommunicator=mpicom, localPet=iam, petcount=petcount, &
       ssiLocalPetCount=stride, rc=rc)
  call ESMF_LogSet(flush=.true.)
  call ESMF_LogWrite("mksurfdata starting", ESMF_LOGMSG_INFO)

  ! Determine root task
  root_task = (iam == 0)

  ! ======================================================================
  ! Initialize PIO
  ! ======================================================================

  ! the following returns pio_iosystem
  call pio_init(iam, mpicom, max(1,petcount/stride), 0, stride, PIO_REARR_SUBSET, pio_iosystem)

  pio_iotype = PIO_IOTYPE_PNETCDF
  pio_ioformat =  PIO_64BIT_DATA

  call ESMF_LogWrite("finished initializing PIO", ESMF_LOGMSG_INFO)

  ! ======================================================================
  ! Read in namelist
  ! ======================================================================

  ! Read input namelist on root_task and broadcast to all pes
  ! root_task is a module variable in mkvarctl
  ! Assumne that input namelist file is 'mksurfdata_in'
  call read_namelist_input(filename='mksurfdata_in')

  ! open output ndiag file
  if (fsurlog == ' ') then
     call shr_sys_abort(' ERROR: must specify fsurlog in namelist')
  end if
  if (root_task) then
     open (newunit=ndiag, file=trim(fsurlog), iostat=ier)
     if (ier /= 0) then
        call shr_sys_abort(' failed to open ndiag file '//trim(fsurlog))
     end if
     write (ndiag,'(a)') 'Attempting to create surface boundary data .....'
     write (ndiag,'(72a1)') ("-",n=1,60)
  else
     ndiag = 6
  end if

  ! Write out namelist input to ndiag
  call check_namelist_input()
  call write_namelist_input()

  ! ======================================================================
  ! Create fsurdat
  ! ======================================================================

  ! Read in model mesh to determine the number of local points
  mesh_model = ESMF_MeshCreate(filename=trim(mksrf_fgrid_mesh), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

  ! Get the number of local destination points on my processor (lsize_o)
  call ESMF_MeshGet(mesh_model, numOwnedElements=lsize_o, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

  ! Initialize urban dimensions (needed to initialize the dimensions in fsurdat)
  call mkurbanInit(mksrf_furban)

  ! Initialize pft/cft dimensions (needed to initialize the dimensions in fsurdat)
  call mkpftInit( )

  ! If fsurdat is blank, then we do not write a surface dataset - but we may still
  ! write a dynamic landuse file. This is useful if we are creating many datasets at
  ! once, and don't want duplicate surface datasets.
  !
  ! TODO(wjs, 2016-01-26) Ideally, we would also avoid doing the processing of
  ! variables that are just needed by the surface dataset (not by the dynamic landuse
  ! file). However, this would require some analysis of the above code, to determine
  ! which processing is needed (directly or indirectly) in order to create a dynamic
  ! landuse file.

  ! Open fsurdat and write out variables
  if (fsurdat == ' ') then
     if (root_task) then
        write (ndiag,'(a)') ' fsurdat is blank: skipping writing surface dataset'
     end if
  else
     if (root_task)then
        write(ndiag,*)
        write(ndiag,'(1x,80a1)') ('=',k=1,80)
     end if

     ! Open file
     ! TODO: what about setting no fill values?
     call mkpio_wopen(trim(fsurdat), clobber=.true., pioid=pioid)

     ! Define dimensions
     call mkfile_define_dims(pioid, nx=mksrf_fgrid_mesh_nx, ny=mksrf_fgrid_mesh_ny, dynlanduse=.false.)

     ! Define global attributes
     call mkfile_define_atts(pioid, dynlanduse = .false.)

     ! Define variables
     call mkfile_define_vars(pioid, dynlanduse = .false.)

     ! End define model
     rcode = pio_enddef(pioid)
  end if

  ! NOTE: do not deallocate pctlak, pctwet, pctglacier and pcturban

  ! -----------------------------------
  ! Write out coordinate variables
  ! -----------------------------------
  rcode = pio_inq_varid(pioid, 'natpft', pio_varid)
  rcode = pio_put_var(pioid, pio_varid, (/(n,n=natpft_lb,natpft_ub)/))
  if (num_cft > 0) then
     rcode = pio_inq_varid(pioid, 'cft', pio_varid)
     rcode = pio_put_var(pioid, pio_varid, (/(n,n=cft_lb,cft_ub)/))
  end if

  ! -----------------------------------
  ! Make lats/lons of model
  ! -----------------------------------
  allocate (lon(lsize_o)) ; lon(:) = spval
  allocate (lat(lsize_o)) ; lat(:) = spval
  if (fsurdat /= ' ') then
     call mkdomain(mesh_model, lon_o=lon, lat_o=lat, pioid_o=pioid, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkdomain')
  end if

  ! -----------------------------------
  ! Make LAI and SAI from 1/2 degree data and write to surface dataset
  ! Write to netcdf file is done inside mklai routine
  ! -----------------------------------
  if (root_task) then
     write(ndiag,'(a)')'calling mklai'
  end if
  call mklai(mksrf_flai_mesh, mksrf_flai, mesh_model, pioid, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mklai')

  ! -----------------------------------
  ! Make PFTs [pctnatpft, pctcft] from dataset [fvegtyp]
  ! Make landfrac_pft and pftdata_mask
  ! -----------------------------------
  ! Determine fractional land from pft dataset
  allocate(pctlnd_pft(lsize_o)); pctlnd_pft(:) = spval
  allocate(pctnatpft(lsize_o)) ;
  allocate(pctcft(lsize_o))    ;
  allocate(pftdata_mask(lsize_o))  ; pftdata_mask(:) = -999
  allocate(landfrac_pft(lsize_o))  ; landfrac_pft(:) = spval
  call mkpft( mksrf_fvegtyp_mesh, mksrf_fvegtyp, mesh_model, &
       pctlnd_o=pctlnd_pft, pctnatpft_o=pctnatpft, pctcft_o=pctcft, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkdomain')

  ! If have pole points on grid - set south pole to glacier
  ! north pole is assumed as non-land
  do n = 1,lsize_o
     if (abs((lat(n) - 90._r8)) < 1.e-6_r8) then
        call pctnatpft(n)%set_pct_l2g(0._r8)
        call pctcft(n)%set_pct_l2g(0._r8)
     end if
     if (pctlnd_pft(n) < 1.e-6_r8) then
        call pctnatpft(n)%set_pct_l2g(0._r8)
        call pctcft(n)%set_pct_l2g(0._r8)
     end if
     landfrac_pft(n) = pctlnd_pft(n)/100._r8
     if (pctlnd_pft(n) < 1.e-6_r8) then
        pftdata_mask(n) = 0
     else
        pftdata_mask(n) = 1
     end if
  end do
  if (fsurdat /= ' ') then
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing land mask from pft dataset"
     call mkfile_output(pioid,  mesh_model, 'PFTDATA_MASK', pftdata_mask, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing land fraction  from pft dataset"
     call mkfile_output(pioid,  mesh_model, 'LANDFRAC_PFT', landfrac_pft, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
     call pio_syncfile(pioid)
  end if

  ! -----------------------------------
  ! Make constant harvesting data at model resolution
  ! -----------------------------------
  ! Note that this call must come after call to mkpftInit - since num_cft is set there
  ! Output data is written in mkharvest
  call mkharvest( mksrf_fhrvtyp_mesh, mksrf_fhrvtyp, mesh_model, pioid_o=pioid, constant=.true., rc=rc )
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkharvest_init')

  ! -----------------------------------
  ! Make inland water [pctlak, pctwet] [flakwat] [fwetlnd]
  ! -----------------------------------
  ! LAKEDEPTH is written out in the subroutine
  ! Need to keep pctlak and pctwet external for use below
  allocate ( pctlak(lsize_o)) ; pctlak(:) = spval
  allocate ( pctwet(lsize_o)) ; pctwet(:) = spval
  call mklakwat(mksrf_flakwat_mesh, mksrf_flakwat, mesh_model, pctlak, pctwet, pioid, fsurdat, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mklatwat')

  ! -----------------------------------
  ! Make glacier fraction [pctgla] from [fglacier] dataset
  ! -----------------------------------
  allocate (pctgla(lsize_o)) ; pctgla(:) = spval
  call mkglacier (mksrf_fglacier_mesh, mksrf_fglacier, mesh_model, glac_o=pctgla, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkglacier')

  ! -----------------------------------
  ! Make glacier region ID [glacier_region] from [fglacierregion] dataset
  ! -----------------------------------
  if (fsurdat /= ' ') then
     ! GLACIER_REGION is written out in the subroutine
     call mkglacierregion(mksrf_fglacierregion_mesh, mksrf_fglacierregion, mesh_model, pioid, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkglacierregion')
  end if

  ! -----------------------------------
  ! Make soil texture [pctsand, pctclay]
  ! -----------------------------------
  if (fsurdat /= ' ') then
     ! mapunits, PCT_SAND and PCT_CLAY are written out in the subroutine
     call mksoiltex( mksrf_fsoitex_mesh, mksrf_fsoitex, mesh_model, pioid, pctlnd_pft, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mksoiltex')
  end if

  ! -----------------------------------
  ! Make soil color classes [soicol] [fsoicol]
  ! -----------------------------------
  if (fsurdat /= ' ') then
     ! SOIL_COLOR and mxsoil_color is written out in the subroutine
     call mksoilcol( mksrf_fsoicol, mksrf_fsoicol_mesh, mesh_model, pctlnd_pft, pioid, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mksoilcol')
  end if

  ! -----------------------------------
  ! Make soil fmax [fmaxsoil]
  ! -----------------------------------
  if (fsurdat /= ' ') then
     ! FMAX is written out in the subroutine
     call mksoilfmax( mksrf_fmax_mesh, mksrf_fmax, mesh_model, pioid, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mksoilfmax')
  end if

  ! -----------------------------------
  ! Make GDP data [gdp] from [gdp]
  ! -----------------------------------
  if (fsurdat /= ' ') then
     ! gdp is written out in the subroutine
     call mkgdp (mksrf_fgdp_mesh, mksrf_fgdp, mesh_model, pioid, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkdomain')
  end if

  ! -----------------------------------
  ! Make peat data [fpeat] from [peatf]
  ! -----------------------------------
  if (fsurdat /= ' ') then
     call mkpeat (mksrf_fpeat_mesh, mksrf_fpeat, mesh_model, pioid, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkpeat')
  end if

  ! -----------------------------------
  ! Make soil depth data [soildepth] from [soildepthf]
  ! -----------------------------------
  if (fsurdat /= ' ') then
     call mksoildepth( mksrf_fsoildepth_mesh, mksrf_fsoildepth, mesh_model, pioid, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mksoildepth')
  end if

  ! -----------------------------------
  ! Make agricultural fire peak month data [abm] from [abm]
  ! -----------------------------------
  if (fsurdat /= ' ') then
     call mkagfirepkmon (mksrf_fabm_mesh, mksrf_fabm, mesh_model, pioid, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkagfirepkmon')
  end if

  ! -----------------------------------
  ! Make urban fraction [pcturb] from [furban] dataset and
  ! -----------------------------------
  allocate (pcturb(lsize_o))                 ; pcturb(:)            = spval
  allocate (urban_classes(lsize_o,numurbl))  ; urban_classes(:,:)   = spval
  allocate (urban_region(lsize_o))           ; urban_region(:)      = -999
  call mkurban(mksrf_furban_mesh, mksrf_furban, mesh_model, pcturb, urban_classes, urban_region, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkurban')
  if (fsurdat /= ' ') then
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out urban region id"
     call mkfile_output(pioid,  mesh_model, 'URBAN_REGION_ID', urban_region, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
     call pio_syncfile(pioid)
  end if
  ! Note that final values of urban are not set and written out until further down

  ! Adjust pcturb
  ! Make elevation [elev] from [ftopo, ffrac] dataset
  ! Used only to screen pcturb, screen pcturb by elevation threshold from elev dataset
  allocate(elev(lsize_o))
  elev(:) = spval
  ! NOTE(wjs, 2016-01-15) This uses the 'TOPO_ICE' variable for historical reasons
  ! (this same dataset used to be used for glacier-related purposes as well).
  ! TODO(wjs, 2016-01-15) A better solution for this urban screening would probably
  ! be to modify the raw urban data; in that case, I believe we could remove furbtopo.
  call mkurban_topo ( mksrf_furbtopo_mesh, mksrf_furbtopo, mesh_model, varname='TOPO_ICE', elev_o=elev, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkurban_topo')
  where (elev > elev_thresh)
     pcturb = 0._r8
  end where
  deallocate(elev)

  ! -----------------------------------
  ! Compute topography statistics [topo_stddev, slope] from [ftopostats]
  ! -----------------------------------
  if (fsurdat /= ' ') then
     call mktopostats ( mksrf_ftopostats_mesh, mksrf_ftopostats, mesh_model, pioid, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mktopostats')
  end if

  ! -----------------------------------
  ! Compute VIC parameters
  ! -----------------------------------
  if (fsurdat /= ' ') then
     if (outnc_vic) then
        call mkVICparams ( mksrf_fvic_mesh, mksrf_fvic, mesh_model, pioid, rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkorganic')
     end if
  end if

  ! -----------------------------------
  ! Make organic matter density [organic] [forganic]
  ! -----------------------------------
  if (fsurdat /= ' ') then
     call mkorganic( mksrf_forganic_mesh, mksrf_forganic, mesh_model, pctlnd_pft, lat, pioid, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkorganic')
  end if

  ! -----------------------------------
  ! Make VOC emission factors for isoprene [ef1_btr,ef1_fet,ef1_fdt,ef1_shr,ef1_grs,ef1_crp]
  ! -----------------------------------
  if (fsurdat /= ' ')  then
     call mkvocef ( mksrf_fvocef_mesh, mksrf_fvocef, mesh_model, pioid, lat, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkvocef')
  end if

  ! -----------------------------------
  ! Adjust pctlak, pctwet, pcturb and pctgla
  ! -----------------------------------
  do n = 1,lsize_o

     ! If have pole points on grid - set south pole to glacier
     ! north pole is assumed as non-land
     if (abs((lat(n) - 90._r8)) < 1.e-6_r8) then
        pctlak(n) = 0._r8
        pctwet(n) = 0._r8
        pcturb(n) = 0._r8
        pctgla(n) = 100._r8
     end if

     ! Truncate all percentage fields on output grid. This is needed to insure that
     ! wt is zero (not a very small number such as 1e-16) where it should be zero
     pctlak(n) = float(nint(pctlak(n)))
     pctwet(n) = float(nint(pctwet(n)))
     pctgla(n) = float(nint(pctgla(n)))

     ! Assume wetland, glacier and/or lake when dataset landmask implies ocean
     if (pctlnd_pft(n) < 1.e-6_r8) then
        if (pctgla(n) < 1.e-6_r8) then
           pctwet(n) = 100._r8 - pctlak(n)
           pctgla(n) = 0._r8
        else
           pctwet(n)  = max(100._r8 - pctgla(n) - pctlak(n), 0.0_r8)
        end if
        pcturb(n) = 0._r8
     end if

     ! Make sure sum of land cover types does not exceed 100. If it does,
     ! subtract excess from most dominant land cover.
     suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n)
     if (suma > 250._r4) then
        write (6,*) subname, ' error: sum of pctlak, pctwet,', &
             'pcturb and pctgla is greater than 250%'
        write (6,*)'n,pctlak,pctwet,pcturb,pctgla= ', &
             n,pctlak(n),pctwet(n),pcturb(n),pctgla(n)
        call shr_sys_abort()
     else if (suma > 100._r4) then
        pctlak(n) = pctlak(n) * 100._r8/suma
        pctwet(n) = pctwet(n) * 100._r8/suma
        pcturb(n) = pcturb(n) * 100._r8/suma
        pctgla(n) = pctgla(n) * 100._r8/suma
     end if
  end do

  ! -----------------------------------
  ! Perform other normalizations
  ! -----------------------------------

  ! Normalize land use and make sure things add up to 100% as well as
  ! checking that things are as they should be.
  call normalize_and_check_landuse(lsize_o)

  ! Write out sum of PFT's
  do k = natpft_lb,natpft_ub
     suma = 0._r8
     do n = 1,lsize_o
        suma = suma + pctnatpft(n)%get_one_pct_p2g(k)
     enddo
     ! TODO: calculate global sum here
     if (root_task) then
        write(ndiag,*) 'sum over domain of pft ',k,suma
     end if
  enddo
  if (root_task) write(ndiag,*)
  do k = cft_lb,cft_ub
     suma = 0._r8
     do n = 1,lsize_o
        suma = suma + pctcft(n)%get_one_pct_p2g(k)
     enddo
     ! TODO: calculate global sum here
     if (root_task) then
        write(6,*) 'sum over domain of cft ',k,suma
     end if
  enddo
  if (root_task) write(ndiag,*)

  ! Make final values of percent urban by class and compute urban parameters
  ! This call needs to occur after all corrections are made to pcturb
  allocate (urban_classes_g(lsize_o,numurbl)); urban_classes_g(:,:) = spval
  call normalize_classes_by_gcell(urban_classes, pcturb, urban_classes_g)
  if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out percnt urban"

  ! Make Urban Parameters from raw input data and write to surface dataset
  ! Write to netcdf file is done inside mkurbanpar routine
  call mkurbanpar(mksrf_furban, pioid, mesh_model, urban_region, urban_classes_g, &
       urban_skip_abort_on_invalid_data_check)
  deallocate(urban_region)

  ! -----------------------------------
  ! Write out PCT_URBAN, PCT_GLACIER, PCT_LAKE and PCT_WETLAND and
  ! PCT_NATVEG, PCT_NAT_PFT, PCT_CROP and PCT_CFT
  ! -----------------------------------

  if (fsurdat /= ' ') then
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out PCT_URBAN"
     call mkfile_output(pioid,  mesh_model,  'PCT_URBAN', urban_classes_g, lev1name='numurbl', rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_URBAN')

     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out PCT_GLACIER"
     call mkfile_output(pioid, mesh_model, 'PCT_GLACIER', pctgla, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in mkfile_output for pctgla')

     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out PCT_LAKE"
     call mkfile_output(pioid,  mesh_model,  'PCT_LAKE', pctlak, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in mkfile_output for pctlak')

     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out PCT_WETLAND"
     call mkfile_output(pioid, mesh_model,  'PCT_WETLAND', pctwet, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in in mkfile_output for pctwet')

     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing PCT_NATVEG"
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing PCT_NATVEG"
     allocate(pctnatveg(lsize_o))
     call get_pct_l2g_array(pctnatpft, pctnatveg)
     call mkfile_output(pioid, mesh_model, 'PCT_NATVEG', pctnatveg, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_NATVEG')

     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing PCT_CROP"
     allocate(pctcrop(lsize_o))
     call get_pct_l2g_array(pctcft, pctcrop)
     call mkfile_output(pioid, mesh_model, 'PCT_CROP', pctcrop, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_CROP')

     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing PCT_NAT_PFT"
     allocate(pct_nat_pft(lsize_o, 0:num_natpft))
     call get_pct_p2l_array(pctnatpft, ndim1=lsize_o, ndim2=num_natpft+1, pct_p2l=pct_nat_pft)
     call mkfile_output(pioid, mesh_model, 'PCT_NAT_PFT', pct_nat_pft, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_NAT_PFT')

     if (num_cft > 0) then
        if (root_task)  write(ndiag, '(a)') trim(subname)//" writing PCT_CFT"
        allocate(pct_cft(lsize_o, num_cft))
        call get_pct_p2l_array(pctcft, ndim1=lsize_o, ndim2=num_cft, pct_p2l=pct_cft)
        call mkfile_output(pioid, mesh_model, 'PCT_CFT', pct_cft, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_CFT')
     end if
  end if

  ! ----------------------------------------------------------------------
  ! Make glacier multiple elevation classes [pctglcmec,topoglcmec] from [fglacier] dataset
  ! ----------------------------------------------------------------------
  ! This call needs to occur after pctgla has been adjusted for the final time
  if (fsurdat /= ' ') then
     call mkglcmecInit (pioid)
     call mkglcmec(mksrf_fglacier_mesh, mksrf_fglacier, mesh_model, pioid, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkglcmec')
  end if

  ! ----------------------------------------------------------------------
  ! Close surface dataset
  ! ----------------------------------------------------------------------

  call pio_closefile(pioid)

  if (root_task) then
     write(ndiag,*)
     write(ndiag,'(a)') 'Successfully created surface data output file = '//trim(fsurdat)
     write(ndiag,'(a)') '   This file contains the land model surface data'
     write(ndiag,*)
  end if

  ! ======================================================================
  ! Create fdyndat if appropriate
  ! ======================================================================

1000 continue

  if (mksrf_fdynuse /= ' ') then

     if (fdyndat == ' ') then
        call shr_sys_abort(' must specify fdyndat in namelist if mksrf_fdynuse is not blank')
     end if

     if (root_task) then
        write(ndiag,*)
        write(ndiag,'(1x,80a1)') ('=',k=1,80)
        write(ndiag,'(1x,80a1)') ('*',k=1,80)
        write(ndiag,'(1x,80a1)') ('=',k=1,80)
        write(ndiag,*)
        write(ndiag,'(a)')'Creating dynamic land use dataset '//trim(fdyndat)
     end if

     allocate(pctcft_max(lsize_o))    ;
     allocate(pctnatpft_max(lsize_o)) ;
     allocate(pctlnd_pft_dyn(lsize_o))

     ! open output file
     call mkpio_wopen(trim(fdyndat), clobber=.true., pioid=pioid)

     ! Define dimensions
     call mkfile_define_dims(pioid, nx=mksrf_fgrid_mesh_nx, ny=mksrf_fgrid_mesh_ny, dynlanduse=.true.)

     ! Define global attributes
     call mkfile_define_atts(pioid, dynlanduse = .true.)

     ! Define variables
     call mkfile_define_vars(pioid, dynlanduse = .true.)

     ! End define model
     rcode = pio_enddef(pioid)

     ! Write out natpft
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out natpft"
     rcode = pio_inq_varid(pioid, 'natpft', pio_varid)
     rcode = pio_put_var(pioid, pio_varid, (/(n,n=natpft_lb,natpft_ub)/))

     ! Write out cft
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out cft"
     rcode = pio_inq_varid(pioid, 'cft', pio_varid)
     rcode = pio_put_var(pioid, pio_varid, (/(n,n=cft_lb,cft_ub)/))

     ! Write out model grid
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out LONGXY"
     call mkfile_output(pioid, mesh_model, 'LONGXY', lon, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out LATIXY"
     call mkfile_output(pioid, mesh_model, 'LATIXY', lat, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')

     ! Write out natpft
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out natpft"
     rcode = pio_inq_varid(pioid, 'natpft', pio_varid)
     rcode = pio_put_var(pioid, pio_varid, (/(n,n=natpft_lb,natpft_ub)/))

     ! Write out cft
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out cft"
     rcode = pio_inq_varid(pioid, 'cft', pio_varid)
     rcode = pio_put_var(pioid, pio_varid, (/(n,n=cft_lb,cft_ub)/))

     ! Write out PFTDATA_MASK
     ! pftdata_mask was calculated ABOVE
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing land mask (calculated in furdata calc)"
     call mkfile_output(pioid,  mesh_model, 'PFTDATA_MASK', pftdata_mask, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')

     ! Write out LANDFRAC_PFT
     ! landfrac_pft was calculated ABOVE
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing land fraction calculated in fsurdata calc)"
     call mkfile_output(pioid,  mesh_model, 'LANDFRAC_PFT', landfrac_pft, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')

     ! -----------------------------------------
     ! Read in each dynamic pft landuse dataset
     ! -----------------------------------------

     ! Open txt file
     if (root_task) then
        write(ndiag,'(a)')' Opening '//trim(mksrf_fdynuse)//' to read dynamic data forcing '
        open (newunit=nfdyn, file=trim(mksrf_fdynuse), form='formatted', iostat=ier)
        if (ier /= 0) then
           call shr_sys_abort(subname//" failed to open file "//trim(mksrf_fdynuse))
        end if
     end if

     pctnatpft_max = pctnatpft
     pctcft_max = pctcft

     end_of_fdynloop = .false.
     ntim = 0
     do

        ! Determine file name - if there are no more files than exit before broadcasting
        if (root_task) then
           read(nfdyn, '(A195,1x,I4)', iostat=ier) string, year
           if (ier /= 0) end_of_fdynloop = .true.
        end if
        call mpi_bcast(end_of_fdynloop, 1, MPI_LOGICAL, 0, mpicom, ier)
        if (end_of_fdynloop) then
           EXIT
        end if
        call mpi_bcast (string, len(string), MPI_CHARACTER, 0, mpicom, ier)
        call mpi_bcast (year, 1, MPI_INTEGER, 0, mpicom, ier)

        ! Intrepret string as a filename with PFT and harvesting values in it

        fname = string
        if (root_task) then
           read(nfdyn, '(A195,1x,I4)', iostat=ier) fhrvname, year2
           write(ndiag,'(a,i8,a)')' input pft dynamic dataset for year ', year,' is : '//trim(fname)
        end if
        call mpi_bcast (fhrvname, len(fhrvname), MPI_CHARACTER, 0, mpicom, ier)
        call mpi_bcast (year2, 1, MPI_INTEGER, 0, mpicom, ier)
        if ( year2 /= year ) then
           if (root_task) then
              write(ndiag,*) subname, ' error: year for harvest not equal to year for PFT files'
           end if
           call shr_sys_abort()
        end if
        ntim = ntim + 1
        if (root_task) then
           write(ndiag,'(a,i8)')subname//' ntime = ',ntim
        end if

        rcode = pio_inq_varid(pioid, 'YEAR', pio_varid)
        rcode = pio_put_var(pioid, pio_varid, (/ntim/), year)
        rcode = pio_inq_varid(pioid, 'time', pio_varid)
        rcode = pio_put_var(pioid, pio_varid, (/ntim/), year)
        ! rcode = pio_inq_varid(pioid, 'input_pftdata_filename', pio_varid)
        ! rcode = pio_put_var(pioid, pio_varid, (/1,ntim/), (/len_trim(string),1/), trim(string))
        call pio_syncfile(pioid)

        ! Create pctpft data at model resolution from file fname
        ! Note that pctlnd_o below is different than the above call and returns pctlnd_pft_dyn

        call mkpft( mksrf_fvegtyp_mesh, fname, mesh_model, &
             pctlnd_o=pctlnd_pft_dyn, pctnatpft_o=pctnatpft, pctcft_o=pctcft, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkpft')
        call pio_syncfile(pioid)

        ! Consistency check on input land fraction
        ! pctlnd_pft was calculated ABOVE
        do n = 1,lsize_o
           if (pctlnd_pft_dyn(n) /= pctlnd_pft(n)) then
              if (root_task) then
                 write(ndiag,*) subname,' error: pctlnd_pft for dynamics data = ',&
                      pctlnd_pft_dyn(n), ' not equal to pctlnd_pft for surface data = ',&
                      pctlnd_pft(n),' at n= ',n
                 if ( trim(fname) == ' ' )then
                    write(ndiag,*) ' PFT string = ',trim(string)
                 else
                    write(ndiag,*) ' PFT file = ', fname
                 end if
              end if
              call shr_sys_abort()
           end if
        end do

        ! Create harvesting data at model resolution
        ! Output data is written in mkharvest
        call mkharvest( mksrf_fhrvtyp_mesh, fhrvname, mesh_model, &
             pioid_o=pioid, constant=.false., ntime=ntim, rc=rc )
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkharvest')
        call pio_syncfile(pioid)

        ! Do landuse changes such as for the poles, etc.
        ! If have pole points on grid - set south pole to glacier
        ! north pole is assumed as non-land
        ! pctlak, pctwet, pcturb and pctgla were calculated ABOVE
        ! pctnatpft and pctcft were calculated ABOVE
        do n = 1,lsize_o
           if (abs(lat(n) - 90._r8) < 1.e-6_r8) then
              pctlak(n) = 0._r8
              pctwet(n) = 0._r8
              pcturb(n) = 0._r8
              pctgla(n) = 100._r8
              call pctnatpft(n)%set_pct_l2g(0._r8)
              call pctcft(n)%set_pct_l2g(0._r8)
           end if
        end do

        ! Normalize land use and make sure things add up to 100% as well as
        ! checking that things are as they should be.
        if (root_task) then
           write(ndiag,*)
           write(ndiag,'(1x,80a1)') ('=',k=1,80)
           write(ndiag,'(a)')' calling normalize_and_check_landuse'
        end if
        call normalize_and_check_landuse(lsize_o)

        ! Given an array of pct_pft_type variables, update all the max_p2l variables.
        call update_max_array(pctnatpft_max, pctnatpft)
        call update_max_array(pctcft_max, pctcft)

        if (root_task)  write(ndiag, '(a,i8)') trim(subname)//" writing PCT_NAT_PFT for year ",year
        rcode = pio_inq_varid(pioid, 'PCT_NAT_PFT', pio_varid)
        call pio_setframe(pioid, pio_varid, int(ntim, kind=Pio_Offset_Kind))
        call get_pct_p2l_array(pctnatpft, ndim1=lsize_o, ndim2=num_natpft+1, pct_p2l=pct_nat_pft)
        call mkfile_output(pioid, mesh_model, 'PCT_NAT_PFT', pct_nat_pft, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_NAT_PFT')
        call pio_syncfile(pioid)

        if (root_task)  write(ndiag, '(a,i8)') trim(subname)//" writing PCT_CROP for year ",year
        rcode = pio_inq_varid(pioid, 'PCT_CROP', pio_varid)
        call pio_setframe(pioid, pio_varid, int(ntim, kind=Pio_Offset_Kind))
        call get_pct_l2g_array(pctcft, pctcrop)
        call mkfile_output(pioid, mesh_model, 'PCT_CROP', pctcrop, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_CROP')
        call pio_syncfile(pioid)

        if (num_cft > 0) then
           if (root_task)  write(ndiag, '(a,i8)') trim(subname)//" writing PCT_CFT for year ",year
           rcode = pio_inq_varid(pioid, 'PCT_CFT', pio_varid)
           call pio_setframe(pioid, pio_varid, int(ntim, kind=Pio_Offset_Kind))
           call get_pct_p2l_array(pctcft, ndim1=lsize_o, ndim2=num_cft, pct_p2l=pct_cft)
           call mkfile_output(pioid, mesh_model, 'PCT_CFT', pct_cft, rc=rc)
           if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_CFT')
           call pio_syncfile(pioid)
        end if

        if (root_task) then
           write(ndiag,'(1x,80a1)') ('=',k=1,80)
           write(ndiag,*)
        end if

     end do   ! end of read loop

     if (root_task)  write(ndiag, '(a,i8)') trim(subname)//" writing PCT_NAT_PFT_MAX "
     call get_pct_p2l_array(pctnatpft_max, ndim1=lsize_o, ndim2=num_natpft+1, pct_p2l=pct_nat_pft)
     call mkfile_output(pioid, mesh_model, 'PCT_NAT_PFT_MAX', pct_nat_pft, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_NAT_PFT')

     if (root_task)  write(ndiag, '(a,i8)') trim(subname)//" writing PCT_CROP_MAX"
     call get_pct_l2g_array(pctcft_max, pctcrop)
     call mkfile_output(pioid, mesh_model, 'PCT_CROP_MAX', pctcrop, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_CROP')

     if (num_cft > 0) then
        if (root_task)  write(ndiag, '(a,i8)') trim(subname)//" writing PCT_CFT_MAX"
        call get_pct_p2l_array(pctcft_max, ndim1=lsize_o, ndim2=num_cft, pct_p2l=pct_cft)
        call mkfile_output(pioid, mesh_model, 'PCT_CFT_MAX', pct_cft, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_CFT')
     end if

     ! Close the file
     call pio_closefile(pioid)

  end if   ! end of if-create dynamic landust dataset

  ! -----------------------------------
  ! Wrap things up
  ! -----------------------------------
  close (ndiag)

  ! TODO: fix the error condition here
  call ESMF_Finalize()

  !-----------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------

    subroutine normalize_and_check_landuse(ns_o)
      !
      ! Normalize land use and make sure things add up to 100% as well as
      ! checking that things are as they should be.
      !
      ! Precondition: pctlak + pctwet + pcturb + pctgla <= 100 (within roundoff)
      !
      use mkpftConstantsMod , only : baregroundindex
      use mkpftUtilsMod     , only : adjust_total_veg_area

      ! input/output variables:
      integer, intent(in) :: ns_o

      ! local variables:
      integer  :: k,n                         ! indices
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

      do n = 1,ns_o

         ! Check preconditions
         if ( pctlak(n) < 0.0_r8 )then
            write(6,*) subname, ' ERROR: pctlak is negative!'
            write(6,*) 'n, pctlak = ', n, pctlak(n)
            call shr_sys_abort()
         end if
         if ( pctwet(n) < 0.0_r8 )then
            write(6,*) subname, ' ERROR: pctwet is negative!'
            write(6,*) 'n, pctwet = ', n, pctwet(n)
            call shr_sys_abort()
         end if
         if ( pcturb(n) < 0.0_r8 )then
            write(6,*) subname, ' ERROR: pcturb is negative!'
            write(6,*) 'n, pcturb = ', n, pcturb(n)
            call shr_sys_abort()
         end if
         if ( pctgla(n) < 0.0_r8 )then
            write(6,*) subname, ' ERROR: pctgla is negative!'
            write(6,*) 'n, pctgla = ', n, pctgla(n)
            call shr_sys_abort()
         end if

         suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n)
         if (suma > (100._r8 + tol_loose)) then
            write(6,*) subname, ' ERROR: pctlak + pctwet + pcturb + pctgla must be'
            write(6,*) '<= 100% before calling this subroutine'
            write(6,*) 'n, pctlak, pctwet, pcturb, pctgla = ', &
                 n, pctlak(n), pctwet(n), pcturb(n), pctgla(n)
            call shr_sys_abort()
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
            write(6,*) subname, ' ERROR in rescaling veg based on (special excluding urban)'
            write(6,*) 'suma = ', suma
            call shr_sys_abort()
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
                     call shr_sys_abort()
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
            call shr_sys_abort()
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
               call shr_sys_abort()
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
            call shr_sys_abort()
         end if

         suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n)
         if (suma < 100._r8-epsilon(suma) .and. suma > (100._r8 - 4._r8*epsilon(suma))) then
            write (6,*) subname, 'n,pctlak,pctwet,pcturb,pctgla,pctnatveg,pctcrop= ', &
                 n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),&
                 pctnatpft(n)%get_pct_l2g(), pctcft(n)%get_pct_l2g()
            call shr_sys_abort()
         end if

         suma = suma + pctnatpft(n)%get_pct_l2g() + pctcft(n)%get_pct_l2g()
         if ( abs(suma-100._r8) > 1.e-10_r8) then
            write (6,*) subname, ' error: sum of pctlak, pctwet,', &
                 'pcturb, pctgla, pctnatveg and pctcrop is NOT equal to 100'
            write (6,*)'n,pctlak,pctwet,pcturb,pctgla,pctnatveg,pctcrop,sum= ', &
                 n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),&
                 pctnatpft(n)%get_pct_l2g(),pctcft(n)%get_pct_l2g(), suma
            call shr_sys_abort()
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
               call shr_sys_abort()
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
               call shr_sys_abort()
            end if
         end do
      end if

      ! Make sure that there is no vegetation outside the pft mask
      do n = 1,ns_o
         if (pftdata_mask(n) == 0 .and. (pctnatpft(n)%get_pct_l2g() > 0 .or. pctcft(n)%get_pct_l2g() > 0)) then
            write (6,*)'vegetation found outside the pft mask at n=',n
            write (6,*)'pctnatveg,pctcrop=', pctnatpft(n)%get_pct_l2g(), pctcft(n)%get_pct_l2g()
            call shr_sys_abort()
         end if
      end do

      ! Make sure that sums at the landunit level all add to 100%
      ! (Note that we don't check pctglcmec here, because it isn't computed at the point
      ! that this subroutine is called -- but the check of sum(pctglcmec) is done in
      ! mkglcmecMod)
      ! (Also note that we don't need to check pctnatpft or pctcft, because a similar check
      ! is done internally by the pct_pft_type routines.)
      do n = 1,ns_o
         if (abs(sum(urban_classes(n,:)) - 100._r8) > 1.e-12_r8) then
            write(6,*) 'sum(urban_classes(n,:)) != 100: ', n, sum(urban_classes(n,:))
            call shr_sys_abort()
         end if
      end do

      if (root_task) then
         if ( nsmall_tot > 0 )then
            write(ndiag,*)'number of small pft = ', nsmall_tot
         end if
      end if

    end subroutine normalize_and_check_landuse

end program mksurfdata
