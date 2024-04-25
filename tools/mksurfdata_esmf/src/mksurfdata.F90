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
  !    mksrf_fpctlak             - Percent lake dataset
  !    mksrf_fpctlak_mesh        - Mesh for mksrf_fpctlak
  !    mksrf_flakdep             - Lake depth dataset
  !    mksrf_flakdep_mesh        - Mesh for mksrf_flakdep
  !    mksrf_fwetlnd             - Wetland water dataset
  !    mksrf_fwetlnd_mesh        - Mesh for mksrf_fwetlnd
  !    mksrf_fmax                - Max fractional saturated area dataset
  !    mksrf_fmax_mesh           - Mesh for mksrf_fmax
  !    mksrf_fsoicol             - Soil color dataset
  !    mksrf_fsoicol_mesh        - Mesh for mksrf_fsoicol
  !    mksrf_fsoitex             - Soil texture dataset in mapunits
  !    mksrf_fsoitex_lookup      - Soil texture lookup for converting mapunits to sand/silt/clay and organic carbon content
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
  !    mksrf_ftopostats_override - Use this file to read in STD_ELEV and SLOPE
  !    mksrf_fvic                - VIC parameters dataset
  !    mksrf_fvic_mesh           - Mesh for mksrf_fvic
  ! ======================================
  ! Optionally specify setting for:
  ! ======================================
  !    mksrf_fdynuse ----- ASCII text file that lists each year of pft, urban, and lake files to use
  !    mksrf_gridtype ---- Type of grid (default is 'global')
  !    outnc_double ------ If output should be in double precision
  !    outnc_large_files - If output should be in NetCDF large file format
  !    outnc_vic --------- Output fields needed for VIC
  !    outnc_3dglc ------- Output 3D glacier fields (normally only needed for comparasion)
  !    nglcec ------------ If you want to change the number of Glacier elevation classes
  !    gitdescribe ------- Description of this version from git
  !    numpft ------------ Iif different than default of 16
  !    urban_skip_abort_on_invalid_data_check--- work around urban bug
  !    no_inlandwet ------ If wetland should be set to 0% over land
  ! ======================================
  ! Note: the folloiwng Optional settings have been REMOVED -
  !  instead should now use tools subset_data and modify_fsurdat
  ! ======================================
  !    all_veg ----------- If entire area is to be vegetated (pft_idx and pft_frc then required)
  !    all_urban --------- If entire area is urban
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
  use mkurbanparMod      , only : mkurbanInit, mkurban, mkurbanpar, mkurban_topo, numurbl, update_max_array_urban
  use mklanwatMod        , only : mkpctlak, mklakdep, mkwetlnd, update_max_array_lake
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
  integer                         :: lsize_o                 ! Size of the local mesh elments
  integer                         :: node_count              ! Number of gridcells on the mesh
  integer                         :: local_nodes(1)          ! Local gridcells on the mesh
  integer                         :: total_nodes(1)          ! Total gridcells on the mesh

  ! error status
  integer                         :: ier,rcode               ! error status

  ! dynamic land use
  integer                         :: nfdyn, nfpio            ! unit numbers
  integer                         :: ntim                    ! time sample for dynamic land use
  integer                         :: year                    ! year for dynamic land use
  integer                         :: year2                   ! year for dynamic land use for harvest file
  real(r8)                        :: suma                    ! local sum for error check
  real(r8)                        :: loc_suma, glob_suma     ! local and global sum for error check with mpi_allreduce

  ! model grid
  real(r8), allocatable           :: lon(:)
  real(r8), allocatable           :: lat(:)

  ! pct vegetation data
  real(r8), allocatable           :: landfrac_pft(:)         ! PFT data: % land per gridcell
  real(r8), allocatable           :: pctlnd_pft(:)           ! PFT data: % of gridcell for PFTs
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
  real(r8), allocatable           :: pctlak_max(:)           ! maximum percent of grid cell that is lake
  real(r8), allocatable           :: pctwet(:)               ! percent of grid cell that is wetland
  real(r8), allocatable           :: pctocn(:)               ! percent of grid cell that is ocean
  real(r8), allocatable           :: pctgla(:)               ! percent of grid cell that is glacier
  integer , allocatable           :: urban_region(:)         ! urban region ID
  real(r8), allocatable           :: pcturb(:)               ! percent of grid cell that is urbanized (total across all urban classes)
  real(r8), allocatable           :: pcturb_max(:,:)         ! maximum percent cover of each urban class, as % of grid cell
  real(r8), allocatable           :: urban_classes(:,:)      ! percent cover of each urban class, as % of total urban area
  real(r8), allocatable           :: urban_classes_g(:,:)    ! percent cover of each urban class, as % of grid cell
  real(r8), allocatable           :: elev(:)                 ! glc elevation (m)
  real(r8), allocatable           :: pctwet_orig(:)          ! percent wetland of gridcell before dynamic land use adjustments
  real(r8), allocatable           :: pctgla_orig(:)          ! percent glacier of gridcell before dynamic land use adjustments

  ! other variables written to file
  real(r8), allocatable           :: landfrac_mksurfdata(:)  ! land fraction used for renormalization of areas

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

  ! character variables
  character(len=CL)               :: default_log_suffix      ! default log file suffix to use for ESMF PET files
  character(len=CL)               :: string                  ! string read in
  character(len=CL)               :: fname
  character(len=*), parameter     :: subname = 'mksrfdata'   ! program name

  character(len=*), parameter :: u_FILE_u = &
       __FILE__

  ! ======================================================================
  ! Initialize MPI get the rank and determine root task
  ! ======================================================================

  call MPI_init(rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
  mpicom = mpi_comm_world

  ! Determine root task
  call MPI_comm_rank(mpicom, iam, rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
  root_task = (iam == 0)
  call MPI_comm_size(mpicom, npes, rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()


  ! ======================================================================
  ! Read in namelist before initializing MPI or ESMF
  ! ======================================================================
  call read_namelist_input()

  ! ======================================================================
  ! Initialize ESMF and get mpicom from ESMF
  ! ======================================================================

  if (create_esmf_pet_files) then
     logkindflag = ESMF_LOGKIND_MULTI
  else
     logkindflag = ESMF_LOGKIND_MULTI_ON_ERROR
  end if
  default_log_suffix = trim(mksrf_grid_name) // '_ESMF_LogFile'
  call ESMF_Initialize(mpiCommunicator=MPICOM, logkindflag=logkindflag, logappendflag=.false., &
       defaultDefaultLogFilename=trim(default_log_suffix), ioUnitLBound=5001, ioUnitUBound=5101, vm=vm, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
  call ESMF_VMGetGlobal(vm, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
  call ESMF_VMGet(vm, mpicommunicator=mpicom, localPet=iam, petcount=petcount, &
       ssiLocalPetCount=stride, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
  call ESMF_LogSet(flush=.true.)
  call ESMF_LogWrite("mksurfdata starting", ESMF_LOGMSG_INFO)
  !
  ! open output ndiag file
  !
  if (root_task) then
     open (newunit=ndiag, file=trim(fsurlog), iostat=ier)
     if (ier /= 0) then
        call shr_sys_abort(' failed to open ndiag file '//trim(fsurlog))
     end if
     write (ndiag,'(a)') 'Attempting to create surface boundary data .....'
     write (ndiag,'(72a1)') ("-",n=1,60)
     flush(ndiag)
  else
     ndiag = 6
  end if
  !
  ! Finish handling of the namelist control variables
  !
  ! Broadcast namelist to all pes
  ! root_task is a module variable in mkvarctl
  call bcast_namelist_input()

  ! Write out namelist input to ndiag
  call check_namelist_input()
  call write_namelist_input()

  ! Some checking
   if (root_task) then
     write(ndiag,'(2(a,I))') ' npes = ', npes, ' grid size = ', grid_size
     flush(ndiag)
  end if
  if (petcount >  grid_size ) then
     call shr_sys_abort(' ERROR: number of tasks exceeds the size of the grid' )
  end if

  ! ======================================================================
  ! Initialize PIO
  ! ======================================================================

  ! the following returns pio_iosystem
  call pio_init(iam, mpicom, max(1,petcount/stride), 0, stride, PIO_REARR_SUBSET, pio_iosystem)

  ! Open txt file
  if (root_task) then
     write(ndiag,*)' Opening file and reading pio_iotype from txt file with the same name'
     flush(ndiag)
     open (newunit=nfpio, file='pio_iotype.txt', status='old', &
           form='formatted', action='read', iostat=ier)
     if (ier /= 0) then
        call shr_sys_abort(subname//" failed to open file pio_iotype.txt")
     end if
     read(nfpio,*)  ! skip file header
     read(nfpio, '(i)', iostat=ier) pio_iotype
     if (ier /= 0) then
        call shr_sys_abort(subname//" failed to read file pio_iotype.txt")
     end if
  end if
  call mpi_bcast(pio_iotype, 1, MPI_INTEGER, 0, mpicom, ier)

  pio_ioformat = PIO_64BIT_DATA

  call ESMF_LogWrite("finished initializing PIO", ESMF_LOGMSG_INFO)

  if (fsurlog == ' ') then
     call shr_sys_abort(' ERROR: must specify fsurlog in namelist')
  end if

  ! ======================================================================
  ! Create fsurdat
  ! ======================================================================

  ! Read in model mesh to determine the number of local points
  call ESMF_LogWrite("MESH creation (if this takes too long [more than an hour] and hangs, you may need more memory...)", ESMF_LOGMSG_INFO)
  mesh_model = ESMF_MeshCreate(filename=trim(mksrf_fgrid_mesh), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

  ! Get the number of local destination points on my processor (lsize_o)
  call ESMF_MeshGet(mesh_model, numOwnedElements=lsize_o, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
  local_nodes(1) = lsize_o
  call ESMF_VMAllReduce(vm, local_nodes, total_nodes, count=1, reduceflag=ESMF_REDUCE_SUM, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
  node_count = total_nodes(1)
  if (node_count /=  grid_size) then
     if (root_task) then
        write (ndiag,'(a, I, a, I)') ' node_count = ', node_count, ' grid_size = ', grid_size
        flush(ndiag)
     end if
     call shr_sys_abort(' ERROR: size of input mesh file does not agree with expected size of nx*ny' )
  end if
  ! Initialize urban dimensions (needed to initialize the dimensions in fsurdat)
  call ESMF_LogWrite("mkurbanInit...")
  call mkurbanInit(mksrf_furban)

  ! Initialize pft/cft dimensions (needed to initialize the dimensions in fsurdat)
  call ESMF_LogWrite("mkpftInit...")
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
        flush(ndiag)
     end if

     call ESMF_LogWrite("mkfile...")
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
  ! Write out natpft, cft, and time
  ! -----------------------------------
  if (fsurdat /= ' ') then
     rcode = pio_inq_varid(pioid, 'natpft', pio_varid)
     rcode = pio_put_var(pioid, pio_varid, (/(n,n=natpft_lb,natpft_ub)/))
     if (num_cft > 0) then
        rcode = pio_inq_varid(pioid, 'cft', pio_varid)
        rcode = pio_put_var(pioid, pio_varid, (/(n,n=cft_lb,cft_ub)/))
     end if
  end if

  ! -----------------------------------
  ! Make lats/lons of model
  ! -----------------------------------
  allocate (lon(lsize_o)) ; lon(:) = spval
  allocate (lat(lsize_o)) ; lat(:) = spval
  call mkdomain(mesh_model, lon_o=lon, lat_o=lat, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkdomain')
  if (fsurdat /= ' ') then
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out model grid"
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out LONGXY"
     if (root_task)  flush(ndiag)
     call mkfile_output(pioid, mesh_model, 'LONGXY', lon, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for LONGXY')
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out LATIXY"
     if (root_task)  flush(ndiag)
     call mkfile_output(pioid, mesh_model, 'LATIXY', lat, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for LATIXY')
     call pio_syncfile(pioid)
  end if

  if (root_task)then
     write(ndiag,*)' Initialization is complete, going on to process types of input files'
     flush(ndiag)
  end if

  ! -----------------------------------
  ! Make LAI and SAI from 1/2 degree data and write to surface dataset
  ! Write to netcdf file is done inside mklai routine
  ! -----------------------------------
  if (fsurdat /= ' ') then
     call mklai(mksrf_flai_mesh, mksrf_flai, mesh_model, pioid, rc=rc)
     flush(ndiag)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mklai')
  end if

  ! -----------------------------------
  ! Make PFTs [pctnatpft, pctcft] from dataset [fvegtyp]
  ! Make landfrac_pft
  ! -----------------------------------
  ! Determine fractional land from pft dataset
  allocate(pctlnd_pft(lsize_o)); pctlnd_pft(:) = spval
  allocate(pctnatpft(lsize_o)) ;
  allocate(pctcft(lsize_o))    ;
  allocate(landfrac_pft(lsize_o))  ; landfrac_pft(:) = spval
  call mkpft( mksrf_fvegtyp_mesh, mksrf_fvegtyp, mesh_model, &
       pctlnd_o=pctlnd_pft, pctnatpft_o=pctnatpft, pctcft_o=pctcft, rc=rc)
  flush(ndiag)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkdomain')

  do n = 1,lsize_o
     landfrac_pft(n) = pctlnd_pft(n)/100._r8
  end do
  if (fsurdat /= ' ') then
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing land fraction  from pft dataset"
     if (root_task)  flush(ndiag)
     call mkfile_output(pioid,  mesh_model, 'LANDFRAC_PFT', landfrac_pft, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
     call pio_syncfile(pioid)
  end if

  ! -----------------------------------
  ! Make constant harvesting data at model resolution
  ! -----------------------------------
  ! Note that this call must come after call to mkpftInit - since num_cft is set there
  ! Output data is written in mkharvest
  if (fsurdat /= ' ') then
     call mkharvest( mksrf_fhrvtyp_mesh, mksrf_fhrvtyp, mesh_model, pioid, &
                     rc=rc )
     flush(ndiag)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkharvest_init')
  end if

  ! -----------------------------------
  ! Make inland water [pctlak, pctwet]
  ! -----------------------------------
  ! LAKEDEPTH is written out in the subroutine
  ! Need to keep pctlak and pctwet external for use below
  allocate ( pctlak(lsize_o)) ; pctlak(:) = spval
  allocate ( pctlak_max(lsize_o)) ; pctlak_max(:) = spval
  call mkpctlak(mksrf_fpctlak_mesh, mksrf_fpctlak, mesh_model, pctlak, pioid, &
                rc=rc)
  flush(ndiag)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkpctlak')
  call mklakdep(mksrf_flakdep_mesh, mksrf_flakdep, mesh_model, pioid, fsurdat, &
                rc=rc)
  flush(ndiag)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mklakdep')

  allocate ( pctwet(lsize_o)) ; pctwet(:) = spval
  allocate ( pctwet_orig(lsize_o)) ; pctwet_orig(:) = spval
  call mkwetlnd(mksrf_fwetlnd_mesh, mksrf_fwetlnd, mesh_model, pctwet, rc=rc)
  flush(ndiag)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkwetlnd')

  ! Initialize pctocn to zero.
  ! Until ctsm5.1 we set pctwet = 100 at ocean points rather than
  ! setting a pctocn. Starting with ctsm5.2, we set pctocn = 100 at
  ! ocean points in subroutine normalize_and_check_landuse.
  ! No regridding required.
  allocate ( pctocn(lsize_o)); pctocn(:) = 0._r8

  ! -----------------------------------
  ! Make glacier fraction [pctgla] from [fglacier] dataset
  ! -----------------------------------
  allocate (pctgla(lsize_o)) ; pctgla(:) = spval
  allocate (pctgla_orig(lsize_o)) ; pctgla_orig(:) = spval
  call mkglacier (mksrf_fglacier_mesh, mksrf_fglacier, mesh_model, glac_o=pctgla, rc=rc)
  flush(ndiag)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkglacier')

  ! -----------------------------------
  ! Make glacier region ID [glacier_region] from [fglacierregion] dataset
  ! -----------------------------------
  if (fsurdat /= ' ') then
     ! GLACIER_REGION is written out in the subroutine
     call mkglacierregion(mksrf_fglacierregion_mesh, mksrf_fglacierregion, mesh_model, pioid, rc)
     flush(ndiag)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkglacierregion')
  end if

  ! -----------------------------------
  ! Make soil texture and organic carbon content [pctsand, pctclay, organic]
  ! -----------------------------------
  if (fsurdat /= ' ') then
     call mksoiltex( mksrf_fsoitex_mesh, file_mapunit_i=mksrf_fsoitex, file_lookup_i=mksrf_fsoitex_lookup, &
          mesh_o=mesh_model, pioid_o=pioid, rc=rc)
     flush(ndiag)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mksoiltex')
  end if

  ! -----------------------------------
  ! Make soil color classes [soicol] [fsoicol]
  ! -----------------------------------
  if (fsurdat /= ' ') then
     ! SOIL_COLOR and mxsoil_color is written out in the subroutine
     call mksoilcol( mksrf_fsoicol, mksrf_fsoicol_mesh, mesh_model, pioid, rc)
     flush(ndiag)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mksoilcol')
  end if

  ! -----------------------------------
  ! Make soil fmax [fmaxsoil]
  ! -----------------------------------
  if (fsurdat /= ' ') then
     ! FMAX is written out in the subroutine
     call mksoilfmax( mksrf_fmax_mesh, mksrf_fmax, mesh_model, pioid, rc)
     flush(ndiag)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mksoilfmax')
  end if

  ! -----------------------------------
  ! Make GDP data [gdp] from [gdp]
  ! -----------------------------------
  if (fsurdat /= ' ') then
     ! gdp is written out in the subroutine
     call mkgdp (mksrf_fgdp_mesh, mksrf_fgdp, mesh_model, pioid, rc=rc)
     flush(ndiag)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkdomain')
  end if

  ! -----------------------------------
  ! Make peat data [fpeat] from [peatf]
  ! -----------------------------------
  if (fsurdat /= ' ') then
     call mkpeat (mksrf_fpeat_mesh, mksrf_fpeat, mesh_model, pioid, rc=rc)
     flush(ndiag)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkpeat')
  end if

  ! -----------------------------------
  ! Make soil depth data [soildepth] from [soildepthf]
  ! -----------------------------------
  if (fsurdat /= ' ') then
     call mksoildepth( mksrf_fsoildepth_mesh, mksrf_fsoildepth, mesh_model, pioid, rc)
     flush(ndiag)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mksoildepth')
  end if

  ! -----------------------------------
  ! Make agricultural fire peak month data [abm] from [abm]
  ! -----------------------------------
  if (fsurdat /= ' ') then
     call mkagfirepkmon (mksrf_fabm_mesh, mksrf_fabm, mesh_model, pioid, rc=rc)
     flush(ndiag)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkagfirepkmon')
  end if

  ! -----------------------------------
  ! Make urban fraction [pcturb] from [furban] dataset and
  ! -----------------------------------
  allocate (pcturb(lsize_o))                 ; pcturb(:)            = spval
  allocate (pcturb_max(lsize_o, numurbl))    ; pcturb_max(:,:)      = spval
  allocate (urban_classes(lsize_o,numurbl))  ; urban_classes(:,:)   = spval
  allocate (urban_region(lsize_o))           ; urban_region(:)      = -999
  call mkurban(mksrf_furban_mesh, mksrf_furban, mesh_model, pcturb, &
               urban_classes, urban_region, rc=rc)
  flush(ndiag)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkurban')
  if (fsurdat /= ' ') then
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out urban region id"
     if (root_task)  flush(ndiag)
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
  flush(ndiag)
  if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkurban_topo')
  where (elev > elev_thresh)
     pcturb = 0._r8
  end where

  ! -----------------------------------
  ! Compute topography statistics [topo_stddev, slope] from [ftopostats]
  ! -----------------------------------
  if (fsurdat /= ' ') then
     call mktopostats ( mksrf_ftopostats_mesh, mksrf_ftopostats, mksrf_ftopostats_override, &
          mesh_model, pioid, rc=rc)
     flush(ndiag)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mktopostats')
  end if

  ! -----------------------------------
  ! Compute VIC parameters
  ! -----------------------------------
  if (fsurdat /= ' ') then
     if (outnc_vic) then
        call mkVICparams ( mksrf_fvic_mesh, mksrf_fvic, mesh_model, pioid, rc)
        flush(ndiag)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkVICparams')
     end if
  end if

  ! -----------------------------------
  ! Make VOC emission factors for isoprene [ef1_btr,ef1_fet,ef1_fdt,ef1_shr,ef1_grs,ef1_crp]
  ! -----------------------------------
  if (fsurdat /= ' ')  then
     call mkvocef ( mksrf_fvocef_mesh, mksrf_fvocef, mesh_model, pioid, lat, rc)
     flush(ndiag)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkvocef')
  end if

  ! -----------------------------------
  ! Save special land unit areas of surface dataset
  ! -----------------------------------

  pctwet_orig(:) = pctwet(:)
  pctgla_orig(:) = pctgla(:)

  ! -----------------------------------
  ! Perform other normalizations
  ! -----------------------------------

  ! Normalize land use and make sure things add up to 100% as well as
  ! checking that things are as they should be.
  allocate(landfrac_mksurfdata(lsize_o))
  call normalize_and_check_landuse(lsize_o)

  ! Write out sum of PFT's
  do k = natpft_lb,natpft_ub
     loc_suma = 0._r8
     do n = 1,lsize_o
        loc_suma = loc_suma + pctnatpft(n)%get_one_pct_p2g(k)
     enddo
     call mpi_reduce(loc_suma, glob_suma, 1, MPI_REAL8, MPI_SUM, 0, mpicom, ier)
     if (root_task) then
        write(ndiag,*) 'sum over domain of pft ',k,glob_suma
        flush(ndiag)
     end if
  enddo
  if (root_task) write(ndiag,*)
  do k = cft_lb,cft_ub
     loc_suma = 0._r8
     do n = 1,lsize_o
        loc_suma = loc_suma + pctcft(n)%get_one_pct_p2g(k)
     enddo
     call mpi_reduce(loc_suma, glob_suma, 1, MPI_REAL8, MPI_SUM, 0, mpicom, ier)
     if (root_task) then
        write(ndiag,*) 'sum over domain of cft ',k,glob_suma
     end if
  enddo
  if (root_task) write(ndiag,*)
  if (root_task) flush(ndiag)

  ! Make final values of percent urban by class and compute urban parameters
  ! This call needs to occur after all corrections are made to pcturb
  allocate (urban_classes_g(lsize_o,numurbl)); urban_classes_g(:,:) = spval
  call normalize_classes_by_gcell(urban_classes, pcturb, urban_classes_g)
  if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out percnt urban"
  if (root_task) flush(ndiag)

  ! Make Urban Parameters from raw input data and write to surface dataset
  ! Write to netcdf file is done inside mkurbanpar routine
  if (fsurdat /= ' ') then
     call mkurbanpar(mksrf_furban, pioid, mesh_model, urban_region, urban_classes_g, &
          urban_skip_abort_on_invalid_data_check)
     flush(ndiag)
  end if

  ! -----------------------------------
  ! Write out PCT_URBAN, PCT_GLACIER, PCT_LAKE and PCT_WETLAND and
  ! PCT_NATVEG, PCT_NAT_PFT, PCT_CROP and PCT_CFT
  ! -----------------------------------

  allocate(pctnatveg(lsize_o))
  allocate(pctcrop(lsize_o))
  allocate(pct_nat_pft(lsize_o, 0:num_natpft))
  if (num_cft > 0) then
     allocate(pct_cft(lsize_o, num_cft))
  end if

  if (fsurdat /= ' ') then
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out PCT_URBAN"
     if (root_task) flush(ndiag)
     call mkfile_output(pioid,  mesh_model,  'PCT_URBAN', urban_classes_g, lev1name='numurbl', rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_URBAN')

     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out PCT_GLACIER"
     if (root_task) flush(ndiag)
     call mkfile_output(pioid, mesh_model, 'PCT_GLACIER', pctgla, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in mkfile_output for pctgla')

     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out PCT_LAKE"
     if (root_task) flush(ndiag)
     call mkfile_output(pioid,  mesh_model,  'PCT_LAKE', pctlak, rc=rc)
     if (root_task) flush(ndiag)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in mkfile_output for pctlak')

     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out PCT_WETLAND"
     if (root_task) flush(ndiag)
     call mkfile_output(pioid, mesh_model,  'PCT_WETLAND', pctwet, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in in mkfile_output for pctwet')

     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out PCT_OCEAN"
     call mkfile_output(pioid, mesh_model,  'PCT_OCEAN', pctocn, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in in mkfile_output for pctocn')

     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing PCT_NATVEG"
     if (root_task) flush(ndiag)
     call get_pct_l2g_array(pctnatpft, pctnatveg)
     call mkfile_output(pioid, mesh_model, 'PCT_NATVEG', pctnatveg, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_NATVEG')

     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing PCT_CROP"
     if (root_task) flush(ndiag)
     call get_pct_l2g_array(pctcft, pctcrop)
     call mkfile_output(pioid, mesh_model, 'PCT_CROP', pctcrop, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_CROP')

     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing PCT_NAT_PFT"
     if (root_task) flush(ndiag)
     if (lsize_o /= 0) then
        call get_pct_p2l_array(pctnatpft, ndim1=lsize_o, ndim2=num_natpft+1, pct_p2l=pct_nat_pft)
     else
        pct_nat_pft(:,:) = 0.
     end if
     call mkfile_output(pioid, mesh_model, 'PCT_NAT_PFT', pct_nat_pft, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_NAT_PFT')

     if (num_cft > 0) then
        if (root_task)  write(ndiag, '(a)') trim(subname)//" writing PCT_CFT"
        if (root_task) flush(ndiag)
        if (lsize_o /= 0) then
           call get_pct_p2l_array(pctcft, ndim1=lsize_o, ndim2=num_cft, pct_p2l=pct_cft)
        else
           pct_cft(:,:) = 0.
        end if
        call mkfile_output(pioid, mesh_model, 'PCT_CFT', pct_cft, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_CFT')
     end if

     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing LANDFRAC_MKSURFDATA"
     if (root_task) flush(ndiag)
     call mkfile_output(pioid, mesh_model, 'LANDFRAC_MKSURFDATA', landfrac_mksurfdata, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for LANDFRAC_MKSURFDATA')
  end if

  ! ----------------------------------------------------------------------
  ! Make glacier multiple elevation classes [pctglcmec,topoglcmec] from [fglacier] dataset
  ! ----------------------------------------------------------------------
  ! This call needs to occur after pctgla has been adjusted for the final time
  if (fsurdat /= ' ') then
     call mkglcmecInit (pioid)
     call mkglcmec(mksrf_fglacier_mesh, mksrf_fglacier, mesh_model, pioid, rc=rc)
     flush(ndiag)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkglcmec')
  end if

  ! ----------------------------------------------------------------------
  ! Close surface dataset
  ! ----------------------------------------------------------------------
  if (fsurdat /= ' ') then
     call pio_closefile(pioid)
     if (root_task) then
        write(ndiag,*)
        write(ndiag,'(a)') 'Successfully created surface data output file = '//trim(fsurdat)
        write(ndiag,'(a)') '   This file contains the land model surface data'
        write(ndiag,*)
        flush(ndiag)
     end if
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
        flush(ndiag)
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

     ! End define mode
     rcode = pio_enddef(pioid)

     ! Write out model grid
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out LONGXY"
     if (root_task)  flush(ndiag)
     call mkfile_output(pioid, mesh_model, 'LONGXY', lon, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out LATIXY"
     if (root_task)  flush(ndiag)
     call mkfile_output(pioid, mesh_model, 'LATIXY', lat, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')

     ! Write out natpft
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out natpft"
     if (root_task)  flush(ndiag)
     rcode = pio_inq_varid(pioid, 'natpft', pio_varid)
     rcode = pio_put_var(pioid, pio_varid, (/(n,n=natpft_lb,natpft_ub)/))

     ! Write out cft
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out cft"
     if (root_task)  flush(ndiag)
     rcode = pio_inq_varid(pioid, 'cft', pio_varid)
     rcode = pio_put_var(pioid, pio_varid, (/(n,n=cft_lb,cft_ub)/))

     ! Write out LANDFRAC_PFT
     ! landfrac_pft was calculated ABOVE
     if (root_task)  write(ndiag, '(a)') trim(subname)//" writing land fraction calculated in fsurdata calc)"
     if (root_task)  flush(ndiag)
     call mkfile_output(pioid,  mesh_model, 'LANDFRAC_PFT', landfrac_pft, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for LANDFRAC_PFT')

     ! -----------------------------------------
     ! Read in each dynamic pft landuse dataset
     ! -----------------------------------------

     ! Open txt file
     if (root_task) then
        write(ndiag,'(a)')' Opening '//trim(mksrf_fdynuse)//' to read dynamic data forcing '
        flush(ndiag)
        open (newunit=nfdyn, file=trim(mksrf_fdynuse), form='formatted', iostat=ier)
        if (ier /= 0) then
           call shr_sys_abort(subname//" failed to open file "//trim(mksrf_fdynuse))
        end if
     end if

     pctnatpft_max = pctnatpft
     pctcft_max = pctcft
     pcturb_max = urban_classes_g
     pctlak_max = pctlak

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
           flush(ndiag)
        end if
        call mpi_bcast (fhrvname, len(fhrvname), MPI_CHARACTER, 0, mpicom, ier)
        call mpi_bcast (year2, 1, MPI_INTEGER, 0, mpicom, ier)
        if ( year2 /= year ) then
           if (root_task) then
              write(ndiag,*) subname, ' error: year for harvest not equal to year for PFT files'
              flush(ndiag)
           end if
           call shr_sys_abort()
        end if
        ! Read input urban data
        if (root_task) then
           read(nfdyn, '(A195,1x,I4)', iostat=ier) furbname, year2
           write(ndiag,*)'input urban dynamic dataset for year ', year2, ' is : ', trim(furbname)
           flush(ndiag)
        end if
        call mpi_bcast (furbname, len(furbname), MPI_CHARACTER, 0, mpicom, ier)
        call mpi_bcast (year2, 1, MPI_INTEGER, 0, mpicom, ier)
        if ( year2 /= year ) then
           if (root_task) then
              write(ndiag,*) subname, ' error: year for urban not equal to year for PFT files'
              flush(ndiag)
           end if
           call shr_sys_abort()
        end if
        ! Read input lake data
        if (root_task) then
           read(nfdyn, '(A195,1x,I4)', iostat=ier) flakname, year2
           write(ndiag,*)'input lake dynamic dataset for year ', year2, ' is : ', trim(flakname)
           flush(ndiag)
        end if
        call mpi_bcast (flakname, len(flakname), MPI_CHARACTER, 0, mpicom, ier)
        call mpi_bcast (year2, 1, MPI_INTEGER, 0, mpicom, ier)
        if ( year2 /= year ) then
           if (root_task) then
              write(ndiag,*) subname, ' error: year for lake not equal to year for PFT files'
              flush(ndiag)
           end if
           call shr_sys_abort()
        end if

        ntim = ntim + 1
        if (root_task) then
           write(ndiag,'(a,i8)')subname//' ntime = ',ntim
           flush(ndiag)
        end if

        rcode = pio_inq_varid(pioid, 'YEAR', pio_varid)
        rcode = pio_put_var(pioid, pio_varid, (/ntim/), year)
        rcode = pio_inq_varid(pioid, 'time', pio_varid)
        rcode = pio_put_var(pioid, pio_varid, (/ntim/), year)
       !rcode = pio_inq_varid(pioid, 'input_pftdata_filename', pio_varid)
       !rcode = pio_put_var(pioid, pio_varid, (/1,ntim/), (/len_trim(string),1/), trim(string))
        call pio_syncfile(pioid)

        ! Create pctpft data at model resolution from file fname
        ! Note that pctlnd_o below is different than the above call and returns pctlnd_pft_dyn
        call mkpft( mksrf_fvegtyp_mesh, fname, mesh_model, &
             pctlnd_o=pctlnd_pft_dyn, pctnatpft_o=pctnatpft, pctcft_o=pctcft, &
             rc=rc)
        flush(ndiag)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkpft')
        call pio_syncfile(pioid)

        ! Consistency check on input land fraction
        ! pctlnd_pft was calculated ABOVE
        do n = 1,lsize_o
           if (pctlnd_pft_dyn(n) /= pctlnd_pft(n)) then
              write(ndiag,*) subname,' error: pctlnd_pft for dynamics data = ',&
                   pctlnd_pft_dyn(n), ' not equal to pctlnd_pft for surface data = ',&
                   pctlnd_pft(n),' at n= ',n
              if ( trim(fname) == ' ' )then
                 write(ndiag,*) ' PFT string = ',trim(string)
              else
                 write(ndiag,*) ' PFT file = ', fname
              end if
              flush(ndiag)
              call shr_sys_abort()
           end if
        end do

        ! Create harvesting data at model resolution
        ! Output data is written in mkharvest
        call mkharvest( mksrf_fhrvtyp_mesh, fhrvname, mesh_model, pioid, &
                        ntime=ntim, rc=rc )
        flush(ndiag)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkharvest')
        call pio_syncfile(pioid)

        ! Create pctlak data at model resolution (use original mapping file from lake data)
        call mkpctlak(mksrf_fpctlak_mesh, flakname, mesh_model, pctlak, pioid, &
                      rc=rc)
        flush(ndiag)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkpctlak')
        call pio_syncfile(pioid)

        call mkurban(mksrf_furban_mesh, furbname, mesh_model, pcturb, &
                     urban_classes, urban_region, rc=rc)
        flush(ndiag)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkurban')
        call pio_syncfile(pioid)
        ! screen pcturb using elevation
        where (elev > elev_thresh)
           pcturb = 0._r8
        end where

        ! For landunits NOT read each year: reset to their pre-adjustment values
        ! in preparation for redoing landunit area normalization
        pctwet(:) = pctwet_orig(:)
        pctgla(:) = pctgla_orig(:)

        ! Normalize land use and make sure things add up to 100% as well as
        ! checking that things are as they should be.
        if (root_task) then
           write(ndiag,*)
           write(ndiag,'(1x,80a1)') ('=',k=1,80)
           write(ndiag,'(a)')' calling normalize_and_check_landuse'
           flush(ndiag)
        end if
        call normalize_and_check_landuse(lsize_o)
        call normalize_classes_by_gcell(urban_classes, pcturb, urban_classes_g)

        ! Given an array of pct_pft_type variables, update all the max_p2l variables.
        call update_max_array(pctnatpft_max, pctnatpft)
        call update_max_array(pctcft_max, pctcft)
        call update_max_array_urban(pcturb_max,urban_classes_g)
        call update_max_array_lake(pctlak_max,pctlak)

        if (root_task)  write(ndiag, '(a,i8)') trim(subname)//" writing PCT_NAT_PFT for year ",year
        if (root_task)  flush(ndiag)
        rcode = pio_inq_varid(pioid, 'PCT_NAT_PFT', pio_varid)
        call pio_setframe(pioid, pio_varid, int(ntim, kind=Pio_Offset_Kind))
        call get_pct_p2l_array(pctnatpft, ndim1=lsize_o, ndim2=num_natpft+1, pct_p2l=pct_nat_pft)
        call mkfile_output(pioid, mesh_model, 'PCT_NAT_PFT', pct_nat_pft, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_NAT_PFT')
        call pio_syncfile(pioid)

        if (root_task)  write(ndiag, '(a,i8)') trim(subname)//" writing PCT_CROP for year ",year
        if (root_task)  flush(ndiag)
        rcode = pio_inq_varid(pioid, 'PCT_CROP', pio_varid)
        call pio_setframe(pioid, pio_varid, int(ntim, kind=Pio_Offset_Kind))
        call get_pct_l2g_array(pctcft, pctcrop)
        call mkfile_output(pioid, mesh_model, 'PCT_CROP', pctcrop, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_CROP')
        call pio_syncfile(pioid)

        if (root_task)  write(ndiag, '(a,i8)') trim(subname)//" writing PCT_URBAN for year ",year
        if (root_task)  flush(ndiag)
        rcode = pio_inq_varid(pioid, 'PCT_URBAN', pio_varid)
        call pio_setframe(pioid, pio_varid, int(ntim, kind=Pio_Offset_Kind))
        call mkfile_output(pioid, mesh_model, 'PCT_URBAN', urban_classes_g, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_URBAN')
        call pio_syncfile(pioid)

        if (root_task)  write(ndiag, '(a,i8)') trim(subname)//" writing PCT_LAKE for year ",year
        if (root_task)  flush(ndiag)
        rcode = pio_inq_varid(pioid, 'PCT_LAKE', pio_varid)
        call pio_setframe(pioid, pio_varid, int(ntim, kind=Pio_Offset_Kind))
        call mkfile_output(pioid, mesh_model, 'PCT_LAKE', pctlak, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_LAKE')
        call pio_syncfile(pioid)

        if (num_cft > 0) then
           if (root_task)  write(ndiag, '(a,i8)') trim(subname)//" writing PCT_CFT for year ",year
           if (root_task)  flush(ndiag)
           rcode = pio_inq_varid(pioid, 'PCT_CFT', pio_varid)
           call pio_setframe(pioid, pio_varid, int(ntim, kind=Pio_Offset_Kind))
           call get_pct_p2l_array(pctcft, ndim1=lsize_o, ndim2=num_cft, pct_p2l=pct_cft)
           call mkfile_output(pioid, mesh_model, 'PCT_CFT', pct_cft, rc=rc)
           if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_CFT')
           call pio_syncfile(pioid)
        end if

        if (root_task)  write(ndiag, '(a,i8)') trim(subname)//" writing LANDFRAC_MKSURFDATA for year ",year
        if (root_task)  flush(ndiag)
        rcode = pio_inq_varid(pioid, 'LANDFRAC_MKSURFDATA', pio_varid)
        call pio_setframe(pioid, pio_varid, int(ntim, kind=Pio_Offset_Kind))
        call mkfile_output(pioid, mesh_model, 'LANDFRAC_MKSURFDATA', landfrac_mksurfdata, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for LANDFRAC_MKSURFDATA')
        call pio_syncfile(pioid)

        if (root_task) then
           write(ndiag,'(1x,80a1)') ('=',k=1,80)
           write(ndiag,*)
           flush(ndiag)
        end if

     end do   ! end of read loop

     if (root_task)  write(ndiag, '(a,i8)') trim(subname)//" writing PCT_NAT_PFT_MAX "
     if (root_task)  flush(ndiag)
     call get_pct_p2l_array(pctnatpft_max, ndim1=lsize_o, ndim2=num_natpft+1, pct_p2l=pct_nat_pft)
     call mkfile_output(pioid, mesh_model, 'PCT_NAT_PFT_MAX', pct_nat_pft, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_NAT_PFT')

     if (root_task)  write(ndiag, '(a,i8)') trim(subname)//" writing PCT_CROP_MAX"
     if (root_task)  flush(ndiag)
     call get_pct_l2g_array(pctcft_max, pctcrop)
     call mkfile_output(pioid, mesh_model, 'PCT_CROP_MAX', pctcrop, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_CROP_MAX')

     if (root_task)  write(ndiag, '(a,i8)') trim(subname)//" writing PCT_URBAN_MAX"
     if (root_task)  flush(ndiag)
     call mkfile_output(pioid, mesh_model, 'PCT_URBAN_MAX', pcturb_max, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_URBAN_MAX')

     if (root_task)  write(ndiag, '(a,i8)') trim(subname)//" writing PCT_LAKE_MAX"
     if (root_task)  flush(ndiag)
     call mkfile_output(pioid, mesh_model, 'PCT_LAKE_MAX', pctlak_max, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_LAKE_MAX')

     if (num_cft > 0) then
        if (root_task)  write(ndiag, '(a,i8)') trim(subname)//" writing PCT_CFT_MAX"
        if (root_task)  flush(ndiag)
        call get_pct_p2l_array(pctcft_max, ndim1=lsize_o, ndim2=num_cft, pct_p2l=pct_cft)
        call mkfile_output(pioid, mesh_model, 'PCT_CFT_MAX', pct_cft, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for PCT_CFT')
     end if

     ! Close the file
     call pio_closefile(pioid)
     if (root_task) then
        write(ndiag,*)
        write(ndiag,'(a)') 'Successfully created landuse timeseries data output file = '//trim(fdyndat)
        write(ndiag,'(a)') '   This file contains the land model time series data'
        write(ndiag,*)
        flush(ndiag)
     end if

  end if   ! end of if-create dynamic landust dataset

  ! -----------------------------------
  ! Wrap things up
  ! -----------------------------------
  if (root_task) then
     write(ndiag,'(a)') 'Successfully ran mksurfdata_esmf'
     close (ndiag)
  end if

  call ESMF_Finalize()

  !-----------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------

    subroutine normalize_and_check_landuse(ns_o)
      !
      ! Normalize land use and make sure things add up to 100% as well as
      ! checking that things are as they should be.
      !
      ! Coming into this subroutine, landunit areas are expressed as percent of the
      ! gridcell and are NOT normalized to sum to 100%. Coming out of this subroutine,
      ! landunit areas are expressed as percent of the land portion of the gridcell and
      ! ARE normalized to sum to 100%.
      !
      ! input/output variables:
      integer, intent(in) :: ns_o

      ! local variables:
      integer  :: k,n                         ! indices
      integer  :: nsmall                      ! number of small PFT values for a single check
      integer  :: nsmall_tot                  ! total number of small PFT values in all grid cells
      real(r8) :: suma                        ! sum for error check
      real(r8) :: pct_land                    ! area considered to be land (% of grid cell)
      real(r8) :: frac_land                   ! area considered to be land (fraction of grid cell)
      real(r8) :: new_total_natveg_pct        ! new % veg (% of grid cell, natural veg)
      real(r8), parameter :: tol_loose = 1.e-4_r8               ! tolerance for some 'loose' error checks
      real(r8), parameter :: toosmallPFT = 1.e-10_r8            ! tolerance for PFT's to ignore
      character(len=32) :: subname = 'normalize_and_check_landuse'  ! subroutine name
      !-----------------------------------------------------------------------

      do n = 1,ns_o

         ! Truncate all percentage fields on output grid. This is needed to
         ! ensure that wt is zero (not a very small number such as
         ! 1e-16) where it really should be zero
         pctlak(n) = float(nint(pctlak(n)))
         pctwet(n) = float(nint(pctwet(n)))
         pctgla(n) = float(nint(pctgla(n)))

         ! Check preconditions
         if ( pctlak(n) < 0.0_r8 )then
            write(ndiag,*) subname, ' ERROR: pctlak is negative!'
            write(ndiag,*) 'n, pctlak = ', n, pctlak(n)
            flush(ndiag)
            call shr_sys_abort()
         end if
         if ( pctwet(n) < 0.0_r8 )then
            write(ndiag,*) subname, ' ERROR: pctwet is negative!'
            write(ndiag,*) 'n, pctwet = ', n, pctwet(n)
            flush(ndiag)
            call shr_sys_abort()
         end if
         if ( pcturb(n) < 0.0_r8 )then
            write(ndiag,*) subname, ' ERROR: pcturb is negative!'
            write(ndiag,*) 'n, pcturb = ', n, pcturb(n)
            flush(ndiag)
            call shr_sys_abort()
         end if
         if ( pctgla(n) < 0.0_r8 )then
            write(ndiag,*) subname, ' ERROR: pctgla is negative!'
            write(ndiag,*) 'n, pctgla = ', n, pctgla(n)
            flush(ndiag)
            call shr_sys_abort()
         end if
         if ( pctcft(n)%get_pct_l2g() < 0.0_r8 )then
            write(ndiag,*) subname, ' ERROR: pctcrop is negative!'
            write(ndiag,*) 'n, pctcrop = ', n, pctcft(n)%get_pct_l2g()
            flush(ndiag)
            call shr_sys_abort()
         end if

         ! Make sure sum of all land cover types except natural vegetation does
         ! not exceed 100. If it does, subtract excess from these land cover
         ! types proportionally.
         suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n) + pctcft(n)%get_pct_l2g()
         if (suma > 100._r4) then
            pctlak(n) = pctlak(n) * 100._r8/suma
            pctwet(n) = pctwet(n) * 100._r8/suma
            pcturb(n) = pcturb(n) * 100._r8/suma
            pctgla(n) = pctgla(n) * 100._r8/suma
            call pctcft(n)%set_pct_l2g(pctcft(n)%get_pct_l2g() * 100._r8/suma)
         end if

         suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n) + pctcft(n)%get_pct_l2g()
         if (suma > (100._r8 + tol_loose)) then
            write(ndiag,*) subname, ' ERROR: pctlak + pctwet + pcturb + pctgla + pctcrop must be'
            write(ndiag,*) '<= 100% before normalizing natural vegetation area'
            write(ndiag,*) 'n, pctlak, pctwet, pcturb, pctgla, pctcrop = ', &
                 n, pctlak(n), pctwet(n), pcturb(n), pctgla(n), pctcft(n)%get_pct_l2g()
            flush(ndiag)
            call shr_sys_abort()
         end if

         ! Determine the percent of each grid cell considered to be land. (See comments in
         ! https://github.com/ESCOMP/CTSM/issues/1716 for details.)
         !
         ! Start by using the land fraction field from the PFT raw data set:
         pct_land = pctlnd_pft(n)
         !
         ! Brief summary of the following: But we don't want to overwrite special
         ! landunits or crop with ocean where these special landunits extend beyond the
         ! PFT data's land fraction.
         !
         ! More details:
         !
         ! In essence, this is saying that we'll let special landunit area grow into the
         ! natveg area before growing into ocean, but we'll have special landunit area
         ! grow into ocean before growing into crop or any other special landunit area.
         ! (This check of special landunit area is particularly important for glaciers,
         ! where we can have floating ice shelves, so we can have a scenario where
         ! pctlnd_pft is 0 but we have non-zero glacier cover and we want the final grid
         ! cell to be glacier-covered.) (We could possibly do better by considering the
         ! land mask from each special landunit raw dataset, and even better by mapping
         ! these various notions of land mask onto each other, but that starts to get
         ! messy, and relies on the trustworthiness of each raw dataset's land mask...
         ! this formulation seems reasonable enough.)
         !
         ! Note that we include pct_crop in the following, but NOT pct_natveg. The
         ! assumption behind that is that pct_crop is more reliable and/or more important,
         ! and should not be overwritten, whereas pct_natveg can be overwritten with a
         ! special landunit. For example, consider a case where the PFT dataset specifies
         ! 40% land, with 20% crop and 10% natveg (so, implicitly, 10% special landunits).
         ! If the only special landunit is glacier and it has 15% cover, then we want to
         ! end up with the glacier overwriting natural veg, so we end up with 20% crop, 5%
         ! natveg, 15% glacier and 60% ocean. However, if glacier has 30% cover, then we
         ! will assume that some of that glacier extends over the ocean rather than
         ! overwriting crop, so we end up with 20% crop, 30% glacier, 50% ocean and 0%
         ! natveg.
         !
         ! Another reason for excluding pct_natveg from the following is more pragmatic:
         ! Typically, we expect the initial sum of all landunit areas to approximately
         ! equal pctlnd_pft. This means that, in a coastal grid cell, if we included all
         ! landunit areas in the following sum, then in many cases we would take the
         ! final pct_land from the landunit sum rather than from pctlnd_pft. But in this
         ! scenario where we take pct_land from the landunit sum, it is likely that
         ! pct_land will vary from year to year on the landuse timeseries file. This
         ! variation from year to year will cause a renormalization of all landunits,
         ! leading to changes in the areas of landunits that should really stay fixed
         ! from one year to the next. By excluding pct_natveg we give more wiggle room:
         ! it will usually be the case that we take the final pct_land from pctlnd_pft,
         ! which stays fixed in time, so that the renormalization stays the same from one
         ! year to the next. The possible downside is that pct_natveg may end up varying
         ! more than is ideal, but this seems better than letting all of the other
         ! landunits vary when they should stay fixed.
         !
         ! Also, note that this landunit sum agrees with the suma sums elsewhere; this
         ! agreement may be important in some cases (so that if we changed the set of
         ! landunits included in the sum here, some other changes may be needed below.)
         pct_land = max(pct_land, suma)
         !
         ! Make sure that we're not ending up with > 100% land area. This can arise if the
         ! sum of special landunits exceeds 100% by roundoff; also, due to rounding
         ! errors, pctlnd_pft may slightly differ from 100% when it should really be
         ! exactly 100%; we want pct_land to be exactly 100% in this case:
         if (pct_land > (100._r8 - 1.e-4_r8)) then
            pct_land = 100._r8
         end if

         if (pct_land < 1.e-6_r8) then
            ! If we have essentially 0 land area, set land area to exactly 0 and put all
            ! area in pctocn. Note that, based on the formulation
            ! for pct_land above, this should only arise if the non-natveg landunits
            ! already have near-zero area (and the natveg landunit should also have
            ! near-zero area in this case, because its area should be no greater than the
            ! land fraction from the PFT raw dataset), so the zeroing of these other
            ! landunits should only result in changes near the roundoff level.
            pct_land = 0._r8
            frac_land = 0._r8
            call pctnatpft(n)%set_pct_l2g(0._r8)
            call pctcft(n)%set_pct_l2g(0._r8)
            pctlak(n) = 0._r8
            pcturb(n) = 0._r8
            pctgla(n) = 0._r8
            pctwet(n) = 0._r8
            pctocn(n) = 100._r8  ! the only asignment of non-zero ocean
         else
            ! Fill the rest of the land area with natveg, then renormalize landunits so
            ! that they are expressed as percent of the land area rather than percent of
            ! the full gridcell area.

            ! First fill the rest of the land area with natveg:
            new_total_natveg_pct = pct_land - suma
            ! Based on the formulation of pct_land above, pct_land is guaranteed to be at
            ! least as large as suma. However, correct it just in case there is rounding
            ! error:
            new_total_natveg_pct = max(new_total_natveg_pct, 0._r8)

            ! Now renormalize landunit areas so they are expressed as percent of the land
            ! area rather than percent of the total gridcell area. Note that we have
            ! already corrected pct_land above so that if it was slightly less than 100
            ! it was set to exactly 100, so we can check for any values less than 100
            ! here (rather than having some tolerance):
            frac_land = pct_land / 100._r8
            if (pct_land < 100._r8) then
               new_total_natveg_pct = new_total_natveg_pct / frac_land
               call pctcft(n)%set_pct_l2g(pctcft(n)%get_pct_l2g() / frac_land)
               pctlak(n) = pctlak(n) / frac_land
               pcturb(n) = pcturb(n) / frac_land
               pctgla(n) = pctgla(n) / frac_land
               pctwet(n) = pctwet(n) / frac_land
            end if

            ! Finally, set the actual pct_natveg:
            call pctnatpft(n)%set_pct_l2g(new_total_natveg_pct)
         end if

         ! Save landfrac for output to file
         landfrac_mksurfdata(n) = frac_land

         ! Confirm that we have done the rescaling correctly: now the sum of all landunits
         ! should be 100% within tol_loose
         suma = pctlak(n) + pctwet(n) + pctgla(n) + pcturb(n) + pctocn(n) +  &
            pctcft(n)%get_pct_l2g() + pctnatpft(n)%get_pct_l2g()
         if (abs(suma - 100._r8) > tol_loose) then
            write(ndiag,*) subname, ' ERROR: landunits do not sum to 100%'
            write(ndiag,*) 'n, suma, pctlak, pctwet, pctgla, pcturb, pctnatveg, pctcrop, pctocn = '
            write(ndiag,*) n, suma, pctlak(n), pctwet(n), pctgla(n), pcturb(n), &
                 pctnatpft(n)%get_pct_l2g(), pctcft(n)%get_pct_l2g(), pctocn(n)
            flush(ndiag)
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

         ! Include pctocn in suma but do not include in the
         ! renormalization. When pctocn /= 0, it is 100, and
         ! all other terms are 0.
         suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n) + pctocn(n) +  &
            pctnatpft(n)%get_pct_l2g() + pctcft(n)%get_pct_l2g()
         if ( abs(suma - 100.0_r8) > 2.0*epsilon(suma) )then
            pctlak(n)    = pctlak(n)    * 100._r8/suma
            pctwet(n)    = pctwet(n)    * 100._r8/suma
            pcturb(n)    = pcturb(n)    * 100._r8/suma
            pctgla(n)    = pctgla(n)    * 100._r8/suma
            call pctnatpft(n)%set_pct_l2g(pctnatpft(n)%get_pct_l2g() * 100._r8/suma)
            call pctcft(n)%set_pct_l2g(pctcft(n)%get_pct_l2g() * 100._r8/suma)
         end if

         ! This roundoff error fix is needed to handle the situation where new_total_natveg_pct
         ! ends up near 0 but not exactly 0 due to roundoff issues. In this situation, we set the
         ! natveg landunit area to exactly 0 and put the remainder into some other landunit. Since
         ! the remainder is very small, it doesn't really matter which other landunit we add it to,
         ! so we just pick some landunit that already has at least 1% cover.
         suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n) + pctcft(n)%get_pct_l2g()
         if ( (suma < 100._r8 .and. suma > (100._r8 - 1.e-6_r8)) .or. &
              (pctnatpft(n)%get_pct_l2g() > 0.0_r8 .and. pctnatpft(n)%get_pct_l2g() <  1.e-6_r8) ) then
            if ( root_task ) then
               write (ndiag,*) 'Special plus crop land units near 100%, but not quite for n,suma =',n,suma
               write (ndiag,*) 'Adjusting special plus crop land units to 100%'
               flush(ndiag)
            end if
            if (pctlak(n) >= 1.0_r8) then
               pctlak(n) = 100._r8 - (pctwet(n) + pcturb(n) + pctgla(n) + pctcft(n)%get_pct_l2g())
            else if (pctwet(n) >= 1.0_r8) then
               pctwet(n) = 100._r8 - (pctlak(n) + pcturb(n) + pctgla(n) + pctcft(n)%get_pct_l2g())
            else if (pcturb(n) >= 1.0_r8) then
               pcturb(n) = 100._r8 - (pctlak(n) + pctwet(n) + pctgla(n) + pctcft(n)%get_pct_l2g())
            else if (pctgla(n) >= 1.0_r8) then
               pctgla(n) = 100._r8 - (pctlak(n) + pctwet(n) + pcturb(n) + pctcft(n)%get_pct_l2g())
            else if (pctcft(n)%get_pct_l2g() >= 1.0_r8) then
               call pctcft(n)%set_pct_l2g(100._r8 - (pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n)))
            else
               write (ndiag,*) subname, 'Error: sum of special plus crop land units nearly 100% but none is >= 1% at ', &
                    'n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),pctnatveg(n),pctcrop(n),suma = ', &
                    n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),&
                    pctnatpft(n)%get_pct_l2g(),pctcft(n)%get_pct_l2g(),suma
               flush(ndiag)
               call shr_sys_abort()
            end if
            call pctnatpft(n)%set_pct_l2g(0._r8)
         end if
         if ( any(pctnatpft(n)%get_pct_p2g() > 0.0_r8 .and. pctnatpft(n)%get_pct_p2g() < toosmallPFT ) .or. &
              any(pctcft(n)%get_pct_p2g()    > 0.0_r8 .and. pctcft(n)%get_pct_p2g()    < toosmallPFT )) then
            write (6,*) 'pctnatpft or pctcft is small at n=', n
            write (6,*) 'pctnatpft%pct_p2l = ', pctnatpft(n)%get_pct_p2l()
            write (6,*) 'pctcft%pct_p2l = ', pctcft(n)%get_pct_p2l()
            write (6,*) 'pctnatpft%pct_l2g = ', pctnatpft(n)%get_pct_l2g()
            write (6,*) 'pctcft%pct_l2g = ', pctcft(n)%get_pct_l2g()
            flush(6)
            call shr_sys_abort()
         end if

         suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n) + pctocn(n) +  &
            pctcft(n)%get_pct_l2g() + pctnatpft(n)%get_pct_l2g()
         if ( abs(suma-100._r8) > 1.e-10_r8) then
            write (6,*) subname, ' error: sum of pctocn, pctlak, pctwet,', &
                 'pcturb, pctgla, pctnatveg and pctcrop is NOT equal to 100'
            write (6,*)'n,pctcon,pctlak,pctwet,pcturb,pctgla,pctnatveg,pctcrop,sum= ', &
                 n,pctocn(n),pctlak(n),pctwet(n),pcturb(n),pctgla(n),&
                 pctnatpft(n)%get_pct_l2g(),pctcft(n)%get_pct_l2g(), suma
            flush(6)
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
            write(ndiag,*) 'sum(urban_classes(n,:)) != 100: ', n, sum(urban_classes(n,:))
            flush(ndiag)
            call shr_sys_abort()
         end if
      end do

      if (root_task) then
         if ( nsmall_tot > 0 )then
            write(ndiag,*)'number of small pft = ', nsmall_tot
            flush(ndiag)
         end if
      end if

    end subroutine normalize_and_check_landuse

end program mksurfdata
