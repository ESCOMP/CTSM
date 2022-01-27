module mkvarctl

  !-----------------------------------------------------------------------
  ! Module containing control variables
  !-----------------------------------------------------------------------

  use ESMF
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_abort

  implicit none
  private

#include <mpif.h>

  ! ---------------------
  ! routines
  ! ---------------------

  public :: read_namelist_input
  public :: check_namelist_input
  public :: write_namelist_input

  ! ---------------------
  ! variables
  ! ---------------------

  integer, public :: ndiag           ! output log unit
  logical, public :: root_task       ! proc 0 logical for printing msgs
  integer, public :: iam             ! processor number
  integer, public :: npes            ! number of processors
  integer, public :: mpicom          ! communicator group

  ! Values from mpif.h that will be used
  public :: MPI_INTEGER
  public :: MPI_REAL8
  public :: MPI_LOGICAL
  public :: MPI_CHARACTER
  public :: MPI_COMM_WORLD
  public :: MPI_MAX
  public :: MPI_SUCCESS

  real(r8), public, parameter :: spval = 1.e36    ! special value
  integer,  public, parameter :: ispval = -9999   ! special value

  character(len=256) , public :: fgrddat  ! grid data file
  character(len=256) , public :: fsurdat  ! output surface data file name (if blank, do not output a surface dataset)
  character(len=256) , public :: fsurlog  ! output surface log file name
  character(len=256) , public :: fdyndat  ! dynamic landuse data file name
  character(len=256) , public :: fhrvname ! generic harvest filename
  logical            , public :: all_veg  ! if gridcell will be 100% vegetated land-cover
  real(r8)           , public :: std_elev = -999.99_r8  ! Standard deviation of elevation (m) to use for ent

  ! NOTE(bja, 2015-01) added to work around a ?bug? causing 1x1_urbanc_alpha to abort. See
  !/glade/p/cesm/cseg/inputdata/lnd/clm2/surfdata_map/README_c141219
  logical :: urban_skip_abort_on_invalid_data_check

  logical, public    :: outnc_large_files     ! output files in 64-bit format for large files
  logical, public    :: outnc_double          ! output ALL data in files as 64-bit
  integer, public    :: outnc_dims = 2        ! only applicable to lat/lon grids
  logical, public    :: outnc_1d              ! true => output file is 1d
  logical, public    :: outnc_vic             ! true => output VIC fields
  logical, public    :: outnc_3dglc           ! true => output 3D glacier fields

  character(len= 32), public :: mksrf_gridnm                 = ' ' ! name of grid to use on output file
  character(len=512), public :: mksrf_gridtype               = ' ' ! land gridtype, global or reg
  character(len=512), public :: mksrf_fgrid_mesh             = ' ' ! land grid file name to use
  integer           , public :: mksrf_fgrid_mesh_nx
  integer           , public :: mksrf_fgrid_mesh_ny

  character(len=512), public :: mksrf_fvegtyp                = ' ' ! vegetation data file name
  character(len=512), public :: mksrf_fhrvtyp                = ' ' ! harvest data file name
  character(len=512), public :: mksrf_fsoitex                = ' ' ! soil texture data file name
  character(len=512), public :: mksrf_forganic               = ' ' ! organic matter data file name
  character(len=512), public :: mksrf_fsoicol                = ' ' ! soil color data file name
  character(len=512), public :: mksrf_fabm                   = ' ' ! ag fire peak month and
  character(len=512), public :: mksrf_fpeat                  = ' ' ! peatlands and
  character(len=512), public :: mksrf_fsoildepth             = ' ' ! soil depth file name
  character(len=512), public :: mksrf_fgdp                   = ' ' ! gdp data file names
  character(len=512), public :: mksrf_flakwat                = ' ' ! inland lake data file name
  character(len=512), public :: mksrf_fwetlnd                = ' ' ! inland wetlands data file name
  character(len=512), public :: mksrf_furban                 = ' ' ! urban data file name
  character(len=512), public :: mksrf_fglacier               = ' ' ! glacier data file name
  character(len=512), public :: mksrf_fglacierregion         = ' ' ! glacier region data file name
  character(len=512), public :: mksrf_furbtopo               = ' ' ! urban topography data file name
  character(len=512), public :: mksrf_fmax                   = ' ' ! fmax data file name
  character(len=512), public :: mksrf_flai                   = ' ' ! lai data filename
  character(len=512), public :: mksrf_fdynuse                = ' ' ! ascii file containing names of dynamic land use files
  character(len=512), public :: mksrf_fvocef                 = ' ' ! VOC Emission Factor data file name
  character(len=512), public :: mksrf_ftopostats             = ' ' ! topography statistics data file name
  character(len=512), public :: mksrf_fvic                   = ' ' ! VIC parameters data file name
  character(len=512), public :: mksrf_irrig                  = ' ' ! TODO: should this namelist be here?

  character(len=512), public :: mksrf_fvegtyp_mesh           = ' ' ! vegetation mesh file name
  character(len=512), public :: mksrf_fhrvtyp_mesh           = ' ' ! harvest mesh file name
  character(len=512), public :: mksrf_fsoitex_mesh           = ' ' ! soil texture mesh file name
  character(len=512), public :: mksrf_forganic_mesh          = ' ' ! organic matter mesh file name
  character(len=512), public :: mksrf_fsoicol_mesh           = ' ' ! soil color mesh file name
  character(len=512), public :: mksrf_fabm_mesh              = ' ' ! ag fire peak month and
  character(len=512), public :: mksrf_fpeat_mesh             = ' ' ! peatlands and
  character(len=512), public :: mksrf_fsoildepth_mesh        = ' ' ! soil depth file name
  character(len=512), public :: mksrf_fgdp_mesh              = ' ' ! gdp mesh file names
  character(len=512), public :: mksrf_flakwat_mesh           = ' ' ! inland lake mesh file name
  character(len=512), public :: mksrf_fwetlnd_mesh           = ' ' ! inland wetlands mesh file name
  character(len=512), public :: mksrf_furban_mesh            = ' ' ! urban mesh file name
  character(len=512), public :: mksrf_fglacier_mesh          = ' ' ! glacier mesh file name
  character(len=512), public :: mksrf_fglacierregion_mesh    = ' ' ! glacier region mesh file name
  character(len=512), public :: mksrf_furbtopo_mesh          = ' ' ! urban topography mesh file name
  character(len=512), public :: mksrf_fmax_mesh              = ' ' ! fmax mesh file name
  character(len=512), public :: mksrf_flai_mesh              = ' ' ! lai mesh filename
  character(len=512), public :: mksrf_fhrv_mesh              = ' ' ! harvest mesh filename
  character(len=512), public :: mksrf_fdynuse_mesh_mesh_mesh = ' ' ! ascii file containing names of dynamic land use files
  character(len=512), public :: mksrf_fvocef_mesh            = ' ' ! VOC Emission Factor mesh file name
  character(len=512), public :: mksrf_ftopostats_mesh        = ' ' ! topography statistics mesh file name
  character(len=512), public :: mksrf_fvic_mesh              = ' ' ! VIC parameters mesh file name
  character(len=512), public :: mksrf_irrig_mesh             = ' ' ! TODO: should this namelist be here?

  integer           , public :: numpft  = 16   ! number of plant types
  !
  ! Variables to override data read in with
  ! (all_urban is mostly for single-point mode, but could be used for sensitivity studies)
  logical,  public :: all_urban              ! output ALL data as 100% covered in urban
  logical,  public :: no_inlandwet           ! set wetland to 0% over land; wetland will only be used for ocean points

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine read_namelist_input()

    ! local variables
    integer :: ier
    integer :: fileunit
    logical :: lexist
    ! ------------------------------------------------------------

    namelist /mksurfdata_input/      &
         mksrf_fvegtyp,             &
         mksrf_fhrvtyp,             &
         mksrf_fsoitex,             &
         mksrf_forganic,            &
         mksrf_fsoicol,             &
         mksrf_fvocef,              &
         mksrf_flakwat,             &
         mksrf_fwetlnd,             &
         mksrf_fglacier,            &
         mksrf_fglacierregion,      &
         mksrf_furbtopo,            &
         mksrf_fmax,                &
         mksrf_furban,              &
         mksrf_flai,                &
         mksrf_fdynuse,             &
         mksrf_fgdp,                &
         mksrf_fpeat,               &
         mksrf_fsoildepth,          &
         mksrf_fabm,                &
         mksrf_ftopostats,          &
         mksrf_fvic,                &
         mksrf_flai_mesh,           &
         mksrf_fhrv_mesh,           &
         mksrf_forganic_mesh,       &
         mksrf_fsoicol_mesh,        &
         mksrf_fsoitex_mesh,        &
         mksrf_fmax_mesh,           &
         mksrf_flakwat_mesh,        &
         mksrf_fwetlnd_mesh,        &
         mksrf_fvocef_mesh,         &
         mksrf_furban_mesh,         &
         mksrf_furbtopo_mesh,       &
         mksrf_fglacier_mesh,       &
         mksrf_fglacierregion_mesh, &
         mksrf_fgdp_mesh,           &
         mksrf_fpeat_mesh,          &
         mksrf_fsoildepth_mesh,     &
         mksrf_fabm_mesh,           &
         mksrf_ftopostats_mesh,     &
         mksrf_fvic_mesh,           &
         mksrf_fvegtyp_mesh,        &
         mksrf_fgrid_mesh,          &
         mksrf_fgrid_mesh_nx,       &
         mksrf_fgrid_mesh_ny,       &
         numpft,                    &
         all_veg,                   &
         all_urban,                 &
         no_inlandwet,              &
#ifdef TODO
         nglcec,                    &
         pft_idx,                   &
         pft_frc,                   &
         gitdescribe,               &
#endif
         outnc_large_files,         &
         outnc_double,              &
         outnc_dims,                &
         outnc_vic,                 &
         outnc_3dglc,               &
         fsurdat,                   &
         fdyndat,                   &
         fsurlog,                   &
         std_elev,                  &
         urban_skip_abort_on_invalid_data_check

    ! Set default namelist values - make these the defaults in gen_mksurfdata_namelist.py
    mksrf_gridtype    = 'global'
    outnc_large_files = .false.
    outnc_double      = .true.
    outnc_vic         = .false.
    outnc_3dglc       = .false.
    all_urban         = .false.
    all_veg           = .false.
    no_inlandwet      = .true.
    urban_skip_abort_on_invalid_data_check = .false.   ! default value for bug work around

    if (root_task) then
       write(6,*) 'Attempting to initialize control settings .....'
    end if

    if (root_task) then
       inquire (file='mksurfdata_in', exist=lexist)
       if (.not. lexist) then
          call shr_sys_abort('mksurfdata_in does not exist')
       end if
       open(newunit=fileunit, status="old", file="mksurfdata_in")
       read(fileunit, nml=mksurfdata_input, iostat=ier)
       if (ier > 0) then
          call shr_sys_abort('error reading in mksurfdata_input namelist from mksurfdata_in')
       end if
       close(fileunit)
    end if

    call mpi_bcast (mksrf_fgrid_mesh, len(mksrf_fgrid_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fvegtyp, len(mksrf_fvegtyp), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fhrvtyp, len(mksrf_fhrvtyp), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fsoitex, len(mksrf_fsoitex), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_forganic, len(mksrf_forganic), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fsoicol, len(mksrf_fsoicol), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fabm, len(mksrf_fabm), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fpeat, len(mksrf_fpeat), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fsoildepth, len(mksrf_fsoildepth), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fgdp, len(mksrf_fgdp), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_flakwat, len(mksrf_flakwat), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fwetlnd, len(mksrf_fwetlnd), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_furban, len(mksrf_furban), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fglacier, len(mksrf_fglacier), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fglacierregion, len(mksrf_fglacierregion), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_furbtopo, len(mksrf_furbtopo), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fmax, len(mksrf_fmax), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_flai, len(mksrf_flai), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fdynuse, len(mksrf_fdynuse), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fvocef, len(mksrf_fvocef), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_ftopostats, len(mksrf_ftopostats), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fvic, len(mksrf_fvic), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fvegtyp_mesh, len(mksrf_fvegtyp_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fhrvtyp_mesh, len(mksrf_fhrvtyp_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fsoitex_mesh, len(mksrf_fsoitex_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_forganic_mesh, len(mksrf_forganic_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fsoicol_mesh, len(mksrf_fsoicol_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fabm_mesh, len(mksrf_fabm_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fpeat_mesh, len(mksrf_fpeat_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fsoildepth_mesh, len(mksrf_fsoildepth_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fgdp_mesh, len(mksrf_fgdp_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_flakwat_mesh, len(mksrf_flakwat_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fwetlnd_mesh, len(mksrf_fwetlnd_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_furban_mesh, len(mksrf_furban_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fglacier_mesh, len(mksrf_fglacier_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fglacierregion_mesh, len(mksrf_fglacierregion_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_furbtopo_mesh, len(mksrf_furbtopo_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fmax_mesh, len(mksrf_fmax_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_flai_mesh, len(mksrf_flai_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fdynuse_mesh_mesh_mesh, len(mksrf_fdynuse_mesh_mesh_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fvocef, len(mksrf_fvocef), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_ftopostats, len(mksrf_ftopostats), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fvic, len(mksrf_fvic), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsurdat, len(fsurdat), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_fgrid_mesh_nx, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fgrid_mesh_ny, 1, MPI_INTEGER, 0, mpicom, ier)

    call mpi_bcast (outnc_dims, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (outnc_large_files, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (outnc_double, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (outnc_1d, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (outnc_vic, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (outnc_3dglc, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (all_urban, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (all_veg, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (no_inlandwet, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (urban_skip_abort_on_invalid_data_check, 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (numpft, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (std_elev, 1, MPI_REAL, 0, mpicom, ier)

  end subroutine read_namelist_input

  !===============================================================
  subroutine check_namelist_input()

    ! error check on namelist input
    if (mksrf_fgrid_mesh /= ' ')then
       fgrddat = mksrf_fgrid_mesh
    else
       call shr_sys_abort(" must specify mksrf_fgrid_mesh")
    endif

    if (trim(mksrf_gridtype) /= 'global' .and. trim(mksrf_gridtype) /= 'regional') then
       call shr_sys_abort(" mksrf_gridtype "//trim(mksrf_gridtype)//" is not supported")
    endif

#ifdef TODO
    if (nglcec <= 0) then
       call shr_sys_abort('nglcec must be at least 1')
    end if
#endif

    if (mksrf_fgrid_mesh_ny == 1) then
       outnc_1d = .true.
       outnc_dims = 1
    else
       outnc_1d = .false.
       outnc_dims = 2
    end if

  end subroutine check_namelist_input

  !===============================================================
  subroutine write_namelist_input()

    ! Note - need to call this after ndiag has been set

    if (root_task) then
       write(ndiag,'(a)')'Input rawdata files and corresponding meshes'
       write(ndiag,'(a)')' PFTs from:                  '//trim(mksrf_fvegtyp)
       write(ndiag,'(a)')' mesh for pft                '//trim(mksrf_fvegtyp_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' inland lake from:           '//trim(mksrf_flakwat)
       write(ndiag,'(a)')' mesh for lake water         '//trim(mksrf_flakwat_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' inland wetland from:        '//trim(mksrf_fwetlnd)
       write(ndiag,'(a)')' mesh for wetland            '//trim(mksrf_fwetlnd_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' soil texture from:          '//trim(mksrf_fsoitex)
       write(ndiag,'(a)')' mesh for soil texture       '//trim(mksrf_fsoitex_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' soil organic from:          '//trim(mksrf_forganic)
       write(ndiag,'(a)')' mesh for soil organic       '//trim(mksrf_forganic_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' soil color from:            '//trim(mksrf_fsoicol)
       write(ndiag,'(a)')' mesh for soil color         '//trim(mksrf_fsoicol_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' fmax from:                  '//trim(mksrf_fmax)
       write(ndiag,'(a)')' mesh for fmax               '//trim(mksrf_fmax_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' soil depth from:            '//trim(mksrf_fsoildepth)
       write(ndiag,'(a)')' mesh for soil depth         '//trim(mksrf_fsoildepth_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' VOC emission factors from:  '//trim(mksrf_fvocef)
       write(ndiag,'(a)')' mesh for VOC pct emis       '//trim(mksrf_fvocef_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' gdp from:                   '//trim(mksrf_fgdp)
       write(ndiag,'(a)')' mesh for gdp                '//trim(mksrf_fgdp_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' peat from:                  '//trim(mksrf_fpeat)
       write(ndiag,'(a)')' mesh for peatlands          '//trim(mksrf_fpeat_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' harvest from:               '//trim(mksrf_fhrvtyp)
       write(ndiag,'(a)')' mesh for harvest            '//trim(mksrf_fhrv_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' topography statistics from: '//trim(mksrf_ftopostats)
       write(ndiag,'(a)')' mesh for topography stats   '//trim(mksrf_ftopostats_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' glaciers from:              '//trim(mksrf_fglacier)
#ifdef TODO
       write(ndiag,'(a)')'            with:            '// nglcec, ' glacier elevation classes'
#endif
       write(ndiag,'(a)')' mesh for glacier            '//trim(mksrf_fglacier_mesh)
       write(ndiag,'(a)')' glacier region ID from:     '//trim(mksrf_fglacierregion)
       write(ndiag,'(a)')' mesh for glacier region     '//trim(mksrf_fglacierregion_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' urban from:                 '//trim(mksrf_furban)
       write(ndiag,'(a)')' mesh for urban              '//trim(mksrf_furban_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' urban topography from:      '//trim(mksrf_furbtopo)
       write(ndiag,'(a)')' mesh for urban topography   '//trim(mksrf_furbtopo_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' mesh for lai/sai            '//trim(mksrf_flai_mesh)
       write(ndiag,'(a)')' mesh for ag fire pk month   '//trim(mksrf_fabm_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' abm from:                   '//trim(mksrf_fabm)
       write(ndiag,*)
       write(ndiag,'(a)')' VIC parameters from:        '//trim(mksrf_fvic)
       write(ndiag,'(a)')' mesh for VIC parameters     '//trim(mksrf_fvic_mesh)
       write(ndiag,*)
       if (mksrf_fdynuse /= ' ') then
          write(ndiag,'(a)')' mksrf_fdynuse = '//trim(mksrf_fdynuse)
       end if
       write(ndiag,*)
       write(ndiag,'(a)')'Model grid configuration variables'  
       write(ndiag,'(a)')' mksrf_fgrid_mesh = '//trim(mksrf_fgrid_mesh)
       write(ndiag,'(a,i8)')' nlon= ',mksrf_fgrid_mesh_nx
       write(ndiag,'(a,i8)')' nlat= ',mksrf_fgrid_mesh_ny
       write(ndiag,'(a)')' mksrf_gridtype = '//trim(mksrf_gridtype)
       write(ndiag,*)
       write(ndiag,'(a)')'Output configuration variables'  
       if (outnc_1d) then
          write(ndiag,'(a)')' output file will be 1d format'
       else
          write(ndiag,'(a)')' fsurdat is 2d lat/lon grid'
       end if
       if ( outnc_large_files ) then
          write(ndiag,'(a)')' Output file in NetCDF 64-bit large_files format'
       end if
       if ( outnc_double )then
          write(ndiag,'(a)')' Output ALL data in file as 64-bit'
       end if
       if ( outnc_vic )then
          write(ndiag,'(a)')' Output VIC fields'
       end if
       if ( outnc_3dglc )then
          write(ndiag,'(a)')' Output optional 3D glacier fields (mostly used for verification of the glacier model)'
       end if
       if ( outnc_3dglc )then
          write(ndiag,'(a)')' Output optional 3D glacier fields (mostly used for verification of the glacier model)'
       end if
       if ( all_urban )then
          write(ndiag,'(a)') ' Output ALL data in file as 100% urban'
       end if
       if ( no_inlandwet )then
          write(ndiag,'(a)') ' Set wetland to 0% over land'
       end if
       if (all_veg) then
          write(ndiag,'(a)') ' Output ALL data in file as 100% vegetated'
       end if
       if (urban_skip_abort_on_invalid_data_check) then
          write(ndiag, '(a)') " WARNING: aborting on invalid data check in urban has been disabled!"
          write(ndiag, '(a)') " WARNING: urban data may be invalid!"
       end if
    end if

  end subroutine write_namelist_input

end module mkvarctl
