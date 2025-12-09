module mkinputMod

  !-----------------------------------------------------------------------
  ! Module containing input namelist settings
  !-----------------------------------------------------------------------

  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_kind_mod , only : CS => shr_kind_CS, CL => shr_kind_CL, CX => shr_kind_CX
  use shr_sys_mod  , only : shr_sys_abort
  use mkvarctl

  implicit none
  private

#include <mpif.h>

  ! ---------------------
  ! routines
  ! ---------------------

  public :: read_namelist_input      ! Read the input control namelist
  public :: bcast_namelist_input     ! Broadcast the namelist to all processors
  public :: check_namelist_input     ! Check the input control namelist for errors
  public :: write_namelist_input     ! Write information on the control namelist

  ! ---------------------
  ! variables
  ! ---------------------

  character(CL) , public    :: fgrddat  ! grid data file
  character(CL) , public    :: fsurdat  ! output surface data file name (if blank, do not output a surface dataset)
  character(CL) , public    :: fsurlog  ! output surface log file name
  character(CL) , public    :: fdyndat  ! dynamic landuse data file name
  character(CL) , public    :: fhrvname ! generic harvest filename
  character(CL) , public    :: furbname ! generic transient urban land cover filename
  character(CL) , public    :: flakname ! generic lake filename

  character(CS) , public    :: mksrf_grid_name              = ' ' ! Name of this grid
  integer       , public    :: grid_size                          ! Number of columne in the grid

  character(CX) , public    :: mksrf_fgrid_mesh             = ' ' ! land grid file name to use
  integer       , public    :: mksrf_fgrid_mesh_nx          = -999
  integer       , public    :: mksrf_fgrid_mesh_ny          = -999

  character(CX) , public    :: mksrf_fvegtyp                = ' ' ! vegetation data file name
  character(CX) , public    :: mksrf_fvegtyp_mesh           = ' ' ! vegetation mesh file name

  character(CX) , public    :: mksrf_fhrvtyp                = ' ' ! harvest data file name
  character(CX) , public    :: mksrf_fhrvtyp_mesh           = ' ' ! harvest mesh file name

  character(CX) , public    :: mksrf_fsoicol                = ' ' ! soil color data file name
  character(CX) , public    :: mksrf_fsoicol_mesh           = ' ' ! soil color mesh file name

  character(CX) , public    :: mksrf_fsoitex                = ' ' ! soil texture mapunit data file name
  character(CX) , public    :: mksrf_fsoitex_lookup         = ' ' ! soil texture lookup data file name
  character(CX) , public    :: mksrf_fsoitex_mesh           = ' ' ! soil texture mesh file name

  character(CX) , public    :: mksrf_fmax                   = ' ' ! fmax data file name
  character(CX) , public    :: mksrf_fmax_mesh              = ' ' ! fmax mesh file name

  character(CX) , public    :: mksrf_fsoildepth             = ' ' ! soil depth file name
  character(CX) , public    :: mksrf_fsoildepth_mesh        = ' ' ! soil depth file name

  character(CX) , public    :: mksrf_fabm                   = ' ' ! ag fire peak month and
  character(CX) , public    :: mksrf_fabm_mesh              = ' ' ! ag fire peak month and

  character(CX) , public    :: mksrf_fpeat                  = ' ' ! peatlands and
  character(CX) , public    :: mksrf_fpeat_mesh             = ' ' ! peatlands and

  character(CX) , public    :: mksrf_fgdp                   = ' ' ! gdp data file names
  character(CX) , public    :: mksrf_fgdp_mesh              = ' ' ! gdp mesh file names

  character(CX) , public    :: mksrf_fpctlak                = ' ' ! percent lake data file name
  character(CX) , public    :: mksrf_fpctlak_mesh           = ' ' ! percent lake file name
  character(CX) , public    :: mksrf_flakdep                = ' ' ! lake depth data file name
  character(CX) , public    :: mksrf_flakdep_mesh           = ' ' ! lake depth file name

  character(CX) , public    :: mksrf_fwetlnd                = ' ' ! inland wetlands data file name
  character(CX) , public    :: mksrf_fwetlnd_mesh           = ' ' ! inland wetlands mesh file name

  character(CX) , public    :: mksrf_furban                 = ' ' ! urban data file name
  character(CX) , public    :: mksrf_furban_mesh            = ' ' ! urban mesh file name

  character(CX) , public    :: mksrf_fglacier               = ' ' ! glacier data file name
  character(CX) , public    :: mksrf_fglacier_mesh          = ' ' ! glacier mesh file name

  character(CX) , public    :: mksrf_fglacierregion         = ' ' ! glacier region data file name
  character(CX) , public    :: mksrf_fglacierregion_mesh    = ' ' ! glacier region mesh file name

  character(CX) , public    :: mksrf_furbtopo               = ' ' ! urban topography data file name
  character(CX) , public    :: mksrf_furbtopo_mesh          = ' ' ! urban topography mesh file name

  character(CX) , public    :: mksrf_flai                   = ' ' ! lai data filename
  character(CX) , public    :: mksrf_flai_mesh              = ' ' ! lai mesh filename

  character(CX) , public    :: mksrf_fdynuse                = ' ' ! ascii file containing names of dynamic land use files
  character(CX) , public    :: mksrf_fdynuse_mesh           = ' ' ! ascii file containing names of dynamic land use files

  character(CX) , public    :: mksrf_fvocef                 = ' ' ! VOC Emission Factor data file name
  character(CX) , public    :: mksrf_fvocef_mesh            = ' ' ! VOC Emission Factor mesh file name

  character(CX) , public    :: mksrf_ftopostats             = ' ' ! topography statistics data file name
  character(CX) , public    :: mksrf_ftopostats_mesh        = ' ' ! topography statistics mesh file name
  character(CX) , public    :: mksrf_ftopostats_override    = ' ' ! read STD_ELEV and SLOPE from this file

  character(CX) , public    :: mksrf_fvic                   = ' ' ! VIC parameters data file name
  character(CX) , public    :: mksrf_fvic_mesh              = ' ' ! VIC parameters mesh file name

  character(CX) , public    :: mksrf_irrig                  = ' ' ! TODO: should this namelist be here?
  character(CX) , public    :: mksrf_irrig_mesh             = ' ' ! TODO: should this namelist be here?

  character(CS) , public    :: gitdescribe                  = ' ' ! Description of model version from git
  character(CS) , public    :: logname                      = ' ' ! user name
  character(CS) , public    :: hostname                     = ' ' ! machine name

  logical       , public    :: create_esmf_pet_files = .false.    ! Always create ESMF PET files (.false. if only on error)
  logical       , public    :: urban_skip_abort_on_invalid_data_check

  character(len=*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine read_namelist_input()

    ! Read in input namelist

    ! local variables
    integer :: ier
    integer :: k
    integer :: fileunit
    logical :: lexist
    character(len=*), parameter :: subname = 'read_namelist_input'
    ! ------------------------------------------------------------

    namelist /mksurfdata_input/      &
         mksrf_grid_name,           &
         mksrf_fvegtyp,             &
         mksrf_fvegtyp_mesh,        &
         mksrf_fhrvtyp,             &
         mksrf_fhrvtyp_mesh,        &
         mksrf_fsoitex,             &
         mksrf_fsoitex_lookup,      &
         mksrf_fsoitex_mesh,        &
         mksrf_fsoicol,             &
         mksrf_fsoicol_mesh,        &
         mksrf_fvocef,              &
         mksrf_fvocef_mesh,         &
         mksrf_fpctlak,             &
         mksrf_fpctlak_mesh,        &
         mksrf_flakdep,             &
         mksrf_flakdep_mesh,        &
         mksrf_fwetlnd,             &
         mksrf_fwetlnd_mesh,        &
         mksrf_fglacier,            &
         mksrf_fglacier_mesh,       &
         mksrf_fglacierregion,      &
         mksrf_fglacierregion_mesh, &
         mksrf_furbtopo,            &
         mksrf_furbtopo_mesh,       &
         mksrf_fmax,                &
         mksrf_fmax_mesh,           &
         mksrf_furban,              &
         mksrf_furban_mesh,         &
         mksrf_flai,                &
         mksrf_flai_mesh,           &
         mksrf_fdynuse,             &
         mksrf_fdynuse_mesh,        &
         mksrf_fgdp,                &
         mksrf_fgdp_mesh,           &
         mksrf_fpeat,               &
         mksrf_fpeat_mesh,          &
         mksrf_fsoildepth,          &
         mksrf_fsoildepth_mesh,     &
         mksrf_fabm,                &
         mksrf_fabm_mesh,           &
         mksrf_ftopostats,          &
         mksrf_ftopostats_mesh,     &
         mksrf_ftopostats_override, &
         mksrf_fvic,                &
         mksrf_fvic_mesh,           &
         mksrf_fgrid_mesh,          &
         mksrf_fgrid_mesh_nx,       &
         mksrf_fgrid_mesh_ny,       &
         numpft,                    &
         no_inlandwet,              &
         nglcec,                    &
         gitdescribe,               &
         logname,                   &
         hostname,                  &
         outnc_large_files,         &
         outnc_double,              &
         outnc_dims,                &
         outnc_vic,                 &
         outnc_3dglc,               &
         fsurdat,                   &
         fdyndat,                   &
         fsurlog,                   &
         std_elev,                  &
         create_esmf_pet_files,     &
         urban_skip_abort_on_invalid_data_check

    ! Set default namelist values - make these the defaults in gen_mksurfdata_namelist.py
    outnc_large_files = .false.
    outnc_double      = .true.
    outnc_vic         = .false.
    outnc_3dglc       = .false.
    no_inlandwet      = .true.
    urban_skip_abort_on_invalid_data_check = .false.   ! default value for bug work around

    if (root_task) then
       write(ndiag,*) 'Attempting to initialize control settings .....'
    end if

    if (root_task) then
       read(5, nml=mksurfdata_input, iostat=ier)
       if (ier > 0) then
          call shr_sys_abort(subname//' error reading in mksurfdata_input namelist from standard input')
       end if
       grid_size = mksrf_fgrid_mesh_nx * mksrf_fgrid_mesh_ny
       if ( mksrf_grid_name == ' ' )then
          write(mksrf_grid_name,'("Cols",I7.7)') grid_size
       end if
    end if

  end subroutine read_namelist_input

  !===============================================================

  subroutine bcast_namelist_input()

    ! Braodcast the namelist to all processors

    ! local variables
    integer :: ier

    call mpi_bcast (mksrf_fgrid_mesh, len(mksrf_fgrid_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fgrid_mesh_nx, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fgrid_mesh_ny, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (grid_size, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (mksrf_grid_name, len(mksrf_grid_name), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_fvegtyp, len(mksrf_fvegtyp), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fvegtyp_mesh, len(mksrf_fvegtyp_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_fhrvtyp, len(mksrf_fhrvtyp), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fhrvtyp_mesh, len(mksrf_fhrvtyp_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_fsoitex, len(mksrf_fsoitex), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fsoitex_lookup, len(mksrf_fsoitex_lookup), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fsoitex_mesh, len(mksrf_fsoitex_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_fmax, len(mksrf_fmax), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fmax_mesh, len(mksrf_fmax_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_fsoicol, len(mksrf_fsoicol), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fsoicol_mesh, len(mksrf_fsoicol_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_fsoildepth, len(mksrf_fsoildepth), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fsoildepth_mesh, len(mksrf_fsoildepth_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_fabm, len(mksrf_fabm), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fabm_mesh, len(mksrf_fabm_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_fpeat, len(mksrf_fpeat), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fpeat_mesh, len(mksrf_fpeat_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_fgdp, len(mksrf_fgdp), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fgdp_mesh, len(mksrf_fgdp_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_fpctlak, len(mksrf_fpctlak), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fpctlak_mesh, len(mksrf_fpctlak_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_flakdep, len(mksrf_flakdep), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_flakdep_mesh, len(mksrf_flakdep_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_fwetlnd, len(mksrf_fwetlnd), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fwetlnd_mesh, len(mksrf_fwetlnd_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_furban, len(mksrf_furban), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_furban_mesh, len(mksrf_furban_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_furbtopo, len(mksrf_furbtopo), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_furbtopo_mesh, len(mksrf_furbtopo_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_fglacier, len(mksrf_fglacier), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fglacier_mesh, len(mksrf_fglacier_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_fglacierregion, len(mksrf_fglacierregion), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fglacierregion_mesh, len(mksrf_fglacierregion_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_flai, len(mksrf_flai), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_flai_mesh, len(mksrf_flai_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_fvocef, len(mksrf_fvocef), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fvocef_mesh, len(mksrf_fvocef_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_ftopostats, len(mksrf_ftopostats), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_ftopostats_mesh, len(mksrf_ftopostats_mesh), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_ftopostats_override, len(mksrf_ftopostats_override), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_fvic, len(mksrf_fvic), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fvic_mesh, len(mksrf_fvic_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (mksrf_fdynuse, len(mksrf_fdynuse), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fdynuse_mesh, len(mksrf_fdynuse_mesh), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (fsurdat, len(fsurdat), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fdyndat, len(fdyndat), MPI_CHARACTER, 0, mpicom, ier)

    call mpi_bcast (outnc_dims, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (outnc_large_files, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (outnc_double, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (outnc_1d, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (outnc_vic, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (outnc_3dglc, 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (urban_skip_abort_on_invalid_data_check, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (create_esmf_pet_files, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (numpft, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (no_inlandwet, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (std_elev, 1, MPI_REAL, 0, mpicom, ier)

    call mpi_bcast (gitdescribe, len(gitdescribe), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (logname, len(logname), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hostname, len(hostname), MPI_CHARACTER, 0, mpicom, ier)

  end subroutine bcast_namelist_input

  !===============================================================
  subroutine check_namelist_input()

    ! error check on namelist input
    if (mksrf_fgrid_mesh /= ' ')then
       fgrddat = mksrf_fgrid_mesh
    else
       call shr_sys_abort(" must specify mksrf_fgrid_mesh")
    endif

    if (nglcec <= 0) then
       call shr_sys_abort('nglcec must be at least 1')
    end if

    ! Reorder user's mksrf_fgrid_mesh_nx and mksrf_fgrid_mesh_ny for 1D
    ! cases if they set model-mesh-nx = 1 instead of model-mesh-ny = 1.
    ! This way the follow-up if-statements works.
    if (mksrf_fgrid_mesh_nx == 1) then
       mksrf_fgrid_mesh_nx = mksrf_fgrid_mesh_ny
       mksrf_fgrid_mesh_ny = 1
       if (root_task) then
          write(ndiag,'(a)') 'WARNING: The code reversed your mksrf_fgrid_mesh_nx and mksrf_fgrid_mesh_ny to '
          write(ndiag,*) 'the expected by the code ', mksrf_fgrid_mesh_nx, mksrf_fgrid_mesh_ny
       end if
    end if

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

    ! Write out namelist input to log file
    ! Note - need to call this after ndiag has been set

    ! local variables
    integer :: k
    ! ------------------------------------------------------------

    if (root_task) then
       write(ndiag,'(a)')'Grid_name:                   '//trim(mksrf_grid_name)
       write(ndiag,'(a)')'Input rawdata files and corresponding meshes'
       write(ndiag,'(a)')' PFTs from:                  '//trim(mksrf_fvegtyp)
       write(ndiag,'(a)')' mesh for pft                '//trim(mksrf_fvegtyp_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' percent lake from:          '//trim(mksrf_fpctlak)
       write(ndiag,'(a)')' mesh for percent lake       '//trim(mksrf_fpctlak_mesh)
       write(ndiag,'(a)')' lake depth from:            '//trim(mksrf_flakdep)
       write(ndiag,'(a)')' mesh for lake depth         '//trim(mksrf_flakdep_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' inland wetland from:        '//trim(mksrf_fwetlnd)
       write(ndiag,'(a)')' mesh for wetland            '//trim(mksrf_fwetlnd_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' soil texture mapunits from: '//trim(mksrf_fsoitex)
       write(ndiag,'(a)')' soil texture (sand/clay, orgc) lookup: '//trim(mksrf_fsoitex_lookup)
       write(ndiag,'(a)')' mesh for soil texture       '//trim(mksrf_fsoitex_mesh)
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
       write(ndiag,'(a)')' mesh for harvest            '//trim(mksrf_fhrvtyp_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' topography statistics from: '//trim(mksrf_ftopostats)
       write(ndiag,'(a)')' mesh for topography stats   '//trim(mksrf_ftopostats_mesh)
       write(ndiag,*)
       write(ndiag,'(a)')' glaciers from:              '//trim(mksrf_fglacier)
       write(ndiag,'(a,i4,a)')'            with:            ',nglcec,' glacier elevation classes'
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
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,'(a)')'Model grid configuration variables'
       write(ndiag,'(a)')' mksrf_fgrid_mesh = '//trim(mksrf_fgrid_mesh)
       write(ndiag,'(a,i8)')' nlon= ',mksrf_fgrid_mesh_nx
       write(ndiag,'(a,i8)')' nlat= ',mksrf_fgrid_mesh_ny
       write(ndiag,'(a,i8)')'Grid_size:                ', grid_size
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)')'Output configuration variables'
       if (outnc_1d) then
          write(ndiag,'(a)')' output file is 1d format'
       else
          write(ndiag,'(a)')' output file is 2d lat/lon format'
       end if
       if ( outnc_large_files ) then
          write(ndiag,'(a)')' Output file in NetCDF 64-bit large_files format'
       end if
       if ( outnc_double )then
          write(ndiag,'(a)')' Output ALL data in file as 64-bit'
       else
          write(ndiag,'(a)')' Output ALL data in file as 32-bit'
       end if
       if ( outnc_vic )then
          write(ndiag,'(a)')' Output VIC fields'
       end if
       if ( outnc_3dglc )then
          write(ndiag,'(a)')' Output optional 3D glacier fields (mostly used for verification of the glacier model)'
       end if
       if (create_esmf_pet_files) then
          write(ndiag,'(a)')' Always output ESMF PET files'
       else
          write(ndiag,'(a)')' Only output ESMF PET files if fatal errors happen in ESMF'
       end if
       if (urban_skip_abort_on_invalid_data_check) then
          write(ndiag, '(a)') " WARNING: aborting on invalid data check in urban has been disabled!"
          write(ndiag, '(a)') " WARNING: urban data may be invalid!"
       end if
       flush(ndiag)
    end if

  end subroutine write_namelist_input

end module mkinputMod
