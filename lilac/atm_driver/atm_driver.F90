program atm_driver

  !----------------------------------------------------------------------------
  ! This is a driver for running lilac with CTSM
  ! There can be no references to ESMF in the driver (the host atmosphere cannot
  ! be required to know or use ESMF)
  !
  ! hierarchy seen here:
  !
  !    atm  driver* (WRF, atm_driver, ...)
  !         |
  !         |
  !    lilac (not an ESMF gridded component!)
  !         |     |________________________.____________.......... gridded components
  !         |                              |                 |
  !   ESMF lilac_atmcap            ESMF CTSM cap     ESMF river cap (Mizzouroute, Mosart)
  !----------------------------------------------------------------------------

  use netcdf      , only : nf90_open, nf90_create, nf90_enddef, nf90_close
  use netcdf      , only : nf90_clobber, nf90_write, nf90_nowrite, nf90_noerr, nf90_double
  use netcdf      , only : nf90_def_dim, nf90_def_var, nf90_put_att, nf90_put_var
  use netcdf      , only : nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, nf90_get_var
  use lilac_mod   , only : lilac_init1, lilac_init2, lilac_run, lilac_final
  use lilac_constants , only : fillvalue => lilac_constants_fillvalue
  use ctsm_LilacCouplingFieldIndices
  use ctsm_LilacCouplingFields, only : lilac_atm2lnd, lilac_lnd2atm
  ! A real atmosphere should not use l2a_fields directly. We use it here just for
  ! convenience of writing every lnd -> atm field to a diagnostic output file.
  use ctsm_LilacCouplingFields, only : l2a_fields

  use shr_cal_mod , only : shr_cal_date2ymd
  use shr_sys_mod , only : shr_sys_abort

  implicit none
#include <mpif.h>

  integer               :: comp_comm
  integer               :: ierr
  real    , allocatable :: centerCoords(:,:)
  real*8  , allocatable :: atm_lons(:), atm_lats(:)
  integer , allocatable :: atm_global_index(:)
  integer               :: mytask, ntasks
  logical               :: masterproc
  integer               :: my_start, my_end
  integer               :: i_local, i_global
  integer               :: nlocal, nglobal
  integer               :: g,i,k         ! indices
  integer               :: fileunit      ! for namelist input
  integer               :: nstep         ! time step counter
  integer               :: atm_nsteps    ! number of time steps of the simulation
  integer               :: nsteps_prev_segs ! number of steps run in previous run segments
  integer               :: atm_nsteps_all_segs ! number of time steps of the simulation, across all run segments (see comment below about atm_ndays_all_segs)
  character(len=512)    :: restart_file  ! local path to lilac restart filename
  integer               :: idfile, varid
  integer               :: atm_restart_ymd
  integer               :: atm_restart_year, atm_restart_mon
  integer               :: atm_restart_day, atm_restart_secs

  ! Namelist and related variables
  character(len=512) :: caseid
  character(len=512) :: atm_mesh_file
  integer            :: atm_global_nx
  integer            :: atm_global_ny
  character(len=128) :: atm_calendar
  integer            :: atm_timestep 
  integer            :: atm_start_year ! (yyyy)
  integer            :: atm_stop_year  ! (yyyy)
  integer            :: atm_start_mon  ! (mm)
  integer            :: atm_stop_mon   ! (mm)
  integer            :: atm_start_day
  integer            :: atm_stop_day
  integer            :: atm_start_secs
  integer            :: atm_stop_secs
  character(len=32)  :: atm_starttype
  ! atm_ndays_all_segs is used for generating the fake data. This should give the total
  ! number of days that will be run across all restart segments. If this isn't exactly
  ! right, it's not a big deal: it just means that the fake data won't be symmetrical in
  ! time. It's just important that we have some rough measure of the run length for the
  ! generation of temporal variability, and that this measure be independent of the
  ! number of restart segments that the run is broken into (so that we can get the same
  ! answers in a restart run as in a straight-through run).
  integer :: atm_ndays_all_segs

  namelist /atm_driver_input/ caseid, atm_mesh_file, atm_global_nx, atm_global_ny, &
       atm_calendar, atm_timestep, &
       atm_start_year, atm_start_mon, atm_start_day, atm_start_secs, &
       atm_stop_year, atm_stop_mon, atm_stop_day, atm_stop_secs, atm_starttype, &
       atm_ndays_all_segs
  !------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Initiallize MPI
  !-----------------------------------------------------------------------------

  call MPI_init(ierr)
  if (ierr .ne. MPI_SUCCESS) then
     print *,'Error starting MPI program. Terminating.'
     call MPI_ABORT(MPI_COMM_WORLD, ierr)
  end if

  comp_comm = MPI_COMM_WORLD
  call MPI_COMM_RANK(comp_comm, mytask, ierr)
  call MPI_COMM_SIZE(comp_comm, ntasks, ierr)
  if (mytask == 0) then
     masterproc = .true.
  else
     masterproc = .false.
  end if

  if (masterproc ) then
     print *, "MPI initialization done ..., ntasks=", ntasks
  end if

  !-----------------------------------------------------------------------------
  ! Read in namelist file ...
  !-----------------------------------------------------------------------------

  if (masterproc) then
     print *,"---------------------------------------"
     print *, "MPI initialized in atm_driver ..."
  end if

  ! The following will read this on all processors - might want to do a read just on the 
  ! master processor and broadcast in the future

  open(newunit=fileunit, status="old", file="atm_driver_in")
  read(fileunit, atm_driver_input, iostat=ierr)
  if (ierr > 0) then
     print *, 'Error on reading atm_driver_in' 
     call MPI_ABORT(MPI_COMM_WORLD, ierr)
  end if
  close(fileunit)

  !-----------------------------------------------------------------------------
  ! Read mesh file to get number of points (n_points)
  ! Also read in global number of lons and lats (needed for lilac history output)
  !-----------------------------------------------------------------------------

  call read_netcdf_mesh(atm_mesh_file, nglobal)
  if (atm_global_nx * atm_global_ny /= nglobal) then
     print *, " atm global nx, ny, nglobal = ",atm_global_nx, atm_global_ny, nglobal
     call shr_sys_abort("Error atm_nx*atm_ny is not equal to nglobal")
  end if
  if (masterproc ) then
     print *, " atm_driver mesh file ",trim(atm_mesh_file)
     print *, " atm global nx = ",atm_global_nx
     print *, " atm global nx = ",atm_global_ny
     print *, " atm number of global points in mesh is:", nglobal
  end if

  !-----------------------------------------------------------------------------
  ! atmosphere domain decomposition
  !
  ! Note that other code in this module relies on this simple decomposition, where we
  ! assign the first points to task 0, then the next points to task 1, etc. Specifically,
  ! code in write_lilac_to_atm_driver_fields relies on this decomposition.
  !-----------------------------------------------------------------------------

  nlocal = nglobal / ntasks

  my_start = nlocal*mytask + min(mytask, mod(nglobal, ntasks)) + 1
  ! The first mod(nglobal,ntasks) of ntasks are the ones that have an extra point
  if (mytask < mod(nglobal, ntasks)) then
     nlocal = nlocal + 1
  end if
  my_end = my_start + nlocal - 1

  allocate(atm_global_index(nlocal))

  i_global = my_start
  do i_local = 1, nlocal
     atm_global_index(i_local) = i_global
     i_global = i_global + 1
  end do

  ! first determine lats and lons
  allocate(atm_lons(nlocal))
  allocate(atm_lats(nlocal))
  do i = 1,nlocal
     i_global = atm_global_index(i)
     atm_lons(i) = centerCoords(1,i_global)
     atm_lats(i) = centerCoords(2,i_global)
  end do

  !------------------------------------------------------------------------
  ! Initialize lilac
  !------------------------------------------------------------------------

  if (masterproc ) then
     print *, " initializing lilac with start type ",trim(atm_starttype)
  end if
  call lilac_init1()
  call lilac_init2( &
       mpicom = comp_comm, &
       atm_global_index = atm_global_index, &
       atm_lons = atm_lons, &
       atm_lats = atm_lats, &
       atm_global_nx = atm_global_nx, &
       atm_global_ny = atm_global_ny, &
       atm_calendar = atm_calendar, &
       atm_timestep = atm_timestep, &
       atm_start_year = atm_start_year, &
       atm_start_mon = atm_start_mon, &
       atm_start_day = atm_start_day, &
       atm_start_secs = atm_start_secs, &
       starttype_in = atm_starttype, &
       fields_needed_from_data = [ &
       ! Deliberately excluding bcphidry to test the logic that says that a field should
       ! only be read from data if explicitly requested by the host atmosphere.
       lilac_a2l_Faxa_bcphodry, lilac_a2l_Faxa_bcphiwet, &
       lilac_a2l_Faxa_ocphidry, lilac_a2l_Faxa_ocphodry, lilac_a2l_Faxa_ocphiwet, &
       lilac_a2l_Faxa_dstwet1, lilac_a2l_Faxa_dstdry1, &
       lilac_a2l_Faxa_dstwet2, lilac_a2l_Faxa_dstdry2, &
       lilac_a2l_Faxa_dstwet3, lilac_a2l_Faxa_dstdry3, &
       lilac_a2l_Faxa_dstwet4, lilac_a2l_Faxa_dstdry4])

  !------------------------------------------------------------------------
  ! Run lilac
  !------------------------------------------------------------------------

  ! Assume that will always run for N days (no partial days)

  if (atm_starttype == 'startup') then
     
     if ( atm_stop_year /= atm_start_year) then
        call shr_sys_abort('not supporting start and stop years to be different')
     else if (atm_stop_mon  /= atm_start_mon) then
        call shr_sys_abort('not supporting start and stop months to be different')
     else if (atm_stop_secs /= 0 .or. atm_start_secs /= 0) then
        call shr_sys_abort('not supporting start and stop secs to be nonzero')
     else
        atm_nsteps = ((atm_stop_day - atm_start_day) * 86400.) / atm_timestep
     end if
     nsteps_prev_segs = 0

  else ! continue

     open(newunit=fileunit, file='rpointer.lilac', form='FORMATTED', status='old',iostat=ierr)
     if (ierr < 0) call shr_sys_abort('Error opening rpointer.lilac')
     read(fileunit,'(a)', iostat=ierr) restart_file
     if (ierr < 0) call shr_sys_abort('Error reading rpointer.lilac')
     close(fileunit)
     if (masterproc) then
        print *,'lilac restart_file = ',trim(restart_file)
     end if

     ierr = nf90_open(restart_file, NF90_NOWRITE, idfile)
     if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_open')

     ierr = nf90_inq_varid(idfile, 'curr_ymd', varid)
     if (ierr /= nf90_NoErr) call shr_sys_abort('ERROR: nf90_inq_varid curr_ymd')
     ierr = nf90_get_var(idfile, varid, atm_restart_ymd)
     if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_get_var curr_ymd')

     ierr = nf90_inq_varid(idfile, 'curr_tod', varid)
     if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_inq_varid curr_tod')
     ierr = nf90_get_var(idfile, varid, atm_restart_secs)
     if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_get_var curr_tod')

     ierr = nf90_close(idfile)
     if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_close')

     if (masterproc) then
        print  *,'restart_ymd = ',atm_restart_ymd
     end if
     call shr_cal_date2ymd(atm_restart_ymd, atm_restart_year, atm_restart_mon, atm_restart_day)

     if ( atm_stop_year /= atm_restart_year .or. atm_restart_year /= atm_start_year) then
        write(6,*)'atm_stop_year, atm_restart_year, atm_start_year = ',&
             atm_stop_year, atm_restart_year, atm_start_year
        call shr_sys_abort('not supporting restart, stop and start years to be different')
     else if (atm_stop_mon /= atm_restart_mon .or. atm_restart_mon /= atm_start_mon) then
        write(6,*)'atm_stop_mon, atm_restart_mon, atm_start_mon = ',&
             atm_stop_mon, atm_restart_mon, atm_start_mon
        call shr_sys_abort('not supporting restart, stop and start months to be different')
     else if (atm_stop_secs /= 0 .or. atm_restart_secs /= 0 .or. atm_start_secs /= 0) then
        write(6,*)'atm_stop_secs, atm_restart_secs, atm_start_secs = ',&
             atm_stop_secs, atm_restart_secs, atm_start_secs
        call shr_sys_abort('not supporting restart, stop or start secs to be nonzero')
     else
        atm_nsteps = ((atm_stop_day - atm_restart_day) * 86400.) / atm_timestep
        ! The following calculation of nsteps_prev_segs is why we need to check the start
        ! time in the above error checks.
        nsteps_prev_segs = ((atm_restart_day - atm_start_day) * 86400.) / atm_timestep
     end if

  end if

  atm_nsteps_all_segs = atm_ndays_all_segs * (86400 / atm_timestep)

  do nstep = 1,atm_nsteps
     ! fill in the dataptr in lilac_coupling_fields
     call atm_driver_to_lilac (atm_lons, atm_lats, nstep, nsteps_prev_segs, atm_nsteps_all_segs)

     if (nstep == atm_nsteps) then
        call lilac_run(write_restarts_now=.true., stop_now=.true.)
     else
        call lilac_run(write_restarts_now=.false., stop_now=.false.)
     end if
  end do

  call write_lilac_to_atm_driver_fields( &
       caseid = caseid, &
       nlocal = nlocal, &
       atm_global_nx = atm_global_nx, &
       atm_global_ny = atm_global_ny, &
       ntasks = ntasks, &
       masterproc = masterproc)

  !------------------------------------------------------------------------
  ! Finalize lilac
  !------------------------------------------------------------------------

  call lilac_final( )

  if (masterproc ) then
     print *,  "======================================="
     print *,  " ............. DONE ..................."
     print *,  "======================================="
  end if

  call MPI_finalize(ierr)

!=======================================================
contains
!=======================================================

  subroutine read_netcdf_mesh(filename, nglobal)

    !  input/output variables
    character(*) , intent(in)  :: filename
    integer      , intent(out) :: nglobal

    !  local Variables
    integer :: idfile
    integer :: ierr
    integer :: dimid_elem
    integer :: dimid_coordDim
    integer :: iddim_elem
    integer :: iddim_coordDim
    integer :: idvar_CenterCoords
    integer :: nelem
    integer :: coordDim
    character (len=100) :: string
    !-----------------------------------------------------------------------------

    ! Open mesh file and get the idfile
    ierr  = nf90_open(filename, NF90_NOWRITE, idfile)
    call nc_check_err(ierr, "opening file", filename)

    ! Get the dimid of  dimensions
    ierr  = nf90_inq_dimid(idfile, 'elementCount', dimid_elem)
    call nc_check_err(ierr, "inq_dimid elementCount", filename)
    ierr  = nf90_inq_dimid(idfile, 'coordDim', dimid_coordDim)
    call nc_check_err(ierr, "coordDim", filename)

    ! Inquire dimensions based on their dimeid(s)
    ierr = nf90_inquire_dimension(idfile, dimid_elem, string, nelem)
    call nc_check_err(ierr, "inq_dim elementCount", filename)
    ierr = nf90_inquire_dimension(idfile, dimid_coordDim, string, coordDim)
    call nc_check_err(ierr, "inq_dim coordDim", filename)

    if (masterproc ) then
       print *,  "======================================="
       print *, "number of elements is : ", nelem
       print *, "coordDim is :", coordDim
       print *,  "======================================="
    end if

    ! Get coordinate values
    allocate (centerCoords(coordDim, nelem))

    ierr = nf90_inq_varid(idfile, 'centerCoords' , idvar_centerCoords)
    call nc_check_err(ierr, "inq_varid centerCoords", filename)
    ierr = nf90_get_var(idfile, idvar_CenterCoords, centerCoords, start=(/1,1/), count=(/coordDim, nelem/))
    call nc_check_err(ierr,"get_var CenterCoords", filename)

    ! Close the file
    ierr = nf90_close(idfile)
    if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_close')

    nglobal = nelem

  end subroutine read_netcdf_mesh

  !========================================================================
  subroutine nc_check_err(ierror, description, filename)
    use netcdf
    integer     , intent(in) :: ierror
    character(*), intent(in) :: description
    character(*), intent(in) :: filename

    if (ierror /= nf90_noerr) then
       write (*,'(6a)') 'ERROR ', trim(description),'. NetCDF file : "', trim(filename),&
            '". Error message:', nf90_strerror(ierror)
    endif
  end subroutine nc_check_err

  !========================================================================
  subroutine atm_driver_to_lilac (lon, lat, nstep, nsteps_prev_segs, atm_nsteps_all_segs)

    ! input/output variables
    real*8, intent(in) :: lon(:)
    real*8, intent(in) :: lat(:)
    integer, intent(in) :: nstep      ! current step number
    integer, intent(in) :: nsteps_prev_segs ! number of time steps in previous run segments
    integer, intent(in) :: atm_nsteps_all_segs ! total number of steps in simulation

    ! local variables
    integer             :: lsize
    real*8              :: time_midpoint
    real*8              :: time_perturbation
    real*8, allocatable :: space_perturbation(:)
    real*8, allocatable :: space_time_perturbation(:)
    real*8, allocatable :: data(:)
    integer             :: i
    integer             :: i_local
    ! --------------------------------------------------------

    lsize = size(lon)
    allocate(space_perturbation(lsize))
    allocate(space_time_perturbation(lsize))
    allocate(data(lsize))

    ! The time perturbation will range from about -0.5 to 0.5
    time_midpoint = atm_nsteps_all_segs / 2.d0
    time_perturbation = 0.5d0 * ((nstep + nsteps_prev_segs) - time_midpoint)/time_midpoint
    space_perturbation(:) = lat(:)*0.01d0 + lon(:)*0.01d0
    space_time_perturbation(:) = time_perturbation + space_perturbation(:)

    ! Only set landfrac in the first time step, similar to what most real atmospheres
    ! will probably do.
    !
    ! We don't have a good way to set a land mask / fraction in this demo driver. Since it
    ! is okay for the atmosphere to call a point ocean when CTSM calls it land, but not
    ! the reverse, here we call all points ocean. In a real atmosphere, the atmosphere
    ! should set landfrac to > 0 for any point for which it needs land input, to ensure
    ! that CTSM is running over all of the necessary points. Note that this landfrac
    ! variable doesn't actually impact the running of CTSM, but it is used for
    ! consistency checking.
    if (nstep == 1) then
       data(:) = 0.d0
       call lilac_atm2lnd(lilac_a2l_Sa_landfrac, data)
    end if

    ! In the following, try to have each field have different values, in order to catch
    ! mis-matches (e.g., if foo and bar were accidentally swapped in CTSM, we couldn't
    ! catch that if they both had the same value).

    ! Sa_z is allowed to be time-constant, but we're keeping it time-varying here in
    ! order to test the ability to have an allowed-to-be-time-constant field actually be
    ! time-varying.
    data(:) = 30.0d0 + space_time_perturbation(:)
    call lilac_atm2lnd(lilac_a2l_Sa_z, data)

    ! Use a time-constant topo field (which may be typical of atmospheres), in order to
    ! test the infrastructure that allows fields to be just set once, in the first time
    ! step.
    if (nstep == 1) then
       data(:) = 10.0d0 + space_perturbation(:)
       call lilac_atm2lnd(lilac_a2l_Sa_topo, data)
    end if

    data(:) = 20.0d0 + space_time_perturbation(:)
    call lilac_atm2lnd(lilac_a2l_Sa_u, data)

    data(:) = 40.0d0 + space_time_perturbation(:)
    call lilac_atm2lnd(lilac_a2l_Sa_v, data)

    data(:) = 280.1d0 + space_time_perturbation(:)
    call lilac_atm2lnd(lilac_a2l_Sa_ptem, data)

    data(:) = 100100.0d0 + space_time_perturbation(:)
    call lilac_atm2lnd(lilac_a2l_Sa_pbot, data)

    data(:) = 280.0d0 + space_time_perturbation(:)
    call lilac_atm2lnd(lilac_a2l_Sa_tbot, data)

    data(:) = 0.0004d0 + space_time_perturbation(:)*1.0d-8
    call lilac_atm2lnd(lilac_a2l_Sa_shum, data)

    data(:) = 200.0d0 + space_time_perturbation(:)
    call lilac_atm2lnd(lilac_a2l_Faxa_lwdn, data)

    data(:) = 1.0d-8 +  space_time_perturbation(:)*1.0d-9
    call lilac_atm2lnd(lilac_a2l_Faxa_rainc, data)

    data(:) = 2.0d-8 +  space_time_perturbation(:)*1.0d-9
    call lilac_atm2lnd(lilac_a2l_Faxa_rainl, data)

    data(:) = 1.0d-9 +  space_time_perturbation(:)*1.0d-10
    call lilac_atm2lnd(lilac_a2l_Faxa_snowc, data)

    data(:) = 2.0d-9 +  space_time_perturbation(:)*1.0d-10
    call lilac_atm2lnd(lilac_a2l_Faxa_snowl, data)

    data(:) = 100.0d0 + space_time_perturbation(:)
    call lilac_atm2lnd(lilac_a2l_Faxa_swndr, data)

    data(:) = 50.0d0 + space_time_perturbation(:)
    call lilac_atm2lnd(lilac_a2l_Faxa_swvdr, data)

    data(:) = 25.0d0 + space_time_perturbation(:)
    call lilac_atm2lnd(lilac_a2l_Faxa_swndf, data)

    data(:) = 45.0d0 + space_time_perturbation(:)
    call lilac_atm2lnd(lilac_a2l_Faxa_swvdf, data)

    ! This field has the potential to be read from data. We're setting it here to provide
    ! a test of the logic that says that a field should only be read from data if
    ! explicitly requested by the host atmosphere.
    data(:) = 1.0d-13 + space_time_perturbation(:)*1.0d-14
    call lilac_atm2lnd(lilac_a2l_Faxa_bcphidry, data)

  end subroutine atm_driver_to_lilac

  !========================================================================
  subroutine write_lilac_to_atm_driver_fields(caseid, nlocal, atm_global_nx, &
       atm_global_ny, ntasks, masterproc)

    ! Fetch lnd2atm fields from LILAC and write them out.
    !
    ! This should only be called once, at the end of the run. (Calling it multiple times
    ! will lead to the output file being overwritten.)

    ! input/output variables
    character(len=*), intent(in) :: caseid
    integer, intent(in) :: nlocal
    integer, intent(in) :: atm_global_nx
    integer, intent(in) :: atm_global_ny
    integer, intent(in) :: ntasks
    logical, intent(in) :: masterproc

    ! local variables
    integer :: nfields
    integer :: ierr
    integer :: ncid
    integer :: dimid_x
    integer :: dimid_y
    integer :: nglobal
    integer :: i
    integer, allocatable :: varids(:)
    character(len=:), allocatable :: field_name
    integer, allocatable :: counts(:)
    integer, allocatable :: displacements(:)
    real*8, allocatable :: data(:)
    real*8, allocatable :: data_global(:)
    real*8, allocatable :: data_2d(:,:)

    ! --------------------------------------------

    ! Implementation note: for convenience and ease of maintenance, we directly leverage
    ! l2a_fields in this subroutine, and loop through all available indices in that
    ! list. A real atmosphere should not use that variable directly; instead, it should
    ! use the indices defined in ctsm_LilacCouplingFieldIndices, similarly to what is done
    ! above in atm_driver_to_lilac.

    nfields = l2a_fields%num_fields()

    ! ------------------------------------------------------------------------
    ! Set up output file
    ! ------------------------------------------------------------------------

    if (masterproc) then
       ! Use an arbitrary time rather than trying to figure out the correct time stamp. This
       ! works because this subroutine is only called once, at the end of the run
       ierr = nf90_create(trim(caseid)//'.clm2.lilac_atm_driver_h0.0001-01.nc', nf90_clobber, ncid)
       if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_create atm driver output file')

       ierr = nf90_def_dim(ncid, 'atm_nx', atm_global_nx, dimid_x)
       if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_def_dim nx atm driver output file')
       ierr = nf90_def_dim(ncid, 'atm_ny', atm_global_ny, dimid_y)
       if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_def_dim ny atm driver output file')

       allocate(varids(nfields))
       do i = 1, nfields
          field_name = l2a_fields%get_fieldname(i)
          ierr = nf90_def_var(ncid, field_name, nf90_double, [dimid_x, dimid_y], varids(i))
          if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_def_var atm driver output file: '//trim(field_name))
          ierr = nf90_put_att(ncid, varids(i), '_FillValue', fillvalue)
          if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_put_att atm driver output file: '//trim(field_name))
       end do

       ierr = nf90_enddef(ncid)
       if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_enddef atm driver output file')
    end if

    ! ------------------------------------------------------------------------
    ! Determine number of points on each processor and set up arrays needed for gathering
    ! data to master proc
    ! ------------------------------------------------------------------------

    allocate(data(nlocal))
    nglobal = atm_global_nx * atm_global_ny
    if (masterproc) then
       allocate(counts(ntasks))
       allocate(displacements(ntasks))
       allocate(data_global(nglobal))
       allocate(data_2d(atm_global_nx, atm_global_ny))
    else
       allocate(counts(1))
       allocate(displacements(1))
       allocate(data_global(1))
    end if

    call mpi_gather(nlocal, 1, mpi_int, counts, 1, mpi_int, 0, mpi_comm_world, ierr)
    if (ierr .ne. MPI_SUCCESS) then
       call shr_sys_abort(' ERROR in mpi_gather for counts')
    end if

    if (masterproc) then
       displacements(1) = 0
       do i = 2, ntasks
          displacements(i) = displacements(i-1) + counts(i-1)
       end do
    end if

    ! ------------------------------------------------------------------------
    ! Retrieve data for each field, gather to master and write to file
    ! ------------------------------------------------------------------------

    do i = 1, nfields
       field_name = l2a_fields%get_fieldname(i)
       ! See implementation note above: typically a host atmosphere should NOT loop
       ! through fields, accessing them anonymously by index as is done here. Instead,
       ! typically the host atmosphere would access specific fields using the indices
       ! defined in ctsm_LilacCouplingFieldIndices.
       call lilac_lnd2atm(i, data)

       ! Because of the way we set up the decomposition, we can use a simple mpi_gatherv
       ! without needing to worry about any rearrangement, and points will appear in the
       ! correct order on the master proc. Specifically, we rely on the fact that the
       ! first points are assigned to task 0, then the next points to task 1, etc.
       call mpi_gatherv(data, size(data), mpi_double, data_global, counts, displacements, &
            mpi_double, 0, mpi_comm_world, ierr)
       if (ierr .ne. MPI_SUCCESS) then
          call shr_sys_abort(' ERROR in mpi_gatherv for ' // trim(field_name))
       end if

       if (masterproc) then
          data_2d = reshape(data_global, [atm_global_nx, atm_global_ny])
          ierr = nf90_put_var(ncid, varids(i), data_2d)
          if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_put_var atm driver output file: '//trim(field_name))
       end if
    end do

    if (masterproc) then
       ierr = nf90_close(ncid)
       if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_close atm driver output file')
    end if

  end subroutine write_lilac_to_atm_driver_fields

end program
