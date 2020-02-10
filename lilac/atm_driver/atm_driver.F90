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
  use netcdf      , only : nf90_def_dim, nf90_def_var, nf90_put_var
  use netcdf      , only : nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, nf90_get_var
  use mpi         , only : MPI_COMM_WORLD, MPI_COMM_NULL, MPI_Init, MPI_FINALIZE, MPI_SUCCESS
  use mpi         , only : MPI_GATHER, MPI_DOUBLE
  use lilac_mod   , only : lilac_init, lilac_run, lilac_final
  use lilac_atmcap, only : lilac_atmcap_atm2lnd, lilac_atmcap_lnd2atm
  use shr_cal_mod , only : shr_cal_date2ymd
  use shr_sys_mod , only : shr_sys_abort

  implicit none

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

  namelist /atm_driver_input/ caseid, atm_mesh_file, atm_global_nx, atm_global_ny, &
       atm_calendar, atm_timestep, &
       atm_start_year, atm_start_mon, atm_start_day, atm_start_secs, &
       atm_stop_year, atm_stop_mon, atm_stop_day, atm_stop_secs, atm_starttype
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
  call lilac_init(comp_comm, atm_global_index, atm_lons, atm_lats, &
       atm_global_nx, atm_global_ny, atm_calendar, atm_timestep, &
       atm_start_year, atm_start_mon, atm_start_day, atm_start_secs, &
       atm_starttype)

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

     if ( atm_stop_year /= atm_restart_year) then
        write(6,*)'atm_stop_year = ',atm_stop_year,'atm_restart_year = ',atm_restart_year 
        call shr_sys_abort('not supporting restart and stop years to be different')
     else if (atm_stop_mon  /= atm_restart_mon) then
        write(6,*)'atm_stop_mon = ',atm_stop_mon,'atm_restart_mon = ',atm_restart_mon 
        call shr_sys_abort('not supporting restart and stop months to be different')
     else if (atm_stop_secs /= 0 .or. atm_restart_secs /= 0) then
        write(6,*)'atm_stop_secs = ',atm_stop_secs,'atm_restart_secs = ',atm_restart_secs 
        call shr_sys_abort('not supporting restart and stop secs to be nonzero')
     else
        atm_nsteps = ((atm_stop_day - atm_restart_day) * 86400.) / atm_timestep
     end if

  end if

  do nstep = 1,atm_nsteps
     ! fill in the dataptr values in atm2lnd type in lilac_atmcap
     call atm_driver_to_lilac (atm_lons, atm_lats, nstep, atm_nsteps)

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
  subroutine atm_driver_to_lilac (lon, lat, nstep, atm_nsteps)

    ! input/output variables
    real*8, intent(in) :: lon(:)
    real*8, intent(in) :: lat(:)
    integer, intent(in) :: nstep      ! current step number
    integer, intent(in) :: atm_nsteps ! total number of steps in simulation

    ! local variables
    integer             :: lsize
    real*8              :: time_midpoint
    real*8              :: time_perturbation
    real*8, allocatable :: space_time_perturbation(:)
    real*8, allocatable :: data(:)
    integer             :: i
    integer             :: i_local
    ! --------------------------------------------------------

    lsize = size(lon)
    allocate(space_time_perturbation(lsize))
    allocate(data(lsize))

    ! The time perturbation will range from about -0.5 to 0.5
    time_midpoint = atm_nsteps / 2.d0
    time_perturbation = 0.5d0 * (nstep - time_midpoint)/time_midpoint
    space_time_perturbation(:) = time_perturbation + lat(:)*0.01d0 + lon(:)*0.01d0

    ! We don't have a good way to set a land mask / fraction in this demo driver. Since it
    ! is okay for the atmosphere to call a point ocean when CTSM calls it land, but not
    ! the reverse, here we call all points ocean. In a real atmosphere, the atmosphere
    ! should set landfrac to > 0 for any point for which it needs land input, to ensure
    ! that CTSM is running over all of the necessary points. Note that this landfrac
    ! variable doesn't actually impact the running of CTSM, but it is used for
    ! consistency checking.
    data(:) = 0.d0
    call lilac_atmcap_atm2lnd('Sa_landfrac', data)

    ! In the following, try to have each field have different values, in order to catch
    ! mis-matches (e.g., if foo and bar were accidentally swapped in CTSM, we couldn't
    ! catch that if they both had the same value).

    data(:) = 30.0d0 + space_time_perturbation(:)
    call lilac_atmcap_atm2lnd('Sa_z', data)

    data(:) = 10.0d0 + space_time_perturbation(:)
    call lilac_atmcap_atm2lnd('Sa_topo', data)

    data(:) = 20.0d0 + space_time_perturbation(:)
    call lilac_atmcap_atm2lnd('Sa_u', data)

    data(:) = 40.0d0 + space_time_perturbation(:)
    call lilac_atmcap_atm2lnd('Sa_v', data)

    data(:) = 280.1d0 + space_time_perturbation(:)
    call lilac_atmcap_atm2lnd('Sa_ptem', data)

    data(:) = 100100.0d0 + space_time_perturbation(:)
    call lilac_atmcap_atm2lnd('Sa_pbot', data)

    data(:) = 280.0d0 + space_time_perturbation(:)
    call lilac_atmcap_atm2lnd('Sa_tbot', data)

    data(:) = 0.0004d0 + space_time_perturbation(:)*1.0e-8
    call lilac_atmcap_atm2lnd('Sa_shum', data)

    data(:) = 200.0d0 + space_time_perturbation(:)
    call lilac_atmcap_atm2lnd('Faxa_lwdn', data)

    data(:) = 1.0d-8 +  space_time_perturbation(:)*1.0e-8
    call lilac_atmcap_atm2lnd('Faxa_rainc', data)

    data(:) = 2.0d-8 +  space_time_perturbation(:)*1.0e-8
    call lilac_atmcap_atm2lnd('Faxa_rainl', data)

    data(:) = 1.0d-9 +  space_time_perturbation(:)*1.0e-9
    call lilac_atmcap_atm2lnd('Faxa_snowc', data)

    data(:) = 2.0d-9 +  space_time_perturbation(:)*1.0e-9
    call lilac_atmcap_atm2lnd('Faxa_snowl', data)

    data(:) = 100.0d0 + space_time_perturbation(:)
    call lilac_atmcap_atm2lnd('Faxa_swndr', data)

    data(:) = 50.0d0 + space_time_perturbation(:)
    call lilac_atmcap_atm2lnd('Faxa_swvdr', data)

    data(:) = 25.0d0 + space_time_perturbation(:)
    call lilac_atmcap_atm2lnd('Faxa_swndf', data)

    data(:) = 45.0d0 + space_time_perturbation(:)
    call lilac_atmcap_atm2lnd('Faxa_swvdf', data)

  end subroutine atm_driver_to_lilac

  !========================================================================
  subroutine write_lilac_to_atm_driver_fields(caseid, nlocal, atm_global_nx, atm_global_ny, masterproc)

    ! Fetch lnd2atm fields from LILAC and write them out.
    !
    ! This should only be called once, at the end of the run. (Calling it multiple times
    ! will lead to the output file being overwritten.)

    ! input/output variables
    character(len=*), intent(in) :: caseid
    integer, intent(in) :: nlocal
    integer, intent(in) :: atm_global_nx
    integer, intent(in) :: atm_global_ny
    logical, intent(in) :: masterproc

    ! local variables
    integer, parameter :: field_name_len = 64
    integer :: ierr
    integer :: ncid
    integer :: dimid_x
    integer :: dimid_y
    integer :: nglobal
    integer :: i
    integer, allocatable :: varids(:)
    character(len=field_name_len) :: field_name
    real*8, allocatable :: data(:)
    real*8, allocatable :: data_global(:)
    real*8, allocatable :: data_2d(:,:)

    character(len=field_name_len), parameter :: fields(23) = [character(len=field_name_len) :: &
         'Sl_t', 'Sl_tref', 'Sl_qref', 'Sl_avsdr', 'Sl_anidr', 'Sl_avsdf', 'Sl_anidf', &
         'Sl_snowh', 'Sl_u10', 'Sl_fv', 'Sl_ram1', 'Sl_z0m', &
         'Fall_taux', 'Fall_tauy', 'Fall_lat', 'Fall_sen', 'Fall_lwup', 'Fall_evap', 'Fall_swnet', &
         'Fall_flxdst1', 'Fall_flxdst2', 'Fall_flxdst3', 'Fall_flxdst4']
    ! --------------------------------------------

    if (masterproc) then
       ! Use an arbitrary time rather than trying to figure out the correct time stamp. This
       ! works because this subroutine is only called once, at the end of the run
       ierr = nf90_create(trim(caseid)//'.atm.h0.0001-01.nc', nf90_clobber, ncid)
       if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_create atm driver output file')

       ierr = nf90_def_dim(ncid, 'atm_nx', atm_global_nx, dimid_x)
       if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_def_dim nx atm driver output file')
       ierr = nf90_def_dim(ncid, 'atm_ny', atm_global_ny, dimid_y)
       if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_def_dim ny atm driver output file')

       allocate(varids(size(fields)))
       do i = 1, size(fields)
          field_name = fields(i)
          ierr = nf90_def_var(ncid, field_name, nf90_double, [dimid_x, dimid_y], varids(i))
          if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_def_var atm driver output file: '//trim(field_name))
       end do

       ierr = nf90_enddef(ncid)
       if (ierr /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_enddef atm driver output file')
    end if

    allocate(data(nlocal))
    nglobal = atm_global_nx * atm_global_ny
    if (masterproc) then
       allocate(data_global(nglobal))
       allocate(data_2d(atm_global_nx, atm_global_ny))
    else
       allocate(data_global(1))
    end if

    do i = 1, size(fields)
       field_name = fields(i)
       call lilac_atmcap_lnd2atm(field_name, data)

       ! Because of the way we set up the decomposition, we can use a simple mpi_gather
       ! without needing to worry about any rearrangement, and points will appear in the
       ! correct order on the master proc. Specifically, we rely on the fact that the
       ! first points are assigned to task 0, then the next points to task 1, etc.
       call mpi_gather(data, size(data), mpi_double, data_global, nglobal, mpi_double, 0, &
            mpi_comm_world, ierr)
       if (ierr .ne. MPI_SUCCESS) then
          call shr_sys_abort(' ERROR in mpi_gather for ' // trim(field_name))
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
