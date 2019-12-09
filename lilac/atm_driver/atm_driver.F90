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
  
  use mpi         , only : MPI_COMM_WORLD, MPI_COMM_NULL, MPI_Init, MPI_FINALIZE, MPI_SUCCESS
  use lilac_mod   , only : lilac_init, lilac_run, lilac_final
  use lilac_atmcap, only : lilac_atmcap_atm2lnd, lilac_atmcap_lnd2atm
  use shr_sys_mod , only : shr_sys_abort

  implicit none

  integer               :: comp_comm
  integer               :: ierr
  real    , allocatable :: centerCoords(:,:)
  real    , allocatable :: atm_lons(:), atm_lats(:)
  integer , allocatable :: atm_global_index(:)
  integer               :: mytask, ntasks
  integer               :: my_start, my_end
  integer               :: i_local, i_global
  integer               :: nlocal, nglobal
  integer               :: g,i,k         ! indices
  integer               :: fileunit      ! for namelist input
  integer               :: nstep         ! time step counter
  integer               :: atm_nsteps    ! number of time steps of the simulation

  ! Namelist and related variables
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

  namelist /atm_driver_input/ atm_mesh_file, atm_global_nx, atm_global_ny, &
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

  if (mytask == 0 ) then
     print *, "MPI initialization done ..., ntasks=", ntasks
  end if

  !-----------------------------------------------------------------------------
  ! Read in namelist file ...
  !-----------------------------------------------------------------------------

  if (mytask == 0) then
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
  if (mytask == 0 ) then
     print *, " atm_driver mesh file ",trim(atm_mesh_file)
     print *, " atm global nx = ",atm_global_nx
     print *, " atm global nx = ",atm_global_ny
     print *, " atm number of global points in mesh is:", nglobal
  end if

  !-----------------------------------------------------------------------------
  ! atmosphere domain decomposition
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
     atm_lons(i) = real(nint(atm_lons(i))) ! rounding to nearest int
     atm_lats(i) = centerCoords(2,i_global)
     atm_lats(i) = real(nint(atm_lats(i))) ! rounding to nearest int
  end do

  !------------------------------------------------------------------------
  ! Initialize lilac
  !------------------------------------------------------------------------

  if (mytask == 0 ) then
     print *, " initializing lilac with start type ",trim(atm_starttype)
  end if
  call lilac_init(comp_comm, atm_global_index, atm_lons, atm_lats, &
       atm_global_nx, atm_global_ny, atm_calendar, atm_timestep, &
       atm_start_year, atm_start_mon, atm_start_day, atm_start_secs, &
       atm_stop_year, atm_stop_mon, atm_stop_day, atm_stop_secs, atm_starttype)

  !------------------------------------------------------------------------
  ! Run lilac
  !------------------------------------------------------------------------

  if ( atm_stop_year == atm_start_year .and. atm_stop_mon == atm_start_mon .and. &
       atm_stop_secs == atm_start_secs) then
     atm_nsteps = ((atm_stop_day - atm_start_day) * 86400.) / atm_timestep
  else
     call shr_sys_abort('not supporting start and stop years,months and secs to be different')
  end if

  do nstep = 1,atm_nsteps
     ! fill in the dataptr values in atm2lnd type in lilac_atmcap
     call atm_driver_to_lilac (atm_lons, atm_lats)

     if (nstep == atm_nsteps) then
        call lilac_run(restart_alarm_is_ringing=.true., stop_alarm_is_ringing=.true.)
     else
        call lilac_run(restart_alarm_is_ringing=.false., stop_alarm_is_ringing=.false.)
     end if
  end do

  !------------------------------------------------------------------------
  ! Finalize lilac
  !------------------------------------------------------------------------

  call lilac_final( )

  if (mytask == 0 ) then
     print *,  "======================================="
     print *,  " ............. DONE ..................."
     print *,  "======================================="
  end if

!=======================================================
contains
!=======================================================

  subroutine read_netcdf_mesh(filename, nglobal)

    use netcdf
    implicit none

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

    if (mytask == 0 ) then
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
  subroutine atm_driver_to_lilac (lon, lat)

    ! input/output variables
    real, intent(in) :: lon(:)
    real, intent(in) :: lat(:)

    ! local variables
    integer             :: lsize
    real*8, allocatable :: data(:)
    integer             :: i
    integer             :: i_local
    ! --------------------------------------------------------

    lsize = size(lon)
    allocate(data(lsize))

    data(:) = 30.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atmcap_atm2lnd('Sa_z', data)

    data(:) = 10.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atmcap_atm2lnd('Sa_topo', data)

    data(:) = 20.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atmcap_atm2lnd('Sa_u', data)

    data(:) = 40.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atmcap_atm2lnd('Sa_v', data)

    data(:) = 280.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atmcap_atm2lnd('Sa_ptem', data)

    data(:) = 100100.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atmcap_atm2lnd('Sa_pbot', data)

    data(:) = 280.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atmcap_atm2lnd('Sa_tbot', data)

    data(:) = 0.0004d0   !+(lat(:)*0.01d0 + lon(:)*0.01d0)*1.0e-8
    call lilac_atmcap_atm2lnd('Sa_shum', data)

    data(:) = 200.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atmcap_atm2lnd('Faxa_lwdn', data)

    data(:) = 0.0d0
    call lilac_atmcap_atm2lnd('Faxa_rainc', data)

    data(:) = 3.0d-8 +  (lat(:)*0.01d0 + lon(:)*0.01d0)*1.0e-8
    call lilac_atmcap_atm2lnd('Faxa_rainl', data)

    data(:) = 1.0d-8 +  (lat(:)*0.01d0 + lon(:)*0.01d0)*1.0e-8
    call lilac_atmcap_atm2lnd('Faxa_snowc', data)

    data(:) = 2.0d-8 +  (lat(:)*0.01d0 + lon(:)*0.01d0)*1.0e-8
    call lilac_atmcap_atm2lnd('Faxa_snowl', data)

    data(:) = 100.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atmcap_atm2lnd('Faxa_swndr', data)

    data(:) = 50.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atmcap_atm2lnd('Faxa_swvdr', data)

    data(:) = 20.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atmcap_atm2lnd('Faxa_swndf', data)

    data(:) = 40.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atmcap_atm2lnd('Faxa_swvdf', data)

  end subroutine atm_driver_to_lilac

  !========================================================================
  subroutine lilac_to_atm_driver

    ! local variables
    integer :: lsize
    real*8, allocatable :: data(:)
    ! --------------------------------------------

    lsize = size(atm_global_index)
    allocate(data(lsize))

    call lilac_atmcap_lnd2atm('Sl_lfrin' , data)
    call lilac_atmcap_lnd2atm('Sl_t'     , data)
    call lilac_atmcap_lnd2atm('Sl_tref'  , data)
    call lilac_atmcap_lnd2atm('Sl_qref'  , data)
    call lilac_atmcap_lnd2atm('Sl_avsdr' , data)
    call lilac_atmcap_lnd2atm('Sl_anidr' , data)
    call lilac_atmcap_lnd2atm('Sl_avsdf' , data)
    call lilac_atmcap_lnd2atm('Sl_anidf' , data)
    call lilac_atmcap_lnd2atm('Sl_snowh' , data)
    call lilac_atmcap_lnd2atm('Sl_u10'   , data)
    call lilac_atmcap_lnd2atm('Sl_fv'    , data)
    call lilac_atmcap_lnd2atm('Sl_ram1'  , data)

  end subroutine lilac_to_atm_driver

end program
