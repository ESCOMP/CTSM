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
  !   ESMF lilac_atmcap            ESMF land cap     ESMF river cap
  !                                        |                 |
  !                                       CTSM          Mizzouroute...  
  !----------------------------------------------------------------------------
  
  use lilac_mod   , only : lilac_init, lilac_run, lilac_final
  use lilac_utils , only : lilac_atm2lnd, lilac_lnd2atm, gindex_atm
  use mpi         , only : MPI_COMM_WORLD, MPI_COMM_NULL, MPI_Init, MPI_FINALIZE, MPI_SUCCESS

  implicit none

  integer                :: comp_comm
  integer                :: ierr
  real    , allocatable  :: centerCoords(:,:)
  real    , allocatable  :: lon(:), lat(:)
  integer                :: mytask, ntasks
  integer                :: my_start, my_end
  integer                :: i_local, i_global
  integer                :: nlocal, nglobal
  integer                :: start_time               !-- start_time    start time
  integer                :: end_time                 !-- end_time      end time
  integer                :: curr_time                !-- cur_time      current time
  integer                :: itime_step               !-- itime_step    counter of time steps
  integer                :: g,i,k                    !-- indices
  character(len=128)     :: filename
  !------------------------------------------------------------------------

  start_time = 1
  end_time   = 48

  !-----------------------------------------------------------------------------
  ! Initiallize MPI
  !-----------------------------------------------------------------------------

  write(*, *) "MPI initialization starts ..."

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
  ! Read mesh file to get number of points (n_points)
  !-----------------------------------------------------------------------------
  filename = '/glade/p/cesmdata/cseg/inputdata/share/meshes/fv4x5_050615_polemod_ESMFmesh.nc'
  call read_netcdf_mesh(filename, nglobal)
  if (mytask == 0 ) then
     print *, "number of global points is is:", nglobal
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

  allocate(gindex_atm(nlocal))

  i_global = my_start
  do i_local = 1, nlocal
     gindex_atm(i_local) = i_global
     i_global = i_global + 1
  end do

  !------------------------------------------------------------------------
  ! Initialize lilac
  !------------------------------------------------------------------------

  call lilac_init(nlocal)

  !------------------------------------------------------------------------
  ! Fill in atm2lnd type pointer data
  !------------------------------------------------------------------------

  ! first determine lats and lons
  allocate(lon(nlocal))
  allocate(lat(nlocal))
  do i = 1,nlocal
     i_global = gindex_atm(i)
     lon(i) = centerCoords(1,i_global)
     lon(i) = real(nint(lon(i))) ! rounding to nearest int
     lat(i) = centerCoords(2,i_global)
     lat(i) = real(nint(lat(i))) ! rounding to nearest int
  end do

  ! now fill in the dataptr values
  call atm_to_lilac (lon, lat)

  !------------------------------------------------------------------------
  ! Run lilac
  !------------------------------------------------------------------------

  itime_step = 1
  do curr_time = start_time, end_time
     call lilac_run( )
     itime_step = itime_step + 1
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

  !=====================
contains
  !=====================

  subroutine read_netcdf_mesh(filename, nglobal)

    use netcdf
    implicit none

    !  input/output variables
    character(*) , intent(in)  :: filename
    integer      , intent(out) :: nglobal

    !  local Variables
    integer :: idfile
    integer :: ierror
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
    ierror  = nf90_open(filename, NF90_NOWRITE, idfile)
    call nc_check_err(ierror, "opening file", filename)

    ! Get the dimid of  dimensions
    ierror  = nf90_inq_dimid(idfile, 'elementCount', dimid_elem)
    call nc_check_err(ierror, "inq_dimid elementCount", filename)
    ierror  = nf90_inq_dimid(idfile, 'coordDim', dimid_coordDim)
    call nc_check_err(ierror, "coordDim", filename)

    ! Inquire dimensions based on their dimeid(s)
    ierror = nf90_inquire_dimension(idfile, dimid_elem, string, nelem)
    call nc_check_err(ierror, "inq_dim elementCount", filename)
    ierror = nf90_inquire_dimension(idfile, dimid_coordDim, string, coordDim)
    call nc_check_err(ierror, "inq_dim coordDim", filename)

    if (mytask == 0 ) then
       print *,  "======================================="
       print *, "number of elements is : ", nelem
       print *, "coordDim is :", coordDim
       print *,  "======================================="
    end if

    ! Get coordinate values
    allocate (centerCoords(coordDim, nelem))

    ierror = nf90_inq_varid(idfile, 'centerCoords' , idvar_centerCoords)
    call nc_check_err(ierror, "inq_varid centerCoords", filename)
    ierror = nf90_get_var(idfile, idvar_CenterCoords, centerCoords, start=(/1,1/), count=(/coordDim, nelem/))
    call nc_check_err(ierror,"get_var CenterCoords", filename)

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
  subroutine atm_to_lilac (lon, lat)

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
    call lilac_atm2lnd('Sa_z', data)

    data(:) = 10.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atm2lnd('Sa_topo', data)

    data(:) = 20.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atm2lnd('Sa_u', data)

    data(:) = 40.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atm2lnd('Sa_v', data)

    data(:) = 280.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atm2lnd('Sa_ptem', data)

    data(:) = 100100.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atm2lnd('Sa_pbot', data)

    data(:) = 280.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atm2lnd('Sa_tbot', data)

    data(:) = 0.0004d0   !+(lat(:)*0.01d0 + lon(:)*0.01d0)*1.0e-8
    call lilac_atm2lnd('Sa_shum', data)

    data(:) = 200.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atm2lnd('Faxa_lwdn', data)

    data(:) = 0.0d0
    call lilac_atm2lnd('Faxa_rainc', data)

    data(:) = 3.0d-8 +  (lat(:)*0.01d0 + lon(:)*0.01d0)*1.0e-8
    call lilac_atm2lnd('Faxa_rainl', data)

    data(:) = 1.0d-8 +  (lat(:)*0.01d0 + lon(:)*0.01d0)*1.0e-8
    call lilac_atm2lnd('Faxa_snowc', data)

    data(:) = 2.0d-8 +  (lat(:)*0.01d0 + lon(:)*0.01d0)*1.0e-8
    call lilac_atm2lnd('Faxa_snowl', data)

    data(:) = 100.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atm2lnd('Faxa_swndr', data)

    data(:) = 50.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atm2lnd('Faxa_swvdr', data)

    data(:) = 20.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atm2lnd('Faxa_swndf', data)

    data(:) = 40.0d0 + lat(:)*0.01d0 + lon(:)*0.01d0
    call lilac_atm2lnd('Faxa_swvdf', data)

  end subroutine atm_to_lilac

  !========================================================================
  subroutine lilac_to_atm ()

    ! local variables
    integer :: lsize
    real*8, allocatable :: data(:)
    ! --------------------------------------------

    lsize = size(gindex_atm)
    allocate(data(lsize))

    call lilac_lnd2atm('Sl_lfrin' , data)
    call lilac_lnd2atm('Sl_t'     , data)
    call lilac_lnd2atm('Sl_tref'  , data)
    call lilac_lnd2atm('Sl_qref'  , data)
    call lilac_lnd2atm('Sl_avsdr' , data)
    call lilac_lnd2atm('Sl_anidr' , data)
    call lilac_lnd2atm('Sl_avsdf' , data)
    call lilac_lnd2atm('Sl_anidf' , data)
    call lilac_lnd2atm('Sl_snowh' , data)
    call lilac_lnd2atm('Sl_u10'   , data)
    call lilac_lnd2atm('Sl_fv'    , data)
    call lilac_lnd2atm('Sl_ram1'  , data)

  end subroutine lilac_to_atm

end program
