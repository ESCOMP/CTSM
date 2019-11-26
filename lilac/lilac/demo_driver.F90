module demo_mod
!----------------------------------------------------------------------------
    use mpi           ,  only : MPI_COMM_WORLD, MPI_COMM_NULL, MPI_Init, MPI_FINALIZE, MPI_SUCCESS
    use spmdMod       ,  only : masterproc
    use lilac_utils , only : atm2lnd_data1d_type , lnd2atm_data1d_type , atm2lnd_data2d_type , atm2lnd_data2d_type , this_clock
    implicit none
    private
    public  :: demo_init
    public  :: read_netcdf_mesh
    integer :: ierr
    integer :: COMP_COMM
    integer :: npts   ! domain global size
    integer :: num_local
    integer                              :: n_points
    real, dimension(:,:), allocatable :: centerCoords
!----------------------------------------------------------------------------
contains
!----------------------------------------------------------------------------
    subroutine demo_init(gindex_atm, atm2lnd, lnd2atm)
        !! TODO: IS THE INTENT CORRECT FOR GINDEX_ATM
        integer , allocatable, intent(inout) :: gindex_atm(:)
        type (atm2lnd_data1d_type), intent(inout)                            :: atm2lnd
        type (lnd2atm_data1d_type), intent(inout)                            :: lnd2atm
        integer                              :: ntasks
        integer                              :: mytask
        character(len=128)                   :: filename
        integer                 :: endc
        !-----------------------------------------------------------------------------
        ! Initiallize MPI
        !-----------------------------------------------------------------------------

        npts = 3312

        write(*, *) "MPI initialization starts ..."

        call MPI_init(ierr)
        COMP_COMM = MPI_COMM_WORLD

        !https://github.com/yudong-tian/LIS-CLM4.5SP/blob/8cec515a628325c73058cfa466db63210cd562ac/xlis-bld/xlis_main.F90
        if (ierr .ne. MPI_SUCCESS) then
            print *,'Error starting MPI program. Terminating.'
            call MPI_ABORT(MPI_COMM_WORLD, ierr)
        end if

        !

        call MPI_COMM_RANK(COMP_COMM, mytask, ierr)
        call MPI_COMM_SIZE(COMP_COMM, ntasks, ierr)

        if (mytask == 0 ) then
            print *, "MPI initialization done ..., ntasks=", ntasks
        end if


        !-----------------------------------------------------------------------------
        ! Read mesh file to get number of points (n_points)
        !-----------------------------------------------------------------------------
        filename = '/glade/p/cesmdata/cseg/inputdata/share/meshes/fv4x5_050615_polemod_ESMFmesh.nc'
        call read_netcdf_mesh(filename, n_points)

        !-----------------------------------------------------------------------------
        ! atmosphere domain decomposition
        !-----------------------------------------------------------------------------

        npts = n_points
        print *, "npts for ", mytask, "is:", npts
        call decompInit_atm( ntasks, mytask, gindex_atm)
        print *, "gindex_atm for ", mytask,"is: ", gindex_atm
        print *, "size gindex_atm for ", mytask,"is: ", size(gindex_atm)

        !-----------------------------------------------------------------------------
        ! allocate and fill in atm2lnd 
        !-----------------------------------------------------------------------------

        endc = npts /ntasks
        call fill_in (atm2lnd, lnd2atm, 1, endc, gindex_atm)
    end subroutine demo_init

    subroutine decompInit_atm( ntasks, mytask,  gindex_atm)

        ! !DESCRIPTION:
        ! !USES:

        ! !ARGUMENTS:
        integer , intent(in)               :: ntasks
        integer , intent(in)               :: mytask
        integer , allocatable, intent(out) :: gindex_atm(:) ! this variable is allocated here, and is assumed to start unallocated
        ! !LOCAL VARIABLES:
        integer                            :: my_start
        integer                            :: my_end
        integer                            :: i_local
        integer                            :: i_global
        !------------------------------------------------------------------------------
        ! create the a global index array for ocean points

        num_local = npts / ntasks

        my_start = num_local*mytask + min(mytask, mod(npts, ntasks)) + 1
        ! The first mod(npts,ntasks) of ntasks are the ones that have an extra point
        if (mytask < mod(npts, ntasks)) then
           num_local = num_local + 1
        end if
        my_end = my_start + num_local - 1

        allocate(gindex_atm(num_local))

        i_global = my_start
        do i_local = 1, num_local
            gindex_atm(i_local) = i_global
            i_global = i_global +1
        end do 

    end subroutine decompInit_atm

    subroutine read_netcdf_mesh(filename, n_points)

        use netcdf
        implicit none

    !
    !  Parameters
    !

    !
    !  Arguments | Global Variables
    !
        character(*)                     , intent(in)  :: filename
        integer , intent(inout)               :: n_points

    !
    !  Local Variables
    !

        integer :: idfile

        integer :: ierror
        integer :: dimid_node
        integer :: dimid_elem
        integer :: dimid_maxnodepe
        integer :: dimid_coordDim

        integer :: iddim_node
        integer :: iddim_elem
        integer :: iddim_maxnodepe
        integer :: iddim_coordDim

        integer :: idvar_nodeCoords
        integer :: idvar_CenterCoords

        character (len=100) :: string


        integer :: nnode
        integer :: nelem
        integer :: maxnodePE
        integer :: coordDim
        real, dimension(:,:), allocatable :: nodeCoords
        !-----------------------------------------------------------------------------
        ! Open mesh file and get the idfile
        ierror  = nf90_open ( filename, NF90_NOWRITE, idfile) ; call nc_check_err(ierror, "opening file", filename)

        ! Get the dimid of  dimensions
        ierror  = nf90_inq_dimid(idfile, 'nodeCount'        , dimid_node ) ; call nc_check_err(ierror, "inq_dimid nodeCount", filename)
        ierror  = nf90_inq_dimid(idfile, 'elementCount'     , dimid_elem ); call nc_check_err(ierror, "inq_dimid elementCount", filename)
        ierror  = nf90_inq_dimid(idfile, 'maxNodePElement'  , dimid_maxnodepe ); call nc_check_err(ierror, "inq_dimid maxNodePElement", filename)
        ierror  = nf90_inq_dimid(idfile, 'coordDim'         , dimid_coordDim  ); call nc_check_err(ierror, "coordDim", filename)

        ! Inquire dimensions based on their dimeid(s)
        ierror = nf90_inquire_dimension(idfile, dimid_node        , string, nnode     ); call nc_check_err(ierror, "inq_dim nodeCount", filename)
        ierror = nf90_inquire_dimension(idfile, dimid_elem        , string, nelem     ); call nc_check_err(ierror, "inq_dim elementCount", filename)
        ierror = nf90_inquire_dimension(idfile, dimid_maxnodepe   , string, maxnodePE ); call nc_check_err(ierror, "inq_dim maxNodePElement", filename)
        ierror = nf90_inquire_dimension(idfile, dimid_coordDim    , string, coordDim  ); call nc_check_err(ierror, "inq_dim coordDim", filename)

        print *,  "======================================="
        print *, "nnode is : ", nnode
        print *, "nelem is : ", nelem
        print *, "coordDim is :", coordDim
        print *,  "======================================="

        allocate (nodeCoords(coordDim, nnode))
        allocate (centerCoords(coordDim, nelem))
        ! Get variable IDs (varid)
        ierror = nf90_inq_varid(idfile, 'nodeCoords'   , idvar_nodeCoords      ); call nc_check_err(ierror, "inq_varid nodeCoords", filename)
        ierror = nf90_inq_varid(idfile, 'centerCoords' , idvar_centerCoords    ); call nc_check_err(ierror, "inq_varid centerCoords", filename)

        ! Get variables values from varids
        ierror = nf90_get_var(idfile, idvar_nodeCoords       , nodeCoords       , start=(/ 1,1/)   , count=(/ coordDim, nnode /)     ); call nc_check_err(ierror,"get_var nodeCoords", filename)
        ierror = nf90_get_var(idfile, idvar_CenterCoords     , centerCoords     , start=(/ 1,1/)   , count=(/ coordDim, nelem /)     ); call nc_check_err(ierror,"get_var CenterCoords", filename)

        !print *, "lons : ",centerCoords(1,:)

        n_points = nelem

    end subroutine read_netcdf_mesh

    subroutine nc_check_err(ierror, description, filename)
    !-------------------------------------------------------------------------------
    !  $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/trunk/src/engines_gpl/wave/packages/data/src/nc_check_err.f90 $
    !!--declarations----------------------------------------------------------------
        use netcdf
        !
        implicit none
    !
    ! Global variables
    !
        integer     , intent(in) :: ierror
        character(*), intent(in) :: description
        character(*), intent(in) :: filename
    !
    ! Local variables
    !
    !
    !     real, parameter :: PI = 3.1415927

    !! executable statements -------------------------------------------------------
    !
        if (ierror /= nf90_noerr) then
            print *,  "ERROR"
            write (*,'(6a)') 'ERROR ', trim(description), '. NetCDF file : "', trim(filename), '". Error message:', nf90_strerror(ierror)
        endif
    end subroutine nc_check_err

    subroutine fill_in (atm2lnd , lnd2atm , begc, endc , gindex_atm) 
        ! !ARGUMENTS:
        type (atm2lnd_data1d_type), intent(inout)                            :: atm2lnd
        type (lnd2atm_data1d_type), intent(inout)                            :: lnd2atm

        integer , intent(in)               :: begc
        integer , intent(in)               :: endc


        real    :: lat
        real    :: lon

        integer , allocatable, intent(in) :: gindex_atm(:)
        !integer :: i 
        integer :: i_local 
        integer :: i_global 
        

        ! tbot is going to be analytical function

        allocate ( atm2lnd%Sa_z       (begc:endc) ) !; atm2lnd%Sa_z       (:) =  30.0d0
        allocate ( atm2lnd%Sa_topo    (begc:endc) ) !; atm2lnd%Sa_topo    (:) =  10.0d0
        allocate ( atm2lnd%Sa_u       (begc:endc) ) !; atm2lnd%Sa_u       (:) =  20.0d0
        allocate ( atm2lnd%Sa_v       (begc:endc) ) !; atm2lnd%Sa_v       (:) =  40.0d0
        allocate ( atm2lnd%Sa_ptem    (begc:endc) ) !; atm2lnd%Sa_ptem    (:) =  280.0d0
        allocate ( atm2lnd%Sa_pbot    (begc:endc) ) !; atm2lnd%Sa_pbot    (:) =  100100.0d0
        allocate ( atm2lnd%Sa_tbot    (begc:endc) ) !; atm2lnd%Sa_tbot    (:) =  280.0
        allocate ( atm2lnd%Sa_shum    (begc:endc) ) !; atm2lnd%Sa_shum    (:) =  0.0004d0

        allocate ( atm2lnd%Faxa_lwdn  (begc:endc) ) !; atm2lnd%Faxa_lwdn  (:) =  200.0d0
        allocate ( atm2lnd%Faxa_rainc (begc:endc) ) !; atm2lnd%Faxa_rainc (:) =  0.0d0
        allocate ( atm2lnd%Faxa_rainl (begc:endc) ) !; atm2lnd%Faxa_rainl (:) =  3.0d-8
        allocate ( atm2lnd%Faxa_snowc (begc:endc) ) !; atm2lnd%Faxa_snowc (:) =  1.0d-8
        allocate ( atm2lnd%Faxa_snowl (begc:endc) ) !; atm2lnd%Faxa_snowl (:) =  2.0d-8

        allocate ( atm2lnd%Faxa_swndr (begc:endc) ) !; atm2lnd%Faxa_swndr (:) =  100.0d0
        allocate ( atm2lnd%Faxa_swvdr (begc:endc) ) !; atm2lnd%Faxa_swvdr (:) =  50.0d0
        allocate ( atm2lnd%Faxa_swndf (begc:endc) ) !; atm2lnd%Faxa_swndf (:) =  20.0d0
        allocate ( atm2lnd%Faxa_swvdf (begc:endc) ) !; atm2lnd%Faxa_swvdf (:) =  40.0d0

        do i_local = begc, endc

            i_global = gindex_atm(i_local)
            lon = centerCoords(1,i_global) 
            lat = centerCoords(2,i_global)

            ! rounding to nearest int 
            lon = real(nint(lon))
            lat = real(nint(lat))
            ! This is i_local
            print *, "i_local is:", i_local, "i_global is :", i_global, "lon:", lon, "lat:", lat
            !atm2lnd%Sa_tbot(i_local) = 280.0d0 + (sin (lat)+ cos(lon))*1.0d0
            !atm2lnd%Sa_tbot(i_local) = 280.0d0 + cos(lon)*1.0d0

            atm2lnd%Sa_z       (i_local) =  30.0d0     + lat *0.01d0             + lon *0.01d0
            atm2lnd%Sa_topo    (i_local) =  10.0d0     + lat *0.01d0             + lon *0.01d0
            atm2lnd%Sa_u       (i_local) =  20.0d0     + lat *0.01d0             + lon *0.01d0
            atm2lnd%Sa_v       (i_local) =  40.0d0     + lat *0.01d0             + lon *0.01d0
            atm2lnd%Sa_ptem    (i_local) =  280.0d0    + lat *0.01d0             + lon *0.01d0
            atm2lnd%Sa_pbot    (i_local) =  100100.0d0 + lat *0.01d0             + lon *0.01d0
            atm2lnd%Sa_tbot    (i_local) = 280.0d0     + lat *0.01d0             + lon *0.01d0
            atm2lnd%Sa_shum    (i_local) =  0.0004d0   !+(lat*0.01d0 + lon*0.01d0)*1.0e-8
            atm2lnd%Faxa_lwdn  (i_local) =  200.0d0    + lat *0.01d0             + lon *0.01d0

            !atm2lnd%Faxa_rainc (i_local) =  0.0d0    +  (lat*0.01d0 + lon*0.01d0)*1.0e-8
            atm2lnd%Faxa_rainl (i_local) =  3.0d-8   +  (lat*0.01d0 + lon*0.01d0)*1.0e-8
            atm2lnd%Faxa_snowc (i_local) =  1.0d-8   +  (lat*0.01d0 + lon*0.01d0)*1.0e-8
            atm2lnd%Faxa_snowl (i_local) =  2.0d-8   +  (lat*0.01d0 + lon*0.01d0)*1.0e-8
            atm2lnd%Faxa_swndr (i_local) =  100.0d0    + lat *0.01d0             + lon *0.01d0
            atm2lnd%Faxa_swvdr (i_local) =  50.0d0     + lat *0.01d0             + lon *0.01d0
            atm2lnd%Faxa_swndf (i_local) =  20.0d0     + lat *0.01d0             + lon *0.01d0
            atm2lnd%Faxa_swvdf (i_local) =  40.0d0     + lat *0.01d0             + lon *0.01d0
            !atm2lnd%Sa_tbot(i) = 280.0                + sin (         lat )*1.0
            !atm2lnd%Sa_tbot(i) = 280.0                + cos(lon)*1.0

            ! radian instead of degrees:
            !lon = lon* PI/180.0
            !lat = lat* PI/180.0
        end do

        !allocating these values from atmosphere for now!
        !allocate ( atm2lnd%Sa_z       (begc:endc) ) ; atm2lnd%Sa_z       (:) =  30.0d0
        !allocate ( atm2lnd%Sa_topo    (begc:endc) ) ; atm2lnd%Sa_topo    (:) =  10.0d0
        !allocate ( atm2lnd%Sa_u       (begc:endc) ) ; atm2lnd%Sa_u       (:) =  20.0d0
        !allocate ( atm2lnd%Sa_v       (begc:endc) ) ; atm2lnd%Sa_v       (:) =  40.0d0
        !allocate ( atm2lnd%Sa_ptem    (begc:endc) ) ; atm2lnd%Sa_ptem    (:) =  280.0d0
        !allocate ( atm2lnd%Sa_pbot    (begc:endc) ) ; atm2lnd%Sa_pbot    (:) =  100100.0d0
        !allocate ( atm2lnd%Sa_tbot    (begc:endc) ) ; atm2lnd%Sa_tbot    (:) =  280.0d0
        !allocate ( atm2lnd%Sa_shum    (begc:endc) ) ; atm2lnd%Sa_shum    (:) =  0.0004d0
        !allocate ( atm2lnd%Faxa_lwdn  (begc:endc) ) ; atm2lnd%Faxa_lwdn  (:) =  200.0d0
        !allocate ( atm2lnd%Faxa_rainc (begc:endc) ) ; atm2lnd%Faxa_rainc (:) =  4.0d-8
        allocate ( atm2lnd%Faxa_rainc (begc:endc) ) ; atm2lnd%Faxa_rainc (:) =  0.0d0
        !allocate ( atm2lnd%Faxa_rainl (begc:endc) ) ; atm2lnd%Faxa_rainl (:) =  3.0d-8
        !allocate ( atm2lnd%Faxa_snowc (begc:endc) ) ; atm2lnd%Faxa_snowc (:) =  1.0d-8
        !allocate ( atm2lnd%Faxa_snowl (begc:endc) ) ; atm2lnd%Faxa_snowl (:) =  2.0d-8

        !allocate ( atm2lnd%Faxa_swndr (begc:endc) ) ; atm2lnd%Faxa_swndr (:) =  100.0d0
        !allocate ( atm2lnd%Faxa_swvdr (begc:endc) ) ; atm2lnd%Faxa_swvdr (:) =  50.0d0
        !allocate ( atm2lnd%Faxa_swndf (begc:endc) ) ; atm2lnd%Faxa_swndf (:) =  20.0d0
        !allocate ( atm2lnd%Faxa_swvdf (begc:endc) ) ; atm2lnd%Faxa_swvdf (:) =  40.0d0
        !allocate ( atm2lnd%Faxa_bcph  (begc:endc) )  ; atm2lnd%Faxa_bcph (:) =  0.0d0


        allocate ( lnd2atm%Sl_lfrin (begc:endc) ) ; lnd2atm%Sl_lfrin (:) =  0
        allocate ( lnd2atm%Sl_t     (begc:endc) ) ; lnd2atm%Sl_t     (:) =  0
        allocate ( lnd2atm%Sl_tref  (begc:endc) ) ; lnd2atm%Sl_tref  (:) =  0
        allocate ( lnd2atm%Sl_qref  (begc:endc) ) ; lnd2atm%Sl_qref  (:) =  0
        allocate ( lnd2atm%Sl_avsdr (begc:endc) ) ; lnd2atm%Sl_avsdr (:) =  0
        allocate ( lnd2atm%Sl_anidr (begc:endc) ) ; lnd2atm%Sl_anidr (:) =  0
        allocate ( lnd2atm%Sl_avsdf (begc:endc) ) ; lnd2atm%Sl_avsdf (:) =  0
        allocate ( lnd2atm%Sl_anidf (begc:endc) ) ; lnd2atm%Sl_anidf (:) =  0
        allocate ( lnd2atm%Sl_snowh (begc:endc) ) ; lnd2atm%Sl_snowh (:) =  0
        allocate ( lnd2atm%Sl_u10   (begc:endc) ) ; lnd2atm%Sl_u10   (:) =  0
        allocate ( lnd2atm%Sl_fv    (begc:endc) ) ; lnd2atm%Sl_fv    (:) =  0
        allocate ( lnd2atm%Sl_ram1  (begc:endc) ) ; lnd2atm%Sl_ram1  (:) =  0
    end subroutine fill_in

end module demo_mod



program demo_lilac_driver

    !----------------------------------------------------------------------------
    !***  All the components are in the hierarchy seen here:
    !
    !           main driver* (WRF)
    !               |
    !               |
    !          lilac (not a gridded component!)
    !               |     |________________________.
    !               |                              |
    !           atmos cap                      land cap ____________.     ......... gridded components
    !               |                              |                |
    !               |                              |             river cap
    !           ocean (MOM, POM)?                  |                |
    !                                              |            Mizzouroute...
    !                                            CTSM
    !
    !
    !----------------------------------------------------------------------------

    ! modules
    use ESMF
    use lilac_mod
    use lilac_utils , only : atm2lnd_data1d_type , lnd2atm_data1d_type , atm2lnd_data2d_type , atm2lnd_data2d_type , this_clock
    use clm_varctl  , only : iulog
    use spmdMod     , only : masterproc
    use demo_mod    , only : demo_init
    use demo_mod    , only : read_netcdf_mesh
    implicit none

    ! TO DO: change the name and the derived data types
    ! data types for 1d arrays for meshes
    type (atm2lnd_data1d_type)                            :: atm2lnd
    type (lnd2atm_data1d_type)                            :: lnd2atm

    type (this_clock)                                     :: this_time

    real    , allocatable                                 :: rand1(:)
    real    , allocatable                                 :: rand2(:)

    integer , allocatable                                 :: seed(:)
    integer                                               :: seed_val, n

    integer                                               :: begc,endc
    integer                                               :: start_time               !-- start_time    start time
    integer                                               :: end_time                 !-- end_time      end time
    integer                                               :: curr_time                !-- cur_time      current time
    integer                                               :: itime_step               !-- itime_step    counter of time steps
    integer                                               :: g,i,k                    !-- indices
    integer, parameter                                    :: debug = 1                !-- internal debug level

    character(len=128)                                    :: fldname

    character(*),parameter                                :: F01 =   "(a,i4,d26.19)"
    character(*),parameter                                :: F02 = "('[demo_driver]',a,i5,2x,d26.19)"
    integer , allocatable :: gindex_atm(:)

    !------------------------------------------------------------------------
    ! real atmosphere:
    begc       = 1
    !endc       = 6912/4/2
    endc        = 3312/4/2/2
    !endc       = 13824
    !endc       = 13968

    start_time = 1
    end_time   = 48
    itime_step = 1

    seed_val   = 0
    n          = endc - begc + 1


    ! making 2 random arrays with a seed.
    call random_seed   (size = n    )
    allocate ( seed                    (n        ) ) ; seed            (:)      =  seed_val
    call random_seed   (put  = seed )

    allocate ( rand1                   (begc:endc) ) ; call random_number (rand1)
    allocate ( rand2                   (begc:endc) ) ; call random_number (rand2)


    !fldname = 'Sa_topo'
    !if (debug > 0) then
    !    do i=begc, endc 
    !        write (iulog,F02)'import: nstep, n, '//trim(fldname)//' = ',i, atm2lnd%Sa_topo(i)
    !    enddo
    ! end if


    !print *, atm2lnd%Sa_topo(1:100)



    !------------------------------------------------------------------------
    ! The newly added demo_init
    ! all allocate will go here:
    !------------------------------------------------------------------------

    call demo_init(gindex_atm, atm2lnd , lnd2atm)

    !------------------------------------------------------------------------
    ! looping over imaginary time ....
    !------------------------------------------------------------------------

    call lilac_init     ( atm2lnd1d = atm2lnd   ,   lnd2atm1d =  lnd2atm , gindex_atm = gindex_atm )
    do curr_time = start_time, end_time
            call lilac_run      ( )
           itime_step = itime_step + 1
    end do
    call lilac_final    ( )
    call ESMF_Finalize  ( )

    print *,  "======================================="
    print *,  " ............. DONE ..................."
    print *,  "======================================="


end program demo_lilac_driver

