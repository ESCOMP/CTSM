module demo_mod
!----------------------------------------------------------------------------
    use mpi           ,  only : MPI_COMM_WORLD, MPI_COMM_NULL, MPI_Init, MPI_FINALIZE, MPI_SUCCESS
    use spmdMod       ,  only : masterproc
    implicit none
    private
    public :: demo_init
    integer                            :: ierr
    integer                                         :: COMP_COMM
    integer                            :: npts   ! domain global size
    integer                            :: num_local
!----------------------------------------------------------------------------
contains
!----------------------------------------------------------------------------
    subroutine demo_init(gindex_atm)
        !! TODO: IS THE INTENT CORRECT FOR GINDEX_ATM
        integer , allocatable, intent(inout)  :: gindex_atm(:)
        integer                :: ntasks
        integer                :: mytask
        !-----------------------------------------------------------------------------
        ! Initiallize MPI
        !-----------------------------------------------------------------------------

        npts = 3312
        ! this is coming from
        ! /glade/work/mvertens/ctsm.nuopc/cime/src/drivers/nuopc/drivers/cime/esmApp.F90
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

        if (masterproc) then
            print *, "MPI initialization done ..., ntasks=", ntasks
        end if
        
        call decompInit_atm( ntasks, mytask, gindex_atm)
        print *, "gindex_atm for ", mytask,"is: ", gindex_atm
        print *, "size gindex_atm for ", mytask,"is: ", size(gindex_atm)
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

    !allocating these values from atmosphere for now!
    allocate ( atm2lnd%Sa_z       (begc:endc) ) ; atm2lnd%Sa_z       (:) =  30.0d0
    allocate ( atm2lnd%Sa_topo    (begc:endc) ) ; atm2lnd%Sa_topo    (:) =  10.0d0
    allocate ( atm2lnd%Sa_u       (begc:endc) ) ; atm2lnd%Sa_u       (:) =  20.0d0
    allocate ( atm2lnd%Sa_v       (begc:endc) ) ; atm2lnd%Sa_v       (:) =  40.0d0
    allocate ( atm2lnd%Sa_ptem    (begc:endc) ) ; atm2lnd%Sa_ptem    (:) =  280.0d0
    allocate ( atm2lnd%Sa_pbot    (begc:endc) ) ; atm2lnd%Sa_pbot    (:) =  100100.0d0
    allocate ( atm2lnd%Sa_tbot    (begc:endc) ) ; atm2lnd%Sa_tbot    (:) =  280.0d0
    allocate ( atm2lnd%Sa_shum    (begc:endc) ) ; atm2lnd%Sa_shum    (:) =  0.0004d0
    allocate ( atm2lnd%Faxa_lwdn  (begc:endc) ) ; atm2lnd%Faxa_lwdn  (:) =  200.0d0
    !allocate ( atm2lnd%Faxa_rainc (begc:endc) ) ; atm2lnd%Faxa_rainc (:) =  4.0d-8
    allocate ( atm2lnd%Faxa_rainc (begc:endc) ) ; atm2lnd%Faxa_rainc (:) =  0.0d0
    allocate ( atm2lnd%Faxa_rainl (begc:endc) ) ; atm2lnd%Faxa_rainl (:) =  3.0d-8
    allocate ( atm2lnd%Faxa_snowc (begc:endc) ) ; atm2lnd%Faxa_snowc (:) =  1.0d-8
    allocate ( atm2lnd%Faxa_snowl (begc:endc) ) ; atm2lnd%Faxa_snowl (:) =  2.0d-8
    allocate ( atm2lnd%Faxa_swndr (begc:endc) ) ; atm2lnd%Faxa_swndr (:) =  100.0d0

    allocate ( atm2lnd%Faxa_swvdr (begc:endc) ) ; atm2lnd%Faxa_swvdr (:) =  50.0d0
    allocate ( atm2lnd%Faxa_swndf (begc:endc) ) ; atm2lnd%Faxa_swndf (:) =  20.0d0
    allocate ( atm2lnd%Faxa_swvdf (begc:endc) ) ; atm2lnd%Faxa_swvdf (:) =  40.0d0


    fldname = 'Sa_topo'
    if (debug > 0) then
        do i=begc, endc 
            write (iulog,F02)'import: nstep, n, '//trim(fldname)//' = ',i, atm2lnd%Sa_topo(i)
        enddo
    end if
    !allocate ( atm2lnd%Faxa_bcph  (begc:endc) )  ; atm2lnd%Faxa_bcph (:) =  0.0d0

    !endc       = 18048 ? should this be the size of the land or atmosphere???


    !print *, atm2lnd%Sa_topo(1:100)


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



    !------------------------------------------------------------------------
    ! The newly added demo_init
    !------------------------------------------------------------------------

    call demo_init(gindex_atm)

    !------------------------------------------------------------------------
    ! looping over imaginary time ....
    !------------------------------------------------------------------------

    call lilac_init     ( atm2lnd1d = atm2lnd   ,   lnd2atm1d =  lnd2atm )
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

