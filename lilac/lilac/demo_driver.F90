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
    !           oceaan (MOM, POM)?                 |                |
    !                                              |            Mizzouroute...
    !                                            CTSM
    !
    !
    !----------------------------------------------------------------------------

    ! modules
    use ESMF
    use lilac_mod
    use lilac_utils, only   : atm2lnd_data1d_type , lnd2atm_data1d_type, atm2lnd_data2d_type, atm2lnd_data2d_type , this_clock

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

    !------------------------------------------------------------------------

    ! real atmosphere:
    begc       = 1
    !endc       = 10
    endc       = 6912/4
    !endc       = 13824
    !endc       = 13968

    start_time = 1
    end_time   = 50
    itime_step = 1

    seed_val   = 0
    n          = endc - begc + 1


    ! making 2 random arrays with a seed.
    call random_seed   (size = n    )
    allocate ( seed                    (n        ) ) ; seed            (:)      =  seed_val
    call random_seed   (put  = seed )

    allocate ( rand1                   (begc:endc) ) ; call random_number (rand1)
    allocate ( rand2                   (begc:endc) ) ; call random_number (rand2)

    !allocating these values of default for now!
    allocate ( atm2lnd%Sa_z       (begc:endc) ) ; atm2lnd%Sa_z       (:) =  30.0
    allocate ( atm2lnd%Sa_topo    (begc:endc) ) ; atm2lnd%Sa_topo    (:) =  10.0
    allocate ( atm2lnd%Sa_u       (begc:endc) ) ; atm2lnd%Sa_u       (:) =  20.0
    allocate ( atm2lnd%Sa_v       (begc:endc) ) ; atm2lnd%Sa_v       (:) =  40.0
    allocate ( atm2lnd%Sa_ptem    (begc:endc) ) ; atm2lnd%Sa_ptem    (:) =  280.0
    allocate ( atm2lnd%Sa_pbot    (begc:endc) ) ; atm2lnd%Sa_pbot    (:) =  100100.0
    allocate ( atm2lnd%Sa_tbot    (begc:endc) ) ; atm2lnd%Sa_tbot    (:) =  280.0
    allocate ( atm2lnd%Sa_shum    (begc:endc) ) ; atm2lnd%Sa_shum    (:) =  0.0004
    allocate ( atm2lnd%Faxa_lwdn  (begc:endc) ) ; atm2lnd%Faxa_lwdn  (:) =  200.0
    allocate ( atm2lnd%Faxa_rainc (begc:endc) ) ; atm2lnd%Faxa_rainc (:) =  4.0e-8
    allocate ( atm2lnd%Faxa_rainl (begc:endc) ) ; atm2lnd%Faxa_rainl (:) =  3.0e-8
    allocate ( atm2lnd%Faxa_snowc (begc:endc) ) ; atm2lnd%Faxa_snowc (:) =  1.0e-8
    allocate ( atm2lnd%Faxa_snowl (begc:endc) ) ; atm2lnd%Faxa_snowl (:) =  2.0e-8
    allocate ( atm2lnd%Faxa_swndr (begc:endc) ) ; atm2lnd%Faxa_swndr (:) =  100.0
    allocate ( atm2lnd%Faxa_swvdr (begc:endc) ) ; atm2lnd%Faxa_swvdr (:) =  90.0
    allocate ( atm2lnd%Faxa_swndf (begc:endc) ) ; atm2lnd%Faxa_swndf (:) =  20.0
    allocate ( atm2lnd%Faxa_swvdf (begc:endc) ) ; atm2lnd%Faxa_swvdf (:) =  40.0

    !endc       = 18048 ? should this be the size of the land or atmosphere???



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
    ! looping over imaginary time ....
    !------------------------------------------------------------------------

    do curr_time = start_time, end_time
        if  (curr_time == start_time) then

            ! Initalization phase
            print *,  "--------------------------"
            print *,  " LILAC Initalization phase"
            print *,  "--------------------------"

            call lilac_init     ( atm2lnd1d = atm2lnd   ,   lnd2atm1d =  lnd2atm )
        else if (curr_time == end_time  ) then
            ! Finalization phase
            call lilac_final    ( )
            call ESMF_Finalize  ( )
        else
            call lilac_run      ( )
        endif
        itime_step = itime_step + 1
    end do

    print *,  "======================================="
    print *,  " ............. DONE ..................."
    print *,  "======================================="


end program demo_lilac_driver

