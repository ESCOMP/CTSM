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
    endc       = 13824

    start_time = 1
    end_time   = 10
    itime_step = 1

    seed_val   = 0
    n          = endc - begc + 1


    ! making 2 random arrays with a seed.
    call random_seed   (size = n    )
    allocate ( seed                    (n        ) ) ; seed            (:)      =  seed_val
    call random_seed   (put  = seed )

    allocate ( rand1                   (begc:endc) ) ; call random_number (rand1)
    allocate ( rand2                   (begc:endc) ) ; call random_number (rand2)


    allocate ( atm2lnd%uwind           (begc:endc) ) ; atm2lnd%uwind   (:)      =  rand1
    allocate ( atm2lnd%vwind           (begc:endc) ) ; atm2lnd%vwind   (:)      =  rand1
    allocate ( atm2lnd%tbot            (begc:endc) ) ; atm2lnd%tbot    (:)      =  rand1
    !endc       = 18048 ? should this be the size of the land or atmosphere???
    allocate ( lnd2atm%lwup            (begc:endc) ) ; lnd2atm%lwup    (:)      =  rand2
    allocate ( lnd2atm%taux            (begc:endc) ) ; lnd2atm%taux    (:)      =  rand2
    allocate ( lnd2atm%tauy            (begc:endc) ) ; lnd2atm%tauy    (:)      =  rand2


    print *,  "======================================="
    print *,  atm2lnd%uwind(1:10)
    print *,  "======================================="

    !------------------------------------------------------------------------
    ! looping over imaginary time ....
    !------------------------------------------------------------------------

    do curr_time = start_time, end_time
        if      (curr_time == start_time) then

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

