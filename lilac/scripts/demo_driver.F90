program demo_lilac_driver

    ! modules
    use ESMF
    use lilac_mod

    use lilac_utils, only   : atm2lnd_data1d_type , lnd2atm_data1d_type, atm2lnd_data2d_type, atm2lnd_data2d_type

    type (atm2lnd_data1d_type)                            :: atm2lnd
    type (lnd2atm_data1d_type)                            :: lnd2atm

    integer                                               :: begc,endc

    real,    dimension(100,100), target                   :: dum_var1
    real,    dimension(4608)                              :: dum_var2

    integer, dimension(:),       allocatable              :: seed
    integer                                               :: seed_val, n

    integer                                               :: start_time               !-- start_time    start time
    integer                                               :: end_time                 !-- end_time      end time
    integer                                               :: curr_time                !-- cur_time      current time
    integer                                               :: itime_step               !-- itime_step    counter of time steps

    !------------------------------------------------------------------------


    begc       = 1
    endc       = 4608

    start_time = 1
    end_time   = 10
    itime_step = 1

    seed_val   = 0
    n          = endc - begc + 1




    call random_seed   (size = n)
    allocate           (seed(n))                     ; seed            (:)      =  seed_val
    call random_seed   (put = seed)
    call random_number (dum_var2)


    allocate( atm2lnd%uwind           (begc:endc) )  ; atm2lnd%uwind   (:)      =  dum_var2
    allocate( atm2lnd%vwind           (begc:endc) )  ; atm2lnd%vwind   (:)      =  dum_var2
    allocate( atm2lnd%tbot            (begc:endc) )  ; atm2lnd%tbot    (:)      =  dum_var2

    allocate( lnd2atm%lwup            (begc:endc) )  ; lnd2atm%lwup    (:)      =  dum_var2
    allocate( lnd2atm%taux            (begc:endc) )  ; lnd2atm%taux    (:)      =  dum_var2
    allocate( lnd2atm%tauy            (begc:endc) )  ; lnd2atm%tauy    (:)      =  dum_var2

    print *,  "======================================="
    print *,  atm2lnd%uwind(1:10)
    print *,  "======================================="


    ! dummy looping over imaginary time ....

    do curr_time = start_time, end_time

        if (curr_time == start_time) then
            ! Initalization phase
            call lilac_init     ( atm2lnd1d = atm2lnd   ,   lnd2atm1d =  lnd2atm )

        else if (curr_time == end_time) then
            !Finalization phase
            call lilac_final    ( )
            call ESMF_Finalize  ( )

        else
            call lilac_run      ( )
        endif

        itime_step = itime_step + 1

    end do

    print *,  "======================================="
    print *,  " ............. DONE ..................."



end program demo_lilac_driver

