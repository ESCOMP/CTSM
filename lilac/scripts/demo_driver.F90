program demo_lilac_driver

    ! modules
    use ESMF
    use Lilac_mod

    use lilac_utils, only   : atm2lnd_data1d_type , lnd2atm_data1d_type, atm2lnd_data2d_type, atm2lnd_data2d_type
    !use lilac_utils, only   : atm2lnd_data2d_type
    implicit none

    type (atm2lnd_data1d_type)                            :: atm2lnd
    type (lnd2atm_data1d_type)                            :: lnd2atm
    integer                                               :: begc,endc
    real, dimension(100,100), target                      :: dum_var1
    real, dimension(100)                                  :: dum_var2

    !------------------------------------------------------------------------


    begc = 1
    endc = 100

    call random_number(dum_var2)

    allocate( atm2lnd%uwind           (begc:endc) )  ; atm2lnd%uwind   (:)      =  dum_var2 
    allocate( atm2lnd%vwind           (begc:endc) )  ; atm2lnd%vwind   (:)      =  dum_var2
    allocate( atm2lnd%tbot            (begc:endc) )  ; atm2lnd%tbot    (:)      =  dum_var2
    allocate( lnd2atm%lwup            (begc:endc) )  ; lnd2atm%lwup    (:)      =  dum_var2
    allocate( lnd2atm%taux            (begc:endc) )  ; lnd2atm%taux    (:)      =  dum_var2
    allocate( lnd2atm%tauy            (begc:endc) )  ; lnd2atm%tauy    (:)      =  dum_var2

    print *,  "======================================="
    print *,  atm2lnd%uwind



    call lilac_init ( atm2lnd1d = atm2lnd   ,   lnd2atm1d =  lnd2atm )

    print *,  "======================================="
    print *,  " ............. DONE ..................."


    call ESMF_Finalize()

end program demo_lilac_driver

