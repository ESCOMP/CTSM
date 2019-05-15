module lilac_mod


    !-----------------------------------------------------------------------
    ! !DESCRIPTION:

    ! !USES
    use ESMF
    use lilac_utils
    use atmos_cap ,  only :         atmos_register
    use lnd_cap   ,  only :         lnd_register
    use cpl_mod   ,  only :         cpl_atm2lnd_register , cpl_lnd2atm_register

    implicit none


    public :: lilac_init
    public :: lilac_run

    character(*) , parameter                          :: modname     = "lilac_mod"
    !type(fld_list_type), public                      :: a2c_fldlist, c2a_fldlist  !defined in atmosphere and land caps....

    !------------------------------------------------------------------------
    ! !Clock, TimeInterval, and Times
    type(ESMF_Clock)                                  :: clock
    type(ESMF_TimeInterval)                           :: timeStep
    type(ESMF_Time)                                   :: startTime
    type(ESMF_Time)                                   :: stopTime
    type(ESMF_Alarm)                                  :: EAlarm_stop, EAlarm_rest
    type(ESMF_Calendar),target                        :: calendar
    integer                                           :: yy,mm,dd,sec
    !  ! Gridded Components and Coupling Components
    type(ESMF_GridComp)                               :: dummy_atmos_comp
    type(ESMF_GridComp)                               :: dummy_land_comp
    type(ESMF_CplComp)                                :: cpl_atm2lnd_comp
    type(ESMF_CplComp)                                :: cpl_lnd2atm_comp
    type(ESMF_State)                                  :: atm2lnd_l_state , atm2lnd_a_state
    type(ESMF_State)                                  :: lnd2atm_a_state, lnd2atm_l_state

    !========================================================================
    contains
    !========================================================================

    subroutine lilac_init( atm2lnd1d, atm2lnd2d, lnd2atm1d, lnd2atm2d)

        use atmos_cap ,  only :          a2c_fldlist, c2a_fldlist
        use lnd_cap   ,  only :          l2c_fldlist, c2l_fldlist

        character(len=*), parameter                      :: subname=trim(modname)//': [lilac_init] '

        ! input/output variables
        type(atm2lnd_data1d_type), intent(in), optional  :: atm2lnd1d
        type(atm2lnd_data2d_type), intent(in), optional  :: atm2lnd2d
        type(lnd2atm_data1d_type), intent(in), optional  :: lnd2atm1d
        type(lnd2atm_data2d_type), intent(in), optional  :: lnd2atm2d

        ! local variables

        type(ESMF_State)                                 :: importState, exportState

        !character(len=*)                                :: atm_mesh_filepath   !!! For now this is hardcoded in the atmos init

        integer                                          :: rc         , urc
        character(len=ESMF_MAXSTR)                       :: gcname1    , gcname2   !    Gridded components names
        character(len=ESMF_MAXSTR)                       :: ccname1    , ccname2   !    Coupling components names
        integer                                          :: a2l_fldnum , l2a_fldnum

        !------------------------------------------------------------------------
        ! Initialize return code
        rc = ESMF_SUCCESS

        a2l_fldnum  = 3
        l2a_fldnum  = 3

        !-------------------------------------------------------------------------
        ! Initialize ESMF, set the default calendar and log type.
        !-------------------------------------------------------------------------
        call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//".........................", ESMF_LOGMSG_INFO)
        call ESMF_LogWrite(subname//"Initializing ESMF        ", ESMF_LOGMSG_INFO)

        print *,  "---------------------------------------"
        print *,  "    Lilac Demo Application Start       "
        print *,  "---------------------------------------"

        !-------------------------------------------------------------------------
        ! Create Field lists -- Basically create a list of fields and add a default
        ! value to them.
        !-------------------------------------------------------------------------

        allocate (a2c_fldlist(a2l_fldnum))
        allocate (c2a_fldlist(l2a_fldnum))

        allocate (l2c_fldlist(l2a_fldnum))
        allocate (c2l_fldlist(a2l_fldnum))
        print *, "creating empty field lists !"

        ! ======================================================================= ! maybe move to create_fldlist?
        !    call create_fldlists(c2a_fldlist, a2c_fldlist, )

        !-------------------------------------------------------------------------
        !            !---- from atm ----! a2c_fldlist & c2l_fldlist
        !-------------------------------------------------------------------------

        a2c_fldlist(1)%stdname      =  'uwind'
        a2c_fldlist(1)%farrayptr1d  => atm2lnd1d%uwind !*** this now sets the module variable memory in atmos_cap.F90
        print *,      a2c_fldlist(1)%stdname
        !print *,      a2c_fldlist(1)%farrayptr1d(:)
        a2c_fldlist(2)%stdname      =  'vwind'
        a2c_fldlist(2)%farrayptr1d  => atm2lnd1d%vwind !*** this now sets the module variable memory in atmos_cap.F90
        print *,      a2c_fldlist(2)%stdname
        !print *,      a2c_fldlist(2)%farrayptr1d(:)
        a2c_fldlist(3)%stdname      =  'tbot'
        a2c_fldlist(3)%farrayptr1d  => atm2lnd1d%vwind
        print *,      a2c_fldlist(3)%stdname
        !print *,      a2c_fldlist(3)%farrayptr1d


        !!! Where should these point to? pointer to an empty array which will be filled in the land....

        ! Similary we need c2a_fldlist
        c2l_fldlist(1)%stdname      =  'uwind'
        print *,      c2l_fldlist(1)%stdname
        c2l_fldlist(2)%stdname      =  'vwind'
        print *,      c2l_fldlist(1)%stdname
        c2l_fldlist(3)%stdname      =  'tbot'
        print *,      c2l_fldlist(1)%stdname

        !-------------------------------------------------------------------------
        !            !---- from land ----! l2c_fldlist & c2a_fldlist
        !-------------------------------------------------------------------------


        l2c_fldlist(1)%stdname      =  'lwup'
        l2c_fldlist(1)%farrayptr1d  => lnd2atm1d%lwup
        print *,      l2c_fldlist(1)%stdname

        l2c_fldlist(2)%stdname      =  'taux'
        print *,      l2c_fldlist(2)%stdname
        l2c_fldlist(2)%farrayptr1d  => lnd2atm1d%taux

        l2c_fldlist(3)%stdname      =  'tauy'
        print *,      l2c_fldlist(3)%stdname
        l2c_fldlist(3)%farrayptr1d  => lnd2atm1d%taux


        c2a_fldlist(1)%stdname      =  'lwup'
        print *,      c2a_fldlist(1)%stdname

        c2a_fldlist(2)%stdname      =  'taux'
        print *,      c2a_fldlist(2)%stdname

        c2a_fldlist(3)%stdname      =  'tauy'
        print *,      c2a_fldlist(3)%stdname

        ! ======================================================================= ! create_fldlist

        !-------------------------------------------------------------------------
        ! Create Gridded Component!  --- dummy atmosphere ( atmos_cap)
        !-------------------------------------------------------------------------
        gcname1 = "Dummy Atmosphere or Atmosphere Cap"
        dummy_atmos_comp = ESMF_GridCompCreate(name=gcname1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"Created "//trim(gcname1)//" component", ESMF_LOGMSG_INFO)
        print *, "Dummy Atmosphere Gridded Component Created!"

        !-------------------------------------------------------------------------
        ! Create Gridded Component!  --- dummy land       ( land cap )
        !-------------------------------------------------------------------------
        gcname2 = "Dummy Land or Land Cap"
        dummy_land_comp = ESMF_GridCompCreate(name=gcname2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"Created "//trim(gcname2)//" component", ESMF_LOGMSG_INFO)
        print *, "Dummy Land  Gridded Component Created!"

        !-------------------------------------------------------------------------
        ! Create Coupling Component! --- Coupler  from atmos  to land
        !-------------------------------------------------------------------------
        ccname1 = "Coupler from atmosphere to land"
        cpl_atm2lnd_comp = ESMF_CplCompCreate(name=ccname1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"Created "//trim(ccname1)//" component", ESMF_LOGMSG_INFO)
        print *, "1st Coupler Component (atmosphere to land ) Created!"

        !-------------------------------------------------------------------------
        ! Create Coupling  Component!  -- Coupler from land to atmos
        !-------------------------------------------------------------------------
        ccname2 = "Coupler from land to atmosphere"
        cpl_lnd2atm_comp = ESMF_CplCompCreate(name=ccname2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"Created "//trim(ccname2)//" component", ESMF_LOGMSG_INFO)
        print *, "2nd Coupler Component (land to atmosphere) Created!"

        ! ========================================================================

        !-------------------------------------------------------------------------
        ! Register section -- set services -- atmos_cap
        !-------------------------------------------------------------------------
        call ESMF_GridCompSetServices(dummy_atmos_comp, userRoutine=atmos_register, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"dummy atmos SetServices finished!", ESMF_LOGMSG_INFO)
        print *, "Dummy Atmosphere Gridded Component SetServices finished!"
        !-------------------------------------------------------------------------
        ! Register section -- set services -- land cap
        !-------------------------------------------------------------------------
        call ESMF_GridCompSetServices(dummy_land_comp, userRoutine=lnd_register, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"land SetServices finished!", ESMF_LOGMSG_INFO)
        print *, "Land  Gridded Component SetServices finished!"
        !-------------------------------------------------------------------------
        ! Register section -- set services -- coupler atmosphere to land
        !-------------------------------------------------------------------------
        call ESMF_CplCompSetServices(cpl_atm2lnd_comp, userRoutine=cpl_atm2lnd_register, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"Coupler from atmosphere to land  SetServices finished!", ESMF_LOGMSG_INFO)
        print *, "Coupler from atmosphere to land SetServices finished!"
        !-------------------------------------------------------------------------
        ! Register section -- set services -- coupler land to atmosphere
        !-------------------------------------------------------------------------
        call ESMF_CplCompSetServices(cpl_lnd2atm_comp, userRoutine=cpl_lnd2atm_register, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"Coupler from land to atmosphere SetServices finished!", ESMF_LOGMSG_INFO)
        print *, "Coupler from land to atmosphere SetServices finished!"

        ! ========================================================================

        !-------------------------------------------------------------------------
        !  Create and initialize a clock!
        ! ????? Should I create a clock here or in driver?
        !-------------------------------------------------------------------------
        calendar = ESMF_CalendarCreate(name='lilac_drv_NOLEAP', calkindflag=ESMF_CALKIND_NOLEAP, rc=rc )
        call ESMF_TimeSet(StartTime, yy=2000, mm=1, dd=1, s=0, calendar=Calendar, rc=rc)
        call ESMF_TimeSet(StopTime , yy=2000, mm=1, dd=10, s=0, calendar=Calendar, rc=rc)
        call ESMF_TimeIntervalSet(TimeStep, s=3600, rc=rc)
        clock = ESMF_ClockCreate(name='lilac_drv_EClock', TimeStep=TimeStep, startTime=StartTime, RefTime=StartTime, stopTime=stopTime, rc=rc)
        !clock = ESMF_ClockCreate(timeStep=timeStep, startTime=startTime, stopTime=stopTime, rc=rc)
        !EClock = ESMF_ClockCreate(name='lilac_drv_EClock', TimeStep=TimeStep, startTime=StartTime, RefTime=StartTime, stopTime=stopTime, rc=rc)

        !-------------------------------------------------------------------------
        ! Create the necessary import and export states used to pass data
        !  between components.
        !-------------------------------------------------------------------------

        ! following 4 states are lilac module variables:
        ! 1- atm2lnd_a_state  2- atm2lnd_l_state 3- lnd2atm_a_state 4-lnd2atm_l_state

        atm2lnd_a_state = ESMF_StateCreate(name=gcname1, stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        atm2lnd_l_state = ESMF_StateCreate(name=gcname1, stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        lnd2atm_a_state = ESMF_StateCreate(name=gcname2, stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        lnd2atm_l_state = ESMF_StateCreate(name=gcname2, stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        call ESMF_LogWrite(subname//"Empty import and export states are created!!", ESMF_LOGMSG_INFO)
        print *, "Empty import and export states are created!!"

        ! returns a valid state_to_lnd_atm and an empty state_from_land_atmgrid

        ! -------------------------------------------------------------------------
        ! Grid Componenet Initialization -- 1- atmos cap   2- lnd cap             !
        !                                   3- cpl_atm2lnd 4- cpl_lnd2atm         !
        ! -------------------------------------------------------------------------

        call ESMF_GridCompInitialize(dummy_atmos_comp, importState=lnd2atm_a_state, exportState=atm2lnd_a_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"atmos_cap or dummy_atmos_comp initialized", ESMF_LOGMSG_INFO)
        print *, "atmos_cap initialize finished, rc =", rc

        call ESMF_GridCompInitialize(dummy_land_comp       , importState=atm2lnd_l_state, exportState=lnd2atm_l_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"lnd_cap or dummy_land_comp initialized", ESMF_LOGMSG_INFO)
        print *, "lnd_cap initialize finished, rc =", rc

        ! All 4 states that are module variables are no longer empty - have been initialized

        call ESMF_CplCompInitialize(cpl_atm2lnd_comp, importState=atm2lnd_a_state, exportState=atm2lnd_l_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"coupler :: cpl_atm2lnd_comp initialized", ESMF_LOGMSG_INFO)
        print *, "coupler :: cpl_atm2lnd_comp initialize finished, rc =", rc

        call ESMF_CplCompInitialize(cpl_lnd2atm_comp, importState=lnd2atm_l_state, exportState=lnd2atm_a_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"coupler :: cpl_lnd2atm_comp initialized", ESMF_LOGMSG_INFO)
        print *, "coupler :: cpl_lnd2atm_comp initialize finished, rc =", rc

    end subroutine lilac_init

    !========================================================================

    subroutine lilac_run( )

        use atmos_cap, only : a2c_fldlist, c2a_fldlist
        use lnd_cap,   only : l2c_fldlist, c2l_fldlist

        character(len=*), parameter                      :: subname=trim(modname)//': [lilac_run] '
        type(ESMF_State)                                 :: importState, exportState

        ! local variables
        integer                                          :: rc, urc
        character(len=ESMF_MAXSTR)                       :: gcname1, gcname2   !    Gridded components names
        character(len=ESMF_MAXSTR)                       :: ccname1, ccname2   !    Coupling components names
        !integer, parameter                              :: fldsMax = 100
        integer                                          :: a2l_fldnum, l2a_fldnum

        ! input/output variables
        !type(atm2lnd_data1d_type), intent(in), optional  :: atm2lnd1d
        !type(atm2lnd_data2d_type), intent(in), optional  :: atm2lnd2d
        !type(lnd2atm_data1d_type), intent(in), optional  :: lnd2atm1d
        !type(lnd2atm_data2d_type), intent(in), optional  :: lnd2atm2d

        !------------------------------------------------------------------------
        ! Initialize return code
        rc = ESMF_SUCCESS


        !-------------------------------------------------------------------------
        ! Initialize ESMF, set the default calendar and log type.
        !-------------------------------------------------------------------------
        call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        print *,  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print *,  "               Lilac Run               "
        print *,  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

        !-------------------------------------------------------------------------
        ! We are running components in this order:
        ! 1- atmos_cap   2- cpl_atm2lnd
        ! 3-   lnd_cap   4- cpl_lnd2atm
        !-------------------------------------------------------------------------
        call ESMF_GridCompRun(dummy_atmos_comp, importState=lnd2atm_a_state, exportState=atm2lnd_a_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"atmos_cap or dummy_atmos_comp is running", ESMF_LOGMSG_INFO)
        print *, "Running atmos_cap gridded component , rc =", rc

        call ESMF_CplCompRun(cpl_atm2lnd_comp, importState=atm2lnd_a_state, exportState=atm2lnd_l_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"running cpl_atm2lnd_comp ", ESMF_LOGMSG_INFO)
        print *, "Running coupler component..... cpl_atm2lnd_comp , rc =", rc

        call ESMF_GridCompRun(dummy_land_comp,  importState=atm2lnd_l_state, exportState=lnd2atm_l_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"lnd_cap or dummy_land_comp is running", ESMF_LOGMSG_INFO)
        print *, "Running lnd_cap gridded component , rc =", rc

        call ESMF_CplCompRun(cpl_lnd2atm_comp, importState=lnd2atm_l_state, exportState=lnd2atm_a_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"running cpl_lnd2atm_comp ", ESMF_LOGMSG_INFO)
        print *, "Running coupler component..... cpl_lnd2atm_comp , rc =", rc

    end subroutine lilac_run


    subroutine lilac_final( )

        use atmos_cap, only : a2c_fldlist, c2a_fldlist
        use lnd_cap,   only : l2c_fldlist, c2l_fldlist


        character(len=*), parameter                      :: subname=trim(modname)//': [lilac_final] '
        type(ESMF_State)                                 :: importState, exportState

        ! local variables
        integer                                          :: rc, urc
        character(len=ESMF_MAXSTR)                       :: gcname1, gcname2   !    Gridded components names
        character(len=ESMF_MAXSTR)                       :: ccname1, ccname2   !    Coupling components names
        !integer, parameter                              :: fldsMax = 100
        integer                                          :: a2l_fldnum, l2a_fldnum

        !------------------------------------------------------------------------
        !------------------------------------------------------------------------
        ! Initialize return code
        rc = ESMF_SUCCESS

        print *,  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print *,  "        Lilac Finalizing               "
        print *,  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        !-------------------------------------------------------------------------
        ! Gridded Component Finalizing!     ---     dummy atmosphere 
        !-------------------------------------------------------------------------
        call ESMF_GridCompFinalize(dummy_atmos_comp, importState=lnd2atm_a_state, exportState=atm2lnd_a_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"atmos_cap or dummy_atmos_comp is running", ESMF_LOGMSG_INFO)
        print *, "Finalizing atmos_cap gridded component , rc =", rc

        !-------------------------------------------------------------------------
        ! Coupler component Finalizing     ---    coupler atmos to land
        !-------------------------------------------------------------------------
        call ESMF_CplCompFinalize(cpl_atm2lnd_comp, importState=atm2lnd_a_state, exportState=atm2lnd_l_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"running cpl_atm2lnd_comp ", ESMF_LOGMSG_INFO)
        print *, "Finalizing coupler component..... cpl_atm2lnd_comp , rc =", rc

        !-------------------------------------------------------------------------
        ! Gridded Component Finalizing!     ---     dummy land 
        !-------------------------------------------------------------------------
        call ESMF_GridCompFinalize(dummy_land_comp,  importState=atm2lnd_l_state, exportState=lnd2atm_l_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"lnd_cap or dummy_land_comp is running", ESMF_LOGMSG_INFO)
        print *, "Finalizing lnd_cap gridded component , rc =", rc

        !-------------------------------------------------------------------------
        ! Coupler component Finalizing     ---    coupler land to atmos
        !-------------------------------------------------------------------------
        call ESMF_CplCompFinalize(cpl_lnd2atm_comp, importState=lnd2atm_l_state, exportState=lnd2atm_a_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"running cpl_lnd2atm_comp ", ESMF_LOGMSG_INFO)
        print *, "Finalizing coupler component..... cpl_lnd2atm_comp , rc =", rc


        ! Then clean them up
        call ESMF_LogWrite(subname//".........................", ESMF_LOGMSG_INFO)
        call ESMF_LogWrite(subname//"destroying all states    ", ESMF_LOGMSG_INFO)

        print *, "ready to destroy all states"
        call ESMF_StateDestroy(atm2lnd_a_state , rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
        call ESMF_StateDestroy(atm2lnd_l_state, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
        call ESMF_StateDestroy(lnd2atm_a_state, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
        call ESMF_StateDestroy(lnd2atm_l_state, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

        call ESMF_LogWrite(subname//"destroying all components    ", ESMF_LOGMSG_INFO)
        print *, "ready to destroy all components"

        call ESMF_GridCompDestroy(dummy_atmos_comp, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
        call ESMF_GridCompDestroy(dummy_land_comp, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
        call ESMF_CplCompDestroy(cpl_atm2lnd_comp, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
        call ESMF_CplCompDestroy(cpl_lnd2atm_comp, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

        call ESMF_LogWrite(subname//".........................", ESMF_LOGMSG_INFO)
        print *, "end of Lilac Finalization routine"

        end subroutine lilac_final



    end module lilac_mod

