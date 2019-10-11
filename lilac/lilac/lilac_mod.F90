module lilac_mod


    !-----------------------------------------------------------------------
    ! !DESCRIPTION:

    ! !USES
    use ESMF
    use lilac_utils   ,  only : fld_list_type, fldsMax, create_fldlists 
    use lilac_utils   ,  only : atm2lnd_data1d_type , lnd2atm_data1d_type
    use lilac_utils   ,  only : atm2lnd_data2d_type , lnd2atm_data2d_type
    use atmos_cap     ,  only :         atmos_register
    !use lnd_shr_methods
    use lnd_comp_esmf ,  only :         lnd_register
    use cpl_mod       ,  only :         cpl_atm2lnd_register , cpl_lnd2atm_register

    use mpi           ,  only : MPI_COMM_WORLD, MPI_COMM_NULL, MPI_Init, MPI_FINALIZE, MPI_SUCCESS
    use shr_pio_mod   ,  only : shr_pio_init1, shr_pio_init2

    use clm_varctl    ,  only : iulog
    use spmdMod       ,  only : masterproc
    implicit none

    !TODO (NS,2019-08-07):
    ! We will move this later to lnd_cap (ctsm_cap) and atmos_cap
    !use atmos_cap     ,  only :         a2l_fldnum
    integer    , public ,  parameter                  :: a2l_fldnum        =  17
    integer    , public ,  parameter                  :: l2a_fldnum        =  12

    public                                            :: lilac_init
    public                                            :: lilac_run

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
    type(ESMF_GridComp)                               :: atmos_gcomp
    type(ESMF_GridComp)                               :: land_gcomp
    type(ESMF_CplComp)                                :: cpl_atm2lnd_comp
    type(ESMF_CplComp)                                :: cpl_lnd2atm_comp
    type(ESMF_State)                                  :: atm2lnd_l_state , atm2lnd_a_state
    type(ESMF_State)                                  :: lnd2atm_a_state, lnd2atm_l_state

    !========================================================================
    contains
    !========================================================================

    subroutine lilac_init( atm2lnd1d, atm2lnd2d, lnd2atm1d, lnd2atm2d)

        use atmos_cap , only :          a2c_fldlist , c2a_fldlist
        use lnd_cap   , only :          l2c_fldlist , c2l_fldlist
 
        character(len=*), parameter                      :: subname=trim(modname)//': [lilac_init] '

        ! input/output variables
        type(atm2lnd_data1d_type), intent(in), optional  :: atm2lnd1d
        type(atm2lnd_data2d_type), intent(in), optional  :: atm2lnd2d
        type(lnd2atm_data1d_type), intent(in), optional  :: lnd2atm1d
        type(lnd2atm_data2d_type), intent(in), optional  :: lnd2atm2d

        ! local variables

        type(ESMF_State)                                 :: importState, exportState

        !character(len=*)                                :: atm_mesh_filepath   !!! For now this is hardcoded in the atmos init

        integer                                          :: rc         , userRC 
        character(len=ESMF_MAXSTR)                       :: gcname1    , gcname2   !    Gridded components names
        character(len=ESMF_MAXSTR)                       :: ccname1    , ccname2   !    Coupling components names


        ! Namelist and related variables
        integer            ::  fileunit
        integer            ::  i_max, j_max
        real(ESMF_KIND_R8) ::  x_min, x_max, y_min, y_max
        integer            ::  s_month, s_day, s_hour, s_min
        integer            ::  e_month, e_day, e_hour, e_min
        namelist /input/ i_max, j_max, x_min, x_max, y_min, y_max, &
                                     s_month, s_day, s_hour, s_min, &
                                     e_month, e_day, e_hour, e_min


        integer                                         :: COMP_COMM
        integer                                         :: ierr
        integer                                         :: ntasks,mytask ! mpicom size and rank
        integer                                         :: ncomps = 1      ! land only
        integer                                         :: n
        integer                                         :: i
        integer, parameter                              :: debug = 1   !-- internal debug level
        !!! above: https://github.com/yudong-tian/LIS-CLM4.5SP/blob/8cec515a628325c73058cfa466db63210cd562ac/pio-xlis-bld/xlis_main.F90


        character(len=128)                                    :: fldname
        integer, parameter     :: begc = 1   !-- internal debug level
        integer, parameter     :: endc = 3312/4/2/2   !-- internal debug level
        character(*),parameter :: F02 =   "('[lilac_mod]',a,i5,2x,d26.19)"
        !------------------------------------------------------------------------
        ! Initialize return code
        rc = ESMF_SUCCESS

        if (masterproc) then
            print *,  "---------------------------------------"
            print *,  "    Lilac Demo Application Start       "
            print *,  "---------------------------------------"
        end if

        !-----------------------------------------------------------------------------
        ! Initiallize MPI
        !-----------------------------------------------------------------------------

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

        !-----------------------------------------------------------------------------
        ! Initialize PIO
        !-----------------------------------------------------------------------------

        ! this is coming from
        ! /glade/work/mvertens/ctsm.nuopc/cime/src/drivers/nuopc/drivers/cime/esmApp.F90
        ! with call shr_pio_init1(8, "drv_in", COMP_COMM)

        ! For planned future use of async io using pio2.  The IO tasks are seperated from the compute tasks here
        ! and COMP_COMM will be MPI_COMM_NULL on the IO tasks which then call shr_pio_init2 and do not return until
        ! the model completes.  All other tasks call ESMF_Initialize.  8 is the maximum number of component models
        ! supported

        call shr_pio_init1(ncomps, "drv_in", COMP_COMM)
        ! NS Question: How many should ncomps (above 1) be??????

        if (COMP_COMM .eq. MPI_COMM_NULL) then
            !call shr_pio_init2(
            call mpi_finalize(ierror=rc)
            stop
        endif

        !-------------------------------------------------------------------------
        ! Initialize ESMF, set the default calendar and log type.
        !-------------------------------------------------------------------------
        call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN,logappendflag=.false., rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogSet(flush=.true.)
        call ESMF_LogWrite(subname//".........................", ESMF_LOGMSG_INFO)
        call ESMF_LogWrite(subname//"Initializing ESMF        ", ESMF_LOGMSG_INFO)

        !-------------------------------------------------------------------------
        ! Read in configuration data -- namelist.input from host atmosphere(wrf)
        !-------------------------------------------------------------------------
        ! Read in namelist file ...
        call ESMF_UtilIOUnitGet(unit=fileunit, rc=rc)     ! get an available Fortran unit number
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        if (masterproc) then
            print *,  "---------------------------------------"
        end if

        open(fileunit, status="old", file="namelist_lilac",  action="read", iostat=rc)

        if (rc .ne. 0) then
            call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_OPEN, msg="Failed to open namelist file 'namelist'",  line=__LINE__, file=__FILE__)
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
        endif
        read(fileunit, input)
        continue
        close(fileunit)

        !-------------------------------------------------------------------------
        ! Create Field lists -- Basically create a list of fields and add a default
        ! value to them.
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        !            !---- from atm ----! a2c_fldlist & c2l_fldlist
        !-------------------------------------------------------------------------
        !allocate (a2c_fldlist(a2l_fldnum))
        !allocate (c2l_fldlist(a2l_fldnum))

        !-------------------------------------------------------------------------
        !            !---- from land ----! l2c_fldlist & c2a_fldlist
        !-------------------------------------------------------------------------
        !allocate (c2a_fldlist(l2a_fldnum))
        !allocate (l2c_fldlist(l2a_fldnum))

        allocate (a2c_fldlist(fldsMax))
        allocate (c2a_fldlist(fldsMax))

        allocate (l2c_fldlist(fldsmax))
        allocate (c2l_fldlist(fldsmax))

        if (masterproc) then
            print *, "creating empty field lists !"
        end if

        call ESMF_LogWrite(subname//"fielldlists are allocated!", ESMF_LOGMSG_INFO)

        ! create field lists  
        call create_fldlists(a2c_fldlist, c2a_fldlist,l2c_fldlist, c2l_fldlist)
        call ESMF_LogWrite(subname//"fielldlists are created!", ESMF_LOGMSG_INFO)

        !-------------------------------------------------------------------------
        !            !---- from atm ----! a2c_fldlist filling the arrayptr..
        !-------------------------------------------------------------------------


        ! FIXME: This should go to the demo_driver or real atmosphere......
        !allocate( a2c_fldlist(fldsmax)%farrayptr1d(1728))
        !do n = 1,a2l_fldnum
        !    print *, " index is ", n
        !    a2c_fldlist(1)%farrayptr1d(:) = 300.0
        !end do 

        a2c_fldlist(1)%farrayptr1d   =>  atm2lnd1d%Sa_z
        a2c_fldlist(2)%farrayptr1d   =>  atm2lnd1d%Sa_topo

        !if (masterproc .and. debug > 0) then
        fldname = 'Sa_topo'
        do i=begc, endc
            write (iulog,F02)'import: nstep, n, '//trim(fldname)//' = ',i, a2c_fldlist(2)%farrayptr1d(i)
        end do
        !end if
        a2c_fldlist(3)%farrayptr1d   =>  atm2lnd1d%Sa_u
        a2c_fldlist(4)%farrayptr1d   =>  atm2lnd1d%Sa_v
        a2c_fldlist(5)%farrayptr1d   =>  atm2lnd1d%Sa_ptem
        a2c_fldlist(6)%farrayptr1d   =>  atm2lnd1d%Sa_pbot
        a2c_fldlist(7)%farrayptr1d   =>  atm2lnd1d%Sa_tbot
        a2c_fldlist(8)%farrayptr1d   =>  atm2lnd1d%Sa_shum

        a2c_fldlist(9)%farrayptr1d   =>  atm2lnd1d%Faxa_lwdn
        a2c_fldlist(10)%farrayptr1d  =>  atm2lnd1d%Faxa_rainc
        a2c_fldlist(11)%farrayptr1d  =>  atm2lnd1d%Faxa_rainl
        a2c_fldlist(12)%farrayptr1d  =>  atm2lnd1d%Faxa_snowc
        a2c_fldlist(13)%farrayptr1d  =>  atm2lnd1d%Faxa_snowl

        a2c_fldlist(14)%farrayptr1d  =>  atm2lnd1d%Faxa_swndr
        a2c_fldlist(15)%farrayptr1d  =>  atm2lnd1d%Faxa_swvdr
        a2c_fldlist(16)%farrayptr1d  =>  atm2lnd1d%Faxa_swndf
        a2c_fldlist(17)%farrayptr1d  =>  atm2lnd1d%Faxa_swvdf
        !-------------------------------------------------------------------------

        ! should I point to zero???

        c2a_fldlist(1)%farrayptr1d   =>  lnd2atm1d%Sl_lfrin
        c2a_fldlist(2)%farrayptr1d   =>  lnd2atm1d%Sl_t
        c2a_fldlist(3)%farrayptr1d   =>  lnd2atm1d%Sl_tref
        c2a_fldlist(4)%farrayptr1d   =>  lnd2atm1d%Sl_qref
        c2a_fldlist(5)%farrayptr1d   =>  lnd2atm1d%Sl_avsdr
        c2a_fldlist(6)%farrayptr1d   =>  lnd2atm1d%Sl_anidr
        c2a_fldlist(7)%farrayptr1d   =>  lnd2atm1d%Sl_avsdf
        c2a_fldlist(8)%farrayptr1d   =>  lnd2atm1d%Sl_anidf

        c2a_fldlist(9)%farrayptr1d   =>  lnd2atm1d%Sl_snowh
        c2a_fldlist(10)%farrayptr1d  =>  lnd2atm1d%Sl_u10
        c2a_fldlist(11)%farrayptr1d  =>  lnd2atm1d%Sl_fv
        c2a_fldlist(12)%farrayptr1d  =>  lnd2atm1d%Sl_ram1



        ! ========================================================================

        !-------------------------------------------------------------------------
        ! Create Gridded Component!  --  atmosphere ( atmos_cap)
        !-------------------------------------------------------------------------
        gcname1 = " Atmosphere or Atmosphere Cap"
        atmos_gcomp = ESMF_GridCompCreate(name=gcname1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"Created "//trim(gcname1)//" component", ESMF_LOGMSG_INFO)
        print *, "Atmosphere Gridded Component Created!"

        !-------------------------------------------------------------------------
        ! Create Gridded Component!  --- CTSM land       ( land_capX )
        !-------------------------------------------------------------------------
        gcname2 = " Land ctsm "
        land_gcomp = ESMF_GridCompCreate(name=gcname2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"Created "//trim(gcname2)//" component", ESMF_LOGMSG_INFO)
        print *, "  Land (ctsm) Gridded Component Created!"

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
        call ESMF_GridCompSetServices(atmos_gcomp, userRoutine=atmos_register, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"  atmos SetServices finished!", ESMF_LOGMSG_INFO)
        print *, "  Atmosphere Gridded Component SetServices finished!"
        !-------------------------------------------------------------------------
        ! Register section -- set services -- land cap
        !-------------------------------------------------------------------------
        call ESMF_GridCompSetServices(land_gcomp, userRoutine=lnd_register, rc=rc)
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
        !  Clock is initialized here from namelist.input from WRF..... still we
        !  are looping over time from host atmosphere
        !-------------------------------------------------------------------------
        calendar = ESMF_CalendarCreate(name='lilac_drv_NOLEAP', calkindflag=ESMF_CALKIND_NOLEAP, rc=rc )
        call ESMF_TimeIntervalSet(TimeStep, s=2, rc=rc) ! time step every 2second
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        !call ESMF_TimeSet(startTime, yy=2003, mm=s_month, dd=s_day, h=s_hour, m=s_min, s=0, rc=rc)
        !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        !call ESMF_TimeSet(stopTime,  yy=2003, mm=e_month, dd=e_day, h=e_hour, m=e_min, s=0, rc=rc)
        !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        !clock = ESMF_ClockCreate(timeStep=timeStep, startTime=startTime, stopTime=stopTime, rc=rc) 
        !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_TimeSet(StartTime, yy=2000, mm=1, dd=1 , s=0, calendar=Calendar, rc=rc)
        call ESMF_TimeSet(StopTime , yy=2000, mm=03, dd=01, s=0, calendar=Calendar, rc=rc)
        !call ESMF_TimeIntervalSet(TimeStep, s=3600, rc=rc)
        call ESMF_TimeIntervalSet(TimeStep, s=1800, rc=rc)
        clock = ESMF_ClockCreate(name='lilac_drv_EClock', TimeStep=TimeStep, startTime=StartTime, RefTime=StartTime, stopTime=stopTime, rc=rc)

        print *,  "---------------------------------------"
        !call ESMF_ClockPrint (clock, rc=rc)
        print *,  "======================================="
        !call ESMF_CalendarPrint ( calendar , rc=rc)
        print *,  "---------------------------------------"

        ! ========================================================================

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

        call ESMF_GridCompInitialize(atmos_gcomp, importState=lnd2atm_a_state, exportState=atm2lnd_a_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"atmos_cap or atmos_gcomp initialized", ESMF_LOGMSG_INFO)

        call ESMF_GridCompInitialize(land_gcomp       , importState=atm2lnd_l_state, exportState=lnd2atm_l_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"lnd_cap or land_gcomp initialized", ESMF_LOGMSG_INFO)

        ! All 4 states that are module variables are no longer empty - have been initialized

        call ESMF_CplCompInitialize(cpl_atm2lnd_comp, importState=atm2lnd_a_state, exportState=atm2lnd_l_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"coupler :: cpl_atm2lnd_comp initialized", ESMF_LOGMSG_INFO)
        print *, "coupler :: cpl_atm2lnd_comp initialize finished" !, rc =", rc

        call ESMF_CplCompInitialize(cpl_lnd2atm_comp, importState=lnd2atm_l_state, exportState=lnd2atm_a_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"coupler :: cpl_lnd2atm_comp initialized", ESMF_LOGMSG_INFO)
        print *, "coupler :: cpl_lnd2atm_comp initialize finished" !, rc =", rc

    end subroutine lilac_init

    !========================================================================

    subroutine lilac_run( )

        use atmos_cap, only : a2c_fldlist, c2a_fldlist
        use lnd_cap,   only : l2c_fldlist, c2l_fldlist

        character(len=*), parameter                      :: subname=trim(modname)//': [lilac_run] '
        type(ESMF_State)                                 :: importState, exportState

        ! local variables
        integer                                          :: rc, userRC
        character(len=ESMF_MAXSTR)                       :: gcname1, gcname2   !    Gridded components names
        character(len=ESMF_MAXSTR)                       :: ccname1, ccname2   !    Coupling components names
        !integer, parameter                              :: fldsMax = 100

        ! input/output variables
        !type(atm2lnd_data1d_type), intent(in), optional  :: atm2lnd1d
        !type(atm2lnd_data2d_type), intent(in), optional  :: atm2lnd2d
        !type(lnd2atm_data1d_type), intent(in), optional  :: lnd2atm1d
        !type(lnd2atm_data2d_type), intent(in), optional  :: lnd2atm2d

        type (ESMF_Clock)                                 :: local_clock

        !------------------------------------------------------------------------
        ! Initialize return code
        rc = ESMF_SUCCESS

        print *,  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print *,  "               Lilac Run               "
        print *,  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

        !-------------------------------------------------------------------------
        ! Create a local clock from the general clock!
        !-------------------------------------------------------------------------

        local_clock = ESMF_ClockCreate(clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        print *, "Run Loop Start time"
        !call ESMF_ClockPrint(local_clock, options="currtime string", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        !-------------------------------------------------------------------------
        ! We are running components in this order:
        ! 1- atmos_cap   2- cpl_atm2lnd
        ! 3-   lnd_cap   4- cpl_lnd2atm
        !-------------------------------------------------------------------------
        ! lilac run the RunComponent phase in a time loop

        !!! if we want to loop through clock in atmos cap. 
        !do while (.NOT. ESMF_ClockIsStopTime(local_clock, rc=rc))
            call ESMF_GridCompRun(atmos_gcomp, importState=lnd2atm_a_state, exportState=atm2lnd_a_state, clock=local_clock, rc=rc, userRC=userRC)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
            if (ESMF_LogFoundError(rcToCheck=userRC, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
            call ESMF_LogWrite(subname//"atmos_cap or atmos_gcomp is running", ESMF_LOGMSG_INFO)
            print *, "Running atmos_cap gridded component , rc =", rc

            call ESMF_CplCompRun(cpl_atm2lnd_comp, importState=atm2lnd_a_state, exportState=atm2lnd_l_state, clock=local_clock, rc=rc , userRC=userRC)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
            if (ESMF_LogFoundError(rcToCheck=userRC, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
            call ESMF_LogWrite(subname//"running cpl_atm2lnd_comp ", ESMF_LOGMSG_INFO)
            print *, "Running coupler component..... cpl_atm2lnd_comp , rc =", rc

            call ESMF_GridCompRun(land_gcomp,  importState=atm2lnd_l_state, exportState=lnd2atm_l_state, clock=local_clock, rc=rc, userRC=userRC)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
            if (ESMF_LogFoundError(rcToCheck=userRC, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
            call ESMF_LogWrite(subname//"lnd_cap or land_gcomp is running", ESMF_LOGMSG_INFO)
            print *, "Running lnd_cap gridded component  , rc =", rc

            call ESMF_CplCompRun(cpl_lnd2atm_comp, importState=lnd2atm_l_state, exportState=lnd2atm_a_state, clock=local_clock, rc=rc, userRC=userRC)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
            if (ESMF_LogFoundError(rcToCheck=userRC, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
            call ESMF_LogWrite(subname//"running cpl_lnd2atm_comp ", ESMF_LOGMSG_INFO)
            print *, "Running coupler component..... cpl_lnd2atm_comp , rc =", rc

            ! Advance the time
            call ESMF_ClockAdvance(local_clock, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
            call ESMF_LogWrite(subname//"time is icremented now... (ClockAdvance)", ESMF_LOGMSG_INFO)
            print *, "time is icremented now... (ClockAdvance) , rc =", rc

        !end do 

    end subroutine lilac_run


    subroutine lilac_final( )

        use atmos_cap, only : a2c_fldlist, c2a_fldlist
        use lnd_cap,   only : l2c_fldlist, c2l_fldlist


        character(len=*), parameter                      :: subname=trim(modname)//': [lilac_final] '
        type(ESMF_State)                                 :: importState, exportState

        ! local variables
        integer                                          :: rc, userRC
        character(len=ESMF_MAXSTR)                       :: gcname1, gcname2   !    Gridded components names
        character(len=ESMF_MAXSTR)                       :: ccname1, ccname2   !    Coupling components names
        !integer, parameter                              :: fldsMax = 100

        !------------------------------------------------------------------------
        !------------------------------------------------------------------------
        ! Initialize return code
        rc = ESMF_SUCCESS

        print *,  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print *,  "        Lilac Finalizing               "
        print *,  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        !-------------------------------------------------------------------------
        ! Gridded Component Finalizing!     ---       atmosphere 
        !-------------------------------------------------------------------------
        call ESMF_GridCompFinalize(atmos_gcomp, importState=lnd2atm_a_state, exportState=atm2lnd_a_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"atmos_cap or atmos_gcomp is running", ESMF_LOGMSG_INFO)
        print *, "Finalizing atmos_cap gridded component , rc =", rc

        !-------------------------------------------------------------------------
        ! Coupler component Finalizing     ---    coupler atmos to land
        !-------------------------------------------------------------------------
        call ESMF_CplCompFinalize(cpl_atm2lnd_comp, importState=atm2lnd_a_state, exportState=atm2lnd_l_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"running cpl_atm2lnd_comp ", ESMF_LOGMSG_INFO)
        print *, "Finalizing coupler component..... cpl_atm2lnd_comp , rc =", rc

        !-------------------------------------------------------------------------
        ! Gridded Component Finalizing!     ---       land 
        !-------------------------------------------------------------------------
        call ESMF_GridCompFinalize(land_gcomp,  importState=atm2lnd_l_state, exportState=lnd2atm_l_state, clock=clock, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"lnd_cap or land_gcomp is running", ESMF_LOGMSG_INFO)
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

        call ESMF_GridCompDestroy(atmos_gcomp, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
        call ESMF_GridCompDestroy(land_gcomp, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
        call ESMF_CplCompDestroy(cpl_atm2lnd_comp, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
        call ESMF_CplCompDestroy(cpl_lnd2atm_comp, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

        call ESMF_LogWrite(subname//".........................", ESMF_LOGMSG_INFO)
        print *, "end of Lilac Finalization routine"

        end subroutine lilac_final



    end module lilac_mod

