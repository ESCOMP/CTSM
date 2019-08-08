module lilac_utils

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    !
    !!! NS: THIS IS FROM JH WORK

    use ESMF

    implicit none

    public fldlist_add , create_fldlists

    integer, parameter              :: fldsMax = 100

    character(*) , parameter                          :: modname     = "lilac_utils"
    ! !PUBLIC TYPES:
    type                            :: fld_list_type
        character(len=128)          :: stdname
        real*8                      :: default_value
        character(len=128)          :: units
        real(ESMF_KIND_R8), pointer :: farrayptr1d(:)  ! this will be filled in by lilac when it gets its data from the host atm
        real(ESMF_KIND_R8), pointer :: farrayptr2d(:,:)  ! this will be filled in by lilac when it gets its data from the host atm
        integer                     :: ungridded_lbound = 0
        integer                     :: ungridded_ubound = 0
    end type                           fld_list_type

    !!! 1d for when we have mesh and 2d for when we have grids....

    type                   :: atm2lnd_data1d_type
        real*8, pointer    :: uwind (:)
        real*8, pointer    :: vwind (:)
        real*8, pointer    :: tbot  (:)
    end type                  atm2lnd_data1d_type

    type                   :: lnd2atm_data1d_type
        real*8, pointer    :: lwup  (:)
        real*8, pointer    :: taux  (:)
        real*8, pointer    :: tauy  (:)
    end type                  lnd2atm_data1d_type

    type                   :: atm2lnd_data2d_type
        real*8, pointer    :: uwind (:,:)
        real*8, pointer    :: vwind (:,:)
        real*8, pointer    :: tbot  (:,:)
    end type                  atm2lnd_data2d_type

    type                   :: lnd2atm_data2d_type
        real*8, pointer    :: lwup  (:,:)
        real*8, pointer    :: taux  (:,:)
        real*8, pointer    :: tauy  (:,:)
     end type                 lnd2atm_data2d_type

    type                   :: this_clock
        integer, pointer   :: yy
        integer, pointer   :: mm
        integer, pointer   :: dd
        integer, pointer   :: hh
        integer, pointer   :: mn
        integer, pointer   :: ss
     end type                 this_clock
    !===============================================================================
    contains
    !===============================================================================

    subroutine fldlist_add(num, fldlist, stdname, default_value, units, ungridded_lbound, ungridded_ubound)
        ! This adds a field to a fieldlist!
        ! input/output variables
        integer,                     intent(inout) :: num
        type(fld_list_type),         intent(inout) :: fldlist(:)
        character(len=*),            intent(in)    :: stdname
        real, optional,              intent(in)    :: default_value
        character(len=*), optional,  intent(in)    :: units
        integer,          optional,  intent(in)    :: ungridded_lbound
        integer,          optional,  intent(in)    :: ungridded_ubound

        ! local variables
        integer :: rc
        character(len=*), parameter                :: subname=trim(modname)//':[fldlist_add]'
        !-------------------------------------------------------------------------------
       call ESMF_LogWrite(subname//"inside fldlist_add!", ESMF_LOGMSG_INFO)

        ! Set up a list of field information
        num = num + 1
        if (num > fldsMax) then
            call ESMF_LogWrite(subname//"?!", ESMF_LOGMSG_INFO)
            call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
                             ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
               return
        endif

        fldlist(num)%stdname = trim(stdname)

        if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
           fldlist(num)%ungridded_lbound = ungridded_lbound
           fldlist(num)%ungridded_ubound = ungridded_ubound
        end if

        if(present(default_value)) then
           fldlist(num)%default_value = default_value
        else
           fldlist(num)%default_value = 0.
        end if
        if(present(units)) then
           fldlist(num)%units = trim(units)
        else
           fldlist(num)%units = ""
        end if

    end subroutine fldlist_add

    !subroutine create_fldlists(a2c_fldlist, c2l_fldlist, l2c_fldlist, c2a_fldlist, rof_prognostic, glc_present )
    subroutine create_fldlists(a2c_fldlist, c2l_fldlist, l2c_fldlist, c2a_fldlist)
        
        ! add all the necessary fields one by one to the fieldlist 
        type(fld_list_type),        intent(inout) :: a2c_fldlist(fldsMax)
        type(fld_list_type),        intent(inout) :: c2a_fldlist(fldsMax)
        type(fld_list_type),        intent(inout) :: l2c_fldlist(fldsMax)
        type(fld_list_type),        intent(inout) :: c2l_fldlist(fldsMax)

        ! I use this as an index!
        integer                                   :: fldsToLnd_num != 0  ! From atmosphere to land (c2l)
        integer                                   :: fldsFrLnd_num != 0  ! From land to atmosphere (l2c)
        integer                                   :: fldsToAtm_num != 0  ! From land to atmosphere (c2a)
        integer                                   :: fldsFrAtm_num != 0  ! From atmosphere to land (a2c)
        integer, parameter                        :: fldsMax = 100


        ! TODO (NS) : Should we move these to the land cap????
        logical                                   :: glc_present    ! .true. => running with a non-stub GLC model
        logical                                   :: rof_prognostic ! .true. => running with a prognostic ROF model

        character(len=*), parameter                :: subname=trim(modname)//':[create_fldlists]'
        ! TODO (NS) : I should add default value and units here.....

        fldsToLnd_num= 0
        fldsFrLnd_num= 0
        fldsToAtm_num= 0
        fldsFrAtm_num= 0

       call ESMF_LogWrite(subname//"is called!", ESMF_LOGMSG_INFO)

        !-------------------------------------------------------------------------
        !            !---- from atm ----! a2c_fldlist &  c2l_fldlist
        !-------------------------------------------------------------------------
        !--------------------------a2c_fldlist------------------------------------
        ! from atm - states
        !call fldlist_add(fldsToLnd_num  , a2c_fldlist , 'Sa_z'         )
        !call fldlist_add(fldsToLnd_num  , a2c_fldlist , 'Sa_topo'      )
        call fldlist_add(fldsToLnd_num  , a2c_fldlist , 'Sa_z'         ,   default_value=30.0      , units='m/s')
        call fldlist_add(fldsToLnd_num  , a2c_fldlist , 'Sa_topo'      ,   default_value=10.0      , units='m')
        call fldlist_add(fldsToLnd_num  , a2c_fldlist , 'Sa_u'           , default_value=0.0      , units='m/s')
        call fldlist_add(fldsToLnd_num  , a2c_fldlist , 'Sa_v'           , default_value=0.0      , units='m/s')
        call fldlist_add(fldsToLnd_num  , a2c_fldlist , 'Sa_ptem'        , default_value=280.0    , units='degK')
        call fldlist_add(fldsToLnd_num  , a2c_fldlist , 'Sa_pbot'        , default_value=100100.0 , units='pa'  )
        call fldlist_add(fldsToLnd_num  , a2c_fldlist , 'Sa_tbot'        , default_value=280.0    , units='degk' )
        call fldlist_add(fldsToLnd_num  , a2c_fldlist , 'Sa_shum'        , default_value=0.0004   , units='kg/kg' )
        !call fldlist_add(fldsToLnd_num , a2c_fldlist , 'Sa_methane'   )

       call ESMF_LogWrite(subname//"from atmosphere states are added!" , ESMF_LOGMSG_INFO)





        ! from atm - fluxes
        call fldlist_add(fldsToLnd_num , a2c_fldlist , 'Faxa_lwdn'  , default_value=200.0  , units='W/m2'   )
        call fldlist_add(fldsToLnd_num , a2c_fldlist , 'Faxa_rainc' , default_value=4.0e-8 , units='kg/m2s' )
        call fldlist_add(fldsToLnd_num , a2c_fldlist , 'Faxa_rainl' , default_value=3.0e-8 , units='kg/m2s' )
        call fldlist_add(fldsToLnd_num , a2c_fldlist , 'Faxa_snowc' , default_value=1.0e-8 , units='kg/m2s' )
        call fldlist_add(fldsToLnd_num , a2c_fldlist , 'Faxa_snowl' , default_value=2.0e-8 , units='kg/m2s' )
        call fldlist_add(fldsToLnd_num , a2c_fldlist , 'Faxa_swndr' , default_value=100.0  , units='W/m2'   )
        call fldlist_add(fldsToLnd_num , a2c_fldlist , 'Faxa_swvdr' , default_value=90.0   , units='W/m2'   )
        call fldlist_add(fldsToLnd_num , a2c_fldlist , 'Faxa_swndf' , default_value=20.0   , units='W/m2'   )
        call fldlist_add(fldsToLnd_num , a2c_fldlist , 'Faxa_swvdf' , default_value=40.0   , units='W/m2'   )

        call ESMF_LogWrite(subname//"from atmosphere fluxes are added!", ESMF_LOGMSG_INFO)

        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_bcphidry')
        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_bcphodry')
        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_bcphiwet')

        !--------------------------c2l_fldlist------------------------------------
        ! from atm - states
        call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Sa_z'         )
        call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Sa_topo'      )
        call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Sa_u'         )
        call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Sa_v'         )
        call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Sa_ptem'      )
        call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Sa_pbot'      )
        call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Sa_tbot'      )
        call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Sa_shum'      )
        !call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Sa_methane'   )
       call ESMF_LogWrite(subname//"from atmosphere states are added!", ESMF_LOGMSG_INFO)

        ! from atm - fluxes
        call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Faxa_lwdn'    )
        call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Faxa_rainc'   )
        call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Faxa_rainl'   )
        call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Faxa_snowc'   )
        call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Faxa_snowl'   )
        call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Faxa_swndr'   )
        call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Faxa_swvdr'   )
        call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Faxa_swndf'   )
        call fldlist_add(fldsToLnd_num, c2l_fldlist, 'Faxa_swvdf'   )
        call ESMF_LogWrite(subname//"from atmosphere fluxes are added!", ESMF_LOGMSG_INFO)

        !-------------------------------------------------------------------------
        !            !---- from lnd ----! l2c_fldlist &  c2a_fldlist
        !-------------------------------------------------------------------------
        !--------------------------l2c_fldlist------------------------------------
        ! export land states
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Sl_lfrin'      )
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Sl_t'          )
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Sl_tref'       )
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Sl_qref'       )
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Sl_avsdr'      )
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Sl_anidr'      )
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Sl_avsdf'      )
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Sl_anidf'      )
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Sl_snowh'      )
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Sl_u10'        )
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Sl_fv'         )
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Sl_ram1'       )
        call ESMF_LogWrite(subname//"l2c: from land states are added!", ESMF_LOGMSG_INFO)

        rof_prognostic = .false.
        ! export fluxes to river
        if (rof_prognostic) then
            call ESMF_LogWrite(subname//"Okay we are in rof_prognostic", ESMF_LOGMSG_INFO)
           call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Flrl_rofsur'   )
            call ESMF_LogWrite(subname//"Okay we are in rof_prognostic 13", ESMF_LOGMSG_INFO)
           call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Flrl_rofgwl'   )
            call ESMF_LogWrite(subname//"Okay we are in rof_prognostic 14", ESMF_LOGMSG_INFO)
           call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Flrl_rofsub'   )
            call ESMF_LogWrite(subname//"Okay we are in rof_prognostic 15", ESMF_LOGMSG_INFO)
           call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Flrl_rofi'     )
            call ESMF_LogWrite(subname//"Okay we are in rof_prognostic 16", ESMF_LOGMSG_INFO)
           call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Flrl_irrig'    )
            call ESMF_LogWrite(subname//"l2c: from land states are added for rof_prognostic!", ESMF_LOGMSG_INFO)
        end if

        ! export fluxes to atm
        call ESMF_LogWrite(subname//"l2c: now adding fluxes to atmosphere!", ESMF_LOGMSG_INFO)
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Fall_taux'     )
        call ESMF_LogWrite(subname//"l2c: Fall_taux!", ESMF_LOGMSG_INFO)
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Fall_tauy'     )
        call ESMF_LogWrite(subname//"l2c: Fall_taux!", ESMF_LOGMSG_INFO)
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Fall_lat'      )
        call ESMF_LogWrite(subname//"l2c: Fall_lat!", ESMF_LOGMSG_INFO)
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Fall_sen'      )
        call ESMF_LogWrite(subname//"l2c: Fall_sen!", ESMF_LOGMSG_INFO)
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Fall_lwup'     )
        call ESMF_LogWrite(subname//"l2c: Fall_lwup!", ESMF_LOGMSG_INFO)
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Fall_evap'     )
        call ESMF_LogWrite(subname//"l2c: Fall_evap!", ESMF_LOGMSG_INFO)
        call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Fall_swnet'    )
        call ESMF_LogWrite(subname//"l2c: Fall_lat!", ESMF_LOGMSG_INFO)
        call ESMF_LogWrite(subname//"l2c: from land fluxes are added!", ESMF_LOGMSG_INFO)

        ! call fldlist_add(fldsFrLnd_num, l2c_fldlist, 'Fall_methane'  )


        !--------------------------c2a_fldlist------------------------------------
        ! export land states
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Sl_lfrin'      )
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Sl_t'          )
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Sl_tref'       )
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Sl_qref'       )
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Sl_avsdr'      )
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Sl_anidr'      )
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Sl_avsdf'      )
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Sl_anidf'      )
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Sl_snowh'      )
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Sl_u10'        )
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Sl_fv'         )
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Sl_ram1'       )


        ! export fluxes to river
        if (rof_prognostic) then
           call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Flrl_rofsur'   )
           call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Flrl_rofgwl'   )
           call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Flrl_rofsub'   )
           call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Flrl_rofi'     )
           call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Flrl_irrig'    )
        end if

        ! export fluxes to atm
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Fall_taux'     )
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Fall_tauy'     )
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Fall_lat'      )
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Fall_sen'      )
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Fall_lwup'     )
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Fall_evap'     )
        call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Fall_swnet'    )

        ! call fldlist_add(fldsFrLnd_num, c2a_fldlist, 'Fall_methane'  )



            !call fldlist_add(fldsToCpl_num, fldsToCpl, 'atmos2lnd_var', default_value=0.0, units='m')
        ! from lnd
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'lnd2atmos_var', default_value=0.0, units='m')


        !  sets the module variable memory in atmos_cap.F9 print *,      a2c_fldlist(1)%stdname
        !!! First from atmosphere to land fields
        ! import fields
        ! call fldlist_add(fldsFrCpl_num, fldsFrCpl, trim(flds_scalar_name))

        !call fldlist_add(fldsToLnd_num, fldsToLnd, trim(flds_scalar_name))

        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_topo')
        !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_u',       default_value=0.0,       units='m/s')
        !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_v',       default_value=0.0,       units='m/s')
        !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_ptem',    default_value=280.0,     units= 'degK')
        !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_pbot',    default_value=100100.0,  units='Pa')
        !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_tbot',    default_value=280.0,     units='degK')
        !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_shum',    default_value=0.0004,    units='kg/kg')
        !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_methane'   )

        !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_lwdn',  default_value=200.0,     units='W/m2')
        !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_rainc', default_value=4.0e-8,    units='kg/m2s')
        !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_rainl', default_value=3.0e-8,    units='kg/m2s')
        !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_snowc', default_value=1.0e-8,    units='kg/m2s')
        !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_snowl', default_value=2.0e-8,    units='kg/m2s')
        !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_swndr', default_value=100.0,     units='W/m2')
        !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_swvdr', default_value=90.0,      units='W/m2')
        !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_swndf', default_value=20.0,      units='W/m2')
        !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_swvdf', default_value=40.0,       units='W/m2')
        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_bcphidry')
        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_bcphodry')
        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_bcphiwet')
        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_ocphidry')
        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_ocphodry')
        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_ocphiwet')
        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstdry1' )
        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstdry2' )
        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstdry3' )
        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstdry4' )
        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstwet1' )
        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstwet2' )
        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstwet3' )
        ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstwet4' )

        ! land states

        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_lfrin'      )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_t'          )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_tref'       )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_qref'       )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_avsdr'      )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_anidr'      )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_avsdf'      )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_anidf'      )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_snowh'      )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_u10'        )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_fv'         )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_ram1'       )

        ! fluxes to atm
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_taux'     )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_tauy'     )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_lat'      )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_sen'      )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_lwup'     )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_evap'     )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_swnet'    )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_flxdst1'  )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_flxdst2'  )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_flxdst3'  )
        !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_flxdst4'  )



        ! more: https://github.com/mvertens/ctsm/blob/ae02ffe25dbc4a85c769c9137b5b3d50f2843e89/src/cpl/nuopc/lnd_import_export.F90#L131
        end subroutine create_fldlists

end module lilac_utils
