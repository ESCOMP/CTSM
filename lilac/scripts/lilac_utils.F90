module lilac_utils
use ESMF
implicit none
!!! NS: THIS IS FROM JH WORK

    integer, parameter     :: fldsMax = 100

    public fldlist_add , create_fldlists

    type fld_list_type
        character(len=128) :: stdname
        real*8             :: default_value
        character(len=128) :: units
        real*8, pointer    :: farrayptr1d(:)  ! this will be filled in by lilac when it gets its data from the host atm
        real*8, pointer    :: farrayptr2d(:,:)  ! this will be filled in by lilac when it gets its data from the host atm
    end type fld_list_type

!===============================================================================
    contains
!===============================================================================

    subroutine fldlist_add(num, fldlist, stdname, default_value, units)
    ! This adds a field to a fieldlist!
        integer,                     intent(inout) :: num
        type(fld_list_type),         intent(inout) :: fldlist(:)
        character(len=*),            intent(in)    :: stdname
        real, optional,              intent(in)    :: default_value
        character(len=*), optional,  intent(in)    :: units

        ! local variables
        integer :: rc
        character(len=*), parameter :: subname='(fldlist_add)'
        !-------------------------------------------------------------------------------

        ! Set up a list of field information
        num = num + 1
        if (num > fldsMax) then
            call ESMF_LogWrite(subname//"?!", ESMF_LOGMSG_INFO)
        endif

        fldlist(num)%stdname = trim(stdname)

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

    subroutine create_fldlists(fldsFrCpl, fldsToCpl, fldsToCpl_num, fldsFrCpl_num)
    ! add all the necessary fields one by one to the fieldlist 
        type(fld_list_type),        intent(inout)      ::  fldsFrCpl(:)
        type(fld_list_type),        intent(inout)      ::  fldsToCpl(:)
        !integer,                    intent(out)        ::  fldsToCpl_num = 0
        !integer,                    intent(out)        ::  fldsFrCpl_num = 0
        integer :: fldsFrCpl_num, fldsToCpl_num

        ! from atm
        call fldlist_add(fldsToCpl_num, fldsToCpl, 'atmos2lnd_var', default_value=0.0, units='m')
        ! from lnd
        call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'lnd2atmos_var', default_value=0.0, units='m')



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
