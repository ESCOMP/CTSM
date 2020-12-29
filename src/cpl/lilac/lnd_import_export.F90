module lnd_import_export

  use ESMF
  use shr_kind_mod          , only : r8 => shr_kind_r8, cx=>shr_kind_cx, cxx=>shr_kind_cxx, cs=>shr_kind_cs
  use shr_sys_mod           , only : shr_sys_abort
  use shr_const_mod         , only : fillvalue=>SHR_CONST_SPVAL
  use clm_varctl            , only : iulog, ndep_from_cpl
  use clm_time_manager      , only : get_nstep
  use clm_instMod           , only : atm2lnd_inst, lnd2atm_inst, water_inst
  use domainMod             , only : ldomain
  use spmdMod               , only : masterproc
  use decompmod             , only : bounds_type
  use lnd2atmType           , only : lnd2atm_type
  use lnd2glcMod            , only : lnd2glc_type
  use atm2lndType           , only : atm2lnd_type
  use lnd_shr_methods       , only : chkerr
  use shr_megan_mod         , only : shr_megan_mechcomps_n  ! TODO: need to add a namelist read here (see https://github.com/ESCOMP/CTSM/issues/926)
  use lnd_import_export_utils, only : derive_quantities, check_for_errors, check_for_nans

  implicit none
  private ! except

  public  :: import_fields
  public  :: export_fields

  private :: state_getimport
  private :: state_setexport
  private :: state_getfldptr

  ! import fields
  character(*), parameter :: Sa_z                = 'Sa_z'
  character(*), parameter :: Sa_topo             = 'Sa_topo'
  character(*), parameter :: Sa_u                = 'Sa_u'
  character(*), parameter :: Sa_v                = 'Sa_v'
  character(*), parameter :: Sa_ptem             = 'Sa_ptem'
  character(*), parameter :: Sa_shum             = 'Sa_shum'
  character(*), parameter :: Sa_pbot             = 'Sa_pbot'
  character(*), parameter :: Sa_tbot             = 'Sa_tbot'
  character(*), parameter :: Sa_methane          = 'Sa_methane'
  character(*), parameter :: Faxa_rainc          = 'Faxa_rainc'
  character(*), parameter :: Faxa_rainl          = 'Faxa_rainl'
  character(*), parameter :: Faxa_snowc          = 'Faxa_snowc'
  character(*), parameter :: Faxa_snowl          = 'Faxa_snowl'
  character(*), parameter :: Faxa_lwdn           = 'Faxa_lwdn'
  character(*), parameter :: Faxa_swvdr          = 'Faxa_swvdr'
  character(*), parameter :: Faxa_swndr          = 'Faxa_swndr'
  character(*), parameter :: Faxa_swvdf          = 'Faxa_swvdf'
  character(*), parameter :: Faxa_swndf          = 'Faxa_swndf'
  character(*), parameter :: Faxa_bcphidry       = 'Faxa_bcphidry'
  character(*), parameter :: Faxa_bcphodry       = 'Faxa_bcphodry'
  character(*), parameter :: Faxa_bcphiwet       = 'Faxa_bcphiwet'
  character(*), parameter :: Faxa_ocphidry       = 'Faxa_ocphidry'
  character(*), parameter :: Faxa_ocphodry       = 'Faxa_ocphodry'
  character(*), parameter :: Faxa_ocphiwet       = 'Faxa_ocphiwet'
  character(*), parameter :: Faxa_dstwet1        = 'Faxa_dstwet1'
  character(*), parameter :: Faxa_dstwet2        = 'Faxa_dstwet2'
  character(*), parameter :: Faxa_dstwet3        = 'Faxa_dstwet3'
  character(*), parameter :: Faxa_dstwet4        = 'Faxa_dstwet4'
  character(*), parameter :: Faxa_dstdry1        = 'Faxa_dstdry1'
  character(*), parameter :: Faxa_dstdry2        = 'Faxa_dstdry2'
  character(*), parameter :: Faxa_dstdry3        = 'Faxa_dstdry3'
  character(*), parameter :: Faxa_dstdry3        = 'Faxa_dstdry4'
  character(*), parameter :: Faxa_ndep           = 'Faxa_ndep'

  ! export fields
  character(*), parameter :: Sl_t           = 'Sl_t'
  character(*), parameter :: Sl_snowh       = 'Sl_snowh'
  character(*), parameter :: Sl_avsdr       = 'Sl_avsdr'
  character(*), parameter :: Sl_anidr       = 'Sl_anidr'
  character(*), parameter :: Sl_avsdf       = 'Sl_avsdf'
  character(*), parameter :: Sl_anidf       = 'Sl_anidf'
  character(*), parameter :: Sl_tref        = 'Sl_tref'
  character(*), parameter :: Sl_qref        = 'Sl_qref'
  character(*), parameter :: Sl_u10         = 'Sl_u10'
  character(*), parameter :: Sl_ram1        = 'Sl_ram1'
  character(*), parameter :: Sl_fv          = 'Sl_fv'
  character(*), parameter :: Sl_z0m         = 'Sl_z0m'
  character(*), parameter :: Sl_soilw       = 'Sl_soilw'
  character(*), parameter :: Sl_ddvel       = 'Sl_ddvel'
  character(*), parameter :: Sl_fztop       = 'Sl_fztop'
  character(*), parameter :: Fall_taux      = 'Fall_taux'
  character(*), parameter :: Fall_tauy      = 'Fall_tauy'
  character(*), parameter :: Fall_lat       = 'Fall_lat'
  character(*), parameter :: Fall_sen       = 'Fall_sen'
  character(*), parameter :: Fall_lwup      = 'Fall_lwup'
  character(*), parameter :: Fall_evap      = 'Fall_evap'
  character(*), parameter :: Fall_swnet     = 'Fall_swnet'
  character(*), parameter :: Fall_flxdst    = 'Fall_flxdst'
  character(*), parameter :: Fall_methane   = 'Fall_methane'
  character(*), parameter :: Fall_voc       = 'Fall_voc'
  character(*), parameter :: Fall_fire      = 'Fall_fire'
  character(*), parameter :: Flrl_rofsur    = 'Flrl_rofsur'
  character(*), parameter :: Flrl_rofsub    = 'Flrl_rofsub'
  character(*), parameter :: Flrl_rofgwl    = 'Flrl_rofgwl'
  character(*), parameter :: Flrl_rofi      = 'Flrl_rofi'
  character(*), parameter :: Flrl_irrig     = 'Flrl_irrig'

  ! from atm->lnd
  integer                :: ndep_nflds               ! number  of nitrogen deposition fields from atm->lnd/ocn

  ! from lnd->atm
  integer                :: drydep_nflds             ! number of dry deposition velocity fields lnd-> atm
  integer                :: emis_nflds               ! number of fire emission fields from lnd-> atm
  integer, parameter     :: debug = 0                ! internal debug level

  character(*),parameter :: F01 = "('(lnd_import_export) ',a,i5,2x,i5,2x,d21.14)"
  character(*),parameter :: F02 = "('(lnd_import_export) ',a,i5,2x,i5,2x,d26.19)"
  character(*),parameter :: u_FILE_u = &
       __FILE__
  character(*),parameter :: modname =  "[lnd_import_export]: "

!===============================================================================
contains
!===============================================================================

  subroutine import_fields( importState, bounds, first_call, rc)

    !---------------------------------------------------------------------------
    ! Convert the input data from lilac to the land model
    !---------------------------------------------------------------------------

    ! input/output variables
    type(ESMF_State)                :: importState
    type(bounds_type) , intent(in)  :: bounds       ! bounds
    logical           , intent(in)  :: first_call   ! true if and only if this is the first time we're calling import_fields from the run method
    integer           , intent(out) :: rc

    ! local variables
    integer                   :: num
    integer                   :: begg, endg                             ! bounds
    integer                   :: g,i,k                                  ! indices
    real(r8)                  :: forc_rainc(bounds%begg:bounds%endg)    ! rainxy Atm flux mm/s
    real(r8)                  :: forc_rainl(bounds%begg:bounds%endg)    ! rainxy Atm flux mm/s
    real(r8)                  :: forc_snowc(bounds%begg:bounds%endg)    ! snowfxy Atm flux  mm/s
    real(r8)                  :: forc_snowl(bounds%begg:bounds%endg)    ! snowfxl Atm flux  mm/s
    real(r8)                  :: qsat_kg_kg                             ! saturation specific humidity (kg/kg)
    real(r8)                  :: forc_noy(bounds%begg:bounds%endg)
    real(r8)                  :: forc_nhx(bounds%begg:bounds%endg)
    character(len=*), parameter :: subname='(lnd_import_export:import_fields)'

    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Set bounds
    begg = bounds%begg; endg=bounds%endg

    if (first_call) then
       ! We only do this for the first call because we assume that the atmosphere's land
       ! mask is constant in time. To allow for a varying land mask in the atmosphere
       ! (doing checking each time), remove this first_call conditional. (It would be
       ! more straightforward to pass this and check it in initialization, but that would
       ! require atm-land communication in initialization, which currently isn't done
       ! with the LILAC coupler.)
       call check_atm_landfrac(importState, bounds, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Note: precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.

    !--------------------------
    ! Required atmosphere input fields
    !--------------------------

    call state_getimport_1d(importState, Sa_z      , atm2lnd_inst%forc_hgt_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Sa_topo   , atm2lnd_inst%forc_topo_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Sa_u      , atm2lnd_inst%forc_u_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Sa_v      , atm2lnd_inst%forc_v_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Sa_shum   , wateratm2lndbulk_inst%forc_q_not_downscaled_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Sa_ptem   , atm2lnd_inst%forc_th_not_downscaled_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Sa_pbot   , atm2lnd_inst%forc_pbot_not_downscaled_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Sa_tbot   , atm2lnd_inst%forc_t_not_downscaled_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_rainc, forc_rainc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_rainl, forc_rainl(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_snowc, forc_snowc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_snowl, forc_snowl(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_lwdn , atm2lnd_inst%forc_lwrad_not_downscaled_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_swvdr, atm2lnd_inst%forc_solad_grc(begg:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_swndr, atm2lnd_inst%forc_solad_grc(begg:,2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_swvdf, atm2lnd_inst%forc_solai_grc(begg:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_swndf, atm2lnd_inst%forc_solai_grc(begg:,2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ! Atmosphere prognostic/prescribed aerosol fields
    call state_getimport_1d(importState, Faxa_bcphidry, atm2lnd_inst%forc_aer_grc(begg:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_bcphodry, atm2lnd_inst%forc_aer_grc(begg:,2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_bcphiwet, atm2lnd_inst%forc_aer_grc(begg:,3), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_ocphidry, atm2lnd_inst%forc_aer_grc(begg:,4), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_ocphodry, atm2lnd_inst%forc_aer_grc(begg:,5), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_ocphiwet, atm2lnd_inst%forc_aer_grc(begg:,6), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_dstwet1, output=atm2lnd_inst%forc_aer_grc(begg:,7), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_dstwet2, output=atm2lnd_inst%forc_aer_grc(begg:,9), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_dstwet3, output=atm2lnd_inst%forc_aer_grc(begg:,11), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_dstwet4, output=atm2lnd_inst%forc_aer_grc(begg:,13), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport_1d(importState, Faxa_dstdry1, output=atm2lnd_inst%forc_aer_grc(begg:,8), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_dstdry2, output=atm2lnd_inst%forc_aer_grc(begg:,10), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_dstdry3, output=atm2lnd_inst%forc_aer_grc(begg:,12), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_dstdry4, output=atm2lnd_inst%forc_aer_grc(begg:,14), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !call state_getimport_1d(importState, Sa_methane, atm2lnd_inst%forc_pch4_grc(begg:), rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! The lilac is sending ndep in units if kgN/m2/s - and ctsm uses units of gN/m2/sec
    ! so the following conversion needs to happen
    ! call state_getimport_1d(importState, Faxa_nhx, output=forc_nhx(begg:), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_getimport_1d(importState, Faxa_nhy, output=forc_nhy(begg:), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! do g = begg,endg
    !    atm2lnd_inst%forc_ndep_grc(g) = (forc_nhx(g) + forc_noy(g))*1000._r8
    ! end do

    !--------------------------
    ! Set force flood back from river to 0
    !--------------------------

    water_inst%wateratm2lndbulk_inst%forc_flood_grc(:) = 0._r8

    do g = begg, endg
       water_inst%wateratm2lndbulk_inst%volr_grc(g) = 0._r8
       water_inst%wateratm2lndbulk_inst%volrmch_grc(g) = 0._r8
    end do

    !--------------------------
    ! Derived quantities for required fields
    ! and corresponding error checks
    !--------------------------

    call derive_quantities(bounds, atm2lnd_inst, water_inst%wateratm2lndbulk_inst, forc_rainc, forc_rainl, forc_snowc, forc_snowl)

    call check_for_errors(bounds, atm2lnd_inst, water_inst%wateratm2lndbulk_inst)

  end subroutine import_fields

  !===============================================================================

  subroutine check_atm_landfrac(importState, bounds, rc)

    ! ------------------------------------------------------------------------
    ! Import Sa_landfrac and check it against CTSM's internal land mask.
    !
    ! We require that CTSM's internal land mask contains all of the atmosphere's land
    ! points (defined as points with landfrac > 0). It is okay for CTSM to include some
    ! points that the atmosphere considers to be ocean, but not the reverse.
    ! ------------------------------------------------------------------------

    ! input/output variables
    type(ESMF_State)                :: importState
    type(bounds_type) , intent(in)  :: bounds
    integer           , intent(out) :: rc

    ! local variables
    real(r8), pointer :: atm_landfrac(:)
    integer :: last_land_index
    integer :: n

    character(len=*), parameter :: subname='(check_atm_landfrac)'
    !---------------------------------------------------------------------------

    ! Implementation notes: The CTSM decomposition is set up so that ocean points appear
    ! at the end of the vectors received from the coupler. Thus, in order to check if
    ! there are any points that the atmosphere considers land but CTSM considers ocean,
    ! it is sufficient to check the points following the typical ending bounds in the
    ! vectors received from the coupler.
    !
    ! Note that we can't use state_getimport here, because that only gets points from
    ! bounds%begg:bounds%endg, whereas we want the points following bounds%endg.

    call state_getfldptr(importState, 'c2l_fb_atm', 'Sa_landfrac', fldptr1d=atm_landfrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    last_land_index = bounds%endg - bounds%begg + 1
    do n = last_land_index + 1, ubound(atm_landfrac, 1)
       if (atm_landfrac(n) > 0._r8) then
          write(iulog,*) 'At point ', n, ' atm landfrac = ', atm_landfrac(n)
          write(iulog,*) 'but CTSM thinks this is ocean.'
          write(iulog,*) "Make sure the mask on CTSM's fatmlndfrc file agrees with the atmosphere's land mask"
          call shr_sys_abort( subname//&
               ' ERROR: atm landfrac > 0 for a point that CTSM thinks is ocean')
       end if
    end do

  end subroutine check_atm_landfrac

  !==============================================================================

  subroutine export_fields(exportState, bounds, rc)

    !-------------------------------
    ! Pack the export state
    !-------------------------------

    ! input/output variables
    type(ESMF_State)                :: exportState
    type(bounds_type) , intent(in)  :: bounds       ! bounds
    integer           , intent(out) :: rc

    ! local variables
    integer                     :: i, g, num
    integer                     :: begg, endg
    real(r8)                    :: array(bounds%begg:bounds%endg)
    character(len=CS)           :: cnum
    character(len=*), parameter :: subname='(lnd_import_export:export_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    begg = bounds%begg
    endg = bounds%endg

    ! -----------------------
    ! output to atm
    ! -----------------------

    call state_setexport_1d(exportState, Sl_t         , lnd2atm_inst%t_rad_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Sl_snowh     , waterlnd2atmbulk_inst%h2osno_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Sl_avsdr     , lnd2atm_inst%albd_grc(begg:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Sl_anidr     , lnd2atm_inst%albd_grc(begg:,2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Sl_avsdf     , lnd2atm_inst%albi_grc(begg:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Sl_anidf     , lnd2atm_inst%albi_grc(begg:,2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Sl_tref      , lnd2atm_inst%t_ref2m_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Sl_qref      , waterlnd2atmbulk_inst%q_ref2m_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Sl_u10       , lnd2atm_inst%u_ref10m_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Fall_taux    , lnd2atm_inst%taux_grc(begg:), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Fall_tauy    , lnd2atm_inst%tauy_grc(begg:), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Fall_lat     , lnd2atm_inst%eflx_lh_tot_grc(begg:), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Fall_sen     , lnd2atm_inst%eflx_sh_tot_grc(begg:), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Fall_lwup    , lnd2atm_inst%eflx_lwrad_out_grc(begg:), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Fall_evap    , waterlnd2atmbulk_inst%qflx_evap_tot_grc(begg:), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Fall_swnet   , lnd2atm_inst%fsa_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Fall_flxdst1 , lnd2atm_inst%flxdst_grc(begg:,1), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Fall_flxdst2 , lnd2atm_inst%flxdst_grc(begg:,2), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Fall_flxdst3 , lnd2atm_inst%flxdst_grc(begg:,3), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Fall_flxdst4 , lnd2atm_inst%flxdst_grc(begg:,4), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Sl_ram1      , lnd2atm_inst%ram1_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Sl_fv        , lnd2atm_inst%fv_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(exportState, Sl_z0m       , lnd2atm_inst%z0m_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! methane
    ! call state_setexport_1d(exportState, Fall_methane , lnd2atm_inst%flux_ch4_grc(begg:), minus=.true., rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! soil water
    ! call state_setexport_1d(exportState, Sl_soilw , water_inst%waterlnd2atmbulk_inst%h2osoi_vol_grc(begg:,1), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! dry dep velocities
    ! do num = 1, drydep_nflds
    !    write(cnum,'(i0)') num
    !    call state_setexport_1d(exportState, trim(Sl_ddvel)//trim(cnum), &
    !         lnd2atm_inst%ddvel_grc(begg:,num), rc=rc)
    ! end do

    ! MEGAN VOC emis fluxes
    ! do num = 1, shr_megan_mechcomps_n
    !    write(cnum,'(i0)') num
    !    call state_setexport_1d(exportState, trim(Fall_voc)//trim(cnum), &
    !         lnd2atm_inst%flxvoc_grc(begg:,num), minus=.true., rc=rc)     
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! end do

    ! fire emis fluxes
    ! do num = 1, emis_nflds
    !    write(cnum,'(i0)') num
    !    call state_setexport_2d(exportState, trim(Fall_fire)//trim(cnum), lnd2atm_inst%fireflx_grc(begg:,num), &
    !         minus = .true., rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! if (emis_nflds > 0) then
    !    call state_setexport_1d(exportState, Sl_fztop, lnd2atm_inst%fireztop_grc(begg:), rc=rc)
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! end if

    ! sign convention is positive downward with hierarchy of atm/glc/lnd/rof/ice/ocn.
    ! i.e. water sent from land to rof is positive

    ! -----------------------
    ! output to river
    ! -----------------------

    ! surface runoff is the sum of qflx_over, qflx_h2osfc_surf
    ! do g = begg,endg
    !    array(g) = water_inst%waterlnd2atmbulk_inst%qflx_rofliq_qsur_grc(g) + &
    !               water_inst%waterlnd2atmbulk_inst%qflx_rofliq_h2osfc_grc(g)
    ! end do

    call state_setexport_1d(exportState, Flrl_rofsur, waterlnd2atmbulk_inst%qflx_rofliq_qsur_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! subsurface runoff is the sum of qflx_drain and qflx_perched_drain
    do g = begg,endg
       array(g) = water_inst%waterlnd2atmbulk_inst%qflx_rofliq_qsub_grc(g) + &
                  water_inst%waterlnd2atmbulk_inst%qflx_rofliq_drain_perched_grc(g)
    end do
    call state_setexport_1d(exportState, Flrl_rofsub, array(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! qgwl sent individually to coupler
    call state_setexport_1d(exportState, Flrl_rofgwl, waterlnd2atmbulk_inst%qflx_rofliq_qgwl_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ice  sent individually to coupler
    call state_setexport_1d(exportState, Flrl_rofi, waterlnd2atmbulk_inst%qflx_rofice_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! irrigation flux to be removed from main channel storage (negative)
    call state_setexport_1d(exportState, Flrl_irrig, waterlnd2atmbulk_inst%qirrig_grc(begg:), &
         minus = .true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine export_fields

  !===============================================================================
  subroutine state_getimport_1d(state, fldname, ctsmdata, rc)

    ! fill in ctsm import data for 1d field

    use ESMF, only : ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT, ESMF_LogFoundError
    use ESMF, only : ESMF_Finalize

    ! input/output variabes
    type(ESMF_State) , intent(in)    :: state
    character(len=*) , intent(in)    :: fldname
    real(r8)         , intent(inout) :: ctsmdata(:)
    integer          , intent(out)   :: rc

    ! local variables
    real(r8), pointer :: fldPtr1d(:)
    integer           :: g
    character(len=*), parameter :: subname='(lnd_import_export:state_getimport_1d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call state_getfldptr(State, trim(fldname), fldptr1d=fldptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = 1,size(ctsmdata)
       ctsmdata(g) = fldptr1d(g)
    end do
    call check_for_nans(ctsmdata, trim(fldname), 1)

  end subroutine state_getimport_1d

  !===============================================================================
  subroutine state_setexport_1d(state, fldname, ctsmdata, minus, rc)

    ! fill in ctsm export data for 1d field

    use ESMF, only : ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT, ESMF_LogFoundError
    use ESMF, only : ESMF_Finalize

    ! input/output variabes
    type(ESMF_State) , intent(in) :: state
    character(len=*) , intent(in) :: fldname
    real(r8)         , intent(in) :: ctsmdata(:)
    logical, optional, intent(in) :: minus
    integer          , intent(out):: rc

    ! local variables
    real(r8), pointer :: fldPtr1d(:)
    integer           :: g
    character(len=*), parameter :: subname='(lnd_export_export:state_setexport_1d)'
    ! ----------------------------------------------

    call state_getfldptr(state, trim(fldname), fldptr1d=fldptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    fldptr1d(:) = 0._r8
    if (present(minus)) then
       do g = 1,size(ctsmdata)
          fldptr1d(g) = -ctsmdata(g)
       end do
    else
       do g = 1,size(ctsmdata)
          fldptr1d(g) = ctsmdata(g)
       end do
    end if
    call check_for_nans(ctsmdata, trim(fldname), 1)

  end subroutine state_setexport_1d

  !===============================================================================
  subroutine state_getfldptr(State, fldname, fldptr1d, fldptr2d, rc)

    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    ! input/output variables
    type(ESMF_State),             intent(in)    :: State
    character(len=*),             intent(in)    :: fldname
    real(R8), pointer, optional , intent(out)   :: fldptr1d(:)
    real(R8), pointer, optional , intent(out)   :: fldptr2d(:,:)
    integer,                      intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    character(len=*), parameter :: subname='(lnd_import_export:state_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (present(fldptr1d)) then
       call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (present(fldptr2d)) then
       call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort("either fldptr1d or fldptr2d must be an input argument")
    end if

  end subroutine state_getfldptr

  !===============================================================================
  logical function fldchk(state, fldname)
    ! ----------------------------------------------
    ! Determine if field with fldname is in the input state
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State), intent(in)  :: state
    character(len=*), intent(in)  :: fldname

    ! local variables
    type(ESMF_StateItem_Flag)   :: itemFlag
    ! ----------------------------------------------
    call ESMF_StateGet(state, trim(fldname), itemFlag)
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then
       fldchk = .true.
    else
       fldchk = .false.
    endif
  end function fldchk

end module lnd_import_export
