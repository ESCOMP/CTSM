module lnd_import_export

  use ESMF
  use shr_kind_mod          , only : r8 => shr_kind_r8, cx=>shr_kind_cx, cxx=>shr_kind_cxx, cs=>shr_kind_cs
  use shr_infnan_mod        , only : isnan => shr_infnan_isnan
  use shr_string_mod        , only : shr_string_listGetName, shr_string_listGetNum
  use shr_sys_mod           , only : shr_sys_abort
  use shr_const_mod         , only : SHR_CONST_TKFRZ, fillvalue=>SHR_CONST_SPVAL
  use clm_varctl            , only : iulog, co2_ppmv, ndep_from_cpl
  use clm_varcon            , only : rair, o2_molar_const
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
  use QSatMod               , only : QSat

  implicit none
  private ! except

  public  :: import_fields
  public  :: export_fields

  private :: state_getimport
  private :: state_setexport
  private :: state_getfldptr
  private :: check_for_nans

  ! from atm->lnd
  integer                :: ndep_nflds               ! number  of nitrogen deposition fields from atm->lnd/ocn

  ! from lnd->atm
  integer                :: drydep_nflds             ! number of dry deposition velocity fields lnd-> atm
  integer                :: emis_nflds               ! number of fire emission fields from lnd-> atm

  integer                :: glc_nec = 10             ! number of glc elevation classes
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
    real(r8)                  :: dum1, dum2, dum3                       ! arguments returned by subr. QSat but not used here
    real(r8)                  :: qsat_kg_kg                             ! saturation specific humidity (kg/kg)
    real(r8)                  :: co2_ppmv_val                           ! temporary
    real(r8)                  :: forc_t                                 ! atmospheric temperature (Kelvin)
    real(r8)                  :: forc_q                                 ! atmospheric specific humidity (kg/kg)
    real(r8)                  :: forc_pbot                              ! atmospheric pressure (Pa)
    real(r8)                  :: forc_rainc(bounds%begg:bounds%endg)    ! rainxy Atm flux mm/s
    real(r8)                  :: forc_rainl(bounds%begg:bounds%endg)    ! rainxy Atm flux mm/s
    real(r8)                  :: forc_snowc(bounds%begg:bounds%endg)    ! snowfxy Atm flux  mm/s
    real(r8)                  :: forc_snowl(bounds%begg:bounds%endg)    ! snowfxl Atm flux  mm/s
    real(r8)                  :: forc_noy(bounds%begg:bounds%endg)
    real(r8)                  :: forc_nhx(bounds%begg:bounds%endg)
    real(r8)                  :: topo_grc(bounds%begg:bounds%endg, 0:glc_nec)
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

    call state_getimport(importState, 'c2l_fb_atm', 'Sa_z', bounds, &
         output=atm2lnd_inst%forc_hgt_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Sa_topo', bounds, &
         output=atm2lnd_inst%forc_topo_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Sa_u', bounds, &
         output=atm2lnd_inst%forc_u_grc, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Sa_v', bounds, &
         output=atm2lnd_inst%forc_v_grc, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Sa_ptem', bounds, &
         output=atm2lnd_inst%forc_th_not_downscaled_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Sa_shum', bounds, &
         output=water_inst%wateratm2lndbulk_inst%forc_q_not_downscaled_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Sa_pbot', bounds, &
         output=atm2lnd_inst%forc_pbot_not_downscaled_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Sa_tbot', bounds, &
         output=atm2lnd_inst%forc_t_not_downscaled_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_rainc', bounds, &
         output=forc_rainc, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_rainl', bounds, &
         output=forc_rainl, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_snowc', bounds, &
         output=forc_snowc, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_snowl', bounds, &
         output=forc_snowl, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_lwdn', bounds, &
         output=atm2lnd_inst%forc_lwrad_not_downscaled_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_swvdr', bounds, &
         output=atm2lnd_inst%forc_solad_grc(:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_swndr', bounds, &
         output=atm2lnd_inst%forc_solad_grc(:,2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_swvdf', bounds, &
         output=atm2lnd_inst%forc_solai_grc(:,1), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_swndf', bounds, &
         output=atm2lnd_inst%forc_solai_grc(:,2), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ! Atmosphere prognostic/prescribed aerosol fields

    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_bcphidry', bounds, &
         output=atm2lnd_inst%forc_aer_grc(:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_bcphodry', bounds, &
         output=atm2lnd_inst%forc_aer_grc(:,2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_bcphiwet', bounds, &
         output=atm2lnd_inst%forc_aer_grc(:,3), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_ocphidry', bounds, &
         output=atm2lnd_inst%forc_aer_grc(:,4), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_ocphodry', bounds, &
         output=atm2lnd_inst%forc_aer_grc(:,5), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_ocphiwet', bounds, &
         output=atm2lnd_inst%forc_aer_grc(:,6), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_dstwet1', bounds, &
         output=atm2lnd_inst%forc_aer_grc(:,7),  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_dstdry1', bounds, &
         output=atm2lnd_inst%forc_aer_grc(:,8),  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_dstwet2', bounds, &
         output=atm2lnd_inst%forc_aer_grc(:,9),  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_dstdry2', bounds, &
         output=atm2lnd_inst%forc_aer_grc(:,10), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_dstwet3', bounds, &
         output=atm2lnd_inst%forc_aer_grc(:,11), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_dstdry3', bounds, &
         output=atm2lnd_inst%forc_aer_grc(:,12), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_dstwet4', bounds, &
         output=atm2lnd_inst%forc_aer_grc(:,13), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'c2l_fb_atm', 'Faxa_dstdry4', bounds, &
         output=atm2lnd_inst%forc_aer_grc(:,14), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! call state_getimport(importState, 'c2l_fb_atm', 'Sa_methane', bounds, output=atm2lnd_inst%forc_pch4_grc, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! The lilac is sending ndep in units if kgN/m2/s - and ctsm uses units of gN/m2/sec
    ! so the following conversion needs to happen
    ! call state_getimport(importState, 'c2l_fb_atm', 'Faxa_nhx', bounds, output=forc_nhx, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_getimport(importState, 'c2l_fb_atm', 'Faxa_noy', bounds, output=forc_noy, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! do g = begg,endg
    !    atm2lnd_inst%forc_ndep_grc(g) = (forc_nhx(g) + forc_noy(g))*1000._r8
    ! end do

    !--------------------------
    ! Set force flood back from river to 0
    !--------------------------

    water_inst%wateratm2lndbulk_inst%forc_flood_grc(:) = 0._r8

    !--------------------------
    ! Derived quantities
    !--------------------------

    do g = begg, endg
       forc_t    = atm2lnd_inst%forc_t_not_downscaled_grc(g)
       forc_q    = water_inst%wateratm2lndbulk_inst%forc_q_not_downscaled_grc(g)
       forc_pbot = atm2lnd_inst%forc_pbot_not_downscaled_grc(g)

       atm2lnd_inst%forc_hgt_u_grc(g) = atm2lnd_inst%forc_hgt_grc(g)    !observational height of wind [m]
       atm2lnd_inst%forc_hgt_t_grc(g) = atm2lnd_inst%forc_hgt_grc(g)    !observational height of temperature [m]
       atm2lnd_inst%forc_hgt_q_grc(g) = atm2lnd_inst%forc_hgt_grc(g)    !observational height of humidity [m]

       atm2lnd_inst%forc_vp_grc(g) = forc_q * forc_pbot  / (0.622_r8 + 0.378_r8 * forc_q)

       atm2lnd_inst%forc_rho_not_downscaled_grc(g) = &
            (forc_pbot - 0.378_r8 * atm2lnd_inst%forc_vp_grc(g)) / (rair * forc_t)

       atm2lnd_inst%forc_po2_grc(g) = o2_molar_const * forc_pbot

       atm2lnd_inst%forc_pco2_grc(g) = co2_ppmv * 1.e-6_r8 * forc_pbot

       atm2lnd_inst%forc_wind_grc(g) = sqrt(atm2lnd_inst%forc_u_grc(g)**2 + atm2lnd_inst%forc_v_grc(g)**2)

       atm2lnd_inst%forc_solar_grc(g) = atm2lnd_inst%forc_solad_grc(g,1) + atm2lnd_inst%forc_solai_grc(g,1) + &
            atm2lnd_inst%forc_solad_grc(g,2) + atm2lnd_inst%forc_solai_grc(g,2)

       water_inst%wateratm2lndbulk_inst%forc_rain_not_downscaled_grc(g)  = forc_rainc(g) + forc_rainl(g)
       water_inst%wateratm2lndbulk_inst%forc_snow_not_downscaled_grc(g)  = forc_snowc(g) + forc_snowl(g)


       call QSat(forc_t, forc_pbot, dum1, dum2, qsat_kg_kg, dum3)

       ! modify specific humidity if precip occurs
       if (1==2) then
          if ((forc_rainc(g)+forc_rainl(g)) > 0._r8) then
             forc_q = 0.95_r8*qsat_kg_kg
             !forc_q = qsat_kg_kg
             water_inst%wateratm2lndbulk_inst%forc_q_not_downscaled_grc(g) = forc_q
          endif
       endif

       water_inst%wateratm2lndbulk_inst%forc_rh_grc(g) = 100.0_r8*(forc_q / qsat_kg_kg)
       water_inst%wateratm2lndbulk_inst%volr_grc(g) = 0._r8
       water_inst%wateratm2lndbulk_inst%volrmch_grc(g) = 0._r8
    end do

    !--------------------------
    ! Error checks
    !--------------------------

    ! Check that solar, specific-humidity and LW downward aren't negative
    do g = begg,endg
       if ( atm2lnd_inst%forc_lwrad_not_downscaled_grc(g) <= 0.0_r8 ) then
          call shr_sys_abort( subname//&
               ' ERROR: Longwave down sent from the atmosphere model is negative or zero' )
       end if
       if ( (atm2lnd_inst%forc_solad_grc(g,1) < 0.0_r8) .or. &
            (atm2lnd_inst%forc_solad_grc(g,2) < 0.0_r8) .or. &
            (atm2lnd_inst%forc_solai_grc(g,1) < 0.0_r8) .or. &
            (atm2lnd_inst%forc_solai_grc(g,2) < 0.0_r8) ) then
          call shr_sys_abort( subname//&
               ' ERROR: One of the solar fields (indirect/diffuse, vis or near-IR)'// &
               ' from the atmosphere model is negative or zero' )
       end if
       if ( water_inst%wateratm2lndbulk_inst%forc_q_not_downscaled_grc(g) < 0.0_r8 )then
          call shr_sys_abort( subname//&
               ' ERROR: Bottom layer specific humidty sent from the atmosphere model is less than zero' )
       end if
    end do

    ! Make sure relative humidity is properly bounded
    ! atm2lnd_inst%forc_rh_grc(g) = min( 100.0_r8, atm2lnd_inst%forc_rh_grc(g) )
    ! atm2lnd_inst%forc_rh_grc(g) = max(   0.0_r8, atm2lnd_inst%forc_rh_grc(g) )

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
    real(r8)                    :: array(bounds%begg:bounds%endg)
    character(len=*), parameter :: subname='(lnd_import_export:export_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! -----------------------
    ! output to atm
    ! -----------------------

    call state_setexport(exportState, 'l2c_fb_atm', 'Sl_t', bounds, &
         input=lnd2atm_inst%t_rad_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Sl_snowh', bounds, &
         input=water_inst%waterlnd2atmbulk_inst%h2osno_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Sl_avsdr', bounds, &
         input=lnd2atm_inst%albd_grc(bounds%begg:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Sl_anidr', bounds, &
         input=lnd2atm_inst%albd_grc(bounds%begg:,2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Sl_avsdf', bounds, &
         input=lnd2atm_inst%albi_grc(bounds%begg:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Sl_anidf', bounds, &
         input=lnd2atm_inst%albi_grc(bounds%begg:,2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Sl_tref', bounds, &
         input=lnd2atm_inst%t_ref2m_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Sl_qref', bounds, &
         input=water_inst%waterlnd2atmbulk_inst%q_ref2m_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Sl_u10', bounds, &
         input=lnd2atm_inst%u_ref10m_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Fall_taux', bounds, &
         input=lnd2atm_inst%taux_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Fall_tauy', bounds, &
         input=lnd2atm_inst%tauy_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Fall_lat', bounds, &
         input=lnd2atm_inst%eflx_lh_tot_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Fall_sen', bounds, &
         input=lnd2atm_inst%eflx_sh_tot_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Fall_lwup', bounds, &
         input=lnd2atm_inst%eflx_lwrad_out_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Fall_evap', bounds, &
         input=water_inst%waterlnd2atmbulk_inst%qflx_evap_tot_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Fall_swnet', bounds, &
         input=lnd2atm_inst%fsa_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Fall_flxdst1', bounds, &
         input=lnd2atm_inst%flxdst_grc(:,1), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'l2c_fb_atm', 'Fall_flxdst2', bounds, &
         input=lnd2atm_inst%flxdst_grc(:,2), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'l2c_fb_atm', 'Fall_flxdst3', bounds, &
         input=lnd2atm_inst%flxdst_grc(:,3), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'l2c_fb_atm', 'Fall_flxdst4', bounds, &
         input=lnd2atm_inst%flxdst_grc(:,4), minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Sl_ram1', bounds, &
         input=lnd2atm_inst%ram1_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Sl_fv', bounds, &
         input=lnd2atm_inst%fv_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'l2c_fb_atm', 'Sl_z0m', bounds, &
         input=lnd2atm_inst%z0m_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! methanem
    ! call state_setexport(exportState, 'l2c_fb_atm', 'Fall_methane', bounds, &
    !    input=lnd2atm_inst%flux_ch4_grc, minus=.true., rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! soil water
    ! call state_setexport(exportState, 'l2c_fb_atm', 'Sl_soilw', bounds, &
    !    input=water_inst%waterlnd2atmbulk_inst%h2osoi_vol_grc(:,1), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! dry dep velocities
    ! do num = 1, drydep_nflds
    !    call state_setexport(exportState, 'l2c_fb_atm', 'Sl_ddvel', bounds, &
    !       input=lnd2atm_inst%ddvel_grc(:,num), ungridded_index=num, rc=rc)
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! end do

    ! MEGAN VOC emis fluxes
    ! do num = 1, shr_megan_mechcomps_n
    !    call state_setexport(exportState, 'l2c_fb_atm', 'Fall_voc', bounds, &
    !       input=lnd2atm_inst%flxvoc_grc(:,num), minus=.true., ungridded_index=num, rc=rc)
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! end do

    ! fire emis fluxes
    ! do num = 1, emis_nflds
    !    call state_setexport(exportState, 'l2c_fb_atm', 'Fall_fire', bounds, &
    !       input=lnd2atm_inst%fireflx_grc(:,num), minus=.true., ungridded_index=num, rc=rc)
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! end do
    ! if (emis_nflds > 0) then
    !    call state_setexport(exportState, 'l2c_fb_atm', 'Sl_fztopo', bounds, input=lnd2atm_inst%fireztop_grc, rc=rc)
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! endif
    ! sign convention is positive downward with hierarchy of atm/glc/lnd/rof/ice/ocn.
    ! i.e. water sent from land to rof is positive

    ! -----------------------
    ! output to river
    ! -----------------------

    ! surface runoff is the sum of qflx_over, qflx_h2osfc_surf
    ! do g = bounds%begg,bounds%endg
    !    array(g) = water_inst%waterlnd2atmbulk_inst%qflx_rofliq_qsur_grc(g) + &
    !               water_inst%waterlnd2atmbulk_inst%qflx_rofliq_h2osfc_grc(g)
    ! end do

    call state_setexport(exportState, 'l2c_fb_rof', 'Flrl_rofsur', bounds, &
         input=water_inst%waterlnd2atmbulk_inst%qflx_rofliq_qsur_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! subsurface runoff is the sum of qflx_drain and qflx_perched_drain
    do g = bounds%begg,bounds%endg
       array(g) = water_inst%waterlnd2atmbulk_inst%qflx_rofliq_qsub_grc(g) + &
                  water_inst%waterlnd2atmbulk_inst%qflx_rofliq_drain_perched_grc(g)
    end do
    call state_setexport(exportState, 'l2c_fb_rof', 'Flrl_rofsub', bounds, &
         input=array, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! qgwl sent individually to coupler
    call state_setexport(exportState, 'l2c_fb_rof', 'Flrl_rofgwl', bounds, &
         input=water_inst%waterlnd2atmbulk_inst%qflx_rofliq_qgwl_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ice  sent individually to coupler
    call state_setexport(exportState, 'l2c_fb_rof', 'Flrl_rofi', bounds, &
         input=water_inst%waterlnd2atmbulk_inst%qflx_rofice_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! irrigation flux to be removed from main channel storage (negative)
    call state_setexport(exportState, 'l2c_fb_rof', 'Flrl_irrig', bounds, &
         input=water_inst%waterlnd2atmbulk_inst%qirrig_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine export_fields

  !===============================================================================

  subroutine state_getimport(state, fb, fldname, bounds, output, ungridded_index, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(in)    :: state
    character(len=*)    , intent(in)    :: fb
    character(len=*)    , intent(in)    :: fldname
    type(bounds_type)   , intent(in)    :: bounds
    real(r8)            , intent(out)   :: output(bounds%begg:bounds%endg)
    integer, optional   , intent(in)    :: ungridded_index
    integer             , intent(out)   :: rc

    ! local variables
    integer                   :: g, i,n
    real(R8), pointer         :: fldptr1d(:)
    real(R8), pointer         :: fldptr2d(:,:)
    character(len=cs)         :: cvalue
    character(len=*), parameter :: subname='(lnd_import_export:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    if (masterproc .and. debug > 0)  then
       write(iulog,F01)' Show me what is in the state? for  '//trim(fldname)
       call ESMF_StatePrint(state, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Get the pointer to data in the field
    if (present(ungridded_index)) then
       write(cvalue,*) ungridded_index
       call ESMF_LogWrite(trim(subname)//": getting import for "//trim(fldname)//" index "//trim(cvalue), &
            ESMF_LOGMSG_INFO)
       call state_getfldptr(state, trim(fb), trim(fldname), fldptr2d=fldptr2d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_LogWrite(trim(subname)//": getting import for "//trim(fldname),ESMF_LOGMSG_INFO)
       call state_getfldptr(state, trim(fb), trim(fldname), fldptr1d=fldptr1d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Fill in output array
    if (present(ungridded_index)) then
       do g = bounds%begg, bounds%endg
          n = g - bounds%begg + 1
          output(g) = fldptr2d(ungridded_index,n)
       end do
    else
       do g = bounds%begg, bounds%endg
          n = g - bounds%begg + 1
          output(g) = fldptr1d(n)
          if (masterproc .and. debug > 0 .and. get_nstep() < 5) then
             write(iulog,F02)' n, g , fldptr1d(n) '//trim(fldname)//' = ',n, g, fldptr1d(n)
          end if
       end do
    end if

    ! Write debug output if appropriate
    if (masterproc .and. debug > 0 .and. get_nstep() < 5) then
       do g = bounds%begg,bounds%endg
          i = 1 + g - bounds%begg
          write(iulog,F02)'import: nstep, n, '//trim(fldname)//' = ',get_nstep(),i,output(g)
       end do
    end if

    ! Check for nans
    call check_for_nans(output, trim(fldname), bounds%begg)

  end subroutine state_getimport

  !===============================================================================

  subroutine state_setexport(state, fb, fldname, bounds, input, minus, ungridded_index, rc)

    ! ----------------------------------------------
    ! Map input array to export state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    character(len=*)    , intent(in)    :: fb
    type(bounds_type)   , intent(in)    :: bounds
    character(len=*)    , intent(in)    :: fldname
    real(r8)            , intent(in)    :: input(bounds%begg:bounds%endg)
    logical, optional   , intent(in)    :: minus
    integer, optional   , intent(in)    :: ungridded_index
    integer             , intent(out)   :: rc

    ! local variables
    logical                     :: l_minus  ! local version of minus
    integer                     :: g, i, n
    real(R8), pointer           :: fldptr1d(:)
    real(R8), pointer           :: fldptr2d(:,:)
    character(len=cs)           :: cvalue
    character(len=*), parameter :: subname='(lnd_import_export:state_setexport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    l_minus = .false.
    if (present(minus)) then
       l_minus = minus
    end if

    ! get field pointer
    if (present(ungridded_index)) then
       call ESMF_LogWrite(trim(subname)//": setting export for "//trim(fldname)//" index "//trim(cvalue), &
            ESMF_LOGMSG_INFO)
       call state_getfldptr(state, trim(fb), trim(fldname), fldptr2d=fldptr2d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_LogWrite(trim(subname)//": setting export for "//trim(fldname), ESMF_LOGMSG_INFO)
       call state_getfldptr(state, trim(fb), trim(fldname), fldptr1d=fldptr1d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! determine output array
    if (present(ungridded_index)) then
       fldptr2d(ungridded_index,:) = fillvalue
       do g = bounds%begg, bounds%endg
          n = g - bounds%begg + 1
          if (l_minus) then
             fldptr2d(ungridded_index,n) = -input(g)
          else
             fldptr2d(ungridded_index,n) = input(g)
          end if
       end do
    else
       fldptr1d(:) = fillvalue
       do g = bounds%begg, bounds%endg
          n = g - bounds%begg + 1
          if (l_minus) then
             fldptr1d(n) = -input(g)
          else
             fldptr1d(n) = input(g)
          end if
       end do
    end if

    ! write debug output if appropriate
    if (masterproc .and. debug > 0 .and. get_nstep() < 5) then
       do g = bounds%begg,bounds%endg
          i = 1 + g - bounds%begg
          write(iulog,F01)'export: nstep, n, '//trim(fldname)//' = ',get_nstep(),i,input(g)
       end do
    end if

    ! check for nans
    call check_for_nans(input, trim(fldname), bounds%begg)

  end subroutine state_setexport

  !===============================================================================

  subroutine state_getfldptr(State, fb, fldname, fldptr1d, fldptr2d, rc)

    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State),             intent(in)    :: State
    character(len=*),             intent(in)    :: fb 
    character(len=*),             intent(in)    :: fldname
    real(R8), pointer, optional , intent(out)   :: fldptr1d(:)
    real(R8), pointer, optional , intent(out)   :: fldptr2d(:,:)
    integer,                      intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    type(ESMF_Mesh)             :: lmesh
    integer                     :: nnodes, nelements
    type(ESMF_FieldBundle)      :: fieldBundle
    character(len=*), parameter :: subname='(lnd_import_export:state_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Get the fieldbundle from the state...
    call ESMF_StateGet(state, trim(fb), fieldBundle, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("ERROR: fb "//trim(fb)//" not found in state")

    ! Get the field from the field bundle
    call ESMF_FieldBundleGet(fieldBundle,fieldName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the status of the field
    call ESMF_FieldGet(lfield, status=status, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then
       call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    else
       call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (nnodes == 0 .and. nelements == 0) then
          call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if

       ! Get the data from the field
       if (present(fldptr1d)) then
          call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (masterproc .and. debug > 0)  then
             write(iulog,F01)' in '//trim(subname)//'fldptr1d for '//trim(fldname)//' is  '
          end if
       else if (present(fldptr2d)) then
          call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call shr_sys_abort("either fldptr1d or fldptr2d must be an input argument")
       end if
    endif  ! status

  end subroutine state_getfldptr

  !===============================================================================

  subroutine check_for_nans(array, fname, begg)

    ! input/output variables
    real(r8)         , intent(in) :: array(:)
    character(len=*) , intent(in) :: fname
    integer          , intent(in) :: begg
    !
    ! local variables
    integer :: i
    !-------------------------------------------------------------------------------

    ! Check if any input from lilac or output to lilac is NaN

    if (any(isnan(array))) then
       write(iulog,*) '# of NaNs = ', count(isnan(array))
       write(iulog,*) 'Which are NaNs = ', isnan(array)
       do i = 1, size(array)
          if (isnan(array(i))) then
             write(iulog,*) "NaN found in field ", trim(fname), ' at gridcell index ',begg+i-1
          end if
       end do
       call shr_sys_abort(' ERROR: One or more of the output from CLM to the coupler are NaN ' )
    end if
  end subroutine check_for_nans

end module lnd_import_export
