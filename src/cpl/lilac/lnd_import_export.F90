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
  use spmdMod               , only : masterproc
  use decompmod             , only : bounds_type
  use lnd2atmType           , only : lnd2atm_type
  use lnd2glcMod            , only : lnd2glc_type
  use atm2lndType           , only : atm2lnd_type
  use domainMod             , only : ldomain
  use shr_megan_mod         , only : shr_megan_mechcomps_n  ! TODO: need to add a namelist read nere
  use lnd_shr_methods       , only : chkerr
  use clm_instMod           , only : atm2lnd_inst, lnd2atm_inst, water_inst

  implicit none
  private ! except

  public  :: import_fields
  public  :: export_fields

  private :: fldlist_add
  private :: state_getimport
  private :: state_setexport
  private :: state_getfldptr
  private :: check_for_nans

  type fld_list_type
     character(len=128) :: stdname
     integer :: ungridded_lbound = 0
     integer :: ungridded_ubound = 0
  end type fld_list_type

  integer, parameter     :: fldsMax = 100
  integer                :: fldsToLnd_num = 0
  integer                :: fldsFrLnd_num = 0
  type (fld_list_type)   :: fldsToLnd(fldsMax)
  type (fld_list_type)   :: fldsFrLnd(fldsMax)
  integer, parameter     :: gridTofieldMap = 2       ! ungridded dimension is innermost

  logical                :: glc_present    = .false. ! .true. => running with a non-stubGLC model
  logical                :: rof_prognostic = .false. ! .true. => running with a prognostic ROF model

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

  subroutine import_fields( gcomp, bounds, rc)

    !---------------------------------------------------------------------------
    ! Convert the input data from the mediator to the land model
    !---------------------------------------------------------------------------

    ! input/output variabes
    type(ESMF_GridComp)             :: gcomp
    type(bounds_type) , intent(in)  :: bounds       ! bounds
    integer           , intent(out) :: rc

    ! local variables
    type(ESMF_State)          :: importState
    integer                   :: num
    integer                   :: begg, endg                             ! bounds
    integer                   :: g,i,k                                  ! indices
    real(r8)                  :: e                                      ! vapor pressure (Pa)
    real(r8)                  :: qsat                                   ! saturation specific humidity (kg/kg)
    real(r8)                  :: co2_ppmv_val                           ! temporary
    real(r8)                  :: esatw                                  ! saturation vapor pressure over water (Pa)
    real(r8)                  :: esati                                  ! saturation vapor pressure over ice (Pa)
    real(r8)                  :: a0,a1,a2,a3,a4,a5,a6                   ! coefficients for esat over water
    real(r8)                  :: b0,b1,b2,b3,b4,b5,b6                   ! coefficients for esat over ice
    real(r8)                  :: tdc, t                                 ! Kelvins to Celcius function and its input
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

    ! Constants to compute vapor pressure
    parameter (a0=6.107799961_r8    , a1=4.436518521e-01_r8, &
         a2=1.428945805e-02_r8, a3=2.650648471e-04_r8, &
         a4=3.031240396e-06_r8, a5=2.034080948e-08_r8, &
         a6=6.136820929e-11_r8)

    parameter (b0=6.109177956_r8    , b1=5.034698970e-01_r8, &
         b2=1.886013408e-02_r8, b3=4.176223716e-04_r8, &
         b4=5.824720280e-06_r8, b5=4.838803174e-08_r8, &
         b6=1.838826904e-10_r8)

    ! function declarations
    tdc(t) = min( 50._r8, max(-50._r8,(t-SHR_CONST_TKFRZ)) )
    esatw(t) = 100._r8*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
    esati(t) = 100._r8*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get import state
    call ESMF_GridCompGet(gcomp, importState=importState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set bounds
    begg = bounds%begg; endg=bounds%endg

    ! Note: precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.

    !--------------------------
    ! Required atmosphere input fields
    !--------------------------

    call state_getimport(importState, 'Sa_z', bounds, output=atm2lnd_inst%forc_hgt_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_topo', bounds, output=atm2lnd_inst%forc_topo_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_u', bounds, output=atm2lnd_inst%forc_u_grc, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_v', bounds, output=atm2lnd_inst%forc_v_grc, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_ptem', bounds, output=atm2lnd_inst%forc_th_not_downscaled_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_shum', bounds, output=water_inst%wateratm2lndbulk_inst%forc_q_not_downscaled_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_pbot', bounds, output=atm2lnd_inst%forc_pbot_not_downscaled_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Sa_tbot', bounds, output=atm2lnd_inst%forc_t_not_downscaled_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_rainc', bounds, output=forc_rainc, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_rainl', bounds, output=forc_rainl, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_snowc', bounds, output=forc_snowc, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_snowl', bounds, output=forc_snowl, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_lwdn', bounds, output=atm2lnd_inst%forc_lwrad_not_downscaled_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_swvdr', bounds, output=atm2lnd_inst%forc_solad_grc(:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_swndr', bounds, output=atm2lnd_inst%forc_solad_grc(:,2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_swvdf', bounds, output=atm2lnd_inst%forc_solai_grc(:,1), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_swndf', bounds, output=atm2lnd_inst%forc_solai_grc(:,2), rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ! Atmosphere prognostic/prescribed aerosol fields

    ! ! bcphidry
    ! call state_getimport(importState, 'Faxa_bcph', bounds, output=atm2lnd_inst%forc_aer_grc(:,1), &
    !      ungridded_index=1, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! ! bcphodry
    ! call state_getimport(importState, 'Faxa_bcph', bounds, output=atm2lnd_inst%forc_aer_grc(:,2), &
    !      ungridded_index=2, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! ! bcphiwet
    ! call state_getimport(importState, 'Faxa_bcph', bounds, output=atm2lnd_inst%forc_aer_grc(:,3), &
    !      ungridded_index=3, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ! ocphidry
    ! call state_getimport(importState, 'Faxa_ocph', bounds, output=atm2lnd_inst%forc_aer_grc(:,4), &
    !      ungridded_index=1, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! ! bcphodry
    ! call state_getimport(importState, 'Faxa_ocph', bounds, output=atm2lnd_inst%forc_aer_grc(:,5), &
    !      ungridded_index=2, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! ! bcphiwet
    ! call state_getimport(importState, 'Faxa_ocph', bounds, output=atm2lnd_inst%forc_aer_grc(:,6), &
    !      ungridded_index=3, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! call state_getimport(importState, 'Faxa_dstwet', bounds, output=atm2lnd_inst%forc_aer_grc(:,7), &
    !      ungridded_index=1, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_getimport(importState, 'Faxa_dstdry', bounds, output=atm2lnd_inst%forc_aer_grc(:,8), &
    !      ungridded_index=1, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! call state_getimport(importState, 'Faxa_dstwet', bounds, output=atm2lnd_inst%forc_aer_grc(:,9), &
    !      ungridded_index=2, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_getimport(importState, 'Faxa_dstdry', bounds, output=atm2lnd_inst%forc_aer_grc(:,10), &
    !      ungridded_index=2, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! call state_getimport(importState, 'Faxa_dstwet', bounds, output=atm2lnd_inst%forc_aer_grc(:,11), &
    !      ungridded_index=3, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_getimport(importState, 'Faxa_dstdry', bounds, output=atm2lnd_inst%forc_aer_grc(:,12), &
    !      ungridded_index=3, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! call state_getimport(importState, 'Faxa_dstwet', bounds, output=atm2lnd_inst%forc_aer_grc(:,13), &
    !      ungridded_index=4, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_getimport(importState, 'Faxa_dstdry', bounds, output=atm2lnd_inst%forc_aer_grc(:,14), &
    !      ungridded_index=4, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! call state_getimport(importState, 'Sa_methane', bounds, output=atm2lnd_inst%forc_pch4_grc, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ! The mediator is sending ndep in units if kgN/m2/s - and ctsm uses units of gN/m2/sec
    ! ! so the following conversion needs to happen

    ! call state_getimport(importState, 'Faxa_nhx', bounds, output=forc_nhx, ungridded_index=1, rc=rc )
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_getimport(importState, 'Faxa_noy', bounds, output=forc_noy, ungridded_index=2, rc=rc )
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


       if (forc_t > SHR_CONST_TKFRZ) then
          e = esatw(tdc(forc_t))
       else
          e = esati(tdc(forc_t))
       end if
       qsat = 0.622_r8*e / (forc_pbot - 0.378_r8*e)

       ! modify specific humidity if precip occurs
       if (1==2) then
          if ((forc_rainc(g)+forc_rainl(g)) > 0._r8) then
             forc_q = 0.95_r8*qsat
             !forc_q = qsat
             water_inst%wateratm2lndbulk_inst%forc_q_not_downscaled_grc(g) = forc_q
          endif
       endif

       water_inst%wateratm2lndbulk_inst%forc_rh_grc(g) = 100.0_r8*(forc_q / qsat)
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

  !==============================================================================

  subroutine export_fields(gcomp, bounds, rc)

    !-------------------------------
    ! Pack the export state
    !-------------------------------

    ! input/output variables
    type(ESMF_GridComp)             :: gcomp
    type(bounds_type) , intent(in)  :: bounds       ! bounds
    integer           , intent(out) :: rc

    ! local variables
    type(ESMF_State)            :: exportState
    integer                     :: i, g, num
    real(r8)                    :: array(bounds%begg:bounds%endg)
    character(len=*), parameter :: subname='(lnd_import_export:export_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get export state (ESMF)
    call ESMF_GridCompGet(gcomp, exportState=exportState, rc=rc) ! do we need the clock now?
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------
    ! output to mediator
    ! -----------------------

    call state_setexport(exportState, 'Sl_lfrin', bounds, input=ldomain%frac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------
    ! output to atm
    ! -----------------------

    call state_setexport(exportState, 'Sl_t', bounds, input=lnd2atm_inst%t_rad_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_snowh', bounds, input=water_inst%waterlnd2atmbulk_inst%h2osno_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_avsdr', bounds, input=lnd2atm_inst%albd_grc(bounds%begg:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_anidr', bounds, input=lnd2atm_inst%albd_grc(bounds%begg:,2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_avsdf', bounds, input=lnd2atm_inst%albi_grc(bounds%begg:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_anidf', bounds, input=lnd2atm_inst%albi_grc(bounds%begg:,2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_tref', bounds, input=lnd2atm_inst%t_ref2m_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_qref', bounds, input=water_inst%waterlnd2atmbulk_inst%q_ref2m_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_u10', bounds, input=lnd2atm_inst%u_ref10m_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_taux', bounds, input=lnd2atm_inst%taux_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_tauy', bounds, input=lnd2atm_inst%tauy_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_lat', bounds, input=lnd2atm_inst%eflx_lh_tot_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_sen', bounds, input=lnd2atm_inst%eflx_sh_tot_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_lwup', bounds, input=lnd2atm_inst%eflx_lwrad_out_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_evap', bounds, input=water_inst%waterlnd2atmbulk_inst%qflx_evap_tot_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_swnet', bounds, input=lnd2atm_inst%fsa_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_flxdst', bounds, input=lnd2atm_inst%flxdst_grc(:,1), &
         minus=.true., ungridded_index=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'Fall_flxdst', bounds, input=lnd2atm_inst%flxdst_grc(:,2), &
         minus=.true., ungridded_index=2, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'Fall_flxdst', bounds, input=lnd2atm_inst%flxdst_grc(:,3), &
         minus=.true., ungridded_index=3, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport(exportState, 'Fall_flxdst', bounds, input=lnd2atm_inst%flxdst_grc(:,4), &
         minus=.true., ungridded_index=4, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Fall_methane', bounds, input=lnd2atm_inst%flux_ch4_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_ram1', bounds, input=lnd2atm_inst%ram1_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_fv', bounds, input=lnd2atm_inst%fv_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Sl_soilw', bounds, input=water_inst%waterlnd2atmbulk_inst%h2osoi_vol_grc(:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! dry dep velocities
    do num = 1, drydep_nflds
       call state_setexport(exportState, 'Sl_ddvel', bounds, input=lnd2atm_inst%ddvel_grc(:,num), &
            ungridded_index=num, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    ! MEGAN VOC emis fluxes
    do num = 1, shr_megan_mechcomps_n
       call state_setexport(exportState, 'Fall_voc', bounds, input=lnd2atm_inst%flxvoc_grc(:,num), minus=.true., &
            ungridded_index=num, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    ! fire emis fluxes
    do num = 1, emis_nflds
       call state_setexport(exportState, 'Fall_fire', bounds, input=lnd2atm_inst%fireflx_grc(:,num), minus=.true., &
            ungridded_index=num, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do
    if (emis_nflds > 0) then
       call state_setexport(exportState, 'Sl_fztopo', bounds, input=lnd2atm_inst%fireztop_grc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
    ! sign convention is positive downward with hierarchy of atm/glc/lnd/rof/ice/ocn.
    ! i.e. water sent from land to rof is positive

    ! -----------------------
    ! output to river
    ! -----------------------

    ! surface runoff is the sum of qflx_over, qflx_h2osfc_surf
    ! do g = bounds%begg,bounds%endg
    !    array(g) = water_inst%waterlnd2atmbulk_inst%qflx_rofliq_qsur_grc(g) + water_inst%waterlnd2atmbulk_inst%qflx_rofliq_h2osfc_grc(g)
    ! end do

    call state_setexport(exportState, 'Flrl_rofsur', bounds, input=water_inst%waterlnd2atmbulk_inst%qflx_rofliq_qsur_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! subsurface runoff is the sum of qflx_drain and qflx_perched_drain
    do g = bounds%begg,bounds%endg
       array(g) = water_inst%waterlnd2atmbulk_inst%qflx_rofliq_qsub_grc(g) + &
                  water_inst%waterlnd2atmbulk_inst%qflx_rofliq_drain_perched_grc(g)
    end do
    call state_setexport(exportState, 'Flrl_rofsub', bounds, input=array, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! qgwl sent individually to coupler
    call state_setexport(exportState, 'Flrl_rofgwl', bounds, input=water_inst%waterlnd2atmbulk_inst%qflx_rofliq_qgwl_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ice  sent individually to coupler
    call state_setexport(exportState, 'Flrl_rofi', bounds, input=water_inst%waterlnd2atmbulk_inst%qflx_rofice_grc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! irrigation flux to be removed from main channel storage (negative)
    call state_setexport(exportState, 'Flrl_irrig', bounds, input=water_inst%waterlnd2atmbulk_inst%qirrig_grc, minus=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine export_fields

  !===============================================================================

  subroutine fldlist_add(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname
    integer,          optional, intent(in)    :: ungridded_lbound
    integer,          optional, intent(in)    :: ungridded_ubound

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(lnd_import_export:fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
       return
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

  end subroutine fldlist_add

  !===============================================================================

  subroutine state_getimport(state, fldname, bounds, output, ungridded_index, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(in)    :: state
    type(bounds_type)   , intent(in)    :: bounds
    character(len=*)    , intent(in)    :: fldname
    real(r8)            , intent(out)   :: output(bounds%begg:bounds%endg)
    integer, optional   , intent(in)    :: ungridded_index
    integer             , intent(out)   :: rc

    ! local variables
    integer                     :: g, i,n

    integer                    :: fieldcount

    real(R8), pointer           :: fldptr1d(:)
    real(R8), pointer           :: fldptr2d(:,:)
    type(ESMF_StateItem_Flag)   :: itemFlag
    character(len=cs)           :: cvalue
    character(len=*), parameter :: subname='(lnd_import_export:state_getimport)'

    type (ESMF_FieldBundle)::  field
    type(ESMF_Field)            :: lfield
    type (ESMF_FieldBundle)::  fieldBundle
    logical                       :: isPresent
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Determine if field with name fldname exists in state

    !call ESMF_StateGet(state, trim(fldname), itemFlag, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! print out what is in our state???
    if (masterproc .and. debug > 0)  then
       write(iulog,F01)' Show me what is in the state? for  '//trim(fldname)
       call ESMF_StatePrint(state, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Determine if fieldbundle  exists in state
    call ESMF_StateGet(state, "c2l_fb", itemFlag, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)

    ! if fieldbundle exists then create output array - else do nothing
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then
       ! Get the field bundle???
       call ESMF_StateGet(state, "c2l_fb", fieldBundle, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_LogWrite(subname//'c2l_fb found and now ... getting '//trim(fldname), ESMF_LOGMSG_INFO)
       call ESMF_FieldBundleGet(fieldBundle,fieldName=trim(fldname), field=lfield,  isPresent=isPresent, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       !call ESMF_FieldBundleGet(fieldBundle,fieldName=trim(fldname), field=field, isPresent=isPresent, rc=rc)
       !call ESMF_FieldBundleGet(fieldBundle,field=field, rc=rc)
       !call ESMF_FieldBundleGet(fieldBundle, fieldCount=fieldCount, rc=rc)


       ! Now for error checking we can put ... if (isPresent...)
       ! get field pointer
       if (present(ungridded_index)) then
          write(cvalue,*) ungridded_index
          call ESMF_LogWrite(trim(subname)//": getting import for "//trim(fldname)//" index "//trim(cvalue), &
               ESMF_LOGMSG_INFO)
          call state_getfldptr(state, trim(fldname), fldptr2d=fldptr2d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call ESMF_LogWrite(trim(subname)//": getting import for "//trim(fldname),ESMF_LOGMSG_INFO)
          call state_getfldptr(state, trim(fldname), fldptr1d=fldptr1d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! determine output array
       if (present(ungridded_index)) then
          if (gridToFieldMap == 1) then
             do g = bounds%begg, bounds%endg
                n = g - bounds%begg + 1
                output(g) = fldptr2d(n,ungridded_index)
             end do
          else if (gridToFieldMap == 2) then
             do g = bounds%begg, bounds%endg
                n = g - bounds%begg + 1
                output(g) = fldptr2d(ungridded_index,n)
             end do
          end if
       else
          do g = bounds%begg, bounds%endg
             n = g - bounds%begg + 1
             output(g) = fldptr1d(n)
             if (masterproc .and. debug > 0 .and. get_nstep() < 5) then
                write(iulog,F02)' n, g , fldptr1d(n) '//trim(fldname)//' = ',n, g, fldptr1d(n)
             end if
          end do
       end if

       ! write debug output if appropriate
       if (masterproc .and. debug > 0 .and. get_nstep() < 5) then
          do g = bounds%begg,bounds%endg
             i = 1 + g - bounds%begg
             write(iulog,F02)'import: nstep, n, '//trim(fldname)//' = ',get_nstep(),i,output(g)
          end do
       end if

       ! check for nans
       call check_for_nans(output, trim(fldname), bounds%begg)
    end if

  end subroutine state_getimport

  !===============================================================================

  subroutine state_setexport(state, fldname, bounds, input, minus, ungridded_index, rc)

    ! ----------------------------------------------
    ! Map input array to export state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    type(bounds_type)   , intent(in)    :: bounds
    character(len=*)    , intent(in)    :: fldname
    real(r8)            , intent(in)    :: input(bounds%begg:bounds%endg)
    logical, optional   , intent(in)    :: minus
    integer, optional   , intent(in)    :: ungridded_index
    integer             , intent(out)   :: rc

    ! local variables
    integer                     :: g, i, n
    real(R8), pointer           :: fldptr1d(:)
    real(R8), pointer           :: fldptr2d(:,:)
    character(len=cs)           :: cvalue
    type(ESMF_StateItem_Flag)   :: itemFlag
    character(len=*), parameter :: subname='(lnd_import_export:state_setexport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine if field with name fldname exists in state
    call ESMF_StateGet(state, trim(fldname), itemFlag, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! if field exists then create output array - else do nothing
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then

       ! get field pointer
       if (present(ungridded_index)) then
          call ESMF_LogWrite(trim(subname)//": setting export for "//trim(fldname)//" index "//trim(cvalue), &
               ESMF_LOGMSG_INFO)
          call state_getfldptr(state, trim(fldname), fldptr2d=fldptr2d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call ESMF_LogWrite(trim(subname)//": setting export for "//trim(fldname), ESMF_LOGMSG_INFO)
          call state_getfldptr(state, trim(fldname), fldptr1d=fldptr1d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! TODO: if fillvalue = shr_const_spval the snowhl sent to the atm will have the spval over some points
       ! rather than 0 - this is very odd and needs to be understood
       ! fldptr(:) = fillvalue

       ! determine output array
       if (present(ungridded_index)) then
          if (gridToFieldMap == 1) then
             fldptr2d(:,ungridded_index) = 0._r8
             do g = bounds%begg, bounds%endg
                n = g - bounds%begg + 1
                fldptr2d(n,ungridded_index) = input(g)
             end do
             if (present(minus)) then
                fldptr2d(:,ungridded_index) = -fldptr2d(:,ungridded_index)
             end if
          else if (gridToFieldMap == 2) then
             fldptr2d(ungridded_index,:) = 0._r8
             do g = bounds%begg, bounds%endg
                n = g - bounds%begg + 1
                fldptr2d(ungridded_index,n) = input(g)
             end do
             if (present(minus)) then
                fldptr2d(ungridded_index,:) = -fldptr2d(ungridded_index,:)
             end if
          end if
       else
          fldptr1d(:) = 0._r8
          do g = bounds%begg, bounds%endg
             n = g - bounds%begg + 1
             fldptr1d(n) = input(g)
          end do
          if (present(minus)) then
             fldptr1d(:) = -fldptr1d(:)
          end if
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
    end if

  end subroutine state_setexport

  !===============================================================================

  subroutine state_getfldptr(State, fldname, fldptr1d, fldptr2d, rc)

    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State),             intent(in)    :: State
    character(len=*),             intent(in)    :: fldname
    real(R8), pointer, optional , intent(out)   :: fldptr1d(:)
    real(R8), pointer, optional , intent(out)   :: fldptr2d(:,:)
    integer,                      intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    type(ESMF_Mesh)             :: lmesh
    integer                     :: nnodes, nelements
    character(len=*), parameter :: subname='(lnd_import_export:state_getfldptr)'

    type(ESMF_StateItem_Flag)   :: itemFlag
    type(ESMF_FieldBundle)      :: fieldBundle
    logical                       :: isPresent
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine if this field bundle exist....
    ! TODO: combine the error checks....

    call ESMF_StateGet(state, "c2l_fb", itemFlag, rc=rc)
    !call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the fieldbundle from state...
    call ESMF_StateGet(state, "c2l_fb", fieldBundle, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return


    call ESMF_FieldBundleGet(fieldBundle,fieldName=trim(fldname), field=lfield,  isPresent=isPresent, rc=rc)
    !call ESMF_FieldBundleGet(fieldBundle,trim(fldname), lfield,  isPresent, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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

       if (present(fldptr1d)) then
          call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (masterproc .and. debug > 0)  then
             write(iulog,F01)' in '//trim(subname)//'fldptr1d for '//trim(fldname)//' is  '
          end if
          !print *, "FLDPTR1D is"
          !print *, FLDPTR1d
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

    ! Check if any input from mediator or output to mediator is NaN

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
