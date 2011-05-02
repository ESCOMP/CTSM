module histFldsMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: histFldsMod
!
! !DESCRIPTION:
! Module containing initialization of clm history fields and files
! This is the module that the user must modify in order to add new
! history fields or modify defaults associated with existing history
! fields.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
  public hist_initFlds ! Build master field list of all possible history
                       ! file fields
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 03/2003
! heald (11/28/06)
!
!EOP
!------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: hist_initFlds
!
! !INTERFACE:
  subroutine hist_initFlds()
!
! !DESCRIPTION:
! Build master field list of all possible fields in a history file.
! Each field has associated with it a ``long\_name'' netcdf attribute that
! describes what the field is, and a ``units'' attribute. A subroutine is
! called to add each field to the masterlist.
!
! !USES:
    use clmtype
    use clm_varcon , only : spval
    use clm_atmlnd , only : clm_a2l, atm_a2l, &
	                    adiag_arain, adiag_asnow, adiag_aflux, adiag_lflux
    use clm_varctl , only : create_glacier_mec_landunit, downscale
#if (defined RTM)
    use RunoffMod  , only : runoff, nt_rtm, rtm_tracers
#endif
    use histFileMod, only : hist_add_subscript, hist_addfld1d, hist_addfld2d, &
                            hist_printflds
    use surfrdMod  , only : crop_prog
#if (defined CASA)
    use CASAMod    , only : nlive, npools, npool_types
#endif
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 03/2003
! Mariana Vertenstein: Updated interface to create history fields 10/2003
!
!EOP
!-----------------------------------------------------------------------

    ! Determine what subscripts to add
    ! (uncomment the following call and modify it appropriately)

    ! call hist_add_subscript(subname='subscript_name', subdim=subscript_dim)

    ! NOTE: make a field not appear on the primary history tape by default -
    ! add the keyword to default='inactive' to the call to addfld_1d or addfld_2d

    ! Snow properties
    ! These will be vertically averaged over the snow profile

    call hist_addfld1d (fname='SNOWDP',  units='m',  &
         avgflag='A', long_name='snow height', &
         ptr_col=clm3%g%l%c%cps%snowdp, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSNO',  units='unitless',  &
         avgflag='A', long_name='fraction of ground covered by snow', &
         ptr_col=clm3%g%l%c%cps%frac_sno, c2l_scale_type='urbanf')

    ! Temperatures

    call hist_addfld1d (fname='TSA', units='K',  &
         avgflag='A', long_name='2m air temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m)

    call hist_addfld1d (fname='TSA_U', units='K',  &
         avgflag='A', long_name='Urban 2m air temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_u, set_nourb=spval)

    call hist_addfld1d (fname='TSA_R', units='K',  &
         avgflag='A', long_name='Rural 2m air temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_r, set_spec=spval)

    call hist_addfld1d(fname='TBUILD', units='K',  &
         avgflag='A', long_name='internal urban building temperature', &
         ptr_lunit=clm3%g%l%lps%t_building, set_nourb=spval, l2g_scale_type='unity')

    call hist_addfld1d (fname='TREFMNAV', units='K',  &
         avgflag='A', long_name='daily minimum of average 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_min)

    call hist_addfld1d (fname='TREFMXAV', units='K',  &
         avgflag='A', long_name='daily maximum of average 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_max)

    call hist_addfld1d (fname='TREFMNAV_U', units='K',  &
         avgflag='A', long_name='Urban daily minimum of average 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_min_u, set_nourb=spval)

    call hist_addfld1d (fname='TREFMXAV_U', units='K',  &
         avgflag='A', long_name='Urban daily maximum of average 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_max_u, set_nourb=spval)

    call hist_addfld1d (fname='TREFMNAV_R', units='K',  &
         avgflag='A', long_name='Rural daily minimum of average 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_min_r, set_spec=spval)

    call hist_addfld1d (fname='TREFMXAV_R', units='K',  &
         avgflag='A', long_name='Rural daily maximum of average 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_max_r, set_spec=spval)

    call hist_addfld1d (fname='TV', units='K',  &
         avgflag='A', long_name='vegetation temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_veg)

    call hist_addfld1d (fname='TV24', units='K',  &
         avgflag='A', long_name='vegetation temperature (last 24hrs)', &
         ptr_pft=clm3%g%l%c%p%pvs%t_veg24, default='inactive')

    call hist_addfld1d (fname='TV240', units='K',  &
         avgflag='A', long_name='vegetation temperature (last 240hrs)', &
         ptr_pft=clm3%g%l%c%p%pvs%t_veg240, default='inactive')

    call hist_addfld1d (fname='TG',  units='K',  &
         avgflag='A', long_name='ground temperature', &
         ptr_col=clm3%g%l%c%ces%t_grnd, c2l_scale_type='urbans')

    call hist_addfld1d (fname='TG_U', units='K',  &
         avgflag='A', long_name='Urban ground temperature', &
         ptr_col=clm3%g%l%c%ces%t_grnd_u, set_nourb=spval, c2l_scale_type='urbans')

    call hist_addfld1d (fname='TG_R', units='K',  &
         avgflag='A', long_name='Rural ground temperature', &
         ptr_col=clm3%g%l%c%ces%t_grnd_r, set_spec=spval)

    call hist_addfld1d (fname='HCSOI',  units='MJ',  &
         avgflag='A', long_name='soil heat content', &
         ptr_col=clm3%g%l%c%ces%hc_soi, set_lake=spval, set_urb=spval)

    call hist_addfld1d (fname='HC',  units='MJ',  &
         avgflag='A', long_name='heat content of soil/snow/lake', &
         ptr_col=clm3%g%l%c%ces%hc_soisno, set_urb=spval)

    call hist_addfld2d (fname='TSOI',  units='K', type2d='levgrnd', &
         avgflag='A', long_name='soil temperature', &
         ptr_col=clm3%g%l%c%ces%t_soisno, c2l_scale_type='urbanh')

    call hist_addfld1d (fname='TSOI_10CM',  units='K', &
         avgflag='A', long_name='soil temperature in top 10cm of soil', &
         ptr_col=clm3%g%l%c%ces%t_soi_10cm, set_urb=spval)

    call hist_addfld2d (fname='TLAKE',  units='K', type2d='levlak', &
         avgflag='A', long_name='lake temperature', &
         ptr_col=clm3%g%l%c%ces%t_lake)

    ! Specific humidity

    call hist_addfld1d (fname='Q2M', units='kg/kg',  &
         avgflag='A', long_name='2m specific humidity', &
         ptr_pft=clm3%g%l%c%p%pes%q_ref2m)

    ! Relative humidity

    call hist_addfld1d (fname='RH2M', units='%',  &
         avgflag='A', long_name='2m relative humidity', &
         ptr_pft=clm3%g%l%c%p%pes%rh_ref2m)

    call hist_addfld1d (fname='RH2M_U', units='%',  &
         avgflag='A', long_name='Urban 2m relative humidity', &
         ptr_pft=clm3%g%l%c%p%pes%rh_ref2m_u, set_nourb=spval)

    call hist_addfld1d (fname='RH2M_R', units='%',  &
         avgflag='A', long_name='Rural 2m specific humidity', &
         ptr_pft=clm3%g%l%c%p%pes%rh_ref2m_r, set_spec=spval)

    ! Wind

    call hist_addfld1d (fname='U10', units='m/s', &
         avgflag='A', long_name='10-m wind', &
         ptr_pft=clm3%g%l%c%p%pps%u10_clm)
    call hist_addfld1d (fname='VA', units='m/s', &
         avgflag='A', long_name='atmospheric wind speed plus convective velocity', &
         ptr_pft=clm3%g%l%c%p%pps%va, default='inactive')

    ! Surface radiation

    call hist_addfld1d (fname='SABV', units='watt/m^2',  &
         avgflag='A', long_name='solar rad absorbed by veg', &
         ptr_pft=clm3%g%l%c%p%pef%sabv, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='SABG', units='watt/m^2',  &
         avgflag='A', long_name='solar rad absorbed by ground', &
         ptr_pft=clm3%g%l%c%p%pef%sabg, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSDSVD', units='watt/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_vis_d)

    call hist_addfld1d (fname='FSDSND', units='watt/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_nir_d)

    call hist_addfld1d (fname='FSDSVI', units='watt/m^2',  &
         avgflag='A', long_name='diffuse vis incident solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_vis_i)

    call hist_addfld1d (fname='FSDSNI', units='watt/m^2',  &
         avgflag='A', long_name='diffuse nir incident solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_nir_i)

    call hist_addfld1d (fname='FSRVD', units='watt/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_vis_d, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSRND', units='watt/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_nir_d, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSRVI', units='watt/m^2',  &
         avgflag='A', long_name='diffuse vis reflected solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_vis_i, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSRNI', units='watt/m^2',  &
         avgflag='A', long_name='diffuse nir reflected solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_nir_i, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSDSVDLN', units='watt/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation at local noon', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_vis_d_ln)

    call hist_addfld1d (fname='FSDSNDLN', units='watt/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation at local noon', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_nir_d_ln)

    call hist_addfld1d (fname='FSRVDLN', units='watt/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation at local noon', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_vis_d_ln, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSRNDLN', units='watt/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation at local noon', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_nir_d_ln, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSA', units='watt/m^2',  &
         avgflag='A', long_name='absorbed solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsa, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSA_U', units='watt/m^2',  &
         avgflag='A', long_name='Urban absorbed solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsa_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='FSA_R', units='watt/m^2',  &
         avgflag='A', long_name='Rural absorbed solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsa_r, set_spec=spval)

    call hist_addfld1d (fname='FSR', units='watt/m^2',  &
         avgflag='A', long_name='reflected solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsr, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='SWup', units='watt/m^2',  &
         avgflag='A', long_name='upwelling shortwave radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsr, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='FIRA', units='watt/m^2',  &
         avgflag='A', long_name='net infrared (longwave) radiation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lwrad_net, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FIRA_U', units='watt/m^2',  &
         avgflag='A', long_name='Urban net infrared (longwave) radiation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lwrad_net_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='FIRA_R', units='watt/m^2',  &
         avgflag='A', long_name='Rural net infrared (longwave) radiation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lwrad_net_r, set_spec=spval)

    call hist_addfld1d (fname='FIRE', units='watt/m^2',  &
         avgflag='A', long_name='emitted infrared (longwave) radiation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lwrad_out, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='LWup', units='watt/m^2',  &
         avgflag='A', long_name='upwelling longwave radiation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lwrad_out, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='BUILDHEAT', units='watt/m^2',  &
         avgflag='A', long_name='heat flux from urban building interior to walls and roof', &
         ptr_col=clm3%g%l%c%cef%eflx_building_heat, set_nourb=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='URBAN_AC', units='watt/m^2',  &
         avgflag='A', long_name='urban air conditioning flux', &
         ptr_col=clm3%g%l%c%cef%eflx_urban_ac, set_nourb=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='URBAN_HEAT', units='watt/m^2',  &
         avgflag='A', long_name='urban heating flux', &
         ptr_col=clm3%g%l%c%cef%eflx_urban_heat, set_nourb=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='TRAFFICFLUX', units='watt/m^2',  &
         avgflag='A', long_name='sensible heat flux from urban traffic', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_traffic_pft, set_nourb=0._r8, c2l_scale_type='urbanf', &
         default='inactive')

    call hist_addfld1d (fname='WASTEHEAT', units='watt/m^2',  &
         avgflag='A', long_name='sensible heat flux from heating/cooling sources of urban waste heat', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_wasteheat_pft, set_nourb=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='HEAT_FROM_AC', units='watt/m^2',  &
         avgflag='A', long_name='sensible heat flux put into canyon due to heat removed from air conditioning', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_heat_from_ac_pft, set_nourb=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='Qanth', units='watt/m^2',  &
         avgflag='A', long_name='anthropogenic heat flux', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_anthro, set_nourb=0._r8, c2l_scale_type='urbanf', &
         default='inactive')

    call hist_addfld1d (fname='Rnet', units='watt/m^2',  &
         avgflag='A', long_name='net radiation', &
         ptr_pft=clm3%g%l%c%p%pef%netrad, c2l_scale_type='urbanf', &
         default='inactive')

    ! Solar zenith angle and solar declination angle

    call hist_addfld1d (fname='COSZEN', units='none', &
         avgflag='A', long_name='cosine of solar zenith angle', &
         ptr_col=clm3%g%l%c%cps%coszen, default='inactive')

    call hist_addfld1d (fname='DECL', units='radians', &
         avgflag='A', long_name='solar declination angle', &
         ptr_col=clm3%g%l%c%cps%decl, default='inactive')

    ! Surface energy fluxes

    call hist_addfld1d (fname='FCTR', units='watt/m^2',  &
         avgflag='A', long_name='canopy transpiration', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lh_vegt, set_lake=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FCEV', units='watt/m^2',  &
         avgflag='A', long_name='canopy evaporation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lh_vege, set_lake=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FGEV', units='watt/m^2',  &
         avgflag='A', long_name='ground evaporation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lh_grnd, c2l_scale_type='urbanf') 

    call hist_addfld1d (fname='FSH_NODYNLNDUSE', units='watt/m^2',  &
         avgflag='A', long_name='sensible heat', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_sh_tot, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSH', units='watt/m^2',  &
         avgflag='A', long_name='sensible heat', &
         ptr_lnd=clm3%g%gef%eflx_sh_totg)

    call hist_addfld1d (fname='FSH_U', units='watt/m^2',  &
         avgflag='A', long_name='Urban sensible heat', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_sh_tot_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='FSH_R', units='watt/m^2',  &
         avgflag='A', long_name='Rural sensible heat', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_sh_tot_r, set_spec=spval)

    call hist_addfld1d (fname='GC_HEAT1',  units='J/m^2',  &
         avgflag='A', long_name='initial gridcell total heat content', &
         ptr_lnd=clm3%g%ges%gc_heat1)

    call hist_addfld1d (fname='GC_HEAT2',  units='J/m^2',  &
         avgflag='A', long_name='post land cover change total heat content', &
         ptr_lnd=clm3%g%ges%gc_heat2, default='inactive')

    call hist_addfld1d (fname='EFLX_DYNBAL',  units='W/m^2',  &
         avgflag='A', long_name='dynamic land cover change conversion energy flux', &
         ptr_lnd=clm3%g%gef%eflx_dynbal)

    call hist_addfld1d (fname='Qh', units='watt/m^2',  &
         avgflag='A', long_name='sensible heat', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_sh_tot, c2l_scale_type='urbanf', &
         default = 'inactive')

    call hist_addfld1d (fname='Qle', units='watt/m^2',  &
         avgflag='A', long_name='total evaporation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lh_tot, c2l_scale_type='urbanf', &
         default = 'inactive')

    call hist_addfld1d (fname='EFLX_LH_TOT_U', units='watt/m^2',  &
         avgflag='A', long_name='Urban total evaporation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lh_tot_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='EFLX_LH_TOT_R', units='watt/m^2',  &
         avgflag='A', long_name='Rural total evaporation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lh_tot_r, set_spec=spval)

    call hist_addfld1d (fname='Qstor', units='watt/m^2',  &
         avgflag='A', long_name='storage heat flux (includes snowmelt)', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_soil_grnd, c2l_scale_type='urbanf', &
         default = 'inactive')

    call hist_addfld1d (fname='FSH_V', units='watt/m^2',  &
         avgflag='A', long_name='sensible heat from veg', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_sh_veg, set_lake=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSH_G', units='watt/m^2',  &
         avgflag='A', long_name='sensible heat from ground', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_sh_grnd, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FGR', units='watt/m^2',  &
         avgflag='A', long_name='heat flux into soil/snow including snow melt', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_soil_grnd, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FGR_U', units='watt/m^2',  &
         avgflag='A', long_name='Urban heat flux into soil/snow including snow melt', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_soil_grnd_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='FGR_R', units='watt/m^2',  &
         avgflag='A', long_name='Rural heat flux into soil/snow including snow melt', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_soil_grnd_r, set_spec=spval)

    call hist_addfld1d (fname='FSM',  units='watt/m^2',  &
         avgflag='A', long_name='snow melt heat flux', &
         ptr_col=clm3%g%l%c%cef%eflx_snomelt, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSM_U',  units='watt/m^2',  &
         avgflag='A', long_name='Urban snow melt heat flux', &
         ptr_col=clm3%g%l%c%cef%eflx_snomelt_u, c2l_scale_type='urbanf', set_nourb=spval)

    call hist_addfld1d (fname='FSM_R',  units='watt/m^2',  &
         avgflag='A', long_name='Rural snow melt heat flux', &
         ptr_col=clm3%g%l%c%cef%eflx_snomelt_r, set_spec=spval)

    call hist_addfld1d (fname='FGR12',  units='watt/m^2',  &
         avgflag='A', long_name='heat flux between soil layers 1 and 2', &
         ptr_col=clm3%g%l%c%cef%eflx_fgr12, set_lake=spval)

    call hist_addfld1d (fname='TAUX', units='kg/m/s^2',  &
         avgflag='A', long_name='zonal surface stress', &
         ptr_pft=clm3%g%l%c%p%pmf%taux)

    call hist_addfld1d (fname='Qtau', units='kg/m/s^2',  &
         avgflag='A', long_name='momentum flux', &
         ptr_pft=clm3%g%l%c%p%pmf%taux, default='inactive')

    call hist_addfld1d (fname='TAUY', units='kg/m/s^2',  &
         avgflag='A', long_name='meridional surface stress', &
         ptr_pft=clm3%g%l%c%p%pmf%tauy)

    ! Vegetation phenology

    call hist_addfld1d (fname='ELAI', units='m^2/m^2', &
          avgflag='A', long_name='exposed one-sided leaf area index', &
         ptr_pft=clm3%g%l%c%p%pps%elai)

    call hist_addfld1d (fname='ESAI', units='m^2/m^2', &
          avgflag='A', long_name='exposed one-sided stem area index', &
         ptr_pft=clm3%g%l%c%p%pps%esai)

    call hist_addfld1d (fname='LAISUN', units='none', &
         avgflag='A', long_name='sunlit projected leaf area index', &
         ptr_pft=clm3%g%l%c%p%pps%laisun, set_urb=0._r8)

    call hist_addfld1d (fname='LAISHA', units='none', &
         avgflag='A', long_name='shaded projected leaf area index', &
         ptr_pft=clm3%g%l%c%p%pps%laisha, set_urb=0._r8)

    call hist_addfld1d (fname='TLAI', units='none', &
         avgflag='A', long_name='total projected leaf area index', &
         ptr_pft=clm3%g%l%c%p%pps%tlai)

    call hist_addfld1d (fname='TSAI', units='none', &
         avgflag='A', long_name='total projected stem area index', &
         ptr_pft=clm3%g%l%c%p%pps%tsai)

    call hist_addfld1d (fname='SLASUN', units='m^2/gC', &
         avgflag='A', long_name='specific leaf area for sunlit canopy, projected area basis', &
         ptr_pft=clm3%g%l%c%p%pps%slasun, set_urb=0._r8, default='inactive')

    call hist_addfld1d (fname='SLASHA', units='m^2/gC', &
         avgflag='A', long_name='specific leaf area for shaded canopy, projected area basis', &
         ptr_pft=clm3%g%l%c%p%pps%slasha, set_urb=0._r8, default='inactive')

    ! Canopy physiology

    call hist_addfld1d (fname='RSSUN', units='s/m',  &
         avgflag='M', long_name='sunlit leaf stomatal resistance', &
         ptr_pft=clm3%g%l%c%p%pps%rssun, set_lake=spval, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='RSSHA', units='s/m',  &
         avgflag='M', long_name='shaded leaf stomatal resistance', &
         ptr_pft=clm3%g%l%c%p%pps%rssha, set_lake=spval, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='BTRAN', units='unitless',  &
         avgflag='A', long_name='transpiration beta factor', &
         ptr_pft=clm3%g%l%c%p%pps%btran, set_lake=spval, set_urb=spval)

    call hist_addfld1d (fname='FPSN', units='umol/m2s',  &
         avgflag='A', long_name='photosynthesis', &
         ptr_pft=clm3%g%l%c%p%pcf%fpsn, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='DSTFLXT', units='kg/m2/s',  &
         avgflag='A', long_name='total surface dust emission', &
         ptr_pft=clm3%g%l%c%p%pdf%flx_mss_vrt_dst_tot, set_lake=0._r8, set_urb=0._r8)
    call hist_addfld1d (fname='DPVLTRB1', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 1', &
         ptr_pft=clm3%g%l%c%p%pdf%vlc_trb_1, default='inactive')
    call hist_addfld1d (fname='DPVLTRB2', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 2', &
         ptr_pft=clm3%g%l%c%p%pdf%vlc_trb_2, default='inactive')
    call hist_addfld1d (fname='DPVLTRB3', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 3', &
         ptr_pft=clm3%g%l%c%p%pdf%vlc_trb_3, default='inactive')
    call hist_addfld1d (fname='DPVLTRB4', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 4', &
         ptr_pft=clm3%g%l%c%p%pdf%vlc_trb_4, default='inactive')

    call hist_addfld1d (fname='VOCFLXT', units='uGC/M2/H',  &
         avgflag='A', long_name='total VOC flux into atmosphere', &
         ptr_pft=clm3%g%l%c%p%pvf%vocflx_tot, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='ISOPRENE', units='uGC/M2/H',  &
         avgflag='A', long_name='isoprene flux', &
         ptr_pft=clm3%g%l%c%p%pvf%vocflx_1, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='MONOTERP', units='uGC/M2/H',  &
         avgflag='A', long_name='monoterpene flux', &
         ptr_pft=clm3%g%l%c%p%pvf%vocflx_2, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='OVOC', units='uGC/M2/H',  &
         avgflag='A', long_name='other VOC flux', &
         ptr_pft=clm3%g%l%c%p%pvf%vocflx_3, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='ORVOC', units='uGC/M2/H',  &
         avgflag='A', long_name='other reactive VOC flux', &
         ptr_pft=clm3%g%l%c%p%pvf%vocflx_4, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='BIOGENCO', units='uGC/M2/H',  &
         avgflag='A', long_name='biogenic CO flux', &
         ptr_pft=clm3%g%l%c%p%pvf%vocflx_5, set_lake=0._r8, set_urb=0._r8)

    call hist_addfld1d (fname='EOPT', units='non',  &
         avgflag='A', long_name='Eopt coefficient for VOC calc', &
         ptr_pft=clm3%g%l%c%p%pvf%Eopt_out, set_lake=0._r8, default='inactive')

    call hist_addfld1d (fname='TOPT', units='non',  &
         avgflag='A', long_name='topt coefficient for VOC calc', &
         ptr_pft=clm3%g%l%c%p%pvf%topt_out, set_lake=0._r8, default='inactive')

    call hist_addfld1d (fname='ALPHA', units='non',  &
         avgflag='A', long_name='alpha coefficient for VOC calc', &
         ptr_pft=clm3%g%l%c%p%pvf%alpha_out, set_lake=0._r8, default='inactive')

    call hist_addfld1d (fname='CP', units='non',  &
         avgflag='A', long_name='cp coefficient for VOC calc', &
         ptr_pft=clm3%g%l%c%p%pvf%cp_out, set_lake=0._r8, default='inactive')

    call hist_addfld1d (fname='FSUN24', units='K',  &
         avgflag='A', long_name='fraction sunlit (last 24hrs)', &
         ptr_pft=clm3%g%l%c%p%pvs%fsun24, default='inactive')

    call hist_addfld1d (fname='FSUN240', units='K',  &
         avgflag='A', long_name='fraction sunlit (last 240hrs)', &
         ptr_pft=clm3%g%l%c%p%pvs%fsun240, default='inactive')

    call hist_addfld1d (fname='FSI24', units='K',  &
         avgflag='A', long_name='indirect radiation (last 24hrs)', &
         ptr_pft=clm3%g%l%c%p%pvs%fsi24, default='inactive')

    call hist_addfld1d (fname='FSI240', units='K',  &
         avgflag='A', long_name='indirect radiation (last 240hrs)', &
         ptr_pft=clm3%g%l%c%p%pvs%fsi240, default='inactive')

    call hist_addfld1d (fname='FSD24', units='K',  &
         avgflag='A', long_name='direct radiation (last 24hrs)', &
         ptr_pft=clm3%g%l%c%p%pvs%fsd24, default='inactive')

    call hist_addfld1d (fname='FSD240', units='K',  &
         avgflag='A', long_name='direct radiation (last 240hrs)', &
         ptr_pft=clm3%g%l%c%p%pvs%fsd240, default='inactive')

    call hist_addfld1d (fname='PAR_sun', units='umol/m2/s', &
         avgflag='A', long_name='sunlit PAR', &
         ptr_pft=clm3%g%l%c%p%pvf%paru_out, set_lake=0._r8, default='inactive')

    call hist_addfld1d (fname='PAR24_sun', units='umol/m2/s', &
         avgflag='A', long_name='sunlit PAR (24 hrs)', &
         ptr_pft=clm3%g%l%c%p%pvf%par24u_out, set_lake=0._r8, default='inactive')

    call hist_addfld1d (fname='PAR240_sun', units='umol/m2/s', &
         avgflag='A', long_name='sunlit PAR (240 hrs)', &
         ptr_pft=clm3%g%l%c%p%pvf%par240u_out, set_lake=0._r8, default='inactive')

    call hist_addfld1d (fname='PAR_shade', units='umol/m2/s', &
         avgflag='A', long_name='shade PAR', &
         ptr_pft=clm3%g%l%c%p%pvf%para_out, set_lake=0._r8, default='inactive')

    call hist_addfld1d (fname='PAR24_shade', units='umol/m2/s', &
         avgflag='A', long_name='shade PAR (24 hrs)', &
         ptr_pft=clm3%g%l%c%p%pvf%par24a_out, set_lake=0._r8, default='inactive')

    call hist_addfld1d (fname='PAR240_shade', units='umol/m2/s', &
         avgflag='A', long_name='shade PAR (240 hrs)', &
         ptr_pft=clm3%g%l%c%p%pvf%par240a_out, set_lake=0._r8, default='inactive')

    call hist_addfld1d (fname='GAMMA', units='non',  &
         avgflag='A', long_name='total gamma for VOC calc', &
         ptr_pft=clm3%g%l%c%p%pvf%gamma_out, set_lake=0._r8, default='inactive')

    call hist_addfld1d (fname='GAMMAL', units='non',  &
         avgflag='A', long_name='gamma L for VOC calc', &
         ptr_pft=clm3%g%l%c%p%pvf%gammaL_out, set_lake=0._r8, default='inactive')

    call hist_addfld1d (fname='GAMMAT', units='non',  &
         avgflag='A', long_name='gamma T for VOC calc', &
         ptr_pft=clm3%g%l%c%p%pvf%gammaT_out, set_lake=0._r8, default='inactive')

    call hist_addfld1d (fname='GAMMAP', units='non',  &
         avgflag='A', long_name='gamma P for VOC calc', &
         ptr_pft=clm3%g%l%c%p%pvf%gammaP_out, set_lake=0._r8, default='inactive')

    call hist_addfld1d (fname='GAMMAA', units='non',  &
         avgflag='A', long_name='gamma A for VOC calc', &
         ptr_pft=clm3%g%l%c%p%pvf%gammaA_out, set_lake=0._r8, default='inactive')

    call hist_addfld1d (fname='GAMMAS', units='non',  &
         avgflag='A', long_name='gamma S for VOC calc', &
         ptr_pft=clm3%g%l%c%p%pvf%gammaS_out, set_lake=0._r8, default='inactive')

    ! Hydrology

    call hist_addfld1d (fname='SoilAlpha',  units='unitless',  &
         avgflag='A', long_name='factor limiting ground evap', &
         ptr_col=clm3%g%l%c%cws%soilalpha, set_urb=spval)

    call hist_addfld1d (fname='SoilAlpha_U',  units='unitless',  &
         avgflag='A', long_name='urban factor limiting ground evap', &
         ptr_col=clm3%g%l%c%cws%soilalpha_u, set_nourb=spval)

    call hist_addfld1d (fname='FCOV',  units='unitless',  &
         avgflag='A', long_name='fractional impermeable area', &
         ptr_col=clm3%g%l%c%cws%fcov, c2l_scale_type='urbanh')
    call hist_addfld1d (fname='FSAT',  units='unitless',  &
         avgflag='A', long_name='fractional area with water table at surface', &
         ptr_col=clm3%g%l%c%cws%fsat, c2l_scale_type='urbanh')
    call hist_addfld1d (fname='ZWT',  units='m',  &
         avgflag='A', long_name='water table depth', &
         ptr_col=clm3%g%l%c%cws%zwt, c2l_scale_type='urbanh')
    !call hist_addfld1d (fname='FROST_TABLE',  units='m',  &
    !     avgflag='A', long_name='frost table depth', &
    !     ptr_col=clm3%g%l%c%cws%frost_table, c2l_scale_type='urbanh')
    !call hist_addfld1d (fname='ZWT_PERCH',  units='m',  &
    !     avgflag='A', long_name='perched water table depth', &
    !     ptr_col=clm3%g%l%c%cws%zwt_perched, c2l_scale_type='urbanh')
    !call hist_addfld1d (fname='QDRAI_PERCH',  units='mm/s',  &
    !     avgflag='A', long_name='perched wt drainage', &
    !     ptr_col=clm3%g%l%c%cwf%qflx_drain_perched, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='WA',  units='mm',  &
         avgflag='A', long_name='water in the unconfined aquifer', &
         ptr_col=clm3%g%l%c%cws%wa, c2l_scale_type='urbanh')

    call hist_addfld1d (fname='WT',  units='mm',  &
         avgflag='A', long_name='total water storage (unsaturated soil water + groundwater)', &
         ptr_col=clm3%g%l%c%cws%wt, c2l_scale_type='urbanh')

    call hist_addfld1d (fname='QCHARGE',  units='mm/s',  &
         avgflag='A', long_name='aquifer recharge rate', &
         ptr_col=clm3%g%l%c%cws%qcharge, c2l_scale_type='urbanh')

    call hist_addfld2d (fname='SMP',  units='mm', type2d='levgrnd',  &
         avgflag='A', long_name='soil matric potential', &
         ptr_col=clm3%g%l%c%cws%smp_l, set_spec=spval, default='inactive')

    call hist_addfld2d (fname='HK',  units='mm/s', type2d='levgrnd',  &
         avgflag='A', long_name='hydraulic conductivity', &
         ptr_col=clm3%g%l%c%cws%hk_l, set_spec=spval, default='inactive')

    call hist_addfld1d (fname='H2OSNO',  units='mm',  &
         avgflag='A', long_name='snow depth (liquid water)', &
         ptr_col=clm3%g%l%c%cws%h2osno, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='H2OCAN', units='mm',  &
         avgflag='A', long_name='intercepted water', &
         ptr_pft=clm3%g%l%c%p%pws%h2ocan, set_lake=0._r8)

    call hist_addfld2d (fname='H2OSOI',  units='mm3/mm3', type2d='levgrnd', &
         avgflag='A', long_name='volumetric soil water', &
         ptr_col=clm3%g%l%c%cws%h2osoi_vol, c2l_scale_type='urbanh')

    call hist_addfld2d (fname='SOILLIQ',  units='kg/m2', type2d='levgrnd', &
         avgflag='A', long_name='soil liquid water', &
         ptr_col=clm3%g%l%c%cws%h2osoi_liq, c2l_scale_type='urbanh')

    call hist_addfld2d (fname='SOILICE',  units='kg/m2', type2d='levgrnd', &
         avgflag='A', long_name='soil ice', &
         ptr_col=clm3%g%l%c%cws%h2osoi_ice, c2l_scale_type='urbanh')

    call hist_addfld1d (fname='SOILWATER_10CM',  units='kg/m2', &
         avgflag='A', long_name='soil liquid water + ice in top 10cm of soil', &
         ptr_col=clm3%g%l%c%cws%h2osoi_liqice_10cm, set_urb=spval)

    call hist_addfld1d (fname='SNOWLIQ',  units='kg/m2',  &
         avgflag='A', long_name='snow liquid water', &
         ptr_col=clm3%g%l%c%cws%snowliq)

    call hist_addfld1d (fname='SNOWICE',  units='kg/m2',  &
         avgflag='A', long_name='snow ice', &
         ptr_col=clm3%g%l%c%cws%snowice)

    call hist_addfld1d (fname='QTOPSOIL',  units='mm/s',  &
         avgflag='A', long_name='water input to surface', &
         ptr_col=clm3%g%l%c%cwf%qflx_top_soil, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='QINFL',  units='mm/s',  &
         avgflag='A', long_name='infiltration', &
         ptr_col=clm3%g%l%c%cwf%qflx_infl, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QOVER',  units='mm/s',  &
         avgflag='A', long_name='surface runoff', &
         ptr_col=clm3%g%l%c%cwf%qflx_surf, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QRGWL',  units='mm/s',  &
         avgflag='A', long_name='surface runoff at glaciers (liquid only), wetlands, lakes', &
         ptr_col=clm3%g%l%c%cwf%qflx_qrgwl, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QSNWCPLIQ', units='mm H2O/s', &
         avgflag='A', long_name='excess rainfall due to snow capping', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_snwcp_liq, default='inactive')

    call hist_addfld1d (fname='QSNWCPICE_NODYNLNDUSE', units='mm H2O/s', &
         avgflag='A', long_name='excess snowfall due to snow capping', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_snwcp_ice)

    call hist_addfld1d (fname='QSNWCPICE',  units='mm/s',  &
         avgflag='A', long_name='excess snowfall due to snow capping', &
         ptr_lnd=clm3%g%gwf%qflx_snwcp_iceg)

    call hist_addfld1d (fname='QDRAI',  units='mm/s',  &
         avgflag='A', long_name='sub-surface drainage', &
         ptr_col=clm3%g%l%c%cwf%qflx_drain, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QRUNOFF_NODYNLNDUSE',  units='mm/s',  &
         avgflag='A', long_name='total liquid runoff (does not include QSNWCPICE)', &
         ptr_col=clm3%g%l%c%cwf%qflx_runoff, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QRUNOFF',  units='mm/s',  &
         avgflag='A', long_name='total liquid runoff (does not include QSNWCPICE)', &
         ptr_lnd=clm3%g%gwf%qflx_runoffg)

    call hist_addfld1d (fname='GC_LIQ1',  units='mm',  &
         avgflag='A', long_name='initial gridcell total liq content', &
         ptr_lnd=clm3%g%gws%gc_liq1)

    call hist_addfld1d (fname='GC_LIQ2',  units='mm',  &  
         avgflag='A', long_name='initial gridcell total liq content', &              
         ptr_lnd=clm3%g%gws%gc_liq2, default='inactive')     

    call hist_addfld1d (fname='QFLX_LIQ_DYNBAL',  units='mm/s',  &  
         avgflag='A', long_name='liq dynamic land cover change conversion runoff flux', &              
         ptr_lnd=clm3%g%gwf%qflx_liq_dynbal)     

    call hist_addfld1d (fname='GC_ICE1',  units='mm',  &  
         avgflag='A', long_name='initial gridcell total ice content', &              
         ptr_lnd=clm3%g%gws%gc_ice1)     

    call hist_addfld1d (fname='GC_ICE2',  units='mm',  &  
         avgflag='A', long_name='post land cover change total ice content', &              
         ptr_lnd=clm3%g%gws%gc_ice2, default='inactive')

    call hist_addfld1d (fname='QFLX_ICE_DYNBAL',  units='mm/s',  &
         avgflag='A', long_name='ice dynamic land cover change conversion runoff flux', &                                   
         ptr_lnd=clm3%g%gwf%qflx_ice_dynbal)

    call hist_addfld1d (fname='QRUNOFF_U', units='mm/s',  &
         avgflag='A', long_name='Urban total runoff', &
         ptr_col=clm3%g%l%c%cwf%qflx_runoff_u, set_nourb=spval, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QRUNOFF_R', units='mm/s',  &
         avgflag='A', long_name='Rural total runoff', &
         ptr_col=clm3%g%l%c%cwf%qflx_runoff_r, set_spec=spval)

    call hist_addfld1d (fname='QINTR', units='mm/s',  &
         avgflag='A', long_name='interception', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_prec_intr, set_lake=0._r8)

    call hist_addfld1d (fname='QDRIP', units='mm/s',  &
         avgflag='A', long_name='throughfall', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_prec_grnd, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QMELT',  units='mm/s',  &
         avgflag='A', long_name='snow melt', &
         ptr_col=clm3%g%l%c%cwf%qflx_snomelt, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QSOIL', units='mm/s',  &
         avgflag='A', long_name='ground evaporation', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_evap_soi, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QVEGE', units='mm/s',  &
         avgflag='A', long_name='canopy evaporation', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_evap_can, set_lake=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QVEGT', units='mm/s',  &
         avgflag='A', long_name='canopy transpiration', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_tran_veg, set_lake=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QIRRIG', units='mm/s', &
         avgflag='A', long_name='water added through irrigation', &
         ptr_col=clm3%g%l%c%cwf%qflx_irrig, set_lake=0._r8)

    if (create_glacier_mec_landunit) then

       call hist_addfld1d (fname='QICE',  units='mm/s',  &
            avgflag='A', long_name='ice growth/melt', &
            ptr_col=clm3%g%l%c%cwf%qflx_glcice, set_noglcmec=spval)

       call hist_addfld1d (fname='QICEYR',  units='mm/s',  &
            avgflag='A', long_name='ice growth/melt', &
            ptr_col=clm3%g%l%c%cwf%qflx_glcice, set_noglcmec=spval)

       call hist_addfld1d (fname='gris_mask',  units='unitless',  &
            avgflag='A', long_name='Greenland mask', &
            ptr_gcell=clm3%g%gris_mask)

       call hist_addfld1d (fname='gris_area',  units='km^2',  &
            avgflag='A', long_name='Greenland ice area', &
            ptr_gcell=clm3%g%gris_area)

       call hist_addfld1d (fname='aais_mask',  units='unitless',  &
            avgflag='A', long_name='Antarctic mask', &
            ptr_gcell=clm3%g%aais_mask)

       call hist_addfld1d (fname='aais_area',  units='km^2',  &
            avgflag='A', long_name='Antarctic ice area', &
            ptr_gcell=clm3%g%aais_area)

   endif

#if (defined RTM)
    ! RTM River Routing

    call hist_addfld1d (fname='QCHANR', units='m3/s',  &
      avgflag='A', long_name='RTM river flow: '//trim(rtm_tracers(1)), &
      ptr_rof=runoff%runofflnd_nt1)

    call hist_addfld1d (fname='QCHANR'//'_'//trim(rtm_tracers(2)), units='m3/s',  &
      avgflag='A', long_name='RTM river flow: '//trim(rtm_tracers(2)), &
      ptr_rof=runoff%runofflnd_nt2)

    call hist_addfld1d (fname='QCHOCNR', units='m3/s', &
      avgflag='A', long_name='RTM river discharge into ocean: '//trim(rtm_tracers(1)), &
      ptr_rof=runoff%runoffocn_nt1)

    call hist_addfld1d (fname='QCHOCNR'//'_'//trim(rtm_tracers(2)), units='m3/s', &
      avgflag='A', long_name='RTM river discharge into ocean: '//trim(rtm_tracers(2)), &
      ptr_rof=runoff%runoffocn_nt2)

    call hist_addfld1d (fname='VOLR', units='m3',  &
         avgflag='A', long_name='RTM storage: '//trim(rtm_tracers(1)), &
         ptr_rof=runoff%volr_nt1)

    call hist_addfld1d (fname='VOLR'//'_'//trim(rtm_tracers(2)), units='m3',  &
         avgflag='A', long_name='RTM storage: '//trim(rtm_tracers(2)), &
         ptr_rof=runoff%volr_nt2, default='inactive')

    call hist_addfld1d (fname='DVOLRDT_LND', units='mm/s',  &
      avgflag='A', long_name='RTM land change in storage: '//trim(rtm_tracers(1)), &
      ptr_rof=runoff%dvolrdtlnd_nt1, default='inactive')

    call hist_addfld1d (fname='DVOLRDT_LND'//'_'//trim(rtm_tracers(2)), units='mm/s',  &
      avgflag='A', long_name='RTM land change in storage: '//trim(rtm_tracers(2)), &
      ptr_rof=runoff%dvolrdtlnd_nt2, default='inactive')

    call hist_addfld1d (fname='DVOLRDT_OCN', units='mm/s',  &
      avgflag='A', long_name='RTM ocean change of storage: '//trim(rtm_tracers(1)), &
      ptr_rof=runoff%dvolrdtocn_nt1, default='inactive')

    call hist_addfld1d (fname='DVOLRDT_OCN'//'_'//trim(rtm_tracers(2)), units='mm/s',  &
      avgflag='A', long_name='RTM ocean change of storage: '//trim(rtm_tracers(2)), &
      ptr_rof=runoff%dvolrdtocn_nt2, default='inactive')
#endif

    ! Water and energy balance checks

    call hist_addfld1d (fname='ERRSOI',  units='watt/m^2',  &
         avgflag='A', long_name='soil/lake energy conservation error', &
         ptr_col=clm3%g%l%c%cebal%errsoi)

    call hist_addfld1d (fname='ERRSEB',  units='watt/m^2',  &
         avgflag='A', long_name='surface energy conservation error', &
         ptr_pft=clm3%g%l%c%p%pebal%errseb)

    call hist_addfld1d (fname='ERRSOL',  units='watt/m^2',  &
         avgflag='A', long_name='solar radiation conservation error', &
         ptr_pft=clm3%g%l%c%p%pebal%errsol, set_urb=spval)

    call hist_addfld1d (fname='ERRH2O', units='mm',  &
         avgflag='A', long_name='total water conservation error', &
         ptr_col=clm3%g%l%c%cwbal%errh2o)

    ! Atmospheric forcing

    if (downscale) then
       call hist_addfld1d (fname='RAINATM', units='mm/s',  &
            avgflag='A', long_name='atmospheric rain forcing', &
            ptr_atm=atm_a2l%forc_rain)
       call hist_addfld1d (fname='SNOWATM', units='mm/s',  &
            avgflag='A', long_name='atmospheric snow forcing', &
            ptr_atm=atm_a2l%forc_snow)
    else
       call hist_addfld1d (fname='RAINATM', units='mm/s',  &
            avgflag='A', long_name='atmospheric rain forcing', &
            ptr_lnd=clm_a2l%forc_rain)
       call hist_addfld1d (fname='SNOWATM', units='mm/s',  &
            avgflag='A', long_name='atmospheric snow forcing', &
            ptr_lnd=clm_a2l%forc_snow)
    end if

    call hist_addfld1d (fname='RAINFM2A', units='mm/s',  &
         avgflag='A', long_name='land rain on atm grid', &
         ptr_atm=adiag_arain)

    call hist_addfld1d (fname='SNOWFM2A', units='mm/s',  &
         avgflag='A', long_name='land snow on atm grid', &
         ptr_atm=adiag_asnow)

    call hist_addfld1d (fname='FLUXFM2A', units='W/m2',  &
         avgflag='A', long_name='heat flux for rain to snow conversion', &
         ptr_atm=adiag_aflux)

    call hist_addfld1d (fname='FLUXFMLND', units='W/m2',  &
         avgflag='A', long_name='heat flux from rain to snow conversion', &
         ptr_lnd=adiag_lflux)

    call hist_addfld1d (fname='RAIN', units='mm/s',  &
         avgflag='A', long_name='atmospheric rain', &
         ptr_lnd=clm_a2l%forc_rain)

    call hist_addfld1d (fname='SNOW', units='mm/s',  &
         avgflag='A', long_name='atmospheric snow', &
         ptr_lnd=clm_a2l%forc_snow)

    call hist_addfld1d (fname='TBOT', units='K',  &
         avgflag='A', long_name='atmospheric air temperature', &
         ptr_lnd=clm_a2l%forc_t)

    call hist_addfld1d (fname='THBOT', units='K',  &
         avgflag='A', long_name='atmospheric air potential temperature', &
         ptr_lnd=clm_a2l%forc_th)

    call hist_addfld1d (fname='WIND', units='m/s',  &
         avgflag='A', long_name='atmospheric wind velocity magnitude', &
         ptr_lnd=clm_a2l%forc_wind)

    call hist_addfld1d (fname='Wind', units='m/s',  &
         avgflag='A', long_name='atmospheric wind velocity magnitude', &
         ptr_gcell=clm_a2l%forc_wind, default = 'inactive')

    call hist_addfld1d (fname='Tair', units='K',  &
         avgflag='A', long_name='atmospheric air temperature', &
         ptr_gcell=clm_a2l%forc_t, default='inactive')

    call hist_addfld1d (fname='PSurf', units='Pa',  &
         avgflag='A', long_name='surface pressure', &
         ptr_gcell=clm_a2l%forc_pbot, default='inactive')

    call hist_addfld1d (fname='Rainf', units='mm/s',  &
         avgflag='A', long_name='atmospheric rain', &
         ptr_gcell=clm_a2l%forc_rain, default='inactive')

    call hist_addfld1d (fname='SWdown', units='watt/m^2',  &
         avgflag='A', long_name='atmospheric incident solar radiation', &
         ptr_gcell=clm_a2l%forc_solar, default='inactive')

    call hist_addfld1d (fname='LWdown', units='watt/m^2',  &
         avgflag='A', long_name='atmospheric longwave radiation', &
         ptr_gcell=clm_a2l%forc_lwrad, default='inactive')

    call hist_addfld1d (fname='RH', units='%',  &
         avgflag='A', long_name='atmospheric relative humidity', &
         ptr_gcell=clm_a2l%forc_rh, default='inactive')

    call hist_addfld1d (fname='QBOT', units='kg/kg',  &
         avgflag='A', long_name='atmospheric specific humidity', &
         ptr_lnd=clm_a2l%forc_q)

    call hist_addfld1d (fname='Qair', units='kg/kg',  &
         avgflag='A', long_name='atmospheric specific humidity', &
         ptr_lnd=clm_a2l%forc_q, default='inactive')

    call hist_addfld1d (fname='ZBOT', units='m',  &
         avgflag='A', long_name='atmospheric reference height', &
         ptr_lnd=clm_a2l%forc_hgt)

    call hist_addfld1d (fname='FLDS', units='watt/m^2',  &
         avgflag='A', long_name='atmospheric longwave radiation', &
         ptr_lnd=clm_a2l%forc_lwrad)

    call hist_addfld1d (fname='FSDS', units='watt/m^2',  &
         avgflag='A', long_name='atmospheric incident solar radiation', &
         ptr_lnd=clm_a2l%forc_solar)

    call hist_addfld1d (fname='PCO2', units='Pa',  &
         avgflag='A', long_name='atmospheric partial pressure of CO2', &
         ptr_lnd=clm_a2l%forc_pco2)

    call hist_addfld1d (fname='PBOT', units='Pa',  &
         avgflag='A', long_name='atmospheric pressure', &
         ptr_lnd=clm_a2l%forc_pbot)

#if (defined CNDV) || (defined CROP)
    call hist_addfld1d (fname='T10', units='K',  &
         avgflag='A', long_name='10-day running mean of 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t10)
#endif

#if (defined CNDV)
    call hist_addfld1d (fname='TDA', units='K',  &
         avgflag='A', long_name='daily average 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pdgvs%t_mo)

    call hist_addfld1d (fname='AGDD', units='K',  &
         avgflag='A', long_name='growing degree-days base 5C', &
         ptr_pft=clm3%g%l%c%p%pdgvs%agdd)
#endif

#if (defined CASA) || (defined CN)
    call hist_addfld2d (fname='SOILPSI', units='MPa', type2d='levgrnd', &
         avgflag='A', long_name='soil water potential in each soil layer', &
         ptr_col=clm3%g%l%c%cps%soilpsi)
#endif

#if (defined CN)
    ! add history fields for all CN variables, always set as default='inactive'

    if ( crop_prog )then

       call hist_addfld1d (fname='A5TMIN', units='K',  &
            avgflag='A', long_name='5-day running mean of min 2-m temperature', &
            ptr_pft=clm3%g%l%c%p%pes%a5tmin, default='inactive')

       call hist_addfld1d (fname='A10TMIN', units='K',  &
            avgflag='A', long_name='10-day running mean of min 2-m temperature', &
            ptr_pft=clm3%g%l%c%p%pes%a10tmin, default='inactive')

    end if
    
    !-------------------------------
    ! C state variables - native to PFT 
    !-------------------------------
    ! add history fields for all CLAMP CN variables

    call hist_addfld1d (fname='WOODC', units='gC/m^2', &
             avgflag='A', long_name='wood C', &
             ptr_pft=clm3%g%l%c%p%pcs%woodc)
    
    call hist_addfld1d (fname='LEAFC', units='gC/m^2', &
         avgflag='A', long_name='leaf C', &
         ptr_pft=clm3%g%l%c%p%pcs%leafc)

    call hist_addfld1d (fname='LEAFC_STORAGE', units='gC/m^2', &
         avgflag='A', long_name='leaf C storage', &
         ptr_pft=clm3%g%l%c%p%pcs%leafc_storage, default='inactive')

    call hist_addfld1d (fname='LEAFC_XFER', units='gC/m^2', &
         avgflag='A', long_name='leaf C transfer', &
         ptr_pft=clm3%g%l%c%p%pcs%leafc_xfer, default='inactive')

    call hist_addfld1d (fname='FROOTC', units='gC/m^2', &
         avgflag='A', long_name='fine root C', &
         ptr_pft=clm3%g%l%c%p%pcs%frootc)

    call hist_addfld1d (fname='FROOTC_STORAGE', units='gC/m^2', &
         avgflag='A', long_name='fine root C storage', &
         ptr_pft=clm3%g%l%c%p%pcs%frootc_storage, default='inactive')

    call hist_addfld1d (fname='FROOTC_XFER', units='gC/m^2', &
         avgflag='A', long_name='fine root C transfer', &
         ptr_pft=clm3%g%l%c%p%pcs%frootc_xfer, default='inactive')

    call hist_addfld1d (fname='LIVESTEMC', units='gC/m^2', &
         avgflag='A', long_name='live stem C', &
         ptr_pft=clm3%g%l%c%p%pcs%livestemc)

    call hist_addfld1d (fname='LIVESTEMC_STORAGE', units='gC/m^2', &
         avgflag='A', long_name='live stem C storage', &
         ptr_pft=clm3%g%l%c%p%pcs%livestemc_storage, default='inactive')

    call hist_addfld1d (fname='LIVESTEMC_XFER', units='gC/m^2', &
         avgflag='A', long_name='live stem C transfer', &
         ptr_pft=clm3%g%l%c%p%pcs%livestemc_xfer, default='inactive')

    call hist_addfld1d (fname='DEADSTEMC', units='gC/m^2', &
         avgflag='A', long_name='dead stem C', &
         ptr_pft=clm3%g%l%c%p%pcs%deadstemc)

    call hist_addfld1d (fname='DEADSTEMC_STORAGE', units='gC/m^2', &
         avgflag='A', long_name='dead stem C storage', &
         ptr_pft=clm3%g%l%c%p%pcs%deadstemc_storage, default='inactive')

    call hist_addfld1d (fname='DEADSTEMC_XFER', units='gC/m^2', &
         avgflag='A', long_name='dead stem C transfer', &
         ptr_pft=clm3%g%l%c%p%pcs%deadstemc_xfer, default='inactive')

    call hist_addfld1d (fname='LIVECROOTC', units='gC/m^2', &
         avgflag='A', long_name='live coarse root C', &
         ptr_pft=clm3%g%l%c%p%pcs%livecrootc)

    call hist_addfld1d (fname='LIVECROOTC_STORAGE', units='gC/m^2', &
         avgflag='A', long_name='live coarse root C storage', &
         ptr_pft=clm3%g%l%c%p%pcs%livecrootc_storage, default='inactive')

    call hist_addfld1d (fname='LIVECROOTC_XFER', units='gC/m^2', &
         avgflag='A', long_name='live coarse root C transfer', &
         ptr_pft=clm3%g%l%c%p%pcs%livecrootc_xfer, default='inactive')

    call hist_addfld1d (fname='DEADCROOTC', units='gC/m^2', &
         avgflag='A', long_name='dead coarse root C', &
         ptr_pft=clm3%g%l%c%p%pcs%deadcrootc)

    call hist_addfld1d (fname='DEADCROOTC_STORAGE', units='gC/m^2', &
         avgflag='A', long_name='dead coarse root C storage', &
         ptr_pft=clm3%g%l%c%p%pcs%deadcrootc_storage,  default='inactive')

    call hist_addfld1d (fname='DEADCROOTC_XFER', units='gC/m^2', &
         avgflag='A', long_name='dead coarse root C transfer', &
         ptr_pft=clm3%g%l%c%p%pcs%deadcrootc_xfer, default='inactive')

    call hist_addfld1d (fname='GRESP_STORAGE', units='gC/m^2', &
         avgflag='A', long_name='growth respiration storage', &
         ptr_pft=clm3%g%l%c%p%pcs%gresp_storage, default='inactive')

    call hist_addfld1d (fname='GRESP_XFER', units='gC/m^2', &
         avgflag='A', long_name='growth respiration transfer', &
         ptr_pft=clm3%g%l%c%p%pcs%gresp_xfer, default='inactive')

    call hist_addfld1d (fname='CPOOL', units='gC/m^2', &
         avgflag='A', long_name='temporary photosynthate C pool', &
         ptr_pft=clm3%g%l%c%p%pcs%cpool)

    call hist_addfld1d (fname='XSMRPOOL', units='gC/m^2', &
         avgflag='A', long_name='temporary photosynthate C pool', &
         ptr_pft=clm3%g%l%c%p%pcs%xsmrpool)

    call hist_addfld1d (fname='PFT_CTRUNC', units='gC/m^2', &
         avgflag='A', long_name='pft-level sink for C truncation', &
         ptr_pft=clm3%g%l%c%p%pcs%pft_ctrunc)

    call hist_addfld1d (fname='DISPVEGC', units='gC/m^2', &
         avgflag='A', long_name='displayed veg carbon, excluding storage and cpool', &
         ptr_pft=clm3%g%l%c%p%pcs%dispvegc)

    call hist_addfld1d (fname='STORVEGC', units='gC/m^2', &
         avgflag='A', long_name='stored vegetation carbon, excluding cpool', &
         ptr_pft=clm3%g%l%c%p%pcs%storvegc)

    call hist_addfld1d (fname='TOTVEGC', units='gC/m^2', &
         avgflag='A', long_name='total vegetation carbon, excluding cpool', &
         ptr_pft=clm3%g%l%c%p%pcs%totvegc)

    call hist_addfld1d (fname='TOTPFTC', units='gC/m^2', &
         avgflag='A', long_name='total pft-level carbon, including cpool', &
         ptr_pft=clm3%g%l%c%p%pcs%totpftc)

#if (defined C13)
    !-------------------------------
    ! C13 state variables - native to PFT 
    !-------------------------------
    
    call hist_addfld1d (fname='C13_LEAFC', units='gC13/m^2', &
         avgflag='A', long_name='C13 leaf C', &
         ptr_pft=clm3%g%l%c%p%pc13s%leafc)

    call hist_addfld1d (fname='C13_LEAFC_STORAGE', units='gC13/m^2', &
         avgflag='A', long_name='C13 leaf C storage', &
         ptr_pft=clm3%g%l%c%p%pc13s%leafc_storage, default='inactive')

    call hist_addfld1d (fname='C13_LEAFC_XFER', units='gC13/m^2', &
         avgflag='A', long_name='C13 leaf C transfer', &
         ptr_pft=clm3%g%l%c%p%pc13s%leafc_xfer, default='inactive')

    call hist_addfld1d (fname='C13_FROOTC', units='gC13/m^2', &
         avgflag='A', long_name='C13 fine root C', &
         ptr_pft=clm3%g%l%c%p%pc13s%frootc)

    call hist_addfld1d (fname='C13_FROOTC_STORAGE', units='gC13/m^2', &
         avgflag='A', long_name='C13 fine root C storage', &
         ptr_pft=clm3%g%l%c%p%pc13s%frootc_storage, default='inactive')

    call hist_addfld1d (fname='C13_FROOTC_XFER', units='gC13/m^2', &
         avgflag='A', long_name='C13 fine root C transfer', &
         ptr_pft=clm3%g%l%c%p%pc13s%frootc_xfer, default='inactive')

    call hist_addfld1d (fname='C13_LIVESTEMC', units='gC13/m^2', &
         avgflag='A', long_name='C13 live stem C', &
         ptr_pft=clm3%g%l%c%p%pc13s%livestemc)

    call hist_addfld1d (fname='C13_LIVESTEMC_STORAGE', units='gC13/m^2', &
         avgflag='A', long_name='C13 live stem C storage', &
         ptr_pft=clm3%g%l%c%p%pc13s%livestemc_storage, default='inactive')

    call hist_addfld1d (fname='C13_LIVESTEMC_XFER', units='gC13/m^2', &
         avgflag='A', long_name='C13 live stem C transfer', &
         ptr_pft=clm3%g%l%c%p%pc13s%livestemc_xfer, default='inactive')

    call hist_addfld1d (fname='C13_DEADSTEMC', units='gC13/m^2', &
         avgflag='A', long_name='C13 dead stem C', &
         ptr_pft=clm3%g%l%c%p%pc13s%deadstemc)

    call hist_addfld1d (fname='C13_DEADSTEMC_STORAGE', units='gC13/m^2', &
         avgflag='A', long_name='C13 dead stem C storage', &
         ptr_pft=clm3%g%l%c%p%pc13s%deadstemc_storage, default='inactive')

    call hist_addfld1d (fname='C13_DEADSTEMC_XFER', units='gC13/m^2', &
         avgflag='A', long_name='C13 dead stem C transfer', &
         ptr_pft=clm3%g%l%c%p%pc13s%deadstemc_xfer, default='inactive')

    call hist_addfld1d (fname='C13_LIVECROOTC', units='gC13/m^2', &
         avgflag='A', long_name='C13 live coarse root C', &
         ptr_pft=clm3%g%l%c%p%pc13s%livecrootc)

    call hist_addfld1d (fname='C13_LIVECROOTC_STORAGE', units='gC13/m^2', &
         avgflag='A', long_name='C13 live coarse root C storage', &
         ptr_pft=clm3%g%l%c%p%pc13s%livecrootc_storage, default='inactive')

    call hist_addfld1d (fname='C13_LIVECROOTC_XFER', units='gC13/m^2', &
         avgflag='A', long_name='C13 live coarse root C transfer', &
         ptr_pft=clm3%g%l%c%p%pc13s%livecrootc_xfer, default='inactive')

    call hist_addfld1d (fname='C13_DEADCROOTC', units='gC13/m^2', &
         avgflag='A', long_name='C13 dead coarse root C', &
         ptr_pft=clm3%g%l%c%p%pc13s%deadcrootc)

    call hist_addfld1d (fname='C13_DEADCROOTC_STORAGE', units='gC13/m^2', &
         avgflag='A', long_name='C13 dead coarse root C storage', &
         ptr_pft=clm3%g%l%c%p%pc13s%deadcrootc_storage,  default='inactive')

    call hist_addfld1d (fname='C13_DEADCROOTC_XFER', units='gC13/m^2', &
         avgflag='A', long_name='C13 dead coarse root C transfer', &
         ptr_pft=clm3%g%l%c%p%pc13s%deadcrootc_xfer, default='inactive')

    call hist_addfld1d (fname='C13_GRESP_STORAGE', units='gC13/m^2', &
         avgflag='A', long_name='C13 growth respiration storage', &
         ptr_pft=clm3%g%l%c%p%pc13s%gresp_storage, default='inactive')

    call hist_addfld1d (fname='C13_GRESP_XFER', units='gC13/m^2', &
         avgflag='A', long_name='C13 growth respiration transfer', &
         ptr_pft=clm3%g%l%c%p%pc13s%gresp_xfer, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL', units='gC13/m^2', &
         avgflag='A', long_name='C13 temporary photosynthate C pool', &
         ptr_pft=clm3%g%l%c%p%pc13s%cpool)

    call hist_addfld1d (fname='C13_XSMRPOOL', units='gC13/m^2', &
         avgflag='A', long_name='C13 temporary photosynthate C pool', &
         ptr_pft=clm3%g%l%c%p%pc13s%xsmrpool)

    call hist_addfld1d (fname='C13_PFT_CTRUNC', units='gC13/m^2', &
         avgflag='A', long_name='C13 pft-level sink for C truncation', &
         ptr_pft=clm3%g%l%c%p%pc13s%pft_ctrunc)

    call hist_addfld1d (fname='C13_DISPVEGC', units='gC13/m^2', &
         avgflag='A', long_name='C13 displayed veg carbon, excluding storage and cpool', &
         ptr_pft=clm3%g%l%c%p%pc13s%dispvegc)

    call hist_addfld1d (fname='C13_STORVEGC', units='gC13/m^2', &
         avgflag='A', long_name='C13 stored vegetation carbon, excluding cpool', &
         ptr_pft=clm3%g%l%c%p%pc13s%storvegc)

    call hist_addfld1d (fname='C13_TOTVEGC', units='gC13/m^2', &
         avgflag='A', long_name='C13 total vegetation carbon, excluding cpool', &
         ptr_pft=clm3%g%l%c%p%pc13s%totvegc)

    call hist_addfld1d (fname='C13_TOTPFTC', units='gC13/m^2', &
         avgflag='A', long_name='C13 total pft-level carbon, including cpool', &
         ptr_pft=clm3%g%l%c%p%pc13s%totpftc)
#endif

    !-------------------------------
    ! C state variables - native to column
    !-------------------------------
     ! add history fields for all CLAMP CN variables
     call hist_addfld1d (fname='SOILC', units='gC/m^2', &
          avgflag='A', long_name='soil C', &
          ptr_col=clm3%g%l%c%ccs%totsomc)

     call hist_addfld1d (fname='LITTERC', units='gC/m^2', &
          avgflag='A', long_name='litter C', &
          ptr_col=clm3%g%l%c%ccs%totlitc)
    
    call hist_addfld1d (fname='CWDC', units='gC/m^2', &
         avgflag='A', long_name='coarse woody debris C', &
         ptr_col=clm3%g%l%c%ccs%cwdc)

    call hist_addfld1d (fname='LITR1C', units='gC/m^2', &
         avgflag='A', long_name='litter labile C', &
         ptr_col=clm3%g%l%c%ccs%litr1c)

    call hist_addfld1d (fname='LITR2C', units='gC/m^2', &
         avgflag='A', long_name='litter cellulose C', &
         ptr_col=clm3%g%l%c%ccs%litr2c)

    call hist_addfld1d (fname='LITR3C', units='gC/m^2', &
         avgflag='A', long_name='litter lignin C', &
         ptr_col=clm3%g%l%c%ccs%litr3c)

    call hist_addfld1d (fname='SOIL1C', units='gC/m^2', &
         avgflag='A', long_name='soil organic matter C (fast pool)', &
         ptr_col=clm3%g%l%c%ccs%soil1c)

    call hist_addfld1d (fname='SOIL2C', units='gC/m^2', &
         avgflag='A', long_name='soil organic matter C (medium pool)', &
         ptr_col=clm3%g%l%c%ccs%soil2c)

    call hist_addfld1d (fname='SOIL3C', units='gC/m^2', &
         avgflag='A', long_name='soil organic matter C (slow pool)', &
         ptr_col=clm3%g%l%c%ccs%soil3c)
    
    call hist_addfld1d (fname='SOIL4C', units='gC/m^2', &
         avgflag='A', long_name='soil organic matter C (slowest pool)', &
         ptr_col=clm3%g%l%c%ccs%soil4c)
    
    call hist_addfld1d (fname='SEEDC', units='gC/m^2', &
         avgflag='A', long_name='pool for seeding new PFTs', &
         ptr_col=clm3%g%l%c%ccs%seedc)
    
    call hist_addfld1d (fname='COL_CTRUNC', units='gC/m^2', &
         avgflag='A', long_name='column-level sink for C truncation', &
         ptr_col=clm3%g%l%c%ccs%col_ctrunc)
    
    call hist_addfld1d (fname='TOTLITC', units='gC/m^2', &
         avgflag='A', long_name='total litter carbon', &
         ptr_col=clm3%g%l%c%ccs%totlitc)
    
    call hist_addfld1d (fname='TOTSOMC', units='gC/m^2', &
         avgflag='A', long_name='total soil organic matter carbon', &
         ptr_col=clm3%g%l%c%ccs%totsomc)
    
    call hist_addfld1d (fname='TOTECOSYSC', units='gC/m^2', &
         avgflag='A', long_name='total ecosystem carbon, incl veg but excl cpool', &
         ptr_col=clm3%g%l%c%ccs%totecosysc)
    
    call hist_addfld1d (fname='TOTCOLC', units='gC/m^2', &
         avgflag='A', long_name='total column carbon, incl veg and cpool', &
         ptr_col=clm3%g%l%c%ccs%totcolc)

    call hist_addfld1d (fname='PROD10C', units='gC/m^2', &
         avgflag='A', long_name='10-yr wood product C', &
         ptr_col=clm3%g%l%c%ccs%prod10c)

    call hist_addfld1d (fname='PROD100C', units='gC/m^2', &
         avgflag='A', long_name='100-yr wood product C', &
         ptr_col=clm3%g%l%c%ccs%prod100c)

    call hist_addfld1d (fname='TOTPRODC', units='gC/m^2', &
         avgflag='A', long_name='total wood product C', &
         ptr_col=clm3%g%l%c%ccs%totprodc)

    
#if (defined C13)
    !-------------------------------
    ! C13 state variables - native to column
    !-------------------------------
    
    call hist_addfld1d (fname='C13_CWDC', units='gC13/m^2', &
         avgflag='A', long_name='C13 coarse woody debris C', &
         ptr_col=clm3%g%l%c%cc13s%cwdc)

    call hist_addfld1d (fname='C13_LITR1C', units='gC13/m^2', &
         avgflag='A', long_name='C13 litter labile C', &
         ptr_col=clm3%g%l%c%cc13s%litr1c)

    call hist_addfld1d (fname='C13_LITR2C', units='gC13/m^2', &
         avgflag='A', long_name='C13 litter cellulose C', &
         ptr_col=clm3%g%l%c%cc13s%litr2c)

    call hist_addfld1d (fname='C13_LITR3C', units='gC13/m^2', &
         avgflag='A', long_name='C13 litter lignin C', &
         ptr_col=clm3%g%l%c%cc13s%litr3c)

    call hist_addfld1d (fname='C13_SOIL1C', units='gC13/m^2', &
         avgflag='A', long_name='C13 soil organic matter C (fast pool)', &
         ptr_col=clm3%g%l%c%cc13s%soil1c)

    call hist_addfld1d (fname='C13_SOIL2C', units='gC13/m^2', &
         avgflag='A', long_name='C13 soil organic matter C (medium pool)', &
         ptr_col=clm3%g%l%c%cc13s%soil2c)

    call hist_addfld1d (fname='C13_SOIL3C', units='gC13/m^2', &
         avgflag='A', long_name='C13 soil organic matter C (slow pool)', &
         ptr_col=clm3%g%l%c%cc13s%soil3c)
    
    call hist_addfld1d (fname='C13_SOIL4C', units='gC13/m^2', &
         avgflag='A', long_name='C13 soil organic matter C (slowest pool)', &
         ptr_col=clm3%g%l%c%cc13s%soil4c)
    
    call hist_addfld1d (fname='C13_SEEDC', units='gC13/m^2', &
         avgflag='A', long_name='C13 pool for seeding new PFTs', &
         ptr_col=clm3%g%l%c%ccs%seedc)
    
    call hist_addfld1d (fname='C13_COL_CTRUNC', units='gC13/m^2', &
         avgflag='A', long_name='C13 column-level sink for C truncation', &
         ptr_col=clm3%g%l%c%cc13s%col_ctrunc)
    
    call hist_addfld1d (fname='C13_TOTLITC', units='gC13/m^2', &
         avgflag='A', long_name='C13 total litter carbon', &
         ptr_col=clm3%g%l%c%cc13s%totlitc)
    
    call hist_addfld1d (fname='C13_TOTSOMC', units='gC13/m^2', &
         avgflag='A', long_name='C13 total soil organic matter carbon', &
         ptr_col=clm3%g%l%c%cc13s%totsomc)
    
    call hist_addfld1d (fname='C13_TOTECOSYSC', units='gC13/m^2', &
         avgflag='A', long_name='C13 total ecosystem carbon, incl veg but excl cpool', &
         ptr_col=clm3%g%l%c%cc13s%totecosysc)
    
    call hist_addfld1d (fname='C13_TOTCOLC', units='gC13/m^2', &
         avgflag='A', long_name='C13 total column carbon, incl veg and cpool', &
         ptr_col=clm3%g%l%c%cc13s%totcolc)

    call hist_addfld1d (fname='C13_PROD10C', units='gC13/m^2', &
         avgflag='A', long_name='C13 10-yr wood product C', &
         ptr_col=clm3%g%l%c%cc13s%prod10c)

    call hist_addfld1d (fname='C13_PROD100C', units='gC13/m^2', &
         avgflag='A', long_name='C13 100-yr wood product C', &
         ptr_col=clm3%g%l%c%cc13s%prod100c)

    call hist_addfld1d (fname='C13_TOTPRODC', units='gC13/m^2', &
         avgflag='A', long_name='C13 total wood product C', &
         ptr_col=clm3%g%l%c%cc13s%totprodc)
#endif

    !-------------------------------
    ! N state variables - native to PFT
    !-------------------------------

    call hist_addfld1d (fname='LEAFN', units='gN/m^2', &
         avgflag='A', long_name='leaf N', &
         ptr_pft=clm3%g%l%c%p%pns%leafn)

    call hist_addfld1d (fname='LEAFN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='leaf N storage', &
         ptr_pft=clm3%g%l%c%p%pns%leafn_storage, default='inactive')

    call hist_addfld1d (fname='LEAFN_XFER', units='gN/m^2', &
         avgflag='A', long_name='leaf N transfer', &
         ptr_pft=clm3%g%l%c%p%pns%leafn_xfer, default='inactive')

    call hist_addfld1d (fname='FROOTN', units='gN/m^2', &
         avgflag='A', long_name='fine root N', &
         ptr_pft=clm3%g%l%c%p%pns%frootn)

    call hist_addfld1d (fname='FROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='fine root N storage', &
         ptr_pft=clm3%g%l%c%p%pns%frootn_storage, default='inactive')

    call hist_addfld1d (fname='FROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='fine root N transfer', &
         ptr_pft=clm3%g%l%c%p%pns%frootn_xfer, default='inactive')

    call hist_addfld1d (fname='LIVESTEMN', units='gN/m^2', &
         avgflag='A', long_name='live stem N', &
         ptr_pft=clm3%g%l%c%p%pns%livestemn)

    call hist_addfld1d (fname='LIVESTEMN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='live stem N storage', &
         ptr_pft=clm3%g%l%c%p%pns%livestemn_storage, default='inactive')

    call hist_addfld1d (fname='LIVESTEMN_XFER', units='gN/m^2', &
         avgflag='A', long_name='live stem N transfer', &
         ptr_pft=clm3%g%l%c%p%pns%livestemn_xfer, default='inactive')

    call hist_addfld1d (fname='DEADSTEMN', units='gN/m^2', &
         avgflag='A', long_name='dead stem N', &
         ptr_pft=clm3%g%l%c%p%pns%deadstemn)

    call hist_addfld1d (fname='DEADSTEMN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='dead stem N storage', &
         ptr_pft=clm3%g%l%c%p%pns%deadstemn_storage, default='inactive')

    call hist_addfld1d (fname='DEADSTEMN_XFER', units='gN/m^2', &
         avgflag='A', long_name='dead stem N transfer', &
         ptr_pft=clm3%g%l%c%p%pns%deadstemn_xfer, default='inactive')

    call hist_addfld1d (fname='LIVECROOTN', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N', &
         ptr_pft=clm3%g%l%c%p%pns%livecrootn)

    call hist_addfld1d (fname='LIVECROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N storage', &
         ptr_pft=clm3%g%l%c%p%pns%livecrootn_storage, default='inactive')

    call hist_addfld1d (fname='LIVECROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N transfer', &
         ptr_pft=clm3%g%l%c%p%pns%livecrootn_xfer, default='inactive')

    call hist_addfld1d (fname='DEADCROOTN', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N', &
         ptr_pft=clm3%g%l%c%p%pns%deadcrootn)

    call hist_addfld1d (fname='DEADCROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N storage', &
         ptr_pft=clm3%g%l%c%p%pns%deadcrootn_storage, default='inactive')

    call hist_addfld1d (fname='DEADCROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N transfer', &
         ptr_pft=clm3%g%l%c%p%pns%deadcrootn_xfer, default='inactive')

    call hist_addfld1d (fname='RETRANSN', units='gN/m^2', &
         avgflag='A', long_name='plant pool of retranslocated N', &
         ptr_pft=clm3%g%l%c%p%pns%retransn)

    call hist_addfld1d (fname='NPOOL', units='gN/m^2', &
         avgflag='A', long_name='temporary plant N pool', &
         ptr_pft=clm3%g%l%c%p%pns%npool, default='inactive')

    call hist_addfld1d (fname='PFT_NTRUNC', units='gN/m^2', &
         avgflag='A', long_name='pft-level sink for N truncation', &
         ptr_pft=clm3%g%l%c%p%pns%pft_ntrunc)

    call hist_addfld1d (fname='DISPVEGN', units='gN/m^2', &
         avgflag='A', long_name='displayed vegetation nitrogen', &
         ptr_pft=clm3%g%l%c%p%pns%dispvegn)

    call hist_addfld1d (fname='STORVEGN', units='gN/m^2', &
         avgflag='A', long_name='stored vegetation nitrogen', &
         ptr_pft=clm3%g%l%c%p%pns%storvegn)

    call hist_addfld1d (fname='TOTVEGN', units='gN/m^2', &
         avgflag='A', long_name='total vegetation nitrogen', &
         ptr_pft=clm3%g%l%c%p%pns%totvegn)

    call hist_addfld1d (fname='TOTPFTN', units='gN/m^2', &
         avgflag='A', long_name='total PFT-level nitrogen', &
         ptr_pft=clm3%g%l%c%p%pns%totpftn)

    !-------------------------------
    ! N state variables - native to column
    !-------------------------------

    call hist_addfld1d (fname='CWDN', units='gN/m^2', &
         avgflag='A', long_name='coarse woody debris N', &
         ptr_col=clm3%g%l%c%cns%cwdn)

    call hist_addfld1d (fname='LITR1N', units='gN/m^2', &
         avgflag='A', long_name='litter labile N', &
         ptr_col=clm3%g%l%c%cns%litr1n)

    call hist_addfld1d (fname='LITR2N', units='gN/m^2', &
         avgflag='A', long_name='litter cellulose N', &
         ptr_col=clm3%g%l%c%cns%litr2n)

    call hist_addfld1d (fname='LITR3N', units='gN/m^2', &
         avgflag='A', long_name='litter lignin N', &
         ptr_col=clm3%g%l%c%cns%litr3n)

    call hist_addfld1d (fname='SOIL1N', units='gN/m^2', &
         avgflag='A', long_name='soil organic matter N (fast pool)', &
         ptr_col=clm3%g%l%c%cns%soil1n)

    call hist_addfld1d (fname='SOIL2N', units='gN/m^2', &
         avgflag='A', long_name='soil organic matter N (medium pool)', &
         ptr_col=clm3%g%l%c%cns%soil2n)

    call hist_addfld1d (fname='SOIL3N', units='gN/m^2', &
         avgflag='A', long_name='soil orgainc matter N (slow pool)', &
         ptr_col=clm3%g%l%c%cns%soil3n)

    call hist_addfld1d (fname='SOIL4N', units='gN/m^2', &
         avgflag='A', long_name='soil orgainc matter N (slowest pool)', &
         ptr_col=clm3%g%l%c%cns%soil4n)

    call hist_addfld1d (fname='SMINN', units='gN/m^2', &
         avgflag='A', long_name='soil mineral N', &
         ptr_col=clm3%g%l%c%cns%sminn)

    call hist_addfld1d (fname='COL_NTRUNC', units='gN/m^2', &
         avgflag='A', long_name='column-level sink for N truncation', &
         ptr_col=clm3%g%l%c%cns%col_ntrunc)

    call hist_addfld1d (fname='TOTLITN', units='gN/m^2', &
         avgflag='A', long_name='total litter N', &
         ptr_col=clm3%g%l%c%cns%totlitn)

    call hist_addfld1d (fname='TOTSOMN', units='gN/m^2', &
         avgflag='A', long_name='total soil organic matter N', &
         ptr_col=clm3%g%l%c%cns%totsomn)

    call hist_addfld1d (fname='TOTECOSYSN', units='gN/m^2', &
         avgflag='A', long_name='total ecosystem N', &
         ptr_col=clm3%g%l%c%cns%totecosysn)

    call hist_addfld1d (fname='TOTCOLN', units='gN/m^2', &
         avgflag='A', long_name='total column-level N', &
         ptr_col=clm3%g%l%c%cns%totcoln)

    call hist_addfld1d (fname='SEEDN', units='gN/m^2', &
         avgflag='A', long_name='pool for seeding new PFTs ', &
         ptr_col=clm3%g%l%c%cns%seedn)

    call hist_addfld1d (fname='PROD10N', units='gN/m^2', &
         avgflag='A', long_name='10-yr wood product N', &
         ptr_col=clm3%g%l%c%cns%prod10n)

    call hist_addfld1d (fname='PROD100N', units='gN/m^2', &
         avgflag='A', long_name='100-yr wood product N', &
         ptr_col=clm3%g%l%c%cns%prod100n)

    call hist_addfld1d (fname='TOTPRODN', units='gN/m^2', &
         avgflag='A', long_name='total wood product N', &
         ptr_col=clm3%g%l%c%cns%totprodn)

    !-------------------------------
    ! C flux variables - native to PFT
    !-------------------------------

     ! add history fields for all CLAMP CN variables

     call hist_addfld1d (fname='WOODC_ALLOC', units='gC/m^2/s', &
          avgflag='A', long_name='wood C allocation', &
          ptr_pft=clm3%g%l%c%p%pcf%woodc_alloc)

     call hist_addfld1d (fname='WOODC_LOSS', units='gC/m^2/s', &
          avgflag='A', long_name='wood C loss', &
          ptr_pft=clm3%g%l%c%p%pcf%woodc_loss)

     call hist_addfld1d (fname='LEAFC_LOSS', units='gC/m^2/s', &
          avgflag='A', long_name='leaf C loss', &
          ptr_pft=clm3%g%l%c%p%pcf%leafc_loss)

     call hist_addfld1d (fname='LEAFC_ALLOC', units='gC/m^2/s', &
          avgflag='A', long_name='leaf C allocation', &
          ptr_pft=clm3%g%l%c%p%pcf%leafc_alloc)

     call hist_addfld1d (fname='FROOTC_LOSS', units='gC/m^2/s', &
          avgflag='A', long_name='fine root C loss', &
          ptr_pft=clm3%g%l%c%p%pcf%frootc_loss)

     call hist_addfld1d (fname='FROOTC_ALLOC', units='gC/m^2/s', &
          avgflag='A', long_name='fine root C allocation', &
          ptr_pft=clm3%g%l%c%p%pcf%frootc_alloc)

    call hist_addfld1d (fname='PSNSUN', units='umolCO2/m^2/s', &
         avgflag='A', long_name='sunlit leaf photosynthesis', &
         ptr_pft=clm3%g%l%c%p%pcf%psnsun)

    call hist_addfld1d (fname='PSNSHA', units='umolCO2/m^2/s', &
         avgflag='A', long_name='shaded leaf photosynthesis', &
         ptr_pft=clm3%g%l%c%p%pcf%psnsha)

    call hist_addfld1d (fname='M_LEAFC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_leafc_to_litter, default='inactive')

    call hist_addfld1d (fname='M_FROOTC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_frootc_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LEAFC_STORAGE_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_leafc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_FROOTC_STORAGE_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_frootc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMC_STORAGE_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livestemc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMC_STORAGE_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadstemc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTC_STORAGE_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livecrootc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTC_STORAGE_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadcrootc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LEAFC_XFER_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_leafc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_FROOTC_XFER_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_frootc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMC_XFER_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livestemc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMC_XFER_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadstemc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTC_XFER_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livecrootc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTC_XFER_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadcrootc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livestemc_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadstemc_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livecrootc_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadcrootc_to_litter, default='inactive')

    call hist_addfld1d (fname='M_GRESP_STORAGE_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='growth respiration storage mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_gresp_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_GRESP_XFER_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='growth respiration transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pcf%m_gresp_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LEAFC_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_leafc_to_fire, default='inactive')

    call hist_addfld1d (fname='M_FROOTC_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_frootc_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LEAFC_STORAGE_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_leafc_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_FROOTC_STORAGE_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_frootc_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMC_STORAGE_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livestemc_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMC_STORAGE_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadstemc_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTC_STORAGE_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livecrootc_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTC_STORAGE_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadcrootc_storage_to_fire,  default='inactive')

    call hist_addfld1d (fname='M_LEAFC_XFER_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_leafc_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_FROOTC_XFER_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_frootc_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMC_XFER_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livestemc_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMC_XFER_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadstemc_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTC_XFER_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livecrootc_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTC_XFER_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadcrootc_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMC_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livestemc_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMC_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadstemc_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMC_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadstemc_to_litter_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTC_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_livecrootc_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTC_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadcrootc_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTC_TO_LITTER_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pcf%m_deadcrootc_to_litter_fire, default='inactive')

    call hist_addfld1d (fname='M_GRESP_STORAGE_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='growth respiration storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_gresp_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_GRESP_XFER_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='growth respiration transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%m_gresp_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='LEAFC_XFER_TO_LEAFC', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%leafc_xfer_to_leafc, default='inactive')

    call hist_addfld1d (fname='FROOTC_XFER_TO_FROOTC', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%frootc_xfer_to_frootc, default='inactive')

    call hist_addfld1d (fname='LIVESTEMC_XFER_TO_LIVESTEMC', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%livestemc_xfer_to_livestemc, default='inactive')

    call hist_addfld1d (fname='DEADSTEMC_XFER_TO_DEADSTEMC', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%deadstemc_xfer_to_deadstemc, default='inactive')

    call hist_addfld1d (fname='LIVECROOTC_XFER_TO_LIVECROOTC', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%livecrootc_xfer_to_livecrootc, default='inactive')

    call hist_addfld1d (fname='DEADCROOTC_XFER_TO_DEADCROOTC', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%deadcrootc_xfer_to_deadcrootc, default='inactive')

    call hist_addfld1d (fname='LEAFC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C litterfall', &
         ptr_pft=clm3%g%l%c%p%pcf%leafc_to_litter, default='inactive')

    call hist_addfld1d (fname='FROOTC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C litterfall', &
         ptr_pft=clm3%g%l%c%p%pcf%frootc_to_litter, default='inactive')

    call hist_addfld1d (fname='LEAF_MR', units='gC/m^2/s', &
         avgflag='A', long_name='leaf maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%leaf_mr, default='inactive')

    call hist_addfld1d (fname='FROOT_MR', units='gC/m^2/s', &
         avgflag='A', long_name='fine root maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%froot_mr, default='inactive')

    call hist_addfld1d (fname='LIVESTEM_MR', units='gC/m^2/s', &
         avgflag='A', long_name='live stem maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%livestem_mr, default='inactive')

    call hist_addfld1d (fname='LIVECROOT_MR', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%livecroot_mr, default='inactive')

    call hist_addfld1d (fname='PSNSUN_TO_CPOOL', units='gC/m^2/s', &
         avgflag='A', long_name='C fixation from sunlit canopy', &
         ptr_pft=clm3%g%l%c%p%pcf%psnsun_to_cpool)

    call hist_addfld1d (fname='PSNSHADE_TO_CPOOL', units='gC/m^2/s', &
         avgflag='A', long_name='C fixation from shaded canopy', &
         ptr_pft=clm3%g%l%c%p%pcf%psnshade_to_cpool)

    call hist_addfld1d (fname='CPOOL_TO_LEAFC', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to leaf C', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_leafc, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_LEAFC_STORAGE', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to leaf C storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_leafc_storage, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_FROOTC', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to fine root C', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_frootc, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_FROOTC_STORAGE', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to fine root C storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_frootc_storage, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_LIVESTEMC', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to live stem C', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_livestemc, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_LIVESTEMC_STORAGE', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to live stem C storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_livestemc_storage, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_DEADSTEMC', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to dead stem C', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_deadstemc, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_DEADSTEMC_STORAGE', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to dead stem C storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_deadstemc_storage, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_LIVECROOTC', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root C', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_livecrootc, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_LIVECROOTC_STORAGE', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root C storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_livecrootc_storage, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_DEADCROOTC', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root C', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_deadcrootc, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_DEADCROOTC_STORAGE', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root C storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_deadcrootc_storage, default='inactive')

    call hist_addfld1d (fname='CPOOL_TO_GRESP_STORAGE', units='gC/m^2/s', &
         avgflag='A', long_name='allocation to growth respiration storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_to_gresp_storage, default='inactive')

    call hist_addfld1d (fname='CPOOL_LEAF_GR', units='gC/m^2/s', &
         avgflag='A', long_name='leaf growth respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_leaf_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_LEAF_STORAGE_GR', units='gC/m^2/s', &
         avgflag='A', long_name='leaf growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_leaf_storage_gr, default='inactive')

    call hist_addfld1d (fname='TRANSFER_LEAF_GR', units='gC/m^2/s', &
         avgflag='A', long_name='leaf growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%transfer_leaf_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_FROOT_GR', units='gC/m^2/s', &
         avgflag='A', long_name='fine root growth respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_froot_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_FROOT_STORAGE_GR', units='gC/m^2/s', &
         avgflag='A', long_name='fine root  growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_froot_storage_gr, default='inactive')

    call hist_addfld1d (fname='TRANSFER_FROOT_GR', units='gC/m^2/s', &
         avgflag='A', long_name='fine root  growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%transfer_froot_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_LIVESTEM_GR', units='gC/m^2/s', &
         avgflag='A', long_name='live stem growth respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_livestem_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_LIVESTEM_STORAGE_GR', units='gC/m^2/s', &
         avgflag='A', long_name='live stem growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_livestem_storage_gr, default='inactive')

    call hist_addfld1d (fname='TRANSFER_LIVESTEM_GR', units='gC/m^2/s', &
         avgflag='A', long_name='live stem growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%transfer_livestem_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_DEADSTEM_GR', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem growth respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_deadstem_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_DEADSTEM_STORAGE_GR', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_deadstem_storage_gr, default='inactive')

    call hist_addfld1d (fname='TRANSFER_DEADSTEM_GR', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%transfer_deadstem_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_LIVECROOT_GR', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root growth respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_livecroot_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_LIVECROOT_STORAGE_GR', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_livecroot_storage_gr, default='inactive')

    call hist_addfld1d (fname='TRANSFER_LIVECROOT_GR', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%transfer_livecroot_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_DEADCROOT_GR', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root growth respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_deadcroot_gr, default='inactive')

    call hist_addfld1d (fname='CPOOL_DEADCROOT_STORAGE_GR', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pcf%cpool_deadcroot_storage_gr, default='inactive')

    call hist_addfld1d (fname='TRANSFER_DEADCROOT_GR', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pcf%transfer_deadcroot_gr, default='inactive')

    call hist_addfld1d (fname='LEAFC_STORAGE_TO_XFER', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pcf%leafc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='FROOTC_STORAGE_TO_XFER', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pcf%frootc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='LIVESTEMC_STORAGE_TO_XFER', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pcf%livestemc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='DEADSTEMC_STORAGE_TO_XFER', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pcf%deadstemc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='LIVECROOTC_STORAGE_TO_XFER', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pcf%livecrootc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='DEADCROOTC_STORAGE_TO_XFER', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pcf%deadcrootc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='GRESP_STORAGE_TO_XFER', units='gC/m^2/s', &
         avgflag='A', long_name='growth respiration shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pcf%gresp_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='LIVESTEMC_TO_DEADSTEMC', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C turnover', &
         ptr_pft=clm3%g%l%c%p%pcf%livestemc_to_deadstemc, default='inactive')

    call hist_addfld1d (fname='LIVECROOTC_TO_DEADCROOTC', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C turnover', &
         ptr_pft=clm3%g%l%c%p%pcf%livecrootc_to_deadcrootc, default='inactive')

    call hist_addfld1d (fname='GPP', units='gC/m^2/s', &
         avgflag='A', long_name='gross primary production', &
         ptr_pft=clm3%g%l%c%p%pcf%gpp)

    call hist_addfld1d (fname='MR', units='gC/m^2/s', &
         avgflag='A', long_name='maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%mr)

    call hist_addfld1d (fname='CURRENT_GR', units='gC/m^2/s', &
         avgflag='A', long_name='growth resp for new growth displayed in this timestep', &
         ptr_pft=clm3%g%l%c%p%pcf%current_gr, default='inactive')

    call hist_addfld1d (fname='TRANSFER_GR', units='gC/m^2/s', &
         avgflag='A', long_name='growth resp for transfer growth displayed in this timestep', &
         ptr_pft=clm3%g%l%c%p%pcf%transfer_gr, default='inactive')

    call hist_addfld1d (fname='STORAGE_GR', units='gC/m^2/s', &
         avgflag='A', long_name='growth resp for growth sent to storage for later display', &
         ptr_pft=clm3%g%l%c%p%pcf%storage_gr, default='inactive')

    call hist_addfld1d (fname='GR', units='gC/m^2/s', &
         avgflag='A', long_name='total growth respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%gr)

    call hist_addfld1d (fname='AR', units='gC/m^2/s', &
         avgflag='A', long_name='autotrophic respiration (MR + GR)', &
         ptr_pft=clm3%g%l%c%p%pcf%ar)

    call hist_addfld1d (fname='RR', units='gC/m^2/s', &
         avgflag='A', long_name='root respiration (fine root MR + total root GR)', &
         ptr_pft=clm3%g%l%c%p%pcf%rr)

    call hist_addfld1d (fname='NPP', units='gC/m^2/s', &
         avgflag='A', long_name='net primary production', &
         ptr_pft=clm3%g%l%c%p%pcf%npp)

    call hist_addfld1d (fname='AGNPP', units='gC/m^2/s', &
         avgflag='A', long_name='aboveground NPP', &
         ptr_pft=clm3%g%l%c%p%pcf%agnpp)

    call hist_addfld1d (fname='BGNPP', units='gC/m^2/s', &
         avgflag='A', long_name='belowground NPP', &
         ptr_pft=clm3%g%l%c%p%pcf%bgnpp)

    call hist_addfld1d (fname='LITFALL', units='gC/m^2/s', &
         avgflag='A', long_name='litterfall (leaves and fine roots)', &
         ptr_pft=clm3%g%l%c%p%pcf%litfall)

    call hist_addfld1d (fname='VEGFIRE', units='gC/m^2/s', &
         avgflag='A', long_name='pft-level fire loss', &
         ptr_pft=clm3%g%l%c%p%pcf%vegfire, default='inactive')

    call hist_addfld1d (fname='WOOD_HARVESTC', units='gC/m^2/s', &
         avgflag='A', long_name='wood harvest (to product pools)', &
         ptr_pft=clm3%g%l%c%p%pcf%wood_harvestc)

    call hist_addfld1d (fname='PFT_FIRE_CLOSS', units='gC/m^2/s', &
         avgflag='A', long_name='total pft-level fire C loss', &
         ptr_pft=clm3%g%l%c%p%pcf%pft_fire_closs)

#if (defined C13)
    !-------------------------------
    ! C13 flux variables - native to PFT
    !-------------------------------

    call hist_addfld1d (fname='C13_PSNSUN', units='umolCO2/m^2/s', &
         avgflag='A', long_name='C13 sunlit leaf photosynthesis', &
         ptr_pft=clm3%g%l%c%p%pc13f%psnsun)

    call hist_addfld1d (fname='C13_PSNSHA', units='umolCO2/m^2/s', &
         avgflag='A', long_name='C13 shaded leaf photosynthesis', &
         ptr_pft=clm3%g%l%c%p%pc13f%psnsha)

    call hist_addfld1d (fname='C13_M_LEAFC_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf C mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_leafc_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_FROOTC_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root C mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_frootc_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_LEAFC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_leafc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_FROOTC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_frootc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVESTEMC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live stem C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_livestemc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADSTEMC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead stem C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_deadstemc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVECROOTC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live coarse root C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_livecrootc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADCROOTC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead coarse root C storage mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_deadcrootc_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_LEAFC_XFER_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_leafc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_FROOTC_XFER_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_frootc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVESTEMC_XFER_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live stem C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_livestemc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADSTEMC_XFER_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead stem C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_deadstemc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVECROOTC_XFER_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live coarse root C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_livecrootc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADCROOTC_XFER_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead coarse root C transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_deadcrootc_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVESTEMC_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live stem C mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_livestemc_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead stem C mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_deadstemc_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVECROOTC_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live coarse root C mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_livecrootc_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead coarse root C mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_deadcrootc_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_GRESP_STORAGE_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 growth respiration storage mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_gresp_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_GRESP_XFER_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 growth respiration transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_gresp_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_M_LEAFC_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf C fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_leafc_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_FROOTC_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root C fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_frootc_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_LEAFC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_leafc_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_FROOTC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_frootc_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVESTEMC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live stem C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_livestemc_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADSTEMC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead stem C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_deadstemc_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVECROOTC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live coarse root C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_livecrootc_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADCROOTC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead coarse root C storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_deadcrootc_storage_to_fire,  default='inactive')

    call hist_addfld1d (fname='C13_M_LEAFC_XFER_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_leafc_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_FROOTC_XFER_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_frootc_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVESTEMC_XFER_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live stem C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_livestemc_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADSTEMC_XFER_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead stem C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_deadstemc_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVECROOTC_XFER_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live coarse root C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_livecrootc_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADCROOTC_XFER_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead coarse root C transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_deadcrootc_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVESTEMC_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live stem C fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_livestemc_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead stem C fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_deadstemc_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_LITTER_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead stem C fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_deadstemc_to_litter_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVECROOTC_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live coarse root C fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_livecrootc_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead coarse root C fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_deadcrootc_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_LITTER_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead coarse root C fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_deadcrootc_to_litter_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_GRESP_STORAGE_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 growth respiration storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_gresp_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_GRESP_XFER_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 growth respiration transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%m_gresp_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_LEAFC_XFER_TO_LEAFC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%leafc_xfer_to_leafc, default='inactive')

    call hist_addfld1d (fname='C13_FROOTC_XFER_TO_FROOTC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%frootc_xfer_to_frootc, default='inactive')

    call hist_addfld1d (fname='C13_LIVESTEMC_XFER_TO_LIVESTEMC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live stem C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%livestemc_xfer_to_livestemc, default='inactive')

    call hist_addfld1d (fname='C13_DEADSTEMC_XFER_TO_DEADSTEMC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead stem C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%deadstemc_xfer_to_deadstemc, default='inactive')

    call hist_addfld1d (fname='C13_LIVECROOTC_XFER_TO_LIVECROOTC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live coarse root C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%livecrootc_xfer_to_livecrootc, default='inactive')

    call hist_addfld1d (fname='C13_DEADCROOTC_XFER_TO_DEADCROOTC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead coarse root C growth from storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%deadcrootc_xfer_to_deadcrootc, default='inactive')

    call hist_addfld1d (fname='C13_LEAFC_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf C litterfall', &
         ptr_pft=clm3%g%l%c%p%pc13f%leafc_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_FROOTC_TO_LITTER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root C litterfall', &
         ptr_pft=clm3%g%l%c%p%pc13f%frootc_to_litter, default='inactive')

    call hist_addfld1d (fname='C13_LEAF_MR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pc13f%leaf_mr, default='inactive')

    call hist_addfld1d (fname='C13_FROOT_MR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pc13f%froot_mr, default='inactive')

    call hist_addfld1d (fname='C13_LIVESTEM_MR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live stem maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pc13f%livestem_mr, default='inactive')

    call hist_addfld1d (fname='C13_LIVECROOT_MR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live coarse root maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pc13f%livecroot_mr, default='inactive')

    call hist_addfld1d (fname='C13_PSNSUN_TO_CPOOL', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 C fixation from sunlit canopy', &
         ptr_pft=clm3%g%l%c%p%pc13f%psnsun_to_cpool)

    call hist_addfld1d (fname='C13_PSNSHADE_TO_CPOOL', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 C fixation from shaded canopy', &
         ptr_pft=clm3%g%l%c%p%pc13f%psnshade_to_cpool)

    call hist_addfld1d (fname='C13_CPOOL_TO_LEAFC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 allocation to leaf C', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_leafc, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_TO_LEAFC_STORAGE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 allocation to leaf C storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_leafc_storage, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_TO_FROOTC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 allocation to fine root C', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_frootc, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_TO_FROOTC_STORAGE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 allocation to fine root C storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_frootc_storage, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_TO_LIVESTEMC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 allocation to live stem C', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_livestemc, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_TO_LIVESTEMC_STORAGE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 allocation to live stem C storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_livestemc_storage, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_TO_DEADSTEMC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 allocation to dead stem C', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_deadstemc, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_TO_DEADSTEMC_STORAGE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 allocation to dead stem C storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_deadstemc_storage, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_TO_LIVECROOTC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 allocation to live coarse root C', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_livecrootc, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_TO_LIVECROOTC_STORAGE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 allocation to live coarse root C storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_livecrootc_storage, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_TO_DEADCROOTC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 allocation to dead coarse root C', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_deadcrootc, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_TO_DEADCROOTC_STORAGE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 allocation to dead coarse root C storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_deadcrootc_storage, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_TO_GRESP_STORAGE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 allocation to growth respiration storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_to_gresp_storage, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_LEAF_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf growth respiration', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_leaf_gr, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_LEAF_STORAGE_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_leaf_storage_gr, default='inactive')

    call hist_addfld1d (fname='C13_TRANSFER_LEAF_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%transfer_leaf_gr, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_FROOT_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root growth respiration', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_froot_gr, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_FROOT_STORAGE_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root  growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_froot_storage_gr, default='inactive')

    call hist_addfld1d (fname='C13_TRANSFER_FROOT_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root  growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%transfer_froot_gr, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_LIVESTEM_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live stem growth respiration', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_livestem_gr, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_LIVESTEM_STORAGE_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live stem growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_livestem_storage_gr, default='inactive')

    call hist_addfld1d (fname='C13_TRANSFER_LIVESTEM_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live stem growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%transfer_livestem_gr, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_DEADSTEM_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead stem growth respiration', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_deadstem_gr, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_DEADSTEM_STORAGE_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead stem growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_deadstem_storage_gr, default='inactive')

    call hist_addfld1d (fname='C13_TRANSFER_DEADSTEM_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead stem growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%transfer_deadstem_gr, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_LIVECROOT_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live coarse root growth respiration', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_livecroot_gr, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_LIVECROOT_STORAGE_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live coarse root growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_livecroot_storage_gr, default='inactive')

    call hist_addfld1d (fname='C13_TRANSFER_LIVECROOT_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live coarse root growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%transfer_livecroot_gr, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_DEADCROOT_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead coarse root growth respiration', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_deadcroot_gr, default='inactive')

    call hist_addfld1d (fname='C13_CPOOL_DEADCROOT_STORAGE_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead coarse root growth respiration to storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%cpool_deadcroot_storage_gr, default='inactive')

    call hist_addfld1d (fname='C13_TRANSFER_DEADCROOT_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead coarse root growth respiration from storage', &
         ptr_pft=clm3%g%l%c%p%pc13f%transfer_deadcroot_gr, default='inactive')

    call hist_addfld1d (fname='C13_LEAFC_STORAGE_TO_XFER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pc13f%leafc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='C13_FROOTC_STORAGE_TO_XFER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pc13f%frootc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='C13_LIVESTEMC_STORAGE_TO_XFER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live stem C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pc13f%livestemc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='C13_DEADSTEMC_STORAGE_TO_XFER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead stem C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pc13f%deadstemc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='C13_LIVECROOTC_STORAGE_TO_XFER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live coarse root C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pc13f%livecrootc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='C13_DEADCROOTC_STORAGE_TO_XFER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead coarse root C shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pc13f%deadcrootc_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='C13_GRESP_STORAGE_TO_XFER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 growth respiration shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pc13f%gresp_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='C13_LIVESTEMC_TO_DEADSTEMC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live stem C turnover', &
         ptr_pft=clm3%g%l%c%p%pc13f%livestemc_to_deadstemc, default='inactive')

    call hist_addfld1d (fname='C13_LIVECROOTC_TO_DEADCROOTC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live coarse root C turnover', &
         ptr_pft=clm3%g%l%c%p%pc13f%livecrootc_to_deadcrootc, default='inactive')

    call hist_addfld1d (fname='C13_GPP', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 gross primary production', &
         ptr_pft=clm3%g%l%c%p%pc13f%gpp)

    call hist_addfld1d (fname='C13_MR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pc13f%mr)

    call hist_addfld1d (fname='C13_CURRENT_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 growth resp for new growth displayed in this timestep', &
         ptr_pft=clm3%g%l%c%p%pc13f%current_gr, default='inactive')

    call hist_addfld1d (fname='C13_TRANSFER_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 growth resp for transfer growth displayed in this timestep', &
         ptr_pft=clm3%g%l%c%p%pc13f%transfer_gr, default='inactive')

    call hist_addfld1d (fname='C13_STORAGE_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 growth resp for growth sent to storage for later display', &
         ptr_pft=clm3%g%l%c%p%pc13f%storage_gr, default='inactive')

    call hist_addfld1d (fname='C13_GR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 total growth respiration', &
         ptr_pft=clm3%g%l%c%p%pc13f%gr)

    call hist_addfld1d (fname='C13_AR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 autotrophic respiration (MR + GR)', &
         ptr_pft=clm3%g%l%c%p%pc13f%ar)

    call hist_addfld1d (fname='C13_RR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 root respiration (fine root MR + total root GR)', &
         ptr_pft=clm3%g%l%c%p%pc13f%rr)

    call hist_addfld1d (fname='C13_NPP', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 net primary production', &
         ptr_pft=clm3%g%l%c%p%pc13f%npp)

    call hist_addfld1d (fname='C13_AGNPP', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 aboveground NPP', &
         ptr_pft=clm3%g%l%c%p%pc13f%agnpp)

    call hist_addfld1d (fname='C13_BGNPP', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 belowground NPP', &
         ptr_pft=clm3%g%l%c%p%pc13f%bgnpp)

    call hist_addfld1d (fname='C13_LITFALL', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 litterfall (leaves and fine roots)', &
         ptr_pft=clm3%g%l%c%p%pc13f%litfall, default='inactive')

    call hist_addfld1d (fname='C13_VEGFIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 pft-level fire loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%vegfire, default='inactive')

    call hist_addfld1d (fname='C13_PFT_FIRE_CLOSS', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 total pft-level fire C loss', &
         ptr_pft=clm3%g%l%c%p%pc13f%pft_fire_closs)
#endif

    !-------------------------------
    ! C flux variables - native to column 
    !-------------------------------
    ! add history fields for all CLAMP CN variables

    call hist_addfld1d (fname='CWDC_HR', units='gC/m^2/s', &
         avgflag='A', long_name='coarse woody debris C heterotrophic respiration', &
         ptr_col=clm3%g%l%c%ccf%cwdc_hr)

    call hist_addfld1d (fname='CWDC_LOSS', units='gC/m^2/s', &
         avgflag='A', long_name='coarse woody debris C loss', &
         ptr_col=clm3%g%l%c%ccf%cwdc_loss)

    call hist_addfld1d (fname='LITTERC_HR', units='gC/m^2/s', &
         avgflag='A', long_name='litter C heterotrophic respiration', &
         ptr_col=clm3%g%l%c%ccf%lithr)

    call hist_addfld1d (fname='LITTERC_LOSS', units='gC/m^2/s', &
         avgflag='A', long_name='litter C loss', &
         ptr_col=clm3%g%l%c%ccf%litterc_loss)

    call hist_addfld1d (fname='SOILC_HR', units='gC/m^2/s', &
         avgflag='A', long_name='soil C heterotrophic respiration', &
         ptr_col=clm3%g%l%c%ccf%somhr)

    call hist_addfld1d (fname='SOILC_LOSS', units='gC/m^2/s', &
         avgflag='A', long_name='soil C loss', &
         ptr_col=clm3%g%l%c%ccf%somhr)

    call hist_addfld1d (fname='M_LEAFC_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%m_leafc_to_litr1c, default='inactive')

    call hist_addfld1d (fname='M_LEAFC_TO_LITR2C', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C mortality to litter 2 C', &
         ptr_col=clm3%g%l%c%ccf%m_leafc_to_litr2c, default='inactive')

    call hist_addfld1d (fname='M_LEAFC_TO_LITR3C', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C mortality to litter 3 C', &
         ptr_col=clm3%g%l%c%ccf%m_leafc_to_litr3c, default='inactive')

    call hist_addfld1d (fname='M_FROOTC_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%m_frootc_to_litr1c, default='inactive')

    call hist_addfld1d (fname='M_FROOTC_TO_LITR2C', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C mortality to litter 2 C', &
         ptr_col=clm3%g%l%c%ccf%m_frootc_to_litr2c, default='inactive')

    call hist_addfld1d (fname='M_FROOTC_TO_LITR3C', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C mortality to litter 3 C', &
         ptr_col=clm3%g%l%c%ccf%m_frootc_to_litr3c, default='inactive')

    call hist_addfld1d (fname='M_LEAFC_STORAGE_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C storage mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%m_leafc_storage_to_litr1c, default='inactive')

    call hist_addfld1d (fname='M_FROOTC_STORAGE_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C storage mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%m_frootc_storage_to_litr1c, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMC_STORAGE_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C storage mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%m_livestemc_storage_to_litr1c, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMC_STORAGE_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C storage mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%m_deadstemc_storage_to_litr1c, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTC_STORAGE_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C storage mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%m_livecrootc_storage_to_litr1c, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTC_STORAGE_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C storage mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%m_deadcrootc_storage_to_litr1c, default='inactive')

    call hist_addfld1d (fname='M_LEAFC_XFER_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C transfer mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%m_leafc_xfer_to_litr1c, default='inactive')

    call hist_addfld1d (fname='M_FROOTC_XFER_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C transfer mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%m_frootc_xfer_to_litr1c, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMC_XFER_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C transfer mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%m_livestemc_xfer_to_litr1c, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMC_XFER_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C transfer mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%m_deadstemc_xfer_to_litr1c, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTC_XFER_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C transfer mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%m_livecrootc_xfer_to_litr1c, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTC_XFER_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C transfer mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%m_deadcrootc_xfer_to_litr1c, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMC_TO_CWDC', units='gC/m^2/s', &
         avgflag='A', long_name='live stem C mortality to coarse woody debris C', &
         ptr_col=clm3%g%l%c%ccf%m_livestemc_to_cwdc, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMC_TO_CWDC', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C mortality to coarse woody debris C', &
         ptr_col=clm3%g%l%c%ccf%m_deadstemc_to_cwdc, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTC_TO_CWDC', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root C mortality to coarse woody debris C', &
         ptr_col=clm3%g%l%c%ccf%m_livecrootc_to_cwdc, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTC_TO_CWDC', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C mortality to coarse woody debris C', &
         ptr_col=clm3%g%l%c%ccf%m_deadcrootc_to_cwdc, default='inactive')

    call hist_addfld1d (fname='M_GRESP_STORAGE_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='growth respiration storage mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%m_gresp_storage_to_litr1c, default='inactive')

    call hist_addfld1d (fname='M_GRESP_XFER_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='growth respiration transfer mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%m_gresp_xfer_to_litr1c, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMC_TO_CWDC_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead stem C to coarse woody debris C by fire', &
         ptr_col=clm3%g%l%c%ccf%m_deadstemc_to_cwdc_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTC_TO_CWDC_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root C to to woody debris C by fire', &
         ptr_col=clm3%g%l%c%ccf%m_deadcrootc_to_cwdc_fire, default='inactive')

    call hist_addfld1d (fname='M_LITR1C_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='litter 1 C fire loss', &
         ptr_col=clm3%g%l%c%ccf%m_litr1c_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LITR2C_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='litter 2 C fire loss', &
         ptr_col=clm3%g%l%c%ccf%m_litr2c_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LITR3C_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='litter 3 C fire loss', &
         ptr_col=clm3%g%l%c%ccf%m_litr3c_to_fire, default='inactive')

    call hist_addfld1d (fname='M_CWDC_TO_FIRE', units='gC/m^2/s', &
         avgflag='A', long_name='coarse woody debris C fire loss', &
         ptr_col=clm3%g%l%c%ccf%m_cwdc_to_fire, default='inactive')

    call hist_addfld1d (fname='LEAFC_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C litterfall to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%leafc_to_litr1c, default='inactive')

    call hist_addfld1d (fname='LEAFC_TO_LITR2C', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C litterfall to litter 2 C', &
         ptr_col=clm3%g%l%c%ccf%leafc_to_litr2c, default='inactive')

    call hist_addfld1d (fname='LEAFC_TO_LITR3C', units='gC/m^2/s', &
         avgflag='A', long_name='leaf C litterfall to litter 3 C', &
         ptr_col=clm3%g%l%c%ccf%leafc_to_litr3c, default='inactive')

    call hist_addfld1d (fname='FROOTC_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C litterfall to litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%frootc_to_litr1c, default='inactive')

    call hist_addfld1d (fname='FROOTC_TO_LITR2C', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C litterfall to litter 2 C', &
         ptr_col=clm3%g%l%c%ccf%frootc_to_litr2c, default='inactive')

    call hist_addfld1d (fname='FROOTC_TO_LITR3C', units='gC/m^2/s', &
         avgflag='A', long_name='fine root C litterfall to litter 3 C', &
         ptr_col=clm3%g%l%c%ccf%frootc_to_litr3c, default='inactive')

    call hist_addfld1d (fname='CWDC_TO_LITR2C', units='gC/m^2/s', &
         avgflag='A', long_name='decomp. of coarse woody debris C to litter 2 C', &
         ptr_col=clm3%g%l%c%ccf%cwdc_to_litr2c, default='inactive')

    call hist_addfld1d (fname='CWDC_TO_LITR3C', units='gC/m^2/s', &
         avgflag='A', long_name='decomp. of coarse woody debris C to litter 3 C', &
         ptr_col=clm3%g%l%c%ccf%cwdc_to_litr3c, default='inactive')

    call hist_addfld1d (fname='LITR1_HR', units='gC/m^2/s', &
         avgflag='A', long_name='het. resp. from litter 1 C', &
         ptr_col=clm3%g%l%c%ccf%litr1_hr, default='inactive')

    call hist_addfld1d (fname='LITR1C_TO_SOIL1C', units='gC/m^2/s', &
         avgflag='A', long_name='decomp. of litter 1 C to SOM 1 C', &
         ptr_col=clm3%g%l%c%ccf%litr1c_to_soil1c)

    call hist_addfld1d (fname='LITR2_HR', units='gC/m^2/s', &
         avgflag='A', long_name='het. resp. from litter 2 C', &
         ptr_col=clm3%g%l%c%ccf%litr2_hr, default='inactive')

    call hist_addfld1d (fname='LITR2C_TO_SOIL2C', units='gC/m^2/s', &
         avgflag='A', long_name='decomp. of litter 2 C to SOM 2 C', &
         ptr_col=clm3%g%l%c%ccf%litr2c_to_soil2c)

    call hist_addfld1d (fname='LITR3_HR', units='gC/m^2/s', &
         avgflag='A', long_name='het. resp. from litter 3 C', &
         ptr_col=clm3%g%l%c%ccf%litr3_hr, default='inactive')

    call hist_addfld1d (fname='LITR3C_TO_SOIL3C', units='gC/m^2/s', &
         avgflag='A', long_name='decomp. of litter 3 C to SOM 3 C', &
         ptr_col=clm3%g%l%c%ccf%litr3c_to_soil3c)

    call hist_addfld1d (fname='SOIL1_HR', units='gC/m^2/s', &
         avgflag='A', long_name='het. resp. from SOM 1 C', &
         ptr_col=clm3%g%l%c%ccf%soil1_hr, default='inactive')

    call hist_addfld1d (fname='SOIL1C_TO_SOIL2C', units='gC/m^2/s', &
         avgflag='A', long_name='decomp. of SOM 1 C to SOM 2 C', &
         ptr_col=clm3%g%l%c%ccf%soil1c_to_soil2c, default='inactive')

    call hist_addfld1d (fname='SOIL2_HR', units='gC/m^2/s', &
         avgflag='A', long_name='het. resp. from SOM 2 C', &
         ptr_col=clm3%g%l%c%ccf%soil2_hr, default='inactive')

    call hist_addfld1d (fname='SOIL2C_TO_SOIL3C', units='gC/m^2/s', &
         avgflag='A', long_name='decomp. of SOM 2 C to SOM 3 C', &
         ptr_col=clm3%g%l%c%ccf%soil2c_to_soil3c, default='inactive')

    call hist_addfld1d (fname='SOIL3_HR', units='gC/m^2/s', &
         avgflag='A', long_name='het. resp. from SOM 3 C', &
         ptr_col=clm3%g%l%c%ccf%soil3_hr, default='inactive')

    call hist_addfld1d (fname='SOIL3C_TO_SOIL4C', units='gC/m^2/s', &
         avgflag='A', long_name='decomp. of SOM 3 C to SOM 4 C', &
         ptr_col=clm3%g%l%c%ccf%soil3c_to_soil4c, default='inactive')

    call hist_addfld1d (fname='SOIL4_HR', units='gC/m^2/s', &
         avgflag='A', long_name='het. resp. from SOM 4 C', &
         ptr_col=clm3%g%l%c%ccf%soil4_hr, default='inactive')

    call hist_addfld1d (fname='LITHR', units='gC/m^2/s', &
         avgflag='A', long_name='litter heterotrophic respiration', &
         ptr_col=clm3%g%l%c%ccf%lithr)

    call hist_addfld1d (fname='SOMHR', units='gC/m^2/s', &
         avgflag='A', long_name='soil organic matter heterotrophic respiration', &
         ptr_col=clm3%g%l%c%ccf%somhr)

    call hist_addfld1d (fname='HR', units='gC/m^2/s', &
         avgflag='A', long_name='total heterotrophic respiration', &
         ptr_col=clm3%g%l%c%ccf%hr)

    call hist_addfld1d (fname='SR', units='gC/m^2/s', &
         avgflag='A', long_name='total soil respiration (HR + root resp)', &
         ptr_col=clm3%g%l%c%ccf%sr)

    call hist_addfld1d (fname='ER', units='gC/m^2/s', &
         avgflag='A', long_name='total ecosystem respiration, autotrophic + heterotrophic', &
         ptr_col=clm3%g%l%c%ccf%er)

    call hist_addfld1d (fname='LITFIRE', units='gC/m^2/s', &
         avgflag='A', long_name='litter fire losses', &
         ptr_col=clm3%g%l%c%ccf%litfire, default='inactive')

    call hist_addfld1d (fname='SOMFIRE', units='gC/m^2/s', &
         avgflag='A', long_name='soil organic matter fire losses', &
         ptr_col=clm3%g%l%c%ccf%somfire, default='inactive')

    call hist_addfld1d (fname='TOTFIRE', units='gC/m^2/s', &
         avgflag='A', long_name='total ecosystem fire losses', &
         ptr_col=clm3%g%l%c%ccf%totfire, default='inactive')

    call hist_addfld1d (fname='NEP', units='gC/m^2/s', &
         avgflag='A', long_name='net ecosystem production, excludes fire, landuse, and harvest flux, positive for sink', &
         ptr_col=clm3%g%l%c%ccf%nep)

    call hist_addfld1d (fname='NBP', units='gC/m^2/s', &
         avgflag='A', long_name='net biome production, includes fire, landuse, and harvest flux, positive for sink', &
         ptr_col=clm3%g%l%c%ccf%nbp)

    call hist_addfld1d (fname='NEE', units='gC/m^2/s', &
         avgflag='A', long_name=&
         'net ecosystem exchange of carbon, includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source', &
         ptr_col=clm3%g%l%c%ccf%nee)

    call hist_addfld1d (fname='COL_FIRE_CLOSS', units='gC/m^2/s', &
         avgflag='A', long_name='total column-level fire C loss', &
         ptr_col=clm3%g%l%c%ccf%col_fire_closs)

    call hist_addfld1d (fname='DWT_SEEDC_TO_LEAF', units='gC/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level leaf', &
         ptr_col=clm3%g%l%c%ccf%dwt_seedc_to_leaf)

    call hist_addfld1d (fname='DWT_SEEDC_TO_DEADSTEM', units='gC/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level deadstem', &
         ptr_col=clm3%g%l%c%ccf%dwt_seedc_to_deadstem)

    call hist_addfld1d (fname='DWT_CONV_CFLUX', units='gC/m^2/s', &
         avgflag='A', long_name='conversion C flux (immediate loss to atm)', &
         ptr_col=clm3%g%l%c%ccf%dwt_conv_cflux)

    call hist_addfld1d (fname='DWT_PROD10C_GAIN', units='gC/m^2/s', &
         avgflag='A', long_name='landcover change-driven addition to 10-yr wood product pool', &
         ptr_col=clm3%g%l%c%ccf%dwt_prod10c_gain)

    call hist_addfld1d (fname='PROD10C_LOSS', units='gC/m^2/s', &
         avgflag='A', long_name='loss from 10-yr wood product pool', &
         ptr_col=clm3%g%l%c%ccf%prod10c_loss)

    call hist_addfld1d (fname='DWT_PROD100C_GAIN', units='gC/m^2/s', &
         avgflag='A', long_name='landcover change-driven addition to 100-yr wood product pool', &
         ptr_col=clm3%g%l%c%ccf%dwt_prod100c_gain)

    call hist_addfld1d (fname='PROD100C_LOSS', units='gC/m^2/s', &
         avgflag='A', long_name='loss from 100-yr wood product pool', &
         ptr_col=clm3%g%l%c%ccf%prod100c_loss)

    call hist_addfld1d (fname='DWT_FROOTC_TO_LITR1C', units='gC/m^2/s', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=clm3%g%l%c%ccf%dwt_frootc_to_litr1c, default='inactive')

    call hist_addfld1d (fname='DWT_FROOTC_TO_LITR2C', units='gC/m^2/s', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=clm3%g%l%c%ccf%dwt_frootc_to_litr2c, default='inactive')

    call hist_addfld1d (fname='DWT_FROOTC_TO_LITR3C', units='gC/m^2/s', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=clm3%g%l%c%ccf%dwt_frootc_to_litr3c, default='inactive')

    call hist_addfld1d (fname='DWT_LIVECROOTC_TO_CWDC', units='gC/m^2/s', &
         avgflag='A', long_name='live coarse root to CWD due to landcover change', &
         ptr_col=clm3%g%l%c%ccf%dwt_livecrootc_to_cwdc, default='inactive')

    call hist_addfld1d (fname='DWT_DEADCROOTC_TO_CWDC', units='gC/m^2/s', &
         avgflag='A', long_name='dead coarse root to CWD due to landcover change', &
         ptr_col=clm3%g%l%c%ccf%dwt_deadcrootc_to_cwdc, default='inactive')

    call hist_addfld1d (fname='DWT_CLOSS', units='gC/m^2/s', &
         avgflag='A', long_name='total carbon loss from land cover conversion', &
         ptr_col=clm3%g%l%c%ccf%dwt_closs)

    call hist_addfld1d (fname='PRODUCT_CLOSS', units='gC/m^2/s', &
         avgflag='A', long_name='total carbon loss from wood product pools', &
         ptr_col=clm3%g%l%c%ccf%product_closs)

    call hist_addfld1d (fname='LAND_USE_FLUX', units='gC/m^2/s', &
         avgflag='A', long_name='total C emitted from land cover conversion and wood product pools', &
         ptr_col=clm3%g%l%c%ccf%landuseflux)

    call hist_addfld1d (fname='LAND_UPTAKE', units='gC/m^2/s', &
         avgflag='A', long_name='NEE minus LAND_USE_FLUX, negative for update', &
         ptr_col=clm3%g%l%c%ccf%landuptake)

#if (defined C13)
    !-------------------------------
    ! C13 flux variables - native to column 
    !-------------------------------

    call hist_addfld1d (fname='C13_M_LEAFC_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf C mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%m_leafc_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_M_LEAFC_TO_LITR2C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf C mortality to litter 2 C', &
         ptr_col=clm3%g%l%c%cc13f%m_leafc_to_litr2c, default='inactive')

    call hist_addfld1d (fname='C13_M_LEAFC_TO_LITR3C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf C mortality to litter 3 C', &
         ptr_col=clm3%g%l%c%cc13f%m_leafc_to_litr3c, default='inactive')

    call hist_addfld1d (fname='C13_M_FROOTC_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root C mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%m_frootc_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_M_FROOTC_TO_LITR2C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root C mortality to litter 2 C', &
         ptr_col=clm3%g%l%c%cc13f%m_frootc_to_litr2c, default='inactive')

    call hist_addfld1d (fname='C13_M_FROOTC_TO_LITR3C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root C mortality to litter 3 C', &
         ptr_col=clm3%g%l%c%cc13f%m_frootc_to_litr3c, default='inactive')

    call hist_addfld1d (fname='C13_M_LEAFC_STORAGE_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf C storage mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%m_leafc_storage_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_M_FROOTC_STORAGE_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root C storage mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%m_frootc_storage_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVESTEMC_STORAGE_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live stem C storage mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%m_livestemc_storage_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADSTEMC_STORAGE_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead stem C storage mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%m_deadstemc_storage_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVECROOTC_STORAGE_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live coarse root C storage mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%m_livecrootc_storage_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADCROOTC_STORAGE_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead coarse root C storage mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%m_deadcrootc_storage_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_M_LEAFC_XFER_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf C transfer mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%m_leafc_xfer_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_M_FROOTC_XFER_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root C transfer mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%m_frootc_xfer_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVESTEMC_XFER_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live stem C transfer mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%m_livestemc_xfer_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADSTEMC_XFER_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead stem C transfer mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%m_deadstemc_xfer_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVECROOTC_XFER_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live coarse root C transfer mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%m_livecrootc_xfer_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADCROOTC_XFER_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead coarse root C transfer mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%m_deadcrootc_xfer_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVESTEMC_TO_CWDC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live stem C mortality to coarse woody debris C', &
         ptr_col=clm3%g%l%c%cc13f%m_livestemc_to_cwdc, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_CWDC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead stem C mortality to coarse woody debris C', &
         ptr_col=clm3%g%l%c%cc13f%m_deadstemc_to_cwdc, default='inactive')

    call hist_addfld1d (fname='C13_M_LIVECROOTC_TO_CWDC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live coarse root C mortality to coarse woody debris C', &
         ptr_col=clm3%g%l%c%cc13f%m_livecrootc_to_cwdc, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_CWDC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead coarse root C mortality to coarse woody debris C', &
         ptr_col=clm3%g%l%c%cc13f%m_deadcrootc_to_cwdc, default='inactive')

    call hist_addfld1d (fname='C13_M_GRESP_STORAGE_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 growth respiration storage mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%m_gresp_storage_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_M_GRESP_XFER_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 growth respiration transfer mortality to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%m_gresp_xfer_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_CWDC_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead stem C to coarse woody debris C by fire', &
         ptr_col=clm3%g%l%c%cc13f%m_deadstemc_to_cwdc_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_CWDC_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead coarse root C to to woody debris C by fire', &
         ptr_col=clm3%g%l%c%cc13f%m_deadcrootc_to_cwdc_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_LITR1C_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 litter 1 C fire loss', &
         ptr_col=clm3%g%l%c%cc13f%m_litr1c_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_LITR2C_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 litter 2 C fire loss', &
         ptr_col=clm3%g%l%c%cc13f%m_litr2c_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_LITR3C_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 litter 3 C fire loss', &
         ptr_col=clm3%g%l%c%cc13f%m_litr3c_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_M_CWDC_TO_FIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 coarse woody debris C fire loss', &
         ptr_col=clm3%g%l%c%cc13f%m_cwdc_to_fire, default='inactive')

    call hist_addfld1d (fname='C13_LEAFC_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf C litterfall to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%leafc_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_LEAFC_TO_LITR2C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf C litterfall to litter 2 C', &
         ptr_col=clm3%g%l%c%cc13f%leafc_to_litr2c, default='inactive')

    call hist_addfld1d (fname='C13_LEAFC_TO_LITR3C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 leaf C litterfall to litter 3 C', &
         ptr_col=clm3%g%l%c%cc13f%leafc_to_litr3c, default='inactive')

    call hist_addfld1d (fname='C13_FROOTC_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root C litterfall to litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%frootc_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_FROOTC_TO_LITR2C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root C litterfall to litter 2 C', &
         ptr_col=clm3%g%l%c%cc13f%frootc_to_litr2c, default='inactive')

    call hist_addfld1d (fname='C13_FROOTC_TO_LITR3C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root C litterfall to litter 3 C', &
         ptr_col=clm3%g%l%c%cc13f%frootc_to_litr3c, default='inactive')

    call hist_addfld1d (fname='C13_CWDC_TO_LITR2C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 decomp. of coarse woody debris C to litter 2 C', &
         ptr_col=clm3%g%l%c%cc13f%cwdc_to_litr2c, default='inactive')

    call hist_addfld1d (fname='C13_CWDC_TO_LITR3C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 decomp. of coarse woody debris C to litter 3 C', &
         ptr_col=clm3%g%l%c%cc13f%cwdc_to_litr3c, default='inactive')

    call hist_addfld1d (fname='C13_LITR1_HR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 het. resp. from litter 1 C', &
         ptr_col=clm3%g%l%c%cc13f%litr1_hr, default='inactive')

    call hist_addfld1d (fname='C13_LITR1C_TO_SOIL1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 decomp. of litter 1 C to SOM 1 C', &
         ptr_col=clm3%g%l%c%cc13f%litr1c_to_soil1c, default='inactive')

    call hist_addfld1d (fname='C13_LITR2_HR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 het. resp. from litter 2 C', &
         ptr_col=clm3%g%l%c%cc13f%litr2_hr, default='inactive')

    call hist_addfld1d (fname='C13_LITR2C_TO_SOIL2C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 decomp. of litter 2 C to SOM 2 C', &
         ptr_col=clm3%g%l%c%cc13f%litr2c_to_soil2c, default='inactive')

    call hist_addfld1d (fname='C13_LITR3_HR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 het. resp. from litter 3 C', &
         ptr_col=clm3%g%l%c%cc13f%litr3_hr, default='inactive')

    call hist_addfld1d (fname='C13_LITR3C_TO_SOIL3C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 decomp. of litter 3 C to SOM 3 C', &
         ptr_col=clm3%g%l%c%cc13f%litr3c_to_soil3c, default='inactive')

    call hist_addfld1d (fname='C13_SOIL1_HR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 het. resp. from SOM 1 C', &
         ptr_col=clm3%g%l%c%cc13f%soil1_hr, default='inactive')

    call hist_addfld1d (fname='C13_SOIL1C_TO_SOIL2C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 decomp. of SOM 1 C to SOM 2 C', &
         ptr_col=clm3%g%l%c%cc13f%soil1c_to_soil2c, default='inactive')

    call hist_addfld1d (fname='C13_SOIL2_HR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 het. resp. from SOM 2 C', &
         ptr_col=clm3%g%l%c%cc13f%soil2_hr, default='inactive')

    call hist_addfld1d (fname='C13_SOIL2C_TO_SOIL3C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 decomp. of SOM 2 C to SOM 3 C', &
         ptr_col=clm3%g%l%c%cc13f%soil2c_to_soil3c, default='inactive')

    call hist_addfld1d (fname='C13_SOIL3_HR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 het. resp. from SOM 3 C', &
         ptr_col=clm3%g%l%c%cc13f%soil3_hr, default='inactive')

    call hist_addfld1d (fname='C13_SOIL3C_TO_SOIL4C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 decomp. of SOM 3 C to SOM 4 C', &
         ptr_col=clm3%g%l%c%cc13f%soil3c_to_soil4c, default='inactive')

    call hist_addfld1d (fname='C13_SOIL4_HR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 het. resp. from SOM 4 C', &
         ptr_col=clm3%g%l%c%cc13f%soil4_hr, default='inactive')

    call hist_addfld1d (fname='C13_LITHR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 litter heterotrophic respiration', &
         ptr_col=clm3%g%l%c%cc13f%lithr)

    call hist_addfld1d (fname='C13_SOMHR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 soil organic matter heterotrophic respiration', &
         ptr_col=clm3%g%l%c%cc13f%somhr)

    call hist_addfld1d (fname='C13_HR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 total heterotrophic respiration', &
         ptr_col=clm3%g%l%c%cc13f%hr)

    call hist_addfld1d (fname='C13_SR', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 total soil respiration (HR + root resp)', &
         ptr_col=clm3%g%l%c%cc13f%sr)

    call hist_addfld1d (fname='C13_ER', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 total ecosystem respiration, autotrophic + heterotrophic', &
         ptr_col=clm3%g%l%c%cc13f%er)

    call hist_addfld1d (fname='C13_LITFIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 litter fire losses', &
         ptr_col=clm3%g%l%c%cc13f%litfire, default='inactive')

    call hist_addfld1d (fname='C13_SOMFIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 soil organic matter fire losses', &
         ptr_col=clm3%g%l%c%cc13f%somfire, default='inactive')

    call hist_addfld1d (fname='C13_TOTFIRE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 total ecosystem fire losses', &
         ptr_col=clm3%g%l%c%cc13f%totfire, default='inactive')

    call hist_addfld1d (fname='C13_NEP', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 net ecosystem production, excludes fire flux, positive for sink', &
         ptr_col=clm3%g%l%c%cc13f%nep)

    call hist_addfld1d (fname='C13_NEE', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 net ecosystem exchange of carbon, includes fire flux, positive for source', &
         ptr_col=clm3%g%l%c%cc13f%nee)

    call hist_addfld1d (fname='C13_COL_FIRE_CLOSS', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 total column-level fire C loss', &
         ptr_col=clm3%g%l%c%cc13f%col_fire_closs)

    call hist_addfld1d (fname='C13_DWT_SEEDC_TO_LEAF', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 seed source to PFT-level leaf', &
         ptr_col=clm3%g%l%c%cc13f%dwt_seedc_to_leaf)

    call hist_addfld1d (fname='C13_DWT_SEEDC_TO_DEADSTEM', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 seed source to PFT-level deadstem', &
         ptr_col=clm3%g%l%c%cc13f%dwt_seedc_to_deadstem)

    call hist_addfld1d (fname='C13_DWT_CONV_CFLUX', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 conversion C flux (immediate loss to atm)', &
         ptr_col=clm3%g%l%c%cc13f%dwt_conv_cflux)

    call hist_addfld1d (fname='C13_DWT_PROD10C_GAIN', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 addition to 10-yr wood product pool', &
         ptr_col=clm3%g%l%c%cc13f%dwt_prod10c_gain)

    call hist_addfld1d (fname='C13_PROD10C_LOSS', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 loss from 10-yr wood product pool', &
         ptr_col=clm3%g%l%c%cc13f%prod10c_loss)

    call hist_addfld1d (fname='C13_DWT_PROD100C_GAIN', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 addition to 100-yr wood product pool', &
         ptr_col=clm3%g%l%c%cc13f%dwt_prod100c_gain)

    call hist_addfld1d (fname='C13_PROD100C_LOSS', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 loss from 100-yr wood product pool', &
         ptr_col=clm3%g%l%c%cc13f%prod100c_loss)

    call hist_addfld1d (fname='C13_DWT_FROOTC_TO_LITR1C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root to litter due to landcover change', &
         ptr_col=clm3%g%l%c%cc13f%dwt_frootc_to_litr1c, default='inactive')

    call hist_addfld1d (fname='C13_DWT_FROOTC_TO_LITR2C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root to litter due to landcover change', &
         ptr_col=clm3%g%l%c%cc13f%dwt_frootc_to_litr2c, default='inactive')

    call hist_addfld1d (fname='C13_DWT_FROOTC_TO_LITR3C', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 fine root to litter due to landcover change', &
         ptr_col=clm3%g%l%c%cc13f%dwt_frootc_to_litr3c, default='inactive')

    call hist_addfld1d (fname='C13_DWT_LIVECROOTC_TO_CWDC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 live coarse root to CWD due to landcover change', &
         ptr_col=clm3%g%l%c%cc13f%dwt_livecrootc_to_cwdc, default='inactive')

    call hist_addfld1d (fname='C13_DWT_DEADCROOTC_TO_CWDC', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 dead coarse root to CWD due to landcover change', &
         ptr_col=clm3%g%l%c%cc13f%dwt_deadcrootc_to_cwdc, default='inactive')

    call hist_addfld1d (fname='C13_DWT_CLOSS', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 total carbon loss from land cover conversion', &
         ptr_col=clm3%g%l%c%cc13f%dwt_closs)

    call hist_addfld1d (fname='C13_PRODUCT_CLOSS', units='gC13/m^2/s', &
         avgflag='A', long_name='C13 total carbon loss from wood product pools', &
         ptr_col=clm3%g%l%c%cc13f%product_closs)
#endif

    !-------------------------------
    ! N flux variables - native to PFT
    !-------------------------------

    call hist_addfld1d (fname='M_LEAFN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_leafn_to_litter, default='inactive')

    call hist_addfld1d (fname='M_FROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_frootn_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LEAFN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N storage mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_leafn_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_FROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N storage mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_frootn_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N storage mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livestemn_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N storage mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadstemn_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N storage mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livecrootn_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N storage mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LEAFN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_leafn_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_FROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_frootn_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livestemn_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadstemn_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livecrootn_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N transfer mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livestemn_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadstemn_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livecrootn_to_litter, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadcrootn_to_litter, default='inactive')

    call hist_addfld1d (fname='M_RETRANSN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='retranslocated N pool mortality', &
         ptr_pft=clm3%g%l%c%p%pnf%m_retransn_to_litter, default='inactive')

    call hist_addfld1d (fname='M_LEAFN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_leafn_to_fire, default='inactive')

    call hist_addfld1d (fname='M_FROOTN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N fire loss ', &
         ptr_pft=clm3%g%l%c%p%pnf%m_frootn_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LEAFN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_leafn_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_FROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_frootn_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livestemn_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadstemn_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livecrootn_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N storage fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LEAFN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_leafn_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_FROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_frootn_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livestemn_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadstemn_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livecrootn_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N transfer fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livestemn_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadstemn_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_TO_LITTER_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadstemn_to_litter_fire, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_livecrootn_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadcrootn_to_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_TO_LITTER_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N fire mortality to litter', &
         ptr_pft=clm3%g%l%c%p%pnf%m_deadcrootn_to_litter_fire, default='inactive')

    call hist_addfld1d (fname='M_RETRANSN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='retranslocated N pool fire loss', &
         ptr_pft=clm3%g%l%c%p%pnf%m_retransn_to_fire, default='inactive')

    call hist_addfld1d (fname='LEAFN_XFER_TO_LEAFN', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N growth from storage', &
         ptr_pft=clm3%g%l%c%p%pnf%leafn_xfer_to_leafn, default='inactive')

    call hist_addfld1d (fname='FROOTN_XFER_TO_FROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N growth from storage', &
         ptr_pft=clm3%g%l%c%p%pnf%frootn_xfer_to_frootn, default='inactive')

    call hist_addfld1d (fname='LIVESTEMN_XFER_TO_LIVESTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N growth from storage', &
         ptr_pft=clm3%g%l%c%p%pnf%livestemn_xfer_to_livestemn, default='inactive')

    call hist_addfld1d (fname='DEADSTEMN_XFER_TO_DEADSTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N growth from storage', &
         ptr_pft=clm3%g%l%c%p%pnf%deadstemn_xfer_to_deadstemn, default='inactive')

    call hist_addfld1d (fname='LIVECROOTN_XFER_TO_LIVECROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N growth from storage', &
         ptr_pft=clm3%g%l%c%p%pnf%livecrootn_xfer_to_livecrootn, default='inactive')

    call hist_addfld1d (fname='DEADCROOTN_XFER_TO_DEADCROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N growth from storage', &
         ptr_pft=clm3%g%l%c%p%pnf%deadcrootn_xfer_to_deadcrootn, default='inactive')

    call hist_addfld1d (fname='LEAFN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N litterfall', &
         ptr_pft=clm3%g%l%c%p%pnf%leafn_to_litter, default='inactive')

    call hist_addfld1d (fname='LEAFN_TO_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N to retranslocated N pool', &
         ptr_pft=clm3%g%l%c%p%pnf%leafn_to_retransn, default='inactive')

    call hist_addfld1d (fname='FROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N litterfall', &
         ptr_pft=clm3%g%l%c%p%pnf%frootn_to_litter, default='inactive')

    call hist_addfld1d (fname='RETRANSN_TO_NPOOL', units='gN/m^2/s', &
         avgflag='A', long_name='deployment of retranslocated N', &
         ptr_pft=clm3%g%l%c%p%pnf%retransn_to_npool)

    call hist_addfld1d (fname='SMINN_TO_NPOOL', units='gN/m^2/s', &
         avgflag='A', long_name='deployment of soil mineral N uptake', &
         ptr_pft=clm3%g%l%c%p%pnf%sminn_to_npool)

    call hist_addfld1d (fname='NPOOL_TO_LEAFN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to leaf N', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_leafn, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_LEAFN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to leaf N storage', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_leafn_storage, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_FROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to fine root N', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_frootn, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_FROOTN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to fine root N storage', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_frootn_storage, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_LIVESTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live stem N', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_livestemn, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_LIVESTEMN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live stem N storage', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_livestemn_storage, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_DEADSTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead stem N', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_deadstemn, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_DEADSTEMN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead stem N storage', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_deadstemn_storage, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_LIVECROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root N', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_livecrootn, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_LIVECROOTN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root N storage', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_livecrootn_storage, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_DEADCROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root N', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_deadcrootn, default='inactive')

    call hist_addfld1d (fname='NPOOL_TO_DEADCROOTN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root N storage', &
         ptr_pft=clm3%g%l%c%p%pnf%npool_to_deadcrootn_storage, default='inactive')

    call hist_addfld1d (fname='LEAFN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pnf%leafn_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='FROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pnf%frootn_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='LIVESTEMN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pnf%livestemn_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='DEADSTEMN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pnf%deadstemn_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='LIVECROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pnf%livecrootn_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='DEADCROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N shift storage to transfer', &
         ptr_pft=clm3%g%l%c%p%pnf%deadcrootn_storage_to_xfer, default='inactive')

    call hist_addfld1d (fname='LIVESTEMN_TO_DEADSTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N turnover', &
         ptr_pft=clm3%g%l%c%p%pnf%livestemn_to_deadstemn, default='inactive')

    call hist_addfld1d (fname='LIVESTEMN_TO_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N to retranslocated N pool', &
         ptr_pft=clm3%g%l%c%p%pnf%livestemn_to_retransn, default='inactive')

    call hist_addfld1d (fname='LIVECROOTN_TO_DEADCROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N turnover', &
         ptr_pft=clm3%g%l%c%p%pnf%livecrootn_to_deadcrootn, default='inactive')

    call hist_addfld1d (fname='LIVECROOTN_TO_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N to retranslocated N pool', &
         ptr_pft=clm3%g%l%c%p%pnf%livecrootn_to_retransn, default='inactive')

    call hist_addfld1d (fname='NDEPLOY', units='gN/m^2/s', &
         avgflag='A', long_name='total N deployed in new growth', &
         ptr_pft=clm3%g%l%c%p%pnf%ndeploy)

    call hist_addfld1d (fname='WOOD_HARVESTN', units='gN/m^2/s', &
         avgflag='A', long_name='wood harvest (to product pools)', &
         ptr_pft=clm3%g%l%c%p%pnf%wood_harvestn)

    call hist_addfld1d (fname='PFT_FIRE_NLOSS', units='gN/m^2/s', &
         avgflag='A', long_name='total pft-level fire N loss', &
         ptr_pft=clm3%g%l%c%p%pnf%pft_fire_nloss)

    !-------------------------------
    ! N flux variables - native to column
    !-------------------------------

    call hist_addfld1d (fname='NDEP_TO_SMINN', units='gN/m^2/s', &
         avgflag='A', long_name='atmospheric N deposition to soil mineral N', &
         ptr_col=clm3%g%l%c%cnf%ndep_to_sminn)

    call hist_addfld1d (fname='NFIX_TO_SMINN', units='gN/m^2/s', &
         avgflag='A', long_name='symbiotic/asymbiotic N fixation to soil mineral N', &
         ptr_col=clm3%g%l%c%cnf%nfix_to_sminn)

    call hist_addfld1d (fname='M_LEAFN_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N mortality to litter 1 N', &
         ptr_col=clm3%g%l%c%cnf%m_leafn_to_litr1n, default='inactive')

    call hist_addfld1d (fname='M_LEAFN_TO_LITR2N', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N mortality to litter 2 N', &
         ptr_col=clm3%g%l%c%cnf%m_leafn_to_litr2n, default='inactive')

    call hist_addfld1d (fname='M_LEAFN_TO_LITR3N', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N mortality to litter 3 N', &
         ptr_col=clm3%g%l%c%cnf%m_leafn_to_litr3n, default='inactive')

    call hist_addfld1d (fname='M_FROOTN_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N mortality to litter 1 N', &
         ptr_col=clm3%g%l%c%cnf%m_frootn_to_litr1n, default='inactive')

    call hist_addfld1d (fname='M_FROOTN_TO_LITR2N', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N mortality to litter 2 N', &
         ptr_col=clm3%g%l%c%cnf%m_frootn_to_litr2n, default='inactive')

    call hist_addfld1d (fname='M_FROOTN_TO_LITR3N', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N mortality to litter 3 N', &
         ptr_col=clm3%g%l%c%cnf%m_frootn_to_litr3n, default='inactive')

    call hist_addfld1d (fname='M_LEAFN_STORAGE_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N storage mortality to litter 1 N', &
         ptr_col=clm3%g%l%c%cnf%m_leafn_storage_to_litr1n, default='inactive')

    call hist_addfld1d (fname='M_FROOTN_STORAGE_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N storage mortality to litter 1 N', &
         ptr_col=clm3%g%l%c%cnf%m_frootn_storage_to_litr1n, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMN_STORAGE_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N storage mortality to litter 1 N', &
         ptr_col=clm3%g%l%c%cnf%m_livestemn_storage_to_litr1n, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_STORAGE_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N storage mortality to litter 1 N', &
         ptr_col=clm3%g%l%c%cnf%m_deadstemn_storage_to_litr1n, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTN_STORAGE_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N storage mortality to litter 1 N', &
         ptr_col=clm3%g%l%c%cnf%m_livecrootn_storage_to_litr1n, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_STORAGE_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N storage mortality to litter 1 N', &
         ptr_col=clm3%g%l%c%cnf%m_deadcrootn_storage_to_litr1n, default='inactive')

    call hist_addfld1d (fname='M_LEAFN_XFER_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N transfer mortality to litter 1 N', &
         ptr_col=clm3%g%l%c%cnf%m_leafn_xfer_to_litr1n, default='inactive')

    call hist_addfld1d (fname='M_FROOTN_XFER_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N transfer mortality to litter 1 N', &
         ptr_col=clm3%g%l%c%cnf%m_frootn_xfer_to_litr1n, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMN_XFER_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N transfer mortality to litter 1 N', &
         ptr_col=clm3%g%l%c%cnf%m_livestemn_xfer_to_litr1n, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_XFER_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N transfer mortality to litter 1 N', &
         ptr_col=clm3%g%l%c%cnf%m_deadstemn_xfer_to_litr1n, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTN_XFER_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N transfer mortality to litter 1 N', &
         ptr_col=clm3%g%l%c%cnf%m_livecrootn_xfer_to_litr1n, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_XFER_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N transfer mortality to litter 1 N', &
         ptr_col=clm3%g%l%c%cnf%m_deadcrootn_xfer_to_litr1n, default='inactive')

    call hist_addfld1d (fname='M_LIVESTEMN_TO_CWDN', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N mortality to coarse woody debris N', &
         ptr_col=clm3%g%l%c%cnf%m_livestemn_to_cwdn, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_TO_CWDN', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N mortality to coarse woody debris N', &
         ptr_col=clm3%g%l%c%cnf%m_deadstemn_to_cwdn, default='inactive')

    call hist_addfld1d (fname='M_LIVECROOTN_TO_CWDN', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N mortality to coarse woody debris N', &
         ptr_col=clm3%g%l%c%cnf%m_livecrootn_to_cwdn, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_TO_CWDN', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N mortality to coarse woody debris N', &
         ptr_col=clm3%g%l%c%cnf%m_deadcrootn_to_cwdn, default='inactive')

    call hist_addfld1d (fname='M_RETRANSN_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='retranslocated N pool mortality to litter 1 N', &
         ptr_col=clm3%g%l%c%cnf%m_retransn_to_litr1n, default='inactive')

    call hist_addfld1d (fname='M_DEADSTEMN_TO_CWDN_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N to coarse woody debris N by fire', &
         ptr_col=clm3%g%l%c%cnf%m_deadstemn_to_cwdn_fire, default='inactive')

    call hist_addfld1d (fname='M_DEADCROOTN_TO_CWDN_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N to to woody debris N by fire', &
         ptr_col=clm3%g%l%c%cnf%m_deadcrootn_to_cwdn_fire, default='inactive')

    call hist_addfld1d (fname='M_LITR1N_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='litter 1 N fire loss', &
         ptr_col=clm3%g%l%c%cnf%m_litr1n_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LITR2N_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='litter 2 N fire loss', &
         ptr_col=clm3%g%l%c%cnf%m_litr2n_to_fire, default='inactive')

    call hist_addfld1d (fname='M_LITR3N_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='litter 3 N fire loss', &
         ptr_col=clm3%g%l%c%cnf%m_litr3n_to_fire, default='inactive')

    call hist_addfld1d (fname='M_CWDN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='coarse woody debris N fire loss', &
         ptr_col=clm3%g%l%c%cnf%m_cwdn_to_fire, default='inactive')

    call hist_addfld1d (fname='LEAFN_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N litterfall to litter 1 N', &
         ptr_col=clm3%g%l%c%cnf%leafn_to_litr1n, default='inactive')

    call hist_addfld1d (fname='LEAFN_TO_LITR2N', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N litterfall to litter 2 N', &
         ptr_col=clm3%g%l%c%cnf%leafn_to_litr2n, default='inactive')

    call hist_addfld1d (fname='LEAFN_TO_LITR3N', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N litterfall to litter 3 N', &
         ptr_col=clm3%g%l%c%cnf%leafn_to_litr3n, default='inactive')

    call hist_addfld1d (fname='FROOTN_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N litterfall to litter 1 N', &
         ptr_col=clm3%g%l%c%cnf%frootn_to_litr1n, default='inactive')

    call hist_addfld1d (fname='FROOTN_TO_LITR2N', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N litterfall to litter 2 N', &
         ptr_col=clm3%g%l%c%cnf%frootn_to_litr2n, default='inactive')

    call hist_addfld1d (fname='FROOTN_TO_LITR3N', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N litterfall to litter 3 N ', &
         ptr_col=clm3%g%l%c%cnf%frootn_to_litr3n, default='inactive')

    call hist_addfld1d (fname='CWDN_TO_LITR2N', units='gN/m^2/s', &
         avgflag='A', long_name='decomp. of coarse woody debris N to litter 2 N', &
         ptr_col=clm3%g%l%c%cnf%cwdn_to_litr2n, default='inactive')

    call hist_addfld1d (fname='CWDN_TO_LITR3N', units='gN/m^2/s', &
         avgflag='A', long_name='decomp. of coarse woody debris N to litter 3 N', &
         ptr_col=clm3%g%l%c%cnf%cwdn_to_litr3n, default='inactive')

    call hist_addfld1d (fname='LITR1N_TO_SOIL1N', units='gN/m^2/s', &
         avgflag='A', long_name='decomp. of litter 1 N to SOM 1 N', &
         ptr_col=clm3%g%l%c%cnf%litr1n_to_soil1n, default='inactive')

    call hist_addfld1d (fname='SMINN_TO_SOIL1N_L1', units='gN/m^2/s', &
         avgflag='A', long_name='mineral N flux for decomp. of litter 1 to SOM 1', &
         ptr_col=clm3%g%l%c%cnf%sminn_to_soil1n_l1, default='inactive')

    call hist_addfld1d (fname='LITR2N_TO_SOIL2N', units='gN/m^2/s', &
         avgflag='A', long_name='decomp. of litter 2 N to SOM 2 N', &
         ptr_col=clm3%g%l%c%cnf%litr2n_to_soil2n, default='inactive')

    call hist_addfld1d (fname='SMINN_TO_SOIL2N_L2', units='gN/m^2/s', &
         avgflag='A', long_name='mineral N flux for decomp. of litter 2 to SOM 2', &
         ptr_col=clm3%g%l%c%cnf%sminn_to_soil2n_l2, default='inactive')

    call hist_addfld1d (fname='LITR3N_TO_SOIL3N', units='gN/m^2/s', &
         avgflag='A', long_name='decomp. of litter 3 N to SOM 3 N', &
         ptr_col=clm3%g%l%c%cnf%litr3n_to_soil3n, default='inactive')

    call hist_addfld1d (fname='SMINN_TO_SOIL3N_L3', units='gN/m^2/s', &
         avgflag='A', long_name='mineral N flux for decomp. of litter 3 to SOM 3', &
         ptr_col=clm3%g%l%c%cnf%sminn_to_soil3n_l3, default='inactive')

    call hist_addfld1d (fname='SOIL1N_TO_SOIL2n', units='gN/m^2/s', &
         avgflag='A', long_name='decomp. of SOM 1 N to SOM 2 N', &
         ptr_col=clm3%g%l%c%cnf%soil1n_to_soil2n, default='inactive')

    call hist_addfld1d (fname='SMINN_TO_SOIL2N_S1', units='gN/m^2/s', &
         avgflag='A', long_name='mineral N flux for decomp. of SOM 1 to SOM 2', &
         ptr_col=clm3%g%l%c%cnf%sminn_to_soil2n_s1, default='inactive')

    call hist_addfld1d (fname='SOIL2N_TO_SOIL3N', units='gN/m^2/s', &
         avgflag='A', long_name='decomp. of SOM 2 N to SOM 3 N', &
         ptr_col=clm3%g%l%c%cnf%soil2n_to_soil3n, default='inactive')

    call hist_addfld1d (fname='SMINN_TO_SOIL3N_S2', units='gN/m^2/s', &
         avgflag='A', long_name='mineral N flux for decomp. of SOM 2 to SOM 3', &
         ptr_col=clm3%g%l%c%cnf%sminn_to_soil3n_s2, default='inactive')

    call hist_addfld1d (fname='SOIL3N_TO_SOIL4N', units='gN/m^2/s', &
         avgflag='A', long_name='decomp. of SOM 3 N to SOM 4 N', &
         ptr_col=clm3%g%l%c%cnf%soil3n_to_soil4n, default='inactive')

    call hist_addfld1d (fname='SMINN_TO_SOIL4N_S3', units='gN/m^2/s', &
         avgflag='A', long_name='mineral N flux for decomp. of SOM 3 to SOM 4', &
         ptr_col=clm3%g%l%c%cnf%sminn_to_soil4n_s3, default='inactive')

    call hist_addfld1d (fname='SOIL4N_TO_SMINN', units='gN/m^2/s', &
         avgflag='A', long_name='N mineralization for decomp. of SOM 4', &
         ptr_col=clm3%g%l%c%cnf%soil4n_to_sminn, default='inactive')

    call hist_addfld1d (fname='SMINN_TO_DENIT_L1S1', units='gN/m^2/s', &
         avgflag='A', long_name='denitrification for decomp. of litter 1 to SOM 1', &
         ptr_col=clm3%g%l%c%cnf%sminn_to_denit_l1s1, default='inactive')

    call hist_addfld1d (fname='SMINN_TO_DENIT_L2S2', units='gN/m^2/s', &
         avgflag='A', long_name='denitrification for decomp. of litter 2 to SOM 2', &
         ptr_col=clm3%g%l%c%cnf%sminn_to_denit_l2s2, default='inactive')

    call hist_addfld1d (fname='SMINN_TO_DENIT_L3S3', units='gN/m^2/s', &
         avgflag='A', long_name='denitrification for decomp. of litter 3 to SOM 3', &
         ptr_col=clm3%g%l%c%cnf%sminn_to_denit_l3s3, default='inactive')

    call hist_addfld1d (fname='SMINN_TO_DENIT_S1S2', units='gN/m^2/s', &
         avgflag='A', long_name='denitrification for decomp. of SOM 1 to SOM 2', &
         ptr_col=clm3%g%l%c%cnf%sminn_to_denit_s1s2, default='inactive')

    call hist_addfld1d (fname='SMINN_TO_DENIT_S2S3', units='gN/m^2/s', &
         avgflag='A', long_name='denitrification for decomp. of SOM 2 to SOM 3', &
         ptr_col=clm3%g%l%c%cnf%sminn_to_denit_s2s3, default='inactive')

    call hist_addfld1d (fname='SMINN_TO_DENIT_S3S4', units='gN/m^2/s', &
         avgflag='A', long_name='denitrification for decomp. of SOM 3 to SOM 4', &
         ptr_col=clm3%g%l%c%cnf%sminn_to_denit_s3s4, default='inactive')

    call hist_addfld1d (fname='SMINN_TO_DENIT_S4', units='gN/m^2/s', &
         avgflag='A', long_name='denitrification for decomp. of SOM 4', &
         ptr_col=clm3%g%l%c%cnf%sminn_to_denit_s4, default='inactive')

    call hist_addfld1d (fname='SMINN_TO_DENIT_EXCESS', units='gN/m^2/s', &
         avgflag='A', long_name='denitrification from excess mineral N pool', &
         ptr_col=clm3%g%l%c%cnf%sminn_to_denit_excess, default='inactive')

    call hist_addfld1d (fname='SMINN_LEACHED', units='gN/m^2/s', &
         avgflag='A', long_name='soil mineral N pool loss to leaching', &
         ptr_col=clm3%g%l%c%cnf%sminn_leached)

    call hist_addfld1d (fname='POTENTIAL_IMMOB', units='gN/m^2/s', &
         avgflag='A', long_name='potential N immobilization', &
         ptr_col=clm3%g%l%c%cnf%potential_immob)

    call hist_addfld1d (fname='ACTUAL_IMMOB', units='gN/m^2/s', &
         avgflag='A', long_name='actual N immobilization', &
         ptr_col=clm3%g%l%c%cnf%actual_immob)

    call hist_addfld1d (fname='SMINN_TO_PLANT', units='gN/m^2/s', &
         avgflag='A', long_name='plant uptake of soil mineral N', &
         ptr_col=clm3%g%l%c%cnf%sminn_to_plant)

    call hist_addfld1d (fname='SUPPLEMENT_TO_SMINN', units='gN/m^2/s', &
         avgflag='A', long_name='supplemental N supply', &
         ptr_col=clm3%g%l%c%cnf%supplement_to_sminn)

    call hist_addfld1d (fname='GROSS_NMIN', units='gN/m^2/s', &
         avgflag='A', long_name='gross rate of N mineralization', &
         ptr_col=clm3%g%l%c%cnf%gross_nmin)

    call hist_addfld1d (fname='NET_NMIN', units='gN/m^2/s', &
         avgflag='A', long_name='net rate of N mineralization', &
         ptr_col=clm3%g%l%c%cnf%net_nmin)

    call hist_addfld1d (fname='DENIT', units='gN/m^2/s', &
         avgflag='A', long_name='total rate of denitrification', &
         ptr_col=clm3%g%l%c%cnf%denit)

    call hist_addfld1d (fname='COL_FIRE_NLOSS', units='gN/m^2/s', &
         avgflag='A', long_name='total column-level fire N loss', &
         ptr_col=clm3%g%l%c%cnf%col_fire_nloss)

    call hist_addfld1d (fname='DWT_SEEDN_TO_LEAF', units='gN/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level leaf', &
         ptr_col=clm3%g%l%c%cnf%dwt_seedn_to_leaf)

    call hist_addfld1d (fname='DWT_SEEDN_TO_DEADSTEM', units='gN/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level deadstem', &
         ptr_col=clm3%g%l%c%cnf%dwt_seedn_to_deadstem)

    call hist_addfld1d (fname='DWT_CONV_NFLUX', units='gN/m^2/s', &
         avgflag='A', long_name='conversion N flux (immediate loss to atm)', &
         ptr_col=clm3%g%l%c%cnf%dwt_conv_nflux)

    call hist_addfld1d (fname='DWT_PROD10N_GAIN', units='gN/m^2/s', &
         avgflag='A', long_name='addition to 10-yr wood product pool', &
         ptr_col=clm3%g%l%c%cnf%dwt_prod10n_gain)

    call hist_addfld1d (fname='PROD10N_LOSS', units='gN/m^2/s', &
         avgflag='A', long_name='loss from 10-yr wood product pool', &
         ptr_col=clm3%g%l%c%cnf%prod10n_loss)

    call hist_addfld1d (fname='DWT_PROD100N_GAIN', units='gN/m^2/s', &
         avgflag='A', long_name='addition to 100-yr wood product pool', &
         ptr_col=clm3%g%l%c%cnf%dwt_prod100n_gain)

    call hist_addfld1d (fname='PROD100N_LOSS', units='gN/m^2/s', &
         avgflag='A', long_name='loss from 100-yr wood product pool', &
         ptr_col=clm3%g%l%c%cnf%prod100n_loss)

    call hist_addfld1d (fname='PRODUCT_NLOSS', units='gN/m^2/s', &
         avgflag='A', long_name='total N loss from wood product pools', &
         ptr_col=clm3%g%l%c%cnf%product_nloss)

    call hist_addfld1d (fname='DWT_FROOTN_TO_LITR1N', units='gN/m^2/s', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=clm3%g%l%c%cnf%dwt_frootn_to_litr1n, default='inactive')

    call hist_addfld1d (fname='DWT_FROOTN_TO_LITR2N', units='gN/m^2/s', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=clm3%g%l%c%cnf%dwt_frootn_to_litr2n, default='inactive')

    call hist_addfld1d (fname='DWT_FROOTN_TO_LITR3N', units='gN/m^2/s', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=clm3%g%l%c%cnf%dwt_frootn_to_litr3n, default='inactive')

    call hist_addfld1d (fname='DWT_LIVECROOTN_TO_CWDN', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root to CWD due to landcover change', &
         ptr_col=clm3%g%l%c%cnf%dwt_livecrootn_to_cwdn, default='inactive')

    call hist_addfld1d (fname='DWT_DEADCROOTN_TO_CWDN', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root to CWD due to landcover change', &
         ptr_col=clm3%g%l%c%cnf%dwt_deadcrootn_to_cwdn, default='inactive')

    call hist_addfld1d (fname='DWT_NLOSS', units='gN/m^2/s', &
         avgflag='A', long_name='total nitrogen loss from landcover conversion', &
         ptr_col=clm3%g%l%c%cnf%dwt_nloss)

    !-------------------------------
    ! PFT ecophysiological variables (pepv) 
    !-------------------------------

    call hist_addfld1d (fname='DORMANT_FLAG', units='none', &
         avgflag='A', long_name='dormancy flag', &
         ptr_pft=clm3%g%l%c%p%pepv%dormant_flag, default='inactive')

    call hist_addfld1d (fname='DAYS_ACTIVE', units='days', &
         avgflag='A', long_name='number of days since last dormancy', &
         ptr_pft=clm3%g%l%c%p%pepv%days_active, default='inactive')

    call hist_addfld1d (fname='ONSET_FLAG', units='none', &
         avgflag='A', long_name='onset flag', &
         ptr_pft=clm3%g%l%c%p%pepv%onset_flag, default='inactive')

    call hist_addfld1d (fname='ONSET_COUNTER', units='days', &
         avgflag='A', long_name='onset days counter', &
         ptr_pft=clm3%g%l%c%p%pepv%onset_counter, default='inactive')

    call hist_addfld1d (fname='ONSET_GDDFLAG', units='none', &
         avgflag='A', long_name='onset flag for growing degree day sum', &
         ptr_pft=clm3%g%l%c%p%pepv%onset_gddflag, default='inactive')

    call hist_addfld1d (fname='ONSET_FDD', units='C degree-days', &
         avgflag='A', long_name='onset freezing degree days counter', &
         ptr_pft=clm3%g%l%c%p%pepv%onset_fdd, default='inactive')

    call hist_addfld1d (fname='ONSET_GDD', units='C degree-days', &
         avgflag='A', long_name='onset growing degree days', &
         ptr_pft=clm3%g%l%c%p%pepv%onset_gdd, default='inactive')

    call hist_addfld1d (fname='ONSET_SWI', units='none', &
         avgflag='A', long_name='onset soil water index', &
         ptr_pft=clm3%g%l%c%p%pepv%onset_swi, default='inactive')

    call hist_addfld1d (fname='OFFSET_FLAG', units='none', &
         avgflag='A', long_name='offset flag', &
         ptr_pft=clm3%g%l%c%p%pepv%offset_flag, default='inactive')

    call hist_addfld1d (fname='OFFSET_COUNTER', units='days', &
         avgflag='A', long_name='offset days counter', &
         ptr_pft=clm3%g%l%c%p%pepv%offset_counter, default='inactive')

    call hist_addfld1d (fname='OFFSET_FDD', units='C degree-days', &
         avgflag='A', long_name='offset freezing degree days counter', &
         ptr_pft=clm3%g%l%c%p%pepv%offset_fdd, default='inactive')

    call hist_addfld1d (fname='OFFSET_SWI', units='none', &
         avgflag='A', long_name='offset soil water index', &
         ptr_pft=clm3%g%l%c%p%pepv%offset_swi, default='inactive')

    call hist_addfld1d (fname='LGSF', units='proportion', &
         avgflag='A', long_name='long growing season factor', &
         ptr_pft=clm3%g%l%c%p%pepv%lgsf, default='inactive')

    call hist_addfld1d (fname='BGLFR', units='1/s', &
         avgflag='A', long_name='background litterfall rate', &
         ptr_pft=clm3%g%l%c%p%pepv%bglfr, default='inactive')

    call hist_addfld1d (fname='BGTR', units='1/s', &
         avgflag='A', long_name='background transfer growth rate', &
         ptr_pft=clm3%g%l%c%p%pepv%bgtr, default='inactive')

    call hist_addfld1d (fname='DAYL',  units='s', &
         avgflag='A', long_name='daylength', &
         ptr_pft=clm3%g%l%c%p%pepv%dayl, default='inactive')

    call hist_addfld1d (fname='PREV_DAYL', units='s', &
         avgflag='A', long_name='daylength from previous timestep', &
         ptr_pft=clm3%g%l%c%p%pepv%prev_dayl, default='inactive')

    call hist_addfld1d (fname='ANNAVG_T2M', units='K', &
         avgflag='A', long_name='annual average 2m air temperature', &
         ptr_pft=clm3%g%l%c%p%pepv%annavg_t2m, default='inactive')

    call hist_addfld1d (fname='TEMPAVG_T2M', units='K', &
         avgflag='A', long_name='temporary average 2m air temperature', &
         ptr_pft=clm3%g%l%c%p%pepv%tempavg_t2m, default='inactive')

    call hist_addfld1d (fname='INIT_GPP', units='gC/m^2/s', &
         avgflag='A', long_name='GPP flux before downregulation', &
         ptr_pft=clm3%g%l%c%p%pepv%gpp, default='inactive')

    call hist_addfld1d (fname='AVAILC', units='gC/m^2/s', &
         avgflag='A', long_name='C flux available for allocation', &
         ptr_pft=clm3%g%l%c%p%pepv%availc, default='inactive')

    call hist_addfld1d (fname='XSMRPOOL_RECOVER', units='gC/m^2/s', &
         avgflag='A', long_name='C flux assigned to recovery of negative xsmrpool', &
         ptr_pft=clm3%g%l%c%p%pepv%xsmrpool_recover)

#if (defined C13)
    call hist_addfld1d (fname='XSMRPOOL_C13RATIO', units='proportion', &
         avgflag='A', long_name='C13/C(12+13) ratio for xsmrpool', &
         ptr_pft=clm3%g%l%c%p%pepv%xsmrpool_c13ratio, default='inactive')
#endif

    call hist_addfld1d (fname='ALLOC_PNOW', units='proportion', &
         avgflag='A', long_name='fraction of current allocation to display as new growth', &
         ptr_pft=clm3%g%l%c%p%pepv%alloc_pnow, default='inactive')

    call hist_addfld1d (fname='C_ALLOMETRY', units='none', &
         avgflag='A', long_name='C allocation index', &
         ptr_pft=clm3%g%l%c%p%pepv%c_allometry, default='inactive')

    call hist_addfld1d (fname='N_ALLOMETRY', units='none', &
         avgflag='A', long_name='N allocation index', &
         ptr_pft=clm3%g%l%c%p%pepv%n_allometry, default='inactive')

    call hist_addfld1d (fname='PLANT_NDEMAND', units='gN/m^2/s', &
         avgflag='A', long_name='N flux required to support initial GPP', &
         ptr_pft=clm3%g%l%c%p%pepv%plant_ndemand)

    call hist_addfld1d (fname='TEMPSUM_POTENTIAL_GPP', units='gC/m^2/yr', &
         avgflag='A', long_name='temporary annual sum of potential GPP', &
         ptr_pft=clm3%g%l%c%p%pepv%tempsum_potential_gpp, default='inactive')

    call hist_addfld1d (fname='ANNSUM_POTENTIAL_GPP', units='gN/m^2/yr', &
         avgflag='A', long_name='annual sum of potential GPP', &
         ptr_pft=clm3%g%l%c%p%pepv%annsum_potential_gpp, default='inactive')

    call hist_addfld1d (fname='TEMPMAX_RETRANSN', units='gN/m^2', &
         avgflag='A', long_name='temporary annual max of retranslocated N pool', &
         ptr_pft=clm3%g%l%c%p%pepv%tempmax_retransn, default='inactive')

    call hist_addfld1d (fname='ANNMAX_RETRANSN', units='gN/m^2', &
         avgflag='A', long_name='annual max of retranslocated N pool', &
         ptr_pft=clm3%g%l%c%p%pepv%annmax_retransn, default='inactive')

    call hist_addfld1d (fname='AVAIL_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='N flux available from retranslocation pool', &
         ptr_pft=clm3%g%l%c%p%pepv%avail_retransn, default='inactive')

    call hist_addfld1d (fname='PLANT_NALLOC', units='gN/m^2/s', &
         avgflag='A', long_name='total allocated N flux', &
         ptr_pft=clm3%g%l%c%p%pepv%plant_nalloc, default='inactive')

    call hist_addfld1d (fname='PLANT_CALLOC', units='gC/m^2/s', &
         avgflag='A', long_name='total allocated C flux', &
         ptr_pft=clm3%g%l%c%p%pepv%plant_calloc, default='inactive')

    call hist_addfld1d (fname='EXCESS_CFLUX', units='gC/m^2/s', &
         avgflag='A', long_name='C flux not allocated due to downregulation', &
         ptr_pft=clm3%g%l%c%p%pepv%excess_cflux, default='inactive')

    call hist_addfld1d (fname='DOWNREG', units='proportion', &
         avgflag='A', long_name='fractional reduction in GPP due to N limitation', &
         ptr_pft=clm3%g%l%c%p%pepv%downreg, default='inactive')

    call hist_addfld1d (fname='PREV_LEAFC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='previous timestep leaf C litterfall flux', &
         ptr_pft=clm3%g%l%c%p%pepv%prev_leafc_to_litter, default='inactive')

    call hist_addfld1d (fname='PREV_FROOTC_TO_LITTER', units='gC/m^2/s', &
         avgflag='A', long_name='previous timestep froot C litterfall flux', &
         ptr_pft=clm3%g%l%c%p%pepv%prev_frootc_to_litter, default='inactive')

    call hist_addfld1d (fname='ANNSUM_NPP', units='gC/m^2/yr', &
         avgflag='A', long_name='annual sum of NPP', &
         ptr_pft=clm3%g%l%c%p%pepv%annsum_npp, default='inactive')

#if (defined C13)
    call hist_addfld1d (fname='RC13_CANAIR', units='proportion', &
         avgflag='A', long_name='C13/C(12+13) for canopy air', &
         ptr_pft=clm3%g%l%c%p%pepv%rc13_canair, default='inactive')

    call hist_addfld1d (fname='RC13_PSNSUN', units='proportion', &
         avgflag='A', long_name='C13/C(12+13) for sunlit photosynthesis', &
         ptr_pft=clm3%g%l%c%p%pepv%rc13_psnsun, default='inactive')

    call hist_addfld1d (fname='RC13_PSNSHA', units='proportion', &
         avgflag='A', long_name='C13/C(12+13) for shaded photosynthesis', &
         ptr_pft=clm3%g%l%c%p%pepv%rc13_psnsha, default='inactive')
#endif

    !-------------------------------
    ! PFT physical state variables not already defined by default
    !-------------------------------

    call hist_addfld1d (fname='EMV', units='proportion', &
         avgflag='A', long_name='vegetation emissivity', &
         ptr_pft=clm3%g%l%c%p%pps%emv, default='inactive')

    call hist_addfld1d (fname='Z0MV', units='m', &
         avgflag='A', long_name='roughness length over vegetation, momentum', &
         ptr_pft=clm3%g%l%c%p%pps%z0mv, default='inactive')

    call hist_addfld1d (fname='Z0HV', units='m', &
         avgflag='A', long_name='roughness length over vegetation, sensible heat', &
         ptr_pft=clm3%g%l%c%p%pps%z0hv, default='inactive')

    call hist_addfld1d (fname='Z0QV', units='m', &
         avgflag='A', long_name='roughness length over vegetation, latent heat', &
         ptr_pft=clm3%g%l%c%p%pps%z0qv, default='inactive')

    call hist_addfld1d (fname='DEWMX', units='mm', &
         avgflag='A', long_name='Maximum allowed dew', &
         ptr_pft=clm3%g%l%c%p%pps%dewmx, default='inactive')

    call hist_addfld1d (fname='LNCSUN', units='gN/m^2', &
         avgflag='A', long_name='leaf N concentration per unit projected LAI', &
         ptr_pft=clm3%g%l%c%p%pps%lncsun, default='inactive')

    call hist_addfld1d (fname='LNCSHA', units='gN/m^2', &
         avgflag='A', long_name='leaf N concentration per unit projected LAI', &
         ptr_pft=clm3%g%l%c%p%pps%lncsha, default='inactive')

    call hist_addfld1d (fname='VCMXSUN', units='umolCO2/m^2/s', &
         avgflag='A', long_name='sunlit leaf Vcmax', &
         ptr_pft=clm3%g%l%c%p%pps%vcmxsun, default='inactive')

    call hist_addfld1d (fname='VCMXSHA', units='umolCO2/m^2/s', &
         avgflag='A', long_name='shaded leaf Vcmax', &
         ptr_pft=clm3%g%l%c%p%pps%vcmxsha, default='inactive')

    call hist_addfld1d (fname='FSUN', units='proportion', &
         avgflag='A', long_name='sunlit fraction of canopy', &
         ptr_pft=clm3%g%l%c%p%pps%fsun, default='inactive')

    call hist_addfld1d (fname='GDIR', units='proportion', &
         avgflag='A', long_name='leaf projection in solar direction', &
         ptr_pft=clm3%g%l%c%p%pps%gdir, default='inactive')

    call hist_addfld1d (fname='CISUN', units='Pa', &
         avgflag='A', long_name='sunlit intracellular CO2', &
         ptr_pft=clm3%g%l%c%p%pps%gdir, default='inactive')

    call hist_addfld1d (fname='CISHA', units='Pa', &
         avgflag='A', long_name='shaded intracellular CO2', &
         ptr_pft=clm3%g%l%c%p%pps%gdir, default='inactive')

#if (defined C13)
    call hist_addfld1d (fname='ALPHAPSNSUN', units='proportion', &
         avgflag='A', long_name='sunlit c13 fractionation', &
         ptr_pft=clm3%g%l%c%p%pps%gdir, default='inactive')

    call hist_addfld1d (fname='ALPHAPSNSHA', units='proportion', &
         avgflag='A', long_name='shaded c13 fractionation', &
         ptr_pft=clm3%g%l%c%p%pps%gdir, default='inactive')
#endif

    call hist_addfld1d (fname='FWET', units='proportion', &
         avgflag='A', long_name='fraction of canopy that is wet', &
         ptr_pft=clm3%g%l%c%p%pps%fwet, default='inactive')

    call hist_addfld1d (fname='FDRY', units='proportion', &
         avgflag='A', long_name='fraction of foliage that is green and dry', &
         ptr_pft=clm3%g%l%c%p%pps%fdry, default='inactive')

    call hist_addfld1d (fname='DT_VEG', units='K', &
         avgflag='A', long_name='change in t_veg, last iteration', &
         ptr_pft=clm3%g%l%c%p%pps%dt_veg, default='inactive')

    call hist_addfld1d (fname='HTOP', units='m', &
         avgflag='A', long_name='canopy top', &
         ptr_pft=clm3%g%l%c%p%pps%htop)

    call hist_addfld1d (fname='HBOT', units='m', &
         avgflag='A', long_name='canopy bottom', &
         ptr_pft=clm3%g%l%c%p%pps%hbot, default='inactive')

    call hist_addfld1d (fname='Z0M', units='m', &
         avgflag='A', long_name='momentum roughness length', &
         ptr_pft=clm3%g%l%c%p%pps%z0m, default='inactive')

    call hist_addfld1d (fname='DISPLA', units='m', &
         avgflag='A', long_name='displacement height', &
         ptr_pft=clm3%g%l%c%p%pps%displa, default='inactive')

    call hist_addfld1d (fname='U10_DUST', units='m/s', &
         avgflag='A', long_name='10-m wind for dust model', &
         ptr_pft=clm3%g%l%c%p%pps%u10, default='inactive')

    call hist_addfld1d (fname='RAM1', units='s/m', &
         avgflag='A', long_name='aerodynamical resistance ', &
         ptr_pft=clm3%g%l%c%p%pps%ram1, default='inactive')

    call hist_addfld1d (fname='FV', units='m/s', &
         avgflag='A', long_name='friction velocity for dust model', &
         ptr_pft=clm3%g%l%c%p%pps%fv, default='inactive')

    call hist_addfld2d (fname='ROOTFR', units='proportion', type2d='levgrnd', &
         avgflag='A', long_name='fraction of roots in each soil layer', &
         ptr_pft=clm3%g%l%c%p%pps%rootfr, default='inactive')
                                                                       
    call hist_addfld2d (fname='ROOTR', units='proportion', type2d='levgrnd', &
         avgflag='A', long_name='effective fraction of roots in each soil layer', &
         ptr_pft=clm3%g%l%c%p%pps%rootr, default='inactive')
                                                                       
    call hist_addfld2d (fname='RRESIS', units='proportion', type2d='levgrnd', &
         avgflag='A', long_name='root resistance in each soil layer', &
         ptr_pft=clm3%g%l%c%p%pps%rresis, default='inactive')
                                                                       
    call hist_addfld2d (fname='ALBD', units='proportion', type2d='numrad', &
         avgflag='A', long_name='surface albedo (direct)', &
         ptr_pft=clm3%g%l%c%p%pps%albd, default='inactive', c2l_scale_type='urbanf')
                                                                        
    call hist_addfld2d (fname='ALBI', units='proportion', type2d='numrad', &
         avgflag='A', long_name='surface albedo (indirect)', &
         ptr_pft=clm3%g%l%c%p%pps%albi, default='inactive', c2l_scale_type='urbanf')
                                                                       
    call hist_addfld2d (fname='FABD', units='proportion', type2d='numrad', &
         avgflag='A', long_name='flux absorbed by veg per unit direct flux', &
         ptr_pft=clm3%g%l%c%p%pps%fabd, default='inactive')
                                                                       
    call hist_addfld2d (fname='FABI', units='proportion', type2d='numrad', &
         avgflag='A', long_name='flux absorbed by veg per unit indirect flux', &
         ptr_pft=clm3%g%l%c%p%pps%fabi, default='inactive')
                                                                       
    call hist_addfld2d (fname='FTDD', units='proportion', type2d='numrad', &
         avgflag='A', long_name='down direct flux below veg per unit dir flx', &
         ptr_pft=clm3%g%l%c%p%pps%ftdd, default='inactive')
                                                                       
    call hist_addfld2d (fname='FTID', units='proportion', type2d='numrad', &
         avgflag='A', long_name='down indirect flux below veg per unit dir flx', &
         ptr_pft=clm3%g%l%c%p%pps%ftid, default='inactive')
                                                                       
    call hist_addfld2d (fname='FTII', units='proportion', type2d='numrad', &
         avgflag='A', long_name='down indirect flux below veg per unit indirect flx', &
         ptr_pft=clm3%g%l%c%p%pps%ftii, default='inactive')
                                                                       
    call hist_addfld2d (fname='OMEGA', units='proportion', type2d='numrad', &
         avgflag='A', long_name='fraction of intercepted radiation that is scattered', &
         ptr_pft=clm3%g%l%c%p%pps%omega, default='inactive')
                                                                       
    call hist_addfld2d (fname='EFF_KID', units='none', type2d='numrad', &
         avgflag='A', long_name='effective extinction coefficient for indirect from direct', &
         ptr_pft=clm3%g%l%c%p%pps%eff_kid, default='inactive')
                                                                       
    call hist_addfld2d (fname='EFF_KII', units='none', type2d='numrad', &
         avgflag='A', long_name='effective extinction coefficient for indirect from indirect', &
         ptr_pft=clm3%g%l%c%p%pps%eff_kii, default='inactive')
                                                                       
    call hist_addfld2d (fname='SUN_FAID', units='proportion', type2d='numrad', &
         avgflag='A', long_name='fraction sun canopy absorbed indirect from direct', &
         ptr_pft=clm3%g%l%c%p%pps%sun_faid, default='inactive')
                                                                       
    call hist_addfld2d (fname='SUN_FAII', units='proportion', type2d='numrad', &
         avgflag='A', long_name='fraction sun canopy absorbed indirect from indirect', &
         ptr_pft=clm3%g%l%c%p%pps%sun_faii, default='inactive')
                                                                       
    call hist_addfld2d (fname='SHA_FAID', units='proportion', type2d='numrad', &
         avgflag='A', long_name='fraction shade canopy absorbed indirect from direct', &
         ptr_pft=clm3%g%l%c%p%pps%sha_faid, default='inactive')
                                                                       
    call hist_addfld2d (fname='SHA_FAII', units='proportion', type2d='numrad', &
         avgflag='A', long_name='fraction shade canopy absorbed indirect from indirect', &
         ptr_pft=clm3%g%l%c%p%pps%sha_faii, default='inactive')

    if ( crop_prog )then

        call hist_addfld1d (fname='GDD0', units='ddays', &
             avgflag='A', long_name='Growing degree days base  0C from planting', &
             ptr_pft=clm3%g%l%c%p%pps%gdd0, default='inactive')

        call hist_addfld1d (fname='GDD8', units='ddays', &
             avgflag='A', long_name='Growing degree days base  8C from planting', &
             ptr_pft=clm3%g%l%c%p%pps%gdd8, default='inactive')

        call hist_addfld1d (fname='GDD10', units='ddays', &
             avgflag='A', long_name='Growing degree days base 10C from planting', &
             ptr_pft=clm3%g%l%c%p%pps%gdd10, default='inactive')

        call hist_addfld1d (fname='GDD020', units='ddays', &
             avgflag='A', long_name='Twenty year average of growing degree days base  0C from planting', &
             ptr_pft=clm3%g%l%c%p%pps%gdd020, default='inactive')

        call hist_addfld1d (fname='GDD820', units='ddays', &
             avgflag='A', long_name='Twenty year average of growing degree days base  8C from planting', &
             ptr_pft=clm3%g%l%c%p%pps%gdd820, default='inactive')

        call hist_addfld1d (fname='GDD1020', units='ddays', &
             avgflag='A', long_name='Twenty year average of growing degree days base 10C from planting', &
             ptr_pft=clm3%g%l%c%p%pps%gdd1020, default='inactive')

        call hist_addfld1d (fname='GDDPLANT', units='ddays', &
             avgflag='A', long_name='Accumulated growing degree days past planting date for crop', &
             ptr_pft=clm3%g%l%c%p%pps%gddplant, default='inactive')

        call hist_addfld1d (fname='GDDHARV', units='ddays', &
             avgflag='A', long_name='Growing degree days (gdd) needed to harvest', &
             ptr_pft=clm3%g%l%c%p%pps%gddmaturity, default='inactive')

        call hist_addfld1d (fname='GDDTSOI', units='ddays', &
             avgflag='A', long_name='Growing degree-days from planting (top two soil layers)', &
             ptr_pft=clm3%g%l%c%p%pps%gddtsoi, default='inactive')

    end if

    !-------------------------------
    ! Column physical state variables not already defined by default
    !-------------------------------

    call hist_addfld1d (fname='EMG', units='proportion', &
         avgflag='A', long_name='ground emissivity', &
         ptr_col=clm3%g%l%c%cps%emg, default='inactive')

    call hist_addfld1d (fname='Z0MG', units='m', &
         avgflag='A', long_name='roughness length over ground, momentum', &
         ptr_col=clm3%g%l%c%cps%z0mg, default='inactive')

    call hist_addfld1d (fname='Z0HG', units='m', &
         avgflag='A', long_name='roughness length over ground, sensible heat', &
         ptr_col=clm3%g%l%c%cps%z0hg, default='inactive')

    call hist_addfld1d (fname='Z0QG', units='m', &
         avgflag='A', long_name='roughness length over ground, latent heat', &
         ptr_col=clm3%g%l%c%cps%z0qg, default='inactive')

    call hist_addfld1d (fname='BETA', units='none', &
         avgflag='A', long_name='coefficient of convective velocity', &
         ptr_col=clm3%g%l%c%cps%beta, default='inactive')

    call hist_addfld1d (fname='ZII', units='m', &
         avgflag='A', long_name='convective boundary height', &
         ptr_col=clm3%g%l%c%cps%zii, default='inactive')

    call hist_addfld1d (fname='WF', units='proportion', &
         avgflag='A', long_name='soil water as frac. of whc for top 0.5 m', &
         ptr_col=clm3%g%l%c%cps%wf, default='inactive')

    call hist_addfld1d (fname='FPI', units='proportion', &
         avgflag='A', long_name='fraction of potential immobilization', &
         ptr_col=clm3%g%l%c%cps%fpi)

    call hist_addfld1d (fname='FPG', units='proportion', &
         avgflag='A', long_name='fraction of potential gpp', &
         ptr_col=clm3%g%l%c%cps%fpg)

    call hist_addfld1d (fname='ANNSUM_COUNTER', units='s', &
         avgflag='A', long_name='seconds since last annual accumulator turnover', &
         ptr_col=clm3%g%l%c%cps%annsum_counter, default='inactive')

    call hist_addfld1d (fname='CANNSUM_NPP', units='gC/m^2/s', &
         avgflag='A', long_name='annual sum of column-level NPP', &
         ptr_col=clm3%g%l%c%cps%cannsum_npp, default='inactive')

    call hist_addfld1d (fname='CANNAVG_T2M', units='K', &
         avgflag='A', long_name='annual average of 2m air temperature', &
         ptr_col=clm3%g%l%c%cps%cannavg_t2m, default='inactive')

    call hist_addfld2d (fname='FRAC_ICEOLD', units='proportion', type2d='levgrnd', &
         avgflag='A', long_name='fraction of ice relative to the tot water', &
         ptr_col=clm3%g%l%c%cps%frac_iceold, default='inactive')

    call hist_addfld2d (fname='EFF_POROSITY', units='proportion', type2d='levgrnd', &
         avgflag='A', long_name='effective porosity = porosity - vol_ice', &
         ptr_col=clm3%g%l%c%cps%eff_porosity, default='inactive')

    call hist_addfld2d (fname='ROOTR_COLUMN', units='proportion', type2d='levgrnd', &
         avgflag='A', long_name='effective fraction of roots in each soil layer', &
         ptr_col=clm3%g%l%c%cps%rootr_column, default='inactive')

    call hist_addfld2d (fname='ALBGRD', units='proportion', type2d='numrad', &
         avgflag='A', long_name='ground albedo (direct)', &
         ptr_col=clm3%g%l%c%cps%albgrd, default='inactive')

    call hist_addfld2d (fname='ALBGRI', units='proportion', type2d='numrad', &
         avgflag='A', long_name='ground albedo (indirect)', &
         ptr_col=clm3%g%l%c%cps%albgri, default='inactive')

    call hist_addfld1d (fname='ME',  units='proportion', &
         avgflag='A', long_name='moisture of extinction', &
         ptr_col=clm3%g%l%c%cps%me, default='inactive')

    call hist_addfld1d (fname='FIRE_PROB',  units='0-1', &
         avgflag='A', long_name='daily fire probability', &
         ptr_col=clm3%g%l%c%cps%fire_prob, default='inactive')

    call hist_addfld1d (fname='MEAN_FIRE_PROB',  units='0-1', &
         avgflag='A', long_name='e-folding mean of daily fire probability', &
         ptr_col=clm3%g%l%c%cps%mean_fire_prob)

    call hist_addfld1d (fname='FIRESEASONL',  units='days', &
         avgflag='A', long_name='annual fire season length', &
         ptr_col=clm3%g%l%c%cps%fireseasonl)

    call hist_addfld1d (fname='FAREA_BURNED',  units='proportion', &
         avgflag='A', long_name='timestep fractional area burned', &
         ptr_col=clm3%g%l%c%cps%farea_burned, default='inactive')

    call hist_addfld1d (fname='ANN_FAREA_BURNED',  units='proportion', &
         avgflag='A', long_name='annual total fractional area burned', &
         ptr_col=clm3%g%l%c%cps%ann_farea_burned)

    !-------------------------------
    ! Energy flux variables not already defined by default - native PFT 
    !-------------------------------

    call hist_addfld1d (fname='PARSUN', units='W/m^2', &
         avgflag='A', long_name='average absorbed PAR for sunlit leaves', &
         ptr_pft=clm3%g%l%c%p%pef%parsun, default='inactive')

    call hist_addfld1d (fname='PARSHA', units='W/m^2', &
         avgflag='A', long_name='average absorbed PAR for shaded leaves', &
         ptr_pft=clm3%g%l%c%p%pef%parsha, default='inactive')

    call hist_addfld1d (fname='DLRAD', units='W/m^2', &
         avgflag='A', long_name='downward longwave radiation below the canopy', &
         ptr_pft=clm3%g%l%c%p%pef%dlrad, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='ULRAD', units='W/m^2', &
         avgflag='A', long_name='upward longwave radiation above the canopy', &
         ptr_pft=clm3%g%l%c%p%pef%ulrad, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='EFLX_LH_TOT', units='W/m^2', &
         avgflag='A', long_name='total latent heat flux [+ to atm]', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lh_tot, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='EFLX_SOIL_GRND', units='W/m^2', &
         avgflag='A', long_name='soil heat flux [+ into soil]', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_soil_grnd, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='CGRND', units='W/m^2/K', &
         avgflag='A', long_name='deriv. of soil energy flux wrt to soil temp', &
         ptr_pft=clm3%g%l%c%p%pef%cgrnd, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='CGRNDL', units='W/m^2/K', &
         avgflag='A', long_name='deriv. of soil latent heat flux wrt soil temp', &
         ptr_pft=clm3%g%l%c%p%pef%cgrndl, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='CGRNDS', units='W/m^2/K', &
         avgflag='A', long_name='deriv. of soil sensible heat flux wrt soil temp', &
         ptr_pft=clm3%g%l%c%p%pef%cgrnds, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='EFLX_GNET', units='W/m^2', &
         avgflag='A', long_name='net heat flux into ground', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_gnet, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='DGNETDT', units='W/m^2/K', &
         avgflag='A', long_name='derivative of net ground heat flux wrt soil temp', &
         ptr_pft=clm3%g%l%c%p%pef%dgnetdT, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld2d (fname='SUN_ADD', units='W/m^2', type2d='numrad', &
         avgflag='A', long_name='sun canopy absorbed direct from direct', &
         ptr_pft=clm3%g%l%c%p%pef%sun_add, default='inactive')

    call hist_addfld2d (fname='TOT_AID', units='W/m^2', type2d='numrad', &
         avgflag='A', long_name='total canopy absorbed indirect from direct', &
         ptr_pft=clm3%g%l%c%p%pef%tot_aid, default='inactive')

    call hist_addfld2d (fname='SUN_AID', units='W/m^2', type2d='numrad', &
         avgflag='A', long_name='sun canopy absorbed indirect from direct', &
         ptr_pft=clm3%g%l%c%p%pef%sun_aid, default='inactive')

    call hist_addfld2d (fname='SUN_AII', units='W/m^2', type2d='numrad', &
         avgflag='A', long_name='sun canopy absorbed indirect from indirect', &
         ptr_pft=clm3%g%l%c%p%pef%sun_aii, default='inactive')

    call hist_addfld2d (fname='SHA_AID', units='W/m^2', type2d='numrad', &
         avgflag='A', long_name='shade canopy absorbed indirect from direct', &
         ptr_pft=clm3%g%l%c%p%pef%sha_aid, default='inactive')

    call hist_addfld2d (fname='SHA_AII', units='W/m^2', type2d='numrad', &
         avgflag='A', long_name='shade canopy absorbed indirect from indirect', &
         ptr_pft=clm3%g%l%c%p%pef%sha_aii, default='inactive')

    call hist_addfld2d (fname='SUN_ATOT', units='W/m^2', type2d='numrad', &
         avgflag='A', long_name='sun canopy total absorbed', &
         ptr_pft=clm3%g%l%c%p%pef%sun_atot, default='inactive')

    call hist_addfld2d (fname='SHA_ATOT', units='W/m^2', type2d='numrad', &
         avgflag='A', long_name='shade canopy total absorbed', &
         ptr_pft=clm3%g%l%c%p%pef%sha_atot, default='inactive')

    call hist_addfld2d (fname='SUN_ALF', units='W/m^2', type2d='numrad', &
         avgflag='A', long_name='sun canopy total absorbed by leaves', &
         ptr_pft=clm3%g%l%c%p%pef%sun_alf, default='inactive')

    call hist_addfld2d (fname='SHA_ALF', units='W/m^2', type2d='numrad', &
         avgflag='A', long_name='shade canopy total absored by leaves', &
         ptr_pft=clm3%g%l%c%p%pef%sha_alf, default='inactive')

    call hist_addfld2d (fname='SUN_APERLAI', units='W/m^2', type2d='numrad', &
         avgflag='A', long_name='sun canopy total absorbed per unit LAI', &
         ptr_pft=clm3%g%l%c%p%pef%sun_aperlai, default='inactive')

    call hist_addfld2d (fname='SHA_APERLAI', units='W/m^2', type2d='numrad', &
         avgflag='A', long_name='shade canopy total absorbed per unit LAI', &
         ptr_pft=clm3%g%l%c%p%pef%sha_aperlai, default='inactive')
                                                                       
    !-------------------------------
    ! Water flux variables not already defined by default  - native PFT
    !-------------------------------

    call hist_addfld1d (fname='QFLX_RAIN_GRND', units='mm H2O/s', &
         avgflag='A', long_name='rain on ground after interception', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_rain_grnd, default='inactive')

    call hist_addfld1d (fname='QFLX_SNOW_GRND', units='mm H2O/s', &
         avgflag='A', long_name='snow on ground after interception', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_snow_grnd, default='inactive')

    call hist_addfld1d (fname='QFLX_EVAP_VEG', units='mm H2O/s', &
         avgflag='A', long_name='vegetation evaporation', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_evap_veg, default='inactive')

    call hist_addfld1d (fname='QFLX_EVAP_TOT', units='mm H2O/s', &
         avgflag='A', long_name='qflx_evap_soi + qflx_evap_can + qflx_tran_veg', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_evap_tot, default='inactive', c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QFLX_DEW_GRND', units='mm H2O/s', &
         avgflag='A', long_name='ground surface dew formation', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_dew_grnd, default='inactive')

    call hist_addfld1d (fname='QFLX_SUB_SNOW', units='mm H2O/s', &
         avgflag='A', long_name='sublimation rate from snow pack', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_sub_snow, default='inactive')

    call hist_addfld1d (fname='QFLX_DEW_SNOW', units='mm H2O/s', &
         avgflag='A', long_name='surface dew added to snow pacK', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_dew_snow, default='inactive')

    
#endif

#if (defined CASA)
    call hist_addfld1d (fname='CO2FLUX', units='g/m2/s',  &
         avgflag='A', long_name='net CO2 flux', &
         ptr_pft=clm3%g%l%c%p%pps%co2flux, set_lake=0._r8)

    call hist_addfld1d (fname='FNPP', units='g/m2/s',  &
         avgflag='A', long_name='net primary production', &
         ptr_pft=clm3%g%l%c%p%pps%fnpp, set_lake=0._r8)

    call hist_addfld1d (fname='PLAI', units='m2/m2',  &
         avgflag='A', long_name='Prognostic Leaf Area Index', &
         ptr_pft=clm3%g%l%c%p%pps%plai, set_lake=0._r8)

    call hist_addfld1d (fname='PET', units='mm h2o/s',  &
         avgflag='A', long_name='Potential Evapotranspiration', &
         ptr_pft=clm3%g%l%c%p%pps%pet, set_lake=0._r8)

    call hist_addfld1d (fname='DEGDAY', units='deg C',  &
         avgflag='A', long_name='Accumulated degree days', &
         ptr_pft=clm3%g%l%c%p%pps%degday, set_lake=0._r8)

    call hist_addfld1d (fname='TDAYAVG', units='deg C',  &
         avgflag='A', long_name='Daily Averaged Temperature', &
         ptr_pft=clm3%g%l%c%p%pps%tdayavg, set_lake=0._r8)

    ! the next three will not give bfb upon restart, due to the frequency of
    ! calculation; however, they do not affect other variables (slevis)

    call hist_addfld1d (fname='STRESST', units=' ',  &
         avgflag='A', long_name='temperature stress function for leaf loss', &
         ptr_pft=clm3%g%l%c%p%pps%stressT, set_lake=0._r8)

    call hist_addfld1d (fname='STRESSW', units=' ',  &
         avgflag='A', long_name='water stress function for leaf loss', &
         ptr_pft=clm3%g%l%c%p%pps%stressW, set_lake=0._r8)

    call hist_addfld1d (fname='STRESSCD', units=' ', &
         avgflag='A', long_name='cold + drought stress function for leaf loss',&
         ptr_pft=clm3%g%l%c%p%pps%stressCD, set_lake=0._r8)

    call hist_addfld1d (fname='LGROW', units=' ',  &
         avgflag='A', long_name='growing season index (0 or 1)', &
         ptr_pft=clm3%g%l%c%p%pps%lgrow, set_lake=0._r8)

    call hist_addfld1d (fname='ISEABEG', units=' ',  &
         avgflag='A', long_name='index for start of growing season', &
         ptr_pft=clm3%g%l%c%p%pps%iseabeg, set_lake=0._r8)

    call hist_addfld1d (fname='NSTEPBEG', units=' ',  &
         avgflag='A', long_name='nstep at start of growing season', &
         ptr_pft=clm3%g%l%c%p%pps%nstepbeg, set_lake=0._r8)

    call hist_addfld1d (fname='BGTEMP', units=' ',  &
         avgflag='A', long_name='temperature dependence', &
         ptr_pft=clm3%g%l%c%p%pps%bgtemp, set_lake=0._r8)

    call hist_addfld1d (fname='BGMOIST', units=' ',  &
         avgflag='A', long_name='moisture dependence', &
         ptr_pft=clm3%g%l%c%p%pps%bgmoist, set_lake=0._r8)

    call hist_addfld1d (fname='EXCESSC', units='g/m2/timestep',  &
         avgflag='A', long_name='excess carbon', &
         ptr_pft=clm3%g%l%c%p%pps%excessC, set_lake=0._r8)

    call hist_addfld1d (fname='CFLUX', units='g/m2/s',  &
         avgflag='A', long_name='total Carbon flux', &
         ptr_pft=clm3%g%l%c%p%pps%Cflux, set_lake=0._r8)

    call hist_addfld1d (fname='XSCPOOL', units='g/m2',  &
         avgflag='A', long_name='total excess Carbon', &
         ptr_pft=clm3%g%l%c%p%pps%XSCpool, set_lake=0._r8)

    call hist_addfld1d (fname='WLIM', units=' ',  &
         avgflag='A', long_name='Water limitation used in bgmoist (atmp factor)', &
         ptr_pft=clm3%g%l%c%p%pps%Wlim, set_lake=0._r8)

    call hist_addfld1d (fname='SOILT', units='deg C',  &
         avgflag='A', long_name='Soil temperature for top 30cm', &
         ptr_pft=clm3%g%l%c%p%pps%soilt, set_lake=0._r8)

    call hist_addfld1d (fname='SMOIST', units='mm3/mm3',  &
         avgflag='A', long_name='Soil moisture for top 30cm', &
         ptr_pft=clm3%g%l%c%p%pps%smoist, set_lake=0._r8)

    call hist_addfld1d (fname='WATOPT', units='mm3/mm3',  &
         avgflag='A', long_name='watopt for entire column', &
         ptr_pft=clm3%g%l%c%p%pps%watoptc, set_lake=0._r8)

    call hist_addfld1d (fname='WATDRY', units='mm3/mm3',  &
         avgflag='A', long_name='watdry for entire column', &
         ptr_pft=clm3%g%l%c%p%pps%watdryc, set_lake=0._r8)

   call hist_addfld2d (fname='LIVEFR', units='g/m2', type2d='nlive', &
        avgflag='A', long_name='live fraction for leaf, wood, froot', &
        ptr_pft=clm3%g%l%c%p%pps%livefr, set_lake=0._r8)

   call hist_addfld2d (fname='CLOSS', units='g/m2/s', type2d='npools', & 
        avgflag='A', long_name='Amt. of Carbon lost by pool', &
        ptr_pft=clm3%g%l%c%p%pps%Closs, set_lake=0._r8)

   call hist_addfld2d (fname='CTRANS', units='g/m2/s', type2d='npool_t', & 
        avgflag='A', long_name='Amt. of Carbon transferred out of pool types', &
        ptr_pft=clm3%g%l%c%p%pps%Ctrans, set_lake=0._r8)

   call hist_addfld2d (fname='RESP_C', units='g/m2/s', type2d='npools', &  
        avgflag='A', long_name='Amt. of Carbon lost to atm by pool', &
        ptr_pft=clm3%g%l%c%p%pps%Resp_C, set_lake=0._r8)

   call hist_addfld2d (fname='TPOOL_C', units='g/m2', type2d='npools', &  
        avgflag='A', long_name='Total Carbon for pool', &
        ptr_pft=clm3%g%l%c%p%pps%Tpool_C, set_lake=0._r8)

   ! Summary variables added for the C-LAMP Experiments

    call hist_addfld1d (fname='AGNPP', units='gC/m2/s',  &
         avgflag='A', long_name='above-ground net primary production', &
         ptr_pft=clm3%g%l%c%p%pps%casa_agnpp, set_lake=0._r8)

    call hist_addfld1d (fname='AR', units='gC/m2/s',  &
         avgflag='A', long_name='autotrophic respiration', &
         ptr_pft=clm3%g%l%c%p%pps%casa_ar, set_lake=0._r8)

    call hist_addfld1d (fname='BGNPP', units='gC/m2/s',  &
         avgflag='A', long_name='below-ground net primary production', &
         ptr_pft=clm3%g%l%c%p%pps%casa_bgnpp, set_lake=0._r8)

    call hist_addfld1d (fname='CWDC', units='gC/m2',  &
         avgflag='A', long_name='coarse woody debris C', &
         ptr_pft=clm3%g%l%c%p%pps%casa_cwdc, set_lake=0._r8)

    call hist_addfld1d (fname='CWDC_HR', units='gC/m2/s',  &
         avgflag='A', long_name='coarse woody debris heterotrophic respiration', &
         ptr_pft=clm3%g%l%c%p%pps%casa_cwdc_hr, set_lake=0._r8)

    call hist_addfld1d (fname='CWDC_LOSS', units='gC/m2/s',  &
         avgflag='A', long_name='coarse woody debris C loss', &
         ptr_pft=clm3%g%l%c%p%pps%casa_cwdc_loss, set_lake=0._r8)

    call hist_addfld1d (fname='FROOTC', units='gC/m2',  &
         avgflag='A', long_name='fine root C', &
         ptr_pft=clm3%g%l%c%p%pps%casa_frootc, set_lake=0._r8)

    call hist_addfld1d (fname='FROOTC_ALLOC', units='gC/m2/s',  &
         avgflag='A', long_name='fine root C allocation', &
         ptr_pft=clm3%g%l%c%p%pps%casa_frootc_alloc, set_lake=0._r8)

    call hist_addfld1d (fname='FROOTC_LOSS', units='gC/m2/s',  &
         avgflag='A', long_name='fine root C loss', &
         ptr_pft=clm3%g%l%c%p%pps%casa_frootc_loss, set_lake=0._r8)

    call hist_addfld1d (fname='GPP', units='gC/m2/s',  &
         avgflag='A', long_name='gross primary production', &
         ptr_pft=clm3%g%l%c%p%pps%casa_gpp, set_lake=0._r8)

    call hist_addfld1d (fname='HR', units='gC/m2/s',  &
         avgflag='A', long_name='total heterotrophic respiration', &
         ptr_pft=clm3%g%l%c%p%pps%casa_hr, set_lake=0._r8)

    call hist_addfld1d (fname='LEAFC', units='gC/m2',  &
         avgflag='A', long_name='leaf C', &
         ptr_pft=clm3%g%l%c%p%pps%casa_leafc, set_lake=0._r8)

    call hist_addfld1d (fname='LEAFC_ALLOC', units='gC/m2/s',  &
         avgflag='A', long_name='leaf C allocation', &
         ptr_pft=clm3%g%l%c%p%pps%casa_leafc_alloc, set_lake=0._r8)

    call hist_addfld1d (fname='LEAFC_LOSS', units='gC/m2/s',  &
         avgflag='A', long_name='leaf C loss', &
         ptr_pft=clm3%g%l%c%p%pps%casa_leafc_loss, set_lake=0._r8)

    call hist_addfld1d (fname='LITTERC', units='gC/m2',  &
         avgflag='A', long_name='litter C', &
         ptr_pft=clm3%g%l%c%p%pps%casa_litterc, set_lake=0._r8)

    call hist_addfld1d (fname='LITTERC_HR', units='gC/m2/s',  &
         avgflag='A', long_name='litter heterotrophic respiration', &
         ptr_pft=clm3%g%l%c%p%pps%casa_litterc_hr, set_lake=0._r8)

    call hist_addfld1d (fname='LITTERC_LOSS', units='gC/m2/s',  &
         avgflag='A', long_name='litter C loss', &
         ptr_pft=clm3%g%l%c%p%pps%casa_litterc_loss, set_lake=0._r8)

    call hist_addfld1d (fname='NEE', units='gC/m2/s',  &
         avgflag='A', long_name='net ecosystem exchange (GPP, Re, and fire)', &
         ptr_pft=clm3%g%l%c%p%pps%casa_nee, set_lake=0._r8)

    call hist_addfld1d (fname='NEP', units='gC/m2/s',  &
         avgflag='A', long_name='net ecosystem production', &
         ptr_pft=clm3%g%l%c%p%pps%casa_nep, set_lake=0._r8)

    call hist_addfld1d (fname='NPP', units='gC/m2/s',  &
         avgflag='A', long_name='net primary production', &
         ptr_pft=clm3%g%l%c%p%pps%casa_npp, set_lake=0._r8)

    call hist_addfld1d (fname='SOILC', units='gC/m2',  &
         avgflag='A', long_name='total soil organic matter C (excluding CWDC and LITTERC)', &
         ptr_pft=clm3%g%l%c%p%pps%casa_soilc, set_lake=0._r8)

    call hist_addfld1d (fname='SOILC_HR', units='gC/m2/s',  &
         avgflag='A', long_name='soil heterotrophic respiration', &
         ptr_pft=clm3%g%l%c%p%pps%casa_soilc_hr, set_lake=0._r8)

    call hist_addfld1d (fname='SOILC_LOSS', units='gC/m2/s',  &
         avgflag='A', long_name='soil C loss', &
         ptr_pft=clm3%g%l%c%p%pps%casa_soilc_loss, set_lake=0._r8)

    call hist_addfld1d (fname='WOODC', units='gC/m2',  &
         avgflag='A', long_name='wood C', &
         ptr_pft=clm3%g%l%c%p%pps%casa_woodc, set_lake=0._r8)

    call hist_addfld1d (fname='WOODC_ALLOC', units='gC/m2/s',  &
         avgflag='A', long_name='wood C allocation', &
         ptr_pft=clm3%g%l%c%p%pps%casa_woodc_alloc, set_lake=0._r8)

    call hist_addfld1d (fname='WOODC_LOSS', units='gC/m2/s',  &
         avgflag='A', long_name='wood C loss', &
         ptr_pft=clm3%g%l%c%p%pps%casa_woodc_loss, set_lake=0._r8)

#endif

    call hist_addfld1d (fname='SNORDSL', units='m^-6', &
         avgflag='A', long_name='top snow layer effective grain radius', &
         ptr_col=clm3%g%l%c%cps%snw_rds_top, set_lake=spval, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='SNOTTOPL', units='K/m', &
         avgflag='A', long_name='snow temperature (top layer)', &
         ptr_col=clm3%g%l%c%cps%snot_top, set_lake=spval, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='SNOdTdzL', units='K/m', &
         avgflag='A', long_name='top snow layer temperature gradient (land)', &
         ptr_col=clm3%g%l%c%cps%dTdz_top, set_lake=spval, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='SNOLIQFL', units='fraction', &
         avgflag='A', long_name='top snow layer liquid water fraction (land)', &
         ptr_col=clm3%g%l%c%cps%sno_liq_top, set_lake=spval, set_urb=spval, &
         default='inactive')

    call hist_addfld1d (fname='SNOFSRVD', units='watt/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation from snow', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_sno_vd, &
         default='inactive')

    call hist_addfld1d (fname='SNOFSRND', units='watt/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation from snow', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_sno_nd, &
         default='inactive')
   
    call hist_addfld1d (fname='SNOFSRVI', units='watt/m^2',  &
         avgflag='A', long_name='diffuse vis reflected solar radiation from snow', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_sno_vi, &
         default='inactive')

    call hist_addfld1d (fname='SNOFSRNI', units='watt/m^2',  &
         avgflag='A', long_name='diffuse nir reflected solar radiation from snow', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_sno_ni, &
         default='inactive')

    call hist_addfld1d (fname='SNOFSDSVD', units='watt/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation on snow', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_sno_vd, &
         default='inactive')
   
    call hist_addfld1d (fname='SNOFSDSND', units='watt/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation on snow', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_sno_nd, &
         default='inactive')
   
    call hist_addfld1d (fname='SNOFSDSVI', units='watt/m^2',  &
         avgflag='A', long_name='diffuse vis incident solar radiation on snow', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_sno_vi, &
         default='inactive')
   
    call hist_addfld1d (fname='SNOFSDSNI', units='watt/m^2',  &
         avgflag='A', long_name='diffuse nir incident solar radiation on snow', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_sno_ni, &
         default='inactive')

    call hist_addfld1d (fname='H2OSNO_TOP', units='kg', &
         avgflag='A', long_name='mass of snow in top snow layer', &
         ptr_col=clm3%g%l%c%cps%h2osno_top, set_lake=spval, set_urb=spval)

    call hist_addfld1d (fname='SNOBCMCL', units='kg/m2', &
         avgflag='A', long_name='mass of BC in snow column', &
         ptr_col=clm3%g%l%c%cps%mss_bc_col, set_lake=spval, set_urb=spval)
    
    call hist_addfld1d (fname='SNOBCMSL', units='kg/m2', &
         avgflag='A', long_name='mass of BC in top snow layer', &
         ptr_col=clm3%g%l%c%cps%mss_bc_top, set_lake=spval, set_urb=spval)

    call hist_addfld1d (fname='SNOOCMCL', units='kg/m2', &
         avgflag='A', long_name='mass of OC in snow column', &
         ptr_col=clm3%g%l%c%cps%mss_oc_col, set_lake=spval, set_urb=spval)
   
    call hist_addfld1d (fname='SNOOCMSL', units='kg/m2', &
         avgflag='A', long_name='mass of OC in top snow layer', &
         ptr_col=clm3%g%l%c%cps%mss_oc_top, set_lake=spval, set_urb=spval)

    call hist_addfld1d (fname='SNODSTMCL', units='kg/m2', &
         avgflag='A', long_name='mass of dust in snow column', &
         ptr_col=clm3%g%l%c%cps%mss_dst_col, set_lake=spval, set_urb=spval)
    
    call hist_addfld1d (fname='SNODSTMSL', units='kg/m2', &
         avgflag='A', long_name='mass of dust in top snow layer', &
         ptr_col=clm3%g%l%c%cps%mss_dst_top, set_lake=spval, set_urb=spval)

    call hist_addfld1d (fname='DSTDEP', units='kg/m^2/s', &
         avgflag='A', long_name='total dust deposition (dry+wet) from atmosphere', &
         ptr_col=clm3%g%l%c%cwf%flx_dst_dep, set_lake=spval, set_urb=spval)

    call hist_addfld1d (fname='BCDEP', units='kg/m^2/s', &
         avgflag='A', long_name='total BC deposition (dry+wet) from atmosphere', &
         ptr_col=clm3%g%l%c%cwf%flx_bc_dep, set_lake=spval, set_urb=spval)
   
    call hist_addfld1d (fname='OCDEP', units='kg/m^2/s', &
         avgflag='A', long_name='total OC deposition (dry+wet) from atmosphere', &
         ptr_col=clm3%g%l%c%cwf%flx_oc_dep, set_lake=spval, set_urb=spval)

#if (defined SNICAR_FRC)
    call hist_addfld1d (fname='SNOAERFRCL', units='watt/m^2', &
         avgflag='A', long_name='surface forcing of all aerosols in snow (land) ', &
         ptr_pft=clm3%g%l%c%p%pef%sfc_frc_aer, set_lake=spval, set_urb=spval)
   
    call hist_addfld1d (fname='SNOAERFRC2L', units='watt/m^2', &
         avgflag='A', long_name='surface forcing of all aerosols in snow, averaged only when snow is present (land)', &
         ptr_pft=clm3%g%l%c%p%pef%sfc_frc_aer_sno, set_lake=spval, set_urb=spval)

    call hist_addfld1d (fname='SNOBCFRCL', units='watt/m^2', &
         avgflag='A', long_name='surface forcing of BC in snow (land) ', &
         ptr_pft=clm3%g%l%c%p%pef%sfc_frc_bc, set_lake=spval, set_urb=spval)
   
    call hist_addfld1d (fname='SNOBCFRC2L', units='watt/m^2', &
         avgflag='A', long_name='surface forcing of BC in snow, averaged only when snow is present (land)', &
         ptr_pft=clm3%g%l%c%p%pef%sfc_frc_bc_sno, set_lake=spval, set_urb=spval)

    call hist_addfld1d (fname='SNOOCFRCL', units='watt/m^2', &
         avgflag='A', long_name='surface forcing of OC in snow (land) ', &
         ptr_pft=clm3%g%l%c%p%pef%sfc_frc_oc, set_lake=spval, set_urb=spval)
   
    call hist_addfld1d (fname='SNOOCFRC2L', units='watt/m^2', &
         avgflag='A', long_name='surface forcing of OC in snow, averaged only when snow is present (land)', &
         ptr_pft=clm3%g%l%c%p%pef%sfc_frc_oc_sno, set_lake=spval, set_urb=spval)

    call hist_addfld1d (fname='SNODSTFRCL', units='watt/m^2', &
         avgflag='A', long_name='surface forcing of dust in snow (land) ', &
         ptr_pft=clm3%g%l%c%p%pef%sfc_frc_dst, set_lake=spval, set_urb=spval)
    
    call hist_addfld1d (fname='SNODSTFRC2L', units='watt/m^2', &
         avgflag='A', long_name='surface forcing of dust in snow, averaged only when snow is present (land)', &
         ptr_pft=clm3%g%l%c%p%pef%sfc_frc_dst_sno, set_lake=spval, set_urb=spval)
#endif

    ! Print masterlist of history fields

    call hist_printflds()

  end subroutine hist_initFlds

end module histFldsMod
