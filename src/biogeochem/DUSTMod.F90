module DUSTMod

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Routines in this module calculate Dust mobilization and dry deposition for dust.
  ! Simulates dust mobilization due to wind from the surface into the 
  ! lowest atmospheric layer. On output flx_mss_vrt_dst(ndst) is the surface dust 
  ! emission (kg/m**2/s) [ + = to atm].
  ! Calculates the turbulent component of dust dry deposition, (the turbulent deposition 
  ! velocity through the lowest atmospheric layer). CAM will calculate the settling 
  ! velocity through the whole atmospheric column. The two calculations will determine 
  ! the dust dry deposition flux to the surface.
  !                              
  ! !USES:
  use shr_kind_mod         , only : r8 => shr_kind_r8 
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use shr_infnan_mod       , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar           , only : dst_src_nbr, ndst, sz_nbr, &
                                    natpft_lb, natpft_ub, natpft_size     ! -dmleung added 24 Jul 2022
  use clm_varcon           , only : grav, spval
  use landunit_varcon      , only : istcrop, istsoil
  use clm_varctl           , only : iulog
  use abortutils           , only : endrun
  use decompMod            , only : bounds_type, subgrid_level_landunit, subgrid_level_patch
  use atm2lndType          , only : atm2lnd_type
  use SoilStateType        , only : soilstate_type
  use CanopyStateType      , only : canopystate_type
  use WaterStateBulkType       , only : waterstatebulk_type
  use WaterDiagnosticBulkType       , only : waterdiagnosticbulk_type
  use FrictionVelocityMod  , only : frictionvel_type
  use LandunitType         , only : lun
  use ColumnType           , only : col
  use PatchType            , only : patch
  !use clm_instur           , only : wt_lunit, wt_nat_patch  ! dmleung added 24 Jul 2022
  !use landunit_varcon      , only : istsoil, istcrop        ! dmleung added 24 Jul 2022 (refering to main/landunit_varcon.F90, for wt_lunit, istsoil=1 is nat veg, istcrop=2 is crop)
  use pftconMod            , only : noveg              ! dmleung added 24 Jul 2022
  use PrigentRoughnessStreamType    , only : prigentroughnessstream_type
  !  
  ! !PUBLIC TYPES
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  !
  public DustEmission   ! Dust mobilization 
  public DustDryDep     ! Turbulent dry deposition for dust
  !
  ! !PUBLIC DATA:
  !
  real(r8) , allocatable :: ovr_src_snk_mss(:,:)
  real(r8) , allocatable :: dmt_vwr(:) ![m] Mass-weighted mean diameter resolved
  real(r8) , allocatable :: stk_crc(:) ![frc] Correction to Stokes settling velocity
  real(r8) tmp1                        !Factor in saltation computation (named as in Charlie's code)
  real(r8) dns_aer                     ![kg m-3] Aerosol density
  !
  ! !PUBLIC DATA TYPES:
  !
  type, public :: dust_type

     real(r8), pointer, PUBLIC  :: flx_mss_vrt_dst_patch     (:,:) ! surface dust emission (kg/m**2/s) [ + = to atm] (ndst) 
     real(r8), pointer, private :: flx_mss_vrt_dst_tot_patch (:)   ! total dust flux into atmosphere
     real(r8), pointer, private :: vlc_trb_patch             (:,:) ! turbulent deposition velocity  (m/s) (ndst) 
     real(r8), pointer, private :: vlc_trb_1_patch           (:)   ! turbulent deposition velocity 1(m/s)
     real(r8), pointer, private :: vlc_trb_2_patch           (:)   ! turbulent deposition velocity 2(m/s)
     real(r8), pointer, private :: vlc_trb_3_patch           (:)   ! turbulent deposition velocity 3(m/s)
     real(r8), pointer, private :: vlc_trb_4_patch           (:)   ! turbulent deposition velocity 4(m/s)
     real(r8), pointer, private :: mbl_bsn_fct_col           (:)   ! basin factor
     !########### added by dmleung 27 Nov 2021 ########################################################################
     real(r8), pointer, private :: dst_emiss_coeff_patch     (:)   ! dust emission coefficient (unitless)
     real(r8), pointer, private :: wnd_frc_thr_patch         (:)   ! wet fluid threshold (m/s)
     real(r8), pointer, private :: wnd_frc_thr_dry_patch     (:)   ! dry fluid threshold (m/s)
     real(r8), pointer, private :: lnd_frc_mble_patch        (:)   ! land mobile fraction -dmleung
     real(r8), pointer, private :: liq_frac_patch            (:)   ! liquid fraction of total water
     real(r8), pointer, private :: wnd_frc_soil_patch        (:)   ! soil wind friction velocity (m/s)
     real(r8), pointer, private :: gwc_patch                 (:)   ! gravimetric water content (kg/kg)
     !########### added by dmleung 2 Dec 2021 #########################################################################
     real(r8), pointer, private :: intrmtncy_fct_patch       (:)   ! intermittency factor, accounting for turbulence shutting down dust emissions (unitless)
     real(r8), pointer, private :: stblty_patch              (:)   ! stability parameter for checking stability condition (stblty < 0 is unstable atmosphere)
     real(r8), pointer, private :: u_mean_slt_patch          (:)   ! wind speed 0.1 m level of dust saltation (m/s)
     real(r8), pointer, private :: u_sd_slt_patch            (:)   ! sd of wind speed 0.1 m level of dust saltation (m/s)
     real(r8), pointer, private :: u_fld_thr_patch           (:)   ! fluid threshold wind speed 0.1 m level of dust saltation (m/s)
     real(r8), pointer, private :: u_impct_thr_patch         (:)   ! impact threshold wind speed at 0.1 m level of dust saltation (m/s)
     real(r8), pointer, private :: thr_crs_rate_patch        (:)   ! threshold crossing rate (unitless)
     real(r8), pointer, private :: prb_crs_fld_thr_patch     (:)   ! probability of wind speed crossing fluid threshold
     real(r8), pointer, private :: prb_crs_impct_thr_patch   (:)   ! probability of wind speed crossing impact threshold
     !########### added by dmleung 20 Dec 2021 ########################################################################
     real(r8), pointer, private :: ssr_patch                 (:)   ! [dimless] integrated shear stress ratiio, defined by Okin (2008) and then integrated by Caroline Pierre et al. (2014)
     real(r8), pointer, private :: lai_patch                 (:)   ! [m2 leaf /m2 land] LAI+SAI for calculating Okin's drag partition, averaged to landunit level
     real(r8), pointer, private :: frc_thr_rghn_fct_patch    (:)   ! [dimless] hybrid drag partition (or called roughness) factor
     !########### added by dmleung 28 Jul 2022 ########################################################################
     real(r8), pointer, private :: wnd_frc_thr_std_patch     (:)   ! standardized fluid threshold friction velocity (m/s)
     !########### added by dmleung 31 Dec 2022 ########################################################################
     type(prigentroughnessstream_type), private :: prigentroughnessstream      ! Prigent roughness stream data
     real(r8), pointer, private :: dpfct_rock_patch          (:)   ! [fraction] rock drag partition factor, time-constant
   contains

     procedure , public  :: Init
     procedure , private :: InitAllocate 
     procedure , private :: InitHistory  
     procedure , private :: InitCold     
     procedure , private :: InitDustVars ! Initialize variables used in subroutine Dust

  end type dust_type
  !------------------------------------------------------------------------

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !------------------------------------------------------------------------
  !##### dmleung edited for initializing stream files 31 Dec 2022  ########
  !subroutine Init(this, bounds, NLFilename, num_nolakep, filter_nolakep)
  subroutine Init(this, bounds, NLFilename)

    class(dust_type) :: this
    type(bounds_type), intent(in) :: bounds 
    character(len=*),  intent(in) :: NLFilename   ! dmleung added 31 Dec 2022
    !integer                , intent(in)    :: num_nolakep                 ! number of column non-lake points in patch filter
    !integer                , intent(in)    :: filter_nolakep(num_nolakep) ! patch filter for non-lake points  

    call this%InitAllocate (bounds)
    call this%InitHistory  (bounds)
    call this%prigentroughnessstream%Init( bounds, NLFilename )  ! dmleung added 31 Dec 2022
    !call this%InitCold     (bounds, num_nolakep, filter_nolakep)
    call this%InitCold     (bounds)                             ! dmleung commented 31 Dec 2022
    call this%InitDustVars (bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class (dust_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp,endp
    integer :: begc,endc
    !------------------------------------------------------------------------

    begp = bounds%begp ; endp = bounds%endp
    begc = bounds%begc ; endc = bounds%endc

    allocate(this%flx_mss_vrt_dst_patch     (begp:endp,1:ndst)) ; this%flx_mss_vrt_dst_patch     (:,:) = nan
    allocate(this%flx_mss_vrt_dst_tot_patch (begp:endp))        ; this%flx_mss_vrt_dst_tot_patch (:)   = nan
    allocate(this%vlc_trb_patch             (begp:endp,1:ndst)) ; this%vlc_trb_patch             (:,:) = nan
    allocate(this%vlc_trb_1_patch           (begp:endp))        ; this%vlc_trb_1_patch           (:)   = nan
    allocate(this%vlc_trb_2_patch           (begp:endp))        ; this%vlc_trb_2_patch           (:)   = nan 
    allocate(this%vlc_trb_3_patch           (begp:endp))        ; this%vlc_trb_3_patch           (:)   = nan
    allocate(this%vlc_trb_4_patch           (begp:endp))        ; this%vlc_trb_4_patch           (:)   = nan
    allocate(this%mbl_bsn_fct_col           (begc:endc))        ; this%mbl_bsn_fct_col     (:)   = nan
    !#### added by dmleung 27 Nov 2021 #####################################
    allocate(this%dst_emiss_coeff_patch     (begp:endp))        ; this%dst_emiss_coeff_patch     (:)   = nan
    allocate(this%wnd_frc_thr_patch         (begp:endp))        ; this%wnd_frc_thr_patch         (:)   = nan
    allocate(this%wnd_frc_thr_dry_patch     (begp:endp))        ; this%wnd_frc_thr_dry_patch     (:)   = nan
    allocate(this%lnd_frc_mble_patch        (begp:endp))        ; this%lnd_frc_mble_patch        (:)   = nan
    allocate(this%wnd_frc_soil_patch        (begp:endp))        ; this%wnd_frc_soil_patch        (:)   = nan
    allocate(this%gwc_patch                 (begp:endp))        ; this%gwc_patch                 (:)   = nan
    allocate(this%liq_frac_patch            (begp:endp))        ; this%liq_frac_patch            (:)   = nan
    !#### added by dmleung 2 Dec 2021 ######################################
    allocate(this%intrmtncy_fct_patch       (begp:endp))        ; this%intrmtncy_fct_patch       (:)   = nan
    allocate(this%stblty_patch              (begp:endp))        ; this%stblty_patch              (:)   = nan
    allocate(this%u_mean_slt_patch          (begp:endp))        ; this%u_mean_slt_patch          (:)   = nan
    allocate(this%u_sd_slt_patch            (begp:endp))        ; this%u_sd_slt_patch            (:)   = nan
    allocate(this%u_fld_thr_patch           (begp:endp))        ; this%u_fld_thr_patch           (:)   = nan
    allocate(this%u_impct_thr_patch         (begp:endp))        ; this%u_impct_thr_patch         (:)   = nan
    allocate(this%thr_crs_rate_patch        (begp:endp))        ; this%thr_crs_rate_patch        (:)   = nan
    allocate(this%prb_crs_fld_thr_patch     (begp:endp))        ; this%prb_crs_fld_thr_patch     (:)   = nan
    allocate(this%prb_crs_impct_thr_patch   (begp:endp))        ; this%prb_crs_impct_thr_patch   (:)   = nan
    !#### added by dmleung 17 Dec 2021 ######################################
    allocate(this%ssr_patch                 (begp:endp))        ; this%ssr_patch                 (:)   = nan
    allocate(this%lai_patch                 (begp:endp))        ; this%lai_patch                 (:)   = nan
    allocate(this%frc_thr_rghn_fct_patch    (begp:endp))        ; this%frc_thr_rghn_fct_patch    (:)   = nan
    !#### added by dmleung 28 Jul 2022 ######################################
    allocate(this%wnd_frc_thr_std_patch     (begp:endp))        ; this%wnd_frc_thr_std_patch     (:)   = nan
    !#### added by dmleung 31 Dec 2022 ######################################
    allocate(this%dpfct_rock_patch          (begp:endp))        ; this%dpfct_rock_patch          (:)   = nan
  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d
    !
    !
    ! !ARGUMENTS:
    class (dust_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp,endp
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    this%flx_mss_vrt_dst_tot_patch(begp:endp) = spval
    call hist_addfld1d (fname='DSTFLXT', units='kg/m2/s',  &
         avgflag='A', long_name='total surface dust emission', &
         ptr_patch=this%flx_mss_vrt_dst_tot_patch, set_lake=0._r8, set_urb=0._r8)

    this%vlc_trb_1_patch(begp:endp) = spval
    call hist_addfld1d (fname='DPVLTRB1', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 1', &
         ptr_patch=this%vlc_trb_1_patch, default='inactive')

    this%vlc_trb_2_patch(begp:endp) = spval
    call hist_addfld1d (fname='DPVLTRB2', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 2', &
         ptr_patch=this%vlc_trb_2_patch, default='inactive')

    this%vlc_trb_3_patch(begp:endp) = spval
    call hist_addfld1d (fname='DPVLTRB3', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 3', &
         ptr_patch=this%vlc_trb_3_patch, default='inactive')

    this%vlc_trb_4_patch(begp:endp) = spval
    call hist_addfld1d (fname='DPVLTRB4', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 4', &
         ptr_patch=this%vlc_trb_4_patch, default='inactive')

    !#####added by dmleung 27 Nov 2021#########################################
    this%dst_emiss_coeff_patch(begp:endp) = spval
    call hist_addfld1d (fname='C_d', units='dimensionless',  &
         avgflag='A', long_name='dust emission coefficient', &
         ptr_patch=this%dst_emiss_coeff_patch, set_lake=0._r8, set_urb=0._r8)
    this%wnd_frc_thr_patch(begp:endp) = spval
    call hist_addfld1d (fname='WND_FRC_FT', units='m/s',  &
         avgflag='A', long_name='fluid threshold friction velocity', &
         ptr_patch=this%wnd_frc_thr_patch, set_lake=0._r8, set_urb=0._r8)
    this%wnd_frc_thr_dry_patch(begp:endp) = spval
    call hist_addfld1d (fname='WND_FRC_FT_DRY', units='m/s',  &
         avgflag='A', long_name='dry fluid threshold friction velocity', &
         ptr_patch=this%wnd_frc_thr_dry_patch, set_lake=0._r8, set_urb=0._r8)
    this%wnd_frc_soil_patch(begp:endp) = spval
    call hist_addfld1d (fname='WND_FRC_SOIL', units='m/s',  &
         avgflag='A', long_name='soil surface wind friction velocity', &
         ptr_patch=this%wnd_frc_soil_patch, set_lake=0._r8, set_urb=0._r8)
    this%lnd_frc_mble_patch(begp:endp) = spval
    call hist_addfld1d (fname='LND_FRC_MBLE', units='dimensionless',  &
         avgflag='A', long_name='land mobile fraction', &
         ptr_patch=this%lnd_frc_mble_patch, set_lake=0._r8, set_urb=0._r8)
    this%gwc_patch(begp:endp) = spval
    call hist_addfld1d (fname='GWC', units='kg/kg',  &
         avgflag='A', long_name='gravimetric water content', &
         ptr_patch=this%gwc_patch, set_lake=0._r8, set_urb=0._r8)
    this%liq_frac_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIQ_FRAC', units='dimensionless',  &
         avgflag='A', long_name='fraction of total water that is liquid', &
         ptr_patch=this%liq_frac_patch, set_lake=0._r8, set_urb=0._r8)
    !#####added by dmleung 2 Dec 2021 #########################################
    this%u_mean_slt_patch(begp:endp) = spval
    call hist_addfld1d (fname='U_S_MEAN', units='m/s',  &
         avgflag='A', long_name='mean wind velocity at saltation level', &
         ptr_patch=this%u_mean_slt_patch, set_lake=0._r8, set_urb=0._r8)
    this%u_sd_slt_patch(begp:endp) = spval
    call hist_addfld1d (fname='U_S_SIGMA', units='m/s',  &
         avgflag='A', long_name='sd of wind velocity at saltation level', &
         ptr_patch=this%u_sd_slt_patch, set_lake=0._r8, set_urb=0._r8)
    this%stblty_patch(begp:endp) = spval
    call hist_addfld1d (fname='ZETAOBU', units='',  &
         avgflag='A', long_name='stability parameter', &
         ptr_patch=this%stblty_patch, set_lake=0._r8, set_urb=0._r8)
    this%u_fld_thr_patch(begp:endp) = spval
    call hist_addfld1d (fname='U_FT', units='m/s',  &
         avgflag='A', long_name='fluid threshold velocity at saltation level', &
         ptr_patch=this%u_fld_thr_patch, set_lake=0._r8, set_urb=0._r8)
    this%u_impct_thr_patch(begp:endp) = spval
    call hist_addfld1d (fname='U_IT', units='m/s',  &
         avgflag='A', long_name='impact threshold velocity at saltation level', &
         ptr_patch=this%u_impct_thr_patch, set_lake=0._r8, set_urb=0._r8)
    this%thr_crs_rate_patch(begp:endp) = spval
    call hist_addfld1d (fname='ALPHA_TC_RATE', units='',  &
         avgflag='A', long_name='threshold crossing rate', &
         ptr_patch=this%thr_crs_rate_patch, set_lake=0._r8, set_urb=0._r8)
    this%prb_crs_fld_thr_patch(begp:endp) = spval
    call hist_addfld1d (fname='P_FT', units='',  &
         avgflag='A', long_name='probability of winds crossing fluid threshold', &
         ptr_patch=this%prb_crs_fld_thr_patch, set_lake=0._r8, set_urb=0._r8)
    this%prb_crs_impct_thr_patch(begp:endp) = spval
    call hist_addfld1d (fname='P_IT', units='',  &
         avgflag='A', long_name='probability of winds crossing impact threshold', &
         ptr_patch=this%prb_crs_impct_thr_patch, set_lake=0._r8, set_urb=0._r8)
    this%intrmtncy_fct_patch(begp:endp) = spval
    call hist_addfld1d (fname='ETA', units='',  &
         avgflag='A', long_name='intermittency factor', &
         ptr_patch=this%intrmtncy_fct_patch, set_lake=0._r8, set_urb=0._r8)
    !#####added by dmleung 20 Dec 2021 ########################################
    this%ssr_patch(begp:endp) = spval
    call hist_addfld1d (fname='SSR', units='m/s',  &
         avgflag='A', long_name='Okin-Pierre vegetation shear stress ratio (drag partition factor)', &
         ptr_patch=this%ssr_patch, set_lake=0._r8, set_urb=0._r8)
    this%lai_patch(begp:endp) = spval
    call hist_addfld1d (fname='LAI', units='m/s',  &
         avgflag='A', long_name='landunit-mean LAI for Okin-Pierre scheme', &
         ptr_patch=this%lai_patch, set_lake=0._r8, set_urb=0._r8)
    this%frc_thr_rghn_fct_patch(begp:endp) = spval
    call hist_addfld1d (fname='FRC_THR_RGHN_FCT', units='dimensionless',  &
         avgflag='A', long_name='hybrid drag partition (or roughness) factor', &
         ptr_patch=this%frc_thr_rghn_fct_patch, set_lake=0._r8, set_urb=0._r8)
    !#####added by dmleung 28 Jul 2022 ########################################
    this%wnd_frc_thr_std_patch(begp:endp) = spval
    call hist_addfld1d (fname='WND_FRC_FT_STD', units='m/s',  &
         avgflag='A', long_name='standardized fluid threshold friction velocity', &
         ptr_patch=this%wnd_frc_thr_std_patch, set_lake=0._r8, set_urb=0._r8)
    !#####added by dmleung 31 Dec 2022 ########################################
    this%dpfct_rock_patch(begp:endp) = spval
    call hist_addfld1d (fname='DPFCT_ROCK', units='m/s',  &
         avgflag='A', long_name='rock drag partition factor', &
         ptr_patch=this%dpfct_rock_patch)
    !##########################################################################

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  !subroutine InitCold(this, bounds, num_nolakep, filter_nolakep)   !dmleung commented 31 Dec 2022
  subroutine InitCold(this, bounds)
    !
    !USES dmleung added 31 Dec 2022
    !use landunit_varcon      , only : istdlak
    !use LandunitType         , only : lun
    !
    ! !ARGUMENTS:
    class(dust_type), intent(inout) :: this  ! dmleung used class instead of type, 31 Dec 2022
    type(bounds_type), intent(in) :: bounds 
    ! dmleung also added the below no-lake filter 31 Dec 2022
    !integer                , intent(in)    :: num_nolakep                 ! number of column non-lake points in patch filter
    !integer                , intent(in)    :: filter_nolakep(num_nolakep) ! patch filter for non-lake points 
    !type(dust_type), intent(inout) :: dust_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c,l
    !-----------------------------------------------------------------------

    ! Set basin factor to 1 for now

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)

       if (.not.lun%lakpoi(l)) then
          this%mbl_bsn_fct_col(c) = 1.0_r8
       end if
    end do

    !associate(                                                                 &
    !     dpfct_rock                 =>   this%dpfct_rock_patch                 &
    !     )
      ! Caulculate Drag Partition factor, dmleung added 31 Dec 2022
      if ( this%prigentroughnessstream%useStreams() )then !if usestreams == true, and it should be always true
         call this%prigentroughnessstream%CalcDragPartition( bounds, &
         !                       num_nolakep, filter_nolakep, this%dpfct_rock_patch(bounds%begp:bounds%endp) )
                                this%dpfct_rock_patch(bounds%begp:bounds%endp) )
      else

         call endrun( "ERROR:: Drag partitioning MUST now use a streams file of aeolian roughness length to calculate, it can no longer read from the fsurdat file" )
      end if
    !end associate

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine DustEmission (bounds, &
       num_nolakep, filter_nolakep, &
       atm2lnd_inst, soilstate_inst, canopystate_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, &
       frictionvel_inst, dust_inst)
    !
    ! !DESCRIPTION: 
    ! Dust mobilization. This code simulates dust mobilization due to wind
    ! from the surface into the lowest atmospheric layer
    ! On output flx_mss_vrt_dst(ndst) is the surface dust emission 
    ! (kg/m**2/s) [ + = to atm]
    ! Source: C. Zender's dust model
    !
    ! !USES
    use shr_const_mod, only : SHR_CONST_RHOFW
    use subgridaveMod, only : p2g
    !use clm_instur           , only : wt_lunit, wt_nat_patch  ! dmleung added 24 Jul 2022
    !use landunit_varcon      , only : istsoil, istcrop        ! dmleung added 24 Jul 2022 (refering to main/landunit_varcon.F90, for wt_lunit, istsoil=1 is nat veg, istcrop=2 is crop)
    use pftconMod            , only : noveg              ! dmleung added 24 Jul 2022
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds                      
    integer                , intent(in)    :: num_nolakep                 ! number of column non-lake points in patch filter
    integer                , intent(in)    :: filter_nolakep(num_nolakep) ! patch filter for non-lake points
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(soilstate_type)   , intent(in)    :: soilstate_inst
    type(canopystate_type) , intent(in)    :: canopystate_inst
    type(waterstatebulk_type)  , intent(in)    :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)  , intent(in)    :: waterdiagnosticbulk_inst
    type(frictionvel_type) , intent(in)    :: frictionvel_inst
    type(dust_type)        , intent(inout) :: dust_inst

    !
    ! !LOCAL VARIABLES
    integer  :: fp,p,c,l,g,m,n      ! indices
    real(r8) :: liqfrac             ! fraction of total water that is liquid
    real(r8) :: wnd_frc_rat         ! [frc] Wind friction threshold over wind friction
    real(r8) :: wnd_frc_slt_dlt     ! [m s-1] Friction velocity increase from saltatn
    real(r8) :: wnd_rfr_dlt         ! [m s-1] Reference windspeed excess over threshld
    real(r8) :: dst_slt_flx_rat_ttl
    real(r8) :: flx_mss_hrz_slt_ttl
    real(r8) :: flx_mss_vrt_dst_ttl(bounds%begp:bounds%endp)
    real(r8) :: frc_thr_wet_fct
    real(r8) :: frc_thr_rgh_fct
    !real(r8) :: wnd_frc_thr_slt    ! dmleung commented and put below 2 Dec 2021
    real(r8) :: wnd_rfr_thr_slt
    real(r8) :: wnd_frc_slt
    real(r8) :: lnd_frc_mbl(bounds%begp:bounds%endp)
    real(r8) :: bd
    real(r8) :: gwc_sfc
    real(r8) :: ttlai(bounds%begp:bounds%endp)
    real(r8) :: tlai_lu(bounds%begl:bounds%endl)
    real(r8) :: sumwt(bounds%begl:bounds%endl) ! sum of weights
    logical  :: found                          ! temporary for error check
    integer  :: index
    !########### added by dmleung 27 Nov 2021 #########################################
    real(r8) :: tmp2   ! calculates the dry fluid threshold using Shao and Lu (2000) scheme; replace the tmp1 (Iversen and White, 1982) that was passed from Dustini to DustEmission; tmp2 will be calculated here  23 May 2020 -dmleung
    real(r8) :: wnd_frc_thr_slt_std ! [m/s] The soil threshold friction speed at standard air density (1.2250 kg/m3) -jfk
    real(r8) :: frag_expt           ! fragmentation exponent, -dmleung 22 Jun 2021
    !########### added by dmleung 2 Dec 2021 for intermittency scheme #################
    real(r8) :: wnd_frc_thr_slt_it  ! [m/s] created for impact threshold friction velocity, dmleung 9 Jun 2021
    real(r8) :: wnd_frc_thr_slt     ! [m/s] used for wet fluid threshold friction velocity, dmleung 9 Jun 2021
    !########### added by dmleung 20 Dec 2021 for drag partition effect #################
    real(r8) :: K_length            ! [dimless] normalized mean interobstacle distance, or called gap length (Okin, 2008)
    !########### added by dmleung 22 Jul 2022 for LUH2 land cover ####################
    real(r8) :: bare_frc            ! LUH2 bare soil land cover fraction
    real(r8) :: veg_frc             ! LUH2 natural vegetation + crop land cover fraction
    !    
    ! constants
    !
    real(r8), parameter :: cst_slt = 2.61_r8           ! [frc] Saltation constant
    real(r8), parameter :: flx_mss_fdg_fct = 5.0e-4_r8 ! [frc] Empir. mass flx tuning eflx_lh_vegt
    !real(r8), parameter :: vai_mbl_thr = 0.3_r8        ! [m2 m-2] VAI threshold quenching dust mobilization
    !####### added by dmleung 27 Nov 2021 ###########################################################################
    character(len=*),parameter :: subname = 'DUSTEmission'
    real(r8), parameter :: vai_mbl_thr = 1.0_r8        ! [m2 m-2] new VAI threshold; dmleung suggests 1 or 0.5, and the default 0.3 seems a bit too small -dmleung 27 Nov 2021
    real(r8), parameter :: Cd0 = 4.4e-5_r8             ! [dimless] proportionality constant in calculation of dust emission coefficient -jfk
    real(r8), parameter :: Ca = 2.7_r8                 ! [dimless] proportionality constant in scaling of dust emission exponent -jfk
    real(r8), parameter :: Ce = 2.0_r8                 ! [dimless] proportionality constant scaling exponential dependence of dust emission coefficient on standardized soil threshold friction speed -jfk
    real(r8), parameter :: C_tune = 0.05_r8             ! [dimless] global tuning constant for vertical dust flux; set to produce ~same global dust flux in control sim (I_2000) as old parameterization -jfk
    real(r8), parameter :: wnd_frc_thr_slt_std_min = 0.16_r8 ! [m/s] minimum standardized soil threshold friction speed -jfk
    real(r8), parameter :: forc_rho_std = 1.2250_r8    ! [kg/m3] density of air at standard pressure (101325) and temperature (293 K) -jfk
    real(r8), parameter :: dns_slt = 2650.0_r8         ! [kg m-3] Density of optimal saltation particles, dml 23 May 2020
    !####### added by dmleung 2 Dec 2021 for intermittency ##########################################################
    real(r8), parameter :: B_it = 0.82_r8              ! [dimless] ratio = u_star_it / u_star_ft0 (may need to change into a fn of moisture later on) -dml
    real(r8), parameter :: k = 0.4_r8                  ! [dimless] von Karman constant -dml
    !####### added by dmleung 2 Dec 2021 for Okin (2008) drag partition for plants ##########################################################
    real(r8), parameter :: f_0 = 0.32_r8               ! [dimless] SSR in the immediate lee of a plant, dimensionless
    real(r8), parameter :: c_e = 4.8_r8                  ! [dimless] e-folding distance velocity recovery, dimensionless
    !################################################################################################################
    !################################################################################################################
    !------------------------------------------------------------------------

    associate(                                                         & 
         forc_rho            => atm2lnd_inst%forc_rho_downscaled_col , & ! Input:  [real(r8) (:)   ]  downscaled density (kg/m**3)                                 
         
         gwc_thr             => soilstate_inst%gwc_thr_col           , & ! Input:  [real(r8) (:)   ]  threshold gravimetric soil moisture based on clay content
         mss_frc_cly_vld     => soilstate_inst%mss_frc_cly_vld_col   , & ! Input:  [real(r8) (:)   ]  [frc] Mass fraction clay limited to 0.20          
         watsat              => soilstate_inst%watsat_col            , & ! Input:  [real(r8) (:,:) ]  saturated volumetric soil water                 
         
         tlai                => canopystate_inst%tlai_patch          , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index, no burying by snow     
         tsai                => canopystate_inst%tsai_patch          , & ! Input:  [real(r8) (:)   ]  one-sided stem area index, no burying by snow     
         
         frac_sno            => waterdiagnosticbulk_inst%frac_sno_col         , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)       
         h2osoi_vol          => waterstatebulk_inst%h2osoi_vol_col       , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat)   
         h2osoi_liq          => waterstatebulk_inst%h2osoi_liq_col       , & ! Input:  [real(r8) (:,:) ]  liquid soil water (kg/m2)                       
         h2osoi_ice          => waterstatebulk_inst%h2osoi_ice_col       , & ! Input:  [real(r8) (:,:) ]  frozen soil water (kg/m2)                       
         
         fv                  => frictionvel_inst%fv_patch            , & ! Input:  [real(r8) (:)   ]  friction velocity (m/s) (for dust model)          
         u10                 => frictionvel_inst%u10_patch           , & ! Input:  [real(r8) (:)   ]  10-m wind (m/s) (created for dust model)          
         
         mbl_bsn_fct         => dust_inst%mbl_bsn_fct_col            , & ! Input:  [real(r8) (:)   ]  basin factor                                      
         flx_mss_vrt_dst     => dust_inst%flx_mss_vrt_dst_patch      , & ! Output: [real(r8) (:,:) ]  surface dust emission (kg/m**2/s)               
         flx_mss_vrt_dst_tot => dust_inst%flx_mss_vrt_dst_tot_patch  , & ! Output: [real(r8) (:)   ]  total dust flux back to atmosphere (pft)
         ! the following are added by dmleung 27 Nov 2021
         dst_emiss_coeff     => dust_inst%dst_emiss_coeff_patch      , & ! Output dust emission coefficient
         wnd_frc_thr         => dust_inst%wnd_frc_thr_patch          , & ! output impact threshold -dmleung
         wnd_frc_thr_dry     => dust_inst%wnd_frc_thr_dry_patch      , & ! output dry threshold
         lnd_frc_mble        => dust_inst%lnd_frc_mble_patch         , & ! -dmleung, 3 Feb 2020
         wnd_frc_soil        => dust_inst%wnd_frc_soil_patch         , & ! soil friction velocity u_*s = (u_*)(f_eff)
         gwc                 => dust_inst%gwc_patch                  , & ! output gravimetric water content
         liq_frac            => dust_inst%liq_frac_patch             , &
         ! added by dmleung 8 Jul 2019, recoded 2 Dec 2021
         intrmtncy_fct       => dust_inst%intrmtncy_fct_patch        , &
         stblty              => dust_inst%stblty_patch               , &
         u_mean_slt          => dust_inst%u_mean_slt_patch           , &
         u_sd_slt            => dust_inst%u_sd_slt_patch             , &
         u_fld_thr           => dust_inst%u_fld_thr_patch            , &
         u_impct_thr         => dust_inst%u_impct_thr_patch          , &
         thr_crs_rate        => dust_inst%thr_crs_rate_patch         , &
         prb_crs_fld_thr     => dust_inst%prb_crs_fld_thr_patch      , &
         prb_crs_impct_thr   => dust_inst%prb_crs_impct_thr_patch    , &
         ! added by dmleung 17 Dec 2021
         roughfct            => soilstate_inst%roughfct_patch        , &  ! dmleung replaced it by dpfct_rock below, 31 Dec 2022
         ! added by dmleung 20 Dec 2021
         ssr                 => dust_inst%ssr_patch                  , &
         lai                 => dust_inst%lai_patch                  , &
         frc_thr_rghn_fct    => dust_inst%frc_thr_rghn_fct_patch     , &
         ! added by dmleung 28 Jul 2022
         wnd_frc_thr_std     => dust_inst%wnd_frc_thr_std_patch      , &
         ! added by dmleung 31 Dec 2022
         dpfct_rock          => dust_inst%dpfct_rock_patch            &  ! dmleung used roughfct (roughness factor) instead of dpfct_rock (rock drag partition factor) here. Could change it back to dpfct_rock here and below later, 31 Dec 2022
         )

      ttlai(bounds%begp : bounds%endp) = 0._r8
      ! make lai average at landunit level
      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         ttlai(p) = tlai(p)+tsai(p)
      enddo

      tlai_lu(bounds%begl : bounds%endl) = spval
      sumwt(bounds%begl : bounds%endl) = 0._r8
      do p = bounds%begp,bounds%endp
         if (ttlai(p) /= spval .and. patch%active(p) .and. patch%wtlunit(p) /= 0._r8) then
            c = patch%column(p)
            l = patch%landunit(p)
            if (sumwt(l) == 0._r8) tlai_lu(l) = 0._r8
            tlai_lu(l) = tlai_lu(l) + ttlai(p) * patch%wtlunit(p)
            sumwt(l) = sumwt(l) + patch%wtlunit(p)
         end if
      end do
      found = .false.
      do l = bounds%begl,bounds%endl
         if (sumwt(l) > 1.0_r8 + 1.e-6_r8) then
            found = .true.
            index = l
            exit
         else if (sumwt(l) /= 0._r8) then
            tlai_lu(l) = tlai_lu(l)/sumwt(l)
         end if
      end do
      if (found) then
         write(iulog,*) subname//':: error: sumwt is greater than 1.0 at l= ',index
         call endrun(subgrid_index=index, subgrid_level=subgrid_level_landunit, msg=errMsg(sourcefile, __LINE__))
      end if

      ! Loop through patches

      ! initialize variables which get passed to the atmosphere
      flx_mss_vrt_dst(bounds%begp:bounds%endp,:)=0._r8

      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         c = patch%column(p)
         l = patch%landunit(p)

         ! the following code from subr. lnd_frc_mbl_get was adapted for lsm use
         ! purpose: return fraction of each gridcell suitable for dust mobilization

         ! the "bare ground" fraction of the current sub-gridscale cell decreases
         ! linearly from 1 to 0 as VAI(=tlai+tsai) increases from 0 to vai_mbl_thr
         ! if ice sheet, wetland, or lake, no dust allowed

         if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
            if (tlai_lu(l) < vai_mbl_thr) then
               lnd_frc_mbl(p) = 1.0_r8 - (tlai_lu(l))/vai_mbl_thr
            else
               lnd_frc_mbl(p) = 0.0_r8
            endif
            lnd_frc_mbl(p) = lnd_frc_mbl(p) * (1.0_r8 - frac_sno(c))
         else          
            lnd_frc_mbl(p) = 0.0_r8   
         end if
      end do

      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         if (lnd_frc_mbl(p)>1.0_r8 .or. lnd_frc_mbl(p)<0.0_r8) then
            write(iulog,*)'Error dstmbl: pft= ',p,' lnd_frc_mbl(p)= ',lnd_frc_mbl(p)
            call endrun(subgrid_index=p, subgrid_level=subgrid_level_patch, msg=errMsg(sourcefile, __LINE__))
         end if
      end do

      ! dmleung add output for bare_frc and veg_frc here if wanted !!!!!!!!!----------------------

      ! reset history output variables before next if-statement to avoid output = inf

      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         flx_mss_vrt_dst_tot(p) = 0.0_r8
         ! the following are added by dmleung 27 Nov 2021
         dst_emiss_coeff(p) = 0.0_r8
         wnd_frc_thr(p) = 0.0_r8
         wnd_frc_thr_dry(p) = 0.0_r8
         lnd_frc_mble(p) = 0.0_r8
         wnd_frc_soil(p) = 0.0_r8
         gwc(p) = 0.0_r8
         liq_frac(p) = 0.0_r8
         ! dmleung's edit, 8 Jul 2019; added by dmleung 2 Dec 2021
         u_mean_slt(p) = 0.0_r8
         u_sd_slt(p) = 0.0_r8
         stblty(p)   = 0.0_r8
         u_fld_thr(p) = 0.0_r8
         u_impct_thr(p) = 0.0_r8
         thr_crs_rate(p) = 0.0_r8
         prb_crs_fld_thr(p) = 0.0_r8
         prb_crs_impct_thr(p) = 0.0_r8
         intrmtncy_fct(p) = 0.0_r8
         ! dmleung's edit, 20 Dec 2021
         ssr(p) = 0.0_r8
         lai(p) = 0.0_r8
         frc_thr_rghn_fct(p) = 0.0_r8
         ! dmleung added 28 Jul 2022
         wnd_frc_thr_std(p) = 0.0_r8
      end do
      do n = 1, ndst
         do fp = 1,num_nolakep
            p = filter_nolakep(fp)
            flx_mss_vrt_dst(p,n) = 0.0_r8
         end do
      end do

      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         c = patch%column(p)
         l = patch%landunit(p)
         g = patch%gridcell(p)

         !################################################################################################
         ! put dust emission calculation here to output threshold friction velocity for the whole globe,
         ! not just when lnd_frc_mbl = 0. Edited by dmleung 27 Nov 2021
         bd = (1._r8-watsat(c,1))*2.7e3_r8      ![kg m-3] Bulk density of dry surface soil
         gwc_sfc = h2osoi_vol(c,1)*SHR_CONST_RHOFW/bd    ![kg kg-1] Gravimetric H2O cont
         if (gwc_sfc > gwc_thr(c)) then
            frc_thr_wet_fct = sqrt(1.0_r8 + 1.21_r8 * (100.0_r8*(gwc_sfc - gwc_thr(c)))**0.68_r8)
         else
            frc_thr_wet_fct = 1.0_r8
         end if

         ! output moisture variables -dmleung, coded Jul 2020, recoded 18 Mar 2021, added to CLM5 27 Nov 2021
         gwc(p) = gwc_sfc     ! output surface gravimetric water content

         ! slevis: adding liqfrac here, because related to effects from soil water
         liqfrac = max( 0.0_r8, min( 1.0_r8, h2osoi_liq(c,1) / (h2osoi_ice(c,1)+h2osoi_liq(c,1)+1.0e-6_r8) ) ) !-dmleung 27 Nov 2021
         ! output liquid fraction -dmleung 27 Nov 2021
         liq_frac(p) = liqfrac

         !#######################################################################################################
         ! calculate Shao & Lu (2000) dust emission threshold scheme here
         ! use tmp1 from DUSTini for Iversen and White I&W (1982) (75 um is optimal); use tmp2 for S&L (2000) (80 um is optimal)
         ! recoded to CLM5 27 Nov 2021
         !#######################################################################################################

         tmp2 = 1.0_r8*sqrt(0.0123_r8 * (dns_slt*grav*130.0e-6_r8 + 1.65e-4_r8/130.0e-6_r8)) ! calculate S&L (2000) scheme here for threshold; gamma = 1.65e-4 following S&L00, D_p = 127 um ~ 130 um following dmleung's dust paper. As this is a global constant, this line can be put outside the loop to save computational power.
         wnd_frc_thr_dry(p) = tmp2 / sqrt(forc_rho(c))    ! output dry fluid threshold
         wnd_frc_thr_slt = tmp2 / sqrt(forc_rho(c)) * frc_thr_wet_fct !* frc_thr_rgh_fct   ! use as threshold in this module
         wnd_frc_thr_slt_it = B_it * tmp2 / sqrt(forc_rho(c)) ! define impact threshold -dml 9 Jun 2021, recoded to CLM5 27 Nov 2021

         !wnd_frc_thr_dry(p) = tmp1 / sqrt(forc_rho(c))    ! output dry fluid threshold         
         !wnd_frc_thr_slt = tmp1 / sqrt(forc_rho(c)) * frc_thr_wet_fct !* frc_thr_rgh_fct   ! use as threshold in this module
         !wnd_frc_thr_slt_it = B_it * tmp1 / sqrt(forc_rho(c)) ! define impact threshold -dmleung 9 Jun 2021, recoded to CLM5 27 Nov 2021
         ! the above formula is true for Iversen and White (1982) and Shao and Lu (2000) scheme -dmleung, 23 Feb 2020, added to CLM5 27 Nov 2021
         wnd_frc_thr(p) = wnd_frc_thr_slt          ! output fluid threshold -dmleung

         ! use emission threshold to calculate standardized threshold and dust emission coefficient dmleung 27 Nov 2021
         wnd_frc_thr_slt_std = wnd_frc_thr_slt * sqrt(forc_rho(c) / forc_rho_std) ! standardized soil threshold friction speed -jfk (defined using fluid threshold
         wnd_frc_thr_std(p) = wnd_frc_thr_slt_std          ! output standardized fluid threshold -dmleung added 28 Jul 2022
         dst_emiss_coeff(p) = Cd0 * exp(-Ce * (wnd_frc_thr_slt_std - wnd_frc_thr_slt_std_min) / wnd_frc_thr_slt_std_min) ! save dust emission coefficient here for all grids, -dml, 1 Mar 2021 

         ! framentation exponent dmleung 27 Nov 2021; moved to this block 23 Dec 2021
         frag_expt = (Ca * (wnd_frc_thr_slt_std - wnd_frc_thr_slt_std_min) / wnd_frc_thr_slt_std_min)  ! fragmentation exponent, defined in Kok et al. (2014a) -dmleung 27 Nov 2021
         if (frag_expt > 5.0_r8) then   ! set fragmentation exponent to be 3 or 5 at maximum, to avoid local AOD blowup
            frag_expt = 5.0_r8
         end if

         !################ drag partition effect, and soil friction velocity############################
         ! subsection on computing vegetation drag partition and hybrid drag partition factors 
         ! in our scheme, drag partition effect is applied on the wind instead of the threshold
         !  -dmleung, 7 Jul 2021 , coded to CLM5 27 Nov 2021
         !##############################################################################################
         ! the following comes from subr. frc_thr_rgh_fct_get
         ! purpose: compute factor by which surface roughness increases threshold
         !          friction velocity (currently a constant)

         if (lnd_frc_mbl(p) > 0.0_r8  .AND. tlai_lu(l)<=1_r8) then
            ! vegetation drag partition equation following Gregory Okin (2008) + Caroline Pierre et al. (2014), dmleung 20 Dec 2021
            lai(p) = tlai_lu(l)+0.1_r8       ! LAI+SAI averaged to landunit level; saved for output
            if (lai(p) > 1_r8) then
               lai(p)  = 1_r8   ! setting LAI ~ 0.1 to be a min value (since the value goes to infinity when LAI=0)
            end if              ! and 1 to be a max value as computing K involves 1 / LAI

            ! calculate Okin's shear stress ratio (which is drag partition factor) using Pierre's equation   
            K_length = 2_r8 * (1_r8/lai(p) - 1_r8)   ! Here LAI has to be non-zero to avoid blowup
            ssr(p) = (K_length+f_0*c_e)/(K_length+c_e)

            !frc_thr_rgh_fct = (rockfrc(p)*(roughfct(p))**3_r8 + (vegefrc(p)+sparfrc(p))*(ssr(p))**3_r8 )**(0.3333_r8)   ! land cover weighted mean using static GLCNMo bare land fraction LC0, dmleung 20 Dec 2021i; dmleung commented 24 Jul 2022

            ! dmleung added calculation of LUH2 bare vs veg fraction within a grid 24 Jul 2022
            !bare_frc = wt_lunit(g,istsoil) * wt_nat_patch(g,noveg) 
            !veg_frc = wt_lunit(g,istsoil) * sum(wt_nat_patch(g,(noveg+1):natpft_ub)) + wt_lunit(g,istcrop)

            !frc_thr_rgh_fct = (bare_frc*(roughfct(p))**3_r8 + veg_frc*(ssr(p))**3_r8 )**(0.3333_r8)   ! land cover weighted mean using LUH2 land cover, dmleung 24 Jul 2022

            ! dmleung added calculation of LUH2 bare vs veg fraction within a grid 6 Oct 2022
            if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
               if (patch%itype(p) == noveg) then
                  !frc_thr_rgh_fct = roughfct(p)  ! dmleung commented out 13 Dec 2022 and added next line
                  frc_thr_rgh_fct = dpfct_rock(p)
               else 
                  frc_thr_rgh_fct = ssr(p)
               end if
            else 
               frc_thr_rgh_fct = 1.0_r8
            end if

            wnd_frc_slt = fv(p) * frc_thr_rgh_fct   ! wnd_frc_slt will be used in the dust emission equation  -dmleung

            frc_thr_rghn_fct(p) = frc_thr_rgh_fct   ! save and output hybrid drag partition factor, dmleung 20 Dec 2021
         else
            wnd_frc_slt = fv(p)                     ! The value here is not important since once lnd_frc_mbl(p) <= 0.0_r8 there will be no emission.
            frc_thr_rghn_fct(p) = 0.0_r8            ! save and output hybrid drag partition factor, dmleung 20 Dec 2021
         end if

         !########## end of drag partition effect #######################################################

         !############ Add Owen effect; if not, comment out this block !-dmleung, 27 Nov 2021 ###########
         ! the following if-block comes from subr. wnd_frc_slt_get 
         ! purpose: compute the saltating friction velocity
         ! theory: saltation roughens the boundary layer, AKA "Owen's effect"

         !if (u10(p) >= wnd_rfr_thr_slt) then
         !   wnd_rfr_dlt = u10(p) - wnd_rfr_thr_slt
         !   wnd_frc_slt_dlt = 0.003_r8 * wnd_rfr_dlt * wnd_rfr_dlt
         !   wnd_frc_slt = wnd_frc_slt + wnd_frc_slt_dlt   ! careful that RHS is now wnd_frc_slt instead of fv(p)
         ! ! because wnd_frc_slt takes drag partition effect into account, but fv(p) doesn't. dmleung 27 Nov 2021
         !end if
         !########## end of Owen effect ################################################################

         ! save soil friction velocity and roughness effect before the if-statement, -dml, 1 Mar 2021, coded to CLM5 27 Nov 2021
         wnd_frc_soil(p) = wnd_frc_slt  ! save soil friction velocity for CLM output, which has drag partition and Owen effect  -dml
         ! save land mobile fraction
         lnd_frc_mble(p) = lnd_frc_mbl(p)  ! save land mobile fraction first, before the if-statement, -dml, 1 Mar 2021
         ! only perform the following calculations if lnd_frc_mbl is non-zero 

         if (lnd_frc_mbl(p) > 0.0_r8) then

            ! the following comes from subr. frc_thr_rgh_fct_get
            ! purpose: compute factor by which surface roughness increases threshold
            !          friction velocity (currently a constant)

            !frc_thr_rgh_fct = 1.0_r8

            ! the following comes from subr. frc_thr_wet_fct_get
            ! purpose: compute factor by which soil moisture increases threshold friction velocity
            ! adjust threshold velocity for inhibition by moisture
            ! modified 4/5/2002 (slevis) to use gravimetric instead of volumetric
            ! water content

            !bd = (1._r8-watsat(c,1))*2.7e3_r8      ![kg m-3] Bulk density of dry surface soil
            !gwc_sfc = h2osoi_vol(c,1)*SHR_CONST_RHOFW/bd    ![kg kg-1] Gravimetric H2O cont
            !if (gwc_sfc > gwc_thr(c)) then
            !   frc_thr_wet_fct = sqrt(1.0_r8 + 1.21_r8 * (100.0_r8*(gwc_sfc - gwc_thr(c)))**0.68_r8)
            !else
            !   frc_thr_wet_fct = 1.0_r8
            !end if

            ! slevis: adding liqfrac here, because related to effects from soil water

            !liqfrac = max( 0.0_r8, min( 1.0_r8, h2osoi_liq(c,1) / (h2osoi_ice(c,1)+h2osoi_liq(c,1)+1.0e-6_r8) ) )

            ! the following lines come from subr. dst_mbl
            ! purpose: adjust threshold friction velocity to acct for moisture and
            !          roughness. The ratio tmp1 / sqrt(forc_rho) comes from
            !          subr. wnd_frc_thr_slt_get which computes dry threshold
            !          friction velocity for saltation

            !wnd_frc_thr_slt = tmp1 / sqrt(forc_rho(c)) * frc_thr_wet_fct * frc_thr_rgh_fct

            ! reset these variables which will be updated in the following if-block

            !wnd_frc_slt = fv(p)
            flx_mss_hrz_slt_ttl = 0.0_r8
            flx_mss_vrt_dst_ttl(p) = 0.0_r8

            ! the following line comes from subr. dst_mbl
            ! purpose: threshold saltation wind speed

            wnd_rfr_thr_slt = u10(p) * wnd_frc_thr_slt / fv(p)     ! keep and use if I want Z03 scheme -dmleung

            ! the following if-block comes from subr. wnd_frc_slt_get 
            ! purpose: compute the saltating friction velocity
            ! theory: saltation roughens the boundary layer, AKA "Owen's effect"

            !if (u10(p) >= wnd_rfr_thr_slt) then
            !   wnd_rfr_dlt = u10(p) - wnd_rfr_thr_slt
            !   wnd_frc_slt_dlt = 0.003_r8 * wnd_rfr_dlt * wnd_rfr_dlt
            !   wnd_frc_slt = fv(p) + wnd_frc_slt_dlt
            !end if

            ! the following comes from subr. flx_mss_hrz_slt_ttl_Whi79_get
            ! purpose: compute vertically integrated streamwise mass flux of particles

            !if (wnd_frc_slt > wnd_frc_thr_slt) then! if want to use fluid threshold for dust emission, uncomment this one,  -dmleung 2 Dec 2021
            if (wnd_frc_slt > wnd_frc_thr_slt_it) then! if want to use impact threshold for dust emission, uncomment this one, -dmleung 2 Dec 2021

               !################### for Zender et al. (2003) scheme -dmleung ###########################
               !################ uncomment the below block if want to use Z03 scheme ###################
               !wnd_frc_rat = wnd_frc_thr_slt / wnd_frc_slt
               !flx_mss_hrz_slt_ttl = cst_slt * forc_rho(c) * (wnd_frc_slt**3.0_r8) * &
               !     (1.0_r8 - wnd_frc_rat) * (1.0_r8 + wnd_frc_rat) * (1.0_r8 + wnd_frc_rat) / grav

               ! the following loop originates from subr. dst_mbl
               ! purpose: apply land sfc and veg limitations and global tuning factor
               ! slevis: multiply flx_mss_hrz_slt_ttl by liqfrac to incude the effect 
               ! of frozen soil

               !flx_mss_hrz_slt_ttl = flx_mss_hrz_slt_ttl * lnd_frc_mbl(p) * mbl_bsn_fct(c) * &
               !     flx_mss_fdg_fct * liqfrac

               ! dmleung moved to this block
               !dst_slt_flx_rat_ttl = 100.0_r8 * exp( log(10.0_r8) * (13.4_r8 * mss_frc_cly_vld(c) - 6.0_r8) )
               !flx_mss_vrt_dst_ttl(p) = flx_mss_hrz_slt_ttl * dst_slt_flx_rat_ttl
               !########################################################################################

               !################### for Kok et al. (2014) scheme -dmleung ##############################
               !################ uncomment the below block if want to use K14 scheme ###################

               ! if want to use fluid threshold for dust emission, uncomment this one, -dmleung 27 Nov 2021
               !flx_mss_vrt_dst_ttl(p) = dst_emiss_coeff(p) * mss_frc_cly_vld(c) * forc_rho(c) * ((wnd_frc_slt**2.0_r8 - wnd_frc_thr_slt**2.0_r8) / wnd_frc_thr_slt_std) * (wnd_frc_slt / wnd_frc_thr_slt)**frag_expt  ! change forc_rho(g) to forc_rho(c) to avoid passing Nan values to the coupler -Longlei ! if want to use fluid threshold for dust emission, uncomment this one, -dml 27 Nov 2021

               ! if want to use impact threshold for dust emission, uncomment this one, -dmleung 2 Dec 2021
               flx_mss_vrt_dst_ttl(p) = dst_emiss_coeff(p) * mss_frc_cly_vld(c) * forc_rho(c) * ((wnd_frc_slt**2.0_r8 - wnd_frc_thr_slt_it**2.0_r8) / wnd_frc_thr_slt_std) * (wnd_frc_slt / wnd_frc_thr_slt_it)**frag_expt  ! if want to use impact threshold for dust emission, uncomment this one, -dml 2 Dec 2021

               ! account for bare soil fraction, frozen soil fraction, and apply global tuning parameter (Kok et al. 2014)
               flx_mss_vrt_dst_ttl(p) = flx_mss_vrt_dst_ttl(p) * lnd_frc_mbl(p) * C_tune * liqfrac
               !########################################################################################
            end if

            ! the following comes from subr. flx_mss_vrt_dst_ttl_MaB95_get
            ! purpose: diagnose total vertical mass flux of dust from vertically
            !          integrated streamwise mass flux

            !dst_slt_flx_rat_ttl = 100.0_r8 * exp( log(10.0_r8) * (13.4_r8 * mss_frc_cly_vld(c) - 6.0_r8) )  ! dmleung commented and moved to the previous block
            !flx_mss_vrt_dst_ttl(p) = flx_mss_hrz_slt_ttl * dst_slt_flx_rat_ttl

            !############## added by dmleung 2 Dec 2021 #############################################
            ! subsection for intermittency factor calculation
            ! need to use with impact threshold and cannot be used with fluid threshold
            ! Danny M. Leung, 24 Jun 2019, readded into CLM5 by dmleung 2 Dec 2021
            ! 2 Dec 2021 note: assume no buoyancy contribution to the wind fluctuation (u_sd_slt), so no obul(p) is needed. It is shown to be important for the wind fluctuations contribute little to the intermittency factor. We might add this back in the future revisions.

            ! mean lowpass-filtered wind speed at 0.1 m saltation height (assuming aerodynamic roughness length = 1e-4 m globally for ease; also assuming neutral condition)
            u_mean_slt(p) = (wnd_frc_slt/k) * log(0.1_r8 / 1e-4_r8)

            ! sd of lowpass-filtered wind speed
            !if (obul(p)==0) then
            !   zetaobu = 0
            !else 
               !zetaobu = zii(p) / obul(p)   ! For now zii is a constant of 1000 m in CLM -dml, 24 Aug 2021
            !   zetaobu = 1000_r8 / obul(p)   ! For now zii is a constant of 1000 m in CLM -dml, 24 Aug 2021
            !end if
            !stblty(p) = zetaobu    ! zetaobu get outputted as the Obukhov stability parameter
            stblty(p) = 0   ! -dmleung 2 Dec 2021: use 0 for now, assuming no buoyancy contribution. Might uncomment the above lines in future revisions.
            if ((12_r8 - 0.5_r8 * stblty(p)) .GE. 0.001_r8) then
               u_sd_slt(p) = wnd_frc_slt * (12_r8 - 0.5_r8 * stblty(p))**0.333_r8
            else
               u_sd_slt(p) = 0.001_r8   ! should have used 0 theoretically; used 0.001 here to avoid undefined values
            end if

            ! threshold velocities
            ! Here wnd_frc_thr_slt is the fluid threshold; wnd_frc_thr_dry(p) is the dry fluid threshold; B_it*wnd_frc_thr_dry(p) is the impact threshold, -dml, 1 Mar 2021
            ! fluid threshold wind at 0.1 m saltation height
            u_fld_thr(p) = (wnd_frc_thr_slt/k) * log(0.1_r8 / 1e-4_r8)
            ! impact threshold wind at 0.1 m saltation height
            u_impct_thr(p) = (wnd_frc_thr_slt_it/k) * log(0.1_r8 / 1e-4_r8)  ! to avoid model error

            ! threshold crossing rate
            thr_crs_rate(p) = (exp((u_fld_thr(p)**2_r8 - u_impct_thr(p)**2_r8 - 2_r8 * u_mean_slt(p) * (u_fld_thr(p) - u_impct_thr(p))) / (2_r8 * u_sd_slt(p)**2_r8)) + 1_r8)**(-1_r8)

            ! probability that lowpass-filtered wind speed does not exceed u_ft
            prb_crs_fld_thr(p) = 0.5_r8 * (1_r8 + erf((u_fld_thr(p) - u_mean_slt(p)) / (1.414_r8 * u_sd_slt(p))))
            ! probability that lowpass-filtered wind speed does not exceed u_it
            prb_crs_impct_thr(p) = 0.5_r8 * (1_r8 + erf((u_impct_thr(p) - u_mean_slt(p)) / (1.414_r8 * u_sd_slt(p))))

            ! intermittency factor (from 0 to 1)
            intrmtncy_fct(p) = 1_r8 - prb_crs_fld_thr(p) + thr_crs_rate(p) * (prb_crs_fld_thr(p) - prb_crs_impct_thr(p))

            ! multiply dust emission flux by intermittency factor
            if (intrmtncy_fct(p) /= intrmtncy_fct(p)) then  ! if intrmtncy_fct(p) is not NaN then multiply by intermittency factor; this statement is needed because dust emission flx_mss_vrt_dst_ttl(p) has to be non NaN (at least zero) to be outputted, dmleung 9 Jun 2021
               flx_mss_vrt_dst_ttl(p) = flx_mss_vrt_dst_ttl(p)  !  -dmleung
            else
               flx_mss_vrt_dst_ttl(p) = flx_mss_vrt_dst_ttl(p) * intrmtncy_fct(p)  ! multiply dust flux by intermittency -dmleung
            end if

            !############### end my intermittency subsection here -dmleung ########################################

         end if   ! lnd_frc_mbl > 0.0

      end do

      ! the following comes from subr. flx_mss_vrt_dst_prt in C. Zender's code
      ! purpose: partition total vertical mass flux of dust into transport bins

      do n = 1, ndst
         do m = 1, dst_src_nbr
            do fp = 1,num_nolakep
               p = filter_nolakep(fp)
               if (lnd_frc_mbl(p) > 0.0_r8) then
                  flx_mss_vrt_dst(p,n) = flx_mss_vrt_dst(p,n) +  ovr_src_snk_mss(m,n) * flx_mss_vrt_dst_ttl(p)
               end if
            end do
         end do
      end do

      do n = 1, ndst
         do fp = 1,num_nolakep
            p = filter_nolakep(fp)
            if (lnd_frc_mbl(p) > 0.0_r8) then
               flx_mss_vrt_dst_tot(p) = flx_mss_vrt_dst_tot(p) + flx_mss_vrt_dst(p,n)
            end if
         end do
      end do

    end associate

  end subroutine DustEmission

   !------------------------------------------------------------------------
  subroutine DustDryDep (bounds, &
       atm2lnd_inst, frictionvel_inst, dust_inst)
    !
    ! !DESCRIPTION: 
    !
    ! Determine Turbulent dry deposition for dust. Calculate the turbulent 
    ! component of dust dry deposition, (the turbulent deposition velocity 
    ! through the lowest atmospheric layer. CAM will calculate the settling 
    ! velocity through the whole atmospheric column. The two calculations 
    ! will determine the dust dry deposition flux to the surface.
    ! Note: Same process should occur over oceans. For the coupled CESM,
    ! we may find it more efficient to let CAM calculate the turbulent dep
    ! velocity over all surfaces. This would require passing the
    ! aerodynamic resistance, ram(1), and the friction velocity, fv, from
    ! the land to the atmosphere component. In that case, dustini need not
    ! calculate particle diamter (dmt_vwr) and particle density (dns_aer).
    ! Source: C. Zender's dry deposition code
    !
    ! !USES
    use shr_const_mod, only : SHR_CONST_PI, SHR_CONST_RDAIR, SHR_CONST_BOLTZ
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds 
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(frictionvel_type) , intent(in)    :: frictionvel_inst
    type(dust_type)        , intent(inout) :: dust_inst
    !
    ! !LOCAL VARIABLES
    integer  :: p,c,g,m,n                             ! indices
    real(r8) :: vsc_dyn_atm(bounds%begp:bounds%endp)  ! [kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm(bounds%begp:bounds%endp)  ! [m2 s-1] Kinematic viscosity of atmosphere
    real(r8) :: shm_nbr_xpn                           ! [frc] Sfc-dep exponent for aerosol-diffusion dependence on Schmidt number
    real(r8) :: shm_nbr                               ! [frc] Schmidt number
    real(r8) :: stk_nbr                               ! [frc] Stokes number
    real(r8) :: mfp_atm                               ! [m] Mean free path of air
    real(r8) :: dff_aer                               ! [m2 s-1] Brownian diffusivity of particle
    real(r8) :: rss_trb                               ! [s m-1] Resistance to turbulent deposition
    real(r8) :: slp_crc(bounds%begp:bounds%endp,ndst) ! [frc] Slip correction factor
    real(r8) :: vlc_grv(bounds%begp:bounds%endp,ndst) ! [m s-1] Settling velocity
    real(r8) :: rss_lmn(bounds%begp:bounds%endp,ndst) ! [s m-1] Quasi-laminar layer resistance
    real(r8) :: tmp                                   ! temporary 
    real(r8), parameter::shm_nbr_xpn_lnd=-2._r8/3._r8 ![frc] shm_nbr_xpn over land
    !------------------------------------------------------------------------

    associate(                                                   & 
         forc_pbot =>    atm2lnd_inst%forc_pbot_downscaled_col , & ! Input:  [real(r8)  (:)   ]  atm pressure (Pa)                                 
         forc_rho  =>    atm2lnd_inst%forc_rho_downscaled_col  , & ! Input:  [real(r8)  (:)   ]  atm density (kg/m**3)                             
         forc_t    =>    atm2lnd_inst%forc_t_downscaled_col    , & ! Input:  [real(r8)  (:)   ]  atm temperature (K)                               
         
         ram1      =>    frictionvel_inst%ram1_patch           , & ! Input:  [real(r8)  (:)   ]  aerodynamical resistance (s/m)                    
         fv        =>    frictionvel_inst%fv_patch             , & ! Input:  [real(r8)  (:)   ]  friction velocity (m/s)                           
         
         vlc_trb   =>    dust_inst%vlc_trb_patch               , & ! Output:  [real(r8) (:,:) ]  Turbulent deposn velocity (m/s)                 
         vlc_trb_1 =>    dust_inst%vlc_trb_1_patch             , & ! Output:  [real(r8) (:)   ]  Turbulent deposition velocity 1                   
         vlc_trb_2 =>    dust_inst%vlc_trb_2_patch             , & ! Output:  [real(r8) (:)   ]  Turbulent deposition velocity 2                   
         vlc_trb_3 =>    dust_inst%vlc_trb_3_patch             , & ! Output:  [real(r8) (:)   ]  Turbulent deposition velocity 3                   
         vlc_trb_4 =>    dust_inst%vlc_trb_4_patch               & ! Output:  [real(r8) (:)   ]  Turbulent deposition velocity 4                   
         )

      do p = bounds%begp,bounds%endp
         if (patch%active(p)) then
            g = patch%gridcell(p)
            c = patch%column(p)

            ! from subroutine dst_dps_dry (consider adding sanity checks from line 212)
            ! when code asks to use midlayer density, pressure, temperature,
            ! I use the data coming in from the atmosphere, ie forc_t, forc_pbot, forc_rho

            ! Quasi-laminar layer resistance: call rss_lmn_get
            ! Size-independent thermokinetic properties

            vsc_dyn_atm(p) = 1.72e-5_r8 * ((forc_t(c)/273.0_r8)**1.5_r8) * 393.0_r8 / &
                 (forc_t(c)+120.0_r8)      ![kg m-1 s-1] RoY94 p. 102
            mfp_atm = 2.0_r8 * vsc_dyn_atm(p) / &   ![m] SeP97 p. 455
                 (forc_pbot(c)*sqrt(8.0_r8/(SHR_CONST_PI*SHR_CONST_RDAIR*forc_t(c))))
            vsc_knm_atm(p) = vsc_dyn_atm(p) / forc_rho(c) ![m2 s-1] Kinematic viscosity of air

            do m = 1, ndst
               slp_crc(p,m) = 1.0_r8 + 2.0_r8 * mfp_atm * &
                    (1.257_r8+0.4_r8*exp(-1.1_r8*dmt_vwr(m)/(2.0_r8*mfp_atm))) / &
                    dmt_vwr(m)   ![frc] Slip correction factor SeP97 p. 464
               vlc_grv(p,m) = (1.0_r8/18.0_r8) * dmt_vwr(m) * dmt_vwr(m) * dns_aer * &
                    grav * slp_crc(p,m) / vsc_dyn_atm(p)   ![m s-1] Stokes' settling velocity SeP97 p. 466
               vlc_grv(p,m) = vlc_grv(p,m) * stk_crc(m)    ![m s-1] Correction to Stokes settling velocity
            end do
         end if
      end do

      do m = 1, ndst
         do p = bounds%begp,bounds%endp
            if (patch%active(p)) then
               g = patch%gridcell(p)
               c = patch%column(p)

               stk_nbr = vlc_grv(p,m) * fv(p) * fv(p) / (grav * vsc_knm_atm(p))  ![frc] SeP97 p.965
               dff_aer = SHR_CONST_BOLTZ * forc_t(c) * slp_crc(p,m) / &          ![m2 s-1]
                    (3.0_r8*SHR_CONST_PI * vsc_dyn_atm(p) * dmt_vwr(m))          !SeP97 p.474
               shm_nbr = vsc_knm_atm(p) / dff_aer                                ![frc] SeP97 p.972
               shm_nbr_xpn = shm_nbr_xpn_lnd                                     ![frc]

               ! fxm: Turning this on dramatically reduces
               ! deposition velocity in low wind regimes
               ! Schmidt number exponent is -2/3 over solid surfaces and
               ! -1/2 over liquid surfaces SlS80 p. 1014
               ! if (oro(i)==0.0) shm_nbr_xpn=shm_nbr_xpn_ocn else shm_nbr_xpn=shm_nbr_xpn_lnd
               ! [frc] Surface-dependent exponent for aerosol-diffusion dependence on Schmidt # 

               tmp = shm_nbr**shm_nbr_xpn + 10.0_r8**(-3.0_r8/stk_nbr)
               rss_lmn(p,m) = 1.0_r8 / (tmp * fv(p)) ![s m-1] SeP97 p.972,965
            end if
         end do
      end do

      ! Lowest layer: Turbulent deposition (CAM will calc. gravitational dep)

      do m = 1, ndst
         do p = bounds%begp,bounds%endp
            if (patch%active(p)) then
               rss_trb = ram1(p) + rss_lmn(p,m) + ram1(p) * rss_lmn(p,m) * vlc_grv(p,m) ![s m-1]
               vlc_trb(p,m) = 1.0_r8 / rss_trb                                          ![m s-1]
            end if
         end do
      end do

      do p = bounds%begp,bounds%endp
         if (patch%active(p)) then
            vlc_trb_1(p) = vlc_trb(p,1)
            vlc_trb_2(p) = vlc_trb(p,2)
            vlc_trb_3(p) = vlc_trb(p,3)
            vlc_trb_4(p) = vlc_trb(p,4)
         end if
      end do

    end associate

  end subroutine DustDryDep

   !------------------------------------------------------------------------
   subroutine InitDustVars(this, bounds)
     !
     ! !DESCRIPTION: 
     !
     ! Compute source efficiency factor from topography
     ! Initialize other variables used in subroutine Dust:
     ! ovr_src_snk_mss(m,n) and tmp1.
     ! Define particle diameter and density needed by atm model
     ! as well as by dry dep model
     ! Source: Paul Ginoux (for source efficiency factor)
     ! Modifications by C. Zender and later by S. Levis
     ! Rest of subroutine from C. Zender's dust model
     !
     ! !USES
     use shr_const_mod , only: SHR_CONST_PI, SHR_CONST_RDAIR
     use shr_spfn_mod  , only: erf => shr_spfn_erf
     use decompMod     , only : get_proc_bounds
     !
     ! !ARGUMENTS:
     class(dust_type)  :: this
     type(bounds_type), intent(in) :: bounds  
     !
     ! !LOCAL VARIABLES
    integer  :: fc,c,l,m,n              ! indices
    real(r8) :: ovr_src_snk_frc
    real(r8) :: sqrt2lngsdi             ! [frc] Factor in erf argument
    real(r8) :: lndmaxjovrdmdni         ! [frc] Factor in erf argument
    real(r8) :: lndminjovrdmdni         ! [frc] Factor in erf argument
    real(r8) :: ryn_nbr_frc_thr_prx_opt ! [frc] Threshold friction Reynolds number approximation for optimal size
    real(r8) :: ryn_nbr_frc_thr_opt_fnc ! [frc] Threshold friction Reynolds factor for saltation calculation
    real(r8) :: icf_fct                 ! Interpartical cohesive forces factor for saltation calc
    real(r8) :: dns_fct                 ! Density ratio factor for saltation calculation
    real(r8) :: dmt_min(ndst)           ! [m] Size grid minimum
    real(r8) :: dmt_max(ndst)           ! [m] Size grid maximum
    real(r8) :: dmt_ctr(ndst)           ! [m] Diameter at bin center
    real(r8) :: dmt_dlt(ndst)           ! [m] Width of size bin
    real(r8) :: slp_crc(ndst)           ! [frc] Slip correction factor
    real(r8) :: vlm_rsl(ndst)           ! [m3 m-3] Volume concentration resolved
    real(r8) :: vlc_stk(ndst)           ! [m s-1] Stokes settling velocity
    real(r8) :: vlc_grv(ndst)           ! [m s-1] Settling velocity
    real(r8) :: ryn_nbr_grv(ndst)       ! [frc] Reynolds number at terminal velocity
    real(r8) :: cff_drg_grv(ndst)       ! [frc] Drag coefficient at terminal velocity
    real(r8) :: tmp                     ! temporary 
    real(r8) :: ln_gsd                  ! [frc] ln(gsd)
    real(r8) :: gsd_anl                 ! [frc] Geometric standard deviation
    real(r8) :: dmt_vma                 ! [m] Mass median diameter analytic She84 p.75 Tabl.1
    real(r8) :: dmt_nma                 ! [m] Number median particle diameter
    real(r8) :: lgn_dst                 ! Lognormal distribution at sz_ctr
    real(r8) :: eps_max                 ! [frc] Relative accuracy for convergence
    real(r8) :: eps_crr                 ! [frc] Current relative accuracy
    real(r8) :: itr_idx                 ! [idx] Counting index
    real(r8) :: dns_mdp                 ! [kg m-3] Midlayer density
    real(r8) :: mfp_atm                 ! [m] Mean free path of air
    real(r8) :: vsc_dyn_atm             ! [kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm             ! [kg m-1 s-1] Kinematic viscosity of air
    real(r8) :: vlc_grv_old             ! [m s-1] Previous gravitational settling velocity
    real(r8) :: series_ratio            ! Factor for logarithmic grid
    real(r8) :: lngsdsqrttwopi_rcp      ! Factor in lognormal distribution
    real(r8) :: sz_min(sz_nbr)          ! [m] Size Bin minima
    real(r8) :: sz_max(sz_nbr)          ! [m] Size Bin maxima
    real(r8) :: sz_ctr(sz_nbr)          ! [m] Size Bin centers
    real(r8) :: sz_dlt(sz_nbr)          ! [m] Size Bin widths
    
    ! constants
    real(r8), allocatable :: dmt_vma_src(:) ! [m] Mass median diameter       BSM96 p. 73 Table 2
    real(r8), allocatable :: gsd_anl_src(:) ! [frc] Geometric std deviation  BSM96 p. 73 Table 2
    real(r8), allocatable :: mss_frc_src(:) ! [frc] Mass fraction            BSM96 p. 73 Table 2

    real(r8) :: dmt_grd(5) =                  &     ! [m] Particle diameter grid
         (/ 0.1e-6_r8, 1.0e-6_r8, 2.5e-6_r8, 5.0e-6_r8, 10.0e-6_r8 /)
    real(r8), parameter :: dmt_slt_opt = 75.0e-6_r8    ! [m] Optim diam for saltation
    real(r8), parameter :: dns_slt = 2650.0_r8         ! [kg m-3] Density of optimal saltation particles
    !------------------------------------------------------------------------

    associate(& 
         mbl_bsn_fct  =>  this%mbl_bsn_fct_col   & ! Output:  [real(r8) (:)] basin factor                                       
         )

      ! allocate module variable
      allocate (ovr_src_snk_mss(dst_src_nbr,ndst))  
      allocate (dmt_vwr(ndst))
      allocate (stk_crc(ndst))

      ! allocate local variable
      allocate (dmt_vma_src(dst_src_nbr))  
      allocate (gsd_anl_src(dst_src_nbr))  
      allocate (mss_frc_src(dst_src_nbr))  

      dmt_vma_src(:) = (/ 0.832e-6_r8 , 4.82e-6_r8 , 19.38e-6_r8 /)        
      gsd_anl_src(:) = (/ 2.10_r8     , 1.90_r8    , 1.60_r8     /)        
      mss_frc_src(:) = (/ 0.036_r8    , 0.957_r8   , 0.007_r8 /)                  

      ! the following comes from (1) szdstlgn.F subroutine ovr_src_snk_frc_get
      !                      and (2) dstszdst.F subroutine dst_szdst_ini
      ! purpose(1): given one set (the "source") of lognormal distributions,
      !             and one set of bin boundaries (the "sink"), compute and return
      !             the overlap factors between the source and sink distributions
      ! purpose(2): set important statistics of size distributions

      do m = 1, dst_src_nbr
         sqrt2lngsdi = sqrt(2.0_r8) * log(gsd_anl_src(m))
         do n = 1, ndst
            lndmaxjovrdmdni = log(dmt_grd(n+1)/dmt_vma_src(m))
            lndminjovrdmdni = log(dmt_grd(n  )/dmt_vma_src(m))
            ovr_src_snk_frc = 0.5_r8 * (erf(lndmaxjovrdmdni/sqrt2lngsdi) - &
                 erf(lndminjovrdmdni/sqrt2lngsdi))
            ovr_src_snk_mss(m,n) = ovr_src_snk_frc * mss_frc_src(m)
         end do
      end do

      ! The following code from subroutine wnd_frc_thr_slt_get was placed 
      ! here because tmp1 needs to be defined just once

      ryn_nbr_frc_thr_prx_opt = 0.38_r8 + 1331.0_r8 * (100.0_r8*dmt_slt_opt)**1.56_r8

      if (ryn_nbr_frc_thr_prx_opt < 0.03_r8) then
         write(iulog,*) 'dstmbl: ryn_nbr_frc_thr_prx_opt < 0.03'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      else if (ryn_nbr_frc_thr_prx_opt < 10.0_r8) then
         ryn_nbr_frc_thr_opt_fnc = -1.0_r8 + 1.928_r8 * (ryn_nbr_frc_thr_prx_opt**0.0922_r8)
         ryn_nbr_frc_thr_opt_fnc = 0.1291_r8 * 0.1291_r8 / ryn_nbr_frc_thr_opt_fnc
      else
         ryn_nbr_frc_thr_opt_fnc = 1.0_r8 - 0.0858_r8 * exp(-0.0617_r8*(ryn_nbr_frc_thr_prx_opt-10.0_r8))
         ryn_nbr_frc_thr_opt_fnc = 0.120_r8 * 0.120_r8 * ryn_nbr_frc_thr_opt_fnc * ryn_nbr_frc_thr_opt_fnc
      end if

      icf_fct = 1.0_r8 + 6.0e-07_r8 / (dns_slt * grav * (dmt_slt_opt**2.5_r8))
      dns_fct = dns_slt * grav * dmt_slt_opt
      tmp1 = sqrt(icf_fct * dns_fct * ryn_nbr_frc_thr_opt_fnc)

      ! Introducing particle diameter. Needed by atm model and by dry dep model.
      ! Taken from Charlie Zender's subroutines dst_psd_ini, dst_sz_rsl,
      ! grd_mk (dstpsd.F90) and subroutine lgn_evl (psdlgn.F90)

      ! Charlie allows logarithmic or linear option for size distribution
      ! however, he hardwires the distribution to logarithmic in his code
      ! therefore, I take his logarithmic code only
      ! furthermore, if dst_nbr == 4, he overrides the automatic grid calculation
      ! he currently works with dst_nbr = 4, so I only take the relevant code
      ! if ndst ever becomes different from 4, must add call grd_mk (dstpsd.F90)
      ! as done in subroutine dst_psd_ini
      ! note that here ndst = dst_nbr

      ! Override automatic grid with preset grid if available

      if (ndst == 4) then
         do n = 1, ndst
            dmt_min(n) = dmt_grd(n)                       ![m] Max diameter in bin
            dmt_max(n) = dmt_grd(n+1)                     ![m] Min diameter in bin
            dmt_ctr(n) = 0.5_r8 * (dmt_min(n)+dmt_max(n)) ![m] Diameter at bin ctr
            dmt_dlt(n) = dmt_max(n)-dmt_min(n)            ![m] Width of size bin
         end do
      else
         write(iulog,*) 'Dustini error: ndst must equal to 4 with current code'
         call endrun(msg=errMsg(sourcefile, __LINE__))
         !see more comments above end if ndst == 4
      end if

      ! Bin physical properties

      gsd_anl = 2.0_r8      ! [frc] Geometric std dev PaG77 p. 2080 Table1
      ln_gsd = log(gsd_anl)
      dns_aer = 2.5e+3_r8   ! [kg m-3] Aerosol density

      ! Set a fundamental statistic for each bin

      dmt_vma = 3.5000e-6_r8 ! [m] Mass median diameter analytic She84 p.75 Table1

      ! Compute analytic size statistics
      ! Convert mass median diameter to number median diameter (call vma2nma)

      dmt_nma = dmt_vma * exp(-3.0_r8*ln_gsd*ln_gsd) ! [m]

      ! Compute resolved size statistics for each size distribution
      ! In C. Zender's code call dst_sz_rsl

      do n = 1, ndst

         series_ratio = (dmt_max(n)/dmt_min(n))**(1.0_r8/sz_nbr)
         sz_min(1) = dmt_min(n)
         do m = 2, sz_nbr                            ! Loop starts at 2
            sz_min(m) = sz_min(m-1) * series_ratio
         end do

         ! Derived grid values
         do m = 1, sz_nbr-1                          ! Loop ends at sz_nbr-1
            sz_max(m) = sz_min(m+1)                  ! [m]
         end do
         sz_max(sz_nbr) = dmt_max(n)                 ! [m]

         ! Final derived grid values
         do m = 1, sz_nbr
            sz_ctr(m) = 0.5_r8 * (sz_min(m)+sz_max(m))
            sz_dlt(m) = sz_max(m)-sz_min(m)
         end do

         lngsdsqrttwopi_rcp = 1.0_r8 / (ln_gsd*sqrt(2.0_r8*SHR_CONST_PI))
         dmt_vwr(n) = 0.0_r8 ! [m] Mass wgted diameter resolved
         vlm_rsl(n) = 0.0_r8 ! [m3 m-3] Volume concentration resolved

         do m = 1, sz_nbr

            ! Evaluate lognormal distribution for these sizes (call lgn_evl)
            tmp = log(sz_ctr(m)/dmt_nma) / ln_gsd
            lgn_dst = lngsdsqrttwopi_rcp * exp(-0.5_r8*tmp*tmp) / sz_ctr(m)

            ! Integrate moments of size distribution
            dmt_vwr(n) = dmt_vwr(n) + sz_ctr(m) *                    &
                 SHR_CONST_PI / 6.0_r8 * (sz_ctr(m)**3.0_r8) * & ![m3] Volume
                 lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn
            vlm_rsl(n) = vlm_rsl(n) +                                &
                 SHR_CONST_PI / 6.0_r8 * (sz_ctr(m)**3.0_r8) * & ![m3] Volume
                 lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn

         end do

         dmt_vwr(n) = dmt_vwr(n) / vlm_rsl(n) ![m] Mass weighted diameter resolved

      end do

      ! calculate correction to Stokes' settling velocity (subroutine stk_crc_get)

      eps_max = 1.0e-4_r8
      dns_mdp = 100000._r8 / (295.0_r8*SHR_CONST_RDAIR) ![kg m-3] const prs_mdp & tpt_vrt

      ! Size-independent thermokinetic properties

      vsc_dyn_atm = 1.72e-5_r8 * ((295.0_r8/273.0_r8)**1.5_r8) * 393.0_r8 / &
           (295.0_r8+120.0_r8)      ![kg m-1 s-1] RoY94 p.102 tpt_mdp=295.0
      mfp_atm = 2.0_r8 * vsc_dyn_atm / &  !SeP97 p. 455 constant prs_mdp, tpt_mdp
           (100000._r8*sqrt(8.0_r8/(SHR_CONST_PI*SHR_CONST_RDAIR*295.0_r8)))
      vsc_knm_atm = vsc_dyn_atm / dns_mdp ![m2 s-1] Kinematic viscosity of air

      do m = 1, ndst
         slp_crc(m) = 1.0_r8 + 2.0_r8 * mfp_atm *                      &
              (1.257_r8+0.4_r8*exp(-1.1_r8*dmt_vwr(m)/(2.0_r8*mfp_atm))) / &
              dmt_vwr(m)                      ! [frc] Slip correction factor SeP97 p.464
         vlc_stk(m) = (1.0_r8/18.0_r8) * dmt_vwr(m) * dmt_vwr(m) * dns_aer * &
              grav * slp_crc(m) / vsc_dyn_atm ! [m s-1] SeP97 p.466
      end do

      ! For Reynolds number flows Re < 0.1 Stokes' velocity is valid for
      ! vlc_grv SeP97 p. 466 (8.42). For larger Re, inertial effects become
      ! important and empirical drag coefficients must be employed
      ! Implicit equation for Re, Cd, and Vt is SeP97 p. 467 (8.44)
      ! Using Stokes' velocity rather than iterative solution with empirical
      ! drag coefficient causes 60% errors for D = 200 um SeP97 p. 468

      ! Iterative solution for drag coefficient, Reynolds number, and terminal veloc
      do m = 1, ndst

         ! Initialize accuracy and counter
         eps_crr = eps_max + 1.0_r8  ![frc] Current relative accuracy
         itr_idx = 0                 ![idx] Counting index

         ! Initial guess for vlc_grv is exact for Re < 0.1
         vlc_grv(m) = vlc_stk(m)     ![m s-1]

         do while(eps_crr > eps_max)

            ! Save terminal velocity for convergence test
            vlc_grv_old = vlc_grv(m) ![m s-1]
            ryn_nbr_grv(m) = vlc_grv(m) * dmt_vwr(m) / vsc_knm_atm !SeP97 p.460

            ! Update drag coefficient based on new Reynolds number
            if (ryn_nbr_grv(m) < 0.1_r8) then
               cff_drg_grv(m) = 24.0_r8 / ryn_nbr_grv(m) !Stokes' law Sep97 p.463 (8.32)
            else if (ryn_nbr_grv(m) < 2.0_r8) then
               cff_drg_grv(m) = (24.0_r8/ryn_nbr_grv(m)) *    &
                    (1.0_r8 + 3.0_r8*ryn_nbr_grv(m)/16.0_r8 + &
                    9.0_r8*ryn_nbr_grv(m)*ryn_nbr_grv(m)*     &
                    log(2.0_r8*ryn_nbr_grv(m))/160.0_r8)        !Sep97 p.463 (8.32)
            else if (ryn_nbr_grv(m) < 500.0_r8) then
               cff_drg_grv(m) = (24.0_r8/ryn_nbr_grv(m)) * &
                    (1.0_r8 + 0.15_r8*ryn_nbr_grv(m)**0.687_r8) !Sep97 p.463 (8.32)
            else if (ryn_nbr_grv(m) < 2.0e5_r8) then
               cff_drg_grv(m) = 0.44_r8                         !Sep97 p.463 (8.32)
            else
               write(iulog,'(a,es9.2)') "ryn_nbr_grv(m) = ",ryn_nbr_grv(m)
               write(iulog,*)'Dustini error: Reynolds number too large in stk_crc_get()'
               call endrun(msg=errMsg(sourcefile, __LINE__))
            end if

            ! Update terminal velocity based on new Reynolds number and drag coeff
            ! [m s-1] Terminal veloc SeP97 p.467 (8.44)

            vlc_grv(m) = sqrt(4.0_r8 * grav * dmt_vwr(m) * slp_crc(m) * dns_aer / &
                 (3.0_r8*cff_drg_grv(m)*dns_mdp))   
            eps_crr = abs((vlc_grv(m)-vlc_grv_old)/vlc_grv(m)) !Relative convergence
            if (itr_idx == 12) then
               ! Numerical pingpong may occur when Re = 0.1, 2.0, or 500.0
               ! due to discontinuities in derivative of drag coefficient
               vlc_grv(m) = 0.5_r8 * (vlc_grv(m)+vlc_grv_old)  ! [m s-1]
            end if
            if (itr_idx > 20) then
               write(iulog,*) 'Dustini error: Terminal velocity not converging ',&
                    ' in stk_crc_get(), breaking loop...'
               goto 100                                        !to next iteration
            end if
            itr_idx = itr_idx + 1

         end do                                                !end while

100      continue   !Label to jump to when iteration does not converge
      end do   !end loop over size

      ! Compute factors to convert Stokes' settling velocities to
      ! actual settling velocities

      do m = 1, ndst
         stk_crc(m) = vlc_grv(m) / vlc_stk(m)
      end do

    end associate 

  end subroutine InitDustVars

end module DUSTMod
