module DustEmisLeung2023

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
  use clm_varpar           , only : dst_src_nbr, ndst
  use clm_varcon           , only : grav, spval
  use landunit_varcon      , only : istcrop, istsoil
  use clm_varctl           , only : iulog
  use abortutils           , only : endrun
  use decompMod            , only : bounds_type, subgrid_level_landunit
  use atm2lndType          , only : atm2lnd_type
  use SoilStateType        , only : soilstate_type
  use CanopyStateType      , only : canopystate_type
  use WaterStateBulkType       , only : waterstatebulk_type
  use WaterDiagnosticBulkType       , only : waterdiagnosticbulk_type
  use FrictionVelocityMod  , only : frictionvel_type
  use LandunitType         , only : lun
  use PatchType            , only : patch
  use DustEmisBase         , only : dust_emis_base_type
  use pftconMod            , only : noveg
  use PrigentRoughnessStreamType, only : prigent_roughness_stream_type
  !
  ! !PUBLIC TYPES
  implicit none
  private
  !
  ! !PRIVATE DATA:
  !
  !
  ! !PUBLIC DATA TYPES:
  !
  type, public, extends(dust_emis_base_type) :: dust_emis_leung2023_type

      real(r8), pointer, private :: dst_emiss_coeff_patch     (:)   ! dust emission coefficient (unitless)
      real(r8), pointer, private :: wnd_frc_thr_patch         (:)   ! wet fluid threshold (m/s)
      real(r8), pointer, private :: wnd_frc_thr_dry_patch     (:)   ! dry fluid threshold (m/s)
      real(r8), pointer, private :: wnd_frc_thr_it_patch     (:)   ! impact threshold (m/s)
      real(r8), pointer, private :: lnd_frc_mble_patch        (:)   ! land mobile fraction
      real(r8), pointer, private :: liq_frac_patch            (:)   ! liquid fraction of total water
      real(r8), pointer, private :: wnd_frc_soil_patch        (:)   ! soil wind friction velocity (m/s)
      real(r8), pointer, private :: gwc_patch                 (:)   ! gravimetric water content (kg/kg)
      real(r8), pointer, private :: intrmtncy_fct_patch       (:)   ! intermittency factor, accounting for turbulence shutting down dust emissions (unitless)
      real(r8), pointer, private :: stblty_patch              (:)   ! stability parameter for checking stability condition (stblty < 0 is unstable atmosphere)
      real(r8), pointer, private :: u_mean_slt_patch          (:)   ! wind speed 0.1 m level of dust saltation (m/s)
      real(r8), pointer, private :: u_sd_slt_patch            (:)   ! sd of wind speed 0.1 m level of dust saltation (m/s)
      real(r8), pointer, private :: u_fld_thr_patch           (:)   ! fluid threshold wind speed 0.1 m level of dust saltation (m/s)
      real(r8), pointer, private :: u_impct_thr_patch         (:)   ! impact threshold wind speed at 0.1 m level of dust saltation (m/s)
      real(r8), pointer, private :: thr_crs_rate_patch        (:)   ! threshold crossing rate (unitless)
      real(r8), pointer, private :: prb_crs_fld_thr_patch     (:)   ! probability of wind speed crossing fluid threshold
      real(r8), pointer, private :: prb_crs_impct_thr_patch   (:)   ! probability of wind speed crossing impact threshold
      real(r8), pointer, private :: ssr_patch                 (:)   ! [dimless] integrated shear stress ratiio, defined by Okin (2008) and then integrated by Caroline Pierre et al. (2014)
      real(r8), pointer, private :: vai_Okin_patch             (:)   ! [m2 leaf /m2 land] LAI+SAI for calculating Okin drag partition
      real(r8), pointer, private :: frc_thr_rghn_fct_patch    (:)   ! [dimless] hybrid drag partition (or called roughness) factor
      real(r8), pointer, private :: wnd_frc_thr_std_patch     (:)   ! standardized fluid threshold friction velocity (m/s)
      type(prigent_roughness_stream_type), private :: prigent_roughness_stream      ! Prigent roughness stream data
      real(r8), pointer, private :: dpfct_rock_patch          (:)   ! [fraction] rock drag partition factor, time-constant

   contains

     procedure , public  :: Init => InitLeung2023
     procedure , public  :: DustEmission    ! Dust mobilization
     procedure , public  :: Clean => CleanLeung2023
     procedure , private :: InitAllocate    ! Allocate data
     procedure , private :: InitHistory     ! History initialization
     procedure , private :: InitCold
     procedure , private :: CalcDragPartition ! Calculate drag partitioning based on Prigent roughness stream

  end type dust_emis_leung2023_type

  interface dust_emis_leung2023_type
     ! initialize a new dust emission object
      module procedure constructor
  end interface dust_emis_leung2023_type
  !------------------------------------------------------------------------

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  type(dust_emis_leung2023_type) function constructor()
  !
  ! Creates a dust emission object for Leung-2023 type
  ! For now this is just a placeholder
  !-----------------------------------------------------------------------

  end function constructor

  !------------------------------------------------------------------------

  subroutine InitLeung2023(this, bounds, NLFilename)

   ! Initialization for this extended class, calling base level initiation and adding to it
    class(dust_emis_leung2023_type) :: this
    type(bounds_type), intent(in) :: bounds
    character(len=*),  intent(in) :: NLFilename

    call this%InitBase(bounds, NLFilename)
    call this%prigent_roughness_stream%Init( bounds, NLFilename )
    call this%InitAllocate (bounds)
    call this%InitHistory  (bounds)
    call this%InitCold     (bounds)

  end subroutine InitLeung2023

  !------------------------------------------------------------------------

  subroutine InitAllocate(this, bounds)
   !
   ! !ARGUMENTS:
   class (dust_emis_leung2023_type) :: this
   type(bounds_type), intent(in) :: bounds
   !
   ! !LOCAL VARIABLES:
   integer :: begc,endc
   integer :: begp,endp
   !------------------------------------------------------------------------

   begc = bounds%begc ; endc = bounds%endc
   begp = bounds%begp ; endp = bounds%endp
   allocate(this%dst_emiss_coeff_patch     (begp:endp))        ; this%dst_emiss_coeff_patch     (:)   = nan
   allocate(this%wnd_frc_thr_patch         (begp:endp))        ; this%wnd_frc_thr_patch         (:)   = nan
   allocate(this%wnd_frc_thr_dry_patch     (begp:endp))        ; this%wnd_frc_thr_dry_patch     (:)   = nan
   allocate(this%wnd_frc_thr_it_patch      (begp:endp))        ; this%wnd_frc_thr_it_patch      (:)   = nan
   allocate(this%lnd_frc_mble_patch        (begp:endp))        ; this%lnd_frc_mble_patch        (:)   = nan
   allocate(this%wnd_frc_soil_patch        (begp:endp))        ; this%wnd_frc_soil_patch        (:)   = nan
   allocate(this%gwc_patch                 (begp:endp))        ; this%gwc_patch                 (:)   = nan
   allocate(this%liq_frac_patch            (begp:endp))        ; this%liq_frac_patch            (:)   = nan
   allocate(this%intrmtncy_fct_patch       (begp:endp))        ; this%intrmtncy_fct_patch       (:)   = nan
   allocate(this%stblty_patch              (begp:endp))        ; this%stblty_patch              (:)   = nan
   allocate(this%u_mean_slt_patch          (begp:endp))        ; this%u_mean_slt_patch          (:)   = nan
   allocate(this%u_sd_slt_patch            (begp:endp))        ; this%u_sd_slt_patch            (:)   = nan
   allocate(this%u_fld_thr_patch           (begp:endp))        ; this%u_fld_thr_patch           (:)   = nan
   allocate(this%u_impct_thr_patch         (begp:endp))        ; this%u_impct_thr_patch         (:)   = nan
   allocate(this%thr_crs_rate_patch        (begp:endp))        ; this%thr_crs_rate_patch        (:)   = nan
   allocate(this%prb_crs_fld_thr_patch     (begp:endp))        ; this%prb_crs_fld_thr_patch     (:)   = nan
   allocate(this%prb_crs_impct_thr_patch   (begp:endp))        ; this%prb_crs_impct_thr_patch   (:)   = nan
   allocate(this%ssr_patch                 (begp:endp))        ; this%ssr_patch                 (:)   = nan
   allocate(this%vai_Okin_patch            (begp:endp))        ; this%vai_Okin_patch            (:)   = nan
   allocate(this%frc_thr_rghn_fct_patch    (begp:endp))        ; this%frc_thr_rghn_fct_patch    (:)   = nan
   allocate(this%wnd_frc_thr_std_patch     (begp:endp))        ; this%wnd_frc_thr_std_patch     (:)   = nan
   allocate(this%dpfct_rock_patch          (begp:endp))        ; this%dpfct_rock_patch          (:)   = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------

  subroutine CleanLeung2023(this)
    !
    ! Deallocation for this extended class, calling base level deallocation and adding to it
    ! !ARGUMENTS:
    class (dust_emis_leung2023_type) :: this
    !
    ! !LOCAL VARIABLES:
    !------------------------------------------------------------------------

    call this%CleanBase()
    call this%prigent_roughness_stream%Clean( )
    deallocate(this%dst_emiss_coeff_patch     )
    deallocate(this%wnd_frc_thr_patch         )
    deallocate(this%wnd_frc_thr_dry_patch     )
    deallocate(this%wnd_frc_thr_it_patch      )
    deallocate(this%lnd_frc_mble_patch        )
    deallocate(this%wnd_frc_soil_patch        )
    deallocate(this%gwc_patch                 )
    deallocate(this%liq_frac_patch            )
    deallocate(this%intrmtncy_fct_patch       )
    deallocate(this%stblty_patch              )
    deallocate(this%u_mean_slt_patch          )
    deallocate(this%u_sd_slt_patch            )
    deallocate(this%u_fld_thr_patch           )
    deallocate(this%u_impct_thr_patch         )
    deallocate(this%thr_crs_rate_patch        )
    deallocate(this%prb_crs_fld_thr_patch     )
    deallocate(this%prb_crs_impct_thr_patch   )
    deallocate(this%ssr_patch                 )
    deallocate(this%vai_Okin_patch            )
    deallocate(this%frc_thr_rghn_fct_patch    )
    deallocate(this%wnd_frc_thr_std_patch     )
    deallocate(this%dpfct_rock_patch          )

  end subroutine CleanLeung2023

  !------------------------------------------------------------------------

  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d
    !
    !
    ! !ARGUMENTS:
    class (dust_emis_leung2023_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begc,endc
    integer :: begp,endp
    !------------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc
    begp = bounds%begp; endp = bounds%endp
    this%dst_emiss_coeff_patch(begp:endp) = spval
    call hist_addfld1d (fname='DUST_EMIS_COEFF', units='dimensionless',  &
         avgflag='A', long_name='soil erodibility or dust emission coefficient for Kok emission scheme', &
         ptr_patch=this%dst_emiss_coeff_patch, set_lake=0._r8, set_urb=0._r8)
    this%wnd_frc_thr_patch(begp:endp) = spval
    call hist_addfld1d (fname='WND_FRC_FT', units='m/s',  &
         avgflag='A', long_name='fluid threshold friction velocity', &
         ptr_patch=this%wnd_frc_thr_patch, set_lake=0._r8, set_urb=0._r8)
    this%wnd_frc_thr_dry_patch(begp:endp) = spval
    call hist_addfld1d (fname='WND_FRC_FT_DRY', units='m/s',  &
         avgflag='A', long_name='dry fluid threshold friction velocity', &
         ptr_patch=this%wnd_frc_thr_dry_patch, set_lake=0._r8, set_urb=0._r8)
    this%wnd_frc_thr_it_patch(begp:endp) = spval
    call hist_addfld1d (fname='WND_FRC_IT', units='m/s',  &
         avgflag='A', long_name='impact threshold friction velocity', &
         ptr_patch=this%wnd_frc_thr_it_patch, set_lake=0._r8, set_urb=0._r8)
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
         avgflag='A', long_name='gravimetric soil moisture at the topmost soil layer', &
         ptr_patch=this%gwc_patch, set_lake=0._r8, set_urb=0._r8)
    this%liq_frac_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIQ_FRAC', units='dimensionless',  &
         avgflag='A', long_name='fraction of total water that is liquid', &
         ptr_patch=this%liq_frac_patch, set_lake=0._r8, set_urb=0._r8)
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
    this%ssr_patch(begp:endp) = spval
    call hist_addfld1d (fname='SSR', units='m/s',  &
         avgflag='A', long_name='Okin-Pierre vegetation shear stress ratio (drag partition factor)', &
         ptr_patch=this%ssr_patch, set_lake=0._r8, set_urb=0._r8)
    this%vai_Okin_patch(begp:endp) = spval
    call hist_addfld1d (fname='VAI_OKIN', units='m/s',  &
         avgflag='A', long_name='vegetation area index used in the Okin-Pierre plant drag partition scheme', &
         ptr_patch=this%vai_Okin_patch, set_lake=0._r8, set_urb=0._r8)
    this%frc_thr_rghn_fct_patch(begp:endp) = spval
    call hist_addfld1d (fname='FRC_THR_RGHN_FCT', units='dimensionless',  &
         avgflag='A', long_name='hybrid drag partition (or roughness) factor', &
         ptr_patch=this%frc_thr_rghn_fct_patch, set_lake=0._r8, set_urb=0._r8)
    this%wnd_frc_thr_std_patch(begp:endp) = spval
    call hist_addfld1d (fname='WND_FRC_FT_STD', units='m/s',  &
         avgflag='A', long_name='standardized fluid threshold friction velocity', &
         ptr_patch=this%wnd_frc_thr_std_patch, set_lake=0._r8, set_urb=0._r8)
    this%dpfct_rock_patch(begp:endp) = spval
    call hist_addfld1d (fname='DPFCT_ROCK', units='m/s',  &
         avgflag='A', long_name='rock drag partition factor', &
         ptr_patch=this%dpfct_rock_patch)

  end subroutine InitHistory

  !-----------------------------------------------------------------------

  subroutine InitCold(this, bounds)
    !
    ! Initialize values from a cold start
    ! !ARGUMENTS:
    class (dust_emis_leung2023_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------
          ! Caulculate Drag Partition factor (Marticorena and Bergametti 1995 formulation)
    if ( this%prigent_roughness_stream%useStreams() )then !if usestreams == true, and it should be always true
      call this%CalcDragPartition( bounds )
    else

      call endrun( "ERROR:: Drag partitioning MUST now use a streams file of aeolian roughness length to calculate, it can no longer read from the fsurdat file" )
    end if

  end subroutine InitCold

  !------------------------------------------------------------------------

  subroutine DustEmission (this, bounds, &
       num_nolakep, filter_nolakep, &
       atm2lnd_inst, soilstate_inst, canopystate_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, &
       frictionvel_inst)
    !
    ! !DESCRIPTION:
    ! Dust mobilization. This code simulates dust mobilization due to wind
    ! from the surface into the lowest atmospheric layer
    ! On output flx_mss_vrt_dst(ndst) is the surface dust emission
    ! (kg/m**2/s) [ + = to atm]
    !
    ! !USES
    use shr_const_mod, only : SHR_CONST_RHOFW
    use subgridaveMod, only : p2g
    !
    ! !ARGUMENTS:
    class (dust_emis_leung2023_type)      :: this
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_nolakep                 ! number of column non-lake points in patch filter
    integer                , intent(in)    :: filter_nolakep(num_nolakep) ! patch filter for non-lake points
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(soilstate_type)   , intent(in)    :: soilstate_inst
    type(canopystate_type) , intent(in)    :: canopystate_inst
    type(waterstatebulk_type)  , intent(in)    :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)  , intent(in)    :: waterdiagnosticbulk_inst
    type(frictionvel_type) , intent(in)    :: frictionvel_inst

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
    real(r8) :: wnd_frc_thr_slt
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
    !
    ! constants
    !
    real(r8), parameter :: cst_slt = 2.61_r8           ! [frc] Saltation constant
    real(r8), parameter :: flx_mss_fdg_fct = 5.0e-4_r8 ! [frc] Empir. mass flx tuning eflx_lh_vegt
    real(r8), parameter :: vai_mbl_thr = 0.3_r8        ! [m2 m-2] VAI threshold quenching dust mobilization
    character(len=*),parameter :: subname = 'DUSTEmission'
    !------------------------------------------------------------------------

    write(iulog,*)
    write(iulog,*)
    write(iulog,*) subname//'::WARNING: CURRENTLY THIS IS JUST THE ZENDER 2003 VERIONS OF DUST EMISSIONS!'
    write(iulog,*)
    write(iulog,*)
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

         flx_mss_vrt_dst     => this%flx_mss_vrt_dst_patch           , & ! Output: [real(r8) (:,:) ]  surface dust emission (kg/m**2/s)
         flx_mss_vrt_dst_tot => this%flx_mss_vrt_dst_tot_patch         & ! Output: [real(r8) (:)   ]  total dust flux back to atmosphere (pft)
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
            write(iulog,*)'Error dstmbl: pft= ',p,' lnd_frc_mbl(p)= ',lnd_frc_mbl(p), &
                           errMsg(sourcefile, __LINE__)
            call endrun("Bad value for dust mobilization fraction")
            return
         end if
      end do

      ! reset history output variables before next if-statement to avoid output = inf

      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         flx_mss_vrt_dst_tot(p) = 0.0_r8
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

         ! only perform the following calculations if lnd_frc_mbl is non-zero

         if (lnd_frc_mbl(p) > 0.0_r8) then

            ! the following comes from subr. frc_thr_rgh_fct_get
            ! purpose: compute factor by which surface roughness increases threshold
            !          friction velocity (currently a constant)

            frc_thr_rgh_fct = 1.0_r8

            ! the following comes from subr. frc_thr_wet_fct_get
            ! purpose: compute factor by which soil moisture increases threshold friction velocity
            ! adjust threshold velocity for inhibition by moisture
            ! modified 4/5/2002 (slevis) to use gravimetric instead of volumetric
            ! water content

            bd = (1._r8-watsat(c,1))*2.7e3_r8      ![kg m-3] Bulk density of dry surface soil
            gwc_sfc = h2osoi_vol(c,1)*SHR_CONST_RHOFW/bd    ![kg kg-1] Gravimetric H2O cont
            if (gwc_sfc > gwc_thr(c)) then
               frc_thr_wet_fct = sqrt(1.0_r8 + 1.21_r8 * (100.0_r8*(gwc_sfc - gwc_thr(c)))**0.68_r8)
            else
               frc_thr_wet_fct = 1.0_r8
            end if

            ! slevis: adding liqfrac here, because related to effects from soil water

            liqfrac = max( 0.0_r8, min( 1.0_r8, h2osoi_liq(c,1) / (h2osoi_ice(c,1)+h2osoi_liq(c,1)+1.0e-6_r8) ) )

            ! the following lines come from subr. dst_mbl
            ! purpose: adjust threshold friction velocity to acct for moisture and
            !          roughness. The ratio saltation_factor / sqrt(forc_rho) comes from
            !          subr. wnd_frc_thr_slt_get which computes dry threshold
            !          friction velocity for saltation

            wnd_frc_thr_slt = this%saltation_factor / sqrt(forc_rho(c)) * frc_thr_wet_fct * frc_thr_rgh_fct

            ! reset these variables which will be updated in the following if-block

            wnd_frc_slt = fv(p)
            flx_mss_hrz_slt_ttl = 0.0_r8
            flx_mss_vrt_dst_ttl(p) = 0.0_r8

            ! the following line comes from subr. dst_mbl
            ! purpose: threshold saltation wind speed

            wnd_rfr_thr_slt = u10(p) * wnd_frc_thr_slt / fv(p)

            ! the following if-block comes from subr. wnd_frc_slt_get
            ! purpose: compute the saltating friction velocity
            ! theory: saltation roughens the boundary layer, AKA "Owen's effect"

            if (u10(p) >= wnd_rfr_thr_slt) then
               wnd_rfr_dlt = u10(p) - wnd_rfr_thr_slt
               wnd_frc_slt_dlt = 0.003_r8 * wnd_rfr_dlt * wnd_rfr_dlt
               wnd_frc_slt = fv(p) + wnd_frc_slt_dlt
            end if

            ! the following comes from subr. flx_mss_hrz_slt_ttl_Whi79_get
            ! purpose: compute vertically integrated streamwise mass flux of particles

            if (wnd_frc_slt > wnd_frc_thr_slt) then
               wnd_frc_rat = wnd_frc_thr_slt / wnd_frc_slt
               flx_mss_hrz_slt_ttl = cst_slt * forc_rho(c) * (wnd_frc_slt**3.0_r8) * &
                    (1.0_r8 - wnd_frc_rat) * (1.0_r8 + wnd_frc_rat) * (1.0_r8 + wnd_frc_rat) / grav

               ! the following loop originates from subr. dst_mbl
               ! purpose: apply land sfc and veg limitations and global tuning factor
               ! slevis: multiply flx_mss_hrz_slt_ttl by liqfrac to incude the effect
               ! of frozen soil

               flx_mss_hrz_slt_ttl = flx_mss_hrz_slt_ttl * lnd_frc_mbl(p) * &
                    flx_mss_fdg_fct * liqfrac
            end if

            ! the following comes from subr. flx_mss_vrt_dst_ttl_MaB95_get
            ! purpose: diagnose total vertical mass flux of dust from vertically
            !          integrated streamwise mass flux

            dst_slt_flx_rat_ttl = 100.0_r8 * exp( log(10.0_r8) * (13.4_r8 * mss_frc_cly_vld(c) - 6.0_r8) )
            flx_mss_vrt_dst_ttl(p) = flx_mss_hrz_slt_ttl * dst_slt_flx_rat_ttl

         end if   ! lnd_frc_mbl > 0.0

      end do

      ! the following comes from subr. flx_mss_vrt_dst_prt in C. Zender's code
      ! purpose: partition total vertical mass flux of dust into transport bins

      do n = 1, ndst
         do m = 1, dst_src_nbr
            do fp = 1,num_nolakep
               p = filter_nolakep(fp)
               if (lnd_frc_mbl(p) > 0.0_r8) then
                  flx_mss_vrt_dst(p,n) = flx_mss_vrt_dst(p,n) +  this%ovr_src_snk_mss(m,n) * flx_mss_vrt_dst_ttl(p)
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
  subroutine CalcDragPartition(this, bounds)
   !
   ! !DESCRIPTION:
   ! Commented below by Danny M. Leung 31 Dec 2022
   ! Calculate the drag partition effect of friction velocity due to surface roughness following
   ! Leung et al. (2023).  This module is used in the dust emission module DUSTMod.F90 for
   ! calculating drag partitioning. The drag partition equation comes from Marticorena and
   ! Bergametti (1995) with constants modified by Darmenova et al. (2009). Here it is assumed
   ! that this equation is used only over arid/desertic regions, such that Catherine Prigent's
   ! roughness measurements represents mostly rocks. For more vegetated areas, the vegetation
   ! roughness and drag partitioning are calculated in the DustEmission subroutine. This
   ! subroutine is used in the InitCold subroutine of DUSTMod.F90.
   !
   ! !USES:
   use PatchType               , only : patch
   use landunit_varcon         , only : istdlak
   use LandunitType            , only : lun
   !
   ! !ARGUMENTS:
   implicit none
   class (dust_emis_leung2023_type) :: this
   type(bounds_type), intent(in) :: bounds
   !
   ! !LOCAL VARIABLES:
   integer  :: g, p, fp, l    ! Indices
   real(r8) :: z0s         ! smooth roughness length (m)

   ! constants
   real(r8), parameter :: D_p = 130e-6_r8           ! [m] Medium soil particle diameter, assuming a global constant
                                                    ! of ~130 um following Leung et al. (2023)
   real(r8), parameter :: X = 10_r8                 ! [m] distance downwind of the roughness element (rock). Assume
                                                    ! estimating roughness effect at a distance of 10 m following Leung et al. (2023)
   real(r8), parameter :: b1 = 0.7_r8               ! [dimless] first fitting coefficient for the drag partition equation by Marticorena and Bergametti (1995), later modified by Darmenova et al. (2009).
   real(r8), parameter :: b2 = 0.8_r8               ! [dimless] second fitting coefficient for the drag partition equation by Marticorena and Bergametti (1995), later modified by Darmenova et al. (2009).
   character(len=*), parameter :: subname = 'PrigentRoughnessStream::CalcDragPartition'
   !---------------------------------------------------------------------

   ! Make sure we've initialized the Prigent roughness streams
   if ( .not. this%prigent_roughness_stream%IsStreamInit() )then
      call endrun(msg=subname//' ERROR: Streams have not been initialized, make sure Init is called first' &
                             //', and streams are on')
   end if

   ! dmleung: this loop calculates the drag partition effect (or roughness effect) of rocks.
   !          We save the drag partition factor as a patch level quantity.
   ! TODO: EBK 02/13/2024: Several magic numbers here that should become parameters so the meaning is preserved
   !z0s = 2_r8/30_r8 * D_p   ! equation for smooth roughness length for soil grain. See Danny M. Leung et al. (2023) and Martina Klose et al. (2021) for instance. 1/15 is a coefficient that relates roughness to soil particle diameter D_p.
   z0s = 2_r8 * D_p / 30_r8 ! equation for smooth roughness length for soil grain. See Danny M. Leung et al. (2023) and Martina Klose et al. (2021) for instance. 1/15 is a coefficient that relates roughness to soil particle diameter D_p.
                            ! Here we assume soil medium size is a global constant, and so is smooth roughness length.
   do p = bounds%begp,bounds%endp
      g = patch%gridcell(p)
      l = patch%landunit(p)
      if (lun%itype(l) /= istdlak) then
         ! Calculating rock drag partition factor using the Marticorena and Bergametti (1995) formulation.
         ! 0.01 is used to convert Prigent's roughness length dataset from centimeter to meter.
         this%dpfct_rock_patch(p) = 1._r8 - ( log(this%prigent_roughness_stream%prigent_rghn(g)*0.01_r8/z0s) &
                            / log(b1 * (X/z0s)**b2 ) )
      end if
   end do

 end subroutine CalcDragPartition

 !==============================================================================

end module DustEmisLeung2023