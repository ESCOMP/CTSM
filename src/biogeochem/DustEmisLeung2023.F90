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
  use shr_infnan_mod       , only : nan => shr_infnan_nan, assignment(=), shr_infnan_isnan
  use clm_varpar           , only : dst_src_nbr, ndst
  use clm_varcon           , only : grav, spval
  use landunit_varcon      , only : istcrop, istsoil
  use clm_varctl           , only : iulog
  use abortutils           , only : endrun
  use decompMod            , only : bounds_type, subgrid_level_landunit
  use atm2lndType          , only : atm2lnd_type
  use SoilStateType        , only : soilstate_type
  use CanopyStateType      , only : canopystate_type
  use WaterStateBulkType   , only : waterstatebulk_type
  use WaterDiagnosticBulkType, only : waterdiagnosticbulk_type
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
      real(r8), pointer, private :: wnd_frc_thr_it_patch      (:)   ! impact threshold (m/s)
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
      real(r8), pointer, private :: vai_Okin_patch            (:)   ! [m2 leaf /m2 land] LAI+SAI for calculating Okin drag partition
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
     ! Public for unit testing
     procedure , public  :: SetDragPartition  ! Set drag partitioning for testing

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
   integer :: begp,endp
   !------------------------------------------------------------------------

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
    integer :: begp,endp
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    this%dst_emiss_coeff_patch(begp:endp) = spval
    call hist_addfld1d (fname='DUST_EMIS_COEFF', units='dimensionless',  &
         avgflag='A', long_name='soil erodibility or dust emission coefficient for Kok emission scheme', &
         ptr_patch=this%dst_emiss_coeff_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%wnd_frc_thr_patch(begp:endp) = spval
    call hist_addfld1d (fname='WND_FRC_FT', units='m/s',  &
         avgflag='A', long_name='fluid threshold friction velocity', &
         ptr_patch=this%wnd_frc_thr_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%wnd_frc_thr_dry_patch(begp:endp) = spval
    call hist_addfld1d (fname='WND_FRC_FT_DRY', units='m/s',  &
         avgflag='A', long_name='dry fluid threshold friction velocity', &
         ptr_patch=this%wnd_frc_thr_dry_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%wnd_frc_thr_it_patch(begp:endp) = spval
    call hist_addfld1d (fname='WND_FRC_IT', units='m/s',  &
         avgflag='A', long_name='impact threshold friction velocity', &
         ptr_patch=this%wnd_frc_thr_it_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%wnd_frc_soil_patch(begp:endp) = spval
    call hist_addfld1d (fname='WND_FRC_SOIL', units='m/s',  &
         avgflag='A', long_name='soil surface wind friction velocity', &
         ptr_patch=this%wnd_frc_soil_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%lnd_frc_mble_patch(begp:endp) = spval
    call hist_addfld1d (fname='LND_FRC_MBLE', units='dimensionless',  &
         avgflag='A', long_name='land mobile fraction', &
         ptr_patch=this%lnd_frc_mble_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%gwc_patch(begp:endp) = spval
    call hist_addfld1d (fname='GWC', units='kg/kg',  &
         avgflag='A', long_name='gravimetric soil moisture at the topmost soil layer', &
         ptr_patch=this%gwc_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%liq_frac_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIQ_FRAC', units='dimensionless',  &
         avgflag='A', long_name='fraction of total water that is liquid', &
         ptr_patch=this%liq_frac_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%u_mean_slt_patch(begp:endp) = spval
    call hist_addfld1d (fname='U_S_MEAN', units='m/s',  &
         avgflag='A', long_name='mean wind velocity at saltation level', &
         ptr_patch=this%u_mean_slt_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%u_sd_slt_patch(begp:endp) = spval
    call hist_addfld1d (fname='U_S_SIGMA', units='m/s',  &
         avgflag='A', long_name='sd of wind velocity at saltation level', &
         ptr_patch=this%u_sd_slt_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%stblty_patch(begp:endp) = spval
    call hist_addfld1d (fname='ZETAOBU', units='',  &
         avgflag='A', long_name='stability parameter', &
         ptr_patch=this%stblty_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%u_fld_thr_patch(begp:endp) = spval
    call hist_addfld1d (fname='U_FT', units='m/s',  &
         avgflag='A', long_name='fluid threshold velocity at saltation level', &
         ptr_patch=this%u_fld_thr_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%u_impct_thr_patch(begp:endp) = spval
    call hist_addfld1d (fname='U_IT', units='m/s',  &
         avgflag='A', long_name='impact threshold velocity at saltation level', &
         ptr_patch=this%u_impct_thr_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%thr_crs_rate_patch(begp:endp) = spval
    call hist_addfld1d (fname='ALPHA_TC_RATE', units='',  &
         avgflag='A', long_name='threshold crossing rate', &
         ptr_patch=this%thr_crs_rate_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%prb_crs_fld_thr_patch(begp:endp) = spval
    call hist_addfld1d (fname='P_FT', units='',  &
         avgflag='A', long_name='probability of winds crossing fluid threshold', &
         ptr_patch=this%prb_crs_fld_thr_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%prb_crs_impct_thr_patch(begp:endp) = spval
    call hist_addfld1d (fname='P_IT', units='',  &
         avgflag='A', long_name='probability of winds crossing impact threshold', &
         ptr_patch=this%prb_crs_impct_thr_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%intrmtncy_fct_patch(begp:endp) = spval
    call hist_addfld1d (fname='ETA', units='',  &
         avgflag='A', long_name='intermittency factor', &
         ptr_patch=this%intrmtncy_fct_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%ssr_patch(begp:endp) = spval
    call hist_addfld1d (fname='SSR', units='m/s',  &
         avgflag='A', long_name='Okin-Pierre vegetation shear stress ratio (drag partition factor)', &
         ptr_patch=this%ssr_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%vai_Okin_patch(begp:endp) = spval
    call hist_addfld1d (fname='VAI_OKIN', units='m/s',  &
         avgflag='A', long_name='vegetation area index used in the Okin-Pierre plant drag partition scheme', &
         ptr_patch=this%vai_Okin_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%frc_thr_rghn_fct_patch(begp:endp) = spval
    call hist_addfld1d (fname='FRC_THR_RGHN_FCT', units='dimensionless',  &
         avgflag='A', long_name='hybrid drag partition (or roughness) factor', &
         ptr_patch=this%frc_thr_rghn_fct_patch, set_lake=0.0_r8, set_urb=0.0_r8)
    this%wnd_frc_thr_std_patch(begp:endp) = spval
    call hist_addfld1d (fname='WND_FRC_FT_STD', units='m/s',  &
         avgflag='A', long_name='standardized fluid threshold friction velocity', &
         ptr_patch=this%wnd_frc_thr_std_patch, set_lake=0.0_r8, set_urb=0.0_r8)
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
    ! Calculate Drag Partition factor (Marticorena and Bergametti 1995 formulation)
    if ( this%prigent_roughness_stream%useStreams() )then !if usestreams == true, and it should be always true
      call this%CalcDragPartition( bounds )
    else

      call endrun( "ERROR:: dus_emis_Leung_2023 requires requires a streams file of aeolian roughness length to calculate drag partitioning" )

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
    real(r8) :: flx_mss_hrz_slt_ttl
    real(r8) :: flx_mss_vrt_dst_ttl(bounds%begp:bounds%endp)
    real(r8) :: frc_thr_wet_fct
    real(r8) :: frc_thr_rgh_fct
    real(r8) :: wnd_frc_slt
    real(r8) :: lnd_frc_mbl(bounds%begp:bounds%endp)
    real(r8) :: bd
    real(r8) :: gwc_sfc
    real(r8) :: ttlai(bounds%begp:bounds%endp)
    real(r8) :: tlai_lu(bounds%begl:bounds%endl)
    real(r8) :: sumwt(bounds%begl:bounds%endl) ! sum of weights
    logical  :: found                          ! temporary for error check
    integer  :: index

    real(r8) :: tmp2                ! calculates the dry fluid threshold using Shao and Lu (2000) scheme;
                                    ! replace the tmp1 (Iversen and White, 1982) that was passed from Dustini to DustEmission;
                                    ! tmp2 will be calculated here
    real(r8) :: wnd_frc_thr_slt_std ! [m/s] The soil threshold friction speed at standard air density (1.2250 kg/m3)
    real(r8) :: frag_expt           ! fragmentation exponent
    real(r8) :: wnd_frc_thr_slt_it  ! [m/s] created for impact threshold friction velocity
    real(r8) :: wnd_frc_thr_slt     ! [m/s] used for wet fluid threshold friction velocity
    real(r8) :: K_length            ! [dimless] normalized mean interobstacle distance, or called gap length (Okin, 2008)
    ! dmleung has these variables and will change them into pointers and prepare for their history outputs. 30 Sep 2024
    real(r8) :: bare_frc            ! LUH2 bare soil land cover fraction
    real(r8) :: natveg_frc          ! LUH2 natural vegetation cover fraction
    real(r8) :: crop_frc            ! LUH2 crop cover fraction.
    !
    ! constants
    !
    real(r8), parameter :: vai_mbl_thr = 0.6_r8        ! [m2 m-2] new VAI threshold; Danny M. Leung suggests something between 0.6 and 1 for tuning. Zender's scheme uses 0.3. Simone Tilmes might want this as a namelist variable for easier CESM tuning. dmleung 30 Sep 2024.

    real(r8), parameter :: Cd0 = 4.4e-5_r8             ! [dimless] proportionality constant in calculation of dust emission coefficient
    real(r8), parameter :: Ca = 2.7_r8                 ! [dimless] proportionality constant in scaling of dust emission exponent
    real(r8), parameter :: Ce = 2.0_r8                 ! [dimless] proportionality constant scaling exponential dependence of dust emission coefficient on standardized soil threshold friction speed
    real(r8), parameter :: C_tune = 0.05_r8            ! [dimless] global tuning constant for vertical dust flux; set to produce ~same global dust flux in control sim (I_2000) as old parameterization
    real(r8), parameter :: wnd_frc_thr_slt_std_min = 0.16_r8 ! [m/s] minimum standardized soil threshold friction speed
    real(r8), parameter :: forc_rho_std = 1.2250_r8    ! [kg/m3] density of air at standard pressure (101325) and temperature (293 K)
    real(r8), parameter :: dns_slt = 2650.0_r8         ! [kg m-3] Density of optimal saltation particles
    real(r8), parameter :: B_it = 0.82_r8              ! [dimless] ratio = u_star_it / u_star_ft0
    real(r8), parameter :: k = 0.4_r8                  ! [dimless] von Karman constant
    real(r8), parameter :: f_0 = 0.32_r8               ! [dimless] SSR in the immediate lee of a plant
    real(r8), parameter :: c_e = 4.8_r8                ! [dimless] e-folding distance velocity recovery
    real(r8), parameter :: D_p = 130e-6_r8             ! [m] Medium soil particle diameter, assuming a global constant of ~130 um following Leung et al. (2023). dmleung 16 Feb 2024
    real(r8), parameter :: gamma_Shao = 1.65e-4_r8     ! [kg s-2] interparticle cohesion: fitting parameter in Shao and Lu (2000) (S&L00). dmleung 16 Feb 2024
    real(r8), parameter :: A_Shao = 0.0123_r8          ! [dimless] coefficient for aerodynamic force: fitting parameter in Shao and Lu (2000). dmleung 16 Feb 2024
    real(r8), parameter :: frag_expt_thr = 2.5_r8      ! [dimless] Maximum value or threshold for fragmentation exponent defined in Leung et al. (2023). Danny M. Leung suggested it to be somewhere between 3 and 5 for tuning. It is used to prevent a local AOD blowup (over Patagonia, Argentina), but one can test larger values and relax the threshold if wanted. dmleung 16 Feb 2024. Update: Simone Tilmes might want this as a namelist variable for easier CESM tuning. 30 Sep 2024.
    real(r8), parameter :: z0a_glob = 1e-4_r8          ! [m] assumed globally constant aeolian roughness length value in Leung et al. (2023), for the log law of the wall for Comola et al. (2019) intermittency scheme. dmleung 20 Feb 2024
    real(r8), parameter :: hgt_sal = 0.1_r8            ! [m] saltation height used by Comola et al. (2019) intermittency scheme for the log law of the wall. dmleung 20 Feb 2024
    real(r8), parameter :: vai0_Okin = 0.1_r8          ! [m2/m2] minimum VAI needed for Okin-Pierre's vegetation drag partition equation. lai=0 in the equation will lead to infinity, so a small value is added into this lai dmleung defined.
    real(r8), parameter :: zii = 1000.0_r8             ! [m] convective boundary layer height added by dmleung 20 Feb 2024, following other CTSM modules (e.g., CanopyFluxesMod). Should we transfer PBL height (PBLH) from CAM?
    real(r8), parameter :: dust_veg_drag_fact = 0.7_r8 ! [dimless] dmleung added a tuning factor for Greg Okin's vegetation drag partition effect. dmleung suggested a smaller vegetation drag partition effect given an increase in vegetation roughness after CTSM switched from using ZengWang2007 to Meier2022. This is simply because the drag partition effect should decrease with increasing roughness, but Okin's scheme is only a function of LAI. One might want to change this factor to 1_r8 when using ZengWang2007. dmleung 30 Sep 2024
    real(r8) :: numer                                  ! Numerator term for threshold crossing rate
    real(r8) :: denom                                  ! Denominator term for threshold crossing rate
    !------------------------------------------------------------------------

    associate(                                                         &
         forc_rho            => atm2lnd_inst%forc_rho_downscaled_col , & ! Input:  [real(r8) (:)   ]  downscaled density (kg/m**3)

         gwc_thr             => soilstate_inst%gwc_thr_col           , & ! Input:  [real(r8) (:)   ]  threshold gravimetric soil moisture based on clay content
         mss_frc_cly_vld     => soilstate_inst%mss_frc_cly_vld_col   , & ! Input:  [real(r8) (:)   ]  [frc] Mass fraction clay limited to 0.20
         watsat              => soilstate_inst%watsat_col            , & ! Input:  [real(r8) (:,:) ]  saturated volumetric soil water

         tlai                => canopystate_inst%tlai_patch          , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index, no burying by snow
         tsai                => canopystate_inst%tsai_patch          , & ! Input:  [real(r8) (:)   ]  one-sided stem area index, no burying by snow

         frac_sno            => waterdiagnosticbulk_inst%frac_sno_col, & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
         h2osoi_vol          => waterstatebulk_inst%h2osoi_vol_col   , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat)
         h2osoi_liq          => waterstatebulk_inst%h2osoi_liq_col   , & ! Input:  [real(r8) (:,:) ]  liquid soil water (kg/m2)
         h2osoi_ice          => waterstatebulk_inst%h2osoi_ice_col   , & ! Input:  [real(r8) (:,:) ]  frozen soil water (kg/m2)

         fv                  => frictionvel_inst%fv_patch            , & ! Input:  [real(r8) (:)   ]  friction velocity (m/s) (for dust model)
         obu                 => frictionvel_inst%obu_patch           , & ! Input:  [real(r8) (:)   ] Monin-Obukhov length from the friction Velocity module          obu                 => frictionvel_inst%obu_patch             & ! Input:  [real(r8) (:)   ] Monin-Obukhov length from the friction Velocity module 

         dpfct_rock          => this%dpfct_rock_patch                , & ! Input:  rock drag partition factor defined in Marticorena and Bergametti 1995. A fraction between 0 and 1.

         flx_mss_vrt_dst     => this%flx_mss_vrt_dst_patch           , & ! Output: [real(r8) (:,:) ]  surface dust emission (kg/m**2/s)
         flx_mss_vrt_dst_tot => this%flx_mss_vrt_dst_tot_patch       , & ! Output: [real(r8) (:)   ]  total dust flux back to atmosphere (pft)
         ! below variables are defined in Kok et al. (2014) or (mostly) Leung et al. (2023) dust emission scheme. dmleung 16 Feb 2024
         dst_emiss_coeff     => this%dst_emiss_coeff_patch           , & ! Output: dust emission coefficient
         wnd_frc_thr         => this%wnd_frc_thr_patch               , & ! Output: fluid threshold
         wnd_frc_thr_dry     => this%wnd_frc_thr_dry_patch           , & ! Output: dry fluid threshold
         wnd_frc_thr_it      => this%wnd_frc_thr_it_patch            , & ! Output: impact threshold
         lnd_frc_mble        => this%lnd_frc_mble_patch              , & ! Output: bare land fraction
         wnd_frc_soil        => this%wnd_frc_soil_patch              , & ! Output: soil friction velocity u_*s = (u_*)(f_eff)
         gwc                 => this%gwc_patch                       , & ! output: gravimetric water content
         liq_frac            => this%liq_frac_patch                  , & ! Output: fraction of liquid moisture
         intrmtncy_fct       => this%intrmtncy_fct_patch             , & ! Output: intermittency factor eta (fraction of time that dust emission is active within a timestep)
         stblty              => this%stblty_patch                    , & ! Output: stability in similarity theory (no need to output)
         u_mean_slt          => this%u_mean_slt_patch                , & ! Output: mean wind speed at 0.1 m height translated from friction velocity using the log law of the wall, assuming neutral condition
         u_sd_slt            => this%u_sd_slt_patch                  , & ! Output: standard deviation of wind speed from similarity theory
         u_fld_thr           => this%u_fld_thr_patch                 , & ! Output: fluid threshold wind speed at 0.1 m height translated from the log law of the wall
         u_impct_thr         => this%u_impct_thr_patch               , & ! Output: impact threshold wind speed at 0.1 m height translated from the log law of the wall
         thr_crs_rate        => this%thr_crs_rate_patch              , & ! Output: threshold crossing rate in Comola 2019 intermittency parameterization
         prb_crs_fld_thr     => this%prb_crs_fld_thr_patch           , & ! Output: probability of instantaneous wind crossing fluid threshold in Comola 2019 intermittency parameterization
         prb_crs_impct_thr   => this%prb_crs_impct_thr_patch         , & ! Output: probability of instantaneous wind crossing impact threshold in Comola 2019 intermittency parameterization
         ssr                 => this%ssr_patch                       , & ! Output: vegetation drag partition factor in Okin 2008 vegetation roughness effect (called shear stress ratio, SSR in Okin 2008)
         vai_Okin            => this%vai_Okin_patch                  , & ! Output: vegetation area index for calculating Okin-Pierre vegetation drag partitioning. vai=0 in the ssr equation will lead to infinity, so a small value is added into this vai dmleung defined. (no need to output) 16 Feb 2024
         frc_thr_rghn_fct    => this%frc_thr_rghn_fct_patch          , & ! Output: hybrid/total drag partition factor considering both rock and vegetation drag partition factors.
         wnd_frc_thr_std     => this%wnd_frc_thr_std_patch             & ! Output: standardized dust emission threshold friction velocity defined in Jasper Kok et al. (2014).
         )

      ttlai(bounds%begp : bounds%endp) = 0.0_r8
      ! make lai average at landunit level
      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         ttlai(p) = tlai(p)+tsai(p)
      enddo

      tlai_lu(bounds%begl : bounds%endl) = spval
      sumwt(bounds%begl : bounds%endl) = 0.0_r8
      do p = bounds%begp,bounds%endp
         if (ttlai(p) /= spval .and. patch%active(p) .and. patch%wtlunit(p) /= 0.0_r8) then
            c = patch%column(p)
            l = patch%landunit(p)
            if (sumwt(l) == 0.0_r8) tlai_lu(l) = 0.0_r8
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
         else if (sumwt(l) /= 0.0_r8) then
            tlai_lu(l) = tlai_lu(l)/sumwt(l)
         end if
      end do
      if (found) then
         write(iulog,*) 'error: sumwt is greater than 1.0 at l= ',index
         call endrun(subgrid_index=index, subgrid_level=subgrid_level_landunit, msg=errMsg(sourcefile, __LINE__))
         return
      end if

      ! Loop through patches

      ! initialize variables which get passed to the atmosphere
      flx_mss_vrt_dst(bounds%begp:bounds%endp,:)=0.0_r8

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

      ! dmleung add output for bare_frc and veg_frc here if wanted !!!!----------------------

      ! reset history output variables before next if-statement to avoid output = inf

      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         flx_mss_vrt_dst_tot(p) = 0.0_r8
         dst_emiss_coeff(p) = 0.0_r8
         wnd_frc_thr(p) = 0.0_r8
         wnd_frc_thr_dry(p) = 0.0_r8
         lnd_frc_mble(p) = 0.0_r8
         wnd_frc_soil(p) = 0.0_r8
         gwc(p) = 0.0_r8
         liq_frac(p) = 0.0_r8
         u_mean_slt(p) = 0.0_r8
         u_sd_slt(p) = 0.0_r8
         stblty(p)   = 0.0_r8
         u_fld_thr(p) = 0.0_r8
         u_impct_thr(p) = 0.0_r8
         thr_crs_rate(p) = 0.0_r8
         prb_crs_fld_thr(p) = 0.0_r8
         prb_crs_impct_thr(p) = 0.0_r8
         intrmtncy_fct(p) = 0.0_r8
         ssr(p) = 0.0_r8
         vai_Okin(p) = 0.0_r8
         frc_thr_rghn_fct(p) = 0.0_r8
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

         !--------------------------------------------------------------------------------------------------
         ! put dust emission calculation here to output threshold friction velocity for the whole globe,
         ! not just when lnd_frc_mbl = 0. Danny M. Leung 27 Nov 2021
          
         !####################################################################################################
         ! calculate soil moisture effect for dust emission threshold
         ! following Fecan, Marticorena et al. (1999)
         ! also see Zender et al. (2003) for DUST emission scheme and Kok et al. (2014b) for K14 emission scheme in CESM
         bd = (1.0_r8-watsat(c,1))*dns_slt      ![kg m-3] Bulk density of dry surface soil (dmleung changed from 2700 to dns_slt, soil particle density, on 16 Feb 2024. Note that dn s_slt=2650 kg m-3 so the value is changed by a tiny bit from 2700 to 2650. dns_slt has been here for many years so dns_slt should be used here instead of explicitly typing the value out. dmleung 16 Feb 2024)

         ! use emission threshold to calculate standardized threshold and dust emission coefficient

         ! Here convert h2osoi_vol (H2OSOI) at the topmost CTSM soil layer from volumetric (m3 water / m3 soil) to gravimetric soil moisture (kg water / kg soil)
         gwc_sfc = h2osoi_vol(c,1)*SHR_CONST_RHOFW/bd    ![kg kg-1] Gravimetric H2O cont
         if (gwc_sfc > gwc_thr(c)) then
            frc_thr_wet_fct = sqrt(1.0_r8 + 1.21_r8 * (100.0_r8*(gwc_sfc - gwc_thr(c)))**0.68_r8)   ! dmleung's comment: this is an empirical equation by Fecan, Marticorena et al. (1999) on relating the soil moisture factor on enhancing dust emission threshold to gravimetric soil moisture. 1.21 and 0.68 are fitting parameters in the regression done by Fecan; 100 is to convert gracimetric soil moisture from fraction (kg water / kg soil) to percentage. Note that gwc_thr was defined in SoilStateInitConst.F90 as a fraction. 1.0_r8 means there is no soil moisture effect on enhancing dust emission threhsold. dmleung 16 Feb 2024.
         else
            frc_thr_wet_fct = 1.0_r8
         end if

         ! output moisture variables
         gwc(p) = gwc_sfc     ! output surface gravimetric water content

         ! slevis: adding liqfrac here, because related to effects from soil water
         liqfrac = max( 0.0_r8, min( 1.0_r8, h2osoi_liq(c,1) / (h2osoi_ice(c,1)+h2osoi_liq(c,1)+1.0e-6_r8) ) )
         ! dmleung: output liquid fraction
         liq_frac(p) = liqfrac

         !#######################################################################################################
         ! calculate Shao & Lu (2000) dust emission threshold scheme here
         ! use tmp1 from DUSTini for Iversen and White I&W (1982) (~75 um is optimal); use tmp2 for S&L (2000) (~80 um is optimal)
         ! see Danny M. Leung et al. (2023)
         !#######################################################################################################
         tmp2 = sqrt(A_Shao * (dns_slt*grav*D_p + gamma_Shao/D_p))  ! calculate S&L (2000) scheme here for threshold; gamma = 1.65e-4 following S&L00, D_p = 127 um ~ 130 um following Leung et al. (2023). dmleung use defined parameters instead of typing numerical values 16 Feb 2024
         wnd_frc_thr_dry(p) = tmp2 / sqrt(forc_rho(c))    ! dry fluid threshold
         wnd_frc_thr_slt = tmp2 / sqrt(forc_rho(c)) * frc_thr_wet_fct !* frc_thr_rgh_fct   ! fluid threshold. dmleung commented out frc_thr_rgh_fct since it is used to modify the wind, not the wind threshold.
         wnd_frc_thr_slt_it = B_it * tmp2 / sqrt(forc_rho(c)) ! define impact threshold

         ! the above formula is true for Iversen and White (1982) and Shao and Lu (2000) scheme
         wnd_frc_thr(p) = wnd_frc_thr_slt          ! output fluid threshold
         wnd_frc_thr_it(p) = wnd_frc_thr_slt_it    ! output impact threshold

         !##############################################################################################
         ! dmleung: here, calculate quantities relevant to the fluid threshold
         ! standardized fluid threshold
         wnd_frc_thr_slt_std = wnd_frc_thr_slt * sqrt(forc_rho(c) / forc_rho_std) ! standardized soil threshold friction speed (defined using fluid threshold)
         wnd_frc_thr_std(p) = wnd_frc_thr_slt_std          ! output standardized fluid threshold
         ! dust emission coefficient or soil erodibility coefficient (this is analogous to the soil erodibility map or prefenertial source filter in Zender; see zendersoilerodstream)
         dst_emiss_coeff(p) = Cd0 * exp(-Ce * (wnd_frc_thr_slt_std - wnd_frc_thr_slt_std_min) / wnd_frc_thr_slt_std_min) ! save dust emission coefficient here for all grids

         ! framentation exponent (dependent on fluid threshold)
         frag_expt = (Ca * (wnd_frc_thr_slt_std - wnd_frc_thr_slt_std_min) / wnd_frc_thr_slt_std_min)  ! fragmentation exponent, defined in Kok et al. (2014a)
         if (frag_expt > frag_expt_thr) then   ! set fragmentation exponent to be 3 or 5 at maximum, to avoid local AOD blowup
            frag_expt = frag_expt_thr
         end if

         !##############################################################################################
         !################ drag partition effect, and soil-surface friction velocity ###################
         ! subsection on computing vegetation drag partition and hybrid drag partition factors
         ! in Leung et al. (2023), drag partition effect is applied on the wind instead of the threshold
         !##############################################################################################
         ! the following comes from subr. frc_thr_rgh_fct_get
         ! purpose: compute factor by which surface roughness increases threshold
         !          friction velocity (currently a constant)

         if (lnd_frc_mbl(p) > 0.0_r8  .AND. tlai_lu(l)<= vai_mbl_thr) then

            vai_Okin(p) = tlai_lu(l)+vai0_Okin       ! LAI+SAI averaged to landunit level; the equation is undefined at lai=0, and LAI in CTSM has some zeros over deserts, so we add in a small number.
            if (vai_Okin(p) > 1.0_r8) then
               vai_Okin(p)  = 1.0_r8   ! setting LAI = 1 to be a max value (since K_length goes to negative when LAI>1) 
            end if


            ! calculate Okin's shear stress ratio (SSR, which is vegetation drag partition factor) using Pierre's equation
            K_length = 2.0_r8 * (1.0_r8/vai_Okin(p) - 1.0_r8)   ! Here VAI has to be non-zero to avoid blowup, and < 1 to avoid -ve K_length. See this equation in Leung et al. (2023). This line is Okin's formulation
            ssr(p) = dust_veg_drag_fact * (K_length+f_0*c_e)/(K_length+c_e) ! see this equation in Caroline Pierre et al. (2014) or Leung et al. (2023). This line is Pierre's formulation. dmleung added a tuning factor for Okin's vegetation drag partition effect (SSR) on 30 Sep 2024.

            ! calculation of the hybrid/total drag partition effect considering both rock and vegetation drag partitioning using LUH2 bare and veg fractions within a grid
            if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
               if (patch%itype(p) == noveg) then ! if bare, uses rock drag partition factor
                  if (shr_infnan_isnan(dpfct_rock(p)) ) then ! dmleung added 24 May 2024: dpfct_rock(p) could be NaN; CLM could run when DEBUG=FALSE in env_build.xml but dies when DEBUG=TRUE (usually when checking if wnd_frc_slt > wnd_frc_thr_slt_it and if numer/denom < 30 below)
                     frc_thr_rgh_fct = 0.001_r8 ! Set drag partition effect to be a very small value (or zero) such that there is no emission whenever dpfct_rock(p) = NaN; dmleung 24 May 2024
                  else
                     frc_thr_rgh_fct = dpfct_rock(p)
                  end if
               else   ! if vegetation, uses vegetation drag partition factor
                  frc_thr_rgh_fct = ssr(p)
               end if
            else
               frc_thr_rgh_fct = 1.0_r8
            end if

            wnd_frc_slt = fv(p) * frc_thr_rgh_fct   ! wnd_frc_slt is the drag-parition-modified wind speed and will be used in the dust emission equation below

            frc_thr_rghn_fct(p) = frc_thr_rgh_fct   ! save and output hybrid drag partition factor

         else    ! for lnd_frc_mbl=0, do not change friction velocity and assume drag partition factor = 0
            wnd_frc_slt = fv(p) * frc_thr_rghn_fct(p) ! The value here is not important since once lnd_frc_mbl(p) <= 0.0_r8 there will be no emission. dmleung added 5 Jul 2024
            frc_thr_rghn_fct(p) = 0.0_r8            ! When LAI > vai_mbl_thr, the drag partition effect is zero. dmleung 16 Feb 2024.
         end if

         !########## end of drag partition effect #######################################################

         ! save soil friction velocity and roughness effect before the if-statement
         wnd_frc_soil(p) = wnd_frc_slt  ! save soil friction velocity for CLM output, which has drag partition and Owen effect
         ! 20 Feb 2024: dmleung notes that Leung does not consider the Owen's effect. This is Jasper Kok's decision. The Owen's effect should be in Zender's DUST emission scheme.

         ! save land mobile fraction
         lnd_frc_mble(p) = lnd_frc_mbl(p)  ! save land mobile fraction first, before the if-statement
         ! only perform the following calculations if lnd_frc_mbl is non-zero

         if (lnd_frc_mbl(p) > 0.0_r8) then  ! if bare land fraction is larger than 0 then calculate the dust emission equation

            ! reset these variables which will be updated in the following if-block

            flx_mss_hrz_slt_ttl = 0.0_r8
            flx_mss_vrt_dst_ttl(p) = 0.0_r8

            ! the following comes from subr. flx_mss_hrz_slt_ttl_Whi79_get
            ! purpose: compute vertically integrated streamwise mass flux of particles

            if (wnd_frc_slt > wnd_frc_thr_slt_it) then! if using Leung's scheme, use impact threshold for dust emission equation; if Zender, uses fluid threshold (wnd_frc_thr_slt) for dust emission equation

               !################### for Leung et al. (2023) ################################################
               ! dmleung: instead of using mss_frc_cly_vld(c) with a range of [0,0.2] , which makes dust too sensitive to input clay surface dataset), for now use 0.1 + mss_frc_cly_vld(c) * 0.1 / 0.20 with a range of [0.1,0.2]. So, instead of scaling dust emission to 1/20 times for El Djouf (western Sahara) because of its 1 % clay fraction, scale its emission to 1/2 times. This reduces the sensitivity of dust emission to the clay input dataset. In particular, because dust emission is a small-scale process and the grid-averaged clay from the new WISE surface dataset over El Djouf is 1 %, much lower than the former FAO estimation of 5 %, dmleung is not sure if the clay value of 1 % suppresses too much of El Djouf's small-scale dust emission process. Similar to Charlie Zender's feeling suspicious about the soil texture datasets (or the dust emission schemes' sandblasting process), dmleung feels the same and for now decides to still allow dust emission  to weakly scale with clay fraction, however limiting the scaling factor to 0.1-0.2. See modification in SoilStateInitTimeConst.F90. dmleung 5 Jul 2024
               flx_mss_vrt_dst_ttl(p) = dst_emiss_coeff(p) * mss_frc_cly_vld(c) * forc_rho(c) * ((wnd_frc_slt**2.0_r8 - wnd_frc_thr_slt_it**2.0_r8) / wnd_frc_thr_slt_it) * (wnd_frc_slt / wnd_frc_thr_slt_it)**frag_expt  ! Leung et al. (2022) uses Kok et al. (2014) dust emission euqation for emission flux

               ! account for bare soil fraction, frozen soil fraction, and apply global tuning parameter (Kok et al. 2014)
               flx_mss_vrt_dst_ttl(p) = flx_mss_vrt_dst_ttl(p) * lnd_frc_mbl(p) * C_tune * liqfrac
               !########################################################################################
            end if

            !############## Danny M. Leung added the intermittency calculation #################################
            ! subsection for intermittency factor calculation (only used by Leung's scheme, not Zender's scheme)
            ! Leung et al. (2023) uses the Comola et al. (2019) intermittency scheme for the calculation of intermittent dust emissions.
            ! This part takes care of the sub-timestep, high-frequency (< 1 minute period) turblent wind fluctuations occuring at the planetary boundary layer (PBL) near surface. Subtimestep wind gusts and episodes are important for generating emissions in marginal dust source regions, such as semiarid areas and high-latitude polar deserts.
            ! 2 Dec 2021: assume no buoyancy contribution to the wind fluctuation (u_sd_slt), so no obul(p) is needed. It is shown to be important for the wind fluctuations contribute little to the intermittency factor. We might add this back in the future revisions.
            ! 20 Feb 2024: dmleung notes that dmleung may revise Comola's scheme in the future to improve Comola's formulation of the statistical parameterization.

            ! mean lowpass-filtered wind speed at hgt_sal = 0.1 m saltation height (assuming aerodynamic roughness length z0a_glob = 1e-4 m globally for ease; also assuming neutral condition)
            u_mean_slt(p) = (wnd_frc_slt/k) * log(hgt_sal / z0a_glob)  ! translating from ustar (velocity scale) to actual wind

            if ( obu(p) == 0.0_r8 )then
               call endrun(msg='Input obu is zero and can NOT be' )
               return
            end if
            stblty(p) = zii / obu(p)   ! -dmleung 20 Feb 2024: use obu from CTSM and PBL height = zii (= 1000_r8) which is default in CTSM. Should we transfer PBL height from CAM?
            if ((12.0_r8 - 0.5_r8 * stblty(p)) .GE. 0.001_r8) then ! should have used 0 theoretically; used 0.001 here to avoid undefined values
               u_sd_slt(p) = wnd_frc_slt * (12.0_r8 - 0.5_r8 * stblty(p))**0.333_r8
            else
               u_sd_slt(p) = wnd_frc_slt * (0.001_r8)**0.333_r8   ! should have used 0 theoretically; used 0.001 here to avoid undefined values
            end if

            ! threshold velocities
            ! Here wnd_frc_thr_slt is the fluid threshold; wnd_frc_thr_dry(p) is the dry fluid threshold; B_it*wnd_frc_thr_dry(p) is the impact threshold
            ! fluid threshold wind at 0.1 m saltation height
            u_fld_thr(p) = (wnd_frc_thr_slt/k) * log(hgt_sal / z0a_glob)  ! assume a globally constant z0a value for the log law of the wall, but it can be z0m from CLM or, better, z0a from Prigent's roughness dataset. Danny M. Leung et al. (2023) chose to assume a global constant z0a = 1e-4 m. dmleung 20 Feb 2024
            ! impact threshold wind at 0.1 m saltation height
            u_impct_thr(p) = (wnd_frc_thr_slt_it/k) * log(hgt_sal / z0a_glob)

            ! threshold crossing rate
            numer = (u_fld_thr(p)**2.0_r8 - u_impct_thr(p)**2.0_r8 - 2.0_r8 * u_mean_slt(p) * (u_fld_thr(p) - u_impct_thr(p)))
            denom = (2.0_r8 * u_sd_slt(p)**2.0_r8) ! note that u_sd_slt should be always positive
            ! Truncate to zero if the expression inside exp is becoming too large
            if ( numer/denom < 30.0_r8 ) then  ! set numer/denom to be < 30 given exp(30) below is already very large; also denom itself should be non-zero and non-negative given the standard deviation (u_sd_slt) of the subtimestep wind fluctuation is non-negative. dmleung 28 May 2024
               thr_crs_rate(p) = (exp((u_fld_thr(p)**2.0_r8 - u_impct_thr(p)**2.0_r8 - 2.0_r8 * u_mean_slt(p) * (u_fld_thr(p) - u_impct_thr(p))) / (2.0_r8 * u_sd_slt(p)**2.0_r8)) + 1.0_r8)**(-1.0_r8)
            else
               thr_crs_rate(p) = 0.0_r8
            end if

            ! probability that lowpass-filtered wind speed does not exceed u_ft
            prb_crs_fld_thr(p) = 0.5_r8 * (1.0_r8 + erf((u_fld_thr(p) - u_mean_slt(p)) / ( sqrt(2.0_r8) * u_sd_slt(p))))
            ! probability that lowpass-filtered wind speed does not exceed u_it
            prb_crs_impct_thr(p) = 0.5_r8 * (1.0_r8 + erf((u_impct_thr(p) - u_mean_slt(p)) / ( sqrt(2.0_r8) * u_sd_slt(p))))

            ! intermittency factor (eta; ranging from 0 to 1)
            intrmtncy_fct(p) = 1.0_r8 - prb_crs_fld_thr(p) + thr_crs_rate(p) * (prb_crs_fld_thr(p) - prb_crs_impct_thr(p))

            ! multiply dust emission flux by intermittency factor
            if ( shr_infnan_isnan(intrmtncy_fct(p)) ) then  ! if intrmtncy_fct(p) is not NaN then multiply by intermittency factor; this statement is needed because dust emission flx_mss_vrt_dst_ttl(p) has to be non NaN (at least zero) to be outputted
               flx_mss_vrt_dst_ttl(p) = flx_mss_vrt_dst_ttl(p)
            else
               flx_mss_vrt_dst_ttl(p) = flx_mss_vrt_dst_ttl(p) * intrmtncy_fct(p)  ! multiply dust flux by intermittency
            end if

            !############ end the intermittency subsection here; only use for Leung's scheme  ##########################

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
   !---------------------------------------------------------------------

   ! Make sure we've initialized the Prigent roughness streams
   if ( .not. this%prigent_roughness_stream%IsStreamInit() )then
      write(iulog,*)'Error : Prigent roughness stream is NOT on: ', errMsg(sourcefile, __LINE__)
      call endrun(msg=' ERROR: Streams have not been initialized, make sure Init is called first' &
                       //', and streams are on')
      return
   end if

   ! dmleung: this loop calculates the drag partition effect (or roughness effect) of rocks.
   !          We save the drag partition factor as a patch level quantity.
   ! TODO: EBK 02/13/2024: Several magic numbers here that should become parameters so the meaning is preserved
   z0s = 2.0_r8 * D_p / 30.0_r8 ! equation for smooth roughness length for soil grain. See Danny M. Leung et al. (2023) and Martina Klose et al. (2021) for instance. 1/15 is a coefficient that relates roughness to soil particle diameter D_p.
                            ! Here we assume soil medium size is a global constant, and so is smooth roughness length.
   do p = bounds%begp,bounds%endp
      g = patch%gridcell(p)
      l = patch%landunit(p)
      if (lun%itype(l) /= istdlak) then
         ! Calculating rock drag partition factor using the Marticorena and Bergametti (1995) formulation.
         ! 0.01 is used to convert Prigent's roughness length dataset from centimeter to meter.
         this%dpfct_rock_patch(p) = 1.0_r8 - ( log(this%prigent_roughness_stream%prigent_rghn(g)*0.01_r8/z0s) &
                            / log(b1 * (X/z0s)**b2 ) )
      end if
   end do

 end subroutine CalcDragPartition


  !------------------------------------------------------------------------

 subroutine SetDragPartition(this, bounds, drag_partition)
   !
   ! !DESCRIPTION:
   ! Set the drag partition for testing
   !
   ! !USES:
   !
   ! !ARGUMENTS:
   implicit none
   class (dust_emis_leung2023_type) :: this
   type(bounds_type), intent(in) :: bounds
   real(r8), intent(in) :: drag_partition
   !
   ! !LOCAL VARIABLES:
   integer  :: p     ! Indices

   !---------------------------------------------------------------------

   do p = bounds%begp,bounds%endp
      this%dpfct_rock_patch(p) = drag_partition
   end do

 end subroutine SetDragPartition

 !==============================================================================

end module DustEmisLeung2023
