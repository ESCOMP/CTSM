module LunaMod

#include "shr_assert.h"
  
  !********************************************************************************************************************************************************************** 
  ! !DESCRIPTION:
  ! Calculates the photosynthetic capacities based on a prescribed leaf nitrogen content, using the LUNA model, developed by Chonggang Xu, Ashehad Ali and Rosie Fisher
  ! Currently only works for C3 plants. See Xu et al 2012; Ali et al 2015a. Ecological Applications. http://dx.doi.org/10.1890/14-2111.1.  and Ali et al 2015b.In Review.
  ! !USES:
  use shr_kind_mod          , only : r8  => shr_kind_r8
  use shr_log_mod           , only : errMsg  => shr_log_errMsg
  use clm_varcon            , only : rgas, tfrz,spval
  use abortutils            , only : endrun
  use clm_varctl            , only : iulog
  use clm_varcon            , only : namep 
  use clm_varpar            , only : nlevcan
  use decompMod             , only : bounds_type
  use pftconMod             , only : pftcon
  use FrictionvelocityMod   , only : frictionvel_type 
  use atm2lndType           , only : atm2lnd_type
  use CanopyStateType       , only : canopystate_type
  use PhotosynthesisMod     , only : photosyns_type
  use TemperatureType       , only : temperature_type
  use PatchType             , only : patch
  use GridcellType          , only : grc     
  use SolarAbsorbedType     , only : solarabs_type
  use SurfaceAlbedoType     , only : surfalb_type
  use WaterDiagnosticBulkType        , only : waterdiagnosticbulk_type
  !use EDPhotosynthesisMod  , only : vcmaxc, jmaxc
  
  
  implicit none
  save
  
  !------------------------------------------------------------------------------
  ! PRIVATE MEMBER FUNCTIONS:
  public  :: LunaReadNML                                   !subroutine to read in the Luna namelist
  public  :: Update_Photosynthesis_Capacity                !subroutine to update the canopy nitrogen profile
  public  :: Acc24_Climate_LUNA                            !subroutine to accumulate 24 hr climates
  public  :: Acc240_Climate_LUNA                           !subroutine to accumulate 10 day climates
  public  :: Clear24_Climate_LUNA                          !subroutine to clear 24 hr climates
  private :: NitrogenAllocation                            !subroutine to update the Vcmax25 and Jmax25 at the leaf level
  private :: NUEref                                        !Calculate the Nitrogen use effieciency based on reference CO2 and leaf temperature
  private :: NUE                                           !Calculate the Nitrogen use effieciency based on current CO2 and leaf temperature
  private :: JmxTLeuning                                   !Calculate the temperature response for Jmax, based on Leunning 2002 Plant, Cell & Environment
  private :: JmxTKattge                                    !Calculate the temperature response for Jmax, based on Kattge and Knorr  2007
  private :: VcmxTLeuning                                  !Calculate the temperature response for Vcmax, based on Leunning 2002 Plant, Cell & Environment
  private :: VcmxTKattge                                   !Calculate the temperature response for Vcmax, based on Kattge and Knorr  2007
  private :: RespTBernacchi                                !Calculate the temperature response for respiration, following Bernacchi PCE 2001
  private :: Photosynthesis_luna                           !calculate the photosynthetic rate for nitrogen allocation
  private :: Quadratic                                     !Calculate the soultion using the quadratic formula

  !------------------------------------------------------------------------------ 
  !Constants  
  real(r8), parameter :: Cv = 1.2e-5_r8 * 3600.0           ! conversion factor from umol CO2 to g carbon
  real(r8), parameter :: Kc25 = 40.49_r8                   ! Mechanis constant of CO2 for rubisco(Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
  real(r8), parameter :: Ko25 = 27840_r8                   ! Mechanis constant of O2 for rubisco(Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
  real(r8), parameter :: Cp25 = 4.275_r8                   ! CO2 compensation point at 25C (Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
  real(r8), parameter :: Fc25 = 294.2_r8                   ! Fc25 = 6.22*47.3 #see Rogers (2014) Photosynthesis Research 
  real(r8), parameter :: Fj25 = 1257.0_r8                  ! Fj25 = 8.06*156 # #see COSTE 2005 and Xu et al 2012
  real(r8), parameter :: NUEr25 = 33.69_r8                 ! nitrogen use efficiency for respiration, see Xu et al 2012
  real(r8), parameter :: Cb = 1.78_r8                      ! nitrogen use effiency for choloraphyll for light capture, see Evans 1989  
  real(r8), parameter :: O2ref = 209460.0_r8                 ! ppm of O2 in the air
  real(r8), parameter :: CO2ref = 380.0_r8                   ! reference CO2 concentration for calculation of reference NUE. 
  real(r8), parameter :: forc_pbot_ref = 101325.0_r8       ! reference air pressure for calculation of reference NUE
  real(r8), parameter :: Q10Enz = 2.0_r8                   ! Q10 value for enzyme decay rate
  real(r8), parameter :: Jmaxb0 = 0.0311_r8                ! the baseline proportion of nitrogen allocated for electron transport (J)     
  real(r8)            :: Jmaxb1 = 0.1_r8                   ! the baseline proportion of nitrogen allocated for electron transport (J)    
  real(r8), parameter :: Wc2Wjb0 = 0.8054_r8               ! the baseline ratio of rubisco limited rate vs light limited photosynthetic rate (Wc:Wj) 
  real(r8), parameter :: relhExp = 6.0999_r8               ! electron transport parameters related to relative humidity
  real(r8), parameter :: Enzyme_turnover_daily = 0.1_r8    ! the daily turnover rate for photosynthetic enzyme at 25oC in view of ~7 days of half-life time for Rubisco (Suzuki et al. 2001)
  real(r8), parameter :: NMCp25 = 0.715_r8                 ! estimated by assuming 80% maintenance respiration is used for photosynthesis enzyme maintenance
  real(r8), parameter :: Trange1 = 5.0_r8                  ! lower temperature limit (oC) for nitrogen optimization  
  real(r8), parameter :: Trange2 = 42.0_r8                 ! upper temperature limit (oC) for nitrogen optimization
  real(r8), parameter :: SNC = 0.004_r8                    ! structural nitrogen concentration (g N g-1 dry mass carbon)
  real(r8), parameter :: mp = 9.0_r8                       ! slope of stomatal conductance; this is used to estimate model parameter, but may need to be updated from the physiology file, 
  real(r8), parameter :: PARLowLim = 200.0_r8              ! minimum photosynthetically active radiation for nitrogen optimization
  real(r8), parameter :: minrelh = 0.25_r8                 ! minimum relative humdity for nitrogen optimization

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------------
  
  contains

  !********************************************************************************************************************************************************************** 
  ! Read in LUNA namelist
  subroutine LunaReadNML( NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for LUNA
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    use abortutils     , only : endrun
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'lunaReadNML'
    character(len=*), parameter :: nmlname = 'luna'
    !-----------------------------------------------------------------------
    namelist /luna/ Jmaxb1

    ! Initialize options to default values, in case they are not specified in
    ! the namelist


    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=luna, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(__FILE__, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(__FILE__, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (Jmaxb1, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=luna)
       write(iulog,*) ' '
    end if

  end subroutine lunaReadNML

  !********************************************************************************************************************************************************************** 
  ! this subroutine updates the photosynthetic capacity as determined by Vcmax25 and Jmax25
  subroutine Update_Photosynthesis_Capacity(bounds, fn, filterp, &
    dayl_factor, atm2lnd_inst, temperature_inst, canopystate_inst, photosyns_inst, &
    surfalb_inst, solarabs_inst, waterdiagnosticbulk_inst, frictionvel_inst)
    !
    ! !DESCRIPTION:
    ! Calculates Nitrogen fractionation within the leaf, based on optimum calculated fractions in rubisco, cholorophyll, 
    ! Respiration and Storage. Based on Xu et al. 2012 and Ali et al 2015.In Review 
    
    !
    ! !REVISION HISTORY:
    ! version 1.0, by Chonggang Xu, Ashehad Ali and Rosie Fisher. July 14  2015.
    ! version 0.1, by Chonggang Xu, Ashehad Ali and Rosie Fisher. October 30 2014. 
    
    ! CALLED FROM:
    ! subroutine CanopyFluxes 
  
    ! !USES:
    use clm_time_manager      , only : get_step_size, is_end_curr_day
    use clm_varpar            , only : nlevsoi, mxpft
    use perf_mod              , only : t_startf, t_stopf
    use clm_varctl            , only : use_cn
    use quadraticMod          , only : quadratic
    use CNSharedParamsMod     , only : CNParamsShareInst
    use shr_infnan_mod, only : isnan => shr_infnan_isnan
        
    implicit none
    
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds                      
    integer                , intent(in)    :: fn                          ! size of pft filter
    integer                , intent(in)    :: filterp(fn)                 ! pft filter
    real(r8)               , intent(in)    :: dayl_factor( bounds%begp: ) ! scalar (0-1) for daylength

                       
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(temperature_type) , intent(inout) :: temperature_inst
    type(canopystate_type) , intent(inout) :: canopystate_inst
    type(photosyns_type)   , intent(inout) :: photosyns_inst
    type(surfalb_type)     , intent(in)    :: surfalb_inst
    type(solarabs_type)    , intent(inout) :: solarabs_inst
    type(waterdiagnosticbulk_type)  , intent(inout) :: waterdiagnosticbulk_inst
    type(frictionvel_type) , intent(inout) :: frictionvel_inst

    
    ! !LOCAL VARIABLES:
    !
    ! local pointers to implicit in variables
   
    integer   :: c,CL,f,g,iv,j,p,ps                                       ! indices
    integer   :: NCL_p                                                    ! number of canopy layers in patch 
    integer   :: ft                                                       ! plant functional type
    integer   :: z                                                        ! the index across leaf layers
    real (r8) :: PNstoreopt,PNlcopt,PNetopt,PNrespopt,PNcbopt             ! the optimal nitrogen allocations
    real (r8) :: PNstoreold,PNlcold,PNetold,PNrespold,PNcbold             ! the previous time step nitrogen allocations
    real (r8) :: delta_fn                                                 ! daily change in nitrogen investiment
    real (r8) :: relCLNCa                                                 ! the relative factor for LNCa due to canopy location and seasonal growth
    real (r8) :: relSLNCa                                                 ! the relative factor for LNCa due to seasonal growth
    real (r8) :: relRad                                                   ! the realtive radiation to the top of the canopy
    real (r8) :: FNCmtar                                                  ! target functional nitrogen content (g N/g leaf c)
    real (r8) :: LMA                                                      ! leaf mass per unit area     (g leaf c/m2 leaf)
    real (r8) :: PARTop                                                   ! photosynthetic active radiation on the top of canopy (umol/m2/s)
    real (r8) :: RadTop                                                   ! short-wave radiation on the top of canopy (w/m2)
    real (r8) :: TRad                                                     ! total short-wave radiation on the top of canopy (w/m2)
    real (r8) :: PARi10                                                   ! 10-day mean photosynthetic active radiation on in the canopy (umol/m2/s) 
    real (r8) :: PARimx10                                                 ! 10-day mean maximum photosynthetic active radiation on in the canopy (umol/m2/s)   
    real (r8) :: tleaf10                                                  ! 10-day mean  leaf temperature (oC)
    real (r8) :: tleafd10                                                 ! 10-day mean daytime leaf temperature (oC)
    real (r8) :: tleafn10                                                 ! 10-day mean nighttime leaf temperature (oC)
    real (r8) :: hourpd                                                   ! hours per day (hours)
    real (r8) :: CO2a10                                                   ! 10-day mean air co2 concentration (pa)
    real (r8) :: O2a10                                                    ! 10-day mean air o2 concentration (pa)
    real (r8) :: max_daily_pchg                                           ! maximum daily percentrage change  for nitrogen allocation
    real (r8) :: max_daily_decay                                          ! maximum daily decay  for nitrogen allocation  
    real (r8) :: radk                                                     ! light extintion factor
    real (r8) :: FNCa                                                     ! leaf functional nitrogen content (g/m2)
    real (r8) :: FNCa_z(1:nlevcan)                                        ! profile of leaf functional nitrogen content (g/m2) 
    real (r8) :: fnps                                                     ! fraction of light absorbed by non-photosynthetic pigments
    real (r8) :: radmax2mean                                              ! ratio of max radiation to mean
    real (r8) :: qabs                                                     ! PAR absorbed by PS II (umol photons/m**2/s)
    real (r8) :: EnzTurnoverTFactor                                       ! temperature adjust factor for enzyme decay
    real (r8) :: vcmax25                                                  ! Predicted vcmax25 from EDN model umol CO2/m**2/s
    real (r8) :: jmax25                                                   ! Predicted jmax25  from EDN model umol electrons/m**2/s
    real (r8) :: dtime                                                    ! stepsize in seconds
    real (r8) :: rb10v                                                    ! 10-day mean boundary layer resistance
    real (r8) :: relh10                                                   ! 10-day mean relative humidity (unitless)
    real (r8) :: tair10                                                   ! 10-day running mean of the 2m temperature (oC)
    real (r8) :: rabsorb                                                  ! ratio of absorbed raditation to the total incident radiation
    real (r8) :: tlaii                                                    ! total leaf area index for a certain canopy layer     
    real (r8) :: SNCa                                                     ! structural leaf nitrogen content (g N/m2 leaf)  
    real (r8) :: vcmx25_opt	                                          ! optimal Vc,max25 (umol CO2/m**2/s) 
    real (r8) :: jmx25_opt	                                          ! optimal Jmax25 (umol electron/m**2/s)        
    real (r8) :: chg                                                      ! change in Vcmax25  or Jmax25     
    real (r8) :: chg_constrn                                              ! constrained change in Vcmax25  or Jmax25
    logical   :: is_end_day                                               ! is end of current day
    !-------------------------------------------------------------------------------------------------------------------------------------------------       
    associate(                                                          &
    c3psn         => pftcon%c3psn                                     , & ! photosynthetic pathway: 0.  =  c4, 1.  =  c3
    slatop        => pftcon%slatop                                    , & ! specific leaf area at top of canopy, projected area basis [m^2/gC]
    leafcn        => pftcon%leafcn                                    , & ! leaf C:N (gC/gN)        
    forc_pbot10   => atm2lnd_inst%forc_pbot240_downscaled_patch       , & ! Input:  [real(r8) (:)   ] 10 day mean atmospheric pressure(Pa)  
    CO2_p240      => atm2lnd_inst%forc_pco2_240_patch                 , & ! Input:  [real(r8) (:)   ] 10-day mean CO2 partial pressure (Pa)
    O2_p240       => atm2lnd_inst%forc_po2_240_patch                  , & ! Input:  [real(r8) (:)   ] 10-day mean O2 partial pressure (Pa)    
    elai          => canopystate_inst%elai_patch                      , & ! Input:  [real(r8) (:)   ] one-sided leaf area index with burying by snow                        
    tlai          => canopystate_inst%tlai_patch                      , & ! Input:  [real(r8) (:)   ] one-sided leaf area index
    tlai_z        => surfalb_inst%tlai_z_patch                        , & ! Input:  [real(r8) (:)   ] one-sided leaf area index
    dayl          => grc%dayl                                         , & ! Input:  [real(r8) (:)   ] daylength (s)
    sabv          => solarabs_inst%sabv_patch                         , & ! Input:  [real(r8) (:)   ] patch solar radiation absorbed by vegetation (W/m**2)
    t_veg         => temperature_inst%t_veg_patch                     , & ! Input:  [real(r8) (:)   ] vegetation temperature (Kelvin) 
    rhol          => pftcon%rhol                                      , & ! Input:  [real(r8) (:)   ] leaf reflectance: 1=vis, 2=nir
    taul          => pftcon%taul                                      , & ! Input:  [real(r8) (:)   ] leaf transmittance: 1=vis, 2=nir
    par240d_z     => solarabs_inst%par240d_z_patch                    , & ! Input:  [real(r8) (:,:) ] 10-day running mean of daytime patch absorbed PAR for leaves in canopy layer (W/m**2) 
    par24d_z      => solarabs_inst%par24d_z_patch                     , & ! Input:  [real(r8) (:,:) ] daily accumulated absorbed PAR for leaves in canopy layer (W/m**2) 
    par240x_z     => solarabs_inst%par240x_z_patch                    , & ! Input:  [real(r8) (:,:) ] 10-day running mean of maximum patch absorbed PAR for leaves in canopy layer (W/m**2) 
    par24x_z      => solarabs_inst%par24x_z_patch                     , & ! Input:  [real(r8) (:,:) ] daily maximum of patch absorbed PAR for leaves in canopy layer (W/m**2) 
    nrad          => surfalb_inst%nrad_patch                          , & ! Input:  [integer  (:)   ] pft number of canopy layers, above snow for radiative transfer
    lnc           => photosyns_inst%lnca_patch                        , & ! Input:  [real(r8) (:)   ] top leaf layer leaf N concentration (gN leaf/m^2)
    t10           => temperature_inst%t_a10_patch                     , & ! Input:  [real(r8) (:)   ] 10-day running mean of the 2 m temperature (K)  
    t_veg_day     => temperature_inst%t_veg_day_patch                 , & ! Input:  [real(r8) (:)   ] daytime mean vegetation temperature (Kelvin)  
    t_veg_night   => temperature_inst%t_veg_night_patch               , & ! Input:  [real(r8) (:)   ] nighttime mean vegetation temperature (Kelvin)
    t_veg10_day   => temperature_inst%t_veg10_day_patch               , & ! Input:  [real(r8) (:)   ] 10-day mean daytime vegetation temperature (Kelvin)  
    t_veg10_night => temperature_inst%t_veg10_night_patch             , & ! Input:  [real(r8) (:)   ] 10-day mean nighttime vegetation temperature (Kelvin)
    rh10_p	  => waterdiagnosticbulk_inst%rh10_af_patch                    , & ! Input:  [real(r8) (:)   ] 10-day mean canopy air relative humidity at the pacth (unitless)
    rb10_p        => frictionvel_inst%rb10_patch                      , & ! Input:  [real(r8) (:)   ] 10-day mean boundary layer resistance at the pacth (s/m)
    gpp_day       => photosyns_inst%fpsn24_patch                      , & ! Input:  [real(r8) (:)   ] patch 24 hours mean gpp(umol CO2/m**2 ground/day) for canopy layer
    vcmx25_z      => photosyns_inst%vcmx25_z_patch                    , & ! Output: [real(r8) (:,:) ] patch leaf Vc,max25 (umol/m2 leaf/s) for canopy layer 
    jmx25_z       => photosyns_inst%jmx25_z_patch                     , & ! Output: [real(r8) (:,:) ] patch leaf Jmax25 (umol electron/m**2/s) for canopy layer
    pnlc_z        => photosyns_inst%pnlc_z_patch                      , & ! Output: [real(r8) (:,:) ] patch proportion of leaf nitrogen allocated for light capture for canopy layer 
    enzs_z        => photosyns_inst%enzs_z_patch                        & ! Output: [real(r8) (:,:) ] enzyme decay status 1.0-fully active; 0-all decayed during stress
    )  
    !----------------------------------------------------------------------------------------------------------------------------------------------------------
    !set timestep

    !Initialize enzyme decay Q10
    dtime        =  get_step_size()

    is_end_day   =  is_end_curr_day()
    fnps         =  0.15_r8
    call t_startf('LUNA')
    do f  =  1,fn
     p  =  filterp(f)
     ft =  patch%itype(p)     
     g  =  patch%gridcell(p)
     c  =  patch%column(p)
     !----------------------------------------------------
     !store the daily mean climate conditions
     if(t_veg_day(p).ne.spval) then    !check whether it is the first day 
         !------------------------------------------
         !get the climate driver    
         CO2a10 = CO2_p240(p)   
         O2a10  = O2_p240(p)
         hourpd = dayl(g) / 3600._r8             
         tleafd10 = t_veg10_day(p) - tfrz
         tleafn10 = t_veg10_night(p) - tfrz
         tleaf10  = (dayl(g)*tleafd10 +(86400._r8-dayl(g)) * tleafd10)/86400._r8 	     
         tair10 = t10(p)- tfrz
         relh10 = min(1.0_r8, rh10_p(p))  
	 rb10v = rb10_p(p)	     
         !--------------------------------------------------------------------
         !calculate the enzyme ternover rate
         EnzTurnoverTFactor = Q10Enz**(0.1_r8*(min(40.0_r8, tleaf10) - 25.0_r8))            
         max_daily_pchg = EnzTurnoverTFactor * Enzyme_turnover_daily
         !-----------------------------------------------------------------
         rabsorb = 1.0_r8-rhol(ft,1)-taul(ft,1)
         !Implemented the nitrogen allocation model
         if(tlai(p) > 0.0_r8 .and. lnc(p) > 0._r8)then   
                RadTop = par240d_z(p,1)/rabsorb
                PARTop = RadTop*4.6    !conversion from w/m2 to umol/m2/s. PAR is still in umol photons, not electrons. Also the par240d_z is only for radiation at visible range. Hence 4.6 not 2.3 multiplier. 
                !-------------------------------------------------------------
                !the nitrogen allocation model, may need to be feed from the parameter file in CLM
                if (nint(c3psn(ft)) == 1)then
                   if(gpp_day(p)>0.0 )then   !only optimize if there is growth and it is C3 plants
                      !-------------------------------------------------------------
                      do z = 1, nrad(p)
                         if(tlai_z(p,z)>0.0_r8)then
                            qabs  = par240d_z(p,z)/rabsorb
                            PARi10 =  qabs * 4.6_r8                 
                         else
                            PARi10  =  0.01_r8
                         endif
                         !-----------------------------------------------------------------------
                         relRad   = PARi10/PARTop
                         relCLNCa = 0.1802_r8*log(relRad)+1.0_r8 !see Ali et al 2015.
                         relCLNCa = max(0.2_r8,relCLNCa)
                         relCLNCa = min(1.0_r8,relCLNCa)
                         relSLNCa = 1.0_r8
                         !------------------------------------------------------------------
                         SNCa     =  1.0_r8/slatop(ft) * SNC
                         if(0.9_r8 * lnc(p)> SNCa)then
                           FNCa_z(z)= relCLNCa*(lnc(p)-SNCa)
                         else
                           FNCa_z(z)= relCLNCa*0.1_r8*lnc(p)
                         endif
                      enddo
                      
                      !----------------------------------------------------------------------
                      !nitrogen allocation model 
                      do z = 1 , nrad(p)
                         
                         !-------------------------------------------------------------------------------------------
                         !for different layers of leaves
                         FNCa     = FNCa_z(z)
                         if(FNCa>15.0_r8) then !boundary condition check for unrealistically high leaf nitrogen content
                             FNCa = 15.0_r8  
                             write(iulog, *) 'Warning: leaf nitrogen content become unrealistically high (>15.0 g N/m2 leaf) ', &
                                  'for patch=', p, 'z=', z, "pft=", ft
                         endif
                         radmax2mean = par240x_z(p,z) / par240d_z(p,z)
                         if(tlai_z(p,z)>0.0_r8)then
                            qabs   = par240d_z(p,z)/rabsorb
                            PARi10 =  qabs * 4.6_r8                 
                         else
                            PARi10 =  0.01_r8
                         endif                         
                         PARimx10  = PARi10*radmax2mean
                         !-----------------------------------------------------------------------------------------------------
 
                         !nitrogen allocastion model-start          
                         PNlcold   = PNlc_z(p,z)
                         PNetold   = 0.0_r8
                         PNrespold = 0.0_r8
                         PNcbold   = 0.0_r8                                     
                         call NitrogenAllocation(FNCa,forc_pbot10(p), relh10, CO2a10, O2a10, PARi10, PARimx10, rb10v, hourpd, &
                              tair10, tleafd10, tleafn10, &
                              Jmaxb0, Jmaxb1, Wc2Wjb0, relhExp, PNlcold, PNetold, PNrespold, &
                              PNcbold, PNstoreopt, PNlcopt, PNetopt, PNrespopt, PNcbopt)
                         vcmx25_opt= PNcbopt * FNCa * Fc25
                         jmx25_opt= PNetopt * FNCa * Fj25
                          
                         chg = vcmx25_opt-vcmx25_z(p, z)
                         chg_constrn = min(abs(chg),vcmx25_z(p, z)*max_daily_pchg)
                         vcmx25_z(p, z)  = vcmx25_z(p, z)+sign(1.0_r8,chg)*chg_constrn
                          
                         chg = jmx25_opt-jmx25_z(p, z)
                         chg_constrn = min(abs(chg),jmx25_z(p, z)*max_daily_pchg)
                         jmx25_z(p, z)  = jmx25_z(p, z)+sign(1.0_r8,chg)*chg_constrn 

                         PNlc_z(p, z)= PNlcopt

                         if(enzs_z(p,z)<1.0) then
                            enzs_z(p,z) = enzs_z(p,z)* (1.0_r8 + max_daily_pchg)
                         endif
                         !nitrogen allocastion model-end  

!DML turn off endrun and instead modify vcmx25_z(p,z) and jmx25_z(p,z) to a reasonable value
                         !-----------------------------------------------------------------------------------------------------  
                         if(isnan(vcmx25_z(p, z)))then
                             write(iulog, *) 'Error: Vc,mx25 is NaN for patch=', &
                                  p, 'z=', z, "pft=", ft
                             write(iulog, *) 'LUNA env:',FNCa,forc_pbot10(p), relh10, CO2a10, O2a10, PARi10, PARimx10, rb10v, &
                                  hourpd, tair10, tleafd10, tleafn10
                             call endrun(msg=errmsg(sourcefile, __LINE__))
                         endif
                         if(vcmx25_z(p, z)>1000._r8 .or. vcmx25_z(p, z)<0._r8)then
                             write(iulog, *) 'Warning: Vc,mx25 become unrealistic (>1000 or negative) for patch=', &
                                  p, 'z=', z, "pft=", ft
                             write(iulog, *) 'LUNA env:',vcmx25_z(p,z),FNCa,forc_pbot10(p), relh10, CO2a10, &
                                  O2a10, PARi10, PARimx10, rb10v, hourpd, tair10, tleafd10, tleafn10
                             vcmx25_z(p,z) = 50._r8
                         endif
                         if(isnan(jmx25_z(p, z)))then
                             write(iulog, *) 'Error: Jmx25 is NaN for patch=', &
                                  p, 'z=', z, "pft=", ft
                             write(iulog, *) 'LUNA env:', FNCa,forc_pbot10(p), relh10, CO2a10, O2a10, PARi10, PARimx10, rb10v, &
                                  hourpd, tair10, tleafd10, tleafn10
                             call endrun(msg=errmsg(sourcefile, __LINE__))
                         endif
                         if(jmx25_z(p, z)>2000._r8 .or.  jmx25_z(p, z)<0._r8)then
                             write(iulog, *) 'Warning: Jmx25 become unrealistic (>2000, or negative) for patch=', &
                                  p, 'z=', z, "pft=", ft
                             write(iulog, *) 'LUNA env:', jmx25_z(p,z),FNCa,forc_pbot10(p), relh10, CO2a10, &
                                  O2a10, PARi10, PARimx10, rb10v, hourpd, tair10, tleafd10, tleafn10
                             jmx25_z(p,z) = 85._r8
                         endif

                      enddo ! finished loop of leaf layers  
                   else !decay during drought or winter
                      max_daily_decay = min(0.5_r8, 0.1_r8 * max_daily_pchg)
                      !assume enzyme turnover under maintenance is 10
                      !times lower than enzyme change under growth
                      do z = 1 , nrad(p)
                         if(enzs_z(p,z)>0.5_r8) then
                            !decay is set at only 50% of original
                            !enzyme in view that plant will need to
                            !maintain their basic functionality
                          enzs_z(p,z) = enzs_z(p,z)* (1.0_r8 - max_daily_decay)
                          jmx25_z(p, z) = jmx25_z(p, z)* (1.0_r8 - max_daily_decay) 
                          vcmx25_z(p, z) = vcmx25_z(p, z)* (1.0_r8 - max_daily_decay) 
                         endif
                      end do              
                   endif !checking for growth                   
                endif !if not C3 plants                   
         else
            do z = 1 , nrad(p)
               jmx25_z(p, z) = 85._r8
               vcmx25_z(p, z) = 50._r8
            end do
         endif !checking for LAI and LNC
     endif !the first daycheck 
    end do !fn loop    
    call t_stopf('LUNA')
  end associate
  
end subroutine Update_Photosynthesis_Capacity



subroutine Acc240_Climate_LUNA(bounds, fn, filterp, oair, cair, &
    rb,rh, temperature_inst, photosyns_inst, &
    surfalb_inst, solarabs_inst, waterdiagnosticbulk_inst, frictionvel_inst)
    !
    ! !DESCRIPTION:
    ! Accumulate the 10-day running mean climates for LUNA model 
    
    !
    ! !REVISION HISTORY:
    ! version 1.0, by Chonggang Xu July 14  2015.
   
    ! CALLED FROM:
    ! subroutine CanopyFluxes 
  
    ! !USES:
    use clm_time_manager      , only : get_step_size, is_end_curr_day
    implicit none
    
      ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds                      
    integer                , intent(in)    :: fn                          ! size of pft filter
    integer                , intent(in)    :: filterp(fn)                 ! pft filter
    real(r8)               , intent(in)    :: oair( bounds%begp: )        ! Atmospheric O2 partial pressure (Pa)
    real(r8)               , intent(in)    :: cair( bounds%begp: )        ! Atmospheric CO2 partial pressure (Pa)
    real(r8)               , intent(in)    :: rb( bounds%begp: )          ! boundary layer resistance (s/m)
    real(r8)               , intent(in)    :: rh( bounds%begp: )          ! canopy air relative humidity
                       
    type(temperature_type) , intent(inout) :: temperature_inst
    type(photosyns_type)   , intent(inout) :: photosyns_inst
    type(surfalb_type)     , intent(in)    :: surfalb_inst
    type(solarabs_type)    , intent(inout) :: solarabs_inst
    type(waterdiagnosticbulk_type)  , intent(inout) :: waterdiagnosticbulk_inst
    type(frictionvel_type) , intent(inout) :: frictionvel_inst
    
    ! !LOCAL VARIABLES:
    !
    ! local pointers to implicit in variables
   
    integer   :: c,f,g,iv,j,p                                             ! indices
    integer   :: ft                                                       ! plant functional type
    integer   :: z                                                        ! the index across leaf layers
    real (r8) :: dtime                                                    ! stepsize in seconds
    real (r8) :: TRad                                                     ! total short-wave radiation on the top of canopy (w/m2)
    real (r8) :: tlaii                                                    ! total leaf area index for a certain canopy layer     
    real (r8) :: t_veg_dayi                                               ! daytime mean vegetation temperature (Kelvin)
    real (r8) :: t_veg_nighti                                             ! nighttime mean vegetation temperature (Kelvin)
    real (r8) :: par24d_z_i(1:nlevcan)                                    ! daytime mean radiation (w/m**2)             
    logical   :: is_end_day                                               ! is end of current day
    !-------------------------------------------------------------------------------------------------------------------------------------------------       
    associate(                                                          &
    par24d_z      => solarabs_inst%par24d_z_patch                     , & ! Input:  [real(r8) (:,:) ] daily accumulated absorbed PAR for leaves in canopy layer (W/m**2) 
    par24x_z      => solarabs_inst%par24x_z_patch                     , & ! Input:  [real(r8) (:,:) ] daily maximum of patch absorbed PAR for leaves in canopy layer (W/m**2) 
    nrad          => surfalb_inst%nrad_patch                          , & ! Input:  [integer  (:)   ] pft number of canopy layers, above snow for radiative transfer
    t_veg_day     => temperature_inst%t_veg_day_patch                 , & ! Input:  [real(r8) (:)   ] daytime accumulative vegetation temperature (Kelvin*nsteps)  
    t_veg_night   => temperature_inst%t_veg_night_patch               , & ! Input:  [real(r8) (:)   ] nighttime accumulative vegetation temperature (Kelvin*nsteps)
    nnightsteps   => temperature_inst%nnightsteps_patch               , & ! Input:  [integer  (:)   ] number of nighttime steps in 24 hours from mid-night, LUNA specific
    ndaysteps     => temperature_inst%ndaysteps_patch                 , & ! Input:  [integer  (:)   ] number of daytime steps in 24 hours from mid-night, LUNA specific
    t_veg10_day   => temperature_inst%t_veg10_day_patch               , & ! Output: [real(r8) (:)   ] 10-day mean vegetation temperature (Kelvin)  
    t_veg10_night => temperature_inst%t_veg10_night_patch             , & ! Output: [real(r8) (:)   ] 10-day mean vegetation temperature (Kelvin)
    rh10_p	  => waterdiagnosticbulk_inst%rh10_af_patch                    , & ! Output: [real(r8) (:)   ] 10-day mean canopy air relative humidity at the pacth (s/m)
    rb10_p        => frictionvel_inst%rb10_patch                      , & ! Output: [real(r8) (:)   ] 10-day mean boundary layer resistance at the pacth (s/m)
    par240d_z     => solarabs_inst%par240d_z_patch                    , & ! Output:  [real(r8) (:,:) ] 10-day running mean of daytime patch absorbed PAR for leaves in canopy layer (W/m**2) 
    par240x_z     => solarabs_inst%par240x_z_patch                      & ! Output:  [real(r8) (:,:) ] 10-day running mean of maximum patch absorbed PAR for leaves in canopy layer (W/m**2)

    )  
    !----------------------------------------------------------------------------------------------------------------------------------------------------------
    !set timestep

    !Initialize enzyme decay Q10
    dtime        =  get_step_size()
    is_end_day   =  is_end_curr_day()
    do f  =  1,fn
      p  =  filterp(f)
      ft =  patch%itype(p)
      g  =  patch%gridcell(p)
      c  =  patch%column(p)
      if(t_veg_day(p).ne.spval) then    !check whether it is the first day            
             !---------------------------------------------------------
             !calculate the 10 day running mean radiations
             if(ndaysteps(p)>0.0) then
                 par24d_z_i=par24d_z(p,:)/(dtime * ndaysteps(p))
             else
                 par24d_z_i = 0._r8
             endif
             if(par240d_z(p,1).eq. spval)then  !first day set as the same of first day environmental conditions
                par240x_z(p,:)= par24x_z(p,:)
                par240d_z(p,:)= par24d_z_i
             else
                par240x_z(p,:)= 0.9_r8 * par240x_z(p,:) + 0.1_r8 * par24x_z(p,:)
                par240d_z(p,:)= 0.9_r8 * par240d_z(p,:) + 0.1_r8 * par24d_z_i
             endif
             !-------------------------------------------------------
             !calculate the 10 day running mean daytime temperature
             if(ndaysteps(p)>0.0)then
                t_veg_dayi    =  t_veg_day(p)   / ndaysteps(p)
             else
                t_veg_dayi    =  t_veg_night(p) / nnightsteps(p)
             endif
             if(t_veg10_day(p).eq. spval)then
                  t_veg10_day(p)  =  t_veg_dayi
             endif
             t_veg10_day(p)  =  0.9_r8 * t_veg10_day(p)+ 0.1_r8 * t_veg_dayi
             !-------------------------------------------------------
             !calculate the 10 day running mean nighttime temperature
             if(nnightsteps(p)>0)then
               t_veg_nighti  =  t_veg_night(p) / nnightsteps(p)
             else
               t_veg_nighti  =  t_veg_day(p)   / ndaysteps(p)
             endif
             if(t_veg10_night(p).eq. spval)then
                  t_veg10_night(p)  =  t_veg_nighti
             endif
             t_veg10_night(p)  =  0.9_r8 * t_veg10_night(p) + 0.1_r8 * t_veg_nighti
             !--------------------------------------------------------------------
             if(rh10_p(p).eq. spval)then
                rh10_p(p)  =  rh(p)
             endif
             rh10_p(p) = 0.9_r8 * rh10_p(p) + 0.1_r8 * min(1.0_r8, rh(p))

             if(rb10_p(p).eq. spval)then
                rb10_p(p)  =  rb(p)
             endif
             rb10_p(p) = 0.9_r8 * rb10_p(p) + 0.1_r8 * rb(p)
       endif !the first day check  
    end do !fn loop    
  end associate
end subroutine Acc240_Climate_LUNA


subroutine Acc24_Climate_LUNA(bounds, fn, filterp, canopystate_inst, photosyns_inst, &
    surfalb_inst, solarabs_inst,temperature_inst)
    !
    ! !DESCRIPTION:
    ! Accumulate the 24 hr climates for LUNA model 
    
    !
    ! !REVISION HISTORY:
    ! version 1.0, by Chonggang Xu July 14  2015.
   
    ! CALLED FROM:
    ! subroutine CanopyFluxes 
  
    ! !USES:
    use clm_time_manager      , only : get_step_size
    implicit none
    
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds                      
    integer                , intent(in)    :: fn                          ! size of pft filter
    integer                , intent(in)    :: filterp(fn)                 ! pft filter
                

    type(canopystate_type) , intent(inout) :: canopystate_inst
    type(photosyns_type)   , intent(inout) :: photosyns_inst
    type(surfalb_type)     , intent(in)    :: surfalb_inst
    type(solarabs_type)    , intent(inout) :: solarabs_inst
    type(temperature_type) , intent(inout) :: temperature_inst
    
    ! !LOCAL VARIABLES:
    !
    ! local pointers to implicit in variables
   
    integer   :: c,f,g,iv,j,p                                             ! indices
    integer   :: ft                                                       ! plant functional type
    integer   :: z                                                        ! the index across leaf layers
    real (r8) :: dtime                                                    ! stepsize in seconds
    real (r8) :: TRad                                                     ! total short-wave radiation on the top of canopy (w/m2)
    real (r8) :: tlaii                                                    ! total leaf area index for a certain canopy layer     

    !-------------------------------------------------------------------------------------------------------------------------------------------------       
    associate(                                                          &
    sabv          => solarabs_inst%sabv_patch                         , & ! Input:  [real(r8) (:)   ] patch solar radiation absorbed by vegetation (W/m**2)
    t_veg         => temperature_inst%t_veg_patch                     , & ! Input:  [real(r8) (:)   ] vegetation temperature (Kelvin) 
    par_sun_z     => solarabs_inst%parsun_z_patch                     , & ! Input:  [real(r8) (:,:) ] par absorbed per unit lai for sunlit canopy layer (w/m**2)
    par_sha_z     => solarabs_inst%parsha_z_patch                     , & ! Input:  [real(r8) (:,:) ] par absorbed per unit lai for shaded canopy layer (w/m**2)
    lai_sun_z     => canopystate_inst%laisun_z_patch                  , & ! Input:  [real(r8) (:,:) ] leaf area index for sunlit canopy layer
    lai_sha_z     => canopystate_inst%laisha_z_patch                  , & ! Input:  [real(r8) (:,:) ] leaf area index for canopy shaded layer
    par24d_z      => solarabs_inst%par24d_z_patch                     , & ! Input:  [real(r8) (:,:) ] daily accumulated absorbed PAR for leaves in canopy layer (W/m**2) 
    par24x_z      => solarabs_inst%par24x_z_patch                     , & ! Input:  [real(r8) (:,:) ] daily maximum of patch absorbed PAR for leaves in canopy layer (W/m**2) 
    nrad          => surfalb_inst%nrad_patch                          , & ! Input:  [integer  (:)   ] pft number of canopy layers, above snow for radiative transfer
    gpp           => photosyns_inst%fpsn_patch                        , & ! Input:  [real(r8) (:)   ] patch instaneous gpp (umol CO2/m**2 ground/s) for canopy layer
    gpp_day       => photosyns_inst%fpsn24_patch                      , & ! Output: [real(r8) (:)   ] patch 24 hours acculative gpp(umol CO2/m**2 ground/day) from mid-night
    t_veg_day     => temperature_inst%t_veg_day_patch                 , & ! Output: [real(r8) (:)   ] daytime mean vegetation temperature (Kelvin)  
    t_veg_night   => temperature_inst%t_veg_night_patch               , & ! Output: [real(r8) (:)   ] nighttime mean vegetation temperature (Kelvin)
    nnightsteps   => temperature_inst%nnightsteps_patch               , & ! Output: [integer  (:)   ] number of nighttime steps in 24 hours from mid-night, LUNA specific
    ndaysteps     => temperature_inst%ndaysteps_patch                   & ! Output: [integer  (:)   ] number of daytime steps in 24 hours from mid-night, LUNA specific
    )  
    !----------------------------------------------------------------------------------------------------------------------------------------------------------
    !set timestep

    !Initialize enzyme decay Q10
    dtime        =  get_step_size()
    do f  =  1,fn
      p  =  filterp(f)
      ft =  patch%itype(p)
      g  =  patch%gridcell(p)
      c  =  patch%column(p)
      !----------------------------------------------------
      !store the daily mean climate conditions
      if(t_veg_day(p).ne.spval) then    !check whether it is the first day 
         if(sabv(p)>0)then
             t_veg_day(p)   = t_veg_day(p)   + t_veg(p)    
             ndaysteps(p)   = ndaysteps(p)   + 1          
         else
             t_veg_night(p) = t_veg_night(p) + t_veg(p) 
             nnightsteps(p) = nnightsteps(p) + 1 
         endif
         do z = 1, nrad(p)
          !average of sunlit and shaded leaves          
          tlaii = lai_sun_z(p,z) + lai_sha_z(p,z)          
          if(tlaii > 0._r8)then
              TRad = (par_sun_z(p,z)*lai_sun_z(p,z)+par_sha_z(p,z)*lai_sha_z(p,z))/tlaii
              TRad = par_sun_z(p,z) !RF & GBB. Make LUNA predict sunlit fraction N fractionation, scale in PhotosynthesisMod. 
              par24d_z(p,z)=   par24d_z(p,z)+ dtime * TRad 
             if(TRad > par24x_z(p,z))then
                par24x_z(p,z) = TRad
             endif
           endif 
         enddo
         gpp_day(p) = gpp_day(p) + dtime * gpp(p) 
      endif !first day check
    end do !fn loop    
  end associate
end subroutine Acc24_Climate_LUNA




subroutine Clear24_Climate_LUNA(bounds, fn, filterp, canopystate_inst, photosyns_inst, &
    surfalb_inst, solarabs_inst,temperature_inst)
    !
    ! !DESCRIPTION:
    ! Zero out the 24 hr climates for LUNA model 
    
    !
    ! !REVISION HISTORY:
    ! version 1.0, by Chonggang Xu July 14  2015.
   
    ! CALLED FROM:
    ! subroutine CanopyFluxes 
  
    ! !USES:
    use clm_time_manager      , only : get_step_size, is_end_curr_day
    implicit none
    
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds                      
    integer                , intent(in)    :: fn                          ! size of pft filter
    integer                , intent(in)    :: filterp(fn)                 ! pft filter
                

    type(canopystate_type) , intent(inout) :: canopystate_inst
    type(photosyns_type)   , intent(inout) :: photosyns_inst
    type(surfalb_type)     , intent(in)    :: surfalb_inst
    type(solarabs_type)    , intent(inout) :: solarabs_inst
    type(temperature_type) , intent(inout) :: temperature_inst
    
    ! !LOCAL VARIABLES:
    !
    ! local pointers to implicit in variables
   
    integer   :: c,f,g,iv,j,p                                             ! indices
    integer   :: ft                                                       ! plant functional type
    integer   :: z                                                        ! the index across leaf layers
    real (r8) :: dtime                                                    ! stepsize in seconds
    logical   :: is_end_day                                               ! is end of current day
    !-------------------------------------------------------------------------------------------------------------------------------------------------       
    associate(                                                          &
    par24d_z      => solarabs_inst%par24d_z_patch                     , & ! Output:  [real(r8) (:,:) ] daily accumulated absorbed PAR for leaves in canopy layer (W/m**2) 
    par24x_z      => solarabs_inst%par24x_z_patch                     , & ! Output:  [real(r8) (:,:) ] daily maximum of patch absorbed PAR for leaves in canopy layer (W/m**2) 
    gpp_day       => photosyns_inst%fpsn24_patch                      , & ! Output: [real(r8) (:)   ] patch 24 hours acculative gpp(umol CO2/m**2 ground/day) from mid-night
    t_veg_day     => temperature_inst%t_veg_day_patch                 , & ! Output: [real(r8) (:)   ] daytime mean vegetation temperature (Kelvin)  
    t_veg_night   => temperature_inst%t_veg_night_patch               , & ! Output: [real(r8) (:)   ] nighttime mean vegetation temperature (Kelvin)
    nnightsteps   => temperature_inst%nnightsteps_patch               , & ! Output: [integer  (:)   ] number of nighttime steps in 24 hours from mid-night, LUNA specific
    ndaysteps     => temperature_inst%ndaysteps_patch                   & ! Output: [integer  (:)   ] number of daytime steps in 24 hours from mid-night, LUNA specific
    )  
    !----------------------------------------------------------------------------------------------------------------------------------------------------------
    !set timestep

    !Initialize enzyme decay Q10
    dtime        =  get_step_size()
    is_end_day   =  is_end_curr_day()
    do f  =  1,fn
      p  =  filterp(f)
      ft =  patch%itype(p)
      g  =  patch%gridcell(p)
      c  =  patch%column(p)
      !------------------------------------------------------------------------------
      !clear out the daily state variables at the begining of simulations 
      t_veg_day(p)                      =  0.0_r8
      t_veg_night(p)                    =  0.0_r8
      par24d_z(p,:)                     =  0.0_r8
      par24x_z(p,:)                     =  0.0_r8
      gpp_day(p)                        =  0.0_r8 
      nnightsteps(p)                    =  0.0_r8
      ndaysteps(p)                      =  0.0_r8

    end do !fn loop    
  end associate
end subroutine Clear24_Climate_LUNA


!************************************************************************************************************************************************
!Use the LUNA model to calculate the Nitrogen partioning 
subroutine NitrogenAllocation(FNCa,forc_pbot10, relh10, CO2a10,O2a10, PARi10,PARimx10,rb10, hourpd, tair10, tleafd10, tleafn10, &
     Jmaxb0, Jmaxb1, Wc2Wjb0, relhExp,&
     PNlcold, PNetold, PNrespold, PNcbold, &
     PNstoreopt, PNlcopt, PNetopt, PNrespopt, PNcbopt)
  implicit none
  real(r8), intent (in) :: FNCa                       !Area based functional nitrogen content (g N/m2 leaf)
  real(r8), intent (in) :: forc_pbot10                !10-day mean air pressure (Pa)
  real(r8), intent (in) :: relh10                     !10-day mean relative humidity (unitless)
  real(r8), intent (in) :: CO2a10                     !10-day meanCO2 concentration in the air (Pa)
  real(r8), intent (in) :: O2a10                      !10-day mean O2 concentration in the air (Pa)
  real(r8), intent (in) :: PARi10                     !10-day mean photosynthetic active radiation on in a canopy (umol/m2/s)
  real(r8), intent (in) :: PARimx10                   !10-day mean 24hr maximum photosynthetic active radiation on in a canopy (umol/m2/s)
  real(r8), intent (in) :: rb10                       !10-day mean boundary layer resistance
  real(r8), intent (in) :: hourpd                     !hours of light in a the day (hrs)
  real(r8), intent (in) :: tair10                     !10-day running mean of the 2m temperature (oC)
  real(r8), intent (in) :: tleafd10                   !10-day running mean of daytime leaf temperature (oC) 
  real(r8), intent (in) :: tleafn10                   !10-day running mean of nighttime leaf temperature (oC) 
  real(r8), intent (in) :: Jmaxb0                     !baseline proportion of nitrogen allocated for electron transport rate (unitless)
  real(r8), intent (in) :: Jmaxb1                     !coefficient determining the response of electron transport rate to light availability (unitless) 
  real(r8), intent (in) :: Wc2Wjb0                    !the baseline ratio of rubisco-limited rate vs light-limited photosynthetic rate (Wc:Wj)
  real(r8), intent (in) :: relhExp                    !specifies the impact of relative humidity on electron transport rate (unitless)
  real(r8), intent (in) :: PNlcold                    !old value of the proportion of nitrogen allocated to light capture (unitless)
  real(r8), intent (in) :: PNetold                    !old value of the proportion of nitrogen allocated to electron transport (unitless)
  real(r8), intent (in) :: PNrespold                  !old value of the proportion of nitrogen allocated to respiration (unitless)
  real(r8), intent (in) :: PNcbold                    !old value of the proportion of nitrogen allocated to carboxylation (unitless)  
  real(r8), intent (out):: PNstoreopt                 !optimal proportion of nitrogen for storage 
  real(r8), intent (out):: PNlcopt                    !optimal proportion of nitrogen for light capture 
  real(r8), intent (out):: PNetopt                    !optimal proportion of nitrogen for electron transport 
  real(r8), intent (out):: PNrespopt                  !optimal proportion of nitrogen for respiration 
  real(r8), intent (out):: PNcbopt                    !optial proportion of nitrogen for carboxyaltion  
 
  !-------------------------------------------------------------------------------------------------------------------------------
  !intermediate variables
  real(r8) :: Carboncost1                             !absolute amount of carbon cost associated with maintenance respiration due to deccrease in light capture nitrogen(g dry mass per day) 
  real(r8) :: Carboncost2                             !absolute amount of carbon cost associated with maintenance respiration due to increase in light capture nitrogen(g dry mass per day) 
  real(r8) :: Carbongain1                             !absolute amount of carbon gain associated with maintenance respiration due to deccrease in light capture nitrogen(g dry mass per day) 
  real(r8) :: Carbongain2                             !absolute amount of carbon gain associated with maintenance respiration due to increase in light capture nitrogen(g dry mass per day) 
  real(r8) :: Fc                                      !the temperature adjustment factor for Vcmax 
  real(r8) :: Fj                                      !the temperature adjustment factor for Jmax 
  real(r8) :: PNlc                                    !the current nitrogen allocation proportion for light capture
  real(r8) :: Jmax                                    !the maximum electron transport rate (umol/m2/s) 
  real(r8) :: JmaxCoef                                !coefficient determining the response of electron transport rate to light availability (unitless) and humidity
  real(r8) :: Jmaxb0act                               !base value of Jmax (umol/m2/s) 
  real(r8) :: JmaxL                                   !the electron transport rate with maximum daily radiation (umol/m2/s)  
  real(r8) :: JmeanL                                  !the electron transport rate with mean radiation (umol/m2/s) 
  real(r8) :: Nstore                                  !absolute amount of nitrogen allocated to storage (gN/m2 leaf)
  real(r8) :: Nresp                                   !absolute amount of nitrogen allocated to respiration (gN/m2 leaf) 
  real(r8) :: Nlc                                     !absolute amount of nitrogen allocated to light capture (gN/m2 leaf) 
  real(r8) :: Net                                     !absolute amount of nitrogen allocated to electron transport (gN/m2 leaf) 
  real(r8) :: Ncb                                     !absolute amount of nitrogen allocated to carboxylation (gN/m2 leaf) 
  real(r8) :: Nresp1                                  !absolute amount of nitrogen allocated to respiration due to increase in light capture nitrogen(gN/m2 leaf)  
  real(r8) :: Nlc1                                    !absolute amount of nitrogen allocated to light capture due to increase in light capture nitrogen(gN/m2 leaf) 
  real(r8) :: Net1                                    !absolute amount of nitrogen allocated to electron transport due to increase in light capture nitrogen(gN/m2 leaf)
  real(r8) :: Ncb1                                    !absolute amount of nitrogen allocated to carboyxlation due to increase in light capture nitrogen(gN/m2 leaf) 
  real(r8) :: Nresp2                                  !absolute amount of nitrogen allocated to respiration due to decrease in light capture nitrogen(gN/m2 leaf) 
  real(r8) :: Nlc2                                    !absolute amount of nitrogen allocated to light capture due to decrease in light capture nitrogen(gN/m2 leaf) 
  real(r8) :: Net2                                    !absolute amount of nitrogen allocated to electron transport due to decrease in light capture nitrogen(gN/m2 leaf) 
  real(r8) :: Ncb2                                    !absolute amount of nitrogen allocated to carboxylation due to increase in light capture nitrogen(gN/m2 leaf) 
  real(r8) :: PSN                                     !g carbon photosynthesized per day per unit(m2) of leaf
  real(r8) :: RESP                                    !g carbon respired per day per unit(m2) of leaf due to increase in light capture nitrogen(gN/m2 leaf) 
  real(r8) :: PSN1                                    !g carbon photosynthesized per day per unit(m2) of leaf due to increase in light capture nitrogen(gN/m2 leaf) 
  real(r8) :: RESP1                                   !g carbon respired per day per unit(m2) of leaf due to decrease in light capture nitrogen(gN/m2 leaf) 
  real(r8) :: PSN2                                    !g carbon photosynthesized per day per unit(m2) of leaf due to decrease in light capture nitrogen(gN/m2 leaf) 
  real(r8) :: RESP2                                   !g carbon respired per day per unit(m2) of leaf
  real(r8) :: Npsntarget                              !absolute amount of target nitrogen for photosynthesis(gN/m2 leaf) 
  real(r8) :: Npsntarget1                             !absolute amount of target nitrogen for photosynthesis due to increase in light capture nitrogen(gN/m2 leaf) 
  real(r8) :: Npsntarget2                             !absolute amount of target nitrogen for photosynthesis due to decrease in light capture nitrogen(gN/m2 leaf) 
  real(r8) :: NUEj                                    !nitrogen use efficiency for electron transport under current environmental conditions 
  real(r8) :: NUEc                                    !nitrogen use efficiency for carboxylation under current environmental conditions  
  real(r8) :: NUEjref                                 !nitrogen use efficiency for electron transport under reference environmental conditions (25oC and 385ppm Co2) 
  real(r8) :: NUEcref                                 !nitrogen use efficiency for carboxylation under reference environmental conditions (25oC and 385ppm Co2) 
  real(r8) :: NUEr                                    !nitrogen use efficiency for respiration 
  real(r8) :: PARi10c                                 !10-day mean constrained photosynthetic active radiation on in a canopy (umol/m2/s)
  real(r8) :: PARimx10c                               !10-day mean constrained 24hr maximum photosynthetic active radiation on in a canopy (umol/m2/s)
  real(r8) :: Kj2Kcref                                !the ratio of rubisco-limited photosynthetic rate (Wc) to light limited photosynthetic rate (Wj)
  real(r8) :: PNlcoldi                                !old value of the proportion of nitrogen allocated to light capture (unitless) 
  real(r8) :: Kj2Kc                                   !the ratio of Wc to Wj under changed conditions 
  real(r8) :: Kc                                      !conversion factors for Vc,max to Wc 
  real(r8) :: Kj                                      !conversion factor for electron transport rate to Wj 
  real(r8) :: theta                                   !efficiency of light energy conversion (unitless) 
  real(r8) :: chg_per_step                            !the nitrogen change per interation
  real(r8) :: Vcmaxnight                              !Vcmax during night (umol/m2/s)
  real(r8) :: ci                                      !inter-cellular CO2 concentration (Pa)
  real(r8) :: theta_cj                                !interpolation coefficient
  real(r8) :: tleafd10c                               !10-day mean daytime leaf temperature, contrained for physiological range (oC)
  real(r8) :: tleafn10c                               !10-day mean leaf temperature for night, constrained for physiological range (oC)
  real(r8) :: Vcmax                                   !the maximum carboxyaltion rate (umol/m2/s) 
  integer  :: KcKjFlag                                !flag to indicate whether to update the Kc and Kj using the photosynthesis subroutine; 0--Kc and Kj need to be calculated; 1--Kc and Kj is prescribed.
  integer  :: jj                                      !index record fo the number of iterations
  integer  :: increase_flag                           !whether to increase or decrease

  call NUEref(NUEjref, NUEcref, Kj2Kcref)
  theta_cj = 0.95_r8
  Nlc = PNlcold * FNCa                                !proportion of light capturing nitrogen in functional nitrogen
  Net = PNetold * FNCa                                !proportion of light harvesting (electron transport) nitrogen in functional nitrogen
  Nresp = PNrespold * FNCa                            !proportion of respirational nitrogen in functional nitrogen
  Ncb = PNcbold * FNCa                                !proportion of carboxylation nitrogen in functional nitrogen
  if (Nlc > FNCa * 0.5_r8) Nlc = 0.5_r8 * FNCa
  chg_per_step = 0.02* FNCa
  PNlc = PNlcold
  PNlcoldi = PNlcold  - 0.001_r8
  PARi10c = max(PARLowLim, PARi10)
  PARimx10c = max(PARLowLim, PARimx10)
  increase_flag = 0
  jj = 1
  tleafd10c = min(max(tleafd10, Trange1), Trange2)    !constrain the physiological range
  tleafn10c = min(max(tleafn10, Trange1), Trange2)    !constrain the physiological range
  ci = 0.7_r8 * CO2a10 
  JmaxCoef = Jmaxb1 * ((hourpd / 12.0_r8)**2.0_r8) * (1.0_r8 - exp(-relhExp * max(relh10 - minrelh, 0.0_r8) / &
       (1.0_r8 - minrelh)))
  do while (PNlcoldi .NE. PNlc .and. jj < 100)      
     Fc = VcmxTKattge(tair10, tleafd10c) * Fc25
     Fj = JmxTKattge(tair10, tleafd10c) * Fj25
     NUEr = Cv * NUEr25 * (RespTBernacchi(tleafd10c) * hourpd + RespTBernacchi(tleafn10c) * (24.0_r8 - hourpd)) !nitrogen use efficiency for respiration (g biomass/m2/day/g N)
     !****************************************************
     !Nitrogen Allocation Scheme: store the initial value
     !****************************************************
     KcKjFlag = 0
     call NUE(O2a10, ci, tair10, tleafd10c, NUEj, NUEc, Kj2Kc)
     call Nitrogen_investments (KcKjFlag,FNCa, Nlc, forc_pbot10, relh10, CO2a10,O2a10, PARi10c, PARimx10c,rb10, hourpd, tair10, &
          tleafd10c,tleafn10c, &
          Kj2Kc, Wc2Wjb0, JmaxCoef, Fc,Fj, NUEc, NUEj, NUEcref, NUEjref, NUEr, Kc, Kj, ci, &
          Vcmax, Jmax,JmeanL,JmaxL, Net, Ncb, Nresp, PSN, RESP)

     Npsntarget = Nlc + Ncb + Net                                                         !target nitrogen allocated to photosynthesis, which may be lower or higher than Npsn_avail
     PNlcoldi = Nlc / FNCa
     Nstore = FNCa - Npsntarget - Nresp
     !------------------------------------------------------------------------------------
     !test the increase of light capture nitrogen
     if (Nstore > 0.0_r8 .and.(increase_flag .eq. 1 .or. jj .eq. 1)) then
        Nlc2 = Nlc + chg_per_step
        if (Nlc2 / FNCa > 0.95_r8) Nlc2 = 0.95_r8 * FNCa
        KcKjFlag = 1
        call Nitrogen_investments (KcKjFlag,FNCa, Nlc2, forc_pbot10, relh10, CO2a10,O2a10, PARi10c, PARimx10c,rb10, hourpd, &
             tair10, tleafd10c,tleafn10c, &
             Kj2Kc, Wc2Wjb0, JmaxCoef, Fc,Fj, NUEc, NUEj, NUEcref, NUEjref,NUEr, Kc, Kj, ci, &
             Vcmax, Jmax,JmeanL,JmaxL, Net2, Ncb2, Nresp2, PSN2, RESP2)

        Npsntarget2 = Nlc2 + Ncb2 + Net2
        !update the nitrogen change
        Carboncost2 = (Npsntarget2 - Npsntarget) * NMCp25 * Cv * (RespTBernacchi(tleafd10c) * hourpd + &
             RespTBernacchi(tleafn10c) * (24.0_r8  - hourpd))
        Carbongain2 =  PSN2 - PSN
        if(Carbongain2 > Carboncost2 .and. (Npsntarget2 + Nresp2 < 0.95_r8 * FNCa))then
           Nlc = Nlc2
           Net = Net2
           Ncb = Ncb2
           Nstore = FNCa - Npsntarget2 - Nresp2 
           if (jj == 1) increase_flag = 1
        end if
     end if
     !------------------------------------------------------------------------------------
     !test the decrease of light capture nitrogen
     if (increase_flag == 0) then  
        if (Nstore < 0.0_r8) then
            Nlc1 = Nlc * 0.8_r8 !bigger step of decrease if it is negative            
        else
            Nlc1 = Nlc - chg_per_step
	end if
        if (Nlc1 < 0.05_r8) Nlc1 = 0.05_r8
        KcKjFlag = 1
        call Nitrogen_investments (KcKjFlag,FNCa, Nlc1,forc_pbot10, relh10, CO2a10,O2a10, PARi10c, PARimx10c,rb10, hourpd, &
             tair10, tleafd10c,tleafn10c, &
             Kj2Kc, Wc2Wjb0, JmaxCoef, Fc,Fj, NUEc, NUEj, NUEcref, NUEjref,NUEr, Kc, Kj, ci,&
             Vcmax, Jmax,JmeanL,JmaxL, Net1, Ncb1, Nresp1, PSN1, RESP1)
        Npsntarget1 = Nlc1 + Ncb1 + Net1
        Carboncost1 = (Npsntarget - Npsntarget1) * NMCp25 * Cv * (RespTBernacchi(tleafd10c) * hourpd + &
             RespTBernacchi(tleafn10c) * (24.0_r8  - hourpd))
        Carbongain1 =  PSN - PSN1
        if((Carbongain1 < Carboncost1 .and. Nlc1 > 0.05_r8) .or. (Npsntarget + Nresp) > 0.95_r8 * FNCa)then
          Nlc = Nlc1 
          Net = Net1   
          Ncb = Ncb1
          Nstore = FNCa - Npsntarget1 - Nresp1  
        end if
     end if
     PNlc = Nlc / FNCa
     jj = jj + 1  
  end do                        
  PNlcopt = Nlc / FNCa
  PNstoreopt = Nstore / FNCa
  PNcbopt = Ncb / FNCa
  PNetopt = Net / FNCa
  PNrespopt = Nresp / FNCa 
  
end subroutine NitrogenAllocation

!*****************************************************************************************************************
!calcualte the nitrogen investment for electron transport, carb10oxylation, respiration given a specified value 
!of nitrogen allocation in light capture [Nlc]. This equation are based on Ali et al 2015b.

subroutine Nitrogen_investments (KcKjFlag, FNCa, Nlc, forc_pbot10, relh10, &
     CO2a10, O2a10, PARi10, PARimx10, rb10, hourpd, tair10, tleafd10, tleafn10, &
     Kj2Kc, Wc2Wjb0, JmaxCoef, Fc, Fj, NUEc, NUEj, NUEcref, NUEjref, NUEr, Kc, &
     Kj, ci, Vcmax, Jmax, JmeanL, JmaxL, Net, Ncb, Nresp, PSN, RESP)
  implicit none
  integer,  intent (in) :: KcKjFlag                   !flag to indicate whether to update the Kc and Kj using the photosynthesis subroutine; 0--Kc and Kj need to be calculated; 1--Kc and Kj is prescribed.
  real(r8), intent (in) :: FNCa                       !Area based functional nitrogen content (g N/m2 leaf)
  real(r8), intent (in) :: Nlc                        !nitrogen content for light capture(g N/m2 leaf)
  real(r8), intent (in) :: forc_pbot10                !10-day mean air pressure (Pa)
  real(r8), intent (in) :: relh10                     !10-day mean relative humidity (unitless)
  real(r8), intent (in) :: CO2a10                     !10-day mean CO2 concentration in the air (Pa)
  real(r8), intent (in) :: O2a10                      !10-day mean O2 concentration in the air (Pa)
  real(r8), intent (in) :: PARi10                     !10-day mean photosynthetic active radiation on in a canopy (umol/m2/s)
  real(r8), intent (in) :: PARimx10                   !10-day mean 24hr maximum photosynthetic active radiation on in a canopy (umol/m2/s)
  real(r8), intent (in) :: rb10                       !10-day mean boundary layer resistance (s/m)
  real(r8), intent (in) :: hourpd                     !hours of light in a the day (hrs)
  real(r8), intent (in) :: tair10                     !10-day running mean of the 2m temperature (oC)
  real(r8), intent (in) :: tleafd10                   !10-day mean daytime leaf temperature (oC) 
  real(r8), intent (in) :: tleafn10                   !10-day mean nighttime leaf temperature (oC) 
  real(r8), intent (in) :: Kj2Kc                      !ratio:  Kj / Kc
  real(r8), intent (in) :: Wc2Wjb0                    !the baseline ratio of rubisco-limited rate vs light-limited photosynthetic rate (Wc:Wj)
  real(r8), intent (in) :: JmaxCoef                   !coefficient determining the response of electron transport rate to light availability (unitless) and humidity
  real(r8), intent (in) :: Fc                         !the temperature adjustment factor for Vcmax 
  real(r8), intent (in) :: Fj                         !the temperature adjustment factor for Jmax 
  real(r8), intent (in) :: NUEc                       !nitrogen use efficiency for carboxylation 
  real(r8), intent (in) :: NUEj                       !nitrogen use efficiency for electron transport
  real(r8), intent (in) :: NUEcref                    !nitrogen use efficiency for carboxylation under reference climates
  real(r8), intent (in) :: NUEjref                    !nitrogen use efficiency for electron transport under reference climates
  real(r8), intent (in) :: NUEr                       !nitrogen use efficiency for respiration
  real(r8), intent (inout) :: Kc                      !conversion factors from Vc,max to Wc 
  real(r8), intent (inout) :: Kj                      !conversion factor from electron transport rate to Wj 
  real(r8), intent (inout) :: ci                      !inter-cellular CO2 concentration (Pa) 
  real(r8), intent (out) :: Vcmax                     !the maximum carboxyaltion rate (umol/m2/s) 
  real(r8), intent (out) :: Jmax                      !the maximum electron transport rate (umol/m2/s) 
  real(r8), intent (out) :: JmaxL                     !the electron transport rate with maximum daily radiation (umol/m2/s)  
  real(r8), intent (out) :: JmeanL                    !the electron transport rate with mean radiation (umol/m2/s) 
  real(r8), intent (out)  :: Net                      !nitrogen content for electron transport(g N/m2 leaf)
  real(r8), intent (out)  :: Ncb                      !nitrogen content for carboxylation(g N/m2 leaf)
  real(r8), intent (out)  :: Nresp                    !nitrogen content for respiration(g N/m2 leaf)
  real(r8), intent (out)  :: PSN                      !daily photosynthetic rate(g C/day/m2 leaf)
  real(r8), intent (out)  :: RESP                     !daily respiration rate(g C/day/m2 leaf)
  !-------------------------------------------------------------------------------------------------------------------------------
  !intermediate variables
  real(r8) :: A                                       !Gross photosynthetic rate (umol CO2/m2/s)
  real(r8) :: Wc2Wj                                   !ratio: Wc/Wj  
  real(r8) :: ELTRNabsorb                             !absorbed electron rate, umol electron/m2 leaf /s
  real(r8) :: Jmaxb0act                               !base value of Jmax (umol/m2/s) 
  real(r8) :: theta_cj                                !interpolation coefficient
  real(r8) :: theta                                   !light absorption rate (0-1)
  real(r8) :: Vcmaxnight                              !Vcmax during night (umol/m2/s)
  real(r8) :: Wc                                      !rubisco-limited photosynthetic rate (umol/m2/s)
  real(r8) :: Wj                                      !light limited photosynthetic rate (umol/m2/s)
  real(r8) :: NUECHG                                  !the nitrogen use efficiency change under current conidtions compared to reference climate conditions (25oC and 385 ppm )
  real(r8), parameter :: leaf_mr_vcm = 0.015_r8       !Scalar constant of leaf respiration with Vcmax (should use parameter in CanopyStateMod)
  
  theta_cj = 0.95_r8
  theta = 0.292_r8 / (1.0_r8 + 0.076_r8 / (Nlc * Cb))
  ELTRNabsorb = theta * PARi10
  Jmaxb0act = Jmaxb0 * FNCa * Fj
  Jmax = Jmaxb0act + JmaxCoef * ELTRNabsorb
  JmaxL = theta * PARimx10 / (sqrt(1.0_r8 + (theta * PARimx10 / Jmax)**2.0_r8))        
  NUEchg = (NUEc / NUEcref) * (NUEjref / NUEj)
  Wc2Wj = Wc2Wjb0 * (NUEchg**0.5_r8)
  Wc2Wj = min(1.0_r8, Wc2Wj)
  Vcmax = Wc2Wj * JmaxL * Kj2Kc
  JmeanL = theta * PARi10 / (sqrt(1.0_r8 + (ELTRNabsorb / Jmax)**2.0_r8))
  if(KcKjFlag.eq.0)then      !update the Kc,Kj, anc ci information
     call Photosynthesis_luna(forc_pbot10, tleafd10, relh10, CO2a10, O2a10,rb10, Vcmax, JmeanL, ci, Kc, Kj, A) 
  else
    Wc = Kc * Vcmax
    Wj = Kj * JmeanL
    A = (1.0_r8 - theta_cj) * max(Wc, Wj) + theta_cj * min(Wc, Wj) 
  endif
  PSN = Cv * A * hourpd
  Vcmaxnight = VcmxTKattge(tair10, tleafn10) / VcmxTKattge(tair10, tleafd10) * Vcmax
  RESP = Cv * leaf_mr_vcm * (Vcmax * hourpd + Vcmaxnight * (24.0_r8 - hourpd))
  Net = Jmax / Fj
  Ncb = Vcmax / Fc
  Nresp = RESP / NUEr

end subroutine Nitrogen_investments



!********************************************************************************************************************
! Calculate the photosynthesis by solving the following 3 equations for 3 unknowns (A, gs, Ci): Farquahr's non-linear equation (A versus Ci), 
! Ball-Berry equation (gs versus A) and the diffusion equation (A = gs * (Ca - Ci). The approach taken is the following; Solve the 3 equations for 
! two phases. First phase is where Rubisco is limiting (Wc <= Wj) and second phase is where light is limiting (Wj > Wc).    

subroutine Photosynthesis_luna(forc_pbot, tleafd, relh, CO2a,O2a, rb, Vcmax, JmeanL, ci, Kc, Kj, A)
  implicit none 
  real(r8), intent (in) :: forc_pbot                  !air presure (Pa)  
  real(r8), intent (in) :: tleafd                     !daytime leaf temperature (oC) 
  real(r8), intent (in) :: relh                       !relative humidity (unitless)
  real(r8), intent (in) :: CO2a                       !atmospheric CO2 partial pressure(Pa)
  real(r8), intent (in) :: O2a                        !atmospheric O2 partial pressure(Pa)
  real(r8), intent (in) :: rb                         !boundary layer resistance (s/m)
  real(r8), intent (in) :: Vcmax                      !maximum carboxylation rate (umol/m2/s)
  real(r8), intent (in) :: JmeanL                     !average electron transport rate (umol/m2/s)
  real(r8), intent (out):: ci                         !inter-cellular CO2 concentration (ppm)
  real(r8), intent (out):: Kc                         !conversion factors for Vc,max to Wc
  real(r8), intent (out):: Kj                         !conversion factors for Jmax to Wj 
  real(r8), intent (out):: A                          !g dry mass photosynthesized per day
    
  !-------------------------------------------------------------------------------------------------------------------------------
  !intermediate variables
  real(r8) :: awc                                     !second deminator term for rubsico limited carboxylation rate based on Farquhar model
  real(r8) :: cf                                      !conversion factor of resistance: m**2/umol -> s/m
  real(r8) :: bp                                      !maximum stomatal resistance
  real(r8) :: mpe                                     !plant functional type dependent parameter for stomatal conductance 
  real(r8) :: rs                                      !stomatal resistance (s/m)
  real(r8) :: r1                                      !root1 of quadratic equations
  real(r8) :: r2                                      !root2 of quadratic equations
  real(r8) :: Wc                                      !rubisco-limited photosynthetic rate (umol/m2/s)
  real(r8) :: Wj                                      !light-limited photosynthetic rate (umol/m2/s)
  real(r8) :: k_o                                     !Michaelis-menten constant for O2 in Farquhar's model
  real(r8) :: k_c                                     !Michaelis-menten constant for CO2 in Farquhar's model
  real(r8) :: CO2c                                    !partial pressure of CO2 (Pa)
  real(r8) :: O2c                                     !partial pressure of oxygen (Pa)
  real(r8) :: c_p                                     !Michaelis-menten constant for Farquhar's model related to rubisco specificity factor
  real(r8) :: tdayk                                   !daytime temperature in Kelvin
  real(r8) :: ciold                                   !old value of inter-cellular CO2 concentration for convergence check
  real(r8) :: bbb                                     !Ball-Berry minimum leaf conductance (umol H20/m2/s) 
  real(r8) :: mbb                                     !Ball-Berry slope of conductance photosynthesis relationship (stressed)
  real(r8) :: gs_mol                                  !leaf stomatal conductance (umol H20/m2/s)
  real(r8) :: gb_mol                                  !leaf boundary layer conductance (umol H20/m2/s)
  real(r8) :: aquad                                   !terms of quadratic equations
  real(r8) :: bquad                                   !terms of quadratic equations
  real(r8) :: cquad                                   !terms of quadratic equations
  real(r8) :: phi                                     !terms of quadratic equations
  real(r8) :: rsmax0                                  !maximum stomata conductance (s/m)
  real(r8) :: tleaf                                   !daytime leaf temperature (oC)
  real(r8) :: tleafk                                  !the temperature of the leaf in Kelvin
  real(r8) :: theta_cj                                !the interpolation coefficient for Wj and Wc
  real(r8) :: relhc                                   !constrained relative humidity (unitless)
  integer  :: i                                       !index record the number of iterations
  
  theta_cj = 0.95_r8
  rsmax0 = 2.0_r8 * 1.0e4_r8
  bp = 2000.0_r8
  tleaf = tleafd
  tleafk = tleaf + tfrz
  aquad = 1.0_r8
  relhc = max(minrelh, relh)
  bbb = 1.0_r8 / bp
  mbb = mp
  CO2c = CO2a 
  O2c = O2a 
  ci = 0.7_r8 * CO2c  
  ciold = ci - 0.02_r8
  cf = forc_pbot / (8.314_r8 * tleafk) * 1.0e6_r8
  gb_mol = cf / rb
  k_c = kc25 * exp((79430.0_r8 / (8.314_r8 * (25.0_r8 + tfrz))) * (1.0_r8 - (tfrz + 25.0_r8) / (tfrz + tleaf)))
  k_o = ko25 * exp((36380.0_r8 / (8.314_r8 * (25.0_r8 + tfrz))) * (1.0_r8 - (tfrz + 25.0_r8) / (tfrz + tleaf)))
  c_p = Cp25 * exp((37830.0_r8 / (8.314_r8 * (25.0_r8 + tfrz))) * (1.0_r8 - (tfrz + 25_r8) / (tfrz + tleaf)))
  awc = k_c * (1.0_r8 + O2c / k_o)
  i = 1
  do while (abs(ci - ciold) > 0.01_r8 .and. i < 100)   ! for RUBISCO limitation
        i = i + 1                
        ciold = ci
        Kc = max(ci - c_p, 0.0_r8) / (ci + awc)
        Wc = Kc * Vcmax
        gs_mol = bbb + mbb * Wc / CO2c * forc_pbot * relhc
        phi = forc_pbot * (1.37_r8 * gs_mol + 1.6_r8 * gb_mol) / (gb_mol * gs_mol)
        bquad = awc - CO2c + phi * Vcmax
        cquad = -(c_p * phi * Vcmax + awc * CO2c)
        call Quadratic(aquad, bquad, cquad, r1, r2)
        ci = max(r1, r2)
        if (ci < 0.0_r8) ci = c_p + 0.5_r8 * ciold
  end do
  Kj = max(ci - c_p, 0.0_r8) / (4.0_r8 * ci + 8.0_r8 * c_p)
  Kc = max(ci - c_p, 0.0_r8) / (ci + awc)
  Wc = Kc * Vcmax
  Wj = Kj * JmeanL
  ciold = ci - 0.02_r8
  if (Wj < Wc)  then               !light limitation
    i = 1
    do while (abs(ci - ciold) > 0.01_r8 .and. i < 100)
         i = i + 1                     
         ciold = ci
         gs_mol = bbb + mbb * Wj / CO2c * forc_pbot * relhc
         phi = forc_pbot * (1.37_r8 * gs_mol + 1.6_r8 * gb_mol) / (gb_mol * gs_mol)
         bquad = 2.0_r8 * c_p - CO2c + phi * JmeanL / 4.0_r8
         cquad = -(c_p * phi * JmeanL / 4.0_r8 + 2.0_r8 * c_p * CO2c)
         call Quadratic(aquad, bquad, cquad, r1, r2)
         ci = max(r1, r2)
         if (ci < 0.0_r8) ci = c_p + 0.5_r8 * ciold
         Kj = max(ci - c_p, 0.0_r8) / (4.0_r8 * ci + 8.0_r8 * c_p)
         Wj = Kj * JmeanL
    end do
    Kj = max(ci - c_p, 0.0_r8) / (4.0_r8 * ci + 8.0_r8 * c_p)
    Kc = max(ci - c_p, 0.0_r8) / (ci + awc)
    Wc = Kc * Vcmax
    Wj = Kj * JmeanL
  end if                
  A = (1.0_r8 - theta_cj) * max(Wc, Wj) + theta_cj * min(Wc, Wj)   !use this instead of the quadratic to avoid values not in the range of wc and wj
  rs = cf / gs_mol
  rs =  min(rsmax0, rs)
                      
end subroutine Photosynthesis_luna



!********************************************************************************************************************************************************************** 
!Calculate the reference nitrogen use effieciency dependence on CO2 and leaf temperature

subroutine NUEref(NUEjref,NUEcref,Kj2Kcref)
  implicit none
  real(r8), intent (out):: NUEjref                    !nitrogen use efficiency for electron transport under refernce environmental conditions (25oC and 385 ppm co2)
  real(r8), intent (out):: NUEcref                    !nitrogen use efficiency for carboxylation under reference environmental conditions  (25oC and 385 ppm co2)
  real(r8), intent (out):: Kj2Kcref                   !the ratio of Wc to Wj under reference (25oC and 385 ppm co2) conditions  
  !---------------------------------------------
  !intermediate variables
  real(r8) :: Fj                                      !the temperature adjust factor for Jmax 
  real(r8) :: Fc                                      !the temperatuer adjust factor for Vcmax 
  real(r8) :: tgrow                                   !10 day mean growth temperature (oC), 24 hour mean temperature
  real(r8) :: tleaf                                   !leaf temperature (oC)
  real(r8) :: CO2c                                    !CO2 concentration (ppm)
  real(r8) :: O2c                                     !O2 concentration (ppm) 
  real(r8) :: k_o                                     !Rubsico O2 specifity
  real(r8) :: k_c                                     !Rubsico CO2 specifity
  real(r8) :: awc                                     !second deminator term for rubsico limited carboxylation rate based on Farquhar model
  real(r8) :: c_p                                     !CO2 compenstation point (Pa)
  real(r8) :: ci                                      !leaf internal [CO2] (Pa)
  real(r8) :: Kc                                      !converstion factor from Vcmax to Wc
  real(r8) :: Kj                                      !converstion factor from J to Wc

  tgrow   = 25.0_r8
  tleaf   = 25.0_r8
  Fc = VcmxTKattge(tgrow, tleaf) * Fc25
  Fj = JmxTKattge(tgrow, tleaf) * Fj25
  CO2c = co2ref * forc_pbot_ref * 1.0e-6_r8 !pa
  O2c = O2ref * forc_pbot_ref * 1.0e-6_r8   !pa
  k_c = Kc25 * exp((79430.0_r8 / (rgas*1.e-3_r8 * (25.0_r8 + tfrz))) * (1.0_r8 - (tfrz + 25.0_r8) / (tfrz + tleaf)))
  k_o = Ko25 * exp((36380.0_r8 / (rgas*1.e-3_r8 * (25.0_r8 + tfrz))) * (1.0_r8 - (tfrz + 25.0_r8) / (tfrz + tleaf)))
  c_p = Cp25 * exp((37830.0_r8 / (rgas*1.e-3_r8 * (25.0_r8 + tfrz))) * (1.0_r8 - (tfrz + 25.0_r8) / (tfrz + tleaf))) 
  awc  = k_c * (1.0_r8+O2c/k_o)
  ci = 0.7_r8 * CO2c
  Kj = max( ci-c_p,0.0_r8 ) / ( 4.0_r8*ci + 8.0_r8*c_p )
  Kc = max( ci-c_p,0.0_r8 ) / (ci+awc)  
  NUEjref  = Kj * Fj
  NUEcref = Kc * Fc
  Kj2Kcref = Kj / Kc
  
end subroutine NUEref


!******************************************************************************************************************** 
!Calculate the Nitrogen use effieciency dependence on CO2 and leaf temperature

subroutine NUE(O2a, ci, tgrow, tleaf, NUEj,NUEc,Kj2Kc)
  implicit none
  real(r8), intent (in) :: o2a                        !air O2 partial presuure (Pa)
  real(r8), intent (in) :: ci                         !leaf inter-cellular [CO2] (PPM)
  real(r8), intent (in) :: tgrow                      !10 day growth temperature (oC), 24 hour mean temperature
  real(r8), intent (in) :: tleaf                      !leaf temperature (oC)
  real(r8), intent (out):: NUEj                       !nitrogen use efficiency for electron transport under refernce environmental conditions (25oC and 385 ppm co2)
  real(r8), intent (out):: NUEc                       !nitrogen use efficiency for carboxylation under reference environmental conditions  (25oC and 385 ppm co2)
  real(r8), intent (out):: Kj2Kc                      !the ratio of Kj to Kc 
  !------------------------------------------------
  !intermediate variables
  real(r8) :: Fj                                      !the temperatuer adjust factor for Jmax 
  real(r8) :: Fc                                      !the temperatuer adjust factor for Vcmax 
  real(r8) :: Kc                                      !conversion factor from Vcmax to Wc 
  real(r8) :: Kj                                      !conversion factor from J to W 
  real(r8) :: k_o                                     !Rubsico O2 specifity
  real(r8) :: k_c                                     !Rubsico CO2 specifity
  real(r8) :: awc                                     !second deminator term for rubsico limited carboxylation rate based on Farquhar model
  real(r8) :: c_p                                     !CO2 compenstation point (Pa)
  
  Fc = VcmxTKattge(tgrow, tleaf) * Fc25
  Fj = JmxTKattge(tgrow, tleaf) * Fj25
  k_c = Kc25 * exp((79430.0_r8 / (rgas*1.e-3_r8 * (25.0_r8 + tfrz))) * (1.0_r8 - (tfrz + 25.0_r8) / (tfrz + tleaf)))
  k_o = Ko25 * exp((36380.0_r8 / (rgas*1.e-3_r8 * (25.0_r8 + tfrz))) * (1.0_r8 - (tfrz + 25.0_r8) / (tfrz + tleaf)))
  c_p = Cp25 * exp((37830.0_r8 / (rgas*1.e-3_r8 * (25.0_r8 + tfrz))) * (1.0_r8 - (tfrz + 25.0_r8) / (tfrz + tleaf)))
  awc = k_c * ( 1.0_r8 + O2a/k_o )
  Kj = max( ci-c_p,0.0_r8 ) / ( 4.0_r8*ci + 8.0_r8*c_p )
  Kc = max( ci-c_p,0.0_r8 ) / ( ci+awc )
  NUEj = Kj * Fj
  NUEc = Kc * Fc  
  Kj2Kc = Kj / Kc
end subroutine NUE

!************************************************************************************************************************************************
!Calculate the temperature response for Vcmax; assuming temperature acclimation as in CLM4.5, based on Kattge and Knorr  2007

real(r8) function VcmxTKattge(tgrow, tleaf)
  implicit none
  real(r8), intent(in):: tgrow !daytime and nightime growth temperature (oC) 
  real(r8), intent(in):: tleaf !leaf temperature (oC)  
  real(r8) :: TlimVcmx         !Vcmax activation energy
  real(r8) :: Vcmxf1           !Vcmax coef1
  real(r8) :: Vcmxf2           !Vcmax coef2
  real(r8) :: Vcmxf3           !Vcmax coef3 

  TlimVcmx = 668.39_r8- 1.07_r8 * (min(max(tgrow, 11.0_r8), 35.0_r8))
  Vcmxf1 = 1.0_r8 + exp((TlimVcmx * (25.0_r8 + tfrz) - 200000.0_r8) / (rgas*1.e-3_r8 * (25.0_r8 + tfrz)))
  Vcmxf2 = exp((72000.0_r8 / (rgas*1.e-3_r8 * (25.0_r8 + tfrz))) * (1.0_r8 - (tfrz+ 25.0_r8) / (tfrz + tleaf)))
  Vcmxf3 = 1.0_r8 + exp((TlimVcmx * (tleaf + tfrz) - 200000.0_r8) / (rgas*1.e-3_r8 * (tleaf + tfrz)))
  VcmxTKattge = Vcmxf1 * Vcmxf2 / Vcmxf3        
        
end function VcmxTKattge

!************************************************************************************************************************************************
!Calculate the temperature response for Jmax; assuming temperature acclimation as in CLM4.5, based on Kattge and Knorr  2007

real(r8) function JmxTKattge(tgrow, tleaf)
  implicit none
  real(r8), intent(in):: tgrow !daytime and nightime growth temperature (oC) 
  real(r8), intent(in):: tleaf !leaf temperature (oC)
  real(r8) :: TlimJmx          !Jmax activation energy
  real(r8) :: Jmxf1            !Jmax coef1
  real(r8) :: Jmxf2            !Jmax coef2
  real(r8) :: Jmxf3            !Jmax coef3

  TlimJmx = 659.7_r8 - 0.75_r8 * (min(max(tgrow, 11.0_r8), 35.0_r8))
  Jmxf1 = 1.0_r8 + exp((TlimJmx * (25.0_r8 + tfrz) - 200000.0_r8) / (rgas*1.e-3_r8 * (25.0_r8 + tfrz)))
  Jmxf2 = exp((50000.0_r8 / (rgas*1.e-3_r8 * (25.0_r8 + tfrz))) * (1._r8 - (tfrz + 25.0_r8) / (tleaf+tfrz)))
  Jmxf3 = 1.0_r8 + exp((TlimJmx * (tleaf + tfrz) - 200000.0_r8) / (rgas*1.e-3_r8 * (tleaf + tfrz)))
  JmxTKattge = Jmxf1 * Jmxf2 / Jmxf3        
        
end function JmxTKattge

!******************************************************************************************************************** 
!Calculate the temperature response for Vcmax; without assuming temperature acclimation and following Leunning 2002 Plant, Cell & Environment 

real(r8) function VcmxTLeuning(tgrow, tleaf)
  implicit none
  real(r8), intent(in) :: tgrow     !daytime and nightime growth temperature (oC) 
  real(r8), intent(in) :: tleaf     !leaf temperature (oC)
  real(r8) :: TlimVcmx              !Vcmax activation energy
  real(r8) :: Vcmxf1                !Vcmax coef1
  real(r8) :: Vcmxf2                !Vcmax coef2
  real(r8) :: Vcmxf3                !Vcmax coef3

  TlimVcmx = 486.0_r8
  Vcmxf1 = 1.0_r8 + exp((TlimVcmx * (25.0_r8 + tfrz) - 149252.0_r8) / (rgas*1.e-3_r8 * (25.0_r8 + tfrz)))
  Vcmxf2 = exp((73637.0_r8 / (rgas*1.e-3_r8 * (25.0_r8 + tfrz))) * (1._r8 - (tfrz + 25.0_r8) / (tfrz + tleaf)))
  Vcmxf3 = 1.0_r8 + exp((TlimVcmx * (tleaf + tfrz) - 149252.0_r8) / (rgas*1.e-3_r8 * (tleaf + tfrz)))
  VcmxTLeuning = Vcmxf1 * Vcmxf2 / Vcmxf3        
        
end function VcmxTLeuning

!******************************************************************************************************************** 
!Calculate the temperature response for Jmax; without assuming temperature acclimation and following Leunning 2002 Plant, Cell & Environment
real(r8) function  JmxTLeuning(tgrow, tleaf)
  implicit none
  real(r8), intent(in):: tgrow  !daytime and nightime growth temperature (oC) 
  real(r8), intent(in):: tleaf  !leaf temperature (oC)  
  real(r8) :: TlimJmx           !Jmax activation energy
  real(r8) :: Jmxf1             !Jmax coef1
  real(r8) :: Jmxf2             !Jmax coef2
  real(r8) :: Jmxf3             !Jmax coef3

  TlimJmx = 495.0_r8
  Jmxf1 = 1.0_r8 + exp((TlimJmx * (25.0_r8 + tfrz) - 152044.0_r8) / (rgas*1.e-3_r8 * (25.0_r8 + tfrz)))
  Jmxf2 = exp((50300.0_r8 / (rgas*1.e-3_r8 * (25.0_r8 + tfrz))) * (1._r8 - (tfrz + 25.0_r8) / (tfrz + tleaf)))
  Jmxf3 = 1.0_r8 + exp((TlimJmx * (tleaf + tfrz) - 152044.0_r8) / (rgas*1.e-3_r8 * (tleaf + tfrz)))
  JmxTLeuning = Jmxf1 * Jmxf2 / Jmxf3        
        
end function JmxTLeuning

!******************************************************************************************************************** 
!Calculate the temperature response for respiration, following Bernacchi PCE 2001

real(r8) function  RespTBernacchi(tleaf)
  implicit none
  real(r8), intent(in):: tleaf  !leaf temperature (oC)
  RespTBernacchi= exp(18.72_r8-46.39_r8/(rgas*1.e-6_r8 *(tleaf+tfrz)))
        
end function RespTBernacchi


!******************************************************************************************************************** 
!Calculate the soultion using the quadratic formula

subroutine  Quadratic(a,b,c,r1,r2) 
  implicit none
  real(r8), intent(in)  :: a   !coefficient a
  real(r8), intent(in)  :: b   !coefficient b
  real(r8), intent(in)  :: c   !coefficient c
  real(r8), intent(out) :: r1  !root one
  real(r8), intent(out) :: r2  !root one
  real(r8)  :: q               ! temporary term for quadratic solution
  
  r1 = 1.0e36_r8
  r2 = 1.0e36_r8
  
  if (a == 0.0_r8) return 

  if (b .GE. 0.0_r8) then 
      q = -0.5_r8 * (b + sqrt(b*b - 4.0_r8*a*c))
  else
      q = -0.5_r8 * (b - sqrt(b*b - 4.0_r8*a*c))
  end if 
  
  r1 = q / a
  
  if (q .NE. 0.0_r8)then
      r2 = c / q
  else 
      r2 = 1.0e36_r8
  end if
        
end subroutine Quadratic
	

end module LunaMod

