module CASAMod

#if (defined CASA)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CASAMod
!
! !DESCRIPTION:
! Terrestrial carbon cycle submodel patterned after the CASA model.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use clmtype
  use clm_atmlnd  , only : clm_a2l
  use clm_varcon  , only : denh2o, hvap, istsoil, tfrz, spval
  use clm_varcon  , only : istcrop
  use clm_varpar  , only : numpft, nlevsoi, nlevgrnd
  use clm_varctl  , only : iulog
  use spmdMod     , only : masterproc
  use CASAPhenologyMod, only : CASAPhenology
!
! !PUBLIC TYPES:
  implicit none
  save

  ! -----------------------------------------------------------------
  ! source file:       casa_params.h
  ! purpose:           CASA V2.1 parameters and variables
  ! modified for LSM/CASA interface by J.John (2001)
  ! -----------------------------------------------------------------

  ! Namelist parameters for CASA

  integer  :: spunup        ! 0=no, 1=yes (used with runtype non Continue only)
  integer  :: lalloc        ! 0=fixed allocation, 1=dynamic allocation
  integer  :: lnpp          ! 1=gpp*gppfact,2=fn(lgrow)*gppfact
  real(r8) :: q10
  character(len=256) :: fcpool          ! Carbon Pool initial state filename
  ! Logical to flag whether C pools have been read in on initial dataset
  logical :: cpool_inic = .false.

  ! Define parameters used in CASA 

  ! Pool Definitions

  integer, parameter :: nlive   = 3
  integer, parameter :: ndead   = 9
  integer, parameter :: npools  = nlive + ndead
  integer, parameter :: LEAF    = 1
  integer, parameter :: WOOD    = 2
  integer, parameter :: FROOT   = 3
  integer, parameter :: SURFMET = 4
  integer, parameter :: SURFSTR = 5
  integer, parameter :: SOILMET = 6
  integer, parameter :: SOILSTR = 7
  integer, parameter :: CWD     = 8
  integer, parameter :: SURFMIC = 9
  integer, parameter :: SOILMIC = 10
  integer, parameter :: SLOW    = 11
  integer, parameter :: PASSIVE = 12
  integer, parameter :: npool_types = 4
  integer, parameter :: LIVE_TYPE = 1
  integer, parameter :: LITTER_TYPE = 2
  integer, parameter :: SOIL_TYPE = 3
  integer, parameter :: CWD_TYPE = 4

  ! Tracer Definitions 

  integer, parameter :: ptrace   = 2
  integer, parameter :: Carbon   = 1
  integer, parameter :: Nitrogen = 2

  ! Respiration definitions

  integer, parameter :: nresp_pools = 14
  integer resp_pool_index(2,nresp_pools)   ! Indices of Respiring Pools in the
  integer pool_type_index(npools)          ! Index of pool type
  ! type definitions for pools in the order specified above
  data pool_type_index/            &
       LIVE_TYPE,   &
       LIVE_TYPE,   &
       LIVE_TYPE,   &
       LITTER_TYPE, &
       LITTER_TYPE, &
       LITTER_TYPE, &
       LITTER_TYPE, &
       CWD_TYPE,    &
       SOIL_TYPE,   &
       SOIL_TYPE,   &
       SOIL_TYPE,   &
       SOIL_TYPE/
  ! order that respiration is called
  data resp_pool_index/            &
       SLOW     ,PASSIVE, &
       SLOW     ,SOILMIC, &
       SURFMET  ,SURFMIC, &
       SURFSTR  ,SURFMIC, &
       SURFSTR  ,SLOW   , &
       SOILMET  ,SOILMIC, &
       SOILSTR  ,SOILMIC, &
       SOILSTR  ,SLOW   , &
       CWD      ,SURFMIC, &
       CWD      ,SLOW   , &
       SURFMIC  ,SLOW   , &
       SOILMIC  ,PASSIVE, &
       SOILMIC  ,SLOW   , &
       PASSIVE  ,SOILMIC/

  ! C:N ratio for pools

  real(r8) CNratio(npools)

  data CNratio/    &
       30.0_r8, &     ! C:N ratio of leaf pool
       130.0_r8, &    ! C:N ratio of wood pool
       55.0_r8, &     ! C:N ratio of froot pool
       30.0_r8, &     ! C:N ratio of surfmet pool
       50.0_r8, &     ! C:N ratio of surfstr pool
       25.0_r8, &     ! C:N ratio of soilmet pool
       50.0_r8, &     ! C:N ratio of soilstr pool
       135.0_r8, &    ! C:N ratio of cwd pool
       12.5_r8, &     ! C:N ratio of surfmic pool
       12.5_r8, &     ! C:N ratio of soilmic pool
       12.5_r8, &     ! C:N ratio of slow pool
       8.5_r8/        ! C:N ratio of passive pool

  ! LSM PFT assigned to CASA veg type
  ! 02/07/08 assign LSM PFT 14 (bare) to CASA type 8 (which does not exist)
  !       LSM:   1  2  3  4  5  6  7  8  9 10 11 12 13 14
  !      CASA:   4  5  1  2  6  7  9 11 10 10 12 12  6  8

  ! Mapping of LSM types to CLM types
  !      CLM :   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
  !      LSM :   1  1  2  3  3  5  4  4  7  8  9 10  6 13 11 11


  ! set values of parameters/constants used to determine PLAI

  ! min/max values of PLAI
  ! July 8 2002 - From Inez
  ! (2)  Leaf and Roots:  Max LAI
  ! Trees - PFT 1,2,3,4,5:  Max lai=6
  ! grass and crops - PFT 6,10, 11, 12, 13:  Max lai=3
  ! shrubs - PFT 7,8,9 - Max lai=1.5
  ! bare:  PFT 14:  max lai=0
  !(3)  Leaf and Roots - min LAI = 0.4 for photosynthesis.
  ! I got this by looking at non-zero gai values in vegconi.F (Figure 2). I am
  ! hopeful that LGROW will turn off photosynthesis during the winter.
  ! CAVEAT:  On page 19 of LSM documentation, it seems that leaf and stem
  ! areas enter into the calculation of albedo.  I am not sure where to set
  ! the min LAI, so that min LAI does not contribute to winter albedos.

  !.. 03/02/28 change plai_min to Dickinson value
  real(r8) plai_min(0:numpft)              ! min value of PLAI
  !.. mod 02/07/17  these are peak gai values from LSM1.1 (see vegconi.F)
  real(r8) plai_max(0:numpft)              ! max value of PLAI

  data plai_min/0.0_r8, 16*0.8_r8/
  data plai_max/0.0_r8, 5.0_r8, 5.0_r8, 2.6_r8, 4.5_r8, 4.5_r8, 3.0_r8, 4.7_r8, 4.7_r8, &
       1.0_r8, 0.9_r8, 1.4_r8, 3.5_r8, 3.5_r8, 3.5_r8, 3.0_r8, 3.0_r8/

  ! sla values below were used in LSM_CASA
  ! see initCasa below for current CLM values

  ! Specific Leaf area from Dickinson et al. (J.Clim, Nov. 1998)
  ! as inferred from summary of observed values by Schulze et al.
  ! Unit is m2 leaf area/kg C
  !     Veg Type                  SLA (m2/kg C)     LSM veg type
  !     crops                      60                  11,12
  !     short grass                40                  6,13
  !     needleleaf evergreen       10                  1
  !     deciduous needleleaf       30                  2
  !     deciduous broadleaf        30                  4,5
  !     broadleaf evergreen        25                  3
  !     evergreen shrubs           25                  7
  !     deciduous shrubs           25                  8,9
  !     tall grass                 35                  6,13
  !     tundra and semidesert      20                  10

  ! check SLA values for LSM veg types 5, 6/13, 8/9, 10, 11/12

  !     data sla/
  !    &         10.0, 30.0, 25.0, 30.0, 30.0, 40.0, 25.0,
  !    &         25.0, 25.0, 20.0, 60.0, 60.0, 35.0,  0.0/

  real(r8) lrage(0:numpft)            
  real(r8) woodage(0:numpft)         
  real(r8) litcn(0:numpft)            
  real(r8) lignin(0:numpft)         

  ! age characteristic parameters
  data lrage/ 0.00_r8,                                            &
       5.00_r8, 5.00_r8, 1.80_r8, 1.80_r8, 1.80_r8, 1.80_r8, 1.20_r8, 1.20_r8,  &
       1.00_r8, 1.00_r8, 2.80_r8, 2.80_r8, 1.50_r8, 1.80_r8, 1.00_r8, 1.00_r8/
  data woodage/0.00_r8,                                           &
       42.00_r8,42.00_r8,27.00_r8,41.00_r8,41.00_r8,25.00_r8,58.00_r8,58.00_r8, &
       5.50_r8, 1.00_r8, 5.50_r8, 5.50_r8, 0.00_r8,25.00_r8, 0.00_r8, 0.00_r8/

  ! litter characteristic parameters - lignin:N, C:N, lignin
  data litcn/  0.0_r8,                                           &
       80.0_r8, 80.0_r8, 50.0_r8, 40.0_r8, 40.0_r8, 50.0_r8, 50.0_r8, 50.0_r8, &
       65.0_r8, 50.0_r8, 50.0_r8, 50.0_r8, 50.0_r8, 50.0_r8, 40.0_r8, 40.0_r8/
  data lignin/ 0.0_r8,                                          &
       0.25_r8, 0.25_r8, 0.20_r8, 0.20_r8, 0.20_r8, 0.15_r8, 0.20_r8, 0.20_r8, &
       0.20_r8, 0.15_r8, 0.15_r8, 0.15_r8, 0.10_r8, 0.15_r8, 0.10_r8, 0.10_r8/

  ! Estimate of lignin content of wood C
  real(r8), parameter :: woodligninfract = 0.40_r8

  ! scaling factors for NPP
  real(r8), parameter :: gppfact = 0.5_r8    ! converts GPP to NPP

  ! 30cm depth for watdry, watopt, smoist, soilt (m)
  real(r8), parameter :: z30       = 0.3_r8  ! 30cm depth for smoist, soilt (m)

  !  set up array of nonwood types (based on pft) used in allocation
  !       lnonwood=1 if nonwoods (ivt = 6,7,8, 11, 12, 13, 14)
  !       lnonwood=0 if wood
  ! Look at LSM Table 2 - PFT composition of surface_types.  Note at
  ! "nonwoods" category.  All grass, shrubs are non-woods.  However, if you
  ! look at stembvt, artic shrub and arctic grass have 0.1 kg/m2.  So I
  ! suggest initial wood biomass = 0 for PFT types 6,7,8, 11, 12, 13, 14.
  ! For CLM, this corresponds to pfts: 9, 10, 13, 14, 15, 16

  integer :: lnonwood(0:numpft)
  data lnonwood/0,                       &
       0, 0, 0, 0, 0, 0, 0, 0,  &
       1, 1, 0, 0, 1, 1, 1, 1/
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: initCASA                 ! initialize the CASA submodel
  public :: CASA_ecosystemDyn        ! the main submodel interface
  public :: CASARest             ! CASA restart
!
! !REVISION HISTORY:
! Ported to the Community Land Model (CLM) by Jasmin John, Sam Levis, and
! Mariana Vertenstein.
! 2004.06.08 Vectorized and reformatted by Forrest Hoffman
!
!-----------------------------------------------------------------------

! !PRIVATE DATA MEMBERS:
!EOP

  real(r8) sla(0:numpft)          ! specific leaf area
  real(r8) leafmin(0:numpft)      ! min leafmass (g C/m2 ground)
  real(r8) leafmax(0:numpft)      ! max leafmass (g C/m2 ground)
  real(r8) solubfract(0:numpft)
  real(r8) lignineffect(0:numpft)
  real(r8) structurallignin(0:numpft)
  real(r8) fact_soilmic(0:numpft)
  real(r8) fact_slow(0:numpft)
  real(r8) fact_passive(0:numpft)
  real(r8) annK(0:numpft,npools)         
  real(r8) kdt(0:numpft,npools)         
!-----------------------------------------------------------------------

  private :: CASAPot_Evptr            ! potential evapotranspiration computation

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initCASA
!
! !INTERFACE:
  subroutine initCASA()
!
! !DESCRIPTION:
! Initialize the CASA submodel.
!
! !USES:
    use fileutils    , only : getfil
    use shr_const_mod, only : SHR_CONST_CDAY
    use decompMod    , only : get_proc_bounds, get_proc_global
    use clm_varctl   , only : nsrest, nsrStartup, nsrContinue
    use clm_varpar   , only : lsmlon, lsmlat, max_pft_per_gcell
    use spmdMod      , only : masterproc
    use clm_time_manager , only : get_step_size
    use pftvarcon    , only : noveg, nc3_nonarctic_grass, nc3crop, nirrig
    use ncdio_pio    
!
! !ARGUMENTS:
    implicit none
! 
! !LOCAL VARIABLES:
    real(r8), parameter :: plai_min_ic = 0.8_r8 ! init plai => Dickinson value
    real(r8), parameter :: secpy = 365._r8*SHR_CONST_CDAY ! no of secs/yr 

    ! local variables

    integer :: g,c,i,j,l,m,n,p,pi ! indices
    character(len=256) :: locfn   ! local file name
    type(file_desc_t)  :: ncid    ! netCDF file id
    integer :: varid         ! netCDF variable id
    integer :: begp, endp    ! per-proc beginning and ending pft indices
    integer :: begc, endc    ! per-proc beginning and ending column indices
    integer :: begl, endl    ! per-proc beginning and ending landunit indices
    integer :: begg, endg    ! per-proc gridcell ending gridcell indices
    integer :: numg          ! total number of gridcells across all processors
    integer :: numl          ! total number of landunits across all processors
    integer :: numc          ! total number of columns across all processors
    integer :: nump          ! total number of pfts across all processors
    logical :: readvar       ! is variable on file
    integer :: ier, ret      ! error return code

    real(r8) dtime           ! land model time step (sec)
    real(r8) lnscl
    real(r8) hardwire_sla(0:numpft)
    real(r8), pointer :: sumwts(:)
    real(r8), pointer :: vege_wts(:)
    real(r8), pointer :: wood_wts(:)
    real(r8), pointer :: vege_scale(:)
    real(r8), pointer :: wood_scale(:)
    real(r8), pointer :: rloc(:)

    ! pointers

    integer , pointer :: pgridcell(:) ! gridcell index of corresponding pft
    integer , pointer :: pcolumn(:)   ! pft's column
    integer , pointer :: plandunit(:) ! landunit index associated with pft
    integer , pointer :: npfts(:)     ! number of pfts on gridcell
    integer , pointer :: pfti(:)      ! initial pft on gridcell
    integer , pointer :: ltype(:)     ! landunit type for corresponding pft
    integer , pointer :: ivt(:)       ! pft vegetation type
    real(r8), pointer :: wtgcell(:)   ! pft weight relative to gridcell
    real(r8), pointer :: XSCpool(:)
    real(r8), pointer :: eff(:,:)
    real(r8), pointer :: frac_donor(:,:)
    real(r8), pointer :: Tpool_C(:,:) ! Total C pool size
    real(r8), pointer :: plai(:)      ! prognostic LAI (m2 leaf/m2 ground)
    real(r8), pointer :: sandfrac(:)
    real(r8), pointer :: clayfrac(:)
    real(r8), pointer :: co2flux(:)   ! net CO2 flux (gC/m2/s) [+ = to atm]
    real(r8), pointer :: fnpp(:)      ! NPP (gC/m2/sec)
    real(r8), pointer :: Resp_C(:,:)  ! could dimension by ndead, but caution!!!
    real(r8), pointer :: Cflux(:)
    real(r8), pointer :: watopt(:)    !optimal soil water content for et for top 30cm (mm3/mm3)
    real(r8), pointer :: watdry(:)    !soil water when et stops for top 30cm (mm3/mm3)
    real(r8), pointer :: sz(:)        !thickness of soil layers contributing to output
    real(r8), pointer :: watoptc(:)   !optimal soil water content for et for entire column (mm3/mm3)
    real(r8), pointer :: watdryc(:)   !soil water when et stops for entire column (mm3/mm3)
    real(r8), pointer :: szc(:)       !thickness of soil layers contributing to output
    real(r8), pointer :: watsat(:,:)  !saturated volumetric soil water content (porosity)
    real(r8), pointer :: sucsat(:,:)  ! minimum soil suction (mm)
    real(r8), pointer :: bsw(:,:)     !Clapp and Hornberger "b" (nlevsoi)  
    real(r8), pointer :: z(:,:)       ! soil layer depth (m)
    real(r8), pointer :: dz(:,:)      ! soil layer thickness (m)
    character(len=32) :: subname='initCasa' ! subroutine name
!
! !CALLED FROM:
! initialize in initializeMod
!
! !REVISION HISTORY:
! 2004.06.08 Vectorized and reformatted by Forrest Hoffman
!
!EOP
!-----------------------------------------------------------------------

    pgridcell  => clm3%g%l%c%p%gridcell
    pcolumn    => clm3%g%l%c%p%column
    plandunit  => clm3%g%l%c%p%landunit
    npfts      => clm3%g%npfts
    pfti       => clm3%g%pfti
    ltype      => clm3%g%l%itype
    ivt        => clm3%g%l%c%p%itype
    wtgcell    => clm3%g%l%c%p%wtgcell
    eff        => clm3%g%l%c%p%pps%eff  
    frac_donor => clm3%g%l%c%p%pps%frac_donor
    XSCpool    => clm3%g%l%c%p%pps%XSCpool
    Tpool_C    => clm3%g%l%c%p%pps%Tpool_C
    plai       => clm3%g%l%c%p%pps%plai
    sandfrac   => clm3%g%l%c%p%pps%sandfrac
    clayfrac   => clm3%g%l%c%p%pps%clayfrac
    co2flux    => clm3%g%l%c%p%pps%co2flux
    fnpp       => clm3%g%l%c%p%pps%fnpp
    Resp_C     => clm3%g%l%c%p%pps%Resp_C
    Cflux      => clm3%g%l%c%p%pps%Cflux  
    watdry     => clm3%g%l%c%p%pps%watdry 
    watopt     => clm3%g%l%c%p%pps%watopt    
    sz         => clm3%g%l%c%p%pps%sz 
    watdryc    => clm3%g%l%c%p%pps%watdryc
    watoptc    => clm3%g%l%c%p%pps%watoptc   
    szc        => clm3%g%l%c%p%pps%szc
    watsat     => clm3%g%l%c%cps%watsat
    sucsat     => clm3%g%l%c%cps%sucsat
    bsw        => clm3%g%l%c%cps%bsw
    z          => clm3%g%l%c%cps%z
    dz         => clm3%g%l%c%cps%dz

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    !iyf ========= set up decay constants ===============================
    ! EXPONENTIAL RATE CONSTANTS (per year) - convert to per second
    ! these are veg dependent (for live pools only), but not time dependent
    ! we have chosen to put inside time loop for clarity

    !iyf turnover times of live pools.  Want annK in sec-1.

    !iyf 02/07/11
    ! stressCD to be added to "annK(m,LEAF)" and "annK(m,FROOT)"
    ! in casa_litterfall.F

    do m = 0, numpft
       if(lrage(m) > 0.0_r8)then
          annK(m,LEAF)  = 1.0_r8/(lrage(m)*secpy)
          annK(m,FROOT) = 1.0_r8/(lrage(m)*secpy)
       else
          annK(m,LEAF)  = 1.0e-40_r8
          annK(m,FROOT) = 1.0e-40_r8
       end if

       if(woodage(m) > 0.0_r8)then
          annK(m,WOOD)  = 1.0_r8/(woodage(m)*secpy)
       else
          annK(m,WOOD)  = 1.0e-40_r8
       end if

       ! mod 03/07/16 multiply dead pool turnover times by 4 to account
       ! for effective turnover time (bgtemp*bgmoist) in casa_bgfluxes.F
       ! mod 03/08/19 multiply passive pool annK by 1/10
       ! want turnover time of passive pool to be 1000 years
       ! mod 03/11/21 remove factor of 4 from annK (as of lbgc.17f)

       !iyf: 1/(turnover times) for dead pools.  Want annK in sec-1.
       annK(m,SURFMET)    =       14.8_r8  /secpy
       annK(m,SURFMIC)    =       6.0_r8   /secpy
       annK(m,SURFSTR)    =       3.9_r8   /secpy
       annK(m,SOILMET)    =       18.5_r8  /secpy
       annK(m,SOILMIC)    =       7.3_r8   /secpy
       annK(m,SOILSTR)    =       4.9_r8   /secpy    ! 4.8 in casa v3.0
       annK(m,CWD)        =       0.2424_r8/secpy
       annK(m,SLOW)       =       0.2_r8   /secpy
       !orig   annK(m,PASSIVE)    =  0.0045/secpy
       !.. 03/04/03 change turnover time of passive pool to 50 years
       annK(m,PASSIVE)    = 0.1_r8 * 0.02_r8  /secpy
    end do

    ! Maximum RATE CONSTANTS FOR EACH POOL SCALED TO LENGTH OF TIME STEP
    ! For small delta_t, kdt   is the same as annK*delta_t
    !iyf:  Consider dM/dt =  - M/tau
    !iyf:  Analytic solution:  M(t) = M0 exp (-t/tau)
    !iyf:  Integrate Flux*dt from t=(n-1)*dt to t=n*dt:
    !iyf:  integral = M(n*dt) - M[(n-1)*dt]
    !iyf:           = - M[(n-1)*dt] {1 - exp [-dt/tau]}
    !iyf:    approx = M[(n-1)*dt] * [dt/tau]
    !iyf:  variable kdt was previously Krate in CASA

    ! NOTE: kdt will be used for WOOD and dead pools only
    !       so no need to worry about adding stressCD to annK here

    dtime = get_step_size()
    do n = 1, npools
       do m = 1, numpft
          kdt(m,n)= 1.0_r8 - (exp(-annK(m,n)*dtime))
       end do
       kdt(0,n)= 0.0_r8
    end do

    do m = 0, numpft

       ! LIGNIN TO NITROGEN SCALAR
       lnscl = litcn(m) * lignin(m) * 2.22_r8

       ! DETERMINE FRACTION OF LITTER THAT WILL BE METABOLIC FROM lignin:N RATIO
       solubfract(m) = 0.85_r8 - (0.018_r8 * lnscl)
       if(solubfract(m) < 0.0_r8)solubfract(m)=0.0_r8

       ! DETERMINE FRACTION OF C IN STRUCTURAL LITTER POOLS FROM LIGNIN
       ! For leaf and root
       structurallignin(m) = (lignin(m) * 0.65_r8 * 2.22_r8)   &
            / (1._r8 - solubfract(m))

       ! DETERMINE EFFECT OF LIGNIN CONTENT ON DECOMPOSITION RATES
       lignineffect(m) = exp(-3.0_r8 * structurallignin(m))

    end do

    ! 01/10/01 see Jim Randerson's version of CASA
    ! assign fact_soilmic, fact_slow, fact_passive for cultivation
    ! used in casa_bgfluxes.F

    do m = 0, numpft
       fact_soilmic(m) = 0.0_r8
       fact_slow(m)    = 0.0_r8
       fact_passive(m) = 0.0_r8
       if(m == nc3crop .or. m == nirrig)then    ! crops (corn, wheat in CLM)
          fact_soilmic(m) = 1.25_r8
          fact_slow(m)    = 1.50_r8
          fact_passive(m) = 1.50_r8
       else
          fact_soilmic(m) = 1.00_r8
          fact_slow(m)    = 1.00_r8
          fact_passive(m) = 1.00_r8
       end if
    end do

    ! sla values below are hardwired in biogeochem/VOCEmissionMod.F90

    hardwire_sla( 0) = 0._r8
    hardwire_sla( 1) = 0.0125_r8 !needleleaf
    hardwire_sla( 2) = 0.0125_r8 !Gordon Bonan suggests NET = 0.0076
    hardwire_sla( 3) = 0.0125_r8 !Gordon Bonan suggests NDT = 0.0200
    hardwire_sla( 4) = 0.0250_r8 !broadleaf
    hardwire_sla( 5) = 0.0250_r8 !Gordon Bonan suggests BET = 0.0178
    hardwire_sla( 6) = 0.0250_r8 !Gordon Bonan suggests BDT = 0.0274
    hardwire_sla( 7) = 0.0250_r8
    hardwire_sla( 8) = 0.0250_r8
    hardwire_sla( 9) = 0.0250_r8
    hardwire_sla(10) = 0.0250_r8
    hardwire_sla(11) = 0.0250_r8
    hardwire_sla(12) = 0.0200_r8 !grass
    hardwire_sla(13) = 0.0200_r8
    hardwire_sla(14) = 0.0200_r8
    hardwire_sla(15) = 0.0200_r8
    hardwire_sla(16) = 0.0200_r8 !numpft = 16

    do m = 0, numpft
       sla(m) = hardwire_sla(m) * 1000.0_r8      ! m2/kg C

       ! 02/07/22 assign min leafmass, max leafmass
       ! used in casa_bgfluxes.F

       if (sla(m) /= 0.0_r8) then
          leafmin(m) = (plai_min(m)/sla(m))*1.e3_r8   ! g C/m2 ground
          leafmax(m) = (plai_max(m)/sla(m))*1.e3_r8   ! g C/m2 ground
       else
          leafmin(m) = 0.0_r8
          leafmax(m) = 0.0_r8
       end if
    end do

    do p = begp,endp
       l = plandunit(p)

       plai(p) = 0.0_r8

       if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
          plai(p) = plai_min_ic

          ! CHECK !!!
          ! soil textures should be checked for soil depth
          ! These are now set in iniTimeConst
          !sandfrac(p) = sand3d(ixy(pgridcell(p)),jxy(pgridcell(p)),1)/100.0
          !clayfrac    = clay3d(ixy(pgridcell(p)),jxy(pgridcell(p)),1)/100.0

          eff(p, 1) =  0.45_r8    ! SLOW,PASSIVE
          eff(p, 2) =  0.45_r8    ! SLOW,SOILMIC
          eff(p, 3) =  0.40_r8    ! SURFMET,SURFMIC
          eff(p, 4) =  0.40_r8    ! SURFSTR,SURFMIC
          eff(p, 5) =  0.70_r8    ! SURFSTR,SLOW
          eff(p, 6) =  0.45_r8    ! SOILMET,SOILMIC
          eff(p, 7) =  0.45_r8    ! SOILSTR,SOILMIC
          eff(p, 8) =  0.70_r8    ! SOILSTR,SLOW
          eff(p, 9) =  0.40_r8    ! CWD,SURFMIC
          eff(p,10) =  0.70_r8    ! CWD,SLOW
          eff(p,11) =  0.40_r8    ! SURFMIC,SLOW
          eff(p,12) =  0.85_r8 - (0.68_r8 * (1.0_r8 - sandfrac(p)))  ! SOILMIC,PASSIVE
          eff(p,13) =  0.85_r8 - (0.68_r8 * (1.0_r8 - sandfrac(p)))  ! SOILMIC,SLOW
          eff(p,14) =  0.45_r8    ! PASSIVE,SOILMIC

          ! EXTRA RESPIRATION TRANSFER EFFICIENCIES

          frac_donor(p, 1) =  0.003_r8 + (0.009_r8*clayfrac(p))
          frac_donor(p, 2) =  1.0_r8 - frac_donor(p,1)
          frac_donor(p, 3) =  1.0_r8
          frac_donor(p, 4) =  1.0_r8 - structurallignin(ivt(p))
          frac_donor(p, 5) =  structurallignin(ivt(p))
          frac_donor(p, 6) =  1.0_r8
          frac_donor(p, 7) =  1.0_r8 - structurallignin(ivt(p))
          frac_donor(p, 8) =  structurallignin(ivt(p))
          frac_donor(p, 9) =  1.0_r8 - woodligninfract
          frac_donor(p,10) =  woodligninfract
          frac_donor(p,11) =  1.0_r8
          frac_donor(p,12) =  0.003_r8 + (0.032_r8*clayfrac(p))
          frac_donor(p,13) =  1.0_r8 - frac_donor(p,12)
          frac_donor(p,14) =  1.0_r8

       end if
    end do


    ! Initialize Pool Sizes to 0.0 or to spun up data (Tpool_C only)
    !!! 04/01/15 JJ
    !!! modify to set Tpool = 0 on initial start (Startup)
    !!! BUT, if SPUNPUP=1 and runtype NOT Continue, then read in initial pool size data 
    !!! ie on initial or branch run, if SPUNUP=1, use initial pool data.

    if (nsrest == nsrStartup .and. .not. cpool_inic) then
       if (masterproc) &
          write(iulog,*)'WARNING: Initializing Tpool_C and XSCpool to 0.0'
       do n = 1, npools
          do p = begp, endp
             Tpool_C(p,n) = 0.0_r8
          end do
       end do
       do p = begp, endp
          XSCpool(p) = 0.0_r8
       end do
    end if 

!!!  read spun up Tpool if available (use only for initial or branch runs)
!!!  don't overwrite restart values 

    if (SPUNUP == 1 .and. nsrest /= nsrContinue) then

       ! Allocate dynamic memory

       allocate(sumwts(begg:endg), vege_wts(begg:endg), wood_wts(begg:endg), stat=ier)
       if (ier /= 0) then
          call endrun('allocation error for sumwts, vege_wts, and wood_wts')
       end if
       allocate(vege_scale(begp:endp), wood_scale(begp:endp), stat=ier)
       if (ier /= 0) then
          call endrun('allocation error for vege_scale and wood_scale')
       end if
       allocate(rloc(begg:endg), stat=ier)
       if (ier /= 0) then
          call endrun('allocation error for rloc')
       end if

       ! Compute vegetated and woody scale factor to adjust initial carbon
       ! pools accounting for bare ground patches.  This will scale up initial
       ! values for vegetated PFTs mixed with bare ground and provide a zero
       ! multiplier for non-vegetated and non-woody patches.
       sumwts(begg:endg) = 0._r8
       vege_wts(begg:endg) = 0._r8
       wood_wts(begg:endg) = 0._r8

       do pi = 1,max_pft_per_gcell
          do g = begg, endg
             if (pi <= npfts(g)) then
                p = pfti(g) + pi - 1
                sumwts(g) = sumwts(g) + wtgcell(p)
                if (ivt(p) > noveg) vege_wts(g) = vege_wts(g) + wtgcell(p)
                if (ivt(p) > noveg .and. ivt(p) < nc3_nonarctic_grass ) wood_wts(g) = wood_wts(g) + wtgcell(p)
             end if
          end do
       end do

       do p = begp,endp
          g = pgridcell(p)
          if (ivt(p) > noveg .and. vege_wts(g) > 0._r8) then
             !vege_scale(p) = sumwts(g) / vege_wts(g)
             vege_scale(p) = 1._r8
          else
             vege_scale(p) = 0._r8
          end if
          if (ivt(p) > noveg .and. ivt(p) < nc3_nonarctic_grass .and. vege_wts(g) > 0._r8) then
             !wood_scale(p) = sumwts(g) / wood_wts(g)
             wood_scale(p) = 1._r8
          else
             wood_scale(p) = 0._r8
          end if
       end do

       ! Read file

       if (masterproc) then
          write(iulog,*)'Reading initial carbon pools from ',locfn
       end if
       call getfil(fcpool, locfn, 0)
       call ncd_pio_openfile (ncid, locfn, 0)
       call check_dim(ncid, 'longitude', lsmlon)
       call check_dim(ncid, 'latitude',  lsmlat)

       ! TPOOL_C_LEAF
       call ncd_io(ncid=ncid, varname='TPOOL_C_LEAF',flag='read', data=rloc, dim1name=nameg, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: TPOOL_C_LEAF NOT on fcpool file' )
       do p = begp,endp
          Tpool_C(p,LEAF) = rloc(pgridcell(p)) * vege_scale(p)
       end do

       ! TPOOL_C_WOOD
       call ncd_io(ncid=ncid, varname='TPOOL_C_WOOD',flag='read', data=rloc, dim1name=nameg, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: TPOOL_C_WOOD NOT on fcpool file' )
       do p = begp,endp
          Tpool_C(p,WOOD) = rloc(pgridcell(p)) * wood_scale(p)
       end do

       ! TPOOL_C_FROOT
       call ncd_io(ncid=ncid, varname='TPOOL_C_FROOT',flag='read', data=rloc, dim1name=nameg, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: TPOOL_C_FROOT NOT on fcpool file' )
       do p = begp,endp
          Tpool_C(p,FROOT) = rloc(pgridcell(p)) * vege_scale(p)
       end do

       ! TPOOL_C_SURFMET
       call ncd_io(ncid=ncid, varname='TPOOL_C_SURFMET',flag='read', data=rloc, dim1name=nameg, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: TPOOL_C_SURFMET NOT on fcpool file' )
       do p = begp,endp
          Tpool_C(p,SURFMET) = rloc(pgridcell(p)) * vege_scale(p)
       end do

       ! TPOOL_C_SURFSTR
       call ncd_io(ncid=ncid, varname='TPOOL_C_SURFSTR',flag='read', data=rloc, dim1name=nameg, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: TPOOL_C_SURFSTR NOT on fcpool file' )
       do p = begp,endp
          Tpool_C(p,SURFSTR) = rloc(pgridcell(p)) * vege_scale(p)
       end do

       ! TPOOL_C_SOILMET
       call ncd_io(ncid=ncid, varname='TPOOL_C_SOILMET',flag='read', data=rloc, dim1name=nameg, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: TPOOL_C_SOILMET NOT on fcpool file' )
       do p = begp, endp
          Tpool_C(p,SOILMET) = rloc(pgridcell(p)) * vege_scale(p)
       end do

       ! TPOOL_C_SOILSTR
       call ncd_io(ncid=ncid, varname='TPOOL_C_SOILSTR',flag='read', data=rloc, dim1name=nameg, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: TPOOL_C_SOILSTR NOT on fcpool file' )
       do p = begp,endp
          Tpool_C(p,SOILSTR) = rloc(pgridcell(p)) * vege_scale(p)
       end do

       ! TPOOL_C_CWD
       call ncd_io(ncid=ncid, varname='TPOOL_C_CWD',flag='read', data=rloc, dim1name=nameg, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: TPOOL_C_CWD NOT on fcpool file' )
       do p = begp,endp
          Tpool_C(p,CWD) = rloc(pgridcell(p)) * wood_scale(p)
       end do

       ! TPOOL_C_SURFMIC
       call ncd_io(ncid=ncid, varname='TPOOL_C_SURFMIC',flag='read', data=rloc, dim1name=nameg, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: TPOOL_C_SURFMIC NOT on fcpool file' )
       do p = begp,endp
          Tpool_C(p,SURFMIC) = rloc(pgridcell(p)) * vege_scale(p)
       end do

       ! TPOOL_C_SOILMIC
       call ncd_io(ncid=ncid, varname='TPOOL_C_SOILMIC',flag='read', data=rloc, dim1name=nameg, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: TPOOL_C_SOILMIC NOT on fcpool file' )
       do p = begp, endp
          Tpool_C(p,SOILMIC) = rloc(pgridcell(p)) * vege_scale(p)
       end do

       ! TPOOL_C_SLOW
       call ncd_io(ncid=ncid, varname='TPOOL_C_SLOW',flag='read', data=rloc, dim1name=nameg, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: TPOOL_C_SLOW  NOT on fcpool file' )
       do p = begp,endp
          Tpool_C(p,SLOW) = rloc(pgridcell(p)) * vege_scale(p)
       end do

       ! TPOOL_C_PASSIVE
       call ncd_io(ncid=ncid, varname='TPOOL_C_PASSIVE',flag='read', data=rloc, dim1name=nameg , readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: TPOOL_C_PASSIVE NOT on fcpool file' )
       do p = begp,endp
          Tpool_C(p,PASSIVE) = rloc(pgridcell(p)) * vege_scale(p)
       end do

       ! Close netcdf file

       call ncd_pio_closefile(ncid)

       ! Deallocate dynamic memory

       deallocate(sumwts, vege_wts, wood_wts)
       deallocate(vege_scale, wood_scale)
       deallocate(rloc)

       call casa_write_cpool()

    end if ! SPUNUP = 1

    ! Initialize all fluxes to 0._r8
    co2flux(begp:endp) = 0._r8
    fnpp(begp:endp) = 0._r8
    Cflux(begp:endp) = 0._r8
    Resp_C(begp:endp,1:npools) = 0._r8

    ! Initialize watdry, watopt, sz and watdryc, watoptc, szc
    do p = begp,endp
       l = plandunit(p)
       if (ltype(l) == istsoil .or. ltype(l) == istcrop) then

          ! top 30cm
          watdry(p)  = 0._r8
          watopt(p)  = 0._r8
          sz(p)      = 0._r8

          ! entire water column
          watdryc(p) = 0._r8
          watoptc(p) = 0._r8
          szc(p)     = 0._r8
       end if
    end do

    ! Compute watdry, watopt, and watdryc, watoptc in mm^3/mm^3
    do j = 1, nlevsoi
       do p = begp,endp
          c = pcolumn(p)
          l = plandunit(p)
          if (ltype(l) == istsoil .or. ltype(l) == istcrop) then

             ! top 30cm
             if (z(c,j)+0.5_r8*dz(c,j) <= z30) then
                watdry(p) = watdry(p) + watsat(c,j) * ((316230._r8/sucsat(c,j)) ** (-1._r8/bsw(c,j))) * dz(c,j)
                watopt(p) = watopt(p) + watsat(c,j) * ((158490._r8/sucsat(c,j)) ** (-1._r8/bsw(c,j))) * dz(c,j)
                sz(p)     = sz(p) + dz(c,j)
             end if

             ! entire water column
             watdryc(p) = watdryc(p) + watsat(c,j) * ((316230._r8/sucsat(c,j)) ** (-1._r8/bsw(c,j))) * dz(c,j)
             watoptc(p) = watoptc(p) + watsat(c,j) * ((158490._r8/sucsat(c,j)) ** (-1._r8/bsw(c,j))) * dz(c,j)
             szc(p)     = szc(p) + dz(c,j)
          end if
       end do
    end do

    do p = begp,endp
       l = plandunit(p)
       if (ltype(l) == istsoil .or. ltype(l) == istcrop) then

          ! top 30cm
          watdry(p)  = watdry(p)/sz(p)
          watopt(p)  = watopt(p)/sz(p)

          ! entire water column
          watdryc(p) = watdryc(p)/szc(p)
          watoptc(p) = watoptc(p)/szc(p)
       end if
    end do

 !  get fcap and wp - used in allocation
 !    wp = permanent wilting point (mm3/mm3)
 !    fcap = field capacity (mm3/mm3)
 !   do p = begp, endp
 !      l = plandunit(p)
 !      if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
 !        sand1 = sand(p)         ! already in pct
 !        clay1 = clay(p)         ! already in pct
 !        sand2 = sand1 * sand1
 !        clay2 = clay1 * clay1
 !        soilalpha = exp(-4.396 - 0.0715*clay1
 !    &             - 4.88e-4*sand2 - 4.285e-5*sand2*clay1)
 !        soilbeta = -3.140 - 0.00222*clay2
 !    &            - 3.484e-5*sand2*clay1
 ! convert to mm3  (x1000)  - CHECK units
 !        if ((soilalpha.ne.0.0).and.(soilbeta.ne.0.0)) then
 !          fcap(p) = ((0.3333/soilalpha)**(1.0/soilbeta))*1000.0
 !          wp(p) = ((15.0/soilalpha)**(1.0/soilbeta))*1000.0
 !        end if
 !      end if
 !   end do

  end subroutine initCASA

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: casa_write_cpool
!
! !INTERFACE:
  subroutine casa_write_cpool()
!
! !DESCRIPTION:
! Writes out a CPOOL_INITIAL file containing the mapped carbon pools.  It
! may be read later in initCASA() during initialization.  This file is not
! saved or shipped to the Mass Storage System since it is really just for
! debugging.  Carbon pool states are saved now in initial datasets.
!
! !USES:
    use decompMod    , only : get_proc_bounds, get_proc_global
    use subgridAveMod, only : p2g
    use domainMod    , only : ldomain, llatlon
    use clm_varpar   , only : lsmlon, lsmlat
    use spmdMod      , only : masterproc
    use ncdio_pio
!
! !ARGUMENTS:
    implicit none
!
! !LOCAL VARIABLES:
    type(file_desc_t):: ncid ! file id
    integer :: omode         ! netCDF mode returned
    integer :: dimid         ! netCDF dimension id
    integer :: begp, endp    ! per-proc beginning and ending pft indices
    integer :: begc, endc    ! per-proc beginning and ending column indices
    integer :: begl, endl    ! per-proc beginning and ending landunit indices
    integer :: begg, endg    ! per-proc gridcell ending gridcell indices
    integer :: numg          ! total number of gridcells across all processors
    integer :: numl          ! total number of landunits across all processors
    integer :: numc          ! total number of columns across all processors
    integer :: nump          ! total number of pfts across all processors
    integer :: ier           ! error flag
    integer :: n,m,g         ! index
    real(r8), pointer :: histi(:,:)
    real(r8), pointer :: histo(:,:)
    real(r8), pointer :: hist1do(:)
    character(len=32) :: subname='casa_write_cpool' ! subroutine name
    ! pointers
    real(r8), pointer :: Tpool_C(:,:) ! Total C pool size
!
! !CALLED FROM:
! initCASA() for debugging.
!
! !REVISION HISTORY:
! 2004.08.17 Created by Forrest Hoffman
!
!
!EOP
!-----------------------------------------------------------------------

    Tpool_C    => clm3%g%l%c%p%pps%Tpool_C

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    call ncd_pio_createfile(ncid, 'CPOOL_INITIAL.nc')

    call ncd_defdim(ncid, 'longitude', lsmlon, dimid)
    call ncd_defdim(ncid, 'latitude' , lsmlat, dimid)

    call ncd_defvar(ncid=ncid, varname='longitude', xtype=ncd_double, &
         dim1name='longitude', long_name='coordinate longitude', &
         units='degrees_east')
    call ncd_defvar(ncid=ncid, varname='latitude', xtype=ncd_double, &
         dim1name='latitude', long_name='coordinate latitude', &
         units='degrees_north')
    call ncd_defvar(ncid=ncid, varname='landfrac', xtype=ncd_double, &
         dim1name='longitude', dim2name='latitude', long_name='land fraction')
    call ncd_defvar(ncid=ncid, varname='landmask', xtype=ncd_int, &
         dim1name='longitude', dim2name='latitude', &
         long_name='land/ocean mask (0.=ocean and 1.=land)')
    call ncd_defvar(ncid=ncid, varname='TPOOL_C_LEAF', &
         xtype=ncd_double, dim1name='longitude', dim2name='latitude', &
         long_name='total C in leaf pool', &
         missing_value=spval, fill_value=spval)
    call ncd_defvar(ncid=ncid, varname='TPOOL_C_WOOD', &
         xtype=ncd_double, dim1name='longitude', dim2name='latitude', &
         long_name='total C in wood pool', &
         missing_value=spval, fill_value=spval)
    call ncd_defvar(ncid=ncid, varname='TPOOL_C_CWD', &
         xtype=ncd_double, dim1name='longitude', dim2name='latitude', &
         long_name='total C in coarse woody debris pool', &
         missing_value=spval, fill_value=spval)
    call ncd_defvar(ncid=ncid, varname='TPOOL_C_FROOT', &
         xtype=ncd_double, dim1name='longitude', dim2name='latitude', &
         long_name='total C in fine root pool', &
         missing_value=spval, fill_value=spval)
    call ncd_defvar(ncid=ncid, varname='TPOOL_C_SURFMET', &
         xtype=ncd_double, dim1name='longitude', dim2name='latitude', &
         long_name='total C in pool', &
         missing_value=spval, fill_value=spval)
    call ncd_defvar(ncid=ncid, varname='TPOOL_C_SURFSTR', &
         xtype=ncd_double, dim1name='longitude', dim2name='latitude', &
         long_name='total C in pool', &
         missing_value=spval, fill_value=spval)
    call ncd_defvar(ncid=ncid, varname='TPOOL_C_SOILMET', &
         xtype=ncd_double, dim1name='longitude', dim2name='latitude', &
         long_name='total C in pool', &
         missing_value=spval, fill_value=spval)
    call ncd_defvar(ncid=ncid, varname='TPOOL_C_SOILSTR', &
         xtype=ncd_double, dim1name='longitude', dim2name='latitude', &
         long_name='total C in pool', &
         missing_value=spval, fill_value=spval)
    call ncd_defvar(ncid=ncid, varname='TPOOL_C_SURFMIC', &
         xtype=ncd_double, dim1name='longitude', dim2name='latitude', &
         long_name='total C in pool', &
         missing_value=spval, fill_value=spval)
    call ncd_defvar(ncid=ncid, varname='TPOOL_C_SOILMIC', &
         xtype=ncd_double, dim1name='longitude', dim2name='latitude', &
         long_name='total C in pool', &
         missing_value=spval, fill_value=spval)
    call ncd_defvar(ncid=ncid, varname='TPOOL_C_SLOW', &
         xtype=ncd_double, dim1name='longitude', dim2name='latitude', &
         long_name='total C in pool', &
         missing_value=spval, fill_value=spval)
    call ncd_defvar(ncid=ncid, varname='TPOOL_C_PASSIVE', &
         xtype=ncd_double, dim1name='longitude', dim2name='latitude', &
         long_name='total C in pool', &
         missing_value=spval, fill_value=spval)

    call ncd_enddef(ncid)

    call ncd_io(varname='longitude', data=llatlon%lonc, ncid=ncid, flag='write')
    call ncd_io(varname='latitude' , data=llatlon%latc, ncid=ncid, flag='write')

    call ncd_io(varname='landfrac', data=ldomain%frac, ncid=ncid, &
         flag='write', dim1name=grlnd)
    call ncd_io(varname='landmask', data=ldomain%mask, ncid=ncid, &
         flag='write', dim1name=grlnd)

    allocate(histi(begp:endp, 1),histo(begg:endg,1),hist1do(begg:endg), stat=ier)
    if (ier /= 0) then
       call endrun( 'Allocation error for TPOOL_C write' )
    end if

    ! TPOOL_C_LEAF
    histi(begp:endp,1) = Tpool_C(begp:endp,LEAF)
    histo(begg:endg,1) = spval
    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, 1, &
       histi, histo, 'unity', 'unity', 'unity')
       hist1do(begg:endg) = histo(begg:endg, 1)
    call ncd_io(flag='write', varname='TPOOL_C_LEAF', dim1name=grlnd, &
       data=hist1do, ncid=ncid)

    ! TPOOL_C_WOOD
    histi(begp:endp,1) = Tpool_C(begp:endp,WOOD)
    histo(begg:endg,1) = spval
    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, 1, &
       histi, histo, 'unity', 'unity', 'unity')
       hist1do(begg:endg) = histo(begg:endg, 1)
    call ncd_io(flag='write', varname='TPOOL_C_WOOD', dim1name=grlnd, &
       data=hist1do, ncid=ncid)

    ! TPOOL_C_CWD
    histi(begp:endp,1) = Tpool_C(begp:endp,CWD)
    histo(begg:endg,1) = spval
    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, 1, &
       histi, histo, 'unity', 'unity', 'unity')
       hist1do(begg:endg) = histo(begg:endg, 1)
    call ncd_io(flag='write', varname='TPOOL_C_CWD', dim1name=grlnd, &
       data=hist1do, ncid=ncid)

    ! TPOOL_C_FROOT
    histi(begp:endp,1) = Tpool_C(begp:endp,FROOT)
    histo(begg:endg,1) = spval
    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, 1, &
       histi, histo, 'unity', 'unity', 'unity')
       hist1do(begg:endg) = histo(begg:endg, 1)
    call ncd_io(flag='write', varname='TPOOL_C_FROOT', dim1name=grlnd, &
       data=hist1do, ncid=ncid)

    ! TPOOL_C_SURFMET
    histi(begp:endp,1) = Tpool_C(begp:endp,SURFMET)
    histo(begg:endg,1) = spval
    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, 1, &
       histi, histo, 'unity', 'unity', 'unity')
       hist1do(begg:endg) = histo(begg:endg, 1)
    call ncd_io(flag='write', varname='TPOOL_C_SURFMET', dim1name=grlnd, &
       data=hist1do, ncid=ncid)

    ! TPOOL_C_SURFSTR
    histi(begp:endp,1) = Tpool_C(begp:endp,SURFSTR)
    histo(begg:endg,1) = spval
    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, 1, &
       histi, histo, 'unity', 'unity', 'unity')
       hist1do(begg:endg) = histo(begg:endg, 1)
    call ncd_io(flag='write', varname='TPOOL_C_SURFSTR', dim1name=grlnd, &
       data=hist1do, ncid=ncid)

    ! TPOOL_C_SOILMET
    histi(begp:endp,1) = Tpool_C(begp:endp,SOILMET)
    histo(begg:endg,1) = spval
    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, 1, &
       histi, histo, 'unity', 'unity', 'unity')
       hist1do(begg:endg) = histo(begg:endg, 1)
    call ncd_io(flag='write', varname='TPOOL_C_SOILMET', dim1name=grlnd, &
       data=hist1do, ncid=ncid)

    ! TPOOL_C_SOILSTR
    histi(begp:endp,1) = Tpool_C(begp:endp,SOILSTR)
    histo(begg:endg,1) = spval
    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, 1, &
       histi, histo, 'unity', 'unity', 'unity')
       hist1do(begg:endg) = histo(begg:endg, 1)
    call ncd_io(flag='write', varname='TPOOL_C_SOILSTR', dim1name=grlnd, &
       data=hist1do, ncid=ncid)

    ! TPOOL_C_SURFMIC
    histi(begp:endp,1) = Tpool_C(begp:endp,SURFMIC)
    histo(begg:endg,1) = spval
    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, 1, &
       histi, histo, 'unity', 'unity', 'unity')
       hist1do(begg:endg) = histo(begg:endg, 1)
    call ncd_io(flag='write', varname='TPOOL_C_SURFMIC', dim1name=grlnd, &
       data=hist1do, ncid=ncid)

    ! TPOOL_C_SOILMIC
    histi(begp:endp,1) = Tpool_C(begp:endp,SOILMIC)
    histo(begg:endg,1) = spval
    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, 1, &
       histi, histo, 'unity', 'unity', 'unity')
       hist1do(begg:endg) = histo(begg:endg, 1)
    call ncd_io(flag='write', varname='TPOOL_C_SOILMIC', dim1name=grlnd, &
       data=hist1do, ncid=ncid)

    ! TPOOL_C_SLOW
    histi(begp:endp,1) = Tpool_C(begp:endp,SLOW)
    histo(begg:endg,1) = spval
    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, 1, &
       histi, histo, 'unity', 'unity', 'unity')
       hist1do(begg:endg) = histo(begg:endg, 1)
    call ncd_io(flag='write', varname='TPOOL_C_SLOW', dim1name=grlnd, &
       data=hist1do, ncid=ncid)

    ! TPOOL_C_PASSIVE
    histi(begp:endp,1) = Tpool_C(begp:endp,PASSIVE)
    histo(begg:endg,1) = spval
    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, 1, &
       histi, histo, 'unity', 'unity', 'unity')
       hist1do(begg:endg) = histo(begg:endg, 1)
    call ncd_io(flag='write', varname='TPOOL_C_PASSIVE', dim1name=grlnd, &
       data=hist1do, ncid=ncid)

    deallocate(hist1do)
    deallocate(histo)
    deallocate(histi)

    call ncd_pio_closefile(ncid)

  end subroutine casa_write_cpool

!===============================================================================
!BOP
!
! !IROUTINE: casa_ecosystemDyn
!
! !INTERFACE:
  subroutine casa_ecosystemDyn(begc,endc, begp,endp, num_soilc,filter_soilc, &
                               num_soilp,filter_soilp, doalb, init)
!
! !DESCRIPTION:
!    Primary calling interface to the CASA submodel.
!
! !CALLED FROM:
!    clm_driver1 and initSurfAlb
!
! !REVISION HISTORY:
!    2004.06.08  F. Hoffman: Vectorized and reformatted
!    2008.11.12  B. Kauffman: include phenology, vegStruct, rename "ecosystemDyn"
! 
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: begc, endc      ! column bounds
    integer, intent(in) :: begp, endp      ! pft bounds
    integer, intent(in) :: num_soilc       ! number of soil points in column filter
    integer, intent(in) :: filter_soilc(:) ! column filter for soil points
    integer, intent(in) :: num_soilp       ! number of soil points in pft filter
    integer, intent(in) :: filter_soilp(:) ! pft filter for soil points
    logical, intent(in), optional :: doalb ! T => run phase call and albedo timestep
    logical, intent(in), optional :: init  ! T => init phase call: for albedo init
!
!EOP

    !--- local variables ---
    integer :: f, p   ! indicies
    logical :: linit  ! local init
    logical :: ldoalb ! local doalb

    ! implicit intent inout
    !============================================================
    real(r8), pointer :: plai(:)     ! prognostic LAI (m2 leaf/m2 ground)

    character(*),parameter :: subname='(casa_ecosystemDyn) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

    linit = .false.
    if (present(init)) lInit = init
    ldoalb = .false.
    if (present(doalb)) ldoalb = doalb

    ! phenology
    call CASAPhenology(begp, endp, num_soilp, filter_soilp)

    if ( linit) then
       if (masterproc) write(iulog,*) subName,"initialize option is ON"

       ! explicitly set prognostic LAI to zero when called from initSurfAlb
       plai       => clm3%g%l%c%p%pps%plai
       do f = 1,num_soilp
          p = filter_soilp(f)

          plai(p) = 0.0_r8
       end do

    else
       !--- potential evapo-transpiration ---
       call CASAPot_Evptr(begp, endp, num_soilp, filter_soilp)

       !--- carbon allocation ---
       call casa_allocate(begp, endp, num_soilp, filter_soilp)

       !--- plant net primary production ---
       call casa_npp(begp,endp, num_soilp,filter_soilp)

       ! bgfluxes (litterfall, soil respiration)
       !  Compute Carbon flux to send to atm
       !  fnpp (gC/m2/sec), Cflux (gC/m2/sec), co2flux (gC/m2/sec)
       call casa_bgfluxes(begp, endp, num_soilp, filter_soilp)

       call CASASummary(begp, endp, num_soilp, filter_soilp )

    end if

  end subroutine casa_ecosystemDyn

!===============================================================================
!BOP
!
! !IROUTINE: CASAPot_Evptr
!
! !INTERFACE:
  subroutine CASAPot_Evptr(lbp, ubp, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Potential Evapotranspiration
! Priestely-Taylor Equation
! Baldocchi et al.  (2000): Climate and vegetation controls on boreal 
! zone energy exchange. Global Change Biology, 6, (Suppl. 1), 69-83.
!
! In Baldocchi et al, PET (equilibrium evaporation) is calculated 
! for time step - as he compares the instantaneous evapotranspiration
! to the eqm evapotranspiration.
! I think that this is simplest, and avoids the definition of 
! an averaging period.  iyf 2002/05/09
!
! The following calculation is done at every land point.
! No explicit discrimination among veg type or soil type
! Local ecology is implictly dealt with in the energy fluxes.
!
!*************************************************************************
! the model partitions the latent heat flux into three components:
!    o fcev: evaporation of intercepted water
!    o fctr: transpiration
!    o fgev: soil evaporation or snow sublimation
!
! the model conserves surface energy fluxes as:
!    o -fsa + fira + fsh + (fcev+fctr+fgev) + fcst + fgr + fsm = 0
!    o fsa + fsr = [solad(1)+solad(2)+solai(1)+solai(2)]
!                = total incident solar
!    o fira = -firgcm + fire
! currently canopy heat storage fcst = 0
!
! ------------------------ code history ---------------------------
!    pot_evptr.F - From Inez
!    modified for LSM/CASA interface by J.John (2002)
! -----------------------------------------------------------------
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp      ! pft bounds
    integer, intent(in) :: num_soilp     ! number of soil points in pft filter
    integer, intent(in) :: filter_soilp(ubp-lbp+1) ! pft filter for soil points
!
! !LOCAL VARIABLES:
    integer  :: f                   ! filter index
    integer  :: g                   ! gridcell index
    integer  :: p                   ! pft index
    real(r8) :: qstar_net, q_grd, flh
    real(r8) :: tdegC, tadd, e_s, s, gamma, Q_E, factor, a_psy
    real(r8) :: fcst                        !canopy heat storage (w/m**2)

    ! inputs:

    integer , pointer :: pgridcell(:)      ! gridcell index of corresponding pft
    real(r8), pointer :: t_ref2m(:)        !2m surface air temperature (K)
    real(r8), pointer :: forc_pbot(:)      !atmospheric pressure (Pa)
    real(r8), pointer :: fsa(:)            !absorbed solar radiation (w/m**2)
    real(r8), pointer :: eflx_lwrad_net(:) !net infrared (longwave) rad (w/m**2) [+ = to atm]
    real(r8), pointer :: eflx_sh_tot(:)    !sensible heat flux (w/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_vege(:)   !veg evaporation heat flux (w/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_grnd(:)   !ground evaporation heat flux (w/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_vegt(:)   !veg transpiration heat flux (w/m**2) [+ to atm]

    !  outputs:

    real(r8), pointer :: pet(:)            !potential evaporation (mm h2o/s)
!
! !CALLED FROM:
! Casa in CASAMod
!
! !REVISION HISTORY:
! 2004.06.08 Vectorized and reformatted by Forrest Hoffman
!
!EOP
!-----------------------------------------------------------------------

    ! implicit intent in
    !============================================================
    pgridcell      => clm3%g%l%c%p%gridcell
    t_ref2m        => clm3%g%l%c%p%pes%t_ref2m 
    eflx_lwrad_net => clm3%g%l%c%p%pef%eflx_lwrad_net  
    eflx_sh_tot    => clm3%g%l%c%p%pef%eflx_sh_tot  
    eflx_lh_vege   => clm3%g%l%c%p%pef%eflx_lh_vege    
    eflx_lh_vegt   => clm3%g%l%c%p%pef%eflx_lh_vegt    
    eflx_lh_grnd   => clm3%g%l%c%p%pef%eflx_lh_grnd    
    fsa            => clm3%g%l%c%p%pef%fsa
    forc_pbot      => clm_a2l%forc_pbot

    ! implicit intent inout
    !============================================================
    pet            => clm3%g%l%c%p%pps%pet

    ! ----------------------------------------------------------------------

    ! reminder:  watt = joule/sec 
    !      latent heat of vaporization = 2.5e6  ! J/kg
    !      density h2o                 = 1.e3   ! kg/m3

    factor = hvap * denh2o * 1.e-3_r8  ! to convert to mm

    ! gamma:  phychromatic constant 
    !      a_psy = 0.000662   for ventilated, u~5 m/s
    !      a_psy = 0.000800 for naturally ventilated, u~1 m/s             

    a_psy = 0.000662_r8                ! unit is s-1

    do f = 1,num_soilp
       p = filter_soilp(f)
       g = pgridcell(p)

       ! set fcst to zero
       fcst = 0.0_r8

       ! net radiation abs by canopy = SW_net - LW_net 
       ! CHECK SIGN CONVENTION

       qstar_net = fsa(p) - eflx_lwrad_net(p)


       ! soil conductive heat flux (W/m2)
       ! CHECK SIGN CONVENTION

       flh  = eflx_lh_grnd(p)+eflx_lh_vege(p)+eflx_lh_vegt(p)  ! latent heat
       q_grd = qstar_net - eflx_sh_tot(p) - flh  ! compare to fgr

       !.......................................................................

       tdegC = t_ref2m(p) - tfrz     ! Kelvin to deg C
       tadd  = tdegC + 237.3_r8     

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
       ! e_s(T) is Clausius Clapeyron relationship
       !       = 0.6108 exp[ (17.27 T )/ (T+237.3) ]  (T is in deg C)
       !
       !**** cross-check:  T=0C, e_s=6.11mb
       !****               T=20C, e_s~20 mb
       !**** see also Andrews page 37, Holton page 484 etc.
       ! unit of e_s is kPa
       !                 
       ! at T= 0C, e_s(T) = 0.6108 kPa = 6.108 mb
       ! at T=20C, e_s(T) = 2.338  kPa = 23.38 mb
       !
       !   s  is slope of Clausius Clapeyron relationship = d/dT (e_s(T)) 
       !    unit is kPa/deg C
       !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

       e_s = 0.6108_r8 * exp(17.27_r8*tdegC/tadd)    ! kPa
       s = 4098._r8 * e_s / (tadd*tadd)           ! kPa/deg C
       ! 
       ! gamma:  phychromatic constant 
       !      a_psy = 0.000662   for ventilated, u~5 m/s
       !      a_psy = 0.000800 for naturally ventilated, u~1 m/s             

       gamma = a_psy * (forc_pbot(g)*1.e-3_r8) ! surface pressure converted to kPa

       !
       ! latent heat of evaporation (in same units W/m2 as energy fluxes) 
       !  According to Baldocchi reference, (s/s+gamma) is a strong function of
       !  temperature.  = 0.32 at -5C and =0.47 at +5C.
       !
       Q_E = (s/(s+gamma))*(qstar_net - q_grd - fcst)

       ! Potential evapotranspiration
       ! Pierre Friedlingstein's allocation code expects PET in (mm h20/s)
       ! Check units
       !   Q_E is in W/m2 = J/s * 1/m2     
       !   factor has units of 1.e-3 * J/m3 

       pet(p) = Q_E/factor         ! mm/sec

    end do

  end subroutine CASAPot_Evptr

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: casa_npp
!
! !INTERFACE:
  subroutine casa_npp(lbp, ubp, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Compute NPP from GPP.
! ------------------------ code history ---------------------------
!    casa_npp.f - get NPP from GPP
!    modified for LSM/CASA interface by J.John (2001)
! -----------------------------------------------------------------
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp   ! pft bounds
    integer, intent(in) :: num_soilp  ! number of soil points in pft filter
    integer, intent(in) :: filter_soilp(ubp-lbp+1) ! pft filter for soil points
!
! !LOCAL VARIABLES:
    integer :: p                      ! pft index
    integer :: f                      ! filter index

    ! pointers
    real(r8), pointer :: fpsn(:)
    real(r8), pointer :: lgrow(:)     ! growing season index (0 or 1)
    real(r8), pointer :: fnpp(:)      ! NPP (g/m2/sec)
!
! !CALLED FROM:
! Casa in CASAMod
!
! !REVISION HISTORY:
! 2004.06.08 Vectorized and reformatted by Forrest Hoffman
!
!EOP
!-----------------------------------------------------------------------

    ! implicit intent in
    !============================================================
    fpsn      => clm3%g%l%c%p%pcf%fpsn
    lgrow     => clm3%g%l%c%p%pps%lgrow

    ! implicit intent out
    !============================================================
    fnpp      => clm3%g%l%c%p%pps%fnpp

    do f = 1,num_soilp
       p = filter_soilp(f)

       !.. LSM overestimates FPSN: use scaled GPP to get npp 
       !.. LNPP = 1 - no lgrow

       fnpp(p) = gppfact * fpsn(p) * 12.0_r8 * 1.e-6_r8   ! umolC/m2/sec to gC/m2/sec

       !.. LNPP = 2 - use lgrow and scaled GPP to get npp

       if (LNPP == 2) then
          if (lgrow(p) /= 1.0_r8) then
             fnpp(p) = 0.0_r8
          end if
       end if

    end do


  end subroutine casa_npp

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: casa_allocate
!
! !INTERFACE:
  subroutine casa_allocate(lbp, ubp, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Provides dynamic or fixed carbon allocation.
! ------------------------ code history ---------------------------
!    allocate.f - Allocation Sub-model                            
!    modified from allocate.c - Allocation Sub-model             
!                                                                
!    Version 2.1                                                 
!                                                                
!    Designed by:             Pierre Freidlingstein             
!    Implemented: 12-12-97    Greg Asner                         
!    Modified: 9-14-98        Greg Asner                         
!    Modified: 08-31-00       Jeff Hicke                         
!    modified for LSM/CASA interface by J.John (2001)
!                                                                
! -----------------------------------------------------------------
!
! ------------------------ notes ----------------------------------
!
! This code allows for either fixed or dynamic allocation using a flag
! Need to check units for all variables (CASA vs LSM)
!
! code only executed for soils (ist = 1)
!
! -----------------------------------------------------------------
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp   ! pft bounds
    integer, intent(in) :: num_soilp  ! number of soil points in pft filter
    integer, intent(in) :: filter_soilp(ubp-lbp+1) ! pft filter for soil points
!
! !LOCAL VARIABLES:
    integer f,c,g,j,l,p

    real(r8) livesum

    real(r8) leaf_fract
    real(r8) wood_fract
    real(r8) root_fract
    real(r8) tfact
    real(r8) pfact
    real(r8) Llim
    real(r8) Nutrient
    real(r8) WorN

    real(r8), parameter :: fixed_leaf = 0.33333333_r8 
    real(r8), parameter :: fixed_stem = 0.33333333_r8 
    real(r8), parameter :: fixed_root = 0.33333333_r8 

    real(r8), parameter :: S0        = 0.30_r8
    real(r8), parameter :: R0        = 0.30_r8
    real(r8), parameter :: Llim_min  = 0.1_r8  ! Light limitation min value
    real(r8), parameter :: Llim_max  = 1.0_r8  ! Light limitation max value
    real(r8), parameter :: Nut_min   = 0.1_r8  ! Nutrient limitation min value
    real(r8), parameter :: Nut_max   = 1.0_r8  ! Nutrient limitation max value
    real(r8), parameter :: pfact_min = 0.5_r8  ! min value of precip
    real(r8), parameter :: pfact_mid = 1.0_r8  ! mid value of precip
    real(r8), parameter :: pfact_max = 2.0_r8  ! max value of precip
    real(r8), parameter :: tfact_min = 0.5_r8  ! min value of temperature
    real(r8), parameter :: tfact_max = 1.0_r8  ! max value of temperature
    real(r8), parameter :: Wlim_min  = 0.1_r8  ! Water limitation min value
    real(r8), parameter :: Wlim_max  = 1.0_r8  ! Water limitation max value

    ! ------------------------ input/output variables -----------------
    ! input

    integer , pointer :: pgridcell(:)    !gridcell index of corresponding pft
    integer , pointer :: pcolumn(:)      !pft's column
    integer , pointer :: ivt(:)          !pft vegetation type

    real(r8), pointer :: btran(:)        ! transpiration factor (0 to 1)
    real(r8), pointer :: tlai(:)         ! leaf area index, one-sided, unadjust for burying by snow
    real(r8), pointer :: forc_rain(:)    ! convective precipitation (mm h2o /s)
    real(r8), pointer :: forc_snow(:)    ! large-scale precipitation (mm h2o /s)

    real(r8), pointer :: z(:,:)          ! soil layer depth (m)
    real(r8), pointer :: dz(:,:)         ! soil layer thickness (m)
    real(r8), pointer :: sz(:)           !thickness of soil layers contributing to output
    real(r8), pointer :: szc(:)          !thickness of soil layers contributing to output
    real(r8), pointer :: t_soisno(:,:)   ! soil temperature (K)
    real(r8), pointer :: h2osoi_liq(:,:) ! liquid water (kg/m2)

    real(r8), pointer :: livefr(:,:)     !live fraction
    real(r8), pointer :: soilt(:)        !soil temp for top 30cm
    real(r8), pointer :: smoist(:)       !soil moisture for top 30cm
    real(r8), pointer :: soiltc(:)       !soil temp for entire column
    real(r8), pointer :: smoistc(:)      !soil moisture for entire column
    real(r8), pointer :: watoptc(:)      !optimal soil water content for et for entire column (mm3/mm3)
    real(r8), pointer :: watdryc(:)      !soil water when et stops for entire column (mm3/mm3)
    real(r8), pointer :: Wlim(:)
!
! !CALLED FROM:
! Casa in CASAMod
!
! !REVISION HISTORY:
! 2004.06.08 Vectorized and reformatted by Forrest Hoffman
!
!EOP
!-----------------------------------------------------------------------

    ! implicit intent in
    !============================================================

    pgridcell    => clm3%g%l%c%p%gridcell
    pcolumn      => clm3%g%l%c%p%column
    ivt          => clm3%g%l%c%p%itype
    dz           => clm3%g%l%c%cps%dz
    z            => clm3%g%l%c%cps%z
    sz           => clm3%g%l%c%p%pps%sz 
    szc          => clm3%g%l%c%p%pps%szc
    h2osoi_liq   => clm3%g%l%c%cws%h2osoi_liq
    t_soisno     => clm3%g%l%c%ces%t_soisno
    forc_rain    => clm_a2l%forc_rain
    forc_snow    => clm_a2l%forc_snow
    btran        => clm3%g%l%c%p%pps%btran 
    tlai         => clm3%g%l%c%p%pps%tlai
    livefr       => clm3%g%l%c%p%pps%livefr
    soilt        => clm3%g%l%c%p%pps%soilt 
    smoist       => clm3%g%l%c%p%pps%smoist 
    soiltc       => clm3%g%l%c%p%pps%soiltc
    smoistc      => clm3%g%l%c%p%pps%smoistc
    watoptc      => clm3%g%l%c%p%pps%watoptc   
    watdryc      => clm3%g%l%c%p%pps%watdryc
    Wlim         => clm3%g%l%c%p%pps%Wlim

    !============================================================

    ! Get avg soil moisture, avg soil temperature over all layers
    ! only for soils (ist = 1)

    do f = 1,num_soilp
       p = filter_soilp(f)
       soilt(p)   = 0._r8
       smoist(p)  = 0._r8
       soiltc(p)  = 0._r8
       smoistc(p) = 0._r8
    end do

!! convert soil temperature to deg C 

    do j = 1, nlevgrnd
       do f = 1,num_soilp
          p = filter_soilp(f)
          c = pcolumn(p)

!! top 30 cm
          if (z(c,j)+0.5_r8*dz(c,j) <= z30) then
             soilt(p)  = soilt(p) + (t_soisno(c,j)-tfrz)*dz(c,j)    
          end if

!! entire column
          soiltc(p)  = soiltc(p) + (t_soisno(c,j)-tfrz)*dz(c,j)    

       end do
    end do

!! convert liquid water to m3/m3 (need mm3/mm3 to match watsat, watdry, watopt)
!! h2osoi_liq /(denh2o*dz) gives m3/m3

    do j = 1, nlevsoi
       do f = 1,num_soilp
          p = filter_soilp(f)
          c = pcolumn(p)

!! top 30 cm
          if (z(c,j)+0.5_r8*dz(c,j) <= z30) then
             smoistc(p) = smoistc(p) + h2osoi_liq(c,j)/denh2o    ! dz cancels
          end if
!! entire column
          smoistc(p) = smoistc(p) + h2osoi_liq(c,j)/denh2o    ! dz cancels

       end do
    end do

    do f = 1,num_soilp
       p = filter_soilp(f)

       soilt(p)   = soilt(p)/sz(p)
       smoist(p)  = smoist(p)/sz(p)
       soiltc(p)  = soiltc(p)/szc(p)
       smoistc(p) = smoistc(p)/szc(p)
    end do

    ! -----------------------------------------------------------------
    ! dynamic allocation; Pierre's model (LALLOC=1)
    ! -----------------------------------------------------------------

    if (LALLOC == 1) then

       do f = 1,num_soilp
          p = filter_soilp(f)
          g = pgridcell(p)

          !  Light limitation calculation 
          !  ----------------------------

          Llim = exp(-0.5_r8 * tlai(p))
          if (Llim <= Llim_min) Llim = Llim_min
          if (Llim >= Llim_max) Llim = Llim_max

          !  Pseudo-nutrient limitation calculation 
          !  --------------------------------------
          !  set pfact = bevap (entire column)

          if (soiltc(p) > 0.0_r8) then
             pfact = min( max(smoistc(p)-watdryc(p),0._r8) / &
                       (watoptc(p)-watdryc(p)), 1._r8 )
          else
             pfact = 0.01_r8
          end if

          tfact = Q10**((soilt(p)-30.0_r8)/10.0_r8)    ! why not soiltc ?
          if (tfact >= tfact_max) tfact = tfact_max
          if (tfact <= tfact_min) tfact = tfact_min

          Nutrient    = pfact*tfact
          if (Nutrient  <= Nut_min) Nutrient    = Nut_min
          if (Nutrient  >= Nut_max) Nutrient    = Nut_max

          !  Water limitation calculation 
          !  ----------------------------
          !  set Wlim = btran (entire column)

         Wlim(p)=btran(p)
         if (Wlim(p)    <= Wlim_min) Wlim(p)    = Wlim_min
         if (Wlim(p)    >= Wlim_max) Wlim(p)    = Wlim_max

         WorN = min(Wlim(p),Nutrient)

!! determine fraction allocated to live pools
         if (lnonwood(ivt(p)) == 0) then
            livefr(p,FROOT) = R0 * 3.0_r8 * Llim/(Llim+2.0_r8*WorN)
            livefr(p,WOOD)  = S0 * 3.0_r8 * WorN/(2.0_r8*Llim+WorN)
            livefr(p,LEAF)  = 1.0_r8 - livefr(p,FROOT) - livefr(p,WOOD)
         else
            livefr(p,FROOT) = R0 * 3.0_r8 * Llim/(Llim+2.0_r8*WorN)
            livefr(p,WOOD)  = 0._r8
            livefr(p,LEAF)  = 1.0_r8 - livefr(p,FROOT)
         end if

       end do  ! end do loop

    end if  ! end dynamic allocation

    ! -----------------------------------------------------------------
    ! fixed allocation
    ! this determines allocation based on vegetation type - if there is
    ! wood in the vegetation class, 1/3 of NPP is given to wood litter -
    ! the remainder is divided evenly between leaf and root material
    ! -----------------------------------------------------------------

    if (LALLOC == 0) then

       do f = 1,num_soilp
          p = filter_soilp(f)

          if (lnonwood(ivt(p)) == 0) then
            wood_fract = fixed_stem      ! (1/3)
         else
            wood_fract = 0._r8
         end if
         leaf_fract = (1.0_r8 - wood_fract)/2.0_r8
         root_fract = (1.0_r8 - wood_fract)/2.0_r8
         livefr(p,LEAF)  = leaf_fract
         livefr(p,WOOD)  = wood_fract
         livefr(p,FROOT) = root_fract
       end do

    end if

  end subroutine casa_allocate

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: casa_bgfluxes
!
! !INTERFACE:
  subroutine casa_bgfluxes(lbp, ubp, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Compute biogeochemical fluxes.
!
! ------------------------ code history ---------------------------
!    bgfluxes_BASIC.c - BG Spin-Up Flux Sub-model VERSION 1.0     
!                                                                
!    Version 2.1                                                
!                                                              
!    Created 7-13-99 by Greg Asner                            
!    Modified for CASA2a on 8-14-00 by Greg Asner            
!    Gleaned from CASA2b                                    
!    modified for LSM/CASA interface by J.John (2001)
!iyf modified 7-18-01 Inez:  all modifications marked by iyf
!                                                          
! -----------------------------------------------------------------
! code only executed for soils (ist = 1)
!
! !USES:
    use clm_time_manager , only : get_step_size
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp   ! pft bounds
    integer, intent(in) :: num_soilp  ! number of soil points in pft filter
    integer, intent(in) :: filter_soilp(ubp-lbp+1) ! pft filter for soil points
!
! !LOCAL VARIABLES:
    integer  f,g,i,j,l,n,p
    integer iptype
    real(r8) leafmass
    real(r8) dtime                  !land model time step (sec)

    ! ------------------------ input/output variables -----------------
    ! input

    integer , pointer :: ivt(:)       !pft vegetation type

    real(r8), pointer :: fnpp(:)      !NPP  (gC/m2/sec)
    real(r8), pointer :: excessC(:)   !excess Carbon (gC/m2/timestep)
    real(r8), pointer :: bgtemp(:)    !temperature dependence for C pools
    real(r8), pointer :: bgmoist(:)   !moisture dependence for C pools
    real(r8), pointer :: soilt(:)     !soil temp for top 30cm
    real(r8), pointer :: smoist(:)    !soil moisture for top 30cm
    real(r8), pointer :: watopt(:)
    real(r8), pointer :: watdry(:)    
    real(r8), pointer :: Wlim(:)
    real(r8), pointer :: livefr(:,:)  !live fraction
    real(r8), pointer :: sandfrac(:)

    ! implicit intent inout
    !============================================================
    real(r8), pointer :: plai(:)     ! prognostic LAI (m2 leaf/m2 ground)
    real(r8), pointer :: Closs(:,:)
    real(r8), pointer :: Ctrans(:,:) ! C transfers out of pool types
    real(r8), pointer :: Resp_C(:,:) ! could dimension by ndead, but caution!!!
    real(r8), pointer :: Tpool_C(:,:)
    real(r8), pointer :: Cflux(:)
    real(r8), pointer :: XSCpool(:)
    real(r8), pointer :: co2flux(:)  !net CO2 flux (gC/m2/s) [+ = to atm]
!
! !CALLED FROM:
! Casa in CASAMod
!
! !REVISION HISTORY:
! 2004.06.08 Vectorized and reformatted by Forrest Hoffman
!
!EOP
!-----------------------------------------------------------------------

    ivt        => clm3%g%l%c%p%itype
    fnpp       => clm3%g%l%c%p%pps%fnpp   
    bgtemp     => clm3%g%l%c%p%pps%bgtemp
    bgmoist    => clm3%g%l%c%p%pps%bgmoist
    excessC    => clm3%g%l%c%p%pps%excessC   
    plai       => clm3%g%l%c%p%pps%plai
    Closs      => clm3%g%l%c%p%pps%Closs
    Ctrans     => clm3%g%l%c%p%pps%Ctrans
    Resp_C     => clm3%g%l%c%p%pps%Resp_C
    Tpool_C    => clm3%g%l%c%p%pps%Tpool_C
    XSCpool    => clm3%g%l%c%p%pps%XSCpool
    Cflux      => clm3%g%l%c%p%pps%Cflux  
    soilt      => clm3%g%l%c%p%pps%soilt 
    smoist     => clm3%g%l%c%p%pps%smoist
    watdry     => clm3%g%l%c%p%pps%watdry
    watopt     => clm3%g%l%c%p%pps%watopt
    Wlim       => clm3%g%l%c%p%pps%Wlim       
    co2flux    => clm3%g%l%c%p%pps%co2flux
    livefr     => clm3%g%l%c%p%pps%livefr
    sandfrac   => clm3%g%l%c%p%pps%sandfrac

    ! -----------------------------------------------------------------

    ! Get step size

    dtime = get_step_size()

    !---Step 1:  growth and allocation -----------------------------

    ! INCREMENT PLANT CARBON POOLS  (Initialized in casainit.F)
    ! Allocation (livefr) may or may not depend on climate
    ! Tpool is unit of gC/m2

    do n = 1, nlive
       do f = 1,num_soilp
          p = filter_soilp(f)
          Tpool_C(p,n) = Tpool_C(p,n) + livefr(p,n) * fnpp(p) * dtime
       end do
    end do

    do f = 1,num_soilp
       p = filter_soilp(f)

       !iyf *****
       !iyf 02/07/22
       !iyf Keep track of fluxes
       !iyf Transfer excess leaf carbon to root
       !jgj 02/09/19
       !>> root pool is getting too big - don't dump excess carbon into root pool
       !>> save excess in separate (phothosynthate) pool for diagnostic purposes
       !>> dump excess into autotrophic respiration  (see below: litterfall)

       excessC(p) = max(Tpool_C(p,LEAF)-leafmax(ivt(p)), 0.0_r8)
       Tpool_C(p,LEAF) = min(Tpool_C(p,LEAF), leafmax(ivt(p)))
       XSCpool(p) = XSCpool(p) + excessC(p)

       !iyf - end1 *****

    end do

    !---Step 1a:  determine prognostic LAI --------------------------------

    ! Leaf mass (gC/m2 ground) to LAI (m2 leaf/m2 ground)
    ! use Specific Leaf Area (m2 leaf/kgC) from Dickinson et al.
    !     SLA is a function of veg type
    ! Tpool is unit of gC/m2

    do f = 1,num_soilp
       p = filter_soilp(f)

       leafmass  = Tpool_C(p,LEAF)*1.e-3_r8 ! (kg C/m2 ground)
       plai(p)   = leafmass * sla(ivt(p))! m2 leaf/m2 ground

       !iyf    upper limit is placed here on leaf mass in this subroutine above
       !iyf    lower limit is placed on CLoss in casa_litterfall.F

       if(plai(p) >= plai_max(ivt(p)))plai(p) = plai_max(ivt(p))
       if(plai(p) <= plai_min(ivt(p)))plai(p) = plai_min(ivt(p))

    end do

    !---Step 2:  death and litterfall -------------------------------
    ! COMPUTE FLUXES FROM LITTERFALL
    !iyf: these use "annK" for the three live pools

    ! 03/03/11 dump excessC into CLOSS101 (g/m2/timestep)

    call casa_litscl(lbp, ubp)
    call casa_litterfall(lbp, ubp, num_soilp, filter_soilp)


    !---Step 3:  respiration ---------------------------------------
    !iyf this should be heterotrophic respiration only
    !.. Initialize respiration fluxes each timestep
    !.. Note: these should be over dead pools only (see resp_pool_index)
    !iyf:  Resp in unit of gC/m2/timestep

    do n = 1, npools
       do f = 1,num_soilp
          p = filter_soilp(f)
          Resp_C(p,n) = 0._r8 ! could dimension by ndead, but caution later!!!
       end do
    end do

    !---Step 3a: TEMPERATURE AND MOISTURE CONSTRAINTS ON DECOMP 
    !  bgmoist currently set to 0.5 for all points
    !  01/07/04 change bgmoist to CLM2/LPJ temperature dependence

    do f = 1,num_soilp
       p = filter_soilp(f)

       ! temperature dependence
       bgtemp(p) = (Q10 ** ((soilt(p) - 30.0_r8) / 10.0_r8))

       !...................................................................
       ! moisture dependence 
       ! this is like Water limitation calculation in casa_allocate.F
       !       moist_resp = 0.25 + 0.75 * clm%wf
       !
       ! there are different parametrizations of moisture dependence
       ! eg bgmoist = 0.5
       ! eg: Pierre Friedlingstein's moisture dependence - fn(Wlim)
       !...................................................................

       ! mod 02/07/17 go back to linear calculation of bgmoist
       ! mimic calculation of bevap in surphy.F to get Wlim 
       ! but use smoist,soilt instead of h2osoi,tsoi and watsat instead of watopt
       !   watdry = water content when evapotranspiration stops = wp
       !   watsat = volumetric soil water content, saturation (porosity) = fcap

       if (soilt(p) > 0.0_r8) then
          Wlim(p) = min( max(smoist(p)-watdry(p),0._r8) /  &
               (watopt(p)-watdry(p)), 1._r8 )
       else
          Wlim(p) = 0.01_r8
       end if

       bgmoist(p) = 0.25_r8 + 0.75_r8*Wlim(p)

       !...................................................................

    end do

    !-------------------------------------------------------------------------

    !---Step 3b: DETERMINE loss of C FROM EACH DEAD POOL (donor) PER TIMESTEP
    !iyf:  Closs is the amount of carbon each pool loses in gC/m2/timestep.
    !iyf:  A fraction of Closs is transferred to another pool (receiver pool), 
    !iyf:  the remainder to the atm (resp).
    !iyf:  The distribution of Closs is done in subroutine casa_respire, Step 3c. 

    do n = nlive+1, npools
       do f = 1,num_soilp
          p = filter_soilp(f)
          Closs(p,n) = Tpool_C(p,n) * kdt(ivt(p),n) * bgtemp(p) * bgmoist(p)
       end do
    end do

    do f = 1,num_soilp
       p = filter_soilp(f)

       !* adjust pools
       Closs(p,SURFSTR) = Closs(p,SURFSTR) * lignineffect(ivt(p))
       Closs(p,SOILSTR) = Closs(p,SOILSTR) * lignineffect(ivt(p))
       Closs(p,SOILMIC) = Closs(p,SOILMIC) &
            * (1.0_r8 - 0.75_r8*(1.0_r8 - sandfrac(p))) * fact_soilmic(ivt(p))
       Closs(p,SLOW) = Closs(p,SLOW) * fact_slow(ivt(p))
       Closs(p,PASSIVE) = Closs(p,PASSIVE) * fact_passive(ivt(p))
    end do

    !iyf:  limits on loss from dead pools.
    do n = nlive+1,npools
       do f = 1,num_soilp
          p = filter_soilp(f)

          !iyf - replace IF test with MIN/MAX functions.
          !iyf - No need to track limits on loss rate
          !iyf (only need to track limits on inventories)

          Closs(p,n)=min(Closs(p,n),Tpool_C(p,n))

      end do
    end do

    !---step 3c:  SOM C DECOMPOSITION
    !iyf:  update Tpool's:  C is transferred between donor and receiver pools 
    !iyf:  Resp is the amount of C to atm, in gC/m2/timestep

    call casa_respire(lbp, ubp, num_soilp, filter_soilp)

    !---Step 4-----------------------------------------------------------
    ! CALCULATE NITROGEN POOLS                                 

    !   do n = 1,npools
    !      do f = 1,num_soilp
    !         p = filter_soilp(f)
    !         Tpool_N(p,n)=Tpool_C(p,n)/CNratio(n)
    !      end do   
    !   end do   

    !--- Step 5 --------------------------------------------------------
    !  Get Total C Fluxes to atm = Sum over respiring (dead) pools/dtlsm
    !  Note: Resp for live pools should be zero
    !iyf:  Cflux in gC/m2/s

    Cflux(lbp:ubp)=0._r8

    do n = nlive+1,npools
       do f = 1,num_soilp
          p = filter_soilp(f)

          Resp_C(p,n)=Resp_C(p,n)/dtime ! g/m2/s instead of g/m2/timestep
          Cflux(p)=Cflux(p)+Resp_C(p,n)

       end do
    end do

    do n = 1,npools
       do f = 1,num_soilp
          p = filter_soilp(f)

          Closs(p,n)=Closs(p,n)/dtime
          !            Nloss(p,n)=Nloss(p,n)/dtime  ! not used

       end do
    end do

    ! Loop over all pool types to adjust units on Ctrans to gC/m2/s
    do iptype = 1, npool_types
       do f = 1,num_soilp
          p = filter_soilp(f)
          Ctrans(p,iptype) = Ctrans(p,iptype) / dtime
       end do
    end do

    do f = 1,num_soilp
       p = filter_soilp(f)

       ! compute co2flux
       co2flux(p) = (-fnpp(p) + Cflux(p))

    end do

  end subroutine casa_bgfluxes

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: casa_litscl
!
! !INTERFACE:
  subroutine casa_litscl(lbp, ubp)
!
! !DESCRIPTION:
! Compute litter fall scalars
!
! ------------------------ code history ---------------------------
!    casa_litscl.F - Litterfall scalars 
!    modified for LSM/CASA interface by J.John (2001)
!                                                                   
! -----------------------------------------------------------------
! Notes:
!
! These need to be changed (dependent on monthly LAI) 
!
! code only executed for soils (ist = 1)
!                                                                     
! -----------------------------------------------------------------
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                    ! pft bounds

!
! !LOCAL VARIABLES:
    integer  l,p

    ! ------------------------ input/output variables -----------------
    !  inputs:  

    integer , pointer :: ltype(:)     ! landunit type for corresponding pft
    integer , pointer :: plandunit(:) ! landunit index associated with each pft
    real(r8), pointer :: litterscalar(:)
    real(r8), pointer :: rootlitscalar(:)
!
! !CALLED FROM:
! casa_bgfluxes in CASAMod
!
! !REVISION HISTORY:
! 2004.06.08 Vectorized and reformatted by Forrest Hoffman
!
!EOP
!-----------------------------------------------------------------------


    ltype         => clm3%g%l%itype
    plandunit     => clm3%g%l%c%p%landunit
    litterscalar  => clm3%g%l%c%p%pps%litterscalar
    rootlitscalar => clm3%g%l%c%p%pps%rootlitscalar


    ! These need to be changed (dependent on monthly LAI) 

    do p = lbp, ubp
       l = plandunit(p)
       if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
          litterscalar(p) = 1.0_r8
          rootlitscalar(p) = 1.0_r8
       else
          litterscalar(p) = 0.0_r8
          rootlitscalar(p) = 0.0_r8
       end if
    end do

  end subroutine casa_litscl

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: casa_litterfall
!
! !INTERFACE:
  subroutine casa_litterfall(lbp, ubp, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Compute dynamic litter fall.
!
! ------------------------ code history ---------------------------
!    litterfall_BASIC.c - Litterfall Spin-up Sub-model VERSION 1.0   
!                                                                   
!    Version 2.1                                                      
!                                                                     
!    Implemented by:    4-22-99     Greg Asner                        
!    Dynamic litterfall by Jim Randerson
!    modified for LSM/CASA interface by J.John (2001)
!                                                                     
! -----------------------------------------------------------------
!
! This routine is used to determine the timing of litterfall
! code only executed for soils (ist = 1)
!
! !USES:
    use clm_time_manager , only : get_step_size, get_nstep
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp   ! pft bounds
    integer, intent(in) :: num_soilp  ! number of soil points in pft filter
    integer, intent(in) :: filter_soilp(ubp-lbp+1) ! pft filter for soil points
!
! !LOCAL VARIABLES:
    integer  f,l,p,n
    integer  nstep
    real(r8) dtime                   !land model time step (sec)

    ! ------------------------ input/output variables -----------------
    !  inputs:  

    integer , pointer :: ivt(:)       !pft vegetation type

    real(r8), pointer :: stressCD(:)  ! cold and drought stress function (sec-1)
    ! add to "annK(m,LEAF)" and "annK(m,FROOT)"
    ! in casa_litterfall.F
    real(r8), pointer :: excessC(:)   ! excess Carbon (gC/m2/timestep)
    real(r8), pointer :: Closs(:,:)   ! C lost to atm
    real(r8), pointer :: Tpool_C(:,:) ! Total C pool size
    real(r8), pointer :: litterscalar(:)
    real(r8), pointer :: rootlitscalar(:)
!
! !CALLED FROM:
! casa_bgfluxes in CASAMod
!
! !REVISION HISTORY:
! 2004.06.08 Vectorized and reformatted by Forrest Hoffman
!
!EOP
!-----------------------------------------------------------------------

    ivt          => clm3%g%l%c%p%itype
    stressCD     => clm3%g%l%c%p%pps%stressCD 
    excessC      => clm3%g%l%c%p%pps%excessC  
    Closs        => clm3%g%l%c%p%pps%Closs
    Tpool_C      => clm3%g%l%c%p%pps%Tpool_C
    litterscalar => clm3%g%l%c%p%pps%litterscalar
    rootlitscalar => clm3%g%l%c%p%pps%rootlitscalar

    ! ----------------------------------------------------------
    ! Get step size

    dtime = get_step_size()
    nstep = get_nstep()

    !  FOLIAGE and ROOT LOSS; Jim's model; Needs to be checked

    do f = 1,num_soilp
       p = filter_soilp(f)

       ! Carbon lost as litter per timestep
       !    dL/dt = -L/tau
       !    Closs = L(n) - L(n+1) = M*delta_t/tau

       !iyf 02/07/11
       ! stressCD to be added to "annK(m,LEAF)" and "annK(m,FROOT)"
       ! iyf 02/07/17 apply stressCD to leaves only (not to roots)
       ! in casa_litterfall.F
       ! kdt is used for WOOD and dead pools only (ie no stressCD needed)

       Closs(p,LEAF) = Tpool_C(p,LEAF)   &
            * ((annK(ivt(p),LEAF)+stressCD(p))*dtime) * litterscalar(p)
       Closs(p,WOOD) = Tpool_C(p,WOOD)   &
            * kdt(ivt(p),WOOD)
       Closs(p,FROOT) = Tpool_C(p,FROOT)   &
            * ((annK(ivt(p),FROOT)           )*dtime) * rootlitscalar(p)

       !.. mod IYF 02/07/22

       !iyf
       !iyf Maintain a minimum Leaf Mass.  Pretend LeafMin=photosynthate to start
       !iyf photosynthesis in the spring.
       !iyf
       If (Tpool_C(p,LEAF) > leafmin(ivt(p))) then
          Closs(p,LEAF) =   &
               MIN(Tpool_C(p,LEAF)-leafmin(ivt(p)),Closs(p,LEAF))
       else
          Closs(p,LEAF) = 0._r8
       end if

    end do

    ! Decrement plant C pools : LEAF, WOOD, FROOT

    do n = 1,nlive
       do f = 1,num_soilp
          p = filter_soilp(f)
          Tpool_C(p,n)  =  Tpool_C(p,n)-Closs(p,n)    
       end do
    end do

    do f = 1,num_soilp
       p = filter_soilp(f)

       !iyf 03/03/25
       ! add excessC to CLOSS101 (g/m2/timestep)
       ! excessC is going from NPP to litter, without going through leaves.
       Closs(p,LEAF) = Closs(p,LEAF) + excessC(p)

       !Increment litter carbon pools 

       Tpool_C(p,CWD)     = Tpool_C(p,CWD)  &
            + Closs(p,WOOD)
       Tpool_C(p,SURFMET) = Tpool_C(p,SURFMET)  &
            + Closs(p,LEAF)  * solubfract(ivt(p))
       Tpool_C(p,SOILMET) = Tpool_C(p,SOILMET) &
            + Closs(p,FROOT) * solubfract(ivt(p))
       Tpool_C(p,SURFSTR) = Tpool_C(p,SURFSTR) &
            + Closs(p,LEAF)   &
            * (1._r8 - solubfract(ivt(p)))
       Tpool_C(p,SOILSTR) = Tpool_C(p,SOILSTR) &
            + Closs(p,FROOT)  &
            * (1._r8 - solubfract(ivt(p)))

    end do

  end subroutine casa_litterfall

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: casa_respire
!
! !INTERFACE:
  subroutine casa_respire(lbp, ubp, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Compute respiration.
!
! ------------------------ code history ---------------------------
!    respire_BASIC.c - Respiration Sub-model VERSION 1.0             
!                                                                    
!    Version 2.1                                                     
!                                                                    
!    Created 7-13-99 by Greg Asner                                   
!    Modified for CASA2a on 8-16-00 by Greg Asner                    
!    Gleaned from CASA2b                                             
!    modified for LSM/CASA interface by J.John (2001)
!                                                                    
! -----------------------------------------------------------------
!
! code only executed for soils (ist = 1)
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp   ! pft bounds
    integer, intent(in) :: num_soilp  ! number of soil points in pft filter
    integer, intent(in) :: filter_soilp(ubp-lbp+1) ! pft filter for soil points
!
! !LOCAL VARIABLES:
    ! ------------------------ input/output variables -----------------

    ! implicit intent in
    !============================================================

    real(r8), pointer :: Closs(:,:)   ! C lost to atm
    real(r8), pointer :: Ctrans(:,:)  ! C transfers out of pool types
    real(r8), pointer :: eff(:,:)
    real(r8), pointer :: frac_donor(:,:)

    ! implicit intent out
    !============================================================

    real(r8), pointer :: Resp_C(:,:)  ! 
    real(r8), pointer :: Tpool_C(:,:) ! Total C pool size

    ! ------------------------ local variables -----------------

    integer f,l,p,n
    integer irtype,iptype
    integer donor_pool
    integer recvr_pool
    integer donor_type
    integer recvr_type

    real(r8) Out

!
! !CALLED FROM:
! casa_bgfluxes in CASAMod
!
! !REVISION HISTORY:
! 2004.06.08 Vectorized and reformatted by Forrest Hoffman
!
!EOP
!-----------------------------------------------------------------------

    Closs      => clm3%g%l%c%p%pps%Closs
    Ctrans     => clm3%g%l%c%p%pps%Ctrans
    Resp_C     => clm3%g%l%c%p%pps%Resp_C
    Tpool_C    => clm3%g%l%c%p%pps%Tpool_C
    eff        => clm3%g%l%c%p%pps%eff
    frac_donor => clm3%g%l%c%p%pps%frac_donor

    ! Loop over all pool types to initialize transfer to zero
    do iptype = 1, npool_types
       do f = 1,num_soilp
          p = filter_soilp(f)
          Ctrans(p,iptype) = 0._r8
       end do
    end do

    ! Loop over all respiring pools
    do irtype = 1, nresp_pools

       donor_pool = resp_pool_index(1,irtype)
       recvr_pool = resp_pool_index(2,irtype)

       donor_type = pool_type_index(donor_pool)
       recvr_type = pool_type_index(recvr_pool)

       ! Loop over pfts
       do f = 1,num_soilp
          p = filter_soilp(f)

          Out  = Closs(p,donor_pool) * frac_donor(p,irtype)

          ! accumulate total pool transfers by pool type
          if (donor_type .ne. recvr_type) then
             Ctrans(p,donor_type) = Ctrans(p,donor_type) + (Out * eff(p,irtype))
          end if

          Tpool_C(p,donor_pool) = Tpool_C(p,donor_pool) &
               - Out
          Tpool_C(p,recvr_pool) = Tpool_C(p,recvr_pool)  &
               + (Out * eff(p,irtype))
          Resp_C(p,donor_pool) =  Resp_C(p,donor_pool) &
               + Out * (1._r8 - eff(p,irtype))

          ! make sure donor pool does not fall below zero
          if (Tpool_C(p,donor_pool) <= 0.0_r8) then
             Tpool_C(p,donor_pool) = 0.0_r8
          end if

       end do      ! end number of points
    end do         ! end number of pools

    ! Stuff the total C lost by live pools into live pool type transfers
    do n = 1, nlive
       do f = 1,num_soilp
          p = filter_soilp(f)
          Ctrans(p,pool_type_index(n)) = Ctrans(p,pool_type_index(n)) + &
                                         Closs(p,n)
       end do
    end do

    ! Add in the respired C to transfers out of pool types
    do n = nlive+1, npools
       do f = 1,num_soilp
          p = filter_soilp(f)
          Ctrans(p,pool_type_index(n)) = Ctrans(p,pool_type_index(n)) + &
                                         Resp_C(p,n)
       end do
    end do

  end subroutine casa_respire

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CASARest
!
! !INTERFACE:
  subroutine CASARest(ncid, flag)
!
! !DESCRIPTION:
! Read/Write CASA information to/from restart file.
!
! !USES:
    use clmtype
    use ncdio_pio
    use clm_time_manager, only : is_restart
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid  !netcdf id
    character(len=*) , intent(in)    :: flag  !flag='read, data= or 'write'
!
! !CALLED FROM:
! restart in restFileMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
! 2004.06.08 Vectorized and reformatted by Forrest Hoffman
!
! !LOCAL VARIABLES:
    integer :: p,j                        ! indices
    logical :: readvar                    ! determine if variable is on initial file
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
!EOP
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! pft type physical state variable - livefr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livefr', xtype=ncd_double,  &
            dim1name='pft', dim2name='nlive',&
            long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livefr', data=pptr%pps%livefr, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft type physical state variable - Tpool_C  (Tracer 1)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='Tpool_C', xtype=ncd_double,  &
            dim1name='pft', dim2name='npools',&
            long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='Tpool_C', data=pptr%pps%Tpool_C, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read') then
          if (.not. readvar) then
             if (is_restart()) then
                call endrun
             else
                if (masterproc) &
                   write(iulog,*)'WARNING: TPOOL_C not contained on initial/restart dataset'
             end if
          else
             if (masterproc) &
                write(iulog,*)'TPOOL_C read from initial/restart dataset'
             cpool_inic = .true.
          end if
       end if
    end if

    ! pft type physical state variable - lgrow
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='lgrow', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='lgrow', data=pptr%pps%lgrow, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft type physical state variable - iseabeg
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='iseabeg', xtype=ncd_double, &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='iseabeg', data=pptr%pps%iseabeg, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft type physical state variable - nstepbeg
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='nstepbeg', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='nstepbeg', data=pptr%pps%nstepbeg, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft type physical state variable - degday
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='degday', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='degday', data=pptr%pps%degday, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft type physical state variable - ndegday
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ndegday', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='ndegday', data=pptr%pps%ndegday, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft type physical state variable - tday
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tday', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='tday', data=pptr%pps%tday, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft type physical state variable - tcount
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tcount', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='tcount', data=pptr%pps%tcount, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft type physical state variable - tdayavg
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tdayavg', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='tdayavg', data=pptr%pps%tdayavg, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft type physical state variable - stressCD
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='stressCD', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='stressCD', data=pptr%pps%stressCD, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

!! added so history output is correct

    ! pft type physical state variable - stressT
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='stressT', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='stressT', data=pptr%pps%stressT, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft type physical state variable - stressW
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='stressW', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='stressW', data=pptr%pps%stressW, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

!! end addition

    ! pft type physical state variable - XSCpool
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='XSCpool', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='XSCpool', data=pptr%pps%XSCpool, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! eflx_lwrad_net => clm3%g%l%c%p%pef%eflx_lwrad_net  
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='eflx_lwrad_net', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='eflx_lwrad_net', data=pptr%pef%eflx_lwrad_net, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! eflx_lh_grnd => clm3%g%l%c%p%pef%eflx_lh_grnd  
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='eflx_lh_grnd', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='eflx_lh_grnd', data=pptr%pef%eflx_lh_grnd, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! eflx_lh_vege => clm3%g%l%c%p%pef%eflx_lh_vege
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='eflx_lh_vege', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='eflx_lh_vege', data=pptr%pef%eflx_lh_vege, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! eflx_lh_vegt => clm3%g%l%c%p%pef%eflx_lh_vegt
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='eflx_lh_vegt', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='eflx_lh_vegt', data=pptr%pef%eflx_lh_vegt, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

  end subroutine CASARest

!===============================================================================
!BOP
!
! !IROUTINE: CASASummary
!
! !INTERFACE:
subroutine CASASummary(lbp, ubp, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Perform pft carbon summary calculations
!
! !USES:
   use clm_time_manager, only: get_step_size
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbp, ubp       ! pft bounds
   integer, intent(in) :: num_soilp      ! number of soil points in pft filter
   integer, intent(in) :: filter_soilp(ubp-lbp+1) ! pft filter for soil points
!
! !CALLED FROM:
! subroutine casa_ecocystemDyn
!
! !REVISION HISTORY:
! 22 Sept 2006: Created by Forrest Hoffman
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
   real(r8), pointer :: fpsn(:)              !photosynthesis (umol CO2 /m**2 /s)
   real(r8), pointer :: Tpool_C(:,:)         ! Total C by pool
   real(r8), pointer :: Resp_C(:,:)          ! Respired C by pool
   real(r8), pointer :: Closs(:,:)           ! Lost C by pool
   real(r8), pointer :: Ctrans(:,:)          ! Transferred C by pool type
   real(r8), pointer :: Cflux(:)             ! C flux
   real(r8), pointer :: livefr(:,:)          ! Live fraction
   real(r8), pointer :: fnpp(:)              ! NPP (gC/m2/sec)
   real(r8), pointer :: co2flux(:)           ! net CO2 flux (gC/m2/s) [+= atm]
!
!
! local pointers to implicit in/out scalars
!
!
! local pointers to implicit out scalars
   real(r8), pointer :: casa_agnpp(:)        ! above-ground net primary production [gC/m2/s]
   real(r8), pointer :: casa_ar(:)           ! autotrophic respiration [gC/m2/s]
   real(r8), pointer :: casa_bgnpp(:)        ! below-ground net primary production [gC/m2/s]
   real(r8), pointer :: casa_cwdc(:)         ! coarse woody debris C [gC/m2]
   real(r8), pointer :: casa_cwdc_hr(:)      ! cwd heterotrophic respiration [gC/m2/s]
   real(r8), pointer :: casa_cwdc_loss(:)    ! cwd C loss [gC/m2/s]
   real(r8), pointer :: casa_frootc(:)       ! fine root C [gC/m2]
   real(r8), pointer :: casa_frootc_alloc(:) ! fine root C allocation [gC/m2/s]
   real(r8), pointer :: casa_frootc_loss(:)  ! fine root C loss [gC/m2/s]
   real(r8), pointer :: casa_gpp(:)          ! gross primary production [gC/m2/s]
   real(r8), pointer :: casa_hr(:)           ! total heterotrophic respiration [gC/m2/s]
   real(r8), pointer :: casa_leafc(:)        ! leaf C [gC/m2]
   real(r8), pointer :: casa_leafc_alloc(:)  ! leaf C allocation [gC/m2/s]
   real(r8), pointer :: casa_leafc_loss(:)   ! leaf C loss [gC/m2/s]
   real(r8), pointer :: casa_litterc(:)      ! total litter C (excluding cwd C) [gC/m2]
   real(r8), pointer :: casa_litterc_hr(:)   ! litter heterotrophic respiration [gC/m2/s]
   real(r8), pointer :: casa_litterc_loss(:) ! litter C loss [gC/m2/s]
   real(r8), pointer :: casa_nee(:)          ! net ecosystem exchange [gC/m2/s]
   real(r8), pointer :: casa_nep(:)          ! net ecosystem production [gC/m2/s]
   real(r8), pointer :: casa_npp(:)          ! net primary production [gC/m2/s]
   real(r8), pointer :: casa_soilc(:)        ! total soil organic matter C (excluding cwd and litter C) [gC/m2]
   real(r8), pointer :: casa_soilc_hr(:)     ! soil heterotrophic respiration [gC/m2/s]
   real(r8), pointer :: casa_soilc_loss(:)   ! total soil organic matter C loss [gC/m2/s]
   real(r8), pointer :: casa_woodc(:)        ! wood C [gC/m2]
   real(r8), pointer :: casa_woodc_alloc(:)  ! wood C allocation [gC/m2/s]
   real(r8), pointer :: casa_woodc_loss(:)   ! wood C loss [gC/m2/s]
!
!
! !OTHER LOCAL VARIABLES:
   integer :: f,p          ! indices
   real(r8):: dtime        ! land model time step (sec)

!EOP
!-----------------------------------------------------------------------

   ! assign local pointers at the pft level
   fpsn                           => clm3%g%l%c%p%pcf%fpsn
   Tpool_C                        => clm3%g%l%c%p%pps%Tpool_C
   Resp_C                         => clm3%g%l%c%p%pps%Resp_C
   Closs                          => clm3%g%l%c%p%pps%Closs
   Ctrans                         => clm3%g%l%c%p%pps%Ctrans
   Cflux                          => clm3%g%l%c%p%pps%Cflux  
   livefr                         => clm3%g%l%c%p%pps%livefr
   fnpp                           => clm3%g%l%c%p%pps%fnpp
   co2flux                        => clm3%g%l%c%p%pps%co2flux

   casa_agnpp                     => clm3%g%l%c%p%pps%casa_agnpp
   casa_ar                        => clm3%g%l%c%p%pps%casa_ar
   casa_bgnpp                     => clm3%g%l%c%p%pps%casa_bgnpp
   casa_cwdc                      => clm3%g%l%c%p%pps%casa_cwdc
   casa_cwdc_hr                   => clm3%g%l%c%p%pps%casa_cwdc_hr
   casa_cwdc_loss                 => clm3%g%l%c%p%pps%casa_cwdc_loss
   casa_frootc                    => clm3%g%l%c%p%pps%casa_frootc
   casa_frootc_alloc              => clm3%g%l%c%p%pps%casa_frootc_alloc
   casa_frootc_loss               => clm3%g%l%c%p%pps%casa_frootc_loss
   casa_gpp                       => clm3%g%l%c%p%pps%casa_gpp
   casa_hr                        => clm3%g%l%c%p%pps%casa_hr
   casa_leafc                     => clm3%g%l%c%p%pps%casa_leafc
   casa_leafc_alloc               => clm3%g%l%c%p%pps%casa_leafc_alloc
   casa_leafc_loss                => clm3%g%l%c%p%pps%casa_leafc_loss
   casa_litterc                   => clm3%g%l%c%p%pps%casa_litterc
   casa_litterc_hr                => clm3%g%l%c%p%pps%casa_litterc_hr
   casa_litterc_loss              => clm3%g%l%c%p%pps%casa_litterc_loss
   casa_nee                       => clm3%g%l%c%p%pps%casa_nee
   casa_nep                       => clm3%g%l%c%p%pps%casa_nep
   casa_npp                       => clm3%g%l%c%p%pps%casa_npp
   casa_soilc                     => clm3%g%l%c%p%pps%casa_soilc
   casa_soilc_hr                  => clm3%g%l%c%p%pps%casa_soilc_hr
   casa_soilc_loss                => clm3%g%l%c%p%pps%casa_soilc_loss
   casa_woodc                     => clm3%g%l%c%p%pps%casa_woodc
   casa_woodc_alloc               => clm3%g%l%c%p%pps%casa_woodc_alloc
   casa_woodc_loss                => clm3%g%l%c%p%pps%casa_woodc_loss

   ! Get step size
   dtime = get_step_size()

   ! pft loop
   do f = 1, num_soilp
      p = filter_soilp(f)

      ! calculate pft-level summary carbon fluxes and states

      ! autotrophic respiration
      casa_ar(p) = gppfact * fpsn(p) * 12.0_r8 * 1.e-6_r8 ! umolC/m2/s to gC/m2/s
      ! gross primary production
      casa_gpp(p) = fpsn(p) * 12.0_r8 * 1.e-6_r8 ! umolC/m2/s to gC/m2/s

      ! total heterotrophic respiration
      casa_hr(p) = Cflux(p)

      ! net ecosystem exchange
      casa_nee(p) = co2flux(p)

      ! net ecosystem production
      casa_nep(p) =  -1._r8 * co2flux(p)

      ! net primary production
      casa_npp(p) = fnpp(p)

      ! above-ground and below-ground NPP
      casa_agnpp(p) = (livefr(p,LEAF)  + 0.8_r8 * livefr(p,WOOD)) * fnpp(p)
      casa_bgnpp(p) = (livefr(p,FROOT) + 0.2_r8 * livefr(p,WOOD)) * fnpp(p)

      ! leaf C
      casa_leafc(p)       = Tpool_C(p,LEAF)
      casa_leafc_alloc(p) = livefr(p,LEAF) * fnpp(p)
      casa_leafc_loss(p)  = Closs(p,LEAF)

      ! wood C
      casa_woodc(p)       = Tpool_C(p,WOOD)
      casa_woodc_alloc(p) = livefr(p,WOOD) * fnpp(p)
      casa_woodc_loss(p)  = Closs(p,WOOD)

      ! fine root C
      casa_frootc(p)       = Tpool_C(p,FROOT)
      casa_frootc_alloc(p) = livefr(p,FROOT) * fnpp(p)
      casa_frootc_loss(p)  = Closs(p,FROOT)

      ! coarse woody debris C
      casa_cwdc(p)      = Tpool_C(p,CWD)
      casa_cwdc_hr(p)   = Resp_C(p,CWD)
      casa_cwdc_loss(p) = Ctrans(p,CWD_TYPE)

      ! litter C
      casa_litterc(p)      = Tpool_C(p,SURFMET) + Tpool_C(p,SURFSTR) + &
                             Tpool_C(p,SOILMET) + Tpool_C(p,SOILSTR)
      casa_litterc_hr(p)   = Resp_C(p,SURFMET) + Resp_C(p,SURFSTR) + &
                             Resp_C(p,SOILMET) + Resp_C(p,SOILSTR)
      casa_litterc_loss(p) = Ctrans(p,LITTER_TYPE)

      ! soil C
      casa_soilc(p)      = Tpool_C(p,SURFMIC) + Tpool_C(p,SOILMIC) + &
                           Tpool_C(p,SLOW)    + Tpool_C(p,PASSIVE)
      casa_soilc_hr(p)   = Resp_C(p,SURFMIC) + Resp_C(p,SOILMIC) + &
                           Resp_C(p,SLOW)    + Resp_C(p,PASSIVE)
      casa_soilc_loss(p) = Ctrans(p,SOIL_TYPE)

   end do  ! end of pfts loop

end subroutine CASASummary
!===============================================================================

#endif

end module CASAMod
