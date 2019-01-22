module CNFUNMod
!--------------------------------------------------------------------
  !---
! ! DESCRIPTION
! ! The FUN model developed by Fisher et al. 2010 and
! ! end Brzostek et al. 2014. Coded by Mingjie Shi 2015.
! ! Coding logic and structure altered by Rosie Fisher. October 2015. 
! ! Critically, this removes the 'FUN-resistors' idea of Brzostek et
  !  al. 2014
! ! and replaces it with uptake that is proportional to the N/C
  !  exchange rate. 
! ! and adjusts the logic so that FUN does not depends upon the
  !  CLM4.0 'FPG' downregulation idea
! ! and instead it takes C spent on N uptake away from growth.
! ! The critical output so fthis code are sminn_to_plant_fun and
  !  npp_Nuptake, which are the N 
! ! available to the plant for grwoth, and the C spent on obtaining
  !  it. 

! !USES: 
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use shr_log_mod                     , only : errMsg => shr_log_errMsg
  use clm_varctl                      , only : iulog
  use PatchType                       , only : patch
  use ColumnType                      , only : col
  use pftconMod                       , only : pftcon, npcropmin
  use decompMod                       , only : bounds_type
  use clm_varctl                      , only : use_nitrif_denitrif,use_flexiblecn
  use abortutils                      , only : endrun
  use CNVegstateType                  , only : cnveg_state_type
  use CNVegCarbonStateType            , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType             , only : cnveg_carbonflux_type
  use CNVegnitrogenstateType          , only : cnveg_nitrogenstate_type
  use CNVegnitrogenfluxType           , only : cnveg_nitrogenflux_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use SoilBiogeochemNitrogenStateType  , only : soilbiogeochem_nitrogenstate_type

  use SoilBiogeochemCarbonFluxType    , only : soilbiogeochem_carbonflux_type
  use WaterStateBulkType                  , only : waterstatebulk_type
  use WaterFluxBulkType                   , only : waterfluxbulk_type
  use TemperatureType                 , only : temperature_type
  use SoilStateType                   , only : soilstate_type
  use CanopyStateType                 , only : canopystate_type
  use perf_mod                        , only : t_startf, t_stopf

  implicit none
  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public:: readParams            ! Read in parameters needed for FUN
  public:: CNFUNInit             ! FUN calculation initialization
  public:: CNFUN                 ! Run FUN
  
  type, private :: params_type
     real(r8) :: ndays_on        ! number of days to complete leaf onset
     real(r8) :: ndays_off       ! number of days to complete leaf offset
  end type params_type   
 
  !
  type(params_type), private :: params_inst  ! params_inst is
  !  populated in readParamsMod
  !
  !
  ! !PRIVATE DATA MEMBERS:
  real(r8) :: dt              ! decomp timestep (seconds)
  real(r8) :: ndays_on        ! number of days to complete onset
  real(r8) :: ndays_off       ! number of days to complete offset
  
  integer, private, parameter :: COST_METHOD = 2 !new way of doing the N uptake
  ! resistances. see teamwork thread on over-cheap uptake in N
  !  resistors. 
  integer,  private, parameter :: nstp            = 2             ! Number of
  !  calculation part
  integer,  private, parameter :: ncost6          = 6             ! Number of
  !  N transport pathways

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

!
!--------------------------------------------------------------------
  !---
 contains
!--------------------------------------------------------------------
   !---
 subroutine readParams ( ncid )
  !
  ! !USES:
  use ncdio_pio , only : file_desc_t,ncd_io

  ! !ARGUMENTS:
  implicit none
  type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
  !
  ! !LOCAL VARIABLES:
  character(len=32)  :: subname = 'CNFUNParamsType'
  character(len=100) :: errCode = '-Error reading in parameters file:'
  logical            :: readv ! has variable been read in or not
  real(r8)           :: tempr ! temporary to read in parameter
  character(len=100) :: tString ! temp. var for reading
!--------------------------------------------------------------------
  !---

  ! read in parameters

    tString='ndays_on'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%ndays_on=tempr

    tString='ndays_off'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%ndays_off=tempr


 end subroutine readParams

!--------------------------------------------------------------------
 !---
 subroutine CNFUNInit (bounds,cnveg_state_inst,cnveg_carbonstate_inst,cnveg_nitrogenstate_inst)
  !
  ! !DESCRIPTION:
  !
  ! !USES:
  use clm_varcon      , only: secspday, fun_period
  use clm_time_manager, only: get_step_size,get_nstep,get_curr_date,get_days_per_year
  !
  ! !ARGUMENTS:
  type(bounds_type)             , intent(in)    :: bounds
  type(cnveg_state_type)        , intent(inout) :: cnveg_state_inst
  type(cnveg_carbonstate_type)  , intent(inout) :: cnveg_carbonstate_inst
  type(cnveg_nitrogenstate_type), intent(inout) :: cnveg_nitrogenstate_inst
  !
  ! !LOCAL VARIABLES:
  real(r8)          :: dayspyr                  ! days per year (days)
  real(r8)          :: timestep_fun             ! Timestep length for
  !  FUN (s)
  real(r8)          :: numofyear                ! number of days per
  !  year
  integer           :: nstep                    ! time step number
  integer           :: nstep_fun                ! Number of
  !  atmospheric timesteps between calls to FUN
  character(len=32) :: subname = 'CNFUNInit'
!--------------------------------------------------------------------
  !---

! Set local pointers
  associate(ivt                    => patch%itype                                          , & ! Input:  [integer  (:)   ]  p
         leafcn                 => pftcon%leafcn                                        , & ! Input:  leaf C:N (gC/gN)
         leafcn_offset          => cnveg_state_inst%leafcn_offset_patch                 , & ! Output:
         !  [real(r8) (:)   ]  Leaf C:N used by FUN  
         leafc_storage_xfer_acc => cnveg_carbonstate_inst%leafc_storage_xfer_acc_patch  , & ! Output: [real(r8) (:)
         !   ]  Accmulated leaf C transfer (gC/m2)
         storage_cdemand        => cnveg_carbonstate_inst%storage_cdemand_patch         , & ! Output: [real(r8) (:)
         !   ]  C use from the C storage pool
         leafn_storage_xfer_acc => cnveg_nitrogenstate_inst%leafn_storage_xfer_acc_patch, & ! Output: [real(r8) (:)
         !  ]  Accmulated leaf N transfer (gC/m2)                 
         storage_ndemand        => cnveg_nitrogenstate_inst%storage_ndemand_patch         & ! Output: [real(r8) (:)
         !  ]  N demand during the offset period 
         )
  !--------------------------------------------------------------------
  !---
  ! Calculate some timestep-related values.
  !--------------------------------------------------------------------
  !---
  ! set time steps
  dt           = real(get_step_size(), r8)
  dayspyr      = get_days_per_year()
  nstep        = get_nstep()
  timestep_fun = real(secspday * fun_period)
  nstep_fun    = int(secspday * dayspyr / dt) 

  ndays_on     = params_inst%ndays_on
  ndays_off    = params_inst%ndays_off

  !--------------------------------------------------------------------
  !---
  ! Decide if FUN will be called on this timestep.
  !--------------------------------------------------------------------
  !---
  numofyear = nstep/nstep_fun
  if (mod(nstep,nstep_fun) == 0) then
     leafcn_offset(bounds%begp:bounds%endp)          = leafcn(ivt(bounds%begp:bounds%endp))
     storage_cdemand(bounds%begp:bounds%endp)        = 0._r8
     storage_ndemand(bounds%begp:bounds%endp)        = 0._r8
     leafn_storage_xfer_acc(bounds%begp:bounds%endp) = 0._r8
     leafc_storage_xfer_acc(bounds%begp:bounds%endp) = 0._r8
  end if  
!--------------------------------------------------------------------
  !---
  end associate
  end subroutine CNFUNInit 
!--------------------------------------------------------------------
  !---

  !--------------------------------------------------------------------
  !---  
  ! Start the CNFUN subroutine
  !--------------------------------------------------------------------
  !---
  subroutine CNFUN(bounds,num_soilc, filter_soilc,num_soilp&
       &,filter_soilp,waterstatebulk_inst, &
       & waterfluxbulk_inst,temperature_inst,soilstate_inst&
       &,cnveg_state_inst,cnveg_carbonstate_inst,&
       & cnveg_carbonflux_inst,cnveg_nitrogenstate_inst&
       &,cnveg_nitrogenflux_inst                ,&
       & soilbiogeochem_nitrogenflux_inst&
       &,soilbiogeochem_carbonflux_inst,canopystate_inst,&
       & soilbiogeochem_nitrogenstate_inst)

! !USES:
   use clm_time_manager, only : get_step_size, get_curr_date, get_days_per_year 
   use clm_varpar      , only : nlevdecomp
   use clm_varcon      , only : secspday, smallValue, fun_period, tfrz, dzsoi_decomp, spval
   use clm_varctl      , only : use_nitrif_denitrif
   use PatchType       , only : patch
   use subgridAveMod   , only : p2c
   use pftconMod       , only : npcropmin
!
! !ARGUMENTS: 
   type(bounds_type)                       , intent(in)    :: bounds
   integer                                 , intent(in)    :: num_soilc             ! number of soil columns in filter
   integer                                 , intent(in)    :: filter_soilc(:)       ! filter for soil columns
   integer                                 , intent(in)    :: num_soilp             ! number of soil patches in filter
   integer                                 , intent(in)    :: filter_soilp(:)       ! filter for soil patches
   type(waterstatebulk_type)                   , intent(in)    :: waterstatebulk_inst
   type(waterfluxbulk_type)                    , intent(in)    :: waterfluxbulk_inst
   type(temperature_type)                  , intent(in)    :: temperature_inst
   type(soilstate_type)                    , intent(in)    :: soilstate_inst
   type(cnveg_state_type)                  , intent(inout) :: cnveg_state_inst
   type(cnveg_carbonstate_type)            , intent(inout) :: cnveg_carbonstate_inst
   type(cnveg_carbonflux_type)             , intent(inout) :: cnveg_carbonflux_inst
   type(cnveg_nitrogenstate_type)          , intent(inout) :: cnveg_nitrogenstate_inst
   type(cnveg_nitrogenflux_type)           , intent(inout) :: cnveg_nitrogenflux_inst
   type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
   type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst 
   type(canopystate_type)                  , intent(inout) :: canopystate_inst  
   type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
  !
  ! !LOCAL VARIABLES:
  ! local pointers to implicit in arrays
  ! 
  !--------------------------------------------------------------------
  ! ------------
  ! Integer parameters
  !--------------------------------------------------------------------
  !-----------
  integer,  parameter :: icostFix        = 1             ! Process
  !  number for fixing.
  integer,  parameter :: icostRetrans    = 2             ! Process
  !  number for retranslocation.
  integer,  parameter :: icostActiveNO3  = 3             ! Process
  !  number for mycorrhizal uptake of NO3.
  integer,  parameter :: icostActiveNH4  = 4             ! Process
  !  number for mycorrhizal uptake of NH4
  integer,  parameter :: icostnonmyc_no3  = 5            ! Process
  !  number for nonmyc uptake of NO3.
  integer,  parameter :: icostnonmyc_nh4  = 6            ! Process
  !  number for nonmyc uptake of NH4.
  real(r8), parameter :: big_cost        = 1000000000._r8! An arbitrarily large cost

  !  array index when plant is fixing
  integer, parameter :: plants_are_fixing = 1
  integer, parameter :: plants_not_fixing = 2

  !  array index for ECM step versus AM step
  integer, parameter :: ecm_step          = 1
  integer, parameter :: am_step           = 2
  !  arbitrary large cost (gC/gN).
  !--------------------------------------------------------------------
  !-----------------------------------------------
  ! Local Real variables.
  !--------------------------------------------------------------------
  !-----------------------------------------------
  real(r8)  :: excess                                                ! excess N taken up by transpiration    (gN/m2) 
  real(r8)  :: steppday                                              ! model time steps in each day          (-)
  real(r8)  :: rootc_dens_step                                       ! root C for each PFT in each soil layer(gC/m2)
  real(r8)  :: retrans_limit1                                        ! a temporary variable for leafn        (gN/m2)
  real(r8)  :: qflx_tran_veg_layer                                   ! transpiration in each soil layer      (mm H2O/S)
  real(r8)  :: dn                                                    ! Increment of N                        (gN/m2)  
  real(r8)  :: dn_retrans                                            ! Increment of N                        (gN/m2)
  real(r8)  :: dnpp                                                  ! Increment of NPP                      (gC/m2)
  real(r8)  :: dnpp_retrans                                          ! Increment of NPP                      (gC/m2)
  real(r8)  :: rootc_dens(bounds%begp:bounds%endp,1:nlevdecomp)      ! the root carbon density               (gC/m2)
  real(r8)  :: rootC(bounds%begp:bounds%endp)                        ! root biomass                          (gC/m2)
  real(r8)  :: permyc(bounds%begp:bounds%endp,1:nstp)                ! the arrary for the ECM and AM ratio   (-) 
  real(r8)  :: kc_active(bounds%begp:bounds%endp,1:nstp)             ! the kc_active parameter               (gC/m2)
  real(r8)  :: kn_active(bounds%begp:bounds%endp,1:nstp)             ! the kn_active parameter               (gC/m2)
  real(r8)  :: availc_pool(bounds%begp:bounds%endp)                  ! The avaible C pool for allocation     (gC/m2)
  real(r8)  :: plantN(bounds%begp:bounds%endp)                       ! Plant N                               (gN/m2)
  real(r8)  :: plant_ndemand_pool(bounds%begp:bounds%endp)           ! The N demand pool                     (gN/m2)
  real(r8)  :: plant_ndemand_pool_step(bounds%begp:bounds%endp,1:nstp)   ! the N demand pool                     (gN/m2)
  real(r8)  :: leafn_step(bounds%begp:bounds%endp,1:nstp)            ! N loss based for deciduous trees      (gN/m2)
  real(r8)  :: leafn_retrans_step(bounds%begp:bounds%endp,1:nstp)    ! N loss based for deciduous trees      (gN/m2) 
  real(r8)  :: litterfall_n(bounds%begp:bounds%endp)                 ! N loss based on the leafc to litter   (gN/m2)  
  real(r8)  :: litterfall_n_step(bounds%begp:bounds%endp,1:nstp)       ! N loss based on the leafc to litter   (gN/m2)
  real(r8)  :: litterfall_c_step(bounds%begp:bounds%endp,1:nstp)       ! N loss based on the leafc to litter   (gN/m2)
  real(r8)  :: tc_soisno(bounds%begc:bounds%endc,1:nlevdecomp)       ! Soil temperature            (degrees Celsius)
  real(r8)  :: npp_remaining(bounds%begp:bounds%endp,1:nstp)         ! A temporary variable for npp_remaining(gC/m2) 
  real(r8)  :: n_passive_step(bounds%begp:bounds%endp,1:nstp)        ! N taken up by transpiration at substep(gN/m2)
  real(r8)  :: n_passive_acc(bounds%begp:bounds%endp)                ! N acquired by passive uptake          (gN/m2)
  real(r8)  :: cost_retran(bounds%begp:bounds%endp,1:nlevdecomp)     ! cost of retran                        (gC/gN)
  real(r8)  :: cost_fix(bounds%begp:bounds%endp,1:nlevdecomp)        ! cost of fixation                      (gC/gN)
  real(r8)  :: cost_resis(bounds%begp:bounds%endp,1:nlevdecomp)      ! cost of resis                         (gC/gN)
  real(r8)  :: cost_res_resis(bounds%begp:bounds%endp,1:nlevdecomp)  ! The cost of resis                     (gN/gC)
  real(r8)  :: n_fix_acc(bounds%begp:bounds%endp,1:nstp)             ! N acquired by fixation                (gN/m2)
  real(r8)  :: n_fix_acc_total(bounds%begp:bounds%endp)              ! N acquired by fixation                (gN/m2)
  real(r8)  :: npp_fix_acc(bounds%begp:bounds%endp,1:nstp)           ! Amount of NPP used by fixation        (gC/m2)
  real(r8)  :: npp_fix_acc_total(bounds%begp:bounds%endp)            ! Amount of NPP used by fixation        (gC/m2)
  real(r8)  :: n_retrans_acc(bounds%begp:bounds%endp,1:nstp)         ! N acquired by retranslocation         (gN/m2)
  real(r8)  :: n_retrans_acc_total(bounds%begp:bounds%endp)          ! N acquired by retranslocation         (gN/m2)
  real(r8)  :: free_nretrans_acc(bounds%begp:bounds%endp,1:nstp)     ! N acquired by retranslocation         (gN/m2)
  real(r8)  :: npp_retrans_acc(bounds%begp:bounds%endp,1:nstp)       ! NPP used for the extraction           (gC/m2)
  real(r8)  :: npp_retrans_acc_total(bounds%begp:bounds%endp)        ! NPP used for the extraction           (gC/m2)
  real(r8)  :: nt_uptake(bounds%begp:bounds%endp,1:nstp)             ! N uptake from retrans, active, and fix(gN/m2)
  real(r8)  :: npp_uptake(bounds%begp:bounds%endp,1:nstp)            ! NPP used by the uptakes               (gC/m2)

  !----------NITRIF_DENITRIF-------------!

  real(r8)  :: sminn_no3_diff                                        ! A temporary limit for N uptake                  (gN/m2)
  real(r8)  :: sminn_nh4_diff                                        ! A temporary limit for N uptake                  (gN/m2)
  real(r8)  :: active_no3_limit1                                     ! A temporary limit for N uptake                  (gN/m2)
  real(r8)  :: active_nh4_limit1                                     ! A temporary limit for N uptake                  (gN/m2)
  real(r8)  :: cost_active_no3(bounds%begp:bounds%endp,1:nlevdecomp) ! cost of mycorrhizal                             (gC/gN)
  real(r8)  :: cost_active_nh4(bounds%begp:bounds%endp,1:nlevdecomp) ! cost of mycorrhizal                             (gC/gN)
  real(r8)  :: cost_nonmyc_no3(bounds%begp:bounds%endp,1:nlevdecomp) ! cost of nonmyc                                  (gC/gN)
  real(r8)  :: cost_nonmyc_nh4(bounds%begp:bounds%endp,1:nlevdecomp) ! cost of nonmyc                                  (gC/gN)

  real(r8)  :: sminn_no3_conc(bounds%begc:bounds%endc,1:nlevdecomp)             ! Concentration of no3 in soil water (gN/gH2O)
  real(r8)  :: sminn_no3_conc_step(bounds%begp:bounds%endp,1:nlevdecomp,1:nstp) ! A temporary variable for soil mineral N (gN/gH2O)
  real(r8)  :: sminn_no3_layer(bounds%begc:bounds%endc,1:nlevdecomp)            ! Available no3 in each soil layer (gN/m2)
  real(r8)  :: sminn_no3_layer_step(bounds%begp:bounds%endp,1:nlevdecomp,1:nstp)! A temporary variable for soil no3 (gN/m2) 
  real(r8)  :: sminn_no3_uptake(bounds%begp:bounds%endp,1:nlevdecomp,1:nstp)    ! A temporary variable for soil mineral N (gN/m2/s)
  real(r8)  :: sminn_nh4_conc(bounds%begc:bounds%endc,1:nlevdecomp)             ! Concentration of nh4 in soil water (gN/gH2O)
  real(r8)  :: sminn_nh4_conc_step(bounds%begp:bounds%endp,1:nlevdecomp,1:nstp) ! A temporary variable for soil mineral N (gN/gH2O)
  real(r8)  :: sminn_nh4_layer(bounds%begc:bounds%endc,1:nlevdecomp)            ! Available nh4 in each soil layer (gN/m2)
  real(r8)  :: sminn_nh4_layer_step(bounds%begp:bounds%endp,1:nlevdecomp,1:nstp)! A temporary variable for soil mineral N (gN/m2)
  real(r8)  :: sminn_nh4_uptake(bounds%begp:bounds%endp,1:nlevdecomp,1:nstp)    ! A temporary variable for soil mineral N (gN/m2/s)

  real(r8)  :: active_no3_uptake1(bounds%begp:bounds%endp,1:nlevdecomp)         ! no3 mycorrhizal uptake (gN/m2)
  real(r8)  :: active_nh4_uptake1(bounds%begp:bounds%endp,1:nlevdecomp)         ! nh4 mycorrhizal uptake (gN/m2) 
  real(r8)  :: nonmyc_no3_uptake1(bounds%begp:bounds%endp,1:nlevdecomp)         ! no3 non-mycorrhizal uptake (gN/m2) 
  real(r8)  :: nonmyc_nh4_uptake1(bounds%begp:bounds%endp,1:nlevdecomp)         ! nh4 non-mycorrhizal uptake (gN/m2) 
  real(r8)  :: active_no3_uptake2(bounds%begp:bounds%endp,1:nlevdecomp)         ! no3 mycorrhizal uptake (gN/m2) 
  real(r8)  :: active_nh4_uptake2(bounds%begp:bounds%endp,1:nlevdecomp)         ! nh4 mycorrhizal uptake (gN/m2) 
  real(r8)  :: nonmyc_no3_uptake2(bounds%begp:bounds%endp,1:nlevdecomp)         ! no3 non-mycorrhizal uptake (gN/m2) 
  real(r8)  :: nonmyc_nh4_uptake2(bounds%begp:bounds%endp,1:nlevdecomp)         ! nh4 non-mycorrhizal uptake (gN/m2) 
  real(r8)  :: n_am_no3_acc(bounds%begp:bounds%endp)                            ! AM no3 uptake (gN/m2)
  real(r8)  :: n_am_nh4_acc(bounds%begp:bounds%endp)                            ! AM nh4 uptake (gN/m2)
  real(r8)  :: n_ecm_no3_acc(bounds%begp:bounds%endp)                           ! ECM no3 uptake (gN/m2)
  real(r8)  :: n_ecm_nh4_acc(bounds%begp:bounds%endp)                           ! ECM nh4 uptake (gN/m2)
  real(r8)  :: n_active_no3_acc(bounds%begp:bounds%endp,1:nstp)                 ! Mycorrhizal no3 uptake (gN/m2)
  real(r8)  :: n_active_nh4_acc(bounds%begp:bounds%endp,1:nstp)                 ! Mycorrhizal nh4 uptake (gN/m2)
  real(r8)  :: n_nonmyc_no3_acc(bounds%begp:bounds%endp,1:nstp)                 ! Non-myc     no3 uptake (gN/m2)
  real(r8)  :: n_nonmyc_nh4_acc(bounds%begp:bounds%endp,1:nstp)                 ! Non-myc     nh4 uptake (gN/m2) 
  real(r8)  :: n_active_no3_acc_total(bounds%begp:bounds%endp)                  ! Mycorrhizal no3 uptake (gN/m2)
  real(r8)  :: n_active_nh4_acc_total(bounds%begp:bounds%endp)                  ! Mycorrhizal no3 uptake (gN/m2)   
     
  real(r8)  :: n_nonmyc_no3_acc_total(bounds%begp:bounds%endp)                  ! Non-myc     no3 uptake (gN/m2)
  real(r8)  :: n_nonmyc_nh4_acc_total(bounds%begp:bounds%endp)                  ! Non-myc     no3 uptake (gN/m2)
  real(r8)  :: npp_active_no3_acc(bounds%begp:bounds%endp,1:nstp)               ! Mycorrhizal no3 uptake used C (gC/m2)
  real(r8)  :: npp_active_nh4_acc(bounds%begp:bounds%endp,1:nstp)               ! Mycorrhizal nh4 uptake used C (gC/m2)
  real(r8)  :: npp_nonmyc_no3_acc(bounds%begp:bounds%endp,1:nstp)               ! Non-myc     no3 uptake used C (gC/m2)
  real(r8)  :: npp_nonmyc_nh4_acc(bounds%begp:bounds%endp,1:nstp)               ! Non-myc     nh4 uptake used C (gC/m2)
  real(r8)  :: npp_active_no3_acc_total(bounds%begp:bounds%endp)                ! Mycorrhizal no3 uptake used C (gC/m2)
  real(r8)  :: npp_active_nh4_acc_total(bounds%begp:bounds%endp)                ! Mycorrhizal nh4 uptake used C (gC/m2)
  real(r8)  :: npp_nonmyc_no3_acc_total(bounds%begp:bounds%endp)                ! Non-myc     no3 uptake used C (gC/m2)
  real(r8)  :: npp_nonmyc_nh4_acc_total(bounds%begp:bounds%endp)                ! Non-myc     nh4 uptake used C (gC/m2)
  real(r8)  :: n_am_no3_retrans(bounds%begp:bounds%endp)                        ! AM no3 uptake for offset (gN/m2)
  real(r8)  :: n_am_nh4_retrans(bounds%begp:bounds%endp)                        ! AM nh4 uptake for offset (gN/m2)
  real(r8)  :: n_ecm_no3_retrans(bounds%begp:bounds%endp)                       ! ECM no3 uptake for offset (gN/m2)
  real(r8)  :: n_ecm_nh4_retrans(bounds%begp:bounds%endp)                       ! ECM nh4 uptake for offset (gN/m2)
  real(r8)  :: n_active_no3_retrans(bounds%begp:bounds%endp,1:nstp)             ! Mycorrhizal no3 for offset (gN/m2)
  real(r8)  :: n_active_nh4_retrans(bounds%begp:bounds%endp,1:nstp)             ! Mycorrhizal nh4 for offset (gN/m2)
  real(r8)  :: n_nonmyc_no3_retrans(bounds%begp:bounds%endp,1:nstp)             ! Non-myc     no3 for offset (gN/m2)
  real(r8)  :: n_nonmyc_nh4_retrans(bounds%begp:bounds%endp,1:nstp)             ! Non-myc     nh4 for offset (gN/m2)
  real(r8)  :: n_active_no3_retrans_total(bounds%begp:bounds%endp)              ! Mycorrhizal no3 for offset (gN/m2)
  real(r8)  :: n_active_nh4_retrans_total(bounds%begp:bounds%endp)              ! Mycorrhizal nh4 for offset (gN/m2)
  real(r8)  :: n_nonmyc_no3_retrans_total(bounds%begp:bounds%endp)              ! Non-myc     no3 for offset (gN/m2)
  real(r8)  :: n_nonmyc_nh4_retrans_total(bounds%begp:bounds%endp)              ! Non-myc     nh4 for offset (gN/m2)
  real(r8)  :: n_passive_no3_vr(bounds%begp:bounds%endp,1:nlevdecomp)           ! Layer passive no3 uptake (gN/m2)
  real(r8)  :: n_passive_nh4_vr(bounds%begp:bounds%endp,1:nlevdecomp)           ! Layer passive nh4 uptake (gN/m2)
  real(r8)  :: n_fix_no3_vr(bounds%begp:bounds%endp,1:nlevdecomp)               ! Layer fixation no3 uptake (gN/m2)
  real(r8)  :: n_fix_nh4_vr(bounds%begp:bounds%endp,1:nlevdecomp)               ! Layer fixation nh4 uptake (gN/m2)
  real(r8)  :: n_active_no3_vr(bounds%begp:bounds%endp,1:nlevdecomp)            ! Layer mycorrhizal no3 uptake (gN/m2)
  real(r8)  :: n_nonmyc_no3_vr(bounds%begp:bounds%endp,1:nlevdecomp)            ! Layer non-myc     no3 uptake (gN/m2)
  real(r8)  :: n_active_nh4_vr(bounds%begp:bounds%endp,1:nlevdecomp)            ! Layer mycorrhizal nh4 uptake (gN/m2)
  real(r8)  :: n_nonmyc_nh4_vr(bounds%begp:bounds%endp,1:nlevdecomp)            ! Layer non-myc     nh4 uptake (gN/m2)
  real(r8)  :: npp_active_no3_retrans(bounds%begp:bounds%endp,1:nstp)           ! Mycorrhizal no3 uptake used C for offset (gN/m2)
  real(r8)  :: npp_active_nh4_retrans(bounds%begp:bounds%endp,1:nstp)           ! Mycorrhizal nh4 uptake used C for offset (gN/m2) 
  real(r8)  :: npp_nonmyc_no3_retrans(bounds%begp:bounds%endp,1:nstp)           ! Non-myc no3 uptake used C for offset (gN/m2)
  real(r8)  :: npp_nonmyc_nh4_retrans(bounds%begp:bounds%endp,1:nstp)           ! Non-myc nh4 uptake used C for offset (gN/m2) 
  real(r8)  :: npp_active_no3_retrans_total(bounds%begp:bounds%endp)            ! Mycorrhizal no3 uptake used C for offset (gN/m2)
  real(r8)  :: npp_active_nh4_retrans_total(bounds%begp:bounds%endp)            ! Mycorrhizal nh4 uptake used C for offset (gN/m2)
  real(r8)  :: npp_nonmyc_no3_retrans_total(bounds%begp:bounds%endp)             ! Non-myc no3 uptake used C for offset (gN/m2)
  real(r8)  :: npp_nonmyc_nh4_retrans_total(bounds%begp:bounds%endp)            ! Non-myc nh4 uptake used C for offset (gN/m2)
  

  real(r8)  :: costNit(1:nlevdecomp,ncost6)                          ! Cost of N via each process                      (gC/gN)

  ! Uptake fluxes for COST_METHOD=2
  ! actual npp to each layer for each uptake process
  real(r8)  ::                   npp_to_fixation(1:nlevdecomp) 
  real(r8)  ::                   npp_to_retrans(1:nlevdecomp)
  real(r8)  ::                   npp_to_active_nh4(1:nlevdecomp)
  real(r8)  ::                   npp_to_nonmyc_nh4(1:nlevdecomp)
  real(r8)  ::                   npp_to_active_no3(1:nlevdecomp)
  real(r8)  ::                   npp_to_nonmyc_no3 (1:nlevdecomp)  

  ! fraction of carbon to each uptake process 
  real(r8)  ::                   npp_frac_to_fixation(1:nlevdecomp) 
  real(r8)  ::                   npp_frac_to_retrans(1:nlevdecomp)
  real(r8)  ::                   npp_frac_to_active_nh4(1:nlevdecomp)
  real(r8)  ::                   npp_frac_to_nonmyc_nh4(1:nlevdecomp)
  real(r8)  ::                   npp_frac_to_active_no3(1:nlevdecomp)
  real(r8)  ::                   npp_frac_to_nonmyc_no3 (1:nlevdecomp)  
   
  ! hypothetical fluxes on N in each layer 
  real(r8)  ::                  n_exch_fixation(1:nlevdecomp)        ! N aquired from one unit of C for fixation (unitless)
  real(r8)  ::                  n_exch_retrans(1:nlevdecomp)         ! N aquired from one unit of C for retrans (unitless)
  real(r8)  ::                  n_exch_active_nh4(1:nlevdecomp)      ! N aquired from one unit of C for act nh4(unitless)
  real(r8)  ::                  n_exch_nonmyc_nh4(1:nlevdecomp)      ! N aquired from one unit of C for nonmy nh4 (unitless) 
  real(r8)  ::                  n_exch_active_no3(1:nlevdecomp)      ! N aquired from one unit of C for act no3 (unitless)
  real(r8)  ::                  n_exch_nonmyc_no3(1:nlevdecomp)      ! N aquired from one unit of C for nonmyc no3 (unitless) 
  
   !actual fluxes of N in each layer
  real(r8)  ::                  n_from_fixation(1:nlevdecomp)        ! N aquired in each layer for fixation       (gN m-2 s-1)
  real(r8)  ::                  n_from_retrans(1:nlevdecomp)         ! N aquired in each layer of C for retrans (gN m-2 s-1)
  real(r8)  ::                  n_from_active_nh4(1:nlevdecomp)      ! N aquired in each layer of C for act nh4 (gN m-2 s-1)
  real(r8)  ::                  n_from_nonmyc_nh4(1:nlevdecomp)      ! N aquired in each layer of C for nonmy nh4 (gN m-2 s-1)
  real(r8)  ::                  n_from_active_no3(1:nlevdecomp)      ! N aquired in each layer of C for act no3 (gN m-2 s-1)
  real(r8)  ::                  n_from_nonmyc_no3(1:nlevdecomp)      ! N aquired in each layer of C for nonmyc no3 (gN m-2 s-1) 

  real(r8)  :: free_Nretrans(bounds%begp:bounds%endp)                     ! the total amount of NO3 and NH4                 (gN/m3/s)

  ! Uptake fluxes for COST_METHOD=2
  !actual fluxes of N in each layer
  real(r8)  ::                  frac_ideal_C_use                     ! How much less C do we use for 'buying' N than that
  !  needed to get to the ideal ratio?  fraction. 

  real(r8)  ::                  N_acquired
  real(r8)  ::                  C_spent
  real(r8)  ::                  leaf_narea ! leaf n per unit leaf
  !  area in gN/m2 (averaged across canopy, which is OK for the cost
  !   calculation)

                       
  real(r8)  ::                  sum_n_acquired                          ! Sum N aquired from one unit of C (unitless)  
  real(r8)  ::                  burned_off_carbon                       ! carbon wasted by poor allocation algorithm. If
  !  this is too big, we need a better iteration. 
  real(r8)   ::                 temp_n_flux  
  real(r8)  ::                  delta_cn                                ! difference between 'ideal' leaf CN ration and
  !  actual leaf C:N ratio. C/N
  real(r8) :: excess_carbon        ! how much carbon goes into the leaf C
  !  pool on account of the flexibleCN modifications.   
  real(r8) :: excess_carbon_acc    ! excess accumulated over layers.
  !  WITHOUT GROWTH RESP
  real(r8) :: fixerfrac            ! what fraction of plants can fix?
  real(r8) :: npp_to_spend         ! how much carbon do we need to get
  !  rid of? 
  real(r8) :: soil_n_extraction    ! calculates total N pullled from
  !  soil
  real(r8) :: total_N_conductance  !inverse of C to of N for whole soil
  ! -leaf pathway
  real(r8) :: total_N_resistance   ! C to of N for whole soil -leaf
  !  pathway
  real(r8) :: free_RT_frac=0.0_r8  !fraction of N retranslocation which is automatic/free.
  !  SHould be made into a PFT parameter. 
  
  real(r8) :: paid_for_n_retrans
  real(r8) :: free_n_retrans
  real(r8) :: total_c_spent_retrans
  real(r8) :: total_c_accounted_retrans

  
  !------end of not_use_nitrif_denitrif------!
  !--------------------------------------------------------------------
  !------------
  ! Local Integer variables
  !--------------------------------------------------------------------
  !------------
  integer   :: fn                                ! number of values
  !  in pft filter
  integer   :: fp                                ! lake filter pft
  !  index
  integer   :: fc                                ! lake filter column
  !  index
  integer   :: p, c                              ! pft index
  integer   :: g, l                              ! indices
  integer   :: j, i, k                           ! soil/snow level
  !  index
  integer   :: istp                              ! Loop counters/work
  integer   :: icost                             ! a local index
  integer   :: fixer                             ! 0 = non-fixer, 1
  ! =fixer 
  logical   :: unmetDemand                       ! True while there
  !  is still demand for N
  logical   :: local_use_flexibleCN              ! local version of use_flexCN
  integer   :: FIX                               ! for loop. 1 for
  !  fixers, 2 for non fixers. This will become redundant with the
  !   'fixer' parameter if it works. 
  
  !--------------------------------------------------------------------
  !---------------------------------
  associate(ivt                    => patch%itype                                          , & ! Input:   [integer  (:) ]  p
         leafcn                 => pftcon%leafcn                                        , & ! Input:   leaf C:N (gC/gN)
         lflitcn                => pftcon%lflitcn                                       , & ! Input:   leaf litter C:N (gC/gN)
         season_decid           => pftcon%season_decid                                  , & ! Input:   binary flag for seasonal
         ! -deciduous leaf habit (0 or 1)
         stress_decid           => pftcon%stress_decid                                  , & ! Input:   binary flag for stress
         ! -deciduous leaf habit (0 or 1)
         a_fix                  => pftcon%a_fix                                         , & ! Input:   A BNF parameter
         b_fix                  => pftcon%b_fix                                         , & ! Input:   A BNF parameter
         c_fix                  => pftcon%c_fix                                         , & ! Input:   A BNF parameter
         s_fix                  => pftcon%s_fix                                         , & ! Input:   A BNF parameter
         akc_active             => pftcon%akc_active                                    , & ! Input:   A mycorrhizal uptake
         !  parameter
         akn_active             => pftcon%akn_active                                    , & ! Input:   A mycorrhizal uptake
         !  parameter
         ekc_active             => pftcon%ekc_active                                    , & ! Input:   A mycorrhizal uptake
         !  parameter
         ekn_active             => pftcon%ekn_active                                    , & ! Input:   A mycorrhizal upatke
         !  parameter
         kc_nonmyc              => pftcon%kc_nonmyc                                     , & ! Input:   A non-mycorrhizal uptake
         !  parameter
         kn_nonmyc              => pftcon%kn_nonmyc                                     , & ! Input:   A non-mycorrhizal uptake
         !  parameter
         perecm                 => pftcon%perecm                                        , & ! Input:   The fraction of ECM
         ! -associated PFT 
         grperc                 => pftcon%grperc                                        , & ! Input:   growth percentage
         fun_cn_flex_a           => pftcon%fun_cn_flex_a                                , & ! Parameter a of FUN-flexcn link code (def 5)
         fun_cn_flex_b           => pftcon%fun_cn_flex_b                                , & ! Parameter b of FUN-flexcn link code (def 200)
         fun_cn_flex_c           => pftcon%fun_cn_flex_c                                , & ! Parameter b of FUN-flexcn link code (def 80)         
         FUN_fracfixers          => pftcon%FUN_fracfixers                               , & ! Fraction of C that can be used for fixation.    
         leafcn_offset          => cnveg_state_inst%leafcn_offset_patch                 , & ! Output:
         !  [real(r8)  (:)]  Leaf C:N used by FUN
         plantCN                => cnveg_state_inst%plantCN_patch                       , & ! Output:  [real(r8)  (:)]  Plant
         !  C:N used by FUN
         onset_flag             => cnveg_state_inst%onset_flag_patch                    , & ! Output:  [real(r8)  (:)]  onset
         !  flag
         offset_flag            => cnveg_state_inst%offset_flag_patch                   , & ! Output:  [real(r8)  (:)]  offset
         !  flag
         availc                 => cnveg_carbonflux_inst%availc_patch                   , & ! Iutput:  [real(r8)  (:)]  C flux
         !  available for allocation (gC/m2/s)
         leafc                  => cnveg_carbonstate_inst%leafc_patch                   , & ! Input:   [real(r8)  (:)]  (gC/m2)
         !  leaf C
         leafc_storage          => cnveg_carbonstate_inst%leafc_storage_patch           , & ! Input:   [real(r8)
         !  (:)]  (gC/m2) leaf C storage
         frootc                 => cnveg_carbonstate_inst%frootc_patch                  , & ! Input:   [real(r8)
         !  (:)]  (gC/m2) fine root C
         frootc_storage         => cnveg_carbonstate_inst%frootc_storage_patch          , & ! Input:   [real(r8)
         !  (:)]  (gC/m2) fine root C storage
         livestemc              => cnveg_carbonstate_inst%livestemc_patch               , & ! Input:   [real(r8)
         !  (:)]  (gC/m2) live stem C
         livecrootc             => cnveg_carbonstate_inst%livecrootc_patch              , & ! Input:   [real(r8)
         !  (:)]  (gC/m2) live coarse root C
         leafc_storage_xfer_acc => cnveg_carbonstate_inst%leafc_storage_xfer_acc_patch  , & ! uutput:  [real(r8)
         !  (:)]  Accmulated leaf C transfer (gC/m2)
         storage_cdemand        => cnveg_carbonstate_inst%storage_cdemand_patch         , & ! Output:  [real(r8)
         !  (:)]  C use f rom the C storage pool
         tlai                   => canopystate_inst%tlai_patch                          , & ! Input:  [real(r8) (:)   ] one
         ! -sided leaf area index
         leafn                  => cnveg_nitrogenstate_inst%leafn_patch                 , & ! Input:   [real(r8)  (:)]
         !   (gN/m2) leaf N
         frootn                 => cnveg_nitrogenstate_inst%frootn_patch                , & ! Input:   [real(r8)  (:)]
         !   (gN/m2) fine root N
         livestemn              => cnveg_nitrogenstate_inst%livestemn_patch             , & ! Input:   [real(r8)  (:)]
         !   (gN/m2) live stem N
         livecrootn             => cnveg_nitrogenstate_inst%livecrootn_patch            , & ! Input:   [real(r8)  (:)]
         !   (gN/m2) live coarse root N
         leafn_storage_xfer_acc => cnveg_nitrogenstate_inst%leafn_storage_xfer_acc_patch, & ! Output:  [real(r8)  (:)]
         !   Accmulated leaf N transfer (gC/m2)
         storage_ndemand        => cnveg_nitrogenstate_inst%storage_ndemand_patch       , & ! Output:  [real(r8)  (:)]
         !   N demand during the offset period
         leafc_to_litter        => cnveg_carbonflux_inst%leafc_to_litter_patch          , & ! Output:  [real(r8)
         !  (:) ]  leaf C litterfall (gC/m2/s)
         leafc_to_litter_fun    => cnveg_carbonflux_inst%leafc_to_litter_fun_patch      , & ! Output:  [real(r8)
         !  (:) ]  leaf C litterfall used by FUN (gC/m2/s)
         prev_leafc_to_litter   => cnveg_carbonflux_inst%prev_leafc_to_litter_patch     , & ! Output: [real(r8) (:)
         !  ] previous timestep leaf C litterfall flux (gC/m2/s)
         leafc_storage_to_xfer  => cnveg_carbonflux_inst%leafc_storage_to_xfer_patch    , & ! Output:  [real(r8)
         !  (:) ] 
         npp_Nactive            => cnveg_carbonflux_inst%npp_Nactive_patch              , & ! Output:  [real(r8)
         !  (:) ]  Mycorrhizal N uptake used C (gC/m2/s)
         npp_Nnonmyc            => cnveg_carbonflux_inst%npp_Nnonmyc_patch              , & ! Output:  [real(r8)
         !  (:) ]  Non-mycorrhizal N uptake use C (gC/m2/s)
         npp_Nam                => cnveg_carbonflux_inst%npp_Nam_patch                  , & ! Output:  [real(r8)
         !  (:) ]  AM uptake use C (gC/m2/s)
         npp_Necm               => cnveg_carbonflux_inst%npp_Necm_patch                 , & ! Output:  [real(r8)
         !  (:) ]  ECM uptake use C (gC/m2/s)
         npp_Nactive_no3        => cnveg_carbonflux_inst%npp_Nactive_no3_patch          , & ! Output:  [real(r8)
         !  (:) ]  Mycorrhizal N uptake used C (gC/m2/s)
         npp_Nnonmyc_no3        => cnveg_carbonflux_inst%npp_Nnonmyc_no3_patch          , & ! Output:  [real(r8)
         !  (:) ]  Non-myco uptake use C (gC/m2/s) rrhizal N uptake
         !   (gN/m2/s)         
         npp_Nam_no3            => cnveg_carbonflux_inst%npp_Nam_no3_patch              , & ! Output:  [real(r8)
         !  (:) ]  AM uptake use C (gC/m2/s)
         npp_Necm_no3           => cnveg_carbonflux_inst%npp_Necm_no3_patch             , & ! Output:  [real(r8)
         !  (:) ]  ECM uptake use C (gC/m2/s)
         npp_Nactive_nh4        => cnveg_carbonflux_inst%npp_Nactive_nh4_patch          , & ! Output:  [real(r8)
         !  (:) ]  Mycorrhizal N uptake used C (gC/m2/s)
         npp_Nnonmyc_nh4        => cnveg_carbonflux_inst%npp_Nnonmyc_nh4_patch          , & ! Output:  [real(r8)
         !  (:) ]  Non-mycorrhizal N uptake used C (gC/m2/s)
         npp_Nam_nh4            => cnveg_carbonflux_inst%npp_Nam_nh4_patch              , & ! Output:  [real(r8)
         !  (:) ]  AM uptake used C(gC/m2/s)
         npp_Necm_nh4           => cnveg_carbonflux_inst%npp_Necm_nh4_patch             , & ! Output:  [real(r8)
         !  (:) ]  ECM uptake used C (gC/m2/s)
         npp_Nfix               => cnveg_carbonflux_inst%npp_Nfix_patch                 , & ! Output:  [real(r8)
         !  (:) ]  Symbiotic BNF used C (gC/m2/s)
         npp_Nretrans           => cnveg_carbonflux_inst%npp_Nretrans_patch             , & ! Output:  [real(r8)
         !  (:) ]  Retranslocation N uptake used C (gC/m2/s)
         npp_Nuptake            => cnveg_carbonflux_inst%npp_Nuptake_patch              , & ! Output:  [real(r8)
         !  (:) ]  Total N uptake of FUN used C (gC/m2/s)
         npp_growth             => cnveg_carbonflux_inst%npp_growth_patch               , & ! Output:  [real(r8)
         !  (:) ]  Total N uptake of FUN used C (gC/m2/s) 
         burnedoff_carbon       => cnveg_carbonflux_inst%npp_burnedoff_patch            , & ! Output:  [real(r8)
         !  (:) ]  C  that cannot be used for N uptake(gC/m2/s)   
         leafc_change           => cnveg_carbonflux_inst%leafc_change_patch             , & ! Output:  [real(r8)
         !  (:) ]  Used C from the leaf (gC/m2/s)
         leafn_storage_to_xfer  => cnveg_nitrogenflux_inst%leafn_storage_to_xfer_patch  , & ! Output:  [real(r8) (:) ]
         plant_ndemand          => cnveg_nitrogenflux_inst%plant_ndemand_patch          , & ! Iutput:  [real(r8) (:)
         !  ]  N flux required to support initial GPP (gN/m2/s)
         plant_ndemand_retrans  => cnveg_nitrogenflux_inst%plant_ndemand_retrans_patch  , & ! Output:  [real(r8) (:)
         !  ]  N demand generated for FUN (gN/m2/s)
         plant_ndemand_season   => cnveg_nitrogenflux_inst%plant_ndemand_season_patch   , & ! Output:  [real(r8) (:)
         !  ]  N demand for seasonal deciduous forest (gN/m2/s)
         plant_ndemand_stress   => cnveg_nitrogenflux_inst%plant_ndemand_stress_patch   , & ! Output:  [real(r8) (:)
         !  ]  N demand for stress deciduous forest   (gN/m2/s)
         Nactive                => cnveg_nitrogenflux_inst%Nactive_patch                , & ! Output:  [real(r8) (:)
         !  ]  Mycorrhizal N uptake (gN/m2/s)
         Nnonmyc                => cnveg_nitrogenflux_inst%Nnonmyc_patch                , & ! Output:  [real(r8) (:)
         !  ]  Non-mycorrhizal N uptake (gN/m2/s)
         Nam                    => cnveg_nitrogenflux_inst%Nam_patch                    , & ! Output:  [real(r8) (:) ]  AM
         !  uptake (gN/m2/s)
         Necm                   => cnveg_nitrogenflux_inst%Necm_patch                   , & ! Output:  [real(r8) (:) ]  ECM
         !  uptake (gN/m2/s)
         Nactive_no3            => cnveg_nitrogenflux_inst%Nactive_no3_patch            , & ! Output:  [real(r8) (:)
         !  ]  Mycorrhizal N uptake (gN/m2/s)
         Nnonmyc_no3            => cnveg_nitrogenflux_inst%Nnonmyc_no3_patch            , & ! Output:  [real(r8) (:)
         !  ]  Non-mycorrhizal N uptake (gN/m2/s)
         Nam_no3                => cnveg_nitrogenflux_inst%Nam_no3_patch                , & ! Output:  [real(r8) (:)
         !  ]  AM uptake (gN/m2/s)
         Necm_no3               => cnveg_nitrogenflux_inst%Necm_no3_patch               , & ! Output:  [real(r8) (:)
         !  ]  ECM uptake (gN/m2/s)
         Nactive_nh4            => cnveg_nitrogenflux_inst%Nactive_nh4_patch            , & ! Output:  [real(r8) (:)
         !  ]  Mycorrhizal N uptake (gN/m2/s)
         Nnonmyc_nh4            => cnveg_nitrogenflux_inst%Nnonmyc_nh4_patch            , & ! Output:  [real(r8) (:)
         !  ]  Non-mycorrhizal N uptake (gN/m2/s)
         Nam_nh4                => cnveg_nitrogenflux_inst%Nam_nh4_patch                , & ! Output:  [real(r8) (:)
         !  ]  AM uptake (gN/m2/s)
         Necm_nh4               => cnveg_nitrogenflux_inst%Necm_nh4_patch               , & ! Output:  [real(r8) (:)
         !  ]  ECM uptake (gN/m2/s)
         Npassive               => cnveg_nitrogenflux_inst%Npassive_patch               , & ! Output:  [real(r8) (:)
         !  ]  Passive N uptake (gN/m2/s)
         Nfix                   => cnveg_nitrogenflux_inst%Nfix_patch                   , & ! Output:  [real(r8) (:) ]
         !  Symbiotic BNF (gN/m2/s)
         cost_nfix              => cnveg_nitrogenflux_inst%cost_Nfix_patch              , & ! Output:  [real(r8) (:)
         !  ]  Cost of fixation gC:gN
         cost_nactive           => cnveg_nitrogenflux_inst%cost_Nactive_patch        , & ! Output:  [real(r8) (:) ]
         !  Cost of active uptake gC:gN         
         cost_nretrans          => cnveg_nitrogenflux_inst%cost_Nretrans_patch      , & ! Output:  [real(r8) (:) ]
         !  Cost of retranslocation gC:gN         
         nuptake_npp_fraction_patch => cnveg_nitrogenflux_inst%nuptake_npp_fraction_patch    , & ! Output:  [real(r8) (:)
         !  ]  frac of NPP in NUPTAKE 
            
         c_allometry            => cnveg_state_inst%c_allometry_patch                   , & ! Output: [real(r8) (:)   ]  C
         !  allocation index (DIM)                
         n_allometry            => cnveg_state_inst%n_allometry_patch                   , & ! Output: [real(r8) (:)   ]  N
         !  allocation index (DIM)                
         leafn_storage          => cnveg_nitrogenstate_inst%leafn_storage_patch         , & ! Input:  [real(r8) (:)
         !  ]  (gN/m2) leaf N store
         nfix_to_sminn          => soilbiogeochem_nitrogenflux_inst%nfix_to_sminn_col   , & ! Output:  [real(r8) (:)]
         !  symbiotic/asymbiotic N fixation to soil mineral N (gN/m2
         !  /s)
         Nretrans               => cnveg_nitrogenflux_inst%Nretrans_patch               , & ! Output:  [real(r8) (:)
         !  ]  Retranslocation N uptake (gN/m2/s)
         Nretrans_season        => cnveg_nitrogenflux_inst%Nretrans_season_patch        , & ! Output:  [real(r8) (:)
         !  ]  Retranslocation N uptake (gN/m2/s)
         Nretrans_stress        => cnveg_nitrogenflux_inst%Nretrans_stress_patch        , & ! Output:  [real(r8) (:)
         !  ]  Retranslocation N uptake (gN/m2/s)
         Nuptake                => cnveg_nitrogenflux_inst%Nuptake_patch                , & ! Output:  [real(r8) (:)
         !  ]  Total N uptake of FUN (gN/m2/s)
         retransn_to_npool      => cnveg_nitrogenflux_inst%retransn_to_npool_patch               , & ! Output: [real(r8)
         !  (:)   ]  deployment of retranslocated N (gN/m2/s)
         free_retransn_to_npool => cnveg_nitrogenflux_inst%free_retransn_to_npool_patch          , & ! Output: [real(r8)
         ! uptake of free N from leaves (needed to allow RT during the night with no NPP
         sminn_to_plant_fun     => cnveg_nitrogenflux_inst%sminn_to_plant_fun_patch              , & ! Output:
         !  [real(r8) (:) ]  Total soil N uptake of FUN (gN/m2/s)
         sminn_to_plant_fun_vr  => cnveg_nitrogenflux_inst%sminn_to_plant_fun_vr_patch           , & ! Output:
         !  [real(r8) (:) ]  Total layer soil N uptake of FUN (gN/m2
         !  /s) 
         sminn_to_plant_fun_no3_vr  => cnveg_nitrogenflux_inst%sminn_to_plant_fun_no3_vr_patch   , & ! Output:  [real(r8)
         !  (:) ]  Total layer no3 uptake of FUN (gN/m2/s)
         sminn_to_plant_fun_nh4_vr  => cnveg_nitrogenflux_inst%sminn_to_plant_fun_nh4_vr_patch   , & ! Output:  [real(r8)
         !  (:) ]  Total layer nh4 uptake of FUN (gN/m2/s)
         sminn_to_plant_vr      => soilbiogeochem_nitrogenflux_inst%sminn_to_plant_vr_col        , & ! Output:  [real(r8) (:
         ! ,:) ]
         smin_no3_to_plant_vr   => soilbiogeochem_nitrogenflux_inst%smin_no3_to_plant_vr_col     , & ! Output:  [real(r8) (:
         ! ,:) ]
         smin_nh4_to_plant_vr   => soilbiogeochem_nitrogenflux_inst%smin_nh4_to_plant_vr_col     , & ! Output:  [real(r8) (:
         ! ,:) ]
         smin_vr_nh4                  => soilbiogeochem_nitrogenstate_inst%smin_nh4_vr_col             , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral
         !  NH4              
         smin_vr_no3                  => soilbiogeochem_nitrogenstate_inst%smin_no3_vr_col             , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral
         !  NO3              
         soilc_change           => cnveg_carbonflux_inst%soilc_change_patch               , & ! Output:  [real(r8)
         !  (:) ]  Used C from the soil (gC/m2/s)
         h2osoi_liq             => waterstatebulk_inst%h2osoi_liq_col                                , & ! Input:   [real(r8) (:,:)]
         !   liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
         qflx_tran_veg          => waterfluxbulk_inst%qflx_tran_veg_patch                            , & ! Input:   [real(r8) (:)  ]
         !   vegetation transpiration (mm H2O/s) (+ = to atm)
         t_soisno               => temperature_inst%t_soisno_col                                 , & ! Input:   [real(r8) (:,:)]
         !   soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         crootfr                => soilstate_inst%crootfr_patch                                    & ! Input:   [real(r8) (:,:)]
         !   fraction of roots for carbon in each soil layer  (nlevgrnd)
         )
  !--------------------------------------------------------------------
  !-----------
  ! Initialize output fluxes, which were also initialized in CNFUNMod.
  !--------------------------------------------------------------------
  !-----------
  local_use_flexibleCN            = use_flexibleCN
  steppday                        = 48._r8
  qflx_tran_veg_layer             = 0._r8
  rootc_dens_step                 = 0._r8
  plant_ndemand_pool              = 0._r8

  call t_startf('CNFUNzeroarrays')
  do fp = 1,num_soilp        ! PFT Starts
     p = filter_soilp(fp)
     availc_pool(p)                  = 0._r8
     rootC(p)                        = 0._r8
     litterfall_n(p)                 = 0._r8
     burnedoff_carbon(p)             = 0._r8
  end do


  do j = 1, nlevdecomp
     do fp = 1,num_soilp        ! PFT Starts
        p = filter_soilp(fp)
        c = patch%column(p)
        rootc_dens(p,j)                 = 0._r8 
        cost_retran(p,j)                = 0._r8
        cost_fix(p,j)                   = 0._r8
        cost_resis(p,j)                 = 0._r8
        cost_res_resis(p,j)             = 0._r8
        cost_active_no3(p,j)            = 0._r8
        cost_active_nh4(p,j)            = 0._r8
        cost_nonmyc_no3(p,j)            = 0._r8
        cost_nonmyc_nh4(p,j)            = 0._r8
   
        sminn_no3_conc(c,j)             = 0._r8
        sminn_no3_layer(c,j)            = 0._r8
        sminn_nh4_conc(c,j)             = 0._r8
        sminn_nh4_layer(c,j)            = 0._r8
     end do
  end do

  do istp = 1, nstp
     do fp = 1,num_soilp        ! PFT Starts
        p = filter_soilp(fp)
        npp_remaining(p,istp)           = 0._r8
        permyc(p,istp)                  = 0._r8
        plant_ndemand_pool_step(p,istp) = 0._r8
        nt_uptake(p,istp)               = 0._r8
        npp_uptake(p,istp)              = 0._r8
        leafn_step(p,istp)              = 0._r8
        leafn_retrans_step(p,istp)      = 0._r8
        litterfall_n_step(p,istp)       = 0._r8
        litterfall_c_step(p,istp)       = 0._r8
     end do
     do j = 1, nlevdecomp
        do fp = 1,num_soilp        ! PFT Starts
           p = filter_soilp(fp)
           sminn_no3_conc_step(p,j,istp)      = 0._r8
           sminn_no3_layer_step(p,j,istp)     = 0._r8
           sminn_no3_uptake(p,j,istp)         = 0._r8
           sminn_nh4_conc_step(p,j,istp)      = 0._r8
           sminn_nh4_layer_step(p,j,istp)     = 0._r8
           sminn_nh4_uptake(p,j,istp)         = 0._r8
        end do
     end do
  end do

  do icost = 1, ncost6
     do j = 1, nlevdecomp
        costNit(j,icost)                    = big_cost 
     end do
  end do
  
  ! Time step of FUN
  dt           =  real(get_step_size(), r8)
  call t_stopf('CNFUNzeroarrays')
  !--------------------------------------------------------------------
  !----------------------------
  ! Calculation starts
  !--------------------------------------------------------------------
  call t_startf('CNFUNcalcs1')
  !----------------------------
  do fp = 1,num_soilp        ! PFT Starts
     p = filter_soilp(fp)

     litterfall_n(p) =  (leafc_to_litter_fun(p) / leafcn_offset(p))  * dt
     rootC(p)        =  frootc(p)

     plantN(p)       =  leafn(p) + frootn(p) + livestemn(p) + livecrootn(p)
     if (n_allometry(p).gt.0._r8) then 
         plantCN(p)  = c_allometry(p)/n_allometry(p) !changed RF.
         ! above code gives CN ratio too low. 
     else
         plantCN(p)  = 0._r8 
     end if
  end do   ! PFT ends
  do istp = 1, nstp
     do fp = 1,num_soilp        ! PFT Starts
        p = filter_soilp(fp)

        if (istp.eq.ecm_step) then
           permyc(p,istp)      = perecm(ivt(p))
           kc_active(p,istp)   = ekc_active(ivt(p))
           kn_active(p,istp)   = ekn_active(ivt(p))
        else
           permyc(p,istp)      = 1._r8 - perecm(ivt(p))
           kc_active(p,istp)   = akc_active(ivt(p))
           kn_active(p,istp)   = akn_active(ivt(p))
        end if

        if(leafc(p)>0.0_r8)then
           ! N available in leaf which fell off in this timestep. Same fraction loss as C.    
           litterfall_c_step(p,istp)         =   dt * permyc(p,istp) * leafc_to_litter_fun(p) 
           litterfall_n_step(p,istp)         =   dt * permyc(p,istp) * leafn(p) * leafc_to_litter_fun(p)/leafc(p) 
        endif 

        if (season_decid(ivt(p)) == 1._r8.or.stress_decid(ivt(p)) == 1._r8) then
          if (offset_flag(p) .ne. 1._r8) then
            litterfall_n_step(p,istp) = 0.0_r8      
            litterfall_c_step(p,istp) = 0.0_r8      
          endif
        endif

     end do
  end do     
      
  do j = 1, nlevdecomp
     do fp = 1,num_soilp        ! PFT Starts
       p = filter_soilp(fp)
       c = patch%column(p)
       sminn_no3_layer(c,j)= smin_no3_to_plant_vr(c,j) * dzsoi_decomp(j) * dt
       sminn_nh4_layer(c,j)= smin_nh4_to_plant_vr(c,j) * dzsoi_decomp(j) * dt
       if (h2osoi_liq(c,j) < smallValue) then
          sminn_no3_layer(c,j) = 0._r8
          sminn_nh4_layer(c,j) = 0._r8
       end if
       sminn_no3_layer(c,j)    = max(sminn_no3_layer(c,j),0._r8)
       sminn_nh4_layer(c,j)    = max(sminn_nh4_layer(c,j),0._r8)
       if (h2osoi_liq(c,j) > smallValue) then
          sminn_no3_conc(c,j)  = sminn_no3_layer(c,j) / (h2osoi_liq(c,j) * 1000._r8) ! (gN/m2)/(gH2O/m2) (coverted from
          !  kg2g)
          sminn_nh4_conc(c,j)  = sminn_nh4_layer(c,j) / (h2osoi_liq(c,j) * 1000._r8) ! (gN/m2)/(gH2O/m2) (coverted from
          !  kg2g)
       else
          sminn_no3_conc(c,j)  = 0._r8
          sminn_nh4_conc(c,j)  = 0._r8
       end if
     end do
  end do

  do istp = 1, nstp
     do j = 1, nlevdecomp
        do fp = 1,num_soilp        ! PFT Starts
           p = filter_soilp(fp)
           c = patch%column(p)

           sminn_no3_layer_step(p,j,istp)  =   sminn_no3_layer(c,j) * permyc(p,istp)
           sminn_nh4_layer_step(p,j,istp)  =   sminn_nh4_layer(c,j) * permyc(p,istp)
           sminn_no3_conc_step(p,j,istp)   =   sminn_no3_conc(c,j)  * permyc(p,istp)
           sminn_nh4_conc_step(p,j,istp)   =   sminn_nh4_conc(c,j)  * permyc(p,istp)
        end do
     end do
  end do
  call t_stopf('CNFUNcalcs1')

  call t_startf('CNFUNzeroarrays2')
  do fp = 1,num_soilp        ! PFT Starts
     p = filter_soilp(fp)
     n_passive_acc(p)                = 0._r8
     n_fix_acc_total(p)              = 0._r8
     n_retrans_acc_total(p)          = 0._r8
     npp_fix_acc_total(p)            = 0._r8
     n_nonmyc_no3_retrans_total(p)   = 0._r8
     n_nonmyc_nh4_retrans_total(p)   = 0._r8
     npp_retrans_acc_total(p)        = 0._r8
     n_am_no3_acc(p)                 = 0._r8
     n_am_nh4_acc(p)                 = 0._r8
     n_am_no3_retrans(p)             = 0._r8
     n_am_nh4_retrans(p)             = 0._r8
     n_ecm_no3_acc(p)                = 0._r8
     n_ecm_nh4_acc(p)                = 0._r8
     n_ecm_no3_retrans(p)            = 0._r8
     n_ecm_nh4_retrans(p)            = 0._r8
     n_active_no3_acc_total(p)       = 0._r8 
     n_active_nh4_acc_total(p)       = 0._r8
     n_active_no3_retrans_total(p)   = 0._r8 
     n_active_nh4_retrans_total(p)   = 0._r8
     n_nonmyc_no3_acc_total(p)       = 0._r8
     n_nonmyc_nh4_acc_total(p)       = 0._r8
     npp_active_no3_acc_total(p)     = 0._r8
     npp_active_nh4_acc_total(p)     = 0._r8
     npp_active_no3_retrans_total(p) = 0._r8
     npp_active_nh4_retrans_total(p) = 0._r8
     npp_nonmyc_no3_acc_total(p)     = 0._r8
     npp_nonmyc_nh4_acc_total(p)     = 0._r8
     npp_nonmyc_no3_retrans_total(p) = 0._r8
     npp_nonmyc_nh4_retrans_total(p) = 0._r8
     free_Nretrans(p)                = 0._r8
  end do

  do j = 1, nlevdecomp
     do fp = 1,num_soilp        ! PFT Starts
        p = filter_soilp(fp)
        n_passive_no3_vr(p,j)           = 0._r8
        n_passive_nh4_vr(p,j)           = 0._r8
        n_active_no3_vr(p,j)            = 0._r8
        n_nonmyc_no3_vr(p,j)            = 0._r8
        n_active_nh4_vr(p,j)            = 0._r8
        n_nonmyc_nh4_vr(p,j)            = 0._r8
     end do
  end do
  do istp = 1, nstp
     do fp = 1,num_soilp        ! PFT Starts
        p = filter_soilp(fp)
        n_passive_step(p,istp)          = 0._r8
        n_fix_acc(p,istp)               = 0._r8
        n_retrans_acc(p,istp)           = 0._r8
        npp_fix_acc(p,istp)             = 0._r8
        npp_retrans_acc(p,istp)         = 0._r8
        n_active_no3_acc(p,istp)        = 0._r8   
        n_active_nh4_acc(p,istp)        = 0._r8  
        n_active_no3_retrans(p,istp)    = 0._r8
        n_active_nh4_retrans(p,istp)    = 0._r8
        n_nonmyc_no3_acc(p,istp)        = 0._r8  
        n_nonmyc_nh4_acc(p,istp)        = 0._r8  
        n_nonmyc_no3_retrans(p,istp)    = 0._r8
        n_nonmyc_nh4_retrans(p,istp)    = 0._r8
        npp_active_no3_acc(p,istp)      = 0._r8
        npp_active_nh4_acc(p,istp)      = 0._r8
        npp_active_no3_retrans(p,istp)  = 0._r8
        npp_active_nh4_retrans(p,istp)  = 0._r8
        npp_nonmyc_no3_acc(p,istp)      = 0._r8
        npp_nonmyc_nh4_acc(p,istp)      = 0._r8
        npp_nonmyc_no3_retrans(p,istp)  = 0._r8
        npp_nonmyc_nh4_retrans(p,istp)  = 0._r8
     end do
  end do

  burned_off_carbon               = 0._r8 
  call t_stopf('CNFUNzeroarrays2')


  call t_startf('CNFUNcalcs')
pft:do fp = 1,num_soilp        ! PFT Starts
      p = filter_soilp(fp)
      c = patch%column(p)
      excess_carbon_acc               = 0.0_r8
      burned_off_carbon               = 0.0_r8
     
      sminn_to_plant_fun_nh4_vr(p,:) = 0._r8
      sminn_to_plant_fun_no3_vr(p,:) = 0._r8      
     
      ! I have turned off this r etranslocation functionality for now. To
      !  be rolled back in to a new version later on once the rest of
      !   th
      ! mode is working OK. RF

      if (season_decid(ivt(p)) == 1._r8.or.stress_decid(ivt(p)) == 1._r8) then
         if (onset_flag(p) == 1._r8) then
            leafc_storage_xfer_acc(p) = leafc_storage_xfer_acc(p) + leafc_storage_to_xfer(p) * dt
            leafn_storage_xfer_acc(p) = leafn_storage_xfer_acc(p) + leafn_storage_to_xfer(p) * dt
         end if
         if (offset_flag(p) == 1._r8) then
            storage_cdemand(p)        = leafc_storage(p)          / (ndays_off * steppday)
            storage_ndemand(p)        = leafn_storage_xfer_acc(p) / (ndays_off * steppday)
            storage_ndemand(p)        = max(storage_ndemand(p),0._r8)
         else
            storage_cdemand(p)        = 0._r8    
            storage_ndemand(p)        = 0._r8   
         end if
      else
          storage_cdemand(p)          = 0._r8
          storage_ndemand(p)          = 0._r8 
      end if   ! end for deciduous

      !---------How much carbon is provided, to be used for either growth
      ! or Nitrogen uptake?-------------------
      availc_pool(p)            =  availc(p)        *  dt

      if (availc_pool(p) > 0._r8) then
         do j = 1, nlevdecomp
            rootc_dens(p,j)     =  crootfr(p,j) * rootC(p)
         end do
      end if

      plant_ndemand_pool(p)     =  plant_ndemand(p) *  dt
      plant_ndemand_pool(p)     =  max(plant_ndemand_pool(p),0._r8)
      plant_ndemand_retrans(p)  =  storage_ndemand(p)

      !--------------------------------------------------------------------
      !----------
stp:  do istp = ecm_step, am_step        ! TWO STEPS
         retrans_limit1              = 0._r8
         dn                          = 0._r8
         dnpp                        = 0._r8
      
         ! zero out all of the fluxes that get accumulated accross ISTP 
         sminn_no3_diff              = 0._r8
         sminn_nh4_diff              = 0._r8
         active_no3_limit1           = 0._r8
         active_nh4_limit1           = 0._r8

         
         n_from_active_no3(:)        = 0.0_r8
         n_from_active_nh4(:)        = 0.0_r8
         n_from_nonmyc_no3(:)        = 0.0_r8
         n_from_nonmyc_nh4(:)        = 0.0_r8
         n_from_fixation(:)          = 0.0_r8
         n_from_retrans(:)           = 0.0_r8
         
         n_active_no3_acc(p,istp)       = 0.0_r8
         n_active_nh4_acc(p,istp)       = 0.0_r8
         n_nonmyc_no3_acc(p,istp)       = 0.0_r8
         n_nonmyc_nh4_acc(p,istp)       = 0.0_r8
         n_fix_acc(p,istp)              = 0.0_r8
         n_retrans_acc(p,istp)          = 0.0_r8
         free_nretrans_acc(p,istp)     = 0.0_r8
 
         npp_active_no3_acc(p,istp)     = 0.0_r8
         npp_active_nh4_acc(p,istp)     = 0.0_r8
         npp_nonmyc_no3_acc(p,istp)     = 0.0_r8
         npp_nonmyc_no3_acc(p,istp)     = 0.0_r8
         npp_fix_acc(p,istp)            = 0.0_r8
         npp_retrans_acc(p,istp)        = 0.0_r8
         
         npp_to_active_no3(:)        = 0.0_r8
         npp_to_active_nh4(:)        = 0.0_r8
         npp_to_nonmyc_no3(:)        = 0.0_r8
         npp_to_nonmyc_nh4(:)        = 0.0_r8
         npp_to_fixation(:)          = 0.0_r8
         npp_to_retrans(:)           = 0.0_r8
     
  
      
         unmetDemand              = .TRUE.
         plant_ndemand_pool_step(p,istp)   = plant_ndemand_pool(p)    * permyc(p,istp) 
         npp_remaining(p,istp)             = availc_pool(p)           * permyc(p,istp)
         
  
         ! if (plant_ndemand_pool_step(p,istp) .gt. 0._r8) then   !
            !  plant_ndemand_pool_step > 0.0
            
            do j = 1, nlevdecomp
               tc_soisno(c,j)          = t_soisno(c,j)  -   tfrz
               if(pftcon%c3psn(patch%itype(p)).eq.1)then
                 fixer=1
               else
                 fixer=0
               endif
               costNit(j,icostFix)     = fun_cost_fix(fixer,a_fix(ivt(p)),b_fix(ivt(p))&
               ,c_fix(ivt(p)) ,big_cost,crootfr(p,j),s_fix(ivt(p)),tc_soisno(c,j))
            end do
            cost_fix(p,1:nlevdecomp)      = costNit(:,icostFix)
            
             
            !--------------------------------------------------------------------
            !------------
            !         If passive uptake is insufficient, consider fixation,
            !          mycorrhizal 
            !         non-mycorrhizal, storage, and retranslocation.
            !--------------------------------------------------------------------
            !------------
            !--------------------------------------------------------------------
            !------------
            !          Costs of active uptake.
            !--------------------------------------------------------------------
            !------------
            !------Mycorrhizal Uptake Cost-----------------!
            do j = 1,nlevdecomp
               rootc_dens_step            = rootc_dens(p,j) *  permyc(p,istp)
               costNit(j,icostActiveNO3)  = fun_cost_active(sminn_no3_layer_step(p,j,istp) &
               ,big_cost,kc_active(p,istp),kn_active(p,istp) ,rootc_dens_step,crootfr(p,j),smallValue)
               costNit(j,icostActiveNH4)  = fun_cost_active(sminn_nh4_layer_step(p,j,istp) &
               ,big_cost,kc_active(p,istp),kn_active(p,istp) ,rootc_dens_step,crootfr(p,j),smallValue)
            end do
            cost_active_no3(p,1:nlevdecomp)  = costNit(:,icostActiveNO3) 
            cost_active_nh4(p,1:nlevdecomp)  = costNit(:,icostActiveNH4)
            

            !------Non-mycorrhizal Uptake Cost-------------!
            do j = 1,nlevdecomp
               rootc_dens_step             = rootc_dens(p,j)  *  permyc(p,istp)
               costNit(j,icostnonmyc_no3)   = fun_cost_nonmyc(sminn_no3_layer_step(p,j,istp) &
               ,big_cost,kc_nonmyc(ivt(p)),kn_nonmyc(ivt(p)) ,rootc_dens_step,crootfr(p,j),smallValue)
               costNit(j,icostnonmyc_nh4)   = fun_cost_nonmyc(sminn_nh4_layer_step(p,j,istp) &
               ,big_cost,kc_nonmyc(ivt(p)),kn_nonmyc(ivt(p)) ,rootc_dens_step,crootfr(p,j),smallValue)
            end do
            cost_nonmyc_no3(p,1:nlevdecomp)   = costNit(:,icostnonmyc_no3)
            cost_nonmyc_nh4(p,1:nlevdecomp)   = costNit(:,icostnonmyc_nh4)
            

            ! Remove C required to pair with N from passive uptake
            !  from the available pool. 
            npp_remaining(p,istp)  =   npp_remaining(p,istp) - n_passive_step(p,istp)*plantCN(p)
              
fix_loop:   do FIX =plants_are_fixing, plants_not_fixing !loop around percentages of fixers and non
               ! fixers, with differnt costs. 
               if(FIX==plants_are_fixing)then ! How much of the carbon in this PFT can in principle be used for fixation? 
                 ! This is analagous to fixing the % of fixers for a given PFT - may not be realistic in the long run
                 ! but prevents wholesale switching to fixer dominance during e.g. CO2 fertilization.  
                 fixerfrac = FUN_fracfixers(ivt(p))
               else
                 fixerfrac = 1.0_r8 - FUN_fracfixers(ivt(p))
               endif 
               npp_to_spend = npp_remaining(p,istp)  * fixerfrac !put parameter here.



               n_from_active_no3(1:nlevdecomp) = 0._r8
               n_from_active_nh4(1:nlevdecomp) = 0._r8
               n_from_nonmyc_no3(1:nlevdecomp) = 0._r8
               n_from_nonmyc_nh4(1:nlevdecomp) = 0._r8
               !--------------------------------------------------------------------
               !-----------
               !           Calculate Integrated Resistance OF WHOLE SOIL COLUMN
               !--------------------------------------------------------------------
               !----------- 

               sum_n_acquired      = 0.0_r8
               total_N_conductance = 0.0_r8
               do j = 1, nlevdecomp
                  !----------!
                  ! Method changed from FUN-resistors method to a method which 
                  ! allocates fluxs based on conductance. rosief
                  !----------!
             
                  ! Sum the conductances             
                  total_N_conductance  = total_N_conductance + 1._r8/ &
                                         cost_active_no3(p,j) + 1._r8/cost_active_nh4(p,j) &
                                         + 1._r8/cost_nonmyc_no3(p,j)      &
                                         + 1._r8/cost_nonmyc_nh4(p,j) 
                  if(FIX==plants_are_fixing)then
                      total_N_conductance  = total_N_conductance  + 1.0_r8 * 1._r8/cost_fix(p,j)
                  end if 
                      
                end do 
             
                do j = 1, nlevdecomp     
                  ! Calculate npp allocation to pathways proportional to their exchange rate (N/C) 
                
                  npp_frac_to_active_nh4(j) = (1._r8/cost_active_nh4(p,j)) / total_N_conductance
                  npp_frac_to_nonmyc_nh4(j) = (1._r8/cost_nonmyc_nh4(p,j)) / total_N_conductance
                  npp_frac_to_active_no3(j) = (1._r8/cost_active_no3(p,j)) / total_N_conductance
                  npp_frac_to_nonmyc_no3(j) = (1._r8/cost_nonmyc_no3(p,j)) / total_N_conductance
                  if(FIX==plants_are_fixing)then
                    npp_frac_to_fixation(j)   = (1.0_r8 * 1._r8/cost_fix(p,j)) / total_N_conductance
                  else
                    npp_frac_to_fixation(j)   = 0.0_r8 
                  end if
                     
                  ! Calculate hypothetical N uptake from each source   
                  if(FIX==plants_are_fixing)then
                    n_exch_fixation(j)   = npp_frac_to_fixation(j)   / cost_fix(p,j)
                  else
                    n_exch_fixation(j)   = 0.0_r8 
                  end if                   
              
                  n_exch_active_nh4(j) = npp_frac_to_active_nh4(j) / cost_active_nh4(p,j) 
                  n_exch_nonmyc_nh4(j) = npp_frac_to_nonmyc_nh4(j) / cost_nonmyc_nh4(p,j) 
                  n_exch_active_no3(j) = npp_frac_to_active_no3(j) / cost_active_no3(p,j) 
                  n_exch_nonmyc_no3(j) = npp_frac_to_nonmyc_no3(j) / cost_nonmyc_no3(p,j) 
               
                  ! Total N aquired from one unit of carbon  (N/C)
                  sum_n_acquired        =  sum_n_acquired  + n_exch_active_nh4(j) +&
                                        n_exch_nonmyc_nh4(j)+ n_exch_active_no3(j) + n_exch_nonmyc_no3(j)
                                          
                  if(FIX==plants_are_fixing)then
                    sum_n_acquired= sum_n_acquired +  n_exch_fixation(j)
                  end if 
                                                                           
               end do !nlevdecomp
            
               total_N_resistance = 1.0_r8/sum_n_acquired

               !-------------------------------------------------------------------------------
               !           Calculate appropriate degree of retranslocation
               !-------------------------------------------------------------------------------
      
               if(leafc(p).gt.0.0_r8.and.litterfall_n_step(p,istp)* fixerfrac>0.0_r8.and.ivt(p) <npcropmin)then
                  call fun_retranslocation(p,dt,npp_to_spend,&
                                litterfall_c_step(p,istp)* fixerfrac,&
                                litterfall_n_step(p,istp)* fixerfrac,&
                                total_n_resistance, total_c_spent_retrans,total_c_accounted_retrans, &
                                free_n_retrans,paid_for_n_retrans, leafcn(ivt(p)), & 
                                grperc(ivt(p)), plantCN(p))
                                 
               else
                   total_c_accounted_retrans = 0.0_r8
                   total_c_spent_retrans     = 0.0_r8
                   total_c_accounted_retrans = 0.0_r8
                   paid_for_n_retrans        = 0.0_r8
                   free_n_retrans            = 0.0_r8
               endif
  
               !---------- add retrans fluxes in to total budgets. --------
               ! remove C from available pool, both directly spent and accounted for by N uptake

               npp_to_spend  = npp_to_spend - total_c_spent_retrans - total_c_accounted_retrans
               npp_retrans_acc(p,istp) = npp_retrans_acc(p,istp) + total_c_spent_retrans 
               ! add to to C spent pool                              
               n_retrans_acc(p,istp)     = n_retrans_acc(p,istp)     + paid_for_n_retrans
               free_nretrans_acc(p,istp) = free_nretrans_acc(p,istp) + free_n_retrans
               ! add N to the acquired from retrans pool 
  
            
               !-------------------------------------------------------------------------------
               !           Spend C on extracting N.
               !-------------------------------------------------------------------------------
               if (plant_ndemand_pool_step(p,istp) .gt. 0._r8) then    ! unmet demand
  
                 if(local_use_flexiblecn)then   
                     if (leafn(p) == 0.0_r8) then   ! to avoid division by zero
                       delta_CN = fun_cn_flex_c(ivt(p))   ! Max CN ratio over standard
                     else
                       delta_CN = (leafc(p)+leafc_storage(p))/(leafn(p)+leafn_storage(p)) - leafcn(ivt(p)) ! leaf CN ratio                                                              
                     end if
                     ! C used for uptake is reduced if the cost of N is very high                
                     frac_ideal_C_use = max(0.0_r8,1.0_r8 - (total_N_resistance-fun_cn_flex_a(ivt(p)))/fun_cn_flex_b(ivt(p)) )
                     ! then, if the plant is very much in need of N, the C used for uptake is increased accordingly.                  
                     if(delta_CN .gt.0.and. frac_ideal_C_use.lt.1.0)then           
                       frac_ideal_C_use = frac_ideal_C_use + (1.0_r8-frac_ideal_C_use)*min(1.0_r8, delta_CN/fun_cn_flex_c(ivt(p)))
                     end if    
                     ! If we have too much N (e.g. from free N retranslocation) then make frac_ideal_c_use even lower.    
                     ! For a CN delta of fun_cn_flex_c, then we reduce C expendiure to the minimum of 0.5. 
                     ! This seems a little intense? 
                     if(delta_CN.lt.0.0)then
                        frac_ideal_C_use = frac_ideal_C_use + 0.5_r8*(1.0_r8*delta_CN/fun_cn_flex_c(ivt(p)))
                     endif 
                     frac_ideal_C_use = max(min(1.0_r8,frac_ideal_C_use),0.5_r8) 
                     ! don't let this go above 1 or below an arbirtray minimum (to prevent zero N uptake). 
                 else
                     frac_ideal_C_use= 1.0_r8
                 end if
               
                 excess_carbon                 = npp_to_spend * (1.0_r8-frac_ideal_c_use)
                 if(excess_carbon*(1.0_r8+grperc(ivt(p))).gt.npp_to_spend)then !prevent negative dnpp
                      excess_carbon =  npp_to_spend/(1.0_r8+grperc(ivt(p)))
                 endif
                 excess_carbon_acc             = excess_carbon_acc + excess_carbon

                 ! spend less C than you have to to meet the target, thus allowing C:N ratios to rise. 
                 npp_to_spend         = npp_to_spend - excess_carbon*(1.0_r8+grperc(ivt(p)))
            
                 ! This is the main equation of FUN, which figures out how much C to spend on uptake to remain at the target CN ratio. 
                 ! nb. This term assumes that cost of N is constant through the timestep, because we don't have 
                 ! a concept of Michealis Menten kinetics. 
                 ! 
                 !Calculate the hypothetical amount of NPP that we should use to extract N over whole profile
                 !This calculation is based on the simulataneous solution of the uptake and extrction N balance. 
                 !It satisfies the criteria (spentC+growthC=availC AND spentC/cost=growthC/plantCN
                 !Had to add growth respiration here to balance carbon pool. 

               
                 dnpp  = npp_to_spend / ( (1.0_r8+grperc(ivt(p)))*(plantCN(p) / total_N_resistance) + 1._r8)  
                 dnpp  = dnpp * frac_ideal_C_use
           
                 !hypothetical amount of N acquired. 
                 dn    = dnpp / total_N_resistance
                 do j = 1,nlevdecomp
                        
                     ! RF How much of this NPP carbon do we allocate to the different pathways? fraction x gC/m2/s?
                     ! Could this code now be put in a matrix? 

                     npp_to_active_nh4(j) = npp_frac_to_active_nh4(j) * dNPP
                     npp_to_nonmyc_nh4(j) = npp_frac_to_nonmyc_nh4(j) * dNPP
                     npp_to_active_no3(j) = npp_frac_to_active_no3(j) * dNPP
                     npp_to_nonmyc_no3(j) = npp_frac_to_nonmyc_no3(j) * dNPP 
                     
                     if(FIX==plants_are_fixing)then
                       npp_to_fixation(j) = npp_frac_to_fixation(j) * dNPP
                     else
                       npp_to_fixation(j) = 0.0_r8
                     end if    

                     n_from_active_nh4(j) = npp_to_active_nh4(j)  / cost_active_nh4(p,j)
                     n_from_nonmyc_nh4(j) = npp_to_nonmyc_nh4(j)  / cost_nonmyc_nh4(p,j)
                     n_from_active_no3(j) = npp_to_active_no3(j)  / cost_active_no3(p,j)
                     n_from_nonmyc_no3(j) = npp_to_nonmyc_no3(j)  / cost_nonmyc_no3(p,j)
                
                     if(FIX==plants_are_fixing)then
                       n_from_fixation(j) = npp_to_fixation(j)    / cost_fix(p,j)
                     else
                       n_from_fixation(j) = 0.0_r8
                     end if
                                                                     
                 end do
              
                 ! did we exceed the limits of uptake for any of these pools?
                 do j = 1,nlevdecomp    
                  
                       ! --------------------ACTIVE UPTAKE NO3 UPTAKE LIMIT------------------------!  
                       active_no3_limit1          = sminn_no3_layer_step(p,j,istp) * fixerfrac 
                  
                        ! trying to remove too much nh4 from soil. 
                        if (n_from_active_no3(j) + n_from_nonmyc_no3(j).gt.active_no3_limit1) then 
                           sminn_no3_diff          = n_from_active_no3(j) + n_from_nonmyc_no3(j) - active_no3_limit1
                           temp_n_flux = n_from_active_no3(j)
                           ! divide discrepancy between sources
                           n_from_active_no3(j)     = n_from_active_no3(j) - sminn_no3_diff &
                                                      * (n_from_active_no3(j) /(n_from_active_no3(j) + n_from_nonmyc_no3(j)))
                           n_from_nonmyc_no3(j)     = n_from_nonmyc_no3(j) - sminn_no3_diff &
                                                 * (n_from_nonmyc_no3(j) /(temp_n_flux + n_from_nonmyc_no3(j)))
                           npp_to_active_no3(j)     = n_from_active_no3(j) * cost_active_no3(p,j) 
                           npp_to_nonmyc_no3(j)     = n_from_nonmyc_no3(j) * cost_nonmyc_no3(p,j)
                                      
                       end if
                  
                       ! --------------------ACTIVE UPTAKE NH4 UPTAKE LIMIT------------------------!  
                       active_nh4_limit1          = sminn_nh4_layer_step(p,j,istp) *fixerfrac
                  
                      
                       ! trying to remove too much nh4 from soil. 
                       if (n_from_active_nh4(j) + n_from_nonmyc_nh4(j).gt.active_nh4_limit1) then 
                    
                          sminn_nh4_diff          = n_from_active_nh4(j) + n_from_nonmyc_nh4(j) - active_nh4_limit1
                          temp_n_flux = n_from_active_nh4(j)
                          ! divide discrepancy between sources
                          n_from_active_nh4(j)  = n_from_active_nh4(j)    - (sminn_nh4_diff &
                                                   * n_from_active_nh4(j) /(n_from_active_nh4(j) + n_from_nonmyc_nh4(j)))
                          n_from_nonmyc_nh4(j)  = n_from_nonmyc_nh4(j)    - (sminn_nh4_diff &
                                                   * n_from_nonmyc_nh4(j) /(temp_n_flux+ n_from_nonmyc_nh4(j)))
                          npp_to_active_nh4(j)  = n_from_active_nh4(j)    * cost_active_nh4(p,j) 
                          npp_to_nonmyc_nh4(j)  = n_from_nonmyc_nh4(j)    * cost_nonmyc_nh4(p,j)    
                                 
                       end if
                                                               
                       ! How much N did we end up with
                       N_acquired                    =  n_from_active_no3(j)+n_from_nonmyc_no3(j) &
                                                       + n_from_active_nh4(j)+n_from_nonmyc_nh4(j)
                  
                  
                       ! How much did it actually cost? 
                       C_spent                       =   npp_to_active_no3(j)+npp_to_nonmyc_no3(j) &
                                                       + npp_to_active_nh4(j)+npp_to_nonmyc_nh4(j)
                                                  
                       if(FIX==plants_are_fixing)then
                          N_acquired = N_acquired + n_from_fixation(j)
                          C_spent    = C_spent + npp_to_fixation(j) 
                       end if
                  
                  

                       ! How much C did we allocate or spend in this layer? 
                       npp_to_spend                 = npp_to_spend   - C_spent - (N_acquired &
                                                       * plantCN(p)*(1.0_r8+ grperc(ivt(p))))
                  
                       ! Accumulate those fluxes
                       nt_uptake(p,istp)             = nt_uptake(p,istp)       + N_acquired
                       npp_uptake(p,istp)            = npp_uptake(p,istp)      + C_spent
                  
                  
                  
                                                                                 
                       !-------------------- N flux accumulation------------!
                       n_active_no3_acc(p,istp)      = n_active_no3_acc(p,istp) + n_from_active_no3(j)
                       n_active_nh4_acc(p,istp)      = n_active_nh4_acc(p,istp) + n_from_active_nh4(j)
                       n_nonmyc_no3_acc(p,istp)      = n_nonmyc_no3_acc(p,istp) + n_from_nonmyc_no3(j)
                       n_nonmyc_nh4_acc(p,istp)      = n_nonmyc_nh4_acc(p,istp) + n_from_nonmyc_nh4(j)                 
            
                       !-------------------- C flux accumulation------------!
                       npp_active_no3_acc(p,istp) = npp_active_no3_acc(p,istp)  + npp_to_active_no3(j)
                       npp_active_nh4_acc(p,istp) = npp_active_nh4_acc(p,istp)  + npp_to_active_nh4(j)
                       npp_nonmyc_no3_acc(p,istp) = npp_nonmyc_no3_acc(p,istp)  + npp_to_nonmyc_no3(j)
                       npp_nonmyc_nh4_acc(p,istp) = npp_nonmyc_nh4_acc(p,istp)  + npp_to_nonmyc_nh4(j)               
                  
                       if(FIX == plants_are_fixing)then
                         n_fix_acc(p,istp)          = n_fix_acc(p,istp)         + n_from_fixation(j)
                         npp_fix_acc(p,istp)        = npp_fix_acc(p,istp)       + npp_to_fixation(j)
                       end if
                  
                 end do    ! j
             
             
                 ! check that we get the right amount of N...
             
            
             
                 ! Occasionally, the algorithm will want to extract a high fraction of NPP from a pool (eg leaves) that                               
                 ! quickly empties. One solution to this is to iterate round all the calculations starting from                                       
                 ! the cost functions. The other is to burn off the extra carbon and hope this doesn't happen very often...                     

                 if (npp_to_spend .ge. 0.0000000000001_r8)then
                       burned_off_carbon =  burned_off_carbon + npp_to_spend
                 end if
           
             
                      
                 ! add vertical fluxes to patch arrays 
                 do j = 1,nlevdecomp
                       n_active_no3_vr(p,j)      =  n_active_no3_vr(p,j)      + n_from_active_no3(j)
                       n_active_nh4_vr(p,j)      =  n_active_nh4_vr(p,j)      + n_from_active_nh4(j)
                       n_nonmyc_no3_vr(p,j)      =  n_nonmyc_no3_vr(p,j)      + n_from_nonmyc_no3(j)
                       n_nonmyc_nh4_vr(p,j)      =  n_nonmyc_nh4_vr(p,j)      + n_from_nonmyc_nh4(j)  
                 end do
               end if !unmet demand`
            
            end do fix_loop ! FIXER. 
             
            if (istp.eq.ecm_step) then
               n_ecm_no3_acc(p)          =  n_active_no3_acc(p,istp)
               n_ecm_nh4_acc(p)          =  n_active_nh4_acc(p,istp)
            else
               n_am_no3_acc(p)           =  n_active_no3_acc(p,istp)
               n_am_nh4_acc(p)           =  n_active_nh4_acc(p,istp)
            end if
            
            ! Accumulate total column N fluxes over istp
            n_active_no3_acc_total(p)    =  n_active_no3_acc_total(p)    + n_active_no3_acc(p,istp)
            n_active_nh4_acc_total(p)    =  n_active_nh4_acc_total(p)    + n_active_nh4_acc(p,istp)
            n_nonmyc_no3_acc_total(p)    =  n_nonmyc_no3_acc_total(p)    + n_nonmyc_no3_acc(p,istp)
            n_nonmyc_nh4_acc_total(p)    =  n_nonmyc_nh4_acc_total(p)    + n_nonmyc_nh4_acc(p,istp)
            n_fix_acc_total(p)           =  n_fix_acc_total(p)           + n_fix_acc(p,istp)
            n_retrans_acc_total(p)       =  n_retrans_acc_total(p)       + n_retrans_acc(p,istp)
            free_nretrans(p)             =  free_nretrans(p)            + free_nretrans_acc(p,istp)
      
            ! Accumulate total column C fluxes over istp
            npp_active_no3_acc_total(p)  =  npp_active_no3_acc_total(p)  + npp_active_no3_acc(p,istp)
            npp_active_nh4_acc_total(p)  =  npp_active_nh4_acc_total(p)  + npp_active_nh4_acc(p,istp)
            npp_nonmyc_no3_acc_total(p)  =  npp_nonmyc_no3_acc_total(p)  + npp_nonmyc_no3_acc(p,istp)
            npp_nonmyc_nh4_acc_total(p)  =  npp_nonmyc_nh4_acc_total(p)  + npp_nonmyc_nh4_acc(p,istp)
            npp_fix_acc_total(p)         =  npp_fix_acc_total(p)         + npp_fix_acc(p,istp)
            npp_retrans_acc_total(p)     =  npp_retrans_acc_total(p)     + npp_retrans_acc(p,istp) 
                   
         !end if  ! plant_ndemand_pool_step > 0._r8
      end do stp ! NSTEP 

   
      !-------------------------------------------------------------------------------
      ! Turn step level quantities back into fluxes per second. 
      !-------------------------------------------------------------------------------

      !---------------------------N fluxes--------------------!
      Npassive(p)               = n_passive_acc(p)/dt
      Nfix(p)                   = n_fix_acc_total(p)/dt                   
      retransn_to_npool(p)      = n_retrans_acc_total(p)/dt
      free_retransn_to_npool(p) = free_nretrans(p)/dt
      ! this is the N that comes off leaves. 
      Nretrans(p)               = retransn_to_npool(p) + free_retransn_to_npool(p)
      
      
      
      
      !Extract active uptake N from soil pools. 
      do j = 1, nlevdecomp
         !RF change. The N fixed doesn't actually come out of the soil mineral pools, it is 'new'... 
         sminn_to_plant_fun_no3_vr(p,j)    = (n_passive_no3_vr(p,j)  + n_active_no3_vr(p,j) &
                                             + n_nonmyc_no3_vr(p,j))/(dzsoi_decomp(j)*dt)
         sminn_to_plant_fun_nh4_vr(p,j)    = (n_passive_nh4_vr(p,j)  + n_active_nh4_vr(p,j) &
                                             + n_nonmyc_nh4_vr(p,j))/(dzsoi_decomp(j)*dt)
         
      end do
      
     
      
      Nactive_no3(p)            = n_active_no3_acc_total(p)/dt   + n_active_no3_retrans_total(p)/dt
      Nactive_nh4(p)            = n_active_nh4_acc_total(p)/dt   + n_active_nh4_retrans_total(p)/dt  
      
      
    
      Necm_no3(p)               = n_ecm_no3_acc(p)/dt            + n_ecm_no3_retrans(p)/dt
      Necm_nh4(p)               = n_ecm_nh4_acc(p)/dt            + n_ecm_nh4_retrans(p)/dt      
      Necm(p)                   = Necm_no3(p) + Necm_nh4(p)
      Nam_no3(p)                = n_am_no3_acc(p)/dt             + n_am_no3_retrans(p)/dt
      Nam_nh4(p)                = n_am_nh4_acc(p)/dt             + n_am_nh4_retrans(p)/dt
      Nam(p)                    = Nam_no3(p) + Nam_nh4(p)
      Nnonmyc_no3(p)            = n_nonmyc_no3_acc_total(p)/dt   + n_nonmyc_no3_retrans_total(p)/dt
      Nnonmyc_nh4(p)            = n_nonmyc_nh4_acc_total(p)/dt   + n_nonmyc_nh4_retrans_total(p)/dt
      Nnonmyc(p)                = Nnonmyc_no3(p) + Nnonmyc_nh4(p)
      plant_ndemand_retrans(p)  = plant_ndemand_retrans(p)/dt
      Nuptake(p)                = Nactive_no3(p) + Nactive_nh4(p) + Nnonmyc_no3(p) &
                                  + Nnonmyc_nh4(p) + Nfix(p) + Npassive(p) + &
                                  retransn_to_npool(p)+free_retransn_to_npool(p) 
      Nactive(p)                = Nactive_no3(p)  + Nactive_nh4(p) + Nnonmyc_no3(p) + Nnonmyc_nh4(p)
                                   
     ! free N goes straight to the npool, not throught Nuptake...
      sminn_to_plant_fun(p)     = Nactive_no3(p) + Nactive_nh4(p) + Nnonmyc_no3(p) + Nnonmyc_nh4(p) + Nfix(p) + Npassive(p)
 
 
      soil_n_extraction = ( sum(n_active_no3_vr(p,1: nlevdecomp))+sum(n_nonmyc_no3_vr(p,1: nlevdecomp))+&
      sum(n_active_nh4_vr(p,1: nlevdecomp)) + sum(n_nonmyc_nh4_vr(p,1: nlevdecomp)))
      
      !---------------------------C fluxes--------------------!

      npp_Nactive_no3(p)        = npp_active_no3_acc_total(p)/dt + npp_active_no3_retrans_total(p)/dt
      npp_Nactive_nh4(p)        = npp_active_nh4_acc_total(p)/dt + npp_active_nh4_retrans_total(p)/dt
      
      npp_Nnonmyc_no3(p)        = npp_nonmyc_no3_acc_total(p)/dt + npp_nonmyc_no3_retrans_total(p)/dt
      npp_Nnonmyc_nh4(p)        = npp_nonmyc_nh4_acc_total(p)/dt + npp_nonmyc_nh4_retrans_total(p)/dt
      npp_Nactive(p)            = npp_Nactive_no3(p) + npp_Nactive_nh4(p) + npp_Nnonmyc_no3(p) + npp_Nnonmyc_nh4(p)
      npp_Nnonmyc(p)            = npp_Nnonmyc_no3(p) + npp_Nnonmyc_nh4(p)              
      npp_Nfix(p)               = npp_fix_acc_total(p)/dt      
      npp_Nretrans(p)           = npp_retrans_acc_total(p)/dt  
     
      !---------------------------Extra Respiration Fluxes--------------------!      
      soilc_change(p)           = (npp_active_no3_acc_total(p)      + npp_active_nh4_acc_total(p) &
                                    + npp_nonmyc_no3_acc_total(p)     &  
                                    + npp_nonmyc_nh4_acc_total(p)     + npp_fix_acc_total(p))/dt      &
                                    + npp_Nretrans(p)
      soilc_change(p)           = soilc_change(p) + burned_off_carbon / dt                 
      burnedoff_carbon(p)       = burned_off_carbon/dt          
      npp_Nuptake(p)            = soilc_change(p)
      ! how much carbon goes to growth of tissues?  
      npp_growth(p)             = (Nuptake(p)- free_retransn_to_npool(p))*plantCN(p)+(excess_carbon_acc/dt) !does not include gresp, since this is calculated from growth


     
      !-----------------------Diagnostic Fluxes------------------------------!
      if(availc(p).gt.0.0_r8)then !what happens in the night? 
        nuptake_npp_fraction_patch(p) = npp_Nuptake(p)/availc(p)
      else
        nuptake_npp_fraction_patch(p) = spval
      endif 
      if(npp_Nfix(p).gt.0.0_r8)then
        cost_nfix(p) = Nfix(p)/npp_Nfix(p)
      else
        cost_nfix(p) = spval
      endif 
      if(npp_Nactive(p).gt.0.0_r8)then
        cost_nactive(p) = Nactive(p)/npp_Nactive(p)
      else
        cost_nactive(p) = spval
      endif 
      if(npp_Nretrans(p).gt.0.0_r8)then
        cost_nretrans(p) = Nretrans(p)/npp_Nretrans(p)
      else
        cost_nretrans(p) = spval
      endif 
       
       
  end do pft ! PFT Ends 

  call t_stopf('CNFUNcalcs')
  
  call p2c(bounds, num_soilc, filter_soilc,                               &
           cnveg_carbonflux_inst%soilc_change_patch(bounds%begp:bounds%endp), &
           soilbiogeochem_carbonflux_inst%soilc_change_col(bounds%begc:bounds%endc))
           
  call p2c(bounds, num_soilc, filter_soilc,                               &
           cnveg_nitrogenflux_inst%Nfix_patch(bounds%begp:bounds%endp), &
           soilbiogeochem_nitrogenflux_inst%nfix_to_sminn_col(bounds%begc:bounds%endc))
    
  end associate
  end subroutine CNFUN
!=========================================================================================
  real(r8) function fun_cost_fix(fixer,a_fix,b_fix,c_fix,big_cost,crootfr,s_fix, tc_soisno)

! Description:
!   Calculate the cost of fixing N by nodules.
! Code Description:
!   This code is written to CLM4CN by Mingjie Shi on 06/27/2013

  implicit none
!--------------------------------------------------------------------------
! Function result.
!--------------------------------------------------------------------------
! real(r8) , intent(out) :: cost_of_n   !!! cost of fixing N (kgC/kgN)
!--------------------------------------------------------------------------
! Scalar arguments with intent(in).
!--------------------------------------------------------------------------
  integer,  intent(in) :: fixer     ! flag indicating if plant is a fixer
                                    ! 1=yes, otherwise no.
  real(r8), intent(in) :: a_fix     ! As in Houlton et al. (Nature) 2008
  real(r8), intent(in) :: b_fix     ! As in Houlton et al. (Nature) 2008
  real(r8), intent(in) :: c_fix     ! As in Houlton et al. (Nature) 2008
  real(r8), intent(in) :: big_cost  ! an arbitrary large cost (gC/gN)
  real(r8), intent(in) :: crootfr   ! fraction of roots for carbon that are in this layer
  real(r8), intent(in) :: s_fix     ! Inverts Houlton et al. 2008 and constrains between 7.5 and 12.5
  real(r8), intent(in) :: tc_soisno ! soil temperature (degrees Celsius)

  if (fixer == 1 .and. crootfr > 1.e-6_r8) then
     fun_cost_fix  = s_fix * (exp(a_fix + b_fix * tc_soisno * (1._r8 - 0.5_r8 * tc_soisno / c_fix)) - 2._r8)
     
     
     ! New term to directly account for Ben Houlton's temperature response function. 
     ! Assumes s_fix is -6.  (RF, Jan 2015)  
     ! 1.25 converts from the Houlton temp response function to a 0-1 limitation factor. 
     ! The cost of N should probably be 6 gC/gN (or 9, including maintenance costs of nodules) 
     ! for 'optimal' temperatures. This cost should increase in a way that mirrors 
     ! Houlton et al's observations of temperautre limitations on the mirboial fixation rates. 
     ! We don't actually simulate the rate of fixation (and assume that N uptake is instantaneous) 
     ! here, so instead the limitation term is here rolled into the cost function.  
     
     ! Here we invert the 'cost' to give the optimal N:C ratio (1/6 gN/gC)  The amount of N 
     ! you get for a given C goes down as it gets colder, so this can be multiplied by 
     ! the temperature function to give a temperature-limited N:C of  f/6. This number 
     ! can then be inverted to give a temperature limited C:N, as 1/(f/6). Which is the 
     ! same as 6/f, given here" 
     fun_cost_fix  = (-1*s_fix) * 1.0_r8 / (1.25_r8* (exp(a_fix + b_fix * tc_soisno * (1._r8 - 0.5_r8 * tc_soisno / c_fix)) ))
  else
     fun_cost_fix = big_cost
  end if    ! ends up with the fixer or non-fixer decision
  
  end function fun_cost_fix
!=========================================================================================
  real(r8) function fun_cost_active(sminn_layer,big_cost,kc_active,kn_active,rootc_dens,crootfr,smallValue)         

! Description:
!    Calculate the cost of active uptake of N frm the soil.
! Code Description:
!   This code is written to CLM4 by Mingjie Shi.

  implicit none
!--------------------------------------------------------------------------
! Function result.
!--------------------------------------------------------------------------
  real(r8), intent(in) :: sminn_layer   !  Amount of N (as NH4 or NO3) in the soil that is available to plants (gN/m2).
  real(r8), intent(in) :: big_cost      !  An arbitrary large cost (gC/gN).
  real(r8), intent(in) :: kc_active     !  Constant for cost of active uptake (gC/m2).
  real(r8), intent(in) :: kn_active     !  Constant for cost of active uptake (gC/m2).
  real(r8), intent(in) :: rootc_dens    !  Root carbon density in layer (gC/m3).
  real(r8), intent(in) :: crootfr        !  Fraction of roots that are in this layer.
  real(r8), intent(in) :: smallValue    !  A small number.

  if (rootc_dens > 1.e-6_r8.and.sminn_layer > smallValue) then
     fun_cost_active =  kn_active/sminn_layer + kc_active/rootc_dens 
  else
!    There are very few roots in this layer. Set a high cost.
     fun_cost_active =  big_cost
  end if
 
  end function fun_cost_active
!=========================================================================================
  real(r8) function fun_cost_nonmyc(sminn_layer,big_cost,kc_nonmyc,kn_nonmyc,rootc_dens,crootfr,smallValue)         

! Description:
!    Calculate the cost of nonmyc uptake of N frm the soil.
! Code Description:
!   This code is written to CLM4 by Mingjie Shi.

  implicit none
!--------------------------------------------------------------------------
! Function result.
!--------------------------------------------------------------------------
  real(r8), intent(in) :: sminn_layer   !  Amount of N (as NH4 or NO3) in the soil that is available to plants (gN/m2).
  real(r8), intent(in) :: big_cost      !  An arbitrary large cost (gC/gN).
  real(r8), intent(in) :: kc_nonmyc     !  Constant for cost of nonmyc uptake (gC/m2).
  real(r8), intent(in) :: kn_nonmyc     !  Constant for cost of nonmyc uptake (gC/m2).
  real(r8), intent(in) :: rootc_dens   !  Root carbon density in layer (gC/m3).
  real(r8), intent(in) :: crootfr        !  Fraction of roots that are in this layer.
  real(r8), intent(in) :: smallValue    !  A small number.

  if (rootc_dens > 1.e-6_r8.and.sminn_layer > smallValue) then
    fun_cost_nonmyc =  kn_nonmyc / sminn_layer + kc_nonmyc / rootc_dens 
  else
!   There are very few roots in this layer. Set a high cost.
    fun_cost_nonmyc = big_cost
  end if

  end function fun_cost_nonmyc

!==========================================================================

 subroutine fun_retranslocation(p,dt,npp_to_spend,total_falling_leaf_c,         &
               total_falling_leaf_n, total_n_resistance, total_c_spent_retrans, &
               total_c_accounted_retrans, free_n_retrans, paid_for_n_retrans,   &
               target_leafcn, grperc, plantCN)
!
! Description:
! This subroutine (should it be a function?) calculates the amount of N absorbed and C spent 
! during retranslocation. 
! Rosie Fisher. April 2016. 
! !USES:
  implicit none 

! !ARGUMENTS:
  real(r8), intent(IN) :: total_falling_leaf_c  ! INPUT  gC/m2/timestep
  real(r8), intent(IN) :: total_falling_leaf_n  ! INPUT  gC/m2/timestep
  real(r8), intent(IN) :: total_n_resistance    ! INPUT  gC/gN
  real(r8), intent(IN) :: npp_to_spend          ! INPUT  gN/m2/timestep
  real(r8), intent(IN) :: target_leafcn         ! INPUT  gC/gN
  real(r8), intent(IN) :: dt                    ! INPUT  seconds
  real(r8), intent(IN) :: grperc                ! INPUT growth respiration fraction
  real(r8), intent(IN) :: plantCN               ! INPUT plant CN ratio
  integer, intent(IN)  :: p                     ! INPUT  patch index

  real(r8), intent(OUT) :: total_c_spent_retrans     ! OUTPUT gC/m2/timestep
  real(r8), intent(OUT) :: total_c_accounted_retrans ! OUTPUT gC/m2/timestep
  real(r8), intent(OUT) :: paid_for_n_retrans        ! OUTPUT gN/m2/timestep
  real(r8), intent(OUT) :: free_n_retrans            ! OUTPUT gN/m2/timestep

  !
  ! !LOCAL VARIABLES:
  real(r8) :: kresorb               ! INTERNAL used factor
  real(r8) :: falling_leaf_c        ! INTERNAL gC/m2/timestep
  real(r8) :: falling_leaf_n        ! INTERNAL gN/m2/timestep
  real(r8) :: falling_leaf_cn       ! INTERNAL gC/gN
  real(r8) :: cost_retrans_temp     ! INTERNAL gC/gN
  real(r8) :: leaf_n_ext            ! INTERNAL gN/m2/timestep
  real(r8) :: c_spent_retrans       ! INTERNAL gC/m2/timestep
  real(r8) :: c_accounted_retrans   ! INTERNAL gC/m2/timestep
  real(r8) :: npp_to_spend_temp     ! INTERNAL gC/m2/timestep
  real(r8) :: max_falling_leaf_cn   ! INTERNAL gC/gN
  real(r8) :: min_falling_leaf_cn   ! INTERNAL gC/gN
  real(r8) :: cost_escalation       ! INTERNAL cost function parameter
  integer  :: iter                  ! INTERNAL
  integer  :: exitloop              ! INTERNAL
  ! ------------------------------------------------------------------------------- 


   ! ------------------ Initialize total fluxes. ------------------!
   total_c_spent_retrans = 0.0_r8
   total_c_accounted_retrans = 0.0_r8
   c_accounted_retrans   = 0.0_r8
   paid_for_n_retrans    = 0.0_r8
   npp_to_spend_temp     = npp_to_spend

   ! ------------------ Initial C and N pools in falling leaves. ------------------!
   falling_leaf_c       =  total_falling_leaf_c      
   falling_leaf_n       =  total_falling_leaf_n 

   !  ------------------ PARAMETERS ------------------ 
   max_falling_leaf_cn = target_leafcn * 3.0_r8 
   min_falling_leaf_cn = target_leafcn * 1.5_r8
   cost_escalation     = 1.3_r8

   !  ------------------ Free uptake ------------------ 
   free_n_retrans  = max(falling_leaf_n -  (falling_leaf_c/min_falling_leaf_cn),0.0_r8)
   falling_leaf_n = falling_leaf_n -  free_n_retrans 

   ! ------------------ Initial CN ratio and costs ------------------!  
   falling_leaf_cn      = falling_leaf_c/falling_leaf_n 
   kresorb =  (1.0_r8/target_leafcn)
   cost_retrans_temp    = kresorb / ((1.0_r8/falling_leaf_cn )**1.3_r8)

   ! ------------------ Iteration loops to figure out extraction limit ------------!
   iter = 0
   exitloop = 0
   do while(exitloop==0.and.cost_retrans_temp .lt. total_n_resistance.and. &
            falling_leaf_n.ge.0.0_r8.and.npp_to_spend.gt.0.0_r8)
      ! ------------------ Spend some C on removing N ------------!
      ! spend enough C to increase leaf C/N by 1 unit. 
      c_spent_retrans   = cost_retrans_temp * (falling_leaf_n - falling_leaf_c / &
                          (falling_leaf_cn + 1.0_r8))
      ! don't spend more C than you have  
      c_spent_retrans   = min(npp_to_spend_temp, c_spent_retrans) 
      ! N extracted, per this amount of C expenditure
      leaf_n_ext        = c_spent_retrans / cost_retrans_temp     
      ! Do not empty N pool 
      leaf_n_ext        = min(falling_leaf_n, leaf_n_ext)    
      !How much C do you need to account for the N that got taken up? 
      c_accounted_retrans = leaf_n_ext * plantCN * (1.0_r8 + grperc)      

      ! ------------------ Update leafCN, recalculate costs ------------!
      falling_leaf_n    = falling_leaf_n - leaf_n_ext          ! remove N from falling leaves pool 
      if(falling_leaf_n.gt.0.0_r8)then
         falling_leaf_cn   = falling_leaf_c/falling_leaf_n     ! C/N ratio
         cost_retrans_temp = kresorb /((1.0_r8/falling_leaf_cn)**1.3_r8) ! cost function. PARAMETER
      else
         exitloop=1
      endif 
 
      ! ------------------ Accumulate total fluxes ------------!
      total_c_spent_retrans     = total_c_spent_retrans + c_spent_retrans 
      total_c_accounted_retrans = total_c_accounted_retrans + c_accounted_retrans 
      paid_for_n_retrans    = paid_for_n_retrans    + leaf_n_ext
      npp_to_spend_temp     = npp_to_spend_temp     - c_spent_retrans  - c_accounted_retrans
      iter = iter+1
   
      ! run out of C or N
      if(npp_to_spend_temp.le.0.0_r8)then
         exitloop=1
         ! if we made a solving error on this (expenditure and n uptake should 
         ! really be solved simultaneously)
         ! then remove the error from the expenditure. This changes the notional cost, 
         ! but only by a bit and prevents cpool errors. 

         total_c_spent_retrans  = total_c_spent_retrans + npp_to_spend_temp 
      endif 
      ! leaf CN is too high
      if(falling_leaf_cn.ge.max_falling_leaf_cn)then
         exitloop=1
      endif
      ! safety check to prevent hanging code
      if(iter.ge.150)then
          exitloop=1
      endif 
   end do

 end subroutine fun_retranslocation

!==========================================================================

end module CNFUNMod 
