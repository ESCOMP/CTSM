module CNAllocationMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains subroutines to calculate allocation of C and N to different
  ! plant components. It also contains subroutines to calculate gpp and maintenance
  ! respiration.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use abortutils           , only : endrun
  use decompMod            , only : bounds_type
  use clm_varcon           , only : secspday
  use clm_varctl           , only : use_c13, use_c14, iulog
  use PatchType            , only : patch
  use pftconMod            , only : pftcon, npcropmin
  use CropType             , only : crop_type
  use CropType             , only : cphase_planted, cphase_leafemerge, cphase_grainfill
  use PhotosynthesisMod    , only : photosyns_type
  use CanopyStateType      , only : canopystate_type
  use CNVegCarbonStateType , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType  , only : cnveg_carbonflux_type
  use CNVegStateType       , only : cnveg_state_type
  use CropReprPoolsMod     , only : nrepr
  use CNPhenologyMod       , only : CropPhase
  use CNSharedParamsMod    , only : use_fun
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams                     ! Read in parameters from file
  public :: calc_gpp_mr_availc             ! Calculate total GPP, various maintenance respiration terms, and total available C for allocation
  public :: calc_crop_allocation_fractions ! Calculate crop allocation fractions to leaf, stem, root and repr
  public :: calc_allometry                 ! Calculate c_allometry and n_allometry terms based on allocation fractions

  ! !PRIVATE MEMBER VARIABLES:
  type, private :: params_type
     real(r8) :: dayscrecover          ! number of days to recover negative cpool
  end type params_type

  type(params_type), private :: params_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine readParams (ncid)
    !
    ! !USES:
    use ncdio_pio   , only : file_desc_t,ncd_io
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNAllocParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in parameter
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    ! read in parameters

    tString='dayscrecover'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%dayscrecover=tempr

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine calc_gpp_mr_availc(bounds, num_soilp, filter_soilp, &
       crop_inst, photosyns_inst, canopystate_inst, &
       cnveg_carbonstate_inst, cnveg_carbonflux_inst, &
       c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst)
    !
    ! !DESCRIPTION:
    ! Calculate total GPP, various maintenance respiration terms, and total available C
    ! for allocation
    !
    ! !USES:
    use CNSharedParamsMod           , only : use_matrixcn
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilp        ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:)  ! filter for soil patches
    type(crop_type)                 , intent(in)    :: crop_inst
    type(photosyns_type)            , intent(in)    :: photosyns_inst
    type(canopystate_type)          , intent(in)    :: canopystate_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: c13_cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: c14_cnveg_carbonflux_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: p, k                ! indices
    integer  :: fp                  ! filter patch index
    real(r8) :: dayscrecover        ! number of days to recover negative cpool
    real(r8) :: mr                  ! maintenance respiration (gC/m2/s)
    real(r8) :: reproductive_mr_tot ! total maintenance respiration from grain components (gC/m2/s)
    real(r8) :: curmr, curmr_ratio  ! xsmrpool temporary variables

    character(len=*), parameter :: subname = 'calc_gpp_mr_availc'
    !-----------------------------------------------------------------------

    ! set number of days to recover negative cpool
    dayscrecover = params_inst%dayscrecover

    associate ( &
         ivt                   => patch%itype                                        ,  & ! Input:  [integer  (:) ]  patch vegetation type
         woody                 => pftcon%woody                                     ,  & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
         croplive              => crop_inst%croplive_patch                          , & ! Input:  [logical  (:)   ]  flag, true if planted, not harvested
         psnsun                => photosyns_inst%psnsun_patch                       , & ! Input:  [real(r8) (:)   ]  sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
         psnsha                => photosyns_inst%psnsha_patch                       , & ! Input:  [real(r8) (:)   ]  shaded leaf-level photosynthesis (umol CO2 /m**2/ s)
         c13_psnsun            => photosyns_inst%c13_psnsun_patch                   , & ! Input:  [real(r8) (:)   ]  sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
         c13_psnsha            => photosyns_inst%c13_psnsha_patch                   , & ! Input:  [real(r8) (:)   ]  shaded leaf-level photosynthesis (umol CO2 /m**2/ s)
         c14_psnsun            => photosyns_inst%c14_psnsun_patch                   , & ! Input:  [real(r8) (:)   ]  sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
         c14_psnsha            => photosyns_inst%c14_psnsha_patch                   , & ! Input:  [real(r8) (:)   ]  shaded leaf-level photosynthesis (umol CO2 /m**2/ s)
         laisun                => canopystate_inst%laisun_patch                     , & ! Input:  [real(r8) (:)   ]  sunlit projected leaf area index
         laisha                => canopystate_inst%laisha_patch                     , & ! Input:  [real(r8) (:)   ]  shaded projected leaf area index
         xsmrpool              => cnveg_carbonstate_inst%xsmrpool_patch             , & ! Input:  [real(r8) (:)   ]  (gC/m2) temporary photosynthate C pool
         leaf_mr               => cnveg_carbonflux_inst%leaf_mr_patch               , & ! Input:  [real(r8) (:)   ]
         froot_mr              => cnveg_carbonflux_inst%froot_mr_patch              , & ! Input:  [real(r8) (:)   ]
         livestem_mr           => cnveg_carbonflux_inst%livestem_mr_patch           , & ! Input:  [real(r8) (:)   ]
         livecroot_mr          => cnveg_carbonflux_inst%livecroot_mr_patch          , & ! Input:  [real(r8) (:)   ]
         reproductive_mr       => cnveg_carbonflux_inst%reproductive_mr_patch       , & ! Input:  [real(r8) (:,:) ]
         psnsun_to_cpool       => cnveg_carbonflux_inst%psnsun_to_cpool_patch       , & ! Output: [real(r8) (:)   ]
         psnshade_to_cpool     => cnveg_carbonflux_inst%psnshade_to_cpool_patch     , & ! Output: [real(r8) (:)   ]
         gpp                   => cnveg_carbonflux_inst%gpp_before_downreg_patch    , & ! Output: [real(r8) (:)   ]  GPP flux before downregulation (gC/m2/s)
         availc                => cnveg_carbonflux_inst%availc_patch                , & ! Output: [real(r8) (:)   ]  C flux available for allocation (gC/m2/s)
         leaf_curmr            => cnveg_carbonflux_inst%leaf_curmr_patch            , & ! Output: [real(r8) (:)   ]
         froot_curmr           => cnveg_carbonflux_inst%froot_curmr_patch           , & ! Output: [real(r8) (:)   ]
         livestem_curmr        => cnveg_carbonflux_inst%livestem_curmr_patch        , & ! Output: [real(r8) (:)   ]
         livecroot_curmr       => cnveg_carbonflux_inst%livecroot_curmr_patch       , & ! Output: [real(r8) (:)   ]
         reproductive_curmr    => cnveg_carbonflux_inst%reproductive_curmr_patch    , & ! Output: [real(r8) (:,:) ]
         leaf_xsmr             => cnveg_carbonflux_inst%leaf_xsmr_patch             , & ! Output: [real(r8) (:)   ]
         froot_xsmr            => cnveg_carbonflux_inst%froot_xsmr_patch            , & ! Output: [real(r8) (:)   ]
         livestem_xsmr         => cnveg_carbonflux_inst%livestem_xsmr_patch         , & ! Output: [real(r8) (:)   ]
         livecroot_xsmr        => cnveg_carbonflux_inst%livecroot_xsmr_patch        , & ! Output: [real(r8) (:)   ]
         reproductive_xsmr     => cnveg_carbonflux_inst%reproductive_xsmr_patch     , & ! Output: [real(r8) (:,:) ]
         cpool_to_xsmrpool     => cnveg_carbonflux_inst%cpool_to_xsmrpool_patch     , & ! Output: [real(r8) (:)   ]
         xsmrpool_recover      => cnveg_carbonflux_inst%xsmrpool_recover_patch        & ! Output: [real(r8) (:)   ]  C flux assigned to recovery of negative cpool (gC/m2/s)
         )

    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! get the time step total gross photosynthesis
       ! this is coming from the canopy fluxes code, and is the
       ! gpp that is used to control stomatal conductance.
       ! For the nitrogen downregulation code, this is assumed
       ! to be the potential gpp, and the actual gpp will be
       ! reduced due to N limitation.

       ! Convert psn from umol/m2/s -> gC/m2/s

       ! The input psn (psnsun and psnsha) are expressed per unit LAI
       ! in the sunlit and shaded canopy, respectively. These need to be
       ! scaled by laisun and laisha to get the total gpp for allocation

       ! Note that no associate statement is used for the isotope carbon fluxes below
       ! since they are not always allocated AND nag compiler will complain if you try to
       ! to have an associate statement with unallocated memory

       psnsun_to_cpool(p)   = psnsun(p) * laisun(p) * 12.011e-6_r8
       psnshade_to_cpool(p) = psnsha(p) * laisha(p) * 12.011e-6_r8

       if ( use_c13 ) then
          c13_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)   = c13_psnsun(p) * laisun(p) * 12.011e-6_r8
          c13_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p) = c13_psnsha(p) * laisha(p) * 12.011e-6_r8
       end if

       if ( use_c14 ) then
          c14_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)   = c14_psnsun(p) * laisun(p) * 12.011e-6_r8
          c14_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p) = c14_psnsha(p) * laisha(p) * 12.011e-6_r8
       end if

       gpp(p) = psnsun_to_cpool(p) + psnshade_to_cpool(p)

       ! get the time step total maintenance respiration
       ! These fluxes should already be in gC/m2/s

       mr = leaf_mr(p) + froot_mr(p)
       if (woody(ivt(p)) == 1.0_r8) then
          mr = mr + livestem_mr(p) + livecroot_mr(p)
       else if (ivt(p) >= npcropmin) then
          if (croplive(p)) then
             reproductive_mr_tot = 0._r8
             do k = 1, nrepr
                reproductive_mr_tot = reproductive_mr_tot + reproductive_mr(p,k)
             end do
             mr = mr + livestem_mr(p) + reproductive_mr_tot
          end if
       end if
       ! For Matrix solution if mr is very small set it to zero
       if(mr < -1.e-15_r8 .and. use_matrixcn) mr = 0.0_r8

       ! carbon flux available for allocation
       availc(p) = gpp(p) - mr

       ! new code added for isotope calculations, 7/1/05, PET
       ! If mr > gpp, then some mr comes from gpp, the rest comes from
       ! cpool (xsmr)
       if (mr > 0._r8 .and. availc(p) < 0._r8) then
          curmr = gpp(p)
          curmr_ratio = curmr / mr
       else
          curmr_ratio = 1._r8
       end if
       leaf_curmr(p)      = leaf_mr(p) * curmr_ratio
       leaf_xsmr(p)       = leaf_mr(p) - leaf_curmr(p)
       froot_curmr(p)     = froot_mr(p) * curmr_ratio
       froot_xsmr(p)      = froot_mr(p) - froot_curmr(p)
       livestem_curmr(p)  = livestem_mr(p) * curmr_ratio
       livestem_xsmr(p)   = livestem_mr(p) - livestem_curmr(p)
       livecroot_curmr(p) = livecroot_mr(p) * curmr_ratio
       livecroot_xsmr(p)  = livecroot_mr(p) - livecroot_curmr(p)
       do k = 1, nrepr
          reproductive_curmr(p,k) = reproductive_mr(p,k) * curmr_ratio
          reproductive_xsmr(p,k)  = reproductive_mr(p,k) - reproductive_curmr(p,k)
       end do

       ! no allocation when available c is negative
       availc(p) = max(availc(p),0.0_r8)

       ! test for an xsmrpool deficit
       if (xsmrpool(p) < 0.0_r8) then
          ! Running a deficit in the xsmrpool, so the first priority is to let
          ! some availc from this timestep accumulate in xsmrpool.
          ! Determine rate of recovery for xsmrpool deficit

          xsmrpool_recover(p) = -xsmrpool(p)/(dayscrecover*secspday)
          if (xsmrpool_recover(p) < availc(p)) then
             ! available carbon reduced by amount for xsmrpool recovery
             availc(p) = availc(p) - xsmrpool_recover(p)
          else
             ! all of the available carbon goes to xsmrpool recovery
             xsmrpool_recover(p) = availc(p)
             availc(p) = 0.0_r8
          end if
          cpool_to_xsmrpool(p) = xsmrpool_recover(p)
       end if

    end do

    end associate

  end subroutine calc_gpp_mr_availc

  !-----------------------------------------------------------------------
  subroutine calc_crop_allocation_fractions(bounds, num_pcropp, filter_pcropp, &
       crop_inst, cnveg_state_inst)
    !
    ! !DESCRIPTION:
    ! Calculate crop allocation fractions to leaf, stem, root and repr, following
    ! AgroIBIS subroutine phenocrop
    !
    ! This sets the following variables in cnveg_state_inst for all patches in the pcrop
    ! filter:
    ! - aleaf
    ! - astem
    ! - aroot
    ! - arepr
    !
    ! And under some conditions it updates the following variables:
    ! - astemi
    ! - aleafi
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_pcropp       ! number of prog crop patches in filter
    integer                         , intent(in)    :: filter_pcropp(:) ! filter for prognostic crop patches
    type(crop_type)                 , intent(in)    :: crop_inst
    type(cnveg_state_type)          , intent(inout) :: cnveg_state_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p, fp, k
    real(r8) :: fleaf                                      ! fraction allocated to leaf
    real(r8) :: crop_phase(bounds%begp:bounds%endp)

    character(len=*), parameter :: subname = 'calc_crop_allocation_fractions'
    !-----------------------------------------------------------------------

    associate( &
         ivt                   => patch%itype                                       , & ! Input:  [integer  (:) ]  patch vegetation type
         arooti                => pftcon%arooti                                     , & ! Input:  parameter used below
         arootf                => pftcon%arootf                                     , & ! Input:  parameter used below
         bfact                 => pftcon%bfact                                      , & ! Input:  parameter used below
         fleafi                => pftcon%fleafi                                     , & ! Input:  parameter used below
         aleaff                => pftcon%aleaff                                     , & ! Input:  parameter used below
         astemf                => pftcon%astemf                                     , & ! Input:  parameter used below
         allconss              => pftcon%allconss                                   , & ! Input:  parameter used below
         allconsl              => pftcon%allconsl                                   , & ! Input:  parameter used below
         declfact              => pftcon%declfact                                   , & ! Input:  parameter used below
         croplive              => crop_inst%croplive_patch                          , & ! Input:  [logical  (:)   ]  flag, true if planted, not harvested
         hui                   => crop_inst%hui_patch                               , & ! Input:  [real(r8) (:)   ]  crop patch heat unit index (growing degree-days); set to 0 at sowing and accumulated until harvest
         peaklai               => cnveg_state_inst%peaklai_patch                    , & ! Input:  [integer  (:)   ]  1: max allowed lai; 0: not at max
         gddmaturity           => cnveg_state_inst%gddmaturity_patch                , & ! Input:  [real(r8) (:)   ]  gdd needed to harvest
         huigrain              => cnveg_state_inst%huigrain_patch                   , & ! Input:  [real(r8) (:)   ]  same to reach vegetative maturity
         aleafi                => cnveg_state_inst%aleafi_patch                     , & ! Output: [real(r8) (:)   ]  saved allocation coefficient from phase 2
         astemi                => cnveg_state_inst%astemi_patch                     , & ! Output: [real(r8) (:)   ]  saved allocation coefficient from phase 2
         aleaf                 => cnveg_state_inst%aleaf_patch                      , & ! Output: [real(r8) (:)   ]  leaf allocation coefficient
         astem                 => cnveg_state_inst%astem_patch                      , & ! Output: [real(r8) (:)   ]  stem allocation coefficient
         aroot                 => cnveg_state_inst%aroot_patch                      , & ! Output: [real(r8) (:)   ]  root allocation coefficient
         arepr                 => cnveg_state_inst%arepr_patch                        & ! Output: [real(r8) (:,:) ]  reproductive allocation coefficient(s)
         )

    call CropPhase(bounds, num_pcropp, filter_pcropp, crop_inst, cnveg_state_inst, &
         crop_phase = crop_phase(bounds%begp:bounds%endp))

    do fp = 1, num_pcropp
       p = filter_pcropp(fp)

       if (croplive(p)) then
          ! same phases appear in subroutine CropPhenology

          ! Phase 1 completed:
          ! ==================
          ! if hui is less than the number of gdd needed for filling of grain
          ! leaf emergence also has to have taken place for lai changes to occur
          ! and carbon assimilation
          ! Next phase: leaf emergence to start of leaf decline

          if (crop_phase(p) == cphase_leafemerge) then

             ! allocation rules for crops based on maturity and linear decrease
             ! of amount allocated to roots over course of the growing season

             do k = 1, nrepr
                arepr(p,k) = 0._r8
             end do
             if (peaklai(p) == 1) then ! lai at maximum allowed
                aleaf(p) = 1.e-5_r8
                astem(p) = 0._r8
                aroot(p) = 1._r8 - aleaf(p)
             else
                aroot(p) = max(0._r8, min(1._r8, arooti(ivt(p)) -   &
                     (arooti(ivt(p)) - arootf(ivt(p))) *  &
                     min(1._r8, hui(p)/gddmaturity(p))))
                fleaf = fleafi(ivt(p)) * (exp(-bfact(ivt(p))) -         &
                     exp(-bfact(ivt(p))*hui(p)/huigrain(p))) / &
                     (exp(-bfact(ivt(p)))-1) ! fraction alloc to leaf (from J Norman alloc curve)
                aleaf(p) = max(1.e-5_r8, (1._r8 - aroot(p)) * fleaf)
                astem(p) = 1._r8 - aleaf(p) - aroot(p)
             end if

             ! AgroIBIS included here an immediate adjustment to aleaf & astem if the
             ! predicted lai from the above allocation coefficients exceeded laimx.
             ! We have decided to live with lais slightly higher than laimx by
             ! enforcing the cap in the following tstep through the peaklai logic above.

             astemi(p) = astem(p) ! save for use by equations after shift
             aleafi(p) = aleaf(p) ! to reproductive phenology stage begins

             ! Phase 2 completed:
             ! ==================
             ! shift allocation either when enough gdd are accumulated or maximum number
             ! of days has elapsed since planting

          else if (crop_phase(p) == cphase_grainfill) then
             aroot(p) = max(0._r8, min(1._r8, arooti(ivt(p)) - &
                  (arooti(ivt(p)) - arootf(ivt(p))) * min(1._r8, hui(p)/gddmaturity(p))))
             if (astemi(p) > astemf(ivt(p))) then
                astem(p) = max(0._r8, max(astemf(ivt(p)), astem(p) * &
                     (1._r8 - min((hui(p)-                 &
                     huigrain(p))/((gddmaturity(p)*declfact(ivt(p)))- &
                     huigrain(p)),1._r8)**allconss(ivt(p)) )))
             end if

             ! If crops have hit peaklai, then set leaf allocation to small value
             if (peaklai(p) == 1) then
                aleaf(p) = 1.e-5_r8
             else if (aleafi(p) > aleaff(ivt(p))) then
                aleaf(p) = max(1.e-5_r8, max(aleaff(ivt(p)), aleaf(p) * &
                     (1._r8 - min((hui(p)-                    &
                     huigrain(p))/((gddmaturity(p)*declfact(ivt(p)))- &
                     huigrain(p)),1._r8)**allconsl(ivt(p)) )))
             end if

             ! For AgroIBIS-based crop model, all repr allocation is assumed to go
             ! into the last reproductive pool. In practice there is only a single
             ! reproductive pool with the AgroIBIS-based crop model, but for
             ! software testing we can have multiple, in which situation we want the
             ! active pool to be the last one.
             do k = 1, nrepr-1
                arepr(p,k) = 0._r8
             end do
             arepr(p,nrepr) = 1._r8 - aroot(p) - astem(p) - aleaf(p)

          else if (crop_phase(p) == cphase_planted) then
             ! pre emergence
             ! allocation coefficients should be irrelevant because crops have no
             ! live carbon pools
             aleaf(p) = 1._r8
             aleafi(p) = 1._r8
             astem(p) = 0._r8
             astemi(p) = 0._r8
             aroot(p) = 0._r8
             do k = 1, nrepr
                arepr(p,k) = 0._r8
             end do

          else
             write(iulog,*) "ERROR in " // subname // ": unexpected crop_phase: ", crop_phase(p)
             call endrun(msg="ERROR: unexpected crop_phase "//errmsg(sourcefile, __LINE__))
          end if

       else   ! .not croplive
          ! allocation coefficients should be irrelevant because crops have no
          ! live carbon pools
          aleaf(p) = 1._r8
          aleafi(p) = 1._r8
          astem(p) = 0._r8
          astemi(p) = 0._r8
          aroot(p) = 0._r8
          do k = 1, nrepr
             arepr(p,k) = 0._r8
          end do
       end if

    end do

    end associate

  end subroutine calc_crop_allocation_fractions

  !-----------------------------------------------------------------------
  subroutine calc_allometry(num_soilp, filter_soilp, &
       cnveg_carbonflux_inst, cnveg_state_inst)
    !
    ! !DESCRIPTION:
    ! Calculate c_allometry and n_allometry terms based on allocation fractions
    !
    ! !ARGUMENTS:
    integer                         , intent(in)    :: num_soilp        ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:)  ! filter for soil patches
    type(cnveg_carbonflux_type)     , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_state_type)          , intent(inout) :: cnveg_state_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p, fp, k
    real(r8):: f1,f2,f3,f4,g1     ! allocation parameters
    real(r8):: g1a                ! g1 included in allocation/allometry
    real(r8):: cnl,cnfr,cnlw,cndw ! C:N ratios for leaf, fine root, and wood
    real(r8):: f5(nrepr)          ! reproductive allocation parameters
    real(r8):: cng                ! C:N ratio for grain (= cnlw for now; slevis)
    real(r8):: f5_tot             ! sum of f5 terms
    real(r8):: f5_n_tot           ! sum of f5 terms converted from C to N

    character(len=*), parameter :: subname = 'calc_allometry'
    !-----------------------------------------------------------------------

    associate( &
         ivt         => patch%itype                            ,  & ! Input:  [integer  (:) ]  patch vegetation type
         woody       => pftcon%woody                           ,  & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
         froot_leaf  => pftcon%froot_leaf                      ,  & ! Input:  allocation parameter: new fine root C per new leaf C (gC/gC)
         croot_stem  => pftcon%croot_stem                      ,  & ! Input:  allocation parameter: new coarse root C per new stem C (gC/gC)
         stem_leaf   => pftcon%stem_leaf                       ,  & ! Input:  allocation parameter: new stem c per new leaf C (gC/gC)
         flivewd     => pftcon%flivewd                         ,  & ! Input:  allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
         leafcn      => pftcon%leafcn                          ,  & ! Input:  leaf C:N (gC/gN)
         frootcn     => pftcon%frootcn                         ,  & ! Input:  fine root C:N (gC/gN)
         livewdcn    => pftcon%livewdcn                        ,  & ! Input:  live wood (phloem and ray parenchyma) C:N (gC/gN)
         deadwdcn    => pftcon%deadwdcn                        ,  & ! Input:  dead wood (xylem and heartwood) C:N (gC/gN)
         graincn     => pftcon%graincn                         ,  & ! Input:  grain C:N (gC/gN)
         grperc      => pftcon%grperc                          , & ! Input:  parameter used below
         annsum_npp  => cnveg_carbonflux_inst%annsum_npp_patch ,  & ! Input:  [real(r8) (:)   ]  annual sum of NPP, for wood allocation
         aleaf       => cnveg_state_inst%aleaf_patch           , & ! Input: [real(r8) (:)   ]  leaf allocation coefficient
         astem       => cnveg_state_inst%astem_patch           , & ! Input: [real(r8) (:)   ]  stem allocation coefficient
         aroot       => cnveg_state_inst%aroot_patch           , & ! Input: [real(r8) (:)   ]  root allocation coefficient
         arepr       => cnveg_state_inst%arepr_patch           , & ! Input: [real(r8) (:,:) ]  reproductive allocation coefficient(s)
         c_allometry => cnveg_state_inst%c_allometry_patch     ,  & ! Output: [real(r8) (:)   ]  C allocation index (DIM)
         n_allometry => cnveg_state_inst%n_allometry_patch        & ! Output: [real(r8) (:)   ]  N allocation index (DIM)
         )

    do fp = 1, num_soilp
       p = filter_soilp(fp)

       f1 = froot_leaf(ivt(p))
       f2 = croot_stem(ivt(p))

       ! modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
       ! constrained so that it does not go lower than 0.2 (under negative annsum_npp)
       ! This variable allocation is only for trees. Shrubs have a constant
       ! allocation as specified in the pft-physiology file.  The value is also used
       ! as a trigger here: -1.0 means to use the dynamic allocation (trees).

       if (stem_leaf(ivt(p)) == -1._r8) then
          f3 = (2.7_r8/(1.0_r8+exp(-0.004_r8*(annsum_npp(p) - 300.0_r8)))) - 0.4_r8
       else
          f3 = stem_leaf(ivt(p))
       end if

       f4   = flivewd(ivt(p))
       if (ivt(p) >= npcropmin) then
          g1 = 0.25_r8
       else
          g1 = grperc(ivt(p))
       end if
       cnl  = leafcn(ivt(p))
       cnfr = frootcn(ivt(p))
       cnlw = livewdcn(ivt(p))
       cndw = deadwdcn(ivt(p))

       ! based on available C, use constant allometric relationships to
       ! determine N requirements
       if (.not. use_fun) then
          g1a = g1
       else
          g1a = 0._r8
       end if
       if (woody(ivt(p)) == 1.0_r8) then
          c_allometry(p) = (1._r8+g1a)*(1._r8+f1+f3*(1._r8+f2))
          n_allometry(p) = 1._r8/cnl + f1/cnfr + (f3*f4*(1._r8+f2))/cnlw + &
               (f3*(1._r8-f4)*(1._r8+f2))/cndw
       else if (ivt(p) >= npcropmin) then ! skip generic crops
          cng = graincn(ivt(p))
          f1 = aroot(p) / aleaf(p)
          f3 = astem(p) / aleaf(p)
          do k = 1, nrepr
             f5(k) = arepr(p,k) / aleaf(p)
          end do
          f5_tot = 0._r8
          f5_n_tot = 0._r8
          do k = 1, nrepr
             f5_tot = f5_tot + f5(k)
             ! Note that currently we use the same C/N ratio for all grain components:
             f5_n_tot = f5_n_tot + f5(k)/cng
          end do
          c_allometry(p) = (1._r8+g1a)*(1._r8+f1+f5_tot+f3*(1._r8+f2))
          n_allometry(p) = 1._r8/cnl + f1/cnfr + f5_n_tot + (f3*f4*(1._r8+f2))/cnlw + &
               (f3*(1._r8-f4)*(1._r8+f2))/cndw
       else
          c_allometry(p) = 1._r8+g1a+f1+f1*g1a
          n_allometry(p) = 1._r8/cnl + f1/cnfr
       end if

    end do

    end associate

  end subroutine calc_allometry


end module CNAllocationMod
