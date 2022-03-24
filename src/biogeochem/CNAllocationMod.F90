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
  use clm_varctl           , only : use_c13, use_c14
  use PatchType            , only : patch
  use pftconMod            , only : pftcon, npcropmin
  use CropType             , only : crop_type
  use PhotosynthesisMod    , only : photosyns_type
  use CanopyStateType      , only : canopystate_type
  use CNVegCarbonStateType , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType  , only : cnveg_carbonflux_type
  use CropReprPoolsMod     , only : nrepr
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams           ! Read in parameters from file
  public :: calc_gpp_mr_availc

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

end module CNAllocationMod
