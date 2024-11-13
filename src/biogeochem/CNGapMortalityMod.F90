module CNGapMortalityMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines used in gap mortality for coupled carbon
  ! nitrogen code.
  !
  ! !USES:
  use shr_kind_mod                   , only : r8 => shr_kind_r8
  use shr_infnan_mod                 , only : nan => shr_infnan_nan, assignment(=)
  use decompMod                      , only : bounds_type
  use abortutils                     , only : endrun
  use shr_log_mod                    , only : errMsg => shr_log_errMsg
  use clm_varpar                     , only : mxpft
  use pftconMod                      , only : pftcon
  use CNDVType                       , only : dgvs_type
  use CNVegCarbonStateType           , only : cnveg_carbonstate_type, spinup_factor_deadwood
  use CNVegCarbonFluxType            , only : cnveg_carbonflux_type
  use CNVegNitrogenStateType         , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType          , only : cnveg_nitrogenflux_type
  use SoilBiogeochemNitrogenFluxType , only : soilbiogeochem_nitrogenflux_type
  use CanopyStateType                , only : canopystate_type
  use ColumnType                     , only : col
  use PatchType                      , only : patch
  use GridcellType                   , only : grc
  use CNSharedParamsMod              , only : use_matrixcn
  use CNVegMatrixMod                 , only : matrix_update_gmc, matrix_update_gmn
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams
  public :: CNGapMortality

  type, private :: params_type
     real(r8):: k_mort ! coeff. of growth efficiency in mortality equation
     real(r8), allocatable :: r_mort(:)  ! Mortality rate (1/year)
  contains
     procedure, private :: allocParams    ! Allocate the parameters
  end type params_type
  !
  type(params_type), private :: params_inst
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: CNGap_PatchToColumn

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine allocParams ( this )
    !
    implicit none

    ! !ARGUMENTS:
    class(params_type) :: this
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'allocParams'
    !-----------------------------------------------------------------------

    ! allocate parameters

    allocate( this%r_mort    (0:mxpft) )          ; this%r_mort(:)   = nan

  end subroutine allocParams

  !-----------------------------------------------------------------------
  subroutine readParams ( ncid )
    !
    ! !DESCRIPTION:
    ! Read in parameters
    !
    ! !USES:
    use ncdio_pio  , only : file_desc_t,ncd_io
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNGapMortParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    real(r8)           :: temp1d(0:mxpft) ! temporary to read in parameter
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    tString='k_mort'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%k_mort=tempr   

    call params_inst%allocParams()

    tString='r_mort'
    call ncd_io(varname=trim(tString),data=temp1d, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%r_mort=temp1d
    
  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine CNGapMortality (bounds, num_soilp, filter_soilp, &
       dgvs_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst,&
       cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, canopystate_inst, &  
       leaf_prof_patch, froot_prof_patch, croot_prof_patch, stem_prof_patch)  
    !
    ! !DESCRIPTION:
    ! Gap-phase mortality routine for coupled carbon-nitrogen code (CN)
    !
    ! !USES:
    use clm_time_manager , only: get_average_days_per_year, get_step_size
    use clm_varpar       , only: nlevdecomp_full
    use clm_varcon       , only: secspday
    use clm_varctl       , only: use_cndv, spinup_state
    use pftconMod        , only: npcropmin
    !
    ! !ARGUMENTS:
    type(bounds_type)                      , intent(in)    :: bounds
    integer                                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                                , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    type(dgvs_type)                        , intent(inout) :: dgvs_inst
    type(cnveg_carbonstate_type)           , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type)         , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)            , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)          , intent(inout) :: cnveg_nitrogenflux_inst
    type(canopystate_type)                 , intent(in)    :: canopystate_inst            
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    real(r8)                               , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
    real(r8)                               , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
    real(r8)                               , intent(in)    :: croot_prof_patch(bounds%begp:,1:)
    real(r8)                               , intent(in)    :: stem_prof_patch(bounds%begp:,1:)
    !
    ! !LOCAL VARIABLES:
    integer :: p                ! patch index
    integer :: fp               ! patch filter index
    real(r8):: dt               ! time step (sec)
    real(r8):: am               ! rate for fractional mortality (1/yr)
    real(r8):: m                ! rate for fractional mortality (1/s)
    real(r8):: mort_max         ! asymptotic max mortality rate (/yr)
    real(r8):: k_mort = 0.3_r8  ! coeff of growth efficiency in mortality equation
    logical,parameter :: matrixcheck_gm = .False. ! If matrix check should be done
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(leaf_prof_patch)   == (/bounds%endp,nlevdecomp_full/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(froot_prof_patch)  == (/bounds%endp,nlevdecomp_full/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(croot_prof_patch)  == (/bounds%endp,nlevdecomp_full/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(stem_prof_patch)   == (/bounds%endp,nlevdecomp_full/)), sourcefile, __LINE__)

    associate(                                         & 
         ivt                      => patch%itype                  , & ! Input:  [integer  (:) ]  patch vegetation type                                

         woody                    => pftcon%woody               , & ! Input:  binary flag for woody lifeform                    
         
         greffic                  => dgvs_inst%greffic_patch    , & ! Input:  [real(r8) (:) ]                                                    
         heatstress               => dgvs_inst%heatstress_patch , & ! Input:  [real(r8) (:) ]    
         
         leafcn                   => pftcon%leafcn               , & ! Input:  [real(r8) (:)]  leaf C:N (gC/gN)                        
         livewdcn                 => pftcon%livewdcn             , & ! Input:  [real(r8) (:)]  live wood (phloem and ray parenchyma) C:N (gC/gN) 
         laisun                   => canopystate_inst%laisun_patch  , & ! Input:  [real(r8) (:)   ]  sunlit projected leaf area index      
         laisha                   => canopystate_inst%laisha_patch  , & ! Input:  [real(r8) (:)   ]  shaded projected leaf area index   
         nind                     => dgvs_inst%nind_patch                           , & ! Output:[real(r8)(:)] number of individuals (#/m2) added by F. Li and S. Levis
         ileaf_to_iout_gmc        => cnveg_carbonflux_inst%ileaf_to_iout_gm         , & ! Input: [integer (:)] Index of gap mortality related C transfer from leaf pool to outside of vegetation pools
         ileafst_to_iout_gmc      => cnveg_carbonflux_inst%ileafst_to_iout_gm       , & ! Input: [integer (:)] Index of gap mortality related C transfer from leaf storage pool to outside of vegetation pools
         ileafxf_to_iout_gmc      => cnveg_carbonflux_inst%ileafxf_to_iout_gm       , & ! Input: [integer (:)] Index of gap mortality related C transfer from leaf transfer pool to outside of vegetation pools
         ifroot_to_iout_gmc       => cnveg_carbonflux_inst%ifroot_to_iout_gm        , & ! Input: [integer (:)] Index of gap mortality related C transfer from fine root pool to outside of vegetation pools
         ifrootst_to_iout_gmc     => cnveg_carbonflux_inst%ifrootst_to_iout_gm      , & ! Input: [integer (:)] Index of gap mortality related C transfer from fine root storage pool to outside of vegetation pools
         ifrootxf_to_iout_gmc     => cnveg_carbonflux_inst%ifrootxf_to_iout_gm      , & ! Input: [integer (:)] Index of gap mortality related C transfer from fine root transfer pool to outside of vegetation pools
         ilivestem_to_iout_gmc    => cnveg_carbonflux_inst%ilivestem_to_iout_gm     , & ! Input: [integer (:)] Index of gap mortality related C transfer from live stem pool to outside of vegetation pools
         ilivestemst_to_iout_gmc  => cnveg_carbonflux_inst%ilivestemst_to_iout_gm   , & ! Input: [integer (:)] Index of gap mortality related C transfer from live stem storage pool to outside of vegetation pools
         ilivestemxf_to_iout_gmc  => cnveg_carbonflux_inst%ilivestemxf_to_iout_gm   , & ! Input: [integer (:)] Index of gap mortality related C transfer from live stem transfer pool to outside of vegetation pools
         ideadstem_to_iout_gmc    => cnveg_carbonflux_inst%ideadstem_to_iout_gm     , & ! Input: [integer (:)] Index of gap mortality related C transfer from dead stem pool to outside of vegetation pools
         ideadstemst_to_iout_gmc  => cnveg_carbonflux_inst%ideadstemst_to_iout_gm   , & ! Input: [integer (:)] Index of gap mortality related C transfer from dead stem storage pool to outside of vegetation pools
         ideadstemxf_to_iout_gmc  => cnveg_carbonflux_inst%ideadstemxf_to_iout_gm   , & ! Input: [integer (:)] Index of gap mortality related C transfer from dead stem transfer pool to outside of vegetation pools
         ilivecroot_to_iout_gmc   => cnveg_carbonflux_inst%ilivecroot_to_iout_gm    , & ! Input: [integer (:)] Index of gap mortality related C transfer from live coarse root pool to outside of vegetation pools
         ilivecrootst_to_iout_gmc => cnveg_carbonflux_inst%ilivecrootst_to_iout_gm  , & ! Input: [integer (:)] Index of gap mortality related C transfer from live coarse root storage pool to outside of vegetation pools
         ilivecrootxf_to_iout_gmc => cnveg_carbonflux_inst%ilivecrootxf_to_iout_gm  , & ! Input: [integer (:)] Index of gap mortality related C transfer from live coarse root transfer pool to outside of vegetation pools
         ideadcroot_to_iout_gmc   => cnveg_carbonflux_inst%ideadcroot_to_iout_gm    , & ! Input: [integer (:)] Index of gap mortality related C transfer from dead coarse root pool to outside of vegetation pools
         ideadcrootst_to_iout_gmc => cnveg_carbonflux_inst%ideadcrootst_to_iout_gm  , & ! Input: [integer (:)] Index of gap mortality related C transfer from dead coarse root storage pool to outside of vegetation pools
         ideadcrootxf_to_iout_gmc => cnveg_carbonflux_inst%ideadcrootxf_to_iout_gm  , & ! Input: [integer (:)] Index of gap mortality related C transfer from dead coarse root transfer pool to outside of vegetation pools
         ileaf_to_iout_gmn        => cnveg_nitrogenflux_inst%ileaf_to_iout_gm       , & ! Input: [integer (:)] Index of gap mortality related N transfer from leaf pool to outside of vegetation pools
         ileafst_to_iout_gmn      => cnveg_nitrogenflux_inst%ileafst_to_iout_gm     , & ! Input: [integer (:)] Index of gap mortality related N transfer from leaf storage pool to outside of vegetation pools
         ileafxf_to_iout_gmn      => cnveg_nitrogenflux_inst%ileafxf_to_iout_gm     , & ! Input: [integer (:)] Index of gap mortality related N transfer from leaf transfer pool to outside of vegetation pools
         ifroot_to_iout_gmn       => cnveg_nitrogenflux_inst%ifroot_to_iout_gm      , & ! Input: [integer (:)] Index of gap mortality related N transfer from fine root pool to outside of vegetation pools
         ifrootst_to_iout_gmn     => cnveg_nitrogenflux_inst%ifrootst_to_iout_gm    , & ! Input: [integer (:)] Index of gap mortality related N transfer from fine root storage pool to outside of vegetation pools
         ifrootxf_to_iout_gmn     => cnveg_nitrogenflux_inst%ifrootxf_to_iout_gm    , & ! Input: [integer (:)] Index of gap mortality related N transfer from fine root transfer pool to outside of vegetation pools
         ilivestem_to_iout_gmn    => cnveg_nitrogenflux_inst%ilivestem_to_iout_gm   , & ! Input: [integer (:)] Index of gap mortality related N transfer from live stem pool to outside of vegetation pools
         ilivestemst_to_iout_gmn  => cnveg_nitrogenflux_inst%ilivestemst_to_iout_gm , & ! Input: [integer (:)] Index of gap mortality related N transfer from live stem storage pool to outside of vegetation pools
         ilivestemxf_to_iout_gmn  => cnveg_nitrogenflux_inst%ilivestemxf_to_iout_gm , & ! Input: [integer (:)] Index of gap mortality related N transfer from live stem transfer pool to outside of vegetation pools
         ideadstem_to_iout_gmn    => cnveg_nitrogenflux_inst%ideadstem_to_iout_gm   , & ! Input: [integer (:)] Index of gap mortality related N transfer from dead stem pool to outside of vegetation pools
         ideadstemst_to_iout_gmn  => cnveg_nitrogenflux_inst%ideadstemst_to_iout_gm , & ! Input: [integer (:)] Index of gap mortality related N transfer from dead stem storage pool to outside of vegetation pools
         ideadstemxf_to_iout_gmn  => cnveg_nitrogenflux_inst%ideadstemxf_to_iout_gm , & ! Input: [integer (:)] Index of gap mortality related N transfer from dead stem transfer pool to outside of vegetation pools
         ilivecroot_to_iout_gmn   => cnveg_nitrogenflux_inst%ilivecroot_to_iout_gm  , & ! Input: [integer (:)] Index of gap mortality related N transfer from live coarse root pool to outside of vegetation pools
         ilivecrootst_to_iout_gmn => cnveg_nitrogenflux_inst%ilivecrootst_to_iout_gm, & ! Input: [integer (:)] Index of gap mortality related N transfer from live coarse root storage pool to outside of vegetation pools
         ilivecrootxf_to_iout_gmn => cnveg_nitrogenflux_inst%ilivecrootxf_to_iout_gm, & ! Input: [integer (:)] Index of gap mortality related N transfer from live coarse root transfer pool to outside of vegetation pools
         ideadcroot_to_iout_gmn   => cnveg_nitrogenflux_inst%ideadcroot_to_iout_gm  , & ! Input: [integer (:)] Index of gap mortality related N transfer from dead coarse root pool to outside of vegetation pools
         ideadcrootst_to_iout_gmn => cnveg_nitrogenflux_inst%ideadcrootst_to_iout_gm, & ! Input: [integer (:)] Index of gap mortality related N transfer from dead coarse root storage pool to outside of vegetation pools
         ideadcrootxf_to_iout_gmn => cnveg_nitrogenflux_inst%ideadcrootxf_to_iout_gm, & ! Input: [integer (:)] Index of gap mortality related N transfer from dead coarse root transfer pool to outside of vegetation pools
         iretransn_to_iout_gmn    => cnveg_nitrogenflux_inst%iretransn_to_iout_gm     & ! Input: [integer (:)] Index of gap mortality related N transfer from retranslocation pool to outside of vegetation pools
         )

      dt = real( get_step_size(), r8 )
      ! set coeff of growth efficiency in mortality equation 
      k_mort = params_inst%k_mort

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         if (use_cndv) then
            ! Stress mortality from lpj's subr Mortality.

            if (woody(ivt(p)) == 1._r8) then

               if (ivt(p) == 8) then
                  mort_max = 0.03_r8 ! BDT boreal
               else
                  mort_max = 0.01_r8 ! original value for all patches
               end if

               ! heatstress and greffic calculated in Establishment once/yr

               ! Mortality rate inversely related to growth efficiency
               ! (Prentice et al 1993)
               am = mort_max / (1._r8 + k_mort * greffic(p))

               ! Mortality rate inversely related to growth efficiency
               ! (Prentice et al 1993)
               am = mort_max / (1._r8 + k_mort * greffic(p))

               am = min(1._r8, am + heatstress(p))
            else ! lpj didn't set this for grasses; cn does
               ! set the mortality rate based on annual rate
               am = params_inst%r_mort(ivt(p))
            end if

         else
            am = params_inst%r_mort(ivt(p)) 
         end if

         m  = am/(get_average_days_per_year() * secspday)

         !------------------------------------------------------
         ! patch-level gap mortality carbon fluxes
         !------------------------------------------------------

         if(.not. use_matrixcn)then
            ! displayed pools
            cnveg_carbonflux_inst%m_leafc_to_litter_patch(p)               = cnveg_carbonstate_inst%leafc_patch(p)               * m
            cnveg_carbonflux_inst%m_frootc_to_litter_patch(p)              = cnveg_carbonstate_inst%frootc_patch(p)              * m
            cnveg_carbonflux_inst%m_livestemc_to_litter_patch(p)           = cnveg_carbonstate_inst%livestemc_patch(p)           * m
            cnveg_carbonflux_inst%m_livecrootc_to_litter_patch(p)          = cnveg_carbonstate_inst%livecrootc_patch(p)          * m
         else
            cnveg_carbonflux_inst%m_leafc_to_litter_patch(p)               = cnveg_carbonstate_inst%leafc_patch(p)       * matrix_update_gmc(p,ileaf_to_iout_gmc,m,dt,cnveg_carbonflux_inst,matrixcheck_gm,.True.)
            cnveg_carbonflux_inst%m_frootc_to_litter_patch(p)              = cnveg_carbonstate_inst%frootc_patch(p)      * matrix_update_gmc(p,ifroot_to_iout_gmc,m,dt,cnveg_carbonflux_inst,.true.,.True.)
            cnveg_carbonflux_inst%m_livestemc_to_litter_patch(p)           = cnveg_carbonstate_inst%livestemc_patch(p)   * matrix_update_gmc(p,ilivestem_to_iout_gmc,m,dt,cnveg_carbonflux_inst,matrixcheck_gm,.True.)
            cnveg_carbonflux_inst%m_livecrootc_to_litter_patch(p)          = cnveg_carbonstate_inst%livecrootc_patch(p)  * matrix_update_gmc(p,ilivecroot_to_iout_gmc,m,dt,cnveg_carbonflux_inst,matrixcheck_gm,.True.)
         end if
         if(.not. use_matrixcn)then
            cnveg_carbonflux_inst%m_deadstemc_to_litter_patch(p)         = cnveg_carbonstate_inst%deadstemc_patch(p)  * m * spinup_factor_deadwood
            cnveg_carbonflux_inst%m_deadcrootc_to_litter_patch(p)        = cnveg_carbonstate_inst%deadcrootc_patch(p) * m * spinup_factor_deadwood
         else
            cnveg_carbonflux_inst%m_deadstemc_to_litter_patch(p)         = cnveg_carbonstate_inst%deadstemc_patch(p)   * matrix_update_gmc(p,ideadstem_to_iout_gmc, &
                                                                             m*spinup_factor_deadwood,dt,cnveg_carbonflux_inst,matrixcheck_gm,.True.)
            cnveg_carbonflux_inst%m_deadcrootc_to_litter_patch(p)        = cnveg_carbonstate_inst%deadcrootc_patch(p)  * matrix_update_gmc(p,ideadcroot_to_iout_gmc, &
                                                                             m*spinup_factor_deadwood,dt,cnveg_carbonflux_inst,matrixcheck_gm,.True.)
         end if !use_matrixcn

         if(.not. use_matrixcn)then
            ! storage pools
            cnveg_carbonflux_inst%m_leafc_storage_to_litter_patch(p)       = cnveg_carbonstate_inst%leafc_storage_patch(p)       * m
            cnveg_carbonflux_inst%m_frootc_storage_to_litter_patch(p)      = cnveg_carbonstate_inst%frootc_storage_patch(p)      * m
            cnveg_carbonflux_inst%m_livestemc_storage_to_litter_patch(p)   = cnveg_carbonstate_inst%livestemc_storage_patch(p)   * m
            cnveg_carbonflux_inst%m_deadstemc_storage_to_litter_patch(p)   = cnveg_carbonstate_inst%deadstemc_storage_patch(p)   * m
            cnveg_carbonflux_inst%m_livecrootc_storage_to_litter_patch(p)  = cnveg_carbonstate_inst%livecrootc_storage_patch(p)  * m
            cnveg_carbonflux_inst%m_deadcrootc_storage_to_litter_patch(p)  = cnveg_carbonstate_inst%deadcrootc_storage_patch(p)  * m
            cnveg_carbonflux_inst%m_gresp_storage_to_litter_patch(p)       = cnveg_carbonstate_inst%gresp_storage_patch(p)       * m
         
            ! transfer pools
            cnveg_carbonflux_inst%m_leafc_xfer_to_litter_patch(p)          = cnveg_carbonstate_inst%leafc_xfer_patch(p)          * m
            cnveg_carbonflux_inst%m_frootc_xfer_to_litter_patch(p)         = cnveg_carbonstate_inst%frootc_xfer_patch(p)         * m
            cnveg_carbonflux_inst%m_livestemc_xfer_to_litter_patch(p)      = cnveg_carbonstate_inst%livestemc_xfer_patch(p)      * m
            cnveg_carbonflux_inst%m_deadstemc_xfer_to_litter_patch(p)      = cnveg_carbonstate_inst%deadstemc_xfer_patch(p)      * m
            cnveg_carbonflux_inst%m_livecrootc_xfer_to_litter_patch(p)     = cnveg_carbonstate_inst%livecrootc_xfer_patch(p)     * m
            cnveg_carbonflux_inst%m_deadcrootc_xfer_to_litter_patch(p)     = cnveg_carbonstate_inst%deadcrootc_xfer_patch(p)     * m
            cnveg_carbonflux_inst%m_gresp_xfer_to_litter_patch(p)          = cnveg_carbonstate_inst%gresp_xfer_patch(p)          * m
         else
           ! NOTE: The non-matrix version of this is in CNCStateUpdate2Mod CStateUpdate2 (EBK 11/25/2019)

         ! storage pools
             cnveg_carbonflux_inst%m_leafc_storage_to_litter_patch(p)         = cnveg_carbonstate_inst%leafc_storage_patch(p)      * matrix_update_gmc(p,ileafst_to_iout_gmc,m,dt,cnveg_carbonflux_inst,matrixcheck_gm,.True.)
             cnveg_carbonflux_inst%m_frootc_storage_to_litter_patch(p)        = cnveg_carbonstate_inst%frootc_storage_patch(p)     * matrix_update_gmc(p,ifrootst_to_iout_gmc,m,dt,cnveg_carbonflux_inst,matrixcheck_gm,.True.)
             cnveg_carbonflux_inst%m_livestemc_storage_to_litter_patch(p)     = cnveg_carbonstate_inst%livestemc_storage_patch(p)  * matrix_update_gmc(p,ilivestemst_to_iout_gmc,m,dt,cnveg_carbonflux_inst,matrixcheck_gm,.True.)
             cnveg_carbonflux_inst%m_deadstemc_storage_to_litter_patch(p)     = cnveg_carbonstate_inst%deadstemc_storage_patch(p)  * matrix_update_gmc(p,ideadstemst_to_iout_gmc,m,dt,cnveg_carbonflux_inst,matrixcheck_gm,.True.)
             cnveg_carbonflux_inst%m_livecrootc_storage_to_litter_patch(p)    = cnveg_carbonstate_inst%livecrootc_storage_patch(p) * matrix_update_gmc(p,ilivecrootst_to_iout_gmc,m,dt,cnveg_carbonflux_inst,matrixcheck_gm,.True.)
             cnveg_carbonflux_inst%m_deadcrootc_storage_to_litter_patch(p)    = cnveg_carbonstate_inst%deadcrootc_storage_patch(p) * matrix_update_gmc(p,ideadcrootst_to_iout_gmc,m,dt,cnveg_carbonflux_inst,matrixcheck_gm,.True.)
            
         ! transfer pools
             cnveg_carbonflux_inst%m_leafc_xfer_to_litter_patch(p)         = cnveg_carbonstate_inst%leafc_xfer_patch(p)         * matrix_update_gmc(p,ileafxf_to_iout_gmc,m,dt,cnveg_carbonflux_inst,matrixcheck_gm,.True.)
             cnveg_carbonflux_inst%m_frootc_xfer_to_litter_patch(p)        = cnveg_carbonstate_inst%frootc_xfer_patch(p)        * matrix_update_gmc(p,ifrootxf_to_iout_gmc,m,dt,cnveg_carbonflux_inst,matrixcheck_gm,.True.)
             cnveg_carbonflux_inst%m_livestemc_xfer_to_litter_patch(p)     = cnveg_carbonstate_inst%livestemc_xfer_patch(p)     * matrix_update_gmc(p,ilivestemxf_to_iout_gmc,m,dt,cnveg_carbonflux_inst,matrixcheck_gm,.True.)
             cnveg_carbonflux_inst%m_deadstemc_xfer_to_litter_patch(p)     = cnveg_carbonstate_inst%deadstemc_xfer_patch(p)     * matrix_update_gmc(p,ideadstemxf_to_iout_gmc,m,dt,cnveg_carbonflux_inst,matrixcheck_gm,.True.)
             cnveg_carbonflux_inst%m_livecrootc_xfer_to_litter_patch(p)    = cnveg_carbonstate_inst%livecrootc_xfer_patch(p)    * matrix_update_gmc(p,ilivecrootxf_to_iout_gmc,m,dt,cnveg_carbonflux_inst,matrixcheck_gm,.True.)
             cnveg_carbonflux_inst%m_deadcrootc_xfer_to_litter_patch(p)    = cnveg_carbonstate_inst%deadcrootc_xfer_patch(p)    * matrix_update_gmc(p,ideadcrootxf_to_iout_gmc,m,dt,cnveg_carbonflux_inst,matrixcheck_gm,.True.)
         end if !use_matrixcn

         !------------------------------------------------------
         ! patch-level gap mortality nitrogen fluxes
         !------------------------------------------------------

         ! displayed pools
         if(.not. use_matrixcn)then
            cnveg_nitrogenflux_inst%m_leafn_to_litter_patch(p)            = cnveg_nitrogenstate_inst%leafn_patch(p)               * m
            cnveg_nitrogenflux_inst%m_frootn_to_litter_patch(p)           = cnveg_nitrogenstate_inst%frootn_patch(p)              * m
            cnveg_nitrogenflux_inst%m_livestemn_to_litter_patch(p)        = cnveg_nitrogenstate_inst%livestemn_patch(p)           * m
            cnveg_nitrogenflux_inst%m_livecrootn_to_litter_patch(p)       = cnveg_nitrogenstate_inst%livecrootn_patch(p)          * m
         else
            cnveg_nitrogenflux_inst%m_leafn_to_litter_patch(p)            = cnveg_nitrogenstate_inst%leafn_patch(p)       * matrix_update_gmn(p,ileaf_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
            cnveg_nitrogenflux_inst%m_frootn_to_litter_patch(p)           = cnveg_nitrogenstate_inst%frootn_patch(p)      * matrix_update_gmn(p,ifroot_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,.true.,.True.)
            cnveg_nitrogenflux_inst%m_livestemn_to_litter_patch(p)        = cnveg_nitrogenstate_inst%livestemn_patch(p)   * matrix_update_gmn(p,ilivestem_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
            cnveg_nitrogenflux_inst%m_livecrootn_to_litter_patch(p)       = cnveg_nitrogenstate_inst%livecrootn_patch(p)  * matrix_update_gmn(p,ilivecroot_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
         end if

         if (spinup_state == 2 .and. .not. use_cndv) then   !accelerate mortality of dead woody pools 
            if(.not. use_matrixcn)then
               cnveg_nitrogenflux_inst%m_deadstemn_to_litter_patch(p)     = cnveg_nitrogenstate_inst%deadstemn_patch(p)  * m * spinup_factor_deadwood
               cnveg_nitrogenflux_inst%m_deadcrootn_to_litter_patch(p)    = cnveg_nitrogenstate_inst%deadcrootn_patch(p) * m * spinup_factor_deadwood
            else
               cnveg_nitrogenflux_inst%m_deadstemn_to_litter_patch(p)     = cnveg_nitrogenstate_inst%deadstemn_patch(p)  * matrix_update_gmn(p,ideadstem_to_iout_gmn , &
                                                                              m*spinup_factor_deadwood,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
               cnveg_nitrogenflux_inst%m_deadcrootn_to_litter_patch(p)    = cnveg_nitrogenstate_inst%deadcrootn_patch(p) * matrix_update_gmn(p,ideadcroot_to_iout_gmn, &
                                                                              m*spinup_factor_deadwood,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
            end if !.not. use_matrixcn
         else
            if (.not. use_matrixcn) then
               cnveg_nitrogenflux_inst%m_deadstemn_to_litter_patch(p)     = cnveg_nitrogenstate_inst%deadstemn_patch(p)           * m 
               cnveg_nitrogenflux_inst%m_deadcrootn_to_litter_patch(p)    = cnveg_nitrogenstate_inst%deadcrootn_patch(p)          * m 
            else
               cnveg_nitrogenflux_inst%m_deadstemn_to_litter_patch(p)     = cnveg_nitrogenstate_inst%deadstemn_patch(p)  * matrix_update_gmn(p,ideadstem_to_iout_gmn ,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
               cnveg_nitrogenflux_inst%m_deadcrootn_to_litter_patch(p)    = cnveg_nitrogenstate_inst%deadcrootn_patch(p) * matrix_update_gmn(p,ideadcroot_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
            end if !use_matrixcn
         end if

         if (ivt(p) < npcropmin) then
            if(.not. use_matrixcn)then
               cnveg_nitrogenflux_inst%m_retransn_to_litter_patch(p) = cnveg_nitrogenstate_inst%retransn_patch(p) * m
            else
               cnveg_nitrogenflux_inst%m_retransn_to_litter_patch(p) = cnveg_nitrogenstate_inst%retransn_patch(p) * matrix_update_gmn(p,iretransn_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
            end if
         end if
            
         if(.not. use_matrixcn)then
            ! storage pools
            cnveg_nitrogenflux_inst%m_leafn_storage_to_litter_patch(p)       = cnveg_nitrogenstate_inst%leafn_storage_patch(p)      * m
            cnveg_nitrogenflux_inst%m_frootn_storage_to_litter_patch(p)      = cnveg_nitrogenstate_inst%frootn_storage_patch(p)     * m
            cnveg_nitrogenflux_inst%m_livestemn_storage_to_litter_patch(p)   = cnveg_nitrogenstate_inst%livestemn_storage_patch(p)  * m
            cnveg_nitrogenflux_inst%m_deadstemn_storage_to_litter_patch(p)   = cnveg_nitrogenstate_inst%deadstemn_storage_patch(p)  * m
            cnveg_nitrogenflux_inst%m_livecrootn_storage_to_litter_patch(p)  = cnveg_nitrogenstate_inst%livecrootn_storage_patch(p) * m
            cnveg_nitrogenflux_inst%m_deadcrootn_storage_to_litter_patch(p)  = cnveg_nitrogenstate_inst%deadcrootn_storage_patch(p) * m

            ! transfer pools
            cnveg_nitrogenflux_inst%m_leafn_xfer_to_litter_patch(p)          = cnveg_nitrogenstate_inst%leafn_xfer_patch(p)         * m
            cnveg_nitrogenflux_inst%m_frootn_xfer_to_litter_patch(p)         = cnveg_nitrogenstate_inst%frootn_xfer_patch(p)        * m
            cnveg_nitrogenflux_inst%m_livestemn_xfer_to_litter_patch(p)      = cnveg_nitrogenstate_inst%livestemn_xfer_patch(p)     * m
            cnveg_nitrogenflux_inst%m_deadstemn_xfer_to_litter_patch(p)      = cnveg_nitrogenstate_inst%deadstemn_xfer_patch(p)     * m
            cnveg_nitrogenflux_inst%m_livecrootn_xfer_to_litter_patch(p)     = cnveg_nitrogenstate_inst%livecrootn_xfer_patch(p)    * m
            cnveg_nitrogenflux_inst%m_deadcrootn_xfer_to_litter_patch(p)     = cnveg_nitrogenstate_inst%deadcrootn_xfer_patch(p)    * m
         else
         ! storage pools
            cnveg_nitrogenflux_inst%m_leafn_storage_to_litter_patch(p)         = cnveg_nitrogenstate_inst%leafn_storage_patch(p)      * matrix_update_gmn(p,ileafst_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
            cnveg_nitrogenflux_inst%m_frootn_storage_to_litter_patch(p)        = cnveg_nitrogenstate_inst%frootn_storage_patch(p)     * matrix_update_gmn(p,ifrootst_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
            cnveg_nitrogenflux_inst%m_livestemn_storage_to_litter_patch(p)     = cnveg_nitrogenstate_inst%livestemn_storage_patch(p)  * matrix_update_gmn(p,ilivestemst_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
            cnveg_nitrogenflux_inst%m_deadstemn_storage_to_litter_patch(p)     = cnveg_nitrogenstate_inst%deadstemn_storage_patch(p)  * matrix_update_gmn(p,ideadstemst_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
            cnveg_nitrogenflux_inst%m_livecrootn_storage_to_litter_patch(p)    = cnveg_nitrogenstate_inst%livecrootn_storage_patch(p) * matrix_update_gmn(p,ilivecrootst_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
            cnveg_nitrogenflux_inst%m_deadcrootn_storage_to_litter_patch(p)    = cnveg_nitrogenstate_inst%deadcrootn_storage_patch(p) * matrix_update_gmn(p,ideadcrootst_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
            
         ! transfer pools
            cnveg_nitrogenflux_inst%m_leafn_xfer_to_litter_patch(p)         = cnveg_nitrogenstate_inst%leafn_xfer_patch(p)      * matrix_update_gmn(p,ileafxf_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
            cnveg_nitrogenflux_inst%m_frootn_xfer_to_litter_patch(p)        = cnveg_nitrogenstate_inst%frootn_xfer_patch(p)     * matrix_update_gmn(p,ifrootxf_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
            cnveg_nitrogenflux_inst%m_livestemn_xfer_to_litter_patch(p)     = cnveg_nitrogenstate_inst%livestemn_xfer_patch(p)  * matrix_update_gmn(p,ilivestemxf_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
            cnveg_nitrogenflux_inst%m_deadstemn_xfer_to_litter_patch(p)     = cnveg_nitrogenstate_inst%deadstemn_xfer_patch(p)  * matrix_update_gmn(p,ideadstemxf_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
            cnveg_nitrogenflux_inst%m_livecrootn_xfer_to_litter_patch(p)    = cnveg_nitrogenstate_inst%livecrootn_xfer_patch(p) * matrix_update_gmn(p,ilivecrootxf_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
            cnveg_nitrogenflux_inst%m_deadcrootn_xfer_to_litter_patch(p)    = cnveg_nitrogenstate_inst%deadcrootn_xfer_patch(p) * matrix_update_gmn(p,ideadcrootxf_to_iout_gmn,m,dt,cnveg_nitrogenflux_inst,matrixcheck_gm,.True.)
         end if !use_matrixcn

         ! added by F. Li and S. Levis
         if (use_cndv) then
            if (woody(ivt(p)) == 1._r8)then
               if (cnveg_carbonstate_inst%livestemc_patch(p) + cnveg_carbonstate_inst%deadstemc_patch(p)> 0._r8)then
                  nind(p)=nind(p)*(1._r8-m)
               else
                  nind(p) = 0._r8 
               end if
            end if
         end if

      end do ! end of patch loop

      ! gather all patch-level litterfall fluxes to the column
      ! for litter C and N inputs

      call CNGap_PatchToColumn(bounds, num_soilp, filter_soilp, &
           cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
           leaf_prof_patch(bounds%begp:bounds%endp, 1:nlevdecomp_full), &
           froot_prof_patch(bounds%begp:bounds%endp, 1:nlevdecomp_full), & 
           croot_prof_patch(bounds%begp:bounds%endp, 1:nlevdecomp_full), &
           stem_prof_patch(bounds%begp:bounds%endp, 1:nlevdecomp_full))
      
    end associate

  end subroutine CNGapMortality

  !-----------------------------------------------------------------------
  subroutine CNGap_PatchToColumn (bounds, num_soilp, filter_soilp, &
       cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
       leaf_prof_patch, froot_prof_patch, croot_prof_patch, stem_prof_patch)
    !
    ! !DESCRIPTION:
    ! gathers all patch-level gap mortality fluxes to the column level and
    ! assigns them to the three litter pools
    !
    ! !USES:
    use clm_varpar , only : maxsoil_patches, nlevdecomp, nlevdecomp_full, i_litr_min, i_litr_max, i_met_lit
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:) ! soil patch filter
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    real(r8)                        , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
    real(r8)                        , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
    real(r8)                        , intent(in)    :: croot_prof_patch(bounds%begp:,1:)
    real(r8)                        , intent(in)    :: stem_prof_patch(bounds%begp:,1:)
    !
    ! !LOCAL VARIABLES:
    integer :: fp,c,p,j,i  ! indices
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(leaf_prof_patch)   == (/bounds%endp,nlevdecomp_full/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(froot_prof_patch)  == (/bounds%endp,nlevdecomp_full/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(croot_prof_patch)  == (/bounds%endp,nlevdecomp_full/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(stem_prof_patch)   == (/bounds%endp,nlevdecomp_full/)), sourcefile, __LINE__)

    associate(                                                                                                 & 
         leaf_prof                           => leaf_prof_patch                                              , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
         froot_prof                          => froot_prof_patch                                             , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
         croot_prof                          => croot_prof_patch                                             , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
         stem_prof                           => stem_prof_patch                                              , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          

         ivt                                 => patch%itype                                                    , & ! Input:  [integer  (:)   ]  patch vegetation type                                
         wtcol                               => patch%wtcol                                                    , & ! Input:  [real(r8) (:)   ]  patch weight relative to column (0-1)               
         
         lf_f                                => pftcon%lf_f                                                  , & ! Input:  [real(r8) (:,:) ]  leaf litter fractions
         fr_f                                => pftcon%fr_f                                                  , & ! Input:  [real(r8) (:,:) ]  fine root litter fractions
         
         m_leafc_to_litter                   => cnveg_carbonflux_inst%m_leafc_to_litter_patch                , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootc_to_litter                  => cnveg_carbonflux_inst%m_frootc_to_litter_patch               , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemc_to_litter               => cnveg_carbonflux_inst%m_livestemc_to_litter_patch            , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemc_to_litter               => cnveg_carbonflux_inst%m_deadstemc_to_litter_patch            , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootc_to_litter              => cnveg_carbonflux_inst%m_livecrootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootc_to_litter              => cnveg_carbonflux_inst%m_deadcrootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_leafc_storage_to_litter           => cnveg_carbonflux_inst%m_leafc_storage_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootc_storage_to_litter          => cnveg_carbonflux_inst%m_frootc_storage_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemc_storage_to_litter       => cnveg_carbonflux_inst%m_livestemc_storage_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemc_storage_to_litter       => cnveg_carbonflux_inst%m_deadstemc_storage_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootc_storage_to_litter      => cnveg_carbonflux_inst%m_livecrootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootc_storage_to_litter      => cnveg_carbonflux_inst%m_deadcrootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                                                    
         m_gresp_storage_to_litter           => cnveg_carbonflux_inst%m_gresp_storage_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
         m_leafc_xfer_to_litter              => cnveg_carbonflux_inst%m_leafc_xfer_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootc_xfer_to_litter             => cnveg_carbonflux_inst%m_frootc_xfer_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemc_xfer_to_litter          => cnveg_carbonflux_inst%m_livestemc_xfer_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemc_xfer_to_litter          => cnveg_carbonflux_inst%m_deadstemc_xfer_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootc_xfer_to_litter         => cnveg_carbonflux_inst%m_livecrootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootc_xfer_to_litter         => cnveg_carbonflux_inst%m_deadcrootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
         m_gresp_xfer_to_litter              => cnveg_carbonflux_inst%m_gresp_xfer_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         gap_mortality_c_to_litr_c           => cnveg_carbonflux_inst%gap_mortality_c_to_litr_c_col          , & ! Output: [real(r8) (:,:,:) ]  C fluxes associated with gap mortality to litter pools (gC/m3/s)
         gap_mortality_c_to_cwdc             => cnveg_carbonflux_inst%gap_mortality_c_to_cwdc_col            , & ! Output: [real(r8) (:,:) ]  C fluxes associated with gap mortality to CWD pool (gC/m3/s)
         
         m_leafn_to_litter                   => cnveg_nitrogenflux_inst%m_leafn_to_litter_patch              , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootn_to_litter                  => cnveg_nitrogenflux_inst%m_frootn_to_litter_patch             , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemn_to_litter               => cnveg_nitrogenflux_inst%m_livestemn_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemn_to_litter               => cnveg_nitrogenflux_inst%m_deadstemn_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootn_to_litter              => cnveg_nitrogenflux_inst%m_livecrootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootn_to_litter              => cnveg_nitrogenflux_inst%m_deadcrootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
         m_retransn_to_litter                => cnveg_nitrogenflux_inst%m_retransn_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_leafn_storage_to_litter           => cnveg_nitrogenflux_inst%m_leafn_storage_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootn_storage_to_litter          => cnveg_nitrogenflux_inst%m_frootn_storage_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemn_storage_to_litter       => cnveg_nitrogenflux_inst%m_livestemn_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemn_storage_to_litter       => cnveg_nitrogenflux_inst%m_deadstemn_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootn_storage_to_litter      => cnveg_nitrogenflux_inst%m_livecrootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootn_storage_to_litter      => cnveg_nitrogenflux_inst%m_deadcrootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
         m_leafn_xfer_to_litter              => cnveg_nitrogenflux_inst%m_leafn_xfer_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootn_xfer_to_litter             => cnveg_nitrogenflux_inst%m_frootn_xfer_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemn_xfer_to_litter          => cnveg_nitrogenflux_inst%m_livestemn_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemn_xfer_to_litter          => cnveg_nitrogenflux_inst%m_deadstemn_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootn_xfer_to_litter         => cnveg_nitrogenflux_inst%m_livecrootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootn_xfer_to_litter         => cnveg_nitrogenflux_inst%m_deadcrootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
         gap_mortality_n_to_litr_n           => cnveg_nitrogenflux_inst%gap_mortality_n_to_litr_n_col        , & ! Output: [real(r8) (:,:,:)]  N fluxes associated with gap mortality to litter pools (gN/m3/s)
         gap_mortality_n_to_cwdn             => cnveg_nitrogenflux_inst%gap_mortality_n_to_cwdn_col            & ! Output: [real(r8) (:,:) ]  N fluxes associated with gap mortality to CWD pool (gN/m3/s)
         )

      do j = 1,nlevdecomp
         do fp = 1,num_soilp
            p = filter_soilp(fp)
            c = patch%column(p)

            do i = i_litr_min, i_litr_max
               gap_mortality_c_to_litr_c(c,j,i) = &
                  gap_mortality_c_to_litr_c(c,j,i) + &
                  ! leaf gap mortality carbon fluxes
                  m_leafc_to_litter(p) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j) + &
                  ! fine root gap mortality carbon fluxes
                  m_frootc_to_litter(p) * fr_f(ivt(p),i) * wtcol(p) * froot_prof(p,j)
            end do

            ! wood gap mortality carbon fluxes
            gap_mortality_c_to_cwdc(c,j)  = gap_mortality_c_to_cwdc(c,j)  + &
                 (m_livestemc_to_litter(p) + m_deadstemc_to_litter(p))  * wtcol(p) * stem_prof(p,j)
            gap_mortality_c_to_cwdc(c,j) = gap_mortality_c_to_cwdc(c,j) + &
                 (m_livecrootc_to_litter(p) + m_deadcrootc_to_litter(p)) * wtcol(p) * croot_prof(p,j)

            ! storage gap mortality carbon fluxes
            ! Metabolic litter is treated differently than other types
            ! of litter, so it gets this additional line after the
            ! most recent loop over all litter types
            gap_mortality_c_to_litr_c(c,j,i_met_lit) = &
               gap_mortality_c_to_litr_c(c,j,i_met_lit) + &
               (m_leafc_storage_to_litter(p) + m_gresp_storage_to_litter(p)) * wtcol(p) * leaf_prof(p,j) + &
               m_frootc_storage_to_litter(p) * wtcol(p) * froot_prof(p,j) + &
               (m_livestemc_storage_to_litter(p) + m_deadstemc_storage_to_litter(p)) * wtcol(p) * stem_prof(p,j) + &
               (m_livecrootc_storage_to_litter(p) + m_deadcrootc_storage_to_litter(p)) * wtcol(p) * croot_prof(p,j) + &

            ! transfer gap mortality carbon fluxes
               (m_leafc_xfer_to_litter(p) + m_gresp_xfer_to_litter(p)) * wtcol(p) * leaf_prof(p,j) + &
               m_frootc_xfer_to_litter(p) * wtcol(p) * froot_prof(p,j) + &
               (m_livestemc_xfer_to_litter(p) + m_deadstemc_xfer_to_litter(p))  * wtcol(p) * stem_prof(p,j) + &
               (m_livecrootc_xfer_to_litter(p) + m_deadcrootc_xfer_to_litter(p)) * wtcol(p) * croot_prof(p,j)

            do i = i_litr_min, i_litr_max
               gap_mortality_n_to_litr_n(c,j,i) = &
                  gap_mortality_n_to_litr_n(c,j,i) + &
                  ! leaf gap mortality nitrogen fluxes
                  m_leafn_to_litter(p) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j) + &
                  ! fine root litter nitrogen fluxes
                  m_frootn_to_litter(p) * fr_f(ivt(p),i) * wtcol(p) * froot_prof(p,j)
            end do

            ! wood gap mortality nitrogen fluxes
            gap_mortality_n_to_cwdn(c,j) = gap_mortality_n_to_cwdn(c,j)  + &
                 (m_livestemn_to_litter(p) + m_deadstemn_to_litter(p))  * wtcol(p) * stem_prof(p,j)
            gap_mortality_n_to_cwdn(c,j) = gap_mortality_n_to_cwdn(c,j) + &
                 (m_livecrootn_to_litter(p) + m_deadcrootn_to_litter(p)) * wtcol(p) * croot_prof(p,j)

            ! Metabolic litter is treated differently than other types
            ! of litter, so it gets this additional line after the
            ! most recent loop over all litter types
            gap_mortality_n_to_litr_n(c,j,i_met_lit) = &
               gap_mortality_n_to_litr_n(c,j,i_met_lit) + &
               ! retranslocated N pool gap mortality fluxes
               m_retransn_to_litter(p) * wtcol(p) * leaf_prof(p,j) + &
               ! storage gap mortality nitrogen fluxes
               m_leafn_storage_to_litter(p) * wtcol(p) * leaf_prof(p,j) + &
               m_frootn_storage_to_litter(p) * wtcol(p) * froot_prof(p,j) + &
               (m_livestemn_storage_to_litter(p) + m_deadstemn_storage_to_litter(p))  * wtcol(p) * stem_prof(p,j) + &
               (m_livecrootn_storage_to_litter(p) + m_deadcrootn_storage_to_litter(p)) * wtcol(p) * croot_prof(p,j) + &
               ! transfer gap mortality nitrogen fluxes
               m_leafn_xfer_to_litter(p) * wtcol(p) * leaf_prof(p,j) + &
               m_frootn_xfer_to_litter(p) * wtcol(p) * froot_prof(p,j) + &
               (m_livestemn_xfer_to_litter(p) + m_deadstemn_xfer_to_litter(p)) * wtcol(p) * stem_prof(p,j) + &
               (m_livecrootn_xfer_to_litter(p) + m_deadcrootn_xfer_to_litter(p)) * wtcol(p) * croot_prof(p,j)

         end do
      end do

    end associate

  end subroutine CNGap_PatchToColumn

end module CNGapMortalityMod
