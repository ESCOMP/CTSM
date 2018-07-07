module CNFireBaseMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! module for fire dynamics 
  ! created in Nov, 2012  and revised in Apr, 2013 by F. Li and S. Levis
  ! based on Li et al. (2012a,b; 2013)
  ! revised in Apr, 2014 according Li et al.(2014)
  ! Fire-related parameters were calibrated or tuned in Apr, 2013 based on the 
  ! 20th Century transient simulations at f19_g16 with (newfire05_clm45sci15_clm4_0_58) 
  ! a CLM4.5 version, Qian et al. (2006) atmospheric forcing, and
  ! climatological lightning data.
  !
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_strdata_mod                    , only : shr_strdata_type, shr_strdata_create, shr_strdata_print
  use shr_strdata_mod                    , only : shr_strdata_advance
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varctl                         , only : iulog
  use spmdMod                            , only : masterproc, mpicom, comp_id
  use fileutils                          , only : getavu, relavu
  use decompMod                          , only : gsmap_lnd_gdc2glo
  use domainMod                          , only : ldomain
  use pftconMod                          , only : noveg, pftcon
  use abortutils                         , only : endrun
  use decompMod                          , only : bounds_type
  use atm2lndType                        , only : atm2lnd_type
  use CNDVType                           , only : dgvs_type
  use CNVegStateType                     , only : cnveg_state_type
  use CNVegCarbonStateType               , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType                , only : cnveg_carbonflux_type
  use CNVegNitrogenStateType             , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType              , only : cnveg_nitrogenflux_type
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use EnergyFluxType                     , only : energyflux_type
  use SaturatedExcessRunoffMod           , only : saturated_excess_runoff_type
  use WaterDiagnosticBulkType                     , only : waterdiagnosticbulk_type
  use GridcellType                       , only : grc
  use ColumnType                         , only : col
  use PatchType                          , only : patch
  use mct_mod
  use CNFireMethodMod                    , only : cnfire_method_type
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: cnfire_base_type

  integer, private, parameter :: num_fp = 2   ! Number of pools relevent for fire
  integer, private, parameter :: lit_fp = 1   ! Pool for liter
  integer, private, parameter :: cwd_fp = 2   ! Pool for CWD Course woody debris
  type, public :: cnfire_const_type
      ! !PRIVATE MEMBER DATA:
      real(r8) :: borealat = 40._r8                    ! Latitude for boreal peat fires
      real(r8) :: lfuel=75._r8                         ! lower threshold of fuel mass (gC/m2) for ignition, Li et al.(2014)
      real(r8) :: ufuel=650._r8                        ! upper threshold of fuel mass(gC/m2) for ignition 
      real(r8) :: g0=0.05_r8                           ! g(W) when W=0 m/s
      real(r8) :: rh_low=30.0_r8                       ! Relative humidty low (%)
      real(r8) :: rh_hgh=80.0_r8                       ! Relative humidty high (%)
      real(r8) :: bt_min=0.3_r8                        ! btran minimum (fraction)
      real(r8) :: bt_max=0.7_r8                        ! btran maximum (fraction)
      real(r8) :: cli_scale=0.035_r8                   ! global constant for deforestation fires (/d)
      real(r8) :: boreal_peatfire_c = 4.2e-5_r8        ! c parameter for boreal peatland fire in Li et. al. (2013) (/hr)
      real(r8) :: pot_hmn_ign_counts_alpha=0.0035_r8   ! Potential human ignition counts (alpha in Li et. al. 2012) (/person/month)
      real(r8) :: non_boreal_peatfire_c  = 0.001_r8    ! c parameter for non-boreal peatland fire in Li et. al. (2013) (/hr)
      real(r8) :: cropfire_a1 = 0.3_r8                 ! a1 parameter for cropland fire in (Li et. al., 2014) (/hr)
      real(r8) :: occur_hi_gdp_tree = 0.39_r8          ! fire occurance for high GDP areas that are tree dominated (fraction)

      real(r8) :: cmb_cmplt_fact(num_fp) = (/ 0.5_r8, 0.25_r8 /) ! combustion completion factor (unitless)
  end type

  !
  type, extends(cnfire_method_type) :: cnfire_base_type
    private
      ! !PRIVATE MEMBER DATA:

      real(r8), public, pointer :: forc_lnfm(:)        ! Lightning frequency
      real(r8), public, pointer :: forc_hdm(:)         ! Human population density

      type(shr_strdata_type) :: sdat_hdm    ! Human population density input data stream
      type(shr_strdata_type) :: sdat_lnfm   ! Lightning input data stream


    contains
      !
      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public :: CNFireInit        ! Initialization of CNFire
      procedure, public :: CNFireReadNML     ! Read in namelist for CNFire
      procedure, public :: CNFireInterp      ! Interpolate fire data
      procedure, public :: CNFireArea        ! Calculate fire area
      procedure, public :: CNFireFluxes      ! Calculate fire fluxes
      !
      ! !PRIVATE MEMBER FUNCTIONS:
      procedure, private :: hdm_init    ! position datasets for dynamic human population density
      procedure, private :: hdm_interp  ! interpolates between two years of human pop. density file data
      procedure, private :: lnfm_init   ! position datasets for Lightning
      procedure, private :: lnfm_interp ! interpolates between two years of Lightning file data
  end type cnfire_base_type
  !-----------------------------------------------------------------------

  type(cnfire_const_type), public, protected :: cnfire_const          ! Fire constants shared by Li versons

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine CNFireInit( this, bounds, NLFilename )
    !
    ! !DESCRIPTION:
    ! Initialize CN Fire module
    ! !USES:
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(cnfire_base_type) :: this
    type(bounds_type), intent(in) :: bounds  
    character(len=*),  intent(in) :: NLFilename
    !-----------------------------------------------------------------------

    if ( this%need_lightning_and_popdens ) then
       ! Allocate lightning forcing data
       allocate( this%forc_lnfm(bounds%begg:bounds%endg) )
       this%forc_lnfm(bounds%begg:) = nan
       ! Allocate pop dens forcing data
       allocate( this%forc_hdm(bounds%begg:bounds%endg) )
       this%forc_hdm(bounds%begg:) = nan

       call this%hdm_init(bounds, NLFilename)
       call this%hdm_interp(bounds)
       call this%lnfm_init(bounds, NLFilename)
       call this%lnfm_interp(bounds)
    end if

  end subroutine CNFireInit

  !-----------------------------------------------------------------------
  subroutine CNFireReadNML( this, NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for CNFire
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    class(cnfire_base_type) :: this
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'CNFireReadNML'
    character(len=*), parameter :: nmlname = 'lifire_inparm'
    !-----------------------------------------------------------------------
    real(r8) :: cli_scale, boreal_peatfire_c, pot_hmn_ign_counts_alpha
    real(r8) :: non_boreal_peatfire_c, cropfire_a1
    real(r8) :: rh_low, rh_hgh, bt_min, bt_max, occur_hi_gdp_tree
    real(r8) :: lfuel, ufuel, cmb_cmplt_fact(num_fp)

    namelist /lifire_inparm/ cli_scale, boreal_peatfire_c, pot_hmn_ign_counts_alpha, &
                             non_boreal_peatfire_c, cropfire_a1,                &
                             rh_low, rh_hgh, bt_min, bt_max, occur_hi_gdp_tree, &
                             lfuel, ufuel, cmb_cmplt_fact

    if ( this%need_lightning_and_popdens ) then
       cli_scale                 = cnfire_const%cli_scale
       boreal_peatfire_c         = cnfire_const%boreal_peatfire_c
       non_boreal_peatfire_c     = cnfire_const%non_boreal_peatfire_c
       pot_hmn_ign_counts_alpha  = cnfire_const%pot_hmn_ign_counts_alpha
       cropfire_a1               = cnfire_const%cropfire_a1
       rh_low                    = cnfire_const%rh_low
       rh_hgh                    = cnfire_const%rh_hgh
       lfuel                     = cnfire_const%lfuel
       ufuel                     = cnfire_const%ufuel
       bt_min                    = cnfire_const%bt_min
       bt_max                    = cnfire_const%bt_max
       occur_hi_gdp_tree         = cnfire_const%occur_hi_gdp_tree
       cmb_cmplt_fact(:)         = cnfire_const%cmb_cmplt_fact(:)
       ! Initialize options to default values, in case they are not specified in
       ! the namelist

       if (masterproc) then
          unitn = getavu()
          write(iulog,*) 'Read in '//nmlname//'  namelist'
          call opnfil (NLFilename, unitn, 'F')
          call shr_nl_find_group_name(unitn, nmlname, status=ierr)
          if (ierr == 0) then
             read(unitn, nml=lifire_inparm, iostat=ierr)
             if (ierr /= 0) then
                call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
             end if
          else
             call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
          call relavu( unitn )
       end if

       call shr_mpi_bcast (cli_scale               , mpicom)
       call shr_mpi_bcast (boreal_peatfire_c       , mpicom)
       call shr_mpi_bcast (pot_hmn_ign_counts_alpha, mpicom)
       call shr_mpi_bcast (non_boreal_peatfire_c   , mpicom)
       call shr_mpi_bcast (cropfire_a1             , mpicom)
       call shr_mpi_bcast (rh_low                  , mpicom)
       call shr_mpi_bcast (rh_hgh                  , mpicom)
       call shr_mpi_bcast (lfuel                   , mpicom)
       call shr_mpi_bcast (ufuel                   , mpicom)
       call shr_mpi_bcast (bt_min                  , mpicom)
       call shr_mpi_bcast (bt_max                  , mpicom)
       call shr_mpi_bcast (occur_hi_gdp_tree       , mpicom)
       call shr_mpi_bcast (cmb_cmplt_fact          , mpicom)

       cnfire_const%cli_scale                 = cli_scale
       cnfire_const%boreal_peatfire_c         = boreal_peatfire_c
       cnfire_const%non_boreal_peatfire_c     = non_boreal_peatfire_c
       cnfire_const%pot_hmn_ign_counts_alpha  = pot_hmn_ign_counts_alpha
       cnfire_const%cropfire_a1               = cropfire_a1
       cnfire_const%rh_low                    = rh_low
       cnfire_const%rh_hgh                    = rh_hgh
       cnfire_const%lfuel                     = lfuel
       cnfire_const%ufuel                     = ufuel
       cnfire_const%bt_min                    = bt_min
       cnfire_const%bt_max                    = bt_max
       cnfire_const%occur_hi_gdp_tree         = occur_hi_gdp_tree
       cnfire_const%cmb_cmplt_fact(:)         = cmb_cmplt_fact(:)

       if (masterproc) then
          write(iulog,*) ' '
          write(iulog,*) nmlname//' settings:'
          write(iulog,nml=lifire_inparm)
          write(iulog,*) ' '
       end if
    end if

  end subroutine CNFireReadNML

  !-----------------------------------------------------------------------
  subroutine CNFireInterp(this,bounds)
    !
    ! !DESCRIPTION:
    ! Interpolate CN Fire datasets
    !
    ! !ARGUMENTS:
    class(cnfire_base_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !-----------------------------------------------------------------------

    if ( this%need_lightning_and_popdens ) then
       call this%hdm_interp(bounds)
       call this%lnfm_interp(bounds)
    end if

  end subroutine CNFireInterp

  !-----------------------------------------------------------------------
  subroutine CNFireArea (this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       atm2lnd_inst, energyflux_inst, saturated_excess_runoff_inst, waterdiagnosticbulk_inst, &
       cnveg_state_inst, cnveg_carbonstate_inst, totlitc_col, decomp_cpools_vr_col, t_soi17cm_col)
    !
    ! !DESCRIPTION:
    ! Computes column-level burned area 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cnfire_base_type)                               :: this
    type(bounds_type)                     , intent(in)    :: bounds 
    integer                               , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                               , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                               , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                               , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(atm2lnd_type)                    , intent(in)    :: atm2lnd_inst
    type(energyflux_type)                 , intent(in)    :: energyflux_inst
    type(saturated_excess_runoff_type)    , intent(in)    :: saturated_excess_runoff_inst
    type(waterdiagnosticbulk_type)                 , intent(in)    :: waterdiagnosticbulk_inst
    type(cnveg_state_type)                , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)          , intent(inout) :: cnveg_carbonstate_inst
    real(r8)                              , intent(in)    :: totlitc_col(bounds%begc:)
    real(r8)                              , intent(in)    :: decomp_cpools_vr_col(bounds%begc:,1:,1:)
    real(r8)                              , intent(in)    :: t_soi17cm_col(bounds%begc:)
    !

    call endrun( 'cnfire_base::CNFireArea: this method MUST be implemented!' )

  end subroutine CNFireArea

  !-----------------------------------------------------------------------
  subroutine CNFireFluxes (this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
      dgvs_inst, cnveg_state_inst,                                                                      &
      cnveg_carbonstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
      leaf_prof_patch, froot_prof_patch, croot_prof_patch, stem_prof_patch, &
      totsomc_col, decomp_cpools_vr_col, decomp_npools_vr_col, somc_fire_col)
   !
   ! !DESCRIPTION:
   ! Fire effects routine for coupled carbon-nitrogen code (CN).
   ! Relies primarily on estimate of fractional area burned, from CNFireArea().
   !
   ! Total fire carbon emissions (g C/m2 land area/yr) 
   !  =avg(COL_FIRE_CLOSS)*seconds_per_year + avg(SOMC_FIRE)*seconds_per_year + 
   !   avg(LF_CONV_CFLUX)*seconds_per_year*min(1.0,avg(LFC2)*seconds_per_year)*0.8
   ! where avg means the temporal average in a year
   ! seconds_per_year is the number of seconds in a year.
   !
   ! !USES:
   use clm_time_manager     , only: get_step_size,get_days_per_year,get_curr_date
   use clm_varpar           , only: max_patch_per_col
   use clm_varctl           , only: use_cndv, spinup_state
   use clm_varcon           , only: secspday
   use pftconMod            , only: nc3crop
   use dynSubgridControlMod , only: run_has_transient_landcover
   use clm_varpar           , only: nlevdecomp_full, ndecomp_pools, nlevdecomp
   !
   ! !ARGUMENTS:
   class(cnfire_base_type)                        :: this
   type(bounds_type)              , intent(in)    :: bounds  
   integer                        , intent(in)    :: num_soilc       ! number of soil columns in filter
   integer                        , intent(in)    :: filter_soilc(:) ! filter for soil columns
   integer                        , intent(in)    :: num_soilp       ! number of soil patches in filter
   integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
   type(dgvs_type)                , intent(inout) :: dgvs_inst
   type(cnveg_state_type)         , intent(inout) :: cnveg_state_inst
   type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst
   type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
   type(cnveg_nitrogenstate_type) , intent(in)    :: cnveg_nitrogenstate_inst
   type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
   real(r8)                       , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
   real(r8)                       , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
   real(r8)                       , intent(in)    :: croot_prof_patch(bounds%begp:,1:)
   real(r8)                       , intent(in)    :: stem_prof_patch(bounds%begp:,1:)
   real(r8)                       , intent(in)    :: totsomc_col(bounds%begc:)                ! (gC/m2) total soil organic matter C
   real(r8)                       , intent(in)    :: decomp_cpools_vr_col(bounds%begc:,1:,1:) ! (gC/m3)  VR decomp. (litter, cwd, soil)
   real(r8)                       , intent(in)    :: decomp_npools_vr_col(bounds%begc:,1:,1:) ! (gC/m3)  VR decomp. (litter, cwd, soil)
   real(r8)                       , intent(out)   :: somc_fire_col(bounds%begc:)              ! (gC/m2/s) fire C emissions due to peat burning
   !
   ! !LOCAL VARIABLES:
   integer :: g,c,p,j,l,pi,kyr, kmo, kda, mcsec   ! indices
   integer :: fp,fc                ! filter indices
   real(r8):: f                    ! rate for fire effects (1/s)
   real(r8):: m                    ! acceleration factor for fuel carbon
   real(r8):: dt                   ! time step variable (s)
   real(r8):: dayspyr              ! days per year
   logical :: transient_landcover  ! whether this run has any prescribed transient landcover
   !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(leaf_prof_patch)      == (/bounds%endp,nlevdecomp_full/))               , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(froot_prof_patch)     == (/bounds%endp,nlevdecomp_full/))               , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(croot_prof_patch)     == (/bounds%endp,nlevdecomp_full/))               , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(stem_prof_patch)      == (/bounds%endp,nlevdecomp_full/))               , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(totsomc_col)          == (/bounds%endc/))                               , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(decomp_cpools_vr_col) == (/bounds%endc,nlevdecomp_full,ndecomp_pools/)) , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(decomp_npools_vr_col) == (/bounds%endc,nlevdecomp_full,ndecomp_pools/)) , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(somc_fire_col)        == (/bounds%endc/))                               , errMsg(sourcefile, __LINE__))

   ! NOTE: VR      = Vertically Resolved
   !       conv.   = conversion
   !       frac.   = fraction
   !       BAF     = Burned Area Fraction
   !       ann.    = annual
   !       GC      = gridcell
   !       dt      = timestep
   !       C       = Carbon
   !       N       = Nitrogen
   !       emis.   = emissions
   !       decomp. = decomposing

    associate(                                                                                                      & 
         croot_prof                          => croot_prof_patch                                                  , & ! Input:  [real(r8) (:,:)   ]  (1/m) profile of coarse roots                   
         stem_prof                           => stem_prof_patch                                                   , & ! Input:  [real(r8) (:,:)   ]  (1/m) profile of stems                          
         froot_prof                          => froot_prof_patch                                                  , & ! Input:  [real(r8) (:,:)   ]  (1/m) profile of fine roots                     
         leaf_prof                           => leaf_prof_patch                                                   , & ! Input:  [real(r8) (:,:)   ]  (1/m) profile of leaves                         
         totsomc                             => totsomc_col                                                       , & ! Input:  [real(r8) (:)     ]  (gC/m2) total soil organic matter C
         decomp_cpools_vr                    => decomp_cpools_vr_col                                              , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  VR decomp. (litter, cwd, soil)
         decomp_npools_vr                    => decomp_npools_vr_col                                              , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  VR decomp. (litter, cwd, soil)
         somc_fire                           => somc_fire_col                                                     , & ! Output: [real(r8) (:)     ]  (gC/m2/s) fire C emissions due to peat burning
         
         is_cwd                              => decomp_cascade_con%is_cwd                                         , & ! Input:  [logical  (:)     ]  TRUE => pool is a cwd pool                         
         is_litter                           => decomp_cascade_con%is_litter                                      , & ! Input:  [logical  (:)     ]  TRUE => pool is a litter pool                      
         
         woody                               => pftcon%woody                                                      , & ! Input:  woody lifeform (1=woody, 0=not woody)             
         cc_leaf                             => pftcon%cc_leaf                                                    , & ! Input: 
         cc_lstem                            => pftcon%cc_lstem                                                   , & ! Input: 
         cc_dstem                            => pftcon%cc_dstem                                                   , & ! Input: 
         cc_other                            => pftcon%cc_other                                                   , & ! Input: 
         fm_leaf                             => pftcon%fm_leaf                                                    , & ! Input: 
         fm_lstem                            => pftcon%fm_lstem                                                   , & ! Input: 
         fm_other                            => pftcon%fm_other                                                   , & ! Input: 
         fm_root                             => pftcon%fm_root                                                    , & ! Input: 
         fm_lroot                            => pftcon%fm_lroot                                                   , & ! Input: 
         fm_droot                            => pftcon%fm_droot                                                   , & ! Input: 
         lf_flab                             => pftcon%lf_flab                                                    , & ! Input: 
         lf_fcel                             => pftcon%lf_fcel                                                    , & ! Input: 
         lf_flig                             => pftcon%lf_flig                                                    , & ! Input: 
         fr_flab                             => pftcon%fr_flab                                                    , & ! Input: 
         fr_fcel                             => pftcon%fr_fcel                                                    , & ! Input: 
         fr_flig                             => pftcon%fr_flig                                                    , & ! Input: 

         cmb_cmplt_fact                      => cnfire_const%cmb_cmplt_fact                                       , & ! Input:  [real(r8) (:)     ]  Combustion completion factor (unitless)
         
         nind                                => dgvs_inst%nind_patch                                              , & ! Input:  [real(r8) (:)     ]  number of individuals (#/m2)                      
         
         cropf_col                           => cnveg_state_inst%cropf_col                                        , & ! Input:  [real(r8) (:)     ]  cropland fraction in veg column                   
         farea_burned                        => cnveg_state_inst%farea_burned_col                                 , & ! Input:  [real(r8) (:)     ]  fractional area burned (/sec)
         fbac1                               => cnveg_state_inst%fbac1_col                                        , & ! Input:  [real(r8) (:)     ]  burned area out of conv. region due to LU fire 
         fbac                                => cnveg_state_inst%fbac_col                                         , & ! Input:  [real(r8) (:)     ]  total burned area out of conversion (/sec)
         baf_crop                            => cnveg_state_inst%baf_crop_col                                     , & ! Input:  [real(r8) (:)     ]  BAF for cropland                                  
         baf_peatf                           => cnveg_state_inst%baf_peatf_col                                    , & ! Input:  [real(r8) (:)     ]  BAF for peatlabd                                  
         trotr1_col                          => cnveg_state_inst%trotr1_col                                       , & ! Input:  [real(r8) (:)     ]  patch weight of BET on the column (0-1)           
         trotr2_col                          => cnveg_state_inst%trotr2_col                                       , & ! Input:  [real(r8) (:)     ]  patch weight of BDT on the column (0-1)           
         dtrotr_col                          => cnveg_state_inst%dtrotr_col                                       , & ! Input:  [real(r8) (:)     ]  ann. decreased frac. coverage of BET+BDT (0-1) on GC
         lfc                                 => cnveg_state_inst%lfc_col                                          , & ! Input:  [real(r8) (:)     ]  conv. area frac. of BET+BDT that haven't burned before
         lfc2                                => cnveg_state_inst%lfc2_col                                         , & ! Output: [real(r8) (:)     ]  conv. area frac. of BET+BDT burned this dt (/sec)
         
         leafcmax                            => cnveg_carbonstate_inst%leafcmax_patch                             , & ! Output: [real(r8) (:)     ]  (gC/m2) ann max leaf C                            
         leafc                               => cnveg_carbonstate_inst%leafc_patch                                , & ! Input:  [real(r8) (:)     ]  (gC/m2) leaf C                                    
         leafc_storage                       => cnveg_carbonstate_inst%leafc_storage_patch                        , & ! Input:  [real(r8) (:)     ]  (gC/m2) leaf C storage                            
         leafc_xfer                          => cnveg_carbonstate_inst%leafc_xfer_patch                           , & ! Input:  [real(r8) (:)     ]  (gC/m2) leaf C transfer                           
         livestemc                           => cnveg_carbonstate_inst%livestemc_patch                            , & ! Input:  [real(r8) (:)     ]  (gC/m2) live stem C                               
         livestemc_storage                   => cnveg_carbonstate_inst%livestemc_storage_patch                    , & ! Input:  [real(r8) (:)     ]  (gC/m2) live stem C storage                       
         livestemc_xfer                      => cnveg_carbonstate_inst%livestemc_xfer_patch                       , & ! Input:  [real(r8) (:)     ]  (gC/m2) live stem C transfer                      
         deadstemc                           => cnveg_carbonstate_inst%deadstemc_patch                            , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead stem C                               
         deadstemc_storage                   => cnveg_carbonstate_inst%deadstemc_storage_patch                    , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead stem C storage                       
         deadstemc_xfer                      => cnveg_carbonstate_inst%deadstemc_xfer_patch                       , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead stem C transfer                      
         frootc                              => cnveg_carbonstate_inst%frootc_patch                               , & ! Input:  [real(r8) (:)     ]  (gC/m2) fine root C                               
         frootc_storage                      => cnveg_carbonstate_inst%frootc_storage_patch                       , & ! Input:  [real(r8) (:)     ]  (gC/m2) fine root C storage                       
         frootc_xfer                         => cnveg_carbonstate_inst%frootc_xfer_patch                          , & ! Input:  [real(r8) (:)     ]  (gC/m2) fine root C transfer                      
         livecrootc                          => cnveg_carbonstate_inst%livecrootc_patch                           , & ! Input:  [real(r8) (:)     ]  (gC/m2) live coarse root C                        
         livecrootc_storage                  => cnveg_carbonstate_inst%livecrootc_storage_patch                   , & ! Input:  [real(r8) (:)     ]  (gC/m2) live coarse root C storage                
         livecrootc_xfer                     => cnveg_carbonstate_inst%livecrootc_xfer_patch                      , & ! Input:  [real(r8) (:)     ]  (gC/m2) live coarse root C transfer               
         deadcrootc                          => cnveg_carbonstate_inst%deadcrootc_patch                           , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead coarse root C                        
         deadcrootc_storage                  => cnveg_carbonstate_inst%deadcrootc_storage_patch                   , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead coarse root C storage                
         deadcrootc_xfer                     => cnveg_carbonstate_inst%deadcrootc_xfer_patch                      , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer               
         gresp_storage                       => cnveg_carbonstate_inst%gresp_storage_patch                        , & ! Input:  [real(r8) (:)     ]  (gC/m2) growth respiration storage                
         gresp_xfer                          => cnveg_carbonstate_inst%gresp_xfer_patch                           , & ! Input:  [real(r8) (:)     ]  (gC/m2) growth respiration transfer               
         
         leafn                               => cnveg_nitrogenstate_inst%leafn_patch                              , & ! Input:  [real(r8) (:)     ]  (gN/m2) leaf N                                    
         leafn_storage                       => cnveg_nitrogenstate_inst%leafn_storage_patch                      , & ! Input:  [real(r8) (:)     ]  (gN/m2) leaf N storage                            
         leafn_xfer                          => cnveg_nitrogenstate_inst%leafn_xfer_patch                         , & ! Input:  [real(r8) (:)     ]  (gN/m2) leaf N transfer                           
         livestemn                           => cnveg_nitrogenstate_inst%livestemn_patch                          , & ! Input:  [real(r8) (:)     ]  (gN/m2) live stem N                               
         livestemn_storage                   => cnveg_nitrogenstate_inst%livestemn_storage_patch                  , & ! Input:  [real(r8) (:)     ]  (gN/m2) live stem N storage                       
         livestemn_xfer                      => cnveg_nitrogenstate_inst%livestemn_xfer_patch                     , & ! Input:  [real(r8) (:)     ]  (gN/m2) live stem N transfer                      
         deadstemn                           => cnveg_nitrogenstate_inst%deadstemn_patch                          , & ! Input:  [real(r8) (:)     ]  (gN/m2) dead stem N                               
         deadstemn_storage                   => cnveg_nitrogenstate_inst%deadstemn_storage_patch                  , & ! Input:  [real(r8) (:)     ]  (gN/m2) dead stem N storage                       
         deadstemn_xfer                      => cnveg_nitrogenstate_inst%deadstemn_xfer_patch                     , & ! Input:  [real(r8) (:)     ]  (gN/m2) dead stem N transfer                      
         frootn                              => cnveg_nitrogenstate_inst%frootn_patch                             , & ! Input:  [real(r8) (:)     ]  (gN/m2) fine root N                               
         frootn_storage                      => cnveg_nitrogenstate_inst%frootn_storage_patch                     , & ! Input:  [real(r8) (:)     ]  (gN/m2) fine root N storage                       
         frootn_xfer                         => cnveg_nitrogenstate_inst%frootn_xfer_patch                        , & ! Input:  [real(r8) (:)     ]  (gN/m2) fine root N transfer                      
         livecrootn                          => cnveg_nitrogenstate_inst%livecrootn_patch                         , & ! Input:  [real(r8) (:)     ]  (gN/m2) live coarse root N                        
         livecrootn_storage                  => cnveg_nitrogenstate_inst%livecrootn_storage_patch                 , & ! Input:  [real(r8) (:)     ]  (gN/m2) live coarse root N storage                
         livecrootn_xfer                     => cnveg_nitrogenstate_inst%livecrootn_xfer_patch                    , & ! Input:  [real(r8) (:)     ]  (gN/m2) live coarse root N transfer               
         deadcrootn                          => cnveg_nitrogenstate_inst%deadcrootn_patch                         , & ! Input:  [real(r8) (:)     ]  (gN/m2) dead coarse root N                        
         deadcrootn_storage                  => cnveg_nitrogenstate_inst%deadcrootn_storage_patch                 , & ! Input:  [real(r8) (:)     ]  (gN/m2) dead coarse root N storage                
         deadcrootn_xfer                     => cnveg_nitrogenstate_inst%deadcrootn_xfer_patch                    , & ! Input:  [real(r8) (:)     ]  (gN/m2) dead coarse root N transfer               
         retransn                            => cnveg_nitrogenstate_inst%retransn_patch                           , & ! Input:  [real(r8) (:)     ]  (gN/m2) plant pool of retranslocated N            
         
         fire_mortality_c_to_cwdc            => cnveg_carbonflux_inst%fire_mortality_c_to_cwdc_col                , & ! Input:  [real(r8) (:,:)   ]  C flux fire mortality to CWD (gC/m3/s)
         m_leafc_to_fire                     => cnveg_carbonflux_inst%m_leafc_to_fire_patch                       , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) fire C emis. from leafc	    
         m_leafc_storage_to_fire             => cnveg_carbonflux_inst%m_leafc_storage_to_fire_patch               , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) fire C emis. from leafc_storage   
         m_leafc_xfer_to_fire                => cnveg_carbonflux_inst%m_leafc_xfer_to_fire_patch                  , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) fire C emis. from leafc_xfer	    
         m_livestemc_to_fire                 => cnveg_carbonflux_inst%m_livestemc_to_fire_patch                   , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) fire C emis. from livestemc	         
         m_livestemc_storage_to_fire         => cnveg_carbonflux_inst%m_livestemc_storage_to_fire_patch           , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. livestemc_storage	       
         m_livestemc_xfer_to_fire            => cnveg_carbonflux_inst%m_livestemc_xfer_to_fire_patch              , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. livestemc_xfer	       
         m_deadstemc_to_fire                 => cnveg_carbonflux_inst%m_deadstemc_to_fire_patch                   , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. deadstemc_xfer	       
         m_deadstemc_storage_to_fire         => cnveg_carbonflux_inst%m_deadstemc_storage_to_fire_patch           , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. deadstemc_storage	       
         m_deadstemc_xfer_to_fire            => cnveg_carbonflux_inst%m_deadstemc_xfer_to_fire_patch              , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. deadstemc_xfer	       
         m_frootc_to_fire                    => cnveg_carbonflux_inst%m_frootc_to_fire_patch                      , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. frootc		       
         m_frootc_storage_to_fire            => cnveg_carbonflux_inst%m_frootc_storage_to_fire_patch              , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. frootc_storage	       
         m_frootc_xfer_to_fire               => cnveg_carbonflux_inst%m_frootc_xfer_to_fire_patch                 , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. frootc_xfer		       
         m_livecrootc_to_fire                => cnveg_carbonflux_inst%m_livecrootc_to_fire_patch                  , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. livecrootc		    	
         m_livecrootc_storage_to_fire        => cnveg_carbonflux_inst%m_livecrootc_storage_to_fire_patch          , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. livecrootc_storage	       
         m_livecrootc_xfer_to_fire           => cnveg_carbonflux_inst%m_livecrootc_xfer_to_fire_patch             , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. livecrootc_xfer	       
         m_deadcrootc_to_fire                => cnveg_carbonflux_inst%m_deadcrootc_to_fire_patch                  , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. deadcrootc		    	
         m_deadcrootc_storage_to_fire        => cnveg_carbonflux_inst%m_deadcrootc_storage_to_fire_patch          , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. deadcrootc_storage	       
         m_deadcrootc_xfer_to_fire           => cnveg_carbonflux_inst%m_deadcrootc_xfer_to_fire_patch             , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. deadcrootc_xfer	       
         m_gresp_storage_to_fire             => cnveg_carbonflux_inst%m_gresp_storage_to_fire_patch               , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. gresp_storage	
         m_gresp_xfer_to_fire                => cnveg_carbonflux_inst%m_gresp_xfer_to_fire_patch                  , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) C emis. gresp_xfer           
         m_leafc_to_litter_fire              => cnveg_carbonflux_inst%m_leafc_to_litter_fire_patch                , & ! Output: [real(r8) (:)     ]                                                    
         m_leafc_storage_to_litter_fire      => cnveg_carbonflux_inst%m_leafc_storage_to_litter_fire_patch        , & ! Output: [real(r8) (:)     ]                                                    
         m_leafc_xfer_to_litter_fire         => cnveg_carbonflux_inst%m_leafc_xfer_to_litter_fire_patch           , & ! Output: [real(r8) (:)     ]                                                    
         m_livestemc_to_litter_fire          => cnveg_carbonflux_inst%m_livestemc_to_litter_fire_patch            , & ! Output: [real(r8) (:)     ]                                                    
         m_livestemc_storage_to_litter_fire  => cnveg_carbonflux_inst%m_livestemc_storage_to_litter_fire_patch    , & ! Output: [real(r8) (:)     ]                                                    
         m_livestemc_xfer_to_litter_fire     => cnveg_carbonflux_inst%m_livestemc_xfer_to_litter_fire_patch       , & ! Output: [real(r8) (:)     ]                                                    
         m_livestemc_to_deadstemc_fire       => cnveg_carbonflux_inst%m_livestemc_to_deadstemc_fire_patch         , & ! Output: [real(r8) (:)     ]                                                    
         m_deadstemc_to_litter_fire          => cnveg_carbonflux_inst%m_deadstemc_to_litter_fire_patch            , & ! Output: [real(r8) (:)     ]                                                    
         m_deadstemc_storage_to_litter_fire  => cnveg_carbonflux_inst%m_deadstemc_storage_to_litter_fire_patch    , & ! Output: [real(r8) (:)     ]                                                    
         m_deadstemc_xfer_to_litter_fire     => cnveg_carbonflux_inst%m_deadstemc_xfer_to_litter_fire_patch       , & ! Output: [real(r8) (:)     ]                                                    
         m_frootc_to_litter_fire             => cnveg_carbonflux_inst%m_frootc_to_litter_fire_patch               , & ! Output: [real(r8) (:)     ]                                                    
         m_frootc_storage_to_litter_fire     => cnveg_carbonflux_inst%m_frootc_storage_to_litter_fire_patch       , & ! Output: [real(r8) (:)     ]                                                    
         m_frootc_xfer_to_litter_fire        => cnveg_carbonflux_inst%m_frootc_xfer_to_litter_fire_patch          , & ! Output: [real(r8) (:)     ]                                                    
         m_livecrootc_to_litter_fire         => cnveg_carbonflux_inst%m_livecrootc_to_litter_fire_patch           , & ! Output: [real(r8) (:)     ]                                                    
         m_livecrootc_storage_to_litter_fire => cnveg_carbonflux_inst%m_livecrootc_storage_to_litter_fire_patch   , & ! Output: [real(r8) (:)     ]                                                    
         m_livecrootc_xfer_to_litter_fire    => cnveg_carbonflux_inst%m_livecrootc_xfer_to_litter_fire_patch      , & ! Output: [real(r8) (:)     ]                                                    
         m_livecrootc_to_deadcrootc_fire     => cnveg_carbonflux_inst%m_livecrootc_to_deadcrootc_fire_patch       , & ! Output: [real(r8) (:)     ]                                                    
         m_deadcrootc_to_litter_fire         => cnveg_carbonflux_inst%m_deadcrootc_to_litter_fire_patch           , & ! Output: [real(r8) (:)     ]                                                    
         m_deadcrootc_storage_to_litter_fire => cnveg_carbonflux_inst%m_deadcrootc_storage_to_litter_fire_patch   , & ! Output: [real(r8) (:)     ]                                                    
         m_deadcrootc_xfer_to_litter_fire    => cnveg_carbonflux_inst%m_deadcrootc_xfer_to_litter_fire_patch      , & ! Output: [real(r8) (:)     ]                                                    
         m_gresp_storage_to_litter_fire      => cnveg_carbonflux_inst%m_gresp_storage_to_litter_fire_patch        , & ! Output: [real(r8) (:)     ]                                                    
         m_gresp_xfer_to_litter_fire         => cnveg_carbonflux_inst%m_gresp_xfer_to_litter_fire_patch           , & ! Output: [real(r8) (:)     ]                                                    
         m_decomp_cpools_to_fire_vr          => cnveg_carbonflux_inst%m_decomp_cpools_to_fire_vr_col              , & ! Output: [real(r8) (:,:,:) ]  (gC/m3/s) VR decomp. C fire loss
         m_c_to_litr_met_fire                => cnveg_carbonflux_inst%m_c_to_litr_met_fire_col                    , & ! Output: [real(r8) (:,:)   ]                                                  
         m_c_to_litr_cel_fire                => cnveg_carbonflux_inst%m_c_to_litr_cel_fire_col                    , & ! Output: [real(r8) (:,:)   ]                                                  
         m_c_to_litr_lig_fire                => cnveg_carbonflux_inst%m_c_to_litr_lig_fire_col                    , & ! Output: [real(r8) (:,:)   ]                                                  
         
         fire_mortality_n_to_cwdn            => cnveg_nitrogenflux_inst%fire_mortality_n_to_cwdn_col              , & ! Input:  [real(r8) (:,:)   ]  N flux fire mortality to CWD (gN/m3/s)
         m_leafn_to_fire                     => cnveg_nitrogenflux_inst%m_leafn_to_fire_patch                     , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. leafn		  
         m_leafn_storage_to_fire             => cnveg_nitrogenflux_inst%m_leafn_storage_to_fire_patch             , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. leafn_storage	  
         m_leafn_xfer_to_fire                => cnveg_nitrogenflux_inst%m_leafn_xfer_to_fire_patch                , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. leafn_xfer	       
         m_livestemn_to_fire                 => cnveg_nitrogenflux_inst%m_livestemn_to_fire_patch                 , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. livestemn	       
         m_livestemn_storage_to_fire         => cnveg_nitrogenflux_inst%m_livestemn_storage_to_fire_patch         , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. livestemn_s	       
         m_livestemn_xfer_to_fire            => cnveg_nitrogenflux_inst%m_livestemn_xfer_to_fire_patch            , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. livestemn_xfer       
         m_deadstemn_to_fire                 => cnveg_nitrogenflux_inst%m_deadstemn_to_fire_patch                 , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. deadstemn	       
         m_deadstemn_storage_to_fire         => cnveg_nitrogenflux_inst%m_deadstemn_storage_to_fire_patch         , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. deadstemn_storage    
         m_deadstemn_xfer_to_fire            => cnveg_nitrogenflux_inst%m_deadstemn_xfer_to_fire_patch            , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. deadstemn_xfer       
         m_frootn_to_fire                    => cnveg_nitrogenflux_inst%m_frootn_to_fire_patch                    , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. frootn	       
         m_frootn_storage_to_fire            => cnveg_nitrogenflux_inst%m_frootn_storage_to_fire_patch            , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. frootn_storage       
         m_frootn_xfer_to_fire               => cnveg_nitrogenflux_inst%m_frootn_xfer_to_fire_patch               , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. frootn_xfer	       
         m_livecrootn_to_fire                => cnveg_nitrogenflux_inst%m_livecrootn_to_fire_patch                , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. m_livecrootn_to_fire 
         m_livecrootn_storage_to_fire        => cnveg_nitrogenflux_inst%m_livecrootn_storage_to_fire_patch        , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. livecrootn_storage   
         m_livecrootn_xfer_to_fire           => cnveg_nitrogenflux_inst%m_livecrootn_xfer_to_fire_patch           , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. livecrootn_xfer      
         m_deadcrootn_to_fire                => cnveg_nitrogenflux_inst%m_deadcrootn_to_fire_patch                , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. deadcrootn	       
         m_deadcrootn_storage_to_fire        => cnveg_nitrogenflux_inst%m_deadcrootn_storage_to_fire_patch        , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. deadcrootn_storage   
         m_deadcrootn_xfer_to_fire           => cnveg_nitrogenflux_inst%m_deadcrootn_xfer_to_fire_patch           , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. deadcrootn_xfer      
         m_retransn_to_fire                  => cnveg_nitrogenflux_inst%m_retransn_to_fire_patch                  , & ! Input:  [real(r8) (:)     ]  (gN/m2/s) N emis. retransn             
         m_leafn_to_litter_fire              => cnveg_nitrogenflux_inst%m_leafn_to_litter_fire_patch              , & ! Output: [real(r8) (:)     ]                                                    
         m_leafn_storage_to_litter_fire      => cnveg_nitrogenflux_inst%m_leafn_storage_to_litter_fire_patch      , & ! Output: [real(r8) (:)     ]                                                    
         m_leafn_xfer_to_litter_fire         => cnveg_nitrogenflux_inst%m_leafn_xfer_to_litter_fire_patch         , & ! Output: [real(r8) (:)     ]                                                    
         m_livestemn_to_litter_fire          => cnveg_nitrogenflux_inst%m_livestemn_to_litter_fire_patch          , & ! Output: [real(r8) (:)     ]                                                    
         m_livestemn_storage_to_litter_fire  => cnveg_nitrogenflux_inst%m_livestemn_storage_to_litter_fire_patch  , & ! Output: [real(r8) (:)     ]                                                    
         m_livestemn_xfer_to_litter_fire     => cnveg_nitrogenflux_inst%m_livestemn_xfer_to_litter_fire_patch     , & ! Output: [real(r8) (:)     ]                                                    
         m_livestemn_to_deadstemn_fire       => cnveg_nitrogenflux_inst%m_livestemn_to_deadstemn_fire_patch       , & ! Output: [real(r8) (:)     ]                                                    
         m_deadstemn_to_litter_fire          => cnveg_nitrogenflux_inst%m_deadstemn_to_litter_fire_patch          , & ! Output: [real(r8) (:)     ]                                                    
         m_deadstemn_storage_to_litter_fire  => cnveg_nitrogenflux_inst%m_deadstemn_storage_to_litter_fire_patch  , & ! Output: [real(r8) (:)     ]                                                    
         m_deadstemn_xfer_to_litter_fire     => cnveg_nitrogenflux_inst%m_deadstemn_xfer_to_litter_fire_patch     , & ! Output: [real(r8) (:)     ]                                                    
         m_frootn_to_litter_fire             => cnveg_nitrogenflux_inst%m_frootn_to_litter_fire_patch             , & ! Output: [real(r8) (:)     ]                                                    
         m_frootn_storage_to_litter_fire     => cnveg_nitrogenflux_inst%m_frootn_storage_to_litter_fire_patch     , & ! Output: [real(r8) (:)     ]                                                    
         m_frootn_xfer_to_litter_fire        => cnveg_nitrogenflux_inst%m_frootn_xfer_to_litter_fire_patch        , & ! Output: [real(r8) (:)     ]                                                    
         m_livecrootn_to_litter_fire         => cnveg_nitrogenflux_inst%m_livecrootn_to_litter_fire_patch         , & ! Output: [real(r8) (:)     ]                                                    
         m_livecrootn_storage_to_litter_fire => cnveg_nitrogenflux_inst%m_livecrootn_storage_to_litter_fire_patch , & ! Output: [real(r8) (:)     ]                                                    
         m_livecrootn_xfer_to_litter_fire    => cnveg_nitrogenflux_inst%m_livecrootn_xfer_to_litter_fire_patch    , & ! Output: [real(r8) (:)     ]                                                    
         m_livecrootn_to_deadcrootn_fire     => cnveg_nitrogenflux_inst%m_livecrootn_to_deadcrootn_fire_patch     , & ! Output: [real(r8) (:)     ]                                                    
         m_deadcrootn_to_litter_fire         => cnveg_nitrogenflux_inst%m_deadcrootn_to_litter_fire_patch         , & ! Output: [real(r8) (:)     ]                                                    
         m_deadcrootn_storage_to_litter_fire => cnveg_nitrogenflux_inst%m_deadcrootn_storage_to_litter_fire_patch , & ! Output: [real(r8) (:)     ]                                                    
         m_deadcrootn_xfer_to_litter_fire    => cnveg_nitrogenflux_inst%m_deadcrootn_xfer_to_litter_fire_patch    , & ! Output: [real(r8) (:)     ]                                                    
         m_retransn_to_litter_fire           => cnveg_nitrogenflux_inst%m_retransn_to_litter_fire_patch           , & ! Output: [real(r8) (:)     ]                                                    
         m_decomp_npools_to_fire_vr          => cnveg_nitrogenflux_inst%m_decomp_npools_to_fire_vr_col            , & ! Output: [real(r8) (:,:,:) ]  VR decomp. N fire loss (gN/m3/s)
         m_n_to_litr_met_fire                => cnveg_nitrogenflux_inst%m_n_to_litr_met_fire_col                  , & ! Output: [real(r8) (:,:)   ]                                                  
         m_n_to_litr_cel_fire                => cnveg_nitrogenflux_inst%m_n_to_litr_cel_fire_col                  , & ! Output: [real(r8) (:,:)   ]                                                  
         m_n_to_litr_lig_fire                => cnveg_nitrogenflux_inst%m_n_to_litr_lig_fire_col                    & ! Output: [real(r8) (:,:)   ]                                                  
         )

     transient_landcover = run_has_transient_landcover()

     ! Get model step size
     ! calculate burned area fraction per sec
     dt = real( get_step_size(), r8 )

     dayspyr = get_days_per_year()
     !
     ! patch loop
     !
     do fp = 1,num_soilp
        p = filter_soilp(fp)
        c = patch%column(p)

        if( patch%itype(p) < nc3crop .and. cropf_col(c) < 1.0_r8)then
           ! For non-crop (bare-soil and natural vegetation)
           if (transient_landcover) then
              f = (fbac(c)-baf_crop(c))/(1.0_r8-cropf_col(c))
           else
              f = (farea_burned(c)-baf_crop(c))/(1.0_r8-cropf_col(c))
           end if
        else
           ! For crops
           if(cropf_col(c) > 0._r8)then
             f = baf_crop(c) /cropf_col(c)
           else
             f = 0._r8
           end if
        end if

        ! apply this rate to the patch state variables to get flux rates
        ! biomass burning
        ! carbon fluxes
        m = 1._r8
        if (spinup_state == 2) then
           m = 10._r8
        end if

        m_leafc_to_fire(p)               =  leafc(p)              * f * cc_leaf(patch%itype(p))
        m_leafc_storage_to_fire(p)       =  leafc_storage(p)      * f * cc_other(patch%itype(p))
        m_leafc_xfer_to_fire(p)          =  leafc_xfer(p)         * f * cc_other(patch%itype(p))
        m_livestemc_to_fire(p)           =  livestemc(p)          * f * cc_lstem(patch%itype(p))
        m_livestemc_storage_to_fire(p)   =  livestemc_storage(p)  * f * cc_other(patch%itype(p))
        m_livestemc_xfer_to_fire(p)      =  livestemc_xfer(p)     * f * cc_other(patch%itype(p))
        m_deadstemc_to_fire(p)           =  deadstemc(p)          * f * cc_dstem(patch%itype(p)) * m
        m_deadstemc_storage_to_fire(p)   =  deadstemc_storage(p)  * f * cc_other(patch%itype(p))
        m_deadstemc_xfer_to_fire(p)      =  deadstemc_xfer(p)     * f * cc_other(patch%itype(p))
        m_frootc_to_fire(p)              =  frootc(p)             * f * 0._r8
        m_frootc_storage_to_fire(p)      =  frootc_storage(p)     * f * cc_other(patch%itype(p)) 
        m_frootc_xfer_to_fire(p)         =  frootc_xfer(p)        * f * cc_other(patch%itype(p))
        m_livecrootc_to_fire(p)          =  livecrootc(p)         * f * 0._r8
        m_livecrootc_storage_to_fire(p)  =  livecrootc_storage(p) * f * cc_other(patch%itype(p)) 
        m_livecrootc_xfer_to_fire(p)     =  livecrootc_xfer(p)    * f * cc_other(patch%itype(p)) 
        m_deadcrootc_to_fire(p)          =  deadcrootc(p)         * f * 0._r8
        m_deadcrootc_storage_to_fire(p)  =  deadcrootc_storage(p) * f*  cc_other(patch%itype(p)) 
        m_deadcrootc_xfer_to_fire(p)     =  deadcrootc_xfer(p)    * f * cc_other(patch%itype(p)) 
        m_gresp_storage_to_fire(p)       =  gresp_storage(p)      * f * cc_other(patch%itype(p))
        m_gresp_xfer_to_fire(p)          =  gresp_xfer(p)         * f * cc_other(patch%itype(p))


        ! nitrogen fluxes
        m_leafn_to_fire(p)               =  leafn(p)              * f * cc_leaf(patch%itype(p))
        m_leafn_storage_to_fire(p)       =  leafn_storage(p)      * f * cc_other(patch%itype(p))
        m_leafn_xfer_to_fire(p)          =  leafn_xfer(p)         * f * cc_other(patch%itype(p))
        m_livestemn_to_fire(p)           =  livestemn(p)          * f * cc_lstem(patch%itype(p))
        m_livestemn_storage_to_fire(p)   =  livestemn_storage(p)  * f * cc_other(patch%itype(p))
        m_livestemn_xfer_to_fire(p)      =  livestemn_xfer(p)     * f * cc_other(patch%itype(p))
        m_deadstemn_to_fire(p)           =  deadstemn(p)          * f * cc_dstem(patch%itype(p)) * m
        m_deadstemn_storage_to_fire(p)   =  deadstemn_storage(p)  * f * cc_other(patch%itype(p))
        m_deadstemn_xfer_to_fire(p)      =  deadstemn_xfer(p)     * f * cc_other(patch%itype(p))
        m_frootn_to_fire(p)              =  frootn(p)             * f * 0._r8
        m_frootn_storage_to_fire(p)      =  frootn_storage(p)     * f * cc_other(patch%itype(p))
        m_frootn_xfer_to_fire(p)         =  frootn_xfer(p)        * f * cc_other(patch%itype(p))
        m_livecrootn_to_fire(p)          =  livecrootn(p)         * f * 0._r8 
        m_livecrootn_storage_to_fire(p)  =  livecrootn_storage(p) * f * cc_other(patch%itype(p)) 
        m_livecrootn_xfer_to_fire(p)     =  livecrootn_xfer(p)    * f * cc_other(patch%itype(p))
        m_deadcrootn_to_fire(p)          =  deadcrootn(p)         * f * 0._r8
        m_deadcrootn_xfer_to_fire(p)     =  deadcrootn_xfer(p)    * f * cc_other(patch%itype(p)) 
        m_deadcrootn_storage_to_fire(p)  =  deadcrootn_storage(p) * f * cc_other(patch%itype(p))
        m_retransn_to_fire(p)            =  retransn(p)           * f * cc_other(patch%itype(p))

        ! mortality due to fire
        ! carbon pools
        m_leafc_to_litter_fire(p)                   =  leafc(p) * f * &
             (1._r8 - cc_leaf(patch%itype(p))) * &
             fm_leaf(patch%itype(p))
        m_leafc_storage_to_litter_fire(p)           =  leafc_storage(p) * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))
        m_leafc_xfer_to_litter_fire(p)              =  leafc_xfer(p) * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))
        ! NOTE: It looks incorrect to use fm_droot here, but it's used to represent fraction of transport from livestem/livecroot to litter
        ! EBK Oct/06/2017 see bug 2516 http://bugs.cgd.ucar.edu/show_bug.cgi?id=2516 (stem and root live or dead assumed to have the same transport)
        m_livestemc_to_litter_fire(p)               =  livestemc(p) * f * &
             (1._r8 - cc_lstem(patch%itype(p))) * &
             fm_droot(patch%itype(p))    
        m_livestemc_storage_to_litter_fire(p)       =  livestemc_storage(p) * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))
        m_livestemc_xfer_to_litter_fire(p)          =  livestemc_xfer(p) * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p)) 
        ! NOTE: It looks incorrect to use fm_droot here, but it's used to represent the fraction of plant-tissue mortality for deadstem/deadcroot
        ! EBK Oct/06/2017 see bug 2516 http://bugs.cgd.ucar.edu/show_bug.cgi?id=2516
        m_livestemc_to_deadstemc_fire(p)            =  livestemc(p) * f * &
             (1._r8 - cc_lstem(patch%itype(p))) * &
             (fm_lstem(patch%itype(p))-fm_droot(patch%itype(p)))
        ! NOTE: It looks incorrect to use fm_droot here, but it's used to represent fraction of transport from deadstem/deadcroot to litter
        ! EBK Oct/06/2017 see bug 2516 http://bugs.cgd.ucar.edu/show_bug.cgi?id=2516 (stem and root live or dead assumed to have the same transport)
        m_deadstemc_to_litter_fire(p)               =  deadstemc(p) * f * m * &
             (1._r8 - cc_dstem(patch%itype(p))) * &
             fm_droot(patch%itype(p))    
        m_deadstemc_storage_to_litter_fire(p)       =  deadstemc_storage(p) * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))
        m_deadstemc_xfer_to_litter_fire(p)          =  deadstemc_xfer(p) * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))
        m_frootc_to_litter_fire(p)                  =  frootc(p)             * f * &
             fm_root(patch%itype(p))
        m_frootc_storage_to_litter_fire(p)          =  frootc_storage(p)     * f * &
             (1._r8- cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))
        m_frootc_xfer_to_litter_fire(p)             =  frootc_xfer(p)        * f * &
             (1._r8- cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))
        ! NOTE: It looks incorrect to use fm_droot here, but it's used to represent fraction of transport from livestem/livecroot to litter
        ! EBK Oct/06/2017 see bug 2516 http://bugs.cgd.ucar.edu/show_bug.cgi?id=2516 (stem and root live or dead assumed to have the same transport)
        m_livecrootc_to_litter_fire(p)              =  livecrootc(p)         * f * &
             fm_droot(patch%itype(p))
        m_livecrootc_storage_to_litter_fire(p)      =  livecrootc_storage(p) * f * &
             (1._r8- cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p)) 
        m_livecrootc_xfer_to_litter_fire(p)         =  livecrootc_xfer(p)    * f * &
             (1._r8- cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p)) 
        m_livecrootc_to_deadcrootc_fire(p)          =  livecrootc(p)         * f * &
             (fm_lroot(patch%itype(p))-fm_droot(patch%itype(p)))
        m_deadcrootc_to_litter_fire(p)              =  deadcrootc(p)         * f * m * &
             fm_droot(patch%itype(p))
        m_deadcrootc_storage_to_litter_fire(p)      =  deadcrootc_storage(p) * f * &
             (1._r8- cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))
        m_deadcrootc_xfer_to_litter_fire(p)         =  deadcrootc_xfer(p)    * f * &
             (1._r8- cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))      
        m_gresp_storage_to_litter_fire(p)           =  gresp_storage(p) * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))  
        m_gresp_xfer_to_litter_fire(p)              =  gresp_xfer(p) * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p)) 


        ! nitrogen pools    
        m_leafn_to_litter_fire(p)                  =  leafn(p) * f * &
             (1._r8 - cc_leaf(patch%itype(p))) * &
             fm_leaf(patch%itype(p))
        m_leafn_storage_to_litter_fire(p)          =  leafn_storage(p) * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))  
        m_leafn_xfer_to_litter_fire(p)             =  leafn_xfer(p) * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))
        ! NOTE: It looks incorrect to use fm_droot here, but it's used to represent fraction of transport from livestem/livecroot to litter
        ! EBK Oct/06/2017 see bug 2516 http://bugs.cgd.ucar.edu/show_bug.cgi?id=2516 (stem and root live or dead assumed to have the same transport)
        m_livestemn_to_litter_fire(p)              =  livestemn(p) * f * &
             (1._r8 - cc_lstem(patch%itype(p))) * &
             fm_droot(patch%itype(p))
        m_livestemn_storage_to_litter_fire(p)      =  livestemn_storage(p) * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))   
        m_livestemn_xfer_to_litter_fire(p)         =  livestemn_xfer(p) * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))
        ! NOTE: It looks incorrect to use fm_droot here, but it's used to represent the fraction of plant-tissue mortality for deadstem/deadcroot
        ! EBK Oct/06/2017 see bug 2516 http://bugs.cgd.ucar.edu/show_bug.cgi?id=2516
        m_livestemn_to_deadstemn_fire(p)           =  livestemn(p) * f * &
             (1._r8 - cc_lstem(patch%itype(p))) * &
             (fm_lstem(patch%itype(p))-fm_droot(patch%itype(p)))
        ! NOTE: It looks incorrect to use fm_droot here, but it's used to represent fraction of transport from deadstem/deadcroot to litter
        ! EBK Oct/06/2017 see bug 2516 http://bugs.cgd.ucar.edu/show_bug.cgi?id=2516 (stem and root live or dead assumed to have the same transport)
        m_deadstemn_to_litter_fire(p)              =  deadstemn(p) * f * m * &
             (1._r8 - cc_dstem(patch%itype(p))) * &
             fm_droot(patch%itype(p))    
        m_deadstemn_storage_to_litter_fire(p)      =  deadstemn_storage(p) * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))
        m_deadstemn_xfer_to_litter_fire(p)         =  deadstemn_xfer(p) * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))
        m_frootn_to_litter_fire(p)                 =  frootn(p)             * f * &
             fm_root(patch%itype(p))
        m_frootn_storage_to_litter_fire(p)         =  frootn_storage(p)     * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))
        m_frootn_xfer_to_litter_fire(p)            =  frootn_xfer(p)        * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))
        ! NOTE: It looks incorrect to use fm_droot here, but it's used to represent fraction of transport from livestem/livecroot to litter
        ! EBK Oct/06/2017 see bug 2516 http://bugs.cgd.ucar.edu/show_bug.cgi?id=2516 (stem and root live or dead assumed to have the same transport)
        m_livecrootn_to_litter_fire(p)             =  livecrootn(p)         * f * &
             fm_droot(patch%itype(p))
        m_livecrootn_storage_to_litter_fire(p)     =  livecrootn_storage(p) * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))
        m_livecrootn_xfer_to_litter_fire(p)        =  livecrootn_xfer(p)    * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p)) 
        m_livecrootn_to_deadcrootn_fire(p)         =  livecrootn(p)         * f * &
             (fm_lroot(patch%itype(p))-fm_droot(patch%itype(p)))
        m_deadcrootn_to_litter_fire(p)             =  deadcrootn(p)         * f * m * &
             fm_droot(patch%itype(p))
        m_deadcrootn_storage_to_litter_fire(p)     =  deadcrootn_storage(p) * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))
        m_deadcrootn_xfer_to_litter_fire(p)        =  deadcrootn_xfer(p)    * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p))
        m_retransn_to_litter_fire(p)               =  retransn(p)           * f * &
             (1._r8 - cc_other(patch%itype(p))) * &
             fm_other(patch%itype(p)) 

        if (use_cndv) then
           if ( woody(patch%itype(p)) == 1._r8 )then
              if ( livestemc(p)+deadstemc(p) > 0._r8 )then
                 nind(p) = nind(p)*(1._r8-1._r8*fm_droot(patch%itype(p))*f) 
              else
                 nind(p) = 0._r8
              end if
           end if
           leafcmax(p) = max(leafc(p)-m_leafc_to_fire(p)*dt, leafcmax(p))
           if (patch%itype(p) == noveg) leafcmax(p) = 0._r8
        end if

     end do  ! end of patches loop  

     ! fire-induced transfer of carbon and nitrogen pools to litter and cwd

     do j = 1,nlevdecomp
        do pi = 1,max_patch_per_col
           do fc = 1,num_soilc
              c = filter_soilc(fc)
              if (pi <=  col%npatches(c)) then
                 p = col%patchi(c) + pi - 1
                 if ( patch%active(p) ) then

                    fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                         m_deadstemc_to_litter_fire(p) * patch%wtcol(p) * stem_prof(p,j)
                    fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                         m_deadcrootc_to_litter_fire(p) * patch%wtcol(p) * croot_prof(p,j)
                    fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                         m_deadstemn_to_litter_fire(p) * patch%wtcol(p) * stem_prof(p,j)
                    fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                         m_deadcrootn_to_litter_fire(p) * patch%wtcol(p) * croot_prof(p,j)


                    fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                         m_livestemc_to_litter_fire(p) * patch%wtcol(p) * stem_prof(p,j)
                    fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                         m_livecrootc_to_litter_fire(p) * patch%wtcol(p) * croot_prof(p,j)
                    fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                         m_livestemn_to_litter_fire(p) * patch%wtcol(p) * stem_prof(p,j)
                    fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                         m_livecrootn_to_litter_fire(p) * patch%wtcol(p) * croot_prof(p,j)


                    m_c_to_litr_met_fire(c,j)=m_c_to_litr_met_fire(c,j) + &
                         ((m_leafc_to_litter_fire(p)*lf_flab(patch%itype(p)) &
                         +m_leafc_storage_to_litter_fire(p) + &
                         m_leafc_xfer_to_litter_fire(p) + &
                         m_gresp_storage_to_litter_fire(p) &
                         +m_gresp_xfer_to_litter_fire(p))*leaf_prof(p,j) + &
                         (m_frootc_to_litter_fire(p)*fr_flab(patch%itype(p)) &
                         +m_frootc_storage_to_litter_fire(p) + &
                         m_frootc_xfer_to_litter_fire(p))*froot_prof(p,j) &
                         +(m_livestemc_storage_to_litter_fire(p) + &
                         m_livestemc_xfer_to_litter_fire(p) &
                         +m_deadstemc_storage_to_litter_fire(p) + &
                         m_deadstemc_xfer_to_litter_fire(p))* stem_prof(p,j)&
                         +(m_livecrootc_storage_to_litter_fire(p) + &
                         m_livecrootc_xfer_to_litter_fire(p) &
                         +m_deadcrootc_storage_to_litter_fire(p) + &
                         m_deadcrootc_xfer_to_litter_fire(p))* croot_prof(p,j))* patch%wtcol(p)    
                    m_c_to_litr_cel_fire(c,j)=m_c_to_litr_cel_fire(c,j) + &
                         (m_leafc_to_litter_fire(p)*lf_fcel(patch%itype(p))*leaf_prof(p,j) + &
                         m_frootc_to_litter_fire(p)*fr_fcel(patch%itype(p))*froot_prof(p,j))* patch%wtcol(p) 
                    m_c_to_litr_lig_fire(c,j)=m_c_to_litr_lig_fire(c,j) + &
                         (m_leafc_to_litter_fire(p)*lf_flig(patch%itype(p))*leaf_prof(p,j) + &
                         m_frootc_to_litter_fire(p)*fr_flig(patch%itype(p))*froot_prof(p,j))* patch%wtcol(p)  

                    m_n_to_litr_met_fire(c,j)=m_n_to_litr_met_fire(c,j) + &
                         ((m_leafn_to_litter_fire(p)*lf_flab(patch%itype(p)) &
                         +m_leafn_storage_to_litter_fire(p) + &
                         m_leafn_xfer_to_litter_fire(p)+m_retransn_to_litter_fire(p)) &
                         *leaf_prof(p,j) +(m_frootn_to_litter_fire(p)*fr_flab(patch%itype(p)) &
                         +m_frootn_storage_to_litter_fire(p) + &
                         m_frootn_xfer_to_litter_fire(p))*froot_prof(p,j) &
                         +(m_livestemn_storage_to_litter_fire(p) + &
                         m_livestemn_xfer_to_litter_fire(p) &
                         +m_deadstemn_storage_to_litter_fire(p) + &
                         m_deadstemn_xfer_to_litter_fire(p))* stem_prof(p,j)&
                         +(m_livecrootn_storage_to_litter_fire(p) + &
                         m_livecrootn_xfer_to_litter_fire(p) &
                         +m_deadcrootn_storage_to_litter_fire(p) + &
                         m_deadcrootn_xfer_to_litter_fire(p))* croot_prof(p,j))* patch%wtcol(p)    
                    m_n_to_litr_cel_fire(c,j)=m_n_to_litr_cel_fire(c,j) + &
                         (m_leafn_to_litter_fire(p)*lf_fcel(patch%itype(p))*leaf_prof(p,j) + &
                         m_frootn_to_litter_fire(p)*fr_fcel(patch%itype(p))*froot_prof(p,j))* patch%wtcol(p) 
                    m_n_to_litr_lig_fire(c,j)=m_n_to_litr_lig_fire(c,j) + &
                         (m_leafn_to_litter_fire(p)*lf_flig(patch%itype(p))*leaf_prof(p,j) + &
                         m_frootn_to_litter_fire(p)*fr_flig(patch%itype(p))*froot_prof(p,j))* patch%wtcol(p) 
                 end if
              end if
           end do
        end do
     end do
     !
     ! vertically-resolved decomposing C/N fire loss   
     ! column loop
     !
     do fc = 1,num_soilc
        c = filter_soilc(fc)

        f = farea_burned(c) 

        do j = 1, nlevdecomp
           ! carbon fluxes
           do l = 1, ndecomp_pools
              if ( is_litter(l) ) then
                 m_decomp_cpools_to_fire_vr(c,j,l) = decomp_cpools_vr(c,j,l) * f * &
                      cmb_cmplt_fact(lit_fp)
              end if
              if ( is_cwd(l) ) then
                 m_decomp_cpools_to_fire_vr(c,j,l) = decomp_cpools_vr(c,j,l) * &
                      (f-baf_crop(c)) * cmb_cmplt_fact(cwd_fp)
              end if
           end do

           ! nitrogen fluxes
           do l = 1, ndecomp_pools
              if ( is_litter(l) ) then
                 m_decomp_npools_to_fire_vr(c,j,l) = decomp_npools_vr(c,j,l) * f * &
                      cmb_cmplt_fact(lit_fp)
              end if
              if ( is_cwd(l) ) then
                 m_decomp_npools_to_fire_vr(c,j,l) = decomp_npools_vr(c,j,l) * &
                      (f-baf_crop(c)) * cmb_cmplt_fact(cwd_fp)
              end if
           end do

        end do
     end do  ! end of column loop

     ! carbon loss due to deforestation fires

     if (transient_landcover) then
        call get_curr_date (kyr, kmo, kda, mcsec)
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           lfc2(c)=0._r8
           if( .not. (kmo == 1 .and. kda == 1 .and. mcsec == 0) )then
              if( trotr1_col(c)+trotr2_col(c) > 0.6_r8 .and. dtrotr_col(c) > 0._r8 .and. &
                   lfc(c) > 0._r8 .and. fbac1(c) == 0._r8) then
                 lfc2(c) = max(0._r8, min(lfc(c), (farea_burned(c)-baf_crop(c) - &
                      baf_peatf(c))/2.0*dt))/(dtrotr_col(c)*dayspyr*secspday/dt)/dt
                 lfc(c)  = lfc(c) - max(0._r8, min(lfc(c), (farea_burned(c)-baf_crop(c) - &
                      baf_peatf(c))*dt/2.0_r8))
              end if
           end if
        end do
     end if
     !
     ! Carbon loss due to peat fires
     !
     ! somc_fire is not connected to clm45 soil carbon pool, ie does not decrease
     ! soil carbon b/c clm45 soil carbon was very low in several peatland grids
     !
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        g = col%gridcell(c)
        if( grc%latdeg(g)  <  cnfire_const%borealat)then
           somc_fire(c)= totsomc(c)*baf_peatf(c)*6.0_r8/33.9_r8
        else
           somc_fire(c)= baf_peatf(c)*2.2e3_r8
        end if
     end do

     ! Fang Li has not added aerosol and trace gas emissions due to fire, yet
     ! They will be added here in proportion to the carbon emission
     ! Emission factors differ for various fire types

   end associate 

  end subroutine CNFireFluxes

  !-----------------------------------------------------------------------
  subroutine hdm_init( this, bounds, NLFilename )
   !
   ! !DESCRIPTION:
   ! Initialize data stream information for population density.
   !
   ! !USES:
   use clm_varctl       , only : inst_name
   use clm_time_manager , only : get_calendar
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use clm_nlUtilsMod   , only : find_nlgroup_name
   use ndepStreamMod    , only : clm_domain_mct
   use histFileMod      , only : hist_addfld1d
   !
   ! !ARGUMENTS:
   implicit none
   class(cnfire_base_type)       :: this
   type(bounds_type), intent(in) :: bounds  
   character(len=*),  intent(in) :: NLFilename   ! Namelist filename
   !
   ! !LOCAL VARIABLES:
   integer            :: stream_year_first_popdens   ! first year in pop. dens. stream to use
   integer            :: stream_year_last_popdens    ! last year in pop. dens. stream to use
   integer            :: model_year_align_popdens    ! align stream_year_first_hdm with 
   integer            :: nu_nml                      ! unit for namelist file
   integer            :: nml_error                   ! namelist i/o error flag
   type(mct_ggrid)    :: dom_clm                     ! domain information 
   character(len=CL)  :: stream_fldFileName_popdens  ! population density streams filename
   character(len=CL)  :: popdensmapalgo = 'bilinear' ! mapping alogrithm for population density
   character(*), parameter :: subName = "('hdmdyn_init')"
   character(*), parameter :: F00 = "('(hdmdyn_init) ',4a)"
   !-----------------------------------------------------------------------

   namelist /popd_streams/          &
        stream_year_first_popdens,  &
        stream_year_last_popdens,   &
        model_year_align_popdens,   &
        popdensmapalgo,             &
        stream_fldFileName_popdens

   ! Default values for namelist
   stream_year_first_popdens  = 1       ! first year in stream to use
   stream_year_last_popdens   = 1       ! last  year in stream to use
   model_year_align_popdens   = 1       ! align stream_year_first_popdens with this model year
   stream_fldFileName_popdens = ' '

   ! Read popd_streams namelist
   if (masterproc) then
      nu_nml = getavu()
      open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call find_nlgroup_name(nu_nml, 'popd_streams', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=popd_streams,iostat=nml_error)
         if (nml_error /= 0) then
            call endrun(msg='ERROR reading popd_streams namelist'//errMsg(sourcefile, __LINE__))
         end if
      end if
      close(nu_nml)
      call relavu( nu_nml )
   endif

   call shr_mpi_bcast(stream_year_first_popdens, mpicom)
   call shr_mpi_bcast(stream_year_last_popdens, mpicom)
   call shr_mpi_bcast(model_year_align_popdens, mpicom)
   call shr_mpi_bcast(stream_fldFileName_popdens, mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'popdens_streams settings:'
      write(iulog,*) '  stream_year_first_popdens  = ',stream_year_first_popdens  
      write(iulog,*) '  stream_year_last_popdens   = ',stream_year_last_popdens   
      write(iulog,*) '  model_year_align_popdens   = ',model_year_align_popdens   
      write(iulog,*) '  stream_fldFileName_popdens = ',stream_fldFileName_popdens
      write(iulog,*) ' '
   endif

   call clm_domain_mct (bounds, dom_clm)

   call shr_strdata_create(this%sdat_hdm,name="clmhdm",     &
        pio_subsystem=pio_subsystem,                   & 
        pio_iotype=shr_pio_getiotype(inst_name),       &
        mpicom=mpicom, compid=comp_id,                 &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,        &
        nxg=ldomain%ni, nyg=ldomain%nj,                &
        yearFirst=stream_year_first_popdens,           &
        yearLast=stream_year_last_popdens,             &
        yearAlign=model_year_align_popdens,            &
        offset=0,                                      &
        domFilePath='',                                &
        domFileName=trim(stream_fldFileName_popdens),  &
        domTvarName='time',                            &
        domXvarName='lon' ,                            &
        domYvarName='lat' ,                            &  
        domAreaName='area',                            &
        domMaskName='mask',                            &
        filePath='',                                   &
        filename=(/trim(stream_fldFileName_popdens)/) , &
        fldListFile='hdm',                             &
        fldListModel='hdm',                            &
        fillalgo='none',                               &
        mapalgo=popdensmapalgo,                        &
        calendar=get_calendar(),                       &
        tintalgo='nearest',                            &
        taxmode='extend'                           )

   if (masterproc) then
      call shr_strdata_print(this%sdat_hdm,'population density data')
   endif

   ! Add history fields
   call hist_addfld1d (fname='HDM', units='counts/km^2',      &
         avgflag='A', long_name='human population density',   &
         ptr_lnd=this%forc_hdm, default='inactive')

  end subroutine hdm_init
  
  !-----------------------------------------------------------------------
  subroutine hdm_interp( this, bounds)
  !
  ! !DESCRIPTION:
  ! Interpolate data stream information for population density.
  !
  ! !USES:
  use clm_time_manager, only : get_curr_date
  !
  ! !ARGUMENTS:
  class(cnfire_base_type)       :: this
  type(bounds_type), intent(in) :: bounds  
  !
  ! !LOCAL VARIABLES:
  integer :: g, ig
  integer :: year    ! year (0, ...) for nstep+1
  integer :: mon     ! month (1, ..., 12) for nstep+1
  integer :: day     ! day of month (1, ..., 31) for nstep+1
  integer :: sec     ! seconds into current date for nstep+1
  integer :: mcdate  ! Current model date (yyyymmdd)
  !-----------------------------------------------------------------------

   call get_curr_date(year, mon, day, sec)
   mcdate = year*10000 + mon*100 + day

   call shr_strdata_advance(this%sdat_hdm, mcdate, sec, mpicom, 'hdmdyn')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      this%forc_hdm(g) = this%sdat_hdm%avs(1)%rAttr(1,ig)
   end do
   
  end subroutine hdm_interp

  !-----------------------------------------------------------------------
  subroutine lnfm_init( this, bounds, NLFilename )
  !
  ! !DESCRIPTION:
  !
  ! Initialize data stream information for Lightning.
  !
  ! !USES:
  use clm_varctl       , only : inst_name
  use clm_time_manager , only : get_calendar
  use ncdio_pio        , only : pio_subsystem
  use shr_pio_mod      , only : shr_pio_getiotype
  use clm_nlUtilsMod   , only : find_nlgroup_name
  use ndepStreamMod    , only : clm_domain_mct
  use histFileMod      , only : hist_addfld1d
  !
  ! !ARGUMENTS:
  implicit none
  class(cnfire_base_type)       :: this
  type(bounds_type), intent(in) :: bounds  
  character(len=*),  intent(in) :: NLFilename
  !
  ! !LOCAL VARIABLES:
  integer            :: stream_year_first_lightng  ! first year in Lightning stream to use
  integer            :: stream_year_last_lightng   ! last year in Lightning stream to use
  integer            :: model_year_align_lightng   ! align stream_year_first_lnfm with 
  integer            :: nu_nml                     ! unit for namelist file
  integer            :: nml_error                  ! namelist i/o error flag
  type(mct_ggrid)    :: dom_clm                    ! domain information 
  character(len=CL)  :: stream_fldFileName_lightng ! lightning stream filename to read
  character(len=CL)  :: lightngmapalgo = 'bilinear'! Mapping alogrithm
  character(*), parameter :: subName = "('lnfmdyn_init')"
  character(*), parameter :: F00 = "('(lnfmdyn_init) ',4a)"
  !-----------------------------------------------------------------------

   namelist /light_streams/         &
        stream_year_first_lightng,  &
        stream_year_last_lightng,   &
        model_year_align_lightng,   &
        lightngmapalgo,             &
        stream_fldFileName_lightng

   ! Default values for namelist
    stream_year_first_lightng  = 1      ! first year in stream to use
    stream_year_last_lightng   = 1      ! last  year in stream to use
    model_year_align_lightng   = 1      ! align stream_year_first_lnfm with this model year
    stream_fldFileName_lightng = ' '

   ! Read light_streams namelist
   if (masterproc) then
      nu_nml = getavu()
      open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call find_nlgroup_name(nu_nml, 'light_streams', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=light_streams,iostat=nml_error)
         if (nml_error /= 0) then
            call endrun(msg='ERROR reading light_streams namelist'//errMsg(sourcefile, __LINE__))
         end if
      end if
      close(nu_nml)
      call relavu( nu_nml )
   endif

   call shr_mpi_bcast(stream_year_first_lightng, mpicom)
   call shr_mpi_bcast(stream_year_last_lightng, mpicom)
   call shr_mpi_bcast(model_year_align_lightng, mpicom)
   call shr_mpi_bcast(stream_fldFileName_lightng, mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'light_stream settings:'
      write(iulog,*) '  stream_year_first_lightng  = ',stream_year_first_lightng  
      write(iulog,*) '  stream_year_last_lightng   = ',stream_year_last_lightng   
      write(iulog,*) '  model_year_align_lightng   = ',model_year_align_lightng   
      write(iulog,*) '  stream_fldFileName_lightng = ',stream_fldFileName_lightng
      write(iulog,*) ' '
   endif

   call clm_domain_mct (bounds, dom_clm)

   call shr_strdata_create(this%sdat_lnfm,name="clmlnfm",  &
        pio_subsystem=pio_subsystem,                  & 
        pio_iotype=shr_pio_getiotype(inst_name),      &
        mpicom=mpicom, compid=comp_id,                &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,       &
        nxg=ldomain%ni, nyg=ldomain%nj,               &
        yearFirst=stream_year_first_lightng,          &
        yearLast=stream_year_last_lightng,            &
        yearAlign=model_year_align_lightng,           &
        offset=0,                                     &
        domFilePath='',                               &
        domFileName=trim(stream_fldFileName_lightng), &
        domTvarName='time',                           &
        domXvarName='lon' ,                           &
        domYvarName='lat' ,                           &  
        domAreaName='area',                           &
        domMaskName='mask',                           &
        filePath='',                                  &
        filename=(/trim(stream_fldFileName_lightng)/),&
        fldListFile='lnfm',                           &
        fldListModel='lnfm',                          &
        fillalgo='none',                              &
        mapalgo=lightngmapalgo,                       &
        calendar=get_calendar(),                      &
        taxmode='cycle'                            )

   if (masterproc) then
      call shr_strdata_print(this%sdat_lnfm,'Lightning data')
   endif

   ! Add history fields
   call hist_addfld1d (fname='LNFM', units='counts/km^2/hr',  &
         avgflag='A', long_name='Lightning frequency',        &
         ptr_lnd=this%forc_lnfm, default='inactive')

  end subroutine lnfm_init
  
  !-----------------------------------------------------------------------
  subroutine lnfm_interp(this, bounds )
  !
  ! !DESCRIPTION:
  ! Interpolate data stream information for Lightning.
  !
  ! !USES:
  use clm_time_manager, only : get_curr_date
  !
  ! !ARGUMENTS:
  class(cnfire_base_type)       :: this
  type(bounds_type), intent(in) :: bounds  
  !
  ! !LOCAL VARIABLES:
  integer :: g, ig
  integer :: year    ! year (0, ...) for nstep+1
  integer :: mon     ! month (1, ..., 12) for nstep+1
  integer :: day     ! day of month (1, ..., 31) for nstep+1
  integer :: sec     ! seconds into current date for nstep+1
  integer :: mcdate  ! Current model date (yyyymmdd)
  !-----------------------------------------------------------------------

   call get_curr_date(year, mon, day, sec)
   mcdate = year*10000 + mon*100 + day

   call shr_strdata_advance(this%sdat_lnfm, mcdate, sec, mpicom, 'lnfmdyn')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      this%forc_lnfm(g) = this%sdat_lnfm%avs(1)%rAttr(1,ig)
   end do
   
  end subroutine lnfm_interp

end module CNFireBaseMod
