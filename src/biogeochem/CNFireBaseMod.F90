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
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varctl                         , only : iulog
  use clm_varpar                         , only : nlevgrnd
  use pftconMod                          , only : noveg, pftcon
  use abortutils                         , only : endrun
  use decompMod                          , only : bounds_type
  use atm2lndType                        , only : atm2lnd_type
  use CNDVType                           , only : dgvs_type
  use CNVegStateType                     , only : cnveg_state_type
  use CNVegCarbonStateType               , only : cnveg_carbonstate_type, spinup_factor_deadwood
  use CNVegCarbonFluxType                , only : cnveg_carbonflux_type
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use CNVegNitrogenStateType             , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType              , only : cnveg_nitrogenflux_type
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use EnergyFluxType                     , only : energyflux_type
  use SaturatedExcessRunoffMod           , only : saturated_excess_runoff_type
  use WaterDiagnosticBulkType            , only : waterdiagnosticbulk_type
  use Wateratm2lndBulkType               , only : wateratm2lndbulk_type
  use WaterStateBulkType                 , only : waterstatebulk_type
  use SoilStateType                      , only : soilstate_type
  use SoilWaterRetentionCurveMod         , only : soil_water_retention_curve_type
  use GridcellType                       , only : grc
  use ColumnType                         , only : col
  use PatchType                          , only : patch
  use FireMethodType                     , only : fire_method_type
  use FireDataBaseType                   , only : fire_base_type
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: cnfire_base_type

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

      real(r8) :: cmb_cmplt_fact_litter = 0.5_r8       ! combustion completion factor for litter (unitless)
      real(r8) :: cmb_cmplt_fact_cwd    = 0.25_r8      ! combustion completion factor for CWD (unitless)
      real(r8) :: max_rh30_affecting_fuel = 90._r8     ! Value above which 30-day running relative humidity has no effect on fuel combustibility (%)
      real(r8) :: defo_fire_precip_thresh_bet = 4.0_r8     ! Max running mean daily precip (mm/d) allowing deforestation fire for broadleaf evergreen trees
      real(r8) :: defo_fire_precip_thresh_bdt = 1.8_r8     ! Max running mean daily precip (mm/d) allowing deforestation fire for broadleaf deciduous trees
      real(r8) :: borpeat_fire_soilmoist_denom = 0.3  ! Denominator of exponential in soil moisture term of equation relating that and temperature to boreal peat fire (unitless)
      real(r8) :: nonborpeat_fire_precip_denom = 1.0  ! Denominator of precipitation in equation relating that to non-boreal peat fire (unitless)
      end type

  type, public :: params_type
     real(r8) :: prh30                ! Factor related to dependence of fuel combustibility on 30-day running mean of relative humidity (unitless)
     real(r8) :: ignition_efficiency  ! Ignition efficiency of cloud-to-ground lightning (unitless)
  end type params_type

  !
  type, abstract, extends(fire_base_type) :: cnfire_base_type
    private
      ! !PRIVATE MEMBER DATA:
      ! !PUBLIC MEMBER DATA (used by extensions of the base class):
      real(r8), public, pointer :: btran2_patch   (:)   ! patch root zone soil wetness factor (0 to 1)

    contains
      !
      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public :: FireInit => CNFireInit        ! Initialization of Fire
      procedure, public :: FireReadNML                   ! Read in namelist for CNFire
      procedure, public :: CNFireReadParams              ! Read in constant parameters from the paramsfile
      procedure, public :: CNFireFluxes                  ! Calculate fire fluxes
      procedure, public :: CNFire_calc_fire_root_wetness_Li2014 ! Calculate CN-fire specific root wetness: original version
      procedure, public :: CNFire_calc_fire_root_wetness_Li2021 ! Calculate CN-fire specific root wetness: 2021 version
      ! !PRIVATE MEMBER FUNCTIONS:
      procedure, private :: InitAllocate                 ! Memory allocation of Fire
      procedure, private :: InitHistory                  ! History file assignment of fire
      !
  end type cnfire_base_type
  !-----------------------------------------------------------------------

  abstract interface
     !-----------------------------------------------------------------------
     function need_lightning_and_popdens_interface(this) result(need_lightning_and_popdens)
       !
       ! !DESCRIPTION:
       ! Returns true if need lightning and popdens, false otherwise
       !
       ! USES
       import :: cnfire_base_type
       !
       ! !ARGUMENTS:
       class(cnfire_base_type), intent(in) :: this
       logical :: need_lightning_and_popdens  ! function result
       !-----------------------------------------------------------------------
     end function need_lightning_and_popdens_interface
  end interface

  type(cnfire_const_type), public, protected :: cnfire_const          ! Fire constants shared by Li versons
  type(params_type)      , public, protected :: cnfire_params         ! Fire parameters shared by Li versions

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine CNFireInit( this, bounds, NLFilename )
    !
    ! !DESCRIPTION:
    ! Initialize CN Fire module
    ! !ARGUMENTS:
    class(cnfire_base_type) :: this
    type(bounds_type), intent(in) :: bounds
    character(len=*),  intent(in) :: NLFilename
    !-----------------------------------------------------------------------
    ! Call the base-class Initialization method
    call this%BaseFireInit( bounds, NLFilename )

    ! Allocate memory
    call this%InitAllocate( bounds )
    ! History file
    call this%InitHistory( bounds )
  end subroutine CNFireInit
  !----------------------------------------------------------------------

  subroutine InitAllocate( this, bounds )
    !
    ! Initiaze memory allocate's
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(cnfire_base_type) :: this
    type(bounds_type), intent(in) :: bounds
    !-----------------------------------------------------------------------
    integer :: begp, endp
    !------------------------------------------------------------------------
    begp = bounds%begp; endp= bounds%endp

    allocate(this%btran2_patch             (begp:endp))             ; this%btran2_patch            (:)   = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory( this, bounds )
    !
    ! Initailizae history variables
    use clm_varcon      , only : spval
    use histFileMod     , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(cnfire_base_type) :: this
    type(bounds_type), intent(in) :: bounds
    !-----------------------------------------------------------------------
    integer :: begp, endp
    !------------------------------------------------------------------------
    begp = bounds%begp; endp= bounds%endp
    this%btran2_patch(begp:endp) = spval
    call hist_addfld1d(fname='BTRAN2', units='unitless',  &
         avgflag='A', long_name='root zone soil wetness factor', &
         ptr_patch=this%btran2_patch, l2g_scale_type='veg')
  end subroutine InitHistory

  !----------------------------------------------------------------------
  subroutine CNFire_calc_fire_root_wetness_Li2014( this, bounds, &
       num_exposedvegp, filter_exposedvegp, num_noexposedvegp, filter_noexposedvegp, &
       waterstatebulk_inst, soilstate_inst, soil_water_retention_curve )
    !
    ! Calculate the root wetness term that will be used by the fire model
    !
    class(cnfire_base_type) :: this
    type(bounds_type)      , intent(in)   :: bounds                         !bounds
    integer                , intent(in)   :: num_exposedvegp                !number of filters
    integer                , intent(in)   :: filter_exposedvegp(:)          !filter array
    integer                , intent(in)   :: num_noexposedvegp       ! number of points in filter_noexposedvegp
    integer                , intent(in)   :: filter_noexposedvegp(:) ! patch filter where frac_veg_nosno is 0 
    type(waterstatebulk_type), intent(in) :: waterstatebulk_inst
    type(soilstate_type)   , intent(in)   :: soilstate_inst
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    ! !LOCAL VARIABLES:
    real(r8) :: smp_node, s_node  !temporary variables
    real(r8) :: smp_node_lf       !temporary variable
    integer :: p, fp, j, c, l      !indices
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(filter_exposedvegp) >= (/num_exposedvegp/)), sourcefile, __LINE__)

    associate(                                                &
         smpso         => pftcon%smpso                      , & ! Input:  soil water potential at full stomatal opening (mm)
         smpsc         => pftcon%smpsc                      , & ! Input:  soil water potential at full stomatal closure (mm)
         watsat        => soilstate_inst%watsat_col         , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation
         btran2        => this%btran2_patch                 , & ! Output: [real(r8) (:)   ]  integrated soil water stress square
         rootfr        => soilstate_inst%rootfr_patch       , & ! Input:  [real(r8) (:,:) ]  fraction of roots in each soil layer
         h2osoi_vol    => waterstatebulk_inst%h2osoi_vol_col  & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3] (porosity)   (constant)
         )

    do fp = 1, num_noexposedvegp
       p = filter_noexposedvegp(fp)
       ! Set for the sake of history diagnostics. The "normal" btran is set to 0 over
       ! this filter, so we do the same for btran2.
       btran2(p) = 0._r8
    end do

    do fp = 1, num_exposedvegp
       p = filter_exposedvegp(fp)
       btran2(p) = 0._r8
    end do
    do j = 1,nlevgrnd
       do fp = 1, num_exposedvegp
          p = filter_exposedvegp(fp)
          c = patch%column(p)
          l = patch%landunit(p)
          s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)

          call soil_water_retention_curve%soil_suction(c, j, s_node, soilstate_inst, smp_node_lf)

          smp_node_lf = max(smpsc(patch%itype(p)), smp_node_lf)
          btran2(p)   = btran2(p) +rootfr(p,j)*max(0._r8,min((smp_node_lf - smpsc(patch%itype(p))) / &
               (smpso(patch%itype(p)) - smpsc(patch%itype(p))), 1._r8))
       end do
    end do

    do fp = 1, num_exposedvegp
       p = filter_exposedvegp(fp)
       if (btran2(p) > 1._r8) then
          btran2(p) = 1._r8
       end if
    end do

    end associate

  end subroutine CNFire_calc_fire_root_wetness_Li2014

  !----------------------------------------------------------------------
  subroutine CNFire_calc_fire_root_wetness_Li2021( this, bounds, &
       num_exposedvegp, filter_exposedvegp, num_noexposedvegp, filter_noexposedvegp, &
       waterstatebulk_inst, soilstate_inst, soil_water_retention_curve )
    !
    ! Calculate the root wetness term that will be used by the fire model
    !
    use pftconMod                 , only : pftcon
    use PatchType                 , only : patch
    class(cnfire_base_type) :: this
    type(bounds_type)      , intent(in)   :: bounds                         !bounds
    integer                , intent(in)   :: num_exposedvegp                !number of filters
    integer                , intent(in)   :: filter_exposedvegp(:)          !filter array
    integer                , intent(in)   :: num_noexposedvegp       ! number of points in filter_noexposedvegp
    integer                , intent(in)   :: filter_noexposedvegp(:) ! patch filter where frac_veg_nosno is 0 
    type(waterstatebulk_type), intent(in) :: waterstatebulk_inst
    type(soilstate_type)   , intent(in)   :: soilstate_inst
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    ! !LOCAL VARIABLES:
    real(r8) :: s_node  !temporary variables
    integer :: p, fp, j, c         !indices
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(filter_exposedvegp) >= (/num_exposedvegp/)), sourcefile, __LINE__)

    associate(                                                &
         watsat        => soilstate_inst%watsat_col         , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation
         btran2        => this%btran2_patch                 , & ! Output: [real(r8) (:)   ]  integrated soil water stress square
         rootfr        => soilstate_inst%rootfr_patch       , & ! Input:  [real(r8) (:,:) ]  fraction of roots in each soil layer
         h2osoi_vol    => waterstatebulk_inst%h2osoi_vol_col  & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3] (porosity)   (constant)
         )

    do fp = 1, num_noexposedvegp
       p = filter_noexposedvegp(fp)
       ! Set for the sake of history diagnostics. The "normal" btran is set to 0 over
       ! this filter, so we do the same for btran2.
       btran2(p) = 0._r8
    end do

    do fp = 1, num_exposedvegp
       p = filter_exposedvegp(fp)
       btran2(p)   = 0._r8
    end do
    do j = 1,nlevgrnd
       do fp = 1, num_exposedvegp
          p = filter_exposedvegp(fp)
          c = patch%column(p)
          s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)

          btran2(p)   = btran2(p) + rootfr(p,j)*s_node
       end do
    end do

    do fp = 1, num_exposedvegp
       p = filter_exposedvegp(fp)
       if (btran2(p) > 1._r8) then
          btran2(p) = 1._r8
       end if
    end do

    end associate

  end subroutine CNFire_calc_fire_root_wetness_Li2021
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  subroutine FireReadNML( this, NLFilename )
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

    character(len=*), parameter :: subname = 'FireReadNML'
    character(len=*), parameter :: nmlname = 'lifire_inparm'
    !-----------------------------------------------------------------------
    real(r8) :: cli_scale, boreal_peatfire_c, pot_hmn_ign_counts_alpha
    real(r8) :: non_boreal_peatfire_c, cropfire_a1
    real(r8) :: rh_low, rh_hgh, bt_min, bt_max, occur_hi_gdp_tree
    real(r8) :: lfuel, ufuel, cmb_cmplt_fact_litter, cmb_cmplt_fact_cwd
    real(r8) :: max_rh30_affecting_fuel
    real(r8) :: defo_fire_precip_thresh_bet, defo_fire_precip_thresh_bdt
    real(r8) :: borpeat_fire_soilmoist_denom, nonborpeat_fire_precip_denom

    namelist /lifire_inparm/ cli_scale, boreal_peatfire_c, pot_hmn_ign_counts_alpha, &
                             non_boreal_peatfire_c, cropfire_a1,                &
                             rh_low, rh_hgh, bt_min, bt_max, occur_hi_gdp_tree, &
                             lfuel, ufuel, cmb_cmplt_fact_litter, cmb_cmplt_fact_cwd, &
                             max_rh30_affecting_fuel, &
                             defo_fire_precip_thresh_bet, defo_fire_precip_thresh_bdt, &
                             borpeat_fire_soilmoist_denom, nonborpeat_fire_precip_denom

    if ( this%need_lightning_and_popdens() ) then
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
       cmb_cmplt_fact_litter     = cnfire_const%cmb_cmplt_fact_litter
       cmb_cmplt_fact_cwd        = cnfire_const%cmb_cmplt_fact_cwd
       max_rh30_affecting_fuel   = cnfire_const%max_rh30_affecting_fuel
       defo_fire_precip_thresh_bet = cnfire_const%defo_fire_precip_thresh_bet
       defo_fire_precip_thresh_bdt = cnfire_const%defo_fire_precip_thresh_bdt
       borpeat_fire_soilmoist_denom = cnfire_const%borpeat_fire_soilmoist_denom
       nonborpeat_fire_precip_denom = cnfire_const%nonborpeat_fire_precip_denom
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
       call shr_mpi_bcast (cmb_cmplt_fact_litter   , mpicom)
       call shr_mpi_bcast (cmb_cmplt_fact_cwd      , mpicom)
       call shr_mpi_bcast (max_rh30_affecting_fuel , mpicom)
       call shr_mpi_bcast (defo_fire_precip_thresh_bet, mpicom)
       call shr_mpi_bcast (defo_fire_precip_thresh_bdt, mpicom)
       call shr_mpi_bcast (borpeat_fire_soilmoist_denom, mpicom)
       call shr_mpi_bcast (nonborpeat_fire_precip_denom, mpicom)

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
       cnfire_const%cmb_cmplt_fact_litter     = cmb_cmplt_fact_litter
       cnfire_const%cmb_cmplt_fact_cwd        = cmb_cmplt_fact_cwd
       cnfire_const%max_rh30_affecting_fuel   = max_rh30_affecting_fuel
       cnfire_const%defo_fire_precip_thresh_bet = defo_fire_precip_thresh_bet
       cnfire_const%defo_fire_precip_thresh_bdt = defo_fire_precip_thresh_bdt
       cnfire_const%borpeat_fire_soilmoist_denom = borpeat_fire_soilmoist_denom
       cnfire_const%nonborpeat_fire_precip_denom = nonborpeat_fire_precip_denom

       if (masterproc) then
          write(iulog,*) ' '
          write(iulog,*) nmlname//' settings:'
          write(iulog,nml=lifire_inparm)
          write(iulog,*) ' '
       end if
    end if

  end subroutine FireReadNML

  !-----------------------------------------------------------------------
  subroutine CNFireFluxes (this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
      num_actfirec, filter_actfirec, num_actfirep, filter_actfirep,                        &
      dgvs_inst, cnveg_state_inst,                                                                      &
      cnveg_carbonstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
      soilbiogeochem_carbonflux_inst,                                       &
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
   use clm_time_manager     , only: get_step_size_real,get_curr_days_per_year,get_curr_date
   use clm_varctl           , only: use_cndv
   use SoilBiogeochemDecompCascadeConType , only : use_soil_matrixcn
   use CNSharedParamsMod    , only: use_matrixcn
   use clm_varcon           , only: secspday
   use pftconMod            , only: nc3crop
   use dynSubgridControlMod , only: run_has_transient_landcover
   use clm_varpar           , only: nlevdecomp_full, ndecomp_pools, nlevdecomp, i_litr_max, i_met_lit
   use clm_varpar           , only: ileaf,ileaf_st,ileaf_xf,ifroot,ifroot_st,ifroot_xf,&
                                    ilivestem,ilivestem_st,ilivestem_xf,&
                                    ideadstem,ideadstem_st,ideadstem_xf,&
                                    ilivecroot,ilivecroot_st,ilivecroot_xf,&
                                    ideadcroot,ideadcroot_st,ideadcroot_xf,iretransn,ioutc,ioutn
   use CNVegMatrixMod       , only: matrix_update_fic, matrix_update_fin
   !
   ! !ARGUMENTS:
   class(cnfire_base_type)                              :: this
   type(bounds_type)                    , intent(in)    :: bounds  
   integer                              , intent(in)    :: num_soilc                               ! number of soil columns in filter
   integer                              , intent(in)    :: filter_soilc(:)                         ! filter for soil columns
   integer                              , intent(in)    :: num_soilp                               ! number of soil patches in filter
   integer                              , intent(in)    :: filter_soilp(:)                         ! filter for soil patches
   integer                              , intent(out)   :: num_actfirep                            ! number of active patches on fire in filter
   integer                              , intent(out)   :: filter_actfirep(:)                      ! filter for soil patches
   integer                              , intent(out)   :: num_actfirec                            ! number of active columns on fire in filter
   integer                              , intent(out)   :: filter_actfirec(:)                      ! filter for soil columns
   type(dgvs_type)                      , intent(inout) :: dgvs_inst
   type(cnveg_state_type)               , intent(inout) :: cnveg_state_inst
   type(soilbiogeochem_carbonflux_type) , intent(inout) :: soilbiogeochem_carbonflux_inst          ! only for matrix_decomp_fire_k: (gC/m3/step) VR deomp. C fire loss in matrix representation
   type(cnveg_carbonstate_type)         , intent(inout) :: cnveg_carbonstate_inst
   type(cnveg_carbonflux_type)          , intent(inout) :: cnveg_carbonflux_inst
   type(cnveg_nitrogenstate_type)       , intent(in)    :: cnveg_nitrogenstate_inst
   type(cnveg_nitrogenflux_type)        , intent(inout) :: cnveg_nitrogenflux_inst
   real(r8)                             , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
   real(r8)                             , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
   real(r8)                             , intent(in)    :: croot_prof_patch(bounds%begp:,1:)
   real(r8)                             , intent(in)    :: stem_prof_patch(bounds%begp:,1:)
   real(r8)                             , intent(in)    :: totsomc_col(bounds%begc:)                ! (gC/m2) total soil organic matter C
   real(r8)                             , intent(in)    :: decomp_cpools_vr_col(bounds%begc:,1:,1:) ! (gC/m3)  VR decomp. (litter, cwd, soil)
   real(r8)                             , intent(in)    :: decomp_npools_vr_col(bounds%begc:,1:,1:) ! (gC/m3)  VR decomp. (litter, cwd, soil)
   real(r8)                             , intent(out)   :: somc_fire_col(bounds%begc:)              ! (gC/m2/s) fire C emissions due to peat burning
   !
   ! !LOCAL VARIABLES:
   integer :: i,g,c,p,j,l,kyr, kmo, kda, mcsec  ! indices
   integer :: fp,fc                ! filter indices
   real(r8):: f                    ! rate for fire effects (1/s)
   real(r8):: m                    ! acceleration factor for fuel carbon
   real(r8):: dt                   ! time step variable (s)
   real(r8):: dayspyr              ! days per year
   logical :: transient_landcover  ! whether this run has any prescribed transient landcover
   !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(leaf_prof_patch)      == (/bounds%endp,nlevdecomp_full/))               , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(froot_prof_patch)     == (/bounds%endp,nlevdecomp_full/))               , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(croot_prof_patch)     == (/bounds%endp,nlevdecomp_full/))               , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(stem_prof_patch)      == (/bounds%endp,nlevdecomp_full/))               , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(totsomc_col)          == (/bounds%endc/))                               , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(decomp_cpools_vr_col) == (/bounds%endc,nlevdecomp_full,ndecomp_pools/)) , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(decomp_npools_vr_col) == (/bounds%endc,nlevdecomp_full,ndecomp_pools/)) , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(somc_fire_col)        == (/bounds%endc/))                               , sourcefile, __LINE__)

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
         lf_f                                => pftcon%lf_f                                                       , & ! Input:
         fr_f                                => pftcon%fr_f                                                       , & ! Input:

         cmb_cmplt_fact_litter               => cnfire_const%cmb_cmplt_fact_litter                                , & ! Input:  [real(r8) (:)     ]  Combustion completion factor for litter (unitless)
         cmb_cmplt_fact_cwd                  => cnfire_const%cmb_cmplt_fact_cwd                                   , & ! Input:  [real(r8) (:)     ]  Combustion completion factor for CWD (unitless)
         
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
         m_c_to_litr_fire                    => cnveg_carbonflux_inst%m_c_to_litr_fire_col                        , & ! Output: [real(r8) (:,:,:) ]
         
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
         m_n_to_litr_fire                    => cnveg_nitrogenflux_inst%m_n_to_litr_fire_col                      , & ! Output: [real(r8) (:,:)   ]                                                  
         ileaf_to_iout_fic                   => cnveg_carbonflux_inst%ileaf_to_iout_fi                            , & ! Input: [integer (:)] Index of fire related C transfer from leaf pool to outside of vegetation pools
         ileafst_to_iout_fic                 => cnveg_carbonflux_inst%ileafst_to_iout_fi                          , & ! Input: [integer (:)] Index of fire related C transfer from leaf storage pool to outside of vegetation pools
         ileafxf_to_iout_fic                 => cnveg_carbonflux_inst%ileafxf_to_iout_fi                          , & ! Input: [integer (:)] Index of fire related C transfer from leaf transfer pool to outside of vegetation pools
         ifroot_to_iout_fic                  => cnveg_carbonflux_inst%ifroot_to_iout_fi                           , & ! Input: [integer (:)] Index of fire related C transfer from fine root pool to outside of vegetation pools
         ifrootst_to_iout_fic                => cnveg_carbonflux_inst%ifrootst_to_iout_fi                         , & ! Input: [integer (:)] Index of fire related C transfer from fine root storage pool to outside of vegetation pools
         ifrootxf_to_iout_fic                => cnveg_carbonflux_inst%ifrootxf_to_iout_fi                         , & ! Input: [integer (:)] Index of fire related C transfer from fine root transfer pool to outside of vegetation pools
         ilivestem_to_iout_fic               => cnveg_carbonflux_inst%ilivestem_to_iout_fi                        , & ! Input: [integer (:)] Index of fire related C transfer from live stem pool to outside of vegetation pools
         ilivestemst_to_iout_fic             => cnveg_carbonflux_inst%ilivestemst_to_iout_fi                      , & ! Input: [integer (:)] Index of fire related C transfer from live stem storage pool to outside of vegetation pools
         ilivestemxf_to_iout_fic             => cnveg_carbonflux_inst%ilivestemxf_to_iout_fi                      , & ! Input: [integer (:)] Index of fire related C transfer from live stem transfer pool to outside of vegetation pools
         ideadstem_to_iout_fic               => cnveg_carbonflux_inst%ideadstem_to_iout_fi                        , & ! Input: [integer (:)] Index of fire related C transfer from dead stem pool to outside of vegetation pools
         ideadstemst_to_iout_fic             => cnveg_carbonflux_inst%ideadstemst_to_iout_fi                      , & ! Input: [integer (:)] Index of fire related C transfer from dead stem storage pool to outside of vegetation pools
         ideadstemxf_to_iout_fic             => cnveg_carbonflux_inst%ideadstemxf_to_iout_fi                      , & ! Input: [integer (:)] Index of fire related C transfer from dead stem transfer pool to outside of vegetation pools
         ilivecroot_to_iout_fic              => cnveg_carbonflux_inst%ilivecroot_to_iout_fi                       , & ! Input: [integer (:)] Index of fire related C transfer from live coarse root pool to outside of vegetation pools
         ilivecrootst_to_iout_fic            => cnveg_carbonflux_inst%ilivecrootst_to_iout_fi                     , & ! Input: [integer (:)] Index of fire related C transfer from live coarse root storage pool to outside of vegetation pools
         ilivecrootxf_to_iout_fic            => cnveg_carbonflux_inst%ilivecrootxf_to_iout_fi                     , & ! Input: [integer (:)] Index of fire related C transfer from live coarse root transfer pool to outside of vegetation pools
         ideadcroot_to_iout_fic              => cnveg_carbonflux_inst%ideadcroot_to_iout_fi                       , & ! Input: [integer (:)] Index of fire related C transfer from dead coarse root pool to outside of vegetation pools
         ideadcrootst_to_iout_fic            => cnveg_carbonflux_inst%ideadcrootst_to_iout_fi                     , & ! Input: [integer (:)] Index of fire related C transfer from dead coarse root storage pool to outside of vegetation pools
         ideadcrootxf_to_iout_fic            => cnveg_carbonflux_inst%ideadcrootxf_to_iout_fi                     , & ! Input: [integer (:)] Index of fire related C transfer from dead coarse root transfer pool to outside of vegetation pools
         ilivestem_to_ideadstem_fic          => cnveg_carbonflux_inst%ilivestem_to_ideadstem_fi                   , & ! Input: [integer (:)] Index of fire related C transfer from live stem pool to dead stem pool
         ilivecroot_to_ideadcroot_fic        => cnveg_carbonflux_inst%ilivecroot_to_ideadcroot_fi                 , & ! Input: [integer (:)] Index of fire related C transfer from live coarse root pool to dead coarse root pool
         ileaf_to_iout_fin                   => cnveg_nitrogenflux_inst%ileaf_to_iout_fi                          , & ! Input: [integer (:)] Index of fire related N transfer from leaf pool to outside of vegetation pools 
         ileafst_to_iout_fin                 => cnveg_nitrogenflux_inst%ileafst_to_iout_fi                        , & ! Input: [integer (:)] Index of fire related N transfer from leaf storage pool to outside of vegetation pools
         ileafxf_to_iout_fin                 => cnveg_nitrogenflux_inst%ileafxf_to_iout_fi                        , & ! Input: [integer (:)] Index of fire related N transfer from leaf transfer pool to outside of vegetation pools
         ifroot_to_iout_fin                  => cnveg_nitrogenflux_inst%ifroot_to_iout_fi                         , & ! Input: [integer (:)] Index of fire related N transfer from fine root pool to outside of vegetation pools
         ifrootst_to_iout_fin                => cnveg_nitrogenflux_inst%ifrootst_to_iout_fi                       , & ! Input: [integer (:)] Index of fire related N transfer from fine root storage pool to outside of vegetation pools
         ifrootxf_to_iout_fin                => cnveg_nitrogenflux_inst%ifrootxf_to_iout_fi                       , & ! Input: [integer (:)] Index of fire related N transfer from fine transfer pool to outside of vegetation pools
         ilivestem_to_iout_fin               => cnveg_nitrogenflux_inst%ilivestem_to_iout_fi                      , & ! Input: [integer (:)] Index of fire related N transfer from live stem pool to outside of vegetation pools
         ilivestemst_to_iout_fin             => cnveg_nitrogenflux_inst%ilivestemst_to_iout_fi                    , & ! Input: [integer (:)] Index of fire related N transfer from live stem storage pool to outside of vegetation pools
         ilivestemxf_to_iout_fin             => cnveg_nitrogenflux_inst%ilivestemxf_to_iout_fi                    , & ! Input: [integer (:)] Index of fire related N transfer from live stem transfer pool to outside of vegetation pools
         ideadstem_to_iout_fin               => cnveg_nitrogenflux_inst%ideadstem_to_iout_fi                      , & ! Input: [integer (:)] Index of fire related N transfer from dead stem pool to outside of vegetation pools
         ideadstemst_to_iout_fin             => cnveg_nitrogenflux_inst%ideadstemst_to_iout_fi                    , & ! Input: [integer (:)] Index of fire related N transfer from dead stem storage pool to outside of vegetation pools
         ideadstemxf_to_iout_fin             => cnveg_nitrogenflux_inst%ideadstemxf_to_iout_fi                    , & ! Input: [integer (:)] Index of fire related N transfer from dead stem transfer pool to outside of vegetation pools
         ilivecroot_to_iout_fin              => cnveg_nitrogenflux_inst%ilivecroot_to_iout_fi                     , & ! Input: [integer (:)] Index of fire related N transfer from live coarse root pool to outside of vegetation pools
         ilivecrootst_to_iout_fin            => cnveg_nitrogenflux_inst%ilivecrootst_to_iout_fi                   , & ! Input: [integer (:)] Index of fire related N transfer from live coarse root storage pool to outside of vegetation pools
         ilivecrootxf_to_iout_fin            => cnveg_nitrogenflux_inst%ilivecrootxf_to_iout_fi                   , & ! Input: [integer (:)] Index of fire related N transfer from live coarse root transfer pool to outside of vegetation pools
         ideadcroot_to_iout_fin              => cnveg_nitrogenflux_inst%ideadcroot_to_iout_fi                     , & ! Input: [integer (:)] Index of fire related N transfer from dead coarse root pool to outside of vegetation pools
         ideadcrootst_to_iout_fin            => cnveg_nitrogenflux_inst%ideadcrootst_to_iout_fi                   , & ! Input: [integer (:)] Index of fire related N transfer from dead coarse root storage pool to outside of vegetation pools
         ideadcrootxf_to_iout_fin            => cnveg_nitrogenflux_inst%ideadcrootxf_to_iout_fi                   , & ! Input: [integer (:)] Index of fire related N transfer from dead coarse root transfer pool to outside of vegetation pools
         ilivestem_to_ideadstem_fin          => cnveg_nitrogenflux_inst%ilivestem_to_ideadstem_fi                 , & ! Input: [integer (:)] Index of fire related N transfer from live stem to dead stem pool
         ilivecroot_to_ideadcroot_fin        => cnveg_nitrogenflux_inst%ilivecroot_to_ideadcroot_fi               , & ! Input: [integer (:)] Index of fire related N transfer from live coarse root pool to dead coarse root pool
         iretransn_to_iout_fin               => cnveg_nitrogenflux_inst%iretransn_to_iout_fi                        & ! Input: [integer (:)] Index of fire related N transfer from retranslocated N pool to outside of vegetation pools
         )

     transient_landcover = run_has_transient_landcover()

     ! Get model step size
     ! calculate burned area fraction per sec
     dt = get_step_size_real()

     dayspyr = get_curr_days_per_year()
     !
     ! patch loop
     !
     num_actfirep = 0
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
        m = spinup_factor_deadwood

        if(f /= 0)then
           num_actfirep = num_actfirep + 1
           filter_actfirep(num_actfirep) = p
        end if
        m_gresp_storage_to_fire(p)       =  gresp_storage(p)      * f * cc_other(patch%itype(p))
        m_gresp_xfer_to_fire(p)          =  gresp_xfer(p)         * f * cc_other(patch%itype(p))
        if ( .not. use_matrixcn )then
           ! NOTE: The non matrix version of this is in CNCStateUpdate3::CStateUpdate3 EBK (11/26/2019)
           !                                        and CNNStateUpdate3::NStateUpdate3
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

        else
           m_leafc_to_fire(p)               =  leafc(p)              * matrix_update_fic(p,ileaf_to_iout_fic        ,f * cc_leaf(patch%itype(p))   ,dt,cnveg_carbonflux_inst,.True.,.True.)
           m_leafc_storage_to_fire(p)       =  leafc_storage(p)      * matrix_update_fic(p,ileafst_to_iout_fic      ,f * cc_other(patch%itype(p))  ,dt,cnveg_carbonflux_inst,.True.,.True.)
           m_leafc_xfer_to_fire(p)          =  leafc_xfer(p)         * matrix_update_fic(p,ileafxf_to_iout_fic      ,f * cc_other(patch%itype(p))  ,dt,cnveg_carbonflux_inst,.True.,.True.)
           m_livestemc_to_fire(p)           =  livestemc(p)          * matrix_update_fic(p,ilivestem_to_iout_fic    ,f * cc_lstem(patch%itype(p))  ,dt,cnveg_carbonflux_inst,.True.,.True.)
           m_livestemc_storage_to_fire(p)   =  livestemc_storage(p)  * matrix_update_fic(p,ilivestemst_to_iout_fic  ,f * cc_other(patch%itype(p))  ,dt,cnveg_carbonflux_inst,.True.,.True.)
           m_livestemc_xfer_to_fire(p)      =  livestemc_xfer(p)     * matrix_update_fic(p,ilivestemxf_to_iout_fic  ,f * cc_other(patch%itype(p))  ,dt,cnveg_carbonflux_inst,.True.,.True.)
           m_deadstemc_to_fire(p)           =  deadstemc(p)          * matrix_update_fic(p,ideadstem_to_iout_fic    ,f * cc_dstem(patch%itype(p))*m,dt,cnveg_carbonflux_inst,.True.,.True.)
           m_deadstemc_storage_to_fire(p)   =  deadstemc_storage(p)  * matrix_update_fic(p,ideadstemst_to_iout_fic  ,f * cc_other(patch%itype(p))  ,dt,cnveg_carbonflux_inst,.True.,.True.)
           m_deadstemc_xfer_to_fire(p)      =  deadstemc_xfer(p)     * matrix_update_fic(p,ideadstemxf_to_iout_fic  ,f * cc_other(patch%itype(p))  ,dt,cnveg_carbonflux_inst,.True.,.True.)
           m_frootc_to_fire(p)              =  frootc(p)             * matrix_update_fic(p,ifroot_to_iout_fic       ,f * 0._r8                     ,dt,cnveg_carbonflux_inst,.True.,.True.)
           m_frootc_storage_to_fire(p)      =  frootc_storage(p)     * matrix_update_fic(p,ifrootst_to_iout_fic     ,f * cc_other(patch%itype(p))  ,dt,cnveg_carbonflux_inst,.True.,.True.)
           m_frootc_xfer_to_fire(p)         =  frootc_xfer(p)        * matrix_update_fic(p,ifrootxf_to_iout_fic     ,f * cc_other(patch%itype(p))  ,dt,cnveg_carbonflux_inst,.True.,.True.)
           m_livecrootc_to_fire(p)          =  livecrootc(p)         * matrix_update_fic(p,ilivecroot_to_iout_fic   ,f * 0._r8                     ,dt,cnveg_carbonflux_inst,.True.,.True.)
           m_livecrootc_storage_to_fire(p)  =  livecrootc_storage(p) * matrix_update_fic(p,ilivecrootst_to_iout_fic ,f * cc_other(patch%itype(p))  ,dt,cnveg_carbonflux_inst,.True.,.True.)
           m_livecrootc_xfer_to_fire(p)     =  livecrootc_xfer(p)    * matrix_update_fic(p,ilivecrootxf_to_iout_fic ,f * cc_other(patch%itype(p))  ,dt,cnveg_carbonflux_inst,.True.,.True.)
           m_deadcrootc_to_fire(p)          =  deadcrootc(p)         * matrix_update_fic(p,ideadcroot_to_iout_fic   ,f * 0._r8                     ,dt,cnveg_carbonflux_inst,.True.,.True.)
           m_deadcrootc_storage_to_fire(p)  =  deadcrootc_storage(p) * matrix_update_fic(p,ideadcrootst_to_iout_fic ,f * cc_other(patch%itype(p))  ,dt,cnveg_carbonflux_inst,.True.,.True.)
           m_deadcrootc_xfer_to_fire(p)     =  deadcrootc_xfer(p)    * matrix_update_fic(p,ideadcrootxf_to_iout_fic ,f * cc_other(patch%itype(p))  ,dt,cnveg_carbonflux_inst,.True.,.True.)

           m_leafn_to_fire(p)               =  leafn(p)              * matrix_update_fin(p,ileaf_to_iout_fin        ,f * cc_leaf(patch%itype(p))   ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_leafn_storage_to_fire(p)       =  leafn_storage(p)      * matrix_update_fin(p,ileafst_to_iout_fin      ,f * cc_other(patch%itype(p))  ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_leafn_xfer_to_fire(p)          =  leafn_xfer(p)         * matrix_update_fin(p,ileafxf_to_iout_fin      ,f * cc_other(patch%itype(p))  ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_livestemn_to_fire(p)           =  livestemn(p)          * matrix_update_fin(p,ilivestem_to_iout_fin    ,f * cc_lstem(patch%itype(p))  ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_livestemn_storage_to_fire(p)   =  livestemn_storage(p)  * matrix_update_fin(p,ilivestemst_to_iout_fin  ,f * cc_other(patch%itype(p))  ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_livestemn_xfer_to_fire(p)      =  livestemn_xfer(p)     * matrix_update_fin(p,ilivestemxf_to_iout_fin  ,f * cc_other(patch%itype(p))  ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_deadstemn_to_fire(p)           =  deadstemn(p)          * matrix_update_fin(p,ideadstem_to_iout_fin    ,f * cc_dstem(patch%itype(p))*m,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_deadstemn_storage_to_fire(p)   =  deadstemn_storage(p)  * matrix_update_fin(p,ideadstemst_to_iout_fin  ,f * cc_other(patch%itype(p))  ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_deadstemn_xfer_to_fire(p)      =  deadstemn_xfer(p)     * matrix_update_fin(p,ideadstemxf_to_iout_fin  ,f * cc_other(patch%itype(p))  ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_frootn_to_fire(p)              =  frootn(p)             * matrix_update_fin(p,ifroot_to_iout_fin       ,f * 0._r8                     ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_frootn_storage_to_fire(p)      =  frootn_storage(p)     * matrix_update_fin(p,ifrootst_to_iout_fin     ,f * cc_other(patch%itype(p))  ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_frootn_xfer_to_fire(p)         =  frootn_xfer(p)        * matrix_update_fin(p,ifrootxf_to_iout_fin     ,f * cc_other(patch%itype(p))  ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_livecrootn_to_fire(p)          =  livecrootn(p)         * matrix_update_fin(p,ilivecroot_to_iout_fin   ,f * 0._r8                     ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_livecrootn_storage_to_fire(p)  =  livecrootn_storage(p) * matrix_update_fin(p,ilivecrootst_to_iout_fin ,f * cc_other(patch%itype(p))  ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_livecrootn_xfer_to_fire(p)     =  livecrootn_xfer(p)    * matrix_update_fin(p,ilivecrootxf_to_iout_fin ,f * cc_other(patch%itype(p))  ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_deadcrootn_to_fire(p)          =  deadcrootn(p)         * matrix_update_fin(p,ideadcroot_to_iout_fin   ,f * 0._r8                     ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_deadcrootn_storage_to_fire(p)  =  deadcrootn_storage(p) * matrix_update_fin(p,ideadcrootst_to_iout_fin ,f * cc_other(patch%itype(p))  ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_deadcrootn_xfer_to_fire(p)     =  deadcrootn_xfer(p)    * matrix_update_fin(p,ideadcrootxf_to_iout_fin ,f * cc_other(patch%itype(p))  ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_retransn_to_fire(p)            =  retransn(p)           * matrix_update_fin(p,iretransn_to_iout_fin    ,f * 0._r8                     ,dt,cnveg_nitrogenflux_inst,.True.,.True.)
        end if
        ! mortality due to fire
        ! carbon pools
        if ( .not. use_matrixcn )then
           ! NOTE: The non matrix version of this is in CNCStateUpdate3::CStateUpdate3 EBK (11/26/2019)
           !                                        and CNNStateUpdate3::NStateUpdate3
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

        else
           m_leafc_to_litter_fire(p)              = leafc(p) * matrix_update_fic(p,ileaf_to_iout_fic, &
                                                  f * (1._r8 - cc_leaf(patch%itype(p)))     * fm_leaf(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_leafc_storage_to_litter_fire(p)      = leafc_storage(p) * matrix_update_fic(p,ileafst_to_iout_fic, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_leafc_xfer_to_litter_fire(p)         = leafc_xfer(p) * matrix_update_fic(p,ileafxf_to_iout_fic, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_livestemc_to_litter_fire(p)          = livestemc(p) * matrix_update_fic(p,ilivestem_to_iout_fic, &
                                                  f * (1._r8 - cc_lstem(patch%itype(p)))    * fm_droot(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_livestemc_storage_to_litter_fire(p)  = livestemc_storage(p) * matrix_update_fic(p,ilivestemst_to_iout_fic, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_livestemc_xfer_to_litter_fire(p)     = livestemc_xfer(p) * matrix_update_fic(p,ilivestemxf_to_iout_fic, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_livestemc_to_deadstemc_fire(p)       = livestemc(p) * matrix_update_fic(p,ilivestem_to_ideadstem_fic,&
                                                  f * (1._r8 - cc_lstem(patch%itype(p)))    * (fm_lstem(patch%itype(p))-fm_droot(patch%itype(p))),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_deadstemc_to_litter_fire(p)          = deadstemc(p) * matrix_update_fic(p,ideadstem_to_iout_fic, &
                                                  f * (1._r8 - cc_dstem(patch%itype(p)))    * fm_droot(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_deadstemc_storage_to_litter_fire(p)  = deadstemc_storage(p) * matrix_update_fic(p,ideadstemst_to_iout_fic, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_deadstemc_xfer_to_litter_fire(p)     = deadstemc_xfer(p) * matrix_update_fic(p,ideadstemxf_to_iout_fic, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_frootc_to_litter_fire(p)             = frootc(p) * matrix_update_fic(p,ifroot_to_iout_fic, &
                                                  f * fm_root(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_frootc_storage_to_litter_fire(p)     = frootc_storage(p) * matrix_update_fic(p,ifrootst_to_iout_fic, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_frootc_xfer_to_litter_fire(p)        = frootc_xfer(p) * matrix_update_fic(p,ifrootxf_to_iout_fic, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_livecrootc_to_litter_fire(p)         = livecrootc(p) * matrix_update_fic(p,ilivecroot_to_iout_fic, &
                                                  f * fm_droot(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_livecrootc_storage_to_litter_fire(p) = livecrootc_storage(p) * matrix_update_fic(p,ilivecrootst_to_iout_fic, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_livecrootc_xfer_to_litter_fire(p)    = livecrootc_xfer(p) * matrix_update_fic(p,ilivecrootxf_to_iout_fic, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_livecrootc_to_deadcrootc_fire(p)     = livecrootc(p) * matrix_update_fic(p,ilivecroot_to_ideadcroot_fic,&
                                                  f * (fm_lroot(patch%itype(p))-fm_droot(patch%itype(p))),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_deadcrootc_to_litter_fire(p)         = deadcrootc(p) * matrix_update_fic(p,ideadcroot_to_iout_fic, &
                                                  f * m * fm_droot(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_deadcrootc_storage_to_litter_fire(p) = deadcrootc_storage(p) * matrix_update_fic(p,ideadcrootst_to_iout_fic, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)
           m_deadcrootc_xfer_to_litter_fire(p)    = deadcrootc_xfer(p) * matrix_update_fic(p,ideadcrootxf_to_iout_fic, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_carbonflux_inst,.True.,.True.)

           m_leafn_to_litter_fire(p)              = leafn(p) * matrix_update_fin(p,ileaf_to_iout_fin, &
                                                  f * (1._r8 - cc_leaf(patch%itype(p)))     * fm_leaf(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_leafn_storage_to_litter_fire(p)      = leafn_storage(p) * matrix_update_fin(p,ileafst_to_iout_fin, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_leafn_xfer_to_litter_fire(p)         = leafn_xfer(p) * matrix_update_fin(p,ileafxf_to_iout_fin, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_livestemn_to_litter_fire(p)          = livestemn(p) * matrix_update_fin(p,ilivestem_to_iout_fin, &
                                                  f * (1._r8 - cc_lstem(patch%itype(p)))    * fm_droot(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_livestemn_storage_to_litter_fire(p)  = livestemn_storage(p) * matrix_update_fin(p,ilivestemst_to_iout_fin, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_livestemn_xfer_to_litter_fire(p)     = livestemn_xfer(p) * matrix_update_fin(p,ilivestemxf_to_iout_fin, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_livestemn_to_deadstemn_fire(p)       = livestemn(p) * matrix_update_fin(p,ilivestem_to_ideadstem_fin,&
                                                  f * (1._r8 - cc_lstem(patch%itype(p)))    * (fm_lstem(patch%itype(p))-fm_droot(patch%itype(p))),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_deadstemn_to_litter_fire(p)          = deadstemn(p) * matrix_update_fin(p,ideadstem_to_iout_fin, &
                                                  f * (1._r8 - cc_dstem(patch%itype(p)))    * fm_droot(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_deadstemn_storage_to_litter_fire(p)  = deadstemn_storage(p) * matrix_update_fin(p,ideadstemst_to_iout_fin, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_deadstemn_xfer_to_litter_fire(p)     = deadstemn_xfer(p) * matrix_update_fin(p,ideadstemxf_to_iout_fin, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_frootn_to_litter_fire(p)             = frootn(p) * matrix_update_fin(p,ifroot_to_iout_fin, &
                                                  f * fm_root(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_frootn_storage_to_litter_fire(p)     = frootn_storage(p) * matrix_update_fin(p,ifrootst_to_iout_fin, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_frootn_xfer_to_litter_fire(p)        = frootn_xfer(p) * matrix_update_fin(p,ifrootxf_to_iout_fin, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_livecrootn_to_litter_fire(p)         = livecrootn(p) * matrix_update_fin(p,ilivecroot_to_iout_fin, &
                                                  f * fm_droot(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_livecrootn_storage_to_litter_fire(p) = livecrootn_storage(p) * matrix_update_fin(p,ilivecrootst_to_iout_fin, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_livecrootn_xfer_to_litter_fire(p)    = livecrootn_xfer(p) * matrix_update_fin(p,ilivecrootxf_to_iout_fin, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_livecrootn_to_deadcrootn_fire(p)     = livecrootn(p) * matrix_update_fin(p,ilivecroot_to_ideadcroot_fin,&
                                                  f * (fm_lroot(patch%itype(p))-fm_droot(patch%itype(p))),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_deadcrootn_to_litter_fire(p)         = deadcrootn(p) * matrix_update_fin(p,ideadcroot_to_iout_fin, &
                                                  f * m * fm_droot(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_deadcrootn_storage_to_litter_fire(p) = deadcrootn_storage(p) * matrix_update_fin(p,ideadcrootst_to_iout_fin, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
           m_deadcrootn_xfer_to_litter_fire(p)    = deadcrootn_xfer(p) * matrix_update_fin(p,ideadcrootxf_to_iout_fin, &
                                                  f * (1._r8 - cc_other(patch%itype(p)))    * fm_other(patch%itype(p)),dt,cnveg_nitrogenflux_inst,.True.,.True.)
        end if

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
        do fp = 1, num_soilp
           p = filter_soilp(fp)
           c = patch%column(p)

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


           m_c_to_litr_fire(c,j,i_met_lit) = &
                m_c_to_litr_fire(c,j,i_met_lit) + &
                ((m_leafc_to_litter_fire(p) * lf_f(patch%itype(p),i_met_lit) &
                +m_leafc_storage_to_litter_fire(p) + &
                m_leafc_xfer_to_litter_fire(p) + &
                m_gresp_storage_to_litter_fire(p) &
                +m_gresp_xfer_to_litter_fire(p))*leaf_prof(p,j) + &
                (m_frootc_to_litter_fire(p) * fr_f(patch%itype(p),i_met_lit) &
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
           ! Here metabolic litter is treated differently than other
           ! types of litter, so it remains outside this litter loop,
           ! in the line above
           do i = i_met_lit+1, i_litr_max
              m_c_to_litr_fire(c,j,i) = m_c_to_litr_fire(c,j,i) + &
                 (m_leafc_to_litter_fire(p) * lf_f(patch%itype(p),i) * leaf_prof(p,j) + &
                 m_frootc_to_litter_fire(p) * fr_f(patch%itype(p),i) * froot_prof(p,j)) * patch%wtcol(p) 
           end do

           m_n_to_litr_fire(c,j,i_met_lit) = &
              m_n_to_litr_fire(c,j,i_met_lit) + &
              ((m_leafn_to_litter_fire(p) * lf_f(patch%itype(p),i_met_lit) + &
                m_leafn_storage_to_litter_fire(p) + &
                m_leafn_xfer_to_litter_fire(p) + &
                m_retransn_to_litter_fire(p)) * leaf_prof(p,j) + &
               (m_frootn_to_litter_fire(p) * fr_f(patch%itype(p),i_met_lit) + &
                m_frootn_storage_to_litter_fire(p) + &
                m_frootn_xfer_to_litter_fire(p)) * froot_prof(p,j) + &
               (m_livestemn_storage_to_litter_fire(p) + &
                m_livestemn_xfer_to_litter_fire(p) + &
                m_deadstemn_storage_to_litter_fire(p) + &
                m_deadstemn_xfer_to_litter_fire(p)) * stem_prof(p,j) + &
               (m_livecrootn_storage_to_litter_fire(p) + &
                m_livecrootn_xfer_to_litter_fire(p) + &
                m_deadcrootn_storage_to_litter_fire(p) + &
                m_deadcrootn_xfer_to_litter_fire(p)) * croot_prof(p,j)) * patch%wtcol(p)
           ! Here metabolic litter is treated differently than other
           ! types of litter, so it remains outside this litter loop,
           ! in the line above
           do i = i_met_lit+1, i_litr_max
              m_n_to_litr_fire(c,j,i) = &
                 m_n_to_litr_fire(c,j,i) + &
                 (m_leafn_to_litter_fire(p) * lf_f(patch%itype(p),i) * leaf_prof(p,j) + &
                  m_frootn_to_litter_fire(p) * fr_f(patch%itype(p),i) * froot_prof(p,j)) * patch%wtcol(p)
           end do
        end do
     end do
     !
     ! vertically-resolved decomposing C/N fire loss   
     ! column loop
     !
     num_actfirec = 0
     do fc = 1,num_soilc
        c = filter_soilc(fc)

        f = farea_burned(c) 

        if(f /= 0 .or. f /= baf_crop(c))then
           num_actfirec = num_actfirec + 1
           filter_actfirec(num_actfirec) = c
        end if
        do j = 1, nlevdecomp
           ! carbon fluxes
           do l = 1, ndecomp_pools
              if ( is_litter(l) ) then
                 m_decomp_cpools_to_fire_vr(c,j,l) = decomp_cpools_vr(c,j,l) * f * &
                      cmb_cmplt_fact_litter
                 if(use_soil_matrixcn)then! matrix is the same for C and N in the fire.
                    associate( &
                       matrix_decomp_fire_k  => soilbiogeochem_carbonflux_inst%matrix_decomp_fire_k_col & ! Output: [real(r8) (:,:)   ]  (gC/m3/step) VR deomp. C fire loss in matrix representation
                    )
                    matrix_decomp_fire_k(c,j+nlevdecomp*(l-1)) = matrix_decomp_fire_k(c,j+nlevdecomp*(l-1)) &
                     - f * cmb_cmplt_fact_litter * dt
                    end associate
                 end if
              end if
              if ( is_cwd(l) ) then
                 m_decomp_cpools_to_fire_vr(c,j,l) = decomp_cpools_vr(c,j,l) * &
                      (f-baf_crop(c)) * cmb_cmplt_fact_cwd
                 if(use_soil_matrixcn)then
                    associate( &
                       matrix_decomp_fire_k  => soilbiogeochem_carbonflux_inst%matrix_decomp_fire_k_col & ! Output: [real(r8) (:,:)   ]  (gC/m3/step) VR deomp. C fire loss in matrix representation
                    )
                    matrix_decomp_fire_k(c,j+nlevdecomp*(l-1)) = matrix_decomp_fire_k(c,j+nlevdecomp*(l-1)) &
                     - (f-baf_crop(c)) * cmb_cmplt_fact_cwd * dt
                    end associate
                 end if
              end if
           end do

           ! nitrogen fluxes
           do l = 1, ndecomp_pools
              if ( is_litter(l) ) then
                 m_decomp_npools_to_fire_vr(c,j,l) = decomp_npools_vr(c,j,l) * f * &
                      cmb_cmplt_fact_litter
              end if
              if ( is_cwd(l) ) then
                 m_decomp_npools_to_fire_vr(c,j,l) = decomp_npools_vr(c,j,l) * &
                      (f-baf_crop(c)) * cmb_cmplt_fact_cwd
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
                      baf_peatf(c))/2.0_r8*dt))/(dtrotr_col(c)*dayspyr*secspday/dt)/dt
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
  subroutine CNFireReadParams( this, ncid )
    !
    ! Read in the constant parameters from the input NetCDF parameter file
    ! !USES:
    use ncdio_pio   , only: file_desc_t
    use paramUtilMod, only: readNcdioScalar
    !
    ! !ARGUMENTS:
    implicit none
    class(cnfire_base_type)         :: this
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'CNFireReadParams'
    !--------------------------------------------------------------------

    ! Factor related to dependence of fuel combustibility on 30-day running mean of relative humidity (unitless)
    call readNcdioScalar(ncid, 'prh30', subname, cnfire_params%prh30)
    ! Ignition efficiency of cloud-to-ground lightning (unitless)
    call readNcdioScalar(ncid, 'ignition_efficiency', subname, cnfire_params%ignition_efficiency)

  end subroutine CNFireReadParams

end module CNFireBaseMod
